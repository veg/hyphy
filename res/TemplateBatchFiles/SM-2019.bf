LoadFunctionLibrary ("libv3/tasks/trees.bf");
LoadFunctionLibrary ("libv3/tasks/alignments.bf");
LoadFunctionLibrary ("libv3/convenience/regexp.bf");
LoadFunctionLibrary ("libv3/convenience/math.bf");
LoadFunctionLibrary ("libv3/IOFunctions.bf");
LoadFunctionLibrary ("libv3/UtilityFunctions.bf");

//#profile START;

sm.analysis_description = {terms.io.info : "This analysis implements canonical and modified versions of the Slatkin-Maddison
phylogeny based test for population segregation. The test estimates the minimum number of migration events using maximum
parsimony, and then evaluating whether or not this number is lower than expected in a panmictic or unstructured population using
permutation tests",
                           terms.io.version : "0.1",
                           terms.io.reference : "*A cladistic measure of gene flow inferred from the phylogenies of alleles* (1989), Genetics 123(3):603-613",
                           terms.io.authors : "Sergei L Kosakovsky Pond",
                           terms.io.contact : "spond@temple.edu",
                           terms.io.requirements : "a phylogenetic tree with leaf names that can be partitioned into sets using regular expressions"
                          };

io.DisplayAnalysisBanner (sm.analysis_description);

KeywordArgument ("tree", "The Newick tree string defining the topology to use for testing");
KeywordArgument ("groups", "The number of compartments to test", "2");


//KeywordArgument ("output", "Write the output JSON to ", "2");



utility.SetEnvVariable("IS_TREE_PRESENT_IN_DATA", FALSE); // for re-entrance
SetDialogPrompt ("Load the phylogeny to test for compartmentalization");
sm.tree = trees.LoadAnnotatedTopology (TRUE);

sm.leaves = utility.Filter (sm.tree[terms.trees.partitioned], "_v_", "_v_==terms.tree_attributes.leaf");
sm.leaf_count = utility.Array1D (sm.leaves);

io.ReportProgressMessageMD ("SM", "data", "Using a tree with **" +
                              sm.leaf_count +
                              "** leaves" +
                              "\n\n\`" + Join (",", utility.Keys (sm.leaves)) + "\`\n");

sm.group_count = io.PromptUser(">How many branch classes are there", 2, 2, sm.leaf_count, TRUE);
sm.regular_expressions = {};

for (sm.i = 0; sm.i < sm.group_count; sm.i += 1) {

    KeywordArgument ("description-" + (sm.i + 1), "Description for sequences in compartment " + (sm.i+1));
    sm.tag = io.PromptUserForString (">Please provide a description for group of branches " + (sm.i+1));
    KeywordArgument ("regexp-" + (sm.i + 1), "Regular expression to select the branches in compartment _" + sm.tag + "_");
    sm.regexp = io.PromptUserForString (">Please provide a regular expression to select the _" + sm.tag + "_ branches");
    sm.regular_expressions [sm.tag] = sm.regexp;
}

sm.regexp_to_name = utility.SwapKeysAndValues (sm.regular_expressions);

assert (Abs (sm.regular_expressions) == Abs (sm.regexp_to_name), "All regular expressions must be distinct");

sm.partitions = regexp.PartitionByRegularExpressions (utility.Keys (sm.leaves), sm.regular_expressions);
sm.partitions_with_labels = {};

sm.partition_counts = {};

utility.ForEachPair (sm.partitions , "_key_", "_value_", '
    if (Abs (_key_)) {
        assert (Abs (_value_) > 0, "Empty branch set for _" + sm.regexp_to_name[_key_] + "_");
        io.ReportProgressMessageMD ("SM", "data", "* Branch set _" + sm.regexp_to_name[_key_] + "_ has " + Abs (_value_) + " branches.");
        sm.partition_counts[(sm.regexp_to_name[_key_])] = Abs(_value_);
        sm.partitions_with_labels[(sm.regexp_to_name[_key_])] = _value_;
    }
');

KeywordArgument ("replicates", "The number of bootstrap replicates", "1000");

sm.replicates = io.PromptUser(">How many bootstrap replicates", 1000, 1, 1000000, TRUE);

sm.bootstrap = trees.BootstrapSupport (sm.tree);

KeywordArgument ("weight", "Probability of branch selection for structured permutation [0-1]; 0 = classical Slatkin-Maddison, 1 = fully structured", "0.2");
sm.default_weighting = io.PromptUser(">Probability of branch selection ", 0.2, 0.0,1.0, FALSE);

sm.bootstrap_weighting = {"default" : sm.default_weighting};


if (utility.Array1D (sm.bootstrap )) {
 
   KeywordArgument ("use-bootstrap", "Use bootstrap weights to respect well supported clades", "Yes");

   if (io.SelectAnOption ({"No" : "Each internal branch has the same probability of being randomized in structured permutations",
                                "Yes"  : "For branches with bootstrap support, determine the probability of randomization based on the bootstrap (higher -> more likely)"},
                                "Use bootstrap weighting"
                                ) == "Yes") {

        sm.bootstrap_weighting * utility.Map (sm.bootstrap, "_bs_", "sm.bootstrap_weight(_bs_)");
   }

}

sm.tree.string = sm.tree[terms.trees.newick_with_lengths];
sm.node_labels = {};



Topology T = sm.tree.string;

utility.ForEachPair (sm.partitions, "_regexp_", "_leaves_",
	'
		tag = sm.regexp_to_name[_regexp_];
		if (Abs (tag)) {
			utility.ForEach (_leaves_, "_leaf_", "sm.node_labels[_leaf_]=tag");
		}
	 '
);



sm.mp = Max (T, {"labels": sm.node_labels});
sm.node_scores = sm.mp ["node-scores"];
sm.score = sm.mp ["score"];
sm.node_scores_replicates = utility.Map (sm.node_scores, "_value_", "{sm.replicates,1}");
sm.node_scores_standard_replicates = utility.Map (sm.node_scores, "_value_", "{sm.replicates,1}");

io.ReportProgressMessageMD('SM',  'result', 'Inferred **' +  sm.score + '** migration events');

sm.resampled_distribution = {sm.replicates , 2};
sm.resampled_p_value = 0;
sm.resampled_sum = 0;


for (sm.k = 0; sm.k < sm.replicates ; ) {
    sm.reshuffled_tree = "" + Random (T, sm.bootstrap_weighting);
    Topology SimT = sm.reshuffled_tree;

    sm.replicate_scores = Max (SimT, {"labels": sm.node_labels});
    sm.resampled_distribution[sm.k][1] = sm.replicate_scores ["score"];
    sm.replicate_scores_standard = Max (T, {"labels": Random (sm.node_labels,0)});
    sm.resampled_distribution[sm.k][0] = sm.replicate_scores_standard["score"];


    utility.ForEachPair (sm.replicate_scores_standard["node-scores"], "_node_", "_value_", '
         (sm.node_scores_standard_replicates [_node_])[sm.k] = _value_;
         (sm.node_scores_replicates [_node_])[sm.k]         = ((sm.replicate_scores["node-scores"])[_node_]);
     '
    );

    if ( sm.resampled_distribution[sm.k][0] <= sm.score) {
        sm.resampled_p_value += 1;
    }

    sm.resampled_sum += sm.resampled_distribution[sm.k][0];
    sm.k+=1;

    io.ReportProgressBar ("SM", "Replicate " +  (sm.k) + "; mean # events " + Format (sm.resampled_sum/sm.k, 6, 3) + "; running p-value <= "
         + Format ( sm.resampled_p_value/(sm.k+1), 6, 3));


}

io.ClearProgressBar();

sm.dist.shuffled     = sm.resampled_distribution[-1][0];
sm.dist.structured   = sm.resampled_distribution[-1][1];

sm.standard_node.p   = utility.MapWithKey (sm.node_scores_standard_replicates, "_key_", "_value_",
'(+(_value_["_MATRIX_ELEMENT_VALUE_<=sm.node_scores[_key_]"])+1)/(sm.replicates+1)'
);


function sm.node_cutoff_p (node, scores) {
    sm.sampled_stats = sm.node_scores_replicates[node] % 0;
    sm.dist.structured.cutoff  = sm.sampled_stats [sm.replicates * 9 $ 10];
    return (+(scores["_MATRIX_ELEMENT_VALUE_<=sm.dist.structured.cutoff"])+1)/(sm.replicates+1);
}

sm.structured_node.p   = utility.MapWithKey (sm.node_scores_standard_replicates, "_key_", "_value_",
    'sm.node_cutoff_p(_key_,_value_)'
);




sm.dist.shuffled     = sm.dist.shuffled  % 0;
sm.dist.structured   = sm.dist.structured  % 0;

sm.dist.structured.cutoff = sm.dist.structured [sm.replicates * 9 $ 10];

sm.structured.p = +sm.dist.shuffled["_MATRIX_ELEMENT_VALUE_<=sm.dist.structured.cutoff"];

sm.sampled_stats = math.GatherDescriptiveStats (sm.resampled_distribution[-1][0]);
sm.sampled_stuctured_stats = math.GatherDescriptiveStats (sm.resampled_distribution[-1][1]);


io.ReportProgressMessageMD('SM',  'result', 'Based on **' + sm.replicates + '** leaf label permutations, the standard (full panmixia) p-value for compartmentalization test was < '
         + Format ( (1+sm.resampled_p_value)/(sm.replicates+1), 6, 3));
io.ReportProgressMessageMD('SM',  'distribution', 'The null distribution of migration events had mean _'
            + Format (sm.sampled_stats[terms.math.mean], 6, 3)
            + '_, median _' + Format (sm.sampled_stats[terms.math.median], 6, 3) +
            '_, and 2.5% - 97.5% range of ' + Format (sm.sampled_stats[terms.math._2.5], 6, 3) +
            '-' + Format (sm.sampled_stats[terms.math._97.5], 6, 3)

);

io.ReportProgressMessageMD('SM',  'reshuffled', 'Based on **' + sm.replicates + '** _structured_ permutations, the p-value for compartmentalization test was < '
         + Format ( (1+sm.structured.p) /(sm.replicates+1), 6, 3));
io.ReportProgressMessageMD('SM',  'reshuffled', '\n> This p-value is derived by comparing the distribution of migration events from the panmictic reshuffle to the 90% percentile of the simulated distribution of expected migrations if leaf labels are permuted partially respecting subtree structure (block permutations), which results in **' + sm.dist.structured.cutoff + '** expected migrations');


sm.branch_count = {2,1};
sm.table_screen_output  = {{"Branch", "Migrations", "Panmictic P", "Structured P"}};

sm.table_output_options = {
                            terms.table_options.header : TRUE,
                            terms.table_options.minimum_column_width: 16,
                            terms.table_options.align : "center"
                           };

utility.ForEachPair (sm.structured_node.p, "_node_", "_pvalue_",
'
    if (_pvalue_ <= 0.05 || sm.standard_node.p[_node_] <= 0.05) {

        if (Max (sm.branch_count,0) == 0) {
            io.ReportProgressMessageMD("SM", "branches", "List of individual branches contributing to compartmentalization signal");
            fprintf (stdout,
                io.FormatTableRow (sm.table_screen_output,sm.table_output_options));
            sm.table_output_options[terms.table_options.header] = FALSE;
        }
        sm.row = {4,1};
        sm.row[0] = _node_;
        sm.row[1] = sm.node_scores[_node_];
        if (sm.standard_node.p[_node_] < 0.05) {
            sm.row[2] = Format (sm.standard_node.p[_node_], 8, 4);
        } else {
            sm.row [2] = "---";
        }
        if (_pvalue_ < 0.05) {
            sm.row[3] = Format (_pvalue_, 8, 4);
        } else {
            sm.row [3] = "---";
        }

        fprintf (stdout, io.FormatTableRow (sm.row,sm.table_output_options));


        sm.branch_count[0] += (_pvalue_ <= 0.05);
        sm.branch_count[1] += (sm.standard_node.p[_node_] <= 0.05);
    }
'
);

function sm.label_with_standard_p_values (node) {
    if (sm.standard_node.p / node) {
        return "" + sm.standard_node.p[node];
    }
    return node;
}

function sm.label_with_p_values (node) {
    if (sm.standard_node.p / node) {
        return "" +sm.structured_node.p[node];
    }
    return node;
}


sm.json = {
    'partitions' : sm.partitions_with_labels,
    'tree' : sm.tree,
    'tree-p-structured' : tree.Annotate ("T", "sm.label_with_p_values", "[]", TRUE),
    'tree-p-panmictic' : tree.Annotate ("T", "sm.label_with_standard_p_values", "[]", TRUE),
    'leaf-count' : sm.leaf_count,
    'partition-counts' : sm.partition_counts,
    'replicates' : sm.replicates,
    'compartments' : sm.group_count,
    'migrations' : sm.score,
    'node-migrations' : sm.node_scores,
    'structured-cutoff': sm.dist.structured.cutoff,
    'node-p-value' : {
        'panmictic'  : sm.standard_node.p,
        'structured' : sm.structured_node.p
    },
    'p-value' : {
        'panmictic'  :  (1+sm.resampled_p_value)/(sm.replicates+1),
        'structured' :  (1+sm.structured.p) /(sm.replicates+1)
    },
    'simulations' : {
        'panmictic' : sm.dist.shuffled,
        'structured' : sm.dist.structured
    },
    "events" : sm.mp["substitutions"]
};


sm.json_file = sm.tree[terms.data.file] + "_SM.json";

KeywordArgument ("output", "Write the JSON file here (default is to save to the same path as the tree file + '_SM.json')", sm.json_file);
sm.json_file = io.PromptUserForFilePath ("Save the resulting JSON file to");
io.SpoolJSON (sm.json, sm.json_file);
io.ReportProgressMessageMD('SM',  'file', 'Saving detailed report as a JSON file to \`' + sm.json_file + '\`');

return sm.json;

//------

lfunction sm.bootstrap_weight (bs) {
    if (bs > 1) {
        bs = bs / 100;
    }
    
    if (bs < 0.7) {
        return ^"sm.default_weighting";
    }
    return Max (^"sm.default_weighting",bs*bs);
}
