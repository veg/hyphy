LoadFunctionLibrary ("libv3/tasks/trees.bf");
LoadFunctionLibrary ("libv3/tasks/alignments.bf");
LoadFunctionLibrary ("libv3/convenience/regexp.bf");
LoadFunctionLibrary ("libv3/convenience/math.bf");
LoadFunctionLibrary ("libv3/IOFunctions.bf");
LoadFunctionLibrary ("libv3/UtilityFunctions.bf");

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
    sm.tag = io.PromptUserForString (">Please provide a description for group of branches " + (sm.i+1));
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


sm.replicates = io.PromptUser(">How many bootstrap replicates", 100, 1, 1000000, TRUE);

sm.bootstrap = trees.BootstrapSupport (sm.tree);
sm.bootstrap_weighting = {"default" : 0.2};

if (utility.Array1D (sm.bootstrap )) {
   if (io.SelectAnOption ({"No" : "Each internal branch has the same probability of being randomized in structured permutations", 
                                "Yes"  : "For branches with bootstrap support, determine the probability of randomization based on the bootstrap (higher -> more likely)"},
                                "Use bootstrap weighting"
                                ) == "Yes") {
                                
        sm.bootstrap_weighting * utility.Map (sm.bootstrap, "_bs_", "Max(0.2,(_bs_>1)*(_bs_/100)^2+(_bs_<=1)*_bs_^2)");
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
sm.score = sm.mp ["score"];

io.ReportProgressMessageMD('SM',  'result', 'Inferred **' +  sm.score + '** migration events');

sm.resampled_distribution = {sm.replicates , 2};
sm.resampled_p_value = 0;
sm.resampled_sum = 0;

sm.shuffling_probability = Max (0.2, 10 / (BranchCount (T) + TipCount (T)));

for (sm.k = 0; sm.k < sm.replicates ; ) {
    //if (sm.method == "Restricted") {
        sm.reshuffled_tree = "" + Random (T, sm.bootstrap_weighting);
        Topology SimT = sm.reshuffled_tree;
        sm.resampled_distribution[sm.k][1] = Max (SimT, {"labels": sm.node_labels})["score"];
        //(trees.ParsimonyLabel ("SimT", sm.node_labels))["score"];
    //} else {    
        sm.resampled_distribution[sm.k][0] =  Max (SimT, {"labels": Random (sm.node_labels,0)})["score"];
        //(trees.ParsimonyLabel ("T", Random (sm.node_labels,0)))["score"];
    //}
    
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

sm.dist.shuffled     = sm.dist.shuffled  % 0;
sm.dist.structured   = sm.dist.structured  % 0;

sm.dist.structured.cutoff = sm.dist.structured [sm.replicates * 90 $ 100];


sm.structured.p = +sm.dist.shuffled["_MATRIX_ELEMENT_VALUE_<= sm.dist.structured.cutoff"];

sm.sampled_stats = math.GatherDescriptiveStats (sm.resampled_distribution[-1][0]);
sm.sampled_stuctured_stats = math.GatherDescriptiveStats (sm.resampled_distribution[-1][1]);


io.ReportProgressMessageMD('SM',  'result', 'Based on **' + sm.replicates + '** leaf label permutations, the standard (full panmixia) p-value for compartmentalization test was < '
         + Format ( sm.resampled_p_value/(sm.replicates+1), 6, 3));
io.ReportProgressMessageMD('SM',  'distribution', 'The null distribution of migration events had mean _' 
            + Format (sm.sampled_stats[terms.math.mean], 6, 3) 
            + '_, median _' + Format (sm.sampled_stats[terms.math.median], 6, 3) +
            '_, and 2.5% - 97.5% range of ' + Format (sm.sampled_stats[terms.math._2.5], 6, 3) +
            '-' + Format (sm.sampled_stats[terms.math._97.5], 6, 3) 
            
);

io.ReportProgressMessageMD('SM',  'reshuffled', 'Based on **' + sm.replicates + '** _structured_ permutations, the p-value for compartmentalization test was < '
         + Format ( sm.structured.p /(sm.replicates+1), 6, 3));
io.ReportProgressMessageMD('SM',  'reshuffled', '\n> This p-value is derived by comparing the distribution of migration events from the panmictic reshuffle to the 90% percentile of the simulated distribution of expected migrations if leaf labels are permuted partially respecting subtree structure (block permutations), which results in **' + sm.dist.structured.cutoff + '** expected migrations');

sm.json_path = sm.tree[terms.data.file] + ".json";
io.ReportProgressMessageMD('SM',  'file', 'Saving detailed report as a JSON file to \`' + sm.json_path + '\`');

sm.json = {
    'partitions' : sm.partitions_with_labels,
    'tree' : sm.tree,
    'leaf-count' : sm.leaf_count,
    'partition-counts' : sm.partition_counts,
    'replicates' : sm.replicates,
    'compartments' : sm.group_count, 
    'migrations' : sm.score,
    'structured-cutoff': sm.dist.structured.cutoff,
    'p-value' : {
        'panmictic'  :  sm.resampled_p_value/(sm.replicates+1), 
        'structured' :  sm.structured.p /(sm.replicates+1)
    },
    'simulations' : {
        'panmictic' : sm.dist.shuffled,
        'structured' : sm.dist.structured   
    },
    "events" : sm.mp["substitutions"]
};

io.SpoolJSON (sm.json, sm.json_path);
return sm.json;
