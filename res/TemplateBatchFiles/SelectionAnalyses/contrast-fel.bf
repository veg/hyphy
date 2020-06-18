RequireVersion("2.5.1");

/*------------------------------------------------------------------------------
    Load library files
*/

LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/stats.bf");
LoadFunctionLibrary("libv3/all-terms.bf");

LoadFunctionLibrary("libv3/tasks/ancestral.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/tasks/mpi.bf");

LoadFunctionLibrary("libv3/models/codon/MG_REV.bf");

LoadFunctionLibrary("libv3/convenience/math.bf");

LoadFunctionLibrary("modules/io_functions.ibf");

//utility.ToggleEnvVariable ("OPTIMIZATION_PRECISION", 1);
//utility.ToggleEnvVariable ("OPTIMIZATION_TIME_HARD_LIMIT", 1);


/*------------------------------------------------------------------------------ Display analysis information
*/

fel.analysis_description = {
    terms.io.info: "Contrast-FEL (Fixed Effects Likelihood) investigates whether or not selective pressures differ between two or more sets of
    branches at a site. Site-specific synonymous (alpha) and non-synonymous (beta, one per branch set) substitution rates are estimated
    and then beta rates are tested for equality at each site. LRT and permutation tests ar used to assess significance at each site, and FDR is applied alignment wide to call sites with different selective profiles",
    terms.io.version: "0.5",
    terms.io.reference: "TBD",
    terms.io.authors: "Sergei L Kosakovsky Pond and Steven Weaver",
    terms.io.contact: "spond@temple.edu",
    terms.io.requirements: "in-frame codon alignment and a phylogenetic tree; only single partition data are supported",
    terms.io.help : "http://www.hyphy.org/methods/other/contrast-fel/"
};

io.DisplayAnalysisBanner(fel.analysis_description);


/*------------------------------------------------------------------------------
    Environment setup
*/

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);

/*------------------------------------------------------------------------------
    Globals
*/



terms.fel.pairwise      = "pairwise";
terms.fel.permutation   = "permutation";
terms.fel.test_keys     = "test keys";
fel.site_alpha          = "Site relative synonymous rate";
fel.site_beta_reference = "Site relative non-synonymous rate (reference branches)";
fel.site_tested_classes = {};
fel.alpha.scaler        = "fel.alpha_scaler";
fel.site.permutations   = 20;


// default cutoff for printing to screen
fel.p_value = 0.1;
fel.q_value = 0.2;
fel.scaler_prefix = "FEL.scaler";


KeywordArgument ("code",      "Which genetic code should be used", "Universal");
KeywordArgument ("alignment", "An in-frame codon alignment in one of the formats supported by HyPhy");
KeywordArgument ("tree",      "A phylogenetic tree (optionally annotated with {})", null, "Please select a tree file for the data:");
KeywordArgument ("branch-set", "The set of branches to use for testing");


// The dictionary of results to be written to JSON at the end of the run
fel.json = {
    terms.json.analysis: fel.analysis_description,
    terms.json.input: {},
    terms.json.fits: {},
    terms.json.timers: {},
};

fel.display_orders =   {terms.original_name: -1,
                        terms.json.nucleotide_gtr: 0,
                        terms.json.global_mg94xrev: 1
                       };



selection.io.startTimer (fel.json [terms.json.timers], "Total time", 0);



namespace fel {
    LoadFunctionLibrary ("modules/shared-load-file.bf");
    load_file ({utility.getGlobalValue("terms.prefix"): "fel", 
                utility.getGlobalValue("terms.settings") : {utility.getGlobalValue("terms.settings.branch_selector") : "fel.select_branches"}});
}

KeywordArgument ("srv", "Include synonymous rate variation in the model", "Yes");
KeywordArgument ("permutations", "Perform permutation significance tests", "Yes");
KeywordArgument ("p-value", "Significance value for site-tests", "0.05");
KeywordArgument ("q-value", "Significance value for FDR reporting", "0.20");


/* Prompt for one-rate or two-rate analysis */
fel.srv = io.SelectAnOption( {{"Yes", "[Recommended] Consider synonymous rate variation (dS varies across sites)."}, {"No", "Ignore synonymous rate variation (dS := 1 at each site)."}},
                                  "Use synonymous rate variation? Strongly recommended YES for selection inference.");

if (fel.srv == "Yes"){
    fel.srv = TRUE
} else {
    fel.srv = FALSE
}

/* Prompt for one-rate or two-rate analysis */
fel.permutations = io.SelectAnOption( {{"Yes", "For sites with significant p-values, perform an additional permutation test (over branch assignments) to assess significance. Adds computational cost, reduces false positves"}, 
                                       {"No", "Do not perform additional tests (faster, higher risk of false positves)"}},
                                  "Perform permutation significance tests");

if (fel.permutations == "Yes"){
    fel.permutations = TRUE
} else {
    fel.permutations = FALSE
}

/* Prompt for p value threshold */
fel.p_value  = io.PromptUser ("\n>Select nominal p-value threshold to use when testing for selection (FDR correction will be performed on all sites)",0.1,0,1,FALSE);
fel.q_value  = io.PromptUser ("\n>Select nominal the q-value threshold to use when testing for selection (FDR reporting)",0.2,0,1,FALSE);

KeywordArgument ("output", "Write the resulting JSON to this file (default is to save to the same path as the alignment file + 'FEL.json')", fel.codon_data_info [terms.json.json]);
fel.codon_data_info [terms.json.json] = io.PromptUserForFilePath ("Save the resulting JSON file to");

io.ReportProgressMessageMD('FEL',  'selector', 'Branches to use as the test set in the FEL-contrast analysis');

fel.branch_sets = {};

utility.ForEachPair (fel.selected_branches[0], "_branch_", "_model_",
"
    utility.EnsureKey (fel.branch_sets, _model_);
    fel.branch_sets[_model_] + _branch_;
");

fel.branch_class_count = utility.Array1D (fel.branch_sets);
fel.scaler_parameter_names = {};

io.ReportProgressMessageMD('FEL', 'selector', "Selected `fel.branch_class_count` sets of branches to test\n");

fel.branch_class_counter = 0;
fel.branches.testable = {};
fel.branches.has_background = FALSE;


utility.ForEachPair (fel.branch_sets, "_group_", "_branches_",
    "
     if (_group_ != terms.tree_attributes.background) {
        fel.site_tested_classes [_group_] = 'Site relative non-synonymous rate (' + _group_ + ' branches)';
        fel.branch_class_counter += 1;
        fel.branches.testable + _group_;
        fel.scaler_parameter_names [_group_] = 'fel.beta_scaler_group_' + fel.branch_class_counter;
        io.ReportProgressMessageMD('FEL',  'selector', '* Selected ' + Abs(_branches_) + ' branches in group _' + _group_ + '_ : \\\`' + Join (', ',_branches_) + '\\\`')
     } else {
        fel.scaler_parameter_names [_group_] = 'fel.beta_scaler_background';
        fel.site_tested_classes [_group_] = 'Site relative non-synonymous rate (reference branches)';
        fel.branches.has_background = TRUE;
        io.ReportProgressMessageMD('FEL',  'selector', '* ' + Abs(_branches_) + ' branches are in the background group : \\\`' + Join (', ',_branches_) + '\\\`')
     }
     "
);



fel.test_count = 1;
if (fel.branch_class_counter > 2) {
    fel.test_count += fel.branch_class_counter * (fel.branch_class_counter-1) / 2;
}

fel.table_headers = {4 + utility.Array1D (fel.scaler_parameter_names) + fel.branch_class_counter + fel.test_count, 2};
fel.table_headers [0][0] = "alpha"; fel.table_headers [0][1] = "Synonymous substitution rate at a site";

fel.k = 1;
fel.rate.names = utility.Keys (fel.scaler_parameter_names);
utility.ForEach (fel.rate.names, "_n_", ,
'
    fel.table_headers [fel.k][0] = "beta (" + _n_ + ")";
    fel.table_headers [fel.k][1] = "Non-synonymous substitution rate at a site for " + _n_ + " branches";
    fel.k += 1;
');

fel.report.rate_count = fel.k;

utility.ForEach (fel.branches.testable, "_n_", ,
'
    fel.table_headers [fel.k][0] = "subs (" + _n_ + ")";
    fel.table_headers [fel.k][1] = "Substitutions mapped to " + _n_ + " branches";
    fel.k += 1;
');

fel.report.subs_kinds = utility.Array1D (fel.branches.testable);

fel.table_headers [fel.k][0] = "P-value (overall)";
fel.table_headers [fel.k+1][0] = "Q-value (overall)";

if (fel.branch_class_counter == 1) {
    fel.table_headers [fel.k][1] = "P-value for the test that " + fel.branches.testable[0] + " branches have different non-synonymous rates than the rest of the branches";
    fel.table_headers [fel.k+1][1] = "Q-value for the test that " + fel.branches.testable[0] + " branches have different non-synonymous rates than the rest of the branches";
} else {
    fel.table_headers [fel.k][1] = "P-value for the test that non-synonymous rates differ between any of the selected groups: " + Join (",", fel.branches.testable);
    fel.table_headers [fel.k+1][1] = "Q-value for the test that non-synonymous rates differ between any of the selected groups: " + Join (",", fel.branches.testable);
}

fel.k += 1;

fel.report.test_count = 1 + (fel.branch_class_counter > 2) * (fel.branch_class_counter * (fel.branch_class_counter-1) / 2);
fel.tests.key = {fel.report.test_count, 4};
fel.tests.key [0][0] = "overall";
fel.tests.key [0][1] = -1;
fel.tests.key [0][3] = fel.report.rate_count + fel.report.subs_kinds;

fel.test.index = fel.test_count > 1;
fel.k += 1;


for (fel.v = 0; fel.v < fel.branch_class_counter; fel.v += 1) {
    for (fel.v2 = fel.v + 1; fel.v2 < fel.branch_class_counter; fel.v2 += 1) {
        fel.table_headers [fel.k][0] = "P-value for " + fel.branches.testable[fel.v] + " vs " + fel.branches.testable[fel.v2];
        fel.tests.key [fel.test.index][0] = fel.branches.testable[fel.v] + " vs " + fel.branches.testable[fel.v2];
        fel.tests.key [fel.test.index][1] = 1 + fel.v;
        fel.tests.key [fel.test.index][2] = 1 + fel.v2;

        if (fel.test_count > 1) {
            fel.tests.key [fel.test.index][3] = fel.report.rate_count + fel.test.index + fel.report.subs_kinds + 1;
            fel.table_headers [fel.k][1] = "P-value for  the test that non-synonymous rates differ between " + fel.branches.testable[fel.v] + " and " + fel.branches.testable[fel.v2] + " branches";
            fel.k += 1;
        }
        fel.test.index += 1;
    }
}

io.ReportProgressMessageMD('FEL',  'tests', '**' + fel.test_count + "** " + io.SingularOrPlural (fel.test_count , "test", "tests")+" will be performed at each site");

fel.report.counts         = {fel.report.test_count,1};

fel.table_headers [fel.k][0] = "Permutation p-value";
fel.table_headers [fel.k][1] = "Label permutation test for significant sites";

fel.k += 1;


fel.table_headers [fel.k][0] = "Total branch length";
fel.table_headers [fel.k][1] = "The total length of branches contributing to inference at this site, and used to scale beta-alpha";

/**
This table is meant for HTML rendering in the results web-app; can use HTML characters, the second column
is 'pop-over' explanation of terms. This is ONLY saved to the JSON file. For Markdown screen output see
the next set of variables.
*/

fel.table_screen_output  = {{"Codon", "alpha", "beta", "substitutions", "test", "LRT p-value", "Permutation p-value"}};

fel.table_output_options = {terms.table_options.header : TRUE, terms.table_options.align : "center",
                            terms.table_options.column_widths : { "0" : 8, "1" : 16, "2" : 30, "3" : 30, "4" : 40, "5" : 10, "6" : 10}};

selection.io.startTimer (fel.json [terms.json.timers], "Model fitting",1);
estimators.fixSubsetOfEstimates(fel.gtr_results, fel.gtr_results[terms.global]);

namespace fel {
    doPartitionedMG ("fel", FALSE);
}


//fel.final_partitioned_mg_results = fel.partitioned_mg_results;

io.ReportProgressMessageMD ("fel", "codon-refit", "Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model");


fel.final_partitioned_mg_results = estimators.FitMGREV (fel.filter_names, fel.trees, fel.codon_data_info [terms.code], {
    terms.run_options.model_type: terms.local,
    terms.run_options.partitioned_omega: fel.selected_branches,
    terms.run_options.retain_lf_object: TRUE
}, fel.partitioned_mg_results);



io.ReportProgressMessageMD("fel", "codon-refit", "* Log(L) = " + Format(fel.final_partitioned_mg_results[terms.fit.log_likelihood],8,2));
fel.global_dnds = selection.io.extract_global_MLE_re (fel.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio);
utility.ForEach (fel.global_dnds, "_value_", 'io.ReportProgressMessageMD ("fel", "codon-refit", "* " + _value_[terms.description] + " = " + Format (_value_[terms.fit.MLE],8,4));');




//Store MG94 to JSON
selection.io.json_store_lf_withEFV (fel.json,
                            terms.json.global_mg94xrev,
                            fel.final_partitioned_mg_results[terms.fit.log_likelihood],
                            fel.final_partitioned_mg_results[terms.parameters],
                            fel.sample_size,
                            utility.ArrayToDict (utility.Map (fel.global_dnds, "_value_", "{'key': _value_[terms.description], 'value' : Eval({{_value_ [terms.fit.MLE],1}})}")),
                            (fel.final_partitioned_mg_results[terms.efv_estimate])["VALUEINDEXORDER"][0],
                            fel.display_orders[terms.json.global_mg94xrev]);


estimators.fixSubsetOfEstimates(fel.final_partitioned_mg_results, fel.final_partitioned_mg_results[terms.global]);

utility.ForEachPair (fel.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(fel.json, terms.json.global_mg94xrev, terms.branch_length, fel.display_orders[terms.json.global_mg94xrev],
                                             _key_,
                                             selection.io.extract_branch_info((fel.final_partitioned_mg_results[terms.branch_length])[_key_], "selection.io.branch.length"));');


selection.io.stopTimer (fel.json [terms.json.timers], "Model fitting");

// define the site-level likelihood function

fel.site.mg_rev = model.generic.DefineModel("models.codon.MG_REV.ModelDescription",
        "fel_mg", {
            "0": parameters.Quote(terms.local),
            "1": fel.codon_data_info[terms.code]
        },
        fel.filter_names,
        None);


fel.site_model_mapping = {"fel_mg" : fel.site.mg_rev};

/* set up the local constraint model */

fel.alpha = model.generic.GetLocalParameter (fel.site.mg_rev, utility.getGlobalValue("terms.parameters.synonymous_rate"));
fel.beta = model.generic.GetLocalParameter (fel.site.mg_rev, utility.getGlobalValue("terms.parameters.nonsynonymous_rate"));
io.CheckAssertion ("None!=fel.alpha && None!=fel.beta", "Could not find expected local synonymous and non-synonymous rate parameters in \`estimators.FitMGREV\`");

selection.io.startTimer (fel.json [terms.json.timers], "FEL analysis", 2);


model.generic.AddGlobal (fel.site.mg_rev, fel.alpha.scaler, fel.site_alpha);


parameters.DeclareGlobal (fel.alpha.scaler, {});
parameters.DeclareGlobal (fel.scaler_parameter_names, {});

utility.ForEachPair (fel.scaler_parameter_names, "_group_", "_name_",
'
    model.generic.AddGlobal (fel.site.mg_rev, _name_ , fel.site_tested_classes [_group_]);
');


/* echo to screen calls */


fel.report.pairwise = {{"" + (1+((fel.filter_specification[fel.report.partition])[terms.data.coverage])[fel.report.site]),
                                    Format(fel.report.row[0],10,3),
                                    Format(fel.report.row[fel.report.rate1_index],10,3) + " : " + Format(fel.report.row[fel.report.rate2_index],10,3),
                                    fel.report.substitutions,
                                    fel.report.test ,
                                    Format(fel.report.row[fel.report.pvalue_index],6,4),
                                    Format(fel.report.row[fel.report.permutation_index],6,4)}};

fel.report.overall = {{"" + (1+((fel.filter_specification[fel.report.partition])[terms.data.coverage])[fel.report.site]),
                                    Format(fel.report.row[0],10,3),
                                    Format(fel.report.rates[0],10,3) + " - " + Format(fel.report.rates[fel.report.rate_count-2],10,3),
                                    fel.report.substitutions,
                                    fel.report.test,
                                    Format(fel.report.row[fel.report.pvalue_index],6,4),
                                    Format(fel.report.row[fel.report.permutation_index],6,4)}};






fel.site_results = {};
fel.partition_index = 0;

fel.report.header_done = FALSE;
fel.table_output_options[terms.table_options.header] = TRUE;
model.ApplyModelToTree( "fel.site_tree", fel.trees[fel.partition_index], {terms.default : fel.site.mg_rev}, None);

fel.case_respecting_node_names = trees.branch_names (fel.site_tree, TRUE);
fel.site_patterns = alignments.Extract_site_patterns ((fel.filter_specification[fel.partition_index])[utility.getGlobalValue("terms.data.name")]);



// apply constraints to the site tree
// alpha = alpha_scaler * branch_length
// beta  = beta_scaler_test * branch_length or beta_nuisance_test * branch_length

utility.ForEach (fel.case_respecting_node_names, "_node_",
        '_node_class_ = (fel.selected_branches[fel.partition_index])[_node_];
         _beta_scaler = fel.scaler_parameter_names[_node_class_];
         fel.apply_proportional_site_constraint ("fel.site_tree", _node_, fel.alpha, fel.beta, fel.alpha.scaler, _beta_scaler, (( fel.final_partitioned_mg_results[terms.branch_length])[fel.partition_index])[_node_]);
    ');



// create the likelihood function for this site
ExecuteCommands (alignments.serialize_site_filter
                                   ((fel.filter_specification[fel.partition_index])[utility.getGlobalValue("terms.data.name")],
                                   ((fel.site_patterns[0])[utility.getGlobalValue("terms.data.sites")])[0],
               ));

__make_filter ("fel.site_filter");

LikelihoodFunction fel.site_likelihood = (fel.site_filter, fel.site_tree);



estimators.ApplyExistingEstimates ("fel.site_likelihood", fel.site_model_mapping, fel.final_partitioned_mg_results,
                                    terms.globals_only);


fel.queue = mpi.CreateQueue ({"LikelihoodFunctions": {{"fel.site_likelihood"}},
                               "Models" : {{"fel.site.mg_rev"}},
                               "Headers" : {{"libv3/all-terms.bf","libv3/tasks/ancestral.bf", "libv3/convenience/math.bf"}},
                               "Functions" : {{"fel.apply_proportional_site_constraint"}},
                               "Variables" : {{"terms.fel.test_keys","fel.permutations", "fel.alpha","fel.beta","fel.alpha.scaler","terms.fel.permutation","fel.final_partitioned_mg_results","fel.srv","fel.site_tested_classes","fel.scaler_parameter_names","fel.branches.testable","fel.branches.has_background","fel.alpha.scaler","terms.fel.pairwise","fel.branch_class_counter","fel.report.test_count", "fel.p_value","fel.site.permutations"}}
                             });


/* run the main loop over all unique site pattern combinations */
utility.ForEachPair (fel.site_patterns, "_pattern_", "_pattern_info_",
    '
        if (_pattern_info_[terms.data.is_constant]) {
            fel.store_results (-1,None,{"0" : "fel.site_likelihood",
                                                            "1" : None,
                                                            "2" : fel.partition_index,
                                                            "3" : _pattern_info_,
                                                            "4" : fel.site_model_mapping,
                                                            "5" : fel.selected_branches[fel.partition_index],
                                                            "6" : fel.scaler_parameter_names,
                                                            "7" : FALSE
                                                             });
        } else {
            mpi.QueueJob (fel.queue, "fel.handle_a_site", {"0" : "fel.site_likelihood",
                                                            "1" : alignments.serialize_site_filter
                                                               ((fel.filter_specification[fel.partition_index])[terms.data.name],
                                                               (_pattern_info_[utility.getGlobalValue("terms.data.sites")])[0]),
                                                            "2" : fel.partition_index,
                                                            "3" : _pattern_info_,
                                                            "4" : fel.site_model_mapping,
                                                            "5" : fel.selected_branches[fel.partition_index],
                                                            "6" : fel.scaler_parameter_names,
                                                            "7" : FALSE
                                                            },
                                                            "fel.store_results");
        }
    '
);

mpi.QueueComplete (fel.queue);
fel.partition_matrix = {Abs (fel.site_results[fel.partition_index]), Rows (fel.table_headers)};

utility.ForEachPair (fel.site_results[fel.partition_index], "_key_", "_value_",
'
    for (fel.index = 0; fel.index < Rows (fel.table_headers); fel.index += 1) {
        fel.partition_matrix [0+_key_][fel.index] = _value_[fel.index];
    }
'
);

fel.site_results[fel.partition_index] = fel.partition_matrix;

/* compute Benjamini-Hochberg */

fel.bh.pv = {};

fel.overall.index = fel.tests.key[0][3];
utility.ForEachPair ((fel.site_results[0])[-1][fel.overall.index], "_idx_", "_value_", "fel.bh.pv[_idx_[0]]=_value_");
fel.bh.pv = math.BenjaminiHochbergFDR(fel.bh.pv);
fel.bh.count = {};

utility.ForEachPair (fel.bh.pv, "_index_", "_value_", 
'
    if (fel.q_value >= _value_) {
        fel.bh.count [1 + _index_] = _value_;
    }
    (fel.site_results[0])[0+_index_][fel.overall.index+1] = _value_;
');

for (fel.k = 0; fel.k < fel.report.test_count; fel.k += 1) {
    io.ReportProgressMessageMD ("fel", "results", "** Found _" + fel.report.counts[fel.k] + "_ " + io.SingularOrPlural (fel.report.counts[fel.k], "site", "sites") + " with different _" + fel.tests.key[fel.k][0] +"_ dN/dS at p <= " + fel.p_value + "**");
}

io.ReportProgressMessageMD ("fel", "FDR", "### False discovery rate correction");
 
if (utility.Array1D (fel.bh.count)) {
    io.ReportProgressMessageMD ("fel", "FDR", "There are " + utility.Array1D (fel.bh.count) + " sites where the overall p-value passes the False Discovery Rate threshold of " + fel.q_value);
    console.log ("");
    fel.bh.count.mx = {Abs (fel.bh.count), 2};
    fel.i = 0;
    utility.ForEachPair (fel.bh.count, "_key_", "_value_", "
        fel.bh.count.mx[fel.i][0] = +_key_;
        fel.bh.count.mx[fel.i][1] = +_value_;
        fel.i += 1;

    ");

    fel.bh.count.mx =  (fel.bh.count.mx%1);

    fel.table_screen_output.qv  = {{"Codon", "q-value"}};
    fel.table_output_options.qv = {terms.table_options.header : TRUE, terms.table_options.align : "center",
                                terms.table_options.column_widths : { "0" : 15, "1" : 22}};

    fprintf (stdout, io.FormatTableRow (fel.table_screen_output.qv,fel.table_output_options.qv));
    fel.table_output_options.qv[terms.table_options.header] = FALSE;

    for (fel.i = 0; fel.i < Rows (fel.bh.count.mx ); fel.i += 1) {
        fprintf (stdout, io.FormatTableRow (
                {{Format (fel.bh.count.mx[fel.i][0],6,0), Format (fel.bh.count.mx[fel.i][1],16,10)}}
                ,fel.table_output_options.qv));
    }

} else {
    io.ReportProgressMessageMD ("fel", "FDR", "There are no sites where the overall p-value passes the False Discovery Rate threshold of " + fel.q_value);
}


fel.json [terms.json.MLE ] = {
                                terms.json.headers   : fel.table_headers,
                                terms.json.content   : fel.site_results 
                              };

selection.io.stopTimer (fel.json [terms.json.timers], "Total time");
selection.io.stopTimer (fel.json [terms.json.timers], "FEL analysis");

io.SpoolJSON (fel.json, fel.codon_data_info[terms.json.json]);


//------------------------------------------------------------------------------
lfunction fel.select_branches(partition_info) {

    io.CheckAssertion("utility.Array1D (`&partition_info`) == 1", "FEL-contrast only works on a single partition dataset");
    available_models = {};
    branch_set = {};


    tree_for_analysis = (partition_info[0])[utility.getGlobalValue("terms.data.tree")];
    utility.ForEach (tree_for_analysis[utility.getGlobalValue("terms.trees.model_map")], "_value_", "`&available_models`[_value_] += 1");
    list_models   = utility.Keys   (available_models); // get keys
    branch_counts = utility.UniqueValues (available_models);
    option_count  = Abs (available_models) + 3;

    //io.CheckAssertion("`&option_count` >= 2", "FEL-contrast requires at least one designated set of branches in the tree.");

    internal_count = utility.Array1D(utility.Filter (tree_for_analysis[^"terms.trees.partitioned"], "_value_", "_value_ == ^'terms.tree_attributes.internal'"));
    leaf_count     = utility.Array1D(utility.Filter (tree_for_analysis[^"terms.trees.partitioned"], "_value_", "_value_ == ^'terms.tree_attributes.leaf'"));
    //leaf_count     = 

    selectTheseForTesting = {
        option_count, 2
    };
    
    selectTheseForTesting [0][0] = "Internal branches";
    selectTheseForTesting [0][1] = "Set of all `internal_count` internal branches in the tree";
    selectTheseForTesting [1][0] = "Terminal branches";
    selectTheseForTesting [1][1] = "Set of all `leaf_count` terminal branches in the tree";
    selectTheseForTesting [2][0] = "Random set of branches";
    selectTheseForTesting [2][1] = "Partition the tree into two random sets of branches (e.g., for null hypothesis testing); this options supersedes all others";

    for (k = 3; k < option_count; k += 1) {
        if (list_models[k-3] != "") {
            selectTheseForTesting[k][0] = list_models[k-3];
            selectTheseForTesting[k][1] = "Set " + list_models[k-3] + " with " + available_models[list_models[k-3]] + " branches";
        } else {
            selectTheseForTesting[k][0] = "Unlabeled branches";
            selectTheseForTesting[k][1] = "Set of " + available_models[list_models[k-3]] + " unlabeled branches";
        }
    }

    ChoiceList(testSet, "Choose sets of branches to compare. If more than one set is chosen, pairwise comparisons will be carried out in addition to a group-level difference test.", 0, NO_SKIP, selectTheseForTesting);
    io.CheckAssertion ("`&testSet[0]` >= 0", "User cancelled branch selection; analysis terminating");

    return_set = {};
    tree_configuration = {};
    tree_for_analysis  = (partition_info[0])[utility.getGlobalValue("terms.data.tree")];
    branch_set_count   = utility.Array1D (testSet);
    test_sets          = {};

    for (k = 0; k < branch_set_count; k+=1) {
        tag_test = selectTheseForTesting [testSet[k]][0];
        if (tag_test == "Unlabeled branches") {
            tag_test = "_unlabeled_";
        }
        test_sets[tag_test] = TRUE;
    }
    
    if (test_sets / 'Random set of branches') {
        KeywordArgument ("random-subset", "How many branches in the random subset to test", (internal_count+leaf_count)$2);
        test_count = io.PromptUser ("How many branches in the test set", (internal_count+leaf_count)$2,1,(internal_count+leaf_count), TRUE);
        branch_names = utility.Keys (tree_for_analysis[^"terms.trees.partitioned"]);
        test_set = Random ({1,internal_count+leaf_count}["_MATRIX_ELEMENT_COLUMN_"],0);
        for (k = 0; k < test_count; k+=1) {
            tree_configuration[branch_names[test_set[k]]] = ^'terms.tree_attributes.test';
        }
        test_count = utility.Array1D (test_set);
        for (;k < test_count; k+=1) {
            tree_configuration[branch_names[test_set[k]]] = ^'terms.tree_attributes.background';
        }
    } else {
        utility.ForEachPair(tree_for_analysis[^"terms.trees.partitioned"], "_key_","_value_", "
            if (`&test_sets` / 'Terminal branches' && _value_ == ^'terms.tree_attributes.leaf') {
                `&tree_configuration`[_key_] = _value_;
            }
            if (`&test_sets` / 'Internal branches' && _value_ == ^'terms.tree_attributes.internal') {
                `&tree_configuration`[_key_] = _value_;
            }
        
        ");   
        utility.ForEachPair (tree_for_analysis[utility.getGlobalValue("terms.trees.model_map")], "_key_", "_value_", "
            if (`&test_sets`[_value_]) {
                io.CheckAssertion ('(`&tree_configuration`/_key_)==FALSE', 'Branch ' + _key_ + ' belongs to multiple sets');
                `&tree_configuration`[_key_] = _value_;
            } else {
                if ((`&tree_configuration`/_key_)==FALSE) {
                    `&tree_configuration`[_key_] = utility.getGlobalValue('terms.tree_attributes.background');
                }
            }
        ");
    }

    return_set + tree_configuration;
    return return_set;
}

//----------------------------------------------------------------------------------------
function fel.apply_proportional_site_constraint (tree_name, node_name, alpha_parameter, beta_parameter, alpha_factor, beta_factor, branch_length) {

    fel.branch_length = (branch_length[terms.parameters.synonymous_rate])[terms.fit.MLE];

    node_name = tree_name + "." + node_name;

    ExecuteCommands ("
        `node_name`.`alpha_parameter` := (`alpha_factor`) * fel.branch_length__;
        `node_name`.`beta_parameter`  := (`beta_factor`)  * fel.branch_length__;
    ");
}
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
lfunction fel.handle_a_site (lf, filter_data, partition_index, pattern_info, model_mapping, branch_sets, parameter_names, permutation) {

    utility.ForEach (utility.Keys (branch_sets), "_node_",
        '_node_class_ = (`&branch_sets`)[_node_];
         _beta_scaler = (`&parameter_names`)[_node_class_];
         fel.apply_proportional_site_constraint ("fel.site_tree", _node_, fel.alpha, fel.beta, fel.alpha.scaler, _beta_scaler, 
            (( fel.final_partitioned_mg_results[terms.branch_length])[partition_index])[_node_]);
    ');

    GetString (lfInfo, ^lf,-1);    
    ExecuteCommands (filter_data);
    __make_filter ((lfInfo["Datafilters"])[0]);

    utility.SetEnvVariable ("USE_LAST_RESULTS", TRUE);

    if (^"fel.srv"){
        ^(^"fel.alpha.scaler") = 1;
    } else
    {
        ^(^"fel.alpha.scaler") := 1;
    }

    // all rates free
    
    start.grid_init = {};
    pnames = {};
    
    for (k, _pname_; in; ^"fel.scaler_parameter_names") {
        ^_pname_ = 1;
        pnames + _pname_;
        start.grid_init [_pname_] = {"0.1" : 1, "1" : 1};
    }
    
    start.grid = {};
    for (k, v; in; utility.CatersianProduct (start.grid_init)) {
        values = {};
        for (vi, vv; in; v) {
            values[pnames[vi]] = +vv;
        }
        start.grid + values;
    }
        
    Optimize (results, ^lf, {
                "OPTIMIZATION_METHOD" : "nedler-mead", 
                "OPTIMIZATION_PRECISION" : 1e-5,
                "OPTIMIZATION_START_GRID" : start.grid    
            }
           );
    
    /**
        infer and report ancestral substitutions
    */
    
    if (!permutation) {
        ancestors   = ancestral._buildAncestralCacheInternal (lf, 0, FALSE, FALSE);
        counts      = ancestral.ComputeSubstitutionCounts (ancestors, null, null, null);
        branch_map  = utility.SwapKeysAndValues (counts["Branches"]);
        counts_by_branch_set = {};
    
    
        utility.ForEach (counts["Branches"], "_name_", '
            `&counts_by_branch_set`[`&branch_sets`[_name_]] += ((`&counts`)["Counts"])[+`&branch_map`[_name_]];
        ');
    }
    
    
            
    snapshot = estimators.TakeLFStateSnapshot (lf);
    alternative = estimators.ExtractMLEs (lf, model_mapping);
    alternative [utility.getGlobalValue("terms.fit.log_likelihood")] = results[1][0];


    sum = + utility.Map (^"fel.scaler_parameter_names", "_pname_",
        '^_pname_'
    );

    testable = utility.Array1D (^"fel.branches.testable");
    denominator = testable;
    if (^"fel.branches.has_background") {
        denominator = testable + 1;
    }

    // baseline NULL (everything = same rate)

    if (^"fel.srv") {
        ^(^"fel.alpha.scaler") = Min (10, (^(^"fel.alpha.scaler")+denominator*sum)/denominator);
    }

    ref_parameter = (^"fel.scaler_parameter_names")[(^"fel.branches.testable")["VALUEINDEXORDER"][0]];

    ^ref_parameter = sum /denominator;

    if (testable == 1) {
        parameters.SetConstraint ((^"fel.scaler_parameter_names")[^"terms.tree_attributes.background"],ref_parameter, "");
    } else {
        utility.ForEach (^"fel.branches.testable", "_gname_",
        '
            //console.log ("REF " + `&ref_parameter`);
            _pname_ =  (^"fel.scaler_parameter_names")[_gname_];
            if (_pname_ != `&ref_parameter`) {
            	//console.log (_pname_ + "=>" + `&ref_parameter`);
                parameters.SetConstraint (_pname_,`&ref_parameter`, "");
            }
        '
        );
    }

    Optimize (results, ^lf);
    Null = estimators.ExtractMLEs (lf, model_mapping);
    Null [utility.getGlobalValue("terms.fit.log_likelihood")] = results[1][0];

    if (testable > 2) {
        pairwise = {};
        for (v = 0; v < testable; v+=1) {
            for (v2 = v + 1; v2 < testable; v2+=1) {
                v1n = (^"fel.branches.testable")[v];
                v2n = (^"fel.branches.testable")[v2];

                estimators.RestoreLFStateFromSnapshot (lf_id, snapshot);
                parameters.SetConstraint ((^"fel.scaler_parameter_names")[v1n],(^"fel.scaler_parameter_names")[v2n], "");
                Optimize (results, ^lf);
                pairwise[v1n + "|" + v2n] = estimators.ExtractMLEs (lf, model_mapping);
                (pairwise[v1n + "|" + v2n])[utility.getGlobalValue("terms.fit.log_likelihood")] = results[1][0];
           }
        }
    } else {
        pairwise = None;
    }
    
    lrt = math.DoLRT (Null[utility.getGlobalValue("terms.fit.log_likelihood")],
                      alternative[utility.getGlobalValue("terms.fit.log_likelihood")],
                      Max (1,^'fel.branch_class_counter' - 1));

    p_values = {};
    p_values ["overall"] = lrt[^'terms.p_value'];
    test_keys = {};

    if (^'fel.report.test_count' > 1) {

        for (v = 0; v < ^'fel.branch_class_counter'; v += 1) {
            for (v2 = v + 1; v2 < ^'fel.branch_class_counter'; v2 += 1) {
                test_key = (^'fel.branches.testable') [v] + "|" + (^'fel.branches.testable') [v2];
                test_keys + test_key;
                lrt = math.DoLRT ((pairwise[test_key])[utility.getGlobalValue("terms.fit.log_likelihood")],
                                   alternative[utility.getGlobalValue("terms.fit.log_likelihood")],
                                  1);
                p_values [test_key] = lrt[^'terms.p_value'];
            }
        }
    }
    

    p_values = math.HolmBonferroniCorrection (p_values);
        
    if (permutation == FALSE && (Min (p_values,0))["value"] <= ^'fel.p_value') {
        result = {
                    utility.getGlobalValue("terms.alternative") : alternative,
                    utility.getGlobalValue("terms.Null"): Null,
                    utility.getGlobalValue("terms.fel.pairwise"): pairwise,
                    utility.getGlobalValue("terms.substitutions") : counts_by_branch_set,
                    utility.getGlobalValue("terms.p_value") : p_values,
                    utility.getGlobalValue("terms.fel.test_keys") : test_keys
                };
                
            
        if (^'fel.permutations' == TRUE) {
            p_min = (Min (p_values,0))["value"];
            //console.log ("Entering permutation mode with p = " + p_min);
            perm_p_values = {};
            for (rep = 0; rep < ^'fel.site.permutations'; rep+=1) {
                this_set = (fel.handle_a_site (lf, filter_data, partition_index, pattern_info, model_mapping, Random(branch_sets,0), parameter_names, 1))[utility.getGlobalValue("terms.p_value")];
                //console.log (this_set);
                perm_p_values + this_set;
                //console.log ("" + rep + " => " + (Min (this_set,0))["value"]);
                if ((Min (this_set,0))["value"] <= p_min) {
                    break;
                }
            }
            result [utility.getGlobalValue("terms.fel.permutation")]= perm_p_values;
        }
        return result;
    } 
    
    return {
                utility.getGlobalValue("terms.alternative") : alternative,
                utility.getGlobalValue("terms.Null"): Null,
                utility.getGlobalValue("terms.fel.pairwise"): pairwise,
                utility.getGlobalValue("terms.substitutions") : counts_by_branch_set,
                utility.getGlobalValue("terms.fel.test_keys") : test_keys,
                utility.getGlobalValue("terms.p_value") : p_values
           };
}


function fel.report.echo (fel.report.site, fel.report.partition, fel.report.row) {

    fel.k.bound   = Rows (fel.tests.key);

    //console.log (fel.report.row);
    //console.log (fel.tests.key);

    for (fel.k = 0; fel.k < fel.k.bound;  fel.k += 1) {
        fel.print_row = None;
        if (fel.report.row[fel.tests.key[fel.k][3]] <= fel.p_value) {
            fel.report.test         = fel.tests.key[fel.k][0];
            fel.report.rate1_index  = fel.tests.key[fel.k][1];
            fel.report.pvalue_index = fel.tests.key[fel.k][3];
            
            fel.report.permutation_index = Columns (fel.report.row) - 2;
            
            if (fel.report.rate1_index >= 0) {
                fel.report.rate2_index = fel.tests.key[fel.k][2];
                fel.print_row = fel.report.pairwise;
                fel.report.substitutions = "" + fel.report.row[-1 + ^"fel.report.rate_count" + fel.report.rate1_index] + ", " +
                                                fel.report.row[-1 + ^"fel.report.rate_count" + fel.report.rate2_index];
                
            } else {
                fel.print_row = fel.report.overall;
                fel.report.rates = Transpose (fel.report.row[{{0,1}}][{{0,fel.report.rate_count-1}}]) % 0;
                fel.report.substitutions = "" + fel.report.row[^"fel.report.rate_count"];
                for (fel.i = 1; fel.i < fel.report.subs_kinds; fel.i +=1) {
                     fel.report.substitutions += ", " + fel.report.row[^"fel.report.rate_count" + fel.i];
                }
                
            }
            fel.report.counts [fel.k] += 1;
        }

        if (None != fel.print_row) {
            if (!fel.report.header_done) {
                io.ReportProgressMessageMD("FEL", "" + fel.report.partition, "For partition " + (fel.report.partition+1) + " these sites are significant at p <=" + fel.p_value + "\n");
                fprintf (stdout,
                    io.FormatTableRow (fel.table_screen_output,fel.table_output_options));
                fel.report.header_done = TRUE;
                fel.table_output_options[terms.table_options.header] = FALSE;
            }
            fprintf (stdout,
                io.FormatTableRow (fel.print_row,fel.table_output_options));
        }
    }
}


lfunction fel.store_results (node, result, arguments) {
    //console.log (^"fel.table_headers");
        
    partition_index = arguments [2];
    pattern_info    = arguments [3];

    array_size  = utility.getGlobalValue ("fel.report.rate_count") + utility.getGlobalValue ("fel.report.subs_kinds") + utility.getGlobalValue ("fel.report.test_count") + 3;
    result_row  = { 1, array_size } ["_MATRIX_ELEMENT_COLUMN_>=^'fel.report.rate_count'&&_MATRIX_ELEMENT_COLUMN_<array_size-1"];


    if (None != result) { // not a constant site

        result_row [0] = estimators.GetGlobalMLE (result[utility.getGlobalValue("terms.alternative")], ^"fel.site_alpha");
        k = 1;
        utility.ForEach (^"fel.rate.names", "_n_",
            '
                `&result_row`[`&k`] = ( estimators.GetGlobalMLE ((`&result`[utility.getGlobalValue("terms.alternative")], fel.site_tested_classes [_n_])));
                 `&k` += 1;
           '
        );

        utility.ForEach (^"fel.branches.testable", "_n_",
            '
                `&result_row`[`&k`] = (`&result`[utility.getGlobalValue("terms.substitutions")])[_n_];
                 `&k` += 1;
           '
        );
        

        
        p_values = result[^"terms.p_value"];
        result_row [k] = p_values ["overall"];
        k += 1;
        result_row [k] = 0; // q-value to be computed later
         
        test_keys = result[^"terms.fel.test_keys"];
        
        for (v = 1; v < ^'fel.report.test_count'; v+=1) {
            k += 1;
            result_row [k] = p_values[test_keys[v-1]];
        }
        

        k += 1;
        
        alternative_lengths = utility.Array1D(result[^"terms.fel.permutation"]);

        if (alternative_lengths) {
            result_row [k] = 1./alternative_lengths;
        } else {
            result_row [k] = -1;
        }

        k += 1;
        
        sum = 0;
    
        alternative_lengths = ((result[utility.getGlobalValue("terms.alternative")])[utility.getGlobalValue("terms.branch_length")])[0];

        utility.ForEach (^"fel.case_respecting_node_names", "_node_",
                '`&sum` += ((`&alternative_lengths`)[_node_])[utility.getGlobalValue("terms.fit.MLE")];');

        result_row [k] = sum;
        
      }
      
    utility.EnsureKey (^"fel.site_results", partition_index);
    utility.ForEach (pattern_info[utility.getGlobalValue("terms.data.sites")], "_fel_result_",
        '
            (fel.site_results[`&partition_index`])[_fel_result_] = `&result_row`;
            fel.report.echo (_fel_result_, `&partition_index`, `&result_row`);
        '
    );


    //assert (0);
}
