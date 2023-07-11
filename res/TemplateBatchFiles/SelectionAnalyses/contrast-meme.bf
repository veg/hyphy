RequireVersion("2.5.50");

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

LoadFunctionLibrary("libv3/models/codon/BS_REL.bf");


//utility.ToggleEnvVariable ("OPTIMIZATION_PRECISION", 1);
//utility.ToggleEnvVariable ("OPTIMIZATION_TIME_HARD_LIMIT", 1);


/*------------------------------------------------------------------------------ Display analysis information
*/

c_meme.analysis_description = {
    terms.io.info: "Contrast-MEME (Mixed Effects Model of Evolution) investigates whether or not selective pressures differ between two or more sets of
    branches at a site. Site-specific synonymous (alpha) and non-synonymous (beta, one per branch set) substitution rates are estimated
    and then beta rates are tested for equality at each site. LRT and permutation tests ar used to assess significance at each site, and FDR is applied alignment wide to call sites with different selective profiles",
    terms.io.version: "0.5",
    terms.io.reference: "TBD",
    terms.io.authors: "Sergei L Kosakovsky Pond and Steven Weaver",
    terms.io.contact: "spond@temple.edu",
    terms.io.requirements: "in-frame codon alignment and a phylogenetic tree; only single partition data are supported",
    terms.io.help : "http://www.hyphy.org/methods/other/contrast-c_meme/"
};

io.DisplayAnalysisBanner(c_meme.analysis_description);


/*------------------------------------------------------------------------------
    Environment setup
*/

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);

/*------------------------------------------------------------------------------
    Globals
*/



terms.c_meme.pairwise      = "pairwise";
terms.c_meme.permutation   = "permutation";
terms.c_meme.test_keys     = "test keys";
c_meme.site_alpha          = "Site relative synonymous rate";
c_meme.site_tested_classes = {};
c_meme.alpha.scaler        = "c_meme.alpha_scaler";
c_meme.site.permutations   = 20;


// default cutoff for printing to screen
c_meme.p_value = 0.1;
c_meme.q_value = 0.2;
c_meme.scaler_prefix = "c_meme.scaler";


KeywordArgument ("code",      "Which genetic code should be used", "Universal");
KeywordArgument ("alignment", "An in-frame codon alignment in one of the formats supported by HyPhy");
KeywordArgument ("tree",      "A phylogenetic tree (optionally annotated with {})", null, "Please select a tree file for the data:");
KeywordArgument ("branch-set", "The set of branches to use for testing");


// The dictionary of results to be written to JSON at the end of the run
c_meme.json = {
    terms.json.analysis: c_meme.analysis_description,
    terms.json.input: {},
    terms.json.fits: {},
    terms.json.timers: {},
};

c_meme.display_orders =   {terms.original_name: -1,
                        terms.json.nucleotide_gtr: 0,
                        terms.json.global_mg94xrev: 1
                       };



selection.io.startTimer (c_meme.json [terms.json.timers], "Total time", 0);



namespace c_meme {
    LoadFunctionLibrary ("modules/shared-load-file.bf");
    load_file ({utility.getGlobalValue("terms.prefix"): "c_meme", 
                utility.getGlobalValue("terms.settings") : {utility.getGlobalValue("terms.settings.branch_selector") : "c_meme.select_branches"}});
}

KeywordArgument ("srv", "Include synonymous rate variation in the model", "Yes");
KeywordArgument ("permutations", "Perform permutation significance tests", "Yes");
KeywordArgument ("p-value", "Significance value for site-tests", "0.05");
KeywordArgument ("q-value", "Significance value for FDR reporting", "0.20");


/* Prompt for one-rate or two-rate analysis */
c_meme.srv = io.SelectAnOption( {{"Yes", "[Recommended] Consider synonymous rate variation (dS varies across sites)."}, {"No", "Ignore synonymous rate variation (dS := 1 at each site)."}},
                                  "Use synonymous rate variation? Strongly recommended YES for selection inference.");

if (c_meme.srv == "Yes"){
    c_meme.srv = TRUE
} else {
    c_meme.srv = FALSE
}

/* Prompt for one-rate or two-rate analysis */
c_meme.permutations = io.SelectAnOption( {{"Yes", "For sites with significant p-values, perform an additional permutation test (over branch assignments) to assess significance. Adds computational cost, reduces false positves"}, 
                                       {"No", "Do not perform additional tests (faster, higher risk of false positves)"}},
                                  "Perform permutation significance tests");

if (c_meme.permutations == "Yes"){
    c_meme.permutations = TRUE
} else {
    c_meme.permutations = FALSE
}

/* Prompt for p value threshold */
c_meme.p_value  = io.PromptUser ("\n>Select nominal p-value threshold to use when testing for selection (FDR correction will be performed on all sites)",0.1,0,1,FALSE);
c_meme.q_value  = io.PromptUser ("\n>Select nominal the q-value threshold to use when testing for selection (FDR reporting)",0.2,0,1,FALSE);

KeywordArgument ("output", "Write the resulting JSON to this file (default is to save to the same path as the alignment file + 'c_meme.json')", c_meme.codon_data_info [terms.json.json]);
c_meme.codon_data_info [terms.json.json] = io.PromptUserForFilePath ("Save the resulting JSON file to");

io.ReportProgressMessageMD('c_meme',  'selector', 'Branches to use as the test set in the c_meme-contrast analysis');

c_meme.branch_sets = {};

utility.ForEachPair (c_meme.selected_branches[0], "_branch_", "_model_",
"
    utility.EnsureKey (c_meme.branch_sets, _model_);
    c_meme.branch_sets[_model_] + _branch_;
");

c_meme.branch_class_count = utility.Array1D (c_meme.branch_sets);
c_meme.scaler_parameter_names = {};

io.ReportProgressMessageMD('c_meme', 'selector', "Selected `c_meme.branch_class_count` sets of branches to test\n");

c_meme.branch_class_counter = 0;
c_meme.branches.testable = {};
c_meme.branches.has_background = FALSE;


utility.ForEachPair (c_meme.branch_sets, "_group_", "_branches_",
    "
     if (_group_ != terms.tree_attributes.background) {
        c_meme.site_tested_classes [_group_] = {
            c_meme.beta1: 'Site relative non-synonymous rate 1 (' + _group_ + ' branches)',
            c_meme.beta2: 'Site relative non-synonymous rate 2 (' + _group_ + ' branches)',
            c_meme.prop: 'Site mixing proportion (' + _group_ + ' branches)'
        };
        c_meme.branch_class_counter += 1;
        c_meme.branches.testable + _group_;
        c_meme.scaler_parameter_names [_group_] = {
            c_meme.beta1 : 'c_meme.beta1_scaler_group_' + c_meme.branch_class_counter,
            c_meme.beta2 : 'c_meme.beta2_scaler_group_' + c_meme.branch_class_counter,
            c_meme.prop : 'c_meme.prop' + c_meme.branch_class_counter
        };
        io.ReportProgressMessageMD('c_meme',  'selector', '* Selected ' + Abs(_branches_) + ' branches in group _' + _group_ + '_ : \\\`' + Join (', ',_branches_) + '\\\`')
     } else {
        c_meme.scaler_parameter_names [_group_] = {
            c_meme.beta1 : 'c_meme.beta1_beta_scaler_background' ,
            c_meme.beta2 : 'c_meme.beta2_beta_scaler_background' ,
            c_meme.prop : 'c_meme.prop_background'  
        };
        c_meme.site_tested_classes [_group_] = {
            c_meme.beta1: 'Site relative non-synonymous rate 1 (reference branches)',
            c_meme.beta2: 'Site relative non-synonymous rate 2 (reference branches)',
            c_meme.prop: 'Site mixing proportion (reference branches)'
        };       
        c_meme.branches.has_background = TRUE;
        io.ReportProgressMessageMD('c_meme',  'selector', '* ' + Abs(_branches_) + ' branches are in the background group : \\\`' + Join (', ',_branches_) + '\\\`')
     }
     "
);


c_meme.test_count = 1;
if (c_meme.branch_class_counter > 2) {
    c_meme.test_count += c_meme.branch_class_counter * (c_meme.branch_class_counter-1) / 2;
}

c_meme.table_headers = {4 + utility.Array1D (c_meme.scaler_parameter_names)*3 + c_meme.branch_class_counter + c_meme.test_count, 2};
c_meme.table_headers [0][0] = "alpha"; c_meme.table_headers [0][1] = "Synonymous substitution rate at a site";

c_meme.k = 1;
c_meme.rate.names = utility.Keys (c_meme.scaler_parameter_names);

for (_n_; in; c_meme.rate.names) {
    c_meme.table_headers [c_meme.k][0] = "beta1 (" + _n_ + ")";
    c_meme.table_headers [c_meme.k][1] = "Non-synonymous substitution rate 1 at a site for " + _n_ + " branches";
    c_meme.table_headers [c_meme.k+1][0] = "beta2 (" + _n_ + ")";
    c_meme.table_headers [c_meme.k+1][1] = "Non-synonymous substitution rate 2 at a site for " + _n_ + " branches";
    c_meme.table_headers [c_meme.k+2][0] = "prop (" + _n_ + ")";
    c_meme.table_headers [c_meme.k+2][1] = "Mixture proportion at a site for " + _n_ + " branches";
    c_meme.k += 3;
}


c_meme.report.rate_count = c_meme.k;

utility.ForEach (c_meme.branches.testable, "_n_", ,
'
    c_meme.table_headers [c_meme.k][0] = "subs (" + _n_ + ")";
    c_meme.table_headers [c_meme.k][1] = "Substitutions mapped to " + _n_ + " branches";
    c_meme.k += 1;
');

c_meme.report.subs_kinds = utility.Array1D (c_meme.branches.testable);

c_meme.table_headers [c_meme.k][0] = "P-value (overall)";
c_meme.table_headers [c_meme.k+1][0] = "Q-value (overall)";

if (c_meme.branch_class_counter == 1) {
    c_meme.table_headers [c_meme.k][1] = "P-value for the test that " + c_meme.branches.testable[0] + " branches have different non-synonymous rates than the rest of the branches";
    c_meme.table_headers [c_meme.k+1][1] = "Q-value for the test that " + c_meme.branches.testable[0] + " branches have different non-synonymous rates than the rest of the branches";
} else {
    c_meme.table_headers [c_meme.k][1] = "P-value for the test that non-synonymous rates differ between any of the selected groups: " + Join (",", c_meme.branches.testable);
    c_meme.table_headers [c_meme.k+1][1] = "Q-value for the test that non-synonymous rates differ between any of the selected groups: " + Join (",", c_meme.branches.testable);
}



c_meme.k += 1;

c_meme.report.test_count = 1 + (c_meme.branch_class_counter > 2) * (c_meme.branch_class_counter * (c_meme.branch_class_counter-1) / 2);
c_meme.tests.key = {c_meme.report.test_count, 4};
c_meme.tests.key [0][0] = "overall";
c_meme.tests.key [0][1] = -1;
c_meme.tests.key [0][3] = c_meme.report.rate_count + c_meme.report.subs_kinds;

c_meme.test.index = c_meme.test_count > 1;
c_meme.k += 1;


for (c_meme.v = 0; c_meme.v < c_meme.branch_class_counter; c_meme.v += 1) {
    for (c_meme.v2 = c_meme.v + 1; c_meme.v2 < c_meme.branch_class_counter; c_meme.v2 += 1) {
        c_meme.table_headers [c_meme.k][0] = "P-value for " + c_meme.branches.testable[c_meme.v] + " vs " + c_meme.branches.testable[c_meme.v2];
        c_meme.tests.key [c_meme.test.index][0] = c_meme.branches.testable[c_meme.v] + " vs " + c_meme.branches.testable[c_meme.v2];
        c_meme.tests.key [c_meme.test.index][1] = 1 + 3*c_meme.v;
        c_meme.tests.key [c_meme.test.index][2] = 1 + 3*c_meme.v2;

        if (c_meme.test_count > 1) {
            c_meme.tests.key [c_meme.test.index][3] = c_meme.report.rate_count + c_meme.test.index + c_meme.report.subs_kinds + 1;
            c_meme.table_headers [c_meme.k][1] = "P-value for  the test that non-synonymous rates differ between " + c_meme.branches.testable[c_meme.v] + " and " + c_meme.branches.testable[c_meme.v2] + " branches";
            c_meme.k += 1;
        }
        c_meme.test.index += 1;
    }
}


io.ReportProgressMessageMD('c_meme',  'tests', '**' + c_meme.test_count + "** " + io.SingularOrPlural (c_meme.test_count , "test", "tests")+" will be performed at each site");

c_meme.report.counts         = {c_meme.report.test_count,1};

c_meme.table_headers [c_meme.k][0] = "Permutation p-value";
c_meme.table_headers [c_meme.k][1] = "Label permutation test for significant sites";

c_meme.k += 1;


c_meme.table_headers [c_meme.k][0] = "Total branch length";
c_meme.table_headers [c_meme.k][1] = "The total length of branches contributing to inference at this site, and used to scale beta-alpha";


/**
This table is meant for HTML rendering in the results web-app; can use HTML characters, the second column
is 'pop-over' explanation of terms. This is ONLY saved to the JSON file. For Markdown screen output see
the next set of variables.
*/

c_meme.table_screen_output  = {{"Codon", "alpha", "beta1", "beta2", "prop", "substitutions", "test", "LRT p-value", "Permutation p-value"}};

c_meme.table_output_options = {terms.table_options.header : TRUE, terms.table_options.align : "center",
                            terms.table_options.column_widths : { "0" : 8, "1" : 10, "2" : 20, "3" : 20, "4" : 20, "5" : 20, "6" : 30, "7" : 10, "8" : 10}};

selection.io.startTimer (c_meme.json [terms.json.timers], "Model fitting",1);

namespace_tag = "c_meme";

namespace c_meme {
    doGTR ("c_meme");
}


estimators.fixSubsetOfEstimates(c_meme.gtr_results, c_meme.gtr_results[terms.global]);

namespace c_meme {
    doPartitionedMG ("c_meme", FALSE);
}


//c_meme.final_partitioned_mg_results = c_meme.partitioned_mg_results;

io.ReportProgressMessageMD ("c_meme", "codon-refit", "Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model");


c_meme.final_partitioned_mg_results = estimators.FitMGREV (c_meme.filter_names, c_meme.trees, c_meme.codon_data_info [terms.code], {
    terms.run_options.model_type: terms.local,
    terms.run_options.partitioned_omega: c_meme.selected_branches,
    terms.run_options.retain_lf_object: TRUE
}, c_meme.partitioned_mg_results);



io.ReportProgressMessageMD("c_meme", "codon-refit", "* Log(L) = " + Format(c_meme.final_partitioned_mg_results[terms.fit.log_likelihood],8,2));
c_meme.global_dnds = selection.io.extract_global_MLE_re (c_meme.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio);
utility.ForEach (c_meme.global_dnds, "_value_", 'io.ReportProgressMessageMD ("c_meme", "codon-refit", "* " + _value_[terms.description] + " = " + Format (_value_[terms.fit.MLE],8,4));');



//Store MG94 to JSON
selection.io.json_store_lf_withEFV (c_meme.json,
                            terms.json.global_mg94xrev,
                            c_meme.final_partitioned_mg_results[terms.fit.log_likelihood],
                            c_meme.final_partitioned_mg_results[terms.parameters],
                            c_meme.sample_size,
                            utility.ArrayToDict (utility.Map (c_meme.global_dnds, "_value_", "{'key': _value_[terms.description], 'value' : Eval({{_value_ [terms.fit.MLE],1}})}")),
                            (c_meme.final_partitioned_mg_results[terms.efv_estimate])["VALUEINDEXORDER"][0],
                            c_meme.display_orders[terms.json.global_mg94xrev]);


estimators.fixSubsetOfEstimates(c_meme.final_partitioned_mg_results, c_meme.final_partitioned_mg_results[terms.global]);

utility.ForEachPair (c_meme.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(c_meme.json, terms.json.global_mg94xrev, terms.branch_length, c_meme.display_orders[terms.json.global_mg94xrev],
                                             _key_,
                                             selection.io.extract_branch_info((c_meme.final_partitioned_mg_results[terms.branch_length])[_key_], "selection.io.branch.length"));');


selection.io.stopTimer (c_meme.json [terms.json.timers], "Model fitting");

// define the site-level likelihood function

c_meme.site.bsrel =  model.generic.DefineMixtureModel("models.codon.BS_REL_Per_Branch_Mixing.ModelDescription",
        "meme_model", {
            "0": parameters.Quote(terms.local),
            "1": c_meme.codon_data_info[terms.code],
            "2": parameters.Quote (2) // the number of rate classes
        },
        c_meme.filter_names,
        None);
        

c_meme.site_model_mapping = {"meme_model" : c_meme.site.bsrel};

/* set up the local constraint model */

c_meme.alpha = model.generic.GetLocalParameter (c_meme.site.bsrel, utility.getGlobalValue("terms.parameters.synonymous_rate"));
c_meme.beta1 = model.generic.GetLocalParameter (c_meme.site.bsrel, terms.AddCategory (utility.getGlobalValue("terms.parameters.nonsynonymous_rate"),1));
c_meme.beta2 = model.generic.GetLocalParameter (c_meme.site.bsrel, terms.AddCategory (utility.getGlobalValue("terms.parameters.nonsynonymous_rate"),2));
c_meme.prop  = model.generic.GetLocalParameter (c_meme.site.bsrel, terms.AddCategory (utility.getGlobalValue("terms.mixture.mixture_aux_weight"),1));

io.CheckAssertion ("None!=c_meme.alpha && None!=c_meme.beta1  && None!=c_meme.beta2  && None!=c_meme.prop ", "Could not find expected local synonymous and non-synonymous rate parameters in the MEME BS-REL model");


selection.io.startTimer (c_meme.json [terms.json.timers], "c_meme analysis", 2);


model.generic.AddGlobal (c_meme.site.bsrel, c_meme.alpha.scaler, c_meme.site_alpha);

parameters.DeclareGlobal (c_meme.alpha.scaler, {});


for (_group_, _name_; in; c_meme.scaler_parameter_names) {
    parameters.DeclareGlobal (_name_, {});
    for (_gr_, _p_; in; _name_) {
        model.generic.AddGlobal (c_meme.site.bsrel, _p_ , ((c_meme.site_tested_classes[_group_]) [_gr_]));
    }
    //console.log (_name_ ["c_meme.prop"]);
    parameters.SetRange (_name_ ["c_meme.prop"], terms.range_almost_01); 
}


/* echo to screen calls */


c_meme.report.pairwise = {{"" + (1+((c_meme.filter_specification[c_meme.report.partition])[terms.data.coverage])[c_meme.report.site]),
                                    Format(c_meme.report.row[0],0,2),
                                    Format(c_meme.report.row[c_meme.report.rate1_index],0,2) + " : " + Format(c_meme.report.row[c_meme.report.rate2_index],0,2),
                                    Format(c_meme.report.row[c_meme.report.rate1_index+1],0,2) + " : " + Format(c_meme.report.row[c_meme.report.rate2_index+1],0,2),
                                    Format(c_meme.report.row[c_meme.report.rate1_index+2],0,2) + " : " + Format(c_meme.report.row[c_meme.report.rate2_index+2],0,2),
                                    c_meme.report.substitutions,
                                    c_meme.report.test ,
                                    Format(c_meme.report.row[c_meme.report.pvalue_index],0,4),
                                    Format(c_meme.report.row[c_meme.report.permutation_index],0,4)}};

c_meme.report.overall = {{"" + (1+((c_meme.filter_specification[c_meme.report.partition])[terms.data.coverage])[c_meme.report.site]),
                                    Format(c_meme.report.row[0],0,2),
ormat(c_meme.report.row[c_meme.report.rate1_index],0,2) + " : " + Format(c_meme.report.row[c_meme.report.rate2_index],0,2),
                                    Format(c_meme.report.row[c_meme.report.rate1_index+1],0,2) + " : " + Format(c_meme.report.row[c_meme.report.rate2_index+1],0,2),
                                    Format(c_meme.report.row[c_meme.report.rate1_index+2],0,2) + " : " + Format(c_meme.report.row[c_meme.report.rate2_index+2],0,2),                                    c_meme.report.substitutions,
                                    c_meme.report.test,
                                    Format(c_meme.report.row[c_meme.report.pvalue_index],0,4),
                                    Format(c_meme.report.row[c_meme.report.permutation_index],0,4)}};






c_meme.site_results = {};
c_meme.partition_index = 0;

c_meme.report.header_done = FALSE;
c_meme.table_output_options[terms.table_options.header] = TRUE;
model.ApplyModelToTree( "c_meme.site_tree", c_meme.trees[c_meme.partition_index], {terms.default : c_meme.site.bsrel}, None);

c_meme.case_respecting_node_names = trees.branch_names (c_meme.site_tree, TRUE);
c_meme.site_patterns = alignments.Extract_site_patterns ((c_meme.filter_specification[c_meme.partition_index])[utility.getGlobalValue("terms.data.name")]);



// apply constraints to the site tree
// alpha = alpha_scaler * branch_length
// beta  = beta_scaler_test * branch_length or beta_nuisance_test * branch_length

c_meme.lengths_by_class = {};

for (_node_; in; c_meme.case_respecting_node_names) {
    _node_class_ = (c_meme.selected_branches[c_meme.partition_index])[_node_];
    _beta_scaler = c_meme.scaler_parameter_names[_node_class_];
    //console.log (_beta_scaler);
    c_meme.lengths_by_class[_node_class_] += c_meme.apply_proportional_site_constraint ("c_meme.site_tree", _node_, c_meme.alpha, c_meme.beta1, c_meme.beta2, c_meme.prop, c_meme.alpha.scaler, _beta_scaler["c_meme.beta1"],  _beta_scaler["c_meme.beta2"],  _beta_scaler["c_meme.prop"], (( c_meme.final_partitioned_mg_results[terms.branch_length])[c_meme.partition_index])[_node_]);
}


c_meme.ignorable = {};
for (c_meme.t; in; c_meme.branches.testable) {
    if ( c_meme.lengths_by_class[c_meme.t] == 0) {
      fprintf(stdout, "\n-------\n", io.FormatLongStringToWidth(
      ">[WARNING] The cumulative branch length in the _" + c_meme.t + "_ class is 0. 
      Rates along these branches are not identifiable; testing will not be formally conducted (all p-values set to 1).", 72),
      "\n-------\n");        
      for (k; in; c_meme.scaler_parameter_names[c_meme.t]) {
        c_meme.ignorable [k] = 1;
      }
    } else {
        for (k; in; c_meme.scaler_parameter_names[c_meme.t]) {
            c_meme.ignorable [k] = 0;
        }    
    }
}

// create the likelihood function for this site
ExecuteCommands (alignments.serialize_site_filter
                                   ((c_meme.filter_specification[c_meme.partition_index])[utility.getGlobalValue("terms.data.name")],
                                   ((c_meme.site_patterns[0])[utility.getGlobalValue("terms.data.sites")])[0],
               ));

__make_filter ("c_meme.site_filter");

LikelihoodFunction c_meme.site_likelihood = (c_meme.site_filter, c_meme.site_tree);



estimators.ApplyExistingEstimates ("c_meme.site_likelihood", c_meme.site_model_mapping, c_meme.final_partitioned_mg_results,
                                    terms.globals_only);


c_meme.queue = mpi.CreateQueue ({"LikelihoodFunctions": {{"c_meme.site_likelihood"}},
                               "Models" : {{"c_meme.site.bsrel"}},
                               "Headers" : {{"libv3/all-terms.bf","libv3/tasks/ancestral.bf", "libv3/convenience/math.bf", "libv3/tasks/estimators.bf"}},
                               "Functions" : {{"c_meme.apply_proportional_site_constraint"}},
                               "Variables" : {{"terms.c_meme.test_keys","c_meme.permutations","c_meme.ignorable","c_meme.alpha","c_meme.beta1","c_meme.beta2","c_meme.prop", "c_meme.alpha.scaler","terms.c_meme.permutation","c_meme.final_partitioned_mg_results","c_meme.srv","c_meme.site_tested_classes","c_meme.scaler_parameter_names","c_meme.branches.testable","c_meme.branches.has_background","c_meme.alpha.scaler","terms.c_meme.pairwise","c_meme.branch_class_counter","c_meme.report.test_count", "c_meme.p_value","c_meme.site.permutations"}}
                             });

c_meme.pattern_count_all = utility.Array1D (c_meme.site_patterns);
c_meme.pattern_count_this = 0;

/* run the main loop over all unique site pattern combinations */
utility.ForEachPair (c_meme.site_patterns, "_pattern_", "_pattern_info_",
    '
        c_meme.pattern_count_this += 1;
        io.ReportProgressBar("", "Working on site pattern " + (c_meme.pattern_count_this) + "/" +  c_meme.pattern_count_all + " in partition " + (1+c_meme.partition_index));
        
        if (_pattern_info_[terms.data.is_constant]) {
            c_meme.store_results (-1,None,{"0" : "c_meme.site_likelihood",
                                                            "1" : None,
                                                            "2" : c_meme.partition_index,
                                                            "3" : _pattern_info_,
                                                            "4" : c_meme.site_model_mapping,
                                                            "5" : c_meme.selected_branches[c_meme.partition_index],
                                                            "6" : c_meme.scaler_parameter_names,
                                                            "7" : FALSE
                                                             });
        } else {
            mpi.QueueJob (c_meme.queue, "c_meme.handle_a_site", {"0" : "c_meme.site_likelihood",
                                                            "1" : alignments.serialize_site_filter
                                                               ((c_meme.filter_specification[c_meme.partition_index])[terms.data.name],
                                                               (_pattern_info_[utility.getGlobalValue("terms.data.sites")])[0]),
                                                            "2" : c_meme.partition_index,
                                                            "3" : _pattern_info_,
                                                            "4" : c_meme.site_model_mapping,
                                                            "5" : c_meme.selected_branches[c_meme.partition_index],
                                                            "6" : c_meme.scaler_parameter_names,
                                                            "7" : FALSE
                                                            },
                                                            "c_meme.store_results");
        }
    '
);

io.ClearProgressBar();


mpi.QueueComplete (c_meme.queue);
c_meme.partition_matrix = {Abs (c_meme.site_results[c_meme.partition_index]), Rows (c_meme.table_headers)};

for (_key_, _value_; in; c_meme.site_results[c_meme.partition_index]) {
    for (c_meme.index = 0; c_meme.index < Rows (c_meme.table_headers); c_meme.index += 1) {
        c_meme.partition_matrix [0+_key_][c_meme.index] = _value_[c_meme.index];
    }
}


c_meme.site_results[c_meme.partition_index] = c_meme.partition_matrix;

/* compute Benjamini-Hochberg */

c_meme.bh.pv = {};

c_meme.overall.index = c_meme.tests.key[0][3];
utility.ForEachPair ((c_meme.site_results[0])[-1][c_meme.overall.index], "_idx_", "_value_", "c_meme.bh.pv[_idx_[0]]=_value_");
c_meme.bh.pv = math.BenjaminiHochbergFDR(c_meme.bh.pv);
c_meme.bh.count = {};

utility.ForEachPair (c_meme.bh.pv, "_index_", "_value_", 
'
    if (c_meme.q_value >= _value_) {
        c_meme.bh.count [1 + _index_] = _value_;
    }
    (c_meme.site_results[0])[0+_index_][c_meme.overall.index+1] = _value_;
');

for (c_meme.k = 0; c_meme.k < c_meme.report.test_count; c_meme.k += 1) {
    io.ReportProgressMessageMD ("c_meme", "results", "** Found _" + c_meme.report.counts[c_meme.k] + "_ " + io.SingularOrPlural (c_meme.report.counts[c_meme.k], "site", "sites") + " with different _" + c_meme.tests.key[c_meme.k][0] +"_ dN/dS at p <= " + c_meme.p_value + "**");
}

io.ReportProgressMessageMD ("c_meme", "FDR", "### False discovery rate correction");
 
if (utility.Array1D (c_meme.bh.count)) {
    io.ReportProgressMessageMD ("c_meme", "FDR", "There are " + utility.Array1D (c_meme.bh.count) + " sites where the overall p-value passes the False Discovery Rate threshold of " + c_meme.q_value);
    console.log ("");
    c_meme.bh.count.mx = {Abs (c_meme.bh.count), 2};
    c_meme.i = 0;
    utility.ForEachPair (c_meme.bh.count, "_key_", "_value_", "
        c_meme.bh.count.mx[c_meme.i][0] = +_key_;
        c_meme.bh.count.mx[c_meme.i][1] = +_value_;
        c_meme.i += 1;

    ");

    c_meme.bh.count.mx =  (c_meme.bh.count.mx%1);

    c_meme.table_screen_output.qv  = {{"Codon", "q-value"}};
    c_meme.table_output_options.qv = {terms.table_options.header : TRUE, terms.table_options.align : "center",
                                terms.table_options.column_widths : { "0" : 15, "1" : 22}};

    fprintf (stdout, io.FormatTableRow (c_meme.table_screen_output.qv,c_meme.table_output_options.qv));
    c_meme.table_output_options.qv[terms.table_options.header] = FALSE;

    for (c_meme.i = 0; c_meme.i < Rows (c_meme.bh.count.mx ); c_meme.i += 1) {
        fprintf (stdout, io.FormatTableRow (
                {{Format (c_meme.bh.count.mx[c_meme.i][0],6,0), Format (c_meme.bh.count.mx[c_meme.i][1],16,10)}}
                ,c_meme.table_output_options.qv));
    }

} else {
    io.ReportProgressMessageMD ("c_meme", "FDR", "There are no sites where the overall p-value passes the False Discovery Rate threshold of " + c_meme.q_value);
}


c_meme.json [terms.json.MLE ] = {
                                terms.json.headers   : c_meme.table_headers,
                                terms.json.content   : c_meme.site_results 
                              };

selection.io.stopTimer (c_meme.json [terms.json.timers], "Total time");
selection.io.stopTimer (c_meme.json [terms.json.timers], "c_meme analysis");

GetString (_hpv,HYPHY_VERSION,0);
c_meme.json[terms.json.runtime] = _hpv;


io.SpoolJSON (c_meme.json, c_meme.codon_data_info[terms.json.json]);


//------------------------------------------------------------------------------
lfunction c_meme.select_branches(partition_info) {

    io.CheckAssertion("utility.Array1D (`&partition_info`) == 1", "c_meme-contrast only works on a single partition dataset");
    available_models = {};
    branch_set = {};


    tree_for_analysis = (partition_info[0])[utility.getGlobalValue("terms.data.tree")];
    utility.ForEach (tree_for_analysis[utility.getGlobalValue("terms.trees.model_map")], "_value_", "`&available_models`[_value_] += 1");
    list_models   = utility.Keys   (available_models); // get keys
    branch_counts = utility.UniqueValues (available_models);
    option_count  = Abs (available_models) + 3;

    //io.CheckAssertion("`&option_count` >= 2", "c_meme-contrast requires at least one designated set of branches in the tree.");

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
function c_meme.apply_proportional_site_constraint (tree_name, node_name, alpha_parameter, beta1_parameter, beta2_parameter, mix_parameter, alpha_factor, beta1_factor, beta2_factor, mix_factor, branch_length) {

    c_meme.branch_length = (branch_length[terms.parameters.synonymous_rate])[terms.fit.MLE];

    node_name = tree_name + "." + node_name;
    
    //console.log (mix_parameter); console.log (mix_factor);

    ExecuteCommands ("
        `node_name`.`alpha_parameter`  := (`alpha_factor`)  * c_meme.branch_length__;
        `node_name`.`beta1_parameter`  := (`beta1_factor`)  * c_meme.branch_length__;
        `node_name`.`beta2_parameter`  := (`beta2_factor`)  * c_meme.branch_length__;
        `node_name`.`mix_parameter`  := (`mix_factor`);

    ");
    
    return c_meme.branch_length;
}
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
lfunction c_meme.handle_a_site (lf, filter_data, partition_index, pattern_info, model_mapping, branch_sets, parameter_names, permutation) {


    lengths_by_class = {};

    for (_node_; in; utility.Keys (branch_sets)) {
        _node_class_ = branch_sets[_node_];
        _beta_scaler = parameter_names[_node_class_];
        lengths_by_class[_node_class_] += c_meme.apply_proportional_site_constraint ("c_meme.site_tree", _node_, c_meme.alpha, c_meme.beta1, c_meme.beta2, c_meme.prop, c_meme.alpha.scaler, _beta_scaler[c_meme.beta1], _beta_scaler[c_meme.beta2], _beta_scaler[c_meme.prop], 
            (( (^"c_meme.final_partitioned_mg_results")[^"terms.branch_length"])[partition_index])[_node_]);    
    }
    
    ignorable = {};
    for (c_meme.t; in; utility.Values (branch_sets)) {
        if ( lengths_by_class[c_meme.t] == 0) {
          for (pn,pv; in; parameter_names[c_meme.t]) {
            ignorable [pv] = 1;
          }
          
        } else {
            for (pn,pv; in; parameter_names[c_meme.t]) {
                ignorable [pv] = 0;
            }
        }
    }
    

    GetString (lfInfo, ^lf,-1);    
    ExecuteCommands (filter_data);
    __make_filter ((lfInfo["Datafilters"])[0]);

    utility.SetEnvVariable ("USE_LAST_RESULTS", TRUE);

    if (^"c_meme.srv"){
        ^(^"c_meme.alpha.scaler") = 1;
    } else
    {
        ^(^"c_meme.alpha.scaler") := 1;
    }

    // all rates free
    
    pnames = {};
    
    ranges = {};
    
    for (k, _pname_; in; ^"c_meme.scaler_parameter_names") {
        for (p,v; in; _pname_) {
            ^v = 0.5;
            pnames + v;
            ranges [v] = {^"terms.lower_bound" : 0, ^"terms.upper_bound" : 1};
        }
    }

    start.grid_init = {};
    
    for (sp; in; estimators.LHC ( ranges, 50))  {
         start.grid_init + sp;
    }
     
    
    //utility.ToggleEnvVariable ("VERBOSITY_LEVEL", 10);
    
    /*Export (lfe, ^lf);
    console.log (lfe);
    assert (0);*/
    
    Optimize (results, ^lf, {
                "OPTIMIZATION_METHOD" : "nedler-mead", 
                "OPTIMIZATION_PRECISION" : 0.01,
                "OPTIMIZATION_START_GRID" : start.grid    
            }
           );
           
    Optimize (results, ^lf, {
            "OPTIMIZATION_PRECISION" : 0.001
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

    testable = utility.Array1D (^"c_meme.branches.testable");
    
    
    denominator = testable;
    if (^"c_meme.branches.has_background") {
        denominator = testable + 1;
    }
    
    if (^"c_meme.srv") {
        ^(^"c_meme.alpha.scaler") = Min (10, (^(^"c_meme.alpha.scaler")+denominator*sum)/denominator);
    }

    ref_group  = (^"c_meme.branches.testable")["VALUEINDEXORDER"][0];
    ref_params = (^"c_meme.scaler_parameter_names")[ref_group];
    
    df = utility.Array1D (ref_params);
   
    for (p_type, ref_parameter; in; ref_params) {
            
        sum = 0;     
        for (g,p; in; ^"c_meme.scaler_parameter_names") {
            sum += ^(p [p_type]);
        }    
        
        
        // baseline NULL (everything = same rate)


        ^ref_parameter = sum /denominator;
   
     
        if (testable == 1) {
            parameters.SetConstraint (((^"c_meme.scaler_parameter_names")[^"terms.tree_attributes.background"])[p_type],ref_parameter, "");
        } else {
            for (_gname_; in; ^"c_meme.branches.testable") {
                _pname_ =  ((^"c_meme.scaler_parameter_names")[_gname_])[p_type];
                if (_pname_ != ref_parameter && (ignorable)[_pname_] == 0) {
                     //console.log (_pname_ + "=" + ref_parameter );
                     parameters.SetConstraint (_pname_,ref_parameter, "");
                }
            }
        }
    }


    Optimize (results, ^lf, {"OPTIMIZATION_METHOD" : "nedler-mead"});
    Null = estimators.ExtractMLEs (lf, model_mapping);
    Null [utility.getGlobalValue("terms.fit.log_likelihood")] = results[1][0];

    if (testable > 2) {
        pairwise = {};
        for (v = 0; v < testable; v+=1) {
            for (v2 = v + 1; v2 < testable; v2+=1) {
                v1n = (^"c_meme.branches.testable")[v];
                v2n = (^"c_meme.branches.testable")[v2];
                 
                estimators.RestoreLFStateFromSnapshot (lf_id, snapshot);
                
                ignorable_pair = FALSE;
                
                for (p,vr; in; (^"c_meme.scaler_parameter_names")[v1n]) {
                    if ((ignorable)[vr] || (ignorable)[((^"c_meme.scaler_parameter_names")[v2n])[p]]) {
                        ignorable_pair = TRUE;
                        break;
                    }
                    parameters.SetConstraint (vr,((^"c_meme.scaler_parameter_names")[v2n])[p], "");
                }
                
                if (ignorable_pair) {
                    (pairwise[v1n + "|" + v2n]) = alternative;
                } else {                
                    Optimize (results, ^lf, {"OPTIMIZATION_METHOD" : "nedler-mead"});
                    pairwise[v1n + "|" + v2n] = estimators.ExtractMLEs (lf, model_mapping);
                    (pairwise[v1n + "|" + v2n])[utility.getGlobalValue("terms.fit.log_likelihood")] = results[1][0];
                }
           }
         }
    } else {
        pairwise = None;
    }
    
    lrt = math.DoLRT (Null[utility.getGlobalValue("terms.fit.log_likelihood")],
                      alternative[utility.getGlobalValue("terms.fit.log_likelihood")],
                      Max (3, 3*(^'c_meme.branch_class_counter' - 1)));
                      
    //console.log (lrt);

    p_values = {};
    p_values ["overall"] = lrt[^'terms.p_value'];
    test_keys = {};
    
    //estimators.RestoreLFStateFromSnapshot (lf_id, snapshot);


    if (^'c_meme.report.test_count' > 1) {

        for (v = 0; v < ^'c_meme.branch_class_counter'; v += 1) {
            for (v2 = v + 1; v2 < ^'c_meme.branch_class_counter'; v2 += 1) {
                test_key = (^'c_meme.branches.testable') [v] + "|" + (^'c_meme.branches.testable') [v2];
                test_keys + test_key;
                lrt = math.DoLRT ((pairwise[test_key])[utility.getGlobalValue("terms.fit.log_likelihood")],
                                   alternative[utility.getGlobalValue("terms.fit.log_likelihood")],
                                  df);
                p_values [test_key] = lrt[^'terms.p_value'];
            }
        }
    }
    
    p_values = math.HolmBonferroniCorrection (p_values);
        
    if (permutation == FALSE && (Min (p_values,0))["value"] <= ^'c_meme.p_value') {
        result = {
                    utility.getGlobalValue("terms.alternative") : alternative,
                    utility.getGlobalValue("terms.Null"): Null,
                    utility.getGlobalValue("terms.c_meme.pairwise"): pairwise,
                    utility.getGlobalValue("terms.substitutions") : counts_by_branch_set,
                    utility.getGlobalValue("terms.p_value") : p_values,
                    utility.getGlobalValue("terms.c_meme.test_keys") : test_keys
                };
                
             
        if (^'c_meme.permutations' == TRUE) {
            p_min = (Min (p_values,0))["value"];
            //console.log ("Entering permutation mode with p = " + p_min);
            perm_p_values = {};
            for (rep = 0; rep < ^'c_meme.site.permutations'; rep+=1) {
                this_set = (c_meme.handle_a_site (lf, filter_data, partition_index, pattern_info, model_mapping, Random(branch_sets,0), parameter_names, 1))[utility.getGlobalValue("terms.p_value")];
                //console.log (this_set);
                perm_p_values + this_set;
                //console.log ("" + rep + " => " + (Min (this_set,0))["value"]);
                if ((Min (this_set,0))["value"] <= p_min) {
                    break;
                }
            }
            result [utility.getGlobalValue("terms.c_meme.permutation")]= perm_p_values;
        }
        return result;
    } 
    
    return {
                utility.getGlobalValue("terms.alternative") : alternative,
                utility.getGlobalValue("terms.Null"): Null,
                utility.getGlobalValue("terms.c_meme.pairwise"): pairwise,
                utility.getGlobalValue("terms.substitutions") : counts_by_branch_set,
                utility.getGlobalValue("terms.c_meme.test_keys") : test_keys,
                utility.getGlobalValue("terms.p_value") : p_values
           };
}


function c_meme.report.echo (c_meme.report.site, c_meme.report.partition, c_meme.report.row) {

    c_meme.k.bound   = Rows (c_meme.tests.key);

    //console.log (c_meme.report.row);
    //console.log (c_meme.tests.key);

    for (c_meme.k = 0; c_meme.k < c_meme.k.bound;  c_meme.k += 1) {
        c_meme.print_row = None;
        if (c_meme.report.row[c_meme.tests.key[c_meme.k][3]] <= c_meme.p_value) {
            c_meme.report.test         = c_meme.tests.key[c_meme.k][0];
            c_meme.report.rate1_index  = +c_meme.tests.key[c_meme.k][1];
            c_meme.report.pvalue_index = +c_meme.tests.key[c_meme.k][3];
            
            c_meme.report.permutation_index = Columns (c_meme.report.row) - 2;
            
            if (c_meme.report.rate1_index >= 0) {
                c_meme.report.rate2_index = +c_meme.tests.key[c_meme.k][2];
                c_meme.print_row = c_meme.report.pairwise;
                c_meme.report.substitutions = "" + c_meme.report.row[0 + ^"c_meme.report.rate_count" + (c_meme.report.rate1_index-1)$3] + ", " +
                                                   c_meme.report.row[0 + ^"c_meme.report.rate_count" + (c_meme.report.rate2_index-1)$3];
                
            } else {
                c_meme.print_row = c_meme.report.overall;
                c_meme.report.rates = Transpose (c_meme.report.row[{{0,1}}][{{0,c_meme.report.rate_count-1}}]) % 0;
                c_meme.report.substitutions = "" + c_meme.report.row[^"c_meme.report.rate_count"];
                for (c_meme.i = 1; c_meme.i < c_meme.report.subs_kinds; c_meme.i +=1) {
                     c_meme.report.substitutions += ", " + c_meme.report.row[^"c_meme.report.rate_count" + c_meme.i];
                }
                
            }
            c_meme.report.counts [c_meme.k] += 1;
        }

        if (None != c_meme.print_row) {
            //console.log (c_meme.print_row);
            //console.log (c_meme.table_output_options);
            io.ClearProgressBar();

            if (!c_meme.report.header_done) {
                io.ReportProgressMessageMD("c_meme", "" + c_meme.report.partition, "For partition " + (c_meme.report.partition+1) + " these sites are significant at p <=" + c_meme.p_value + "\n");
                fprintf (stdout,
                    io.FormatTableRow (c_meme.table_screen_output,c_meme.table_output_options));
                c_meme.report.header_done = TRUE;
                c_meme.table_output_options[terms.table_options.header] = FALSE;
            }
            fprintf (stdout,
                io.FormatTableRow (c_meme.print_row,c_meme.table_output_options));
        }
    }
}


lfunction c_meme.store_results (node, result, arguments) {



    //console.log (^"c_meme.table_headers");
        
    partition_index = arguments [2];
    pattern_info    = arguments [3];

    array_size  = utility.getGlobalValue ("c_meme.report.rate_count") + utility.getGlobalValue ("c_meme.report.subs_kinds") + utility.getGlobalValue ("c_meme.report.test_count") + 3;
    result_row  = { 1, array_size } ["_MATRIX_ELEMENT_COLUMN_>=^'c_meme.report.rate_count'&&_MATRIX_ELEMENT_COLUMN_<array_size-1"];

   
    if (None != result) { // not a constant site

        result_row [0] = estimators.GetGlobalMLE (result[utility.getGlobalValue("terms.alternative")], ^"c_meme.site_alpha");
        k = 1;
        
        for (g; in; ^"c_meme.rate.names") {
            for (r, p; in; (^"c_meme.scaler_parameter_names")[g]) {
                result_row[k] = estimators.GetGlobalMLE (result[utility.getGlobalValue("terms.alternative")], ((^"c_meme.site_tested_classes")[g])[r]);
                k += 1;
            }
        }
        
        
        for (rn; in; ^"c_meme.branches.testable") {
            result_row[k] = (result[utility.getGlobalValue("terms.substitutions")])[rn];
            k += 1;
        }
        
        
        p_values = result[^"terms.p_value"];
        result_row [k] = p_values ["overall"];
        k += 1;
        result_row [k] = 0; // q-value to be computed later
         
        test_keys = result[^"terms.c_meme.test_keys"];
        
        for (v = 1; v < ^'c_meme.report.test_count'; v+=1) {
            k += 1;
            result_row [k] = p_values[test_keys[v-1]];
        }
        
        k += 1;
        
        alternative_lengths = utility.Array1D(result[^"terms.c_meme.permutation"]);

        if (alternative_lengths) {
            result_row [k] = 1./alternative_lengths;
        } else {
            result_row [k] = -1;
        }

        k += 1;
        
        sum = 0;
    
        alternative_lengths = ((result[utility.getGlobalValue("terms.alternative")])[utility.getGlobalValue("terms.branch_length")])[0];

        for (_node_;in;^"c_meme.case_respecting_node_names") {
            sum += (alternative_lengths[_node_])[utility.getGlobalValue("terms.fit.MLE")];
        }

        result_row [k] = sum;
        
      }
      
      
    utility.EnsureKey (^"c_meme.site_results", partition_index);
    
    for (_fel_result_; in; pattern_info[utility.getGlobalValue("terms.data.sites")]) {
       ((^"c_meme.site_results")[partition_index])[_fel_result_] = result_row;
       c_meme.report.echo (_fel_result_, partition_index, result_row);
    }
    
    //console.log (result_row);

    //assert (0);
}
