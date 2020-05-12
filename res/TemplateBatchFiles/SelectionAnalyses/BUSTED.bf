RequireVersion ("2.5.12");


LoadFunctionLibrary("libv3/all-terms.bf"); // must be loaded before CF3x4
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/models/codon.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/tasks/genetic_code.bf");
LoadFunctionLibrary("modules/io_functions.ibf");
LoadFunctionLibrary("modules/selection_lib.ibf");
LoadFunctionLibrary("libv3/models/codon/BS_REL.bf");
LoadFunctionLibrary("libv3/convenience/math.bf");
LoadFunctionLibrary("libv3/models/rate_variation.bf");


utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);
utility.SetEnvVariable ("LF_SMOOTHING_SCALER", 0.1);


busted.analysis_description = {
                               terms.io.info : 
"BUSTED (branch-site unrestricted statistical test of episodic diversification) uses a random effects branch-site model fitted 
jointly to all or a subset of tree branches in order to test for alignment-wide evidence of episodic diversifying selection. 
Assuming there is evidence of positive selection (i.e. there is an omega > 1),  BUSTED will also perform a quick evidence-ratio 
style analysis to explore which individual sites may have been subject to selection. v2.0 adds support for synonymous rate variation, 
and relaxes the test statistic to 0.5 (chi^2_0 + chi^2_2). Version 2.1 adds a grid search for the initial starting point.
Version 2.2 changes the grid search to LHC, and adds an initial search phase to use adaptive Nedler-Mead. Version 3.0 implements the option
for branch-site variation in synonymous substitution rates. Version 3.1 adds HMM auto-correlation option for SRV, and binds SRV distributions for multiple branch sets
",
                               terms.io.version : "3.1",
                               terms.io.reference : "*Gene-wide identification of episodic selection*, Mol Biol Evol. 32(5):1365-71",
                               terms.io.authors : "Sergei L Kosakovsky Pond",
                               terms.io.contact : "spond@temple.edu",
                               terms.io.requirements : "in-frame codon alignment and a phylogenetic tree (optionally annotated with {})"
                              };

io.DisplayAnalysisBanner (busted.analysis_description);

busted.FG = "Test";
busted.BG = "Background";
busted.SRV = "Synonymous site-to-site rates";
busted.background = "background";
busted.unconstrained = "unconstrained";
busted.constrained = "constrained";
busted.optimized_null = "optimized null";
busted.MG94 = terms.json.mg94xrev_sep_rates;

busted.json.background = busted.background;
busted.json.site_logl  = "Site Log Likelihood";
busted.json.evidence_ratios  = "Evidence Ratios";
busted.json.srv_posteriors  = "Synonymous site-posteriors";
busted.json.srv_viterbi = "Viterbi synonymous rate path";
busted.rate_classes = 3;
busted.synonymous_rate_classes = 3;
busted.initial_grid.N = 250;


busted.json    = { terms.json.analysis: busted.analysis_description,
                   terms.json.input: {},
                   busted.json.background: {},
                   terms.json.fits : {},
                   terms.json.timers : {},
                   busted.json.site_logl : {},
                   busted.json.evidence_ratios: {},
                   busted.json.site_logl : {}
                  };


busted.display_orders = {terms.original_name: -1,
                         terms.json.nucleotide_gtr: 0,
                         busted.MG94: 1,
                         busted.unconstrained: 2,
                         busted.constrained: 3                        
                        };


selection.io.startTimer (busted.json [terms.json.timers], "Overall", 0);

KeywordArgument ("code",      "Which genetic code should be used", "Universal");
    /**
        keyword, description (for inline documentation and help messages), default value
    */
KeywordArgument ("alignment", "An in-frame codon alignment in one of the formats supported by HyPhy");
    /**
        keyword, description (for inline documentation and help messages), no default value,
        meaning that it will be required
    */

KeywordArgument ("tree",      "A phylogenetic tree (optionally annotated with {})", null, "Please select a tree file for the data:");
    /** the use of null as the default argument means that the default expectation is for the 
        argument to be missing, i.e. the tree is expected to be in the file
        the fourth, optional argument, can match this keyword with the dialog prompt / choice list title,
        meaning that it can only be consumed when this dialog prompt / choice list is invoked
        This allows handling some branching logic conditionals
    */
    
KeywordArgument ("branches",  "Branches to test", "All");
KeywordArgument ("srv", "Include synonymous rate variation in the model", "Yes");
KeywordArgument ("rates", "The number omega rate classes to include in the model [1-10, default 3]", busted.rate_classes);

namespace busted {
    LoadFunctionLibrary ("modules/shared-load-file.bf");
    load_file ("busted");
}


busted.do_srv = io.SelectAnOption ({"Yes" : "Allow synonymous substitution rates to vary from site to site (but not from branch to branch)", 
                                    "Branch-site" : "Allow synonymous substitution rates to vary using general branch site models",
                                    "HMM" : "Allow synonymous substitution rates to vary from site to site (but not from branch to branch), and use a hidden markov model for spatial correlation on synonymous rates",
                                    "No"  : "Synonymous substitution rates are constant across sites. This is the 'classic' behavior, i.e. the original published test"},
                                    "Synonymous rate variation"
                                    );
                                    
if (busted.do_srv == "Branch-site") {
    busted.do_bs_srv = TRUE;
    busted.do_srv = TRUE;
    (busted.json[terms.json.analysis])[terms.settings] = "Branch-site synonymous rate variation";
} else {
    if (busted.do_srv == "Yes") {
        busted.do_bs_srv = FALSE;
        busted.do_srv = TRUE;
        busted.do_srv_hmm  = FALSE; 
    } else {
        if (busted.do_srv == "HMM") {
         busted.do_srv      = TRUE;
         busted.do_bs_srv   = FALSE;
         busted.do_srv_hmm  = TRUE; 
        } else {
            busted.do_srv = FALSE;  
        } 
    }
}                       
                                    

busted.rate_classes = io.PromptUser ("The number omega rate classes to include in the model", busted.rate_classes, 1, 10, TRUE);

if (busted.do_srv) {
    KeywordArgument ("syn-rates", "The number alpha rate classes to include in the model [1-10, default 3]", busted.synonymous_rate_classes);
    busted.synonymous_rate_classes = io.PromptUser ("The number omega rate classes to include in the model", busted.synonymous_rate_classes, 1, 10, TRUE);
}

KeywordArgument ("grid-size", "The number of points in the initial distributional guess for likelihood fitting", 250);
busted.initial_grid.N = io.PromptUser ("The number of points in the initial distributional guess for likelihood fitting", 250, 1, 10000, TRUE);
KeywordArgument ("starting-points", "The number of initial random guesses to seed rate values optimization", 1);
busted.N.initial_guesses = io.PromptUser ("The number of initial random guesses to 'seed' rate values optimization", 1, 1, busted.initial_grid.N$10, TRUE);
                                    
KeywordArgument ("output", "Write the resulting JSON to this file (default is to save to the same path as the alignment file + 'BUSTED.json')", busted.codon_data_info [terms.json.json]);
busted.codon_data_info [terms.json.json] = io.PromptUserForFilePath ("Save the resulting JSON file to");
    
io.ReportProgressMessageMD('BUSTED',  'selector', 'Branches to test for selection in the BUSTED analysis');

utility.ForEachPair (busted.selected_branches, "_partition_", "_selection_",
    "_selection_ = utility.Filter (_selection_, '_value_', '_value_ == terms.tree_attributes.test');
     io.ReportProgressMessageMD('BUSTED',  'selector', '* Selected ' + Abs(_selection_) + ' branches to test in the BUSTED analysis: \\\`' + Join (', ',utility.Keys(_selection_)) + '\\\`')");

// check to see if there are any background branches
busted.has_background = FALSE;

utility.ForEachPair (busted.selected_branches, "_partition_", "_selection_",
    "_selection_ = utility.Filter (_selection_, '_value_', '_value_ != terms.tree_attributes.test');
     if (utility.Array1D (_selection_)) { busted.has_background = TRUE;} ");
busted.json[busted.json.background] =  busted.has_background;

selection.io.startTimer (busted.json [terms.json.timers], "Preliminary model fitting", 1);

namespace busted {
    doGTR ("busted");
}


estimators.fixSubsetOfEstimates(busted.gtr_results, busted.gtr_results[terms.global]);

namespace busted {
    scaler_prefix = "busted.scaler";
    doPartitionedMG ("busted", FALSE);
}


io.ReportProgressMessageMD ("BUSTED", "codon-refit", "Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model");

busted.final_partitioned_mg_results = estimators.FitMGREV (busted.filter_names, busted.trees, busted.codon_data_info [terms.code], {
    terms.run_options.model_type: terms.local,
    terms.run_options.partitioned_omega: busted.selected_branches,
}, busted.partitioned_mg_results);


io.ReportProgressMessageMD("BUSTED", "codon-refit", "* " + selection.io.report_fit (busted.final_partitioned_mg_results, 0, busted.codon_data_info[terms.data.sample_size]));
busted.global_dnds = selection.io.extract_global_MLE_re (busted.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio);

utility.ForEach (busted.global_dnds, "_value_", 'io.ReportProgressMessageMD ("BUSTED", "codon-refit", "* " + _value_[terms.description] + " = " + Format (_value_[terms.fit.MLE],8,4));');

selection.io.stopTimer (busted.json [terms.json.timers], "Preliminary model fitting");

//Store MG94 to JSON
selection.io.json_store_lf_withEFV (busted.json,
                            busted.MG94,
                            busted.final_partitioned_mg_results[terms.fit.log_likelihood],
                            busted.final_partitioned_mg_results[terms.parameters],
                            busted.sample_size,
                            utility.ArrayToDict (utility.Map (busted.global_dnds, "_value_", "{'key': _value_[terms.description], 'value' : Eval({{_value_ [terms.fit.MLE],1}})}")),
                            (busted.final_partitioned_mg_results[terms.efv_estimate])["VALUEINDEXORDER"][0],
                            busted.display_orders[busted.MG94]);
                            
utility.ForEachPair (busted.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(busted.json, busted.MG94, terms.branch_length, busted.display_orders[busted.MG94],
                                             _key_,
                                             selection.io.extract_branch_info((busted.final_partitioned_mg_results[terms.branch_length])[_key_], "selection.io.branch.length"));');


busted.model_generator = "models.codon.BS_REL.ModelDescription";


if (busted.do_srv) {

    if (busted.do_bs_srv) {
        busted.model_generator = "busted.model.with.GDD";
        busted.model_generator = "models.codon.BS_REL_SRV.ModelDescription";
        busted.rate_class_arguments = {{busted.synonymous_rate_classes__,busted.rate_classes__}};
    } else {
    
        lfunction busted.model.with.GDD (type, code, rates) {        
            def = models.codon.BS_REL.ModelDescription (type, code, rates);
            options = {utility.getGlobalValue("terms.rate_variation.bins") : utility.getGlobalValue("busted.synonymous_rate_classes"),
                       utility.getGlobalValue("terms._namespace") : "busted._shared_srv"};

            if (utility.getGlobalValue("busted.do_srv_hmm")) {
                    options [utility.getGlobalValue ("terms.rate_variation.HMM")] = "equal";
            }
            def [utility.getGlobalValue("terms.model.rate_variation")] = rate_variation.types.GDD.factory (options);
            return def;
        }
        
        busted.model_generator = "busted.model.with.GDD";
        busted.rate_class_arguments = busted.rate_classes;
    }
} else {
    busted.rate_class_arguments = busted.rate_classes;
}

busted.test.bsrel_model =  model.generic.DefineMixtureModel(busted.model_generator,
        "busted.test", {
            "0": parameters.Quote(terms.global),
            "1": busted.codon_data_info[terms.code],
            "2": parameters.Quote (busted.rate_class_arguments) // the number of rate classes
        },
        busted.filter_names,
        None);
        
busted.background.bsrel_model =  model.generic.DefineMixtureModel(busted.model_generator,
        "busted.background", {
            "0": parameters.Quote(terms.global),
            "1": busted.codon_data_info[terms.code],
            "2": parameters.Quote (busted.rate_class_arguments) // the number of rate classes
        },
        busted.filter_names,
        None);


models.BindGlobalParameters ({"0" : busted.test.bsrel_model, "1" : busted.background.bsrel_model}, terms.nucleotideRate("[ACGT]","[ACGT]"));


/*if (busted.do_srv) {
    if (busted.do_bs_srv) {
        models.BindGlobalParameters ({"0" : busted.test.bsrel_model, "1" : busted.background.bsrel_model}, "Mean scaler variable for");
        models.BindGlobalParameters ({"0" : busted.test.bsrel_model, "1" : busted.background.bsrel_model}, terms.AddCategory (utility.getGlobalValue("terms.mixture.mixture_aux_weight"), "SRV [0-9]+"));
    } else {
        models.BindGlobalParameters ({"0" : busted.test.bsrel_model, "1" : busted.background.bsrel_model}, "GDD rate category");
        models.BindGlobalParameters ({"0" : busted.test.bsrel_model, "1" : busted.background.bsrel_model}, utility.getGlobalValue("terms.mixture.mixture_aux_weight") + " for GDD category ");
    }
}*/

busted.distribution = models.codon.BS_REL.ExtractMixtureDistribution(busted.test.bsrel_model);

// set up parameter constraints

for (busted.i = 1; busted.i < busted.rate_classes; busted.i += 1) {
    parameters.SetRange (model.generic.GetGlobalParameter (busted.test.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,busted.i)), terms.range01);
    parameters.SetRange (model.generic.GetGlobalParameter (busted.background.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,busted.i)), terms.range01);
}


parameters.SetRange (model.generic.GetGlobalParameter (busted.test.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,busted.rate_classes)), terms.range_gte1);

/* create an initial grid to initialize the optimization */
busted.initial_grid = {};

/* if populated, use this as a baseline to generate the distributions from */
busted.initial_grid_presets = {"0" : 0.1};
busted.initial_ranges = {};

busted.init_grid_setup (busted.distribution);

/** setup parameter optimization groups */

PARAMETER_GROUPING = {};
PARAMETER_GROUPING + busted.distribution["rates"];
PARAMETER_GROUPING + busted.distribution["weights"];


if (busted.has_background) {
    busted.model_object_map = { "busted.background" : busted.background.bsrel_model,
                                "busted.test" :       busted.test.bsrel_model };
    busted.background_distribution = models.codon.BS_REL.ExtractMixtureDistribution(busted.background.bsrel_model);
    busted.init_grid_setup (busted.background_distribution);

    PARAMETER_GROUPING = {};
    PARAMETER_GROUPING + busted.background_distribution["rates"];
    PARAMETER_GROUPING + busted.background_distribution["weights"];
    
} else {
    busted.model_object_map = { "busted.test" :       busted.test.bsrel_model };
}

if (busted.do_srv)  {

    if (busted.do_bs_srv) {
        busted.srv_rate_regex   = "Mean scaler variable for";
        busted.srv_weight_regex = terms.AddCategory (utility.getGlobalValue("terms.mixture.mixture_aux_weight"), "SRV [0-9]+");
        
        busted.srv_rate_reporting = regexp.PartitionByRegularExpressions(utility.Keys ((busted.test.bsrel_model[terms.parameters])[terms.global]), 
            {"0" : "^" + utility.getGlobalValue('terms.parameters.synonymous_rate'), 
             "1" :  busted.srv_weight_regex});

        busted.srv_rate_reporting = {
          'rates' : utility.UniqueValues (utility.Map ( busted.srv_rate_reporting ["^" + utility.getGlobalValue('terms.parameters.synonymous_rate') ]  , "_value_", '((busted.test.bsrel_model[terms.parameters])[terms.global])[_value_]')),
           'weights' : utility.UniqueValues (utility.Map (busted.srv_rate_reporting [busted.srv_weight_regex ]  , "_value_", '((busted.test.bsrel_model[terms.parameters])[terms.global])[_value_]'))
        };
    } else {
        busted.srv_rate_regex  = "GDD rate category [0-9]+";
        busted.srv_weight_regex = "Mixture auxiliary weight for GDD category [0-9]+";
    }
    busted.srv_distribution = regexp.PartitionByRegularExpressions(utility.Keys ((busted.test.bsrel_model[terms.parameters])[terms.global]), {"0" : busted.srv_rate_regex, "1" : busted.srv_weight_regex});

    
    busted.srv_distribution = {
        'rates' : utility.UniqueValues (utility.Map (busted.srv_distribution [busted.srv_rate_regex ]  , "_value_", '((busted.test.bsrel_model[terms.parameters])[terms.global])[_value_]')),
        'weights' : utility.UniqueValues (utility.Map (busted.srv_distribution [busted.srv_weight_regex ]  , "_value_", '((busted.test.bsrel_model[terms.parameters])[terms.global])[_value_]'))
    };
    
 
    PARAMETER_GROUPING + busted.srv_distribution["rates"];
    PARAMETER_GROUPING + busted.srv_distribution["weights"];

    busted.init_grid_setup (busted.srv_distribution);
    
}


busted.initial.test_mean    = ((selection.io.extract_global_MLE_re (busted.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio + ".+test.+"))["0"])[terms.fit.MLE];
busted.initial_grid         = estimators.LHC (busted.initial_ranges,busted.initial_grid.N);

//estimators.CreateInitialGrid (busted.initial_grid, busted.initial_grid.N, busted.initial_grid_presets);

busted.initial_grid = utility.Map (busted.initial_grid, "_v_", 
    'busted._renormalize (_v_, "busted.distribution", busted.initial.test_mean)'
);

if (busted.has_background) { //GDD rate category
    busted.initial.background_mean    = ((selection.io.extract_global_MLE_re (busted.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio + ".+background.+"))["0"])[terms.fit.MLE];
    busted.initial_grid = utility.Map (busted.initial_grid, "_v_", 
        'busted._renormalize (_v_, "busted.background_distribution", busted.initial.background_mean)'
    );
}


busted.model_map = {};

for (busted.partition_index = 0; busted.partition_index < busted.partition_count; busted.partition_index += 1) {
    selection.io.json_store_branch_attribute(busted.json, terms.original_name, terms.json.node_label, 0,
                                             busted.partition_index,
                                             busted.name_mapping);

    busted.model_map + { "busted.test" : utility.Filter (busted.selected_branches[busted.partition_index], '_value_', '_value_ == terms.tree_attributes.test'),
					     "busted.background" : utility.Filter (busted.selected_branches[busted.partition_index], '_value_', '_value_ != terms.tree_attributes.test')};
}

utility.SetEnvVariable ("ASSUME_REVERSIBLE_MODELS", TRUE);

selection.io.startTimer (busted.json [terms.json.timers], "Unconstrained BUSTED model fitting", 2);

io.ReportProgressMessageMD ("BUSTED", "main", "Performing the full (dN/dS > 1 allowed) branch-site model fit");

/** 
    perform the initial fit using a constrained branch length model and 
    a low precision Nedler-Mead algorithm pass to find good initial values 
    for the rate distribution parameters
*/

  
busted.nm.precision = -0.00025*busted.final_partitioned_mg_results[terms.fit.log_likelihood];

parameters.DeclareGlobalWithRanges ("busted.bl.scaler", 1, 0, 1000);
busted.grid_search.results =  estimators.FitLF (busted.filter_names, busted.trees, busted.model_map, busted.final_partitioned_mg_results, busted.model_object_map, {
    "retain-lf-object": TRUE,
    terms.run_options.proportional_branch_length_scaler : 
                                            {"0" : "busted.bl.scaler"},
                                            
    terms.run_options.optimization_settings : 
        {
            "OPTIMIZATION_METHOD" : "nedler-mead",
            "MAXIMUM_OPTIMIZATION_ITERATIONS" : 500,
            "OPTIMIZATION_PRECISION" : busted.nm.precision
        } ,
                                     
    terms.search_grid     : busted.initial_grid,
    terms.search_restarts : busted.N.initial_guesses
});
                                
                
busted.full_model =  estimators.FitLF (busted.filter_names, busted.trees, busted.model_map, busted.grid_search.results, busted.model_object_map, {
        "retain-lf-object": TRUE,
        terms.run_options.optimization_settings : 
            {
                "OPTIMIZATION_METHOD" : "hybrid"
            } 
                                    
    });
    


KeywordArgument ("save-fit", "Save BUSTED model fit to this file (default is not to save)", "/dev/null");
io.SpoolLFToPath(busted.full_model[terms.likelihood_function], io.PromptUserForFilePath ("Save BUSTED model fit to this file ['/dev/null' to skip]"));

io.ReportProgressMessageMD("BUSTED", "main", "* " + selection.io.report_fit (busted.full_model, 9, busted.codon_data_info[terms.data.sample_size]));

if (busted.do_srv_hmm) {
    busted.hmm_lambda = selection.io.extract_global_MLE (busted.full_model, terms.rate_variation.hmm_lambda);
    io.ReportProgressMessageMD("BUSTED", "main", "* HMM switching rate = " +  Format (busted.hmm_lambda, 8, 3));

}


io.ReportProgressMessageMD("BUSTED", "main", "* For *test* branches, the following rate distribution for branch-site combinations was inferred");

selection.io.stopTimer (busted.json [terms.json.timers], "Unconstrained BUSTED model fitting");

busted.inferred_test_distribution = parameters.GetStickBreakingDistribution (busted.distribution) % 0;
selection.io.report_dnds (busted.inferred_test_distribution);



busted.distribution_for_json = {busted.FG : utility.Map (utility.Range (busted.rate_classes, 0, 1),
                                                         "_index_",
                                                         "{terms.json.omega_ratio : busted.inferred_test_distribution [_index_][0],
                                                           terms.json.proportion : busted.inferred_test_distribution [_index_][1]}")
                                };

if (busted.do_srv_hmm) {
    busted.distribution_for_json [terms.rate_variation.hmm_lambda] = busted.hmm_lambda;
}

if (busted.has_background) {
    io.ReportProgressMessageMD("BUSTED", "main", "* For *background* branches, the following rate distribution for branch-site combinations was inferred");
    busted.inferred_background_distribution = parameters.GetStickBreakingDistribution (busted.background_distribution) % 0;
    selection.io.report_dnds (busted.inferred_background_distribution);
    busted.distribution_for_json [busted.BG] = utility.Map (utility.Range (busted.rate_classes, 0, 1),
                                                         "_index_",
                                                         "{terms.json.omega_ratio : busted.inferred_background_distribution [_index_][0],
                                                           terms.json.proportion : busted.inferred_background_distribution [_index_][1]}");
}

if (busted.do_srv) {
    
    if (busted.do_bs_srv) {
        busted.srv_info = parameters.GetStickBreakingDistribution ( busted.srv_rate_reporting) % 0
    } else {
        busted.srv_info = Transpose((rate_variation.extract_category_information(busted.test.bsrel_model))["VALUEINDEXORDER"][0])%0;
    }
    io.ReportProgressMessageMD("BUSTED", "main", "* The following rate distribution for site-to-site **synonymous** rate variation was inferred");
    selection.io.report_distribution (busted.srv_info);

    busted.distribution_for_json [busted.SRV] = (utility.Map (utility.Range (busted.synonymous_rate_classes, 0, 1),
                                                             "_index_",
                                                             "{terms.json.rate :busted.srv_info [_index_][0],
                                                               terms.json.proportion : busted.srv_info [_index_][1]}"));
                                                               
    ConstructCategoryMatrix (busted.cmx, ^(busted.full_model[terms.likelihood_function]));
    ConstructCategoryMatrix (busted.cmx_weights, ^(busted.full_model[terms.likelihood_function]), WEIGHTS);
    busted.cmx_weighted         = (busted.cmx_weights[-1][0]) $ busted.cmx; // taking the 1st column fixes a bug with multiple partitions 
    busted.column_weights       = {1, Rows (busted.cmx_weights)}["1"] * busted.cmx_weighted;
    busted.column_weights       = busted.column_weights["1/_MATRIX_ELEMENT_VALUE_"];
    (busted.json [busted.json.srv_posteriors]) =  busted.cmx_weighted $ busted.column_weights;
    
    
    if (busted.do_srv_hmm ) {
        ConstructCategoryMatrix (busted.cmx_viterbi, ^(busted.full_model[terms.likelihood_function]), SHORT);
        (busted.json [busted.json.srv_viterbi]) = busted.cmx_viterbi;
        io.ReportProgressMessageMD("BUSTED", "main", "* The following switch points for synonymous rates were inferred");
        selection.io.report_viterbi_path (busted.cmx_viterbi);
        
    }

}



selection.io.json_store_lf (busted.json,
                            "Unconstrained model",
                            busted.full_model[terms.fit.log_likelihood],
                            busted.full_model[terms.parameters] + 9 , // +9 comes from CF3x4
                            busted.codon_data_info[terms.data.sample_size],
                            busted.distribution_for_json,
                            busted.display_orders[busted.unconstrained]);


busted.run_test = busted.inferred_test_distribution [busted.rate_classes-1][0] > 1 && busted.inferred_test_distribution [busted.rate_classes-1][1] > 0;

utility.ForEachPair (busted.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(busted.json, busted.unconstrained, terms.branch_length, 0,
                                             _key_,
                                             selection.io.extract_branch_info((busted.full_model[terms.branch_length])[_key_], "selection.io.branch.length"));');


(busted.json [busted.json.site_logl])[busted.unconstrained] = busted.ComputeSiteLikelihoods (busted.full_model[terms.likelihood_function]);

if (!busted.run_test) {
    io.ReportProgressMessageMD ("BUSTED", "Results", "No evidence for episodic diversifying positive selection under the unconstrained model, skipping constrained model fitting");
    busted.json [terms.json.test_results] = busted.ComputeLRT (0, 0);
} else {
    selection.io.startTimer (busted.json [terms.json.timers], "Constrained BUSTED model fitting", 3);


    io.ReportProgressMessageMD ("BUSTED", "test", "Performing the constrained (dN/dS > 1 not allowed) model fit");
    parameters.SetConstraint (model.generic.GetGlobalParameter (busted.test.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,busted.rate_classes)), terms.parameters.one, terms.global);
    (busted.json [busted.json.site_logl])[busted.constrained] = busted.ComputeSiteLikelihoods (busted.full_model[terms.likelihood_function]);
    busted.null_results = estimators.FitExistingLF (busted.full_model[terms.likelihood_function], busted.model_object_map);
    (busted.json [busted.json.site_logl])[busted.optimized_null] = busted.ComputeSiteLikelihoods (busted.full_model[terms.likelihood_function]);
    io.ReportProgressMessageMD ("BUSTED", "test", "* " + selection.io.report_fit (busted.null_results, 9, busted.codon_data_info[terms.data.sample_size]));
    busted.LRT = busted.ComputeLRT (busted.full_model[terms.fit.log_likelihood], busted.null_results[terms.fit.log_likelihood]);
    busted.json [terms.json.test_results] = busted.LRT;


    selection.io.stopTimer (busted.json [terms.json.timers], "Constrained BUSTED model fitting");

    utility.ForEachPair (busted.filter_specification, "_key_", "_value_",
        'selection.io.json_store_branch_attribute(busted.json, busted.constrained, terms.branch_length, 1,
                                                 _key_,
                                                 selection.io.extract_branch_info((busted.null_results[terms.branch_length])[_key_], "selection.io.branch.length"));');

    io.ReportProgressMessageMD("BUSTED", "test", "* For *test* branches under the null (no dN/dS > 1 model), the following rate distribution for branch-site combinations was inferred");
    busted.inferred_test_distribution = parameters.GetStickBreakingDistribution (busted.distribution) % 0;
    selection.io.report_dnds (parameters.GetStickBreakingDistribution (busted.distribution) % 0);

    busted.distribution_for_json = {busted.FG : utility.Map (utility.Range (busted.rate_classes, 0, 1),
                                                         "_index_",
                                                         "{terms.json.omega_ratio : busted.inferred_test_distribution [_index_][0],
                                                           terms.json.proportion : busted.inferred_test_distribution [_index_][1]}")};

    if (busted.has_background) {
        busted.inferred_background_distribution = parameters.GetStickBreakingDistribution (busted.background_distribution) % 0;
        //selection.io.report_dnds (busted.inferred_background_distribution);
        busted.distribution_for_json [busted.BG] = utility.Map (utility.Range (busted.rate_classes, 0, 1),
                                                             "_index_",
                                                             "{terms.json.omega_ratio : busted.inferred_background_distribution [_index_][0],
                                                               terms.json.proportion : busted.inferred_background_distribution [_index_][1]}");
    }

    if (busted.do_srv) {
        if (busted.do_bs_srv) {
            busted.srv_info = parameters.GetStickBreakingDistribution ( busted.srv_rate_reporting) % 0
        } else {
            busted.srv_info = Transpose((rate_variation.extract_category_information(busted.test.bsrel_model))["VALUEINDEXORDER"][0])%0;
        }
        io.ReportProgressMessageMD("BUSTED", "main", "* The following rate distribution for site-to-site **synonymous** rate variation was inferred");
        selection.io.report_distribution (busted.srv_info);

         busted.distribution_for_json [busted.SRV] = (utility.Map (utility.Range (busted.synonymous_rate_classes, 0, 1),
                                                                 "_index_",
                                                                 "{terms.json.rate :busted.srv_info [_index_][0],
                                                                   terms.json.proportion : busted.srv_info [_index_][1]}"));

    }
    
    selection.io.json_store_lf (busted.json,
                            "Constrained model",
                            busted.null_results[terms.fit.log_likelihood],
                            busted.null_results[terms.parameters] + 9 , // +9 comes from CF3x4
                            busted.codon_data_info[terms.data.sample_size],
                            busted.distribution_for_json,
                            busted.display_orders[busted.constrained]);

    (busted.json [busted.json.evidence_ratios])[busted.constrained] = busted.EvidenceRatios ( (busted.json [busted.json.site_logl])[busted.unconstrained],  (busted.json [busted.json.site_logl])[busted.constrained]);
    (busted.json [busted.json.evidence_ratios ])[busted.optimized_null] = busted.EvidenceRatios ( (busted.json [busted.json.site_logl])[busted.unconstrained],  (busted.json [busted.json.site_logl])[busted.optimized_null]);
}


console.log ("----\n## Branch-site unrestricted statistical test of episodic diversification [BUSTED]");
console.log ( "Likelihood ratio test for episodic diversifying positive selection, **p = " + Format ((busted.json [terms.json.test_results])[terms.p_value], 8, 4) + "**.");

selection.io.stopTimer (busted.json [terms.json.timers], "Overall");

io.SpoolJSON (busted.json, busted.codon_data_info [terms.json.json]);

return busted.json;

/// HELPERS

//------------------------------------------------------------------------------
lfunction busted.ComputeSiteLikelihoods (id) {
    ConstructCategoryMatrix (sl, ^id, SITE_LOG_LIKELIHOODS);
   return sl;
}
//------------------------------------------------------------------------------

function busted.ComputeLRT (ha, h0) {
    return {terms.LRT : 2*(ha-h0),
            terms.p_value : 0.5*(1-CChi2 (2*(ha-h0),2))};
}

//------------------------------------------------------------------------------
lfunction busted.EvidenceRatios (ha, h0) {
    return ha["Exp(_MATRIX_ELEMENT_VALUE_-h0[_MATRIX_ELEMENT_COLUMN_])"];
}
//------------------------------------------------------------------------------

lfunction busted.DistributionGuess (mean) {
    /*guess = {{Random (0,0.25),Random (0.5,0.7)}
             {Random (0.25,0.6),Random (0.1,0.2)}
             {Random (2,20),Random (0.01, 0.1)}};
    */

    guess = {{0,0.7}{0.25,0.2}{10,0.1}};

    norm = + guess[-1][1];
    guess_mean = 1/(+(guess [-1][0] $ guess [-1][1]))/norm;
    return guess["_MATRIX_ELEMENT_VALUE_*(guess_mean*(_MATRIX_ELEMENT_COLUMN_==0)+(_MATRIX_ELEMENT_COLUMN_==1)*(1/norm))"];
}

//------------------------------------------------------------------------------

// renormalize the grid to have the same mean as the initial omega value


lfunction busted._renormalize (v, distro, mean) {
    parameters.SetValues (v);
    m = parameters.GetStickBreakingDistribution (^distro);
    d = Rows (m);
    m = +(m[-1][0] $ m[-1][1]); // current mean
    for (i = 0; i < d; i+=1) {
        (v[((^"busted.distribution")["rates"])[i]])[^"terms.fit.MLE"] = (v[((^"busted.distribution")["rates"])[i]])[^"terms.fit.MLE"] / m * mean;
    }
    return v;
    
}

//------------------------------------------------------------------------------

function busted.init_grid_setup (omega_distro) {
    utility.ForEachPair (omega_distro[terms.parameters.rates], "_index_", "_name_", 
        '
            if (_index_[0] < busted.rate_classes - 1) { // not the last rate
                busted.initial_grid  [_name_] = {
                    {
                        0.01, 0.1, 0.25, 0.75
                    }
                }["_MATRIX_ELEMENT_VALUE_^(busted.rate_classes-_index_[0]-1)"];
                busted.initial_grid_presets [_name_] = 0;
                busted.initial_ranges [_name_] = {
                    terms.lower_bound : 0,
                    terms.upper_bound : 1
                };
            }  else {
                busted.initial_grid  [_name_] = {
                    {
                        1, 1.5, 2, 4, 10
                    }
                };
                busted.initial_ranges [_name_] = {
                    terms.lower_bound : 1,
                    terms.upper_bound : 10
                };
                busted.initial_grid_presets [_name_] = 2;
            }
        '
    );


    utility.ForEachPair (omega_distro[terms.parameters.weights], "_index_", "_name_", 
        '
            busted.initial_grid  [_name_] = {
                {
                    0.2, 0.5, 0.7, 0.8, 0.9
                }
            };
            busted.initial_grid_presets [_name_] = 3;
            busted.initial_ranges [_name_] = {
                terms.lower_bound : 0,
                terms.upper_bound : 1
            };
        '
    );

}
//------------------------------------------------------------------------------

