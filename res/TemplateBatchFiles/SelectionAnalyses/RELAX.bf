RequireVersion ("2.5.82");

LoadFunctionLibrary("libv3/all-terms.bf"); // must be loaded before CF3x4

// namespace 'utility' for convenience functions
LoadFunctionLibrary("libv3/UtilityFunctions.bf");

// namespace 'io' for interactive/datamonkey i/o functions
LoadFunctionLibrary("libv3/IOFunctions.bf");

LoadFunctionLibrary ("libv3/models/codon/MG_REV.bf");

// namespace 'estimators' for various estimator related functions
LoadFunctionLibrary("libv3/tasks/estimators.bf");

// namespace 'estimators' for various estimator related functions
LoadFunctionLibrary("libv3/tasks/alignments.bf");

// namespace 'estimators' for various estimator related functions
LoadFunctionLibrary("libv3/tasks/trees.bf");

LoadFunctionLibrary("modules/io_functions.ibf");
LoadFunctionLibrary("modules/selection_lib.ibf");
LoadFunctionLibrary("libv3/models/codon/BS_REL.bf");
LoadFunctionLibrary("libv3/convenience/math.bf");
LoadFunctionLibrary("libv3/tasks/mpi.bf");
LoadFunctionLibrary("libv3/models/rate_variation.bf");



utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);
utility.SetEnvVariable ("ASSUME_REVERSIBLE_MODELS", TRUE);
utility.SetEnvVariable ("USE_MEMORY_SAVING_DATA_STRUCTURES", 1e8);
utility.SetEnvVariable ("ENFORCE_CONSTRAINT_VIOLATIONS", TRUE);

//relax.OPTIMIZATION_LOGS = 1;

/*------------------------------------------------------------------------------*/

relax.analysis_description = {
                               terms.io.info : "RELAX (a random effects test of selection relaxation) uses a random effects branch-site model framework to test whether a set of 'Test' branches evolves under relaxed selection relative to a set of 'Reference' branches (R), as measured by the relaxation parameter (K).
                                                Version 2.1 adds a check for stability in K estimates to try to mitigate convergence problems. 
                                                Version 3 provides support for >2 branch sets.
                                                Version 3.1 adds LHC + Nedler-Mead initial fit phase and keyword support.
                                                Version 3.1.1 adds some bug fixes for better convergence.
                                                Version 4.0 adds support for synonymous rate variation.
                                                Version 4.1 adds further support for multiple hit models.
                                                Version 4.1.1 adds reporting for convergence diagnostics.
                                                Version 4.5 adds support for multiple datasets for joint testing.
                                                Version 4.6 adds support for site-level EBF calculation.",
                               terms.io.version : "4.5",
                               terms.io.reference : "RELAX: Detecting Relaxed Selection in a Phylogenetic Framework (2015). Mol Biol Evol 32 (3): 820-832",
                               terms.io.authors : "Sergei L Kosakovsky Pond, Ben Murrell, Steven Weaver and Temple iGEM / UCSD viral evolution g",
                               terms.io.contact : "spond@temple.edu",
                               terms.io.requirements : "in-frame codon alignment and a phylogenetic tree, with at least two groups of branches defined using the {} notation (one group can be defined as all unlabeled branches)",
                               terms.settings: {}
                              };

relax.json    = { terms.json.analysis: relax.analysis_description,
                  terms.json.input: {},
                  terms.json.fits : {},
                  terms.json.timers : {},
                  terms.json.test_results : {},
                  };

relax.relaxation_parameter        = "relax.K";
relax.rate_classes                = 3;
relax.synonymous_rate_classes     = 3;
relax.SRV                         = "Synonymous site-to-site rates";
relax.json.srv_posteriors  = "Synonymous site-posteriors";
relax.json.srv_viterbi = "Viterbi synonymous rate path";
relax.kGroupMode                   = "Group mode";


relax.initial_ranges              = {};
relax.initial_grid.N              = 500;


relax.MG94_name = terms.json.mg94xrev_sep_rates;
relax.general_descriptive_name = "General descriptive";
relax.alternative_name = "RELAX alternative";
relax.null_name = "RELAX null";
relax.partitioned_descriptive_name = "RELAX partitioned descriptive";

relax.numbers_of_tested_groups = 2; 
/* 
   whether or not this is a group-set run, i.e. there are 3 or more sets of branches,
   then run the group test, i.e. K_1 = K_2 == .. K_{N-1} vs unconstrained alternative
*/


terms.relax.k          = "relaxation or intensification parameter";
terms.relax.k_range    = {
        terms.lower_bound: "0",
        terms.upper_bound: "50"
    };

terms.relax.k_range1    = {
        terms.lower_bound: "1",
        terms.upper_bound: "50"
    };

terms.relax.t_range    = {
        terms.lower_bound: "0",
        terms.upper_bound: "1e26"
    };


relax.p_threshold = 0.05;

relax.test_branches_name = "Test";
relax.reference_branches_name = "Reference";
relax.unclassified_branches_name = "Unclassified";

relax.display_orders = {terms.original_name: -1,
                        terms.json.nucleotide_gtr: 0,
                        relax.MG94_name: 1,
                        relax.general_descriptive_name: 4,
                        relax.alternative_name: 2,
                        relax.null_name: 3,
                        relax.partitioned_descriptive_name: 5
                       };


/*------------------------------------------------------------------------------*/

KeywordArgument ("multiple-files", "Use multiple files as input", "No");

relax.multiple_files = io.SelectAnOption ({
                                        {"No", "Load data and trees from a single file"}
                                        {"Yes", "Prompt for a file list and load files (and trees) from it"}
                                  }, "Use multiple files as input") == "Yes";


    
if (relax.multiple_files) {
    KeywordArgument ("filelist", "A line list of file paths for the alignments to include in this analysis");
    KeywordArgument ("code",      "Which genetic code should be used", "Universal");
    /**
        keyword, description (for inline documentation and help messages), default value
    */
} else {
    KeywordArgument ("code",      "Which genetic code should be used", "Universal");
    /**
        keyword, description (for inline documentation and help messages), default value
    */
    KeywordArgument ("alignment", "An in-frame codon alignment in one of the formats supported by HyPhy");
        /**
            keyword, description (for inline documentation and help messages), no default value,
            meaning that it will be required
        */
    
}

KeywordArgument ("tree",      "A phylogenetic tree (optionally annotated with {})", null, "Please select a tree file for the data:");
    /** the use of null as the default argument means that the default expectation is for the 
        argument to be missing, i.e. the tree is expected to be in the file
        the fourth, optional argument, can match this keyword with the dialog prompt / choice list title,
        meaning that it can only be consumed when this dialog prompt / choice list is invoked
        This allows handling some branching logic conditionals
    */


KeywordArgument ("mode",      "Run mode", "Classic mode", "Group test mode");
KeywordArgument ("test",  "Branches to use as the test set", null, "Choose the set of branches to use as the _test_ set");
KeywordArgument ("reference",  "Branches to use as the reference set", null, "Choose the set of branches to use as the _reference_ set");

io.DisplayAnalysisBanner ( relax.analysis_description );

selection.io.startTimer (relax.json [terms.json.timers], "Overall", 0);

namespace relax {
    LoadFunctionLibrary ("modules/shared-load-file.bf");
    load_file ({
                utility.getGlobalValue("terms.multiple_files") : multiple_files,
                utility.getGlobalValue("terms.prefix"): "relax", 
                utility.getGlobalValue("terms.settings") : {utility.getGlobalValue("terms.settings.branch_selector") : "relax.select_branches"}
               });
    LoadFunctionLibrary ("modules/grid_compute.ibf");
}

KeywordArgument ("output", "Write the resulting JSON to this file (default is to save to the same path as the alignment file + 'RELAX.json')", relax.codon_data_info [terms.json.json]);

relax.codon_data_info [terms.json.json] = io.PromptUserForFilePath ("Save the resulting JSON file to");

KeywordArgument ("grid-size", "The number of points in the initial distributional guess for likelihood fitting", 250);
relax.initial_grid.N = io.PromptUser ("The number of points in the initial distributional guess for likelihood fitting", 250, 1, 10000, TRUE);

KeywordArgument ("starting-points", "The number of initial random guesses to seed rate values optimization", 1);
relax.N.initial_guesses = io.PromptUser ("The number of initial random guesses to 'seed' rate values optimization", 1, 1, relax.initial_grid.N$10, TRUE);

io.ReportProgressMessageMD('RELAX',  'selector', 'Branch sets for RELAX analysis');

KeywordArgument ("multiple-hits",  "Include support for multiple nucleotide substitutions", "None");

relax.multi_hit = io.SelectAnOption ({
                                        {"Double", "Include branch-specific rates for double nucleotide substitutions"}
                                        {"Double+Triple", "Include branch-specific rates for double and triple nucleotide substitutions"}
                                        {"None", "[Default] Use standard models which permit only single nucleotide changes to occur instantly"}
                                  }, "Include support for multiple nucleotide substitutions");




if (relax.analysis_run_mode != relax.kGroupMode) {
	relax.branch_sets = {relax.test_branches_name : "test",
						 relax.reference_branches_name : "reference"}; 
} else {
	relax.branch_sets = {};
	utility.ForEachPair (relax.selected_branches[0], "_group_", "_value_","
		if (_value_ != relax.unclassified_branches_name) {
			if (utility.Has (relax.branch_sets, _value_, 'Number') == FALSE) {
				relax.branch_sets [_value_] = Abs (relax.branch_sets);
			}
		}
	");
}

relax.has_unclassified = FALSE;
relax.group_choices = {};

function relax.echo_group (group, description) {
	utility.ForEachPair (relax.selected_branches, "_partition_", "_selection_",
		"_selection_ = utility.Filter (_selection_, '_selector_', '_selector_ == group');
		 relax.group_choices[group] = 'Set ' + group + ' with ' + utility.Array1D (_selection_) + ' branches';
		 io.ReportProgressMessageMD('RELAX',  'selector', '\n* Selected ' + Abs(_selection_) + ' branches as the _`group`_ set: \\\`' + Join (', ',utility.Keys(_selection_)) + '\\\`')");

};



relax.branch_sets ["relax.echo_group"][""];

if (relax.analysis_run_mode == relax.kGroupMode) {
    KeywordArgument ("reference-group",  "Branches to use as the reference group", null);
	relax.reference_set_name = io.SelectAnOption (relax.group_choices, "Choose the set of branches to use as the _reference_ group");
	if (relax.branch_sets [relax.reference_set_name] > 0) {
		relax.k = utility.Keys (relax.branch_sets);
		for (relax.i = 0; relax.i < relax.numbers_of_tested_groups; relax.i += 1) {
			if (relax.k [relax.i] != relax.reference_set_name ) {
				relax.branch_sets [relax.k [relax.i]] = relax.branch_sets [relax.reference_set_name];
				break;
			}
		}
		relax.branch_sets [relax.reference_set_name] = 0;
	}	
}


utility.ForEachPair (relax.selected_branches, "_partition_", "_selection_",
    "_selection_ = utility.Filter (_selection_, '_value_', '_value_ == relax.unclassified_branches_name');
     relax.has_unclassified = Abs(_selection_) > 0;
     if (relax.has_unclassified) {
        io.ReportProgressMessageMD('RELAX',  'selector', '* ' + Abs(_selection_) + ' branches are in the unclassified (nuisance) set: \\\`' + Join (', ',utility.Keys(_selection_)) + '\\\`')
     }");

KeywordArgument ("rates", "The number omega rate classes to include in the model [2-10, default 3]", relax.rate_classes);
relax.rate_classes = io.PromptUser ("The number omega rate classes to include in the model", relax.rate_classes, 2, 10, TRUE);

KeywordArgument ("models", "Which version of the test to run (All or Minimal)", "All");
relax.model_set = io.SelectAnOption ({
                                        {"All", "[Default] Fit descriptive models AND run the relax test (4 models)"}
                                        {"Minimal", "Run only the RELAX test (2 models)"}
                                    }, "RELAX analysis type");
 
KeywordArgument ("srv", "Include synonymous rate variation", "No");          
                          
relax.do_srv = io.SelectAnOption (
                                {
                                    "Yes" : "Allow synonymous substitution rates to vary from site to site (but not from branch to branch)", 
                                    "Branch-site" : "Allow synonymous substitution rates to vary using general branch site models",
                                    "HMM" : "Allow synonymous substitution rates to vary from site to site (but not from branch to branch), and use a hidden markov model for spatial correlation on synonymous rates",
                                    "No"  : "Synonymous substitution rates are constant across sites. This is the 'classic' behavior, i.e. the original published test"
                                },
                                "Synonymous rate variation"
);

if (relax.do_srv == "Branch-site") {
    relax.do_bs_srv = TRUE;
    relax.do_srv = TRUE;
    selection.io.json_store_setting (relax.json, "SRV type", "Branch-site synonymous rate variation");
} else {
    if (relax.do_srv == "Yes") {
        relax.do_bs_srv = FALSE;
        relax.do_srv = TRUE;
        relax.do_srv_hmm  = FALSE; 
    } else {
        if (relax.do_srv == "HMM") {
         relax.do_srv      = TRUE;
         relax.do_bs_srv   = FALSE;
         relax.do_srv_hmm  = TRUE; 
        } else {
            relax.do_srv = FALSE;  
        } 
    }
}    

if (relax.do_srv) {
    KeywordArgument ("syn-rates", "The number alpha rate classes to include in the model [1-10, default 3]", relax.synonymous_rate_classes);
    relax.synonymous_rate_classes = io.PromptUser ("The number omega rate classes to include in the model", relax.synonymous_rate_classes, 1, 10, TRUE);
}  

selection.io.startTimer (relax.json [terms.json.timers], "Preliminary model fitting", 1);

namespace_tag = "relax";

namespace relax {
    doGTR ("relax");
}

estimators.fixSubsetOfEstimates(relax.gtr_results, relax.gtr_results[terms.global]);

namespace relax {
    scaler_prefix = "relax.scaler";
    doPartitionedMG ("relax", FALSE);
}


relax.run_full_mg94 = TRUE;
    
if (Type (relax.save_intermediate_fits) == "AssociativeList") {
    if (None != relax.save_intermediate_fits[^"terms.data.value"]) {
        if (utility.Has (relax.save_intermediate_fits[^"terms.data.value"], "Full-MG94", "AssociativeList")) {
            relax.final_partitioned_mg_results = (relax.save_intermediate_fits[^"terms.data.value"])["Full-MG94"];
            relax.run_full_mg94 = FALSE;
        }        
    }
}

io.ReportProgressMessageMD ("RELAX", "codon-refit", "Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model");

if (relax.run_full_mg94) {        
    relax.final_partitioned_mg_results = estimators.FitMGREV (relax.filter_names, relax.trees, relax.codon_data_info [terms.code], {
        terms.run_options.model_type: terms.local,
        terms.run_options.partitioned_omega: relax.selected_branches,
    }, relax.partitioned_mg_results);

    if (Type (relax.save_intermediate_fits) == "AssociativeList") {
        (relax.save_intermediate_fits[^"terms.data.value"])["Full-MG94"] = relax.final_partitioned_mg_results;        
        Export (lfe, ^relax.final_partitioned_mg_results[^"terms.likelihood_function"]);
        (relax.save_intermediate_fits[^"terms.data.value"])["Full-MG94-LF"] = lfe;
        io.SpoolJSON (relax.save_intermediate_fits[^"terms.data.value"],relax.save_intermediate_fits[^"terms.data.file"]);      
    }
}


io.ReportProgressMessageMD("RELAX", "codon-refit", "* " + selection.io.report_fit (relax.final_partitioned_mg_results, 0, relax.codon_data_info[terms.data.sample_size]));

io.ReportProgressMessageMD ("RELAX", "codon-refit", "* " + selection.io.report_fit_secondary_stats (relax.final_partitioned_mg_results));
        
selection.io.report_fit_secondary_stats (relax.final_partitioned_mg_results);

relax.global_dnds  = selection.io.extract_global_MLE_re (relax.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio);
relax.report_dnds = {};

utility.ForEach (relax.global_dnds, "_value_", '
    io.ReportProgressMessageMD ("RELAX", "codon-refit", "* " + _value_[terms.description] + " = " + Format (_value_[terms.fit.MLE],8,4));
    relax.report_dnds [(regexp.FindSubexpressions (_value_[terms.description], "^" + terms.parameters.omega_ratio + ".+\\*(.+)\\*$"))[1]] = {"0" : {terms.json.omega_ratio : _value_[terms.fit.MLE], terms.json.proportion : 1}};
');


//Store MG94 to JSON
selection.io.json_store_lf_withEFV (relax.json,
                            relax.MG94_name,
                            relax.final_partitioned_mg_results[terms.fit.log_likelihood],
                            relax.final_partitioned_mg_results[terms.parameters],
                            relax.sample_size,
                            utility.ArrayToDict (utility.Map (relax.global_dnds, "_value_", "{'key': _value_[terms.description], 'value' : Eval({{_value_ [terms.fit.MLE],1}})}")),
                            (relax.final_partitioned_mg_results[terms.efv_estimate])["VALUEINDEXORDER"][0],
                            relax.display_orders[relax.MG94_name]);
                            
utility.ForEachPair (relax.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(relax.json,relax.MG94_name, terms.branch_length, relax.display_orders[relax.MG94_name],
                                             _key_,
                                             selection.io.extract_branch_info((relax.final_partitioned_mg_results[terms.branch_length])[_key_], "selection.io.branch.length"));');



if (relax.multi_hit == "None") {
    relax.model_generator = "models.codon.BS_REL.ModelDescription";

    if (relax.do_srv) {
        if (relax.do_bs_srv) {
            relax.model_generator = "models.codon.BS_REL_SRV.ModelDescription";
            relax.rate_class_arguments = {{relax.synonymous_rate_classes__,relax.rate_classes__}};
        } else {
    
            lfunction relax.model.with.GDD (type, code, rates) {        
                def = models.codon.BS_REL.ModelDescription (type, code, rates);
                options = {utility.getGlobalValue("terms.rate_variation.bins") : utility.getGlobalValue("relax.synonymous_rate_classes"),
                           utility.getGlobalValue("terms._namespace") : "relax._shared_srv"};

                if (utility.getGlobalValue("relax.do_srv_hmm")) {
                        options [utility.getGlobalValue ("terms.rate_variation.HMM")] = "equal";
                }
                def [utility.getGlobalValue("terms.model.rate_variation")] = rate_variation.types.GDD.factory (options);
                return def;
            }
        
            relax.model_generator = "relax.model.with.GDD";
            relax.rate_class_arguments = relax.rate_classes;
        }
    } else {
        relax.rate_class_arguments = relax.rate_classes;
    }
} else {

    lfunction relax.model.BS_REL_MH (type, code, rates) {        
        def = models.codon.BS_REL.ModelDescription (type, code, rates);
        def [utility.getGlobalValue("terms.model.defineQ")] = "models.codon.BS_REL._DefineQ.MH";
        if (^'relax.multi_hit' == "Double") {
            def [utility.getGlobalValue("terms.model.rate_generator")] = "function rate_generator (fromChar, toChar, namespace, model_type, model) {
               return models.codon.MG_REV_MH._GenerateRate_generic (fromChar, toChar, namespace, model_type, model[utility.getGlobalValue('terms.translation_table')],
                'alpha', utility.getGlobalValue('terms.parameters.synonymous_rate'),
                'beta_' + component, terms.AddCategory (utility.getGlobalValue('terms.parameters.nonsynonymous_rate'), component),
                'omega' + component, terms.AddCategory (utility.getGlobalValue('terms.parameters.omega_ratio'), component),
                'delta', utility.getGlobalValue('terms.parameters.multiple_hit_rate')
               );}";
            }  else {
                def [utility.getGlobalValue("terms.model.rate_generator")] = "function rate_generator (fromChar, toChar, namespace, model_type, model) {
                       return models.codon.MG_REV_TRIP._GenerateRate_generic (fromChar, toChar, namespace, model_type, model[utility.getGlobalValue('terms.translation_table')],
                        'alpha', utility.getGlobalValue('terms.parameters.synonymous_rate'),
                        'beta_' + component, terms.AddCategory (utility.getGlobalValue('terms.parameters.nonsynonymous_rate'), component),
                        'omega' + component, terms.AddCategory (utility.getGlobalValue('terms.parameters.omega_ratio'), component),
                        'delta', utility.getGlobalValue('terms.parameters.multiple_hit_rate'),
                        'psi', utility.getGlobalValue('terms.parameters.triple_hit_rate'),
                        'psi', utility.getGlobalValue('terms.parameters.triple_hit_rate')
                       );}"                
        }       
        return def;
    }
    
    relax.model_generator = "relax.model.BS_REL_MH";
 
    if (relax.do_srv) {
        if (relax.do_bs_srv) {
            relax.model_generator = "models.codon.BS_REL_SRV.ModelDescription";
            relax.rate_class_arguments = {{relax.synonymous_rate_classes__,relax.rate_classes__}};
        } else {
    
            lfunction relax.model.MH.with.GDD (type, code, rates) {        
                def = relax.model.BS_REL_MH (type, code, rates);
                options = {utility.getGlobalValue("terms.rate_variation.bins") : utility.getGlobalValue("relax.synonymous_rate_classes"),
                           utility.getGlobalValue("terms._namespace") : "relax._shared_srv"};

                if (utility.getGlobalValue("relax.do_srv_hmm")) {
                        options [utility.getGlobalValue ("terms.rate_variation.HMM")] = "equal";
                }
                def [utility.getGlobalValue("terms.model.rate_variation")] = rate_variation.types.GDD.factory (options);
                return def;
            }
        
            relax.model_generator = "relax.model.MH.with.GDD";
            relax.rate_class_arguments = relax.rate_classes;
        }
    } else {
        relax.rate_class_arguments = relax.rate_classes;
    }
}

selection.io.stopTimer (relax.json [terms.json.timers], "Preliminary model fitting");
parameters.DeclareGlobalWithRanges (relax.relaxation_parameter, 1, 0, 50);

PARAMETER_GROUPING = {};
PARAMETER_GROUPING + relax.distribution["rates"];
PARAMETER_GROUPING + relax.distribution["weights"];

if (relax.model_set == "All") { // run all the models

    relax.ge_guess = None;

    while (1) {

        relax.ge.bsrel_model =  model.generic.DefineMixtureModel("relax.BS_REL.ModelDescription",
                "relax.ge", {
                    "0": parameters.Quote(terms.local),
                    "1": relax.codon_data_info[terms.code],
                    "2": parameters.Quote (relax.rate_classes) // the number of rate classes
                },
                relax.filter_names,
                None);
                

        relax.distribution          = models.codon.BS_REL.ExtractMixtureDistribution(relax.ge.bsrel_model);
        relax.weight_multipliers    = parameters.helper.stick_breaking (utility.SwapKeysAndValues(utility.MatrixToDict(relax.distribution["weights"])),None);
        relax.constrain_parameters   = parameters.ConstrainMeanOfSet(relax.distribution["rates"],relax.weight_multipliers,1,"relax");
 
        
        relax.i = 0;
        for (key, value; in; relax.constrain_parameters[terms.global]){
            model.generic.AddGlobal (relax.ge.bsrel_model, value, key);
            relax.i += 1;
            if (relax.i < relax.rate_classes) {
                parameters.SetRange (value, terms.range_almost_01);
            } else {
                parameters.SetRange (value, terms.range_gte1);
            }
            
        }
        
        
        relax.distribution["rates"] = Transpose (utility.Values (relax.constrain_parameters[terms.global]));
        
        for (relax.i = 1; relax.i < relax.rate_classes; relax.i += 1) {
            parameters.SetRange (model.generic.GetGlobalParameter (relax.ge.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,relax.i)), terms.range_almost_01);
        }
        parameters.SetRange (model.generic.GetGlobalParameter (relax.ge.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,relax.rate_classes)), terms.range_gte1);
        
        // constrain the mean of this distribution to 1
        

        relax.model_object_map = { "relax.ge" :       relax.ge.bsrel_model };

        io.ReportProgressMessageMD ("RELAX", "gd", "Fitting the general descriptive (separate k per branch) model");
        selection.io.startTimer (relax.json [terms.json.timers], "General descriptive model fitting", 2);

        //utility.SetEnvVariable ("VERBOSITY_LEVEL",10);
        
        relax.ge_model_map = {};
        for (relax.index, relax.junk ; in; relax.filter_names) {
             relax.ge_model_map [relax.index] = {"DEFAULT" : "relax.ge"};
        }



        if (Type (relax.ge_guess) != "Matrix") {
        
        
            // first time in 
            relax.initial.test_mean    =
            
            math.Mean ( 
                utility.Map (selection.io.extract_global_MLE_re (relax.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio + ".+"), "_v_", "_v_[terms.fit.MLE]"));
                
                
            relax.init_grid_setup       (relax.distribution);
            relax.initial_grid         = estimators.LHC (relax.initial_ranges,relax.initial_grid.N);
            relax.initial_grid = utility.Map (relax.initial_grid, "_v_", 
                'relax._renormalize (_v_, "relax.distribution", relax.initial.test_mean)'
            );
            relax.nm.precision = -0.00025*relax.final_partitioned_mg_results[terms.fit.log_likelihood];
            
            
            relax.nm.bl_scalers = {};
            
            
            for (relax.index, relax.junk ; in; relax.filter_names) {
                relax.scaler.id = "relax.bl.scaler." + relax.index;
                relax.nm.bl_scalers[relax.index] = relax.scaler.id ;
                parameters.DeclareGlobalWithRanges (relax.scaler.id, 1, 0, 1000);
                for (i, v; in; relax.initial_grid) {
                    v[relax.scaler.id] = {terms.id : relax.scaler.id, terms.fit.MLE : Random (2,4)};
                }
            }            
            
                           
            relax.grid_search.results =  estimators.FitLF (relax.filter_names, relax.trees, relax.ge_model_map,
                relax.final_partitioned_mg_results,
                relax.model_object_map, 
                {
                    "retain-lf-object": TRUE,
                    terms.run_options.apply_user_constraints: "relax.init.k",
                    terms.run_options.proportional_branch_length_scaler : 
                                                            relax.nm.bl_scalers,
                    
                    terms.run_options.optimization_settings : 
                        {
                            "OPTIMIZATION_METHOD" : "nedler-mead",
                            "MAXIMUM_OPTIMIZATION_ITERATIONS" : 500,
                            "OPTIMIZATION_PRECISION" : relax.nm.precision
                        } ,
             
                    terms.search_grid : relax.initial_grid,
                    terms.search_restarts : relax.N.initial_guesses,
                    terms.run_options.optimization_log : relax.optimization_log_file (".GE-1-log.json")
                }
            );
            
 
            relax.general_descriptive.fit =  estimators.FitLF (relax.filter_names,
                                        relax.trees,
                                        relax.ge_model_map,
                                        relax.grid_search.results,
                                        relax.model_object_map,
                                        {
                                            terms.run_options.apply_user_constraints: "relax.init.k",
                                            terms.run_options.retain_lf_object : TRUE,
                                            terms.run_options.optimization_log :  relax.optimization_log_file(".GE-2-log.json")

                                        });
                                        
           
      } else {
            parameters.SetStickBreakingDistribution (relax.distribution, relax.ge_guess);
            relax.general_descriptive.fit =  estimators.FitLF (relax.filter_names,
                                            relax.trees,
                                            relax.ge_model_map,
                                            relax.final_partitioned_mg_results,
                                            relax.model_object_map,
                                            {
                                                terms.run_options.apply_user_constraints: "relax.init.k",
                                                terms.run_options.retain_lf_object : TRUE,
                                                terms.run_options.run_options.optimization_log :  relax.optimization_log_file(".GE-log.json")

                                            });
       }
       

        estimators.TraverseLocalParameters (relax.general_descriptive.fit [terms.likelihood_function], relax.model_object_map, "relax.set.k2");
        relax.general_descriptive.fit = estimators.FitExistingLF (relax.general_descriptive.fit [terms.likelihood_function], relax.model_object_map);

        
        selection.io.stopTimer (relax.json [terms.json.timers], "General descriptive model fitting");

        io.ReportProgressMessageMD("RELAX", "ge", "* " + selection.io.report_fit (relax.general_descriptive.fit, 9, relax.codon_data_info[terms.data.sample_size]));
        io.ReportProgressMessageMD("RELAX", "ge", "* The following baseline rate distribution for branch-site combinations was inferred");
        relax.inferred_ge_distribution = parameters.GetStickBreakingDistribution (models.codon.BS_REL.ExtractMixtureDistributionFromFit (relax.ge.bsrel_model, relax.general_descriptive.fit)) % 0;

        selection.io.report_dnds (relax.inferred_ge_distribution);

        if (relax.rate_classes > 2) {
            if (relax.inferred_ge_distribution[0][1] < 1e-5 || relax.inferred_ge_distribution[1][1] < 1e-5) {
                io.ReportProgressMessageMD("RELAX", "ge", "\n ### Because some of the rate classes were collapsed to 0, the model is likely overparameterized. RELAX will reduce the number of site rate classes by one and repeat the fit now.\n----\n");
                relax.rate_classes = relax.rate_classes - 1;
                relax.ge_guess = {relax.rate_classes, 2};
                relax.shift    = 0;
                //console.log (relax.inferred_ge_distribution);
                for (relax.i = 0; relax.i < relax.rate_classes; relax.i += 1) {
                    if (relax.inferred_ge_distribution[relax.i][1] < 1e-5 && relax.shift == 0) {
                        relax.shift += 1;
                        continue;
                    }
                    relax.ge_guess[relax.i][0] = relax.inferred_ge_distribution[relax.i + relax.shift][0];
                    relax.ge_guess[relax.i][1] = relax.inferred_ge_distribution[relax.i + relax.shift][1];
                }
                continue;
            }
        }


        relax.distribution_for_json = {'Shared' : utility.Map (utility.Range (relax.rate_classes, 0, 1),
                                                             "_index_",
                                                             "{terms.json.omega_ratio : relax.inferred_ge_distribution [_index_][0],
                                                               terms.json.proportion  : relax.inferred_ge_distribution [_index_][1]}")
                                       };
                                       
                                       
        relax.report_multi_hit (relax.general_descriptive.fit, relax.distribution_for_json, "RELAX", "gd-mh");
        selection.io.json_store_lf (relax.json,
                                    relax.general_descriptive_name,
                                    relax.general_descriptive.fit[terms.fit.log_likelihood],
                                    relax.general_descriptive.fit[terms.parameters] + 9 , // +9 comes from CF3x4
                                    relax.codon_data_info[terms.data.sample_size],
                                    relax.distribution_for_json,
                                    relax.display_orders[relax.general_descriptive_name]
                                );

        selection.io.json_store_branch_attribute(relax.json, relax.general_descriptive_name, terms.branch_length, relax.display_orders[relax.general_descriptive_name],
                                                     0,
                                                     selection.io.extract_branch_info((relax.general_descriptive.fit[terms.branch_length])[0], "selection.io.branch.length"));

        relax.k_estimates = selection.io.extract_branch_info((relax.general_descriptive.fit[terms.branch_length])[0], "relax.extract.k");

        relax.k_stats = math.GatherDescriptiveStats (utility.Map (utility.UniqueValues (relax.k_estimates), "_value_", "0+_value_"));

        io.ReportProgressMessageMD("RELAX", "ge", "* Branch-level `terms.relax.k` distribution has mean " + Format (relax.k_stats[terms.math.mean], 5,2) + ", median " +
                                                     Format (relax.k_stats[terms.math.median], 5,2) + ", and 95% of the weight in " + Format (relax.k_stats[terms.math._2.5], 5,2) + " - " + Format (relax.k_stats[terms.math._97.5], 5,2));


        selection.io.json_store_branch_attribute(relax.json, "k (general descriptive)", terms.json.branch_label, relax.display_orders[relax.general_descriptive_name],
                                                     0,
                                                     relax.k_estimates);


        for (relax.i = 1; relax.i <= relax.rate_classes; relax.i += 1) {
            //console.log (model.generic.GetGlobalParameter (relax.ge.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,relax.i)));
            parameters.RemoveConstraint (model.generic.GetGlobalParameter (relax.ge.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,relax.i)));
        }

        break;
    }

} else {
    // MINIMAL models branch
    relax.general_descriptive.fit = relax.final_partitioned_mg_results;
}

/* now fit the two main models for RELAX */

if (relax.analysis_run_mode != relax.kGroupMode) {
	io.ReportProgressMessageMD ("RELAX", "alt", "Fitting the alternative model to test K != 1");
} else {
	io.ReportProgressMessageMD ("RELAX", "alt", "Fitting the alternative model with individual K parameters for " + relax.numbers_of_tested_groups + " branch groups");
}


selection.io.startTimer (relax.json [terms.json.timers], "RELAX alternative model fitting", 3);

relax.model_object_map = {};
relax.model_to_relax_parameter = {};

if (relax.analysis_run_mode != relax.kGroupMode) {
	relax.model_object_map  		["relax.reference"] = None;
	relax.reference_model_namespace 		= "relax.reference";
	relax.model_object_map  		["relax.test"] = None;
	relax.model_to_relax_parameter  ["relax.test"] = relax.relaxation_parameter;
	relax.model_namespaces = {{"relax.reference","relax.test"}};
	relax.model_to_group_name = {"relax.reference" : relax.reference_branches_name, "relax.test" : relax.test_branches_name};
	relax.model_for_srv      = "relax.test";

} else {
	relax.model_to_group_name = {};
	relax.model_namespaces = {};
	utility.ForEachPair  (relax.branch_sets, "_key_", "_value_", "
		relax.k = 'relax.model_' + _value_ ;
		relax.model_object_map[relax.k] = None;
		relax.model_to_group_name [relax.k] = _key_;
		relax.model_to_relax_parameter  [relax.k] = relax.relaxation_parameter + '_' + _value_;
		relax.model_namespaces[_value_] = relax.k;
	");
	relax.reference_model_namespace = relax.model_namespaces[0];
	relax.model_for_srv      = relax.model_namespaces[0];
}

 

//console.log (relax.model_namespaces);
//explicit loop to avoid re-entrance errors 

relax.relax_parameter_terms = {};
relax.bound_weights         = {};

for (relax.k = 0; relax.k < relax.numbers_of_tested_groups; relax.k += 1) {
	relax.model_nmsp = relax.model_namespaces[relax.k ];
	relax.model_object_map[relax.model_nmsp] = model.generic.DefineMixtureModel(relax.model_generator,
        relax.model_nmsp, {
            "0": parameters.Quote(terms.global),
            "1": relax.codon_data_info[terms.code],
            "2": parameters.Quote (relax.rate_class_arguments) // the number of rate classes
        },
        relax.filter_names,
        None);
        
	for (relax.i = 1; relax.i < relax.rate_classes; relax.i += 1) {
		parameters.SetRange (model.generic.GetGlobalParameter (relax.model_object_map[relax.model_nmsp] , terms.AddCategory (terms.parameters.omega_ratio,relax.i)), terms.range01);
	}
	parameters.SetRange (model.generic.GetGlobalParameter (relax.model_object_map[relax.model_nmsp] , terms.AddCategory (terms.parameters.omega_ratio,relax.i)), terms.range_gte1);

    if (relax.k > 0) { // not the reference group
    	 parameters.DeclareGlobalWithRanges (relax.model_to_relax_parameter  [relax.model_nmsp], 1, 0, 50);
    	 if (relax.k > 1) {
    	 	relax.relax_parameter_terms [relax.k] = terms.relax.k + " for " + relax.model_to_group_name[relax.model_nmsp];
    	 	//console.log (relax.relax_parameter_terms [relax.k]);
    	 	model.generic.AddGlobal (relax.model_object_map[relax.model_nmsp], relax.model_to_relax_parameter  [relax.model_nmsp], relax.relax_parameter_terms [relax.k]);
    	 } else {
    	 	model.generic.AddGlobal (relax.model_object_map[relax.model_nmsp], relax.model_to_relax_parameter  [relax.model_nmsp], terms.relax.k);
    	 }
    	 relax.bound_weights * models.BindGlobalParameters ({"0" : relax.model_object_map[relax.reference_model_namespace], "1" : relax.model_object_map[relax.model_nmsp]}, terms.mixture.mixture_aux_weight + " for category");
		 models.BindGlobalParameters ({"0" : relax.model_object_map[relax.reference_model_namespace], "1" : relax.model_object_map[relax.model_nmsp]}, terms.nucleotideRate("[ACGT]","[ACGT]"));
		 for (relax.i = 1; relax.i <= relax.rate_classes; relax.i += 1) {
			parameters.SetConstraint (model.generic.GetGlobalParameter (relax.model_object_map[relax.model_nmsp] , terms.AddCategory (terms.parameters.omega_ratio,relax.i)),
									  model.generic.GetGlobalParameter (relax.model_object_map[relax.reference_model_namespace] , terms.AddCategory (terms.parameters.omega_ratio,relax.i)) + "^" + relax.model_to_relax_parameter  [relax.model_nmsp],
									  terms.global);
		}
		if (relax.do_bs_srv) {
		    for (description,var_id; in;  ((relax.model_object_map[relax.reference_model_namespace])[terms.parameters])[terms.global]) {
		        
		        if (None != regexp.Find (description, terms.parameters.synonymous_rate + " for category")) {
		            var_id_2 = utility.GetByKey (((relax.model_object_map[relax.model_nmsp])[terms.parameters])[terms.global], description, "String");
		            if (None != var_id_2) {
		                parameters.SetConstraint (var_id, var_id_2, terms.global);                
		             }
		        }

                if (None != regexp.Find (description, terms.mixture.mixture_aux_weight + " for category SRV ")) {
		            var_id_2 = utility.GetByKey (((relax.model_object_map[relax.model_nmsp])[terms.parameters])[terms.global], description, "String");
		            if (None != var_id_2) {
		                if (parameters.IsIndependent (var_id_2)) {
		                    parameters.SetConstraint (var_id, var_id_2, terms.global);                
		                }
		             }
		        }		        
 		    }
		    
		    
		}
   }
}

relax.model_map = {};
for (relax.index, relax.junk ; in; relax.filter_names) {
    relax.model_map [relax.index]= utility.Map (relax.model_to_group_name, "_groupid_", 'utility.Filter (relax.selected_branches[relax.index], "_branchgroup_", "_branchgroup_ == \'" + _groupid_ + "\'")');
}


// constrain the proportions to be the same

relax.do_lhc = FALSE;

if (relax.do_srv)  {

    if (relax.do_bs_srv) {
        relax.srv_rate_regex   = "Mean scaler variable for";
        relax.srv_weight_regex = terms.AddCategory (utility.getGlobalValue("terms.mixture.mixture_aux_weight"), "SRV [0-9]+");
        
        relax.srv_rate_reporting = regexp.PartitionByRegularExpressions(utility.Keys (((relax.model_object_map[relax.model_for_srv])[terms.parameters])[terms.global]), 
            {"0" : "^" + utility.getGlobalValue('terms.parameters.synonymous_rate'), 
             "1" :  relax.srv_weight_regex});


        relax.srv_rate_reporting = {
          'rates' : utility.UniqueValues (utility.Map ( relax.srv_rate_reporting ["^" + utility.getGlobalValue('terms.parameters.synonymous_rate') ]  , "_value_", '(((relax.model_object_map[relax.model_for_srv])[terms.parameters])[terms.global])[_value_]')),
           'weights' : utility.UniqueValues (utility.Map (relax.srv_rate_reporting [relax.srv_weight_regex ]  , "_value_", '(((relax.model_object_map[relax.model_for_srv])[terms.parameters])[terms.global])[_value_]'))
        };
        
        
        
    } else {
        relax.srv_rate_regex  = "GDD rate category [0-9]+";
        relax.srv_weight_regex = "Mixture auxiliary weight for GDD category [0-9]+";
    }
    relax.srv_distribution = regexp.PartitionByRegularExpressions(utility.Keys (((relax.model_object_map[relax.model_for_srv])[terms.parameters])[terms.global]), {"0" : relax.srv_rate_regex, "1" : relax.srv_weight_regex});

    
    relax.srv_distribution = {
        'rates' : utility.UniqueValues (utility.Map (relax.srv_distribution [relax.srv_rate_regex ]  , "_value_", '(((relax.model_object_map[relax.model_for_srv])[terms.parameters])[terms.global])[_value_]')),
        'weights' : utility.UniqueValues (utility.Map (relax.srv_distribution [relax.srv_weight_regex ]  , "_value_", '(((relax.model_object_map[relax.model_for_srv])[terms.parameters])[terms.global])[_value_]'))
    };
    
 
    PARAMETER_GROUPING + relax.srv_distribution["rates"];
    PARAMETER_GROUPING + relax.srv_distribution["weights"];

    relax.init_grid_setup (relax.srv_distribution);
    
}

if (relax.model_set != "All") {
    relax.do_lhc = TRUE;
    relax.distribution = models.codon.BS_REL.ExtractMixtureDistribution(relax.model_object_map[relax.reference_model_namespace]);
    relax.init_grid_setup        (relax.distribution);
    relax.initial_ranges[relax.relaxation_parameter] = {terms.lower_bound : 0.10,
                    terms.upper_bound : 2.0};
     
                    
}

if (relax.has_unclassified) {

    relax.unclassified.bsrel_model =  model.generic.DefineMixtureModel(relax.model_generator,
        "relax.unclassified", {
            "0": parameters.Quote(terms.global),
            "1": relax.codon_data_info[terms.code],
            "2": parameters.Quote (relax.rate_class_arguments) // the number of rate classes
        },
        relax.filter_names,
        None);

    for (relax.i = 1; relax.i < relax.rate_classes-1; relax.i += 1) {
        parameters.SetRange (model.generic.GetGlobalParameter (relax.unclassified.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,relax.i)), terms.range01);
    }

    if (relax.do_lhc) {
        relax.distribution_uc = models.codon.BS_REL.ExtractMixtureDistribution(relax.unclassified.bsrel_model);
        relax.init_grid_setup  (relax.distribution_uc);
    }
    

    parameters.SetRange (model.generic.GetGlobalParameter (relax.unclassified.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,relax.rate_classes)), terms.range_gte1);

    relax.model_object_map ["relax.unclassified"] = relax.unclassified.bsrel_model;
    
    for (relax.index, relax.junk ; in; relax.filter_names) {
        (relax.model_map[relax.index]) ["relax.unclassified"] = utility.Filter (relax.selected_branches[relax.index], '_value_', '_value_ == relax.unclassified_branches_name');
     }
        
    models.BindGlobalParameters ({"0" : relax.model_object_map[relax.reference_model_namespace], "1" : relax.unclassified.bsrel_model}, terms.nucleotideRate("[ACGT]","[ACGT]"));
}

if (relax.do_lhc) {
    relax.initial_grid         = estimators.LHC (relax.initial_ranges,relax.initial_grid.N);
    relax.initial.test_mean    = ((selection.io.extract_global_MLE_re (relax.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio + ".+`relax.reference_branches_name`.+"))["0"])[terms.fit.MLE];

 
    relax.initial_grid = utility.Map (relax.initial_grid, "_v_", 
        'relax._renormalize_with_weights (_v_, "relax.distribution", relax.initial.test_mean)'
    );
    
    //console.log (relax.initial_grid[0]);
    if (relax.has_unclassified) {
        relax.initial.unc_mean    = ((selection.io.extract_global_MLE_re (relax.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio + ".+`relax.unclassified_branches_name`.+"))["0"])[terms.fit.MLE];
 
        relax.initial_grid = utility.Map (relax.initial_grid, "_v_", 
            'relax._renormalize_with_weights (_v_, "relax.distribution_uc", relax.initial.unc_mean)'
        );
    }
    //console.log (relax.initial_grid[0]);

}


//------------------------------------



function relax.report_multi_class_rates (model_fit, distributions) {
    io.ReportProgressMessageMD("RELAX", "alt", "* The following rate distribution was inferred for **" + relax.model_to_group_name[relax.reference_model_namespace] + "** (reference) branches");
    relax.inferred_distribution_ref = parameters.GetStickBreakingDistribution (models.codon.BS_REL.ExtractMixtureDistribution (relax.model_object_map[relax.reference_model_namespace])) % 0;
    selection.io.report_dnds (relax.inferred_distribution_ref);

    relax.distribution_for_json = {};
    relax.distribution_for_json [relax.model_to_group_name[relax.reference_model_namespace]] = utility.Map (utility.Range (relax.rate_classes, 0, 1),
                                                     "_index_",
                                                     "{terms.json.omega_ratio : relax.inferred_distribution_ref [_index_][0],
                                                       terms.json.proportion  : relax.inferred_distribution_ref [_index_][1]}");

    if (None != model_fit || distributions) {

        if (None != model_fit) {
            relax.fitted.K = {};
            relax.fitted.K [relax.model_to_group_name[relax.reference_model_namespace]] = 1;
        }


        for (relax.k = 1; relax.k < relax.numbers_of_tested_groups; relax.k += 1) {
            relax.model_nmsp = relax.model_namespaces[relax.k ];
            relax.branch_set = relax.model_to_group_name[relax.model_nmsp];
            if (None != model_fit) {
                if (relax.k > 1) {
                    relax.fitted.K_group = estimators.GetGlobalMLE (model_fit,relax.relax_parameter_terms[relax.k]);
                } else {
                    relax.fitted.K_group = estimators.GetGlobalMLE (model_fit,terms.relax.k);
                }
                relax.fitted.K [relax.model_to_group_name[relax.model_nmsp]] = relax.fitted.K_group ;
        
                io.ReportProgressMessageMD("RELAX", "alt", "* Relaxation/intensification parameter (K) for branch set `relax.branch_set` = " + Format(relax.fitted.K_group ,8,2));
            }
            io.ReportProgressMessageMD("RELAX", "alt", "* The following rate distribution was inferred for **`relax.branch_set`** branches");
            relax.inferred_distribution = parameters.GetStickBreakingDistribution (models.codon.BS_REL.ExtractMixtureDistribution (relax.model_object_map[relax.model_nmsp])) % 0;
            selection.io.report_dnds (relax.inferred_distribution);
            relax.distribution_for_json [relax.model_to_group_name[relax.model_nmsp]] = utility.Map (utility.Range (relax.rate_classes, 0, 1),
                                                         "_index_",
                                                         "{terms.json.omega_ratio : relax.inferred_distribution [_index_][0],
                                                           terms.json.proportion  : relax.inferred_distribution [_index_][1]}");
        }
        
        
        
    }
}

//------------------------------------

function relax._report_srv (relax_model_fit, is_null) {

    if (relax.do_srv_hmm) {
        relax.hmm_lambda = selection.io.extract_global_MLE (relax_model_fit, terms.rate_variation.hmm_lambda);
        if (is_null) {
            io.ReportProgressMessageMD("relax", "main", "* HMM switching rate = " +  Format (relax.hmm_lambda, 8, 3));
            relax.distribution_for_json [terms.rate_variation.hmm_lambda] = relax.hmm_lambda;    
        } else {
            relax.hmm_lambda.CI = parameters.GetProfileCI(((relax_model_fit[terms.global])[terms.rate_variation.hmm_lambda])[terms.id],
                                        relax_model_fit[terms.likelihood_function], 0.95);
 
            io.ReportProgressMessageMD ("relax", "main", "* HMM switching rate = " + Format (relax.hmm_lambda,8,4) + 
                            " (95% profile CI " + Format ((relax.hmm_lambda.CI )[terms.lower_bound],8,4) + "-" + Format ((relax.hmm_lambda.CI )[terms.upper_bound],8,4) + ")");
                    
            relax.distribution_for_json [terms.rate_variation.hmm_lambda] = relax.hmm_lambda.CI;    
                        
        }
        
    }

    if (relax.do_srv) {
        if (relax.do_bs_srv) {
            relax.srv_info = parameters.GetStickBreakingDistribution ( relax.srv_rate_reporting) % 0;
        } else {
            relax.srv_info = Transpose((rate_variation.extract_category_information((relax.model_object_map[relax.model_for_srv])))["VALUEINDEXORDER"][0])%0;
        }
        io.ReportProgressMessageMD("RELAX", "alt", "* The following rate distribution for site-to-site **synonymous** rate variation was inferred");
        selection.io.report_distribution (relax.srv_info);

        relax.distribution_for_json [relax.SRV] = (utility.Map (utility.Range (relax.synonymous_rate_classes, 0, 1),
                                                                 "_index_",
                                                                 "{terms.json.rate :relax.srv_info [_index_][0],
                                                                   terms.json.proportion : relax.srv_info [_index_][1]}"));
                                   
        if (is_null == FALSE) {                        
            ConstructCategoryMatrix (relax.cmx, ^(relax_model_fit[terms.likelihood_function]));
            ConstructCategoryMatrix (relax.cmx_weights, ^(relax_model_fit[terms.likelihood_function]), WEIGHTS);
            relax.cmx_weighted         = (relax.cmx_weights[-1][0]) $ relax.cmx; // taking the 1st column fixes a bug with multiple partitions 
            relax.column_weights       = {1, Rows (relax.cmx_weights)}["1"] * relax.cmx_weighted;
            relax.column_weights       = relax.column_weights["1/_MATRIX_ELEMENT_VALUE_"];
            (relax.json [relax.json.srv_posteriors]) =  relax.cmx_weighted $ relax.column_weights;
    	}

        if (relax.do_srv_hmm && is_null == FALSE ) {
            ConstructCategoryMatrix (relax.cmx_viterbi, ^(relax_model_fit[terms.likelihood_function]), SHORT);
            (relax.json [relax.json.srv_viterbi]) = relax.cmx_viterbi;
            io.ReportProgressMessageMD("relax", "main", "* The following switch points for synonymous rates were inferred");
            selection.io.report_viterbi_path (relax.cmx_viterbi);
    
        }
    }
    
    if (relax.do_srv_hmm) {
        relax.distribution_for_json [terms.rate_variation.hmm_lambda] = relax.hmm_lambda;
    }
}

//------------------------------------

function relax.FitMainTestPair (prompt) {
     //_varname_ = model.generic.GetGlobalParameter (relax.model_object_map[relax.model_namespaces[1]] , terms.AddCategory (terms.parameters.omega_ratio,1));

     
     if (relax.do_lhc) {

        relax.nm.precision = -0.00025*relax.final_partitioned_mg_results[terms.fit.log_likelihood];
        //parameters.DeclareGlobalWithRanges ("relax.bl.scaler", 1, 0, 1000);
        
        relax.main.bl_scalers = {};
            
        //utility.SetEnvVariable ("VERBOSITY_LEVEL",10);
            
        for (relax.index, relax.junk ; in; relax.filter_names) {
            relax.scaler.id = "relax.bl.scaler." + relax.index;
            relax.main.bl_scalers[relax.index] = relax.scaler.id ;
            parameters.DeclareGlobalWithRanges (relax.scaler.id, 1, 0, 1000);
        }       
        

        
        relax.general_descriptive.fit =  estimators.FitLF (relax.filter_names, relax.trees, relax.model_map,
                                    relax.general_descriptive.fit,
                                    relax.model_object_map, 
                                    {
                                        "retain-lf-object": TRUE,
                                        terms.run_options.proportional_branch_length_scaler : 
                                                                                relax.main.bl_scalers,
                                        
                                        terms.run_options.optimization_settings : 
                                            {
                                                "OPTIMIZATION_METHOD" : "nedler-mead",
                                                "MAXIMUM_OPTIMIZATION_ITERATIONS" : 500,
                                                "OPTIMIZATION_PRECISION" : relax.nm.precision
                                            } ,
                                 
                                        terms.search_grid : relax.initial_grid,
                                        terms.search_restarts : relax.N.initial_guesses
                        
                                    }
        );
        
    }
    
    relax.general_descriptive.fit = Eval (relax.general_descriptive.fit);
    
	relax.alternative_model.fit =  estimators.FitLF (relax.filter_names, 
	                                                 relax.trees, 
	                                                 relax.model_map, 
	                                                 relax.general_descriptive.fit, 
	                                                 relax.model_object_map, 
	                                                 {
	                                                    terms.run_options.retain_lf_object: TRUE, 
	                                                    terms.run_options.optimization_log :  relax.optimization_log_file ( "MainALT-log.json")});
	                                                    
	io.ReportProgressMessageMD("RELAX", "alt", "* " + selection.io.report_fit (relax.alternative_model.fit, 9, relax.codon_data_info[terms.data.sample_size]));

    if (prompt) {
        KeywordArgument ("save-fit", "Save RELAX alternative model fit to this file (default is not to save)", "/dev/null");
        relax.save_fit_path = io.PromptUserForFilePath ("Save RELAX model fit to this file ['/dev/null' to skip]");
    } 
    
    
    relax.stashLF = estimators.TakeLFStateSnapshot (relax.alternative_model.fit[terms.likelihood_function]);

    io.SpoolLFToPath(relax.alternative_model.fit[terms.likelihood_function], relax.save_fit_path);

	if (relax.numbers_of_tested_groups == 2 && relax.analysis_run_mode != relax.kGroupMode) {

		relax.fitted.K = estimators.GetGlobalMLE (relax.alternative_model.fit,terms.relax.k);
		io.ReportProgressMessageMD("RELAX", "alt", "* Relaxation/intensification parameter (K) = " + Format(relax.fitted.K,8,2));
		io.ReportProgressMessageMD("RELAX", "alt", "* The following rate distribution was inferred for **test** branches");
		relax.inferred_distribution = parameters.GetStickBreakingDistribution (models.codon.BS_REL.ExtractMixtureDistribution (relax.model_object_map ["relax.test"])) % 0;
		selection.io.report_dnds (relax.inferred_distribution);


		io.ReportProgressMessageMD("RELAX", "alt", "* The following rate distribution was inferred for **reference** branches");
		relax.inferred_distribution_ref = parameters.GetStickBreakingDistribution (models.codon.BS_REL.ExtractMixtureDistribution (relax.model_object_map ["relax.reference"])) % 0;
		selection.io.report_dnds (relax.inferred_distribution_ref);

        relax.take1_snapshot = estimators.TakeLFStateSnapshot(relax.alternative_model.fit[terms.likelihood_function]);

		relax.lf.raw = relax.ComputeOnGrid  ( 
											  relax.alternative_model.fit[terms.likelihood_function],
											  relax.grid.MatrixToDict ({200,1}["_MATRIX_ELEMENT_ROW_*0.025"]),
											  "relax.pass1.evaluator",
											  "relax.pass1.result_handler"
											);
											

		// FIND the difference between K < 1 and K > 1

		relax.best_samples = {{-1e100,-1e100}};

		for (relax.k = 0; relax.k < 40; relax.k += 1) {
			relax.best_samples[0] = Max (relax.best_samples[0], relax.lf.raw[relax.k]);
		}

		for (relax.k = 40; relax.k < 200; relax.k += 1) {
			relax.best_samples[1] = Max (relax.best_samples[1], relax.lf.raw[relax.k]);
		}


		if (Abs (relax.best_samples[1] - relax.best_samples[0]) < 5.) { // could be diagnostic of convergence problems
		    selection.io.json_store_setting  (relax.json, "convergence-flat-surface", relax.best_samples[1] - relax.best_samples[0]);

		
			io.ReportProgressMessageMD("RELAX", "alt-2", "* Potential convergence issues due to flat likelihood surfaces; checking to see whether K > 1 or K < 1 is robustly inferred");
			
			if (relax.fitted.K > 1) {
				parameters.SetRange (model.generic.GetGlobalParameter (relax.model_object_map ["relax.test"] , terms.relax.k), terms.range01);
			} else {
				parameters.SetRange (model.generic.GetGlobalParameter (relax.model_object_map ["relax.test"] , terms.relax.k), terms.relax.k_range1);
			}
			
			//assert (__SIGTRAP__);
			
						
			relax.alternative_model.fit.take2 =  estimators.FitLF (relax.filter_names, relax.trees, relax.model_map,
																   relax.alternative_model.fit ,
																   relax.model_object_map,
																   {
																    terms.run_options.retain_lf_object: TRUE,
																    terms.run_options.optimization_log :  relax.optimization_log_file("MainALT-redo-log.json")}
																   );

            io.ReportProgressMessageMD("RELAX", "alt2", "* Attempt to fit an alternative direction of K");
 	        io.ReportProgressMessageMD("RELAX", "alt2", "* " + selection.io.report_fit (relax.alternative_model.fit.take2, 9, relax.codon_data_info[terms.data.sample_size]));
           io.ReportProgressMessageMD("RELAX", "alt2", "* Relaxation/intensification parameter (K) = " + Format(estimators.GetGlobalMLE (relax.alternative_model.fit.take2,terms.relax.k),8,2));
            io.ReportProgressMessageMD("RELAX", "alt2", "* The following rate distribution was inferred for **test** branches");
            relax.inferred_distribution2 = parameters.GetStickBreakingDistribution (models.codon.BS_REL.ExtractMixtureDistribution (relax.model_object_map ["relax.test"])) % 0;
            selection.io.report_dnds (relax.inferred_distribution2);
    
    
            io.ReportProgressMessageMD("RELAX", "alt2", "* The following rate distribution was inferred for **reference** branches");
            relax.inferred_distribution_ref2 = parameters.GetStickBreakingDistribution (models.codon.BS_REL.ExtractMixtureDistribution (relax.model_object_map ["relax.reference"])) % 0;
            selection.io.report_dnds (relax.inferred_distribution_ref2);

            selection.io.json_store_setting  (relax.json, "convergence-unstable-alernative", {
                    {estimators.GetGlobalMLE (relax.alternative_model.fit.take2,terms.relax.k), relax.alternative_model.fit.take2 [terms.fit.log_likelihood]}
                });

			if (relax.alternative_model.fit.take2 [terms.fit.log_likelihood] > relax.alternative_model.fit [terms.fit.log_likelihood]) {
			    relax.stash_fitted_K = relax.fitted.K;

				io.ReportProgressMessageMD("RELAX", "alt-2", "\n### Potential for highly unreliable K inference due to multiple local maxima in the likelihood function, treat results with caution ");
				io.ReportProgressMessageMD("RELAX", "alt-2", "> Relaxation parameter reset to opposite mode of evolution from that obtained in the initial optimization.");
				io.ReportProgressMessageMD("RELAX", "alt-2", "* " + selection.io.report_fit (relax.alternative_model.fit.take2, 9, relax.codon_data_info[terms.data.sample_size]));
				relax.fitted.K = estimators.GetGlobalMLE (relax.alternative_model.fit.take2,terms.relax.k);
				io.ReportProgressMessageMD("RELAX", "alt-2", "* Relaxation/intensification parameter (K) = " + Format(relax.fitted.K,8,2));
				io.ReportProgressMessageMD("RELAX", "alt-2", "* The following rate distribution was inferred for **test** branches");
				relax.inferred_distribution = parameters.GetStickBreakingDistribution (models.codon.BS_REL.ExtractMixtureDistribution (relax.model_object_map ["relax.test"])) % 0;
				selection.io.report_dnds (relax.inferred_distribution);

                selection.io.json_store_setting  (relax.json, "convergence-unstable", {
                    {relax.fitted.K, relax.alternative_model.fit.take2 [terms.fit.log_likelihood]}
                    {relax.stash_fitted_K, relax.alternative_model.fit [terms.fit.log_likelihood]}
                });

				io.ReportProgressMessageMD("RELAX", "alt-2", "* The following rate distribution was inferred for **reference** branches");
				relax.inferred_distribution_ref = parameters.GetStickBreakingDistribution (models.codon.BS_REL.ExtractMixtureDistribution (relax.model_object_map ["relax.reference"])) % 0;
				selection.io.report_dnds (relax.inferred_distribution_ref);

				relax.alternative_model.fit = relax.alternative_model.fit.take2;
                io.SpoolLFToPath(relax.alternative_model.fit.take2[terms.likelihood_function], relax.save_fit_path);
                relax.stashLF = estimators.TakeLFStateSnapshot (relax.alternative_model.fit[terms.likelihood_function]);

			} else {
			    estimators.RestoreLFStateFromSnapshot(relax.alternative_model.fit[terms.likelihood_function], relax.take1_snapshot);
			}

            DeleteObject (relax.alternative_model.fit.take2[terms.likelihood_function]);
            DeleteObject (relax.alternative_model.fit.take2);
 
			parameters.SetRange (model.generic.GetGlobalParameter (relax.model_object_map ["relax.test"] , terms.relax.k), terms.relax.k_range);


		}
	
		relax.distribution_for_json = {relax.test_branches_name : utility.Map (utility.Range (relax.rate_classes, 0, 1),
														 "_index_",
														 "{terms.json.omega_ratio : relax.inferred_distribution [_index_][0],
														   terms.json.proportion  : relax.inferred_distribution [_index_][1]}"),

									   relax.reference_branches_name : utility.Map (utility.Range (relax.rate_classes, 0, 1),
														 "_index_",
														 "{terms.json.omega_ratio : relax.inferred_distribution_ref [_index_][0],
														   terms.json.proportion  : relax.inferred_distribution_ref [_index_][1]}")
								   };

	} else {
		relax.report_multi_class_rates (relax.alternative_model.fit, TRUE);
	}
	
	relax.report_multi_hit  (relax.alternative_model.fit, relax.distribution_for_json, "RELAX", "alt-mh");
    relax._report_srv (relax.alternative_model.fit, FALSE);
        
       
	selection.io.json_store_lf (relax.json,
								relax.alternative_name,
								relax.alternative_model.fit[terms.fit.log_likelihood],
								relax.alternative_model.fit[terms.parameters] + 9 , // +9 comes from CF3x4
								relax.codon_data_info[terms.data.sample_size],
								relax.distribution_for_json,
								relax.display_orders[relax.alternative_name]
							);

	selection.io.json_store_branch_attribute(relax.json, relax.alternative_name, terms.branch_length, relax.display_orders[relax.alternative_name],
												 0,
												 selection.io.extract_branch_info((relax.alternative_model.fit[terms.branch_length])[0], "selection.io.branch.length"));

	selection.io.stopTimer (relax.json [terms.json.timers], "RELAX alternative model fitting");

	// NULL MODEL

	selection.io.startTimer (relax.json [terms.json.timers], "RELAX null model fitting", 4);

	io.ReportProgressMessageMD ("RELAX", "null", "Fitting the null (K := 1) model");
    
    
	for (relax.k = 1; relax.k < relax.numbers_of_tested_groups; relax.k += 1) {
		relax.model_nmsp = relax.model_namespaces[relax.k ];
		if (relax.k > 1) {
			parameters.SetConstraint (model.generic.GetGlobalParameter (relax.model_object_map[relax.model_nmsp] , relax.relax_parameter_terms[relax.k]), terms.parameters.one, terms.global);
		} else {
			parameters.SetConstraint (model.generic.GetGlobalParameter (relax.model_object_map[relax.model_nmsp] , terms.relax.k), terms.parameters.one, terms.global);
		}
	}
	
   /// START EBF CALCULATION
   ///=====================================================================================================================================
	
  utility.ToggleEnvVariable ("KEEP_OPTIMAL_ORDER", TRUE);
        
        relax.er_report.tagged_branches = {};
        relax.er_report.tagged_sites    = {};
        
        for (_partition_, _selection_; in; relax.selected_branches) {
        
            busted.branch_level_ER = {};
            busted.branch_site_level_ER = {};
            busted.current_weights = {};
            
            for (_key_, _value_; in; busted.mixture_weights_parameters) {
                busted.current_weights  [_value_] = Eval (_key_);
            }
            
            busted.tested_branches = {};
            busted.full_to_short = {};
            for (_b_,_class_; in; _selection_) {
                if (_class_ == terms.tree_attributes.test) {
                     busted.node_name = busted.tree_ids[+_partition_] + '.' + _b_;
                     busted.tested_branches [busted.node_name] = 1;
                     busted.full_to_short  [busted.node_name] = _b_;
                } else {
                    busted.branch_level_ER  [_b_] = None;
                }
            }
             
            busted.tested_branches = Rows (busted.tested_branches);    
            
            utility.ToggleEnvVariable ("SET_MODEL_KEEP_LOCALS", TRUE);
            SetParameter ( busted.tested_branches , MODEL, ^busted.er_model);
            utility.ToggleEnvVariable ("SET_MODEL_KEEP_LOCALS", None);
            
           
            for (_b_; in;  busted.tested_branches) {
                for (_key_, _value_; in; busted.mixture_weights_parameters) {
                    ^(_b_ + "." + _value_ ) =  busted.current_weights [_value_];
                }
            }
            
            busted.weight_matrix = {1,busted.rate_classes-1};
            
            for (_key_, _value_; in; busted.mixture_weights_parameters) {
                ^_value_ =^_key_;
            }
            
            for (busted.i = 1; busted.i < busted.rate_classes; busted.i += 1)  {
                busted.weight_matrix[busted.i-1] = Eval (busted.branch_weight_prefix + (busted.i));
            }
            
            busted.weight_matrix = parameters.GetStickBreakingDistributionWeigths (busted.weight_matrix);
            
            io.ReportProgressBar("", "Computing branch and site level posterior probabilities");
            LFCompute (^(busted.full_model[terms.likelihood_function]),LF_START_COMPUTE);
            
            //GetString (exp_counter, MATRIX_EXPONENTIALS_COMPUTED,0);
            //console.log ("\n\n" + exp_counter + "\n\n");
            
            for (_b_; in;  busted.tested_branches) {
    
                io.ReportProgressBar("", "Profiling partition " + (+_partition_ + 1) + "/branch " + busted.full_to_short [_b_]);
                
                // set all weight parameters to 0; this will yield the weight for the N-1 class
                
                busted.branchPP = {busted.rate_classes,1};
                for (_key_, _value_; in; busted.mixture_weights_parameters) {
                    ^(_b_ + "." + _value_ ) =  0;
                }
                
                
                LFCompute (^(busted.full_model[terms.likelihood_function]),logl);
                ConstructCategoryMatrix (busted.siteLL, ^(busted.full_model[terms.likelihood_function]) , SITE_LOG_LIKELIHOODS, {{+_partition_}});
                busted.all_siteLL = {busted.rate_classes, Columns (busted.siteLL)};
                busted.all_siteLL [busted.rate_classes-1][0] = busted.siteLL;
                busted.branchPP [busted.rate_classes-1] = logl;
                
                for (busted.i = 1; busted.i < busted.rate_classes; busted.i += 1) {
                    if (busted.i > 1) {
                        ^(_b_ + "." + busted.branch_weight_prefix + (busted.i-1)) = 0;
                    }
                    ^(_b_ + "." + busted.branch_weight_prefix + (busted.i)) = 1;
                    LFCompute (^(busted.full_model[terms.likelihood_function]),logl);
                    ConstructCategoryMatrix (busted.siteLL, ^(busted.full_model[terms.likelihood_function]) , SITE_LOG_LIKELIHOODS, {{+_partition_}});
                    busted.all_siteLL [busted.i-1][0] = busted.siteLL;
                    busted.branchPP [busted.i-1] = logl;
                }
                
                for (_key_, _value_; in; busted.mixture_weights_parameters) {
                    ^(_b_ + "." + _value_ ) =  busted.current_weights [_value_];
                }
                
                busted.branch_level_ER  [busted.full_to_short [_b_]]      = busted.mixture_site_logl (busted.branchPP,busted.weight_matrix );
                busted.branch_site_level_ER  [busted.full_to_short [_b_]] = busted.mixture_site_logl (busted.all_siteLL,busted.weight_matrix );
                
                busted.branch_site_level_BF = (busted.mixture_site_BF (busted.branch_site_level_ER  [busted.full_to_short [_b_]], busted.weight_matrix))[busted.rate_classes - 1][-1];
                busted.tagged_sites         = busted.branch_site_level_BF ["_MATRIX_ELEMENT_VALUE_>busted.site_BF_reporting"];
                busted.tagged_site_count    = +busted.tagged_sites;
                
                if (busted.tagged_site_count ) {
                    busted.er_report.tagged_branches [_partition_ + ":" + _b_] = busted.tagged_site_count;
                    for (_site_;in;(busted.tagged_sites["(1+_MATRIX_ELEMENT_COLUMN_)*_MATRIX_ELEMENT_VALUE_"])[busted.tagged_sites]) {
                        _site_tag_ = _partition_ + ":" + _site_;
                        if ( busted.er_report.tagged_sites / _site_tag_ == FALSE) {
                             busted.er_report.tagged_sites [_site_tag_] = {};
                        }
                        busted.er_report.tagged_sites [_site_tag_] + _b_;
                    }
                }
                //GetString (exp_counter2, MATRIX_EXPONENTIALS_COMPUTED,0);
                //console.log ("\n\n" + (exp_counter2-exp_counter) + "\n\n");
                exp_counter = exp_counter2;
            
            }	   
            
            //GetString (exp_counter2, MATRIX_EXPONENTIALS_COMPUTED,0);
            //console.log ("\n\n" + (exp_counter2-exp_counter) + "\n\n");
            
            LFCompute (^(busted.full_model[terms.likelihood_function]),LF_DONE_COMPUTE);
            SetParameter ( busted.tested_branches , MODEL, busted.test);
    
            estimators.RestoreLFStateFromSnapshot (busted.full_model[terms.likelihood_function], busted.stashLF);
             
            
            selection.io.json_store_branch_attribute(busted.json, busted.ER, terms.json.branch_annotations, busted.display_orders[busted.ER],
                                                 _partition_,
                                                 busted.branch_level_ER);
                                                 
            selection.io.json_store_branch_attribute(busted.json, busted.bsER, terms.json.branch_annotations, busted.display_orders[busted.bsER],
                                                 _partition_,
                                                 busted.branch_site_level_ER);
                                                 
            //console.log (busted.er_report.tagged_branches);
            //console.log (busted.er_report.tagged_sites);
       }
       io.ClearProgressBar();
       utility.ToggleEnvVariable ("KEEP_OPTIMAL_ORDER", None);
    ///=====================================================================================================================================
	/// END EBF CALCULATION

	relax.null_model.fit = estimators.FitExistingLF (relax.alternative_model.fit[terms.likelihood_function], relax.model_object_map);
	io.ReportProgressMessageMD ("RELAX", "null", "* " + selection.io.report_fit (relax.null_model.fit, 9, relax.codon_data_info[terms.data.sample_size]));
	relax.LRT = math.DoLRT (relax.null_model.fit[terms.fit.log_likelihood], relax.alternative_model.fit[terms.fit.log_likelihood],  relax.numbers_of_tested_groups-1);

	if (relax.numbers_of_tested_groups == 2 && relax.analysis_run_mode != relax.kGroupMode) {
		io.ReportProgressMessageMD("RELAX", "null", "* The following rate distribution for test/reference branches was inferred");
		relax.inferred_distribution = parameters.GetStickBreakingDistribution (models.codon.BS_REL.ExtractMixtureDistributionFromFit (relax.model_object_map ["relax.test"], relax.null_model.fit)) % 0;
		selection.io.report_dnds (relax.inferred_distribution);

		relax.distribution_for_json = {relax.test_branches_name : utility.Map (utility.Range (relax.rate_classes, 0, 1),
															 "_index_",
															 "{terms.json.omega_ratio : relax.inferred_distribution [_index_][0],
															   terms.json.proportion  : relax.inferred_distribution [_index_][1]}")};

		relax.distribution_for_json   [relax.reference_branches_name] =   relax.distribution_for_json   [relax.test_branches_name];
	} else {
		relax.report_multi_class_rates (None, FALSE);
	}
	
	relax.report_multi_hit  (relax.null_model.fit, relax.distribution_for_json, "RELAX", "alt-mh-null");
    relax._report_srv (relax.null_model.fit, TRUE);

	selection.io.json_store_lf (relax.json,
								relax.null_name ,
								relax.null_model.fit[terms.fit.log_likelihood],
								relax.null_model.fit[terms.parameters] + 9 , // +9 comes from CF3x4
								relax.codon_data_info[terms.data.sample_size],
								relax.distribution_for_json,
								relax.display_orders[relax.null_name]
							);

	selection.io.json_store_branch_attribute(relax.json, relax.null_name, terms.branch_length, relax.display_orders[relax.null_name],
												 0,
												 selection.io.extract_branch_info((relax.null_model.fit[terms.branch_length])[0], "selection.io.branch.length"));


	console.log ("----\n## Test for relaxation (or intensification) of selection [RELAX]");
	console.log ( "Likelihood ratio test **p = " + Format (relax.LRT[terms.p_value], 8, 4) + "**.");


	if (relax.LRT[terms.p_value] <= relax.p_threshold) {
		if ( relax.numbers_of_tested_groups == 2 && relax.analysis_run_mode != relax.kGroupMode) {
			if (relax.fitted.K > 1) {
				console.log (">Evidence for *intensification of selection* among **test** branches _relative_ to the **reference** branches at P<="+ relax.p_threshold);
			} else {
				console.log (">Evidence for *relaxation of selection* among **test** branches _relative_ to the **reference** branches at P<="+ relax.p_threshold);
			}
		} else {
			console.log (">Evidence for differences of selective pressures among the test groups at P<="+ relax.p_threshold);
	
		}
	} else {
		if ( relax.numbers_of_tested_groups == 2 && relax.analysis_run_mode != relax.kGroupMode) {
			console.log (">No significant evidence for relaxation (or intensification) of selection among **test** branches _relative_ to the **reference** branches at P<="+ relax.p_threshold);
		} else {
			console.log (">>No significant evidence for differences of selective pressures among the test groups at P<="+ relax.p_threshold);
	
		}
	}

	relax.json [terms.json.test_results] = relax.LRT;
	(relax.json [terms.json.test_results])[terms.relax.k] = relax.fitted.K;

	console.log ("----\n");

	selection.io.stopTimer (relax.json [terms.json.timers], "RELAX null model fitting");
}


//------------------------------------


relax.convergence_loop = TRUE;
relax.loop_passes = 0;

do {
	relax.loop_passes += 1;
	relax.FitMainTestPair (relax.loop_passes == 1);
	relax.do_lhc = FALSE;
	

	if (relax.LRT [terms.LRT] < 0) {
		io.ReportProgressMessageMD("RELAX", "refit", "* Detected convergence issues (negative LRT). Refitting the alterative/null model pair from a new starting point");
		selection.io.json_store_setting  (relax.json, "convergence-negative-lrt", relax.loop_passes);

		relax.general_descriptive.fit = relax.null_model.fit; // reset initial conditions
		for (relax.k = 1; relax.k < relax.numbers_of_tested_groups; relax.k += 1) { // remove constraints on K
			relax.model_nmsp = relax.model_namespaces[relax.k ];
			if (relax.k > 1) {
				parameters.RemoveConstraint (model.generic.GetGlobalParameter (relax.model_object_map[relax.model_nmsp] , relax.relax_parameter_terms[relax.k]));
			} else {
				parameters.RemoveConstraint (model.generic.GetGlobalParameter (relax.model_object_map[relax.model_nmsp] , terms.relax.k));
			}
		}
	} else  {
		relax.convergence_loop = FALSE;
	}	
} while (relax.convergence_loop && relax.loop_passes <= 5);

if (relax.model_set == "All") {
    selection.io.startTimer (relax.json [terms.json.timers], "RELAX partitioned descriptive", 5);

    io.ReportProgressMessageMD ("RELAX", "pe", "Fitting the partitioned descriptive model (completely separate rate distributions for branch sets)");
    parameters.RemoveConstraint (utility.Keys (relax.bound_weights));
	for (relax.k = 1; relax.k < relax.numbers_of_tested_groups; relax.k += 1) {
		relax.model_nmsp = relax.model_namespaces[relax.k ];
		for (relax.i = 1; relax.i <= relax.rate_classes; relax.i += 1) {
        	parameters.RemoveConstraint (model.generic.GetGlobalParameter (relax.model_object_map [relax.model_nmsp] , terms.AddCategory (terms.parameters.omega_ratio,relax.i)));
        }
    }
    relax.pe.fit = estimators.FitExistingLF (relax.alternative_model.fit[terms.likelihood_function], relax.model_object_map);
    io.ReportProgressMessageMD ("RELAX", "pe", "* " + selection.io.report_fit (relax.pe.fit, 9, relax.codon_data_info[terms.data.sample_size]));
    if (relax.numbers_of_tested_groups == 2 && relax.analysis_run_mode != relax.kGroupMode) {
		io.ReportProgressMessageMD ("RELAX", "pe", "* The following rate distribution was inferred for *test* branches ");
		relax.test.inferred_distribution = parameters.GetStickBreakingDistribution (models.codon.BS_REL.ExtractMixtureDistribution(relax.model_object_map ["relax.test"])) % 0;
		selection.io.report_dnds (relax.test.inferred_distribution);
		io.ReportProgressMessageMD("RELAX", "pe", "* The following rate distribution was inferred for *reference* branches ");
		relax.reference.inferred_distribution =  parameters.GetStickBreakingDistribution (models.codon.BS_REL.ExtractMixtureDistribution(relax.model_object_map ["relax.reference"])) % 0;
		selection.io.report_dnds (relax.reference.inferred_distribution);

		relax.distribution_for_json = {relax.test_branches_name : utility.Map (utility.Range (relax.rate_classes, 0, 1),
															 "_index_",
															 "{terms.json.omega_ratio : relax.test.inferred_distribution [_index_][0],
															   terms.json.proportion  : relax.test.inferred_distribution [_index_][1]}"),

										relax.reference_branches_name : utility.Map (utility.Range (relax.rate_classes, 0, 1),
															 "_index_",
															 "{terms.json.omega_ratio : relax.reference.inferred_distribution [_index_][0],
															   terms.json.proportion  : relax.reference.inferred_distribution [_index_][1]}")
									   };
	} else {
		relax.report_multi_class_rates (None, TRUE);
	}
	relax.report_multi_hit (relax.pe.fit, relax.distribution_for_json, "RELAX", "pe-mh");

    selection.io.json_store_lf (relax.json,
                                relax.partitioned_descriptive_name,
                                relax.pe.fit[terms.fit.log_likelihood],
                                relax.pe.fit[terms.parameters] + 9 , // +9 comes from CF3x4
                                relax.codon_data_info[terms.data.sample_size],
                                relax.distribution_for_json,
                                relax.display_orders[relax.partitioned_descriptive_name]
                            );

    selection.io.json_store_branch_attribute(relax.json, relax.partitioned_descriptive_name, terms.branch_length, relax.display_orders[relax.partitioned_descriptive_name],
                                                 0,
                                                 selection.io.extract_branch_info((relax.pe.fit[terms.branch_length])[0], "selection.io.branch.length"));


    selection.io.stopTimer (relax.json [terms.json.timers], "RELAX partitioned descriptive");
}

selection.io.stopTimer (relax.json [terms.json.timers], "Overall");

GetString (_hpv,HYPHY_VERSION,0);
relax.json[terms.json.runtime] = _hpv;


io.SpoolJSON (relax.json, relax.codon_data_info [terms.json.json]);

return relax.json;


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

lfunction relax.extract.k(branch_info) {
    return (branch_info[utility.getGlobalValue("terms.relax.k")])[utility.getGlobalValue("terms.fit.MLE")];
}

//------------------------------------------------------------------------------

lfunction relax.set.k (tree_name, node_name, model_description, ignore) {
    if (utility.Has (model_description [utility.getGlobalValue ("terms.local")], utility.getGlobalValue ("terms.relax.k"), "String")) {
        k = (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.relax.k")];
        t = (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.parameters.synonymous_rate")];
        parameters.SetConstraint (tree_name + "." + node_name + "." + k, utility.getGlobalValue ("relax.relaxation_parameter"), "");
        parameters.SetRange (tree_name + "." + node_name + "." + k, utility.getGlobalValue ("terms.relax.k_range"));
        parameters.SetRange (tree_name + "." + node_name + "." + t, utility.getGlobalValue ("terms.relax.t_range"));
    }
    return tree_name + "." + node_name + "." + k;
}

//------------------------------------------------------------------------------

lfunction relax.set.k2 (tree_name, node_name, model_description, ignore) {
    if (utility.Has (model_description [utility.getGlobalValue ("terms.local")], utility.getGlobalValue ("terms.relax.k"), "String")) {
        k = (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.relax.k")];
        parameters.RemoveConstraint (tree_name + "." + node_name + "." + k);
    }
    return tree_name + "." + node_name + "." + k;
}


//------------------------------------------------------------------------------

lfunction relax.init.k (lf_id, components, data_filter, tree, model_map, initial_values, model_objects) {

    parameter_set = estimators.TraverseLocalParameters (lf_id, model_objects, "relax.set.k");


    return 0;
}

//------------------------------------------------------------------------------

lfunction relax.BS_REL.ModelDescription (type, code, components) {
    model = models.codon.BS_REL.ModelDescription(utility.getGlobalValue ('terms.global'), code, components);
    model [utility.getGlobalValue("terms.model.defineQ")] = "relax.BS_REL._DefineQ";
    return model;
}

//------------------------------------------------------------------------------


lfunction relax.DistributionGuess (mean) {
    rc = utility.getGlobalValue ("relax.rate_classes");

    guess = {rc,2};

    guess[rc-1][0] = 5;
    guess[rc-1][1] = 0.1;

    for (k = 0; k < rc - 1; k += 1) {
        guess[k][0] = 0.1 ^ (1 / (1 + k));
        guess[k][1] = (0.9) / (rc-1) ;
    }

    norm = + guess[-1][1];
    guess_mean = 1/(+(guess [-1][0] $ guess [-1][1]))/norm;
    return guess["_MATRIX_ELEMENT_VALUE_*(guess_mean*(_MATRIX_ELEMENT_COLUMN_==0)+(_MATRIX_ELEMENT_COLUMN_==1)*(1/norm))"];
}

//------------------------------------------------------------------------------

lfunction relax.BS_REL._GenerateRate.MH (fromChar, toChar, namespace, model_type, _tt, alpha, alpha_term, beta, beta_term, omega, omega_term, delta, delta_term, psi, psi_term, psi_s, psi_s_term) {

    p = {};
    _GenerateRate.diff = models.codon.diff.complete(fromChar, toChar);
    diff_count = utility.Array1D (_GenerateRate.diff);

    if (diff_count 
    ) {

        p[model_type] = {};
        p[utility.getGlobalValue("terms.global")] = {};

        nuc_rate = "";

        for (i = 0; i < diff_count; i += 1) {
            if ((_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.from")] > (_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.to")]) {
                nuc_p = "theta_" + (_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.to")] + (_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.from")];
            } else {
                nuc_p = "theta_" + (_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.from")] +(_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.to")];
            }
            nuc_p = parameters.ApplyNameSpace(nuc_p, namespace);
            (p[utility.getGlobalValue("terms.global")])[terms.nucleotideRateReversible((_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.from")], (_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.to")])] = nuc_p;

            nuc_rate = parameters.AppendMultiplicativeTerm (nuc_rate, nuc_p);
       }


        rate_entry = nuc_rate;

        if (_tt[fromChar] != _tt[toChar]) {
            if (model_type == utility.getGlobalValue("terms.global")) {
                aa_rate = parameters.ApplyNameSpace(omega, namespace);
                (p[model_type])[omega_term] = aa_rate;
                utility.EnsureKey (p, utility.getGlobalValue("terms.local"));
                 (p[utility.getGlobalValue("terms.local")])[utility.getGlobalValue ("terms.relax.k")] = "k";
                 aa_rate += "^k";
            } else {
                aa_rate = beta;
                (p[model_type])[beta_term] = aa_rate;
            }
            rate_entry += "*" + aa_rate;
        } else {
            if (model_type == utility.getGlobalValue("terms.local")) {
                (p[model_type])[alpha_term] = alpha;
                rate_entry += "*" + alpha;
            } else {
                p[utility.getGlobalValue("terms.model.rate_entry")] = nuc_rate;
            }


        }

        if (diff_count == 2) {
            if (model_type == utility.getGlobalValue("terms.global")) {
                delta_rate = parameters.ApplyNameSpace(delta, namespace);
                (p[model_type])[delta_term] = delta_rate;
                rate_entry += "*" + delta_rate;
            } else {
                (p[model_type])[delta_term] = delta;
                rate_entry += "*" + delta;
            }
        }

        
        if (diff_count == 3 && ^'relax.multi_hit' == 'Double+Triple') {
            if (_tt[fromChar] != _tt[toChar]) {
                if (model_type == utility.getGlobalValue("terms.global")) {
                    psi_rate = parameters.ApplyNameSpace(psi, namespace);
                    (p[model_type])[psi_term] = psi_rate;
                    rate_entry += "*" + psi_rate;
                } else {
                    (p[model_type])[psi_term] = psi;
                    rate_entry += "*" + psi;
                }
            } else {
               //console.log (fromChar + " <-> " + toChar + " (" + _tt[fromChar] + ")");
               if (model_type == utility.getGlobalValue("terms.global")) {
                    psi_rate = parameters.ApplyNameSpace(psi_s, namespace);
                    (p[model_type])[psi_s_term] = psi_rate;
                    rate_entry += "*" + psi_rate;
                } else {
                    (p[model_type])[psi_s_term] = psi_s;
                    rate_entry += "*" + psi_s;
                }
            
            }
        }


        p[utility.getGlobalValue("terms.model.rate_entry")] = rate_entry;
    }
    return p;
}

//------------------------------------------------------------------------------


lfunction relax.BS_REL._GenerateRate (fromChar, toChar, namespace, model_type, _tt, alpha, alpha_term, beta, beta_term, omega, omega_term) {
    p = {};
    diff = models.codon.diff(fromChar, toChar);

    if (None != diff) {
        p[model_type] = {};
        p[utility.getGlobalValue("terms.global")] = {};

        if (diff[utility.getGlobalValue("terms.diff.from")] > diff[utility.getGlobalValue("terms.diff.to")]) {
            nuc_rate = "theta_" + diff[utility.getGlobalValue("terms.diff.to")] + diff[utility.getGlobalValue("terms.diff.from")];
        } else {
            nuc_rate = "theta_" + diff[utility.getGlobalValue("terms.diff.from")] + diff[utility.getGlobalValue("terms.diff.to")];
        }
        nuc_rate = parameters.ApplyNameSpace(nuc_rate, namespace);
        (p[utility.getGlobalValue("terms.global")])[terms.nucleotideRate(diff[utility.getGlobalValue("terms.diff.from")], diff[utility.getGlobalValue("terms.diff.to")])] = nuc_rate;

        if (_tt[fromChar] != _tt[toChar]) {
            if (model_type == utility.getGlobalValue("terms.global")) {
                aa_rate = parameters.ApplyNameSpace(omega, namespace);
                (p[model_type])[omega_term] = aa_rate;
                utility.EnsureKey (p, utility.getGlobalValue("terms.local"));
                 (p[utility.getGlobalValue("terms.local")])[utility.getGlobalValue ("terms.relax.k")] = "k";
                 aa_rate += "^k";
            } else {
                aa_rate = beta;
                (p[model_type])[beta_term] = aa_rate;
            }
            p[utility.getGlobalValue("terms.model.rate_entry")] = nuc_rate + "*" + aa_rate;
        } else {
            if (model_type == utility.getGlobalValue("terms.local")) {
                (p[model_type])[alpha_term] = alpha;
                p[utility.getGlobalValue("terms.model.rate_entry")] = nuc_rate + "*" + alpha;
            } else {
                p[utility.getGlobalValue("terms.model.rate_entry")] = nuc_rate;
            }
        }
    }
    return p;
}

//------------------------------------------------------------------------------

lfunction relax.BS_REL._DefineQ (bs_rel, namespace) {
    rate_matrices = {};

    bs_rel [utility.getGlobalValue("terms.model.q_ij")] = &rate_generator;
    bs_rel [utility.getGlobalValue("terms.mixture.mixture_components")] = {};

    _aux = parameters.GenerateSequentialNames (namespace + ".bsrel_mixture_aux", bs_rel[utility.getGlobalValue("terms.model.components")] - 1, "_");
    _wts = parameters.helper.stick_breaking (_aux, None);
    mixture = {};

    for (component = 1; component <= bs_rel[utility.getGlobalValue("terms.model.components")]; component += 1) {
       key = "component_" + component;
   
       if (^'relax.multi_hit' == 'None') {
           ExecuteCommands ("
            function rate_generator (fromChar, toChar, namespace, model_type, model) {
               return relax.BS_REL._GenerateRate (fromChar, toChar, namespace, model_type, model[utility.getGlobalValue('terms.translation_table')],
                    'alpha', utility.getGlobalValue('terms.parameters.synonymous_rate'),
                    'beta_`component`', terms.AddCategory (utility.getGlobalValue('terms.parameters.nonsynonymous_rate'), component),
                    'omega`component`', terms.AddCategory (utility.getGlobalValue('terms.parameters.omega_ratio'), component));
                }"
           );
        } else {
            ExecuteCommands ("
            function rate_generator (fromChar, toChar, namespace, model_type, model) {
               return relax.BS_REL._GenerateRate.MH (fromChar, toChar, namespace, model_type, model[utility.getGlobalValue('terms.translation_table')],
                    'alpha', utility.getGlobalValue('terms.parameters.synonymous_rate'),
                    'beta_`component`', terms.AddCategory (utility.getGlobalValue('terms.parameters.nonsynonymous_rate'), component),
                    'omega`component`', terms.AddCategory (utility.getGlobalValue('terms.parameters.omega_ratio'), component),
                    'delta', utility.getGlobalValue('terms.parameters.multiple_hit_rate'),
                    'psi', utility.getGlobalValue('terms.parameters.triple_hit_rate'),
                    'psi', utility.getGlobalValue('terms.parameters.triple_hit_rate'));
                }"
           );
        }
       
       if ( component < bs_rel[utility.getGlobalValue("terms.model.components")]) {
            model.generic.AddGlobal ( bs_rel, _aux[component-1], terms.AddCategory (utility.getGlobalValue("terms.mixture.mixture_aux_weight"), component ));
            parameters.DeclareGlobalWithRanges (_aux[component-1], 0.5, 0, 1);
       } else {

       }
       models.codon.generic.DefineQMatrix(bs_rel, namespace);
       rate_matrices [key] = bs_rel[utility.getGlobalValue("terms.model.rate_matrix")];
       (bs_rel [^'terms.mixture.mixture_components'])[key] = _wts [component-1];
    }


    bs_rel[utility.getGlobalValue("terms.model.rate_matrix")] = rate_matrices;
    parameters.SetConstraint(((bs_rel[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.global")])[terms.nucleotideRate("A", "G")], "1", "");
    
    return bs_rel;
}



//------------------------------------------------------------------------------
lfunction relax.select_branches(partition_info) {

    kGroupMode = ^"relax.kGroupMode";

    tree_count = utility.Array1D (partition_info);
    //io.CheckAssertion(" == 1", "RELAX only works on a single partition dataset");
    available_models = {};

    return_set         = {};
    tree_configuration = {};

    for (i,p; in; partition_info) {
        tree_for_analysis = p[utility.getGlobalValue("terms.data.tree")];
        tree_configuration [i] = {};
        if (+i == 0) {
            for (_value_; in; tree_for_analysis[utility.getGlobalValue("terms.trees.model_map")]) {
                available_models [_value_] = 1;
            }
        } else {
            partition_models = {};
            for (_value_; in; tree_for_analysis[utility.getGlobalValue("terms.trees.model_map")]) {
                partition_models [_value_] = 1;
            }
            for (k, v; in; partition_models) {
                available_models[k] += 1;
            }
        }
    }
    
    list_models = {};
    for (i,p; in; available_models) {
        if (p == tree_count) {
            list_models[i] = 1;
        }
    }
    
     
    list_models   = utility.sortStrings(utility.Keys(list_models)); // get keys
    option_count  = utility.Array1D (list_models);
    can_run_group_mode = (option_count == Abs (available_models));
    
    io.CheckAssertion	("`&option_count` >= 2", "RELAX requires at least one designated set of branches in the tree.");
    
    nontrivial_groups = option_count;
    
    if (utility.Has (available_models, "", "Number")) {
    	nontrivial_groups += -1;
    }
        
    run_mode = None;   
         
	if (nontrivial_groups >= 2 && can_run_group_mode) { // could run as a group set
		run_mode = io.SelectAnOption ({
			{kGroupMode, "Run the test for equality of selective regimes among  " + nontrivial_groups + " groups of branches"}
			{"Classic mode", "Select one test and one reference group of branches, with the rest of the branches treated as unclassified"}
		}, "Group test mode");
		
		if (run_mode == kGroupMode) {
			 utility.SetEnvVariable ("relax.numbers_of_tested_groups", nontrivial_groups);
			 
			 for (i,p; in; partition_info) {		
			     tree_for_analysis = p[utility.getGlobalValue("terms.data.tree")];	 
                 for (_key_, _value_; in; tree_for_analysis[utility.getGlobalValue("terms.trees.model_map")]) {
                    if ('' == _value_ ) {
                        (tree_configuration[i])[_key_] = utility.getGlobalValue('relax.unclassified_branches_name');
                    } else {
                        (tree_configuration[i])[_key_] = _value_;
                    }
                 }			 
                }
			 utility.SetEnvVariable ("relax.analysis_run_mode", kGroupMode);
		}
	}
	
	// console.log ("**** `option_count` ****");
	
	if (run_mode != kGroupMode) {

		selectTheseForTesting = {
			option_count, 2
		};

		for (k = 0; k < option_count; k += 1) {
			if (list_models[k] != "") {
				selectTheseForTesting[k][0] = list_models[k];
				selectTheseForTesting[k][1] = "Set " + list_models[k] + " with " + available_models[list_models[k]] + " branches";
			} else {
				selectTheseForTesting[k][0] = "Unlabeled branches";
				selectTheseForTesting[k][1] = "Set of " + available_models[list_models[k]] + " unlabeled branches";
			}
		}

		ChoiceList(testSet, "Choose the set of branches to use as the _test_ set", 1, NO_SKIP, selectTheseForTesting);
		io.CheckAssertion ("`&testSet` >= 0", "User cancelled branch selection; analysis terminating");
		if (option_count > 2) {
			ChoiceList(referenceSet, "Choose the set of branches to use as the _reference_ set", 1, testSet, selectTheseForTesting);
			io.CheckAssertion ("`&referenceSet` >= 0", "User cancelled branch selection; analysis terminating");
		} else {
			referenceSet = 1-testSet;
		}


		tag_test = selectTheseForTesting [testSet][0];
		if (tag_test == "Unlabeled branches") {
			tag_test = "";
		}
		tag_reference = selectTheseForTesting [referenceSet][0];
		if (tag_reference == "Unlabeled branches") {
			tag_reference = "";
		}


        for (i,p; in; partition_info) {		
		    tree_for_analysis = p[utility.getGlobalValue("terms.data.tree")];	 
		    for (_key_,_value_; in; tree_for_analysis[utility.getGlobalValue("terms.trees.model_map")]) {
		        if (tag_test == _value_ ) {
				    (tree_configuration[i])[_key_] = utility.getGlobalValue('relax.test_branches_name');
                } else {
                    if (tag_reference == _value_ ) {
                        (tree_configuration[i])[_key_] = utility.getGlobalValue('relax.reference_branches_name');
                    } else {
                        (tree_configuration[i])[_key_] = utility.getGlobalValue('relax.unclassified_branches_name');
                    }
                }
            }
		}
	}

    return tree_configuration;
}

//------------------------------------------------------------------------------

lfunction relax.grid.MatrixToDict (grid) {
    return utility.Map (utility.MatrixToListOfRows (grid), "_value_",
                                                                '{  terms.relax.k : {
                                                                            terms.id : relax.relaxation_parameter,
                                                                            terms.fit.MLE : _value_[1]
                                                                        }

                                                                 }');
}

//------------------------------------------------------------------------------

function relax.init_grid_setup (omega_distro) {
    utility.ForEachPair (omega_distro[terms.parameters.rates], "_index_", "_name_", 
        '
            if (_index_[0] < relax.rate_classes - 1) { // not the last rate
                  relax.initial_ranges [_name_] = {
                    terms.lower_bound : 0,
                    terms.upper_bound : 1
                };
            }  else {
                relax.initial_ranges [_name_] = {
                    terms.lower_bound : 1,
                    terms.upper_bound : 5
                };
            }
        '
    );


    utility.ForEachPair (omega_distro[terms.parameters.weights], "_index_", "_name_", 
        '
             relax.initial_ranges [_name_] = {
                terms.lower_bound : 0,
                terms.upper_bound : 1
            };
        '
    );

}

//------------------------------------------------------------------------------

function relax.init_grid_setup_scaled (omega_distro) {
    utility.ForEachPair (omega_distro[terms.parameters.rates], "_index_", "_name_", 
        '
              relax.initial_ranges [_name_] = {
                terms.lower_bound : 0,
                terms.upper_bound : 1
             };
            
        '
    );


    utility.ForEachPair (omega_distro[terms.parameters.weights], "_index_", "_name_", 
        '
             relax.initial_ranges [_name_] = {
                terms.lower_bound : 0,
                terms.upper_bound : 1
            };
        '
    );

}

//------------------------------------------------------------------------------

lfunction relax._renormalize (v, distro, mean) {

    parameters.SetValues (v);
    m = parameters.GetStickBreakingDistribution (^distro);
    d = Rows (m);
    m = +(m[-1][0] $ m[-1][1]); // current mean
    for (i = 0; i < d; i+=1) {
        (v[((^distro)["rates"])[i]])[^"terms.fit.MLE"] = (v[((^distro)["rates"])[i]])[^"terms.fit.MLE"] / m * mean;
    }
    return v;
    
}

//------------------------------------------------------------------------------

lfunction relax._renormalize_with_weights (v, distro, mean) {


    parameters.SetValues (v);
    m = parameters.GetStickBreakingDistribution (^distro);
    //console.log (v);
    //console.log (m);
    //console.log (mean);
    d = Rows (m);
    mean = Max (mean, 1e-3);
    
    over_one = m[d-1][0] * m[d-1][1];
    
    if (over_one >= mean*0.95) {
       //console.log ("OVERAGE");
       new_weight = mean * Random (0.9, 0.95) / m[d-1][0];
       diff = (m[d-1][1] - new_weight)/(d-1);
       for (k = 0; k < d-1; k += 1) {
            m[k][1] += diff;
       }
       m[d-1][1] = new_weight;
    }
    
    over_one = m[d-1][0] * m[d-1][1];
    under_one = (+(m[-1][0] $ m[-1][1])) / (1-m[d-1][1]); // current mean
    
    for (i = 0; i < d-1; i+=1) {
        m[i][0] = m[i][0] * mean / under_one;
    }
  
    under_one = +(m[-1][0] $ m[-1][1]);
    
    if (under_one > 0) {
        for (i = 0; i < d; i+=1) {
            m[i][0] = m[i][0] * mean / under_one;
        }
    }
        
    m = m%0;
    
    wts = parameters.SetStickBreakingDistributionWeigths (m[-1][1]);
 

    for (i = 0; i < d; i+=1) {
        (v[((^distro)["rates"])[i]])[^"terms.fit.MLE"] = m[i][0];
    }
    for (i = 0; i < d-1; i+=1) {
        (v[((^distro)["weights"])[i]])[^"terms.fit.MLE"] = wts[i];
    }
    
    //console.log (v);

    //assert (0);
    return v;
    
}


//------------------------------------------------------------------------------
function relax.optimization_log_file (extension) {
    if (relax.OPTIMIZATION_LOGS) {
        return relax.codon_data_info [terms.json.json] + extension;
    }   
    return None;
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

lfunction relax.report_multi_hit (model_fit, json, l1, l2) {
    if (^'relax.multi_hit' != "None") {
        io.ReportProgressMessageMD(l1, l2, 'Partition-level rates for multiple-hit substitutions\n');
        params = selection.io.extract_global_MLE_re (model_fit, ^'terms.parameters.multiple_hit_rate');
        
        double_rates = {};
        for (mle; in; params) {
            exp = regexp.FindSubexpressions (mle[^'terms.description'], "\[relax\.(.+)\]");
            if (None == exp) {
                exp = 'reference';
            } else {
                exp = exp[1];
            }
            double_rates [exp] = mle[^'terms.fit.MLE'];
            io.ReportProgressMessageMD(l1, l2, '* 2H rate for **' + exp + '**: ' + Format (double_rates [exp], 8, 4));
            
        }
        
        json[^'terms.parameters.multiple_hit_rate'] = double_rates; 
        
        if (^'relax.multi_hit' != 'Double') {
            params = selection.io.extract_global_MLE_re (model_fit, ^'terms.parameters.triple_hit_rate');
        
            double_rates = {};
            for (mle; in; params) {
                exp = regexp.FindSubexpressions (mle[^'terms.description'], "\[relax\.(.+)\]");
                if (None == exp) {
                    exp = 'reference';
                } else {
                    exp = exp[1];
                }
                double_rates [exp] = mle[^'terms.fit.MLE'];
                io.ReportProgressMessageMD(l1, l2, '* 3H rate for **' + exp + '**: ' + Format (double_rates [exp], 8, 4));
    
            }

            json[^'terms.parameters.triple_hit_rate'] = double_rates;            
        }       
    }
}