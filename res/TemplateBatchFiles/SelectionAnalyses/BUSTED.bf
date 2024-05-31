RequireVersion ("2.5.51");


LoadFunctionLibrary("libv3/all-terms.bf"); // must be loaded before CF3x4
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/ancestral.bf");
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
//utility.SetEnvVariable ("LF_SMOOTHING_SCALER", 0.1);


busted.analysis_description = {

                               terms.io.info : "BUSTED (branch-site unrestricted statistical test of episodic diversification) uses a random effects branch-site model fitted  jointly to all or a subset of tree branches in order to test for alignment-wide evidence of episodic diversifying selection. Assuming there is evidence of positive selection (i.e. there is an omega > 1),  BUSTED will also perform a quick evidence-ratio style analysis to explore which individual sites may have been subject to selection. v2.0 adds support for synonymous rate variation, and relaxes the test statistic to 0.5 (chi^2_0 + chi^2_2).
                               
Version 2.1 adds a grid search for the initial starting point.

Version 2.2 changes the grid search to LHC, and adds an initial search phase to use adaptive Nedler-Mead. 

Version 3.0 implements the option for branch-site variation in synonymous substitution rates. 

Version 3.1 adds HMM auto-correlation option for SRV, and binds SRV distributions for multiple branch sets.

Version 4.0 adds support for multiple hits (MH), ancestral state reconstruction saved to JSON, and profiling of branch-site level support for selection / multiple hits.

Version 4.2 adds calculation of MH-attributable fractions of substitutions.

Version 4.5 adds an 'error absorption' component [experimental]
",
                               terms.io.version : "4.5",
                               terms.io.reference : "*Gene-wide identification of episodic selection*, Mol Biol Evol. 32(5):1365-71, *Synonymous Site-to-Site Substitution Rate Variation Dramatically Inflates False Positive Rates of Selection Analyses: Ignore at Your Own Peril*, Mol Biol Evol. 37(8):2430-2439",
                               terms.io.authors : "Sergei L Kosakovsky Pond",
                               terms.io.contact : "spond@temple.edu",
                               terms.settings: {},
                               terms.io.requirements : "in-frame codon alignment and a phylogenetic tree (optionally annotated with {})"
                              };


io.DisplayAnalysisBanner (busted.analysis_description);

busted.FG = "Test";
busted.BG = "Background";
busted.SRV = "Synonymous site-to-site rates";
busted.background = "background";
busted.unconstrained = "unconstrained";
busted.constrained   = "constrained";
busted.optimized_null = "optimized null";
busted.ER = "Posterior prob omega class";
busted.bsER = "Posterior prob omega class by site";
busted.ER2H = "Evidence ratio for 2H";
busted.ER3H = "Evidence ratio for 3H";
busted.ER23H = "Evidence ratio for 2H+3H";
busted.MG94 = terms.json.mg94xrev_sep_rates;

busted.json.background = busted.background;
busted.json.site_logl  = "Site Log Likelihood";
busted.json.evidence_ratios  = "Evidence Ratios";
busted.json.srv_posteriors  = "Synonymous site-posteriors";
busted.json.srv_viterbi = "Viterbi synonymous rate path";
busted.rate_classes = 3;
busted.synonymous_rate_classes = 3;
busted.initial_grid.N = 250;
busted.site_BF_reporting = 100;


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
                         busted.constrained: 3 ,
                         busted.ER : 4,
                         busted.bsER : 5,
                         busted.ER2H : 6,
                         busted.ER3H : 7,
                         busted.ER23H : 8                    
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
                                    
selection.io.json_store_setting  (busted.json, "srv", busted.do_srv);

                                    
if (busted.do_srv == "Branch-site") {
    busted.do_bs_srv = TRUE;
    busted.do_srv = TRUE;
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

KeywordArgument ("multiple-hits",  "Include support for multiple nucleotide substitutions", "None");

busted.multi_hit = io.SelectAnOption ({
                                        {"Double", "Include branch-specific rates for double nucleotide substitutions"}
                                        {"Double+Triple", "Include branch-specific rates for double and triple nucleotide substitutions"}
                                        {"None", "[Default] Use standard models which permit only single nucleotide changes to occur instantly"}
                                  }, "Include support for multiple nucleotide substitutions");


selection.io.json_store_setting  (busted.json, "multiple-hit", busted.multi_hit);

if (busted.do_srv) {
    KeywordArgument ("syn-rates", "The number alpha rate classes to include in the model [1-10, default 3]", busted.synonymous_rate_classes);
    busted.synonymous_rate_classes = io.PromptUser ("The number alpha rate classes to include in the model [1-10, default 3]", busted.synonymous_rate_classes, 1, 10, TRUE);
}

KeywordArgument ("error-sink",  "[Advanced experimental setting] Include a rate class to capture misalignment artifacts", "No");

busted.error_sink = io.SelectAnOption ({
                                        {"No", "Standard dN/dS model (all rates are used for inference)"}
                                        {"Yes", "The dN/dS model has an additional class (high dN/dS, low proportion), which 'absorbs' alignment errors"}
                                  }, "[Advanced experimental setting] Include a rate class to capture misalignment artifacts") == "Yes";

selection.io.json_store_setting  (busted.json, "error-sink", busted.error_sink);

if (busted.error_sink) {
    busted.rate_classes += 1;
}

KeywordArgument ("grid-size", "The number of points in the initial distributional guess for likelihood fitting", 250);
busted.initial_grid.N = io.PromptUser ("The number of points in the initial distributional guess for likelihood fitting", 250, 1, 10000, TRUE);
KeywordArgument ("starting-points", "The number of initial random guesses to seed rate values optimization", 1);
busted.N.initial_guesses = io.PromptUser ("The number of initial random guesses to 'seed' rate values optimization", 1, 1, busted.initial_grid.N$10, TRUE);
                                    
KeywordArgument ("output", "Write the resulting JSON to this file (default is to save to the same path as the alignment file + 'BUSTED.json')", busted.codon_data_info [terms.json.json]);
busted.codon_data_info [terms.json.json] = io.PromptUserForFilePath ("Save the resulting JSON file to");
    
io.ReportProgressMessageMD('BUSTED',  'selector', 'Branches to test for selection in the BUSTED analysis');

for (_partition_, _selection_; in; busted.selected_branches) {
    _selection_ = utility.Filter (_selection_, '_value_', '_value_ == terms.tree_attributes.test');
    io.ReportProgressMessageMD('BUSTED',  'selector', '* Selected ' + Abs(_selection_) + ' branches to test in the BUSTED analysis: \`' + Join (', ',utility.Keys(_selection_)) + '\`');
}


// check to see if there are any background branches
busted.has_background = FALSE;

for (_partition_, _selection_; in; busted.selected_branches) {
     _selection_ = utility.Filter (_selection_, '_value_', '_value_ != terms.tree_attributes.test');
     if (utility.Array1D (_selection_)) { busted.has_background = TRUE; break;}
}

busted.json[busted.json.background] =  busted.has_background;
selection.io.startTimer (busted.json [terms.json.timers], "Preliminary model fitting", 1);

namespace_tag = "busted";

namespace busted {
    doGTR ("busted");
}


estimators.fixSubsetOfEstimates(busted.gtr_results, busted.gtr_results[terms.global]);

namespace busted {
    scaler_prefix = "busted.scaler";
    doPartitionedMG ("busted", FALSE);
}


busted.run_full_mg94 = TRUE;
    
if (Type (busted.save_intermediate_fits) == "AssociativeList") {
    if (None != busted.save_intermediate_fits[^"terms.data.value"]) {
        if (utility.Has (busted.save_intermediate_fits[^"terms.data.value"], "Full-MG94", "AssociativeList")) {
            busted.final_partitioned_mg_results = (busted.save_intermediate_fits[^"terms.data.value"])["Full-MG94"];
            busted.run_full_mg94 = FALSE;
        }        
    }
}

io.ReportProgressMessageMD ("BUSTED", "codon-refit", "Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model");

if (busted.run_full_mg94) {        
    busted.final_partitioned_mg_results = estimators.FitMGREV (busted.filter_names, busted.trees, busted.codon_data_info [terms.code], {
            terms.run_options.model_type: terms.local,
            terms.run_options.partitioned_omega: busted.selected_branches,
            terms.run_options.apply_user_constraints: busted.zero_branch_length_constrain,
            terms.run_options.optimization_settings: {
                "OPTIMIZATION_METHOD" : "hyrbid"
            }
        }, busted.partitioned_mg_results);

    if (Type (busted.save_intermediate_fits) == "AssociativeList") {
        (busted.save_intermediate_fits[^"terms.data.value"])["Full-MG94"] = busted.final_partitioned_mg_results;        
        Export (lfe, ^busted.final_partitioned_mg_results[^"terms.likelihood_function"]);
        (busted.save_intermediate_fits[^"terms.data.value"])["Full-MG94-LF"] = lfe;
        io.SpoolJSON (busted.save_intermediate_fits[^"terms.data.value"],busted.save_intermediate_fits[^"terms.data.file"]);      
    }
}



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



if (busted.multi_hit == "None") {
    busted.model_generator = "models.codon.BS_REL.ModelDescription";
   
} else {
   lfunction busted.model.BS_REL_MH (type, code, rates) {        
        def = models.codon.BS_REL.ModelDescription (type, code, rates);
        def [utility.getGlobalValue("terms.model.defineQ")] = "models.codon.BS_REL._DefineQ.MH";
        if (^'busted.multi_hit' == "Double") {
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
    busted.model_generator = "busted.model.BS_REL_MH";
}

busted.baseline_model_generator = busted.model_generator;
if (busted.do_srv) {
    if (busted.do_bs_srv) {
        busted.model_generator = "models.codon.BS_REL_SRV.ModelDescription";
        busted.rate_class_arguments = {{busted.synonymous_rate_classes__,busted.rate_classes__}};
        assert (busted.multi_hit == "None", "Multiple hit and branch site SRV combination is currently not supported");
        assert (busted.error_sink  == FALSE, "Error sink and branch site SRV combination is currently not supported");
    } else {
    
        lfunction busted.model.with.GDD (type, code, rates) {        
            def = Call (^"busted.baseline_model_generator", type, code, rates);
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
        
        

busted.mixture_weights_parameters = {};
busted.branch_weight_prefix = "branch_level_weight_";
busted.g2l = {};

for (busted.i = 1; busted.i < busted.rate_classes; busted.i += 1) {
    busted.global_parameter = model.generic.GetGlobalParameter (busted.test.bsrel_model , terms.AddCategory (utility.getGlobalValue("terms.mixture.mixture_aux_weight"), busted.i));
    busted.mixture_weights_parameters [busted.global_parameter] = busted.branch_weight_prefix  + busted.i;
    busted.g2l [busted.global_parameter] = 1;
}

        
function busted.export_er_model () {        
    Export (busted.er_omegaN, busted.test, {"RETURN_MODEL_NAME" : TRUE, 
                                            "UNIQUE_COMPONENT_NAMES" : TRUE, 
                                            "VARIABLE_SUBSTITUTIONS" : busted.mixture_weights_parameters,
                                            "GLOBAL_TO_LOCAL" :  busted.g2l
                                            });
    ExecuteCommands (busted.er_omegaN);  
    //console.log (busted.er_omegaN);    branch_level_er_calculator = 1;
}

busted.er_model = busted.export_er_model();

busted._mh_parameters = {};

if (busted.multi_hit != "None") {
    busted.double_hit_parameter = model.generic.GetGlobalParameter (busted.test.bsrel_model , ^'terms.parameters.multiple_hit_rate');
    busted.triple_hit_parameter = model.generic.GetGlobalParameter (busted.test.bsrel_model , ^'terms.parameters.triple_hit_rate');
    busted._mh_to_local   = {};
    busted._mh_parameter_map = {};
    if (None != busted.double_hit_parameter) {
        busted._mh_parameters [busted.double_hit_parameter] = "branch_level_dh";
        busted._mh_to_local [busted.double_hit_parameter] = 1;
        busted._mh_parameter_map [^'terms.parameters.multiple_hit_rate'] = busted._mh_parameters [busted.double_hit_parameter];
        
    }
    if (None != busted.triple_hit_parameter) {
        busted._mh_parameters [busted.triple_hit_parameter] = "branch_level_th";
        busted._mh_to_local [busted.triple_hit_parameter] = 1;
        busted._mh_parameter_map [^'terms.parameters.triple_hit_rate'] = busted._mh_parameters [busted.triple_hit_parameter];
    }
    
    function busted.export_mh_er_model () {        
        Export (busted.er_mh, busted.test, {"RETURN_MODEL_NAME" : TRUE, 
                                                "UNIQUE_COMPONENT_NAMES" : TRUE, 
                                                "VARIABLE_SUBSTITUTIONS" : busted._mh_parameters,
                                                "GLOBAL_TO_LOCAL" : busted._mh_to_local 
                                                });
                                                
        ExecuteCommands (busted.er_mh);  
    }    
    
    busted.mh_er_model = busted.export_mh_er_model();
}


busted.background.bsrel_model =  model.generic.DefineMixtureModel(busted.model_generator,
        "busted.background", {
            "0": parameters.Quote(terms.global),
            "1": busted.codon_data_info[terms.code],
            "2": parameters.Quote (busted.rate_class_arguments) // the number of rate classes
        },
        busted.filter_names,
        None);


models.BindGlobalParameters ({"0" : busted.test.bsrel_model, "1" : busted.background.bsrel_model}, terms.nucleotideRate("[ACGT]","[ACGT]"));

if (busted.do_srv) {
    if (busted.do_bs_srv) {
        for (description,var_id; in;  (busted.background.bsrel_model[terms.parameters])[terms.global]) {
            if (regexp.Find (description, terms.parameters.synonymous_rate + " for category")) {
                var_id_2 = utility.GetByKey ((busted.test.bsrel_model[terms.parameters])[terms.global], description, "String");
                if (None != var_id_2) {
                   parameters.SetConstraint (var_id, var_id_2, terms.global);
                }
	        }
		}

        models.BindGlobalParameters ({"0" : busted.test.bsrel_model, "1" : busted.background.bsrel_model}, terms.AddCategory (utility.getGlobalValue("terms.mixture.mixture_aux_weight"), "SRV [0-9]+"));
    } 
}

if (busted.multi_hit != "None" && busted.has_background) {
     models.BindGlobalParameters ({"0" : busted.test.bsrel_model, "1" : busted.background.bsrel_model}, ^'terms.parameters.multiple_hit_rate');
     models.BindGlobalParameters ({"0" : busted.test.bsrel_model, "1" : busted.background.bsrel_model}, ^'terms.parameters.triple_hit_rate');
}

busted.distribution = models.codon.BS_REL.ExtractMixtureDistribution(busted.test.bsrel_model);


// set up parameter constraints

for (busted.i = 1; busted.i < busted.rate_classes; busted.i += 1) {
    parameters.SetRange (model.generic.GetGlobalParameter (busted.test.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,busted.i)), terms.range01);
    parameters.SetRange (model.generic.GetGlobalParameter (busted.background.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,busted.i)), terms.range01);
}

if (busted.error_sink) {
     parameters.SetRange (model.generic.GetGlobalParameter (busted.test.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,1)), terms.range_high);
     parameters.SetRange (model.generic.GetGlobalParameter (busted.background.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,1)), terms.range_high);
     parameters.SetRange  (model.generic.GetGlobalParameter (busted.test.bsrel_model , terms.AddCategory (utility.getGlobalValue("terms.mixture.mixture_aux_weight"), 1)),terms.range_small_fraction);
     parameters.SetRange  (model.generic.GetGlobalParameter (busted.background.bsrel_model , terms.AddCategory (utility.getGlobalValue("terms.mixture.mixture_aux_weight"), 1)),terms.range_small_fraction);
}

busted.omega_eds_parameter = terms.AddCategory (terms.parameters.omega_ratio,busted.rate_classes);
busted.omega_eds_w_parameter = terms.AddCategory (utility.getGlobalValue("terms.mixture.mixture_aux_weight"), busted.rate_classes-1);

parameters.SetRange (model.generic.GetGlobalParameter (busted.test.bsrel_model , busted.omega_eds_parameter), terms.range_gte1);

busted.branch_length_string = busted.test.bsrel_model [terms.model.branch_length_string];
busted.model_parameters = busted.test.bsrel_model[terms.parameters];


/* create an initial grid to initialize the optimization */
busted.initial_grid = {};

/* if populated, use this as a baseline to generate the distributions from */
busted.initial_grid_presets = {"0" : 0.1};
busted.initial_ranges = {};

busted.init_grid_setup (busted.distribution, busted.error_sink, busted.rate_classes);

/** setup parameter optimization groups */

PARAMETER_GROUPING = {};
//PARAMETER_GROUPING + busted.distribution["rates"];
//PARAMETER_GROUPING + busted.distribution["weights"];

PARAMETER_GROUPING + utility.Concat (busted.distribution["rates"], busted.distribution["weights"]);

if (busted.has_background) {
    busted.model_object_map = { "busted.background" : busted.background.bsrel_model,
                                "busted.test" :       busted.test.bsrel_model };
    busted.background_distribution = models.codon.BS_REL.ExtractMixtureDistribution(busted.background.bsrel_model);
    busted.init_grid_setup (busted.background_distribution, busted.error_sink, busted.rate_classes);

    //PARAMETER_GROUPING + busted.background_distribution["rates"];
    //PARAMETER_GROUPING + busted.background_distribution["weights"];
    PARAMETER_GROUPING + utility.Concat (busted.background_distribution["rates"], busted.background_distribution["weights"]);
    
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
    
 
    //PARAMETER_GROUPING + busted.srv_distribution["rates"];
    //PARAMETER_GROUPING + busted.srv_distribution["weights"];
    PARAMETER_GROUPING + utility.Concat (busted.srv_distribution["rates"],busted.srv_distribution["weights"]);

    busted.init_grid_setup (busted.srv_distribution, FALSE, busted.synonymous_rate_classes);    
}

if (busted.do_srv_hmm) {
    busted.hmm_lambda = model.generic.GetGlobalParameter (busted.test.bsrel_model, terms.rate_variation.hmm_lambda);
    //parameters.SetConstraint (((busted.test.bsrel_model[terms.parameters])[terms.global])[terms.rate_variation.hmm_lambda],"1e-6", "");
    busted.initial_ranges [busted.hmm_lambda] = {
                terms.lower_bound : 0.05,
                terms.upper_bound : 0.95
            };
}

busted.global_scaler_list = {};

for (busted.partition_index = 0; busted.partition_index < busted.partition_count; busted.partition_index += 1) {
    busted.global_scaler_list [busted.partition_index] = "busted.bl.scaler_" + busted.partition_index;
    parameters.DeclareGlobalWithRanges (busted.global_scaler_list [busted.partition_index], 1, 0, 1000);
   // busted.initial_ranges      [busted.global_scaler_list [busted.partition_index]] = {
   //             terms.lower_bound : 1.0,
   //             terms.upper_bound : 3.0
   //         };
}


busted.initial.test_mean    = ((selection.io.extract_global_MLE_re (busted.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio + ".+test.+"))["0"])[terms.fit.MLE];

busted.initial_grid         = estimators.LHC (busted.initial_ranges,busted.initial_grid.N);

busted.initial_grid = utility.Map (busted.initial_grid, "_v_", 
    'busted._renormalize_with_weights (_v_, "busted.distribution", busted.initial.test_mean, busted.error_sink)'
);

if (busted.has_background) { //has background
    busted.initial.background_mean    = ((selection.io.extract_global_MLE_re (busted.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio + ".+background.+"))["0"])[terms.fit.MLE];
    
    busted.initial_grid = utility.Map (busted.initial_grid, "_v_", 
        'busted._renormalize_with_weights (_v_, "busted.background_distribution", busted.initial.background_mean, busted.error_sink)'
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

  
busted.nm.precision = Max (-0.0001*busted.final_partitioned_mg_results[terms.fit.log_likelihood],0.5);

debug.checkpoint = utility.GetEnvVariable ("DEBUG_CHECKPOINT");

busted.converged = FALSE;
busted.restart_optimization = None;

busted.use_cached_full_model = FALSE;
    
if (Type (busted.save_intermediate_fits) == "AssociativeList") {
    if (None != busted.save_intermediate_fits[^"terms.data.value"]) {
        if (utility.Has (busted.save_intermediate_fits[^"terms.data.value"], "BUSTED-alt", "AssociativeList")) {
            busted.full_model_res =  busted.final_partitioned_mg_results;
            busted.full_model_res * (busted.save_intermediate_fits[^"terms.data.value"])["BUSTED-alt"];
            // check to see which *global* variables have no initial values
            busted.missing = {};
            busted.cache_set = {};
            
            busted.use_model_prefix = {};
            busted.variable_overlap = {};
            
            for (model_name, m; in; busted.model_object_map) {
                for (d, p; in; (m[terms.parameters])[terms.global]) {
                    busted.variable_overlap[d] += 1;
                    busted.check_id = "[`model_name`] `d`";
                    if (busted.full_model_res[terms.global] / busted.check_id != 0) {
                        busted.use_model_prefix [model_name] += 1;
                    } 
                }
            }
            
            function busted._get_cache_name (model_name, parameter) {
                if (busted.variable_overlap[parameter] > 1) {
                    if (busted.use_model_prefix [model_name] > 0) {
                        return "[" + model_name + "] " + parameter;
                    }
                }
                return parameter;
            }
            
            for (model_name, m; in; busted.model_object_map) {
                for (d, p; in; (m[terms.parameters])[terms.global]) {
                    busted.check_id = busted._get_cache_name (model_name, d);
                    
                    if (busted.full_model_res[terms.global] / busted.check_id == 0) {
                        busted.missing[busted.check_id] = p;
                    } else {
                        busted.cache_set[p] = ((busted.full_model_res[terms.global])[busted.check_id])[terms.fit.MLE];
                    }
                }
            }
            
             
            if (Abs (busted.missing) > 0) {
                if (busted.error_sink) {
                    if (Abs (busted.missing) == 2 * (1 + busted.has_background)) {
                        if (busted.missing /  busted._get_cache_name ("busted.test", busted.omega_eds_parameter) && busted.missing/busted._get_cache_name ("busted.test", busted.omega_eds_w_parameter)) {

                            if (busted.has_background) {
                                busted.use_cached_full_model = busted.missing /  busted._get_cache_name ("busted.background", busted.omega_eds_parameter) && busted.missing/busted._get_cache_name ("busted.background", busted.omega_eds_w_parameter);
                                busted.use_cached_full_model = TRUE;
                            }  else {
                                busted.use_cached_full_model = TRUE;
                            }                         

                            if (busted.use_cached_full_model) {
                                for (d, k; in; busted.missing) {
                                    (busted.full_model_res[terms.global])[d] = {
                                        terms.id : k,
                                        terms.fit.MLE : 0.
                                    };
                                }
                                
                                
                                function busted._swap_estimates (model_name) {
                                    for (k = busted.rate_classes; k > 1; k+=(-1)) {
                                        busted.cat_k   = busted._get_cache_name (model_name, terms.AddCategory (terms.parameters.omega_ratio,k));
                                        busted.cat_k_1 =  busted._get_cache_name (model_name, terms.AddCategory (terms.parameters.omega_ratio,k-1));
                                        
                                    
                                        ((busted.full_model_res[terms.global])[busted.cat_k])[terms.fit.MLE] = 
                                            ((busted.full_model_res[terms.global])[busted.cat_k_1])[terms.fit.MLE];
                                            
                                        busted.cache_set [((busted.full_model_res[terms.global])[busted.cat_k])[terms.id]] =
                                            busted.cache_set [((busted.full_model_res[terms.global])[busted.cat_k_1])[terms.id]];     
                                    }
                                    
                                    for (k = busted.rate_classes-1; k > 1; k+=(-1)) {
                                        busted.cat_k   = busted._get_cache_name (model_name, terms.AddCategory (terms.mixture.mixture_aux_weight,k));
                                        busted.cat_k_1 =  busted._get_cache_name (model_name, terms.AddCategory (terms.mixture.mixture_aux_weight,k-1));
                                    
                                        ((busted.full_model_res[terms.global])[busted.cat_k])[terms.fit.MLE] = 
                                            ((busted.full_model_res[terms.global])[busted.cat_k_1])[terms.fit.MLE];
                                            
                                        busted.cache_set [((busted.full_model_res[terms.global])[busted.cat_k])[terms.id]] =
                                            busted.cache_set [((busted.full_model_res[terms.global])[busted.cat_k_1])[terms.id]];     
                                    }
                                }
                                
                                function busted._delete_cat_1 (model_name) {
                                    busted.cat_k = busted._get_cache_name (model_name, terms.AddCategory (terms.parameters.omega_ratio,1));
                                    busted.cat_k_1 = busted._get_cache_name (model_name, terms.AddCategory (terms.mixture.mixture_aux_weight,1));
                                    
                                    busted.cache_set - ((busted.full_model_res[terms.global])[busted.cat_k])[terms.id];
                                    busted.cache_set - ((busted.full_model_res[terms.global])[busted.cat_k_1])[terms.id];
        
                                    //busted.full_model_res[terms.global] - busted.cat_k;
                                    //busted.full_model_res[terms.global] - busted.cat_k_1;
                                }
                                
                                busted._swap_estimates ("busted.test");
                                busted._delete_cat_1  ("busted.test");
                                
                                if (busted.has_background) {
                                    busted._swap_estimates ("busted.background");
                                    busted._delete_cat_1  ("busted.background");
                                }
     
                             }
                               
                        }
 
                }
              }
            } else {
                busted.use_cached_full_model = TRUE;
            }
        }
    }
}

while (!busted.converged) {
       
    if (Type (debug.checkpoint) != "String") {
    
        // constrain nucleotide rate parameters
        // copy global nucleotide parameter estimates to busted.test.bsrel_model)
        
        if (None == busted.restart_optimization) {
        
            if (busted.use_cached_full_model) {
                for (k,d;in;busted.initial_grid) {
                    if (+k < 5 ) {
                        for (v,mle; in; d) {
                                if (busted.cache_set / v) {
                                    ((busted.initial_grid[k])[v])[terms.fit.MLE] = busted.cache_set[v];
                                } else {
                                    ((busted.initial_grid[k])[v])[terms.fit.MLE] = Random (0.0,1e-6);
                                }
                            }
            
                    } else {
                        if (Random (0,1) < 0.6) {
                            for (v,mle; in; d) {
                                if (busted.cache_set / v) {
                                    ((busted.initial_grid[k])[v])[terms.fit.MLE] = busted.cache_set[v];
                                } 
                            }
                        }
                    }
                }
                                
                busted.full_model_grid =  estimators.FitLF (busted.filter_names, busted.trees, busted.model_map, busted.full_model_res, busted.model_object_map, {
                        "retain-lf-object": TRUE,
                        terms.run_options.optimization_settings : 
                            {
                                "OPTIMIZATION_METHOD" : "nedler-mead",
                                "MAXIMUM_OPTIMIZATION_ITERATIONS" : 500,
                                "OPTIMIZATION_PRECISION" : busted.nm.precision
                            } ,
                                                     
                        terms.search_grid     : busted.initial_grid,
                        terms.search_restarts : busted.N.initial_guesses
                                                
                });
                
                busted.full_model =  estimators.FitLF (busted.filter_names, busted.trees, busted.model_map, busted.full_model_grid, busted.model_object_map, {
                        "retain-lf-object": TRUE,
                        terms.run_options.optimization_settings : 
                            {
                                "OPTIMIZATION_METHOD" : "hybrid",
                                //"OPTIMIZATION_PRECISION" : 1.
                            } 
                                                
                    });
                
                
                busted.use_cached_full_model = FALSE;
            } else {
            
                busted.tmp_fixed = models.FixParameterSetRegExpFromReference (terms.nucleotideRatePrefix,busted.test.bsrel_model, busted.final_partitioned_mg_results[terms.global]);
                
                  
                busted.grid_search.results =  estimators.FitLF (busted.filter_names, busted.trees, busted.model_map, busted.final_partitioned_mg_results, busted.model_object_map,                    
                    {
                        terms.run_options.retain_lf_object: TRUE,
                        terms.run_options.proportional_branch_length_scaler : 
                                                                busted.global_scaler_list,
                                                            
                        terms.run_options.optimization_settings : 
                            {
                                "OPTIMIZATION_METHOD" : "nedler-mead",
                                "MAXIMUM_OPTIMIZATION_ITERATIONS" : 500,
                                "OPTIMIZATION_PRECISION" : busted.nm.precision
                            } ,
                                                     
                        terms.search_grid     : busted.initial_grid,
                        terms.search_restarts : busted.N.initial_guesses
                    }
                );
                
               
                
                
                parameters.RemoveConstraint (busted.tmp_fixed );
                //console.log (busted.tmp_fixed);
                PARAMETER_GROUPING + Columns (busted.tmp_fixed);
            
                //PRODUCE_OPTIMIZATION_LOG        = 1;
                
                                                            
                busted.full_model =  estimators.FitLF (busted.filter_names, busted.trees, busted.model_map, busted.grid_search.results, busted.model_object_map, {
                        "retain-lf-object": TRUE,
                        terms.run_options.optimization_settings : 
                            {
                                "OPTIMIZATION_METHOD" : "hybrid",
                                //"OPTIMIZATION_PRECISION" : 1.
                            } 
                                                
                    });
            }
       } else {
            busted.full_model =  estimators.FitLF (busted.filter_names, busted.trees, busted.model_map, busted.restart_optimization, busted.model_object_map, {
                    "retain-lf-object": TRUE,
                    terms.run_options.optimization_settings : 
                        {
                            "OPTIMIZATION_METHOD" : "hybrid",
                            //"OPTIMIZATION_PRECISION" : 1.
                        } 
                                            
                });
       }
            
            
       if (Type (busted.save_intermediate_fits) == "AssociativeList") {  
            if (!busted.error_sink) {
               (busted.save_intermediate_fits[^"terms.data.value"])["BUSTED-alt"] = busted.full_model;
               io.SpoolJSON (busted.save_intermediate_fits[^"terms.data.value"],busted.save_intermediate_fits[^"terms.data.file"]);   
            }  
       }       
     
    } else {
        ExecuteAFile (debug.checkpoint);
        GetString (lf, LikelihoodFunction, 0);
        busted.full_model = estimators.ExtractMLEs (lf, busted.model_object_map);
        busted.full_model[terms.likelihood_function] = lf;
        busted.full_model [terms.fit.log_likelihood] = estimators.ComputeLF(lf);
        
    }
    
    // recover ancestral states
    
    
    busted.json [^"terms.substitutions"] = {};
    
    for (_partition_, _selection_; in; busted.selected_branches) {
        busted.ancestral_info = ancestral.build (busted.full_model[terms.likelihood_function],0+_partition_,FALSE);
        (busted.json [^"terms.substitutions"] )[_partition_] = ancestral.ComputeCompressedSubstitutions (busted.ancestral_info);
        DeleteObject (busted.ancestral_info);
    }
    
    
    
    KeywordArgument ("save-fit", "Save BUSTED model fit to this file (default is not to save)", "/dev/null");
    io.SpoolLFToPath(busted.full_model[terms.likelihood_function], io.PromptUserForFilePath ("Save BUSTED model fit to this file ['/dev/null' to skip]"));
    
    
    io.ReportProgressMessageMD("BUSTED", "main", "* " + selection.io.report_fit (busted.full_model, 9, busted.codon_data_info[terms.data.sample_size]));
    
    //VERBOSITY_LEVEL = 101;
    
    
    io.ReportProgressMessageMD("BUSTED", "main", "* For *test* branches, the following rate distribution for branch-site combinations was inferred");
    
    selection.io.stopTimer (busted.json [terms.json.timers], "Unconstrained BUSTED model fitting");
    
    busted.inferred_test_distribution_raw = parameters.GetStickBreakingDistribution (busted.distribution);
    busted.inferred_test_distribution = busted.inferred_test_distribution_raw  % 0;
    
    busted.has_selection_support = busted.inferred_test_distribution_raw[Rows(busted.inferred_test_distribution_raw)-1][-1];
    busted.has_selection_support = busted.has_selection_support[0] * busted.has_selection_support[1] > 1e-8;
    
    busted.distribution_for_json = {busted.FG : utility.Map (utility.Range (busted.rate_classes, 0, 1),
                                                             "_index_",
                                                             "{terms.json.omega_ratio : busted.inferred_test_distribution_raw [_index_][0],
                                                               terms.json.proportion : busted.inferred_test_distribution_raw [_index_][1]}")
                                    };
                                    
    
    
    
    busted.tree_ids = estimators.LFObjectGetTrees (busted.full_model[terms.likelihood_function]);
    busted.EFV_ids = estimators.LFObjectGetEFV (busted.full_model[terms.likelihood_function]);
    
    (busted.json [busted.json.site_logl])[busted.unconstrained] = busted.ComputeSiteLikelihoods (busted.full_model[terms.likelihood_function]);
    
    
                                    
    busted.report_multi_hit  (busted.full_model, busted.distribution_for_json, "MultiHit", "alt-mh", busted.branch_length_string, busted.model_parameters);
    
    if (busted.error_sink) {
        selection.io.report_dnds_with_sink (busted.inferred_test_distribution, busted.inferred_test_distribution_raw[0][0], busted.inferred_test_distribution_raw[0][1]);
    } else {
        selection.io.report_dnds (busted.inferred_test_distribution);
    }
    
    if (busted.has_background) {
        io.ReportProgressMessageMD("BUSTED", "main", "* For *background* branches, the following rate distribution for branch-site combinations was inferred");
        busted.inferred_background_distribution_raw = parameters.GetStickBreakingDistribution (busted.background_distribution);
        busted.inferred_background_distribution  =  busted.inferred_background_distribution_raw % 0;
        if (busted.error_sink) {
            selection.io.report_dnds_with_sink (busted.inferred_background_distribution, busted.inferred_background_distribution_raw[0][0], busted.inferred_background_distribution_raw[0][1]);
        } else {
            selection.io.report_dnds (busted.inferred_background_distribution);
        }
        
        busted.distribution_for_json [busted.BG] = utility.Map (utility.Range (busted.rate_classes, 0, 1),
                                                             "_index_",
                                                             "{terms.json.omega_ratio : busted.inferred_background_distribution_raw [_index_][0],
                                                               terms.json.proportion : busted.inferred_background_distribution_raw [_index_][1]}");
    }
    
    if (busted.do_srv) {
        
        if (busted.do_srv_hmm) {
            busted.hmm_lambda = selection.io.extract_global_MLE (busted.full_model, terms.rate_variation.hmm_lambda);
            busted.hmm_lambda.CI = parameters.GetProfileCI(((busted.full_model[terms.global])[terms.rate_variation.hmm_lambda])[terms.id],
                                        busted.full_model[terms.likelihood_function], 0.95);
            io.ReportProgressMessageMD("BUSTED", "main", "* HMM switching rate = " +  Format (busted.hmm_lambda, 8, 3));
    
            io.ReportProgressMessageMD ("BUSTED", "main", "* HMM switching rate = " + Format (busted.hmm_lambda,8,4) + 
                            " (95% profile CI " + Format ((busted.hmm_lambda.CI )[terms.lower_bound],8,4) + "-" + Format ((busted.hmm_lambda.CI )[terms.upper_bound],8,4) + ")");
                        
            busted.distribution_for_json [terms.rate_variation.hmm_lambda] = busted.hmm_lambda.CI;
        }
    
    
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
    
    busted.stashLF = estimators.TakeLFStateSnapshot (busted.full_model[terms.likelihood_function]);
    
    if (busted.has_selection_support || busted.error_sink) {
    
        utility.ToggleEnvVariable ("KEEP_OPTIMAL_ORDER", TRUE);
        
        busted.er_report.tagged_branches = {};
        busted.er_report.tagged_sites    = {};
        
        for (_partition_, _selection_; in; busted.selected_branches) {
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
        
    } 
    
    
    if (utility.Array1D (busted._mh_parameters)) {
    
        utility.ToggleEnvVariable ("KEEP_OPTIMAL_ORDER", TRUE);
        
        busted.do2H = (^busted.double_hit_parameter) > 0;
        busted.do3H = FALSE;
        
        if (busted._mh_parameter_map / ^'terms.parameters.triple_hit_rate') {
            busted.do3H = (^busted.triple_hit_parameter) > 0;
        }
        
        if (busted.do2H || busted.do3H ) {
        
            for (_partition_, _selection_; in; busted.selected_branches) {
                busted.branch_site_level_ER = {};
            
            
                busted.tested_branches = {};
                busted.full_to_short = {};
                for (_b_,_class_; in; _selection_) {
                    if (_class_ == terms.tree_attributes.test) {
                         busted.node_name = busted.tree_ids[+_partition_] + '.' + _b_;
                         busted.tested_branches [busted.node_name] = 1;
                         busted.full_to_short  [busted.node_name] = _b_;
                    } 
                }
             
                busted.tested_branches = Rows (busted.tested_branches);
                utility.ToggleEnvVariable ("SET_MODEL_KEEP_LOCALS", TRUE);
                SetParameter ( busted.tested_branches , MODEL, ^busted.mh_er_model);
                utility.ToggleEnvVariable ("SET_MODEL_KEEP_LOCALS", None);
            
           
                for (_b_; in;  busted.tested_branches) {
                    for (_key_, _value_; in; busted._mh_parameters) {
                        ^(_b_ + "." + _value_ ) = ^_key_;
                    }
                }
            
            
                io.ReportProgressBar("", "Computing ER for multiple-hits at branch/site level");
    
    
                LFCompute (^(busted.full_model[terms.likelihood_function]),LF_START_COMPUTE);
                
                busted.branch2H = {};
                busted.branch3H = {};
                busted.branch23H = {};
                
            
                for (_b_; in;  busted.tested_branches) {
    
                    io.ReportProgressBar("", "Profiling partition " + (+_partition_ + 1) + "/branch " + busted.full_to_short [_b_]);
                
                    // set all weight parameters to 0; this will yield the weight for the N-1 class
                
                
                    LFCompute (^(busted.full_model[terms.likelihood_function]),logl);
                    ConstructCategoryMatrix (busted.siteLL, ^(busted.full_model[terms.likelihood_function]) , SITE_LOG_LIKELIHOODS, {{+_partition_}});
                
                
                    if (busted.do2H) {
                        busted.er_parameter =  (_b_ + "." + busted._mh_parameter_map[^'terms.parameters.multiple_hit_rate']);
                        ^busted.er_parameter  = 0;
                        ConstructCategoryMatrix (busted.siteLLConstrained, ^(busted.full_model[terms.likelihood_function]) , SITE_LOG_LIKELIHOODS, {{+_partition_}});
                        ^busted.er_parameter  = ^busted.double_hit_parameter;
                        busted.branch2H [busted.full_to_short [_b_]] = busted.FilteredEvidenceRatios (busted.siteLL, busted.siteLLConstrained, 1);
                    }
                    if (busted.do3H) {
                        busted.er_parameter =  (_b_ + "." + busted._mh_parameter_map[^'terms.parameters.triple_hit_rate']);
                        ^busted.er_parameter  = 0;
                        ConstructCategoryMatrix (busted.siteLLConstrained, ^(busted.full_model[terms.likelihood_function]) , SITE_LOG_LIKELIHOODS, {{+_partition_}});
                        busted.branch3H [busted.full_to_short [_b_]] = busted.FilteredEvidenceRatios (busted.siteLL, busted.siteLLConstrained, 1);
                        if (busted.do2H) {
                            busted.er2_parameter =  (_b_ + "." + busted._mh_parameter_map[^'terms.parameters.multiple_hit_rate']);
                            ^busted.er2_parameter  = 0;
                            ConstructCategoryMatrix (busted.siteLLConstrained, ^(busted.full_model[terms.likelihood_function]) , SITE_LOG_LIKELIHOODS, {{+_partition_}});
                            ^busted.er2_parameter  = ^busted.double_hit_parameter;
                            busted.branch23H [busted.full_to_short [_b_]] = busted.FilteredEvidenceRatios (busted.siteLL, busted.siteLLConstrained, 1);
                        }
                    
                        ^busted.er_parameter  = ^busted.triple_hit_parameter;
                    }
                }	   
            
            
                LFCompute (^(busted.full_model[terms.likelihood_function]),LF_DONE_COMPUTE);
    
                SetParameter ( busted.tested_branches , MODEL, busted.test);
                estimators.RestoreLFStateFromSnapshot (busted.full_model[terms.likelihood_function], busted.stashLF);
                
                if (busted.do2H) {
                    selection.io.json_store_branch_attribute(busted.json, busted.ER2H, terms.json.branch_annotations, busted.display_orders[busted.ER],
                                                         _partition_,
                                                         busted.branch2H);
                }
                if (busted.do3H) {
                    selection.io.json_store_branch_attribute(busted.json, busted.ER3H, terms.json.branch_annotations, busted.display_orders[busted.ER],
                                                         _partition_,
                                                         busted.branch3H);
                }
                if (busted.do2H && busted.do3H) {
                    selection.io.json_store_branch_attribute(busted.json, busted.ER23H, terms.json.branch_annotations, busted.display_orders[busted.ER],
                                                         _partition_,
                                                         busted.branch23H);
                }
                                                 
            }
       }
       io.ClearProgressBar();
       utility.ToggleEnvVariable ("KEEP_OPTIMAL_ORDER", None);
        
    } 
    
    
    
    selection.io.json_store_lf (busted.json,
                                "Unconstrained model",
                                busted.full_model[terms.fit.log_likelihood],
                                busted.full_model[terms.parameters] + 9 , // +9 comes from CF3x4
                                busted.codon_data_info[terms.data.sample_size],
                                busted.distribution_for_json,
                                busted.display_orders[busted.unconstrained]);
    
    
    if (busted.error_sink) {
        busted.run_test = busted.inferred_test_distribution_raw [busted.rate_classes-1][0] > 1 && busted.inferred_test_distribution_raw [busted.rate_classes-1][1] > 0;
    } else {
        busted.run_test = busted.inferred_test_distribution [busted.rate_classes-1][0] > 1 && busted.inferred_test_distribution [busted.rate_classes-1][1] > 0;
    }
    
    for (_key_, _value_; in; busted.filter_specification) {
        selection.io.json_store_branch_attribute(busted.json, busted.unconstrained, terms.branch_length, 0,
                                                 _key_,
                                                 selection.io.extract_branch_info((busted.full_model[terms.branch_length])[_key_], "selection.io.branch.length"));
    }
    
    busted.max_site = None;
    
    if (!busted.run_test) {
        io.ReportProgressMessageMD ("BUSTED", "Results", "No evidence for episodic diversifying positive selection under the unconstrained model, skipping constrained model fitting");
        busted.json [terms.json.test_results] = busted.ComputeLRT (0, 0);
        busted.converged = TRUE;
    } else {
        selection.io.startTimer (busted.json [terms.json.timers], "Constrained BUSTED model fitting", 3);
    
    
        io.ReportProgressMessageMD ("BUSTED", "test", "Performing the constrained (dN/dS > 1 not allowed) model fit");
        parameters.SetConstraint (model.generic.GetGlobalParameter (busted.test.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,busted.rate_classes)), terms.parameters.one, terms.global);
        (busted.json [busted.json.site_logl])[busted.constrained] = busted.ComputeSiteLikelihoods (busted.full_model[terms.likelihood_function]);
        busted.null_results = estimators.FitExistingLF (busted.full_model[terms.likelihood_function], busted.model_object_map);
        (busted.json [busted.json.site_logl])[busted.optimized_null] = busted.ComputeSiteLikelihoods (busted.full_model[terms.likelihood_function]);
        
        busted.site_ll_diff = Transpose (((busted.json [busted.json.site_logl])[busted.unconstrained] - (busted.json [busted.json.site_logl])[busted.optimized_null])*2);
        busted.site_ll_sum = +busted.site_ll_diff;
        
        busted.site_ll_diff = ({Rows (busted.site_ll_diff), 2})["(_MATRIX_ELEMENT_COLUMN_==0)*busted.site_ll_diff[_MATRIX_ELEMENT_ROW_]+(_MATRIX_ELEMENT_COLUMN_==1)_MATRIX_ELEMENT_ROW_"] % 0;
        
        busted.max_site = busted.site_ll_diff [Rows (busted.site_ll_diff)-1][-1];
        
        if (busted.max_site[0] >= 0.8 * busted.site_ll_sum ) {
            busted.max_site = busted.max_site [1] + 1;
        } else {
            busted.max_site = None;
        }
        
        io.ReportProgressMessageMD ("BUSTED", "test", "* " + selection.io.report_fit (busted.null_results, 9, busted.codon_data_info[terms.data.sample_size]));
        busted.LRT = busted.ComputeLRT (busted.full_model[terms.fit.log_likelihood], busted.null_results[terms.fit.log_likelihood]);
        busted.json [terms.json.test_results] = busted.LRT;
    
    
        selection.io.stopTimer (busted.json [terms.json.timers], "Constrained BUSTED model fitting");
    
        utility.ForEachPair (busted.filter_specification, "_key_", "_value_",
            'selection.io.json_store_branch_attribute(busted.json, busted.constrained, terms.branch_length, 1,
                                                     _key_,
                                                     selection.io.extract_branch_info((busted.null_results[terms.branch_length])[_key_], "selection.io.branch.length"));');
    
        io.ReportProgressMessageMD("BUSTED", "test", "* For *test* branches under the null (no dN/dS > 1 model), the following rate distribution for branch-site combinations was inferred");
        busted.inferred_test_distribution_raw = parameters.GetStickBreakingDistribution (busted.distribution);
        busted.inferred_test_distribution = busted.inferred_test_distribution_raw  % 0;
    
        busted.distribution_for_json = {busted.FG : utility.Map (utility.Range (busted.rate_classes, 0, 1),
                                                             "_index_",
                                                             "{terms.json.omega_ratio : busted.inferred_test_distribution_raw [_index_][0],
                                                               terms.json.proportion : busted.inferred_test_distribution_raw [_index_][1]}")};
    
        busted.report_multi_hit  (busted.null_results, busted.distribution_for_json, "MultiHit", "null-mh",busted.branch_length_string, busted.model_parameters);
        busted.null_distro_raw = parameters.GetStickBreakingDistribution (busted.distribution);
        
        if (busted.error_sink) {
            selection.io.report_dnds_with_sink (busted.null_distro_raw % 0, busted.null_distro_raw[0][0], busted.null_distro_raw[0][1]);
        } else {
            selection.io.report_dnds (busted.null_distro_raw % 0);
        }
    
    
        if (busted.has_background) {
            busted.inferred_background_distribution_raw = parameters.GetStickBreakingDistribution (busted.background_distribution);
            //selection.io.report_dnds (busted.inferred_background_distribution);
            busted.distribution_for_json [busted.BG] = utility.Map (utility.Range (busted.rate_classes, 0, 1),
                                                                 "_index_",
                                                                 "{terms.json.omega_ratio : busted.inferred_background_distribution_raw [_index_][0],
                                                                   terms.json.proportion : busted.inferred_background_distribution_raw [_index_][1]}");
        }
    
        if (busted.do_srv) {
            if (busted.do_bs_srv) {
                busted.srv_info = parameters.GetStickBreakingDistribution ( busted.srv_rate_reporting) % 0;
            } else {
                busted.srv_info = Transpose((rate_variation.extract_category_information(busted.test.bsrel_model))["VALUEINDEXORDER"][0])%0;
            }
            if (busted.do_srv_hmm) {
                busted.hmm_lambda = selection.io.extract_global_MLE (busted.full_model, terms.rate_variation.hmm_lambda);
                io.ReportProgressMessageMD("BUSTED", "main", "* HMM switching rate = " +  Format (busted.hmm_lambda, 8, 3));
                busted.distribution_for_json [terms.rate_variation.hmm_lambda] = busted.hmm_lambda;
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
        
        if (busted.LRT[terms.LRT] < -0.01) {
            parameters.RemoveConstraint (model.generic.GetGlobalParameter (busted.test.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,busted.rate_classes)));
            console.log ("----\n## Negative test LRT (convergence problem). Refitting the alternative model using the null model as a start.\n----\n");
            busted.restart_optimization = busted.null_results;
        } else {
            busted.converged = TRUE;
        }
    }
}


console.log ("----\n## Branch-site unrestricted statistical test of episodic diversification [BUSTED]");
console.log ( "Likelihood ratio test for episodic diversifying positive selection, **p = " + Format ((busted.json [terms.json.test_results])[terms.p_value], 8, 4) + "**.");

if ((busted.json [terms.json.test_results])[terms.p_value] < 0.1) {
    if (busted.max_site) {
          fprintf(stdout, "\n-------\n", io.FormatLongStringToWidth(
                  ">[WARNING] Most of the statistical signal for episodic diversifying selection in this alignment is derived from a single codon site (** " + busted.max_site  + "**).
                  This could be a sign of possible data quality issues, or outsized influence of a few substitutions, especially if they involve replacing multiple nucleotides along a short branch.
                  You may want to examine the alignment at this site using BUSTED visualization tools, performing model-averaged inference, or rerunning the alignment with data at that site masked to confirm robustness of the result.", 72), "\n");

    }
}

selection.io.stopTimer (busted.json [terms.json.timers], "Overall");

GetString (_hpv,HYPHY_VERSION,0);
busted.json[terms.json.runtime] = _hpv;

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
lfunction busted.FilteredEvidenceRatios (ha, h0, level) {
    h0 = ha["Exp(_MATRIX_ELEMENT_VALUE_-h0[_MATRIX_ELEMENT_COLUMN_])"];
    selected_indices = h0["_MATRIX_ELEMENT_VALUE_>level"];
    if (+selected_indices) {
        m1 = (h0["_MATRIX_ELEMENT_COLUMN_*selected_indices[_MATRIX_ELEMENT_COLUMN_]"])[selected_indices];
        m2 = h0[selected_indices];
        return ({
            Columns (m1), 2
        }["m1[_MATRIX_ELEMENT_ROW_]*(_MATRIX_ELEMENT_COLUMN_==0)+m2[_MATRIX_ELEMENT_ROW_]*(_MATRIX_ELEMENT_COLUMN_==1)"]);
    }
    return {{}};
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


lfunction busted._renormalize (v, distro, mean, skip_first) {
    parameters.SetValues (v);
    m = parameters.GetStickBreakingDistribution (^distro);
    d = Rows (m);
    if (skip_first) {
        m0 = m[0][0]*m[0][1];
    } else {
        m0 = 0;
    }
    m = +(m[-1][0] $ m[-1][1]) -m0; // current mean
    
    for (i = (skip_first != 0); i < d; i+=1) {
        (v[((^distro)["rates"])[i]])[^"terms.fit.MLE"] = (v[((^distro)["rates"])[i]])[^"terms.fit.MLE"] / m * mean;
    }
    return v;
    
}

//------------------------------------------------------------------------------

lfunction busted._renormalize_with_weights (v, distro, mean, skip_first) {
    parameters.SetValues (v);
    m = parameters.GetStickBreakingDistribution (^distro);
    d = Rows (m);

    mean = Max (mean, 1e-3);
   
    if (skip_first) {
        m0 = m[0][0]*m[0][1];
    } else {
        m0 = 0;
    }
    
    
    over_one = m[d-1][0] * m[d-1][1];
    
    if (over_one >= mean*0.95) {
       //console.log ("OVERAGE");
       new_weight = mean * Random (0.9, 0.95) / m[d-1][0];
       diff = (m[d-1][1] - new_weight)/(d-1-(skip_first != 0));
       for (k = (skip_first != 0); k < d-1; k += 1) {
            m[k][1] += diff;
       }
       m[d-1][1] = new_weight;
    }
    
    
    over_one = m[d-1][0] * m[d-1][1];
    
    under_one = (+(m[-1][0] $ m[-1][1]) - m0) / (1-m[d-1][1]); // current mean
        
    for (i = (skip_first != 0); i < d-1; i+=1) {
        m[i][0] = m[i][0] * mean / under_one;
    }

    if (skip_first) {
        m_rest = m [{{1,0}}][{{d-1,1}}];
        under_one = +(m_rest[-1][0] $ m_rest[-1][1]);
    } else {
        under_one = +(m[-1][0] $ m[-1][1]);
    }
    
    
    for (i = (skip_first != 0); i < d; i+=1) {
        m[i][0] = m[i][0] * mean / under_one;
    }
    
   
   
    if (skip_first) {
        m_rest = m [{{1,0}}][{{d-1,1}}];
        m_rest = m_rest % 0;
        for (i = 1; i < d; i+=1) {
            m[i][0] = m_rest[i-1][0];
            m[i][1] = m_rest[i-1][1];
        }
    } else {
        m = m%0;
    }
    

    wts = parameters.SetStickBreakingDistributionWeigths (m[-1][1]);


    for (i = (skip_first != 0); i < d; i+=1) {
        (v[((^distro)["rates"])[i]])[^"terms.fit.MLE"] = m[i][0];
    }
    for (i = (skip_first != 0); i < d-1; i+=1) {
        (v[((^distro)["weights"])[i]])[^"terms.fit.MLE"] = wts[i];
    }
    return v;
    
}

//------------------------------------------------------------------------------

lfunction busted.get_multi_hit (model_fit) {
    params = selection.io.extract_global_MLE_re (model_fit, '^' + ^'terms.parameters.multiple_hit_rate');
    for (k,v; in; selection.io.extract_global_MLE_re (model_fit, '^' + ^'terms.parameters.triple_hit_rate')) {
        params + v;
    }
    return params;
}

//------------------------------------------------------------------------------

lfunction busted.report_multi_hit (model_fit, json, l1, l2, bl, mdl) {

    if (^'busted.multi_hit' != "None") {
        io.ReportProgressMessageMD(l1, l2, 'Partition-level rates for multiple-hit substitutions');
        params = busted.get_multi_hit (model_fit);
        
        fracs = busted.compute_mh_fractions (bl, mdl, model_fit[^"terms.global"]);
                
        for (mle; in; params) {
            io.ReportProgressMessageMD(l1, l2, '* ' + mle[^'terms.description'] + ' : ' + Format (mle[^'terms.fit.MLE'], 8, 4));
            io.ReportProgressMessageMD(l1, l2, '* Corresponding fraction of substitutions : ' + Format (fracs[mle[^'terms.description']]*100, 6, 3) + "%");
            json[mle[^'terms.description']] = mle[^'terms.fit.MLE']; 
            json["Fraction of subs " + mle[^'terms.description']] = fracs[mle[^'terms.description']]; 
        }
    }
}

//------------------------------------------------------------------------------

lfunction busted.mixture_site_logl (ll, p) {
    res = ll["0"];
    S = Columns (ll);
    norm      = {1,S}["Max(ll[-1][_MATRIX_ELEMENT_COLUMN_],0)"];
    scaled_ll = ll["Exp(_MATRIX_ELEMENT_VALUE_-norm[_MATRIX_ELEMENT_COLUMN_])*p[_MATRIX_ELEMENT_ROW_]"];
    
    norm      = p["1"]*scaled_ll;
    ll = scaled_ll["_MATRIX_ELEMENT_VALUE_/norm[_MATRIX_ELEMENT_COLUMN_]"];
    return ll;
}

//------------------------------------------------------------------------------


lfunction busted.mixture_site_BF (ll, p) {
    prior_odds = p["_MATRIX_ELEMENT_VALUE_/(1-_MATRIX_ELEMENT_VALUE_)"];
    ll = ll["(_MATRIX_ELEMENT_VALUE_/(1-_MATRIX_ELEMENT_VALUE_))/prior_odds[_MATRIX_ELEMENT_ROW_]"];
    return ll;
}

//------------------------------------------------------------------------------

function busted.init_grid_setup (omega_distro, error_sink, rate_count) {
   for (_index_,_name_; in; omega_distro[terms.parameters.rates]) {
        if (_index_ < rate_count - 1) { // not the last rate
            busted.initial_grid  [_name_] = {
                {
                    0.01, 0.1, 0.25, 0.75
                }
            }["_MATRIX_ELEMENT_VALUE_^(rate_count-_index_-1)"];
            
            busted.initial_grid_presets [_name_] = 0;
            
            if (error_sink && _index_ == 0) {
                 busted.initial_grid  [_name_] = {{100,500,1000,5000}};
                 busted.initial_ranges [_name_] = {
                    terms.lower_bound : 100,
                    terms.upper_bound : 1000
                };
            } else {                
                busted.initial_ranges [_name_] = {
                    terms.lower_bound : 0,
                    terms.upper_bound : 1
                };
            }
        }  else {
            busted.initial_grid  [_name_] = {
                {
                    1, 1.5, 2, 4, 10
                }
            };
            busted.initial_ranges [_name_] = {
                terms.lower_bound : 1,
                terms.upper_bound : 100,
                terms.range_transform : "Sqrt(`terms.range_transform_variable`)"
            };
            busted.initial_grid_presets [_name_] = 2;
        }
    }

 
   for (_index_, _name_; in; omega_distro[terms.parameters.weights]) {
        busted.initial_grid  [_name_] = {
            {
                0.2, 0.5, 0.7, 0.8, 0.9
            }
        };
        busted.initial_grid_presets [_name_] = 3;
        if (error_sink && _index_ == 0) {
            busted.initial_grid  [_name_] = {
                {
                    0, 0.001, 0.0025, 0.025
                }
            };        
             busted.initial_ranges [_name_] = {
                terms.lower_bound : 0,
                terms.upper_bound : Sqrt (0.005),
                terms.range_transform : "`terms.range_transform_variable`^2"
            };
        } else { 
            busted.initial_ranges [_name_] = {
                terms.lower_bound : 1e-6,
                terms.upper_bound : 0.999,
                terms.range_transform : "Sqrt(`terms.range_transform_variable`)"
            };
        }
    }

}
//------------------------------------------------------------------------------

lfunction busted.compute_mh_fractions (bl, mdl, mles) {
    busted.parameter_substitutions = {};
    
    if (utility.Has (mdl, ^"terms.category", "AssociativeList")) {
        for (k, v; in; mdl[^"terms.category"]) {
             busted.parameter_substitutions [k] = 1;
    
        }
    }

    for (k, v; in; mdl[^"terms.local"]) {
        busted.parameter_substitutions [v] = 1;
    }

    busted.mh_param_mapping = {};

    for (k, v; in; mles) {
        if (k != ^"terms.parameters.multiple_hit_rate" && k != ^"terms.parameters.triple_hit_rate") {
            busted.parameter_substitutions [v[^"terms.id"]] = v[^"terms.fit.MLE"];
        } else {
            busted.mh_param_mapping [k] = v;
        }
    }

    for (k,v; in; busted.mh_param_mapping) {
         busted.parameter_substitutions [v[^"terms.id"]] = v[^"terms.fit.MLE"];
    }


    total_bl = +Simplify (bl, busted.parameter_substitutions);
    bls = {};

    for (k,v; in; busted.mh_param_mapping) {
         busted.parameter_substitutions [v[^"terms.id"]] = 0;
         blr = +Simplify (bl, busted.parameter_substitutions);
         busted.parameter_substitutions [v[^"terms.id"]] = v[^"terms.fit.MLE"];
         bls[k] = (total_bl-blr)/total_bl;
    }
    return bls;
    
}
