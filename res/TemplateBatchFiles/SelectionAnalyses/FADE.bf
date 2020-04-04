RequireVersion ("2.4.0");


LoadFunctionLibrary     ("libv3/all-terms.bf");
LoadFunctionLibrary     ("libv3/UtilityFunctions.bf");
LoadFunctionLibrary     ("libv3/IOFunctions.bf");
LoadFunctionLibrary     ("libv3/tasks/estimators.bf");
LoadFunctionLibrary     ("libv3/tasks/alignments.bf");
LoadFunctionLibrary     ("libv3/tasks/ancestral.bf");
LoadFunctionLibrary     ("libv3/tasks/trees.bf");
LoadFunctionLibrary     ("modules/io_functions.ibf");
LoadFunctionLibrary     ("modules/selection_lib.ibf");
LoadFunctionLibrary     ("libv3/convenience/math.bf");
LoadFunctionLibrary     ("libv3/models/protein.bf");
LoadFunctionLibrary     ("libv3/models/protein/empirical.bf");
LoadFunctionLibrary     ("libv3/models/protein/REV.bf");
LoadFunctionLibrary     ("libv3/tasks/mpi.bf");
LoadFunctionLibrary     ("libv3/stats.bf");
LoadFunctionLibrary ("libv3/convenience/random.bf");



utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);
utility.SetEnvVariable ("ACCEPT_ROOTED_TREES", TRUE);

namespace terms.fade {
    grid        = "grid";
    posterior   = "posterior";
    bias        = "FADE bias";
    rate        = "FADE site rate";

    namespace json {
        site_annotations = "site annotations";
        headers = "headers";
    }

    namespace cache {
        baseline     = "baseline";
        grid         = "grid";
        conditionals = "conditionals";
        settings     = "settings";
        branches     = "branches";
        model       = "model";
        model_generator   = "generator";
        posterior    = "posterior";
        root         = "root";
        substitutions = "substitutions";
        composition = "composition";
        mcmc        = "chain samples";
     }
    namespace settings {
        grid_points = "grid points";
        branches    = "branches";
    }

    namespace methods {
        MH  = "Metropolis-Hastings";
        VB0 = "Variational-Bayes";
        CG  = "Collapsed-Gibbs";
    }

};

fade.display_orders = {
                       terms.original_name: -1
                      };


fade.prompts = {
    "model" :       TRUE,
    "branches" :    TRUE,
    "grid" :        TRUE,
    "method" :      TRUE,
    "chain" :       TRUE
};

fade.run_settings = {
    "grid size" : 20,
    "chains" : 5,
    "chain-length" : 2e6,
    "burn-in" : 1e6,
    "samples" : 100,
    "concentration" : 0.5,
    "bayes factor" : 100
};

fade.table_output_options = {terms.table_options.header : TRUE, terms.table_options.minimum_column_width: 12, terms.table_options.align : "center"};


fade.analysis_description = {terms.io.info :

                            "FADE (FUBAR Approach to Directional Evolution) is a fast method to test whether or not a subset of sites in a protein alignment
                            evolve towards a particular residue along a subset of branches at accelerated rates compared to reference model.
                            FADE uses a random effects model and latent Dirichlet allocation (LDA) - inspired approximation methods to allocate sites to rate classes.",

                           terms.io.version : "0.2",
                           terms.io.reference : "TBD",
                           terms.io.authors : "Sergei L Kosakovsky Pond",
                           terms.io.contact : "spond@temple.edu",
                           terms.io.requirements : "A protein alignment and a **rooted** phylogenetic tree (optionally annotated with {})"
                          };


io.DisplayAnalysisBanner (fade.analysis_description);

fade.json = {
    terms.json.analysis: fade.analysis_description,
    terms.json.fits: {},
    terms.json.timers: {}
};

selection.io.startTimer (fade.json [terms.json.timers], "Overall", 0);

fade.parameter.bias = "FADE.bias";
fade.parameter.rate = "FADE.rate";

namespace fade {
    LoadFunctionLibrary ("modules/shared-load-file.bf");
    LoadFunctionLibrary ("modules/grid_compute.ibf");
}


// =========== LOAD DATA AND SET UP CACHES

//*******************************************************************************************//
KeywordArgument ("alignment",      "Protein alignment to screen for directional selection");
//*******************************************************************************************//

SetDialogPrompt ("Specify a protein multiple sequence alignment file");
fade.alignment_info                  = alignments.ReadProteinDataSet ("fade.dataset", None);

//*******************************************************************************************//
KeywordArgument ("output",   "Save FADE results (JSON) to [default is alignment+.FADE.json]", fade.alignment_info[terms.data.file] + ".FADE.json");
//*******************************************************************************************//
fade.alignment_info[terms.json.json] =  io.PromptUserForFilePath ("Save FADE results (JSON) to");



/** Input attribute to JSON **/

selection.io.json_store_key_value_pair (fade.json, terms.json.input, terms.json.file, fade.alignment_info [terms.data.file]);
selection.io.json_store_key_value_pair (fade.json, terms.json.input, terms.json.sequences, fade.alignment_info [terms.data.sequences]);
selection.io.json_store_key_value_pair (fade.json, terms.json.input, terms.json.sites, fade.alignment_info [terms.data.sites]);

fade.path.base = (fade.json [terms.json.input])[terms.json.file];

//*******************************************************************************************//
KeywordArgument ("cache",   "Save FADE cache to [default is alignment+.FADE.cache]", fade.path.base + ".FADE.cache");
//*******************************************************************************************//
fade.path.cache = io.PromptUserForString ("Save FADE cache to");
io.ReportProgressBar  ("init", "Loading existing cache files");
fade.cache = io.LoadCacheFromFile (fade.path.cache);
io.ClearProgressBar  ();

console.log ( "> FADE will write cache and result files to _`fade.path.cache`_ and _`fade.alignment_info[terms.json.json]`_, respectively \n\n");



fade.alignment_sample_size = fade.alignment_info[terms.data.sequences] * fade.alignment_info[terms.data.sites];

alignments.EnsureMapping ("fade.dataset", fade.alignment_info);

//*******************************************************************************************//
KeywordArgument ("tree",      "A rooted phylogenetic tree", null, "Please select a tree file for the data:");
//*******************************************************************************************//

fade.partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions (fade.alignment_info[utility.getGlobalValue("terms.data.partitions")],
                                                                              fade.alignment_info[utility.getGlobalValue("terms.data.name_mapping")]
                                                                             );
    
                                                                             
fade.partition_count = Abs (fade.partitions_and_trees);

selection.io.json_store_key_value_pair (fade.json, terms.json.input, terms.json.partition_count,fade.partition_count);

fade.filter_specification = alignments.DefineFiltersForPartitions (fade.partitions_and_trees, "fade.dataset" , "fade.filter.", fade.alignment_info);

io.ReportProgressMessageMD ("FADE", "Data", "Loaded **" +
                            fade.alignment_info [terms.data.sequences] + "** sequences, **" +
                            fade.alignment_info [terms.data.sites] + "** sites, and **" + fade.partition_count + "** partitions from \`" + fade.alignment_info [terms.data.file] + "\`");

if (utility.Has (fade.cache, terms.fade.cache.root, "AssociativeList")) {
    fade.roots = fade.cache [terms.fade.cache.root];
} else {
    fade.roots = {};
}

if (utility.Has (fade.roots,index,'String')) {
 (fade.partitions_and_trees[index])[terms.data.tree] = (trees.RootTree (_partition_[terms.data.tree], fade.roots[index]))[terms.data.tree];
}
assert (((fade.partitions_and_trees[index])[terms.data.tree])[terms.trees.rooted], "Input tree MUST be rooted");


fade.name_mapping = fade.alignment_info[utility.getGlobalValue("terms.data.name_mapping")];



if (utility.Has (fade.cache, terms.fade.cache.settings, "AssociativeList")) {
    fade.run_settings = fade.cache [terms.fade.cache.settings];
    fade.prompts["grid"]   = FALSE;
    fade.prompts["method"] = FALSE;
    fade.prompts["chain"]  = FALSE;
} else {
    fade.cache = {};
    fade.RunPrompts (fade.prompts);
    fade.cache [terms.fade.cache.settings] = fade.run_settings;
    io.WriteCacheToFile (fade.path.cache, fade.cache);
}

fade.cache [terms.fade.cache.root] = fade.roots;

if (fade.prompted_for_roots) {
    io.WriteCacheToFile (fade.path.cache, fade.cache);
}

// ===========  DEFINE TEST BRANCH SETS AND STORE TREE INFO   ====================== //

if (utility.Has (fade.cache, terms.fade.cache.branches, "AssociativeList")) {
    fade.selected_branches = fade.cache[terms.fade.cache.branches];
    fade.prompts ["branches"] = FALSE;
} else {
    fade.RunPrompts (fade.prompts);
    io.WriteCacheToFile (fade.path.cache, fade.cache);
}

fade.store_tree_information();

if (utility.Has (fade.cache, terms.fade.cache.model, "String") && utility.Has (fade.cache, terms.fade.cache.model_generator, "String")) {
    fade.baseline_model = fade.cache[terms.fade.cache.model];
    fade.generator = fade.cache[terms.fade.cache.model_generator];
    fade.prompts ["model"] = FALSE;
} else {
    fade.RunPrompts (fade.prompts);
    io.WriteCacheToFile (fade.path.cache, fade.cache);
}

fade.rebuild_lf = FALSE;

if (utility.Has (fade.cache, terms.fade.cache.baseline, "AssociativeList")) {
    io.ReportProgressMessageMD  ("FADE", "baseline", "Loaded baseline model fit from cache");
    fade.baseline_fit = fade.cache [terms.fade.cache.baseline];
    fade.rebuild_lf = TRUE;

} else {
    selection.io.startTimer (fade.json [terms.json.timers], "Baseline Fit", 1);
    io.ReportProgressMessageMD  ("FADE", "baseline", "Fitting the baseline (`fade.baseline_model`) model to obtain branch lengths and rate matrix estimates");

    lfunction fade.generator.MLE (type) {
        model = Call (^"fade.generator", type);
        model [utility.getGlobalValue("terms.model.frequency_estimator")] = "frequencies.mle";
        return model;
    }
    

    fade.baseline_fit = estimators.FitSingleModel_Ext (
                                                          fade.filter_names,
                                                          fade.trees,
                                                          "fade.generator.MLE" ,
                                                          parameters.helper.tree_lengths_to_initial_values (fade.trees, None),
                                                          {terms.run_options.retain_lf_object: TRUE}
                                                   );
    fade.cache [terms.fade.cache.baseline] =  fade.baseline_fit;
    io.WriteCacheToFile (fade.path.cache, fade.cache);
    selection.io.stopTimer (fade.json [terms.json.timers], "Baseline Fit");

}

if (utility.Has (fade.run_settings, "method", "String")) {
    fade.prompts["method"] = FALSE;
} else {
    fade.RunPrompts(fade.prompts);
    fade.cache [terms.fade.cache.settings] = fade.run_settings;
}


// ===========  FIT BASELINE MODEL


io.ReportProgressMessageMD ("FADE", "baseline", ">Fitted an alignment-wide model. " + selection.io.report_fit (fade.baseline_fit, 0, fade.alignment_sample_size ) +  "\n\nTotal tree lengths by partition\n");
utility.ForEachPair (fade.baseline_fit[terms.branch_length], "_part_", "_value_",
'
    io.ReportProgressMessageMD ("FADE", "baseline", "Partition " + (1+_part_) + ". " + Format (+(utility.Map (_value_, "_data_",
    "
        _data_ [terms.fit.MLE]
    "))
    ,6,3) + " subs/site."
    )
'
);

//Normalize EFV for saving to JSON
fade.baseline_efv = {20, 1};
for (i = 0; i < 20; i+=1)
{
    fade.efv_search = terms.characterFrequency(models.protein.alphabet[i]);  
    fade.baseline_efv[i] = ((fade.baseline_fit[terms.global])[fade.efv_search])[terms.fit.MLE];

}
fade.baseline_efv = fade.baseline_efv * (1/(+fade.baseline_efv));




selection.io.json_store_lf_withEFV(fade.json, fade.baseline_model,fade.baseline_fit[terms.fit.log_likelihood],
                            fade.baseline_fit[terms.parameters],
                            fade.alignment_sample_size, 
                            None, 
                            fade.baseline_efv,
                            0);

utility.ForEachPair (fade.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(fade.json, fade.baseline_model, terms.branch_length, 0,
                                     _key_,
                                     selection.io.extract_branch_info((fade.baseline_fit[terms.branch_length])[_key_], "selection.io.branch.length"));');




// ===========  PERFORM ANCESTRAL RECONSTRUCTION AND COUNT CHARACTERS
// TBD

if (utility.Has (fade.cache, terms.fade.cache.substitutions, "AssociativeList") == FALSE) {
    if (fade.rebuild_lf) {
        estimators.CreateLFObject (                   "fade.ancestral.rebuild",
                                                          fade.filter_names,
                                                          fade.trees,
                                                          "fade.generator.MLE" ,
                                                           fade.baseline_fit,
                                                          {terms.run_options.retain_lf_object: TRUE},
                                                          None
                                                   );

        fade.baseline_fit[terms.likelihood_function] = "fade.ancestral.rebuild.likelihoodFunction";
    }

    utility.EnsureKey (fade.cache, terms.fade.cache.substitutions);
    utility.EnsureKey (fade.cache, terms.fade.cache.composition);

    for (fade.i = 0; fade.i < fade.partition_count; fade.i += 1) {
        fade.ancestral_cache = ancestral.build (fade.baseline_fit[terms.likelihood_function], fade.i, None);
        fade.branch_filter = utility.Filter (fade.selected_branches[fade.i], "_class_", "_class_ == terms.tree_attributes.test");
        fade.partition_sites = utility.Array1D ((fade.filter_specification[fade.i])[terms.data.coverage]);

        fade.partition_substitutions = {};
        fade.partition_composition = {};

        for (fade.site = 0; fade.site < fade.partition_sites; fade.site += 1) {
            fade.partition_substitutions [fade.site] = (ancestral.ComputeSubstitutionBySite (fade.ancestral_cache, fade.site, fade.branch_filter))[terms.substitutions];
            fade.partition_composition [fade.site] = (ancestral.ComputeSiteComposition (fade.ancestral_cache, fade.site, fade.branch_filter))[terms.data.composition];
        }

        DeleteObject (fade.ancestral_cache);
        (fade.cache [terms.fade.cache.substitutions])[fade.i] = fade.partition_substitutions;
        (fade.cache [terms.fade.cache.composition])[fade.i] = fade.partition_composition;
    }
    io.WriteCacheToFile (fade.path.cache, fade.cache);
}

lfunction fade.rate.modifier (fromChar, toChar, namespace, model_type, model) {
    baseline = Call (^"fade.baseline_model.rate", fromChar,toChar, namespace, model_type, model);
    utility.EnsureKey (baseline, model_type);
    selection.io.json_store_key_value_pair (baseline, model_type, utility.getGlobalValue("terms.fade.bias"), utility.getGlobalValue("fade.parameter.bias"));
    selection.io.json_store_key_value_pair (baseline, model_type, utility.getGlobalValue("terms.fade.rate"), utility.getGlobalValue("fade.parameter.rate"));
    baseline [utility.getGlobalValue("terms.model.rate_entry")] = parameters.AppendMultiplicativeTerm (baseline [utility.getGlobalValue("terms.model.rate_entry")], utility.getGlobalValue("fade.parameter.rate"));
    if (toChar == model["fade.residue_bias"]) {
        baseline [utility.getGlobalValue("terms.model.rate_entry")] =
            parameters.AppendMultiplicativeTerm ( baseline [utility.getGlobalValue("terms.model.rate_entry")],
                                                  "`utility.getGlobalValue("fade.parameter.bias")`/(1-Exp (-`utility.getGlobalValue("fade.parameter.bias")`))");
     } else {
        if (fromChar == model["fade.residue_bias"]) {
            baseline [utility.getGlobalValue("terms.model.rate_entry")] = 
                parameters.AppendMultiplicativeTerm ( baseline [utility.getGlobalValue("terms.model.rate_entry")],
                                                  "`utility.getGlobalValue("fade.parameter.bias")`/(Exp (`utility.getGlobalValue("fade.parameter.bias")`)-1)");
        }
    }
    return baseline;
}

lfunction fade.biased.model.generator (type, residue) {
    model = Call (^"fade.generator", type);
    utility.setGlobalValue("fade.baseline_model.rate", model[utility.getGlobalValue ("terms.model.q_ij")]);
    model[utility.getGlobalValue ("terms.model.q_ij")] = "fade.rate.modifier";
    model["fade.residue_bias"] = residue;
    return model;
}

//============================================================================================================
// conditional likelihood calculations
//============================================================================================================

fade.alphabet       = "ACDEFGHIKLMNPQRSTVWY";

if (utility.Has (fade.cache, terms.fade.cache.grid, "Matrix") == FALSE) {
    fade.cache[terms.fade.cache.grid]    = fade.DefineGrid (fade.run_settings["grid size"]);
    io.WriteCacheToFile (fade.path.cache, fade.cache);
}


fade.sites_found_summary = {};

if (utility.Has (fade.cache, terms.fade.cache.conditionals, "AssociativeList") == FALSE) {
    fade.cache [terms.fade.cache.conditionals] = {};
}

if (utility.Has (fade.cache, terms.fade.cache.posterior, "AssociativeList") == FALSE) {
    fade.cache [terms.fade.cache.posterior] = {};
}

utility.EnsureKey (fade.cache, terms.fade.cache.mcmc);

fade.site_results = {};
fade.site_annotations = {};
//fade.report.posteriors = {};

namespace fade {

    site.composition.string := fade.CompositionString (((cache [^"terms.fade.cache.composition"])[partition_index])[s]);
    site.substitution.string := fade.SubstitutionHistory (((cache [^"terms.fade.cache.substitutions"])[partition_index])[s]);

    site_annotation_headers = {
                                    "Composition" : "Amino acid composition of site",
                                    "Substitutions" : "Substitution history on selected branches"
                                  };

    if (run_settings["method"] == ^"terms.fade.methods.MH") {
        table_headers = {{"rate", "Mean posterior relative rate at a site"}
                         {"bias", "Mean posterior bias parameter at a site"}
                         {"Prob[bias>0]", "Posterior probability of substitution bias towards `bias.residue`"}
                         {"BayesFactor[bias>0]", "Empiricial Bayes Factor for substitution bias towards `bias.residue`"}
                         {"PSRF", "Potential scale reduction factor - an MCMC mixing measure"}
                         {"Neff", "Estimated effective sample site for Prob [bias>0]"}};

        table_screen_output  = {{"Site", "Partition", "target", "rate", "bias", "N.eff", "Bayes Factor",site_annotation_headers["Composition"], site_annotation_headers["Substitutions"]}};
        report.biased_site = {{"" + (1+filter_info[s]),
                                        partition_index + 1,
                                        bias.residue,
                                        Format(partition_results[s][0],8,2),
                                        Format(partition_results[s][1],8,2),
                                        Format(partition_results[s][5],8,2),
                                        Format(partition_results[s][3],8,2),
                                        site.composition.string,
                                        site.substitution.string}};
    } else {
        table_headers = {{"rate", "Mean posterior relative rate at a site"}
                         {"bias", "Mean posterior bias parameter at a site"}
                         {"Prob[bias>0]", "Posterior probability of substitution bias"}
                         {"BayesFactor[bias>0]", "Empiricial Bayes Factor for substitution bias"}
                         };



         table_screen_output  = {{"Site", "Partition", "target", "rate", "bias", "Bayes Factor", site_annotation_headers["Composition"], site_annotation_headers["Substitutions"]}};
         report.biased_site = {{"" + (1+filter_info[s]),
                                        partition_index + 1,
                                        bias.residue,
                                        Format(partition_results[s][0],8,2),
                                        Format(partition_results[s][1],8,2),
                                        Format(partition_results[s][3],8,2),
                                        site.composition.string,
                                        site.substitution.string
                                        }};
    }

    for (partition_index = 0; partition_index < partition_count; partition_index += 1) {
        filter_info = (filter_specification [partition_index])[utility.getGlobalValue ("terms.data.coverage")];
        sites_in_partition = utility.Array1D (filter_info);
        site_annotations_p = {sites_in_partition, 2};
        site_annotations_p [0] = "";

        for (s = 0; s < sites_in_partition; s += 1) {
            site_annotations_p [s][0] = site.composition.string;
            site_annotations_p [s][1] = site.substitution.string;
        }

        site_annotations [partition_index] = site_annotations_p;

    }
    s = 0; partition_index = 0;
}




for (fade.residue = 0; fade.residue < 20; fade.residue += 1) {


    fade.bias.residue = fade.alphabet[fade.residue];
    selection.io.startTimer (fade.json [terms.json.timers], "Residue `fade.bias.residue` analysis", 2 + fade.residue);

    if (utility.Has (fade.cache [terms.fade.cache.conditionals], fade.bias.residue, "AssociativeList")) {
        fade.conditionals =  (fade.cache [terms.fade.cache.conditionals])[fade.bias.residue];
        io.ReportProgressBar    ("fade", "[`fade.bias.residue`] Loaded the phylogenetic likelihood function on the grid");
    } else {
        io.ReportProgressBar    ("fade", "[`fade.bias.residue`] Computing the phylogenetic likelihood function on the grid");
        fade.model.baseline = model.generic.DefineModel(fade.generator,
            "fade.baseline_model", {
                "0": "terms.global"
             },
             fade.filter_names,
             None);

        fade.model.biased = model.generic.DefineModel("fade.biased.model.generator",
            "fade.biased_model", {
                "0": "terms.global",
                "1": parameters.Quote (fade.bias.residue)
            },
            fade.filter_names,
            None);


        fade.parameter.scalers = {
            terms.fade.bias : fade.parameter.bias,
            terms.fade.rate : fade.parameter.rate
        };

        utility.Extend ((fade.biased [terms.parameters])[terms.global], fade.parameter.scalers);
        parameters.DeclareGlobalWithRanges (fade.parameter.rate, 1, 0, 100);
        parameters.DeclareGlobalWithRanges (fade.parameter.bias, 1e-10, 1e-10, 100);

        fade.model_id_to_object = {
            "fade.biased_model": fade.model.biased,
            "fade.baseline_model": fade.model.baseline
        };

        fade.trees.names = utility.Map (utility.Range (fade.partition_count, 1, 1), "_index_", "'fade.grid_tree_' + _index_");
        fade.lf.components = {fade.partition_count * 2, 1};
        
        
        utility.ForEachPair (fade.filter_names, "_index_", "_filter_",
            '
                fade.model_assignment = {
                    "fade.baseline_model" : utility.Filter (fade.selected_branches[_index_], "_value_", "_value_ == terms.tree_attributes.background"),
                    "fade.biased_model" : utility.Filter (fade.selected_branches[_index_], "_value_", "_value_ == terms.tree_attributes.test"),
                };


                fade.lf.components [2*(0+_index_)] = _filter_;
                fade.lf.components [2*(0+_index_) + 1] = fade.trees.names[_index_];
                model.ApplyModelToTree(fade.trees.names [_index_], fade.trees[_index_], None, fade.model_assignment);
            '
        );



        LikelihoodFunction fade.lf = (fade.lf.components);
        estimators.ApplyExistingEstimates  ("fade.lf", fade.model_id_to_object, fade.baseline_fit, None);
    
        fade.conditionals.raw = fade.ComputeOnGrid  ("fade.lf",
                             fade.grid.MatrixToDict (fade.cache[terms.fade.cache.grid]),
                            "fade.pass2.evaluator",
                            "fade.pass1.result_handler");
  

        (fade.cache [terms.fade.cache.conditionals])[fade.bias.residue] = fade.ConvertToConditionals (fade.conditionals.raw);
        io.WriteCacheToFile (fade.path.cache, fade.cache);
   }


    if (fade.run_settings["method"] == ^"terms.fade.methods.VB0") {
        if (utility.Has (fade.cache [terms.fade.cache.posterior], fade.bias.residue, "Matrix")) {
            io.ReportProgressBar    ("fade", "[`fade.bias.residue`] Loaded posterior means for grid loadings");
        } else {
            io.ReportProgressBar    ("fade", "[`fade.bias.residue`] Estimating posterior means for grid loadings ");
            (fade.cache[terms.fade.cache.posterior])[fade.bias.residue] = fade.RunVariationalBayes  (fade.run_settings,
                                                              fade.cache[terms.fade.cache.grid],
                                                              (fade.cache [terms.fade.cache.conditionals])[fade.bias.residue],
                                                               None
                                                              );
        }
    } else {
        if (fade.run_settings["method"] == terms.fade.methods.MH) {
          if (utility.Has (fade.cache [terms.fade.cache.mcmc], fade.bias.residue, "AssociativeList")) {
                io.ReportProgressBar    ("fade", "[`fade.bias.residue`] Loaded posterior sample for grid loadings");
          } else {
                (fade.cache[terms.fade.cache.mcmc])[fade.bias.residue] = fade.RunMCMC  (fade.run_settings,
                                                                  fade.cache[terms.fade.cache.grid],
                                                                  (fade.cache [terms.fade.cache.conditionals])[fade.bias.residue],
                                                                  "fade.pass1.result_handler",
                                                                  None);
                 io.WriteCacheToFile (fade.path.cache, fade.cache);
          }
        } else {
            if (utility.Has (fade.cache [terms.fade.cache.mcmc], fade.bias.residue, "AssociativeList")) {
                io.ReportProgressBar    ("fade", "[`fade.bias.residue`] Loaded posterior sample for grid loadings");
            } else {
                (fade.cache[terms.fade.cache.mcmc])[fade.bias.residue] = fade.RunCollapsedGibbs  (fade.run_settings,
                                                                        fade.cache[terms.fade.cache.grid],
                                                                        (fade.cache [terms.fade.cache.conditionals])[fade.bias.residue],
                                                                         None
                                                                      );
                io.WriteCacheToFile (fade.path.cache, fade.cache);
            }
        }

    }

    io.ClearProgressBar ();

    namespace fade {
        sites   = (json [utility.getGlobalValue ("terms.json.input")])[utility.getGlobalValue ("terms.json.sites")];
        grid_points = Rows (cache['grid']);
        bias_present_stencil = {grid_points,sites} ["(cache['grid'])[_MATRIX_ELEMENT_ROW_][1]>0."];

        rates  = Transpose ((cache['grid'])[-1][0]);
        biases = Transpose ((cache['grid'])[-1][1]);
                
 

        if (run_settings["method"] != ^"terms.fade.methods.VB0") {

            samples = run_settings["samples"];
            chains  = run_settings["chains"];

            results.log_L   = {1,samples};
            results.samples = {samples,grid_points};

            per_chain      = samples $ chains;


            from = 0;
            to   = per_chain;


            posterior_mean_rates                 = {};
            posterior_mean_biases                = {};
            denominators                         = {};
            posteriors                           = {};
            biased_ks                            = {};


            for (chain_id = 0; chain_id < chains; chain_id += 1) {
                io.ReportProgressBar                  ("PROCESSING", "Samples from chain " + (chain_id + 1));


                grid_samples = (((cache[utility.getGlobalValue("terms.fade.cache.mcmc")])[bias.residue])[chain_id])["weights"];
                grid_samples_T = Transpose (grid_samples);
                P_ks = grid_samples *
                       ((cache[utility.getGlobalValue("terms.fade.cache.conditionals")])[bias.residue])["conditionals"];

                denominators[chain_id] = P_ks;


                posterior_mean_rates[chain_id]      =            (grid_samples $ rates *
                                                                 ((cache[utility.getGlobalValue("terms.fade.cache.conditionals")])[bias.residue])["conditionals"]) / P_ks;

                posterior_mean_biases[chain_id]     =            (grid_samples $ biases *
                                                                 ((cache[utility.getGlobalValue("terms.fade.cache.conditionals")])[bias.residue])["conditionals"]) / P_ks;


                biased_ks[chain_id]  = grid_samples *
                               (bias_present_stencil $ ((cache[utility.getGlobalValue("terms.fade.cache.conditionals")])[bias.residue])["conditionals"]) / P_ks;

                if (run_settings["method"] == ^"terms.fade.methods.MH") {
                    logL_samples = (((cache[utility.getGlobalValue("terms.fade.cache.mcmc")])[bias.residue])[chain_id])["likelihoods"];
                    draw_from_this_chain                =         Random ({1,samples}["_MATRIX_ELEMENT_COLUMN_"], 0);

                    for (i = from; i < to; i += 1) {
                        draw_this_index = draw_from_this_chain[i];
                        results.log_L [i] = logL_samples [draw_this_index];
                        for (r = 0; r < grid_points; r += 1) {
                            results.samples [i][r] = grid_samples[draw_this_index][r];
                        }
                    }
                } else {
                    results.samples = grid_samples;
                }

                from = to;
                if (chain_id == chains - 2) {
                    to = samples;
                } else {
                    to += per_chain;
                }
            }
            io.ClearProgressBar                   ();
            posterior_mean_over_grid                = {grid_points,1}["(+results.samples[-1][_MATRIX_ELEMENT_ROW_])/samples"];

        } else {
             posterior_mean_over_grid                 = (cache[^"terms.fade.cache.posterior"])[bias.residue];
             posterior_mean_over_grid_T               = Transpose (posterior_mean_over_grid);
             cache[terms.fade.cache.posterior]       = posterior_mean_over_grid;
 
             P_ks = posterior_mean_over_grid_T * ((cache[utility.getGlobalValue("terms.fade.cache.conditionals")])[bias.residue])["conditionals"];


             posterior_mean_rates     =             (posterior_mean_over_grid_T $ rates *
                                                              ((cache[utility.getGlobalValue("terms.fade.cache.conditionals")])[bias.residue])["conditionals"]) / P_ks;
             posterior_mean_biases     =            (posterior_mean_over_grid_T $ biases *
                                                              ((cache[utility.getGlobalValue("terms.fade.cache.conditionals")])[bias.residue])["conditionals"]) / P_ks;

             biased_ks  = posterior_mean_over_grid_T *
                               (bias_present_stencil $ ((cache[utility.getGlobalValue("terms.fade.cache.conditionals")])[bias.residue])["conditionals"]) / P_ks;

       }

        prior_weight_bias                        = +posterior_mean_over_grid ["_MATRIX_ELEMENT_VALUE_*((cache['grid'])[_MATRIX_ELEMENT_ROW_][1]>0.)"];

        headers.printed = FALSE;

        i = 0;
        s = 0; // reset to ensure good re-entrant behavior

        report.sites_found = {};
        //report.posteriors[bias.residue] = {};
        site_results[bias.residue] = {};

        chain_iterator = utility.Range (chains, 0, 1);

        for (partition_index = 0; partition_index < partition_count; partition_index += 1) {
            filter_info = (filter_specification [partition_index])[utility.getGlobalValue ("terms.data.coverage")];
            sites_in_partition = utility.Array1D (filter_info);

            partition_posteriors = {};
            if (run_settings["method"] != utility.getGlobalValue ("terms.fade.methods.MH")) {
                partition_results    = {sites_in_partition, 4};
            } else {
                partition_results    = {sites_in_partition, 6};
           }



            for (s = 0; s < sites_in_partition; s += 1) {

                pp = posterior_mean_over_grid $ (((cache[utility.getGlobalValue("terms.fade.cache.conditionals")])[bias.residue])["conditionals"])[-1][i];
                partition_posteriors [s] = Transpose (pp * (1/(+pp)));

                if (run_settings["method"] != utility.getGlobalValue ("terms.fade.methods.VB0")) {
                    partition_results[s][0] = fade.ComputeRandNeff (
                        utility.Map (chain_iterator, "_value_", "((`&posterior_mean_rates`)[_value_])[-1][`&i`]")
                    )[0];
                    partition_results[s][1] = fade.ComputeRandNeff (
                        utility.Map (chain_iterator, "_value_", "((`&posterior_mean_biases`)[_value_])[-1][`&i`]")
                    )[0];

                    biased_posterior = fade.ComputeRandNeff (
                        utility.Map (chain_iterator, "_value_", "((`&biased_ks`)[_value_])[-1][`&i`]")
                    );

                    partition_results[s][2] = biased_posterior[0];
                    partition_results[s][3] = stats.BayesFactor (prior_weight_bias, biased_posterior[0]) ;

                    if (run_settings["method"] == utility.getGlobalValue ("terms.fade.methods.MH")) {
                        partition_results[s][4] = biased_posterior[1];
                        partition_results[s][5] = biased_posterior[2];
                    }

                } else {

                    if (run_settings["method"] == utility.getGlobalValue ("terms.fade.methods.VB0")) {
                        partition_results [s][0] = posterior_mean_rates[i];
                        partition_results [s][1] = posterior_mean_biases[i];
                        partition_results [s][2] = biased_ks[i];
                        partition_results [s][3] = stats.BayesFactor (prior_weight_bias, biased_ks[i]);
                    }
                }


                if (partition_results[s][3] >= run_settings["bayes factor"]) {
                    if (Abs(report.sites_found) == 0 && table_output_options[^"terms.table_options.header"]) {
                        fprintf (stdout, "\n", io.FormatTableRow (table_screen_output,table_output_options));
                        table_output_options[^"terms.table_options.header"] = FALSE;
                    }
                    fprintf (stdout, io.FormatTableRow (report.biased_site,table_output_options));
                    report.sites_found + (1-partition_results [s][3]);
                }

                i+=1;
            }

            (site_results[bias.residue]) [partition_index] = partition_results;
            //(report.posteriors[bias.residue]) [partition_index] = partition_posteriors;
            s = 0; // for re-entrancy
        }

        partition_index = 0;
        sites_found = Abs(report.sites_found);
        sites_found_summary [bias.residue] = Abs(report.sites_found);

     }
     selection.io.stopTimer (fade.json [terms.json.timers], "Residue `fade.bias.residue` analysis");

}



// ===========  ANALYSIS SUMMARY ========

fade.json [terms.fade.cache.settings] = fade.run_settings;
fade.json [terms.fade.json.site_annotations] = {
    terms.fade.json.headers : fade.site_annotation_headers,
    terms.fade.json.site_annotations : fade.site_annotations
};
fade.json [terms.fit.MLE] = {terms.json.headers   : fade.table_headers,
                               terms.json.content : fade.site_results };

//fade.json [terms.fade.posterior] = fade.report.posteriors;

console.log ("----\n## FADE analysis summary. Evidence for directional selection evaluated using empirical Bayes factor threshold of " + fade.run_settings["bayes factor"]);


utility.ForEachPair (fade.sites_found_summary, "_residue_", "_count_",
'
    if (_count_ == 0) {
        console.log ("* No sites are evolving directionally towards " + _residue_);
    } else {
        console.log ("* " + _count_ + " " + io.SingularOrPlural (_count_, "site is", "sites are") + " evolving directionally towards " + _residue_);
    }

');


selection.io.stopTimer (fade.json [terms.json.timers], "Overall");

io.SpoolJSON (fade.json, fade.alignment_info[terms.json.json]);

// HELPER FUNCTIONS GO HERE
//----------------------------------------------------------------------------

function     fade.RunPrompts (prompts) {
    if (prompts["branches"]) {
        //*******************************************************************************************//
        KeywordArgument ("branches", "The branches to test", "All");
        //*******************************************************************************************//
        fade.selected_branches = selection.io.defineBranchSets ( fade.partitions_and_trees );
        fade.cache [terms.fade.cache.branches] = fade.selected_branches;
        prompts["branches"] = FALSE;
    }

     if (prompts["grid"]) {
        //*******************************************************************************************//
        KeywordArgument ("grid", "The number of grid points", "20");
        //*******************************************************************************************//
        fade.run_settings["grid size"] = io.PromptUser ("> Number of grid points per dimension (total number is D^2)",fade.run_settings["grid size"],5,50,TRUE);
        prompts["grid"] = FALSE;
    }



    if (prompts["model"]) {
        utility.Extend (models.protein.empirical_models, {"GTR" : "General time reversible model (189 estimated parameters)."});
        //*******************************************************************************************//
        KeywordArgument ("model", "The substitution model to use", "GTR");
        //*******************************************************************************************//
        fade.baseline_model         = io.SelectAnOption (models.protein.empirical_models, "Baseline substitution model");
        fade.generator              = (utility.Extend (models.protein.empirical.plusF_generators , {"GTR" : "models.protein.REV.ModelDescription"}))[fade.baseline_model ];
        fade.cache[terms.fade.cache.model] = fade.baseline_model;
        fade.cache[terms.fade.cache.model_generator] = fade.generator;
        prompts["model"] = FALSE;
    }


     if (prompts["method"]) {
        //*******************************************************************************************//
        KeywordArgument ("method", "Inference method to use", "`terms.fade.methods.VB0`");
        //*******************************************************************************************//
        fade.run_settings["method"] = io.SelectAnOption  ({
                                                                terms.fade.methods.VB0 : "0-th order Variational Bayes approximations (fastest, recommended default)",
                                                                terms.fade.methods.CG : "Collapsed Gibbs sampler (intermediate speed)",
                                                                terms.fade.methods.MH : "Full Metropolis-Hastings MCMC algorithm (slowest, original 2013 paper implementation)"
                                                            }, "Posterior estimation method");
        prompts["method"] = FALSE;
     }

    if (prompts["chain"]) {
        if (fade.run_settings["method"] ==  terms.fade.methods.MH) {
            //*******************************************************************************************//
            KeywordArgument ("chains", "How many MCMC chains to run", fade.run_settings["chains"]);
            //*******************************************************************************************//
            fade.run_settings["chains"] = io.PromptUser ("> Number of MCMC chains to run",fade.run_settings["chains"],2,20,TRUE);
        } else {
            fade.run_settings["chains"] = 1;
        }
        if (fade.run_settings["method"] !=  terms.fade.methods.VB0) {
            KeywordArgument ("chain-length", "MCMC chain length", fade.run_settings["chain-length"]);
            fade.run_settings["chain-length"] = io.PromptUser ("> The length of each chain",fade.run_settings["chain-length"],5e3,5e7,TRUE);
            KeywordArgument ("burn-in", "MCMC chain burn in", fade.run_settings["chain-length"]$2);
            fade.run_settings["burn-in"] = io.PromptUser ("> Use this many samples as burn-in",fade.run_settings["chain-length"]$2,fade.run_settings["chain-length"]$20,fade.run_settings["chain-length"]*95$100,TRUE);
            KeywordArgument ("samples", "MCMC samples to draw", fade.run_settings["samples"]);
            fade.run_settings["samples"] = io.PromptUser ("> How many samples should be drawn from each chain",fade.run_settings["samples"],50,fade.run_settings["chain-length"]-fade.run_settings["burn-in"],TRUE);
        }
        KeywordArgument ("concentration_parameter", "The concentration parameter of the Dirichlet prior", fade.run_settings["concentration"]);
        fade.run_settings["concentration"] = io.PromptUser  ("> The concentration parameter of the Dirichlet prior",fade.run_settings["concentration"],0.001,1,FALSE);
        prompts["chain"] = FALSE;
    }
}

//------------------------------------------------------------------------------------------------//

lfunction fade.DefineGrid (one_d_points) {
    // only one point for rate = 0, because bias is not identifiable if rate = 0

    one_d_points    = Max (one_d_points, 5);

    alphaBetaGrid = {one_d_points^2,2}; // (alpha, beta) pair

    oneDGridRate      = {one_d_points,1};
    oneDGridBias      = {one_d_points,1};

    below1_frac         = 0.7;
    below1  = ((one_d_points)*below1_frac+0.5)$1;
    above1  = (one_d_points-1)*(1-below1_frac)$1;

    if (below1 + above1 != one_d_points) {
        above1 = one_d_points - below1;
    }

    _neg_step = 1/(below1);
    _neg_stepP1 = 1/(below1+1);

    for (_k = 0; _k < below1; _k += 1) {
        oneDGridBias [_k][0] =  _neg_step * (_k);
        oneDGridRate [_k][0] =  _neg_stepP1 * (_k + 1);
    }

    oneDGridRate [below1-1][0] = 1;
    oneDGridBias [below1-1][0] = 1;

    _pos_step = 49^(1/3)/above1;
    for (_k = 1; _k <= above1; _k += 1) {
        oneDGridBias [below1+_k-1][0] = 1+(_pos_step*_k)^3;
        oneDGridRate [below1+_k-1][0] = 1+(_pos_step*_k)^3;
    }

    _p = 0;

    for (_r = 0; _r < one_d_points; _r += 1) {
        for (_c = 0; _c < one_d_points; _c += 1) {
           alphaBetaGrid[_p][0] = oneDGridRate[_r];
           alphaBetaGrid[_p][1] = oneDGridBias[_c];
           _p += 1;
        }
    }
    alphaBetaGrid[0][0] = 0; alphaBetaGrid[0][1] = 0;
    alphaBetaGrid[1][1] = 0;

    return alphaBetaGrid;
}

//------------------------------------------------------------------------------------------------//

lfunction fade.SubstitutionHistory (subs) {
    result = "";
    result * 128;
    keys = utility.sortStrings (utility.Keys (subs));

    for (i = 0; i < Abs (subs); i+=1) {
        source  = keys[i];
        targets  = subs[source];
        if (i > 0) {
            result * ", ";
        }
        result * (source + "->");
        keys2 =  utility.sortStrings (utility.Keys (targets));
        for (k = 0; k < Abs (targets); k+=1) {
           result * (keys2[k] + "(" + Abs(targets[keys2[k]]) + ")");
        }
    }

    result * 0;
    return result;

}

//------------------------------------------------------------------------------------------------//

lfunction fade.CompositionString (composition) {
    result = "";
    result * 128;
    keys = utility.sortStrings (utility.Keys (composition));

    for (i = 0; i < Abs (composition); i+=1) {
        residue  = keys[i];
        if (i) {
            result * ",";
        }
        result * (residue + composition [residue]);
    }

    result * 0;
    return result;

}


//------------------------------------------------------------------------------------------------//


lfunction fade.grid.MatrixToDict (grid) {
    return utility.Map (utility.MatrixToListOfRows (grid), "_value_",
                                                                '{  terms.fade.bias : {
                                                                            terms.id : fade.parameter.scalers [ terms.fade.bias ],
                                                                            terms.fit.MLE : _value_[1]
                                                                        },
                                                                    terms.fade.rate :  {
                                                                            terms.id : fade.parameter.scalers [ terms.fade.rate ],
                                                                            terms.fit.MLE : _value_[0]
                                                                        }
                                                                 }');
}
