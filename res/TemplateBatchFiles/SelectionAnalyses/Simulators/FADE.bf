RequireVersion ("2.3.12");


LoadFunctionLibrary     ("libv3/all-terms.bf");
LoadFunctionLibrary     ("libv3/UtilityFunctions.bf");
LoadFunctionLibrary     ("libv3/IOFunctions.bf");
LoadFunctionLibrary     ("libv3/tasks/estimators.bf");
LoadFunctionLibrary     ("libv3/tasks/alignments.bf");
LoadFunctionLibrary     ("libv3/tasks/ancestral.bf");
LoadFunctionLibrary     ("libv3/tasks/trees.bf");
LoadFunctionLibrary     ("../modules/io_functions.ibf");
LoadFunctionLibrary     ("../modules/selection_lib.ibf");
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
    mode = "mode";
    regimes = "regimes";
    bias = "substitution bias";
    rate = "rate multiplier";
};


fade.parameter.bias = "FADE.bias";
fade.parameter.rate = "FADE.rate";
fade.tree.name      = "FADE.simulated_tree";

fade.settings = {};
fade.alphabet           = "ACDEFGHIKLMNPQRSTVWY";
fade.alphabet.matrix    = {1,20};
fade.simulation.matrix  = {2,20};
fade.evolutionary_modes = {"Null" : "Evolution under baseline model"};

for (r = 0; r < Abs (fade.alphabet); r += 1) {
    fade.evolutionary_modes [fade.alphabet[r]] = "Directional evolution towards `fade.alphabet[r]`";
    fade.alphabet.matrix [r] = fade.alphabet[r];
    fade.simulation.matrix[0][r] = fade.alphabet[r];
}

fade.simulation.matrix[1][0] = "1";

fade.analysis_description = {terms.io.info :
                            "A companion data simulator for FADE",
                            terms.io.version :          "0.1",
                            terms.io.reference :        "TBD",
                            terms.io.authors :          "Sergei L Kosakovsky Pond",
                            terms.io.contact :          "spond@temple.edu",
                            terms.io.requirements :     "A **rooted** phylogenetic tree with branch lengths (optionally annotated with {} to define a branch partition set)"
                          };

io.DisplayAnalysisBanner (fade.analysis_description);


// =========== LOAD DATA AND SET UP CACHES

SetDialogPrompt ("Specify a rooted tree to use for data simulations");
fade.baseline.tree = trees.LoadAnnotatedTopology (TRUE);
assert (trees.HasBranchLengths(fade.baseline.tree), "Input tree MUST have branch lengths");
assert (fade.baseline.tree[terms.trees.rooted], "Input tree MUST be rooted");


fade.replicates = io.PromptUser ("How many replicate datasets be simulated", 100, 1, 10000, true);
fade.sites_class_count = io.PromptUser ("How many types of sites will be simulated", 2, 1, 10000, true);
fade.site_classes      = {};

for (k = 0; k < fade.sites_class_count; k += 1) {
    this_class = {
                    terms.data.sites : io.PromptUser ("How many sites are in class " + (k+1), 100, 1, 10000, true),
                    terms.fade.mode : io.SelectAnOption (fade.evolutionary_modes, "Evolutionary regime for site class " + (k+1)),
                    fade.parameter.rate : io.PromptUser ("Relative overall rate for site class " + (k+1) + " (1 = average)", 1, 0, 1000, false)
                 };

    if (this_class [terms.fade.mode] != "Null") {
        this_class [fade.parameter.bias] = io.PromptUser ("Substitution bias for site class " + (k+1) + " (0 = no bias)", 1, 0, 1000, false);
    }
    fade.site_classes [k] = this_class;
}


fade.selected_branches = (selection.io.defineBranchSets ( {"0" : { terms.data.tree : fade.baseline.tree}} ))[0];

fade.settings [terms.data.tree]    = fade.baseline.tree[terms.trees.newick_with_lengths];
fade.settings [terms.json.tested]  = fade.selected_branches;
fade.settings [terms.fade.regimes] = fade.site_classes;
fade.settings [terms.replicates]   = fade.replicates;

utility.Extend               (models.protein.empirical_models, {"GTR" : "General time reversible model (189 estimated parameters)."});
fade.baseline_model         = io.SelectAnOption (models.protein.empirical_models, "Baseline substitution model");
fade.generator              = (utility.Extend (models.protein.empirical.plusF_generators , {"GTR" : "models.protein.REV.ModelDescription"}))[fade.baseline_model ];
fade.branch_lengths         = parameters.helper.tree_lengths_to_initial_values ({"0" : fade.baseline.tree}, None);

fade.settings [terms.model]          = fade.baseline_model;
fade.settings [terms.fade.generator] = fade.generator;


lfunction fade.rate.modifier (fromChar, toChar, namespace, model_type, model) {

    baseline = Call (^"fade.baseline_model.rate", fromChar,toChar, namespace, model_type, model);
    utility.EnsureKey (baseline, model_type);
    selection.io.json_store_key_value_pair (baseline, model_type, utility.getGlobalValue("terms.fade.bias"), utility.getGlobalValue("fade.parameter.bias"));
    selection.io.json_store_key_value_pair (baseline, model_type, utility.getGlobalValue("terms.fade.rate"), utility.getGlobalValue("fade.parameter.rate"));
    baseline [utility.getGlobalValue("terms.model.rate_entry")] = parameters.AppendMultiplicativeTerm (baseline [utility.getGlobalValue("terms.model.rate_entry")], utility.getGlobalValue("fade.parameter.rate"));
    if ( Type (model["fade.residue_bias"]) == "String") {
        if (toChar == model["fade.residue_bias"]) {
            baseline [utility.getGlobalValue("terms.model.rate_entry")] =
                parameters.AppendMultiplicativeTerm ( baseline [utility.getGlobalValue("terms.model.rate_entry")],
                                                      "`utility.getGlobalValue("fade.parameter.bias")`/(1-Exp (-`utility.getGlobalValue("fade.parameter.bias")`))");
         } else {
            if (fromChar == model["fade.residue_bias"]) {
                parameters.AppendMultiplicativeTerm ( baseline [utility.getGlobalValue("terms.model.rate_entry")],
                                                      "`utility.getGlobalValue("fade.parameter.bias")`/(Exp (`utility.getGlobalValue("fade.parameter.bias")`-1))");
            }
        }
    }
    return baseline;
}



lfunction fade.biased.model.generator (type, residue) {
    model = Call (^"fade.generator", type);
    utility.setGlobalValue("fade.baseline_model.rate", model[utility.getGlobalValue ("terms.model.q_ij")]);
    model[utility.getGlobalValue ("terms.model.q_ij")] = "fade.rate.modifier";
    model["fade.residue_bias"] = residue;
    model[utility.getGlobalValue ("terms.alphabet")] = utility.getGlobalValue ("fade.alphabet.matrix");
    return model;
}


fade.bias.residue = "F";

fade.model.biased = model.generic.DefineModel("fade.biased.model.generator",
            "fade.biased_model", {
                "0": "terms.global",
                "1": parameters.Quote (fade.bias.residue)
            },
            None,
            "frequencies.equal");


fade.model.baseline = model.generic.DefineModel("fade.biased.model.generator",
            "fade.baseline_model", {
                "0": "terms.global",
                "1": None
            },
            None,
            "frequencies.equal");

fade.model_id_to_object = {
            "fade.biased_model": fade.model.biased,
            "fade.baseline_model": fade.model.baseline
        };


fade.model_assignment = {
    "fade.baseline_model" : utility.Filter (fade.selected_branches, "_value_", "_value_ == terms.tree_attributes.background"),
    "fade.biased_model" : utility.Filter (fade.selected_branches, "_value_", "_value_ == terms.tree_attributes.test"),
};



parameters.DeclareGlobalWithRanges (fade.parameter.rate, 1, 0, 100);
parameters.DeclareGlobalWithRanges (fade.parameter.bias, 1e-10, 1e-10, 100);

model.ApplyModelToTree(fade.tree.name, fade.baseline.tree, None, fade.model_assignment);
//function estimators.ApplyExistingEstimatesToTree (_tree_name, model_descriptions, initial_values, _application_type, keep_track_of_proportional_scalers) {


parameters.SetValue (fade.parameter.bias, 1e-10);
parameters.SetValue (fade.parameter.rate, 1);

estimators.ApplyExistingEstimatesToTree (fade.tree.name, fade.model_id_to_object, (fade.branch_lengths[terms.branch_length])[0], None, {});

fprintf (stdout, Format (^fade.tree.name, 1, 1), "\n");

parameters.SetValue (fade.parameter.bias, 10);
parameters.SetValue (fade.parameter.rate, 5);

fprintf (stdout, Format (^fade.tree.name, 1, 1), "\n");


fade.sim_frequencies = fade.model.biased[terms.efv_estimate];

DataSet simulated_block = Simulate (^fade.tree.name, fade.sim_frequencies , fade.simulation.matrix, 200);
DataSetFilter simulated_block_filter = CreateFilter (simulated_block, 1);

fprintf (stdout, simulated_block_filter, "\n");

return 0;



namespace fade {

    site.composition.string := fade.CompositionString (((cache [^"terms.fade.cache.composition"])[partition_index])[s]);
    site.substitution.string := fade.SubstitutionHistory (((cache [^"terms.fade.cache.substitutions"])[partition_index])[s]);

    site_annotation_headers = {
                                    "Composition" : "Aminoacid composition of site",
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
        fade.selected_branches = selection.io.defineBranchSets ( fade.partitions_and_trees );
        fade.cache [terms.fade.cache.branches] = fade.selected_branches;
        prompts["branches"] = FALSE;
    }

     if (prompts["grid"]) {
        fade.run_settings["grid size"] = io.PromptUser ("> Number of grid points per dimension (total number is D^2)",fade.run_settings["grid size"],5,50,TRUE);
        prompts["grid"] = FALSE;
    }



    if (prompts["model"]) {
        utility.Extend (models.protein.empirical_models, {"GTR" : "General time reversible model (189 estimated parameters)."});
        fade.baseline_model         = io.SelectAnOption (models.protein.empirical_models, "Baseline substitution model");
        fade.generator              = (utility.Extend (models.protein.empirical.plusF_generators , {"GTR" : "models.protein.REV.ModelDescription"}))[fade.baseline_model ];
        fade.cache[terms.fade.cache.model] = fade.baseline_model;
        fade.cache[terms.fade.cache.model_generator] = fade.generator;
        prompts["model"] = FALSE;
    }


     if (prompts["method"]) {
        fade.run_settings["method"] = io.SelectAnOption  ({
                                                                terms.fade.methods.MH : "Full Metropolis-Hastings MCMC algorithm (slowest, original 2013 paper implementation)",
                                                                terms.fade.methods.CG : "Collapsed Gibbs sampler (intermediate speed)",
                                                                terms.fade.methods.VB0 : "0-th order Variational Bayes approximations (fastest, recommended default)"
                                                            }, "Posterior estimation method");
        prompts["method"] = FALSE;
     }

    if (prompts["chain"]) {
        if (fade.run_settings["method"] ==  terms.fade.methods.MH) {
            fade.run_settings["chains"] = io.PromptUser ("> Number of MCMC chains to run",fade.run_settings["chains"],2,20,TRUE);
        } else {
            fade.run_settings["chains"] = 1;
        }
        if (fade.run_settings["method"] !=  terms.fade.methods.VB0) {
            fade.run_settings["chain-length"] = io.PromptUser ("> The length of each chain",fade.run_settings["chain-length"],5e3,5e7,TRUE);
            fade.run_settings["burn-in"] = io.PromptUser ("> Use this many samples as burn-in",fade.run_settings["chain-length"]$2,fade.run_settings["chain-length"]$20,fade.run_settings["chain-length"]*95$100,TRUE);
            fade.run_settings["samples"] = io.PromptUser ("> How many samples should be drawn from each chain",fade.run_settings["samples"],50,fade.run_settings["chain-length"]-fade.run_settings["burn-in"],TRUE);
        }
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
