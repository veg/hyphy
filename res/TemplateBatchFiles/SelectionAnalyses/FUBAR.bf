RequireVersion  ("2.4.0");

LoadFunctionLibrary("libv3/all-terms.bf");
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary ("libv3/models/codon/MG_REV.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/tasks/mpi.bf");

LoadFunctionLibrary("modules/io_functions.ibf");
LoadFunctionLibrary("modules/selection_lib.ibf");
LoadFunctionLibrary("libv3/models/codon/BS_REL.bf");
LoadFunctionLibrary("libv3/convenience/math.bf");
LoadFunctionLibrary ("libv3/convenience/random.bf");


utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);

/*------------------------------------------------------------------------------*/

namespace terms.fubar {
    grid        = "grid";
    posterior   = "posterior";

    namespace cache {
        gtr          = "gtr-fit";
        grid         = "grid";
        conditionals = "conditionals";
        mcmc         = "mcmc";
        posterior    = "posterior";
        settings     = "settings";
     }
    namespace settings {
        grid_points = "grid points";
    }
    namespace branch_length_scaler {
        beta = "fubar beta scaler";
        alpha = "fubar alpha scaler";
    }

    namespace methods {
        MH  = "Metropolis-Hastings";
        VB0 = "Variational-Bayes";
        CG  = "Collapsed-Gibbs";
    }
};

fubar.json = {
                    terms.json.analysis:    {},
                    terms.json.input:       {},
                    terms.json.fits :       {},
                    terms.json.timers :     {},
                    terms.fit.MLE :         {},
                    terms.settings :        {},
                    terms.fubar.grid :      {},
                    terms.fubar.posterior  :      {}
                };

fubar.display_orders = {terms.original_name: -1,
                        terms.json.nucleotide_gtr: 0
                        };


fubar.prompts = {
    "grid" : TRUE,
    "chain" : TRUE,
    "method" : TRUE,
    "non-zero" : TRUE
};

fubar.run_settings = {
    "grid size" : 20,
    "chains" : 5,
    "chain-length" : 2e6,
    "burn-in" : 1e6,
    "samples" : 100,
    "concentration" : 0.5,
    "posterior" : 0.9,
    "non-zero" : FALSE
};



/*------------------------------------------------------------------------------*/

fubar.analysis_description = {terms.io.info :

    "Perform a Fast Unbiased AppRoximate Bayesian (FUBAR)
    analysis of a coding sequence alignment to determine whether some
    sites have been subject to pervasive purifying or diversifying
    selection.

    v2.1 introduces two more methods for estimating the posterior distribution of grid weights:
    collapsed Gibbs MCMC (faster) and 0-th order Variation Bayes approximation (fastest).

    Please note that a FUBAR analysis generates a cache and a results JSON file in the same directory as
    directory as the original alignment. HyPhy needs to have write
    privileges to this directory. For example if the original file is in
    /home/sergei/FUBAR/data/pol.nex then at the end of a FUBAR run, there
    will also exist FUBAR-generated files /home/sergei/FUBAR/data/pol.nex.FUBAR.json,
    /home/sergei/FUBAR/data/pol.nex.fubrar.cache. They also
    provide checkpointing so that a partially completed analysis can be
    restarted.",

           terms.io.version : "2.2",
           terms.io.reference : "FUBAR: a fast, unconstrained bayesian approximation for inferring selection (2013), Mol Biol Evol. 30(5):1196-205",
           terms.io.authors : "Sergei L Kosakovsky Pond",
           terms.io.contact : "spond@temple.edu",
           terms.io.requirements : "in-frame codon alignment (possibly partitioned) and a phylogenetic tree (one per partition)"
          };

io.DisplayAnalysisBanner ( fubar.analysis_description );

fubar.json[terms.json.analysis] = fubar.analysis_description;


/*------------------------------------------------------------------------------
    Key word arguments
*/
KeywordArgument ("code", "Which genetic code should be used", "Universal");
KeywordArgument ("alignment", "An in-frame codon alignment in one of the formats supported by HyPhy");
KeywordArgument ("tree", "A phylogenetic tree (optionally annotated with {})", null, "Please select a tree file for the data:");

// Additional Keyword Arguments ("output" and "cache") called below after namespace fubar.
// Additional Keyword Arguments ("non-zero" "grid", "method", "chain", "chain-length", "burn-in", "samples" and "concentration") called from within the "fubar.RunPrompts" function.


namespace fubar {
    LoadFunctionLibrary ("modules/shared-load-file.bf");
    LoadFunctionLibrary ("modules/grid_compute.ibf");
    load_file ({utility.getGlobalValue("terms.prefix"): "fubar", utility.getGlobalValue("terms.settings") : {utility.getGlobalValue("terms.settings.branch_selector") : "selection.io.SelectAllBranches"}});

}

selection.io.startTimer (fubar.json [terms.json.timers], "Overall", 0);

fubar.path.base = (fubar.json [terms.json.input])[terms.json.file];

KeywordArgument ("cache",   "Save FUBAR cache to [default is alignment+.FUBAR.cache]", fubar.path.base + ".FUBAR.cache");
KeywordArgument ("output",   "Save FUBAR results (JSON) to [default is alignment+.FUBAR.json]", fubar.codon_data_info[terms.data.file] + ".FUBAR.json");


/*------------------------------------------------------------------------------
    Continued Analysis Setup
*/

fubar.path.cache = io.PromptUserForString ("Save FUBAR cache to");
fubar.cache = io.LoadCacheFromFile (fubar.path.cache);

fubar.table_output_options = {terms.table_options.header : TRUE, terms.table_options.minimum_column_width: 16, terms.table_options.align : "center"};

console.log ( "> FUBAR will write cache and result files to _`fubar.path.base`.FUBAR.cache_ and _`fubar.path.base`.FUBAR.json_, respectively \n\n");

//----------------------------------------------------------------------------
// PHASE 1: nucleotide fit
//----------------------------------------------------------------------------


if (utility.Has (fubar.cache, terms.fubar.cache.settings, "AssociativeList")) {
    fubar.run_settings = fubar.cache [terms.fubar.cache.settings];
} else {
    fubar.cache = {};
    fubar.RunPrompts (fubar.prompts);
    fubar.cache [terms.fubar.cache.settings] = fubar.run_settings;
    io.WriteCacheToFile (fubar.path.cache, fubar.cache);
}


selection.io.startTimer (fubar.json [terms.json.timers], "Nucleotide Fit", 1);
if (utility.Has (fubar.cache, terms.fubar.cache.gtr, "AssociativeList")) {
     fubar.gtr_results = fubar.cache[terms.fubar.cache.gtr];
     io.ReportProgressMessageMD("fubar", "nuc-fit", "_Loaded cached GTR fit_ " + selection.io.report_fit (fubar.gtr_results, 0, fubar.codon_data_info[terms.data.sample_size]));

    /* Add nucleotide fit to the JSON */
    gtr_rates = utility.Map(
                    utility.Map (fubar.gtr_results[utility.getGlobalValue("terms.global")], "_value_", '   {terms.fit.MLE : _value_[terms.fit.MLE]}'),
                    "_value_",
                    "_value_[terms.fit.MLE]");
    efv = (fubar.gtr_results[utility.getGlobalValue("terms.efv_estimate")])["VALUEINDEXORDER"][0];
    selection.io.json_store_lf_withEFV (fubar.json,
                                utility.getGlobalValue ("terms.json.nucleotide_gtr"),
                                fubar.gtr_results[utility.getGlobalValue ("terms.fit.log_likelihood")],
                                fubar.gtr_results[utility.getGlobalValue ("terms.parameters")] ,
                                fubar.codon_data_info[utility.getGlobalValue ("terms.data.sample_size")],
                                gtr_rates,
                                efv,
                                fubar.display_orders[utility.getGlobalValue ("terms.json.nucleotide_gtr")]);

    utility.ForEachPair (fubar.filter_specification, "_key_", "_value_",
        'selection.io.json_store_branch_attribute(fubar.json, terms.json.nucleotide_gtr, terms.branch_length, fubar.display_orders[terms.json.nucleotide_gtr],
                                         _key_,
                                         selection.io.extract_branch_info((fubar.gtr_results[terms.branch_length])[_key_], "selection.io.branch.length"));');


}
else {
    fubar.RunPrompts (fubar.prompts);
    namespace fubar {
        doGTR ("fubar");
    }
    fubar.cache [terms.fubar.cache.gtr] = fubar.gtr_results;
    io.WriteCacheToFile (fubar.path.cache, fubar.cache);
}


utility.ForEachPair (fubar.gtr_results [terms.branch_length], "_key_", "_value_",
                        '
                             _value_lengths_ = utility.Map (_value_, "_info_", "_info_[utility.getGlobalValue (\\"terms.fit.MLE\\")]");
                             io.ReportProgressMessageMD ("fubar", "nuc-fit", "* Tree length (expected substitutions/site) for partition " + (1 + _key_) + " : " + Format (+_value_lengths_, 8, 3));
                        ');
selection.io.stopTimer (fubar.json [terms.json.timers], "Nucleotide Fit");


//----------------------------------------------------------------------------
// PHASE 2: compute LF on the grid
//----------------------------------------------------------------------------
selection.io.startTimer (fubar.json [terms.json.timers], "Grid Calculations", 2);

if (utility.Has (fubar.cache, terms.fubar.cache.grid, "Matrix") && utility.Has (fubar.cache, terms.fubar.cache.conditionals, "AssociativeList")) {
    fubar.prompts["grid"] = FALSE;
    fubar.grid.matrix = fubar.cache [terms.fubar.cache.grid];
    fubar.conditionals.raw =  fubar.cache [terms.fubar.cache.conditionals];
    fubar.run_settings["grid size"] = Sqrt(Rows ((fubar.cache [terms.fubar.cache.conditionals])["conditionals"]));
    io.ReportProgressMessageMD                  ("fubar", "codon_fit", "_Loaded cached phylogenetic likelihood function on the grid_ of dimension `(fubar.run_settings['grid size'])`x`(fubar.run_settings['grid size'])` points");
} else {
    // define FUBAR model

    fubar.RunPrompts (fubar.prompts);
    fubar.grid.matrix                            = fubar.DefineAlphaBetaGrid (fubar.run_settings["grid size"], fubar.run_settings["non-zero"]);

    io.ReportProgressMessageMD                  ("fubar", "codon_fit", "Computing the phylogenetic likelihood function on the grid ");
    io.ReportProgressMessageMD                  ("fubar", "codon_fit", "* Determining appropriate tree scaling based on the best score from a  `fubar.run_settings['grid size']` x `fubar.run_settings['grid size']` rate grid");

    fubar.mg_rev = model.generic.DefineModel("models.codon.MG_REV.ModelDescription",
        "fubar.codon_model", {
            "0": "terms.local",
            "1": fubar.codon_data_info [terms.code]
        },
        fubar.filter_names,
        None);

    fubar.mg_rev [utility.getGlobalValue("terms.model.set_branch_length")] = "fubar.scalers.SetBranchLength";

    fubar.parameter.scalers = {
        terms.fubar.branch_length_scaler.alpha : "fubar.scaler.alpha",
        terms.fubar.branch_length_scaler.beta : "fubar.scaler.beta"
    };

    utility.Extend ((fubar.mg_rev [terms.parameters])[terms.global], fubar.parameter.scalers);
    parameters.DeclareGlobal (fubar.parameter.scalers, {});

    fubar.model_assignment = {
        "default": fubar.mg_rev
    };

    fubar.model_id_to_object = {
        "fubar.codon_model": fubar.mg_rev
    };

    fubar.trees.names = utility.Map (utility.Range (fubar.partition_count, 1, 1), "_index_", "'fubar.codon_tree_' + _index_");
    fubar.lf.components = {fubar.partition_count * 2, 1};
    utility.ForEachPair (fubar.filter_names, "_index_", "_filter_",
        '
            fubar.lf.components [2*(0+_index_)] = _filter_;
            fubar.lf.components [2*(0+_index_) + 1] = fubar.trees.names[_index_];
            model.ApplyModelToTree(fubar.trees.names [_index_], fubar.trees[_index_], fubar.model_assignment, None);
        ');


    LikelihoodFunction fubar.lf.codon = (fubar.lf.components);


    estimators.ApplyExistingEstimates  ("fubar.lf.codon", fubar.model_id_to_object, fubar.gtr_results, None);
    estimators.TraverseLocalParameters ("fubar.lf.codon", fubar.model_id_to_object, "fubar.scalers.Constrain");

    fubar.pass1 = Max (fubar.ComputeOnGrid  ("fubar.lf.codon",
                         fubar.grid.MatrixToDict (fubar.grid.matrix),
                        "fubar.pass1.evaluator",
                        "fubar.pass1.result_handler"),
                        1);

    fubar.best_scaler = fubar.grid.matrix[0 + fubar.pass1["key"]][-1];

    parameters.SetValues ((fubar.grid.MatrixToDict (fubar.best_scaler) )[0]);

    io.ReportProgressMessageMD                  ("fubar", "codon_fit", "* Best scaling achieved for \n\t* synonymous rate = " + Format (fubar.best_scaler[0], 6, 3) +
                                                                                                   "\n\t* non-synonymous rate = " + Format (fubar.best_scaler[1], 6, 3));

    fubar.pass1.mle = estimators.ExtractMLEs( "fubar.lf.codon", fubar.model_id_to_object);

    estimators.TraverseLocalParameters ("fubar.lf.codon", fubar.model_id_to_object, "fubar.scalers.Unconstrain");
    estimators.ApplyExistingEstimates  ("fubar.lf.codon", fubar.model_id_to_object, fubar.pass1.mle , None);
    estimators.TraverseLocalParameters ("fubar.lf.codon", fubar.model_id_to_object, "fubar.scalers.Constrain");

    io.ReportProgressMessageMD                  ("fubar", "codon_fit", "* Computing conditional site likelihoods on a `fubar.run_settings["grid size"]` x `fubar.run_settings["grid size"]` rate grid");

    fubar.conditionals.raw = fubar.ComputeOnGrid  ("fubar.lf.codon",
                         fubar.grid.MatrixToDict (fubar.grid.matrix),
                        "fubar.pass2.evaluator",
                        "fubar.pass1.result_handler");


    fubar.cache[terms.fubar.cache.conditionals] = fubar.ConvertToConditionals (fubar.conditionals.raw);

    fubar.cache[terms.fubar.cache.grid] = fubar.grid.matrix;
    fubar.cache - terms.fubar.cache.mcmc; // overwrite old MCMC cache
    io.WriteCacheToFile (fubar.path.cache, fubar.cache);


}
selection.io.stopTimer (fubar.json [terms.json.timers], "Grid Calculations");

//----------------------------------------------------------------------------
// PHASE 3: Posterior estimation
//----------------------------------------------------------------------------

if (utility.Has (fubar.run_settings, "method", "String")) {
    fubar.prompts["method"] = FALSE;
} else {
    fubar.RunPrompts(fubar.prompts);
    fubar.cache [terms.fubar.cache.settings] = fubar.run_settings;
}


selection.io.startTimer (fubar.json [terms.json.timers], "Posterior estimation", 3);

if (fubar.run_settings["method"] != terms.fubar.methods.VB0) {
    if (utility.Has (fubar.cache, terms.fubar.cache.mcmc, "AssociativeList")) {
        io.ReportProgressMessageMD                  ("fubar", "mcmc", "_Loaded cached results from `Abs(fubar.cache[terms.fubar.cache.mcmc])` chains_");
    } else {
        fubar.RunPrompts (fubar.prompts);
        if (fubar.run_settings["method"] == terms.fubar.methods.MH) {
            fubar.cache[terms.fubar.cache.mcmc] = fubar.RunMCMC  (fubar.run_settings,
                                                                  fubar.cache[terms.fubar.cache.grid],
                                                                  fubar.cache[terms.fubar.cache.conditionals],
                                                                  "fubar.pass1.result_handler",
                                                                  "fubar"
                                                                  );
        } else {
            fubar.cache[terms.fubar.cache.mcmc] = fubar.RunCollapsedGibbs  (fubar.run_settings,
                                                                  fubar.cache[terms.fubar.cache.grid],
                                                                  fubar.cache[terms.fubar.cache.conditionals],
                                                                  "fubar"
                                                                  );
        }
        io.WriteCacheToFile (fubar.path.cache, fubar.cache);
    }
} else {
    if (utility.Has (fubar.cache, terms.fubar.cache.posterior, "Matrix")) { // for VB this is going to be the set of posterior weights
        io.ReportProgressMessageMD                  ("fubar", "mcmc", "_Loaded cached estimated posterior probabilities on the grid");
    } else {
        fubar.RunPrompts (fubar.prompts);
        fubar.cache[terms.fubar.cache.posterior] = fubar.RunVariationalBayes  (fubar.run_settings,
                                                              fubar.cache[terms.fubar.cache.grid],
                                                              fubar.cache[terms.fubar.cache.conditionals],
                                                              "fubar"
                                                              );
        io.WriteCacheToFile (fubar.path.cache, fubar.cache);
    }
}

selection.io.stopTimer (fubar.json [terms.json.timers], "Posterior estimation");

//----------------------------------------------------------------------------
// PHASE 4: Processing
//----------------------------------------------------------------------------


io.ReportProgressMessageMD                  ("fubar", "reporting", "Tabulating site-level results");

namespace fubar {
    sites   = (json [utility.getGlobalValue ("terms.json.input")])[utility.getGlobalValue ("terms.json.sites")];
    grid_points = Rows (cache['grid']);
    positive_selection_stencil = {grid_points,sites} ["(cache['grid'])[_MATRIX_ELEMENT_ROW_][0]<(cache['grid'])[_MATRIX_ELEMENT_ROW_][1]"];
    negative_selection_stencil = {grid_points,sites} ["(cache['grid'])[_MATRIX_ELEMENT_ROW_][0]>(cache['grid'])[_MATRIX_ELEMENT_ROW_][1]"];
    diag_alpha                 = {grid_points,grid_points}["(cache['grid'])[_MATRIX_ELEMENT_ROW_][0]*(_MATRIX_ELEMENT_ROW_==_MATRIX_ELEMENT_COLUMN_)"];
    diag_beta                  = {grid_points,grid_points}["(cache['grid'])[_MATRIX_ELEMENT_ROW_][1]*(_MATRIX_ELEMENT_ROW_==_MATRIX_ELEMENT_COLUMN_)"];

    if (run_settings["method"] == ^"terms.fubar.methods.MH") {
        table_headers = {{"alpha", "Mean posterior synonymous substitution rate at a site"}
                         {"beta", "Mean posterior non-synonymous substitution rate at a site"}
                         {"beta-alpha", "Mean posterior beta-alpha"}
                         {"Prob[alpha>beta]", "Posterior probability of negative selection at a site"}
                         {"Prob[alpha<beta]", "Posterior probability of positive selection at a site"}
                         {"BayesFactor[alpha<beta]", "Empiricial Bayes Factor for positive selection at a site"}
                         {"PSRF", "Potential scale reduction factor - an MCMC mixing measure"}
                         {"Neff", "Estimated effective sample site for Prob [alpha<beta]"}};

        table_screen_output  = {{"Codon", "Partition", "alpha", "beta", "N.eff", "Posterior prob for positive selection"}};
        report.positive_site = {{"" + (1+filter_info[s]),
                                        partition_index + 1,
                                        Format(partition_results[s][0],10,3),
                                        Format(partition_results[s][1],10,3),
                                        Format(partition_results[s][7],10,3),
                                        "Pos. posterior = " + Format(partition_results[s][4],6,4)}};
    } else {
        table_headers = {{"alpha", "Mean posterior synonymous substitution rate at a site"}
                         {"beta", "Mean posterior non-synonymous substitution rate at a site"}
                         {"beta-alpha", "Mean posterior beta-alpha"}
                         {"Prob[alpha>beta]", "Posterior probability of negative selection at a site"}
                         {"Prob[alpha<beta]", "Posterior probability of positive selection at a site"}
                         {"BayesFactor[alpha<beta]", "Empiricial Bayes Factor for positive selection at a site"}};

         table_screen_output  = {{"Codon", "Partition", "alpha", "beta", "Posterior prob for positive selection"}};
         report.positive_site = {{"" + (1+filter_info[s]),
                                        partition_index + 1,
                                        Format(partition_results[s][0],10,3),
                                        Format(partition_results[s][1],10,3),
                                        "Pos. posterior = " + Format(partition_results[s][4],6,4)}};
    }

    if (run_settings["method"] != ^"terms.fubar.methods.VB0") {
        samples = run_settings["samples"];
        chains  = run_settings["chains"];

        results.log_L   = {1,samples};
        results.samples = {samples,grid_points};


        per_chain      = samples $ chains;


        from = 0;
        to   = per_chain;



        positive_ks                 = {};
        negative_ks                 = {};
        posterior_mean_alpha        = {};
        posterior_mean_beta         = {};
        denominators                = {};
        posteriors                  = {};

        for (chain_id = 0; chain_id < chains; chain_id += 1) {
            io.ReportProgressBar                  ("PROCESSING", "Samples from chain " + (chain_id + 1));

            /* now, for each posterior sample k of grid weights, i.e. (alpha_i, beta_i) -> weight_ik
               and for each site, s, we compute (up to a factor C), which will be divided out in the Bayes' formula computaion, the
                    P_ks = \sum_i Prob (site s | (alpha_i, beta_i)) Prob (alpha_i, beta_i) ~ Prob (site s) for sample k

                this matrix will have dimension (# samples x # sites)
            */

            P_ks = ((cache[utility.getGlobalValue("terms.fubar.cache.mcmc")])[chain_id])["weights"] * (cache[utility.getGlobalValue("terms.fubar.cache.conditionals")])["conditionals"];
            denominators[chain_id] = P_ks;

            /*
                the next two matrices store posterior calculations for

                positive_ks = P {dN > dS @ site s for grid sample k}
                negative_ks = P {dN < dS @ site s for grid sample k}
            */

            grid_samples = ((cache[utility.getGlobalValue("terms.fubar.cache.mcmc")])[chain_id])["weights"];

            if (run_settings["method"] == ^"terms.fubar.methods.MH") {
                 logL_samples = ((cache[utility.getGlobalValue("terms.fubar.cache.mcmc")])[chain_id])["likelihoods"];
            }

            positive_ks [chain_id] = grid_samples *
                          (positive_selection_stencil $ (cache[utility.getGlobalValue("terms.fubar.cache.conditionals")])["conditionals"]) / P_ks;

            negative_ks [chain_id] = grid_samples *
                          (negative_selection_stencil $ (cache[utility.getGlobalValue("terms.fubar.cache.conditionals")])["conditionals"]) / P_ks;


            /*
                the next two matrices compute posterior mean alpha and beta values for site s and grid sample k

            */

            posterior_mean_alpha [chain_id]     =            (grid_samples * diag_alpha *
                                                          (cache[utility.getGlobalValue("terms.fubar.cache.conditionals")])["conditionals"]) / P_ks;
            posterior_mean_beta [chain_id]      =            (grid_samples * diag_beta *
                                                          (cache[utility.getGlobalValue("terms.fubar.cache.conditionals")])["conditionals"]) / P_ks;

            if (run_settings["method"] == ^"terms.fubar.methods.MH") {
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
        posterior_mean_over_grid                 = cache[^"terms.fubar.cache.posterior"];
        posterior_mean_over_grid_T               = Transpose (posterior_mean_over_grid);
        cache[terms.fubar.cache.posterior]       = posterior_mean_over_grid;


         P_ks = posterior_mean_over_grid_T * (cache[utility.getGlobalValue("terms.fubar.cache.conditionals")])["conditionals"];


         posterior_mean_alpha     =             (posterior_mean_over_grid_T * diag_alpha *
                                                          (cache[utility.getGlobalValue("terms.fubar.cache.conditionals")])["conditionals"]) / P_ks;
         posterior_mean_beta       =            (posterior_mean_over_grid_T * diag_beta *
                                                          (cache[utility.getGlobalValue("terms.fubar.cache.conditionals")])["conditionals"]) / P_ks;

         positive_ks  = posterior_mean_over_grid_T *
                           (positive_selection_stencil $ (cache[utility.getGlobalValue("terms.fubar.cache.conditionals")])["conditionals"]) / P_ks;

         negative_ks  = posterior_mean_over_grid_T *
                          (negative_selection_stencil $ (cache[utility.getGlobalValue("terms.fubar.cache.conditionals")])["conditionals"]) / P_ks;

   }

    /* compute the posterior mean of point loadings */


    weight_non_positive         = +posterior_mean_over_grid [posterior_mean_over_grid ["(cache['grid'])[_MATRIX_ELEMENT_ROW_][0]>=(cache['grid'])[_MATRIX_ELEMENT_ROW_][1]"]];
    /* total posterior weight assigned to negative or neutral selective regimes */


    headers.printed = FALSE;

    site_results = {};
    i = 0;
    s = 0; // reset to ensure good re-entrant behavior

    report.sites_found = {};
    report.posteriors  = {};

    chain_iterator = utility.Range (chains, 0, 1);

    for (partition_index = 0; partition_index < partition_count; partition_index += 1) {
        filter_info = (filter_specification [partition_index])[utility.getGlobalValue ("terms.data.coverage")];
        sites_in_partition = utility.Array1D (filter_info);

        partition_results    = {sites_in_partition, 8};
        partition_posteriors = {};


        for (s = 0; s < sites_in_partition; s += 1) {

            pp = posterior_mean_over_grid $ ((cache[utility.getGlobalValue("terms.fubar.cache.conditionals")])["conditionals"])[-1][i];
            partition_posteriors [s] = Transpose (pp * (1/(+pp)));

            if (run_settings["method"] != utility.getGlobalValue ("terms.fubar.methods.VB0")) {
                partition_results[s][0] = fubar.ComputeRandNeff (
                    utility.Map (chain_iterator, "_value_", "((`&posterior_mean_alpha`)[_value_])[-1][`&i`]")
                )[0];
                partition_results[s][1] = fubar.ComputeRandNeff (
                    utility.Map (chain_iterator, "_value_", "((`&posterior_mean_beta`)[_value_])[-1][`&i`]")
                )[0];
                partition_results[s][3] = fubar.ComputeRandNeff (
                    utility.Map (chain_iterator, "_value_", "((`&negative_ks`)[_value_])[-1][`&i`]")
                )[0];

                pos_selection = fubar.ComputeRandNeff (
                    utility.Map (chain_iterator, "_value_", "((`&positive_ks`)[_value_])[-1][`&i`]")
                );
                partition_results[s][4] = pos_selection[0];

                if (run_settings["method"] == utility.getGlobalValue ("terms.fubar.methods.MH")) {
                    partition_results[s][6] = pos_selection[1];
                    partition_results[s][7] = pos_selection[2];
                }

            } else {
                if (run_settings["method"] == utility.getGlobalValue ("terms.fubar.methods.VB0")) {
                    partition_results [s][0] = posterior_mean_alpha[i];
                    partition_results [s][1] = posterior_mean_beta[i];
                    partition_results [s][3] = negative_ks[i];
                    partition_results [s][4] = positive_ks[i];
               }
            }


            partition_results[s][2] = partition_results[s][1] - partition_results[s][0];
            if (weight_non_positive > 0 && weight_non_positive < 1) {
                partition_results [s][5] = partition_results[s][4] / (1-partition_results[s][4]) / (1-weight_non_positive) * weight_non_positive;
            } else {
                partition_results[s][5] = 1;
            }

            if (partition_results[s][4] >= run_settings["posterior"]) {
                if (Abs(report.sites_found) == 0) {
                    fprintf (stdout, io.FormatTableRow (table_screen_output,table_output_options));
                    table_output_options[^"terms.table_options.header"] = FALSE;
                }
                fprintf (stdout, io.FormatTableRow (report.positive_site,table_output_options));
                report.sites_found + (1-partition_results[s][4]);
            }

            i+=1;
        }
        site_results [partition_index] = partition_results;
        report.posteriors [partition_index] = partition_posteriors;
    }

    sites_found = Abs(report.sites_found);
    if (sites_found == 0) {
        console.log ("----\n## FUBAR inferred no sites under subject to positive selection at posterior probability >= " + run_settings["posterior"]);
    } else {
        console.log ("----\n## FUBAR inferred `sites_found` sites subject to diversifying positive selection at posterior probability >= " + run_settings["posterior"]);
        mean_fp = (+report.sites_found);
        ci = fubar.ComputeENFP_CI(report.sites_found, 0.05);
        console.log ("Of these, " + Format (mean_fp, 5, 2) + " are expected to be false positives (95% confidence interval of " + Join ("-", ci) + " )");
     }
 };


fubar.json [terms.fubar.cache.settings] = fubar.run_settings;
fubar.json [terms.fit.MLE] = {terms.json.headers   : fubar.table_headers,
                               terms.json.content : fubar.site_results };
fubar.json [terms.fubar.posterior] = fubar.report.posteriors;

namespace fubar {
    grid_with_weights = {grid_points, 3};
    for (i = 0; i < grid_points; i+=1) {
        grid_with_weights [i][0] = grid.matrix[i][0];
        grid_with_weights [i][1] = grid.matrix[i][1];
        grid_with_weights [i][2] = posterior_mean_over_grid[i];
    }
}

fubar.json [terms.fubar.grid]           = fubar.grid_with_weights;
selection.io.stopTimer (fubar.json [terms.json.timers], "Overall");

io.SpoolJSON (fubar.json, fubar.codon_data_info[terms.json.json]);

return fubar.json;

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

lfunction fubar.ComputeENFP_CI (p_i,sig_level) {
    N = Abs (p_i);

    PDF_old = {{1-p_i[0],p_i[0]}};
    PDF = PDF_old;

    for (i = 1; i < N; i+=1) {
        PDF = {1,i+2};
        PDF[0] = PDF_old[0] * (1-p_i[i]);
        for (X = 1; X < i+1; X+=1) {
            PDF[X] = PDF_old[X] * (1-p_i[i]) + PDF_old[X-1] * (p_i[i]);
        }
        PDF[i+1] = PDF_old[i] * p_i[i];
        PDF_old = PDF;
    }

    sum = PDF[0];
    _idx = 0;
    while (sum < sig_level/2) {
        _idx+=1;
        sum += PDF[_idx];
    }
    lb = _idx;

    while (sum < 1-sig_level/2) {
        _idx+=1;
        sum += PDF[_idx];
    }

    return {{lb__, _idx__}}
}




//----------------------------------------------------------------------------

function     fubar.RunPrompts (prompts) {

 
     if (prompts["grid"]) {
        KeywordArgument ("grid", "The number of grid points", "20");
        fubar.run_settings["grid size"] = io.PromptUser ("> Number of grid points per dimension (total number is D^2)",fubar.run_settings["grid size"],5,50,TRUE);
        prompts["grid"] = FALSE;
    }

    if (prompts["method"]) {
        // TODO: can we replace the default with `variational-bayes`?
        KeywordArgument ("method", "Inference method to use", "`terms.fubar.methods.VB0`");
        fubar.run_settings["method"] = io.SelectAnOption  ({
                                                                terms.fubar.methods.MH : "Full Metropolis-Hastings MCMC algorithm (slowest, original 2013 paper implementation)",
                                                                terms.fubar.methods.CG : "Collapsed Gibbs sampler (intermediate speed)",
                                                                terms.fubar.methods.VB0 : "0-th order Variational Bayes approximations (fastest, recommended default)"
                                                            }, "Posterior estimation method");
        prompts["method"] = FALSE;
     }

    if (prompts["chain"]) {
        if (fubar.run_settings["method"] ==  terms.fubar.methods.MH) {
            KeywordArgument ("chains", "How many MCMC chains to run", fubar.run_settings["chains"]);
            fubar.run_settings["chains"] = io.PromptUser ("> Number of MCMC chains to run",fubar.run_settings["chains"],2,20,TRUE);
        } else {
            fubar.run_settings["chains"] = 1;
        }
        if (fubar.run_settings["method"] !=  terms.fubar.methods.VB0) {
            KeywordArgument ("chain-length", "MCMC chain length", fubar.run_settings["chain-length"]);
            fubar.run_settings["chain-length"] = io.PromptUser ("> The length of each chain",fubar.run_settings["chain-length"],5e3,5e7,TRUE);
            KeywordArgument ("burn-in", "MCMC chain burn in", fubar.run_settings["chain-length"]$2);
            fubar.run_settings["burn-in"] = io.PromptUser ("> Use this many samples as burn-in",fubar.run_settings["chain-length"]$2,fubar.run_settings["chain-length"]$20,fubar.run_settings["chain-length"]*95$100,TRUE);
            KeywordArgument ("samples", "MCMC samples to draw", fubar.run_settings["samples"]);
            fubar.run_settings["samples"] = io.PromptUser ("> How many samples should be drawn from each chain",fubar.run_settings["samples"],50,fubar.run_settings["chain-length"]-fubar.run_settings["burn-in"],TRUE);
        }
        KeywordArgument ("concentration_parameter", "The concentration parameter of the Dirichlet prior", fubar.run_settings["concentration"]);
        fubar.run_settings["concentration"] = io.PromptUser  ("> The concentration parameter of the Dirichlet prior",fubar.run_settings["concentration"],0.001,1,FALSE);
        prompts["chain"] = FALSE;
    }

   if (prompts["non-zero"]) {
    
        KeywordArgument ("non-zero", "Enforce non-zero synonymous rates on the grid to enable dN/dS calculations", "No");
        
        fubar.run_settings["non-zero"] = io.SelectAnOption  ({
                                                                "Yes" : "The grid will exclude zero synonymous rates, resulting in finite dN/dS values",
                                                                "No"  : "The grid will include zero synonymous rates, resulting in infinite dN/dS values"
                                                             }, "Enforce non-zero synonymous rates on the grid") == "Yes";
                                                             
        prompts["non-zero"] = FALSE;
     }

}

//------------------------------------------------------------------------------


lfunction fubar.grid.MatrixToDict (grid) {
    return utility.Map (utility.MatrixToListOfRows (grid), "_value_",
                                                                '{  terms.fubar.branch_length_scaler.alpha : {
                                                                            terms.id : fubar.parameter.scalers [ terms.fubar.branch_length_scaler.alpha ],
                                                                            terms.fit.MLE : _value_[0]
                                                                        },
                                                                    terms.fubar.branch_length_scaler.beta :  {
                                                                            terms.id : fubar.parameter.scalers [ terms.fubar.branch_length_scaler.beta ],
                                                                            terms.fit.MLE : _value_[1]
                                                                        }
                                                                 }');
}

//------------------------------------------------------------------------------


lfunction fubar.scalers.SetBranchLength (model, value, parameter) {
    io.CheckAssertion ('Type (`&value`) == "Number"', "Internal error in fubar.scalers.SetBranchLength");
    local = (model[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.local")];
    alpha = local [utility.getGlobalValue("terms.parameters.synonymous_rate")];
    beta = local [utility.getGlobalValue("terms.parameters.nonsynonymous_rate")];
    parameters.SetConstraint (beta, alpha, "");
    utility.ExecuteInGlobalNamespace ("FindRoot (`&s`,(" + model[utility.getGlobalValue("terms.model.branch_length_string")] + ")-" + 3*value + "," + alpha + ",0,10000)");
    parameters.RemoveConstraint (beta);
    parameters.SetValue ("`parameter`.`alpha`", s);
    parameters.SetValue ("`parameter`.`beta`", ^beta);
    return 1;
}

//------------------------------------------------------------------------------

lfunction fubar.scalers.Constrain (tree_name, node_name, model_description) {
    parameters.SetProprtionalConstraint (tree_name + "." + node_name + "." + (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.parameters.synonymous_rate")],
                              (model_description[utility.getGlobalValue ("terms.global")])[utility.getGlobalValue ("terms.fubar.branch_length_scaler.alpha")]);

    parameters.SetProprtionalConstraint (tree_name + "." + node_name + "." + (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.parameters.nonsynonymous_rate")],
                              (model_description[utility.getGlobalValue ("terms.global")])[utility.getGlobalValue ("terms.fubar.branch_length_scaler.beta")]);

}

//------------------------------------------------------------------------------

lfunction fubar.scalers.Unconstrain (tree_name, node_name, model_description) {
    parameters.RemoveConstraint (tree_name + "." + node_name + "." + (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.parameters.synonymous_rate")]);
    parameters.RemoveConstraint (tree_name + "." + node_name + "." + (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.parameters.nonsynonymous_rate")]);
}



//------------------------------------------------------------------------------------------------//

lfunction fubar.DefineAlphaBetaGrid (one_d_points, non_zero) {


    one_d_points    = Max (one_d_points, 5);

    alphaBetaGrid = {one_d_points^2,2}; // (alpha, beta) pair
    oneDGrid      = {one_d_points,1};

    neg_sel         = 0.7;
    neg_sel_points  = ((one_d_points)*neg_sel+0.5)$1;
    pos_sel_points  = (one_d_points-1)*(1-neg_sel)$1;

    if (neg_sel_points + pos_sel_points != one_d_points) {
        pos_sel_points = one_d_points - neg_sel_points;
    }
    _neg_step = 1/neg_sel_points;
    for (_k = 0; _k < neg_sel_points; _k += 1) {
        oneDGrid [_k][0] =  _neg_step * _k;
    }
    oneDGrid [neg_sel_points-1][0] = 1;
    _pos_step = 49^(1/3)/pos_sel_points;
    for (_k = 1; _k <= pos_sel_points; _k += 1) {
        oneDGrid [neg_sel_points+_k-1][0] = 1+(_pos_step*_k)^3;
    }

    _p = 0;
    if (non_zero) {
       min_value = Max (1e-3, oneDGrid[0]);
        for (_r = 0; _r < one_d_points; _r += 1) {
            for (_c = 0; _c < one_d_points; _c += 1) {
               alphaBetaGrid[_p][0] = Max (oneDGrid[_r], min_value);
               alphaBetaGrid[_p][1] = oneDGrid[_c];
               _p += 1;
            }
        }
        
    } else {
        for (_r = 0; _r < one_d_points; _r += 1) {
            for (_c = 0; _c < one_d_points; _c += 1) {
               alphaBetaGrid[_p][0] = oneDGrid[_r];
               alphaBetaGrid[_p][1] = oneDGrid[_c];
               _p += 1;
            }
        }       
    }

    return alphaBetaGrid;
}

