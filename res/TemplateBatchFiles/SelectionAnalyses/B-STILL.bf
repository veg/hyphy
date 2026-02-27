RequireVersion  ("2.5.93");

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

lfunction fubar.compute_ebf (p, prior) {
    if (prior > 0 && prior < 1) {
        if (p >= 1) { return 1e10; }
        return (p / (1-p)) / (prior / (1-prior));
    }
    return 0;
}

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
        sub_scaler   = "sub-scaler";
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
    "non-zero" : TRUE,
    "ebf" : TRUE,
    "radius-threshold" : TRUE
};

fubar.run_settings = {
    "grid size" : 20,
    "chains" : 5,
    "chain-length" : 2e6,
    "burn-in" : 1e6,
    "samples" : 100,
    "concentration" : 0.5,
    "posterior" : 0.9,
    "non-zero" : FALSE,
    "ebf" : 10,
    "radius-threshold" : 0.5
};



/*------------------------------------------------------------------------------*/

fubar.analysis_description = {terms.io.info :

    "Perform a B-STILL (Bayesian Significance Test of Invariant Low Likelihoods) 
    analysis to detect invariant sites (alpha=beta=0) and quantify their 
    posterior probabilities and Empirical Bayes Factors. This is a modified 
    version of the standard FUBAR analysis that uses a denser grid around 
    zero and reports the probability of a site being effectively invariant.",

           terms.io.version : "1.0 (B-STILL)",
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

namespace fubar {
    LoadFunctionLibrary ("modules/shared-load-file.bf");
    LoadFunctionLibrary ("modules/grid_compute.ibf");
    load_file ({utility.getGlobalValue("terms.prefix"): "fubar", utility.getGlobalValue("terms.settings") : {utility.getGlobalValue("terms.settings.branch_selector") : "selection.io.SelectAllBranches"}});

}

selection.io.startTimer (fubar.json [terms.json.timers], "Overall", 0);

fubar.path.base = (fubar.json [terms.json.input])[terms.json.file];

KeywordArgument ("cache",   "Save FUBAR cache to [default is alignment+.FUBAR-inv.cache]", fubar.path.base + ".FUBAR-inv.cache");


/*------------------------------------------------------------------------------
    Continued Analysis Setup
*/

fubar.path.cache = io.PromptUserForString ("Save FUBAR cache to");
fubar.cache = io.LoadCacheFromFile (fubar.path.cache);

KeywordArgument ("output",   "Save FUBAR results (JSON) to [default is alignment+.FUBAR-inv.json]", fubar.codon_data_info[terms.data.file] + ".FUBAR-inv.json");
fubar.codon_data_info [terms.json.json] = io.PromptUserForFilePath ("Save the resulting JSON file to");

fubar.table_output_options = {terms.table_options.header : TRUE, terms.table_options.minimum_column_width: 16, terms.table_options.align : "center"};

console.log ( "> FUBAR will write cache and result files to _`fubar.path.cache`_ and _`fubar.codon_data_info[terms.json.json]`_, respectively 

");

//----------------------------------------------------------------------------
// PHASE 1: nucleotide fit
//----------------------------------------------------------------------------


if (utility.Has (fubar.cache, terms.fubar.cache.settings, "AssociativeList")) {
    fubar.run_settings = fubar.cache [terms.fubar.cache.settings];
    fubar.RunPrompts ({"ebf" : TRUE, "radius-threshold" : TRUE});
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


lfunction fubar.report_partition_length (key, value, arguments) {
    _sum = 0;
    _mle_term = ^"terms.fit.MLE";
    
    for (_n, _v; in; value) {
         _sum += _v[_mle_term];
    }
    
    io.ReportProgressMessageMD ("fubar", "nuc-fit", "* Tree length (expected substitutions/site) for partition " + (1 + key) + " : " + Format (_sum, 8, 3));
    
    return _sum;
}


fubar.tree_scaler = 0;

for (_key_, _value_; in; fubar.gtr_results [terms.branch_length]) {
    fubar.tree_scaler = Max (fubar.tree_scaler, fubar.report_partition_length (_key_, _value_, None));
}


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
    fubar.single_sub_scale = fubar.cache[terms.fubar.cache.sub_scaler];

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
    
    alpha = 1;
    beta  = 1;
    

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

    io.ReportProgressMessageMD                  ("fubar", "codon_fit", "* Best scaling achieved for 
	* synonymous rate = " + Format (fubar.best_scaler[0], 6, 3) +
                                                                                                   "
	* non-synonymous rate = " + Format (fubar.best_scaler[1], 6, 3));

    fubar.pass1.mle = estimators.ExtractMLEs( "fubar.lf.codon", fubar.model_id_to_object);

    estimators.TraverseLocalParameters ("fubar.lf.codon", fubar.model_id_to_object, "fubar.scalers.Unconstrain");
    estimators.ApplyExistingEstimates  ("fubar.lf.codon", fubar.model_id_to_object, fubar.pass1.mle , None);
    estimators.TraverseLocalParameters ("fubar.lf.codon", fubar.model_id_to_object, "fubar.scalers.Constrain");

    io.ReportProgressMessageMD                  ("fubar", "codon_fit", "* Computing conditional site likelihoods on a `fubar.run_settings["grid size"]` x `fubar.run_settings["grid size"]` rate grid");

    fubar.conditionals.raw = fubar.ComputeOnGrid  ("fubar.lf.codon",
                         fubar.grid.MatrixToDict (fubar.grid.matrix),
                        "fubar.pass2.evaluator",
                        "fubar.pass1.result_handler");

    fubar.branch_length_expression = fubar.mg_rev[terms.model.branch_length_string];
    
    fubar.subs = {};
    
    for (p, v; in; fubar.pass1.mle[terms.global]) {
        fubar.subs[v[terms.id]] = v[terms.fit.MLE];
    }
    
    
    fubar.single_sub_scale = "3*(" + Simplify (fubar.branch_length_expression,fubar.subs) + ")";
    fubar.cache[terms.fubar.cache.conditionals] = fubar.ConvertToConditionals (fubar.conditionals.raw);

    fubar.cache[terms.fubar.cache.grid] = fubar.grid.matrix;
    fubar.cache[terms.fubar.cache.sub_scaler] = fubar.single_sub_scale;
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
    
    // Invariance stencils
    invariant_stencil  = {grid_points,sites} ["(cache['grid'])[_MATRIX_ELEMENT_ROW_][0]==0 && (cache['grid'])[_MATRIX_ELEMENT_ROW_][1]==0"];
    
    alpha_zero_stencil = {grid_points,sites} ["(cache['grid'])[_MATRIX_ELEMENT_ROW_][0]==0"];
    beta_zero_stencil  = {grid_points,sites} ["(cache['grid'])[_MATRIX_ELEMENT_ROW_][1]==0"];

    proximal_radius = run_settings["radius-threshold"];
    
    function check_radius (x, y) {
        alpha = x;
        beta = y;
        _rv = Eval (single_sub_scale);
        return _rv;
    }
	
    proximal_stencil   = {grid_points,sites};
    
    
    for (r = 0; r < grid_points; r+=1) {
        if (check_radius ((cache['grid'])[r][0], (cache['grid'])[r][1]) <= proximal_radius) {
            for (c = 0; c < sites; c+=1) {
                proximal_stencil[r][c] = 1;
            }
        }
    }
    
    
    diag_alpha                 = {grid_points,grid_points}["(cache['grid'])[_MATRIX_ELEMENT_ROW_][0]*(_MATRIX_ELEMENT_ROW_==_MATRIX_ELEMENT_COLUMN_)"];
    diag_beta                  = {grid_points,grid_points}["(cache['grid'])[_MATRIX_ELEMENT_ROW_][1]*(_MATRIX_ELEMENT_ROW_==_MATRIX_ELEMENT_COLUMN_)"];

    table_headers = {{"alpha", "Mean posterior synonymous substitution rate at a site"}
                     {"beta", "Mean posterior non-synonymous substitution rate at a site"}
                     {"Prob[alpha=beta=0]", "Posterior probability of alpha=beta=0"}
                     {"Prob[alpha=0]", "Posterior probability of alpha=0"}
                     {"Prob[beta=0]", "Posterior probability of beta=0"}
                     {"Prob[alpha,beta~0]", "Posterior probability of alpha and beta within a radius of " + proximal_radius + " of 0"}
                     {"Prob[alpha<beta]", "Posterior probability of positive selection at a site"}
                     {"PSRF", "Potential scale reduction factor - an MCMC mixing measure"}
                     {"Neff", "Estimated effective sample site for Prob [alpha<beta]"}
                     {"EBF[alpha=beta=0]", "Empirical Bayes Factor for alpha=beta=0"}
                     {"EBF[alpha=0]", "Empirical Bayes Factor for alpha=0"}
                     {"EBF[beta=0]", "Empirical Bayes Factor for beta=0"}
                     {"EBF[alpha,beta~0]", "Empirical Bayes Factor for alpha and beta within a radius of " + proximal_radius + " of 0"}};

    /*
    if (run_settings["method"] == ^"terms.fubar.methods.MH") {
        table_headers + {{"PSRF", "Potential scale reduction factor - an MCMC mixing measure"}};
        table_headers + {{"Neff", "Estimated effective sample site for Prob [alpha<beta]"}};
    }
    */

    table_screen_output  = {{"Codon", "Partition", "alpha", "beta", "P[a,b~0]", "EBF[a,b~0]"}};

    if (run_settings["method"] != ^"terms.fubar.methods.VB0") {
        samples = run_settings["samples"];
        chains  = run_settings["chains"];

        results.samples = {samples,grid_points};
        per_chain      = samples $ chains;

        from = 0;
        to   = per_chain;

        positive_ks                 = {};
        invariant_ks                = {};
        alpha_zero_ks               = {};
        beta_zero_ks                = {};
        proximal_ks                 = {};
        posterior_mean_alpha        = {};
        posterior_mean_beta         = {};
        denominators                = {};

        for (chain_id = 0; chain_id < chains; chain_id += 1) {
            io.ReportProgressBar                  ("PROCESSING", "Samples from chain " + (chain_id + 1));

            grid_samples = ((cache[utility.getGlobalValue("terms.fubar.cache.mcmc")])[chain_id])["weights"];
            P_ks = grid_samples * (cache[utility.getGlobalValue("terms.fubar.cache.conditionals")])["conditionals"];
            denominators[chain_id] = P_ks;

            positive_ks [chain_id] = grid_samples *
                          (positive_selection_stencil $ (cache[utility.getGlobalValue("terms.fubar.cache.conditionals")])["conditionals"]) / P_ks;
            
            invariant_ks [chain_id] = grid_samples *
                          (invariant_stencil $ (cache[utility.getGlobalValue("terms.fubar.cache.conditionals")])["conditionals"]) / P_ks;

            alpha_zero_ks [chain_id] = grid_samples *
                          (alpha_zero_stencil $ (cache[utility.getGlobalValue("terms.fubar.cache.conditionals")])["conditionals"]) / P_ks;

            beta_zero_ks [chain_id] = grid_samples *
                          (beta_zero_stencil $ (cache[utility.getGlobalValue("terms.fubar.cache.conditionals")])["conditionals"]) / P_ks;

            proximal_ks [chain_id] = grid_samples *
                          (proximal_stencil $ (cache[utility.getGlobalValue("terms.fubar.cache.conditionals")])["conditionals"]) / P_ks;

            posterior_mean_alpha [chain_id]     =            (grid_samples * diag_alpha *
                                                          (cache[utility.getGlobalValue("terms.fubar.cache.conditionals")])["conditionals"]) / P_ks;
            posterior_mean_beta [chain_id]      =            (grid_samples * diag_beta *
                                                          (cache[utility.getGlobalValue("terms.fubar.cache.conditionals")])["conditionals"]) / P_ks;

            for (i = from; i < to; i += 1) {
                for (r = 0; r < grid_points; r += 1) {
                    results.samples [i][r] = grid_samples[i-from][r];
                }
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
        
        P_ks = posterior_mean_over_grid_T * (cache[utility.getGlobalValue("terms.fubar.cache.conditionals")])["conditionals"];

        posterior_mean_alpha     =             (posterior_mean_over_grid_T * diag_alpha *
                                                          (cache[utility.getGlobalValue("terms.fubar.cache.conditionals")])["conditionals"]) / P_ks;
        posterior_mean_beta       =            (posterior_mean_over_grid_T * diag_beta *
                                                          (cache[utility.getGlobalValue("terms.fubar.cache.conditionals")])["conditionals"]) / P_ks;

        positive_ks  = posterior_mean_over_grid_T *
                           (positive_selection_stencil $ (cache[utility.getGlobalValue("terms.fubar.cache.conditionals")])["conditionals"]) / P_ks;

        invariant_ks = posterior_mean_over_grid_T *
                           (invariant_stencil $ (cache[utility.getGlobalValue("terms.fubar.cache.conditionals")])["conditionals"]) / P_ks;

        alpha_zero_ks = posterior_mean_over_grid_T *
                           (alpha_zero_stencil $ (cache[utility.getGlobalValue("terms.fubar.cache.conditionals")])["conditionals"]) / P_ks;

                 beta_zero_ks = posterior_mean_over_grid_T *
                                   (beta_zero_stencil $ (cache[utility.getGlobalValue("terms.fubar.cache.conditionals")])["conditionals"]) / P_ks;
        
                 proximal_ks = posterior_mean_over_grid_T *
                                   (proximal_stencil $ (cache[utility.getGlobalValue("terms.fubar.cache.conditionals")])["conditionals"]) / P_ks;
           }

    /* Compute Priors for EBF */
    _grid = cache["grid"];
    _inv_vec  = {grid_points, 1} ["_grid[_MATRIX_ELEMENT_ROW_][0]==0 && _grid[_MATRIX_ELEMENT_ROW_][1]==0"];
    _a0_vec   = {grid_points, 1} ["_grid[_MATRIX_ELEMENT_ROW_][0]==0"];
    _b0_vec   = {grid_points, 1} ["_grid[_MATRIX_ELEMENT_ROW_][1]==0"];
    _prox_vec = {grid_points, 1} ["_grid[_MATRIX_ELEMENT_ROW_][0]^2 + _grid[_MATRIX_ELEMENT_ROW_][1]^2 < proximal_radius^2"];

    priors = {
        "inv"  : +(posterior_mean_over_grid $ _inv_vec),
        "a0"   : +(posterior_mean_over_grid $ _a0_vec),
        "b0"   : +(posterior_mean_over_grid $ _b0_vec),
        "prox" : +(posterior_mean_over_grid $ _prox_vec)
    };

    site_results = {};
    i = 0;
    s = 0;

    report.posteriors  = {};
    report.sites_found = 0;
    chain_iterator = utility.Range (chains, 0, 1);

    for (partition_index = 0; partition_index < partition_count; partition_index += 1) {
        filter_info = (filter_specification [partition_index])[utility.getGlobalValue ("terms.data.coverage")];
        sites_in_partition = utility.Array1D (filter_info);

        partition_results    = {sites_in_partition, 14};
        partition_posteriors = {};

        for (s = 0; s < sites_in_partition; s += 1) {

             pp = posterior_mean_over_grid $ ((cache[utility.getGlobalValue("terms.fubar.cache.conditionals")])["conditionals"])[-1][i];
             partition_posteriors [s] = Transpose (pp * (1/(+pp)));
             if (run_settings["method"] != utility.getGlobalValue ("terms.fubar.methods.VB0")) {
                partition_results[s][0] = ComputeRandNeff (
                    utility.Map (chain_iterator, "_value_", "((`&posterior_mean_alpha`)[_value_])[-1][`&i`]")
                )[0];
                partition_results[s][1] = ComputeRandNeff (
                    utility.Map (chain_iterator, "_value_", "((`&posterior_mean_beta`)[_value_])[-1][`&i`]")
                )[0];
                
                partition_results[s][2] = ComputeRandNeff (
                    utility.Map (chain_iterator, "_value_", "((`&invariant_ks`)[_value_])[-1][`&i`]")
                )[0];
                
                partition_results[s][3] = ComputeRandNeff (
                    utility.Map (chain_iterator, "_value_", "((`&alpha_zero_ks`)[_value_])[-1][`&i`]")
                )[0];
                
                partition_results[s][4] = ComputeRandNeff (
                    utility.Map (chain_iterator, "_value_", "((`&beta_zero_ks`)[_value_])[-1][`&i`]")
                )[0];

                partition_results[s][5] = ComputeRandNeff (
                    utility.Map (chain_iterator, "_value_", "((`&proximal_ks`)[_value_])[-1][`&i`]")
                )[0];

                pos_selection = ComputeRandNeff (
                    utility.Map (chain_iterator, "_value_", "((`&positive_ks`)[_value_])[-1][`&i`]")
                );
                partition_results[s][6] = pos_selection[0];

                if (run_settings["method"] == utility.getGlobalValue ("terms.fubar.methods.MH")) {
                    partition_results[s][7] = pos_selection[1];
                    partition_results[s][8] = pos_selection[2];
                }

            } else {
                partition_results [s][0] = posterior_mean_alpha[i];
                partition_results [s][1] = posterior_mean_beta[i];
                partition_results [s][2] = invariant_ks[i];
                partition_results [s][3] = alpha_zero_ks[i];
                partition_results [s][4] = beta_zero_ks[i];
                partition_results [s][5] = proximal_ks[i];
                partition_results [s][6] = positive_ks[i];
            }

            partition_results[s][9]  = compute_ebf(partition_results[s][2], priors["inv"]);
            partition_results[s][10] = compute_ebf(partition_results[s][3], priors["a0"]);
            partition_results[s][11] = compute_ebf(partition_results[s][4], priors["b0"]);
            partition_results[s][12] = compute_ebf(partition_results[s][5], priors["prox"]);

            if (partition_results[s][12] >= run_settings["ebf"]) {
                if (report.sites_found == 0) {
                     fprintf (stdout, io.FormatTableRow (table_screen_output,table_output_options));
                     table_output_options[^"terms.table_options.header"] = FALSE;
                }
                
                report.row = {{"" + (1+filter_info[s]),
                                            partition_index + 1,
                                            Format(partition_results[s][0],10,3),
                                            Format(partition_results[s][1],10,3),
                                            Format(partition_results[s][5],6,4),
                                            Format(partition_results[s][12],8,2)}};
                                            
                fprintf (stdout, io.FormatTableRow (report.row,table_output_options));
                report.sites_found += 1;
            }

            i+=1;
        }
        site_results [partition_index] = partition_results;
        report.posteriors [partition_index] = partition_posteriors;

    }
    

    if (report.sites_found == 0) {
        console.log ("----\n## B-STILL inferred no sites under proximal constraint at EBF >= " + run_settings["ebf"]);
    } else {
        console.log ("----\n## B-STILL inferred " + report.sites_found + " sites under proximal constraint at EBF >= " + run_settings["ebf"]);
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

GetString (_hpv,HYPHY_VERSION,0);
fubar.json[terms.json.runtime] = _hpv;

io.SpoolJSON (fubar.json, fubar.codon_data_info[terms.json.json]);

return fubar.json;

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

function     fubar.RunPrompts (prompts) {

     if (prompts["grid"]) {
        KeywordArgument ("grid", "The number of grid points", "20");
        fubar.run_settings["grid size"] = io.PromptUser ("> Number of grid points per dimension (total number is D^2)",fubar.run_settings["grid size"],5,50,TRUE);
        prompts["grid"] = FALSE;
    }

    if (prompts["method"]) {
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

    if (prompts["ebf"]) {
        KeywordArgument ("ebf", "EBF threshold for reporting proximal invariance", fubar.run_settings["ebf"]);
        fubar.run_settings["ebf"] = io.PromptUser ("> EBF threshold for reporting proximal invariance", fubar.run_settings["ebf"], 0, 1e10, FALSE);
        prompts["ebf"] = FALSE;
    }

    if (prompts["radius-threshold"]) {
        KeywordArgument ("radius-threshold", "Expected substitution multiplier for defining 'near-zero' Selective Regime", fubar.run_settings["radius-threshold"]);
        fubar.run_settings["radius-threshold"] = io.PromptUser ("> Expected substitution multiplier for defining 'near-zero' Selective Regime", fubar.run_settings["radius-threshold"], 0, 10, FALSE);
        prompts["radius-threshold"] = FALSE;
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

lfunction fubar.scalers.Constrain (tree_name, node_name, model_description, ignore) {
    parameters.SetProprtionalConstraint (tree_name + "." + node_name + "." + (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.parameters.synonymous_rate")],
                              (model_description[utility.getGlobalValue ("terms.global")])[utility.getGlobalValue ("terms.fubar.branch_length_scaler.alpha")]);

    parameters.SetProprtionalConstraint (tree_name + "." + node_name + "." + (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.parameters.nonsynonymous_rate")],
                              (model_description[utility.getGlobalValue ("terms.global")])[utility.getGlobalValue ("terms.fubar.branch_length_scaler.beta")]);

}

//------------------------------------------------------------------------------

lfunction fubar.scalers.Unconstrain (tree_name, node_name, model_description, ignore) {
    parameters.RemoveConstraint (tree_name + "." + node_name + "." + (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.parameters.synonymous_rate")]);
    parameters.RemoveConstraint (tree_name + "." + node_name + "." + (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.parameters.nonsynonymous_rate")]);
}

//------------------------------------------------------------------------------------------------//

lfunction fubar.DefineAlphaBetaGrid (one_d_points, non_zero) {

    one_d_points    = Max (one_d_points, 5);
    alphaBetaGrid = {one_d_points^2,2};
    oneDGrid      = {one_d_points,1};

    neg_sel         = 0.7;
    neg_sel_points  = ((one_d_points)*neg_sel+0.5)$1;
    pos_sel_points  = one_d_points - neg_sel_points;

    // Modified for denser grid around zero using x^2
    for (_k = 0; _k < neg_sel_points; _k += 1) {
        oneDGrid [_k][0] =  (_k / (neg_sel_points - 1))^2;
    }

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
