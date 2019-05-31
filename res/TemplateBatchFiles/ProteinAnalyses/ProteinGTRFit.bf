RequireVersion("2.3.3");
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/stats.bf");
LoadFunctionLibrary("libv3/all-terms.bf");

LoadFunctionLibrary("libv3/tasks/ancestral.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/tasks/mpi.bf");
LoadFunctionLibrary("libv3/convenience/math.bf");

LoadFunctionLibrary("libv3/models/rate_variation.bf");

LoadFunctionLibrary("libv3/models/protein/empirical.bf");
LoadFunctionLibrary("libv3/models/protein/REV.bf");
LoadFunctionLibrary("plusF_helper.ibf");
LoadFunctionLibrary("ProteinGTRFit_helper.ibf");


/*------------------------------------------------------------------------------*/

utility.ToggleEnvVariable ("NORMALIZE_SEQUENCE_NAMES", 1);
utility.SetEnvVariable ("OPTIMIZATION_PRECISION", 0.1); // for testing purposes.
utility.ToggleEnvVariable ("PRODUCE_OPTIMIZATION_LOG", 1); // for testing purposes

io.DisplayAnalysisBanner({
    terms.io.info: "Fit a general time reversible model to a collection
    of training protein sequence alignments. Generate substitution
    and scoring matrices following the procedures described in Nickle et al 2007",
    terms.io.version: "0.01",
    terms.io.reference: "Nickle DC, Heath L, Jensen MA, Gilbert PB, Mullins JI, Kosakovsky Pond SL (2007) HIV-Specific Probabilistic Models of Protein Evolution. PLoS ONE 2(6): e503. doi:10.1371/journal.pone.0000503",
    terms.io.authors: "Sergei L Kosakovsky Pond and Stephanie J Spielman",
    terms.io.contact: "{spond,stephanie.spielman}@temple.edu"
});


protein_gtr.filename_to_index = terms.data.filename_to_index;
protein_gtr.logl = terms.fit.log_likelihood;
protein_gtr.phase = terms.fit.phase;
protein_gtr.json.information = "information";
protein_gtr.baseline_phase = "Baseline Phase";

protein_gtr.options.convergence_type = "convergence type";
protein_gtr.options.tolerance        = "tolerance";
protein_gtr.options.baseline_model   = "baseline model";
protein_gtr.options.rate_variation   = "use rate variation";

/********************************************** MENU PROMPTS ********************************************************/
/********************************************************************************************************************/

// Load file containing paths to alignments for fitting and assess whether to start from scratch or resume a cached analysis

SetDialogPrompt ("Supply a list of files to include in the analysis (one per line)");
fscanf (PROMPT_FOR_FILE, "Lines", protein_gtr.file_list);
protein_gtr.listfile = utility.getGlobalValue("LAST_FILE_PATH");
protein_gtr.cache_file = protein_gtr.listfile  + ".cache";
protein_gtr.json_file = protein_gtr.listfile  + ".json";
protein_gtr.file_list = io.validate_a_list_of_files (protein_gtr.file_list);
protein_gtr.file_list_count = Abs (protein_gtr.file_list);

// Populate analysis_results and important variables from cache, or prompt for all variables.
if (io.FileExists(protein_gtr.cache_file)) {
    protein_gtr.use_cache = io.SelectAnOption ({{"YES", "Resume analysis using the detected cache file."}, {"NO", "Launch a new analysis and *overwrite* the detected cache file."}}, "A cache file of a prior analysis on this list of files was detected. Would you like to use it?");
}

// Load all information from cache
if (protein_gtr.use_cache == "YES" || protein_gtr.use_cache == 1){
    protein_gtr.analysis_results = io.LoadCacheFromFile (protein_gtr.cache_file);
    protein_gtr.load_cached_options();  // Load previously selected options into namespace

}
// Prompt for all information
else {

    protein_gtr.analysis_results = {};

    // Prompt for convergence assessment type
    protein_gtr.convergence_type = io.SelectAnOption( protein_gtr.convergence_options, "Select a convergence criterion.");

    if (protein_gtr.convergence_type == "LogL"){
        protein_gtr.tolerance = io.PromptUser ("\n>Provide a tolerance level for convergence assessment (Default 0.01)",0.01,0,1,FALSE); // default, lower, upper, is_integer
    }
    else {
        protein_gtr.tolerance = io.PromptUser ("\n>Provide a tolerance level for convergence assessment (Default 0.001)",0.01,0,1,FALSE); // default, lower, upper, is_integer
    }

    // Prompt for baseline AA model
    protein_gtr.baseline_model  = io.SelectAnOption (models.protein.empirical_models,
                                                     "Select an empirical protein model to use for optimizing the provided branch lengths (we recommend LG):");

    // Prompt for rate variation
    protein_gtr.use_rate_variation = io.SelectAnOption( protein_gtr.rate_variation_options, "Would you like to optimize branch lengths with rate variation (uses a four-category gamma)?");

    protein_gtr.save_options();

}



protein_gtr.baseline_Rij = plusF_helper.Rij_options[protein_gtr.baseline_model]; // Defined in helper
if (protein_gtr.use_rate_variation == "Yes"){
    protein_gtr.final_baseline_model = "protein_gtr.plusF.ModelDescription.withGamma";
    protein_gtr.rev_model_branch_lengths = "protein_gtr.REV.ModelDescription.withGamma";
} else{
    protein_gtr.final_baseline_model = "protein_gtr.plusF.ModelDescription";
    protein_gtr.rev_model_branch_lengths = "protein_gtr.REV.ModelDescription";
}
/********************************************************************************************************************/


protein_gtr.queue = mpi.CreateQueue ({  "Headers"   : utility.GetListOfLoadedModules ("libv3/") ,
                                        "Functions" :
                                        {
                                            {"protein_gtr.REV.ModelDescription",
                                             "protein_gtr.REV.ModelDescription.withGamma",
                                             "protein_gtr.REV.ModelDescription.freqs",
                                             "protein_gtr.plusF.ModelDescription",
                                             "protein_gtr.plusF._GenerateRate",
                                             "protein_gtr.plusF.frequencies"
                                            }
                                        },

                                        "Variables" : {{
                                            "protein_gtr.shared_EFV",
                                            "protein_gtr.final_baseline_model",
                                            "protein_gtr.rev_model_branch_lengths",
                                            "protein_gtr.baseline_Rij"

                                        }}
                                     });

io.ReportProgressMessageMD ("Protein GTR Fitter", "Initial branch length fit", "Initial branch length fit");

protein_gtr.fit_phase = 0;
protein_gtr.scores = {};
protein_gtr.phase_key = protein_gtr.baseline_phase;



/*************************** STEP ONE ***************************
Perform an initial fit of Baseline model+4G to the data (or load cached fit.)
*****************************************************************/
console.log("\n\n[PHASE 1] Performing initial branch length optimization using " + protein_gtr.baseline_model);

for (file_index = 0; file_index < protein_gtr.file_list_count; file_index += 1) {

    cached_file = protein_gtr.file_list[file_index] + "' " + (file_index+1) + "/" + protein_gtr.file_list_count;

    if (utility.Has (protein_gtr.analysis_results, {{ protein_gtr.file_list[file_index], protein_gtr.phase_key}}, None)) {

        thisKey = (protein_gtr.analysis_results[protein_gtr.file_list[file_index]])[protein_gtr.baseline_phase];
        io.ReportProgressMessageMD ("Protein GTR Fitter", " * Initial branch length fit",
                                    "Loaded cached results for '" + cached_file + ". Log(L) = " + thisKey[terms.fit.log_likelihood]);

    } else {
        io.ReportProgressMessageMD ("Protein GTR Fitter", " * Initial branch length fit",
                                    "Dispatching file '" + cached_file);


         mpi.QueueJob (protein_gtr.queue, "protein_gtr.fitBaselineToFile", {"0" : protein_gtr.file_list[file_index]},
                                                            "protein_gtr.handle_baseline_callback");

    }
}
mpi.QueueComplete (protein_gtr.queue);


// Sum of the logL from fitted baseline model across each data set
protein_gtr.baseline_fit_logL = math.Sum (utility.Map (utility.Filter (protein_gtr.analysis_results, "_value_", "_value_/protein_gtr.phase_key"), "_value_", "(_value_[protein_gtr.phase_key])[terms.fit.log_likelihood]"));

io.ReportProgressMessageMD ("Protein GTR Fitter", " * Initial branch length fit",
                            "Overall Log(L) = " + protein_gtr.baseline_fit_logL);


//io.WriteCacheToFile (^"protein_gtr.cache_file", ^"protein_gtr.analysis_results");


/*************************** STEP TWO ***************************
          Perform an initial GTR fit on the data
*****************************************************************/
console.log("\n\n[PHASE 2] Performing initial REV fit to the data");

result_key = "REV-Phase-" + protein_gtr.fit_phase;

if (utility.Has (protein_gtr.analysis_results, result_key, None)) {
    io.ReportProgressMessageMD ("Protein GTR Fitter", result_key,
                                "Loaded cached results for '" + result_key + "'. Log(L) = " + (protein_gtr.analysis_results[result_key])[terms.fit.log_likelihood] );
    protein_gtr.current_gtr_fit = protein_gtr.analysis_results [result_key];
} else {
    protein_gtr.current_gtr_fit = protein_gtr.fitGTRtoFileList (utility.Map (utility.Filter (protein_gtr.analysis_results, "_value_", "_value_/'" + protein_gtr.phase_key + "'"), "_value_", "_value_['" + protein_gtr.phase_key + "']"), protein_gtr.current_gtr_fit, result_key, FALSE);
}

// Record logL
protein_gtr.scores + (protein_gtr.analysis_results[result_key])[terms.fit.log_likelihood];

/* Now that the initial GTR fit has been performed, we toggle between a GTR fit and a branch length fit under the estimated GTR parameters */
protein_gtr.shared_EFV = (utility.UniqueValues (protein_gtr.current_gtr_fit [terms.efv_estimate]))[0];
if (Type (protein_gtr.shared_EFV) == "String") {
    protein_gtr.shared_EFV = Eval (protein_gtr.shared_EFV);
}


/********************************** STEP THREE ******************************************
  Iteratively optimize branch lengths with previous REV fit, and re-optimize REV
*****************************************************************************************/

console.log("\n\n[PHASE 3] Iteratively optimizing branch lengths and fitting REV model until convergence.");
for (;;) {


    // Optimize branch lengths, with 4-category gamma. Record logL
    protein_gtr.phase_results = protein_gtr.run_gtr_iteration_branch_lengths();

    // Commented out below because this is never actually used in the analysis, and it is always cached anyways
    // protein_gtr.scores + protein_gtr.phase_results[terms.fit.log_likelihood];

    result_key = "REV-Phase-" + protein_gtr.fit_phase;

    if (utility.Has (protein_gtr.analysis_results, result_key, None)) {
        io.ReportProgressMessageMD ("Protein GTR Fitter", result_key,
                                    "Loaded cached results for '" + result_key + "'. Log(L) = " + (protein_gtr.analysis_results[result_key])[terms.fit.log_likelihood] );
        protein_gtr.current_gtr_fit = protein_gtr.analysis_results [result_key];
    } else {
        protein_gtr.current_gtr_fit = protein_gtr.fitGTRtoFileList (utility.Map (utility.Filter (protein_gtr.analysis_results, "_value_", "_value_/protein_gtr.phase_results['phase']"), "_value_", "_value_[protein_gtr.phase_results['phase']]"), protein_gtr.current_gtr_fit, result_key, FALSE);
    }

    protein_gtr.scores + (protein_gtr.analysis_results[result_key])[terms.fit.log_likelihood];

    // LogL
    if (protein_gtr.convergence_type == "LogL"){
        console.log("\n\n[PHASE 3] Delta log-L = " + (protein_gtr.scores[Abs(protein_gtr.scores)-1] - protein_gtr.scores[Abs(protein_gtr.scores)-2]));
        if (protein_gtr.scores[Abs(protein_gtr.scores)-1] - protein_gtr.scores[Abs(protein_gtr.scores)-2] < protein_gtr.tolerance){
            break;
        }
    }
    // RMSE
    else {
        previous_Q = (protein_gtr.analysis_results["REV-Phase-" + (protein_gtr.fit_phase-1)])[terms.global]; // isolate Q from previous phase
        current_Q = (protein_gtr.analysis_results[result_key])[terms.global];                                // isolate Q from current phase

        // Calculate RMSE between previous, current fitted Q's
        rmse = 0;
        N = 0;
        for (l1 = 0; l1 < 20; l1 += 1) {
            for (l2 = l1 + 1; l2 < 20; l2 += 1) {

                previous = (previous_Q[ terms.aminoacidRate (models.protein.alphabet[l1],models.protein.alphabet[l2]) ])[terms.fit.MLE];
                current = (current_Q[ terms.aminoacidRate (models.protein.alphabet[l1],models.protein.alphabet[l2]) ])[terms.fit.MLE];

                rmse += (previous - current)^2;
                N += 1;
             }
        }
        rmse = Sqrt( rmse/N );
        console.log("\n\n[PHASE 3] RMSE = " + rmse);
        if (rmse < protein_gtr.tolerance) {
            break;
        }
    }

    // UNCLEAR IF TO DO: Before going to next iteration which begins w/ a branch length fit, we need to perform normalization. HyPhy might already do this.

}




/********************************** STEP FOUR ******************************************
  Perform a final optimization on the REV matrix while also optimizing the frequencies
***************************************************************************************/


console.log("\n\n[PHASE 4] Convergence achieved. Optimizing final model.");
// do a final tune-up by reoptimizing everything
result_key = "REV-Final";
if (utility.Has (protein_gtr.analysis_results, result_key, None)) {
    io.ReportProgressMessageMD ("Protein GTR Fitter", result_key,
                                "Loaded cached results for '" + result_key + "'. Log(L) = " + (protein_gtr.analysis_results[result_key])[terms.fit.log_likelihood] );
    protein_gtr.current_gtr_fit = protein_gtr.analysis_results [result_key];
} else {
    protein_gtr.current_gtr_fit = protein_gtr.fitGTRtoFileList (utility.Map (utility.Filter (protein_gtr.analysis_results, "_value_", "_value_/protein_gtr.phase_results['phase']"), "_value_", "_value_[protein_gtr.phase_results['phase']]"), protein_gtr.current_gtr_fit, result_key, TRUE);
}


/* Save the JSON */
loaded_cache = io.LoadCacheFromFile(protein_gtr.cache_file);
io.SpoolJSON(loaded_cache, protein_gtr.json_file);
