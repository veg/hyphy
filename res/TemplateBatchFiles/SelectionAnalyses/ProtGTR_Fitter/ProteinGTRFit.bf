RequireVersion("2.31");

LoadFunctionLibrary("libv3/UtilityFunctions.bf"); // namespace 'utility' for convenience functions
LoadFunctionLibrary("libv3/tasks/alignments.bf"); // namespace 'alignments' for alignment-related functions
LoadFunctionLibrary("libv3/tasks/trees.bf"); // namespace 'trees' for trees-related functions
LoadFunctionLibrary("libv3/IOFunctions.bf"); // namespace 'io' for interactive/datamonkey i/o functions
LoadFunctionLibrary("libv3/tasks/estimators.bf"); // namespace 'estimators' for various estimator related functions
LoadFunctionLibrary("libv3/tasks/ancestral.bf"); // namespace 'ancestral' for ancestral reconstruction functions
LoadFunctionLibrary("libv3/stats.bf"); // namespace 'stats' for various descriptive stats functions
LoadFunctionLibrary("libv3/terms-json.bf");
LoadFunctionLibrary("libv3/tasks/mpi.bf");
LoadFunctionLibrary("libv3/convenience/math.bf");

LoadFunctionLibrary("libv3/models/rate_variation.bf"); // rate variation

// Protein models
LoadFunctionLibrary("libv3/models/protein/empirical.bf");
LoadFunctionLibrary("libv3/models/protein/REV.bf");

LoadFunctionLibrary("ProteinGTRFit_helper.ibf"); // Moved all functions from this file into loaded file, for clarity.
/*------------------------------------------------------------------------------*/

utility.ToggleEnvVariable ("NORMALIZE_SEQUENCE_NAMES", 1);

io.DisplayAnalysisBanner({
    "info": "Fit a general time reversible model to a collection
    of training protein sequence alignments. Generate substitution
    and scoring matrices following the procedures described in Nickle et al 2007",
    "version": "0.01",
    "reference": "Nickle DC, Heath L, Jensen MA, Gilbert PB, Mullins JI, Kosakovsky Pond SL (2007) HIV-Specific Probabilistic Models of Protein Evolution. PLoS ONE 2(6): e503. doi:10.1371/journal.pone.0000503",
    "authors": "Sergei L Kosakovsky Pond and Stephanie J Spielman",
    "contact": "{spond,stephanie.spielman}@temple.edu"
});



OPTIMIZATION_PRECISION = 0.1; // for testing purposes.
protein_gtr.tolerance=0.01; // should be different if logL or RMSE

// Load file containing paths to alignments for fitting and assess whether to start from scratch or resume a cached analysis
SetDialogPrompt ("Supply a list of files to include in the analysis (one per line)");
fscanf (PROMPT_FOR_FILE, "Lines", protein_gtr.file_list);
protein_gtr.cache_file = utility.getGlobalValue("LAST_FILE_PATH") + ".cache";
protein_gtr.file_list = io.validate_a_list_of_files (protein_gtr.file_list);
protein_gtr.file_list_count = Abs (protein_gtr.file_list);

protein_gtr.analysis_results = {}; //Overwrite with cache below as needed.
if (protein_gtr.cache_file) { 
    protein_gtr.use_cache = io.SelectAnOption ({{"YES", "Resume analysis using the detected cache file."}, {"NO", "Launch a new analysis and *overwrite* the detected cache file."}}, "A cache file of a prior analysis on this list of files was detected. Would you like to use it?");
    if (protein_gtr.use_cache == "YES"){
        protein_gtr.analysis_results = io.LoadCacheFromFile (protein_gtr.cache_file);
     }
        
}
console.log(protein_gtr.analysis_results);

return 0;
//protein_gtr.use_cache = io.SelectAnOption ({{"YES", "Resume analysis using the detected cache file."}, {"NO", "Launch a new analysis and *overwrite* the detected cache file."}}, "A cache file of a prior analysis on this list of files was detected. Would you like to use it?");
/*


if (protein_gtr.cache == "YES"){
    
    SetDialogPrompt ("Supply a list of files to include in the analysis (one per line)");
    fscanf (PROMPT_FOR_FILE, "Lines", protein_gtr.file_list);
    protein_gtr.file_list = io.validate_a_list_of_files (protein_gtr.file_list);
    protein_gtr.file_list_count = Abs (protein_gtr.file_list);
    protein_gtr.analysis_results = {};


else{

    SetDialogPrompt ("Supply the full path to your cached analysis.");
    fscanf (PROMPT_FOR_FILE, "String", protein_gtr.cache_file);
    protein_gtr.file_list = io.validate_a_list_of_files (protein_gtr.file_list);
    protein_gtr.file_list_count = Abs (protein_gtr.file_list);
    protein_gtr.analysis_results = io.LoadCacheFromFile (protein_gtr.cache_file);



protein_gtr.analysis_results = io.LoadCacheFromFile (protein_gtr.cache_file); // either loaded or an empty dictionary

*/

// TO DO: Add option for user to select baseline model for initial branch length fit. This may need to reference a new dictionary of models that will be created in protein.bf.
// TO DO: Related to above TO DO, WAG is currently hardcoded in various places throughout this and the helper batchfile, in particular as the Phase name. Bad news. Either have a generic name ("empirical" or something) or populate it with the chosen model string (JTT, WAG, LG, etc.)
// TO DO: Add option for user to select convergence criterion.
// TO DO: Prompt users if they would like to load the cache file or not


protein_gtr.baseline_model  = io.SelectAnOption (models.protein.empirical_options,
                                                 "Select an empirical protein model to use for optimizing the provided branch lengths:");
/*
if (protein_gtr.baseline_model == "WAG"){
    
}

console.log(protein_gtr.baseline_model);
protein_gtr.convergence_type  = io.SelectAnOption (models.protein.empirical_options,
                                                 "Select an empirical protein model to use for optimizing the provided branch lengths:");
*/




// To store a prior run, and load as needed













return 0;



OPTIMIZATION_PRECISION = 0.1; // for testing purposes.
protein_gtr.tolerance=0.01; 



// To store a prior run, and load as needed
protein_gtr.cache_file = utility.getGlobalValue("LAST_FILE_PATH") + ".cache";
protein_gtr.analysis_results = io.LoadCacheFromFile (protein_gtr.cache_file); // either loaded or an empty dictionary


protein_gtr.queue = mpi.CreateQueue ({"Headers" : utility.GetListOfLoadedModules ()});

io.ReportProgressMessageMD ("Protein GTR Fitter", "Initial branch length fit", "Initial branch length fit");

protein_gtr.fit_phase = 0;
protein_gtr.scores = {};
protein_gtr.phase_key = "WAG-Phase0";

/*************************** STEP ONE ***************************
Perform an initial fit of WAG+4G to the data (or load cached fit.) 
*****************************************************************/
for (file_index = 0; file_index < protein_gtr.file_list_count; file_index += 1) {

    cached_file = protein_gtr.file_list[file_index] + "' " + (file_index+1) + "/" + protein_gtr.file_list_count;

    // This function checks if a key (2nd argument) is found in a given dictionary/matrix (1st arg). 3rd optional argument will check type of key, if specified.
    if (utility.Has (protein_gtr.analysis_results, {{ protein_gtr.file_list[file_index], protein_gtr.phase_key}}, None)) {

        logL = ((protein_gtr.analysis_results[protein_gtr.file_list[file_index]])[protein_gtr.phase_key])["LogL"];
        io.ReportProgressMessageMD ("Protein GTR Fitter", " * Initial branch length fit",
                                    "Loaded cached results for '" + cached_file + ". Log(L) = " + logL);

    } else {
        io.ReportProgressMessageMD ("Protein GTR Fitter", " * Initial branch length fit",
                                    "Dispatching file '" + cached_file);
        // Constant
        /*
        mpi.QueueJob (protein_gtr.queue, "protein_gtr.fitWAGtoFile", {"0" : protein_gtr.file_list[file_index]},
                                                            "protein_gtr.handle_wag_callback");
        */
        
        // Four category Gamma
         mpi.QueueJob (protein_gtr.queue, "protein_gtr.fitWAGwithGammatoFile", {"0" : protein_gtr.file_list[file_index]},
                                                            "protein_gtr.handle_wag_callback");
       
    }
}
mpi.QueueComplete (protein_gtr.queue);

// Sum of the logL from fitted WAG model across each data set
protein_gtr.run_gtr_iteration_branch_lengths.logL = math.Sum (utility.Map (utility.Filter (protein_gtr.analysis_results, "_value_", "_value_/protein_gtr.phase_key"), "_value_", "(_value_[protein_gtr.phase_key])['LogL']"));

io.ReportProgressMessageMD ("Protein GTR Fitter", " * Initial branch length fit",
                            "Overall Log(L) = " + protein_gtr.run_gtr_iteration_branch_lengths.logL);


// Adding to the dictionary adds a new key:value. Key starts from 0 by default, and value here will be the summed logL from initial WAG fits across datasets
protein_gtr.scores + protein_gtr.run_gtr_iteration_branch_lengths.logL;



return 0;

/*************************** STEP TWO ***************************
          Perform an initial GTR fit on the data
*****************************************************************/
fprintf(stdout, "============ Phase 2 ==========");
// Check cache for the REV phase and load, or fit GTR. 
result_key = "REV-Phase-" + protein_gtr.fit_phase;
if (utility.Has (protein_gtr.analysis_results, result_key, None)) {
    io.ReportProgressMessageMD ("Protein GTR Fitter", result_key,
                                "Loaded cached results for '" + result_key + "'. Log(L) = " + (protein_gtr.analysis_results[result_key])["LogL"] );
    protein_gtr.current_gtr_fit = protein_gtr.analysis_results [result_key];
} else {
    
    current_results = utility.Map (utility.Filter (protein_gtr.analysis_results, "_value_", "_value_/'WAG-Phase0'"), "_value_", "_value_['WAG-Phase0']"); // input arguments : current_results, previous_values, phase, final. Pretty sure "previous values" are initial values for fit.

    protein_gtr.current_gtr_fit = protein_gtr.fitGTRtoFileList (current_results, protein_gtr.current_gtr_fit, result_key, FALSE);
}

// Record sum(logL) from this phase
protein_gtr.scores + (protein_gtr.analysis_results[result_key])["LogL"];



/* Now that the initial GTR fit has been performed, we toggle between a GTR fit and a branch length fit under the estimated GTR parameters */
protein_gtr.shared_EFV = (utility.Values (protein_gtr.current_gtr_fit [terms.efv_estimate]))[0];
if (Type (protein_gtr.shared_EFV) == "String") {
    protein_gtr.shared_EFV = Eval (protein_gtr.shared_EFV);
}


for (;;) {
    //console.log("sans bl score");
    //console.log(protein_gtr.scores);
    
    // Optimize branch lengths, with 4-category gamma. Record logL
    protein_gtr.phase_results = protein_gtr.run_gtr_iteration_branch_lengths(); 

    protein_gtr.scores + protein_gtr.phase_results["LogL"]; // WAIT, are too many logL's being saved? We should only be saving these between Q inferences, not between branch length optimizations.
    
    //console.log("with bl score");
    //console.log(protein_gtr.scores);

    result_key = "REV-Phase-" + protein_gtr.fit_phase;

    if (utility.Has (protein_gtr.analysis_results, result_key, None)) {
        io.ReportProgressMessageMD ("Protein GTR Fitter", result_key,
                                    "Loaded cached results for '" + result_key + "'. Log(L) = " + (protein_gtr.analysis_results[result_key])["LogL"] );
        protein_gtr.current_gtr_fit = protein_gtr.analysis_results [result_key];
    } else {
        protein_gtr.current_gtr_fit = protein_gtr.fitGTRtoFileList (utility.Map (utility.Filter (protein_gtr.analysis_results, "_value_", "_value_/protein_gtr.phase_results['phase']"), "_value_", "_value_[protein_gtr.phase_results['phase']]"), protein_gtr.current_gtr_fit, result_key, FALSE);
    }
    protein_gtr.scores + (protein_gtr.analysis_results[result_key])["LogL"];

    console.log("with rev score");
    console.log(protein_gtr.scores);
    console.log("convergence 1");
    console.log( protein_gtr.scores[Abs(protein_gtr.scores)-1] );
    console.log("convergence 2");
    console.log( protein_gtr.scores[Abs(protein_gtr.scores)-2] );
    return 0;

    // Convergence check
    
    /*
    // Here, convergence is assessed via logL.
    //-------------------------------------------------------------------------------------------------------------------------------------------------------//
    // TODO: IS THIS A BUG? SHOULD BE -1, -3 (not -1, -2)
    if (protein_gtr.scores[Abs(protein_gtr.scores)-1] - protein_gtr.scores[Abs(protein_gtr.scores)-2] < protein_gtr.tolerance) {
        break;
    } 
    //-------------------------------------------------------------------------------------------------------------------------------------------------------//
    */
    
    // Here, convergence is assessed via RMSE
    //-------------------------------------------------------------------------------------------------------------------------------------------------------//

    previous_Q = (protein_gtr.analysis_results["REV-Phase-" + (protein_gtr.fit_phase-1)])["global"]; // isolate Q from previous phase
    current_Q = (protein_gtr.analysis_results[result_key])["global"];                                // isolate Q from current phase

    // Calculate RMSE between previous, current fitted Q's
    rmse = 0;
    N = 0;    
    for (l1 = 0; l1 < 20; l1 += 1) {
        for (l2 = l1 + 1; l2 < 20; l2 += 1) {
                
            previous = (previous_Q[ terms.aminoacidRate (models.protein.alphabet[l1],models.protein.alphabet[l2]) ])[terms.MLE];
            current = (current_Q[ terms.aminoacidRate (models.protein.alphabet[l1],models.protein.alphabet[l2]) ])[terms.MLE];
        
            rmse += (previous - current)^2;      
            N += 1;
         }
    }
    rmse = Sqrt( rmse/N );
    if (rmse < protein_gtr.tolerance) {
        break;
    }
    //-------------------------------------------------------------------------------------------------------------------------------------------------------//
    
    
}





result_key = "REV-Final";

if (utility.Has (protein_gtr.analysis_results, result_key, None)) {
    io.ReportProgressMessageMD ("Protein GTR Fitter", result_key,
                                "Loaded cached results for '" + result_key + "'. Log(L) = " + (protein_gtr.analysis_results[result_key])["LogL"] );
    protein_gtr.current_gtr_fit = protein_gtr.analysis_results [result_key];
} else {
    protein_gtr.current_gtr_fit = protein_gtr.fitGTRtoFileList (utility.Map (utility.Filter (protein_gtr.analysis_results, "_value_", "_value_/protein_gtr.phase_results['phase']"), "_value_", "_value_[protein_gtr.phase_results['phase']]"), protein_gtr.current_gtr_fit, result_key, TRUE);
}


// do a final tune-up by reoptimizing everything


