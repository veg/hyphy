RequireVersion("2.31");

// namespace 'utility' for convenience functions
LoadFunctionLibrary("libv3/UtilityFunctions.bf");

// namespace 'alignments' for alignment-related functions
LoadFunctionLibrary("libv3/tasks/alignments.bf");

// namespace 'trees' for trees-related functions
LoadFunctionLibrary("libv3/tasks/trees.bf");

// namespace 'io' for interactive/datamonkey i/o functions
LoadFunctionLibrary("libv3/IOFunctions.bf");

// namespace 'estimators' for various estimator related functions
LoadFunctionLibrary("libv3/tasks/estimators.bf");

// namespace 'ancestral' for ancestral reconstruction functions
LoadFunctionLibrary("libv3/tasks/ancestral.bf");

// namespace 'stats' for various descriptive stats functions
LoadFunctionLibrary("libv3/stats.bf");

LoadFunctionLibrary("libv3/models/protein/empirical.bf");
LoadFunctionLibrary("libv3/models/protein/REV.bf");

LoadFunctionLibrary("libv3/terms-json.bf");
LoadFunctionLibrary("libv3/tasks/mpi.bf");
LoadFunctionLibrary("libv3/convenience/math.bf");

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

// TO DO: Add option for user to select baseline model for initial branch length fit.
// TO DO: Incorporate rate heterogeneity for branch length fits.



SetDialogPrompt ("Supply a list of files to include in the analysis (one per line)");
fscanf (PROMPT_FOR_FILE, "Lines", protein_gtr.file_list);


// To store a prior run, and load as needed
protein_gtr.cache_file = utility.getGlobalValue("LAST_FILE_PATH") + ".cache";
protein_gtr.analysis_results = io.LoadCacheFromFile (protein_gtr.cache_file); // either loaded or an empty dictionary
protein_gtr.file_list = io.validate_a_list_of_files (protein_gtr.file_list);
protein_gtr.file_list_count = Abs (protein_gtr.file_list);

protein_gtr.queue = mpi.CreateQueue ({"Headers" : utility.GetListOfLoadedModules ()});

io.ReportProgressMessageMD ("Protein GTR Fitter", "Initial branch length fit", "Initial branch length fit");

protein_gtr.fit_phase = 0;
protein_gtr.scores = {};
protein_gtr.phase_key = "WAG-Phase0";

/*************************** STEP ONE ***************************
Perform an initial fit of WAG to the data (or load cached fit.) 
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

        mpi.QueueJob (protein_gtr.queue, "protein_gtr.fitWAGtoFile", {"0" : protein_gtr.file_list[file_index]},
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
    // input arguments : current_results, previous_values, phase, final. Pretty sure "previous values" are initial values for fit.
    current_results = utility.Map (utility.Filter (protein_gtr.analysis_results, "_value_", "_value_/'WAG-Phase0'"), "_value_", "_value_['WAG-Phase0']");

    protein_gtr.current_gtr_fit = protein_gtr.fitGTRtoFileList (current_results, protein_gtr.current_gtr_fit, result_key, FALSE);
}

// Record sum(logL) from this phase
protein_gtr.scores + (protein_gtr.analysis_results[result_key])["LogL"];

fprintf(stdout, protein_gtr.analysis_results);
return 0;



/* Now that the initial GTR fit has been performed, we toggle between a GTR fit and a branch length fit under the estimated GTR parameters */
protein_gtr.shared_EFV = (utility.Values (protein_gtr.current_gtr_fit [terms.efv_estimate]))[0];
if (Type (protein_gtr.shared_EFV) == "String") {
    protein_gtr.shared_EFV = Eval (protein_gtr.shared_EFV);
}


// TO DO: Change stopping criterion to matrix convergence based on some tolerance, rather than LogL convergence
TOLERANCE=0.1;
for (;;) {
    protein_gtr.phase_results = protein_gtr.run_gtr_iteration_branch_lengths();
    protein_gtr.scores + protein_gtr.phase_results["LogL"];

    result_key = "REV-Phase-" + protein_gtr.fit_phase;

    if (utility.Has (protein_gtr.analysis_results, result_key, None)) {
        io.ReportProgressMessageMD ("Protein GTR Fitter", result_key,
                                    "Loaded cached results for '" + result_key + "'. Log(L) = " + (protein_gtr.analysis_results[result_key])["LogL"] );
        protein_gtr.current_gtr_fit = protein_gtr.analysis_results [result_key];
    } else {
        protein_gtr.current_gtr_fit = protein_gtr.fitGTRtoFileList (utility.Map (utility.Filter (protein_gtr.analysis_results, "_value_", "_value_/protein_gtr.phase_results['phase']"), "_value_", "_value_[protein_gtr.phase_results['phase']]"), protein_gtr.current_gtr_fit, result_key, FALSE);
    }
    protein_gtr.scores + (protein_gtr.analysis_results[result_key])["LogL"];
    if (protein_gtr.scores[Abs(protein_gtr.scores)-1] - protein_gtr.scores[Abs(protein_gtr.scores)-2] < TOLERANCE) {
        break;
    }
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


