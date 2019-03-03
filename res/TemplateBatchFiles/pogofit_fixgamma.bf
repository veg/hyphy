RequireVersion("2.3.9");
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
LoadFunctionLibrary("libv3/models/protein.bf");
LoadFunctionLibrary("pogofit_helper_fixgamma.bf"); // Functions, model definitions used for this batchfile.


/*------------------------------------------------------------------------------*/

// do not uncomment this line; debugging purposes ONLY.
//utility.ToggleEnvVariable ("OPTIMIZATION_TIME_HARD_LIMIT", 2);

utility.ToggleEnvVariable ("NORMALIZE_SEQUENCE_NAMES", 1);

pogofit.analysis_banner = {
    terms.io.info: "PogoFit, *P*rotein *G*TR *Fit*ter: Fit a general time reversible (GTR) model to a collection of training protein sequence alignments.",
    terms.io.version: "0.01",
    terms.io.reference: "TBD",
    terms.io.authors: "Sergei L Kosakovsky Pond and Stephanie J Spielman",
    terms.io.contact: "spond@temple.edu; spielman@rowan.edu",
    terms.io.requirements: "All alignments must be in HyPhy-format: Each file must contain a protein multiple sequence alignment and newick phylogeny. NEXUS input is not accepted."
};
io.DisplayAnalysisBanner(pogofit.analysis_banner);

pogofit.baseline_phase   = "Baseline Fit";
pogofit.final_phase      = "GTR Fit";

pogofit.options.imputation          = "Impute zero rates";
pogofit.options.dataset_information = "Dataset information";
pogofit.options.number_of_datasets  = "Number of datasets";
pogofit.options.baseline_model      = "Baseline model";

pogofit.single   = "Single";
pogofit.multiple = "Multiple";

pogofit.impute    = "Yes";
pogofit.no_impute = "No";

pogofit.output_hyphy    = "HyPhy";
pogofit.output_paml     = "PAML";
pogofit.output_raxml    = "RAxML";
pogofit.output_all      = "All";
pogofit.hyphy_model_ext = ".POGOFIT.fitted_model";
pogofit.paml_model_ext  = ".POGOFIT.paml";
pogofit.raxml_model_ext = ".POGOFIT.raxml";

pogofit.analysis_results = {};
pogofit.analysis_results[terms.json.analysis] = pogofit.analysis_banner;

pogofit.timers = {};

/********************************************** MENU PROMPTS ********************************************************/
/********************************************************************************************************************/


// Prompt for number of files to analyze, and read file list accordingly //
pogofit.one_or_many  = io.SelectAnOption ({{pogofit.multiple, "Infer a protein model from multiple training datasets (this is more common)."}, 
                                               {pogofit.single, "Infer a protein model from a single training datasets."}}, 
                                                "How many datasets will be used to fit the protein model?");                                    
/****************** Read input data ***********************/
if (pogofit.one_or_many == pogofit.single)
{
    pogofit.input_file = io.PromptUserForString ("Provide the filename of the alignment to analyze");
    pogofit.file_list           = {{pogofit.input_file}};
}

if (pogofit.one_or_many == pogofit.multiple)
{
    SetDialogPrompt ("Supply a list of files to include in the analysis (one per line)");
    fscanf (PROMPT_FOR_FILE, "Lines", pogofit.file_list);
    pogofit.input_file  = utility.getGlobalValue("LAST_FILE_PATH");
}


pogofit.output_model_prefix = pogofit.input_file;
pogofit.json_file           = pogofit.input_file  + ".POGOFIT.json";
pogofit.file_list           = io.validate_a_list_of_files (pogofit.file_list);
pogofit.file_list_count     = Abs (pogofit.file_list);
pogofit.index_to_filename   = utility.SwapKeysAndValues(pogofit.file_list);
/*****************************************************************************/


// Prompt for baseline model, used to optimize branch lengths and as intial GTR matrix //
pogofit.baseline_model  = io.SelectAnOption (models.protein.empirical_models,
                                                "Select a protein model to use for the baseline:");
pogofit.baseline_model_name = pogofit.baseline_model + "+F, with 4 category Gamma rates";
pogofit.baseline_model_desc = "pogofit.Baseline.ModelDescription.withGamma";
pogofit.initial_rates       = Eval("models.protein." + pogofit.baseline_model + ".Rij");
pogofit.rev_model           = "models.protein.REV.ModelDescription.withGamma";
                   

// Prompt for zero-rate imputation //
pogofit.imputation  = io.SelectAnOption ({{pogofit.impute, "Impute zero rates as in Nickle et al. 2007 (Recommended)"},
                                          {pogofit.no_impute, "Leave zero rates at zero"}},
                                           "Impute zero rates for final model files (*excluding* JSON)?:");


// Prompt for output format(s) //
pogofit.output_format  = io.SelectAnOption ({ {pogofit.output_hyphy, "HyPhy-formatted model (extension `.fitted_model`)"},
                                              {pogofit.output_paml, "PAML-formatted model (extension `.paml`)"},
                                              {pogofit.output_raxml, "RAXML-formatted model (extension `.raxml`)"},
                                              {pogofit.output_all, "Output all file formats"}},
                                               "Select an output format for the fitted model:");

pogofit.options  = {pogofit.options.baseline_model: pogofit.baseline_model,
                    pogofit.options.imputation: pogofit.imputation};
pogofit.input    = {terms.json.file: pogofit.input_file,
                     pogofit.options.number_of_datasets: pogofit.file_list_count,
                     pogofit.options.dataset_information: {}};


/********************************************************************************************************************/
/********************************************* ANALYSIS BEGINS HERE *************************************************/
/********************************************************************************************************************/


pogofit.startTimer (pogofit.timers, "Total time");


/******************************************* STEP ONE *******************************************************
        Perform an initial fit of the Baseline model+F+4G for full data. Estimates shared alpha as well.
*************************************************************************************************************/
console.log("\n\n[PHASE 1] Performing initial branch length optimization using " + pogofit.baseline_model);

pogofit.startTimer (pogofit.timers, pogofit.baseline_phase);

pogofit.baseline_fit = pogofit.fitBaselineTogether();

pogofit.stopTimer (pogofit.timers, pogofit.baseline_phase);
/*************************************************************************************************************/
/*************************************************************************************************************/



/******************************************* STEP TWO *******************************************************
        Fit a full GTR model to all dataset(s) jointly, using the baseline model as initial rates
*************************************************************************************************************/
console.log("\n\n[PHASE 2] Optimizing protein model");

pogofit.startTimer (pogofit.timers, pogofit.final_phase);

pogofit.gtr_fit = pogofit.fitGTR_fixalpha(pogofit.baseline_fit);
                                                                                                                              
pogofit.stopTimer (pogofit.timers, pogofit.final_phase);
/*************************************************************************************************************/
/*************************************************************************************************************/

console.log("\n\nOptimization complete! Saving results.");



/*********************** Save custom model to file(s) as specified **************************/
pogofit.final_efv = (pogofit.gtr_fit[terms.efv_estimate])["VALUEINDEXORDER"][0];

if (pogofit.imputation == pogofit.impute) 
{
    // Tree length for each alignment
    pogofit.tree_lengths = {};
    utility.ForEachPair (pogofit.gtr_fit[terms.branch_length], "_part_", "_value_",
    '
        pogofit.tree_lengths[_part_] = math.Sum(utility.Map (_value_, "_data_",
        "
            _data_ [terms.fit.MLE]
        " 
        ))
    '
    );


    // Site counts for each alignment
    pogofit.site_counts = {};
    utility.ForEachPair( (pogofit.analysis_results[terms.json.input])[pogofit.options.dataset_information], "_key_", "_value_",
    '
        pogofit.site_counts[_key_] = _value_[terms.json.sites];
    '
    );
    pogofit.final_rij = pogofit.extract_rates_imputation();
}
else 
{
    pogofit.final_rij = pogofit.extract_rates();
}
pogofit.write_model_to_file();

/************************************* Save analysis JSON ***********************************/
pogofit.stopTimer (pogofit.timers, "Total time");
pogofit.analysis_results[terms.json.input]    = pogofit.input;
pogofit.analysis_results[terms.json.timers]   = pogofit.timers;
io.SpoolJSON(pogofit.analysis_results, pogofit.json_file);

