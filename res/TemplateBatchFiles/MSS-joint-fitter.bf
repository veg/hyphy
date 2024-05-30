
/* 1. SETUP
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/* 1a. Initial Setup
------------------------------------------------------------------------------*/
RequireVersion ("2.5.48");

LoadFunctionLibrary ("libv3/all-terms.bf");
LoadFunctionLibrary ("libv3/convenience/regexp.bf");
LoadFunctionLibrary ("libv3/convenience/math.bf");
LoadFunctionLibrary ("libv3/IOFunctions.bf");
LoadFunctionLibrary ("libv3/UtilityFunctions.bf");
LoadFunctionLibrary ("libv3/tasks/trees.bf");
LoadFunctionLibrary ("libv3/tasks/alignments.bf");
LoadFunctionLibrary ("libv3/tasks/estimators.bf");
LoadFunctionLibrary ("SelectionAnalyses/modules/io_functions.ibf");
LoadFunctionLibrary ("libv3/tasks/mpi.bf");
LoadFunctionLibrary ("libv3/convenience/random.bf");
LoadFunctionLibrary("libv3/models/codon/MSS.bf");



mss_selector.analysisDescription = {terms.io.info : "mss_joint_fitter : Performs a joint MSS model fit to several genes jointly.",
                           terms.io.version : "0.0.1",
                           terms.io.reference : "TBD",
                           terms.io.authors : "Sergei L Kosakovsky Pond",
                           terms.io.contact : "spond@temple.edu",
                           terms.io.requirements : "A collection of coding alignments with trees in the corresponding alignment files "
                          };


mss.json = {   terms.json.analysis: mss_selector.analysisDescription,
                terms.json.input: {},
            };
            

utility.SetEnvVariable ("OPTIMIZE_SUMMATION_ORDER_PARTITION", 100); 
// don't spend too much time optimizing column ordering.

/* 1b. User Input and data load
------------------------------------------------------------------------------*/
io.DisplayAnalysisBanner (mss_selector.analysisDescription);

KeywordArgument ("filelist","List of files to include in this analysis");
mss_selector.file_list = io.get_a_list_of_files(io.PromptUserForFilePathRead ("List of files to include in this analysis"));

mss_selector.file_count = utility.Array1D (mss_selector.file_list);
io.CheckAssertion("mss_selector.file_count >= 1", "A non-empty file list is required");

io.ReportProgressMessageMD("mss", "data" , "* Loaded a list with **" + mss_selector.file_count  + "** files");

KeywordArgument ("code",        "Which genetic code should be used", "Universal");  
mss.genetic_code = alignments.LoadGeneticCode (None);

KeywordArgument ("omega", "How should alignment-level omega be treated?", "Fix");

mss.omega_option = io.SelectAnOption ({
                                        {"Fix", "[Default] Fix omega estimates at those obtained with the standard MG94xREV model"}
                                        {"Fit", "Fit omega (per alignment) together with the MSS model"}
                                    }, "How should alignment-level omega be treated?");


mss.file_records = {};
mss.file_info = {};
mss.tree_info = {};
mss.file_prefix = {};
mss.fits = {};
mss.filters = {};
mss.parameters = 0;
mss.baselineLL = 0;
mss.is_bic = TRUE;

KeywordArgument ("output", "Write the resulting JSON to this file",None); 
mss.json = io.ReadFromOrCreate ("Save the resulting JSON file to", mss.json);

mss.json - terms.data.value;

mss_selector.file_list = utility.DictToArray (mss_selector.file_list);

mss_selector.queue = mpi.CreateQueue (
                            {
                            "Headers" : {{"libv3/UtilityFunctions.bf","libv3/tasks/alignments.bf",
                            "libv3/tasks/trees.bf","SelectionAnalyses/modules/io_functions.ibf",
                            "libv3/tasks/estimators.bf","libv3/models/codon/MSS.bf"}},
                            "Variables" : {{}}
                            }
                        );

function mss.handle_initial_fit (filepath, namespace, genetic_code, run_fit) {

     utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", "/dev/null");
     ExecuteCommands ('
        `namespace`.json = {};
         namespace `namespace` {
            scaler_prefix = "MSS.scaler";
            KeywordArgument ("code",        "Which genetic code should be used", "`genetic_code`");  
            KeywordArgument ("alignment",   "Load this alignment", "`filepath`");  
            utility.ToggleEnvVariable ("ALWAYS_RELOAD_FUNCTION_LIBRARIES", TRUE);
            LoadFunctionLibrary ("SelectionAnalyses/modules/shared-load-file.bf");
            utility.ToggleEnvVariable ("ALWAYS_RELOAD_FUNCTION_LIBRARIES", None);

            load_file ({utility.getGlobalValue("terms.prefix"): "`namespace`", 
                        utility.getGlobalValue("terms.settings") : {utility.getGlobalValue("terms.settings.branch_selector") : "selection.io.SelectAllBranches"}});
        
            if (^"`&run_fit`") {
                doGTR ("`namespace`");
                doPartitionedMG ("`namespace`", FALSE);
            }
    
        };
    ');
    utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", None);
    return {
        "path" : filepath,
        "pt" : Eval (namespace + ".trees"),
        "init" : Eval (namespace + ".partitioned_mg_results")
    };
}

lfunction mss.store_initial_fit (node, result, arguments) {

    mss_selector.print_row = {
            5,
            1
    };

    (^"mss.fits") [result["path"]] = result["init"];
    (^"mss.tree_info") [result["path"]] = result ["pt"];
    
    ^"mss.parameters" += (result["init"])[^"terms.parameters"];
    ^"mss.baselineLL" += (result["init"])[^"terms.fit.log_likelihood"];
    //^"mss.sample_size" += (result["init"])[^"terms.data.sample_size"];
    mss_selector.print_row [0] = result["path"];
    info = (^"mss.file_records")[result["path"]];
    mss_selector.print_row [1] = info [^"terms.data.sequences"];
    mss_selector.print_row [2] = info [^"terms.data.sites"];
    
    mss.T = 0;
    for (mss.b,mss.bi; in;  ((result["init"])[^"terms.branch_length"])["0"]) {
        mss.T += mss.bi [^"terms.fit.MLE"];
    }
    mss_selector.print_row [3] = Format (mss.T, 8, 3);
    mss_selector.print_row [4] = Format ((result["init"])[^"terms.fit.log_likelihood"], 12, 4);
    fprintf(stdout, io.FormatTableRow(mss_selector.print_row, ^"mss_selector.table_output_options"));
    return TRUE;
}


mss_selector.table_output_options = {
        utility.getGlobalValue("terms.table_options.header"): 1,
        utility.getGlobalValue("terms.table_options.minimum_column_width"): 16,
        utility.getGlobalValue("terms.table_options.align"): "left",
        utility.getGlobalValue("terms.table_options.column_widths"): {
            "0": 120,
            "1" :10,
            "2" :10,
            "3" :10,
            "4" :15
        }
    };
    
mss_selector.header = {
        {"Filepath"}
        {"Sequences"}
        {"Sites"}
        {"TreeLength"},
        {"Log(L)"}
    };

KeywordArgument ("model", "Substitution model to use","SynREV"); 
    
mss.model_option = io.SelectAnOption ({
                                        {"SynREV", "[Default] SynREV model (one rate per aa)"}
                                        {"SynREVCodon", "SynREVCodon (one rate per codon pair)"}
                                    }, "Which model?");

  
    
ExecuteCommands ( "mss.codon_classes = model.codon.MSS.prompt_and_define (terms.global, mss.genetic_code[terms.code])", 
                                       {"--mss-type" : mss.model_option}
                                    );
                                        
io.ReportProgressMessageMD("mss", "fit0" , "Individual file statistics and simple model fits\n");


fprintf(stdout, io.FormatTableRow(mss_selector.header, mss_selector.table_output_options));
mss_selector.table_output_options[utility.getGlobalValue("terms.table_options.header")] = FALSE;

mss.path_ordering = {};
mss.filter_names = {};
mss.trees = {};
mss.order2path = {};

                                    

for (mss.counter, mss_selector.path; in; mss_selector.file_list) {
     //io.ReportProgressBar("", "Loading datafile  " + (+mss.counter+1) + "/" +  mss_selector.file_count);
    
     mss.namespace = "mss_" + mss.counter;
     mss.handle_initial_fit (mss_selector.path, mss.namespace, mss.genetic_code[terms.id], FALSE); 
     mss.file_records [mss_selector.path] = ^"`mss.namespace`.codon_data_info";
     mss.sample_size += (^"`mss.namespace`.codon_data_info")[terms.data.sample_size];

   
     mpi.QueueJob (mss_selector.queue, "mss.handle_initial_fit", {"0" : mss_selector.path,
                                                     "1" : mss.namespace,
                                                     "2" : mss.genetic_code[terms.id],
                                                     "3" : TRUE},
                                                     "mss.store_initial_fit");
    
    
    
    mss.file_info [mss_selector.path]    = ^"`mss.namespace`.json";
    mss.filters[mss_selector.path] = ^"`mss.namespace`.filter_names";
    mss.file_prefix  [mss_selector.path] = mss.namespace;
    mss.order2path [Abs (mss.path_ordering)] = mss_selector.path;
    mss.path_ordering [mss_selector.path] = Abs (mss.path_ordering);
    mss.filter_names + (^"`mss.namespace`.filter_names")[0];
}   



//selection.io.getIC(fit[terms.fit.log_likelihood], fit[terms.parameters] + xp, ss)

//io.ClearProgressBar();
mpi.QueueComplete (mss_selector.queue);

mss.baseline_AIC = mss.getIC (mss.baselineLL, mss.parameters, mss.sample_size, mss.is_bic);

io.ReportProgressMessageMD("mss", "fit0done" , "Baseline model fit");
io.ReportProgressMessageMD("mss", "fit0done", "**LogL = " + mss.baselineLL + "**"  + ", **AIC-c = " + mss.baseline_AIC + "** \n");

mss.json ["baseline"] = {terms.json.log_likelihood : mss.baselineLL, terms.json.AICc : mss.baseline_AIC};

function mss.MSS_generator (type) {
    mss.codon_classes [utility.getGlobalValue("terms.model.MSS.normalize")] = TRUE;
    model = Call ("models.codon.MSS.ModelDescription",type, mss.genetic_code[terms.code], mss.codon_classes);
    return model;
}

mss.tree_objects  = {};
mss.fit_objects   = {};
mss.lf_components = {   
                            2*mss_selector.file_count, 
                            1
                    };

mss.lf_id = "MSS.COMPOSITE.LF";
mss.model_map_overall = {};


//#profile START;
mss.constraints = {};

for (mss.counter, mss_selector.path; in; mss_selector.file_list) {
     mss.prefix = mss.file_prefix  [mss_selector.path];
    
     console.log (">" + mss.counter + " / " +  mss_selector.path);
    
     mss.model_name = "`mss.prefix`.model";
     mss.tree_name = "`mss.prefix`.tree";
     mss.lf_components[2*mss.counter] = mss.filter_names[mss.counter];
     mss.lf_components[2*mss.counter+1] = mss.tree_name;
     utility.ExecuteInGlobalNamespace ("`mss.model_name ` = 0");
     ^(mss.model_name) = model.generic.DefineModel("mss.MSS_generator", mss.prefix + ".model_object", {
            "0": "terms.global"
        }, mss.filter_names[mss.counter], None);
        
     mss.model_map_overall [mss.prefix + ".model_object"] = ^(mss.model_name);
        
     mss.model_map = { mss.prefix + ".model_object" :  ^(mss.model_name)};
        
     model.ApplyModelToTree(mss.lf_components[2 * mss.counter + 1], (mss.tree_info[mss_selector.path])[0], {
            "default": ^(mss.model_name)
        }, None);
        
        
     initial_values = mss.fits[mss_selector.path];
    
     for (mss.p; in; estimators.GetGlobalMLE_RegExp( initial_values, terms.parameters.omega_ratio)) {
        (initial_values[terms.global])[terms.parameters.omega_ratio] = mss.p;
     }       
     
     mss.model_map["estimators.SetCategory"][""];
     estimators.ApplyExistingEstimates.set_globals = {};
     mss.model_map["estimators.SetGlobals"][""];
     mss.bl_scaler = "mss.bl_scaler_" + mss.counter; 
     
     parameters.DeclareGlobal (mss.bl_scaler, None);
     
     mss.constraint_count =  estimators.ApplyExistingEstimatesToTree (mss.lf_components[2 * mss.counter + 1], 
                                              mss.model_map, 
                                              (initial_values [terms.branch_length])[0],
                                              mss.bl_scaler,
                                              {});
                                           
     if (mss.omega_option == "Fix") {
        models.FixParameterSetRegExp (terms.parameters.omega_ratio, ^(mss.model_name));
        mss.constraint_count += 1;
     }
     models.FixParameterSetRegExp (terms.nucleotideRatePrefix, ^(mss.model_name));
      
     
    if (mss.counter == 0) {
        mss.reference_set   =  ^(mss.model_name);
        mss.mss_rate_list   =  model.GetParameters_RegExp( ^(mss.model_name),"^" + terms.MeanScaler(""));
        mss.model_dimension = utility.Array1D (mss.mss_rate_list);
        mss.scaler_prefix = 'mss.scaler_parameter_';
        mss.scaler_mapping = {};
        mss.scaler_index = { mss.model_dimension , 1};
        for (mss.i, mss.v; in; mss.mss_rate_list) {
            mss.scaler_mapping [mss.i] = Abs (mss.scaler_mapping);
            mss.scaler_index   [ mss.scaler_mapping [mss.i]] = mss.v;
        }
        for (mss.i2 = 0; mss.i2 < mss.model_dimension; mss.i2 += 1) {
            parameters.DeclareGlobal (mss.scaler_prefix + mss.i2, None);
        }
        mss.json ["mapping"] = mss.scaler_index;
    } else {
        //utility.getGlobalValue ("terms.parameters.synonymous_rate");
        mss.constraints  * models.BindGlobalParametersDeferred ({"0" : mss.reference_set , "1" : ^(mss.model_name)}, "^" + terms.MeanScaler(""));
     }
}

parameters.BatchApplyConstraints (mss.constraints);

//#profile _hyphy_profile_dump;
//utility.FinishAndPrintProfile (_hyphy_profile_dump);
utility.SetEnvVariable ("AUTO_PARALLELIZE_OPTIMIZE", 1);

utility.ExecuteInGlobalNamespace ("LikelihoodFunction `mss.lf_id` = (`&mss.lf_components`)");
 
VERBOSITY_LEVEL                 = 1;
USE_LAST_RESULTS                = 1;
Optimize                        (res, ^mss.lf_id);

mss.joint_AIC = mss.getIC (res[1][0], mss.parameters + res[1][1], mss.sample_size, mss.is_bic);

mss.json ["joint"] = {terms.json.log_likelihood : res[1][0], terms.json.AICc : mss.joint_AIC};


io.ReportProgressMessageMD("mss", "fitAlldone" , "Joint model fit");
io.ReportProgressMessageMD("mss", "fitAlldone", "**LogL = " + res[1][0] + "**"  + ", **AIC-c = " + mss.joint_AIC + "** \n");

mss.json["joint-model"] = estimators.ExtractMLEsOptions (mss.lf_id, mss.model_map_overall, {terms.globals_only : TRUE});
//estimators.RemoveBranchLengthConstraints (mss.json["joint-model"]);

function pfilter (key, value) {
    if (key[0] != "[" || (key $ "ratio")[0] >= 0)  {
        mss.filtered [key] = value;
    }
}

mss.filtered = {};
((mss.json["joint-model"])[terms.global])["pfilter"][""];
(mss.json ["joint-model"])[terms.global] = mss.filtered;


io.SpoolJSON (mss.json, mss.json[terms.data.file]);  

KeywordArgument ("save-fit", "Write the resulting model fit file to this (large!) file", "/dev/null");
io.SpoolLFToPath(mss.lf_id, io.PromptUserForFilePath ("Save model fit to this file ['/dev/null' to skip]"));

//------------------------------------------------------------------------------------------------------------------------

lfunction mss.getIC(logl, params, samples, is_bic) {
    if (is_bic) {
        return -2 * logl + Log (samples) * params;
    }
    return -2 * logl + 2 * samples / (samples - params - 1) * params;
}
