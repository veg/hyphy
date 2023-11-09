
/* 1. SETUP
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/* 1a. Initial Setup
------------------------------------------------------------------------------*/
RequireVersion ("2.5.35");

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

namespace terms.mss_selector {
    mapping = "mapping";
    masterList = "masterList";
    files = "files";
};


mss_selector.analysisDescription = {terms.io.info : "mss_selector : Genetic Algorithms for Recombination Detection. Implements a heuristic
    approach to screening alignments of sequences for recombination, by using the CHC genetic algorithm to search for phylogenetic
    incongruence among different partitions of the data. The number of partitions is determined using a step-up procedure, while the
    placement of breakpoints is searched for with the GA. The best fitting model (based on c-AIC) is returned; and additional post-hoc
    tests run to distinguish topological incongruence from rate-variation. v0.2 adds  and spooling results to JSON after each breakpoint search conclusion",
                           terms.io.version : "0.0.1",
                           terms.io.reference : "**Perform a Genetic Algorithm search model selection for an MSS model on a collection of alignments",
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

KeywordArgument ("ic",        "Which IC should be used for scoring", "AIC-c");  
mss.ic_score = io.SelectAnOption  ({"AIC-c" : "Small Sample AIC score", 
                                   "BIC" : "Bayesian Information Criterion (more conservative)"}, "Which IC should be used for scoring");

mss.is_bic = mss.ic_score == "BIC";

KeywordArgument ("classes",        "How many rate classes should be considered", "2");  
mss.rate_classes = io.PromptUser ("How many rate classes should be considered?", 2, 2, 6, TRUE);


mss.file_records = {};
mss.file_info = {};
mss.tree_info = {};
mss.file_prefix = {};
mss.fits = {};
mss.filters = {};
mss.parameters = 0;
mss.baselineLL = 0;
mss.sample_size = 0;

KeywordArgument ("output", "Write the resulting JSON to this file",None); 
mss.json = io.ReadFromOrCreate ("Save the resulting JSON file to", mss.json);

//console.log (mss.json[terms.data.file]);


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
    
    
ExecuteCommands ( "mss.codon_classes = model.codon.MSS.prompt_and_define (terms.global, mss.genetic_code[terms.code])", 
                                       {"--mss-type" : "SynREVFull"}
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

function mss.MSS_generator (type) {

    model = Call ("models.codon.MSS.ModelDescription",type, mss.genetic_code[terms.code], mss.codon_classes);
    return model;
}

mss.tree_objects = {};
mss.fit_objects = {};
mss.lf_components = {   
                            2*mss_selector.file_count, 
                            1
                    };

mss.lf_id = "MSS.COMPOSITE.LF";

for (mss.counter, mss_selector.path; in; mss_selector.file_list) {
     mss.prefix = mss.file_prefix  [mss_selector.path];
    
     mss.model_name = "`mss.prefix`.model";
     mss.tree_name = "`mss.prefix`.tree";
     mss.lf_components[2*mss.counter] = mss.filter_names[mss.counter];
     mss.lf_components[2*mss.counter+1] = mss.tree_name;
     utility.ExecuteInGlobalNamespace ("`mss.model_name ` = 0");
     ^(mss.model_name) = model.generic.DefineModel("mss.MSS_generator", mss.prefix + ".model_object", {
            "0": "terms.global"
        }, mss.filter_names[mss.counter], None);
        
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
     mss.constraint_count =  estimators.ApplyExistingEstimatesToTree (mss.lf_components[2 * mss.counter + 1], 
                                              mss.model_map, 
                                              (initial_values [terms.branch_length])[0],
                                              terms.model.branch_length_constrain,
                                              TRUE);
                                              
     models.FixParameterSetRegExp (terms.parameters.omega_ratio, ^(mss.model_name));
     models.FixParameterSetRegExp (terms.nucleotideRatePrefix, ^(mss.model_name));
      
     
    if (mss.counter == 0) {
        mss.reference_set   =  ^(mss.model_name);
        mss.mss_rate_list   =  model.GetParameters_RegExp( ^(mss.model_name),terms.parameters.synonymous_rate + terms.model.MSS.syn_rate_within);
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
        mss.json [terms.mss_selector.mapping] = mss.scaler_index;
    } else {
        //utility.getGlobalValue ("terms.parameters.synonymous_rate");
        models.BindGlobalParameters ({"0" : mss.reference_set , "1" : ^(mss.model_name)}, terms.parameters.synonymous_rate + terms.model.MSS.syn_rate_within);
    }
}



utility.ExecuteInGlobalNamespace ("LikelihoodFunction `mss.lf_id` = (`&mss.lf_components`)");
 
namespace mss {

    // GA.1: Setup global parameters
    masterList     = {};
    masterRateList = {};
    populationSize = 32;  

    mutationRate = 0.2; // the GARD paper said "15% of randomly selected bits were toggled"...
    rateOfMutationsTharAreSmallShifts = 0.5; // some mutations are a new random break point; some are small shifts of the break point to an adjacent location.
    maxGenerationsAllowedAtStagnant_cAIC = 100;
    maxGenerationsAllowedWithNoNewModelsAdded = maxGenerationsAllowedAtStagnant_cAIC $ 4;

    maxFailedAttemptsToMakeNewModel = 7;
    cAIC_diversityThreshold   = 1e-6;
    bestOverallModelSoFar = {model_dimension, 1};
    numberOfClassesBeingEvaluated = rate_classes;
    cAIC_improvementThreshold = 0.01;
    parentModels = mss.GA.initializeModels(numberOfClassesBeingEvaluated, populationSize, model_dimension, None);
    
    null_model = Eval((Rows (parentModels))[0])*0;   
    masterList[null_model] = baseline_AIC;
    masterRateList[null_model] = {{baseline_AIC__,baselineLL__,1}};
    
    terminationCondition = FALSE;
    generation = 0;
    bestOverall_cAIC_soFar = baseline_AIC;
    previousBest_cAIC = ^"math.Infinity";
    startTime = Time(1);
    
    while(terminationCondition == FALSE) {
        // GA.2.b.1 Produce the next generation of models with recombination.
        childModels = GA.recombineModels(parentModels, populationSize, numberOfClassesBeingEvaluated);
        //console.log ("\n" + Abs (childModels));

        // GA.2.b.2 Select the fittest models.
        //console.log ("\n\nPARENT MODELS");
        //console.log (parentModels);
        interGenerationalModels = parentModels;
        interGenerationalModels * childModels;
        
        
        mss.GA.evaluateModels(interGenerationalModels);
        
        //console.log (masterList);
        //console.log (interGenerationalModels);
        
        selectedModels = mss.GA.selectModels(interGenerationalModels, populationSize);

        // GA.2.b.3 Evaluate convergence for this modelSet initialization
        // If converged, produce a new parent modelSet (keeping the current fitest)
        if (selectedModels == parentModels) {
            generationsNoNewModelsAdded += 1;
        } else {
            generationsNoNewModelsAdded = 0;
        }

        currentBest_individual = Min(interGenerationalModels,1);
        currentBest_cAIC  = currentBest_individual["value"];
        currentBest_model = Eval(currentBest_individual["key"]);

        if ( (math.minNormalizedRange(selectedModels) < cAIC_diversityThreshold) || (generationsNoNewModelsAdded > maxGenerationsAllowedWithNoNewModelsAdded) ) {
            //console.log ("\n\nMUTATE\n\n" );
            parentModels = mss.GA.generateNewGenerationOfModelsByMutatingModelSet(selectedModels, numberOfClassesBeingEvaluated, mutationRate, rateOfMutationsTharAreSmallShifts);
            if (Abs(parentModels) == 1) {
                //console.log ("REINIT PARENT MODELS");
                parentModels = mss.GA.initializeModels(numberOfClassesBeingEvaluated, populationSize - 1, model_dimension, bestOverallModelSoFar);
                parentModels [currentBest_model] = currentBest_cAIC;
            }
            differenceThreshold = numberOfClassesBeingEvaluated / 4;
        } else {
            parentModels = selectedModels;
        }

        // GA.2.b.4 Evaluate convergence for this number of break points
        // If converged, move on to n+1 break points or end analysis

        if (previousBest_cAIC - currentBest_cAIC < cAIC_improvementThreshold) {
            generationsAtCurrentBest_cAIC += 1;
            if (generationsAtCurrentBest_cAIC >= maxGenerationsAllowedAtStagnant_cAIC) {
                terminationCondition = TRUE;
            }
        } else {
            //console.log ("RESETTING generationsAtCurrentBest_cAIC" + previousBest_cAIC + ", " + currentBest_cAIC + ", " + cAIC_improvementThreshold + "\n" + currentBest_model + "\n");
            generationsAtCurrentBest_cAIC = 0;
        }


        if (previousBest_cAIC > 0 && previousBest_cAIC < currentBest_cAIC) {
            io.CheckAssertion("mss.previousBest_cAIC >= mss.currentBest_cAIC", "Internal error in MSS-Select -- c-AIC INCREASED between two consecutive generations");
        }
        previousBest_cAIC = currentBest_cAIC;
        generation += 1;
        
        if (bestOverall_cAIC_soFar-currentBest_cAIC > 0) {
             io.ReportProgressBar ("MSS-GA", Format (numberOfClassesBeingEvaluated,3,0) + ' rates [ generation ' +  Format (generation, 6, 0) + ", total models " + Format (Abs (masterList), 8, 0) + ", " 
                                    + Format (generationsAtCurrentBest_cAIC/ maxGenerationsAllowedAtStagnant_cAIC*100, 4, 0) +"% converged, " + Format (Abs (masterList)/(Time(1)-startTime),5,2) + "/sec]"
                                    + ". Min (c-AIC) = " + Format (currentBest_cAIC, 12,4) + " [ delta = " + Format (bestOverall_cAIC_soFar-currentBest_cAIC, 8, 2) + "], rate profile " + Join ("", currentBest_model));

        } else {
            io.ReportProgressBar ("MSS-GA", Format (numberOfClassesBeingEvaluated,3,0) + ' rates [ generation ' +  Format (generation, 6, 0) + ", total models " + Format (Abs (masterList), 8, 0) + ", " 
            + Format (generationsAtCurrentBest_cAIC/ maxGenerationsAllowedAtStagnant_cAIC*100, 4, 0) +"% converged, " + Format (Abs (masterList)/(Time(1)-startTime),5,2) + "/sec]"
                                    + ". Min (c-AIC) = " + Format (currentBest_cAIC, 12,4) + " [no improvement], rate profile " + Join ("", currentBest_model));
        }
        
    }
    io.ClearProgressBar();
            
}

mss.json [terms.mss_selector.masterList]  = mss.masterRateList;
mss.json [terms.mss_selector.files] = mss_selector.file_list;

io.SpoolJSON (mss.json, mss.json[terms.data.file]);
 
return 0; 

// ---- Genetic Algorithm (GA) Functions ----

/**
 * @name gard.GA.initializeModels
 * @param {Number} numberOfClasses
 * @param {Number} populationSize
 * @param {Number} dimension
 * @returns a {Dictonary} initializedModels
 */
lfunction mss.GA.initializeModels (numberOfClasses, populationSize, dimension, seed) {

    initializedModels = {};
    modelNumber       = 0;

    if (null != seed) {
       for (; modelNumber < populationSize$2; modelNumber += 1) {
         do {
            model = seed;
            for (i,v; in; model) {
                if (Random (0,1) < 0.2) {
                    model[i] = numberOfClasses - 1;
                }
            }
            model = mss.GA.normalize_and_validate (model, numberOfClasses);
         } while (None == model);

        initializedModels[breakPoints] = ^"math.Infinity";
       }
    }
    for (; modelNumber < populationSize; modelNumber += 1) {
        initializedModels [mss.GA.generate_model (dimension, numberOfClasses)] = ^"math.Infinity";
    }
    return initializedModels;
}

//-------------------------------------------------------------------------------------------
lfunction mss.GA.fit_model (model, lfid, xp, ss) {
    
    if (^"mss.masterRateList" / model) {
        //console.log ("\n**CACHED**\n");
        return (^"mss.masterRateList")[model];
    }
    
    for (i,v; in; model) {
       parameters.SetConstraint ((^"mss.scaler_index")[i], ^"mss.scaler_prefix" + v, "");
    }
    parameters.SetConstraint (^"mss.scaler_prefix" + 0, "1", "");
    utility.SetEnvVariable ("AUTO_PARALLELIZE_OPTIMIZE", 1.); 
    Optimize (res, ^lfid);
    return_expr = {1,res[1][1] + 3};
    return_expr [0] = mss.getIC (res[1][0], xp + res[1][1] ,ss, ^"mss.is_bic"); 
    return_expr [1] = res[1][0]; 
    
    for (i = 0; i < res[1][1] + 1; i+=1) {
        return_expr[2+i] = Eval (^"mss.scaler_prefix" + i);
    }
    
    (^"mss.masterList")[model] = return_expr[0];
    (^"mss.masterRateList")[model] = return_expr;
    
    
    return return_expr;
}

//-------------------------------------------------------------------------------------------

lfunction mss.GA.normalize_and_validate (model, rates) {
    rate_counter = {rates,2}["_MATRIX_ELEMENT_COLUMN_*100000"];
    
    for (i, k; in; model) {
        rate_counter [k][0] += 1;
        rate_counter [k][1] = Min (i, rate_counter [k][1]);
    }
    
    if (Min (rate_counter[-1][0], 0) > 0) { 
        rate_counter = (rate_counter["(_MATRIX_ELEMENT_COLUMN_==0)*_MATRIX_ELEMENT_ROW_+(_MATRIX_ELEMENT_COLUMN_==1)*_MATRIX_ELEMENT_VALUE_"])%1;
        rate_counter = (rate_counter["(_MATRIX_ELEMENT_COLUMN_==1)*_MATRIX_ELEMENT_ROW_+(_MATRIX_ELEMENT_COLUMN_==0)*_MATRIX_ELEMENT_VALUE_"])%0;
        return model ["rate_counter[_MATRIX_ELEMENT_VALUE_][1]"];           
    } 
    return None;
}

//-------------------------------------------------------------------------------------------

lfunction mss.GA.generate_model (dimension, rates) {
    while (TRUE) {
        model = mss.GA.normalize_and_validate({1, dimension}["Random (0, rates)$1"], rates);
        if (None != model) {
            return model;
        }
        
    }
    return None;
}

/**
 * @name mss.GA.evaluateModels
 * @param {Dictonary} models (some models may have cAIC scores already)
 * @returns nothing... the gard.GA.storeMultiBreakPointModelResults function gets called which updates the intergenerationalModels object
 */
function mss.GA.evaluateModels (models) {
    
    for (modelID, cAIC; in; models)  {
        if (cAIC == math.Infinity) {
            models[modelID] = (mss.GA.fit_model ( Eval (modelID), mss.lf_id, mss.parameters, mss.sample_size))[0];
            console.log ("\t" + Join("", Eval (modelID)) + " IC = " + Format (models[modelID], 15, 4) + ", delta = " + Format ((-(models[modelID]) + mss.bestOverall_cAIC_soFar), 8, 2));        }
    }
    
}

//-------------------------------------------------------------------------------------------

/**
 * @name gard.GA.recombineModels
 Given a set of models create a new generation of models by iteratively recombining two random parent models.
 The child model will have a random subset of the breakpoints from the parents

 * @param {Matrix} parentModels
 * @param {Number} populationSize
 * @returns a {Dictonary} childModels
 */
lfunction mss.GA.recombineModels (parentModels, populationSize, rateCount) {
    parentModelIds = utility.Keys(parentModels);
    numberOfParentModels = Abs(parentModels);
    firstModel = gard.Helper.convertMatrixStringToMatrix(parentModelIds[0]);

    childModels = {};
    for(modelIndex=0; modelIndex<populationSize; modelIndex+=1) {
        failedAttempts = 0;
        model = None;
        while(None == model && failedAttempts < ^"mss.maxFailedAttemptsToMakeNewModel") {
            parentModel1 = parentModelIds[Random(0, numberOfParentModels-.0001)$1];
            parentModel2 = parentModelIds[Random(0, numberOfParentModels-.0001)$1];

            while(parentModel1 == parentModel2) {
                parentModel2 = parentModelIds[Random(0, numberOfParentModels-.0001)$1];
            }
            
            parentModel1 = Eval(parentModel1);
            parentModel2 = Eval(parentModel2);

            model = parentModel1;
            for (i, v; in; model) {
                if (Random (0,1) < 0.5) {
                    model[i] = parentModel2[i];
                }   
            }
            
            model = mss.GA.normalize_and_validate (model, rateCount);
            if (None == model /*|| ^"mss.masterList" / model*/) {
                failedAttempts += 1;
                model = None;
            }

        }
        if (None != model) {
            childModels[model] = ^"math.Infinity";
        }

    }
    return childModels;
}

/**
 * @name gard.GA.selectModels
 * @param {Dictionary} evaluatedModels
 * @param {Number} numberOfModelsToKeep
 * @returns a {Matrix} selectedModels
 */
lfunction mss.GA.selectModels (evaluatedModels, numberOfModelsToKeep) {
/**
    evaluatedModels: model -> score
*/
    // sort model scores
    N = utility.Array1D (evaluatedModels);
    idx2score = {};
    modelIds = utility.Keys(evaluatedModels);
    for (i = 0; i < N; i+=1) {
        idx2score[i] = evaluatedModels[modelIds[i]];
    }
    idx2score = utility.DictToSortedArray (idx2score);
    selectedModels = {};
    numberOfModelsToKeep = Min (numberOfModelsToKeep, N);
    for(i=0; i<numberOfModelsToKeep; i+=1) {
        selectedModels[modelIds[idx2score[i][1]]] = idx2score[i][0];
    }
    return selectedModels;
}

/**
 * @name gard.GA.generateNewGenerationOfModelsByMutatingModelSet
 Keeps the current best performing model
 Generates new models for the rest of the population by randomly mutating some of the parent models break points
 // TODO: decide if we want some logic in here to more meaningfully mutate.
 * @param {Dictonary} parentModels
 * @param {Number} numberOfPotentialBreakPoints
 * @param {Number} mutationRate
 * @param {Number} rateOfMutationsThatAreSmallShifts
 * @returns a {Dictonary} the new generation of models
 */
 
lfunction mss.GA.generateNewGenerationOfModelsByMutatingModelSet(parentModels, numberOfRateClasses, mutationRate, smallShift) {



    modelIds = utility.Keys(parentModels);
    populationSize = Columns(modelIds);
    firstModel =Eval (modelIds[0]);
    numberOfBreakPoints = utility.Array1D(firstModel);
 
    // Generate a new set of models
    nextGenOfModels = {};
    for(ind=0; ind<populationSize-1; ind += 1) {
        model = None;
        while(None == model && failedAttempts < ^"mss.maxFailedAttemptsToMakeNewModel") {
            model = Eval(modelIds[ind]);
            if (Random (0,1) < smallShift) {
              i2 = Random (0, utility.Array1D (model))$1;
              do {
                v = Random (0, numberOfRateClasses)$1;
              } while (v == model[i2]);
              model[i2] = v;
            } else {
              for (i, v; in; model) {
                if (Random (0,1) < mutationRate) {
                    model[i] = Random (0,numberOfRateClasses)$1;
                }   
              }     
            }
           
            model = mss.GA.normalize_and_validate (model, numberOfRateClasses);
            if (None == model) {
                failedAttempts += 1;
                model = None;
            }           

        }
        if (None != model) {
            nextGenOfModels[model] = ^"math.Infinity";
        }

    }

    // Add the most fit model from the previous set

    //console.log ("STORING CURRENT BEST INDIVIDUAL");
    cAIC_scores = utility.Values(parentModels);
    //console.log (parentModels);
    lowest_cAIC_index = Min(cAIC_scores,1)[1];
    //console.log (Min(cAIC_scores,1)[1]);
    selectedModelId = modelIds[lowest_cAIC_index];
    //console.log (selectedModelId);
    nextGenOfModels[selectedModelId] = cAIC_scores[lowest_cAIC_index];
    //console.log (nextGenOfModels[selectedModelId]);
    return nextGenOfModels;
}


//------------------------------------------------------------------------------------------------------------------------

lfunction mss.getIC(logl, params, samples, is_bic) {
    if (is_bic) {
        return -2 * logl + Log (samples) * params;
    }
    return -2 * logl + 2 * samples / (samples - params - 1) * params;
}


