/**
    ---- OUTLINE ----
    1. SETUP
        1a. Initial Setup
        1b. User Input
        1c. Load and Filter Data
        1d. Check that there are enough sites to permit cAIC
        1e. Baseline fit on entire alignment
        1f. Set up some book keeping
    2. MAIN ANALYSIS
        2a. Evaluation of single break points with brute force
        2b. Evaluation of multiple break points with genetic algorithm
    3. POST PROCESSING

    FUNCTIONS
        - General GARD Functions
        - Genetic Algorithm Functions
        - Helper Functions
*/

/* 1. SETUP
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/* 1a. Initial Setup
------------------------------------------------------------------------------*/
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


gard.analysis_description = {terms.io.info : "GARD : Genetic Algorithms for Recombination Detection. Implements a heuristic
    approach to screening alignments of sequences for recombination, by using the CHC genetic algorithm to search for phylogenetic
    incongruence among different partitions of the data. The number of partitions is determined using a step-up procedure, while the
    placement of breakpoints is searched for with the GA. The best fitting model (based on c-AIC) is returned; and additional post-hoc
    tests run to distinguish topological incongruence from rate-variation.",
                           terms.io.version : "0.1",
                           terms.io.reference : "**Automated Phylogenetic Detection of Recombination Using a Genetic Algorithm**, _Mol Biol Evol 23(10), 1891–1901",
                           terms.io.authors : "Sergei L Kosakovsky Pond",
                           terms.io.contact : "spond@temple.edu",
                           terms.io.requirements : "A sequence alignment."
                          };

namespace terms.gard {
    nucleotide = "Nucleotide";
    protein    = "Protein";
    codon      = "Codon";
};


/* 1b. User Input
------------------------------------------------------------------------------*/
//TODO: io.DisplayAnalysisBanner (gard.analysis_description);

KeywordArgument ("type",        "The type of data to perform screening on", "Nucleotide");
KeywordArgument ("code",        "Genetic code to use (for codon alignments)", "Universal", "Choose Genetic Code");
KeywordArgument ("alignment",   "Sequence alignment to screen for recombination");

gard.data_type = io.SelectAnOption  ({terms.gard.nucleotide : "A nucleotide (DNA/RNA) alignment",
                                      terms.gard.protein : "A protein alignment",
                                      terms.gard.codon : "An in-frame codon alignment"},
                                      "The type of data to perform screening on");


/* 1c. Load and Filter Data
------------------------------------------------------------------------------*/
if (gard.data_type == terms.gard.nucleotide) {
    LoadFunctionLibrary ("libv3/models/DNA/GTR.bf");
    gard.model.generator = "models.DNA.GTR.ModelDescription";
    gard.alignment = alignments.ReadNucleotideDataSet ("gard.sequences", null);
    DataSetFilter gard.filter = CreateFilter (gard.sequences, 1);
} else {
    // TODO: implement these branches
    if (gard.data_type == terms.gard.protein) {
        gard.alignment = alignments.ReadProteinDataSet ("gard.sequences", null);
        DataSetFilter gard.filter = CreateFilter (gard.sequences, 1);
    } else {
        gard.alignment = alignments.LoadGeneticCodeAndAlignment ("gard.sequences", "gard.filter", null);
    }
}

// Define model to be used in each fit
gard.model = model.generic.DefineModel (gard.model.generator, "gard.overall_model", {"0" : "terms.global"}, "gard.filter", null);

gard.numSites               = gard.alignment[terms.data.sites];
gard.numSeqs                = gard.alignment[terms.data.sequences];

// Get a matrix of the variable sites {"0": site index, "1": site index ...}
gard.variableSiteMap = {};
utility.ForEach (alignments.Extract_site_patterns ("gard.filter"), "_pattern_", "
    if (_pattern_[terms.data.is_constant]==FALSE) {
        utility.ForEachPair (_pattern_[terms.data.sites], '_key_', '_value_',
        '
            gard.variableSiteMap + _value_;
        ');
    }
");
gard.variableSiteMap = Transpose (utility.DictToArray (gard.variableSiteMap)) % 0; // sort by 1st column
gard.variableSites = Rows (gard.variableSiteMap);
gard.numberOfPotentialBreakPoints = gard.variableSites - 1;

io.ReportProgressMessage ("", ">Loaded a `gard.data_type` multiple sequence alignment with **`gard.numSeqs`** sequences, **`gard.numSites`** sites (`gard.variableSites` of which are variable)\n \`" +
                                gard.alignment[utility.getGlobalValue("terms.data.file")] + "\`");

gard.baselineParameters     = (gard.model[terms.parameters])[terms.model.empirical] + // empirical parameters
                              utility.Array1D ((gard.model[terms.parameters])[terms.global]) + // global parameters
                              utility.Array1D ((gard.model[terms.parameters])[terms.local]) * (2*gard.numSeqs-5) * 2; // local parameters X number of branches x at least 2 partitions

gard.min_expected_sites     = (gard.baselineParameters + 2);
gard.minPartitionSize       = 2*gard.numSeqs-3;
/** gard.minPartitionSize:
    a simple heuristic that requires that there c-AIC makes sense for each individual
    partition i.e. there are at least 2 more sites than branch lengths
**/

io.ReportProgressMessage ("", ">Minimum size of a partition is set to be `gard.minPartitionSize` sites");

/* 1d. Check that there are enough sites to permit cAIC
------------------------------------------------------------------------------*/
io.CheckAssertion("gard.min_expected_sites <= gard.numSites", "The alignment is too short to permit c-AIC based model comparison. Need at least `gard.min_expected_sites` sites for `gard.numSeqs` sequences to fit a two-partiton model.");


/* 1e. Baseline fit on entire alignment
------------------------------------------------------------------------------*/
io.ReportProgressMessageMD('GARD', 'baseline-fit', 'Fitting the baseline (single-partition; no breakpoints) model');
// Infer NJ tree, estimting rate parameters (branch-lengths, frequencies and substitution rate matrix)

gard.base_likelihood_info = gard.fit_partitioned_model (null, gard.model, null);
gard.initial_values = gard.base_likelihood_info;
gard.globalParameterCount = utility.Array1D (estimators.fixSubsetOfEstimates(gard.initial_values, gard.initial_values[terms.global]));
(gard.initial_values [terms.branch_length])[0] = null;
(gard.initial_values [terms.branch_length])[1] = null;

gard.baseline_cAIC = math.GetIC (gard.base_likelihood_info[terms.fit.log_likelihood], gard.base_likelihood_info[terms.parameters], gard.numSites);
io.ReportProgressMessageMD("GARD", "baseline-fit", "* " + selection.io.report_fit (gard.base_likelihood_info, 0, gard.numSites));


/* 1f. Set up some book keeping
------------------------------------------------------------------------------*/
// Setup mpi queue
gard.createLikelihoodFunctionForExport ("gard.export", gard.model);
gard.queue = mpi.CreateQueue (
                            {
                            "LikelihoodFunctions" : {{"gard.export"}},
                            "Headers" : {{"libv3/all-terms.bf"}},
                            "Variables" : {{"gard.globalParameterCount", "gard.numSites"}}
                            }
                        );

// Setup master list of evaluated models.
gard.masterList = {};
/** gard.masterList:
    "model string": "model fitness (cAIC score)"
    model string is in the format: "{\n{bp1} \n{bp2} ...\n}"
    model fitness is either infinity (if not evaluated) or the numeric cAIC score
**/

gard.masterList [{{}}] = gard.baseline_cAIC;
gard.bestOverall_cAIC_soFar = gard.baseline_cAIC;


/* 2. MAIN ANALYSIS
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/* 2a. Evaluation of single break points with brute force
------------------------------------------------------------------------------*/
io.ReportProgressMessageMD('GARD', 'single-breakpoint', 'Performing an exhaustive single breakpoint analysis');

namespace gard {

    singleBreakPoint_best_cAIC = ^"math.Infinity";

    for (breakPointIndex = 0; breakPointIndex < variableSites - 1; breakPointIndex += 1) {
        siteIndex = variableSiteMap [breakPointIndex];

        io.ReportProgressBar ("GARD", "Breakpoint " +  (1+breakPointIndex) + " of " + (variableSites-1));


        if (gard.validatePartititon ({{siteIndex}}, minPartitionSize, numSites) == FALSE)  {
            continue;
        }

        mpi.QueueJob (queue, "gard.obtainModel_cAIC", {"0" : {{siteIndex__}},
                                                     "1" : model,
                                                     "2" : base_likelihood_info},
                                                     "gard.storeSingleBreakPointModelResults");
                                                    
    }

    mpi.QueueComplete (queue);

    io.ClearProgressBar();
    io.ReportProgressMessageMD('GARD', 'single-breakpoint', 'Done with single breakpoint analysis.');
    io.ReportProgressMessageMD('GARD', 'single-breakpoint', ("   Best sinlge break point location: " + singleBreakPoint_best_location));
    io.ReportProgressMessageMD('GARD', 'single-breakpoint', ("   c-AIC  = " + singleBreakPoint_best_cAIC));

}

// Evaluate if the analyis should continue to the multi break point stage
if (gard.singleBreakPoint_best_cAIC < gard.bestOverall_cAIC_soFar) {
    gard.bestOverall_cAIC_soFar = gard.singleBreakPoint_best_cAIC;
} else {
    gard.concludeAnalysis();
    return 0;
}

/* 2b. Evaluation of multiple break points with genetic algorithm
------------------------------------------------------------------------------*/
io.ReportProgressMessageMD('GARD', 'multi-breakpoint', 'Performing multi breakpoint analysis using a genetic algorithm');

namespace gard {
    
    addingBreakPointsImproves_cAIC = TRUE;
    numberOfBreakPointsBeingEvaluated = 1;
    while(addingBreakPointsImproves_cAIC) {
    
        numberOfBreakPointsBeingEvaluated+=1;
        populationSize = 20; // the GARD paper used: (numberOfMpiNodes*2 - 2) with 17 mpi nodes
        mutationRate = 0.75; // the GARD paper said "15% of randomly selected bits were toggled"...
        maxFailedAttemptsToMakeNewModel = 5;
        cAIC_diversityThreshold = 0.001;
        cAIC_improvementThreshold = 0.1; // I think this was basically 0 in the gard paper
        maxGenerationsAllowedWithNoNewModelsAdded = 10; // Not in the GARD paper. use 10?
        maxGenerationsAllowedAtStagnent_cAIC = 50; // Set to 100 in the GARD paper

        parentModels = gard.GA.initializeModels(numberOfBreakPointsBeingEvaluated, populationSize, numberOfPotentialBreakPoints);
        
        terminationCondition = FALSE;
        while(terminationCondition == FALSE) {
            io.ReportProgressBar ("GARD", 'Genetic algorithm for ' + numberOfBreakPointsBeingEvaluated + ' breakpoint analysis: generation ' +  (1+generation));

            // 1. Produce the next generation of models with recombination. 
            childModels = gard.GA.recombineModels(parentModels, populationSize, differenceThreshold);
            
            // 2. Select the fittest models. 
            interGenerationalModels = parentModels;
            interGenerationalModels * childModels;
            gard.GA.evaluateModels(interGenerationalModels, model, null, numSites);
            selectedModels = gard.GA.selectModels(interGenerationalModels, populationSize);

            // 3. Evaluate convergence for this modelSet initialization 
            // If converged, produce a new parent modelSet (keeping the current fitest)
            if (gard.GA.modelSetsAreTheSame(selectedModels, parentModels)) {
                generationsNoNewModelsAdded += 1;
            } else {
                generationsNoNewModelsAdded = 0;
            }
            if ( (math.minNormalizedRange(selectedModels) < cAIC_diversityThreshold) || (generationsNoNewModelsAdded > maxGenerationsAllowedWithNoNewModelsAdded) ) {
                parentModels = gard.GA.generateNewGenerationOfModelsByMutatingModelSet(parentModels, numberOfPotentialBreakPoints, mutationRate);
                if (Abs(parentModels) == 1) {
                    parentModels = gard.GA.initializeModels(numberOfBreakPointsBeingEvaluated, populationSize, numberOfPotentialBreakPoints);
                }
                differenceThreshold = numberOfBreakPointsBeingEvaluated / 4;
            } else {
                parentModels = selectedModels;
            }

            // 4. Evaluate convergence for this number of break points
            // If converged, move on to n+1 break points or end analysis
            currentBest_cAIC = Min(utility.Values(interGenerationalModels),0);
            if (previousBest_cAIC - currentBest_cAIC < cAIC_improvementThreshold) {
                generationsAtCurrentBest_cAIC += 1;
                if (generationsAtCurrentBest_cAIC >= maxGenerationsAllowedAtStagnent_cAIC) {
                    terminationCondition = TRUE;
                }
            } else {
                generationsAtCurrentBest_cAIC = 0;
            }
            previousBest_cAIC = currentBest_cAIC;

            generation += 1;
        }

        io.ClearProgressBar();
        // Report status of n-break point analysis.
        statusString = ("Done with " + numberOfBreakPointsBeingEvaluated + " breakpoint analysis.\n" + gard.GA.getMultiBreakPointStatusString(selectedModels));
        io.ReportProgressMessageMD('GARD', 'multi-breakpoint', statusString);

        // Evaluate if the analysis should conclude or if the GA should run with breakpoints + 1
        if (previousBest_cAIC >= bestOverall_cAIC_soFar) {
            addingBreakPointsImproves_cAIC = FALSE;
            gard.concludeAnalysis();
        }
    }

}


/* 3. POST PROCESSING
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

//--------------------------------------------------------------------------------------------------------------------

/* FUNCTIONS
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

// ---- General GARD Functions ----

/**
 * @name gard.fit_partitioned_model
 * Given a list of partitions, specified as increasing breakpoint locations,
   fit the specified model to said partitions, using neighbor joining trees on each partition
   return LogL and IC values
   // TODO: currently just returning the log_likelihood info and then calculating the cAIC from that...
   // If we were to change the function to execute how it is documented we wouldn't have access to the additional info (trees, parameteres, etc.)
   // Not sure if we will need that moving forward but we can revisit this.

 * @param {Matrix} breakPoints : sorted, 0-based breakpoints, e.g.
    {{100,200}} -> 3 partitions : 0-100, 101-200, 201-end
 * @param {Dict} model : an instantiated model to be used for all partitions
 * @param {Dict/null} initial_values : if provided, use as initial values

 * @returns a {Dictionary} :
    terms.fit.log_likelihood -> log likelihood
    terms.fit.AICc -> small sample AIC

 */

lfunction gard.fit_partitioned_model (breakPoints, model, initial_values) {

    current_index = 0;
    current_start = 0;
    breakPoints_count = utility.Array1D (breakPoints);
    part_count = breakPoints_count + 1;
    lf_components = {2 * (part_count), 1};
    trees = {};

    for (p = 0; p < part_count; p += 1) {
        last_partition = p >= breakPoints_count;
        if (!last_partition) {
            current_end = breakPoints[p];
        } else {
            current_end = ^"gard.filter.sites" - 1;
        }
        lf_components [2*p] = "gard.filter.part_" + p;
        lf_components [2*p+1] = "gard.tree.part_" + p;
        DataSetFilter ^(lf_components[2*p]) = CreateFilter (^"gard.filter", 1, "" + current_start + "-" + current_end);
        trees[p] = trees.ExtractTreeInfo (tree.infer.NJ ( lf_components[2*p], null));
        model.ApplyModelToTree(lf_components[2 * p + 1], trees[p], {
            "default": model
        }, None);
        // Increment the current starting point
        if (!last_partition) {
            current_start = breakPoints[p];
        }
    }

    lf_id = &likelihoodFunction;
    utility.ExecuteInGlobalNamespace ("LikelihoodFunction `lf_id` = (`&lf_components`)");

    model_objects = {
        "gard.overall_model" : model
    };

    df = 0;
    if (Type(initial_values) == "AssociativeList") {
        utility.ToggleEnvVariable("USE_LAST_RESULTS", 1);
        df = estimators.ApplyExistingEstimates(&likelihoodFunction, model_objects, initial_values, None);
    }

    res = estimators.FitExistingLF (&likelihoodFunction, model_objects);

    DeleteObject (likelihoodFunction, :shallow);

    res[^"terms.parameters"] += df + (model[^"terms.parameters"])[^"terms.model.empirical"];
    return res;

}

/**
 * @name gard.obtainModel_cAIC
 * Wrap a call to gard.fit_partitioned_model to compute c-AIC

   @returns c-AIC

 */
lfunction gard.obtainModel_cAIC (breakPoints, model, initial_values) {
    /*if (utility.Has (^"gard.masterList", breakPoints, "Number")) {
        return (^"gard.masterList")[breakPoints];
    }*/
    fit = gard.fit_partitioned_model (breakPoints, model, initial_values);
    c_aic = math.GetIC (fit[^"terms.fit.log_likelihood"], fit[^"terms.parameters"] + ^"gard.globalParameterCount", ^"gard.numSites");
    //(^"gard.masterList")[breakPoints] = c_aic;
    return c_aic;
}

function gard.storeSingleBreakPointModelResults (node, result, arguments) {
    result = Eval (result);
    gard.masterList[arguments[0]] = result;
    if (result < gard.singleBreakPoint_best_cAIC) {
        gard.singleBreakPoint_best_cAIC = result;
        gard.singleBreakPoint_best_location = Eval(arguments[0])[0];
    }
}

lfunction gard.createLikelihoodFunctionForExport (id,  model) {
    overall_tree = trees.ExtractTreeInfo (tree.infer.NJ ( "gard.filter", null));
    model.ApplyModelToTree("gard._ignore_tree",overall_tree, {
        "default": model
    }, None);

    utility.ExecuteInGlobalNamespace ("LikelihoodFunction `id` = (gard.filter, gard._ignore_tree)");

}

/**
    definition is assumed to contain a SORTED list of 0-based
    breakpoint locations

    If any of the individual partitions is < min_size, the function
    returns FALSE, otherwise it returns TRUE
*/
lfunction gard.validatePartititon (definition, min_size, total_sites) {
    lastBP = 0;
    bpCount = utility.Array1D (definition);
    for (i = 0; i < bpCount; i+=1) {
        if (definition[i] - lastBP + 1 < min_size) {
            return FALSE;
        }
        lastBP = definition[i];
    }
    if (total_sites - lastBP < min_size) {
        return FALSE;
    }
    return TRUE;
}

/**
   @name gard.modelIsNotInMasterList
 

   @params {Matrix || Dictonary} breakPoints
   @returns {Bolean}

 */
lfunction gard.modelIsNotInMasterList(masterList, breakPoints) {
    if (Type(breakPoints) == "Matrix") {
        breakPointString = '';
        Eval("breakPoints = breakPointString + breakPoints"); 
    }
    if (utility.KeyExists(masterList, breakPoints)) {
        return FALSE;
    } else {
        return TRUE;
    }
}

lfunction gard.concludeAnalysis() {
    console.log('The cAIC score was not improved by adding additional break points; the anlysis is complete.');
    console.log('TODO: any book keeping to close the analysis can go here');
}


// ---- Genetic Algorithm (GA) Functions ----

/**
 * @name gard.GA.initializeModels
 * @param {Number} numberOfBreakPoints
 * @param {Number} populationSize
 * @param {Number} numberOfPotentialBreakPoints
 * @returns a {Dictonary} initializedModels
 */
function gard.GA.initializeModels (numberOfBreakPoints, populationSize, numberOfPotentialBreakPoints) {

    initializedModels = {};
    for (modelNumber=0; modelNumber < populationSize; modelNumber += 1) {
        // TODO : make sure that feasible solutions exist, i.e. that gard.validatePartititon doesn't just keep rejecting all models
        do {
            breakPoints = ({1,numberOfBreakPoints} ["Random(0, numberOfPotentialBreakPoints)$1"]) % 0 ;
        } while (gard.validatePartititon (breakPoints, gard.minPartitionSize, gard.numSites) == FALSE);

        initializedModels[breakPoints] = math.Infinity;
    }
    return initializedModels;
}

/**
 * @name gard.GA.recombineModels
 Given a set of models create a new generation of models by iteratively recombining two random parent models.
 The child model will have a random subset of the breakpoints from the parents

 * @param {Matrix} parentModels
 * @param {Number} populationSize
 * @returns a {Dictonary} childModels
 */
function gard.GA.recombineModels (parentModels, populationSize, differenceThreshold) {
    parentModelIds = utility.Keys(parentModels);
    numberOfParentModels = Columns(parentModelIds);
    firstModel = gard.Helper.convertMatrixStringToMatrix(parentModelIds[0]);
    numberOfBreakPoints = Columns(firstModel);

    childModels = {};
    for(modelIndex=0; modelIndex<populationSize; modelIndex+=1) {

        modelIsValid = FALSE;
        failedAttempts = 0;
        while(modelIsValid == FALSE && failedAttempts < gard.maxFailedAttemptsToMakeNewModel) {

            parentModel1 = gard.Helper.convertMatrixStringToMatrix(parentModelIds[Random(0, numberOfParentModels-.0001)$1]);
            parentModel2 = gard.Helper.convertMatrixStringToMatrix(parentModelIds[Random(0, numberOfParentModels-.0001)$1]);

            while(gard.GA.modelsDifferentEnough(parentModel1, parentModel2, differenceThreshold) == FALSE) {
                parentModel2 = gard.Helper.convertMatrixStringToMatrix(parentModelIds[Random(0, numberOfParentModels-.0001)$1]);
            }

            breakPoints = {1, numberOfBreakPoints};
            for (breakPointNumber=0; breakPointNumber < numberOfBreakPoints; breakPointNumber=breakPointNumber+1) {
                if (random.TRUE_or_FALSE()) {
                    parentModel = parentModel1;
                } else {
                    parentModel = parentModel2;
                }
                breakPoint = parentModel[Random(0, numberOfBreakPoints-.0001)$1];
                
                breakPoints[breakPointNumber] = breakPoint;
            }

            breakPoints = gard.Helper.sortedMatrix(breakPoints);
            if ( (gard.modelIsNotInMasterList(gard.masterList, breakPoints)) && (gard.validatePartititon(breakPoints, gard.minPartitionSize, gard.numSites)) ) { // TODO: also check that no partition is too small
                modelIsValid = TRUE;
            } else {
                failedAttempts+=1;
            }

        }
        if (modelIsValid) {
            childModels[breakPoints] = math.Infinity;
        }

    }

    return childModels;
}

/**
 * @name gard.GA.evaluateModels
 * @param {Dictonary} models (some models may have cAIC scores already)
 * @param {Model} evolutionaryModel; to be used in fit_partitioned_model (i.e. gard.model)
 * @param {??} initial_values; The initial values for fit_partitioned_model (can be null)
 * @param {Number} numSites; the number of sites in the alignment, to be used in the cAIC calculation
 * @returns nothing... the gard.GA.storeMultiBreakPointModelResults function gets called which updates the intergenerationalModels object
 */
function gard.GA.evaluateModels (models, evolutionaryModel, initialValues, numSites) {
    modelIds = utility.Keys(models);
    numberOfModels = Columns(modelIds);

    for(modelIndex=0; modelIndex<numberOfModels; modelIndex=modelIndex+1) {
        modelId = modelIds[modelIndex];
        cAIC = models[modelId];

        if (cAIC == math.Infinity) {
            breakPoints = gard.Helper.convertMatrixStringToMatrix(modelId);
            mpi.QueueJob (queue, "gard.obtainModel_cAIC", {"0" : breakPoints__,
                                                     "1" : evolutionaryModel,
                                                     "2" : base_likelihood_info},
                                                     "gard.GA.storeMultiBreakPointModelResults");
        }
        mpi.QueueComplete (queue);

    }

}

/**
 * @name gard.GA.selectModels
 * @param {Dictionary} evaluatedModels
 * @param {Number} numberOfModelsToKeep
 * @returns a {Matrix} selectedModels
 */
function gard.GA.selectModels (evaluatedModels, numberOfModelsToKeep) {
    modelIds = utility.Keys(evaluatedModels);
    cAIC_scores = utility.Values(evaluatedModels);

    selectedModels = {};
    for(i=0; i<numberOfModelsToKeep; i+=1) {
        lowest_cAIC_index = Min(cAIC_scores,1)[1];
        selectedModelId = modelIds[lowest_cAIC_index];
        selectedModels[selectedModelId] = cAIC_scores[lowest_cAIC_index];

        cAIC_scores[lowest_cAIC_index] = math.Infinity;
    }

    return selectedModels;
}

/**
 * @name gard.GA.modelsDifferentEnough
 * @param {Matrix} model1
 * @param {Matrix} model2
 * @param {Number} differenceThreshold
 * @returns a {Bolean} modelsDifferentEnough
 */
function gard.GA.modelsDifferentEnough (model1, model2, differenceThreshold) {
    // TODO: Actually compare the models, not just make sure they aren't identical
    if (model1 == model2) {
        modelsDifferentEnough = 0;
    } else {
        modelsDifferentEnough = 1;
    }
    return modelsDifferentEnough;
}

/**
 * @name gard.GA.modelSetsAreTheSame
 * @param {Dictionary} modelSet1
 * @param {Dictionary} modelSet2
 * @returns a {Bolean} 
 */
function gard.GA.modelSetsAreTheSame(modelSet1, modelSet2) {
    if(Abs(modelSet1) != Abs(modelSet2)){
        return FALSE;
    }

    IdSet1 = utility.Keys(modelSet1);
    for(i=0; i<Abs(modelSet1); i=i+1) {
        if(utility.KeyExists(modelSet2, IdSet1[i]) == FALSE) {
            return FALSE;
        }
    }

    return TRUE;
}

/**
 * @name gard.GA.generateNewGenerationOfModelsByMutatingModelSet
 Keeps the current best performing model
 Generates new models for the rest of the population by randomly mutating some of the parent models break points
 // TODO: decide if we want some logic in here to more meaningfully mutate.
 * @param {Dictonary} parentModels
 * @param {Number} numberOfPotentialBreakPoints
 * @returns a {Dictonary} the new generation of models
 */
function gard.GA.generateNewGenerationOfModelsByMutatingModelSet(parentModels, numberOfPotentialBreakPoints, mutationRate) {
    modelIds = utility.Keys(parentModels);
    populationSize = Columns(modelIds);
    firstModel = gard.Helper.convertMatrixStringToMatrix(modelIds[0]);
    numberOfBreakPoints = Columns(firstModel);

    // Generate a new set of models
    nextGenOfModels = {};
    for(i=0; i<populationSize-1; i=i+1) {
        
        modelIsValid = FALSE;
        failedAttempts = 0;
        while(modelIsValid == FALSE && failedAttempts < gard.maxFailedAttemptsToMakeNewModel) {
            parentModel = gard.Helper.convertMatrixStringToMatrix(modelIds[i]);
            breakPoints = {1,numberOfBreakPoints};
            for(breakPointIndex=0; breakPointIndex<numberOfBreakPoints; breakPointIndex=breakPointIndex+1) {
                if(Random(0,1) < mutationRate) {
                    breakPoints[breakPointIndex] = parentModel[breakPointIndex];
                } else {
                    breakPoints[breakPointIndex] = Random(0,numberOfPotentialBreakPoints)$1;
                }
            }
            if ((gard.modelIsNotInMasterList(gard.masterList, breakPoints)) && (gard.validatePartititon(breakPoints, gard.minPartitionSize, gard.numSites)) ) { // TODO: also check that no partition is too small
                modelIsValid = TRUE;
            } else {
                failedAttempts+=1;
            }

        }

        if (modelIsValid) {
            nextGenOfModels[breakPoints] = math.Infinity;
        }
    }

    // Add the most fit model from the previous set
    cAIC_scores = utility.Values(parentModels);
    lowest_cAIC_index = Min(cAIC_scores,1)[1];
    selectedModelId = modelIds[lowest_cAIC_index];
    nextGenOfModels[selectedModelId] = cAIC_scores[lowest_cAIC_index];

    return nextGenOfModels;
}

/**
 * @name gard.GA.getMultiBreakPointStatusString
 get the string to print out after the GA has finished n-breakpoint analysis
 * @param {Dictonary} evaluatedModels: the set of models at the end of the analysis
 * @returns a {String} the markdown string to log to console
 */
function gard.GA.getMultiBreakPointStatusString(evaluatedModels) {
    best_cAIC_index = Min(utility.Values(evaluatedModels),1)[1];
    bestModel = gard.Helper.convertMatrixStringToMatrix(utility.Keys(evaluatedModels)[best_cAIC_index]);
    best_cAIC = Min(utility.Values(evaluatedModels),1)[0];
    bestBreakPoints={1,Columns(bestModel)};
    bestBreakPointsString = '';
    for(i=0; i<Columns(bestModel); i=i+1) {
        bestBreakPoints[i] = gard.variableSiteMap[bestModel[i]];
        Eval("bestBreakPointsString += gard.variableSiteMap[bestModel[i]]");
        if (i<Columns(bestModel)-1) {
            bestBreakPointsString += ', ';
        }
    }
    statusString = "    Best break point locations: " + bestBreakPointsString + "\n" +
                   "    c-AIC = " + best_cAIC;
    
    return statusString;
}

/**
 * @name gard.GA.storeMultiBreakPointModelResults
 // The call back function passed to mpi.QueueJob in the gard.GA.evaluateModels function
 * @param {Dictonary} evaluatedModels: the set of models at the end of the analysis
 * @returns a {String} the markdown string to log to console
 */
function gard.GA.storeMultiBreakPointModelResults(node, result, arguments) {
    result = Eval (result);
    gard.masterList[arguments[0]] = result;
    gard.interGenerationalModels[arguments[0]] = result;
}


// ---- Helper Functions ----
// TODO: replace these helper functions with a more HBL way of dealing with matrices

// Return the nth matrix in a nxm matrix.
function gard.Helper.getNthMatrix(setOfMatrices, n) {
    nthMatrix = {1, Columns(setOfMatrices)};
    for(i=0; i<Columns(setOfMatrices); i=i+1) {
        nthMatrix[i] = setOfMatrices[n][i];
    }
    return nthMatrix;
}

function gard.Helper.convertMatrixStringToMatrix(matrixString){
    stringToEval = 'matrix = ' + matrixString;
    Eval(stringToEval);
    return matrix;
}

//sort numeric 1xn matrix into assending order.
function gard.Helper.sortedMatrix(matrix) {
    sortedMatrix = {1,Columns(matrix)};

    for (i=0; i<Columns(matrix); i=i+1) {
        minElementLeft=Min(matrix, 1);
        sortedMatrix[i] = minElementLeft[0];
        matrix[minElementLeft[1]] = math.Infinity;
    }

    return sortedMatrix;
}
