/**
    ---- OUTLINE ----
    1. SETUP
        1a. Initial Setup
        1b. User Input
        1c. Load and Filter Data
        1d. Baseline fit on entire alignment
        1e. Checks
        1f. Set up some book keeping
    2. MAIN ANALYSIS
        2a. Evaluation of single break points with brute force
        2b. Evaluation of multiple break points with genetic algorithm
    3. POST PROCESSING

    GARD FUNCTIONS
*/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    1. SETUP
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/*------------------------------------------------------------------------------
    1a. Initial Setup
*/
LoadFunctionLibrary ("libv3/tasks/trees.bf");
LoadFunctionLibrary ("libv3/tasks/alignments.bf");
LoadFunctionLibrary ("libv3/tasks/estimators.bf");
LoadFunctionLibrary ("libv3/convenience/regexp.bf");
LoadFunctionLibrary ("libv3/convenience/math.bf");
LoadFunctionLibrary ("libv3/IOFunctions.bf");
LoadFunctionLibrary ("libv3/UtilityFunctions.bf");


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

namespace gard.parameters {

};

gard.master_list = {}; // "model string" -> "model fitness"


/*------------------------------------------------------------------------------
    1b. User Input
*/
// TODO: turn this back on... io.DisplayAnalysisBanner (gard.analysis_description);

KeywordArgument ("type",        "The type of data to perform screening on", "Nucleotide");
KeywordArgument ("code",        "Genetic code to use (for codon alignments)", "Universal", "Choose Genetic Code");
KeywordArgument ("alignment",   "Sequence alignment to screen for recombination");

gard.data_type = io.SelectAnOption  ({terms.gard.nucleotide : "A nucleotide (DNA/RNA) alignment",
                                      terms.gard.protein : "A protein alignment",
                                      terms.gard.codon : "An in-frame codon alignment"},
                                      "The type of data to perform screening on");


/*------------------------------------------------------------------------------
    1c. Load and Filter Data
*/
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


/*------------------------------------------------------------------------------
    1d. Baseline fit on entire alignment
*/
console.log ( "\n> Fitting a baseline model...\n");
// Define model to be used in each fit
gard.model = model.generic.DefineModel (gard.model.generator, "gard.overall_model", {"0" : "terms.global"}, "gard.filter", null); 

// Infer NJ tree, estimting rate parameters (branch-lengths, frequencies and substitution rate matrix)
base_likelihood_info = gard.fit_partitioned_model (null, gard.model, null);

// Calculate c-AIC
// TODO: confirm that the calculation of c-AIC is correct.
base_logLikelihood = base_likelihood_info["LogL"];
baseParams = base_likelihood_info["parameters"];
numSites = gard.alignment['sites'];
numSeqs = gard.alignment['sequences'];
base_cAIC = gard.calculate_cAIC(base_logLikelihood, baseParams, numSites);

console.log("Done with single partition analysis.");
console.log("   Log(L) = " + base_logLikelihood);
console.log("   c-AIC  = " + base_cAIC);


/*------------------------------------------------------------------------------
    1e. Checks
*/
// Too few sites for c-AIC inference?.
// from paper: "Note that the use of AICc sensibly requires that there be more observations (alignment columns) 
// than the number of estimated model parameters. The formal requirement for this setting is, consequently, 
// N > 2(2S − 3) + bp"
// TODO: Confirm that the check below is correct. It seems different than the old implementation.
if (numSites < 2 * (2 * numSeqs - 2) + baseParams) {
    console.log("ERROR: Too few sites for c-AIC inference.\n");
    return 0;
}


/*------------------------------------------------------------------------------
    1f. Set up some book keeping
*/
minimumNumberOfVariableSitesInPartitoin = 2; // TODO: what should this number be? and should it be used at all? (not yet used in the GA stage)

// Get a list of the variable/informative sites (this will be the list of potential breakpoints).
variableSitesList = gard.getVariableSitesMatrix();
numberOfPotentialBreakPoints = Abs(variableSitesList);
bppSize = (Log(numberOfPotentialBreakPoints)/Log(2)+1)$1; // number of binary digits to represent a potential break point.
console.log("\nThere are " + numberOfPotentialBreakPoints + " potential breakpoints. Bit size of the sample is " + bppSize + '\n');


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    2. MAIN ANALYSIS
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/*------------------------------------------------------------------------------
    2a. Evaluation of single break points with brute force
*/
console.log ( "\n> Fitting a two partition model...\n");
// TODO: Improve this section
// 
// EVALUATION WITH Flu.fasta dataset:
//                          2.3.11          2.4.0-alpha
//  No BP cAIC              8887.4          8878.7              GOOD   
//  Best single BP          696             696                 GOOD
//  Best single BP cAIC     8764.7          8603                NOT-GOOD

// Don't evaluate break points near the begining or end of the alignment that create partitions that are two small
firstSinglePartitionBreakPointIndex = minimumNumberOfVariableSitesInPartitoin - 1;
lastSingelPartitionBreakPointIndex = numberOfPotentialBreakPoints - minimumNumberOfVariableSitesInPartitoin;

singleBreakPoint_best_cAIC = 10^25;
for (breakPointIndex = firstSinglePartitionBreakPointIndex; breakPointIndex < 5 /*lastSingelPartitionBreakPointIndex*/; breakPointIndex=breakPointIndex+1) {
    siteIndex = variableSitesList[breakPointIndex];
    singleBreakPoint_LikelihoodInfo = gard.fit_partitioned_model ({{siteIndex}}, gard.model, null);
    singleBreakPoint_logL = singleBreakPoint_LikelihoodInfo["LogL"];
    singleBreakPoint_numberParams = singleBreakPoint_LikelihoodInfo["parameters"];
    singleBreakPoint_cAIC = gard.calculate_cAIC(singleBreakPoint_logL, singeleBreakPoint_numberParams, numSites);

    if (singleBreakPoint_cAIC < singleBreakPoint_best_cAIC) {
        singleBreakPoint_best_cAIC = singleBreakPoint_cAIC;
        singleBreakPoint_best_location = siteIndex;
        singleBreakPoint_best_LogL = singleBreakPoint_logL;
    }
}

console.log("Done with two partition analysis.");
console.log("   Best sinlge break point location: " + singleBreakPoint_best_location);
console.log("   Log(L) = " + singleBreakPoint_best_LogL);
console.log("   c-AIC  = " + singleBreakPoint_best_cAIC);


/*------------------------------------------------------------------------------
    2b. Evaluation of multiple break points with genetic algorithm
*/
/** Notes on implementing the genetic algorithm

    ---- Pseudocode from CHC paper ----
    Code:
        t =0;
        d = L/4;
        initialize P(t);
        evaluate structures in P(t);
        while termination condition not satisfied do
        begin
            t = t+1;
            select_r C(t) from P(t-1);
            recombine structures in C(t) forming C'(t);
            evaluate structures in C'(t);
            select_s P(t) from C'(t) and P(t-1);
            if P(t) equals P(t-1)
                d--;
            if d < 0
            begin
                diverge P(t);
                d = r * (1.0 -r) * L;
            end
        end

    Procedures:

        select_r
            copy all members of P(t-1) to C(t) in random order;

        select_s
            form P(t) from P(t-1)
                by replacing the worst members of P(t-1) with the best members of C'(t)
                until no remaning member of C'(t) is any better than any remaining member of P(t-1);

        recombine
            for each of the M/2 pairs of structures in C(t) 
            begin
                determine the Hamming_distance
                if (Hamming_distance/2) > d
                    swap half the differing bits at random;
                else
                    delete the pair of structures from C(t)
            end

        diverge
            replace P(t) with M copies of the best member of P(t-1);
            for all but one member of P(t)
            begin
                flip r * L bits at random;
                evaluate structure;
            end

    Variables:
        M: populationSize
        L: string length
        t: generation
        d: difference threshold
        r: divergence rate

*/


console.log ( "\n> [Work in progress] fitting multi break point models with a genetic algorithm...\n");
numberOfBreakPointsBeingEvaluated = 2; // TODO: Iterate over increasing number of breakPoints until converged
populationSize = 4; // TODO: change to something like 50... set low for debuging
generation = 0;
differenceThreshold = numberOfBreakPointsBeingEvaluated / 4;

parentModels = gard.GA.initializeModels(numberOfBreakPointsBeingEvaluated, populationSize, numberOfPotentialBreakPoints);
terminationCondition = 0;
while(terminationCondition == 0) {
    childModels = gard.GA.recombineModels(parentModels, populationSize, differenceThreshold);
    interGenerationalModels = parentModels;
    interGenerationalModels * childModels;
    evaluatedModels = gard.GA.evaluateModels(interGenerationalModels, gard.model, null, numSites);
    selectedModels = gard.GA.selectModels(evaluatedModels, populationSize, numberOfBreakPointsBeingEvaluated);
    console.log("\n\n --- LEFT OFF HERE ----\n\n");
    return;

    if (gard.GA.modelSetsAreTheSame(selectedModels, parentModels)){
        differenceThreshold = gard.GA.decrementDifferenceThreshold(differenceThreshold);
    }

    if (differenceThreshold < 0) {
        parentModels = gard.GA.generateNewGenerationOfModelsByMutatingModelSet(parentModels);
        differenceThreshold = gard.GA.reinitializeDifferenceThreshold(numberOfBreakPointsBeingEvaluated);
    }

    terminationCondition = gard.GA.evaluateConvergence();
    generation = generation + 1;

}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    3. POST PROCESSING
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

//--------------------------------------------------------------------------------------------------------------------

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    GARD FUNCTIONS
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

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

    return estimators.FitExistingLF (&likelihoodFunction, model_objects);
    
}


/**
 * @name gard.calculate_cAIC   
 * @param {Number} logLikelihood
 * @param {Number} numParameters
 * @param {Number} numObservations (will be the number of sites)
 * @returns a {Number} the AICc score
 */
lfunction gard.calculate_cAIC (logLikelihood, numParameters, numObservations) {
    AIC = 2*numParameters - 2*logLikelihood;
    correction = ( 2*numParameters*numParameters + 2*numParameters ) / ( numObservations - numParameters - 1 );
    return AIC + correction;
}


/**
 * @name gard.getVariableSitesMatrix   
 * @param {DataSetFilter} 
 TODO: The old implementation had a zero as the last argument for HarvestFrequencies but this caused issues at sites that had gaps where it would
        only read the frequencies above the first gap... changing the last argument to a 1 seems to have fixed this but need to confirm.
 
 Also TODO-minor: I can't pass in the datasetfilter (I tired dereferencing, executing in global namespace, passing in as string, etc. but didn't work)
        There's the same issue with gard.fit_partitioned_model...
 * @returns a {Associative List} The variable sites in the form of {"0": siteIndex, "1": siteIndex, "2": siteIndex...}
 */
function gard.getVariableSitesMatrix () {
    breakPoints = {};
    breakPointNumber = 0;
    for (siteIndex=0; siteIndex<gard.filter.sites; siteIndex=siteIndex+1) {
        filterString =Format(siteIndex,20,0);
        DataSetFilter siteFilter = CreateFilter (gard.filter, 1, filterString);
        HarvestFrequencies (siteCharacterFrequencies, siteFilter, 1, 1, 1);
        
        maxFrequency = Abs(siteCharacterFrequencies);
        variableSite = maxFrequency<1;
        if (variableSite){
            breakPoints[breakPointNumber] = siteIndex;
            breakPointNumber = breakPointNumber+1;
        }

    }
    return breakPoints;
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
    for (modelNumber=0; modelNumber < populationSize; modelNumber=modelNumber+1) {
        breakPoints = {1,numberOfBreakPoints};
        for(breakPointIndex=0; breakPointIndex<numberOfBreakPoints; breakPointIndex=breakPointIndex+1) {
            breakPoints[breakPointIndex] = Random(0, numberOfPotentialBreakPoints)$1;
        }
        breakPoints = gard.Helper.sortedMatrix(breakPoints);
        modelId = Join(",", breakPoints);
        initializedModels[modelId] = {'breakPoints': breakPoints, 'cAIC': 'notEvaluated'};
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
    firstModelObject = parentModels[parentModelIds[0]];
    numberOfBreakPoints = Columns(firstModelObject['breakPoints']);

    childModels = {};
    for(modelIndex=0; modelIndex<populationSize; modelIndex=modelIndex+1) {
        parentModel1 = parentModels[parentModelIds[Random(0, populationSize)$1]];
        parentModel2 = parentModels[parentModelIds[Random(0, populationSize)$1]];

        while(gard.GA.modelsDifferentEnough(parentModel1, parentModel2, differenceThreshold) == 0) {
            parentModel2 = parentModels[parentModelIds[Random(0, populationSize)$1]];
        }

        breakPoints = {1, numberOfBreakPoints};
        for (breakPointNumber=0; breakPointNumber < numberOfBreakPoints; breakPointNumber=breakPointNumber+1) {
            // TODO-Major: actually recombine... I'm just randomly generating new modesl for now.
            breakPoint = Random(0, numberOfPotentialBreakPoints)$1;
            breakPoints[breakPointNumber] = breakPoint;
        }
        breakPoints = gard.Helper.sortedMatrix(breakPoints);

        modelId = Join(",", breakPoints);
        childModels[modelId] = {'breakPoints': breakPoints, "cAIC": 'notEvaluated'};

    }

    return childModels;
}


/**
 * @name gard.GA.evaluateModels
 * @param {Dictonary} models (some models may have cAIC scores already)
 * @param {Model} evolutionaryModel; to be used in fit_partitioned_model (i.e. gard.model)
 * @param {??} initial_values; The initial values for fit_partitioned_model (can be null)
 * @param {Number} numSites; the number of sites in the alignment, to be used in the cAIC calculation
 * @returns a {Dictionary} models (with cAIC for each model)
 */
function gard.GA.evaluateModels (models, evolutionaryModel, initial_values, numSites) {
    modelIds = utility.Keys(models);
    numberOfModels = Columns(modelIds);
    firstModelObject = models[modelIds[0]];
    numberOfBreakPoints = Columns(firstModelObject['breakPoints']);

    for(modelIndex=0; modelIndex<numberOfModels; modelIndex=modelIndex+1) {
        modelId = modelIds[modelIndex];
        modelObject = models[modelId];
        cAIC = modelObject['cAIC'];

        if (cAIC == 'notEvaluated') {
            breakPoints = modelObject['breakPoints'];
            modelLikelihoodInfo = gard.fit_partitioned_model(breakPoints, evolutionaryModel, initial_values);
            model_cAIC = gard.calculate_cAIC(modelLikelihoodInfo["LogL"], modelLikelihoodInfo["parameters"], numSites);
            
            models[modelId] = {'breakPoints': breakPoints, 'cAIC': model_cAIC};
        }

    }

    return models;


    /*
    numberOfModels = Rows(models);
    numberOfBreakPoints = Columns(models);

    modelPerformance = {};
    for(modelIndex=0; modelIndex<numberOfModels; modelIndex=modelIndex+1) {
        model = gard.Helper.getNthMatrix(models, modelIndex);
        individualModelLikelihoodInfo = gard.fit_partitioned_model(model, evolutionaryModel, initial_values);
        individualModel_cAIC = gard.calculate_cAIC(individualModelLikelihoodInfo["LogL"], individualModelLikelihoodInfo["parameters"], numSites);
        
        modelId = Join(",", model);
        modelPerformance[modelId] = individualModel_cAIC;
    }

    return modelPerformance;
    */

}


/**
 * @name gard.GA.selectModels
 * @param {Dictionary} evaluatedModels
 * @param {Number} numberOfModelsToKeep
 * @returns a {Matrix} selectedModels
 */
function gard.GA.selectModels (evaluatedModels, numberOfModelsToKeep, numberOfBreakPoints) {
    modelIds = utility.Keys(evaluatedModels);
    numberOfModels = Columns(modelIds);
    firstModelObject = models[modelIds[0]];
    numberOfBreakPoints = Columns(firstModelObject['breakPoints']);
    veryLargeNumber = 100000000;
    
    cAIC_scores = {1,numberOfModels};
    for(modelIndex=0; modelIndex<numberOfModels; modelIndex=modelIndex+1) {
        modelObject = evaluatedModels[modelIds[modelIndex]];
        cAIC_scores[modelIndex] = modelObject["cAIC"];    
    }

    selectedModels = {};
    for(i=0; i<numberOfModelsToKeep; i=i+1) {
        lowest_cAIC_index = Min(cAIC_scores,1)[1];

        selectedModelId = modelIds[lowest_cAIC_index];
        selectedModelObject = evaluatedModels[selectedModelId];
        selectedModels[selectedModelId] = selectedModelObject;

        cAIC_scores[lowest_cAIC_index] = veryLargeNumber;   
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


// ---- HELPER FUNCTIONS ----
// TODO: replace these helper functions with a more HBL way of dealing with matrices

//sort numeric 1xn matrix into assending order.
function gard.Helper.sortedMatrix(matrix) {
    veryLargeNumber = 100000000;
    sortedMatrix = {1,Columns(matrix)};

    for (i=0; i<Columns(matrix); i=i+1) {
        minElementLeft=Min(matrix, 1);
        sortedMatrix[i] = minElementLeft[0];
        matrix[minElementLeft[1]] = veryLargeNumber;
    }

    return sortedMatrix;
}

// Return the nth matrix in a nxm matrix.
function gard.Helper.getNthMatrix(setOfMatrices, n) {
    nthMatrix = {1, Columns(setOfMatrices)};
    for(i=0; i<Columns(setOfMatrices); i=i+1) {
        nthMatrix[i] = setOfMatrices[n][i];
    }
    return nthMatrix;
}

// Return the values of a dictionary as a matrix; not just the unique values and not as strings
function gard.Helper.getValues(dictionary) {
    keys = utility.Keys(dictionary);
    values = {1,Abs(dictionary)};
    for(i=0; i<Abs(dictionary); i=i+1) {
        values[i] = dictionary[keys[i]];
    }
    return values;
}