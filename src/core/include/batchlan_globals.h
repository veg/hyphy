/*
 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (spond@ucsd.edu)
 Art FY Poon    (apoon@cfenet.ubc.ca)
 Steven Weaver (sweaver@ucsd.edu)
 
 Module Developers:
 Lance Hepler (nlhepler@gmail.com)
 Martin Smith (martin.audacis@gmail.com)
 
 Significant contributions from:
 Spencer V Muse (muse@stat.ncsu.edu)
 Simon DW Frost (sdf22@cam.ac.uk)
 
 Permission is hereby granted, free of charge, to any person obtaining a
 copy of this software and associated documentation files (the
 "Software"), to deal in the Software without restriction, including
 without limitation the rights to use, copy, modify, merge, publish,
 distribute, sublicense, and/or sell copies of the Software, and to
 permit persons to whom the Software is furnished to do so, subject to
 the following conditions:
 
 The above copyright notice and this permission notice shall be included
 in all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
 */

#ifndef     __BATCHLANGLOBAL__
#define     __BATCHLANGLOBAL__


extern _List
dataSetList,
dataSetNamesList,
likeFuncList,   // list of all datasets
dataSetFilterList,
dataSetFilterNamesList,
likeFuncNamesList, // list of all dataset filters
pathNames,
theModelList,
allowedFormats,
batchLanguageFunctions,
batchLanguageFunctionNames,
batchLanguageFunctionParameterLists,
compiledFormulaeParameters,
modelNames,
executionStack,
standardLibraryPaths,
standardLibraryExtensions,
loadedLibraryPathsBackend;

extern _SimpleList
returnlist,
batchLanguageFunctionParameters,
batchLanguageFunctionClassification,
modelMatrixIndices,
modelTypeList,
modelFrequenciesIndices,
listOfCompiledFormulae;

extern _String  blFor,               // moved
         blWhile,         // moved
         blFunction,      // moved
         blFFunction,     // moved
         blLFunction,     // moved
         blReturn,        // moved
         blReturn2,            // moved
         blIf,            // moved
         blElse,           // moved
         blDo,            // moved
         blBreak,         // moved
         blContinue,          // moved
         blInclude,           // moved
         blDataSet,           // moved
         blDataSetFilter,
         blConstructCM,
         blTree,
         blLF,
         blLF3,
         blMolClock,
         blfprintf,
         blGetString,
         blfscanf,
         blsscanf,
         blExport,
         blReplicate,
         blImport,
         blCategory,
         blClearConstraints,
         blSetDialogPrompt,
         blModel,
         blChoiceList,
         blOpenDataPanel,
         blGetInformation,
         blExecuteCommands,
         blExecuteAFile,
         blLoadFunctionLibrary,
         blOpenWindow,
         blSpawnLF,
         blDifferentiate,
         blFindRoot,
         blMPIReceive,
         blMPISend,
         blGetDataInfo,
         blStateCounter,
         blIntegrate,
         blLFCompute,
         blGetURL,
         blDoSQL,
         blTopology,
         blAlignSequences,
         blGetNeutralNull,
         blHBLProfile,
         blDeleteObject,
         blRequireVersion,
         blSCFG,
         blNN,
         blBGM,
         blSimulateDataSet,
         blAssert;

extern _String globalPolynomialCap,
    enforceGlobalPolynomialCap,
    dropPolynomialTerms,
    maxPolyTermsPerVariable,
    maxPolyExpIterates,
    polyExpPrecision,
    systemVariableDump,
    selfDump,
    printDigitsSpec,
    explicitFormMExp,
    multByFrequencies,
    getDString,
    useLastFString,
    getFString,
    tempFString,
    defFileString,
    useLastModel,
    VerbosityLevelString,
    hasEndBeenReached,
    clearFile,
    keepFileOpen,
    closeFile,
    useLastDefinedMatrix,
    MessageLogging,
    selectionStrings,
    useNoModel,
    stdoutDestination,
    messageLogDestination,
    lastModelParameterList,
    dataPanelSourcePath,
    windowTypeTree,
    windowTypeClose,
    windowTypeTable,
    windowTypeDistribTable,
    windowTypeDatabase,
    screenWidthVar,
    screenHeightVar,
    useNexusFileData,
    mpiMLELFValue,
    lf2SendBack,
    pcAmbiguitiesResolve,
    pcAmbiguitiesAverage,
    pcAmbiguitiesSkip,
    lfStartCompute,
    lfDoneCompute,
    getURLFileFlag,
    versionString,
    timeStamp,
    simulationFilter,
    prefixDS,
    prefixDF,
    prefixLF,
    replaceTreeStructure,
    hyphyBaseDirectory,
    hyphyLibDirectory,
    platformDirectorySeparator,
    covarianceParameterList,
    matrixEvalCount,
    scfgCorpus,
    _hyLastExecutionError,
    _hyExecutionErrorMode,
    bgmData,
    bgmScores,
    bgmGraph,
    bgmNodeOrder,
    bgmConstraintMx,
    bgmParameters,
    pathToCurrentBF,
    hfCountGap,
    gdiDFAtomSize,
    statusBarProgressValue,
    statusBarUpdateString,
    marginalAncestors,
    doLeavesAncestors,
    blScanfRewind,
    blFprintfRedirect,
    blFprintfDevNull,
    getDataInfoReturnsOnlyTheIndex,
    alwaysReloadLibraries,
    dialogPrompt,
    baseDirectory,
    lastModelUsed,
    libDirectory,
    defFileNameValue;

extern long scanfLastReadPosition;

extern _String  sqlOpen,
                sqlClose,
                sqlRowData,
                sqlColNames,
                seqAlignMap,
                seqAlignScore,
                seqAlignScoreCodon2,
                seqAlignScoreCodon1,
                seqAlignGapChar,
                seqAlignGapOpen,
                seqAlignGapExtend,
                seqAlignGapOpen2,
                seqAlignGapExtend2,
                seqAlignFrameShift,
                seqAlignGapLocal,
                seqAlignGapAffine,
                seqAlignCodonAlign,
                seqAlignGapLinearSpace,
                seqAlignGapCodon3x1,
                seqAlignGapCodon3x2,
                seqAlignGapCodon3x4,
                seqAlignGapCodon3x5,
                seqAlignDoLocal,
                completeFlag,
                conditionalWeights,
                siteProbabilities,
                lastSetOfConstraints,
                deferConstrainAssignment,
                assertionBehavior,
                _hyStatusConditionProbsMatrix,
                isDynamicGraph;

extern _SimpleList sqlDatabases,
            _HY_HBLCommandHelperAux;

#endif
