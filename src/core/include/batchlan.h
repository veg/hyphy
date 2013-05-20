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

#ifndef __BATCHLANGUAGE__
#define __BATCHLANGUAGE__

#include "legacy_parser.h"
#include "site.h"
#include "trie.h"

#include "associativelist.h"
#include "growingvector.h"

#include <stdio.h>

#define HY_BL_ERROR_HANDLING_DEFAULT 0
#define HY_BL_ERROR_HANDLING_SOFT 1

//____________________________________________________________________________________

class _ElementaryCommand;

//____________________________________________________________________________________
_ElementaryCommand *makeNewCommand(long);

//____________________________________________________________________________________

#ifdef __HYPHYMPI__
#include <mpi.h>

extern _String mpiNodeID, mpiNodeCount;

#define HYPHY_MPI_SIZE_TAG 111
#define HYPHY_MPI_STRING_TAG 112
#define HYPHY_MPI_DONE_TAG 113
#define HYPHY_MPI_VARS_TAG 114
#define HYPHY_MPI_DATA_TAG 115
#define HYPHY_MPI_DIE_TAG 666

void ReportMPIError(int, bool);
void MPISendString(_String &, long, bool = false);
_String *MPIRecvString(long, long &);

#endif

//____________________________________________________________________________________
extern _List batchLanguageFunctions, batchLanguageFunctionNames,
    batchLanguageFunctionParameterLists, dataSetList, dataSetNamesList,
    likeFuncList, dataSetFilterList, dataSetFilterNamesList, templateModelList,
    scfgNamesList, scfgList, bgmNamesList, // modified by afyp
    bgmList, likeFuncNamesList, modelNames, executionStack,
    compiledFormulaeParameters, standardLibraryPaths, standardLibraryExtensions;

extern _SimpleList batchLanguageFunctionParameters,
    batchLanguageFunctionClassification, modelMatrixIndices, modelTypeList,
    // SLKP: 20100313 this list stores 0 for  normal (rate-matrix based models),
    //       vs expression based matrices, for which the dimension is stored.
    modelFrequenciesIndices, listOfCompiledFormulae;

extern _String getDString, getFString, tempFString, baseDirectory, libDirectory,
    useLastFString, mpiMLELFValue, lf2SendBack, hyphyBaseDirectory,
    hyphyLibDirectory, platformDirectorySeparator, defFileNameValue,
    defFileString, blConstructCM, blFprintfRedirect, blFprintfDevNull,
    globalPolynomialCap, enforceGlobalPolynomialCap, dropPolynomialTerms,
    maxPolyTermsPerVariable, maxPolyExpIterates, polyExpPrecision,
    systemVariableDump, selfDump, printDigitsSpec, explicitFormMExp,
    multByFrequencies, getDString, useLastFString, getFString, defFileString,
    useLastModel, VerbosityLevelString, hasEndBeenReached, clearFile,
    keepFileOpen, closeFile, useLastDefinedMatrix, MessageLogging,
    selectionStrings, useNoModel, stdoutDestination, messageLogDestination,
    lastModelParameterList, dataPanelSourcePath, windowTypeTree,
    windowTypeClose, windowTypeTable, windowTypeDistribTable,
    windowTypeDatabase, screenWidthVar, screenHeightVar, useNexusFileData,
    mpiMLELFValue, lf2SendBack, pcAmbiguitiesResolve, pcAmbiguitiesAverage,
    pcAmbiguitiesSkip, lfStartCompute, lfDoneCompute, getURLFileFlag,
    versionString, timeStamp, simulationFilter, prefixDS, prefixDF, prefixLF,
    replaceTreeStructure, hyphyBaseDirectory, platformDirectorySeparator,
    covarianceParameterList, matrixEvalCount, scfgCorpus, bgmData, bgmWeights,
    pathToCurrentBF, statusBarUpdateString, statusBarProgressValue,
    errorReportFormatExpression, errorReportFormatExpressionStr,
    errorReportFormatExpressionStack, errorReportFormatExpressionStdin,
    lastModelUsed, deferConstrainAssignment, bgmData, bgmScores, bgmGraph,
    bgmNodeOrder, bgmConstraintMx, bgmParameters, assertionBehavior,
    dialogPrompt, _hyLastExecutionError, _hyExecutionErrorMode,
#ifdef __HYPHYMPI__
    mpiNodeID, mpiNodeCount, mpiLastSentMsg,
#endif
    hfCountGap;

extern _ExecutionList *currentExecutionList;

extern _AVLList loadedLibraryPaths;
extern _AVLListX _HY_HBLCommandHelper, _HY_GetStringGlobalTypes;

extern _Trie _HY_ValidHBLExpressions, _HY_HBL_Namespaces;

extern long globalRandSeed, matrixExpCount;

long FindDataSetName(_String &);
long FindDataSetFilterName(_String &);
long FindSCFGName(_String &);
long FindBFFunctionName(_String &, _VariableContainer * = nil);

long FindBgmName(_String &);
// added by afyp, March 18, 2007

long FindLikeFuncName(_String &, bool = false);
long FindModelName(_String &);
void ScanModelForVariables(long modelID, _AVLList &theReceptacle, bool inclG,
                           long modelID2, bool inclCat);
/* 20100316 SLKP:
    factored out a function call to scan a particular model
    for variables to permit the use of explicit (formula-based) model
    definitions
 */

void ReadBatchFile(_String &, _ExecutionList &);
_String ReturnDialogInput(bool dispPath = false);
_String ReturnFileDialogInput(void);
_String *ProcessCommandArgument(_String *);
_String WriteFileDialogInput(void);
_Parameter ProcessNumericArgument(_String *, _VariableContainer *,
                                  _ExecutionList * = nil);
_String ProcessLiteralArgument(_String *, _VariableContainer *,
                               _ExecutionList * = nil);
_AssociativeList *ProcessDictionaryArgument(_String *data,
                                            _VariableContainer *theP,
                                            _ExecutionList * = nil);

_String GetStringFromFormula(_String *, _VariableContainer *);
void ExecuteBLString(_String &, _VariableContainer *);

void SerializeModel(_String &, long, _AVLList * = nil, bool = false);
bool Get_a_URL(_String &, _String * = nil);

long AddFilterToList(_String &, _DataSetFilter *, bool = false);
long AddDataSetToList(_String &, _DataSet *);
void SetDataFilterParameters(_String &, _DataSetFilter *, bool);
void KillDataFilterRecord(long, bool = false);
void KillLFRecord(long, bool = true);
void KillDataSetRecord(long);
void KillModelRecord(long);
void KillExplicitModelFormulae(void);
bool PushFilePath(_String &, bool = true);
void PopFilePath(void);
_Matrix *CheckMatrixArg(_String *, bool);
_AssociativeList *CheckAssociativeListArg(_String *);
void RetrieveModelComponents(long, _Matrix *&, _Matrix *&, bool &);
void RetrieveModelComponents(long, _Variable *&, _Variable *&, bool &);
bool IsModelReversible(long);
bool IsModelOfExplicitForm(long);

void ReadModelList(void);
_String ProcessStringArgument(_String *data);
_String *_HBLObjectNameByType(const long type, const long index,
                              bool correct_for_empties = true);
_String _hblCommandAccessor(_ExecutionList *, long);
_String _HYGenerateANameSpace(void);

_PMathObj ProcessAnArgumentByType(_String *, _VariableContainer *, long,
                                  _ExecutionList * = nil);

void _HBL_Init_Const_Arrays(void);

void ReturnCurrentCallStack(_List &, _List &);

/**
    An accessor function which attempts to retrieve a reference to a HyPhy Batch
    Language Object
    by name. A list of acceptable object classes can be specified in the type
    parameter. Note that
    types will be searched in the following order:

        HY_BL_DATASET,HY_BL_DATASET_FILTER,HY_BL_LIKELIHOOD_FUNCTION,HY_BL_SCFG,HY_BL_BGM,HY_BL_MODEL,HY_BL_HBL_FUNCTION

    i.e. if there is a dataset named 'foo' and a likelihood function named
    'foo', then the dataset
    will be returned.



    @param   name provides a string with the name of the object to be retrieved.
    @param   type [in] which types of objects will be searched.
                 [out] which type of object was retrieved (HY_BL_NOT_DEFINED if
    not found)
    @param   index (if not nil) will receive the index of the found object in
    the corresponding array
    @param   errMsg if set to True, will cause the function to report an error
    if no object of corresponding type could be found
    @param   tryLiteralLookup if set to True, will cause the function to, upon a
    failed lookup, to also try interpreting name as a string variable ID
    @return  pointer to the retrieved object or nil if not found
    @author  SLKP
    @version 20120324
*/

BaseRef _HYRetrieveBLObjectByName(_String &name, long &type, long *index = nil,
                                  bool errMsg = false,
                                  bool tryLiteralLookup = false);

_String _HYHBLTypeToText(long type);
_String _HYStandardDirectory(const unsigned long);

_HBLCommandExtras *_hyInitCommandExtras(const long = 0, const long = 0,
                                        const _String = empty, const char = ';',
                                        const bool = true, const bool = false,
                                        const bool = false,
                                        _SimpleList * = nil);

extern bool numericalParameterSuccessFlag;
extern _Parameter messageLogFlag;

bool RecurseDownTheTree(_SimpleList &, _List &, _List &, _List &,
                        _SimpleList &);

#endif
