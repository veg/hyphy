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

#include "likefunc.h"
#include "batchlan.h"
#include "string.h"
#include "ctype.h"
#include "polynoml.h"
#include "time.h"
#include "scfg.h"
#include "HYNetInterface.h"
#include "hy_globals.h"
#include "executionlist.h"
#include "datasetfilter.h"
#include "thetree.h"

#include "bayesgraph.h"

#if not defined __HEADLESS__ &&defined __HYPHYQT__
//#include "HYSharedMain.h"
#include "hyphy_qt_helpers.h"
#endif

//______________________________________________________________________________
// global variables
_List dataSetList, dataSetNamesList, likeFuncList, // list of all datasets
    dataSetFilterList, dataSetFilterNamesList,
    likeFuncNamesList, // list of all dataset filters
    pathNames, theModelList, allowedFormats, batchLanguageFunctions,
    batchLanguageFunctionNames, batchLanguageFunctionParameterLists,
    compiledFormulaeParameters, modelNames, executionStack,
    standardLibraryPaths, standardLibraryExtensions, loadedLibraryPathsBackend;

// retrieval functions

_SimpleList returnlist, batchLanguageFunctionParameters,
    batchLanguageFunctionClassification, modelMatrixIndices, modelTypeList,
    modelFrequenciesIndices, listOfCompiledFormulae;

_String globalPolynomialCap("GLOBAL_POLYNOMIAL_CAP"),
    enforceGlobalPolynomialCap("ENFORCE_GLOBAL_POLYNOMIAL_CAP"),
    dropPolynomialTerms("DROP_POLYNOMIAL_TERMS"),
    maxPolyTermsPerVariable("MAX_POLY_TERMS_PER_VARIABLE"),
    maxPolyExpIterates("MAX_POLYNOMIAL_EXP_ITERATES"),
    polyExpPrecision("POLYNOMIAL_EXP_PRECISION"),
    systemVariableDump("LIST_ALL_VARIABLES"), selfDump("PRINT_SELF"),
    printDigitsSpec("PRINT_DIGITS"),
    explicitFormMExp("EXPLICIT_FORM_MATRIX_EXPONENTIAL"),
    multByFrequencies("MULTIPLY_BY_FREQUENCIES"),
    getDString("PROMPT_FOR_STRING"), useLastFString("LAST_FILE_PATH"),
    getFString("PROMPT_FOR_FILE"), tempFString("TEMP_FILE_NAME"),
    defFileString("DEFAULT_FILE_SAVE_NAME"), useLastModel("USE_LAST_MODEL"),
    VerbosityLevelString("VERBOSITY_LEVEL"), hasEndBeenReached("END_OF_FILE"),
    clearFile("CLEAR_FILE"), keepFileOpen("KEEP_OPEN"), closeFile("CLOSE_FILE"),
    useLastDefinedMatrix("USE_LAST_DEFINED_MATRIX"),
    MessageLogging("MESSAGE_LOGGING"), selectionStrings("SELECTION_STRINGS"),
    useNoModel("USE_NO_MODEL"), stdoutDestination("stdout"),
    messageLogDestination("MESSAGE_LOG"),
    lastModelParameterList("LAST_MODEL_PARAMETER_LIST"),
    dataPanelSourcePath("DATA_PANEL_SOURCE_PATH"), windowTypeTree("TREEWINDOW"),
    windowTypeClose("CLOSEWINDOW"), windowTypeTable("CHARTWINDOW"),
    windowTypeDistribTable("DISTRIBUTIONWINDOW"),
    windowTypeDatabase("DATABASEWINDOW"), screenWidthVar("SCREEN_WIDTH"),
    screenHeightVar("SCREEN_HEIGHT"), useNexusFileData("USE_NEXUS_FILE_DATA"),
    mpiMLELFValue("MPI_MLE_LF_VALUE"),
    lf2SendBack("LIKE_FUNC_NAME_TO_SEND_BACK"),
    pcAmbiguitiesResolve("RESOLVE_AMBIGUITIES"),
    pcAmbiguitiesAverage("AVERAGE_AMBIGUITIES"),
    pcAmbiguitiesSkip("SKIP_AMBIGUITIES"), lfStartCompute("LF_START_COMPUTE"),
    lfDoneCompute("LF_DONE_COMPUTE"), getURLFileFlag("SAVE_TO_FILE"),
    versionString("HYPHY_VERSION"), timeStamp("TIME_STAMP"),
    simulationFilter("_SIM_INTERNAL_FILTER_"), prefixDS("DataSet_"),
    prefixDF("Partition_"), prefixLF("LF_"),
    replaceTreeStructure("REPLACE_TREE_STRUCTURE"),
    hyphyBaseDirectory("HYPHY_BASE_DIRECTORY"),
    hyphyLibDirectory("HYPHY_LIB_DIRECTORY"),
    platformDirectorySeparator("DIRECTORY_SEPARATOR"),
    covarianceParameterList("COVARIANCE_PARAMETER"),
    matrixEvalCount("MATRIX_EXPONENTIATION_COUNTS"),
    scfgCorpus("SCFG_STRING_CORPUS"),
    _hyLastExecutionError("LAST_HBL_EXECUTION_ERROR"),
    _hyExecutionErrorMode("HBL_EXECUTION_ERROR_HANDLING"),
    bgmData("BGM_DATA_MATRIX"), bgmScores("BGM_SCORE_CACHE"),
    bgmGraph("BGM_GRAPH_MATRIX"), bgmNodeOrder("BGM_NODE_ORDER"),
    bgmConstraintMx("BGM_CONSTRAINT_MATRIX"),
    bgmParameters("BGM_NETWORK_PARAMETERS"),
    pathToCurrentBF("PATH_TO_CURRENT_BF"),
    hfCountGap("COUNT_GAPS_IN_FREQUENCIES"), gdiDFAtomSize("ATOM_SIZE"),
    statusBarProgressValue("STATUS_BAR_PROGRESS_VALUE"),
    statusBarUpdateString("STATUS_BAR_STATUS_STRING"),
    marginalAncestors("MARGINAL"), doLeavesAncestors("DOLEAVES"),
    blScanfRewind("REWIND"), blFprintfRedirect("GLOBAL_FPRINTF_REDIRECT"),
    blFprintfDevNull("/dev/null"),
    getDataInfoReturnsOnlyTheIndex("GET_DATA_INFO_RETURNS_ONLY_THE_INDEX"),
    alwaysReloadLibraries("ALWAYS_RELOAD_FUNCTION_LIBRARIES"), dialogPrompt,
    baseDirectory, lastModelUsed, libDirectory, defFileNameValue;

#ifdef __HYPHYMPI__

_String mpiNodeID("MPI_NODE_ID"), mpiNodeCount("MPI_NODE_COUNT"),
    mpiLastSentMsg("MPI_LAST_SENT_MSG");

void ReportMPIError(int, bool);

#endif

_Parameter explicitFormMatrixExponential = 0.0, messageLogFlag = 1.0;
long scanfLastReadPosition = 0;
extern _String MATRIX_AGREEMENT, ANAL_COMP_FLAG;
extern _Parameter toPolyOrNot, toMorNot2M, ANALYTIC_COMPUTATION_FLAG;
extern _SimpleList freeSlots;
_AVLList loadedLibraryPaths (&loadedLibraryPathsBackend);

//______________________________________________________________________________
// Function prototypes

#ifdef __HYPHYMPI__

//______________________________________________________________________________
void ReportMPIError(int code, bool send) {
  if (code != MPI_SUCCESS) {
    _String errMsg = "MPI Error while ";
    if (send) {
      errMsg = errMsg & "sending";
    } else {
      errMsg = errMsg & "receiving";
    }
    errMsg = errMsg & _String(" code:") & (long) code;
    FlagError(errMsg);
  }
}

#define MPI_SEND_CHUNK 0xFFFFFFL

//______________________________________________________________________________
void MPISendString(_String &theMessage, long destID, bool isError) {
  long messageLength = theMessage.sLength, transferCount = 0;

  if (isError) {
    messageLength = -messageLength;
  }

  ReportMPIError(MPI_Send(&messageLength, 1, MPI_LONG, destID,
                          HYPHY_MPI_SIZE_TAG, MPI_COMM_WORLD),
                 true);

  if (messageLength == 0) {
    return;
  }

  if (isError) {
    messageLength = -messageLength;
  }

  while (messageLength - transferCount > MPI_SEND_CHUNK) {
    ReportMPIError(
        MPI_Send(theMessage.sData + transferCount, MPI_SEND_CHUNK, MPI_CHAR,
                 destID, HYPHY_MPI_STRING_TAG, MPI_COMM_WORLD),
        true);
    transferCount += MPI_SEND_CHUNK;
  }

  if (messageLength - transferCount) {
    ReportMPIError(MPI_Send(theMessage.sData + transferCount,
                            messageLength - transferCount, MPI_CHAR, destID,
                            HYPHY_MPI_STRING_TAG, MPI_COMM_WORLD),
                   true);
  }

  //ReportMPIError(MPI_Send(&messageLength, 1, MPI_LONG, destID,
  //HYPHY_MPI_DONE_TAG, MPI_COMM_WORLD),true);

  _FString *sentVal = new _FString;
  sentVal->theString = (_String *)theMessage.makeDynamic();
  _Variable *mpiMsgVar = CheckReceptacle(&mpiLastSentMsg, empty, false);
  mpiMsgVar->SetValue(sentVal, false);
  //setParameter (mpiLastSentMsg, &sentVal);

}

//______________________________________________________________________________
_String *MPIRecvString(long senderT, long &senderID) {

  _String *theMessage = nil;
  long messageLength = 0, transferCount = 0;
  int actualReceived = 0;
  bool isError = false;

  if (senderT < 0) {
    senderT = MPI_ANY_SOURCE;
  }

  MPI_Status status;

  ReportMPIError(MPI_Recv(&messageLength, 1, MPI_LONG, senderT, HYPHY_MPI_SIZE_TAG, MPI_COMM_WORLD, &status), false);

  if (messageLength < 0) {
    isError = true;
    messageLength = -messageLength;
  }

  if (!isError) {
    //MPI_Get_count (&status,MPI_CHAR,&actualReceived);

    if (messageLength == 0) {
      return new _String;
    }

    theMessage = new _String(messageLength, false);

    senderT = senderID = status.MPI_SOURCE;

    while (messageLength - transferCount > MPI_SEND_CHUNK) {
      ReportMPIError(
          MPI_Recv(theMessage->sData + transferCount, MPI_SEND_CHUNK, MPI_CHAR,
                   senderT, HYPHY_MPI_STRING_TAG, MPI_COMM_WORLD, &status),
          false);
      MPI_Get_count(&status, MPI_CHAR, &actualReceived);
      if (actualReceived != MPI_SEND_CHUNK) {
        WarnError("Failed in MPIRecvString - some data was not properly received\n");
      }
      //return    nil;
      transferCount += MPI_SEND_CHUNK;
    }

    if (messageLength - transferCount) {
      ReportMPIError(MPI_Recv(theMessage->sData + transferCount,messageLength - transferCount, MPI_CHAR, senderT,HYPHY_MPI_STRING_TAG, MPI_COMM_WORLD, &status),false);
      MPI_Get_count(&status, MPI_CHAR, &actualReceived);
      if (actualReceived != messageLength - transferCount) {
        WarnError("Failed in MPIRecvString - some data was not properly received\n");
      }
      //return    nil;
    }
    //ReportMPIError(MPI_Recv(&messageLength, 1, MPI_LONG, senderT,
    //HYPHY_MPI_DONE_TAG, MPI_COMM_WORLD,&status),false);

    if (isError) {
      FlagError(theMessage);
    }
  }
  return theMessage;
}
#endif

//______________________________________________________________________________
_String GetStringFromFormula(_String *data, _VariableContainer *theP) {
  _Formula nameForm(*data, theP);
  _PMathObj formRes = nameForm.Compute();

  if (formRes && formRes->ObjectClass() == STRING) {
    data = ((_FString *)formRes)->theString;
  }

  return *data;
}

//______________________________________________________________________________
_String *ProcessCommandArgument(_String *data) {
  if (data->sLength > 1 && data->sData[data->sLength - 1] == '&') {
    _String argName(*data, 0, data->sLength - 2);
    _FString *theVar =
        (_FString *)FetchObjectFromVariableByType(&argName, STRING);

    if (theVar) {
      return theVar->theString;
    }

    WarnError(_String("Reference argument \"") & *data &
              "\" is not a valid string variable.");
    return nil;
  }
  return (_String *)data;
}

//______________________________________________________________________________
bool numericalParameterSuccessFlag = true;

_Parameter ProcessNumericArgument(_String *data, _VariableContainer *theP,
                                  _ExecutionList *currentProgram) {
  _String errMsg;
  _Formula nameForm(*data, theP, currentProgram ? &errMsg : nil);

  if (!errMsg.sLength) {
    _PMathObj formRes = nameForm.Compute();
    numericalParameterSuccessFlag = true;
    if (formRes && formRes->ObjectClass() == NUMBER) {
      return formRes->Value();
    } else {
      if (formRes && formRes->ObjectClass() == STRING) {
        return _String((_String *)((_FString *)formRes)->toStr()).toNum();
      } else {
        errMsg = (_String("'") & *data &
                  "' was expected to be a numerical argument.");
      }
    }
  }

  if (currentProgram) {
    currentProgram->ReportAnExecutionError(errMsg);
  } else {
    WarnError(errMsg);
  }

  numericalParameterSuccessFlag = false;
  return 0.0;
}

//______________________________________________________________________________
_PMathObj ProcessAnArgumentByType(_String *expression, _VariableContainer *theP,
                                  long objectType,
                                  _ExecutionList *currentProgram) {
  _String errMsg;

  _Formula expressionProcessor(*expression, theP,
                               currentProgram ? &errMsg : nil);

  if (errMsg.sLength) {
    currentProgram->ReportAnExecutionError(errMsg);
  } else {
    _PMathObj expressionResult = expressionProcessor.Compute(0, theP);
    if (expressionResult && expressionResult->ObjectClass() == objectType) {
      return (_PMathObj) expressionResult->makeDynamic();
    }
  }

  return nil;
}

//______________________________________________________________________________
_String ProcessLiteralArgument(_String *data, _VariableContainer *theP,
                               _ExecutionList *currentProgram) {
  _String errMsg;

  _Formula expressionProcessor(*data, theP, currentProgram ? &errMsg : nil);

  if (errMsg.sLength) {
    currentProgram->ReportAnExecutionError(errMsg);
  } else {
    _PMathObj expressionResult = expressionProcessor.Compute(0, theP);
    if (expressionResult && expressionResult->ObjectClass() == STRING) {
      return *((_FString *)expressionResult)->theString;
    }
  }

  return empty;
}

//______________________________________________________________________________
_AssociativeList *ProcessDictionaryArgument(_String *data,
                                            _VariableContainer *theP,
                                            _ExecutionList *currentProgram) {
  _String errMsg;
  _Formula nameForm(*data, theP, currentProgram ? &errMsg : nil);
  if (errMsg.sLength && currentProgram) {
    currentProgram->ReportAnExecutionError(errMsg);
  } else {
    _PMathObj formRes = nameForm.Compute();
    if (formRes && formRes->ObjectClass() == ASSOCIATIVE_LIST) {
      formRes->AddAReference();
      return (_AssociativeList *)formRes;
    }
  }
  return nil;
}

//______________________________________________________________________________
long FindDataSetName(_String &s) { return dataSetNamesList.Find(&s); }

//______________________________________________________________________________
long FindDataSetFilterName(_String &s) {
  return dataSetFilterNamesList.Find(&s);
}

//______________________________________________________________________________
long FindLikeFuncName(_String &s, bool tryAsAString) {
  long try1 = likeFuncNamesList.Find(&s);
  if (try1 < 0 && tryAsAString) {
    _String s2(ProcessLiteralArgument(&s, nil));
    try1 = likeFuncNamesList.Find(&s2);
  }
  return try1;
}

//______________________________________________________________________________
long FindModelName(_String &s) {
  if (s.Equal(&useLastModel)) {
    return lastMatrixDeclared;
  }
  return modelNames.Find(&s);
}

//______________________________________________________________________________
_LikelihoodFunction *FindLikeFuncByName(_String &s) {
  long i = FindLikeFuncName(s);
  if (i >= 0) {
    return (_LikelihoodFunction *)likeFuncList(i);
  }
  return nil;
}

//______________________________________________________________________________
long FindSCFGName(_String &s) { return scfgNamesList.Find(&s); }

//______________________________________________________________________________
long FindBFFunctionName(_String &s, _VariableContainer *theP) {
  if (theP) {
    _String testName = *(theP->GetName()) & '.' & s;

    long cutAt = testName.sLength - s.sLength - 2;
    do {
      long idx = batchLanguageFunctionNames.Find(&testName);
      if (idx >= 0) {
        s = testName;
        return idx;
      }
      testName.Trim(0, cutAt);
      cutAt = testName.FindBackwards('.', 0, -1) - 1;
    } while (cutAt >= 0);
  }

  return batchLanguageFunctionNames.Find(&s);
}

//______________________________________________________________________________
long FindBgmName(_String &s) { return bgmNamesList.Find(&s); }


//______________________________________________________________________________
long AddFilterToList(_String &partName, _DataSetFilter *theFilter, bool addP) {
  FindUnusedObjectName(prefixDF, partName, dataSetFilterNamesList);
  long k;

  for (k = 0; k < dataSetFilterNamesList.lLength; k++)
    if (((_String *)dataSetFilterNamesList(k))->sLength == 0) {
      break;
    }

  if (addP) {
    SetDataFilterParameters(partName, theFilter, true);
  }

  if (k == dataSetFilterNamesList.lLength) {
    dataSetFilterList << theFilter;
    DeleteObject(theFilter);
    dataSetFilterNamesList &&&partName;
    return dataSetFilterNamesList.lLength - 1;
  }
  dataSetFilterList.lData[k] = (long) theFilter;
  dataSetFilterNamesList.Replace(k, &partName, true);
  return k;
}

//______________________________________________________________________________
long AddDataSetToList(_String &theName, _DataSet *theDS) {
  FindUnusedObjectName(prefixDS, theName, dataSetNamesList);
  long k = dataSetNamesList.Find(&empty);
  if (k == -1) {
    dataSetList.AppendNewInstance(theDS);
    dataSetNamesList &&&theName;
    k = dataSetNamesList.lLength - 1;
  } else {
    dataSetNamesList.Replace(k, &theName, true);
    dataSetList.lData[k] = (long) theDS;
  }
  return k;
}

//______________________________________________________________________________
void KillDataFilterRecord(long dfID, bool addP) {
  if (addP) {
    SetDataFilterParameters(*(_String *)(dataSetFilterNamesList(dfID)), nil, false);
  }

  if (dfID < dataSetFilterList.lLength - 1) {
    DeleteObject(dataSetFilterList(dfID));
    dataSetFilterList.lData[dfID] = 0;
    dataSetFilterNamesList.Replace(dfID, &empty, true);
  } else {
    dataSetFilterList.Delete(dfID);
    dataSetFilterNamesList.Delete(dfID);
    if (dfID)
      while (((_String *)dataSetFilterNamesList(--dfID))->sLength == 0) {
        dataSetFilterList.Delete(dfID);
        dataSetFilterNamesList.Delete(dfID);
        if (dfID == 0) {
          break;
        }
      }
  }
}

//______________________________________________________________________________
void KillLFRecord(long lfID, bool completeKill) {
  /* compile the list of variables which will no longer be referenced */

  if (lfID >= 0) {
    _LikelihoodFunction *me = (_LikelihoodFunction *)likeFuncList(lfID);

    if (completeKill) {
      _SimpleList wastedVars, otherVars, myVars, otherModels, wastedModels;

      long k;

      myVars << me->GetIndependentVars();
      myVars << me->GetDependentVars();

      for (k = 0; k < likeFuncList.lLength; k++)
        if (k != lfID) {
          if (((_String *)likeFuncNamesList(k))->sLength) {
            _LikelihoodFunction *lf = (_LikelihoodFunction *)likeFuncList(k);
            otherVars << lf->GetIndependentVars();
            otherVars << lf->GetDependentVars();
            for (long kk = lf->GetTheTrees().lLength - 1; kk >= 0; kk--) {
              _TheTree *thisTree =
                  (_TheTree *)LocateVar(lf->GetTheTrees().lData[kk]);
              thisTree->CompileListOfModels(otherModels);
            }
          }
        }

      otherVars.Sort();
      otherModels.Sort();

      for (k = 0; k < myVars.lLength; k++)
        if (otherVars.BinaryFind(myVars.lData[k]) < 0) {
          wastedVars << myVars.lData[k];
        }

      myVars.Clear();

      for (k = me->GetTheTrees().lLength - 1; k >= 0; k--) {
        _TheTree *thisTree = (_TheTree *)LocateVar(me->GetTheTrees().lData[k]);
        thisTree->CompileListOfModels(myVars);
        _CalcNode *tNode = thisTree->DepthWiseTraversal(true);
        while (tNode) {
          tNode->SetValue(new _Constant(tNode->BranchLength()), false);
          tNode = thisTree->DepthWiseTraversal();
        }
        thisTree->RemoveModel();
      }

      for (k = 0; k < myVars.lLength; k++)
        if (otherModels.BinaryFind(myVars.lData[k]) < 0) {
          KillModelRecord(myVars.lData[k]);
        }

      for (k = 0; k < wastedVars.lLength; k++) {
        DeleteVariable(*LocateVar(wastedVars.lData[k])->GetName());
      }

    }

    if (lfID < likeFuncList.lLength - 1) {
      DeleteObject(likeFuncList(lfID));
      likeFuncList.lData[lfID] = nil;
      likeFuncNamesList.Replace(lfID, &empty, true);
    } else {
      likeFuncList.Delete(lfID);
      likeFuncNamesList.Delete(lfID);
      if (lfID)
        while (((_String *)likeFuncNamesList(--lfID))->sLength == 0) {
          likeFuncList.Delete(lfID);
          likeFuncNamesList.Delete(lfID);
          if (lfID == 0) {
            break;
          }
        }
    }
  }
}

//______________________________________________________________________________
void KillLFRecordFull(long lfID) {
  _LikelihoodFunction *lf = (_LikelihoodFunction *)likeFuncList(lfID);

  long k;
  //for (k=lf->GetTheFilters().lLength-1; k>=0; k--)
  //  KillDataFilterRecord (lf->GetTheFilters().lData[k]);

  _SimpleList l;

  lf->GetGlobalVars(l);

  for (k = 0; k < l.lLength; k++) {
    DeleteVariable(*LocateVar(l.lData[k])->GetName());
  }

  l.Clear();

  for (k = lf->GetTheTrees().lLength - 1; k >= 0; k--) {
    _TheTree *thisTree = (_TheTree *)LocateVar(lf->GetTheTrees().lData[k]);
    thisTree->CompileListOfModels(l);
    DeleteVariable(*thisTree->GetName());
  }

  for (k = 0; k < l.lLength; k++) {
    KillModelRecord(l.lData[k]);
  }

  KillLFRecord(lfID);
}

//______________________________________________________________________________
void KillDataSetRecord(long dsID) {
  if (dsID < dataSetList.lLength - 1) {
    DeleteObject(dataSetList(dsID));
    dataSetList.lData[dsID] = 0;
    dataSetNamesList.Replace(dsID, &empty, true);
  } else {
    dataSetList.Delete(dsID);
    dataSetNamesList.Delete(dsID);
    if (dsID)
      while (((_String *)dataSetNamesList(--dsID))->sLength == 0) {
        dataSetList.Delete(dsID);
        dataSetNamesList.Delete(dsID);
        if (dsID == 0) {
          break;
        }
      }
  }
}

//______________________________________________________________________________
void KillExplicitModelFormulae(void) {
  for (long i = 0; i < modelTypeList.lLength; i++)
    if (modelTypeList.lData[i]) {
      delete (_Formula *)(modelMatrixIndices.lData[i]);
    }
}

//______________________________________________________________________________
void KillModelRecord(long mdID) {
  if (lastMatrixDeclared == mdID) {
    lastMatrixDeclared = -1;
  }

  // SLKP 20110816
  // can't delete matrices before checking that no other model depends on the

  if (modelTypeList.lData[mdID]) {
    delete (_Formula *)(modelMatrixIndices.lData[mdID]);
  } else {
    _Variable *modelMatrix = nil, *freqMatrix = nil;

    bool multByFreqs = false;

    _SimpleList saveTheseVariablesAux;
    _AVLList saveTheseVariables(&saveTheseVariablesAux);

    for (long k = 0; k < modelNames.lLength; k++) {
      if (k != mdID && ((_String *)modelNames(k))->sLength) {
        if (modelTypeList.lData[k]) {
          _SimpleList dependantMatrices;
          ((_Formula *)(modelMatrixIndices.lData[k]))
              ->ScanFForType(dependantMatrices, MATRIX);
          for (long k2 = 0; k2 < dependantMatrices.lLength; k2++) {
            saveTheseVariables.Insert((BaseRef) dependantMatrices.lData[k2]);
          }
        } else {
          RetrieveModelComponents(k, modelMatrix, freqMatrix, multByFreqs);

          if (modelMatrix) {
            saveTheseVariables.Insert((BaseRef) modelMatrix->GetIndex());
          }
          if (freqMatrix) {
            saveTheseVariables.Insert((BaseRef) freqMatrix->GetIndex());
          }
        }
      }
    }

    RetrieveModelComponents(mdID, modelMatrix, freqMatrix, multByFreqs);
    if (modelMatrix &&
        saveTheseVariables.Find((BaseRef) modelMatrix->GetIndex()) < 0) {
      DeleteVariable(*modelMatrix->GetName());
    }
    if (freqMatrix &&
        saveTheseVariables.Find((BaseRef) freqMatrix->GetIndex()) < 0) {
      DeleteVariable(*freqMatrix->GetName());
    }
  }

  if (mdID < modelNames.lLength - 1) {
    modelMatrixIndices.lData[mdID] = -1;
    modelTypeList.lData[mdID] = 0;
    modelFrequenciesIndices.lData[mdID] = -1;
    modelNames.Replace(mdID, &empty, true);
  } else {
    modelNames.Delete(mdID);
    modelMatrixIndices.Delete(modelMatrixIndices.lLength - 1);
    modelFrequenciesIndices.Delete(modelFrequenciesIndices.lLength - 1);
    modelTypeList.Delete(modelTypeList.lLength - 1);

    if (mdID)
      while (((_String *)modelNames(--mdID))->sLength == 0) {
        modelNames.Delete(mdID);
        modelMatrixIndices.Delete(mdID);
        modelFrequenciesIndices.Delete(mdID);
        modelTypeList.Delete(mdID);
        if (mdID == 0) {
          break;
        }
      }
  }
}

_String blFor("for("),         // moved
    blWhile("while("),         // moved
    blFunction("function "),   // moved
    blFFunction("ffunction "), // moved
    blLFunction("lfunction "), // moved
    blReturn("return "),       // moved
    blReturn2("return("),      // moved
    blIf("if("),               // moved
    blElse("else"),            // moved
    blDo("do{"),               // moved
    blBreak("break;"),         // moved
    blContinue("continue;"),   // moved
    blInclude("#include"),     // moved
    blDataSet("DataSet "),     // moved
    blDataSetFilter("DataSetFilter "),
    blConstructCM("ConstructCategoryMatrix("), blTree("Tree "),
    blLF("LikelihoodFunction "), blLF3("LikelihoodFunction3 "),
    blMolClock("MolecularClock("), blfprintf("fprintf("),
    blGetString("GetString("), blfscanf("fscanf("), blsscanf("sscanf("),
    blExport("Export("), blReplicate("ReplicateConstraint("),
    blImport("Import"), blCategory("category "),
    blClearConstraints("ClearConstraints("),
    blSetDialogPrompt("SetDialogPrompt("), blModel("Model "),
    blChoiceList("ChoiceList("), blOpenDataPanel("OpenDataPanel("),
    blGetInformation("GetInformation("), blExecuteCommands("ExecuteCommands("),
    blExecuteAFile("ExecuteAFile("),
    blLoadFunctionLibrary("LoadFunctionLibrary("), blOpenWindow("OpenWindow("),
    blSpawnLF("SpawnLikelihoodFunction("), blDifferentiate("Differentiate("),
    blFindRoot("FindRoot("), blMPIReceive("MPIReceive("), blMPISend("MPISend("),
    blGetDataInfo("GetDataInfo("), blStateCounter("StateCounter("),
    blIntegrate("Integrate("), blLFCompute("LFCompute("), blGetURL("GetURL("),
    blDoSQL("DoSQL("), blTopology("Topology "),
    blAlignSequences("AlignSequences("), blGetNeutralNull("GetNeutralNull("),
    blHBLProfile("#profile"), blDeleteObject("DeleteObject("),
    blRequireVersion("RequireVersion("), blSCFG("SCFG "), blNN("NeuralNet "),
    blBGM("BayesianGraphicalModel "), blSimulateDataSet("SimulateDataSet"),
    blAssert("assert(");

//______________________________________________________________________________
_String _hblCommandAccessor(_ExecutionList *theList, long index) {
  if (theList) {
    if (index >= 0) {
      if (index < theList->lLength) {
        _ElementaryCommand *aCommand =
            (_ElementaryCommand *)theList->GetItem(index);
        return _String((_String *)aCommand->toStr());
      } else {
        return "<END EXECUTION>";
      }
    }
  }
  return _String("command index ") & index;
}

//______________________________________________________________________________
_ElementaryCommand *makeNewCommand(long ccode) {
  _ElementaryCommand *newC = new _ElementaryCommand(ccode);
  checkPointer(newC);
  return newC;
}

//______________________________________________________________________________
void ReadBatchFile(_String &fName, _ExecutionList &target)
    // read/parse a file into an execution list
    // THE function!!!
    {

  fName.ProcessFileName(target.nameSpacePrefix);

  if (terminateExecution) {
    return;
  }
  /*#else
      _Variable optprec (optimizationPrecision);
      _Constant precvalue (0.01);
      FetchVar(LocateVarByName (optimizationPrecision))->SetValue(&precvalue);
  #endif*/

  FILE *f = doFileOpen(fName.getStr(), "rb");
  SetStatusLine("Parsing File");
  if (!f) {
    WarnError(_String("Could not read batch file '") & fName &
              "'.\nPath stack: " & (_String *)pathNames.toStr());
  } else {
    _String beavis(f);

    if (beavis.beginswith("#NEXUS", false)) {
      ReadDataSetFile(f, 1, nil, &fName);
    } else {
      target.BuildList(beavis);
      target.sourceFile = fName;
    }
    fclose(f);
  }
}

//______________________________________________________________________________
void SetDataFilterParameters(_String &parName, _DataSetFilter *thedf, bool setOrKill) {
  _String varName(parName & ".species");
  _Variable *receptacleVar = nil;

  if (setOrKill) {
    setParameter(varName, thedf->NumberSpecies());
  } else {
    DeleteVariable(varName);
  }

  varName = parName & ".sites";
  if (setOrKill) {
    setParameter(varName, thedf->GetFullLengthSpecies() / thedf->GetUnitLength());
  } else {
    DeleteVariable(varName);
  }

  varName = parName & ".unique_sites";
  if (setOrKill) {
    setParameter(varName, thedf->NumberDistinctSites());
  } else {
    DeleteVariable(varName);
  }

  varName = parName & ".site_freqs";
  _Parameter sizeCutoff;
  if (setOrKill) {
    checkParameter(defaultLargeFileCutoff, sizeCutoff, 100000.);
    if (thedf->theFrequencies.lLength < sizeCutoff) {
      receptacleVar = CheckReceptacle(&varName, empty, false);
      receptacleVar->SetValue(new _Matrix(thedf->theFrequencies), false);
    }
  } else {
    DeleteVariable(varName);
  }

  varName = parName & ".site_map";
  if (setOrKill) {
    if (thedf->theOriginalOrder.lLength < sizeCutoff) {
      receptacleVar = CheckReceptacle(&varName, empty, false);
      receptacleVar->SetValue(new _Matrix(thedf->theOriginalOrder), false);
    }
  } else {
    DeleteVariable(varName);
  }

  varName = parName & ".sequence_map";
  if (setOrKill) {
    if (thedf->theOriginalOrder.lLength < sizeCutoff) {
      receptacleVar = CheckReceptacle(&varName, empty, false);
      receptacleVar->SetValue(new _Matrix(thedf->theNodeMap), false);
    }
  } else {
    DeleteVariable(varName);
  }
}

//______________________________________________________________________________
void SerializeModel(_String &rec, long theModel, _AVLList *alreadyDone,
                    bool completeExport) {
  bool mByF = true, do2 = false;

  _Variable *tV = nil, *tV2 = nil;

  _Formula *theExp = nil;
  _SimpleList matrices;

  if (modelTypeList.lData[theModel]) {
    theExp = (_Formula *)modelMatrixIndices.lData[theModel];
    theExp->ScanFForType(matrices, MATRIX);

    for (long mi = 0; mi < matrices.countitems(); mi++) {
      if (alreadyDone &&
          alreadyDone->Insert((BaseRef) matrices.lData[mi]) < 0) {
        matrices.Delete(mi);
        mi--;
      }
    }
  } else {
    if (!alreadyDone || alreadyDone->Find((BaseRef) modelMatrixIndices.lData[theModel]) < 0) {
      if (alreadyDone) {
        alreadyDone->Insert((BaseRef) modelMatrixIndices.lData[theModel]);
      }
      matrices << modelMatrixIndices.lData[theModel];
    }
    tV = LocateVar(modelMatrixIndices.lData[theModel]);
  }

  long freqID = modelFrequenciesIndices.lData[theModel];

  if (freqID >= 0) {
    tV2 = LocateVar(freqID);
  } else {
    mByF = false;
    tV2 = LocateVar(-freqID - 1);
  }

  if (!alreadyDone || alreadyDone->Find((BaseRef) tV2->GetAVariable()) < 0) {
    if (alreadyDone) {
      alreadyDone->Insert((BaseRef) tV2->GetAVariable());
    }
    do2 = true;
  }

  if (completeExport && (matrices.lLength || do2 || theExp)) {
    _SimpleList vl, ind, dep, cat;

    _AVLList vlst(&vl);

    if (theExp) {
      theExp->ScanFForVariables(vlst, true, false, true);
    }

    for (long mi = 0; mi < matrices.lLength; mi++) {
      LocateVar(matrices.lData[mi])->ScanForVariables(vlst, true);
    }

    if (do2) {
      tV2->ScanForVariables(vlst, true);
    }

    vlst.ReorderList();
    SplitVariablesIntoClasses(vl, ind, dep, cat);
    _String glVars(128L, true), locVars(128L, true);
    ExportIndVariables(glVars, locVars, &ind);
    ExportDepVariables(glVars, locVars, &dep);
    glVars.Finalize();
    locVars.Finalize();
    rec << glVars;
    rec << locVars;
    ExportCatVariables(rec, &cat);
  }

  if (matrices.lLength) {
    for (long k = 0; k < matrices.lLength; k++) {
      _Variable *tV = LocateVar(matrices.lData[k]);
      ((_Matrix *)tV->GetValue())->Serialize(rec, *tV->GetName());
      rec << '\n';
    }
  }

  if (do2) {
    ((_Matrix *)tV2->GetValue())->Serialize(rec, *tV2->GetName());
  }

  rec << "\nModel ";
  rec << *((_String *)modelNames(theModel));
  rec << "=(";
  if (theExp) {
    rec << '"';
    rec << _String((_String *)(theExp->toStr()));
    rec << '"';
  } else {
    rec << *tV->GetName();
  }
  rec << ',';
  rec << *tV2->GetName();
  if (theExp) {
    rec << ',';
    rec << explicitFormMExp;
  } else if (!mByF) {
    rec << ",0";
  }
  rec << ");\n";
}

//______________________________________________________________________________
_List   *_HYFetchFunctionParameters (const unsigned long func_id) {
    if (func_id < batchLanguageFunctionParameterLists.lLength) 
        return (_List*)batchLanguageFunctionParameterLists.GetItem(func_id);
    return nil;
}

//______________________________________________________________________________
_ExecutionList   *_HYFetchFunctionBody (const unsigned long func_id) {
    if (func_id < batchLanguageFunctions.lLength) 
        return (_ExecutionList*)batchLanguageFunctions.GetItem(func_id);
    return nil;
}
