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

#include "bayesgraph.h"


//#include "profiler.h"
#ifndef __HEADLESS__
#ifdef __HYPHY_GTK__
#include "HYConsoleWindow.h"
#include "HYDialogs.h"
#include "HYUtils.h"
#include "HYTreePanel.h"
#include "HYDataPanel.h"
#include "HYChartWindow.h"
#include "HYDBWindow.h"
#endif

#ifdef __MAC__
#include <Aliases.h>
//#include <StandardFile.h>
#include <Files.h>
#include "HYDataPanel.h"
#include "HYDialogs.h"
#include "HYTreePanel.h"
#include "HYDataPanel.h"
#include "HYUtils.h"
#include "HYChartWindow.h"
#include "HYDBWindow.h"
#include "timer.h"
void    GetFullPathName (FSSpec& theReply, _String& feedback);
_String MacSimpleFileOpen (void);
_String MacSimpleFileSave (void);
#endif

#ifdef __WINDOZE__
#include "HYDataPanel.h"
#include "HYDialogs.h"
#include "HYTreePanel.h"
#include "HYDataPanel.h"
#include "HYUtils.h"
#include "HYChartWindow.h"
extern  long    lastFileTypeSelection;
#include "HYUtils.h"
#include "HYDBWindow.h"
#endif

#ifndef __UNIX__
#include "HYEventTypes.h"
#endif
#endif

#if defined __HYPHYQT__
#include "HYSharedMain.h"
#include "hyphy_qt_helpers.h"
#endif


#ifdef      __MACPROFILE__
#include "profiler.h"
#endif

#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif


//____________________________________________________________________________________
// global variables

_List
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
openFileHandlesBackend,
standardLibraryPaths,
standardLibraryExtensions,
loadedLibraryPathsBackend;


#ifdef __MAC__
_String volumeName;
#endif


// retrieval functions

_SimpleList
returnlist,
batchLanguageFunctionParameters,
batchLanguageFunctionClassification,
modelMatrixIndices,
modelTypeList,
modelFrequenciesIndices,
listOfCompiledFormulae;

_String
globalPolynomialCap             ("GLOBAL_POLYNOMIAL_CAP"),
                                enforceGlobalPolynomialCap      ("ENFORCE_GLOBAL_POLYNOMIAL_CAP"),
                                dropPolynomialTerms             ("DROP_POLYNOMIAL_TERMS"),
                                maxPolyTermsPerVariable         ("MAX_POLY_TERMS_PER_VARIABLE"),
                                maxPolyExpIterates              ("MAX_POLYNOMIAL_EXP_ITERATES"),
                                polyExpPrecision                ("POLYNOMIAL_EXP_PRECISION"),
                                systemVariableDump              ("LIST_ALL_VARIABLES"),
                                selfDump                        ("PRINT_SELF"),
                                printDigitsSpec                 ("PRINT_DIGITS"),
                                explicitFormMExp                ("EXPLICIT_FORM_MATRIX_EXPONENTIAL"),
                                multByFrequencies               ("MULTIPLY_BY_FREQUENCIES"),
                                getDString                      ("PROMPT_FOR_STRING"),
                                useLastFString                  ("LAST_FILE_PATH"),
                                getFString                      ("PROMPT_FOR_FILE"),
                                tempFString                     ("TEMP_FILE_NAME"),
                                defFileString                   ("DEFAULT_FILE_SAVE_NAME"),
                                useLastModel                    ("USE_LAST_MODEL"),
                                VerbosityLevelString            ("VERBOSITY_LEVEL"),
                                hasEndBeenReached               ("END_OF_FILE"),
                                clearFile                       ("CLEAR_FILE"),
                                keepFileOpen                    ("KEEP_OPEN"),
                                closeFile                       ("CLOSE_FILE"),
                                useLastDefinedMatrix            ("USE_LAST_DEFINED_MATRIX"),
                                MessageLogging                  ("MESSAGE_LOGGING"),
                                selectionStrings                ("SELECTION_STRINGS"),
                                useNoModel                      ("USE_NO_MODEL"),
                                stdoutDestination               ("stdout"),
                                messageLogDestination           ("MESSAGE_LOG"),
                                lastModelParameterList          ("LAST_MODEL_PARAMETER_LIST"),
                                dataPanelSourcePath             ("DATA_PANEL_SOURCE_PATH"),
                                windowTypeTree                  ("TREEWINDOW"),
                                windowTypeClose                 ("CLOSEWINDOW"),
                                windowTypeTable                 ("CHARTWINDOW"),
                                windowTypeDistribTable          ("DISTRIBUTIONWINDOW"),
                                windowTypeDatabase              ("DATABASEWINDOW"),
                                screenWidthVar                  ("SCREEN_WIDTH"),
                                screenHeightVar                 ("SCREEN_HEIGHT"),
                                useNexusFileData                ("USE_NEXUS_FILE_DATA"),
                                mpiMLELFValue                   ("MPI_MLE_LF_VALUE"),
                                lf2SendBack                     ("LIKE_FUNC_NAME_TO_SEND_BACK"),
                                pcAmbiguitiesResolve            ("RESOLVE_AMBIGUITIES"),
                                pcAmbiguitiesAverage            ("AVERAGE_AMBIGUITIES"),
                                pcAmbiguitiesSkip               ("SKIP_AMBIGUITIES"),
                                lfStartCompute                  ("LF_START_COMPUTE"),
                                lfDoneCompute                   ("LF_DONE_COMPUTE"),
                                getURLFileFlag                  ("SAVE_TO_FILE"),
                                versionString                   ("HYPHY_VERSION"),
                                timeStamp                       ("TIME_STAMP"),
                                simulationFilter                ("_SIM_INTERNAL_FILTER_"),
                                prefixDS                        ("DataSet_"),
                                prefixDF                        ("Partition_"),
                                prefixLF                        ("LF_"),
                                replaceTreeStructure            ("REPLACE_TREE_STRUCTURE"),
                                hyphyBaseDirectory              ("HYPHY_BASE_DIRECTORY"),
                                hyphyLibDirectory               ("HYPHY_LIB_DIRECTORY"),
                                platformDirectorySeparator      ("DIRECTORY_SEPARATOR"),
                                covarianceParameterList         ("COVARIANCE_PARAMETER"),
                                matrixEvalCount                 ("MATRIX_EXPONENTIATION_COUNTS"),
                                scfgCorpus                      ("SCFG_STRING_CORPUS"),
                                _hyLastExecutionError           ("LAST_HBL_EXECUTION_ERROR"),
                                _hyExecutionErrorMode           ("HBL_EXECUTION_ERROR_HANDLING"),

                                bgmData                         ("BGM_DATA_MATRIX"),
                                bgmScores                       ("BGM_SCORE_CACHE"),
                                bgmGraph                        ("BGM_GRAPH_MATRIX"),
                                bgmNodeOrder                    ("BGM_NODE_ORDER"),

                                bgmConstraintMx                 ("BGM_CONSTRAINT_MATRIX"),
                                bgmParameters                   ("BGM_NETWORK_PARAMETERS"),

                                pathToCurrentBF                 ("PATH_TO_CURRENT_BF"),
                                hfCountGap                      ("COUNT_GAPS_IN_FREQUENCIES"),
                                gdiDFAtomSize                   ("ATOM_SIZE"),
                                statusBarProgressValue          ("STATUS_BAR_PROGRESS_VALUE"),
                                statusBarUpdateString           ("STATUS_BAR_STATUS_STRING"),
                                marginalAncestors               ("MARGINAL"),
                                doLeavesAncestors               ("DOLEAVES"),
                                blScanfRewind                   ("REWIND"),
                                blFprintfRedirect               ("GLOBAL_FPRINTF_REDIRECT"),
                                blFprintfDevNull                ("/dev/null"),
                                getDataInfoReturnsOnlyTheIndex  ("GET_DATA_INFO_RETURNS_ONLY_THE_INDEX"),
                                dialogPrompt,
                                baseDirectory,
                                lastModelUsed,
                                libDirectory,
                                scanfLastFilePath,
                                defFileNameValue;




#ifdef      __HYPHYMPI__

_String     mpiNodeID                       ("MPI_NODE_ID"),
            mpiNodeCount                    ("MPI_NODE_COUNT"),
            mpiLastSentMsg                  ("MPI_LAST_SENT_MSG");

void        ReportMPIError                  (int, bool);

#endif

bool        isInFunction = false;

_Parameter  explicitFormMatrixExponential = 0.0,
            messageLogFlag                = 1.0;

long        scanfLastReadPosition         = 0;

extern      _String             MATRIX_AGREEMENT,
            ANAL_COMP_FLAG;


extern      _Parameter          toPolyOrNot,
            toMorNot2M,
            ANALYTIC_COMPUTATION_FLAG;

extern      _SimpleList         freeSlots;


_AVLListX   openFileHandles     (&openFileHandlesBackend);
_AVLList    loadedLibraryPaths  (&loadedLibraryPathsBackend);

_ExecutionList
*currentExecutionList = nil;


//____________________________________________________________________________________
// Function prototypes

#ifdef      __HYPHYMPI__

//____________________________________________________________________________________

void    ReportMPIError      (int code, bool send)
{
    if (code != MPI_SUCCESS) {
        _String errMsg = "MPI Error while ";
        if (send) {
            errMsg = errMsg & "sending";
        } else {
            errMsg = errMsg & "receiving";
        }

        errMsg = errMsg & _String(" code:") & (long)code;
        FlagError (errMsg);
    }
}

#define MPI_SEND_CHUNK 0xFFFFFFL

//____________________________________________________________________________________

void    MPISendString       (_String& theMessage, long destID, bool isError)
{

    long    messageLength = theMessage.sLength,
            transferCount = 0;

    if (isError) {
        messageLength = -messageLength;
    }

    ReportMPIError(MPI_Send(&messageLength, 1, MPI_LONG, destID, HYPHY_MPI_SIZE_TAG, MPI_COMM_WORLD),true);

    if (messageLength == 0) {
        return;
    }

    if (isError) {
        messageLength = -messageLength;
    }

    while (messageLength-transferCount>MPI_SEND_CHUNK) {
        ReportMPIError(MPI_Send(theMessage.sData+transferCount, MPI_SEND_CHUNK, MPI_CHAR, destID, HYPHY_MPI_STRING_TAG, MPI_COMM_WORLD),true);
        transferCount += MPI_SEND_CHUNK;
    }

    if (messageLength-transferCount) {
        ReportMPIError(MPI_Send(theMessage.sData+transferCount, messageLength-transferCount, MPI_CHAR, destID, HYPHY_MPI_STRING_TAG, MPI_COMM_WORLD),true);
    }

    //ReportMPIError(MPI_Send(&messageLength, 1, MPI_LONG, destID, HYPHY_MPI_DONE_TAG, MPI_COMM_WORLD),true);

    _FString*    sentVal = new _FString;
    sentVal->theString = (_String*)theMessage.makeDynamic();
    _Variable *   mpiMsgVar = CheckReceptacle (&mpiLastSentMsg, empty, false);
    mpiMsgVar->SetValue (sentVal, false);
    //setParameter (mpiLastSentMsg, &sentVal);

}

//____________________________________________________________________________________
_String*    MPIRecvString       (long senderT, long& senderID)
{
    _String*    theMessage = nil;
    long        messageLength = 0,
                transferCount = 0;

    int         actualReceived = 0;
    bool        isError       = false;

    if  (senderT<0) {
        senderT = MPI_ANY_SOURCE;
    }

    MPI_Status  status;

    ReportMPIError(MPI_Recv(&messageLength, 1, MPI_LONG, senderT, HYPHY_MPI_SIZE_TAG, MPI_COMM_WORLD,&status),false);

    if (messageLength < 0) {
        isError = true;
        messageLength = -messageLength;
    }

    if (!isError) {
        //MPI_Get_count (&status,MPI_CHAR,&actualReceived);

        if (messageLength==0) {
            return new _String;
        }

        theMessage = new _String (messageLength, false);

        senderT = senderID = status.MPI_SOURCE;

        while (messageLength-transferCount>MPI_SEND_CHUNK) {
            ReportMPIError(MPI_Recv(theMessage->sData+transferCount, MPI_SEND_CHUNK, MPI_CHAR, senderT, HYPHY_MPI_STRING_TAG, MPI_COMM_WORLD,&status),false);
            MPI_Get_count (&status,MPI_CHAR,&actualReceived);
            if (actualReceived!=MPI_SEND_CHUNK) {
                WarnError ("Failed in MPIRecvString - some data was not properly received\n");
            }
            //return    nil;
            transferCount += MPI_SEND_CHUNK;
        }

        if (messageLength-transferCount) {
            ReportMPIError(MPI_Recv(theMessage->sData+transferCount, messageLength-transferCount, MPI_CHAR, senderT, HYPHY_MPI_STRING_TAG, MPI_COMM_WORLD,&status),false);
            MPI_Get_count (&status,MPI_CHAR,&actualReceived);
            if (actualReceived!=messageLength-transferCount) {
                WarnError ("Failed in MPIRecvString - some data was not properly received\n");
            }
            //return    nil;
        }
        //ReportMPIError(MPI_Recv(&messageLength, 1, MPI_LONG, senderT, HYPHY_MPI_DONE_TAG, MPI_COMM_WORLD,&status),false);

        if (isError) {
            FlagError (theMessage);
        }
    }
    return theMessage;
}
#endif


//____________________________________________________________________________________

_String GetStringFromFormula (_String* data,_VariableContainer* theP)
{
    _Formula  nameForm (*data,theP);
    _PMathObj formRes = nameForm.Compute();

    if (formRes&& formRes->ObjectClass()==STRING) {
        data = ((_FString*)formRes)->theString;
    }

    return *data;
}

//____________________________________________________________________________________

_String*    ProcessCommandArgument (_String* data)
{
    if ( data->sLength>1 && data->sData[data->sLength-1]=='&' ) {
        _String argName    (*data,0,data->sLength-2);
        _FString* theVar = (_FString*)FetchObjectFromVariableByType(&argName,STRING);
        if (theVar) {
            return theVar->theString;
        }

        WarnError (_String("Reference argument \"")&*data &"\" is not a valid string variable.");
        return nil;
    }
    return (_String*)data;
}

//____________________________________________________________________________________

bool    numericalParameterSuccessFlag = true;

_Parameter  ProcessNumericArgument (_String* data, _VariableContainer* theP, _ExecutionList* currentProgram) {
    _String   errMsg;
    _Formula  nameForm (*data,theP, currentProgram?&errMsg:nil);
     
    if (errMsg.sLength && currentProgram) {
        currentProgram->ReportAnExecutionError (errMsg);
    }
    else {
        _PMathObj formRes = nameForm.Compute();
        numericalParameterSuccessFlag = true;
        if (formRes&& formRes->ObjectClass()==NUMBER) {
            return formRes->Value();
        } else {
            if (formRes&& formRes->ObjectClass()==STRING) {
                return _String((_String*)((_FString*)formRes)->toStr()).toNum();
            } else {
                WarnError (_String("'") & *data & "' was expected to be a numerical argument.");
            }
        }
    }
    numericalParameterSuccessFlag = false;
    return 0.0;
}

//____________________________________________________________________________________

_PMathObj   ProcessAnArgumentByType (_String* expression, _VariableContainer* theP, long objectType)
{
    _Formula  expressionProcessor (*expression, theP);
    _PMathObj expressionResult = expressionProcessor.Compute();
    if (expressionResult && expressionResult->ObjectClass()==objectType) {
        return (_PMathObj)expressionResult->makeDynamic();
    }

    return nil;
}


//____________________________________________________________________________________

_String ProcessLiteralArgument (_String* data, _VariableContainer* theP, _ExecutionList* currentProgram)
{
    _String   errMsg;
    _Formula  nameForm (*data,theP, currentProgram?&errMsg:nil);
    if (errMsg.sLength && currentProgram) {
        currentProgram->ReportAnExecutionError (errMsg);
    } else {
        _PMathObj formRes = nameForm.Compute();
        if (formRes && formRes->ObjectClass()==STRING) {
            return *((_FString*)formRes)->theString;
        }
    }

    return empty;
}

//____________________________________________________________________________________

_AssociativeList*   ProcessDictionaryArgument (_String* data, _VariableContainer* theP, _ExecutionList* currentProgram)
{
    _String   errMsg;
    _Formula  nameForm (*data,theP, currentProgram?&errMsg:nil);
    if (errMsg.sLength && currentProgram) {
        currentProgram->ReportAnExecutionError (errMsg);
    } else {
        _PMathObj formRes = nameForm.Compute();
        if (formRes && formRes->ObjectClass()==ASSOCIATIVE_LIST) {
            formRes->AddAReference();
            return (_AssociativeList*)formRes;
        }
    }
    return nil;
}




//____________________________________________________________________________________
long    FindDataSetName (_String&s)
{
    return dataSetNamesList.Find (&s);
}
//____________________________________________________________________________________
long    FindDataSetFilterName (_String&s)
{
    return dataSetFilterNamesList.Find (&s);
}
//____________________________________________________________________________________
long    FindLikeFuncName (_String&s, bool tryAsAString)
{
    long try1 = likeFuncNamesList.Find (&s);
    if (try1 < 0 && tryAsAString) {
        _String s2 (ProcessLiteralArgument(&s, nil));
        try1 = likeFuncNamesList.Find(&s2);
    }
    return try1;
}

//____________________________________________________________________________________
long    FindModelName (_String&s)
{
    if (s.Equal (&useLastModel)) {
        return lastMatrixDeclared;
    }

    return modelNames.Find (&s);
}

//____________________________________________________________________________________
_LikelihoodFunction*    FindLikeFuncByName (_String&s)
{
    long i = FindLikeFuncName(s);
    if (i>=0) {
        return (_LikelihoodFunction*)likeFuncList (i);
    }
    return nil;
}

//____________________________________________________________________________________
long    FindSCFGName (_String&s)
{
    return scfgNamesList.Find (&s);
}

//____________________________________________________________________________________
long    FindBFFunctionName (_String&s, _VariableContainer* theP)
{
    if (theP) {
        _String testName = *(theP->GetName()) & '.' & s;

        long cutAt = testName.sLength - s.sLength - 2;
        do {
            long idx = batchLanguageFunctionNames.Find (&testName);
            if (idx >= 0) {
                s = testName;
                return idx;
            }
            testName.Trim (0,cutAt);
            cutAt = testName.FindBackwards('.',0,-1)-1;
        } while (cutAt >= 0);
    }

    return batchLanguageFunctionNames.Find (&s);
}


//____________________________________________________________________________________
long    FindBgmName (_String&s)
{
    return bgmNamesList.Find (&s);
}




//__________________________________________________________

long  AddFilterToList (_String& partName,_DataSetFilter* theFilter, bool addP)
{
    FindUnusedObjectName (prefixDF,partName,dataSetFilterNamesList);
    long k;

    for (k=0; k<dataSetFilterNamesList.lLength; k++)
        if (((_String*)dataSetFilterNamesList(k))->sLength==0) {
            break;
        }

    if (addP) {
        SetDataFilterParameters (partName, theFilter, true);
    }

    if (k==dataSetFilterNamesList.lLength) {
        dataSetFilterList<<theFilter;
        DeleteObject (theFilter);
        dataSetFilterNamesList&& & partName;
        return dataSetFilterNamesList.lLength-1;
    }
    dataSetFilterList.lData[k]=(long)theFilter;
    dataSetFilterNamesList.Replace(k,&partName,true);
    return k;
}

//__________________________________________________________
long  AddDataSetToList (_String& theName,_DataSet* theDS)
{
    FindUnusedObjectName (prefixDS,theName,dataSetNamesList);
    long k = dataSetNamesList.Find (&empty);
    if (k==-1) {
        dataSetList.AppendNewInstance (theDS);
        dataSetNamesList&& & theName;
        k = dataSetNamesList.lLength-1;
    } else {
        dataSetNamesList.Replace (k, &theName, true);
        dataSetList.lData[k]     = (long)theDS;
    }
    return k;
}


//__________________________________________________________

void KillDataFilterRecord (long dfID, bool addP)
{
    if (addP) {
        SetDataFilterParameters (*(_String*)(dataSetFilterNamesList(dfID)), nil, false);
    }

    if (dfID<dataSetFilterList.lLength-1) {
        DeleteObject(dataSetFilterList(dfID));
        dataSetFilterList.lData [dfID] = 0;
        dataSetFilterNamesList.Replace(dfID,&empty,true);
    } else {
        dataSetFilterList.Delete(dfID);
        dataSetFilterNamesList.Delete(dfID);
        if (dfID)
            while (((_String*)dataSetFilterNamesList (--dfID))->sLength==0) {
                dataSetFilterList.Delete(dfID);
                dataSetFilterNamesList.Delete(dfID);
                if (dfID==0) {
                    break;
                }
            }
    }
}

//__________________________________________________________

void KillLFRecord (long lfID, bool completeKill)
{
    /* compile the list of variables which will no longer be referenced */

    if (lfID>=0) {
        _LikelihoodFunction *me = (_LikelihoodFunction*)likeFuncList (lfID);

        if (completeKill) {
            _SimpleList         wastedVars,
                                otherVars,
                                myVars,
                                otherModels,
                                wastedModels;

            long                k;


            myVars  << me->GetIndependentVars();
            myVars  << me->GetDependentVars();

            for (k=0; k<likeFuncList.lLength; k++)
                if (k!=lfID) {
                    if (((_String*)likeFuncNamesList(k))->sLength) {
                        _LikelihoodFunction *lf = (_LikelihoodFunction*)likeFuncList (k);
                        otherVars << lf->GetIndependentVars();
                        otherVars << lf->GetDependentVars();
                        for (long kk=lf->GetTheTrees().lLength-1; kk>=0; kk--) {
                            _TheTree * thisTree = (_TheTree*)LocateVar(lf->GetTheTrees().lData[kk]);
                            thisTree->CompileListOfModels (otherModels);
                        }
                    }
                }

            otherVars.Sort();
            otherModels.Sort();

            for (k=0; k<myVars.lLength; k++)
                if (otherVars.BinaryFind(myVars.lData[k])<0) {
                    wastedVars << myVars.lData[k];
                }

            myVars.Clear();

            for (k=me->GetTheTrees().lLength-1; k>=0; k--) {
                _TheTree * thisTree = (_TheTree*)LocateVar(me->GetTheTrees().lData[k]);
                thisTree->CompileListOfModels (myVars);
                _CalcNode * tNode = thisTree->DepthWiseTraversal (true);
                while (tNode) {
                    tNode->SetValue (new _Constant (tNode->BranchLength()),false);
                    tNode = thisTree->DepthWiseTraversal();
                }
                thisTree->RemoveModel();
            }

            for (k=0; k<myVars.lLength; k++)
                if (otherModels.BinaryFind (myVars.lData[k])<0) {
                    KillModelRecord (myVars.lData[k]);
                }

            for (k=0; k<wastedVars.lLength; k++) {
                DeleteVariable (*LocateVar(wastedVars.lData[k])->GetName());
            }

        }

        if (lfID<likeFuncList.lLength-1) {
            DeleteObject(likeFuncList(lfID));
            likeFuncList.lData[lfID] = nil;
            likeFuncNamesList.Replace(lfID,&empty,true);
        } else {
            likeFuncList.Delete(lfID);
            likeFuncNamesList.Delete(lfID);
            if (lfID)
                while (((_String*)likeFuncNamesList (--lfID))->sLength==0) {
                    likeFuncList.Delete(lfID);
                    likeFuncNamesList.Delete(lfID);
                    if (lfID==0) {
                        break;
                    }
                }
        }
    }
}

//__________________________________________________________

void KillLFRecordFull (long lfID)
{
    _LikelihoodFunction* lf = (_LikelihoodFunction*) likeFuncList (lfID);

    long    k;
    //for (k=lf->GetTheFilters().lLength-1; k>=0; k--)
    //  KillDataFilterRecord (lf->GetTheFilters().lData[k]);

    _SimpleList l;

    lf->GetGlobalVars (l);

    for (k=0; k<l.lLength; k++) {
        DeleteVariable (*LocateVar(l.lData[k])->GetName());
    }

    l.Clear ();

    for (k=lf->GetTheTrees().lLength-1; k>=0; k--) {
        _TheTree * thisTree = (_TheTree*)LocateVar(lf->GetTheTrees().lData[k]);
        thisTree->CompileListOfModels (l);
        DeleteVariable (*thisTree->GetName());
    }

    for (k=0; k<l.lLength; k++) {
        KillModelRecord (l.lData[k]);
    }

    KillLFRecord (lfID);
}

//__________________________________________________________

void KillDataSetRecord (long dsID)
{
    if (dsID<dataSetList.lLength-1) {
        DeleteObject(dataSetList(dsID));
        dataSetList.lData[dsID] = 0;
        dataSetNamesList.Replace(dsID,&empty,true);
    } else {
        dataSetList.Delete(dsID);
        dataSetNamesList.Delete(dsID);
        if (dsID)
            while (((_String*)dataSetNamesList (--dsID))->sLength==0) {
                dataSetList.Delete(dsID);
                dataSetNamesList.Delete(dsID);
                if (dsID==0) {
                    break;
                }
            }
    }
}

//__________________________________________________________

void    KillExplicitModelFormulae (void)
{
    for (long i = 0; i < modelTypeList.lLength; i++)
        if (modelTypeList.lData[i]) {
            delete (_Formula*)(modelMatrixIndices.lData[i]);
        }
}


//__________________________________________________________

void KillModelRecord (long mdID)
{
    if (lastMatrixDeclared == mdID) {
        lastMatrixDeclared = -1;
    }

    // SLKP 20110816
    // can't delete matrices before checking that no other model depends on the


    if (modelTypeList.lData[mdID]) {
        delete (_Formula*)(modelMatrixIndices.lData[mdID]);
    } else {
        _Variable * modelMatrix = nil,
                    * freqMatrix  = nil;

        bool      multByFreqs   = false;



        _SimpleList  saveTheseVariablesAux;
        _AVLList     saveTheseVariables (&saveTheseVariablesAux);

        for (long k = 0; k < modelNames.lLength; k++) {
            if (k != mdID && ((_String*)modelNames(k))->sLength) {
                if (modelTypeList.lData[k]) {
                    _SimpleList dependantMatrices;
                    ((_Formula*)(modelMatrixIndices.lData[k]))->ScanFForType (dependantMatrices, MATRIX);
                    for (long k2 = 0; k2 < dependantMatrices.lLength; k2++) {
                        saveTheseVariables.Insert((BaseRef)dependantMatrices.lData[k2]);
                    }
                } else {
                    RetrieveModelComponents(k, modelMatrix, freqMatrix, multByFreqs);

                    if (modelMatrix) {
                        saveTheseVariables.Insert((BaseRef)modelMatrix->GetIndex());
                    }
                    if (freqMatrix) {
                        saveTheseVariables.Insert((BaseRef)freqMatrix->GetIndex());
                    }
                }
            }
        }

        RetrieveModelComponents(mdID, modelMatrix, freqMatrix, multByFreqs);
        if (modelMatrix && saveTheseVariables.Find ((BaseRef)modelMatrix->GetIndex()) < 0) {
            DeleteVariable (*modelMatrix->GetName());
        }
        if (freqMatrix && saveTheseVariables.Find ((BaseRef)freqMatrix->GetIndex()) < 0) {
            DeleteVariable (*freqMatrix->GetName());
        }

        //DeleteVariable (*LocateVar(mID)->GetName());
    }




    if (mdID<modelNames.lLength-1) {
        modelMatrixIndices.lData[mdID] = -1;
        modelTypeList.lData[mdID] = 0;
        modelFrequenciesIndices.lData[mdID] = -1;
        modelNames.Replace(mdID,&empty,true);
    } else {
        modelNames.Delete(mdID);
        modelMatrixIndices.Delete (modelMatrixIndices.lLength-1);
        modelFrequenciesIndices.Delete (modelFrequenciesIndices.lLength-1);
        modelTypeList.Delete (modelTypeList.lLength-1);

        if (mdID)
            while (((_String*)modelNames (--mdID))->sLength==0) {
                modelNames.Delete(mdID);
                modelMatrixIndices.Delete (mdID);
                modelFrequenciesIndices.Delete (mdID);
                modelTypeList.Delete (mdID);
                if (mdID == 0) {
                    break;
                }
            }
    }
}

//____________________________________________________________________________________
//____________________________________________________________________________________
_ExecutionList::_ExecutionList ()
{
    result         = nil;
    currentCommand = 0;
    cli            = nil;
    profileCounter = nil;
    stdinRedirect  = nil;
    stdinRedirectAux = nil;
    doProfile      = 0;
    nameSpacePrefix = nil;
    if (currentExecutionList) {
        errorHandlingMode  = currentExecutionList->errorHandlingMode;
    } else {
        errorHandlingMode = HY_BL_ERROR_HANDLING_DEFAULT;
    }

} // doesn't do much

//____________________________________________________________________________________
_ExecutionList::_ExecutionList (_String& source, _String* namespaceID, bool copySource, bool* successFlag)
{
    currentCommand = 0;
    result         = nil;
    cli            = nil;
    profileCounter = nil;
    doProfile      = 0;
    stdinRedirect  = nil;
    stdinRedirectAux = nil;
    nameSpacePrefix = nil;
    
    if (namespaceID) {
        SetNameSpace (*namespaceID);
    }
    if (copySource) {
        sourceText.Duplicate (&source);
    }
    if (currentExecutionList) {
        errorHandlingMode  = currentExecutionList->errorHandlingMode;
    } else {
        errorHandlingMode = HY_BL_ERROR_HANDLING_DEFAULT;
    }
    bool result = BuildList (source, nil, false, true);
    if (successFlag) {
        *successFlag = result;
    }
}

//____________________________________________________________________________________

_ExecutionList::~_ExecutionList (void)
{
    if (cli) {
        delete [] cli->values;
        delete [] cli->stack;
        delete cli;
        cli = nil;
    }

    if (profileCounter) {
        DeleteObject (profileCounter);
        profileCounter = nil;
    }

    DeleteObject (stdinRedirect);
    DeleteObject (stdinRedirectAux);
    DeleteObject (nameSpacePrefix);

    ResetFormulae();
    DeleteObject (result);
}

//____________________________________________________________________________________

BaseRef     _ExecutionList::makeDynamic (void)
{
    _ExecutionList * Res = (_ExecutionList*)checkPointer(new _ExecutionList);

    memcpy ((char*)Res, (char*)this, sizeof (_ExecutionList));

    Res->nInstances         = 1;
    Res->Duplicate          (this);
    Res->cli                = nil;
    Res->profileCounter     = nil;
    Res->doProfile          = doProfile;
    Res->errorHandlingMode  = errorHandlingMode;

    if(result) {
        Res->result = (_PMathObj)result->makeDynamic();
    }

    return Res;
}

//____________________________________________________________________________________

void        _ExecutionList::Duplicate   (BaseRef source)
{
    _List::Duplicate    (source);

    _ExecutionList* s = (_ExecutionList*)source;

    if (s->result) {
        s->result=(_PMathObj)result->makeDynamic();
    }

    errorHandlingMode  = s->errorHandlingMode;
}


//____________________________________________________________________________________
void    _ExecutionList::ReportAnExecutionError (_String errMsg, bool doCurrentCommand, bool appendToExisting) {
    if (doCurrentCommand) {
        _ElementaryCommand *theCommand = FetchLastCommand();
        if (theCommand) {
            errMsg = errMsg & " in call to " & _HY_ValidHBLExpressions.RetrieveKeyByPayload(theCommand->GetCode());
        }
    }
    switch (errorHandlingMode) {
        case HY_BL_ERROR_HANDLING_SOFT:
            if (appendToExisting) {
              _FString * existing = (_FString*) FetchObjectFromVariableByType(&_hyLastExecutionError, STRING);
              if (existing) {
                errMsg = *existing->theString & '\n' & errMsg;
              }
            }
            setParameter(_hyLastExecutionError, new _FString (errMsg, false), false);
            
            break;
        default: 
            WarnError (errMsg);
    }
}

//____________________________________________________________________________________
_String*    _ExecutionList::FetchFromStdinRedirect (void)
// grab a string from the front of the input queue
// complain if nothing is left
{
    if (!stdinRedirect) {
        WarnError ("No input buffer was given for a redirected standard input read.");
        return new _String;
    }
    long d = stdinRedirect->First();
    if (d<0) {
        WarnError ("Ran out of input in buffer during a redirected standard input read.");
        return new _String;
    }
    _String *sendBack = (_String*)stdinRedirect->GetXtra (d);
    sendBack->nInstances++;
    stdinRedirect->Delete ((*(_List*)stdinRedirect->dataList)(d),true);
    return sendBack;
}

//____________________________________________________________________________________

_String       _ExecutionList::GetFileName     (void)  {
    if (sourceFile.sLength) {
        return sourceFile;
    } else {
        if (pathNames.lLength)
            return *(_String*)pathNames.GetElement (-1);
    }
    return empty;
}
// doesn't do much
//____________________________________________________________________________________

_PMathObj       _ExecutionList::Execute     (void)      // run this execution list
{

    setParameter(_hyLastExecutionError, new _MathObject, false);
    
    _ExecutionList*      stashCEL = currentExecutionList;
    callPoints << currentCommand;
    executionStack       << this;

    _String             dd (GetPlatformDirectoryChar());

    _FString            bp  (baseDirectory, false),
                        lp  (libDirectory, false),
                        ds  (dd),
                        cfp (pathNames.lLength?*(_String*)pathNames(pathNames.lLength-1):empty),
                        * stashed = (_FString*)FetchObjectFromVariableByType (&pathToCurrentBF, STRING);

    setParameter        (platformDirectorySeparator, &ds);
    setParameter        (hyphyBaseDirectory, &bp);
    setParameter        (hyphyLibDirectory, &lp);

    if (stashed) {
        stashed = (_FString*)stashed->makeDynamic();
    }
    setParameter        (pathToCurrentBF,&cfp);

    DeleteObject        (result);
    result               = nil;
    currentExecutionList = this;
    currentCommand       = 0;

    terminateExecution  = false;
    skipWarningMessages = false;

    while (currentCommand<lLength) {
        if (doProfile == 1 && profileCounter) {
            long        instCounter = currentCommand;
            _Parameter  timeDiff    = 0.0;

            TimerDifferenceFunction (false);
            (((_ElementaryCommand**)lData)[currentCommand])->Execute(*this);
            timeDiff   = TimerDifferenceFunction(true);

            if (profileCounter) {
                profileCounter->theData[instCounter*2]   += timeDiff;
                profileCounter->theData[instCounter*2+1] += 1.0;
            }
        } else {
            (((_ElementaryCommand**)lData)[currentCommand])->Execute(*this);
        }

        if (terminateExecution) {
            break;
        }
    }
    currentCommand = callPoints.lData[callPoints.lLength-1];
    callPoints.Delete (callPoints.lLength-1);
    currentExecutionList = stashCEL;

    if (stashed) {
        setParameter        (pathToCurrentBF,stashed,false);
    }

    executionStack.Delete (executionStack.lLength-1);

    return result;
}

//____________________________________________________________________________________

long        _ExecutionList::ExecuteAndClean     (long g, _String* fName)        // run this execution list
{
    long    f = -1;
    Execute ();

    if (fName && !terminateExecution) {
        f = batchLanguageFunctionNames.Find (fName);
    }

    terminateExecution      = false;
    skipWarningMessages     = false;

    while (g<batchLanguageFunctionNames.lLength) {
        batchLanguageFunctionNames.Delete           (g);
        batchLanguageFunctionParameters.Delete      (g);
        batchLanguageFunctions.Delete               (g);
        batchLanguageFunctionClassification.Delete  (g);
        batchLanguageFunctionParameterLists.Delete  (g);
    }
    return f;
}

//____________________________________________________________________________________

bool        _ExecutionList::TryToMakeSimple     (void)
{
    _SimpleList     varList,
                    formulaeToConvert,
                    parseCodes;

    long            stackDepth  = 0;

    bool            status      = true;

    for (long k = 0; k<lLength && status; k++) {
        _ElementaryCommand * aStatement = (_ElementaryCommand*)(*this)(k);
        switch (aStatement->code) {
        case 0: {
            _String * formulaString = (_String*)aStatement->parameters(0);

            if (formulaString->sData[formulaString->sLength-1]!='}') {
                _Formula *f  = new _Formula,
                *f2 = new _Formula;

                checkPointer ((BaseRef)(f&&f2));

                long          varRef,
                              parseCode = Parse(f,*formulaString,varRef,nameSpacePrefix,f2);

                if (parseCode == HY_FORMULA_EXPRESSION || parseCode == HY_FORMULA_VARIABLE_VALUE_ASSIGNMENT) {
                    if (f->AmISimple(stackDepth,varList)) {
                        aStatement->simpleParameters<<parseCode;
                        aStatement->simpleParameters<<(long)f;
                        aStatement->simpleParameters<<(long)f2;
                        aStatement->simpleParameters<<varRef;

                        formulaeToConvert << (long)f;

                        if (HY_FORMULA_VARIABLE_VALUE_ASSIGNMENT) {
                            parseCodes        << varRef;
                        } else {
                            parseCodes        << -1;
                        }
                        break;
                    }
                }

                delete f;
                delete f2;
            }
            status = false;
            break;
        }

        case 4:
            parseCodes        << -1;
            if (aStatement->simpleParameters.lLength == 3 || aStatement->parameters.lLength) {
                if (aStatement->parameters.lLength) {
                    _Formula f;
                    long     varRef;
                    //printf ("Namespace: %x\nCode: %s\n", chain.nameSpacePrefix, ((_String*)parameters(0))->sData);

                    long status = Parse (&f, *(_String*)aStatement->parameters(0), varRef, nameSpacePrefix,nil);

                    //printf ("Print formula: %s\n", _String((_String*)f.toStr()).sData);

                    if (status== HY_FORMULA_EXPRESSION) {
                        aStatement->simpleParameters<<long(f.makeDynamic());
                    }
                }


                _Formula *cf  = ((_Formula*)aStatement->simpleParameters(2));
                if (cf->AmISimple(stackDepth,varList)) {
                    formulaeToConvert << (long)cf;
                } else {
                    status = false;
                }
            }
            break;

        default:
            status = false;
        }
        if (status == false) {
            ReportWarning (_String ("Failed to compile an execution list: offending command was\n") & _String (((_String*)aStatement->toStr())));
        }
    }

    if (status) {
        cli = new _CELInternals;
        checkPointer (cli);
        checkPointer (cli->values = new _SimpleFormulaDatum[varList.lLength+1]);
        checkPointer (cli->stack  = new _SimpleFormulaDatum[stackDepth+1]);

        _SimpleList  avlData;
        _AVLListX    avlList (&avlData);

        for (long fi = 0; fi < formulaeToConvert.lLength; fi++) {
            ((_Formula*)formulaeToConvert(fi))->ConvertToSimple (varList);
        }

        for (long vi = 0; vi < varList.lLength; vi++) {
            avlList.Insert ((BaseRef)varList.lData[vi], vi);
        }

        for (long ri = 0; ri<parseCodes.lLength; ri++) {
            if (parseCodes.lData[ri] < 0) {
                cli->storeResults << -1;
            } else {
                cli->storeResults << avlList.GetXtra (avlList.Find ((BaseRef) parseCodes.lData[ri]));
            }
        }
        cli->varList.Duplicate(&varList);
    }

    return status;
}

//____________________________________________________________________________________

void        _ExecutionList::ExecuteSimple       (void)
{
    PopulateArraysForASimpleFormula (cli->varList, cli->values);
    Execute();

    for (long vi2 = 0; vi2 < cli->varList.lLength; vi2++) {
        _Variable * mv = LocateVar(cli->varList.lData[vi2]);
        if (mv->ObjectClass() == NUMBER) {
            mv->SetValue (new _Constant (cli->values[vi2].value),false);
        }
    }
}

//____________________________________________________________________________________

void        _ExecutionList::ResetFormulae       (void)      // run this execution list
{
    currentCommand = 0;
    while (currentCommand<lLength) {
        _ElementaryCommand* thisCommand = ((_ElementaryCommand**)lData)[currentCommand];
        if (thisCommand->code==0) {
            if (thisCommand->simpleParameters.lLength) {
                //printf ("[ResetFormulae] %s\n", thisCommand->sData);
                _Formula* f = (_Formula*)
                              thisCommand->simpleParameters.lData[1],
                              *f2 = (_Formula*)
                                    thisCommand->simpleParameters.lData[2] ;
                if (f) {
                    delete f;
                }
                if (f2) {
                    delete f2;
                }
                thisCommand->simpleParameters.Clear();
                long k = listOfCompiledFormulae.Find((long)thisCommand);
                if (k >= 0) {
                    listOfCompiledFormulae.Delete(k);
                    //printf ("[ResetFormulae:listOfCompiledFormulae %d]\n",k);
                    compiledFormulaeParameters.Delete(k);
                    //printf ("[ResetFormulae:compiledFormulaeParameters %d]\n",k);
                }
            }
        } else {
            if (thisCommand->code==4) {
                if (thisCommand->parameters.lLength && thisCommand->simpleParameters.lLength == 3) {
                    _Formula* f = (_Formula*)thisCommand->simpleParameters.lData[2];
                    if (f) {
                        delete f;
                    }
                    thisCommand->simpleParameters.Delete (2);
                }
            }
        }
        currentCommand++;
    }
}
//____________________________________________________________________________________

BaseRef  _ExecutionList::toStr (void)
{
    _String *result = new _String (1,true),
    step ("\n\nStep"),
    dot (".");

    for (long i=0; i<countitems(); i++) {
        (*result) << &step;
        _String lineNumber (i);
        (*result)<< &lineNumber;
        (*result)<< '.';
        result->AppendNewInstance ((_String*)(*this)(i)->toStr());
    }
    result->Finalize();
    return result;
}

//____________________________________________________________________________________

void     _ExecutionList::ResetNameSpace (void)
{
    DeleteObject (nameSpacePrefix);
    nameSpacePrefix = nil;
}

void     _ExecutionList::SetNameSpace (_String nID)
{
    ResetNameSpace ();
    nameSpacePrefix = new _VariableContainer(nID);
    checkPointer(nameSpacePrefix);
}

_String*     _ExecutionList::GetNameSpace ()
{
    if (nameSpacePrefix) {
        return nameSpacePrefix->GetName();
    }
    return nil;
}

_String  _ExecutionList::AddNameSpaceToID (_String& theID)
{
    if (nameSpacePrefix) {
        return (*nameSpacePrefix->GetName())&'.'&theID;
    }
    return theID;
}

_String  _ExecutionList::TrimNameSpaceFromID (_String& theID)
{
    if (nameSpacePrefix) {
        if (theID.startswith(*nameSpacePrefix->GetName())) {
            return theID.Cut(nameSpacePrefix->GetName()->sLength+1,-1);
        }
    }
    return theID;
}


_String  blFor                  ("for("),               // moved
         blWhile                    ("while("),         // moved
         blFunction                 ("function "),      // moved
         blFFunction                ("ffunction "),     // moved
         blReturn                   ("return "),        // moved
         blReturn2              ("return("),            // moved
         blIf                       ("if("),            // moved
         blElse                     ("else"),           // moved
         blDo                       ("do{"),            // moved
         blBreak                    ("break;"),         // moved
         blContinue             ("continue;"),          // moved
         blInclude              ("#include"),           // moved
         blDataSet              ("DataSet "),           // moved
         blDataSetFilter            ("DataSetFilter "),
         blConstructCM          ("ConstructCategoryMatrix("),
         blTree                     ("Tree "),
         blLF                       ("LikelihoodFunction "),
         blLF3                  ("LikelihoodFunction3 "),
         blMolClock                 ("MolecularClock("),
         blfprintf              ("fprintf("),
         blGetString                ("GetString("),
         blfscanf                   ("fscanf("),
         blsscanf                   ("sscanf("),
         blExport                   ("Export("),
         blReplicate                ("ReplicateConstraint("),
         blImport                   ("Import"),
         blCategory             ("category "),
         blClearConstraints         ("ClearConstraints("),
         blSetDialogPrompt      ("SetDialogPrompt("),
         blModel                    ("Model "),
         blChoiceList               ("ChoiceList("),
         blOpenDataPanel            ("OpenDataPanel("),
         blGetInformation           ("GetInformation("),
         blExecuteCommands      ("ExecuteCommands("),
         blExecuteAFile         ("ExecuteAFile("),
         blLoadFunctionLibrary      ("LoadFunctionLibrary("),
         blOpenWindow               ("OpenWindow("),
         blSpawnLF                  ("SpawnLikelihoodFunction("),
         blDifferentiate            ("Differentiate("),
         blFindRoot             ("FindRoot("),
         blMPIReceive               ("MPIReceive("),
         blMPISend                  ("MPISend("),
         blGetDataInfo              ("GetDataInfo("),
         blStateCounter             ("StateCounter("),
         blIntegrate                ("Integrate("),
         blLFCompute                ("LFCompute("),
         blGetURL                   ("GetURL("),
         blDoSQL                    ("DoSQL("),
         blTopology                 ("Topology "),
         blAlignSequences           ("AlignSequences("),
         blGetNeutralNull           ("GetNeutralNull("),
         blHBLProfile               ("#profile"),
         blDeleteObject         ("DeleteObject("),
         blRequireVersion           ("RequireVersion("),
         blSCFG                     ("SCFG "),
         blNN                       ("NeuralNet "),
         blBGM                      ("BayesianGraphicalModel "),
         blSimulateDataSet          ("SimulateDataSet"),
         blAssert                   ("assert(");



bool        _ExecutionList::BuildList   (_String& s, _SimpleList* bc, bool processed, bool empty_is_success)
{
    if (terminateExecution) {
        return false;
    }

    char * savePointer = s.sData;
    
    _SimpleList          triePath;

    while (s.Length()) { // repeat while there is stuff left in the buffer
        _String currentLine (_ElementaryCommand::FindNextCommand (s,true));

        if (currentLine.getChar(0)=='}') {
            currentLine.Trim(1,-1);
        }

        if (!currentLine.Length()) {
            continue;
        }
        
        triePath.Clear(false);
        long prefixTreeCode = _HY_ValidHBLExpressions.Find (currentLine, &triePath, true);
        
        _List *pieces = nil;
        _HBLCommandExtras *commandExtraInfo = nil;
        
        if (prefixTreeCode != HY_TRIE_NOTFOUND) {
            prefixTreeCode = _HY_ValidHBLExpressions.GetValue(prefixTreeCode);
            long commandExtra = _HY_HBLCommandHelper.FindLong (prefixTreeCode);
            if (commandExtra >= 0) { // pre-trim all strings as needed
                commandExtraInfo = (_HBLCommandExtras*)_HY_HBLCommandHelper.GetXtra (commandExtra);
                if (commandExtraInfo->extract_conditions.lLength > 0) {
                    pieces = new _List;
                    long upto = _ElementaryCommand::ExtractConditions (currentLine, commandExtraInfo->cut_string,*pieces,commandExtraInfo->extract_condition_separator),
                         condition_index_match = commandExtraInfo->extract_conditions.Find(pieces->lLength);
                    if (condition_index_match < 0) {
                        // try to see if the command accepts a variable number of arguments (at least X)
                       _String parseFail;
                       if (commandExtraInfo->extract_conditions.lLength == 1 && commandExtraInfo->extract_conditions.lData[0] < 0) {
                            if (pieces->lLength < -commandExtraInfo->extract_conditions.lData[0]) {
                                 parseFail = _String("Incorrect number of arguments (") & (long) pieces->lLength & ") supplied: expected at least " & _String (-commandExtraInfo->extract_conditions.lData[0]) & ", while processing '"& currentLine.Cut (0, upto) & "'. ";
                             }
                        } else {
                            parseFail = _String("Incorrect number of arguments (") & (long) pieces->lLength & ") supplied: expected one of " & _String ((_String*)commandExtraInfo->extract_conditions.toStr()) & ", while processing '"& currentLine.Cut (0, upto) & "'. ";
                        }
                        if (parseFail.sLength) {
                            if (currentExecutionList) {
                                currentExecutionList->ReportAnExecutionError(parseFail, false, true);
                            } else {
                                acknError (parseFail);
                            }
                            DeleteObject (pieces);
                            return false;  
                        }
                    }
                    if (commandExtraInfo->do_trim) {
                        currentLine.Trim (upto, -1);
                    }
                }
            }
        }
        
        bool handled = false;
               
        switch (prefixTreeCode) {
            case HY_HBL_COMMAND_FOR:
                _ElementaryCommand::BuildFor (currentLine, *this, *pieces);
                handled = true;
                break;
            case HY_HBL_COMMAND_WHILE:
                _ElementaryCommand::BuildWhile (currentLine, *this, *pieces);
                handled = true;
                break;
            case HY_HBL_COMMAND_BREAK:
            case HY_HBL_COMMAND_CONTINUE:
                if (bc) {
                    AppendNewInstance(new _ElementaryCommand);
                    (*bc) << ((prefixTreeCode == HY_HBL_COMMAND_BREAK) ? (countitems()-1) : (-(long)countitems()+1));
                } else {
                    WarnError (currentLine & " only makes sense in the context of a loop.");
                    return false;
                }
                handled = true;
                break;
            case HY_HBL_COMMAND_SET_DIALOG_PROMPT:
            case HY_HBL_COMMAND_HARVEST_FREQUENCIES:
            case HY_HBL_COMMAND_OPTIMIZE:
            case HY_HBL_COMMAND_COVARIANCE_MATRIX:
            case HY_HBL_COMMAND_LFCOMPUTE:
            case HY_HBL_COMMAND_SELECT_TEMPLATE_MODEL:
            case HY_HBL_COMMAND_USE_MODEL:
            case HY_HBL_COMMAND_SET_PARAMETER:
            case HY_HBL_COMMAND_ASSERT:
            case HY_HBL_COMMAND_REQUIRE_VERSION:
            case HY_HBL_COMMAND_DELETE_OBJECT:
            case HY_HBL_COMMAND_CLEAR_CONSTRAINTS:
            case HY_HBL_COMMAND_MOLECULAR_CLOCK:
            case HY_HBL_COMMAND_GET_URL:
            case HY_HBL_COMMAND_GET_STRING:
            case HY_HBL_COMMAND_EXPORT:
            case HY_HBL_COMMAND_DIFFERENTIATE:
                _ElementaryCommand::ExtractValidateAddHBLCommand (currentLine, prefixTreeCode, pieces, commandExtraInfo, *this);
                handled = true;
                break;
                
        }
        
        if (handled)
            DeleteObject (pieces);
        
        // 20111212: this horrendous switch statement should be replaced with a 
        // prefix tree lookup 

        if (!handled) {
            if (currentLine.startswith (blFunction)||currentLine.startswith (blFFunction)) { // function declaration
                _ElementaryCommand::ConstructFunction (currentLine, *this);
            } else if (currentLine.startswith (blReturn) || currentLine.startswith (blReturn2)) { // function return statement
                _ElementaryCommand::ConstructReturn (currentLine, *this);
            } else if (currentLine.startswith (blIf)) { // if-then-else statement
                _ElementaryCommand::BuildIfThenElse (currentLine, *this, bc);
            } else if (currentLine.startswith (blElse)) { // else clause of an if-then-else statement
                if (lastif.countitems()) {
                    long    temp = countitems(),
                            lc   = lastif.countitems(),
                            lif  = lastif.lData[lc-1];

                    _ElementaryCommand      * stuff = new _ElementaryCommand ();
                    stuff->MakeJumpCommand  (nil,0,0,*this);
                    AppendNewInstance       (stuff);
                    currentLine.Trim        (4,-1);

                    long  index         = currentLine.Length()-1,
                          scopeIn     = 0;

                    while (currentLine.sData[scopeIn]=='{' && currentLine.sData[index]=='}') {
                        scopeIn++;
                        index--;
                    }

                    if (scopeIn) {
                        currentLine.Trim (scopeIn,index);
                    }

                    BuildList (currentLine,bc,true);

                    if (lif<0 || lif>=lLength) {
                        WarnError ("'else' w/o an if to latch on to...");
                        return false;
                    }

                    ((_ElementaryCommand*)((*this)(lif)))->MakeJumpCommand(nil,-1,temp+1,*this);
                    ((_ElementaryCommand*)(*this)(temp))->simpleParameters[0]=countitems();

                    while (lastif.countitems()>=lc) {
                        lastif.Delete(lastif.countitems()-1);
                    }
                } else {
                    WarnError ("'else' w/o an if to latch on to...");
                    return false;
                }

            } else if (currentLine.startswith (blDo)) { // do {} while statement
                _ElementaryCommand::BuildDoWhile (currentLine, *this);
            }  else if (currentLine.startswith (blInclude)) { // #include
                _ElementaryCommand::ProcessInclude (currentLine, *this);
            } else if (currentLine.startswith (blDataSet)) { // data set definition
                _ElementaryCommand::ConstructDataSet (currentLine, *this);
            } else if (currentLine.startswith (blDataSetFilter)) { // data set filter definition
                _ElementaryCommand::ConstructDataSetFilter (currentLine, *this);
            } else if (currentLine.startswith (blConstructCM)) { // construct category assignments matrix
                _ElementaryCommand::ConstructCategoryMatrix (currentLine, *this);
            } else if (currentLine.startswith (blTree) || currentLine.startswith (blTopology)) { // tree definition
                _ElementaryCommand::ConstructTree (currentLine, *this);
            } else if (currentLine.startswith (blLF) || currentLine.startswith (blLF3)) { // LF definition
                _ElementaryCommand::ConstructLF (currentLine, *this);
            } else if (currentLine.startswith (blfprintf)) { // fpintf call
                _ElementaryCommand::ConstructFprintf (currentLine, *this);
            } else if (currentLine.startswith (blfscanf) || currentLine.startswith (blsscanf)) { // fscanf call
                _ElementaryCommand::ConstructFscanf (currentLine, *this);
            } else if (currentLine.startswith (blReplicate)) { // replicate constraint statement
                _ElementaryCommand::ConstructReplicateConstraint (currentLine, *this);
            } else if (currentLine.startswith (blCategory)) { // category variable declaration
                _ElementaryCommand::ConstructCategory (currentLine, *this);
            } else if (currentLine.startswith (blGetNeutralNull)) { // select a template model
                _ElementaryCommand::ConstructGetNeutralNull (currentLine, *this);
            } else if (currentLine.startswith (blModel)) { // Model declaration
                _ElementaryCommand::ConstructModel (currentLine, *this);
            } else if (currentLine.startswith (blChoiceList)) { // choice list
                _ElementaryCommand::ConstructChoiceList (currentLine, *this);
            } else if (currentLine.startswith (blOpenDataPanel)) { // open data panel window
                _ElementaryCommand::ConstructOpenDataPanel (currentLine, *this);
            } else if (currentLine.startswith (blGetInformation)) { // get information
                _ElementaryCommand::ConstructGetInformation (currentLine, *this);
            } else if (currentLine.startswith (blExecuteCommands) || currentLine.startswith (blExecuteAFile) || currentLine.startswith (blLoadFunctionLibrary))
                // execute commands
            {
                _ElementaryCommand::ConstructExecuteCommands (currentLine, *this);
            } else if (currentLine.startswith (blOpenWindow)) { // execute commands
                _ElementaryCommand::ConstructOpenWindow (currentLine, *this);
            } else if (currentLine.startswith (blSpawnLF)) { // execute commands
                _ElementaryCommand::ConstructSpawnLF (currentLine, *this);
            } else if (currentLine.startswith (blFindRoot)||currentLine.startswith (blIntegrate))
                // find a root of an expression in an interval
                // or an integral
            {
                _ElementaryCommand::ConstructFindRoot (currentLine, *this);
            } else if (currentLine.startswith (blMPISend)) { // MPI Send
                _ElementaryCommand::ConstructMPISend (currentLine, *this);
            } else if (currentLine.startswith (blMPIReceive)) { // MPI Receive
                _ElementaryCommand::ConstructMPIReceive (currentLine, *this);
            } else if (currentLine.startswith (blGetDataInfo)) { // Get Data Info
                _ElementaryCommand::ConstructGetDataInfo (currentLine, *this);
            } else if (currentLine.startswith (blStateCounter)) { // Get Data Info
                _ElementaryCommand::ConstructStateCounter (currentLine, *this);
            } else if (currentLine.startswith (blDoSQL)) { // Do SQL
                _ElementaryCommand::ConstructDoSQL (currentLine, *this);
            } else if (currentLine.startswith (blAlignSequences)) { // Do AlignSequences
                _ElementaryCommand::ConstructAlignSequences (currentLine, *this);
            } else if (currentLine.startswith (blHBLProfile)) { // #profile
                _ElementaryCommand::ConstructProfileStatement (currentLine, *this);
            } else if (currentLine.startswith (blSCFG)) { // SCFG definition
                _ElementaryCommand::ConstructSCFG (currentLine, *this);
            } else if (currentLine.startswith (blNN)) { // Neural Net definition
                _ElementaryCommand::ConstructNN (currentLine, *this);
            } else if (currentLine.startswith (blBGM)) {    // Bayesian Graphical Model definition
                _ElementaryCommand::ConstructBGM (currentLine, *this);
            } 
            // plain ol' formula - parse it as such!
            else {
                _String checker (currentLine);
                if (_ElementaryCommand::FindNextCommand (checker).Length()==currentLine.Length()) {
                    if (currentLine.Length()>1)
                        while (currentLine[currentLine.Length()-1]==';') {
                            currentLine.Trim (0,currentLine.Length()-2);
                        }
                    else {
                        continue;
                    }
                    _ElementaryCommand* oddCommand = new _ElementaryCommand(currentLine);
                    oddCommand->code = 0;
                    oddCommand->parameters&&(&currentLine);
                    AppendNewInstance (oddCommand);
                } else {
                    while (currentLine.Length()) {
                        _String part (_ElementaryCommand::FindNextCommand (currentLine));
                        BuildList (part,bc,processed);
                    }
                }
            }
            
            /*if (currentLine.sLength > 1 || currentLine.sLength == 1 && currentLine.getChar(0) != ';'){
                WarnError (_String ("Missing semicolon before ") & currentLine);
                return false;
            }*/
        }
    }
    s.sData = savePointer;
    s.DuplicateErasing (&empty);
    return empty_is_success || countitems();
}

//____________________________________________________________________________________

_ElementaryCommand::_ElementaryCommand (void)
{
    code = -1;
}

//____________________________________________________________________________________

_ElementaryCommand::_ElementaryCommand (long ccode)
{
    code = ccode;
}

//____________________________________________________________________________________

_ElementaryCommand::_ElementaryCommand (_String& s)
{
    code = -1;
    _String::Duplicate (&s);
}

//____________________________________________________________________________________
_ElementaryCommand::~_ElementaryCommand (void)
{
    if (nInstances==1) {
        if (code==4) {
            if (simpleParameters.lLength>2) {
                _Formula* f = (_Formula*)simpleParameters(2);
                delete (f);
            }
        } else if (code==0) {
            if (simpleParameters.lLength) {
                //long k = listOfCompiledFormulae.Find ((long)this);
                //listOfCompiledFormulae.Delete (k);
                //compiledFormulaeParameters.Delete (k);

                _Formula* f = (_Formula*)simpleParameters(2);
                delete (f);
                f = (_Formula*)simpleParameters(1);
                delete(f);
                simpleParameters.Clear();
            }
        } else if ((code==6)||(code==9)) {
            for (long i = 0; i<simpleParameters.lLength; i++) {
                _Formula* f = (_Formula*)simpleParameters(i);
                delete (f);
            }
        }
    }

}

//____________________________________________________________________________________
BaseRef   _ElementaryCommand::makeDynamic (void)
{
    _ElementaryCommand * nec = new _ElementaryCommand;
    if (!nec) {
        return nil;
    }
    nec->code = code;

    //memcpy ((char*)Res, (char*)this, sizeof (_ElementaryCommand));
    //if (code == 4 || code==6 || code==9 || code==0)
    //nInstances++;

    nec->nInstances = 1;
    nec->Duplicate (this);
    return nec;
}

//____________________________________________________________________________________

void      _ElementaryCommand::Duplicate (BaseRef source)
{
    _ElementaryCommand* sec = (_ElementaryCommand*)source;

    //if (code == 0)
    //_String::DuplicateErasing (source);
    //else
    _String::Duplicate (source);

    parameters.Duplicate(&sec->parameters);
    if (code != 0) {
        simpleParameters.Duplicate(&sec->simpleParameters);
    }
}

//____________________________________________________________________________________

_String _hblCommandAccessor (_ExecutionList* theList, long index) {
    if (theList) {
        _ElementaryCommand * aCommand = (_ElementaryCommand*)theList->GetItem (index);
        if (aCommand) {
            return _String ((_String*)aCommand->toStr());
        }
    }
    return _String ("command index ") & index;
}

//____________________________________________________________________________________

BaseRef   _ElementaryCommand::toStr      (void)
{
    _String result, *converted = nil;
    long k;
    switch (code) {

    case 0: // formula reparser
        converted = (_String*)((parameters(0))->toStr());
        result    = *converted;
        break;

    case 4:

        result = "Branch ";
        if (simpleParameters.countitems()==3) {
            converted = (_String*)((_Formula*)simpleParameters(2))->toStr();
            result = result& "under condition '"& *converted&"'\n\tto "&
                        _hblCommandAccessor (currentExecutionList,simpleParameters(0))&
                        "\n\telse "&
                        _hblCommandAccessor (currentExecutionList,simpleParameters(1));
        } else {
            result = result&"to "& _hblCommandAccessor (currentExecutionList,simpleParameters(0));
        }

        break;

    case 5: // data set contruction
        converted = (_String*)parameters(0)->toStr();
        result = _String("Read Data Set ")&(*converted);
        DeleteObject(converted);
        converted = (_String*)parameters(1)->toStr();
        result = result& _String(" from file ")&(*converted);

        break;

    case 6:  // data set filter construction
    case HY_HBL_COMMAND_HARVEST_FREQUENCIES:  // or HarvestFrequencies
        if (code == HY_HBL_COMMAND_HARVEST_FREQUENCIES) {
            result = "Harvest Frequencies into matrix ";
            converted = (_String*)(parameters(0)->toStr());
            result = result&*converted;
            DeleteObject (converted);
            converted = (_String*)(parameters(1)->toStr());
            result = result&" from DataSet "&*converted;
        } else {
            result = "Build Data Set Filter ";
            converted = (_String*)(parameters(0)->toStr());
            result = result&*converted;
            DeleteObject (converted);
            converted = (_String*)(parameters(1)->toStr());
            result = result&(_String)(" from DataSet ")&*converted;
        }

        DeleteObject (converted);

        if (parameters.lLength > 2) {
            result = result & " with unit size = " & *(_String*)parameters(2);
        }

        if (code == 9) {

            result = result & " with atom size = " & *(_String*)parameters(3);
            result = result & " with position specific flag = " & *(_String*)parameters(4);
        }

        result = result&" Partition:";
        for (k = code==6?3:5; k<parameters.countitems(); k++) {
            result = result & "," & *(_String*)parameters(k);
        }

        converted = nil;

        break;

    case 7: // build a tree
    case 54: // build a tree

        converted = (_String*)parameters(0)->toStr();
        if (code == 7) {
            result = _String("Build Tree ")&*converted;
        } else {
            result = _String("Build Topology ")&*converted;
        }

        DeleteObject (converted);

        converted = (_String*)parameters(1)->toStr();
        result = result & " from string " &*converted;
        break;

    case 8: // print stuff to file (or stdout)

        converted = (_String*)parameters(0)->toStr();
        result = _String("Print to File ")& (*converted) & ":";
        DeleteObject (converted);

        converted = nil;
        for (k = 1; k<parameters.countitems(); k++) {
            result = result&*((_String*)parameters(k))&"\n";
        }


        break;

    case HY_HBL_COMMAND_OPTIMIZE: // optimize the likelihood function
    case HY_HBL_COMMAND_COVARIANCE_MATRIX: // compute the covariance matrix
        converted = (_String*)parameters(0)->toStr();
        if (code == HY_HBL_COMMAND_OPTIMIZE) {
            result = _String("Optimize storing into, ");
        } else {
            result = _String("Calculate the Covariance Matrix storing into, ");
        }

        result = result& (*converted) & ", the following likelihood function:";
        DeleteObject (converted);

        converted = nil;
        for (k = 1; k<parameters.countitems(); k++) {
            result = result&*((_String*)parameters(k))&" ; ";
        }

        break;

    case 11: // build the likelihood function

        converted = (_String*)parameters(0)->toStr();
        result = _String("Construct the following likelihood function:");
        DeleteObject (converted);

        converted = nil;
        for (k = 1; k<parameters.countitems(); k++) {
            result = result&*((_String*)parameters(k))&" ; ";
        }

        break;

    case 12: // data set simulation
        converted = (_String*)parameters(0)->toStr();
        result = _String("Simulate Data Set ")&(*converted);
        DeleteObject(converted);
        converted = (_String*)parameters(1)->toStr();
        result = result& _String(" from the likelihood function ")&(*converted);

        break;

    case 13: // a function
        converted = (_String*)parameters(0)->toStr();
        result = _String("Function ")&(*converted);
        DeleteObject(converted);
        converted = new _String (simpleParameters(0));
        checkPointer(converted);
        result = result& _String(" of ")&(*converted)&_String(" parameters:{");
        DeleteObject (converted);
        converted = (_String*)parameters(simpleParameters(0)-1)->toStr();
        result = result&_String("}");


        break;

    case 14: // return statement
        result = _String("A return statement with:");
        converted = new _String (simpleParameters(0));
        checkPointer(converted);
        result = result& *converted;

        break;

    case 16: // data set merger
        converted = (_String*)parameters(0)->toStr();
        result = _String("Build dataset")&(*converted)&_String(" by ");
        if (abs(simpleParameters(0)==1)) {
            result = result & _String (" concatenating ");
        } else {
            result = result & _String (" combining ");
        }
        if (simpleParameters (0)<0) {
            result = result & _String ("(deleting arguments upon completion)");
        }
        for (k=1; k<parameters.countitems(); k++) {
            result = result & *((_String*)parameters(k));
            if (k<parameters.countitems()-1) {
                result = result & _String(",");
            }
        }
        break;

    case HY_HBL_COMMAND_EXPORT:
        converted = (_String*)parameters(1)->toStr();
        result = _String("Export ")&(*converted);
        DeleteObject(converted);
        converted = (_String*)parameters(0)->toStr();
        checkPointer(converted);
        result = result& _String(" to ")& *converted;
        break;

    case HY_HBL_COMMAND_MOLECULAR_CLOCK: // a call to MolecularClock
        converted = (_String*)parameters(0)->toStr();
        result    = ",";
        result = _String("Molecular clock imposed starting at ")&(*converted)
               & " on variables " & _String((_String*)parameters.Join (&result, 1, -1));
        break;

    case 20: // category variable construction
        converted = (_String*)parameters.toStr();
        result = _String("Category variable: ")&(*converted);
        break;

    case 21: // optimize the likelihood function

    {
        converted = (_String*)parameters(0)->toStr();
        result = _String("Construct the category matrix, ")& (*converted) & ", for the likelihood function:";
        DeleteObject (converted);
        result = result&*((_String*)parameters(1))&" in ";
        converted = nil;
        _String testAgainst ("COMPLETE");
        if (parameters.countitems()>2) {
            if (((_String*)parameters(2))->Equal(&testAgainst)) {
                result = result & "complete format";
                break;
            }
        }
        result = result & "short format";


        break;
    }
    case HY_HBL_COMMAND_CLEAR_CONSTRAINTS: { // clear constraints
        converted = (_String*)parameters.toStr();
        result = _String("Clear contstraints on: ")&(*converted);
        break;
    }
    case HY_HBL_COMMAND_SET_DIALOG_PROMPT: { // set dialog prompt
        converted = (_String*)parameters.toStr();
        result = _String("Set dialog prompt to: ")&(*converted);
        break;
    }
    case 24: { // select standard model
        converted = (_String*)parameters.toStr();
        result = _String("Select Template Model for ")&(*converted);
        break;
    }
    case 25: // fscanf
    case 56: { // sscanf
        converted = (_String*)parameters(0);
        result = _String ("Read the following data from ")&*converted&" : ";
        long shift = 1;
        for (long p = 0; p<simpleParameters.lLength; p++) {
            long theFormat = simpleParameters(p);
            if (theFormat < 0) {
                shift -- ;
                result = result & " REWIND THE FILE";
            } else {
                result = result&*(_String*)allowedFormats (theFormat)
                         &" into "&*(_String*)parameters(p+shift)&";";
            }
        }
        converted = nil;
        break;
    }
    case 26: { // replicate constraint
        converted = (_String*)parameters(0);
        result = _String ("Replicate the following constraint ")&*converted&" using these variables:";
        for (long p = 1; p<parameters.lLength; p++) {
            result = result&*(_String*)parameters(p)&";";
        }
        converted = nil;
        break;
    }

    case HY_HBL_COMMAND_USE_MODEL: { // use matrix
        converted = (_String*)parameters(0);
        result = _String ("Use the matrix ")&*converted;
        converted = nil;
        break;
    }

    case 31: { // define a model
        converted = (_String*)parameters(0);
        result = _String ("Define the model ")&*converted;
        converted = (_String*)parameters(1);
        result = result& _String (" using transition matrix '")&*converted&"'";
        if (parameters.lLength>2) {
            converted = (_String*)parameters(2);
            result = result& _String (" and equilibrium frequencies '")&*converted&"'";
        }
        converted = nil;
        break;
    }

    case 32: { // choice list
        converted = (_String*)parameters(1)->toStr();
        result = _String ("Choice List for ")&*converted;
        DeleteObject (converted);
        converted = (_String*)parameters(3)->toStr();
        result = result& _String (" with choice list:")&*converted;
        DeleteObject (converted);
        converted = (_String*)parameters(0);
        result = result& _String (". Store result in ")&*converted;
        converted = nil;
        break;
    }

    case HY_HBL_COMMAND_GET_STRING: { // get string from object
        converted = (_String*)parameters(2)->toStr();
        result = _String ("Get string ")&*converted;
        DeleteObject (converted);
        converted = (_String*)parameters(1)->toStr();
        result = result& _String (" from object:")&*converted;
        DeleteObject (converted);
        converted = (_String*)parameters(0);
        result = result& _String (". Store result in ")&*converted;
        converted = nil;
        break;
    }
    case HY_HBL_COMMAND_SET_PARAMETER: { // set parameter value
        converted = (_String*)parameters(1)->toStr();
        result = _String ("Set parameter ")&*converted;
        DeleteObject (converted);
        converted = (_String*)parameters(0)->toStr();
        result = result& _String (" of ")&*converted;
        DeleteObject (converted);
        converted = (_String*)parameters(2);
        result = result& _String (" to ")&*converted;
        converted = nil;
        break;
    }
    case 36: { // open data panel
        converted = (_String*)parameters(0)->toStr();
        result = _String ("Open data window for the data set ")&*converted;
        DeleteObject (converted);
        converted = (_String*)parameters(1)->toStr();
        result = result& _String (" list species ")&*converted;
        DeleteObject (converted);
        converted = (_String*)parameters(2);
        result = result& _String (" with window settings ")&*converted;
        converted = (_String*)parameters(3);
        result = result& _String (" with partition settings ")&*converted;
        converted = nil;
        break;
    }
    case 37: { // open data panel
        converted = (_String*)parameters(0)->toStr();
        result = _String ("Get Information for ")&*converted;
        DeleteObject (converted);
        converted = (_String*)parameters(1)->toStr();
        result = result& _String (" storing in ")&*converted;
        break;
    }

    case 38: { // reconsruct ancestors
        converted = (_String*)parameters(0)->toStr();
        result = _String("Reconstruct Ancestors into ")&(*converted);
        DeleteObject(converted);
        converted = (_String*)parameters(1)->toStr();
        result = result& _String(" from the likelihood function ")&(*converted);
        break;
    }

    case 39: // execute commands
    case 62: // execute a file
    case 66: { // load function library
        converted = (_String*)parameters(0)->toStr();
        if (code == 39) {
            result = _String("ExecuteCommands in string ");
        } else if (code == 62) {
            result = _String("ExecuteAFile from file ");
        } else {
            result = _String("LoadFunctionLibrary from file");
        }

        result = result & *converted & _String(" using basepath ") & ((_String*)parameters(1)->toStr()) & '.';
        if (simpleParameters.lLength) {
            result = result & "\nCompiled.";
        } else if (parameters.lLength>2) {
            _String inputName ((_String*)parameters(2)->toStr());
            result = result & " reading input from " & inputName;
            _AssociativeList *inputValues = (_AssociativeList *)FetchObjectFromVariableByType(&inputName, ASSOCIATIVE_LIST);
            if (inputValues) {
                result = result & '\n' & _String ((_String*)inputValues->toStr());
            }

            if (parameters.lLength > 3) {
                result = result & " using name space prefix " & _String ((_String*)parameters(3)->toStr());
            }
        }
        break;
    }

    case 40: { // open window
        converted = (_String*)parameters(0)->toStr();
        result = _String("Open window of type ")&(*converted);
        DeleteObject(converted);
        converted = (_String*)parameters(1)->toStr();
        result = result& _String(" with parameters from ")&(*converted);
        break;
    }

    case 41: { // spawn LF
        converted = (_String*)parameters(0)->toStr();
        result = _String("Spawn Likelihood Function ")&(*converted);
        DeleteObject(converted);
        converted = (_String*)parameters(1)->toStr();
        result = result& _String(" with tree ")&(*converted);
        DeleteObject(converted);
        converted = (_String*)parameters(2)->toStr();
        result = result& _String(" from ")&(*converted);
        DeleteObject(converted);
        converted = (_String*)parameters(3)->toStr();
        result = result& _String(" using sequences ids in ")&(*converted);
        break;
    }

    case HY_HBL_COMMAND_DIFFERENTIATE: { // Differentiate
        converted = (_String*)parameters(1)->toStr();
        result = _String("Differentiate '")&(*converted);
        DeleteObject(converted);
        converted = (_String*)parameters(2)->toStr();
        result = result& _String("' on ")&(*converted);
        DeleteObject(converted);
        if (parameters.lLength==4) {
            converted = (_String*)parameters(3)->toStr();
            result = result& _String(" ")&(*converted) & " times ";
            DeleteObject(converted);
        }
        converted = (_String*)parameters(2)->toStr();
        result = result& _String(" storing result in ")&(*converted);
        break;
    }
    case 43: // Find Root
    case 48: { // Integrate
        converted = (_String*)parameters(1)->toStr();
        if (code == 43) {
            result = _String("Find the root of ")&(*converted);
        } else {
            result = _String("Integrate ")&(*converted);
        }
        DeleteObject(converted);
        converted = (_String*)parameters(2)->toStr();
        result = result& _String(" in ")&(*converted);
        DeleteObject(converted);
        converted = (_String*)parameters(3)->toStr();
        result = result& _String(" between ")&(*converted);
        DeleteObject(converted);
        converted = (_String*)parameters(4)->toStr();
        result = result& _String(" and ")&(*converted);
        DeleteObject(converted);
        converted = (_String*)parameters(0)->toStr();
        result = result& _String(" storing the result in ")&(*converted);
        break;
    }
    case 44: { // MPISend
        converted = (_String*)parameters(1)->toStr();
        result = _String("MPI Send ")&(*converted);
        DeleteObject(converted);
        if (parameters.lLength>2) {
            converted = (_String*)parameters(2)->toStr();
            result = _String(" with input from ")&(*converted);
            DeleteObject(converted);
        }
        converted = (_String*)parameters(0)->toStr();
        result = result& _String(" to MPI Node ")&(*converted);
        break;
    }
    case 45: { // MPIReceive
        converted = (_String*)parameters(0)->toStr();
        result = _String("MPI Receive from ")&(*converted);
        DeleteObject(converted);
        converted = (_String*)parameters(1)->toStr();
        result = result& _String(" storing actual sender node into ")&(*converted);
        DeleteObject(converted);
        converted = (_String*)parameters(2)->toStr();
        result = result& _String(" and storing the string result into ")&(*converted);
        break;
    }
    case 46: { //GetDataInfo
        converted = (_String*)parameters(1)->toStr();
        if (parameters.lLength>4) {
            result = _String("Compute pairwise differences from ")&(*converted);
        } else {
            result = _String("GetDataInfo from ")&(*converted);
        }
        DeleteObject(converted);
        if (parameters.lLength>2) {
            converted = (_String*)parameters(2)->toStr();
            result = result& _String(" for sequence ")&(*converted);
            DeleteObject(converted);
            if (parameters.lLength>3) {
                converted = (_String*)parameters(3)->toStr();
                if (parameters.lLength>4) {
                    result = result& _String(" and site ")&(*converted);
                } else {
                    result = result& _String(" and sequence ")&(*converted);
                }
            }
            DeleteObject(converted);
        }
        converted = (_String*)parameters(0)->toStr();
        result = result& _String(" and storing the result in ")&(*converted);
        break;
    }
    case 47: { //GetDataInfo
        converted = (_String*)parameters(0)->toStr();
        result = _String("StateCounter on ")&(*converted);
        DeleteObject(converted);
        converted = (_String*)parameters(1)->toStr();
        result = result& _String(" using the callback: ")&(*converted);
        break;
    }
    case HY_HBL_COMMAND_LFCOMPUTE: { //Compute LF
        converted = (_String*)parameters(0)->toStr();
        result = blLFCompute&(*converted);
        DeleteObject(converted);
        converted = (_String*)parameters(1)->toStr();
        result = result& ',' &(*converted) & ')';
        break;
    }
    case HY_HBL_COMMAND_GET_URL: { //GetURL
        converted = (_String*)parameters(0)->toStr();
        result = blGetURL&(*converted);
        DeleteObject(converted);
        converted = (_String*)parameters(1)->toStr();
        result = result& ',' &(*converted);
        if (parameters.lLength>2) {
            DeleteObject(converted);
            converted = (_String*)parameters(2)->toStr();
            result = result& ',' &(*converted);
        }
        result = result & ')';
        break;
    }
    case 52: { //Simulate
        converted = (_String*)parameters(0)->toStr();
        result = blDataSet &(*converted) & "=Simulate(";
        DeleteObject(converted);
        converted = (_String*)parameters(1)->toStr();
        result = result &(*converted);
        for (long i=2; i<parameters.lLength; i++) {
            DeleteObject(converted);
            converted = (_String*)parameters(i)->toStr();
            result = result &','&(*converted);
        }
        result = result & ')';
        break;
    }
    case 53: //DoSQL
    case 55: //AlignSequences
    case 57: { //GetNeutralNull
        converted = (_String*)parameters(0)->toStr();
        result = (code==53?blDoSQL:(code==55?blAlignSequences:blGetNeutralNull)) &(*converted) & '(';
        for (long i=1; i<parameters.lLength; i++) {
            DeleteObject(converted);
            converted = (_String*)parameters(i)->toStr();
            result = result &','&(*converted);
        }
        result = result & ')';
        break;
    }
    case 58: {
        converted = (_String*)parameters(0)->toStr();
        result = blHBLProfile & " " & *converted;
        break;
    }
    case HY_HBL_COMMAND_DELETE_OBJECT: {
        converted = (_String*)parameters.toStr();
        result = blDeleteObject & '(' & *converted & ')';
        break;
    }
    case HY_HBL_COMMAND_REQUIRE_VERSION: {
        converted = (_String*)parameters(0)->toStr();
        result = blRequireVersion & '(' & *converted & ')';
        break;
    }
    case 61:
    case 63: {
        converted = (_String*)parameters(0)->toStr();
        result = (code==61?blSCFG:blNN) & *converted & "=(";
        for (long i=1; i<parameters.lLength; i++) {
            DeleteObject(converted);
            converted = (_String*)parameters(i)->toStr();
            if (i>1) {
                result = result &','&(*converted);
            } else {
                result = result &(*converted);
            }
        }
        result = result &')';
        break;
    }
    case 64:
        converted = (_String *) parameters(0)->toStr();
        result = blBGM & *converted & "=(";
        for (long p = 1; p < parameters.lLength; p++) {
            DeleteObject (converted);
            converted = (_String *) parameters(p)->toStr();
            if (p > 1) {
                result = result & ',' & (*converted);
            } else {
                result = result & (*converted);    // first argument
            }
        }
        result = result & ')';
        break;
    case HY_HBL_COMMAND_ASSERT: {
        converted = (_String*)parameters(0)->toStr();
        result = _String ("Assert ") & "'" & *converted & "'";
        break;
    }
    }

    DeleteObject (converted);
    return result.makeDynamic();
}

//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase0 (_ExecutionList& chain)
{
    chain.currentCommand++;

    if (chain.cli) {
        _Parameter result = ((_Formula*)simpleParameters.lData[1])->ComputeSimple (chain.cli->stack, chain.cli->values);
        long sti = chain.cli->storeResults.lData[chain.currentCommand-1];
        if (sti>=0) {
            chain.cli->values[sti].value = result;
        }
        return;
    }

    if (!simpleParameters.lLength) { // not compiled yet
        _Formula f,
                 f2;

        _String* theFla     = (_String*)parameters(0),
                 errMsg;

        bool    doNotCompileThisFormula = false;

        long    varRef,
                parseCode = Parse(&f,(*theFla),varRef,chain.nameSpacePrefix,&f2,nil,&doNotCompileThisFormula);

        if (parseCode != HY_FORMULA_FAILED ) {
            if (doNotCompileThisFormula == false) { // not a matrix constant
                simpleParameters    <<parseCode;
                simpleParameters    <<long (f.makeDynamic());
                simpleParameters    <<long (f2.makeDynamic());
                simpleParameters    <<varRef;

                _SimpleList*        varList = new _SimpleList;
                _AVLList            varListA (varList);
                f.ScanFForVariables (varListA, true, true, true, true);
                f2.ScanFForVariables(varListA, true, true);
                varListA.ReorderList();
                listOfCompiledFormulae<<(long)this;
                compiledFormulaeParameters.AppendNewInstance(varList);
            } else {
                ExecuteFormula(&f,&f2,parseCode,varRef,chain.nameSpacePrefix);
                return;
            }
        } else {
            return;
        }
    }

    ExecuteFormula ((_Formula*)simpleParameters.lData[1],(_Formula*)simpleParameters.lData[2],simpleParameters.lData[0],simpleParameters.lData[3], chain.nameSpacePrefix);

    if (terminateExecution) {
        WarnError (_String("Problem occurred in line: ")&*this);
        return;
    }
}


//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase4 (_ExecutionList& chain)
{
    chain.currentCommand++;
    if (simpleParameters.lLength==2) {

    }

    if (simpleParameters.lLength==3 || parameters.lLength) {
        if ( parameters.lLength && simpleParameters.lLength < 3) {
            _Formula f;
            long     varRef;
            //printf ("Namespace: %x\nCode: %s\n", chain.nameSpacePrefix, ((_String*)parameters(0))->sData);

            long status = Parse (&f, *(_String*)parameters(0), varRef, chain.nameSpacePrefix,nil);

            //printf ("Print formula: %s\n", _String((_String*)f.toStr()).sData);

            if (status== HY_FORMULA_EXPRESSION) {
                simpleParameters<<long(f.makeDynamic());
            } else {
                return;
            }
        }

        if (chain.cli) {
            if ( ((_Formula*)simpleParameters(2))->ComputeSimple(chain.cli->stack, chain.cli->values)==0.0) {
                chain.currentCommand = simpleParameters.lData[1];
                return;
            }
        } else {
            _PMathObj result = ((_Formula*)simpleParameters(2))->Compute();
            if (!result) {
                WarnError ("Condition Evaluation Failed");
                return ;
            }

            if (terminateExecution) {
                subNumericValues = 2;
                _String       *s = (_String*)((_Formula*)simpleParameters(2))->toStr();
                subNumericValues = 0;
                _String     err  = _String("Failed while evaluating: ") & _String((_String*)((_Formula*)simpleParameters(2))->toStr()) & " - " & *s;
                DeleteObject (s);
                WarnError    (err);
                return;
            }

            bool conditionFalse = false;

            switch (result->ObjectClass()) {
            case NUMBER:
                conditionFalse = result->Value()==0.0;
                break;
            case STRING:
                conditionFalse = ((_FString*)result)->IsEmpty();
                break;
            default:
                WarnError ("Condition evaluation result be be a number or a string");
                return;

            }

            if (conditionFalse) {
                chain.currentCommand = simpleParameters.lData[1];
                return;
            }
        }
    }
    chain.currentCommand = simpleParameters.lData[0];

    if (chain.currentCommand == -1) {
        terminateExecution   = true;
        chain.currentCommand = chain.lLength;
    }
}

//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase5 (_ExecutionList& chain)
{
    chain.currentCommand++;
    FILE*    df;
    _String  fName = *(_String*)parameters(1);
    _DataSet*ds;

    if (simpleParameters.lLength == 1) {
        fName = GetStringFromFormula ((_String*)parameters(1),chain.nameSpacePrefix);
        ds = ReadDataSetFile (nil,0,&fName,nil,chain.nameSpacePrefix?chain.nameSpacePrefix->GetName():nil);
    } else {
        if (fName.Equal (&useNexusFileData)) {
            if (!lastNexusDataMatrix) {
                _String errMsg = useNexusFileData & " was used in ReadDataFile, and no NEXUS data matrix was available.";
                acknError (errMsg);
                return;
            }
            ds = lastNexusDataMatrix;
        } else {
#if defined  __WINDOZE__ && ! defined __HEADLESS__
            lastFileTypeSelection = 1;
#endif

            fName.ProcessFileName(false,false,(Ptr)chain.nameSpacePrefix);
            if (terminateExecution) {
                return;
            }
            SetStatusLine ("Loading Data");

            df = doFileOpen (fName.getStr(),"rb");
            if (df==nil) {
                // try reading this file as a string formula
                fName = GetStringFromFormula ((_String*)parameters(1),chain.nameSpacePrefix);
                fName.ProcessFileName(false,false,(Ptr)chain.nameSpacePrefix);

                if (terminateExecution) {
                    return;
                }

                df = doFileOpen (fName.getStr(),"rb");
                if (df==nil) {
                    _String errMsg ("Could not find source dataset file:");
                    errMsg = errMsg & *(_String*)parameters(1) & " Path stack: " & _String((_String*)pathNames.toStr());
                    WarnError (errMsg);
                    return;
                }
            }
            ds = ReadDataSetFile (df,0,nil,nil,chain.nameSpacePrefix?chain.nameSpacePrefix->GetName():nil);
            fclose (df);
        }
    }


    // 20110802: need to check that this data set is not empty

    if (ds->NoOfSpecies() && ds->NoOfColumns()) {
        _String  * dsID = new _String (chain.AddNameSpaceToID(*(_String*)parameters(0)));
        StoreADataSet (ds, dsID);
        DeleteObject  (dsID);
    } else {
        DeleteObject (ds);
        WarnError    ("The format of the sequence file has not been recognized and may be invalid");
    }

    //StoreADataSet (ds, (_String*)parameters(0));
}

//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase8 (_ExecutionList& chain)
{
    chain.currentCommand++;

    _String* targetName = (_String*)parameters(0),
             fnm;


    FILE*    dest = nil;

    int     out2    = 0;
    bool    doClose = true,
            skipFilePathEval = false;

    if (targetName->Equal(&stdoutDestination)) {
        _FString * redirect = (_FString*)FetchObjectFromVariableByType (&blFprintfRedirect, STRING);
        if (redirect && redirect->theString->sLength) {
            if (redirect->theString->Equal (&blFprintfDevNull)) {
                return;
            } else {
                skipFilePathEval = true;
                targetName = redirect->theString;
            }
        } else {
            out2=1;
        }
    }

    checkParameter (printDigitsSpec,printDigits,0);

    if (!out2) {
        fnm = *targetName;
        if (fnm.Equal(&messageLogDestination)) {
            if ((dest = globalMessageFile) == nil) {
                return;
            }
        } else {
            if (skipFilePathEval == false && !fnm.IsALiteralArgument()) {
                fnm = GetStringFromFormula (&fnm,chain.nameSpacePrefix);
            }

            fnm.ProcessFileName(true,false,(Ptr)chain.nameSpacePrefix);
            if (terminateExecution) {
                return;
            }

            long k  = openFileHandles.Find (&fnm);
            doClose = k<0;

            if (!doClose) {
                dest = (FILE*)openFileHandles.GetXtra (k);
            } else {
                if ((dest = doFileOpen (fnm.getStr(), "a")) == nil) {
                    WarnError (_String  ("Could not create/open output file at path '") & fnm & "'.");
                    return;
                }
            }
        }
    }

    long literalshift = simpleParameters.lLength?0:-1;
    for (long i = 1; i<parameters.lLength; i++) {
        if (literalshift>=0) {
            if (i==simpleParameters(literalshift)) {
                if (!out2) {
                    fprintf (dest,"%s", ((_String*)parameters(i))->getStr());
                } else {
                    StringToConsole (*(_String*)parameters(i));
                }

                if (literalshift<simpleParameters.lLength-1) {
                    literalshift++;
                } else {
                    literalshift = -1;
                }
                continue;
            }
        }

        _String *varname = ProcessCommandArgument((_String*)parameters(i));

        if (!varname) {
            return;
        }

        BaseRef    thePrintObject   =   nil;
        _Formula   f;

        if (varname->Equal(&clearFile)) {
            if (!out2 && dest) {
                fclose (dest);
                dest = doFileOpen (fnm.getStr(), "w");
                long k = openFileHandles.Find (&fnm);
                if (k>=0) {
                    openFileHandles.SetXtra(k, (long)dest);
                }
            }
        } else if (varname->Equal(&keepFileOpen) && !out2) {
            if (openFileHandles.Find (&fnm) < 0) {
                openFileHandles.Insert (fnm.makeDynamic(), (long)dest);
            }

            doClose = false;
        } else if (varname->Equal(&closeFile) && !out2) {
            openFileHandles.Delete (&fnm, true);
            doClose = true;
        } else if (varname->Equal(&systemVariableDump)) {
            thePrintObject=&variableNames;
        } else if (varname->Equal(&selfDump)) {
            thePrintObject=&chain;
        } else {
            // check for possible string reference

            _String    temp    = ProcessStringArgument (varname),
                       nmspace;

            if (temp.sLength > 0) {
                nmspace = AppendContainerName(temp,chain.nameSpacePrefix);
                if (temp.IsValidIdentifier()) {
                    thePrintObject = FetchObjectFromVariableByType (&nmspace,HY_ANY_OBJECT);
                }
            } else {
                nmspace = AppendContainerName(*varname,chain.nameSpacePrefix);
            }


            if (thePrintObject == nil) {
                long typeFlag = HY_BL_ANY;

                thePrintObject = _HYRetrieveBLObjectByName (nmspace, typeFlag);

                if (!thePrintObject) {
                    long    varRef = -1;
                    _String argCopy = *varname;
                    if (Parse (&f,argCopy, varRef, chain.nameSpacePrefix,nil) == HY_FORMULA_EXPRESSION) {
                        thePrintObject = f.Compute();
                    }
                }
            }
        }

        if (thePrintObject) {
            if (!out2) {
                thePrintObject->toFileStr (dest);
            } else {
                _String outS ((_String*)thePrintObject->toStr());
                StringToConsole (outS);
            }
        }
    }
#if !defined __UNIX__ || defined __HEADLESS__ || defined __HYPHYQT__
    if (!dest) {
        yieldCPUTime();
    }
#endif
    if (dest && dest!=globalMessageFile && doClose) {
        fclose (dest);
    }
}

//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase11 (_ExecutionList& chain)
/*
 code cleanup SLKP 20090316
*/

{
    chain.currentCommand++;

    _String  parm,
             errMsg;

    bool     explicitFreqs = simpleParameters.lLength,
             // if false then the likelihood function will be of the form (filter1,tree1,filter2,tree2,...,filterN,treeN)
             // if true then the likelihood function will be of the form  (filter1,tree1,freq1,filter2,tree2,freq2,...,filterN,treeN,freqN)
             assumeList    = parameters.lLength > 2;
    // if there is only one parameter to the function constructor, it is assumed to be a string matrix
    // otherwise it is expected to be a collection of literals

    _List    * likelihoodFunctionSpec   = nil,
               passThisToLFConstructor;

    if (assumeList) {
        likelihoodFunctionSpec = new _List (parameters, 1, - 1);
    } else {
        _Matrix * matrixOfStrings = (_Matrix*)ProcessAnArgumentByType ((_String*)parameters(1), chain.nameSpacePrefix, MATRIX);
        if (matrixOfStrings && matrixOfStrings->IsAStringMatrix()) {
            likelihoodFunctionSpec = new _List;
            matrixOfStrings->FillInList(*likelihoodFunctionSpec);
            if (likelihoodFunctionSpec->lLength == 0) {
                DeleteObject (likelihoodFunctionSpec);
                likelihoodFunctionSpec = nil;
            }
        }
        if (likelihoodFunctionSpec == nil) {
            WarnError (_String("Not a valid string matrix object passed to a _LikelihoodFunction constructor: ") & *(_String*)parameters(1));
            return;
        }
    }

    long i       = 0,
         stepper = explicitFreqs?3:2;

    for (; i<=(long)likelihoodFunctionSpec->lLength-stepper; i+=stepper) {
        _String     *dataset = (_String*)(*likelihoodFunctionSpec)(i),
                     *tree   = (_String*)(*likelihoodFunctionSpec)(i+1),
                      *freq    = explicitFreqs?(_String*)(*likelihoodFunctionSpec)(i+2):nil;

        if(FindDataSetFilterName(AppendContainerName(*dataset,chain.nameSpacePrefix))!=-1) {
            _TheTree*   thisTree = (_TheTree*)FetchObjectFromVariableByType(&AppendContainerName(*tree,chain.nameSpacePrefix),TREE);
            if (thisTree) {
                _CalcNode*  thisNode = thisTree->DepthWiseTraversal(true);
                if (!freq) { // no explicit frequency parameter; grab one from the tree
                    long        theFreqID       = -1,
                                theModelID     = -1,
                                finalFreqID        = -1;
                    bool        done = false;

                    while (1) {
                        if ((theModelID     = thisNode->GetModelIndex()) == HY_NO_MODEL) { // this node has no model
                            done = false;
                            break;
                        }
                        theFreqID   = modelFrequenciesIndices.lData[theModelID];
                        thisNode    = thisTree->DepthWiseTraversal();
                        while(thisNode) {
                            theModelID      = thisNode->GetModelIndex();
                            if (theModelID == HY_NO_MODEL) { // no model
                                done = false;
                                break;
                            }
                            if (modelFrequenciesIndices.lData[theModelID]!=theFreqID) {
                                done = true;
                                break;
                            }
                            if (thisTree->IsCurrentNodeTheRoot()) {
                                break;
                            }
                            thisNode = thisTree->DepthWiseTraversal();
                        }
                        if (theFreqID<0) {
                            finalFreqID = -theFreqID-1;
                        } else {
                            finalFreqID = theFreqID;
                        }
                        break;
                    }

                    if (finalFreqID>=0) {
                        _String freqID = chain.TrimNameSpaceFromID(*LocateVar(finalFreqID)->GetName());
                        passThisToLFConstructor &&  dataset;
                        passThisToLFConstructor &&  tree;
                        passThisToLFConstructor &&  freqID;
                        continue;
                    } else {
                        if (!done) {
                            errMsg = (((_String)("LF: Not a well-defined tree/model combination: ")&*tree));
                        } else {
                            errMsg = (((_String)("LF: All models in the tree: ")&*tree&_String(" must share the same frequencies vector")));
                        }
                    }
                } else {
                    if (FetchObjectFromVariableByType(&AppendContainerName(*freq,chain.nameSpacePrefix),MATRIX)) {
                        passThisToLFConstructor &&  dataset;
                        passThisToLFConstructor &&  tree;
                        passThisToLFConstructor &&  freq;
                        continue;
                    }
                    errMsg = (((_String)("LF: Not a valid frequency matrix ID: ")& *freq));
                }
            } else {
                errMsg = (((_String)("LF: Not a valid tree ID: ")& *tree));
            }

        } else {
            errMsg = (((_String)("LF: Not a valid dataset filter: ")& *dataset));
        }

        if (errMsg.sLength) {
            DeleteObject (likelihoodFunctionSpec);
            WarnError    (errMsg);
            return;
        }
    }

    if (i==likelihoodFunctionSpec->lLength-1) { // computing template
        passThisToLFConstructor && *((_String*)(*likelihoodFunctionSpec)(i));
    }


    DeleteObject (likelihoodFunctionSpec);


    _String lfID                  = chain.AddNameSpaceToID(*(_String*)parameters(0)); // the ID of the likelihood function
    long    likeFuncObjectID      = FindLikeFuncName (lfID);
    if (likeFuncObjectID==-1)
        // not an existing LF ID
    {
        _LikelihoodFunction * lkf = new _LikelihoodFunction ();
        if (! lkf->Construct(passThisToLFConstructor,chain.nameSpacePrefix))
            // constructor failed
        {
            DeleteObject (lkf);
        } else {
            likeFuncObjectID = likeFuncNamesList.Find(&empty);
            // see if there are any vacated spots in the list

            if (likeFuncObjectID < 0) {
                likeFuncList << lkf;
                likeFuncNamesList&&(&lfID);
                DeleteObject (lkf);
            } else {
                likeFuncNamesList.Replace(likeFuncObjectID,&lfID,true);
                likeFuncList.lData[likeFuncObjectID] = (long)lkf;
            }
        }
    } else {
        _LikelihoodFunction* lkf = (_LikelihoodFunction*)likeFuncList(likeFuncObjectID);
        if (!lkf->Construct(passThisToLFConstructor,chain.nameSpacePrefix)) {
            KillLFRecord (likeFuncObjectID,false);
        }
    }

}



//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase12 (_ExecutionList& chain)
{
    chain.currentCommand++;
    SetStatusLine ("Simulating Data");

    _String  likefID        = chain.AddNameSpaceToID(*(_String*)parameters(1)),
             tempString     = ProcessStringArgument (&likefID),
             errMsg;

    if (tempString.sLength) {
        likefID = tempString;
    }

    long f  = FindLikeFuncName (likefID),
         s2 = FindSCFGName     (likefID);

    if (f==-1 && s2==-1) {
        WarnError (_String("Likelihood Function (or SCFG)")&likefID& " has not been initialized" );
        return ;
    }

    if (f>=0) {
        _DataSet  * ds = new _DataSet;
        checkPointer (ds);

        _List     theExclusions;

        if (parameters.lLength>2) // there is a list of exclusions there
            // ';'-sep for different partititons
            // ','-sep for states in a given partition
        {
            // SLKP mod 20070622 to allow string expressions as well
            _String theExc (ProcessLiteralArgument((_String*)parameters(2),chain.nameSpacePrefix));
            if (theExc.sLength) {
                long f = theExc.Find(';'),
                     g = 0;

                while(1) {
                    _String subExc (theExc,g,(f==-1)?(-1):(f-1));
                    long    h = subExc.Find(','),
                            l = 0;
                    _List   myExc;

                    while(1) {
                        _String excludeMe (subExc,l,(h==-1)?(-1):(h-1));
                        myExc && & excludeMe;
                        if (h==-1) {
                            break;
                        }
                        l = h+1;
                        h = subExc.Find(',',h+1,-1);
                    }
                    theExclusions&& &myExc;
                    if (f==-1) {
                        break;
                    }
                    g = f+1;
                    f = theExc.Find(';',f+1,-1);
                }
            }

        }

        _Matrix  *   catValues  = nil,
                     *   catNames   = nil;

        _Variable*   catValVar  = nil,
                     *   catNameVar = nil;

        if (parameters.lLength>3)
            // a matrix to store simulated category values
        {
            _String  matrixName (chain.AddNameSpaceToID(*(_String*)parameters(3)));
            if (!(catValVar = CheckReceptacle(&matrixName,blSimulateDataSet,true))) {
                return;
            } else {
                checkPointer (catValues = new _Matrix (1,1,false,true));
            }
        }

        if (parameters.lLength>4)
            // a matrix to store simulated category values
        {
            _String  matrixName (chain.AddNameSpaceToID(*(_String*)parameters(4)));
            if (!(catNameVar = CheckReceptacle(&matrixName,blSimulateDataSet,true))) {
                return;
            } else {
                checkPointer (catNames = new _Matrix (1,1,false,true));
            }
        }

        _String * resultingDSName = new _String (chain.AddNameSpaceToID(*(_String*)parameters(0)));

        if (!resultingDSName->IsValidIdentifier(true)) {
            errMsg = *resultingDSName & " is not a valid receptacle identifier in call to " & blSimulateDataSet;
            DeleteObject (resultingDSName);
            WarnError (errMsg);
            return;
        }

        ((_LikelihoodFunction*)likeFuncList(f))->Simulate(*ds,theExclusions,catValues,catNames);

        if (catValues) {
            catValVar->SetValue(catValues,false);
        }
        if (catNames) {
            catNameVar->SetValue(catNames,false);
        }

        StoreADataSet (ds, resultingDSName);
        DeleteObject  (resultingDSName);
    } else {
        _String newCorpus = chain.AddNameSpaceToID(*(_String*)parameters(0));
        CheckReceptacleAndStore (&newCorpus," SimulateDataSet (SCFG)", true, new _FString(((Scfg*)scfgList (s2))->SpawnRandomString()), false);
    }
}

//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase38 (_ExecutionList& chain, bool sample)
{
    chain.currentCommand++;
    SetStatusLine ("Reconstructing Ancestors");

    _String *likef          = (_String*)parameters(1),
             tempString      = ProcessStringArgument (likef),
             errMsg;

    if (tempString.sLength) {
        likef = &tempString;
    }

    _String name2lookup = AppendContainerName(*likef,chain.nameSpacePrefix);
    long    objectID    = FindLikeFuncName (name2lookup);
    if (objectID >= 0) {
        _DataSet     * ds               = (_DataSet*) checkPointer(new _DataSet);
        _String      * dsName           = new _String (AppendContainerName(*(_String*)parameters(0),chain.nameSpacePrefix));
        _LikelihoodFunction *lf         = ((_LikelihoodFunction*)likeFuncList(objectID));

        _Matrix * partitionList         = nil;
        if (parameters.lLength>2) {
            _String  secondArg = *(_String*)parameters(2);
            partitionList = (_Matrix*)ProcessAnArgumentByType (&secondArg, chain.nameSpacePrefix, MATRIX);
        }
        _SimpleList                     partsToDo;
        if (lf->ProcessPartitionList(partsToDo, partitionList, " ancestral reconstruction")) {
            lf->ReconstructAncestors(*ds, partsToDo, *dsName,  sample, simpleParameters.Find(-1) >= 0, simpleParameters.Find(-2) >= 0 );
        }
        StoreADataSet  (ds, dsName);
        DeleteObject   (dsName);
    } else {
        objectID    =   FindSCFGName       (name2lookup);
        if (objectID>=0)
            /* reconstruct best parse tree for corpus using SCFG */
        {
            CheckReceptacleAndStore (&AppendContainerName(*(_String*)parameters(0),chain.nameSpacePrefix)," ReconstructAncestors (SCFG)", true, new _FString( ((Scfg*)scfgList (objectID))->BestParseTree() ), false);
        } else {
            errMsg =  (((_String)("Likelihood Function/SCFG")&*likef&_String(" has not been initialized")));
            WarnError (errMsg);
            return ;
        }
    }
}

//____________________________________________________________________________________


void      _ElementaryCommand::ExecuteCase39 (_ExecutionList& chain)
{
    chain.currentCommand++;

    _String *commands,
            theCommand,
            *namespc = nil;

    if (code == 39) {
        commands = ProcessCommandArgument((_String*)parameters(0));
    } else {
        _String filePath = GetStringFromFormula((_String*)parameters(0),chain.nameSpacePrefix),
                originalPath = filePath;


        FILE * commandSource = nil;
        if (code == 66) {
            bool hasExtension    = filePath.FindBackwards (".",0,-1) > 0;

            for (unsigned long p = 0; !commandSource && p < standardLibraryPaths.lLength; p++) {
                for (unsigned long e = 0; !commandSource && e < standardLibraryExtensions.lLength; e++) {
                    _String tryPath = *((_String*)standardLibraryPaths(p)) & filePath & *((_String*)standardLibraryExtensions(e));

                    // printf ("%s\n", tryPath.sData);

                    if (loadedLibraryPaths.Find(&tryPath) >= 0 && parameters.lLength == 2) {
                        ReportWarning (_String("Already loaded '") & originalPath & "' from " & tryPath);
                        return;
                    }
                    if ((commandSource = doFileOpen (tryPath.getStr(), "rb"))) {
                        filePath = tryPath;
                        break;
                    }
                    if (hasExtension) {
                        break;
                    }
                }
            }

        }

        if (commandSource == nil) {
            filePath.ProcessFileName (false,false,(Ptr)chain.nameSpacePrefix);
            if ((commandSource = doFileOpen (filePath.getStr(), "rb")) == nil) {
                WarnError (_String("Could not read command file in ExecuteAFile.\nOriginal path: '") &
                                    originalPath & "'.\nExpanded path: '" & filePath & "'");
                return;
            }
        }

        if (code == 66 && commandSource) {
            ReportWarning (_String("Loaded '") & originalPath & "' from " & filePath);
            loadedLibraryPaths.Insert (filePath.makeDynamic(),0,false,true);
        }

        commands = new _String (commandSource);
        if (fclose       (commandSource) ) { // failed to fclose
            DeleteObject (commands);
            WarnError (_String("Internal error: failed in a call to fclose ") & filePath);
        }
        PushFilePath (filePath);
    }

    if (!commands) {
        return;
    }

    if (code == 39) {
        theCommand = ProcessLiteralArgument (commands,chain.nameSpacePrefix);
    } else {
        theCommand = commands;
    }

    if (theCommand.sLength == 0) {
        WarnError (_String("Invalid string argument '") & *commands & "' in call to ExecuteCommands/ExecuteAFile.");
        return ;
    }

    if (code == 39 && ((_String*)parameters(1))->sLength) {
        pathNames << (_String*)parameters(1);
    }

    _AVLListXL * inArg    = nil;
    _List      * inArgAux = nil;

    if (parameters.lLength >= 3)
        // stdin redirect (and/or name space prefix)
    {
        _PMathObj inAVL = ProcessDictionaryArgument ((_String*)parameters(2),chain.nameSpacePrefix);

        if (!inAVL) {
            if (parameters.lLength == 3) {
                WarnError (_String("Not a valid associative array index passed as input redirect argument to ExecuteCommands/ExecuteAFile: )") & *(_String*)parameters(2));
                return;
            }
        } else {
            _AssociativeList * stdinRedirect = (_AssociativeList*)inAVL;

            checkPointer (inArgAux = new _List);
            checkPointer (inArg    = new _AVLListXL (inArgAux));

            _List        *stdKeys = stdinRedirect->GetKeys();

            for (long kid = 0; kid < stdKeys->lLength; kid++) {
                _String  * aKey         = (_String*)(*stdKeys) (kid);
                if (aKey) {
                    _FString * aString      = (_FString*)stdinRedirect->GetByKey (*aKey, STRING);
                    if (!aString) {
                        WarnError    (_String("All entries in the associative array used as input redirect argument to ExecuteCommands/ExecuteAFile must be strings. The following key was not: ") & *aKey);
                        DeleteObject (inAVL);
                        return;
                    }
                    inArg -> Insert (aKey->makeDynamic(),(long)new _String (*aString->theString),false);
                }
            }
        }

        DeleteObject (inAVL);

        if (parameters.lLength > 3) {
            _String nameSpaceID = ProcessLiteralArgument((_String*)parameters(3),chain.nameSpacePrefix);
            if (!nameSpaceID.IsValidIdentifier(true)) {
                WarnError (_String("Invalid namespace ID in call to ExecuteCommands/ExecuteAFile: ") & *(_String*)parameters(3));
                return;
            }
            namespc = new _String (nameSpaceID);
        }
    }

    if (parameters.lLength <4 && chain.nameSpacePrefix) {
        namespc = new _String (*chain.nameSpacePrefix->GetName());
    }

    if (theCommand.beginswith ("#NEXUS")) {
        ReadDataSetFile (nil,1,&theCommand,nil,namespc);
    } else {
        bool result = false;
        _ExecutionList exc (theCommand,namespc, false, &result);
        
        if (!result) {
            chain.ReportAnExecutionError("Encountered an error while parsing HBL", false, true);
        } else {

            exc.stdinRedirectAux = inArgAux?inArgAux:chain.stdinRedirectAux;
            exc.stdinRedirect    = inArg?inArg:chain.stdinRedirect;

            if (simpleParameters.lLength && exc.TryToMakeSimple()) {
                ReportWarning ("Successfully compiled an execution list.");
                exc.ExecuteSimple ();
            } else {
                exc.Execute();
            }

            exc.stdinRedirectAux = nil;
            exc.stdinRedirect    = nil;
            if (exc.result) {
                DeleteObject (chain.result);
                chain.result = exc.result;
                exc.result = nil;
            }
        }
    }

    if (inArg) {
        DeleteObject            (inArg);
        DeleteObject            (inArgAux);
    }

    DeleteObject (namespc);

    if (code == 62) {
        PopFilePath ();
    } else if (((_String*)parameters(1))->sLength) {
        pathNames.Delete (pathNames.lLength-1);
    }
}

void      _ElementaryCommand::ExecuteCase40 (_ExecutionList& chain)
{
    chain.currentCommand++;
    _String errMsg;

#if !defined __UNIX__ && !defined __HEADLESS__

    bool    done = false;

    _String *windowType     = ((_String*)parameters(0)),
             *positionString = nil;

    _HYWindow*
    resultWindow = nil;

    long    positionSpec[] = {-1,-1,-1,-1};

    if (parameters.lLength==3) {
        positionString = ((_String*)parameters(2));
        _List * positionList = positionString->Tokenize(";");

        if ((positionList->lLength>=2)&&(positionList->lLength<=4)) {
            _HYRect      screenDim = GetScreenDimensions();
            setParameter (screenWidthVar , screenDim.right );
            setParameter (screenHeightVar, screenDim.bottom);

            for (long i = 0; i<positionList->lLength; i++) {
                _String* thisSpec = (_String*)(*positionList)(i);
                if (thisSpec->sLength) {
                    _Formula  fla    (*thisSpec,nil,false);
                    _PMathObj flav = fla.Compute();
                    if (flav&&(flav->ObjectClass()==NUMBER)) {
                        positionSpec[i] = flav->Value();
                    }
                }
            }
        }
        DeleteObject (positionList);
    }

    if (windowType->Equal (&windowTypeClose)) {
        _String     wName   = GetStringFromFormula ((_String*)parameters(1),chain.nameSpacePrefix);
        long        wID     = FindWindowByName (wName);
        if (wID >= 0) {
            postWindowCloseEvent (((_HYWindow*)windowObjectRefs (wID))->GetID());
        }
        done = true;
    } else if (windowType->Equal (&windowTypeTree)||windowType->Equal (&windowTypeTable)||windowType->Equal (&windowTypeDistribTable)||windowType->Equal (&windowTypeDatabase)) {
        _String  * windowOptions = ((_String*)parameters(1));
        _Matrix  * theOptions = (_Matrix*)FetchObjectFromVariableByType(&AppendContainerName(*windowOptions,chain.nameSpacePrefix),MATRIX);
        if (!theOptions) {
            if ((windowOptions->sLength>3)&&(windowOptions->sData[0]=='{')
                    &&(windowOptions->sData[windowOptions->sLength-1]=='}')) {
                theOptions = new _Matrix (*windowOptions);
            } else {
                WarnError (*windowOptions & " is not a valid inline matrix specification/variable reference in call to OpenWindow.");
                return;
            }
        } else {
            theOptions->nInstances++;
        }

        if (theOptions->MatrixType()!=2||theOptions->GetVDim()>1||theOptions->GetHDim()<1) {
            DeleteObject (theOptions);
            WarnError (*windowOptions & " is not a valid options matrix in call to OpenWindow.");
            return;
        }

        if (windowType->Equal (&windowTypeTree)) {
            // row 1: a tree var reference
            _PMathObj arg = theOptions->GetFormula(0,0)->Compute();
            if (!(arg&&(arg->ObjectClass () == STRING))) {
                DeleteObject (theOptions);
                errMsg = "The first entry in the specification matrix must be a tree variable ID string in call to OpenWindow.";
                WarnError (errMsg);
                return;
            }

            windowOptions = ((_FString*)arg)->theString;

            long f = LocateVarByName (AppendContainerName(*windowOptions,chain.nameSpacePrefix));
            if (f>=0) {
                _Variable* treeVar = FetchVar (f);
                if (treeVar->ObjectClass()==TREE) {
                    _String windowName = _String ("Tree ")&*treeVar->GetName();
                    long fw = FindWindowByName (windowName);
                    _HYTreePanel * myTP;
                    if (fw>=0) {
                        myTP = (_HYTreePanel*)windowObjectRefs (fw);
                        myTP->SetVariableReference   (*treeVar->GetName());
                        myTP->BuildTree (true);
                    } else {
                        myTP = new _HYTreePanel (*treeVar->GetName(), *treeVar->GetName());
                        checkPointer (myTP);
                    }
                    resultWindow = myTP;
                    for (long k=1; k<theOptions->GetHDim(); k++) {
                        _PMathObj arg = theOptions->GetFormula(k,0)->Compute();
                        if (!(arg&&(arg->ObjectClass () == STRING))) {
                            DeleteObject (theOptions);
                            errMsg = "All entries in the specification matrix must be strings in call to OpenWindow.";
                            WarnError (errMsg);
                            return;
                        }
                        windowOptions = ((_FString*)arg)->theString;
                        if (windowOptions->sLength) {
                            switch (k) {
                            case 1: { // view options
                                unsigned int treeFlags = windowOptions->toNum();
                                if ((treeFlags&HY_TREEPANEL_SCALE_TO_WINDOW)!=(myTP->GetFlags()&HY_TREEPANEL_SCALE_TO_WINDOW)) {
                                    myTP->ToggleScaleOption ();
                                }
                                myTP->SetFlags (treeFlags);
                                break;
                            }
                            case 2: { // window position spec
                                PositionWindow (myTP,windowOptions);
                                break;
                            }
                            case 3: { // scale variable
                                myTP->SetScaleVariable (*windowOptions);
                                break;
                            }
                            case 4: { // select branches
                                _List * nodes2Select = windowOptions->Tokenize(",");
                                if (nodes2Select->lLength) {
                                    myTP->SelectRangeAndScroll (*nodes2Select);
                                }
                                DeleteObject (nodes2Select);
                                break;
                            }
                            }
                        }
                    }
                } else {
                    f = -1;
                }
            }
            if (f<0) {
                errMsg = *windowOptions & " is not the ID of an existing tree in call to OpenWindow.";
                WarnError (errMsg);
            }

            DeleteObject (theOptions);
            done = true;
        } else if (windowType->Equal (&windowTypeDatabase)) {
            _PMathObj arg = theOptions->GetFormula(0,0)->Compute();
            if (!(arg&&(arg->ObjectClass () == STRING))) {
                DeleteObject (theOptions);
                errMsg = "The first entry in the specification matrix must be a file name string in call to OpenWindow.";
                WarnError (errMsg);
                return;
            }

            windowOptions = ((_FString*)arg)->theString;

            _String windowName = _String ("SQLite DB: ")&*windowOptions;
            long fw = FindWindowByName (windowName);
            _HYDBWindow * myDW;
            if (fw>=0) {
                myDW = (_HYDBWindow*)windowObjectRefs (fw);
                myDW->SetDB  (windowOptions);
            } else {
                myDW = new _HYDBWindow (windowName, windowOptions);
                checkPointer (myDW);
            }
            resultWindow = myDW;

            DeleteObject (theOptions);
            done = true;
        } else {
            bool  isD = windowType->Equal (&windowTypeDistribTable);

            if (theOptions->GetHDim()<3) {
                DeleteObject (theOptions);

                errMsg = "The matrix specification for a chart window must have at least 3 entries in call to OpenWindow.";
                WarnError (errMsg);
                return;

            }
            _PMathObj arg  = theOptions->GetFormula(1,0)->Compute(),
                      arg2 = theOptions->GetFormula(2,0)->Compute(),
                      arg3 = theOptions->GetFormula(0,0)->Compute();

            if (!(arg&&(arg->ObjectClass () == STRING)&&arg2&&(arg2->ObjectClass () == STRING)&&arg3&&(arg3->ObjectClass () == STRING))) {
                DeleteObject (theOptions);

                errMsg = "The first two entries in the specification matrix must be strings with matrix ID (e.g. \"matrixID\") in call to OpenWindow.";
                WarnError (errMsg);
                return;
            }

            windowOptions = ((_FString*)arg)->theString;

            _String      *options2 = ((_FString*)arg2)->theString;
            long f = LocateVarByName            (AppendContainerName(*windowOptions,chain.nameSpacePrefix)),
                 f2 = LocateVarByName   (AppendContainerName(*options2,chain.nameSpacePrefix));

            if ((f>=0)&&(f2>=0)) {
                _Variable* labelMatrix = FetchVar (f),
                           * dataMatrix  = FetchVar (f2);

                if ((labelMatrix->ObjectClass()==MATRIX)&&(dataMatrix->ObjectClass()==MATRIX)) {
                    _Matrix*lMatrix = (_Matrix*)labelMatrix->GetValue(),
                            *dMatrix = (_Matrix*)dataMatrix->GetValue();

                    _String windowName = *((_FString*)arg3)->theString;
                    long fw = FindWindowByName (windowName);
                    _List   columnHeaders;

                    for (f2=0; f2<lMatrix->GetVDim(); f2++) {
                        _Formula  *fff= lMatrix->GetFormula(0,f2);
                        if (fff) {
                            _PMathObj arg = fff->Compute();
                            if (arg->ObjectClass()==STRING) {
                                columnHeaders && ((_FString*)arg)->theString;
                                continue;
                            }
                        }
                        columnHeaders && & empty;

                    }

                    if ((dMatrix->GetVDim()!=columnHeaders.lLength)&&(dMatrix->GetVDim()!=columnHeaders.lLength-1)) {
                        errMsg = "The number of columns in the data matrix must match the dimension of the header matrix in call to OpenWindow (CHART_WINDOW).";
                        WarnError (errMsg);
                        return;
                    }

                    if ((dMatrix->GetHDim()==0)||(lMatrix->GetVDim()==0)) {
                        errMsg = "Empty data matrix (or label matrix) in call to OpenWindow (CHART_WINDOW).";
                        WarnError (errMsg);
                        return;
                    }

                    _String      mL[14];

                    for (long k=3; k<theOptions->GetHDim(); k++) {
                        _PMathObj arg = theOptions->GetFormula(k,0)->Compute();
                        if (!(arg&&(arg->ObjectClass () == STRING))) {
                            DeleteObject (theOptions);
                            errMsg = "All entries in the specification matrix must be strings in call to OpenWindow.";
                            WarnError (errMsg);
                            return;
                        }
                        windowOptions = ((_FString*)arg)->theString;
                        mL [k-3] = *windowOptions;

                        if (k==16) {
                            break;
                        }
                    }

                    _HYChartWindow * myTP;

                    if (isD) {
                        _List   *varInfo = mL[13].Tokenize (";"),
                                 cInfo,
                                 dInfo;


                        bool      firstDerived = false;
                        _String   errMsg;

                        for (long k=0; k < varInfo->lLength; k++) {
                            _List *item = ((_String*)(*varInfo)(k))->Tokenize (":");
                            if (item->lLength == 1) {
                                if (cInfo.lLength) {
                                    dInfo << (*item)(0);
                                } else {
                                    errMsg = "Derived variables appear before the atom variables";
                                }
                                firstDerived = true;
                            } else {
                                if (firstDerived) {
                                    errMsg = "Some atom variables appear after derived variables";
                                } else if ((item->lLength>=3) && (item->lLength % 2 == 1)) {
                                    _String * vName = (_String*)(*item)(0);
                                    if (vName->IsValidIdentifier()) {
                                        long      catCount = (item->lLength - 1)/2;

                                        _Matrix   wts (1,catCount,false,true),
                                                  prb (1,catCount,false,true);

                                        for (long k2 = 0, k3 = catCount+1; k2 < catCount; k2++, k3++) {
                                            wts.theData[k2] = ((_String*)(*item)(k2+1))->toNum();
                                            prb.theData[k2] = ((_String*)(*item)(k3))->toNum();
                                        }

                                        _List anItem;
                                        anItem << vName;
                                        anItem && & wts;
                                        anItem && & prb;
                                        cInfo && & anItem;
                                    } else {
                                        errMsg = *vName & " is an invalid atom variable identifier";
                                    }
                                } else {
                                    errMsg = "Invalid item count in atom variable specification";
                                }
                            }

                            if (errMsg.sLength) {
                                DeleteObject (item);
                                DeleteObject (varInfo);
                                errMsg = errMsg  & " in call to OpenWindow";
                                ProblemReport (errMsg);
                                return;
                            }
                        }

                        DeleteObject (varInfo);

                        if (fw>=0) {
                            myTP = (_HYChartWindow*)windowObjectRefs (fw);
                            myTP->SetTable   (columnHeaders,*dMatrix);
                            ((_HYDistributionChartWindow*)myTP)->SetAtoms (*dMatrix,cInfo);
                        } else {
                            myTP = new _HYDistributionChartWindow (windowName,columnHeaders, *dMatrix, cInfo, nil);
                            checkPointer (myTP);
                        }

                        for (long dc = 0; dc < dInfo.lLength; dc++) {
                            ((_HYDistributionChartWindow*)myTP)->AddVariable ((_String*)dInfo(dc));
                        }
                    } else {
                        if (fw>=0) {
                            myTP = (_HYChartWindow*)windowObjectRefs (fw);
                            myTP->SetTable   (columnHeaders,*dMatrix);
                        } else {
                            myTP = new _HYChartWindow (windowName,columnHeaders,*dMatrix,nil);
                            checkPointer (myTP);
                        }
                    }

                    resultWindow = myTP;

                    myTP->SetChartType (mL[0],mL[1],mL[2],false);
                    myTP->ToggleSuspend (true);

                    long  kk;

                    _Parameter  uMin = 0.0,
                                uMax = uMin;

                    _List * ubounds = mL[12].Tokenize(",");
                    if (ubounds->lLength == 3) {
                        mL[12] = *(_String*)(*ubounds)(0);
                        uMin = ((_String*)(*ubounds)(1))->toNum();
                        uMax = ((_String*)(*ubounds)(2))->toNum();
                    }
                    DeleteObject (ubounds);

                    if ((kk = mL[8].Find(';'))>0) {
                        myTP->SetLabels    (mL[3],mL[4],mL[5],mL[6].toNum(),mL[7],mL[8].Cut (0,kk-1).toNum()-1,mL[8].Cut (kk+1,-1).toNum()-1,mL[12].toNum(),uMin,uMax);
                    } else {
                        myTP->SetLabels    (mL[3],mL[4],mL[5],mL[6].toNum(),mL[7],-1,-1,mL[12].toNum(),uMin,uMax);
                    }

                    _List * optList = mL[9].Tokenize (";");

                    if (optList->lLength == 3) {
                        myTP->SetProjection (((_String*)(*optList)(0))->toNum(),((_String*)(*optList)(1))->toNum(),((_String*)(*optList)(2))->toNum());
                    }

                    DeleteObject (optList);

                    optList = mL[10].Tokenize (";");
                    myTP->SetFonts (optList);
                    DeleteObject (optList);

                    optList = mL[11].Tokenize (";");
                    myTP->SetColors (optList);
                    DeleteObject (optList);

                    myTP->ToggleSuspend (false);
                    myTP->DrawChart     ();
                } else {
                    f = -1;
                }
            }

            if (f<0 || f2<0) {
                errMsg = *windowOptions & " and " & *options2 & " must both refer to existing matrices in call to OpenWindow.";
                WarnError (errMsg);
            }

            DeleteObject (theOptions);
            done = true;
        }
    }

    if (resultWindow) {
        if ((positionSpec[0]>-1)&&(positionSpec[1]>-1)) {
            resultWindow->SetWindowRectangle (0,0,positionSpec[1],positionSpec[0]);
        }
        if ((positionSpec[2]>-1)&&(positionSpec[3]>-1)) {
            resultWindow->SetPosition (positionSpec[2],positionSpec[3]);
        }

        resultWindow->BringToFront ();
        handleGUI(true);
    }

    if (!done) {
        errMsg =  *windowType& " is not a valid window type in call to OpenWindow.";
        WarnError (errMsg);
    }
#endif
}

void      _ElementaryCommand::ExecuteCase41 (_ExecutionList& chain)
{
    chain.currentCommand++;

#if !defined    __UNIX__ && ! defined __HEADLESS__
    _String     lfID,
                dsWindow,
                errMsg;

    long        treeID = -1;

    _SimpleList subsetSpec;

    lfID = ProcessLiteralArgument ((_String*)parameters(0),chain.nameSpacePrefix);

    if (!lfID.IsValidIdentifier()) {
        errMsg = lfID & " is not a valid likelihood function identifier.";
        WarnError (errMsg);
        return;
    }

    dsWindow = ProcessLiteralArgument ((_String*)parameters(1),chain.nameSpacePrefix);

    treeID =  LocateVarByName(dsWindow);

    if (treeID<0 || FetchVar(treeID)->ObjectClass()!=TREE) {
        errMsg = dsWindow & " is not the name of an existing tree.";
        WarnError (errMsg);
        return;
    } else {
        treeID = variableNames.GetXtra (treeID);
    }

    dsWindow = ProcessLiteralArgument ((_String*)parameters(2),chain.nameSpacePrefix);
    long     k = FindWindowByName (dsWindow);

    if ((k<0)||(((_HYWindow*)windowObjectRefs (k))->WindowKind () != HY_WINDOW_KIND_DATAPANEL)) {
        errMsg = dsWindow & " is not the name of an open data panel window.";
        WarnError (errMsg);
        return;
    }

    _HYDataPanel * theDP  = (_HYDataPanel*)windowObjectRefs (k);
    _Matrix * theMatrix = (_Matrix*)FetchObjectFromVariableByType(&AppendContainerName(*(_String*)parameters(3),chain.nameSpacePrefix),MATRIX);

    if (theMatrix) {
        for (long i1 = 0; i1< theMatrix->GetHDim(); i1++)
            for (long i2 = 0; i2 < theMatrix->GetVDim (); i2++) {
                subsetSpec << (*theMatrix)(i1,i2);
            }

        theDP->BuildLikelihoodFunction (&lfID,& subsetSpec, treeID);
        return;
    }

    errMsg = (*(_String*)parameters(3)) & " is not the name of an existing matrix.";
    WarnError (errMsg);

#endif
}


//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase25 (_ExecutionList& chain, bool issscanf)
{
    chain.currentCommand++;
    // first of all obtain the string to be parsed
    // either read the file into a string or get a string from standard input
    _String     currentParameter = *(_String*)parameters(0),
                *data = nil;

    long        p,
                p2 = 0,
                r,
                q,
                v,
                t,
                shifter = simpleParameters.lData[0] < 0;

    bool        skipDataDelete = false;
    _Variable*  iseof          = CheckReceptacle (&hasEndBeenReached,empty,false);

    if (currentParameter==_String("stdin")) { //
        if (chain.stdinRedirect) {
            data = chain.FetchFromStdinRedirect ();
        } else {
            if (!CheckEqual(iseof->Compute()->Value(),0) && currentParameter.Equal (&scanfLastFilePath)) {
                WarnError ("Ran out of standard input\n");
                return;
            }
            checkPointer ((Ptr)(data = new _String (StringFromConsole())));
        }
    } else {
        if (issscanf) {
            currentParameter = chain.AddNameSpaceToID(currentParameter);
            _FString * sscanfData = (_FString*)FetchObjectFromVariableByType(&currentParameter,STRING);
            if (!sscanfData) {
                WarnError         (currentParameter& " does not refer to a string variable in call to sscanf");
                return;
            }
            data = sscanfData->theString;
            skipDataDelete = true;

            if (iseof->Compute()->Value() > 0.) {
                scanfLastFilePath = empty;
            }

            if (!currentParameter.Equal (&scanfLastFilePath) || shifter) {
                scanfLastFilePath     = currentParameter;
                p = scanfLastReadPosition = 0;
            } else {
                p = p2 = scanfLastReadPosition;
                if (p>=data->sLength) {
                    iseof->SetValue (new _Constant (1.0), false);
                    return;
                }
            }
        } else {
            FILE*   inputBuffer;
            if (currentParameter.Find('"')==-1) {
                currentParameter = GetStringFromFormula (&currentParameter,chain.nameSpacePrefix);
            }

            currentParameter.ProcessFileName(false,false,(Ptr)chain.nameSpacePrefix);
            if (terminateExecution) {
                return;
            }
            inputBuffer = doFileOpen (currentParameter.getStr(), "rb");
            if (!inputBuffer) {
                WarnError         (currentParameter& " could not be opened for reading by fscanf. Path stack: " & _String((_String*)pathNames.toStr()));
                return;
            }

            if (iseof->Compute()->Value()>0) {
                scanfLastFilePath = empty;
            }

            if (!currentParameter.Equal (&scanfLastFilePath) || shifter) {
                scanfLastFilePath = currentParameter;
                scanfLastReadPosition = 0;
            }

            fseek (inputBuffer,0,SEEK_END);
            p    = ftell (inputBuffer);
            p   -= scanfLastReadPosition;

            if (p<=0) {
                iseof->SetValue (new _Constant (1.0), false);
                fclose(inputBuffer);
                return;
            }

            data = (_String*)checkPointer(new _String ((unsigned long)p));
            rewind (inputBuffer);
            fseek  (inputBuffer, scanfLastReadPosition, SEEK_SET);
            fread  (data->sData, 1, p, inputBuffer);
            fclose (inputBuffer);
        }
    }
    // now that the string has been read in, read in all the terms, ignoring all the characters in between
    if (!skipDataDelete) {
        p = 0;    // will be used to keep track of the position in the string
    }

    q = 0;
    r = shifter;

    while (r<simpleParameters.lLength && p<data->sLength) {
        _String *currentParameter = ProcessCommandArgument((_String*)parameters(r+1-shifter)); // name of the receptacle
        if (!currentParameter) {
            DeleteObject (data);
            return;
        }
        if (!currentParameter->IsValidIdentifier()) {
            WarnError (_String ('\\') & *currentParameter & "\" is not a valid identifier in call to fscanf.");
            DeleteObject (data);
            return;
        }
        _String namespacedParameter (chain.AddNameSpaceToID(*currentParameter));

        v = LocateVarByName (namespacedParameter);
        if (v<0) {
            if (simpleParameters.lData[r]!=2) {
                v = CheckReceptacle(&namespacedParameter,empty,false)->GetAVariable();
            }
        } else {
            if (simpleParameters.lData[r]==2)
                if (FetchVar(v)->ObjectClass()==TREE) {
                    DeleteVariable(*FetchVar(v)->GetName());
                }
        }


        _Variable * theReceptacle = FetchVar(v); //this will return nil for TREE

        if (simpleParameters.lData[r]==0) { // number
            q = p;
            while (!(((('0'<=data->sData[q])&&(data->sData[q]<='9'))||(data->sData[q]=='-')||(data->sData[q]=='+')||(data->sData[q]=='.')||(data->sData[q]=='e')||(data->sData[q]=='E')))
                    &&(q<data->sLength)) {
                q++;
            }
            p = q;
            while (((('0'<=data->sData[q])&&(data->sData[q]<='9'))||(data->sData[q]=='-')||(data->sData[q]=='+')||(data->sData[q]=='.')||(data->sData[q]=='e')||(data->sData[q]=='E'))
                    &&(q<data->sLength)) {
                q++;
            }
            theReceptacle->SetValue (new _Constant (data->Cut (p, q-1).toNum()), false);
            while ((q<data->sLength-1)&&(isspace (data->sData[q+1]))) {
                q++;
            }
        } else {
            if (simpleParameters.lData[r]==3) { // string
                q=0;
                bool  startFound=false;
                while (q+p<data->sLength) {
                    char c = data->sData[q+p];
                    if (!startFound) {
                        if (!isspace(c)) {
                            p+=q;
                            startFound = true;
                            q=0;
                        }
                    } else if (c=='\n' || c=='\r' || c=='\t') {
                        break;
                    }
                    q++;
                }
                if (startFound) {
                    theReceptacle->SetValue (new _FString (new _String(*data,p,q+p-1)),false);
                } else {
                    theReceptacle->SetValue (new _FString, false);
                }

                p+=q;
                r++;
                continue;
            } else if (simpleParameters.lData[r]==5) { // raw
                theReceptacle->SetValue (new _FString (new _String (*data,p,-1)), false);
                p = data->sLength;
                r++;
                continue;
            } else {
                if (simpleParameters.lData[r]==6) { // lines
                    _String  inData  (*data,p,-1);

                    _List     lines;

                    long      lastP = 0,
                              loopP = 0;

                    for (loopP = 0; loopP < inData.sLength; loopP ++) {
                        if (inData.sData[loopP] == '\r' || inData.sData[loopP] == '\n') {
                            if (lastP<loopP) {
                                lines.AppendNewInstance (new _String (inData,lastP, loopP-1));
                            } else {
                                lines && & empty;
                            }

                            lastP = loopP+1;

                            if (lastP < inData.sLength && (inData.sData[lastP] == '\r' || inData.sData[lastP] == '\n') && (inData.sData[lastP] != inData.sData[lastP-1])) {
                                lastP++;
                            }

                            loopP = lastP-1;
                        }
                    }

                    if (lastP < inData.sLength && lastP<loopP) {
                        lines.AppendNewInstance (new _String (inData,lastP, loopP-1));
                    } else if (lines.lLength == 0) {
                        lines && & empty;
                    }

                    theReceptacle->SetValue (new _Matrix (lines), false);
                    p = data->sLength;
                    r++;
                    continue;
                } else {
                    char delimiter1 = (simpleParameters.lData[r]==2)?'(':'{',
                         delimiter2 = delimiter1=='{'?'}':')';
                    q = data->Find (delimiter1,p,-1);
                    if (q==-1) {
                        break;
                    }
                    p = q;
                    t = 0;
                    do {
                        if (data->sData[q]==delimiter1) {
                            t++;
                        } else if (data->sData[q]==delimiter2) {
                            t--;
                        }
                        q++;
                    } while (t&&(q<data->sLength));
                    if (t) {
                        break;
                    }
                    if ((simpleParameters.lData[r]==1)||(simpleParameters.lData[r]==4)) { // matrix
                        _String localData (*data,p,q-1);
                        _Matrix *newMatrixValue = new _Matrix (localData,simpleParameters.lData[r]==4);
                        checkPointer (newMatrixValue);
                        theReceptacle->SetValue (newMatrixValue, false);
                    } else {
                        long  varID = LocateVarByName (namespacedParameter);
                        if (varID>=0)
                            if (FetchVar(varID)->ObjectClass()==TREE) {
                                DeleteVariable(*FetchVar(varID)->GetName());
                            }
                        _String  treeString (*data, p, q-1);
                        _TheTree dummyTree (namespacedParameter, treeString);
                    }
                }
            }
        }
        p = q+1;
        r++;
    }

    if (r<simpleParameters.lLength) {
        ReportWarning ("fscanf could not read all the parameters requested.");
        iseof->SetValue (new _Constant (1.0), false);
    } else {
        iseof->SetValue (new _Constant (0.0), false);
    }

    if (skipDataDelete) {
        scanfLastReadPosition += p-p2;
    } else {
        scanfLastReadPosition += p;
        DeleteObject (data);
    }
}
//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase31 (_ExecutionList& chain)
// 20100312 SLKP: added matrix-expression based model
// definitions
{
    chain.currentCommand++;
    // first check to see if matrix parameters here are valid

    bool     usingLastDefMatrix = false,
             doExpressionBased  = false;

    _Formula *isExpressionBased  = nil;

    _String* parameterName,
             errMsg,
             arg0 = chain.AddNameSpaceToID(*(_String*)parameters(0));

    long     f,
             f2=-1,
             matrixDim,
             f3,
             multFreqs = 1;



    if (parameters.lLength>3) {
        parameterName = (_String*)parameters.lData[3];
        if (parameterName->Equal(&ModelTrainNNFlag)) {
            _String arg1 = chain.AddNameSpaceToID(*(_String*)parameters(1));
            TrainModelNN (&arg0,&arg1);
            return;
        } else

            if (parameterName->Equal(&explicitFormMExp)) {
                doExpressionBased = true;
                multFreqs         = 0;
            } else {
                multFreqs = ProcessNumericArgument (parameterName,chain.nameSpacePrefix);
            }
    }

    _Matrix*  checkMatrix = nil;

    parameterName = (_String*)parameters.lData[1];

    if (parameterName->Equal (&useLastDefinedMatrix)) {
        if (lastMatrixDeclared<0) {
            errMsg = "First Call to Model. USE_LAST_DEFINED_MATRIX is meaningless.";
            acknError (errMsg);
            return;
        }
        f3 = lastMatrixDeclared;
        f  = modelMatrixIndices[f3];
        usingLastDefMatrix = true;
    } else {
        if (doExpressionBased) {
            _String matrixExpression (ProcessLiteralArgument((_String*)parameters.lData[1],chain.nameSpacePrefix)),
                    defErrMsg = _String ("The expression for the explicit matrix exponential passed to Model must be a valid matrix-valued HyPhy formula that is not an assignment.") & ':' & matrixExpression;
            // try to parse the expression, confirm that it is a square  matrix,
            // and that it is a valid transition matrix
            long                varRef = 0;
            isExpressionBased = (_Formula*)checkPointer(new _Formula);
            long parseCode = Parse(isExpressionBased,matrixExpression,varRef,chain.nameSpacePrefix);
            if (parseCode != HY_FORMULA_EXPRESSION || isExpressionBased->ObjectClass()!= MATRIX ) {
                WarnError (defErrMsg );
                return;
            }
            
            //for (unsigned long k = 0; k < isExpressionBased
            
            checkMatrix = (_Matrix*)isExpressionBased->Compute();


        } else {
            parameterName = (_String*)parameters.lData[1];

            _String augName (chain.AddNameSpaceToID(*parameterName));
            f = LocateVarByName (augName);

            if (f<0) {
                WarnError (*parameterName & " has not been defined prior to the call to Model = ...");
                return;
            }

            _Variable* checkVar = usingLastDefMatrix?LocateVar(f):FetchVar (f);
            if (checkVar->ObjectClass()!=MATRIX) {
                WarnError (*parameterName & " must refer to a matrix in the call to Model = ...");
                return;
            }
            checkMatrix = (_Matrix*)checkVar->GetValue();
        }
    }


    matrixDim = checkMatrix->GetHDim();
    if ( matrixDim!=checkMatrix->GetVDim() || matrixDim<2 ) {
        WarnError (*parameterName & " must be a square matrix of dimension>=2 in the call to Model = ...");
        return;
    }

    // so far so good

    parameterName = (_String*)parameters.lData[2]; // this is the frequency matrix (if there is one!)
    _String         freqNameTag (chain.AddNameSpaceToID(*parameterName));

    f2 = LocateVarByName (freqNameTag);
    if (f2<0) {
        WarnError(*parameterName & " has not been defined prior to the call to Model = ...");
        return;
    }
    _Variable * checkVar = FetchVar (f2);
    if (checkVar->ObjectClass()!=MATRIX) {
        WarnError (*parameterName & " must refer to a column/row vector in the call to Model = ...");
        return;
    }
    checkMatrix = (_Matrix*)checkVar->GetValue();
    if (checkMatrix->GetVDim()==1) {
        if (checkMatrix->GetHDim()!=matrixDim) {
            WarnError (*parameterName & " must be a column vector of the same dimension as the model matrix in the call to Model = ...");
            return;
        }
    } else if (checkMatrix->GetHDim()==1) {
        if (checkMatrix->GetVDim()!=matrixDim) {
            WarnError ( *parameterName & " must be a row vector of the same dimension as the model matrix in the call to Model = ...");
            return;
        }
        errMsg = *parameterName & " has been transposed to the default column vector setting ";
        checkMatrix->Transpose();
        ReportWarning (errMsg);
    } else {
        WarnError (*parameterName & " must refer to a column/row vector in the call to Model = ...");
        return;
    }

    if (usingLastDefMatrix) {
        if (modelFrequenciesIndices[f3]<0) {
            f2 = -f2-1;
        }
    } else if (multFreqs == 0) { // optional flag present
        f2 = -f2-1;
    }

    long existingIndex = modelNames.Find(&arg0);

    if (existingIndex == -1) { // name not found
        lastMatrixDeclared = modelNames.Find (&empty);

        if (lastMatrixDeclared>=0) {
            modelNames.Replace (lastMatrixDeclared,&arg0,true);
            modelTypeList.lData[lastMatrixDeclared] = isExpressionBased?matrixDim:0;
            if (isExpressionBased) {
                modelMatrixIndices.lData[lastMatrixDeclared] = (long)isExpressionBased;
            } else {
                modelMatrixIndices.lData[lastMatrixDeclared] = (usingLastDefMatrix?f:variableNames.GetXtra(f));
            }

            if (f2>=0) {
                modelFrequenciesIndices.lData[lastMatrixDeclared] = variableNames.GetXtra(f2);
            } else {
                modelFrequenciesIndices.lData[lastMatrixDeclared] = -variableNames.GetXtra(-f2-1)-1;
            }
        } else {
            modelNames && & arg0;
            modelTypeList << (isExpressionBased?matrixDim:0);
            if (isExpressionBased) {
                modelMatrixIndices << (long)isExpressionBased;
            } else {
                modelMatrixIndices << (usingLastDefMatrix?f:variableNames.GetXtra(f));
            }
            if (f2>=0) {
                modelFrequenciesIndices << variableNames.GetXtra(f2);
            } else {
                modelFrequenciesIndices << -variableNames.GetXtra(-f2-1)-1;
            }
            lastMatrixDeclared = modelNames.lLength-1;
        }
    } else {
        modelNames.Replace(existingIndex,&arg0,true);
        if (modelTypeList.lData[existingIndex]) {
            delete ((_Formula*)modelMatrixIndices[existingIndex]);
        }

        modelTypeList.lData[existingIndex] = isExpressionBased?matrixDim:0;
        if (isExpressionBased) {
            modelMatrixIndices[existingIndex] = (long)isExpressionBased;
        } else {
            modelMatrixIndices[existingIndex] = usingLastDefMatrix?f:variableNames.GetXtra(f);
        }


        if (f2>=0) {
            modelFrequenciesIndices[existingIndex] = variableNames.GetXtra(f2);
        } else {
            modelFrequenciesIndices[existingIndex] = -variableNames.GetXtra(-f2-1)-1;
        }

        lastMatrixDeclared = existingIndex;
    }
}

//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase32 (_ExecutionList& chain)
{
    chain.currentCommand++;
    // first check to see if matrix parameters here are valid
    long        f = LocateVarByName (AppendContainerName(*(_String*)parameters (3),chain.nameSpacePrefix)),
                fixedLength = ProcessNumericArgument((_String*)parameters(2),chain.nameSpacePrefix);


    _String     saveTheArg;
    _SimpleList sel,
                exclusions;

    _Variable* holder;

    if (fixedLength<0) {
        fixedLength = 1;
        saveTheArg = *(_String*)parameters(2) & " should represent a non-negative integer in call to ChoiceList. The value was reset to 1";
        ReportWarning (saveTheArg);
    }

    if (f>=0) {
        holder = FetchVar(f);
        if (holder->ObjectClass()==NUMBER) {
            if ((f = holder->Value())>=0) {
                exclusions<<f;
            }
        } else if (holder->ObjectClass()==MATRIX) {
            _Matrix* theExcl = (_Matrix*)holder->GetValue()->Compute();
            for (long k=theExcl->GetHDim()*theExcl->GetVDim()-1; k>=0; k--) {
                f = (*theExcl)[k];
                if (f>=0) {
                    exclusions<<f;
                }
            }
            exclusions.Sort();
        }
    }

    holder = CheckReceptacle (&AppendContainerName(*(_String*)parameters (0),chain.nameSpacePrefix), "Choice List", true);
    holder->SetBounds (-2.0, holder->GetUpperBound());

    bool    validChoices = simpleParameters.lData[0] == 0;

    if (simpleParameters.lData[0])
        // some data structure present - process accordingly
    {
        saveTheArg = *(_String*)parameters(4);
        // see if there is a "standard argument"
        _List choices;
        if (saveTheArg == _String("LikelihoodFunction")) {
            parameters.Delete(4);
            for (f=0; f<likeFuncList.lLength; f++) {
                if (exclusions.BinaryFind(f)>=0) {
                    continue;
                }

                if (likeFuncList.lData[f]) {
                    _List thisPair;
                    thisPair << likeFuncNamesList(f);
                    _String likeFuncDesc ("Likelihood Function \"");
                    likeFuncDesc = likeFuncDesc&*(_String*)likeFuncNamesList(f)&("\".");

                    thisPair && &likeFuncDesc;
                    choices&& &thisPair;
                }
            }
            validChoices = true;
            parameters&& & choices;
        } else {
            _String nmspName = AppendContainerName(saveTheArg,chain.nameSpacePrefix);
            f = FindDataSetFilterName (nmspName);
            if (f>=0) {
                parameters.Delete(4);
                _DataSetFilter *theFilter = (_DataSetFilter*)dataSetFilterList (f);
                for (f = 0; f<theFilter->NumberSpecies(); f++) {
                    if (exclusions.BinaryFind(f)>=0) {
                        continue;
                    }

                    _List thisPair;
                    thisPair<< theFilter->GetData()->GetNames() (f);
                    _String spNumber ("Taxon ");
                    spNumber = spNumber & (f+1) & '(' & *(_String*)theFilter->GetData()->GetNames() (f) & ')';
                    thisPair && &spNumber;
                    choices&& &thisPair;
                }
                validChoices = true;
                parameters&& & choices;
            } else {
                f = FindDataSetName (nmspName);
                if (f>=0) {
                    parameters.Delete(4);
                    _DataSet *theSet = (_DataSet*)dataSetList (f);
                    for (f = 0; f<theSet->NoOfSpecies(); f++) {
                        if (exclusions.BinaryFind(f)>=0) {
                            continue;
                        }
                        _List thisPair;
                        thisPair<< theSet->GetNames() (f);
                        _String spNumber ("Taxon ");
                        spNumber = spNumber & (f+1) & '(' & *(_String*)theSet->GetNames() (f) & ')';
                        thisPair && &spNumber;
                        choices&& &thisPair;
                    }
                    validChoices = true;
                    parameters&& & choices;
                } else {
                    if (saveTheArg==lastModelParameterList) {
                        f = lastMatrixDeclared;
                    } else {
                        f = modelNames.Find(&nmspName);
                    }

                    if (f>=0) {
                        parameters.Delete(4);
                        _Variable *theSet = LocateVar (modelMatrixIndices.lData[f]);
                        _SimpleList modelParms;
                        _String     ts ("All Parameters");
                        _List       tl;
                        tl && &ts;
                        ts = "All local model parameters are constrained";
                        tl && &ts;
                        choices && &tl;
                        _AVLList modelParmsA (&modelParms);
                        theSet->ScanForVariables(modelParmsA,false);
                        modelParmsA.ReorderList();
                        for (f = 0; f<modelParms.lLength; f++) {
                            if (exclusions.BinaryFind(f)>=0) {
                                continue;
                            }

                            _List thisPair;
                            thisPair<< LocateVar(modelParms.lData[f])->GetName();
                            _String spNumber ("Constrain parameter ");
                            spNumber = spNumber & *LocateVar(modelParms.lData[f])->GetName();
                            thisPair && &spNumber;
                            choices&& &thisPair;
                        }
                        validChoices = true;
                        parameters&& & choices;
                    } else {
                        f = LocateVarByName (nmspName);
                        if (f>=0) {
                            _Variable * theV = FetchVar(f);
                            if (theV->ObjectClass() == MATRIX) {
                                _Matrix * vM = (_Matrix*)theV->GetValue();
                                if (vM->IsAStringMatrix() && (vM->GetVDim () == 2)) {
                                    parameters.Delete(4);
                                    for (f = 0; f<vM->GetHDim(); f++) {
                                        if (exclusions.BinaryFind (f) < 0) {
                                            _Formula   *f1 = vM->GetFormula (f,0),
                                                        *f2 = vM->GetFormula (f,1);

                                            if (f1&&f2) {
                                                _PMathObj p1 = f1->Compute(),
                                                          p2 = f2->Compute();

                                                if (p1&&p2&&(p1->ObjectClass() == STRING)&&(p2->ObjectClass() == STRING)) {
                                                    _List thisPair;
                                                    thisPair << ((_FString*)p1)->theString;
                                                    thisPair << ((_FString*)p2)->theString;
                                                    choices&& &thisPair;
                                                }
                                            }
                                        }
                                    }
                                    validChoices = true;
                                    parameters&& & choices;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (validChoices) {
        long choice = -1;
        _List* theChoices = (_List*)parameters(4);
        if (fixedLength>theChoices->lLength) {
            _String e = "List of selections is too short in ChoiceList";
            acknError (e);
        } else {
            if (chain.stdinRedirect) {
                if (fixedLength == 1) {
                    _String buffer (chain.FetchFromStdinRedirect());
                    for (choice = 0; choice<theChoices->lLength; choice++)
                        if (buffer.Equal ((_String*)(*(_List*)(*theChoices)(choice))(0))) {
                            break;
                        }
                    if (choice == theChoices->lLength) {
                        choice = -1;
                        WarnError (_String("Not a valid option: '") & buffer & "' passed to Choice List '" & ((_String*)parameters(1))->sData & "' using redirected stdin input");
                        return;
                    }
                } else {
                    if (fixedLength>0) {
                        while (sel.lLength<fixedLength) {
                            _String buffer (chain.FetchFromStdinRedirect());
                            for (choice = 0; choice<theChoices->lLength; choice++)
                                if (buffer.Equal ((_String*)(*(_List*)(*theChoices)(choice))(0))) {
                                    break;
                                }
                            if (choice<theChoices->lLength && sel.Find(choice)==-1) {
                                sel<<choice-1;
                            } else {
                                break;
                            }
                        }
                        if (sel.lLength<fixedLength) {
                            WarnError ("Failed to make the required number of choices in ChoiceList using redirected stdin input.");
                            return;
                        }
                    } else
                        while (1) {
                            _String buffer (chain.FetchFromStdinRedirect());
                            if (buffer.sLength) {
                                for (choice = 0; choice<theChoices->lLength; choice++)
                                    if (buffer.Equal ((_String*)(*(_List*)(*theChoices)(choice))(0))) {
                                        break;
                                    }

                                if (choice<theChoices->lLength && sel.Find(choice)==-1) {
                                    sel<<choice;
                                } else {
                                    WarnError (_String("Not a valid (or duplicate) option: '") & buffer & "' passed to ChoiceList (with multiple selections) '" & ((_String*)parameters(1))->sData & "' using redirected stdin input");
                                    return;
                                }
                            } else {
                                break;
                            }
                        }
                }
            } else {
#ifdef  __HEADLESS__
                WarnError ("Unhandled request for data from standard input in ChoiceList in headless HyPhy");
                return;
#else
#if defined __HYPHYQT__
                SetStatusLine ("Waiting for user selection.");
                _String* param = (_String*)parameters(1);

                _SimpleList std(2,0,1),
                            all(theChoices->lLength,0,1);

                choice = HandleListSelection (*theChoices,std, all, *param, sel,fixedLength,(Ptr)_hyPrimaryConsoleWindow);
#else
                _String* param = (_String*)parameters(1);
                printf ("\n\n\t\t\t+");

                for (f = 1; f<param->sLength+1; f++) {
                    printf ("-");
                }

                printf ("+\n\t\t\t|%s|\n\t\t\t+",param->getStr());

                for (f = 1; f<param->sLength+1; f++) {
                    printf ("-");
                }

                printf ("+\n\n");


                long  loopits = 1;

                if (fixedLength == 1) {
                    while (choice == -1) {
                        for (choice = 0; choice<theChoices->lLength; choice++) {
                            printf ("\n\t(%ld):[%s] %s",choice+1,((_String*)(*(_List*)(*theChoices)(choice))(0))->getStr(),((_String*)(*(_List*)(*theChoices)(choice))(1))->getStr());
                        }

                        printf ("\n\n Please choose an option (or press q to cancel selection):");
                        _String buffer (StringFromConsole());
                        if (buffer.sData[0] == 'q' || buffer.sData[0] =='Q') {
                            choice = -1;
                            break;
                        }
                        choice = buffer.toNum();
                        if (choice<1 || choice>theChoices->lLength) {
                            choice = -1;
                            if (loopits++ > 10) {
                                FlagError ("Failed to make a valid selection in ChoiceList after 10 tries");
                                return;
                            }
                        } else {
                            choice--;
                        }
                    }
                } else {
                    if (fixedLength>0)
                        while (sel.lLength<fixedLength) {
                            for (choice = 0; choice<theChoices->lLength; choice++) {
                                if (sel.Find(choice)==-1) {
                                    printf ("\n\t(%ld):%s",choice+1,((_String*)(*(_List*)(*theChoices)(choice))(1))->getStr());
                                }
                            }
                            printf ("\n\n Please choose option %ld of %ld (or press q to cancel selection):",sel.lLength+1,fixedLength);
                            _String buffer (StringFromConsole());
                            if (buffer.sData[0] == 'q' || buffer.sData[0] =='Q') {
                                choice = -1;
                                break;
                            }
                            choice = buffer.toNum();
                            if ((choice>=1)&&(choice<=theChoices->lLength)) {
                                if (sel.Find(choice)==-1) {
                                    sel<<choice-1;
                                }
                            } else {
                                if (loopits++ > 10) {
                                    FlagError ("Failed to make a valid selection in ChoiceList after 10 tries");
                                    return;
                                }
                            }
                        }
                    else
                        while (1) {
                            for (choice = 0; choice<theChoices->lLength; choice++) {
                                if (sel.Find(choice)==-1) {
                                    printf ("\n\t(%ld):[%s] %s",choice+1,((_String*)(*(_List*)(*theChoices)(choice))(0))->getStr(),((_String*)(*(_List*)(*theChoices)(choice))(1))->getStr());
                                }
                            }
                            printf ("\n\n Please choose option %ld, enter d to complete selection, enter q to cancel:",sel.lLength+1);
                            _String buffer (StringFromConsole());
                            if (buffer.sData[0] == 'q' || buffer.sData[0] =='Q') {
                                choice = -1;
                                break;
                            }
                            if (buffer.sData[0] == 'd' || buffer.sData[0] =='D') {
                                break;
                            }

                            choice = buffer.toNum();
                            if ((choice>=1)&&(choice<=theChoices->lLength)) {
                                if (sel.Find(choice)==-1) {
                                    sel<<choice-1;
                                }
                            } else {
                                if (loopits++ > 10) {
                                    FlagError ("Failed to make a valid selection in ChoiceList after 10 tries");
                                    return;
                                }
                            }
                        }
                }
#endif
#endif
            }

            _Variable* sStrV = CheckReceptacle(&selectionStrings,empty,false);

            if (fixedLength == 1) {
                if (choice>=0) {
                    _FString  choiceString (*(_String*) ((_List*)(*theChoices)(choice))->lData[0]);
                    sStrV->SetValue (&choiceString);
                    for (long k=0; k<exclusions.lLength; k++) {
                        if (choice>=exclusions.lData[k]) {
                            choice++;
                        } else {
                            break;
                        }
                    }
                }
                holder->SetValue (new _Constant (choice), false);
            } else {
                if (fixedLength == 0) {
                    fixedLength = sel.lLength;
                    if (fixedLength == 0) {
                        fixedLength = 1;
                    }
                }
                sel.Sort();
                _Matrix   selVector (1,fixedLength,false,true);
                _Matrix   selMatrix (1,fixedLength,false,true);
                if (choice == -1) {
                    selVector[0]=-1;
                } else {
                    for (f=0; f<fixedLength; f++) {
                        choice=sel.lData[f];

                        _FString  *choiceString = new _FString ((*(_String*) ((_List*)(*theChoices)(choice))->lData[0]));
                        _Formula  sf (choiceString);
                        selMatrix.MStore(0,f,sf);
                        for (long k=0; k<exclusions.lLength; k++) {
                            if (choice>=exclusions.lData[k]) {
                                choice++;
                            } else {
                                break;
                            }
                        }
                        selVector[f]=choice;
                        //DeleteObject (choiceString);
                    }
                    sStrV->SetValue (&selMatrix);
                }
                holder->SetValue (&selVector);
            }

            if (choice<0) {
                terminateExecution = true;
            }
        }
    } else {
        WarnError  ("List of selections is invalid in ChoiceList");
    }

    if (simpleParameters.lData[0]) {
        parameters.Delete(4);
        parameters&& &saveTheArg;
    }
}


void      _ElementaryCommand::ExecuteCase36 (_ExecutionList& chain)
{
    chain.currentCommand++;
    _String *currentArgument = (_String*)parameters(0),
             errMsg,
             result;

    long    f = dataSetNamesList.Find(&AppendContainerName(*currentArgument,chain.nameSpacePrefix)),
            s,
            k,
            m;

    if (f<0) {
        ReportWarning (*currentArgument & " is not a valid data set in call to OpenDataPanel");
        return;
    }
    _DataSet* theDS = (_DataSet*)dataSetList(f);

    // process species list
    result = ProcessLiteralArgument ((_String*)parameters(1),chain.nameSpacePrefix);
    if (result.sLength) {
        result.Insert ('"',0);
        result.Insert ('"',-1);
    }
    _SimpleList   speciesList;
    theDS->ProcessPartition (result,speciesList,true);

    // check the validity of list entries
    s = theDS->NoOfSpecies();

    for (m=speciesList.lLength-1; m>=0; m--) {
        k = speciesList.lData[m];
        if (k<0 || k>=s) {
            speciesList.Delete(m);
            m--;
        }
    }

    if (speciesList.lLength==s) {
        speciesList.Clear();
    }

#if !defined __UNIX__ && !defined __HEADLESS__

    _HYDataPanel*  newDP = new _HYDataPanel (empty,empty);
    if (speciesList.lLength) {
        newDP->SetDataSetReference (*(_String*)dataSetNamesList(f),&speciesList);
    } else {
        newDP->SetDataSetReference (*(_String*)dataSetNamesList(f),nil);
    }
    result = dataPanelSourcePath;
    result = ProcessLiteralArgument (&result,chain.nameSpacePrefix);
    if (result.sLength) {
        newDP->SetFilePath (result);
    }
    currentArgument = (_String*)parameters(3);
    newDP->RestorePartInfo(currentArgument);
    currentArgument = (_String*)parameters(2);
    newDP->RestorePanelSettings(currentArgument);
    currentArgument = (_String*)parameters(1);
    newDP->SetSavePath (chain.sourceFile);
    newDP->BringToFront();
    if (parameters.lLength>4) {
        newDP->BuildLikelihoodFunction ((_String*)parameters(4));
        newDP->RestoreSavedLFs ();
    }
#endif
}


//____________________________________________________________________________________
// GetInformation()
void      _ElementaryCommand::ExecuteCase37 (_ExecutionList& chain)
{
    chain.currentCommand++;

    _String matrixName = chain.AddNameSpaceToID(*(_String*)parameters(0)),
            *objectName = (_String*)parameters(1);


    long    sID;
    if (parameters.lLength > 2) {
        sID = ProcessNumericArgument ((_String*)parameters(2), chain.nameSpacePrefix);
    }


    _Matrix *result = nil;

    // object is a non-empty string
    if (objectName->sLength > 2 && objectName->sData[0] == '"' && objectName->sData[objectName->sLength-1] == '"')
        // regular expression
    {
        _String regExp = GetStringFromFormula (objectName,chain.nameSpacePrefix);
        int errNo = 0;
        Ptr regex = PrepRegExp (&regExp, errNo, true);
        if (regex) {
            _List       matches;

            _SimpleList tcache;
            long        iv,
                        k = variableNames.Traverser (tcache, iv, variableNames.GetRoot());

            for (; k>=0; k = variableNames.Traverser (tcache, iv)) {
                _String* vName = (_String*)variableNames.Retrieve (k);
                _SimpleList mtch;
                vName->RegExpMatch (regex,mtch);
                if (mtch.lLength) {
                    matches << vName;
                }

            }

            if (matches.lLength) {
                result = new _Matrix (matches);
            }

            FlushRegExp (regex);
        } else {
            WarnError (GetRegExpError (errNo));
        }
    } else {    // object is not a string, is some kind of variable
        _String objectNameID = chain.AddNameSpaceToID(*objectName);
        long    f = LocateVarByName (objectNameID);
        if      (f>=0) {    // it's a numeric variable
            _Variable* theObject = FetchVar(f);
            if (theObject->ObjectClass()==STRING) {
                objectNameID = _String((_String*)theObject->Compute()->toStr());
                theObject    = FetchVar (LocateVarByName (objectNameID));
            }
            if (theObject) {
                if (theObject->IsCategory()) {
                    _CategoryVariable * thisCV = (_CategoryVariable*)theObject;
                    thisCV->Refresh();

                    _Matrix *values  = thisCV->GetValues(),
                             *weights = thisCV->GetWeights(!thisCV->IsUncorrelated());

                    f = values->GetHDim()*values->GetVDim();
                    checkPointer (result = new _Matrix (2,f,false,true));

                    for (long k = 0; k<f; k++) {
                        result->theData[k]   = values->theData[k];
                        result->theData[f+k] = weights->theData[k];
                    }
                } else {
                    if (theObject->ObjectClass()==TREE_NODE) {
                        _CalcNode* theNode = (_CalcNode*)theObject;
                        if (theNode->GetModelIndex() != HY_NO_MODEL) {
                            checkPointer(result = new _Matrix);
                            theNode->RecomputeMatrix (0,1,result);
                        }
                    }

                    if ((!result)&& theObject->ObjectClass()==NUMBER) {
                        checkPointer(result = new _Matrix (1,3,false,true));
                        result->theData[0]=theObject->Compute()->Value();
                        result->theData[1]=theObject->GetLowerBound();
                        result->theData[2]=theObject->GetUpperBound();
                    }
                }
            }
        } else {
            f = likeFuncNamesList.Find (&objectNameID);
            if (f>=0) {     // it's a likelihood function
                _LikelihoodFunction * lf = (_LikelihoodFunction*)likeFuncList (f);
                f = lf->GetCategoryVars().lLength;
                if (f==0) {
                    f++;
                }

                _List        catVars;

                for (long k=0; k<lf->GetCategoryVars().lLength; k++) {
                    _String varName = *LocateVar(lf->GetCategoryVars().lData[k])->GetName();
                    catVars && & varName;
                }

                result = (_Matrix*) checkPointer(new _Matrix (catVars));
            } else {
				if ((f = dataSetFilterNamesList.Find (&objectNameID))>=0)
					// return a vector of strings - each with actual characters of the corresponding sequence
				{
					_DataSetFilter* daFilter = (_DataSetFilter*)dataSetFilterList (f);
					result = daFilter->GetFilterCharacters();
				} else {
					// it's a tree node with a rate matrix assigned
					f = FindModelName (objectNameID);
					if (f>=0)
						// for models, return the list of variables in the model
					{
						_SimpleList modelParms;
						_AVLList    modelParmsA (&modelParms);

						LocateVar (modelMatrixIndices.lData[f])->ScanForVariables(modelParmsA,false);
						_List       modelPNames;

						for (unsigned long vi=0; vi<modelParms.lLength; vi++) {
							modelPNames << LocateVar(modelParms.lData[vi])->GetName();
						}

						result = new _Matrix (modelPNames);
					}
				}
            }
        }
    }

    if (!result) {
        result = new _Matrix (0,0,false,false);
    }

    CheckReceptacleAndStore (&matrixName, empty, true, result, true);
    DeleteObject (result);

}


//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase43 (_ExecutionList& chain)
{
    chain.currentCommand++;

    _String *currentArgument = (_String*)parameters(0),
             result;

    _Variable * theReceptacle = CheckReceptacle(&AppendContainerName(*currentArgument,chain.nameSpacePrefix),code==43?blFindRoot:blIntegrate,true);

    if (theReceptacle) {
        _String     exprString =  *(_String*)parameters(1);
        _Formula    theExpression (exprString);

        currentArgument =   (_String*)parameters(2);
        long f = LocateVarByName (AppendContainerName(*currentArgument,chain.nameSpacePrefix));

        if (f<0) {
            ReportWarning (*currentArgument & " is not an existing variable to solve for in call to FindRoot/Integrate.");
            return;
        }

        if (terminateExecution) {
            return;
        }

        _Formula * dF = (code==43)?theExpression.Differentiate (*(_String*)parameters(2),false):nil;

        _Parameter    lb = ProcessNumericArgument ((_String*)parameters(3),chain.nameSpacePrefix),
                      ub = ProcessNumericArgument ((_String*)parameters(4),chain.nameSpacePrefix);

        if (ub<=lb && code==48) {
            ReportWarning (_String ('[') & lb & ',' & ub & "] is not a valid search interval in call to FindRoot/Integrate");
            return;
        }

        if (code==43) {
            if (dF) {
                theReceptacle->SetValue (new _Constant (theExpression.Newton (*dF,FetchVar (f), 0.0, lb, ub)),false);
            } else {
                theReceptacle->SetValue (new _Constant (theExpression.Brent (FetchVar(f), lb, ub)), false);
            }
        } else {
            theReceptacle->SetValue (new _Constant (theExpression.Integral (FetchVar (f), lb, ub, ub-lb>1e10)), false);
        }

        if (dF) {
            delete (dF);
        }
    }
}

//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase44 (_ExecutionList& chain)
{
    chain.currentCommand++;

#ifdef __HYPHYMPI__
    _String *arg1 = (_String*)parameters(0),
             *arg2 = (_String*)parameters(1),
              *arg3 = parameters.lLength>2?(_String*)parameters(2):nil,
               *theMessage = nil;


    _Parameter      nodeCount;
    checkParameter (mpiNodeCount,nodeCount,1);

    long            destID = ProcessNumericArgument (arg1,chain.nameSpacePrefix),
                    g;

    if (!numericalParameterSuccessFlag || destID<0 || destID>=nodeCount) {
        WarnError (*arg1 & " is not a valid MPI node ID in call to MPISend.");
        return;
    }

    if (arg3) {
        _AssociativeList * ar = (_AssociativeList *)FetchObjectFromVariableByType (&AppendContainerName(*arg3,chain.nameSpacePrefix), ASSOCIATIVE_LIST);
        if (!ar) {
            WarnError (*arg3 & " is not a valid associative array for input options in call to MPISend.");
            return;
        }
        theMessage = new _String (256L, true);
        checkPointer (theMessage);
        _String arrayID ("_HYPHY_MPI_INPUT_ARRAY_");
        (*theMessage) << arrayID;
        (*theMessage) << '=';
        arg3 = ar->Serialize (arrayID);
        (*theMessage) << arg3;
        DeleteObject (arg3);
        (*theMessage) << ';';
        arrayID = *arg2;
        arrayID.ProcessFileName(false,true,(Ptr)chain.nameSpacePrefix);
        (*theMessage) << "\nExecuteAFile (\"";
        (*theMessage) << arrayID;
        (*theMessage) << "\",_HYPHY_MPI_INPUT_ARRAY_);";
        theMessage->Finalize();
    } else if ((g=FindLikeFuncName(AppendContainerName(*arg2,chain.nameSpacePrefix)))>=0) {
        checkPointer (theMessage = new _String(1024L,true));
        ((_LikelihoodFunction*)likeFuncList(g))->SerializeLF(*theMessage,_hyphyLFSerializeModeOptimize);
        theMessage->Finalize();
    } else {
        theMessage = new _String (ProcessLiteralArgument (arg2,chain.nameSpacePrefix));
    }

    if (theMessage == nil || theMessage->sLength==0) {
        WarnError (*arg2 & " is not a valid (or is an empty) string (LF ID) in call to MPISend.");
    } else {
        MPISendString (*theMessage, destID);
    }

    DeleteObject (theMessage);

#else
    WarnError ("MPISend can't be used by non-MPI versions of HyPhy.");
#endif

}

//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase45 (_ExecutionList& chain)
{
    chain.currentCommand++;

#ifdef __HYPHYMPI__
    _String *arg1 = (_String*)parameters(0),
             *arg2 = (_String*)parameters(1),
              *arg3 = (_String*)parameters(2);

    _Parameter      nodeCount;
    checkParameter (mpiNodeCount,nodeCount,1);

    long            srcT = ProcessNumericArgument (arg1,chain.nameSpacePrefix),
                    srcID,
                    g;

    if ((!numericalParameterSuccessFlag)||(srcT<-1)||(srcT>=nodeCount)) {
        WarnError (*arg1 & " is not a valid MPI node ID in call to MPIReceive.");
        return;
    }

    _Variable* idVar = CheckReceptacle (&AppendContainerName(*arg2,chain.nameSpacePrefix),"MPIReceive"),
               * mVar  = CheckReceptacle (&AppendContainerName(*arg3,chain.nameSpacePrefix),"MPIReceive");

    if (!(idVar&&mVar)) {
        return;
    }

    _FString* theMV = new _FString (MPIRecvString (srcT,srcID));
    checkPointer (theMV);
    idVar->SetValue (new _Constant (srcID),false);
    mVar->SetValue (theMV, false);
#else
    WarnError ("MPIReceive can't be used by non-MPI versions of HyPhy.");
#endif

}

//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase46 (_ExecutionList& chain)
{
    chain.currentCommand++;

    _String *arg1 = (_String*)parameters(1),
             *arg2 = (_String*)parameters(0),
              errMsg;

    long    k = dataSetFilterNamesList.Find (&AppendContainerName(*arg1,chain.nameSpacePrefix));

    if (k<0) {
        errMsg = *arg1 & " is not a defined data set filter ID ";
    } else {
        _DataSetFilter * dsf   = (_DataSetFilter*)dataSetFilterList (k);
        _Variable *      stVar = CheckReceptacle(&AppendContainerName(*arg2,chain.nameSpacePrefix),"GetDataInfo");

        if (stVar) {
            if (parameters.lLength == 2) {
                _Matrix * res = new _Matrix (1,dsf->duplicateMap.lLength, false, true);
                checkPointer (res);
                for (k = 0; k<dsf->duplicateMap.lLength; k++) {
                    res->theData[k] = dsf->duplicateMap.lData[k];
                }
                stVar->SetValue (res,false);
            } else {
                if (parameters.lLength == 3) {
                    _String checker = ProcessLiteralArgument ((_String*)parameters(2),chain.nameSpacePrefix);
                    if (checker == _String ("CHARACTERS")) {
                        _List   characters;
                        k       = dsf->GetDimension(true);
                        long fd = dsf->GetUnitLength();
                        for (long idx = 0; idx < k; idx++) {
                            characters.AppendNewInstance(new _String (dsf->ConvertCodeToLetters (dsf->CorrectCode(idx), fd)));
                        }

                        stVar->SetValue (new _Matrix (characters), false);
                    } else if (checker == _String ("PARAMETERS")) {
                        _AssociativeList * parameterInfo = new _AssociativeList;
                        parameterInfo->MStore ("ATOM_SIZE",             new _Constant (dsf->GetUnitLength()), false);
                        parameterInfo->MStore ("EXCLUSIONS",            new _FString  (dsf->GetExclusions()), false);
                        parameterInfo->MStore ("SITES_STRING",          new _FString  ((_String*)dsf->theOriginalOrder.ListToPartitionString()), false);
                        parameterInfo->MStore ("SEQUENCES_STRING",      new _FString  ((_String*)dsf->theNodeMap.ListToPartitionString()), false);
                        stVar->SetValue (parameterInfo,false);

                    } else if (checker == _String ("CONSENSUS")) {
                        stVar->SetValue (new _FString (new _String(dsf->GenerateConsensusString())), false);
                    } else {
                        long seqID = ProcessNumericArgument ((_String*)parameters(2),chain.nameSpacePrefix);
                        if (seqID>=0 && seqID < dsf->NumberSpecies()) {
                            stVar->SetValue (new _FString (dsf->GetSequenceCharacters(seqID)),false);
                        } else {
                            // 20110916 SLKP : the option for filtering duplicate sequences
                            if (seqID >= -4 && seqID <= -1) {
                                _SimpleList indices, map, counts;
                                long uniqueSequences = dsf->FindUniqueSequences(indices, map, counts, -seqID - 1);
                                _AssociativeList * parameterInfo = new _AssociativeList;
                                parameterInfo->MStore ("UNIQUE_SEQUENCES",             new _Constant (uniqueSequences), false);
                                parameterInfo->MStore ("UNIQUE_INDICES",            new _Matrix (indices), false);
                                parameterInfo->MStore ("SEQUENCE_MAP",          new _Matrix (map), false);
                                parameterInfo->MStore ("UNIQUE_COUNTS",      new _Matrix  (counts), false);
                                stVar->SetValue (parameterInfo,false);
                            }
                        }
                    }
                } else {
                    long seq  = ProcessNumericArgument ((_String*)parameters(2),chain.nameSpacePrefix),
                         site = ProcessNumericArgument ((_String*)parameters(3),chain.nameSpacePrefix);

                    if (parameters.lLength == 4) {
                        if ((seq>=0)&&(site>=0)&&(seq<dsf->NumberSpecies())&&(site<dsf->NumberDistinctSites())) {
                            _Matrix             * res = (_Matrix*)checkPointer(new _Matrix (dsf->GetDimension (true), 1, false, true));

                            _Parameter          onlyTheIndex = 0.0;
                            checkParameter      (getDataInfoReturnsOnlyTheIndex,onlyTheIndex,0.0);

                            long                theValue = dsf->Translate2Frequencies ((*dsf)(site,seq), res->theData,  true);

                            if (onlyTheIndex > 0.5) {
                                stVar->SetValue (new _Constant (theValue),false);
                                DeleteObject     (res);
                            } else {
                                stVar->SetValue (res,false);
                            }
                        } else {
                            errMsg = _String (seq) & "," & _String (site) & " is an invalid site index ";
                        }
                    } else {
                        if ((seq>=0)&&(site>=0)&&(seq<dsf->NumberSpecies())&&(site<dsf->NumberSpecies())) {
                            _String* resFlag = (_String*)parameters(4);
                            _Matrix * res;

                            if (pcAmbiguitiesAverage.Equal (resFlag)) {
                                res = dsf->ComputePairwiseDifferences (seq,site,1);
                            } else if (pcAmbiguitiesResolve.Equal (resFlag)) {
                                res = dsf->ComputePairwiseDifferences (seq,site,2);
                            } else if (pcAmbiguitiesSkip.Equal (resFlag)) {
                                res = dsf->ComputePairwiseDifferences (seq,site,3);
                            } else {
                                res = dsf->ComputePairwiseDifferences (seq,site,0);
                            }

                            stVar->SetValue (res,false);
                        } else {
                            errMsg = _String (seq) & "," & _String (site) & " is an invalid sequence pair specification.";
                        }

                    }
                }
            }
        }
    }

    if (errMsg.sLength) {
        errMsg = errMsg & " in call to GetDataInfo ";
        WarnError (errMsg);
    }
}

//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase47 (_ExecutionList& chain)
{
    chain.currentCommand++;

    _String *arg1 = (_String*)parameters(0),
             *arg2 = (_String*)parameters(1),
              errMsg;

    long    k = FindLikeFuncName(AppendContainerName(*arg1, chain.nameSpacePrefix));

    if (k<0) {
        _String  litArg = ProcessLiteralArgument (arg1,chain.nameSpacePrefix);
        k = FindLikeFuncName (litArg);
        if (k<0) {
            errMsg = *arg1 & " is not a defined likelihood function ID ";
        }
    }

    if (errMsg.sLength == 0) {
        _LikelihoodFunction * lf   = (_LikelihoodFunction *) likeFuncList (k);
        _String         callBack   = ProcessLiteralArgument (arg2,chain.nameSpacePrefix);

        k = batchLanguageFunctionNames.Find (&callBack);

        if (k<0) {
            errMsg = *arg2 & " is not a defined user batch language function ";
        } else {
            if (batchLanguageFunctionParameters.lData[k]!=2) {
                errMsg = *arg2 & " callback function must depend on 2 parameters ";
            } else {
                lf->StateCounter (k);
            }
        }
    }

    if (errMsg.sLength) {
        errMsg = errMsg & " in call to StateCounter.";
        WarnError (errMsg);
    }
}



//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase52 (_ExecutionList& chain)
{

    chain.currentCommand++;

    _String         errMsg ;

    // check validity of an alphabet

    long           siteCount  = ProcessNumericArgument ((_String*)parameters (4),chain.nameSpacePrefix);
    _String        givenState;

    if (siteCount < 1) {
        givenState = ProcessLiteralArgument((_String*)parameters (4),chain.nameSpacePrefix);
        siteCount = givenState.sLength;
    }

    if (siteCount < 1) {
        errMsg = *(_String*)parameters (4) & " must either evaluate to a positive integer or be a non-empty string of root states";
        WarnError (errMsg);
        return;
    }

    _Variable   *  alphabet = FetchVar (LocateVarByName (AppendContainerName(*(_String*)parameters (3),chain.nameSpacePrefix))),
                *  treeVar  = FetchVar (LocateVarByName (AppendContainerName(*(_String*)parameters (1),chain.nameSpacePrefix))),
                *  freqVar  = FetchVar (LocateVarByName (AppendContainerName(*(_String*)parameters (2),chain.nameSpacePrefix)));


    if (alphabet&&treeVar&&freqVar) {
        if (alphabet->ObjectClass() == MATRIX) {
            _Matrix * alphabetMatrix = (_Matrix*)alphabet->GetValue();

            if (alphabetMatrix->IsAStringMatrix() && alphabetMatrix->GetHDim() == 2 && alphabetMatrix->GetVDim () > 1) {
                _String baseSet;

                for (long k=0; k < alphabetMatrix->GetVDim (); k++) {
                    _FString * aState = (_FString*)alphabetMatrix->GetFormula(0,k)->Compute();
                    if (aState) {
                        if (aState->theString->sLength == 1) {
                            char c = aState->theString->sData[0];
                            if (baseSet.Find(c) == -1) {
                                baseSet = baseSet & c;
                            } else {
                                break;
                            }
                        } else {
                            break;
                        }
                    } else {
                        break;
                    }
                }

                if (baseSet.sLength == alphabetMatrix->GetVDim ()) {
                    long unitSize = ((_FString*)alphabetMatrix->GetFormula(1,0)->Compute())->theString->toNum();

                    if (unitSize >= 1) {
                        _Formula* exclusionFormula = alphabetMatrix->GetFormula(1,1);
                        _String* theExclusions = &empty;
                        
                        if (exclusionFormula)
                            theExclusions = ((_FString*)exclusionFormula->Compute())->theString;

                        if (treeVar->ObjectClass() == TREE) {
                            if (freqVar->ObjectClass() == MATRIX) {
                                _TheTree * spawningTree = (_TheTree*)treeVar;

                                if (!(parameters.lLength>6 && (spawningTree->CountTreeCategories()>1))) {

                                    if (givenState.sLength>1)
                                        // root state
                                    {
                                        if ((givenState.sLength >= unitSize)&&(givenState.sLength % unitSize == 0)) {
                                            siteCount                   = givenState.sLength/unitSize;
                                        } else {
                                            errMsg = "Root state string is either too short or has length which is not divisible by the unit size";
                                        }
                                    }
                                    if (errMsg.sLength == 0) {
                                        _TranslationTable newTT (baseSet);
                                        _DataSet * ds = (_DataSet*)checkPointer(new _DataSet);

                                        if (! newTT.IsStandardNucleotide() ) {
                                            ds->SetTranslationTable (&newTT);    // mod 20060113 to properly deal with non-standard alphabets
                                        }
                                        // make a dummy
                                        spawningTree->AddNodeNamesToDS (ds,true,false,1);

                                        char    c = baseSet.sData[0];
                                        long    s = ds->GetNames().lLength;

                                        if (s<2) {
                                            _String rt ("Root");
                                            ds->GetNames().InsertElement (&rt,0,true);
                                            s ++;
                                        }

                                        unsigned long ssi = _String::storageIncrement;

                                        if (s>ssi) {
                                            _String::storageIncrement = s;
                                        }

                                        ds->AddSite(c);
                                        for (long u = 1; u < s; u++) {
                                            ds->Write2Site(0,c);
                                        }
                                        ds->Finalize();

                                        _String::storageIncrement = ssi;
                                        ds->SetNoSpecies (s);

                                        _SimpleList * theMap = & ds->GetTheMap();
                                        theMap->RequestSpace (siteCount*unitSize);
                                        for (long filler = 0; filler < siteCount*unitSize; filler++) {
                                            theMap->lData[filler] = 0;
                                        }

                                        theMap->lLength = siteCount*unitSize;

                                        _DataSetFilter* newFilter = new _DataSetFilter();
                                        checkPointer   (newFilter);
                                        _SimpleList     h,v;

                                        newFilter->SetFilter     (ds,unitSize,h,v,false);
                                        newFilter->SetExclusions (theExclusions,true);
                                        newFilter->SetupConversion ();

                                        /*char buffer[255];
                                        sprintf (buffer,"%d %d\n",siteCount, newFilter->GetFullLengthSpecies(),unitSize);
                                        BufferToConsole (buffer);
                                        */
                                        _Matrix*   rootStates = nil;
                                        if (givenState.sLength>=unitSize) {
                                            rootStates                  = new _Matrix (1,siteCount,false,true);
                                            checkPointer                (rootStates);
                                            _Parameter*  holder         = new _Parameter [newFilter->GetDimension(false)];
                                            checkPointer                (holder);

                                            for (long cc = 0; cc < siteCount; cc++) {
                                                _String aState (givenState.Cut(cc*unitSize,(cc+1)*unitSize-1));
                                                long    stateV = newFilter->Translate2Frequencies (aState,holder,false);
                                                if (stateV<0) {
                                                    errMsg = aState & " found in the root state string at position " & cc*unitSize & " is an invalid state";
                                                    break;
                                                } else {
                                                    rootStates->theData[cc] = stateV;
                                                }
                                            }

                                            delete [] holder;
                                        }
                                        if (errMsg.sLength == 0) {

                                            long       filterID = AddFilterToList (simulationFilter,newFilter);

                                            spawningTree->SetUp();
                                            spawningTree->InitializeTreeFrequencies((_Matrix*)freqVar->Compute(),true);
                                            errMsg = *(_String*)dataSetFilterNamesList(filterID) & ',' & *spawningTree->GetName() & ',' & *freqVar->GetName();


                                            _LikelihoodFunction lf (errMsg, nil);

                                            if (terminateExecution) {
                                                return;
                                            }

                                            bool    doInternals = false;

                                            if (parameters.lLength>5) {
                                                doInternals = (ProcessNumericArgument ((_String*)parameters (5),chain.nameSpacePrefix)>0.5);
                                            }


                                            _String spoolFile;

                                            FILE*   mainFile = nil;

                                            errMsg = empty;

                                            if (parameters.lLength > 6) {
                                                spoolFile = ProcessLiteralArgument ((_String*)parameters (6),chain.nameSpacePrefix);
                                                spoolFile.ProcessFileName();
                                                mainFile = doFileOpen (spoolFile.sData,"w");
                                                if (!mainFile) {
                                                    errMsg = _String("Failed to open ") & spoolFile & " for writing";
                                                }
                                                if (doInternals) {
                                                    spoolFile = spoolFile & ".anc";
                                                }
                                            }

                                            if (errMsg.sLength == 0) {
                                                _DataSet    * simDataSet;

                                                if (mainFile) {
                                                    simDataSet = new _DataSet (mainFile);
                                                } else {
                                                    simDataSet = new _DataSet (siteCount);
                                                }

                                                checkPointer (simDataSet);

                                                _List exclusions;

                                                _String *simName = new _String(AppendContainerName(*(_String*)parameters (0),chain.nameSpacePrefix));
                                                _String mxName = *simName & ".rates";
                                                setParameter (mxName, 0.0);
                                                _Variable *catValVar        = FetchVar (LocateVarByName (mxName));
                                                _Matrix*   catValues        = new _Matrix (1,1,false,true);
                                                checkPointer    (catValues);

                                                mxName = *simName & ".rateVars";
                                                setParameter (mxName, 0.0);
                                                _Variable * catNameVar  = FetchVar (LocateVarByName (mxName));
                                                _Matrix* catNames       = new _Matrix (1,1,false,true);

                                                SetStatusLine ("Simulating Data");
                                                lf.Simulate (*simDataSet, exclusions, catValues, catNames, rootStates, doInternals?(mainFile?&spoolFile:&empty):nil);
                                                SetStatusLine ("Idle");

                                                catValVar->SetValue(catValues, false);
                                                catNameVar->SetValue(catNames, false);

                                                StoreADataSet (simDataSet, simName);
                                                DeleteObject (simName);
                                                KillDataFilterRecord (filterID);
                                                errMsg = empty;
                                            }
                                        }
                                        DeleteObject   (ds);
                                        if (rootStates) {
                                            DeleteObject (rootStates);
                                        }

                                        if (errMsg.sLength == 0) {
                                            return;
                                        }
                                    }
                                } else {
                                    errMsg = "Can't use spool to file option in Simulate when the tree depends on category variables.";
                                }
                            } else {
                                errMsg = *(_String*)parameters (2) & " must be an existing matrix";
                            }
                        } else {
                            errMsg = *(_String*)parameters (1) & " must be an existing tree";
                        }
                    } else {
                        errMsg = "Invalid unit length specification (must be >=1)";
                    }
                } else {
                    errMsg = "Invalid alphabet character specification";
                }
            }
        }
        if (errMsg.sLength == 0) {
            errMsg = _String("Alphabet specification variable ") & *(_String*)parameters (3) & " must be a string matrix with 2 rows and at least 2 columns";
        }
    } else {
        long i = 0;
        if (!alphabet) {
            i = 3;
        } else {
            if (!treeVar) {
                i = 1;
            } else if (!freqVar) {
                i = 2;
            }
        }

        errMsg = _String("Variable ") & *(_String*)parameters (i) & " has not been defined";

    }

    if (errMsg.sLength) {
        errMsg = errMsg & " in Simulate.";
        WarnError (errMsg);
    }

}



//____________________________________________________________________________________

bool      _ElementaryCommand::Execute    (_ExecutionList& chain) // perform this command in a given list
{
    _String errMsg;

    switch (code) {

    case 0: // formula reparser
        ExecuteCase0 (chain);
        break;


    case 4:
        ExecuteCase4 (chain);
        break;


    case 5: // data set contruction

        ExecuteCase5 (chain);
        break;

    case 6:  // data set filter construction
    case 27: // Permute
    case 28: // Bootstrap
        ExecuteDataFilterCases(chain);
        break;


    case 7: { // build a tree
        chain.currentCommand++;

        _String treeIdent   = chain.AddNameSpaceToID(*(_String*)parameters(0)),
                treeString  = *(_String*)parameters(1);

        SetStatusLine (_String("Constructing Tree ")&treeIdent);
        long  varID = LocateVarByName (treeIdent);

        _Parameter rtv = 0.0; // mod 11/19/2003
        checkParameter (replaceTreeStructure, rtv, 0.0); // mod 11/19/2003

        _SimpleList   leftOverVars; // mod 02/03/2003
        if (varID>=0)
            if (FetchVar(varID)->ObjectClass()==TREE) {
                if (rtv>0.5) {
                    DeleteVariable(*FetchVar(varID)->GetName());    // mod 11/19/2003
                } else {
                    DeleteTreeVariable(*FetchVar(varID)->GetName(),leftOverVars,true);    // mod 02/03/2003
                }
            }

        treeString.ProcessParameter();

        _TheTree * tr = nil;

        if (treeString.getChar(0)!='(') {
            _Formula  nameForm (treeString,chain.nameSpacePrefix);
            _PMathObj formRes = nameForm.Compute();
            if (formRes) {
                if (formRes->ObjectClass () == STRING) {
                    tr = new _TheTree (treeIdent,*((_FString*)formRes)->theString,false);
                } else if (formRes->ObjectClass () == TOPOLOGY) {
                    tr = new _TheTree (treeIdent,(_TreeTopology*)nameForm.Compute());
                }
            }
        } else {
            tr = new _TheTree (treeIdent,treeString,false);
        }

        if (!tr) {
            WarnError ("Illegal right hand side in call to Tree id = ...; it must be a string, a Newick tree spec or a topology");
            return false;
        }
        //_TheTree tr (*(_String*)parameters(0),treeStr);


        /*for (varID = 0; varID < leftOverVars.lLength; varID++)
        {
            _Variable* theVar = LocateVar (leftOverVars.lData[varID]);
            if (theVar)
                printf ("%d = %s\n", leftOverVars.lData[varID], theVar->GetName()->getStr());
            else
                printf ("%d Deleted!!!\n", leftOverVars.lData[varID], theVar->GetName()->getStr());

        }*/

        if (leftOverVars.lLength) { // mod 02/03/2003 - the entire "if" block
            _SimpleList indep, dep, holder;
            {
                _AVLList    indepA (&indep),
                            depA   (&dep);

                tr->ScanForVariables (indepA,depA);
                //tr.ScanForVariables (indepA,depA);
                indepA.ReorderList();
                depA.ReorderList();
            }

            //indep.Sort();
            //dep.Sort();

            holder.Union (indep,dep);
            leftOverVars.Sort ();
            indep.Subtract (leftOverVars,holder);

            /* the bit with freeSlots is here b/c
               some nodes variables may have been deleted during the unroot
               in the tree constructor and we don't want to delete them twice,
               do we? 08/22/2003 */

            dep.Clear();
            dep.Duplicate (&freeSlots);
            dep.Sort ();
            holder.Subtract (indep,dep);
            for (varID = holder.lLength-1; varID >=0 ; varID--) {
                DeleteVariable (*LocateVar (holder.lData[varID])->GetName());
            }

            tr->Clear();

        }
        SetStatusLine ("Idle");

    }
    break;

    case 8: // print stuff to file (or stdout)
        ExecuteCase8(chain);
        break;

    case HY_HBL_COMMAND_HARVEST_FREQUENCIES: { // or HarvestFrequencies
        return HandleHarvestFrequencies(chain);
    }

    case HY_HBL_COMMAND_OPTIMIZE: // optimize the likelihood function
    case HY_HBL_COMMAND_COVARIANCE_MATRIX: {
        return HandleOptimizeCovarianceMatrix (chain, code == HY_HBL_COMMAND_OPTIMIZE);
    }

    case 11: // build the likelihood function

        ExecuteCase11 (chain);
        break;

    case 12: // data set contruction by simulation

        ExecuteCase12 (chain);
        break;

    case 14:

    {
        if (parameters.lLength) {
            DeleteObject (chain.result);
            _Formula returnValue (*(_String*)parameters(0),chain.nameSpacePrefix);
            chain.result = returnValue.Compute();
            if (chain.result) {
                chain.result = (_PMathObj) chain.result->makeDynamic();
            }
        }
        chain.currentCommand = simpleParameters(0);
        if (chain.currentCommand<0) {
            chain.currentCommand = 0x7fffffff;
        }
    }

    break;

    case 16: { // data set merger operation
        chain.currentCommand++;
        SetStatusLine ("Merging Datasets");
        _SimpleList     dsIndex;
        for (long di=1; di<parameters.lLength; di++) {
            _String  dsname = chain.AddNameSpaceToID(*(_String*)parameters(di));
            long f = FindDataSetName (dsname);
            if (f==-1) {
                WarnError (((_String)("Identifier ")&dsname&_String(" doesn't correspond to a valid dataset.")));
                return false;
            } else {
                dsIndex<<f;
            }
        }

        _DataSet*  mergeResult = (simpleParameters(0)==1 || simpleParameters(0)==-1)?_DataSet::Concatenate(dsIndex):_DataSet::Combine(dsIndex);
        // xlc mod 03/08/2005
        _String  * resultName = new _String (chain.AddNameSpaceToID(*(_String*)parameters(0)));

        if (StoreADataSet (mergeResult, resultName) && simpleParameters(0)<0)
            // purge all the datasets except the resulting one
        {
            long newSetID = FindDataSetName (*resultName);
            for (long di=0; di<dsIndex.lLength; di++)
                if (dsIndex.lData[di] != newSetID) {
                    KillDataSetRecord(dsIndex.lData[di]);
                }
        }

        DeleteObject (resultName);

    }
    break;

    case HY_HBL_COMMAND_EXPORT: // matrix export operation
        HandleExport (chain);
        break;

    case 18: // matrix import operation

    {
        bool importResult = true;
        chain.currentCommand++;
        _String  fName (*(_String*)parameters(1));
        fName.ProcessFileName();
        if (terminateExecution) {
            return false;
        }
        FILE*   theDump = doFileOpen (fName.getStr(),"rb");
        if (!theDump) {
            WarnError (((_String)("File ")&fName&_String(" couldn't be open for reading.")));
            return false;
        }

        fName = chain.AddNameSpaceToID(*(_String*)parameters(0));
        _Variable * result  = CheckReceptacle(&fName,blImport.Cut(0,blImport.sLength-2),true);
        if (result) {
            _Matrix   * storage = new _Matrix (1,1,false,true);
            result->SetValue(storage,false);
            lastMatrixDeclared = result->GetAVariable();
            if (!storage->ImportMatrixExp(theDump)) {
                WarnError ("Matrix import failed - the file has an invalid format.");
                importResult = false;
                DeleteObject(storage);
            }
        } else {
            importResult = false;
        }
        fclose (theDump);
        return importResult;
    }
    break;

    case HY_HBL_COMMAND_MOLECULAR_CLOCK: // molecular_clock constraint
        HandleMolecularClock(chain);
        break;

    case 20: // category variable construction

    {
        chain.currentCommand++;
        _String cName = chain.AddNameSpaceToID (*(_String*)parameters(0));
        _List parms (parameters);
        parms.Delete(0);
        _CategoryVariable newCat(cName,&parms,chain.nameSpacePrefix);
        ReplaceVar(&newCat);
    }
    break;

    case 21: // construct the category matrix
        ExecuteCase21 (chain);
        break;

    case HY_HBL_COMMAND_CLEAR_CONSTRAINTS:  // clear constraints
        HandleClearConstraints(chain);
        break;

    case HY_HBL_COMMAND_SET_DIALOG_PROMPT: { // set dialog prompt
        chain.currentCommand++;
        dialogPrompt = ProcessLiteralArgument((_String*)parameters(0),chain.nameSpacePrefix);
    }
    break;

    case HY_HBL_COMMAND_SELECT_TEMPLATE_MODEL: // prompt for a model file
        return HandleSelectTemplateModel(chain);

    case 25: // fscanf
    case 56: // sscanf

        ExecuteCase25 (chain,code == 56);
        break;

    case 26: // replicate constraint
        ExecuteCase26 (chain);
        break;

    case HY_HBL_COMMAND_USE_MODEL:
        return HandleUseModel(chain);

    case 31:
        ExecuteCase31 (chain);
        break;

    case 32:
        ExecuteCase32 (chain);
        break;

    case HY_HBL_COMMAND_GET_STRING:
        HandleGetString (chain);
        break;

    case HY_HBL_COMMAND_SET_PARAMETER:
        return HandleSetParameter(chain);

    case 36:
        ExecuteCase36 (chain);
        break;

    case 37:
        ExecuteCase37 (chain);
        break;

    case 38:
        ExecuteCase38 (chain, false);
        break;

    case 39:
    case 62:
    case 66:
        ExecuteCase39 (chain);
        break;

    case 40:
        ExecuteCase40 (chain);
        break;

    case 41:
        ExecuteCase41 (chain);
        break;

    case HY_HBL_COMMAND_DIFFERENTIATE:
        return HandleDifferentiate (chain);
        break;

    case 43:
        ExecuteCase43 (chain);
        break;

    case 44:
        ExecuteCase44 (chain);
        break;

    case 45:
        ExecuteCase45 (chain);
        break;

    case 46:
        ExecuteCase46 (chain);
        break;

    case 47:
        ExecuteCase47 (chain);
        break;

    case 48:
        ExecuteCase43 (chain);
        break;

    case HY_HBL_COMMAND_LFCOMPUTE:
        return HandleComputeLFFunction(chain);

    case 50:
        ExecuteCase38 (chain, true);
        break;

    case HY_HBL_COMMAND_GET_URL:
        return HandleGetURL (chain);

    case 52:
        ExecuteCase52 (chain);
        break;

    case 53:
        ExecuteCase53 (chain);
        break;

    case 54:
        ExecuteCase54 (chain);
        break;

    case 55:
        ExecuteCase55 (chain);
        break;

    case 57:
        ExecuteCase57 (chain);
        break;

    case 58:
        ExecuteCase58 (chain);
        break;

    case HY_HBL_COMMAND_DELETE_OBJECT:
        HandleDeleteObject (chain);
        break;

    case HY_HBL_COMMAND_REQUIRE_VERSION:
        HandleRequireVersion(chain);
        break;

    case 61:
        ExecuteCase61 (chain);
        break;

    case 63:
        ExecuteCase63 (chain);
        break;

    case 64:
        ExecuteCase64 (chain);
        break;

    case HY_HBL_COMMAND_ASSERT:
        HandleAssert (chain);
        break;

    default:
        chain.currentCommand++;
    }

    return true;
}

//____________________________________________________________________________________


_String   _ElementaryCommand::FindNextCommand  (_String& input, bool useSoftTrim)
{

    bool    isString  = false,
            skipping  = false;

    char    isComment = 0;


    long    scopeIn     = 0,
            matrixScope = 0,
            parenIn  = 0,
            bracketIn   = 0,
            index,
            saveSI = _String::storageIncrement;

    _SimpleList isDoWhileLoop;

    if (input.sLength/4 > saveSI) {
        _String::storageIncrement = input.sLength/4;
    }

    _String result (128L,true);

    char    lastChar = 0;

    index = input.Length();

    if (!index) {
        result.Finalize();
        return empty;
    }

    // non printable characters at the end ?
    while (index>=0 && !isprint(input[--index])) ;
    input.Trim (0,index, useSoftTrim);

    for (index = 0; index<input.Length(); index++) {
        char c = input.sData[index];

        if (!isString && c=='\t') {
            c = ' ';
        }

        // check for comments
        if (isComment) {
            if (isComment == 1) {
                if (c=='/' && input.sData[index-1]=='*') {
                    isComment = 0;
                }
            } else if (c == '\r' || c == '\n') {
                isComment = 0;
            }

            lastChar  = 0;
            continue;
        } else {
            if (!isString && c=='/') {
                switch (input.getChar(index+1)) {
                case '*':
                    isComment = 1;
                    break;
                case '/':
                    isComment = 2;
                }

                if (isComment) {
                    lastChar  = 0;
                    index++;
                    continue;
                }
            }
        }


        // skip spaces
        if (!isString && isspace(c)) {
            if (index >= 6 && input.getChar(index-1) == 'n'
                    && input.getChar(index-2) == 'r'
                    && input.getChar(index-3) == 'u'
                    && input.getChar(index-4) == 't'
                    && input.getChar(index-5) == 'e'
                    && input.getChar(index-6) == 'r') {
                if (index == 6 || (index > 6 && !(isalnum(input.getChar(index-7)) || input.getChar(index-7) == '_'))) {
                    result<<' ';
                }
            }


            skipping = true;
            continue;
        }

        if (skipping&&( isalpha(c) || c=='_') && (isalnum(lastChar) || lastChar=='_')) {
            result<<' ';
        }

        skipping = false;

        result<<c;

        if (isString && c == '\\') {
            result<< input.getChar(++index);
            continue;
        }

        // are we inside a string literal?

        if (c=='"') {
            isString = !isString;
            lastChar = 0;
            continue;
        }

        if (isString) {
            continue;
        }

        // maybe we are done?

        if (c==';' && scopeIn == 0 && matrixScope == 0 && parenIn <= 0 && bracketIn <= 0) {
            break;
        }

        // check to see whether we are defining a matrix

        if (c=='(') {
            parenIn++;
            lastChar = 0;
            continue;
        }

        if (c==')') {
            parenIn--;
            if (parenIn < 0) {
                WarnError (_String("Too many closing ')' near '") & input.Cut (MAX(0,index-32),index) & "'.");
                input = empty;
                return empty;
            }
            lastChar = 0;
            continue;
        }

        if (c=='[') {
            bracketIn++;
            lastChar = 0;
            continue;
        }

        if (c==']') {
            bracketIn--;
            lastChar = 0;
            continue;
        }


        if (c=='{') {
            if (matrixScope) {
                matrixScope++;
            } else if (lastChar == '=') { // a matrix def
                matrixScope++;
            } else {
                scopeIn++;
                if (index>=2) {
                    long t = input.FirstNonSpaceIndex (0, index-1, -1);
                    if (t>=1) {
                        if (input.getChar(t)=='o' && input.getChar(t-1)=='d') {
                            isDoWhileLoop << scopeIn-1;
                            //printf ("%d\n%s\n\n", isDoWhileLoop, input.Cut (t,-1).sData);
                        }
                    }
                }
            }
            lastChar = 0;
            continue;
        }

        if (c=='}') {
            if (matrixScope) {
                matrixScope--;
            } else {
                scopeIn--;
                if (!parenIn && !bracketIn) {
                    if (scopeIn >=0 && isDoWhileLoop.lLength && isDoWhileLoop.lData[isDoWhileLoop.lLength-1] == scopeIn) {
                        isDoWhileLoop.Delete (isDoWhileLoop.lLength-1);
                    } else if (scopeIn == 0) {
                        break;
                    }
                }

            }
            lastChar = 0;
            continue;
        }

        lastChar = c;
    }

    result.Finalize();
    _String::storageIncrement = saveSI;

    if (scopeIn||isString||isComment == 1||parenIn||matrixScope) {
        if (result!='}') {
            WarnError (_String("Expression appears to be incomplete/syntax error. Scope: ") &scopeIn & ", paretheses depth: "
                       & parenIn & ", matrix scope: " & matrixScope & '.' & matrixScope & '.' & (isString?"In a literal. ":empty) &
                       (isComment == 1? "In a comment ":empty) & '\n' & input);
            input = empty;
            return empty;
        } else {
            result = empty;
        }
    }

    lastChar=0;
    while (result.getChar(lastChar)=='{') {
        lastChar++;
    }

    if (lastChar) {
        long index2 = result.sLength-1;

        while (result[index2]=='}') {
            index2--;
        }

        if (result.sLength-index2-1 <lastChar) {
            ReportWarning ((_String)("Expression appears to be incomplete/syntax error and will be ignored:")&input);
            result.DuplicateErasing (&empty);
        } else {
            result.Trim(lastChar,result.sLength-1-lastChar);
        }
    }

    if (index<input.Length()-1) {
        input.Trim (index+1,-1, useSoftTrim);
    } else if (useSoftTrim) {
        input.sLength = 0;
    } else {
        input.DuplicateErasing (&empty);
    }


    return result;
}
//____________________________________________________________________________________

long _ElementaryCommand::ExtractConditions (_String& source, long startwith, _List& receptacle, char delimeter, bool includeEmptyConditions)
{
    long    parenLevel      = 1,
            lastsemi        = startwith,
            index,
            quote         = 0,
            curlyLevel        = 0;

    for (index = startwith; index<source.sLength; index++) {
        char c = source.sData[index];
        if (quote==0) {
            if (c=='(') {
                parenLevel++;
                continue;
            }
            if (c=='{') {
                curlyLevel++;
                continue;
            }
            if (c=='}') {
                curlyLevel--;
                continue;
            }
            if (c==')') {
                parenLevel--;
                if (!parenLevel) {
                    break;
                }
                continue;
            }
        }
        if (c=='"') {
            if (index == startwith || source.sData[index-1] != '\\') {
                quote+=quote?-1:1;
            }
            continue;
        }
        if (c==delimeter) {
            if (parenLevel>1 || quote || curlyLevel) {
                continue;
            }

            _String *term = (_String*)checkPointer(new _String (source,lastsemi,index-1));
            receptacle.AppendNewInstance (term);
            lastsemi = index+1;
            continue;
        }
    }

    if (includeEmptyConditions || lastsemi <= index-1) {
        receptacle.AppendNewInstance (new _String(source,lastsemi,index-1));
    }
    return index+1;
}

//____________________________________________________________________________________


bool       _ElementaryCommand::MakeGeneralizedLoop  (_String*p1, _String*p2, _String*p3 , bool fb, _String& source, _ExecutionList&target)
{

    // extract the for enclosure
    long  beginning = target.lLength,
          forreturn = target.lLength,
          index;

    int   success = 1;
    bool  hasIncrement = false;

    _SimpleList bc;

    while (1) {

        if (p1 && p1->Length()) { // initialization stage
            forreturn++;
            success *= target.BuildList (*p1, nil, true); // add init step
        }

        // append condition now

        if (!success) {
            break;
        }

        if (fb)
            if (p2 && p2->Length()) { // condition stage
                _ElementaryCommand condition (*p2);
                target&&(&condition);
            }

        if (source.getChar(0)=='{') {
            source.Trim(1,-1);
        }

        if ((success *= target.BuildList (source, &bc)) == 0) { // construct the main body
            break;
        }

        if (p3 && p3->Length()) { // increment stage
            success *= target.BuildList (*p3, nil,true); // add increment step
            hasIncrement = true;
        }

        if (!success) {
            break;
        }

        if (fb) {
            _ElementaryCommand loopback;
            success*=loopback.MakeJumpCommand (nil,forreturn,0,target);
            target&&(&loopback);
            if (p2 && p2->Length()) {
                success*=((_ElementaryCommand*)(target(forreturn)))->MakeJumpCommand (p2, forreturn+1, target.lLength,target);
            }
        } else {
            if (p2) {
                _ElementaryCommand* loopback = new _ElementaryCommand ;
                checkPointer (loopback);
                success*=loopback->MakeJumpCommand (p2,forreturn,target.lLength+1,target);
                target.AppendNewInstance(loopback);
            }
        }
        break;

    }

    if (!success) { // clean up
        for (index = beginning; index<target.lLength; index++) {
            target.Delete (beginning);
        }
        return false;
    } else {
        // write out the breaks and continues
        for (index = 0; index<bc.lLength; index++) {
            long loc = bc(index);
            if (loc>0) { // break
                ((_ElementaryCommand*)(target(loc)))->MakeJumpCommand (nil, target.lLength, 0,target);
            } else { // continue
                ((_ElementaryCommand*)(target(-loc)))->MakeJumpCommand (nil, target.lLength-(hasIncrement?2:1), 0,target);
            }
        }
    }

    return true;
}

//____________________________________________________________________________________


bool       _ElementaryCommand::BuildFor (_String&source, _ExecutionList&target,  _List & pieces)

// the for loop becomes this:
// initialize
// if (condition) then
// else go to end+2
// end
// increment
// goto if(condition)

{
    return MakeGeneralizedLoop ((_String*)pieces(0),(_String*)pieces(1),(_String*)pieces(2),true,source,target);
}

//____________________________________________________________________________________

bool    _ElementaryCommand::BuildWhile          (_String&source, _ExecutionList&target,  _List &pieces)
{
    return MakeGeneralizedLoop (nil,(_String*)pieces(0),nil,true,source,target);
}

//____________________________________________________________________________________

bool    _ElementaryCommand::BuildIfThenElse (_String&source, _ExecutionList&target, _SimpleList* bc)
{
    _List   pieces;
    long    upto = ExtractConditions (source,3,pieces),
            beginning = target.lLength;
    target.lastif << target.lLength;
    int     success = 1,
            intIfs = target.lastif.lLength;


    {
        if (pieces.lLength!=1) {
            WarnError ("'if' header makes no sense");
        }

        source.Trim (upto,-1);
        target.AppendNewInstance (new _ElementaryCommand);

        _String nextCommand (FindNextCommand(source));
        success *= target.BuildList (nextCommand, bc, true);

    }

    if (!success) { // clean up
        for (unsigned long index = beginning; index<target.lLength; index++) {
            target.Delete (beginning);
        }
        return false;
    } else {
        _ElementaryCommand* ec=(_ElementaryCommand*)(target(beginning));
        ((_ElementaryCommand*)(target(beginning)))->MakeJumpCommand (((_String*)pieces(0)), beginning+1, (ec->simpleParameters.lLength<2)?target.lLength:ec->simpleParameters(1),target);
    }

    while (target.lastif.lLength>intIfs) {
        target.lastif.Delete(target.lastif.lLength-1);
    }

    return target.BuildList(source,bc,true);
}



//____________________________________________________________________________________
bool    _ElementaryCommand::BuildDoWhile            (_String&source, _ExecutionList&target)
{
    long upto = source.FindBackwards(_String('}'), 0, -1);
    if (upto >= 0) {
        _String clipped (source, upto+1, -1);
        if (clipped.beginswith (blWhile)) {
            source.Trim (blDo.sLength,upto);
            _List pieces;
            ExtractConditions (clipped,blWhile.sLength,pieces);
            if (pieces.lLength != 1) {
                WarnError ("Malformed while clause in a do-while loop");
                return false;
            }

            if (!MakeGeneralizedLoop (nil,(_String*)pieces(0),nil,false,source,target)) {
                return false;
            }

            return true;
        }
    }
    WarnError ("Could not find a matching 'while' in the definition of a do-while loop");

    return false;
}

//____________________________________________________________________________________

bool    _ElementaryCommand::ProcessInclude      (_String&source, _ExecutionList&target)
{

    _String         fileName (source,blInclude.sLength,source.sLength-2);
    fileName = ProcessLiteralArgument (&fileName,target.nameSpacePrefix);
    if (fileName.sLength == 0) {
        WarnError (_String("#include missing a meaningful filename. Check that there is a ';' at the end of the statement. Had ")& source.Cut(8,source.sLength-2));
        return false;
    }

    fileName.ProcessFileName(false,false,(Ptr)target.nameSpacePrefix);
    if (terminateExecution) {
        return false;
    }

    PushFilePath  (fileName);
    ReadBatchFile (fileName, target);
    PopFilePath   ();

    return true;
}

//____________________________________________________________________________________

_ElementaryCommand* makeNewCommand (long ccode)
{
    _ElementaryCommand * newC = new _ElementaryCommand (ccode);
    checkPointer        (newC);
    return               newC;
}

//____________________________________________________________________________________

void _ElementaryCommand::addAndClean (_ExecutionList&target,_List* parList, long parFrom)
{
    if (parList)
        for (long k = parFrom; k<parList->lLength; k++) {
            parameters && (*parList) (k);
        }
    target << this;
    DeleteObject (this);
}

//____________________________________________________________________________________


bool    _ElementaryCommand::ConstructDataSet (_String&source, _ExecutionList&target)
// DataSet    dataSetid = ReadDataFile ("..");
// or
// DataSet    dataSetid = SimulateDataSet (likeFunc);
// or
// DataSet    dataSetid = Concatenate (<purge>,list of DataSets);
// or
// DataSet    dataSetid = Combine (<purge>,list of DataSets);
// or
// DataSet    dataSetid = ReconstructAncestors (lf)
// or
// DataSet    dataSetid = SampleAncestors (lf)
// or
// DataSet    dataSetid = Simulate (tree, freqs, alphabet, <store internal nodes, root vector>)
// or
// DataSet    dataSetid = ReadFromString (string);


{
    // first we must segment out the data set name
    // then the ReadDataFile command
    // then the data set file name

    // look for the data set name first

    long    mark1 = source.FirstSpaceIndex(0,-1,1),
            mark2 = source.Find ('=', mark1, -1);

    _String dsID (source,mark1+1,mark2-1);

    if (mark1==-1 || mark2==-1 || dsID.Length()==0) {
        WarnErrorWhileParsing ("DataSet declaration missing a valid identifier", source);
        return false;
    }

    // now look for the opening paren

    mark1 = source.Find ('(',mark2,-1);

    _ElementaryCommand dsc;
    _String            oper (source,mark2+1,mark1-1);

    if (oper ==  _String("ReadDataFile") || oper == _String ("ReadFromString")) { // a switch statement if more than 1
        _List pieces;
        mark2 = ExtractConditions (source,mark1+1,pieces,',');
        if (pieces.lLength!=1) {
            WarnErrorWhileParsing ("DataSet declaration missing a valid filename", source);
            return false;
        }

        _ElementaryCommand * dsc = makeNewCommand (5);

        dsc->parameters&&(&dsID);
        dsc->parameters&&(pieces(0));

        if (oper == _String ("ReadFromString")) {
            dsc->simpleParameters << 1;
        }

        dsc->addAndClean (target);
        return true;
    } else if (oper.Equal(&blSimulateDataSet)) {
        _List pieces;
        mark2 = ExtractConditions (source,mark1+1,pieces,',');
        if ( pieces.lLength>4 || pieces.lLength==0 ) {
            WarnErrorWhileParsing (blSimulateDataSet & "expects 1-4 parameters: likelihood function ident (needed), a list of excluded states, a matrix to store random rates in, and a matrix to store the order of random rates in (last 3 - optional).",
                                   source);
            return false;
        }

        dsc.code = 12;
        dsc.parameters&&(&dsID);
        dsc.parameters&&(pieces(0));
        for (mark2 = 1; mark2 < pieces.lLength; mark2++) {
            dsc.parameters&&(pieces(mark2));
        }

        target&&(&dsc);
        return true;
    } else if ( oper ==  _String("Concatenate") || oper ==  _String("Combine")) {
        _List pieces;
        mark2 = ExtractConditions (source,mark1+1,pieces,',');
        if (pieces.lLength==0) {
            WarnErrorWhileParsing("DataSet merging operation missing a valid list of arguments.",source);
            return false;
        }


        dsc.code = 16;
        dsc.parameters&&(&dsID);

        long i=0;

        dsc.simpleParameters<<((oper==_String("Concatenate"))?1:2);

        _String purge ("purge");
        if (purge.Equal ((_String*)pieces(0))) {
            dsc.simpleParameters[0]*=-1;
            i++;
        }

        for (; i<pieces.lLength; i++) {
            dsc.parameters<<pieces (i);
        }

        if (dsc.parameters.lLength<=1) {
            WarnErrorWhileParsing("DataSet merging operation missing a valid list of arguments.",source);
            return false;
        }

        target&&(&dsc);
        return true;

    } else {
        if (oper ==  _String("ReconstructAncestors") || oper ==  _String("SampleAncestors")) {
            _List pieces;
            mark2 = ExtractConditions (source,mark1+1,pieces,',');
            if (pieces.lLength>3 || pieces.lLength==0) {
                WarnErrorWhileParsing("ReconstructAncestors and SampleAncestors expects 1-4 parameters: likelihood function ident (mandatory), an matrix expression to specify the list of partition(s) to reconstruct/sample from (optional), and, for ReconstructAncestors, an optional MARGINAL flag, plus an optional DOLEAVES flag.",
                                      source);
                return false;
            }

            dsc.code                    = (oper == _String("ReconstructAncestors"))?38:50;
            dsc.parameters              &&(&dsID);
            dsc.parameters              << pieces(0);
            for (long optP = 1; optP < pieces.lLength; optP++)
                if (((_String*)pieces(optP))->Equal(&marginalAncestors)) {
                    dsc.simpleParameters << -1;
                } else if (((_String*)pieces(optP))->Equal(&doLeavesAncestors)) {
                    dsc.simpleParameters << -2;
                } else {
                    dsc.parameters  << pieces(optP);
                }

            target&&(&dsc);
            return true;
        } else if (oper ==  _String("Simulate")) {
            _List pieces;
            mark2 = ExtractConditions (source,mark1+1,pieces,',');
            if ((pieces.lLength>7)||(pieces.lLength<4)) {
                WarnErrorWhileParsing ("Simulate expects 4-6 parameters: tree with attached models, equilibrium frequencies, character map, number of sites|root sequence, <save internal node sequences>, <file name for direct storage>",
                                       source);
                return false;
            }

            dsc.code = 52;
            dsc.parameters&&(&dsID);

            for (mark2 = 0; mark2 < pieces.lLength; mark2++) {
                dsc.parameters&&(pieces(mark2));
            }

            target&&(&dsc);
            return true;
        } else {
            WarnErrorWhileParsing ("Expected DataSet ident = ReadDataFile(filename); or DataSet ident = SimulateDataSet (LikelihoodFunction); or DataSet ident = Combine (list of DataSets); or DataSet ident = Concatenate (list of DataSets); or DataSet ident = ReconstructAnscetors (likelihood function); or DataSet ident = SampleAnscetors (likelihood function) or DataSet	  dataSetid = ReadFromString (string);",
                                   source);
        }
    }

    return false;
}
//____________________________________________________________________________________

bool    _ElementaryCommand::ConstructCategory (_String&source, _ExecutionList&target)
// category <id> = (number of int, weights, method for representation, density, cumulative, left bound, right bound);
{

    long    mark1 = source.FirstSpaceIndex  (0,-1,1),
            mark2 = source.Find             ('=', mark1, -1);

    _String catID (source,mark1+1,mark2-1);

    if (mark1==-1 || mark2==-1 || catID.Length()==0 ) {
        _String errMsg ("Category variable declaration missing a valid identifier");
        WarnError (errMsg);
        return false;
    }

    // now look for the opening paren

    mark1 = source.Find ('(',mark2,-1);

    if (mark1!=-1) {
        mark2 = source.FindBackwards(')',mark1+1,-1);
        if (mark2!=-1) {
            source = source.Cut (mark1+1,mark2-1);
            _List args;
            mark2 = ExtractConditions (source,0,args,',');
            if (args.lLength>=7) {
                _ElementaryCommand * cv = new _ElementaryCommand (20);
                checkPointer (cv);
                cv->parameters&&(&catID);
                cv->addAndClean(target,&args,0);
                return true;
            }
        }
    }
    _String errMsg ("Expected: category <id> = (number of intervals, weights, method for representation, density, cumulative, left bound, right bound,<optional mean cumulative function>,<optional hidden markov matrix>);");
    WarnError (errMsg);
    return false;
}

//____________________________________________________________________________________

bool    _ElementaryCommand::ConstructStateCounter (_String&source, _ExecutionList&target)
// UseMatrix (matrixIdent)
{
    _List args;
    ExtractConditions (source,blStateCounter.sLength,args,',');
    if (args.lLength!=2) {
        WarnError ("Expected: StateCounter(likefuncID, callback function ID)");
        return false;
    }
    _ElementaryCommand * sc = new _ElementaryCommand(47);
    sc->addAndClean (target,&args,0);
    return true;
}


//____________________________________________________________________________________
bool    _ElementaryCommand::ConstructChoiceList(_String&source, _ExecutionList&target)
{
    _List args;
    ExtractConditions (source,blChoiceList.sLength,args,',');
    if (args.lLength<5) {
        WarnError  ("ChoiceList needs at least 5 arguments");
        return false;
    }
    _ElementaryCommand *cv = new _ElementaryCommand (32);

    cv->parameters<<args(0);
    ((_String*)args.lData[1])->StripQuotes();
    cv->parameters<<args(1);
    cv->parameters<<args(2);
    cv->parameters<<args(3);

    if  (args.lLength>5) {
        _List   choices;
        for (long k = 4; k<args.lLength-1; k+=2) {
            ((_String*)args.lData[k])->StripQuotes();
            ((_String*)args.lData[k+1])->StripQuotes();
            _List thisChoice;
            thisChoice<< args(k);
            thisChoice<< args(k+1);
            choices&& & thisChoice;
        }
        cv->parameters && & choices;
        cv->simpleParameters<<0;
    } else {
        cv->parameters<< args(4);
        cv->simpleParameters<<1;
    }

    cv->addAndClean(target,nil,0);
    return true;
}
//____________________________________________________________________________________

bool    _ElementaryCommand::ConstructReplicateConstraint (_String&source, _ExecutionList&target)
// ReplicateConstraint ("constraint to be replicated in terms of this1,...,thisn and"
//                       list of n variables to put in place of this1, this2, ... thisn);
// this1 .. etc are all expected to be either trees of nodes of trees with wildcards.
{
    _List args;
    ExtractConditions (source,20,args,',');
    if (args.lLength<2) {
        _String errMsg ("Expected: ReplicateConstraint (\"constraint to be replicated in terms of this1,...,thisn and wildcard *\", list of n variables to put in place of this1, this2, ... thisn);");
        acknError (errMsg);
        return false;
    }
    /*_String *theConstraint = (_String*)args(0), thisString;
    long k = 0;
    theConstraint->StripQuotes();
    do
    {
        k++;
        thisString  = _String("this")&_String(k);
    }
    while (theConstraint->Find(thisString)!=-1);

    if (args.lLength!=k)
    {
        _String errMsg ("Replicate constraint could not match the number of 'this' arguments with actual variables");
        acknError (errMsg);
        return false;
    }*/

    _ElementaryCommand cv;
    cv.code = 26;
    for (long k=0; k<args.lLength; k++) {
        cv.parameters << args(k);
    }
    target&& &cv;
    return true;
}

//____________________________________________________________________________________

bool    _ElementaryCommand::ConstructTree (_String&source, _ExecutionList&target)
// Tree   treeid = (...) or Topology = (...);
{
    long    mark1 = source.FirstSpaceIndex(0,-1,1), mark2, mark3;
    mark2 = source.Find ('=', mark1, -1);
    mark3 = mark2;

    if ((mark1==-1)||(mark2==-1)||(mark1+1>mark2-1)) {
        _String errMsg ("Tree declaration missing a valid identifier");
        acknError (errMsg);
        return false;
    }

    _String dsID = source.Cut (mark1+1,mark2-1);
    // now look for the opening paren

    mark1 = source.Find ('(',mark2,-1);
    mark2 = source.FindBackwards(')',mark1,-1);
    if ((mark1==-1)||(mark2==-1)||(mark2<mark1)) {
        if (source.Find(getDString)==-1) {
            mark1 = mark3+1;
            mark2 = source.Find (';',mark3,-1)-1;
        } else {
            source = getDString;
            mark1 = 0;
            mark2 = -1;
        }
    }

    _ElementaryCommand * dsc = new _ElementaryCommand(source.startswith(blTree)?7:54);
    checkPointer     (dsc);
    dsc->parameters&&(&dsID);
    dsc->parameters.AppendNewInstance(new _String(source,mark1,mark2));
    dsc->addAndClean(target,nil,0);
    return true;
}



//____________________________________________________________________________________

bool    _ElementaryCommand::ConstructDataSetFilter (_String&source, _ExecutionList&target)
// DataSetFilter      dataSetFilterid = CreateFilter (datasetid;unit;vertical partition; horizontal partition; alphabet exclusions);

{
    // first we must segment out the data set name

    long    mark1 = source.FirstSpaceIndex(0,-1,1),
            mark2 = source.Find ('=', mark1, -1);

    _String dsID    (source,mark1+1,mark2-1),
            command;

    if ( mark1==-1 || mark2==-1 || dsID.Length()==0) {
        _String errMsg ("DataSetFilter declaration missing a valid identifier");
        acknError (errMsg);
        return false;
    }

    // now look for the opening paren

    mark1 = source.Find ('(',mark2,-1);
    command = source.Cut (mark2+1,mark1-1);

    _ElementaryCommand *dsf;

    _List pieces;

    if (command == _String("CreateFilter")) {
        dsf = new _ElementaryCommand(6);
    } else if (command == _String("Permute")) {
        dsf = new _ElementaryCommand(27);
    } else if (command == _String("Bootstrap")) {
        dsf = new _ElementaryCommand(28);
    } else {
        _String errMsg ("Expected: DataSetFilter	  dataSetFilterid = CreateFilter (datasetid,unit,vertical partition,horizontal partition,alphabet exclusions); or Permute/Bootstrap (dataset/filter,<atom>,<column partition>)");
        acknError (errMsg);
        return false;
    }


    mark2 = ExtractConditions (source,mark1+1,pieces,',');
    if (!(pieces.lLength>=2 || (pieces.lLength == 1 && dsf->code == 6))) {
        _String errMsg ("Parameter(s) missing in DataSetFilter definition.");
        acknError (errMsg);
        return false;
    }

    dsf->parameters&&(&dsID);
    dsf->addAndClean (target,&pieces);
    return true;
}

//____________________________________________________________________________________

bool    _ElementaryCommand::ConstructModel (_String&source, _ExecutionList&target)

// Model ID = (inst transition matrix ident, equilibrium frequencies ident, <multiply by frequencies>);
// if the third parameter is explicitFormMExp, then inst transition matrix ident is expected to be an explicit matrix exponential
// EXPRESSION

{
    // first we must segment out the data set name

    long    mark1 = source.FirstSpaceIndex(0,-1,1),
            mark2 = source.Find ('=', mark1, -1);

    _String modelID (source,mark1+1,mark2-1);

    if (mark1==-1 || mark2==-1 || !modelID.IsValidIdentifier()) {
        _String errMsg ("Model declaration missing a valid identifier.");
        acknError (errMsg);
        return false;
    }

    // now look for the opening paren
    mark1 = source.Find ('(',mark2,-1);
    _List pieces;
    mark2 = ExtractConditions (source,mark1+1,pieces,',');

    if (pieces.lLength<2) {
        _String errMsg ("Parameter(s) missing in Model definition. Must have a matrix and a compatible eqiulibrium frequencies vector.");
        acknError (errMsg);
        return false;
    } else {
        if (pieces.lLength>3) {
            _String errMsg ("Too many parameters (3 max) in Model definition");
            acknError (errMsg);
            return false;
        }
    }

    _ElementaryCommand * model = new _ElementaryCommand(31);
    checkPointer (model);
    model->parameters&&(&modelID);
    model->addAndClean (target,&pieces,0);
    return true;

}


//____________________________________________________________________________________

bool    _ElementaryCommand::ConstructFprintf (_String&source, _ExecutionList&target)

{

    _ElementaryCommand  *fpr = (_ElementaryCommand*)checkPointer(new _ElementaryCommand (8));

    long     lastStart = 8;
    bool     done      = false;

    _String  comma (",");

    while (!done) {
        long lastEnd = source.FindTerminator(lastStart,comma);
        if (lastEnd < 0) {
            lastEnd = source.sLength-2;
            done = true;
        }
        _String *thisArgument = new _String (source, lastStart, lastEnd-1);

        if (fpr->parameters.lLength && thisArgument->IsALiteralArgument(true)) {
            fpr->simpleParameters << fpr->parameters.lLength;
            _FString converted (*thisArgument, true);
            fpr->parameters << converted.theString;
            DeleteObject (thisArgument);
        } else {
            fpr->parameters.AppendNewInstance (thisArgument);
        }
        lastStart = lastEnd + 1;
    }

    fpr->addAndClean(target, nil, 0);
    return true;
}

//____________________________________________________________________________________

bool    _ElementaryCommand::ConstructFscanf (_String&source, _ExecutionList&target)
// syntax:
// fscanf (stdin or "file name" or PROMPT_FOR_FILE, "argument descriptor", list of arguments to be read);
// argument descriptor is a comma separated list of one of the three constants
// "number", "matrix", "tree"
// list of arguments to be read specifies which variables will receive the values

{
    if (!allowedFormats.lLength) {
        allowedFormats.AppendNewInstance (new _String ("Number"));
        allowedFormats.AppendNewInstance (new _String ("Matrix"));
        allowedFormats.AppendNewInstance (new _String ("Tree"));
        allowedFormats.AppendNewInstance (new _String ("String"));
        allowedFormats.AppendNewInstance (new _String ("NMatrix"));
        allowedFormats.AppendNewInstance (new _String ("Raw"));
        allowedFormats.AppendNewInstance (new _String ("Lines"));
    }

    _ElementaryCommand  *fscan = new _ElementaryCommand (source.startswith (blsscanf)?56:25);
    _List               arguments, argDesc;
    long                f,p, shifter = 0;


    ExtractConditions   (source,7,arguments,',');
    if (arguments.lLength<3) {
        WarnError (_String("Too few arguments in call to fscanf or sscanf"));
        DeleteObject (fscan);
        return false;
    }
    fscan->parameters<<arguments(0);


    if (((_String*)arguments(1))->Equal (&blScanfRewind)) {
        fscan->simpleParameters << -1;
        shifter = 1;
    }

    ((_String*)arguments(1+shifter))->StripQuotes();
    ExtractConditions   (*((_String*)arguments(1+shifter)),0,argDesc,',');

    for (f = 0; f<argDesc.lLength; f++) {
        p = allowedFormats.Find(argDesc(f));
        if (p==-1) {
            WarnError ( *((_String*)argDesc(f))&" is not a valid type descriptor for fscanf. Allowed ones are:"& _String((_String*)allowedFormats.toStr()));
            DeleteObject (fscan);
            return false;
        } else {
            fscan->simpleParameters<<p;
        }
    }

    if (arguments.lLength!=fscan->simpleParameters.lLength+2) {
        WarnError (_String("fscanf passed ")&_String((long)(fscan->simpleParameters.lLength-shifter))&" parameter type descriptors and "
                   &_String((long)(arguments.lLength-2-shifter))& " actual arguments");
        DeleteObject (fscan);
        return false;
    }

    for (f = 2+shifter; f<arguments.lLength; f++) {
        _String* thisArg = (_String*)arguments(f);
        if (thisArg->IsValidIdentifier()) {
            fscan->parameters<< thisArg;
        } else {
            WarnError (_String("fscanf passed an invalid variable identifier: ")&*thisArg);
            DeleteObject (fscan);
            return false;
        }
    }

    fscan->addAndClean(target,nil,0);
    return true;
}

//____________________________________________________________________________________

bool      _ElementaryCommand::MakeJumpCommand       (_String* source,   long branch1, long branch2, _ExecutionList& parentList)
{
    long oldFla = 0;
    code        = 4;

    if (simpleParameters.lLength==3) {
        if (source) {
            _Formula* f = (_Formula*)simpleParameters(2);
            delete (f);
        } else {
            oldFla = simpleParameters(2);
        }
    }

    if (branch1==-1) {
        if (simpleParameters.lLength == 0) {
            WarnError ("An if-then-else scoping error. Check opening and closing brackets and double else's.");
            return false;
        }
        branch1 = simpleParameters[0];
    }

    simpleParameters.Clear();
    simpleParameters<<branch1;
    simpleParameters<<branch2;
    if (source) {
        /*_Formula f;
        long status = Parse (&f, *source, parentList.nameSpacePrefix,nil);
        if (status==-1)
            simpleParameters<<long(f.makeDynamic());*/
        parameters && source;
    } else if (oldFla) {
        simpleParameters<<oldFla;
    }

    return true;
}


//____________________________________________________________________________________
bool    _ElementaryCommand::ConstructFindRoot (_String&source, _ExecutionList&target)
// syntax: FindRoot (receptacle, expression, variable, left bound, right bound)
// or: Integrate (receptacle, expression, variable, left bound, right bound
{
    _List pieces;
    long    mark1 = source.Find ('(');
    _String oper (source,0,mark1);
    source.Trim(ExtractConditions (source,mark1+1,pieces,','),-1);
    if (pieces.lLength!=5) {
        WarnError ("Expected: FindRoot|Integrate (receptacle, expression, variable, left bound, right bound).");
        return false;
    }
    _ElementaryCommand * fri = new _ElementaryCommand(oper.Equal (&blFindRoot)?43:48);
    fri->addAndClean(target,&pieces,0);
    return true;
}

//____________________________________________________________________________________
bool    _ElementaryCommand::ConstructMPISend (_String&source, _ExecutionList&target)
// syntax: MPISend (numeric node ID, string with HBL code <or> a LF ID, <input redirect target>);
{

    _List pieces;
    ExtractConditions (source, blMPISend.sLength ,pieces,',');
    if (pieces.lLength!=2 && pieces.lLength!=3) {
        WarnError ("Expected: MPISend (numeric node ID, string with HBL code <or> a LF ID).");
        return false;
    }
    _ElementaryCommand * mpiSend = makeNewCommand (44);
    mpiSend->addAndClean (target, &pieces, 0);
    return true;
}

//____________________________________________________________________________________
bool    _ElementaryCommand::ConstructMPIReceive (_String&source, _ExecutionList&target)
// syntax: MPIReceive (can receive from node, received from node, receptacle for the string result);
{
    _List pieces;
    ExtractConditions (source, blMPIReceive.sLength ,pieces,',');
    if (pieces.lLength !=3 ) {
        WarnError ("Expected: MPIReceive (can receive from node, received from node, receptacle for the string result).");
        return false;
    }

    _ElementaryCommand * mpiRecv= makeNewCommand (45);
    mpiRecv->addAndClean (target, &pieces, 0);
    return true;
}

//____________________________________________________________________________________
bool    _ElementaryCommand::ConstructCategoryMatrix (_String&source, _ExecutionList&target)
// syntax: ConstructCategoryMatrix (receptacle, likelihood function, [COMPLETE/SHORT, which partitions to include -- a matrix agrument] )
{
    _List pieces;
    ExtractConditions (source,blConstructCM.sLength,pieces,',');
    if (pieces.lLength<2) {
        WarnError ("Expected: ConstructCategoryMatrix (receptacle, likelihood function,COMPLETE/SHORT/WEIGHTS [optional; default is COMPLETE], [optional matrix argument with partitions to include; default is to include all]");
        return false;
    }

    _ElementaryCommand * constuctCatMatrix = makeNewCommand(21);
    constuctCatMatrix->addAndClean (target, &pieces, 0);
    return true;
}


//____________________________________________________________________________________
bool    _ElementaryCommand::ConstructExecuteCommands (_String&source, _ExecutionList&target)
{

    _List pieces;

    char  execAFile = source.startswith (blExecuteAFile)?1:(source.startswith(blLoadFunctionLibrary)?2:0);
    long  code = 39;

    switch (execAFile) {
    case 0:
        ExtractConditions (source,blExecuteCommands.sLength,pieces,',');
        break;

    case 1:
        ExtractConditions (source,blExecuteAFile.sLength,pieces,',');
        code = 62;
        break;

    case 2:
        ExtractConditions (source,blLoadFunctionLibrary.sLength,pieces,',');
        code = 66;
        break;
    }

    if (pieces.lLength < 1 || pieces.lLength > 3) {
        WarnError ("Expected: ExecuteCommands (identifier, <compiled|(input redirect<,string prefix>)>) or ExecuteAFile (path name, <compiled|(input redirect<,string prefix>)> or LoadFunctionLibrary (path name, <compiled|(input redirect<,string prefix>)>)");
        return false;
    }

    _ElementaryCommand * exc = (_ElementaryCommand *)checkPointer(new _ElementaryCommand (code));

    exc->parameters<<pieces(0);

    if (pathNames.lLength) {
        exc->parameters && pathNames (pathNames.lLength-1);
    } else {
        exc->parameters && & empty;
    }

    if (pieces.lLength >1) {
        if (*(_String*)pieces(1) == _String("compiled")) {
            exc->simpleParameters << 1;
        } else {
            exc->parameters << pieces(1);
            if (pieces.lLength > 2) {
                exc->parameters << pieces(2);
            }
        }
    }

    exc->addAndClean (target);
    return true;
}

//____________________________________________________________________________________
bool    _ElementaryCommand::ConstructOpenWindow (_String&source, _ExecutionList&target)
{

    _List pieces;
    ExtractConditions (source,blOpenWindow.sLength,pieces,',');
    if (pieces.lLength<2 || pieces.lLength>3) {
        WarnError ("Expected: OpenWindow (window type,parameter matrix,<position string>)");
        return false;
    }

    if (pieces.lLength==3) {
        ((_String*)pieces(2))->StripQuotes();
    }


    _ElementaryCommand * exc = new _ElementaryCommand (40);
    exc->addAndClean(target, &pieces, 0);
    return true;
}

//____________________________________________________________________________________
bool    _ElementaryCommand::ConstructSpawnLF (_String&source, _ExecutionList&target)
{
    _List pieces;
    ExtractConditions (source,blSpawnLF.sLength,pieces,',');
    if (pieces.lLength!=4) {
        _String errMsg ("Expected: SpawnLikelihoodFunction (likeFuncID, treeID, window ID, subset matrix)");
        acknError (errMsg);
        return false;
    }

    _ElementaryCommand * exc = new _ElementaryCommand (41);
    exc->addAndClean(target,&pieces,0);
    return true;
}


//____________________________________________________________________________________
bool    _ElementaryCommand::ConstructGetDataInfo (_String&source, _ExecutionList&target)
// syntax: GetDataInfo(matrixID, dataFilterID,<sequence ref, site ref | sequence 1 , sequence 2, DISTANCES>)
{
    _List pieces;
    ExtractConditions (source,blGetDataInfo.sLength,pieces,',');
    if (pieces.lLength<2 ||pieces.lLength>5) {
        WarnError ("Expected: syntax: GetDataInfo(matrix ID, dataFilterID,<sequence ref, site ref | sequence 1 , sequence 2, DISTANCES>)");
        return false;
    }
    _ElementaryCommand * gdi = new _ElementaryCommand(46);
    gdi->addAndClean(target,&pieces,0);
    return true;
}

//____________________________________________________________________________________
bool    _ElementaryCommand::ConstructOpenDataPanel (_String&source, _ExecutionList&target)
// syntax: OpenDataPanel(dataSetID,"species order","display settings","partition settings")
{
    _List pieces;
    ExtractConditions (source,blOpenDataPanel.sLength,pieces,',');
    if (pieces.lLength!=4 && pieces.lLength!=5) {
        ReportWarning("Expected: syntax: OpenDataPanel(dataSetID,\"species order\",\"display settings\",\"partition settings\"),[likefunc ID]");
        return false;
    }

    _ElementaryCommand * sp = new _ElementaryCommand (36);
    sp->addAndClean(target,&pieces,0);
    return true;
}

//____________________________________________________________________________________
bool    _ElementaryCommand::ConstructGetInformation (_String&source, _ExecutionList&target)
// syntax: GetInformation(object,receptacle)
{

    _List pieces;
    ExtractConditions (source,blGetInformation.sLength,pieces,',');

    if (pieces.lLength < 2) {
        _String errMsg ("Expected at least 2 arguments: GetInformation(object,receptacle,...);");
        WarnError (errMsg);
        return false;
    }

    else {
        _String *s0 = (_String*)pieces(0),
                 *s1 = (_String*)pieces(1);

        if (!(s0->IsValidIdentifier()&&((s1->sLength>2&&s1->getChar(0)=='"'&s1->getChar(s1->sLength-1)=='"') || s1->IsValidIdentifier()))) {
            WarnError (_String ("Both ") & *s0 & " and " & *s1 & " must be valid identifiers in call to GetInformation.");
            return     false;
        }
    }

    _ElementaryCommand * sp = makeNewCommand(37);
    sp->addAndClean (target, &pieces, 0);
    return true;
}

//____________________________________________________________________________________
bool    _ElementaryCommand::ConstructLF (_String&source, _ExecutionList&target)
// syntax: LikelihoodFunction id = (filter1, tree1, ..., filterN, treeN, optional compute template)
// or LikelihoodFunction3 id = (filter1, tree1, freq1, ... filterN, treeN, freqN, optional compute template)
{
    long    mark1 = source.FirstSpaceIndex(0,-1,1),
            mark2 = source.Find ('=', mark1, -1);

    if ( mark1==-1 || mark2==-1 || mark1+1>mark2-1 ) {
        _String errMsg ("Likelihood function declaration missing a valid identifier");
        acknError (errMsg);
        return false;
    }

    _String lfID (source,mark1+1,mark2-1);
    // now look for the opening paren

    _List pieces;
    mark1 = source.Find ('(',mark2,-1);
    mark2 = source.FindBackwards(')',mark1,-1);
    ExtractConditions (source,mark1+1,pieces,',');

    if ( mark1==-1 || mark2==-1 || mark2<mark1 ) {
        WarnError ("Expected: Likelihood Function ident = (tree1, datasetfilter1,...)");
        return false;
    }

    _ElementaryCommand*  dsc = (_ElementaryCommand*)checkPointer(new _ElementaryCommand (11));
    dsc->parameters&&(&lfID);

    if (source.startswith(blLF3)) {
        dsc->simpleParameters << 1;
    }

    dsc->addAndClean(target,&pieces,0);
    return true;
}



//____________________________________________________________________________________
bool    _ElementaryCommand::ConstructFunction (_String&source, _ExecutionList& chain)
// syntax: function <ident> (comma separated list of parameters) {body}
{
    if (isInFunction) {
        WarnError ("Nested function declarations are not allowed");
        return false;
    }

    isInFunction = true;

    bool    isFFunction = source.beginswith (blFFunction);

    long    mark1 = source.FirstNonSpaceIndex(isFFunction?blFFunction.sLength:blFunction.sLength,-1,1),
            mark2 = source.Find ('(', mark1, -1);


    if ( mark1==-1 || mark2==-1 || mark1+1>mark2-1) {
        WarnError      ("Function declaration missing a valid function identifier or parameter list.");
        isInFunction = false;
        return false;
    }

    _String*    funcID  = (_String*)checkPointer(new _String(source.Cut (mark1,mark2-1)));
    *funcID = chain.AddNameSpaceToID (*funcID);

    // now look for the opening paren

    if ((mark1=batchLanguageFunctionNames.Find(funcID))!=-1) {
        ReportWarning (_String("Overwritten previously defined function:'") & *funcID & '\'');
    }

    _List pieces;

    long upto = ExtractConditions (source,mark2+1,pieces,',',false);


    if (upto==source.sLength || source[upto]!='{' || source[source.sLength-1]!='}') {
        WarnError (_String("Function declaration is missing a valid function body."));
        isInFunction= false;
        return false;
    }

    for (long k = 0; k < pieces.lLength; k++) {
        pieces.Replace (k,new _String(chain.AddNameSpaceToID (*(_String*)pieces(k))),false);
    }

    _String          sfunctionBody (source, upto+1,source.Length()-2);
    _ExecutionList * functionBody = new _ExecutionList (sfunctionBody,chain.GetNameSpace(),true);

    //  take care of all the return statements
    while (returnlist.lLength) {
        ((_ElementaryCommand*)(*functionBody)(returnlist(0)))->simpleParameters<<functionBody->lLength;
        returnlist.Delete(0);
    }


    if (mark1>=0) {
        batchLanguageFunctions.Replace (mark1, functionBody, false);
        batchLanguageFunctionNames.Replace (mark1, funcID, false);
        batchLanguageFunctionParameterLists.Replace (mark1, &pieces, true);
        batchLanguageFunctionParameters.lData[mark1] = pieces.lLength;
        batchLanguageFunctionClassification.lData[mark1] = isFFunction? BL_FUNCTION_NORMAL_UPDATE :  BL_FUNCTION_ALWAYS_UPDATE;
    } else {
        batchLanguageFunctions.AppendNewInstance(functionBody);
        batchLanguageFunctionNames.AppendNewInstance(funcID);
        batchLanguageFunctionParameterLists &&(&pieces);
        batchLanguageFunctionParameters     <<pieces.lLength;
        batchLanguageFunctionClassification <<(isFFunction? BL_FUNCTION_NORMAL_UPDATE :  BL_FUNCTION_ALWAYS_UPDATE);
    }


    isInFunction = false;
    return true;
}

//____________________________________________________________________________________
bool    _ElementaryCommand::ConstructReturn (_String&source, _ExecutionList&target)
// syntax: return <statement>
{

    /*if (!isInFunction)
    {
        _ElementaryCommand exit;
        exit.code = 4;
        exit.simpleParameters<<-1;
        target&&(&exit);
        return true;
    }*/

    long    mark1 = source.FirstNonSpaceIndex(blReturn.sLength,-1,1);

    _ElementaryCommand ret;

    ret.code = 14;

    if (mark1!=-1) {
        _String cut_s;
        if (source.sData[source.sLength-1]==';') {
            cut_s = source.Cut (mark1,source.sLength-2);
        } else {
            cut_s = source.Cut (mark1,-1);
        }

        ret.parameters&&(&cut_s);
    }

    if (isInFunction) {
        returnlist<<target.lLength;
    } else {
        ret.simpleParameters << -1;
    }

    target&&(&ret);
    return true;
}


//____________________________________________________________________________________

void    ReadBatchFile (_String& fName, _ExecutionList& target)
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

    FILE            *f = doFileOpen (fName.getStr(), "rb");
    SetStatusLine   ("Parsing File");
    if (!f) {
        WarnError (_String("Could not read batch file '") & fName & "'.\nPath stack: " & (_String*)pathNames.toStr());
    } else {
        _String beavis (f);

        if (beavis.beginswith("#NEXUS",false)) {
            ReadDataSetFile (f,1,nil,&fName);
        } else {
            target.BuildList (beavis);
            target.sourceFile = fName;
        }
        fclose (f);
    }
}

//____________________________________________________________________________________

void        SetDataFilterParameters (_String& parName, _DataSetFilter* thedf, bool setOrKill)
{
    _String     varName (parName&".species");
    _Variable*  receptacleVar = nil;

    if (setOrKill) {
        setParameter (varName, thedf->NumberSpecies());
    } else {
        DeleteVariable (varName);
    }

    varName = parName&".sites";
    if (setOrKill) {
        setParameter (varName, thedf->GetFullLengthSpecies()/thedf->GetUnitLength());
    } else {
        DeleteVariable (varName);
    }

    varName = parName&".unique_sites";
    if (setOrKill) {
        setParameter (varName,thedf->NumberDistinctSites());
    } else {
        DeleteVariable (varName);
    }

    varName = parName&".site_freqs";
    _Parameter      sizeCutoff;
    if (setOrKill) {
        checkParameter  (defaultLargeFileCutoff,sizeCutoff, 100000.);

        if (thedf->theFrequencies.lLength < sizeCutoff) {
            receptacleVar = CheckReceptacle (&varName, empty, false);
            receptacleVar->SetValue (new _Matrix(thedf->theFrequencies),false);
        }
    } else {
        DeleteVariable (varName);
    }

    varName = parName&".site_map";
    if (setOrKill) {
        if (thedf->theOriginalOrder.lLength < sizeCutoff) {
            receptacleVar = CheckReceptacle (&varName, empty, false);
            receptacleVar->SetValue (new _Matrix(thedf->theOriginalOrder),false);
        }
    } else {
        DeleteVariable (varName);
    }


    varName = parName&".sequence_map";
    if (setOrKill) {
        if (thedf->theOriginalOrder.lLength < sizeCutoff) {
            receptacleVar = CheckReceptacle (&varName, empty, false);
            receptacleVar->SetValue (new _Matrix(thedf->theNodeMap),false);
       }
    } else {
        DeleteVariable (varName);
    }
}

//____________________________________________________________________________________
void    SerializeModel  (_String& rec, long theModel, _AVLList* alreadyDone, bool completeExport)
{
    bool        mByF = true,
                do2  = false;

    _Variable   * tV  = nil,
                  * tV2 = nil;

    _Formula    * theExp  = nil;
    _SimpleList   matrices;

    if (modelTypeList.lData[theModel]) {
        theExp = (_Formula*)modelMatrixIndices.lData[theModel];
        theExp->ScanFForType(matrices, MATRIX);

        for (long mi = 0; mi < matrices.countitems(); mi++) {
            if (alreadyDone && alreadyDone->Insert ((BaseRef)matrices.lData[mi]) < 0) {
                matrices.Delete(mi);
                mi--;
            }
        }
    } else {
        if (!alreadyDone || alreadyDone->Find ((BaseRef)modelMatrixIndices.lData[theModel]) < 0) {
            if (alreadyDone) {
                alreadyDone->Insert ((BaseRef)modelMatrixIndices.lData[theModel]);
            }
            matrices << modelMatrixIndices.lData[theModel];
        }
        tV = LocateVar(modelMatrixIndices.lData[theModel]);
    }

    long freqID = modelFrequenciesIndices.lData[theModel];

    if (freqID>=0) {
        tV2 = LocateVar(freqID);
    } else {
        mByF = false;
        tV2 = LocateVar(-freqID-1);
    }

    if (!alreadyDone || alreadyDone->Find ((BaseRef)tV2->GetAVariable()) < 0) {
        if (alreadyDone) {
            alreadyDone->Insert ((BaseRef)tV2->GetAVariable());
        }
        do2 = true;
    }

    if (completeExport && (matrices.lLength || do2 || theExp)) {
        _SimpleList    vl,
                       ind,
                       dep,
                       cat;

        _AVLList vlst (&vl);

        if (theExp) {
            theExp->ScanFForVariables(vlst, true, false, true);
        }

        for (long mi = 0; mi < matrices.lLength; mi++) {
            LocateVar(matrices.lData[mi])->ScanForVariables (vlst,true);
        }

        if (do2) {
            tV2->ScanForVariables (vlst,true);
        }
        vlst.ReorderList ();
        SplitVariablesIntoClasses (vl,ind,dep,cat);

        _String glVars (128L, true),
                locVars(128L, true);

        ExportIndVariables (glVars,locVars, &ind);
        ExportDepVariables (glVars,locVars, &dep);
        glVars.Finalize();
        locVars.Finalize();
        rec<<glVars;
        rec<<locVars;
        ExportCatVariables (rec,&cat);
    }

    if (matrices.lLength) {
        for (long k = 0; k < matrices.lLength; k++) {
            _Variable *tV = LocateVar (matrices.lData[k]);
            ((_Matrix*)   tV->GetValue())->Serialize (rec,*tV->GetName());
            rec << '\n';
        }
    }

    if (do2) {
        ((_Matrix*)   tV2->GetValue())->Serialize (rec,*tV2->GetName());
    }

    rec << "\nModel ";
    rec << *((_String*)modelNames (theModel));
    rec << "=(";
    if (theExp) {
        rec << '"';
        rec << _String((_String*)(theExp->toStr()));
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
