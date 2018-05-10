/*
 
 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
 Art FY Poon    (apoon@cfenet.ubc.ca)
 Steven Weaver (sweaver@temple.edu)
 
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
#include "bayesgraph.h"
#include "function_templates.h"
#include "avllistx.h"
#include "global_object_lists.h"
#include "time_difference.h"




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


#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

using namespace hyphy_global_objects;


//____________________________________________________________________________________
// global variables


_List
dataSetList,
dataSetNamesList,
likeFuncList,   // list of all datasets
likeFuncNamesList, // list of all dataset filters
pathNames,
theModelList,
allowedFormats,
batchLanguageFunctions,
batchLanguageFunctionNames,
_batchLanguageFunctionNamesIndexed,
batchLanguageFunctionParameterLists,
batchLanguageFunctionParameterTypes,
compiledFormulaeParameters,
modelNames,
executionStack,
standardLibraryPaths,
standardLibraryExtensions,
loadedLibraryPathsBackend;

_AVLListX batchLanguageFunctionNamesIndexed (&_batchLanguageFunctionNamesIndexed);

#ifdef __MAC__
_String volumeName;
#endif


// retrieval functions

_SimpleList
returnlist,
batchLanguageFunctionClassification,
modelMatrixIndices,
modelTypeList,
modelFrequenciesIndices,
listOfCompiledFormulae;

_String
markdownOutput                  ("MARKDOWN_OUTPUT"),
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
                                listLoadedLibraries             ("LIST_OF_LOADED_LIBRARIES"),
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
                                alwaysReloadLibraries           ("ALWAYS_RELOAD_FUNCTION_LIBRARIES"),
                                dialogPrompt,
                                baseDirectory,
                                lastModelUsed,
                                libDirectory,
                                scanfLastFilePath,
                                defFileNameValue;


//____________________________________________________________________________________


_String  blFor                  ("for("),               // moved
blWhile                    ("while("),         // moved
blFunction                 ("function "),      // moved
blFFunction                ("ffunction "),     // moved
blLFunction                ("lfunction "),     // moved
blNameSpace                ("namespace "),
blReturn                   ("return "),        // moved
blReturnPrefix             ("return"),
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
blBGM                      ("BayesianGraphicalModel "),
blSimulateDataSet          ("SimulateDataSet"),
blAssert                   ("assert(");

_Trie    _HY_HBL_KeywordsPreserveSpaces  ;

#ifdef      __HYPHYMPI__

_String     mpiNodeID                       ("MPI_NODE_ID"),
            mpiNodeCount                    ("MPI_NODE_COUNT"),
            mpiLastSentMsg                  ("MPI_LAST_SENT_MSG");

void        ReportMPIError                  (int, bool);

#endif

_hy_nested_check  isInFunction = _HY_NO_FUNCTION;

_Parameter  explicitFormMatrixExponential = 0.0,
            messageLogFlag                = 1.0;

long        scanfLastReadPosition         = 0;

extern      _String             MATRIX_AGREEMENT,
            ANAL_COMP_FLAG;


extern      _SimpleList         freeSlots;


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
    _Variable *   mpiMsgVar = CheckReceptacle (&mpiLastSentMsg, emptyString, false);
    mpiMsgVar->SetValue (sentVal, false);
    //setParameter (mpiLastSentMsg, &sentVal);

}

//____________________________________________________________________________________
_String*    MPIRecvString       (long senderT, long& senderID) {
    _String*    theMessage = nil;
    long        messageLength = 0,
                transferCount = 0;

    int         actualReceived = 0;
    bool        isError       = false;

    if  (senderT<0) {
        senderT = MPI_ANY_SOURCE;
    }

    MPI_Status  status;
  
    // nonagressive polling mode
  
    int message_received = 0;
    while (! message_received) {
      MPI_Iprobe (senderT, HYPHY_MPI_SIZE_TAG, MPI_COMM_WORLD, &message_received, MPI_STATUS_IGNORE);
      usleep (100);
    }

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

const _String GetStringFromFormula (_String* data,_VariableContainer* theP) {
  
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

_Parameter  ProcessNumericArgument (_String* data, _VariableContainer const* theP, _ExecutionList* currentProgram) {
    _String   errMsg;
    _Formula  nameForm (*data,theP, currentProgram?&errMsg:nil);
     
    if (!errMsg.sLength) {
        _PMathObj formRes = nameForm.Compute();
        numericalParameterSuccessFlag = true;
        if (formRes&& formRes->ObjectClass()==NUMBER) {
            return formRes->Value();
        } else {
            if (formRes&& formRes->ObjectClass()==STRING) {
                return _String((_String*)((_FString*)formRes)->toStr()).toNum();
            } else {
                errMsg = (_String("'") & *data & "' was expected to be a numerical argument.");
            }
        }
    }
    
    if (currentProgram) {
        currentProgram->ReportAnExecutionError (errMsg);
    } else {
        WarnError(errMsg);
    }
    
    numericalParameterSuccessFlag = false;
    return 0.0;
}

//____________________________________________________________________________________

_PMathObj   ProcessAnArgumentByType (_String const* expression, _VariableContainer const* theP, long objectType, _ExecutionList* currentProgram)
{
    _String   errMsg;

  
    _Formula  expressionProcessor (*expression, theP, currentProgram?&errMsg:nil);
  
    if (errMsg.sLength && currentProgram) {
        currentProgram->ReportAnExecutionError (errMsg);
    }
    else {
        _PMathObj expressionResult = expressionProcessor.Compute(0,theP);
        if (expressionResult && expressionResult->ObjectClass()==objectType) {
          expressionResult->AddAReference();
          return expressionResult;
        }
    }
    
    return nil;
}


//____________________________________________________________________________________

const _String ProcessLiteralArgument (_String const* data, _VariableContainer const* theP, _ExecutionList* currentProgram) {
  //NLToConsole(); BufferToConsole("ProcessLiteralArgument:"); StringToConsole(*data); NLToConsole();
   _PMathObj getString = ProcessAnArgumentByType (data, theP, STRING, currentProgram);
  
    if (getString) {
      _String result (*((_FString*)getString)->theString);
      DeleteObject(getString);
      return result;
    }
    
    return emptyString;
}

//____________________________________________________________________________________

_AssociativeList*   ProcessDictionaryArgument (_String* data, _VariableContainer* theP, _ExecutionList* currentProgram) {
  return (_AssociativeList* )ProcessAnArgumentByType (data, theP, ASSOCIATIVE_LIST, currentProgram);
}




//____________________________________________________________________________________
long    FindDataSetName (_String const&s)
{
    return dataSetNamesList.FindObject (&s);
}


//____________________________________________________________________________________
long    FindLikeFuncName (_String const&s, bool tryAsAString)
{
    long try1 = likeFuncNamesList.FindObject (&s);
    if (try1 < 0 && tryAsAString) {
        _String s2 (ProcessLiteralArgument(&s, nil));
        try1 = likeFuncNamesList.FindObject(&s2);
    }
    return try1;
}

//____________________________________________________________________________________
long    FindModelName (_String const &s)
{
    if (s.Equal (&useLastModel)) {
        return lastMatrixDeclared;
    }

    return modelNames.FindObject (&s);
}

//____________________________________________________________________________________
_LikelihoodFunction*    FindLikeFuncByName (_String const&s)
{
    long i = FindLikeFuncName(s);
    if (i>=0) {
        return (_LikelihoodFunction*)likeFuncList (i);
    }
    return nil;
}

//____________________________________________________________________________________
long    FindSCFGName (_String const&s)
{
    return scfgNamesList.FindObject (&s);
}

//____________________________________________________________________________________
const _String&    GetBFFunctionNameByIndex  (long idx) {
  return *GetObjectNameByType (HY_BL_HBL_FUNCTION, idx, false);
}

//____________________________________________________________________________________
long   GetBFFunctionArgumentCount  (long idx) {
  return ((_List*)batchLanguageFunctionParameterLists.Element (idx))->countitems();
}

//____________________________________________________________________________________
_List&   GetBFFunctionArgumentList  (long idx) {
  return *(_List*)batchLanguageFunctionParameterLists.Element (idx);
}

//____________________________________________________________________________________
_SimpleList&   GetBFFunctionArgumentTypes  (long idx) {
  return *(_SimpleList*)batchLanguageFunctionParameterTypes.Element (idx);
}

//____________________________________________________________________________________
_ExecutionList&   GetBFFunctionBody  (long idx) {
  return *(_ExecutionList*)batchLanguageFunctions.Element (idx);
}

//____________________________________________________________________________________
_HY_BL_FUNCTION_TYPE   GetBFFunctionType  (long idx) {
  return (_HY_BL_FUNCTION_TYPE) batchLanguageFunctionClassification.Element (idx);
}

//____________________________________________________________________________________
_String const ExportBFFunction (long idx, bool recursive) {
  
  
  _String bf (8192UL, true);
  if (IsBFFunctionIndexValid(idx)) {
  
    _String hbf_name = GetBFFunctionNameByIndex (idx);
    _ExecutionList * body = &GetBFFunctionBody(idx);
    
    if (body->enclosingNamespace.sLength) {
      bf << "namespace " << body->enclosingNamespace << " {\n";
    }

    switch (GetBFFunctionType (idx)) {
      case BL_FUNCTION_SKIP_UPDATE:
        bf << blFFunction;
        break;
      case BL_FUNCTION_LOCAL:
        bf << blLFunction;
        break;
      default:
        bf << blFunction;
    }
    
    
    bf << hbf_name;
    bf << '(';
    
    long argument_count = GetBFFunctionArgumentCount (idx);
    _List * argument_list = &GetBFFunctionArgumentList (idx);
    for (long argument_id = 0; argument_id < argument_count; argument_id ++) {
      if (argument_id) {
        bf << ',';
      }
      
      
      bf << body->TrimNameSpaceFromID(*(_String*)argument_list->Element (argument_id));
      if (GetBFFunctionArgumentTypes (idx).GetElement(argument_id) == BL_FUNCTION_ARGUMENT_REFERENCE) {
        bf << '&';
      }
    }
    bf << ") {\n";
    bf << body->sourceText;
    bf << "\n}";
 
    if (body->enclosingNamespace.sLength) {
      bf << "\n}";
    }

    if (recursive) {
      _List      hbl_functions;
      _AVLListX other_functions (&hbl_functions);
      
      other_functions.Insert (new _String (hbf_name), HY_BL_HBL_FUNCTION, false, false);
      
      body->BuildListOfDependancies (other_functions, true);
      
      for (long i = 0; i < hbl_functions.lLength; i++) {
        _String * a_name = (_String*)hbl_functions (i);
        if (! a_name -> Equal( &hbf_name)) {
          bf << "\n/*----- Called function '";
          bf << *a_name;
          bf << "' ------*/\n";
          bf << ExportBFFunction (FindBFFunctionName(*a_name), false);
          bf << "\n\n";
        }
      }
    }
  }
  
  bf.Finalize();
  return bf;
  
}

//____________________________________________________________________________________
void ClearBFFunctionLists (long start_here) {
  if (start_here > 0L && start_here < batchLanguageFunctionNames.countitems()) {
    
    _SimpleList delete_me (batchLanguageFunctionNames.countitems()-start_here, start_here, 1L);
    
    for (long k = 0; k < delete_me.countitems(); k++) {
      batchLanguageFunctionNamesIndexed.Delete (batchLanguageFunctionNames.GetItem (delete_me.Get (k)));
    }
    
    batchLanguageFunctionNames.DeleteList           (delete_me);
    batchLanguageFunctions.DeleteList               (delete_me);
    batchLanguageFunctionClassification.DeleteList  (delete_me);
    batchLanguageFunctionParameterLists.DeleteList  (delete_me);
    batchLanguageFunctionParameterTypes.DeleteList  (delete_me);
    
    
    
  } /*else {
    batchLanguageFunctionNames.Clear();
    batchLanguageFunctions.Clear();
    batchLanguageFunctionClassification.Clear();
    batchLanguageFunctionParameterLists.Clear();
    batchLanguageFunctionParameterTypes.Clear();
  }*/
}

//____________________________________________________________________________________
bool IsBFFunctionIndexValid (long index) {
  if (index >= 0L && index < batchLanguageFunctionNames.countitems()) {
    return batchLanguageFunctions.GetItem(index) != nil;
  }
  return false;
}

//____________________________________________________________________________________
long GetBFFunctionCount (void) {
  return batchLanguageFunctions.countitems();
}

//____________________________________________________________________________________
long    FindBFFunctionName (_String const&s, _VariableContainer const* theP) {
    if (theP) {
        _String prefix = *(theP->GetName());

        //ReportWarning (_String ("Looking for ") & s.Enquote() & " in " & prefix.Enquote());

        while (1) {
            _String test_id = prefix & '.' & s;
            long idx = batchLanguageFunctionNamesIndexed.Find (&test_id);
            if (idx >= 0) {
                return batchLanguageFunctionNamesIndexed.GetXtra(idx);
                //s = test_id;
                //return idx;
            }
            long cut_at = prefix.FindBackwards ('.', 0, -1);
            if (cut_at > 0) {
              prefix.Trim (0, cut_at - 1);
            } else {
              break;
            }
        };
    }

  
    //ReportWarning (_String ("Looking for ") & s.Enquote() & " in global context");
    return batchLanguageFunctionNamesIndexed.FindAndGetXtra(&s,-1);
}


//____________________________________________________________________________________
long    FindBgmName (_String const&s) {
    return bgmNamesList.FindObject (&s);
}



//__________________________________________________________
long  AddDataSetToList (_String& theName,_DataSet* theDS) {
    theName = GenerateUniqueObjectIDByType(theName, HY_BL_DATASET);
    long k = dataSetNamesList.FindObject (&emptyString);
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

void KillLFRecord (long lfID, bool completeKill) {
    /* compile the list of variables which will no longer be referenced */

  
  
    if (lfID>=0) {
        //printf ("\n****\nKillLFRecord\n%s\n****\n", (char const*) * (_String*)likeFuncNamesList.GetItem (lfID));
        _LikelihoodFunction *me = (_LikelihoodFunction*)likeFuncList (lfID);

        if (completeKill) {
            _SimpleList         wastedVars,
                                otherVars,
                                myVars,
                                otherModels,
                                wastedModels;



            myVars  << me->GetIndependentVars();
            myVars  << me->GetDependentVars();

          
          
            for (unsigned long k=0UL; k<likeFuncList.lLength; k++) {
                  if (k!=lfID) {
                      if (((_String*)likeFuncNamesList(k))->sLength) {
                          _LikelihoodFunction *lf = (_LikelihoodFunction*)likeFuncList (k);
                          otherVars << lf->GetIndependentVars();
                          otherVars << lf->GetDependentVars();
                        
                          unsigned long component_count = lf->CountObjects(kLFCountPartitions);
                        
                        for (long tree_index = 0UL; tree_index < component_count; tree_index++) {
                          lf->GetIthTree(tree_index)->CompileListOfModels(otherModels);
                        }
                        
                      }
                  }
            }

            myVars.Sort ();
            otherVars.Sort();
            otherModels.Sort();
          
            wastedVars.Subtract(myVars, otherVars);

            for (unsigned long k=0UL; k<myVars.lLength; k++)
                if (otherVars.BinaryFind(myVars.lData[k])<0) {
                    wastedVars << myVars.lData[k];
                }

            myVars.Clear();
          
            unsigned long component_count = me->CountObjects(kLFCountPartitions);

            for (long tree_index = 0UL; tree_index < component_count; tree_index++) {
                _TheTree* thisTree = me->GetIthTree(tree_index);
                thisTree->CompileListOfModels (myVars);
                _TreeIterator ti (thisTree, _HY_TREE_TRAVERSAL_POSTORDER);
                while (_CalcNode* tNode = ti.Next()) {
                    tNode->SetValue (new _Constant (tNode->ComputeBranchLength()),false);
                }
                thisTree->RemoveModel();
            }

            for (unsigned long k=0UL; k<myVars.lLength; k++)
                if (otherModels.BinaryFind (myVars.lData[k])<0) {
                    KillModelRecord (myVars.lData[k]);
                }

            for (unsigned long k=0UL; k<wastedVars.lLength; k++) {
                //printf ("Deleting %ld (%s)\n", wastedVars.lData[k],  ->GetName()->getStr());
                _Variable * check_me = LocateVar(wastedVars.lData[k]);
                if (check_me) {
                  DeleteVariable (*check_me->GetName());
                }
            }

        }
        // unregister event listeners to release reference counts
        
        me->UnregisterListeners();
        


        if (lfID<likeFuncList.lLength-1) {
            DeleteObject(likeFuncList(lfID));
            likeFuncList.lData[lfID] = 0L;
            likeFuncNamesList.Replace(lfID,&emptyString,true);
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

void KillLFRecordFull (long lfID) {
    _LikelihoodFunction* lf = (_LikelihoodFunction*) likeFuncList (lfID);

    _SimpleList l;
    lf->GetGlobalVars (l);
  
    for (unsigned long k=0UL; k<l.lLength; k++) {
        DeleteVariable (*LocateVar(l.lData[k])->GetName());
    }

    l.Clear ();
  
    unsigned long partition_count = lf->CountObjects(kLFCountPartitions);

    for (unsigned long k=0UL; k<partition_count; k++) {
        _TheTree * thisTree = lf->GetIthTree(k);
        thisTree->CompileListOfModels (l);
        DeleteVariable (*thisTree->GetName());
    }

    for (unsigned long k=0UL; k<l.lLength; k++) {
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
        dataSetNamesList.Replace(dsID,&emptyString,true);
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
    }

    if (mdID<modelNames.lLength-1) {
        modelMatrixIndices.lData[mdID] = -1;
        modelTypeList.lData[mdID] = 0;
        modelFrequenciesIndices.lData[mdID] = -1;
        modelNames.Replace(mdID,&emptyString,true);
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
_ExecutionList::_ExecutionList () {
    Init();
} // doesn't do much

//____________________________________________________________________________________
_ExecutionList::_ExecutionList (_String& source, _String* namespaceID , bool copySource, bool* successFlag) {
    Init (namespaceID);
    
    if (copySource) {
        sourceText.Duplicate (&source);
    }

    bool result = BuildList (source, nil, false, true);
    if (successFlag) {
        *successFlag = result;
    }
}

//____________________________________________________________________________________
void _ExecutionList::Init (_String* namespaceID) {
    result              = nil;
    currentCommand      = 0;
    cli                 = nil;
    profileCounter      = nil;
    stdinRedirect       = nil;
    stdinRedirectAux    = nil;
    doProfile           = 0;
    nameSpacePrefix     = nil;
    
    if (currentExecutionList) {
        errorHandlingMode  = currentExecutionList->errorHandlingMode;
        errorState         = currentExecutionList->errorState;
    } else {
        errorHandlingMode = HY_BL_ERROR_HANDLING_DEFAULT;
        errorState = false;
    }

    if (namespaceID) {
        SetNameSpace (*namespaceID);
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
    _ExecutionList * Res = new _ExecutionList;

    memcpy ((char*)Res, (char*)this, sizeof (_ExecutionList));

    Res->nInstances         = 1;
    Res->Duplicate          (this);
    Res->cli                = nil;
    Res->profileCounter     = nil;
    Res->doProfile          = doProfile;
    Res->errorHandlingMode  = errorHandlingMode;
    Res->errorState         = errorState;

    if(result) {
        Res->result = (_PMathObj)result->makeDynamic();
    }

    return Res;
}

//____________________________________________________________________________________

void        _ExecutionList::Duplicate   (BaseRef source) {
    _List::Duplicate    (source);

    _ExecutionList* s = (_ExecutionList*)source;

    if (s->result) {
        s->result=(_PMathObj)result->makeDynamic();
    }

    errorHandlingMode  = s->errorHandlingMode;
    errorState         = s->errorState;
}


//____________________________________________________________________________________
void    _ExecutionList::ReportAnExecutionError (_String errMsg, bool doCurrentCommand, bool appendToExisting) {
    if (doCurrentCommand) {
        _ElementaryCommand *theCommand = FetchLastCommand();
        if (theCommand) {
            errMsg = errMsg & " in call to " & _HY_ValidHBLExpressions.RetrieveKeyByPayload(theCommand->GetCode());
        }
    }
    errorState = true;
    switch (errorHandlingMode) {
        case HY_BL_ERROR_HANDLING_SOFT:
            if (appendToExisting) {
              _FString * existing = (_FString*) FetchObjectFromVariableByType(&_hyLastExecutionError, STRING);
              if (existing) {
                errMsg = *existing->theString & '\n' & errMsg;
              }
            }
            setParameter(_hyLastExecutionError, new _FString (errMsg, false), nil, false);
            
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

_String  const     _ExecutionList::GetFileName     (void)  const {
    if (sourceFile.sLength) {
        return sourceFile;
    } else {
        _String const *top_path = PeekFilePath();
        if (top_path)
          return *top_path;
    }
    return emptyString;
}

//____________________________________________________________________________________

void _ExecutionList::BuildListOfDependancies   (_AVLListX & collection, bool recursive) {
  for (unsigned long step = 0UL; step < lLength; step++) {
    ((_ElementaryCommand*)GetItem(step))->BuildListOfDependancies (collection, recursive, *this);
  }
}

//____________________________________________________________________________________

_PMathObj       _ExecutionList::Execute     (_ExecutionList* parent) {

  //setParameter(_hyLastExecutionError, new _MathObject, nil, false);
  
  
    _ExecutionList*      stashCEL = currentExecutionList;
    callPoints << currentCommand;
    executionStack       << this;

    if (parent && stdinRedirect == nil) {
      stdinRedirect = parent->stdinRedirect;
      stdinRedirectAux = parent->stdinRedirectAux;
    } else {
      parent = nil;
    }
 
    _FString            cfp (PeekFilePath() ? *PeekFilePath () :emptyString),
                      * stashed = (_FString*)FetchObjectFromVariableByType (&pathToCurrentBF, STRING);

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

            TimeDifference timer;
            (((_ElementaryCommand**)lData)[currentCommand])->Execute(*this);
            timeDiff   = timer.TimeSinceStart();
          

          if (profileCounter) {
            // a call to _hyphy_profile_dump can set this to NULL
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
        setParameter        (pathToCurrentBF,stashed,nil, false);
    }

    executionStack.Delete (executionStack.lLength-1);
    if (result == nil) {
        result = new _MathObject();
    }

    if (parent) {
      stdinRedirect = nil;
      stdinRedirectAux = nil;
    }

    return result;
}

//____________________________________________________________________________________

long        _ExecutionList::ExecuteAndClean     (long g, _String* fName)        // run this execution list
{
    long    f = -1;
    Execute ();

    if (fName && !terminateExecution) {
        f = batchLanguageFunctionNamesIndexed.Find (fName);
        if (f >= 0) {
          f = batchLanguageFunctionNamesIndexed.GetXtra (f);
        }
    }

    terminateExecution      = false;
    skipWarningMessages     = false;

    ClearBFFunctionLists    (g);
  
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

    for (unsigned long k = 0; k<lLength && status; k++) {
        _ElementaryCommand * aStatement = (_ElementaryCommand*)(*this)(k);
        switch (aStatement->code) {
        case 0: {
            _String * formulaString = (_String*)aStatement->parameters(0);

            if (formulaString->sData[formulaString->sLength-1]!='}') {
                _Formula *f  = new _Formula,
                *f2 = new _Formula;

                checkPointer ((BaseRef)(f&&f2));

                _FormulaParsingContext fpc (nil, nameSpacePrefix);

                long          parseCode = Parse(f,*formulaString,fpc,f2);

                if (parseCode == HY_FORMULA_EXPRESSION || parseCode == HY_FORMULA_VARIABLE_VALUE_ASSIGNMENT || parseCode == HY_FORMULA_FORMULA_VALUE_ASSIGNMENT) {
                  
                    if (f->AmISimple(stackDepth,varList)) {
                        try {
                          if (parseCode == HY_FORMULA_FORMULA_VALUE_ASSIGNMENT) {
                            if (!f2->AmISimple(stackDepth, varList)) throw 0;
                            long assignment_length = f->NumberOperations();
                            if (assignment_length < 3) throw 0;
                            _Variable * mx = f->GetIthTerm(0)->RetrieveVar();
                            if (! mx) throw 0;
                            f->GetIthTerm (0)->SetAVariable(mx->GetAVariable());
                            _Operation * last = f->GetIthTerm(assignment_length-1);
                            if (! (last->TheCode() == HY_OP_CODE_MCOORD && last->GetNoTerms() == 2)) throw 0;
                            
                            
                            f2->GetList() << f->GetList();
                            f->Clear();
                            
                            _Formula *t = f2;
                            f2 = f;
                            f  = t;
          
                          }
                          
                        } catch (int e) {
                          status = false;
                          break;
                        }
                        aStatement->simpleParameters<<parseCode;
                        aStatement->simpleParameters<<(long)f;
                        aStatement->simpleParameters<<(long)f2;
                      
                      
                        aStatement->simpleParameters<<fpc.assignmentRefID();
                      

                        formulaeToConvert << (long)f;
                      

                        if (parseCode == HY_FORMULA_VARIABLE_VALUE_ASSIGNMENT) {
                            if (varList.Find (fpc.assignmentRefID()) < 0) {
                              varList << fpc.assignmentRefID();
                            }
                            parseCodes        << fpc.assignmentRefID();
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
                    _FormulaParsingContext fpc (nil, nameSpacePrefix);

                    long status = Parse (&f, *(_String*)aStatement->parameters(0), fpc, nil);

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

        for (unsigned long fi = 0; fi < formulaeToConvert.lLength; fi++) {
            ((_Formula*)formulaeToConvert(fi))->ConvertToSimple (varList);
        }

        for (unsigned long vi = 0; vi < varList.lLength; vi++) {
            avlList.Insert ((BaseRef)varList.lData[vi], vi);
        }

        for (unsigned long ri = 0; ri<parseCodes.lLength; ri++) {
            if (parseCodes.lData[ri] < 0) {
                cli->storeResults << -1;
            } else {
                cli->storeResults << avlList.GetXtra (avlList.Find ((BaseRef) parseCodes.lData[ri]));
            }
            //printf ("\n%ld\n",  cli->storeResults.lData[ri]);
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
    currentCommand = 0L;
    _SimpleList to_delete_aux;
    _AVLList to_delete (&to_delete_aux);
    while (currentCommand<lLength) {
        _ElementaryCommand* thisCommand = ((_ElementaryCommand**)lData)[currentCommand];
        if (thisCommand->DecompileFormulae()) {
          to_delete.Insert(thisCommand);
        }
        currentCommand++;
    }
  
    if (to_delete.countitems()) {
      _SimpleList batch_delete;
      for (unsigned long i = 0; i < listOfCompiledFormulae.lLength; i++) {
        if (to_delete.Find ((BaseRef)listOfCompiledFormulae.Element(i)) >= 0) {
          batch_delete << i;
        }
      }
      listOfCompiledFormulae.DeleteList(batch_delete);
      compiledFormulaeParameters.DeleteList(batch_delete);
    }
}
//____________________________________________________________________________________

BaseRef  _ExecutionList::toStr (unsigned long)
{
    _String *result = new _String (1,true),
    step ("\n\nStep "),
    dot (".");
  
    _ExecutionList* stash = currentExecutionList;
  
    currentExecutionList = this;

    for (unsigned long i=0UL; i<countitems(); i++) {
        (*result) << &step;
        (*result)<< _String((long)i);
        (*result)<< '.';
        result->AppendNewInstance ((_String*)GetItem(i)->toStr());
    }
    result->Finalize();
  
    currentExecutionList = stash;
    return result;
}

//____________________________________________________________________________________

void     _ExecutionList::ResetNameSpace (void)
{
    DeleteObject (nameSpacePrefix);
    nameSpacePrefix = nil;
}

//____________________________________________________________________________________

void     _ExecutionList::SetNameSpace (_String nID) {
    ResetNameSpace ();
    nameSpacePrefix = new _VariableContainer(nID);
}

//____________________________________________________________________________________

_String*     _ExecutionList::GetNameSpace () {
    if (nameSpacePrefix) {
        return nameSpacePrefix->GetName();
    }
    return nil;
}

//____________________________________________________________________________________

_String  _ExecutionList::AddNameSpaceToID (_String& theID, _String * extra) {
    _String name_space;
            
    if (extra && extra->sLength) {
        if (nameSpacePrefix) {
            name_space = (*nameSpacePrefix->GetName())&'.'& *extra;
        } else {
            name_space = *extra;
        }
    } else {
        if (nameSpacePrefix) {
            name_space = (*nameSpacePrefix->GetName());        
        }
    }
            
    return AppendContainerName (theID, &name_space);
}

//____________________________________________________________________________________

_String  _ExecutionList::TrimNameSpaceFromID (_String& theID)
{
    if (nameSpacePrefix) {
        if (theID.startswith(*nameSpacePrefix->GetName())) {
            return theID.Cut(nameSpacePrefix->GetName()->sLength+1,-1);
        }
    }
    return theID;
}



/* 
 
  holds all the expressions that require that spaces between them and the next expressions be 
  maintained, like
 
  return expr
  DataSet expr = 
  DateSetFilter expr =
 
  if (expr) is an identifier, then the spaces will be maintained, otherwise they will 
  be squished, causing incorrect behavior (like DataSet(expr) will gets parsed as a formula)
 
  initialized in _HBL_Init_Const_Arrays

*/
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
        long prefixTreeCode = _HY_ValidHBLExpressions.FindKey (currentLine, &triePath, true);
        
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
                _ElementaryCommand::BuildFor (currentLine, *this, pieces);
                handled = true;
                break;
            case HY_HBL_COMMAND_WHILE:
                _ElementaryCommand::BuildWhile (currentLine, *this, pieces);
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
                currentLine = emptyString;
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
            case HY_HBL_COMMAND_FPRINTF:
                _ElementaryCommand::ExtractValidateAddHBLCommand (currentLine, prefixTreeCode, pieces, commandExtraInfo, *this);
                handled = true;
                break;
                
        }
        
        if (handled) {
              DeleteObject (pieces);
              if (currentLine.Length() > 1UL) {
                WarnError (currentLine.Enquote() & " contained syntax errors, possibly a missing semicolon. " & s.Enquote());
              }
        }
        
        // 20111212: this horrendous switch statement should be replaced with a 
        // prefix tree lookup 

        if (!handled) {
            if (currentLine.startswith (blFunction)||currentLine.startswith (blFFunction)||currentLine.startswith (blLFunction) || currentLine.startswith (blNameSpace)) { // function declaration
                _ElementaryCommand::ConstructFunction (currentLine, *this);
            } else if (currentLine.startswith_noident (blReturnPrefix)) { // function return statement
                                                                          //StringToConsole(currentLine); NLToConsole();
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
    s.DuplicateErasing (&emptyString);
    return empty_is_success || countitems();
}

//____________________________________________________________________________________

_ElementaryCommand::_ElementaryCommand (void) {
    code = -1;
}

//____________________________________________________________________________________

_ElementaryCommand::_ElementaryCommand (long ccode) {
    code = ccode;
}

//____________________________________________________________________________________

_ElementaryCommand::_ElementaryCommand (_String& s) {
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

const _String _hblCommandAccessor (_ExecutionList* theList, long index) {
    if (theList) {
        if (index >= 0) {
            if (index < theList->lLength) {
                _ElementaryCommand * aCommand = (_ElementaryCommand*)theList->GetItem (index);
                return _String ((_String*)aCommand->toStr());
            } else {
              return _String("<END EXECUTION>");
            }
        }
    }
    return _String ("command index ") & index;
}

//____________________________________________________________________________________

BaseRef   _ElementaryCommand::toStr      (unsigned long)
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
        if (simpleParameters.countitems()==3 || parameters.countitems() == 1) {
            converted = (_String*)parameters.GetItem(0)->toStr();
            result = result& "under condition '"& *converted&"'\n\tto\n\t\t"&
                        _hblCommandAccessor (currentExecutionList,simpleParameters(0))&
                        "\n\telse\n\t\t"&
                        _hblCommandAccessor (currentExecutionList,simpleParameters(1));
        } else {
            result = result&"to Step "& simpleParameters(0);
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

    case HY_HBL_COMMAND_FPRINTF: // print stuff to file (or stdout)

        converted = (_String*)parameters(0)->toStr();
        result = _String("fprintf(")& (*converted);
        DeleteObject (converted);

        converted = nil;
        for (k = 1; k<parameters.countitems(); k++) {
            result = result&*","& *((_String*)parameters(k));
        }
        result = result & ")";


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
        if (labs(simpleParameters(0))==1) {
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
    case 61: {
        converted = (_String*)parameters(0)->toStr();
        result = blSCFG & *converted & "=(";
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
        
      case HY_HBL_COMMAND_NESTED_LIST: {
        converted = (_String*)parameters(0)->toStr();
        result = _String("Call a nested list (via namespace):\n ") & *converted;
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
  
    _String * errMsg = nil;
  
    try {

      if (chain.cli) {
          _Parameter result = ((_Formula*)simpleParameters.lData[1])->ComputeSimple (chain.cli->stack, chain.cli->values);
          long sti = chain.cli->storeResults.lData[chain.currentCommand-1];
          if (sti>=0) {
            //printf ("%ld, %g\n", sti, result);
              chain.cli->values[sti].value = result;
          }
          return;
      }

      if (!simpleParameters.lLength) { // not compiled yet
          _Formula f,
                   f2;

          _String* theFla     = (_String*)parameters(0);
        
          _FormulaParsingContext fpc (nil, chain.nameSpacePrefix);

          long     parseCode = Parse(&f,(*theFla),fpc,&f2);

          if (parseCode != HY_FORMULA_FAILED ) {
              if (fpc.isVolatile() == false) { // not a matrix constant
                  simpleParameters    <<parseCode;
                  simpleParameters    <<long (f.makeDynamic());
                  simpleParameters    <<long (f2.makeDynamic());
                  simpleParameters    <<fpc.assignmentRefID   ();
                  simpleParameters    <<fpc.assignmentRefType ();
                  appendCompiledFormulae (&f, &f2);
                
              } else {
                  ExecuteFormula(&f,&f2,parseCode,fpc.assignmentRefID(),chain.nameSpacePrefix,fpc.assignmentRefType());
                  if (terminateExecution) {
                    errMsg = new _String ("Error computing the compiled statement: ");
                    throw 0;
                  }
                  return;
              }
          } else {
            errMsg = new _String ("Error compiling the statement: ");
            throw 0;
          }
      }

      ExecuteFormula ((_Formula*)simpleParameters.lData[1],(_Formula*)simpleParameters.lData[2],simpleParameters.lData[0],simpleParameters.lData[3], chain.nameSpacePrefix, simpleParameters.lData[4]);
      
      if (terminateExecution) {
        errMsg = new _String ("Error computing the interpreted statement: ");
        throw 0;
      }
      
    } catch (int e) {
      if (errMsg) {
        WarnError (_String(errMsg) & *this);
      }
    }
}


//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase4 (_ExecutionList& chain)
{
    chain.currentCommand++;
    /*if (simpleParameters.lLength==2) {

    }*/
  
    _Formula * expression = nil;
    _String  * errMsg = nil;

    try {
      if (simpleParameters.lLength==3 || parameters.lLength) {
        
        
          if ( parameters.lLength && simpleParameters.lLength < 3) {
              expression = new _Formula;
              //printf ("Namespace: %x\nCode: %s\n", chain.nameSpacePrefix, ((_String*)parameters(0))->sData);

              _FormulaParsingContext fpc (nil,  chain.nameSpacePrefix);
              long status = Parse (expression, *(_String*)parameters(0), fpc, nil);

              //printf ("Print formula: %s\n", _String((_String*)expression->toStr()).sData);

              if (status== HY_FORMULA_EXPRESSION) {
                if (fpc.isVolatile() == false) {
                    simpleParameters << (long)expression;
                    appendCompiledFormulae (expression);
                    expression = nil;
                }
              } else {
                  errMsg = new _String (" is not a valid conditional expression");
                  throw (0);
              }
          }

          if (chain.cli) {
              if ( ((_Formula*)simpleParameters(2))->ComputeSimple(chain.cli->stack, chain.cli->values)==0.0) {
                  chain.currentCommand = simpleParameters.lData[1];
                  return;
              }
          } else {
              _PMathObj result;
              if (expression) {
                  //printf ("\n*** Interpreted condition\n");
                result = expression->Compute();
              } else {
                  //printf ("\n*** Compiled condition\n");
                result = ((_Formula*)simpleParameters(2))->Compute();
              }

              // printf ("\n*** %s\n", ((_String*)result->toStr())->sData);

            if (terminateExecution && !result) {
                  subNumericValues = 2;
                  _String       *s = (_String*)((_Formula*)simpleParameters(2))->toStr();
                  subNumericValues = 0;
                  errMsg  = new _String(_String("Failed while evaluating: ") & _String((_String*)((_Formula*)simpleParameters(2))->toStr()) & " which expanded to  " & s);
                  throw (1);
               }

              bool conditionFalse = false;

              switch (result->ObjectClass()) {
                case NUMBER:
                    conditionFalse = result->Value()==0.0;
                    break;
                case STRING:
                    conditionFalse = ((_FString*)result)->IsEmpty();
                    break;
                case HY_UNDEFINED:
                    conditionFalse = true;
                    break;
                default:
                    errMsg = new _String(_String(" did not evaluate to a number, a string, or a null (") &  (_String*)result->toStr() & ")");
                    throw (0);
              }
            
              if (expression) {
                delete expression;
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
    catch (int e) {
      if (expression) {
        delete expression;
      }
      if (errMsg) {
        if (e == 0) {
          WarnError (_String ("'") & *(_String*)parameters(0) & "'" & errMsg);
        } else {
          WarnError    (errMsg);
        }
        // note that errMsg will be deleted by _String (*_String) constructors
      }
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
            fName.ProcessFileName(false,true,(Ptr)chain.nameSpacePrefix);
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
                     WarnError (_String ("Could not find source dataset file ") & ((_String*)parameters(1))->Enquote('"')
                                & " (resolved to '" & fName & "')\nPath stack:\n\t" & GetPathStack ("\n\t"));
                    return;
                }
            }
            ds = ReadDataSetFile (df,0,nil,nil,chain.nameSpacePrefix?chain.nameSpacePrefix->GetName():nil);
            fclose (df);
        }
    }


    // 20110802: need to check that this data set is not emptyString

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
            DeleteObject (matrixOfStrings);
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

        if(GetDataFilter (AppendContainerName(*dataset,chain.nameSpacePrefix)) != nil) {
            _TheTree*   thisTree = (_TheTree*)FetchObjectFromVariableByType(&AppendContainerName(*tree,chain.nameSpacePrefix),TREE);
            if (thisTree) {
                _TreeIterator ti (thisTree, _HY_TREE_TRAVERSAL_POSTORDER);
                if (!freq) { // no explicit frequency parameter; grab one from the tree
                    long        theFreqID       = -1,
                                theModelID     = -1,
                                finalFreqID        = -1;
                    bool        done = false;

                    while (1) {
                        _CalcNode *thisNode = ti.Next();
                      
                        if ((theModelID     = thisNode->GetModelIndex()) == HY_NO_MODEL) { // this node has no model
                            done = false;
                            break;
                        }
                        theFreqID   = modelFrequenciesIndices.lData[theModelID];
                      
                        while((thisNode = ti.Next()) && !ti.IsAtRoot()) {
                            theModelID      = thisNode->GetModelIndex();
                            if (theModelID == HY_NO_MODEL) { // no model
                                done = false;
                                break;
                            }
                            if (modelFrequenciesIndices.lData[theModelID]!=theFreqID) {
                                done = true;
                                break;
                            }
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
                errMsg = (((_String)("LF: Not a valid tree ID: `")& *tree & "`"));
            }

        } else {
            errMsg = (((_String)("LF: Not a valid dataset filter `")& *dataset & "`"));
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
            likeFuncObjectID = likeFuncNamesList.FindObject(&emptyString);
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

        

      
        if (parameters.lLength>3) {
            // a matrix to store simulated category values
            _String  matrixName (chain.AddNameSpaceToID(*(_String*)parameters(3)));
          
            if (!(catValVar = CheckReceptacle(&matrixName,blSimulateDataSet,true))) {
                return;
            } else {
                checkPointer (catValues = new _Matrix (1,1,false,true));
            }
        }

        if (parameters.lLength>4) {
            // a matrix to store simulated category values
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
        DeleteObject   (partitionList);
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


void      _ElementaryCommand::ExecuteCase39 (_ExecutionList& chain) {
    chain.currentCommand++;

    _String *commands,
            theCommand,
            *namespc = nil;
  
    _AVLListXL * inArg    = nil;
    _List      * inArgAux = nil;

    bool    pop_path = false;
  
    try {
      if (code == 39) {
          commands = ProcessCommandArgument((_String*)parameters(0));
      } else {
          _String filePath = GetStringFromFormula((_String*)parameters(0),chain.nameSpacePrefix),
                  originalPath = filePath;
        

          FILE * commandSource = nil;
          
          _Parameter reload = 0.;
          checkParameter(alwaysReloadLibraries, reload, 0.);

          if (code == 66) {
              bool hasExtension    = filePath.FindBackwards (".",0,-1) > 0;

              for (unsigned long p = 0; !commandSource && p < standardLibraryPaths.lLength; p++) {
                  for (unsigned long e = 0; !commandSource && e < standardLibraryExtensions.lLength; e++) {
                      _String tryPath = *((_String*)standardLibraryPaths(p)) & filePath & *((_String*)standardLibraryExtensions(e));

                      // printf ("%s\n", tryPath.sData);

                      tryPath.ProcessFileName (false, false, (Ptr)chain.nameSpacePrefix);
                    
                      if (loadedLibraryPaths.Find(&tryPath) >= 0 && parameters.lLength == 2 && reload < 0.5) {
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
   
              if (code == 66 && loadedLibraryPaths.Find(&filePath) >= 0 && parameters.lLength == 2 && reload < 0.5) {
                  ReportWarning (_String("Already loaded '") & originalPath & "' from " & filePath);
                  return;
              }
              
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
          pop_path = true;
          PushFilePath (filePath);
      }

      if (!commands) {
          throw (1);
      }

      if (code == 39) {
          theCommand = ProcessLiteralArgument (commands,chain.nameSpacePrefix);
      } else {
          theCommand = commands;
      }

      if (theCommand.sLength == 0) {
          WarnError (_String("Invalid string argument '") & *commands & "' in call to ExecuteCommands/ExecuteAFile.");
          throw (1);
      }

      if (code == 39 && ((_String*)parameters(1))->sLength) {
        pop_path = true;
        PushFilePath (*(_String*)parameters(1), false, false);
      }


      if (parameters.lLength >= 3)
          // stdin redirect (and/or name space prefix)
      {

         _PMathObj inAVL = ProcessDictionaryArgument ((_String*)parameters(2),chain.nameSpacePrefix);

          if (!inAVL) {
              if (parameters.lLength == 3) {
                  WarnError (_String("Not a valid associative array index passed as input redirect argument to ExecuteCommands/ExecuteAFile: )") & *(_String*)parameters(2));
                  throw (1);
              }
          } else {
              _AssociativeList * stdinRedirect = (_AssociativeList*)inAVL;

              inArgAux = new _List;
              inArg    = new _AVLListXL (inArgAux);

              _List        *stdKeys = stdinRedirect->GetKeys();

              for (long kid = 0; kid < stdKeys->lLength; kid++) {
                  _String  * aKey         = (_String*)(*stdKeys) (kid);
                  if (aKey) {
                      _FString * aString      = (_FString*)stdinRedirect->GetByKey (*aKey, STRING);
                      if (!aString) {
                          WarnError    (_String("All entries in the associative array used as input redirect argument to ExecuteCommands/ExecuteAFile must be strings. The following key was not: ") & *aKey);
                          DeleteObject (inAVL);
                          throw (1);
                      }
                      inArg -> Insert (aKey->makeDynamic(),(long)new _String (*aString->theString),false);
                  }
              }
          }

          DeleteObject (inAVL);

          if (parameters.lLength > 3) {
              _String nameSpaceID = ProcessLiteralArgument((_String*)parameters(3),chain.nameSpacePrefix);
              if (nameSpaceID.sLength > 0UL) {
                if (!nameSpaceID.IsValidIdentifier(true)) {
                    WarnError (_String("Invalid namespace ID in call to ExecuteCommands/ExecuteAFile: ") & *(_String*)parameters(3));
                    throw (1);
                }
                namespc = new _String (nameSpaceID);
              }
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
                  ReportWarning (_String ("Successfully compiled an execution list.\n") & _String ((_String*)exc.toStr()) );
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
    }  catch (int e) {
      
    }
  

    DeleteObject            (inArg);
    DeleteObject            (inArgAux);
    DeleteObject            (namespc);

    if (pop_path) {
      PopFilePath();
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
                        columnHeaders && & emptyString;

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

void      _ElementaryCommand::ExecuteCase25 (_ExecutionList& chain, bool issscanf) {
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
  
    _Variable*  iseof          = CheckReceptacle (&hasEndBeenReached,emptyString,false);
  

    if (currentParameter==_String("stdin")) { //
        if (chain.stdinRedirect) {
            data = chain.FetchFromStdinRedirect ();
            // echo the input if there is no fprintf redirect in effect
            _FString * redirect = (_FString*)FetchObjectFromVariableByType (&blFprintfRedirect, STRING);
            if (! (redirect && redirect->theString->sLength)) {
              StringToConsole (*data); NLToConsole();
            }
          
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
                scanfLastFilePath = emptyString;
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
            currentParameter = GetStringFromFormula (&currentParameter,chain.nameSpacePrefix);
            currentParameter.ProcessFileName(false,false,(Ptr)chain.nameSpacePrefix);
            if (terminateExecution) {
                return;
            }
            inputBuffer = doFileOpen (currentParameter.getStr(), "rb");
            if (!inputBuffer) {
                WarnError         (currentParameter& " could not be opened for reading by fscanf. Path stack:\n\t" & GetPathStack("\n\t"));
                return;
            }

            if (iseof->Compute()->Value()>0) {
                scanfLastFilePath = emptyString;
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
                v = CheckReceptacle(&namespacedParameter,emptyString,false)->GetAVariable();
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
                                lines && & emptyString;
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
                        lines && & emptyString;
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
             f2=-1L,
             matrixDim,
             f3,
             multFreqs = 1;



    if (parameters.lLength>3) {
        parameterName = (_String*)parameters.lData[3];
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
    }
  
  
    if (doExpressionBased) {
        _String matrixExpression (ProcessLiteralArgument((_String*)parameters.lData[1],chain.nameSpacePrefix)),
                defErrMsg = _String ("The expression for the explicit matrix exponential passed to Model must be a valid matrix-valued HyPhy formula that is not an assignment") & ':' & matrixExpression;
        // try to parse the expression, confirm that it is a square  matrix,
        // and that it is a valid transition matrix
        isExpressionBased = new _Formula;
        _FormulaParsingContext fpc (nil, chain.nameSpacePrefix);
        matrixExpression =  _ElementaryCommand::FindNextCommand (matrixExpression);
        long parseCode = Parse(isExpressionBased,matrixExpression,fpc, nil);
        if (parseCode != HY_FORMULA_EXPRESSION || isExpressionBased->ObjectClass()!= MATRIX ) {
            WarnError (defErrMsg & " parse code = " & parseCode & " " & (parseCode == HY_FORMULA_EXPRESSION ? (_String(", object type code ") & _String((long) isExpressionBased->ObjectClass())) : emptyString ));
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
  

  
    // so far so good
    matrixDim = checkMatrix->GetHDim();
    if ( matrixDim!=checkMatrix->GetVDim() || matrixDim<2 ) {
      WarnError (*parameterName & " must be a square matrix of dimension>=2 in the call to Model = ...");
      return;
    }

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
  
   if (checkMatrix->GetVDim()==1UL) {
        if (checkMatrix->GetHDim()!=matrixDim) {
            WarnError (*parameterName & " must be a column vector of the same dimension as the model matrix in the call to Model = ...");
            return;
        }
    } else if (checkMatrix->GetHDim()==1UL) {
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

    long existingIndex = modelNames.FindObject(&arg0);

    if (existingIndex == -1) { // name not found
        lastMatrixDeclared = modelNames.FindObject (&emptyString);

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


    _String     saveTheArg,
                dialog_title (ProcessLiteralArgument((_String*)parameters(1),chain.nameSpacePrefix));
    _SimpleList sel,
                exclusions;

    long         do_markdown;
    checkParameter (markdownOutput,do_markdown,0L);
  
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
            f = FindDataFilter (nmspName);
            if (f>=0) {
                parameters.Delete(4);
              
                _DataSetFilter const *theFilter = GetDataFilter (f);
                _DataSet *linked_set = theFilter->GetData();
              
                for (unsigned long species_index  = 0; species_index < theFilter->NumberSpecies(); species_index ++) {
                    if (exclusions.BinaryFind(species_index) >= 0) {
                        continue;
                    }
                  
                  choices < &((*new _List)
                              << linked_set->GetSequenceName(species_index)
                              < new _String (_String ("Taxon ") & (species_index + 1) & '(' & *linked_set->GetSequenceName(species_index) & ')'));
                }
              
                validChoices = true;
                parameters&& & choices;
            } else {
                f = FindDataSetName (nmspName);
                if (f>=0) {
                    parameters.Delete(4);
                    _DataSet *linked_set = (_DataSet*)dataSetList (f);
                    for (unsigned long species_index  = 0; species_index < linked_set->NoOfSpecies(); species_index ++) {
                        if (exclusions.BinaryFind(species_index) >= 0) {
                            continue;
                        }

                        choices < &((*new _List)
                                  << linked_set->GetSequenceName(species_index)
                                  < new _String (_String ("Taxon ") & (species_index + 1) & '(' & *linked_set->GetSequenceName(species_index) & ')'));

                     }
                  
                    validChoices = true;
                    parameters&& & choices;
                } else {
                    if (saveTheArg==lastModelParameterList) {
                        f = lastMatrixDeclared;
                    } else {
                        f = modelNames.FindObject(&nmspName);
                    }

                    if (f>=0) {
                        parameters.Delete(4);
                      
                        _Variable *theSet = LocateVar (modelMatrixIndices.lData[f]);
                        _SimpleList modelParms;
                      
                        choices << &((*new _List) < "All Parameters" < "All local model parameters are constrained");
                      
                        _AVLList modelParmsA (&modelParms);
                        theSet->ScanForVariables(modelParmsA,false);
                        modelParmsA.ReorderList();
                        for (f = 0; f<modelParms.lLength; f++) {
                            if (exclusions.BinaryFind(f)>=0) {
                                continue;
                            }

                            choices << &((*new _List)
                                         << LocateVar(modelParms.lData[f])->GetName()
                                         < new _String (_String ("Constrain parameter ") & *LocateVar(modelParms.lData[f])->GetName()));
                          
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
                    for (choice = 0; choice<theChoices->lLength; choice++) {
                        if (buffer.Equal ((_String*)(*(_List*)(*theChoices)(choice))(0))) {
                            break;
                        }
                    }
                    if (choice == theChoices->lLength) {
                        WarnError (_String("Not a valid option: '") & buffer & "' passed to Choice List '" & dialog_title & "' using redirected stdin input");
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
                                for (choice = 0; choice<theChoices->lLength; choice++) {
                                    if (buffer.Equal((_String*)(*(_List*)(*theChoices)(choice))(0))) {
                                        break;
                                    }
                                }
                                if (choice<theChoices->lLength && sel.Find(choice)==-1) {
                                    sel<<choice;
                                } else {
                                    WarnError (_String("Not a valid (or duplicate) option: '") & buffer & "' passed to ChoiceList (with multiple selections) '" & dialog_title & "' using redirected stdin input");
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
              
              
              
                if (do_markdown) {
                  printf ("\n\n####%s\n", dialog_title.getStr());
                } else {
                
                  printf ("\n\n\t\t\t+");

                  for (f = 1; f<dialog_title.sLength+1; f++) {
                      printf ("-");
                  }

                  printf ("+\n\t\t\t|%s|\n\t\t\t+",(const char*)dialog_title);

                  for (f = 1; f<dialog_title.sLength+1; f++) {
                      printf ("-");
                  }

                  printf ("+\n\n");
                }


                long  loopits = 1;

                if (fixedLength == 1) {
                    while (choice == -1) {
                        for (choice = 0; choice<theChoices->lLength; choice++) {
                          if (do_markdown) {
                            printf ("\n%ld. [**%s**] %s", choice+1, ((_String*)(*(_List*)(*theChoices)(choice))(0))->getStr(),((_String*)(*(_List*)(*theChoices)(choice))(1))->getStr());
                          } else {
                            printf ("\n\t(%ld):[%s] %s",choice+1,((_String*)(*(_List*)(*theChoices)(choice))(0))->getStr(),((_String*)(*(_List*)(*theChoices)(choice))(1))->getStr());
                          }
                        }
                      
                        printf ("\n\n%sPlease choose an option (or press q to cancel selection):", do_markdown ? ">" : "");
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
                            printf ("\n\n%sPlease choose option %ld of %ld (or press q to cancel selection):",do_markdown ? ">" : "", sel.lLength+1,fixedLength);
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
                            printf ("\n\n%sPlease choose option %ld, enter d to complete selection, enter q to cancel:",do_markdown ? ">" : "",sel.lLength+1);
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
            }

            _Variable* sStrV = CheckReceptacle(&selectionStrings,emptyString,false);

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

    long    f = dataSetNamesList.FindObject(&AppendContainerName(*currentArgument,chain.nameSpacePrefix)),
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

    _HYDataPanel*  newDP = new _HYDataPanel (emptyString,emptyString);
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
void      _ElementaryCommand::ExecuteCase37 (_ExecutionList& chain) {
  chain.currentCommand++;
  
  _String matrixName = chain.AddNameSpaceToID(*(_String*)parameters(0)),
         *objectName = (_String*)parameters(1);
  
  
  _Matrix *result = nil;
  
  // object is a non-emptyString string
  if (objectName->sLength > 2 && objectName->sData[0] == '"' && objectName->sData[objectName->sLength-1] == '"') {
    // regular expression
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
    _String objectNameID = AppendContainerName(*objectName,chain.nameSpacePrefix);
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
          result = new _Matrix (2,f,false,true);
          
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
          } else {
            if (theObject->ObjectClass() == TOPOLOGY || theObject->ObjectClass() == TREE) {
              
              _List* map = ((_TreeTopology*)theObject)->MapNodesToModels ();
              _AssociativeList* return_this = new _AssociativeList();
              
              for (unsigned long i = 0; i < map->lLength; i++) {
                _List * nodeInfo = (_List*) map->GetItem(i);
                return_this->MStore(*(_String*)nodeInfo->GetItem(0), *(_String*)nodeInfo->GetItem (1));
              }
              result = (_Matrix*) return_this;
              DeleteObject (map);
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
      f = likeFuncNamesList.FindObject (&objectNameID);
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
        if ((f = FindDataFilter(objectNameID))>=0)
          // return a vector of strings - each with actual characters of the corresponding sequence
        {
          _DataSetFilter const * daFilter = GetDataFilter (f);
          result = daFilter->GetFilterCharacters();
        } else {
          // it could be a model
          f = FindModelName (objectNameID);
          if (f>=0) {
            // for models, return the list of variables in the model
            _SimpleList modelParms;
            _AVLList    modelParmsA (&modelParms);
            
            
              if (IsModelOfExplicitForm (f)) {
                  ((_Formula*)modelMatrixIndices.lData[f])->ScanFForVariables(modelParmsA,false);
              } else {
                  LocateVar (modelMatrixIndices.lData[f])->ScanForVariables(modelParmsA,false);
            
              }
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
  
  CheckReceptacleAndStore (&matrixName, emptyString, true, result, false);
  
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


    long      nodeCount;
    checkParameter (mpiNodeCount,nodeCount,1L);

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
        theMessage->AppendNewInstance(ar->Serialize(0UL));
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
        WarnError (*arg2 & " is not a valid (or is an emptyString) string (LF ID) in call to MPISend.");
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

    long      nodeCount;
    checkParameter (mpiNodeCount,nodeCount,1L);

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

void      _ElementaryCommand::ExecuteCase46 (_ExecutionList& chain) {
  chain.currentCommand++;
  
  _String *arg1 = (_String*)parameters(1),
          *arg2 = (_String*)parameters(0),
          errMsg;
  
  const _String source_name = AppendContainerName(*arg1,chain.nameSpacePrefix);
  
  long            object_type = HY_BL_DATASET|HY_BL_DATASET_FILTER;
  BaseRefConst    source_object = _HYRetrieveBLObjectByName (source_name, object_type,nil,false);
  
    //_DataSetFilter const * dsf = GetDataFilter    (filter_name);
    //_DataSet       const * ds  =
  
  if (source_object == nil) {
    errMsg = source_name.Enquote('\'') & " is not a defined data set / filter ID ";
  } else {
    
    const _String receptacle_name = AppendContainerName(*arg2,chain.nameSpacePrefix);
    _Variable *      stVar = CheckReceptacle(&receptacle_name,"GetDataInfo");
    
    if (stVar) {
      
      _DataSetFilter const * filter_source  = object_type == HY_BL_DATASET_FILTER ? (_DataSetFilter const *)source_object : nil;
      _DataSet       const * dataset_source = filter_source ? nil : (_DataSet const *)source_object;
      
      switch (parameters.lLength) {
        case 2UL: { // get site->pattern map
          if (filter_source) {
            stVar->SetValue (new _Matrix (filter_source->duplicateMap),false);
          } else {
            stVar->SetValue (new _Matrix (dataset_source->DuplicateMap()),false);
          }
        }
          break;
          
        case 3UL: { // data parameters, or sequence string
          _String argument = ProcessLiteralArgument ((_String*)parameters(2),chain.nameSpacePrefix);
          if (argument == _String("CHARACTERS")) {
            _List characters;
            if (filter_source) {
              
              unsigned long character_count = filter_source->GetDimension(true),
              fd = filter_source->GetUnitLength();
              
              for (unsigned long idx = 0UL; idx < character_count; idx++) {
                characters < new _String (filter_source->ConvertCodeToLetters (filter_source->CorrectCode(idx), fd));
              }
            } else {
              _String alphabet_string = dataset_source->GetTT () ? dataset_source->GetTT ()->GetAlphabetString() : emptyString;
              for (unsigned long idx = 0UL; idx < alphabet_string.sLength; idx++) {
                characters < new _String (alphabet_string (idx));
              }
            }
            stVar->SetValue (new _Matrix (characters), false);
          } else if (argument == _String ("PARAMETERS")) { // argument == _String("CHARACTERS")
            if (filter_source) {
              _AssociativeList * parameterInfo = new _AssociativeList;
              
              (*parameterInfo) < (_associative_list_key_value){"ATOM_SIZE", new _Constant (filter_source->GetUnitLength())}
                               < (_associative_list_key_value){"EXCLUSIONS", new _FString  (filter_source->GetExclusions())}
                               < (_associative_list_key_value){"SITES_STRING", new _FString  ((_String*)filter_source->theOriginalOrder.ListToPartitionString())}
                               < (_associative_list_key_value){"SEQUENCES_STRING", new _FString  ((_String*)filter_source->theNodeMap.ListToPartitionString())};
              
              stVar->SetValue (parameterInfo,false);
              
            } else {
              errMsg = argument.Enquote('\'') & " is a supported argument for the dataset source";
            }
          } else if (argument == _String ("CONSENSUS")) { // argument == _String("PARAMETERS")
            if (filter_source) {
              stVar->SetValue (new _FString (new _String(filter_source->GenerateConsensusString())), false);
            } else {
              _DataSetFilter temp;
              _SimpleList l1, l2;
              temp.SetFilter (dataset_source, 1, l1, l2, false);
              stVar->SetValue (new _FString (new _String(temp.GenerateConsensusString())), false);
            }
          } else { // argument == _String("CONSENSUS")
            long seqID = ProcessNumericArgument ((_String*)parameters(2),chain.nameSpacePrefix);
            
            if (filter_source) {
              if (seqID>=0 && seqID < filter_source->NumberSpecies()) {
                stVar->SetValue (new _FString (filter_source->GetSequenceCharacters(seqID)),false);
              } else  if (seqID >= -4 && seqID <= -1) {
                _SimpleList indices, map, counts;
                long uniqueSequences = filter_source->FindUniqueSequences(indices, map, counts, -seqID - 1);
                _AssociativeList * parameterInfo = new _AssociativeList;
                parameterInfo->MStore ("UNIQUE_SEQUENCES",             new _Constant (uniqueSequences), false);
                parameterInfo->MStore ("UNIQUE_INDICES",            new _Matrix (indices), false);
                parameterInfo->MStore ("SEQUENCE_MAP",          new _Matrix (map), false);
                parameterInfo->MStore ("UNIQUE_COUNTS",      new _Matrix  (counts), false);
                stVar->SetValue (parameterInfo,false);
              }
            } else { // filter_source
              if (seqID>=0 && seqID < dataset_source->NoOfSpecies()) {
                stVar->SetValue (new _FString (dataset_source->GetSequenceCharacters(seqID)),false);
              }
            }
          } // else numeric cases
        }
          break;
          
        case 4UL : {
          if (filter_source) {
            long seq  = ProcessNumericArgument ((_String*)parameters(2),chain.nameSpacePrefix),
            site = ProcessNumericArgument ((_String*)parameters(3),chain.nameSpacePrefix);
            
            if (site >=0 && site<filter_source->GetPatternCount()) {
              if ( seq>=0 && seq<filter_source->NumberSpecies()) {
                _Matrix             * res = new _Matrix (filter_source->GetDimension (true), 1, false, true);
                
                _Parameter          onlyTheIndex = 0.0;
                checkParameter      (getDataInfoReturnsOnlyTheIndex,onlyTheIndex,0.0);
                
                
                _String             character (filter_source->RetrieveState(site, seq));
                long                theValue = filter_source->Translate2Frequencies (character, res->theData,  true);
                
                if (onlyTheIndex > 0.5) {
                  stVar->SetValue (new _Constant (theValue),false);
                  DeleteObject     (res);
                } else {
                  stVar->SetValue (res,false);
                }
              } else {
                _Parameter          count_gaps = 0.0;
                checkParameter      (hfCountGap,count_gaps,1.0);
                
                
                _Matrix * accumulator = new _Matrix (filter_source->GetDimension (true), 1, false, true),
                * storage     = new _Matrix (filter_source->GetDimension (true), 1, false, true);
                
                
                
                _String *buffer = filter_source->MakeSiteBuffer();
                
                for (long species_index = filter_source->NumberSpecies()-1; species_index >= 0; species_index --) {
                  filter_source->RetrieveState(site,species_index,*buffer, false);
                  filter_source->Translate2Frequencies (*buffer, storage->theData,  count_gaps >= 0.5);
                  *accumulator += *storage;
                }
                DeleteObject (storage);
                stVar -> SetValue (accumulator, false);
                
                DeleteObject (buffer);
                
              }
            } else {
              errMsg =  _String (site) & " is an invalid site index";
            }
          } else {
            errMsg = "This set of options is not supported for DataSet arguments";
          }
        }
          break;
          
        case 5UL: {
          if (filter_source) {
            long seq1  = ProcessNumericArgument ((_String*)parameters(2),chain.nameSpacePrefix),
            seq2  = ProcessNumericArgument ((_String*)parameters(3),chain.nameSpacePrefix);
            if ( seq1>=0 && seq2 >=0 && seq1< filter_source->NumberSpecies() && seq2 <filter_source->NumberSpecies()) {
              _String* resFlag = (_String*)parameters(4);
              _Matrix * res;
              
              if (pcAmbiguitiesAverage.Equal (resFlag)) {
                res = filter_source->ComputePairwiseDifferences (seq1,seq2,kAmbiguityHandlingAverageFrequencyAware);
              } else if (pcAmbiguitiesResolve.Equal (resFlag)) {
                res = filter_source->ComputePairwiseDifferences (seq1,seq2,kAmbiguityHandlingResolve);
              } else if (pcAmbiguitiesSkip.Equal (resFlag)) {
                res = filter_source->ComputePairwiseDifferences (seq1,seq2,kAmbiguityHandlingSkip);
              } else {
                res = filter_source->ComputePairwiseDifferences (seq1,seq2,kAmbiguityHandlingResolveFrequencyAware);
              }
              
              stVar->SetValue (res,false);
            } else {
              errMsg = _String (seq1).Enquote() & "," & _String (seq2).Enquote() & " is an invalid sequence pair specification.";
            }
          } else {
            errMsg = "This set of options is not supported for DataSet arguments";
          }
        }
        break;
      } // switch
    } else {
      errMsg = receptacle_name.Enquote() & " is not a valid receptacle identifier";
    }
          
          
  }
  if (errMsg.sLength) {
    errMsg = errMsg & " in call to GetDataInfo ";
    WarnError (errMsg);
  }
}

//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase47 (_ExecutionList& chain) {
    chain.currentCommand++;

    _String *arg1 = (_String*)parameters(0),
             *arg2 = (_String*)parameters(1),
              errMsg;

    long    k = FindLikeFuncName(AppendContainerName(*arg1, chain.nameSpacePrefix));

    if (k<0L) {
        _String  litArg = ProcessLiteralArgument (arg1,chain.nameSpacePrefix);
        k = FindLikeFuncName (litArg);
        if (k<0) {
            errMsg = *arg1 & " is not a defined likelihood function ID ";
        }
    }

    if (errMsg.sLength == 0UL) {
        _LikelihoodFunction * lf   = (_LikelihoodFunction *) likeFuncList (k);
        _String         callBack   = ProcessLiteralArgument (arg2,chain.nameSpacePrefix);
        k = FindBFFunctionName (callBack);

        if (k<0) {
            errMsg = arg2->Enquote() & " is not a defined user batch language function ";
        } else {
            if (GetBFFunctionArgumentCount(k)!=2L) {
                errMsg = arg2->Enquote() & " callback function must depend on 2 parameters ";
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

void      _ElementaryCommand::ExecuteCase52 (_ExecutionList& chain) {

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
        errMsg = *(_String*)parameters (4) & " must either evaluate to a positive integer or be a non-emptyString string of root states";
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
                        _String* theExclusions = &emptyString;
                        
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
                                            ds->InsertName (_String ("Root"),0L);
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
                                        _SimpleList     h,v;

                                        newFilter->SetFilter     (ds,unitSize,h,v,false);
                                        newFilter->SetExclusions (theExclusions,true);
                                        newFilter->SetupConversion ();

                                        /*char buffer[255];
                                        snprintf (buffer, sizeof(buffer),"%d %d\n",siteCount, newFilter->GetSiteCount(),unitSize);
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

                                            long       filterID = StoreDataFilter (simulationFilter, newFilter);
                                          
                                            spawningTree->SetUp();
                                            spawningTree->InitializeTreeFrequencies((_Matrix*)freqVar->Compute(),true);
                                          
                                            _String filter_specification = *GetFilterName (filterID) & spawningTree->GetName()->Enquote(',') & *freqVar->GetName();
                                          
                                            _LikelihoodFunction lf (filter_specification, nil);

                                            if (terminateExecution) {
                                                return;
                                            }

                                            bool    doInternals = false;

                                            if (parameters.lLength>5) {
                                                doInternals = (ProcessNumericArgument ((_String*)parameters (5),chain.nameSpacePrefix)>0.5);
                                            }


                                            _String spoolFile;

                                            FILE*   mainFile = nil;

                                            errMsg = emptyString;

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
                                                lf.Simulate (*simDataSet, exclusions, catValues, catNames, rootStates, doInternals?(mainFile?&spoolFile:&emptyString):nil);
                                                SetStatusLine ("Idle");

                                                catValVar->SetValue(catValues, false);
                                                catNameVar->SetValue(catNames, false);

                                                StoreADataSet (simDataSet, simName);
                                                DeleteObject (simName);
                                                DeleteDataFilter (filterID);
                                                errMsg = emptyString;
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

bool      _ElementaryCommand::Execute    (_ExecutionList& chain) {
  
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
                    tr = new _TheTree (treeIdent,(_TreeTopology*)formRes);
                } else if (formRes->ObjectClass () == TREE) {
                    for (unsigned long i = 0; i < leftOverVars.lLength; i++) {
                        //printf ("%s\n", LocateVar(leftOverVars.lData[i])->GetName()->sData);
                        DeleteVariable(leftOverVars.lData[i], true);
                    }
                    leftOverVars.Clear();
                    tr = new _TheTree (treeIdent,(_TheTree*)formRes);
                }
            }
        } else {
            tr = new _TheTree (treeIdent,treeString,false);
        }

        if (!tr) {
            WarnError ("Illegal right hand side in call to Tree id = ...; it must be a string, a Newick tree spec or a topology");
            return false;
        }
      
        if (leftOverVars.lLength) { // mod 02/03/2003 - the entire "if" block
            _SimpleList indep, dep, holder;
            {
                _AVLList    indepA (&indep),
                            depA   (&dep);

                tr->ScanContainerForVariables (indepA,depA);
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

    case HY_HBL_COMMAND_FPRINTF: { // print stuff to file (or stdout)
        return HandleFprintf(chain);
    }

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

    case 14: {
      // a return statement
    
      if (parameters.lLength) {
        
        _Formula * expression = nil;
        _String  * errMsg     = nil;
        try {
          
          
          if (simpleParameters.lLength < 2) {
            
            expression = new _Formula;
            //printf ("Namespace: %x\nCode: %s\n", chain.nameSpacePrefix, ((_String*)parameters(0))->sData);
            
            _FormulaParsingContext fpc (nil,  chain.nameSpacePrefix);
            long status = Parse (expression, *(_String*)parameters(0), fpc, nil);
            
            if (status== HY_FORMULA_EXPRESSION) {
              if (fpc.isVolatile() == false) {
                simpleParameters<<(long)expression;
                appendCompiledFormulae (expression);
                expression = nil;
              }
            } else {
                errMsg = new _String ("Invalid return statement");
                throw 0;
            }
          }
          
          _PMathObj ret_val = nil;
          // important to store the return value in a local variable
          // because chain.result may be overwritten by recursive calls to
          // this function
         
          if (expression) {
            //printf ("Return interpreted\n");
            ret_val = expression->Compute();
          }
          else{
            //printf ("Return compiled %d\n", ((_Formula*)simpleParameters(1))->GetList().lLength);
            ret_val = ((_Formula*)simpleParameters(1))->Compute();
          }
          
          DeleteObject (chain.result);
          
          chain.result = ret_val;
          if (ret_val) {
            chain.result->AddAReference();
          }
          
          if (expression) {
            delete (expression);
          }
        }
        catch (int e) {
          if (expression)
            delete expression;
          if (errMsg)
            WarnError (errMsg);
          return false;
        }
      }
      
      chain.currentCommand = simpleParameters(0);
      if (chain.currentCommand<0) {
        chain.currentCommand = 0x7fffffff;
      }
      break;
   }
      

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

        if (StoreADataSet (mergeResult, resultName) && simpleParameters(0)<0) {
            // purge all the datasets except the resulting one
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

    case HY_HBL_COMMAND_NESTED_LIST:
      chain.currentCommand++;
      {
        ((_ExecutionList*)parameters.GetItem(0))->Execute(&chain);
      }
      break;
      
    default:
        chain.currentCommand++;
    }

    return true;
}

//____________________________________________________________________________________


const _String   _ElementaryCommand::FindNextCommand  (_String& input, bool useSoftTrim)
{

    long    index     = input.Length();
  
    if (index == 0L) {
      return emptyString;
    }
    
    bool    isStringDouble  = false,
            isStringSingle  = false,
            skipping  = false;

    char    isComment = 0;


    long    scopeIn     = 0,
            matrixScope = 0,
            parenIn  = 0,
            bracketIn   = 0,
            saveSI = _String::storageIncrement;

    _SimpleList isDoWhileLoop;

    if ((input.sLength >> 2) > saveSI) {
        _String::storageIncrement = (input.sLength >> 2);
    }

    _String result (128L,true);

    char    lastChar = 0;

 
    // non printable characters at the end ?
    while (index>=0 && !isprint(input[--index])) ;
    input.Trim (0,index, useSoftTrim);

    for (index = 0L; index<input.Length(); index++) {
        char c = input.sData[index];

        if (!(isStringDouble || isStringSingle) && c=='\t') {
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
            if (!(isStringDouble || isStringSingle) && c=='/') {
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


        // skip spaces, except for special cases, like return and data set filters
      
        if (!(isStringDouble || isStringSingle) && isspace(c)) {
          
          // skip/compress spaces, unless we are in a higher level HBL statement
          // where spaces can't be compressed
          // examples include
          // DataSet|DataSetFilter|return|LikelihoodFunction (something)
          // need to maintain spaces for this to work appropriately
          
          
            /*if (index >= 6 && input.getChar(index-1) == 'n'
                    && input.getChar(index-2) == 'r'
                    && input.getChar(index-3) == 'u'
                    && input.getChar(index-4) == 't'
                    && input.getChar(index-5) == 'e'
                    && input.getChar(index-6) == 'r') {
                if (index == 6 || (index > 6 && !(isalnum(input.getChar(index-7)) || input.getChar(index-7) == '_'))) {
                    result<<' ';
                }
            }*/
          
            if (!skipping && index > 0) {
              _String lookback = input.Cut (MAX (0, index - 20), index-1);
              long trie_match = _HY_HBL_KeywordsPreserveSpaces.FindKey(lookback.Flip(), nil, true);
              if (trie_match != HY_TRIE_NOTFOUND) {
                long matched_length = _HY_HBL_KeywordsPreserveSpaces.GetValue(trie_match);
                if (matched_length == index || !(isalnum(input.getChar(index-matched_length-1)) || input.getChar(index-matched_length-1) == '_')) {
                  result << ' ';
                }
              }
            }


            skipping = true;
            continue;
        }

        if (skipping&&( isalpha(c) || c=='_') && (isalnum(lastChar) || lastChar=='_')) {
          // this is meant to determine that we are at the beginning of a new ident-like
          // token and insert a space
            result<<' ';
        }

        skipping = false;

        result<<c;

        if ((isStringDouble || isStringSingle) && c == '\\') {
            result<< input.getChar(++index);
            continue;
        }

        // are we inside a string literal?

        if (c=='"') {
          if (!isStringSingle) {
            isStringDouble = !isStringDouble;
            lastChar = 0;
            continue;
          }
        } else {
          if (c == '\'') {
            if (!isStringDouble) {
              isStringSingle = !isStringSingle;
              lastChar = 0;
              continue;
            }
          }
        }

        if (isStringDouble || isStringSingle) {
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
                input = emptyString;
                return emptyString;
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

    if (scopeIn||isStringDouble||isStringSingle||isComment == 1||parenIn||matrixScope) {
        if (result!='}') {
            WarnError (_String("Expression appears to be incomplete/syntax error. Scope: ") &scopeIn & ", paretheses depth: "
                       & parenIn & ", matrix scope: " & matrixScope & '.' & (isStringDouble?" In a \"\" literal. ":emptyString)
                       & (isStringSingle?" In a '' literal. ":emptyString) &
                       (isComment == 1? " In a comment ":emptyString) & '\n' & input);
            input = emptyString;
            return emptyString;
        } else {
            result = emptyString;
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
            result.DuplicateErasing (&emptyString);
        } else {
            result.Trim(lastChar,result.sLength-1-lastChar);
        }
    }

    if (index<input.Length()-1) {
        input.Trim (index+1,-1, useSoftTrim);
    } else if (useSoftTrim) {
        input.sLength = 0;
    } else {
        input.DuplicateErasing (&emptyString);
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


bool       _ElementaryCommand::BuildFor (_String&source, _ExecutionList&target,  _List * pieces)

// the for loop becomes this:
// initialize
// if (condition) then
// else go to end+2
// end
// increment
// goto if(condition)

{
  if (pieces)
    return MakeGeneralizedLoop ((_String*)pieces->GetItem(0),(_String*)pieces->GetItem(1),(_String*)pieces->GetItem(2),true,source,target);
  else
    return MakeGeneralizedLoop (nil,nil,nil,true,source,target);
}

//____________________________________________________________________________________

bool    _ElementaryCommand::BuildWhile          (_String&source, _ExecutionList&target,  _List * pieces)
{
    if (pieces)
      return MakeGeneralizedLoop (nil,(_String*)pieces->GetItem(0),nil,true,source,target);
    else
      return MakeGeneralizedLoop (nil,nil,nil,true,source,target);
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

    long    mark1 = source.FirstNonSpaceFollowingSpace(),
            mark2 = source.FindTerminator(mark1, '='); ;

 
    if (mark1==-1 || mark2==-1 || mark2 - 1 <= mark1 ) {
        WarnErrorWhileParsing ("DataSet declaration missing a valid identifier", source);
        return false;
    }

    _String dsID (source,mark1,mark2-1);
    // now look for the opening paren

    mark1 = source.Find ('(',mark2,-1);

    _ElementaryCommand dsc;
    _String            oper (source,mark2+1,mark1-1);

    if (oper ==  _String("ReadDataFile") || oper == _String ("ReadFromString")) { // a switch statement if more than 1
        _List pieces;
        ExtractConditions (source,mark1+1,pieces,',');
        if (pieces.lLength!=1UL) {
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
        ExtractConditions (source,mark1+1,pieces,',');
        if ( pieces.lLength>4UL || pieces.lLength==0UL ) {
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
        ExtractConditions (source,mark1+1,pieces,',');
        if (pieces.lLength==0UL) {
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
            ExtractConditions (source,mark1+1,pieces,',');
            if (pieces.lLength>3UL || pieces.lLength==0UL) {
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
            ExtractConditions (source,mark1+1,pieces,',');
            if ((pieces.lLength>7)||(pieces.lLength<4UL)) {
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
            ExtractConditions (source,0,args,',');
            if (args.lLength>=7UL) {
                _ElementaryCommand * cv = new _ElementaryCommand (20);
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
bool    _ElementaryCommand::ConstructChoiceList(_String&source, _ExecutionList&target) {
    _List args;
  
  
    ExtractConditions (source,blChoiceList.sLength,args,',');
    if (args.lLength<5UL) {
        WarnError  ("ChoiceList needs at least 5 arguments");
        return false;
    }
    _ElementaryCommand *cv = new _ElementaryCommand (32);

    cv->parameters<<args(0);
    //((_String*)args.lData[1])->StripQuotes();
    cv->parameters<<args(1)
                  <<args(2)
                  <<args(3);

    if  (args.lLength>5UL) {
        _List * choices = new _List;
        for (long k = 4L; k<args.lLength-1; k+=2) {
            ((_String*)args.lData[k])->StripQuotes();
            ((_String*)args.lData[k+1])->StripQuotes();
            _List * thisChoice = new _List;
            *thisChoice << args(k);
            *thisChoice << args(k+1);
            *choices < thisChoice;
        }
        cv->parameters < choices;
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
  
    ExtractConditions (source,blReplicate.sLength,args,',');
    if (args.lLength<2) {
        _String errMsg ("Expected: ReplicateConstraint (\"constraint to be replicated in terms of this1,...,thisn and wildcard *\", list of n variables to put in place of this1, this2, ... thisn);");
        acknError (errMsg);
        return false;
    }
  
    _ElementaryCommand cv;
    cv.code = 26;
    cv.parameters << args;
    target&& &cv;
    return true;
}

//____________________________________________________________________________________

bool    _ElementaryCommand::ConstructTree (_String&source, _ExecutionList&target)
// Tree   treeid = (...) or Topology = (...);
{
    long    mark1 = source.FirstSpaceIndex(0,-1,1);
    if (mark1 > 0) {
      mark1 = source.FirstNonSpaceIndex (mark1 + 1, -1);
    }
    
  
    long    mark2 = source.FindTerminator(mark1, "=");
    long    mark3 = mark2;

    if ( mark1 < 0 || mark2 < 0 || mark2 - mark1 < 1) {
        acknError ("Tree declaration missing a valid identifier");
        return false;
    }

    _String dsID = source.Cut (mark1,mark2-1);
    // now look for the opening paren
  
    //(long& from, char open, char close, bool respectQuote, bool respectEscape)
    mark3 = source.ExtractEnclosedExpression (mark1, '(', ')', true, true);
  

    if (mark1 < 0 || mark3 < 0 || mark3 <= mark1) {
        if (source.Find(getDString)==-1) {
            mark1 = mark2+1;
            mark3 = source.FindTerminator (mark1,";")-1;
        } else {
            source = getDString;
            mark1 = 0;
            mark3 = -1;
        }
    }

    _ElementaryCommand * dsc = new _ElementaryCommand(source.startswith(blTree)?7:54);
 
    dsc->parameters&&(&dsID);
    dsc->parameters.AppendNewInstance(new _String(source,mark1,mark3));
  
    dsc->addAndClean(target,nil,0);
    return true;
}



//____________________________________________________________________________________

bool    _ElementaryCommand::ConstructDataSetFilter (_String&source, _ExecutionList&target)
// DataSetFilter      dataSetFilterid = CreateFilter (datasetid;unit;vertical partition; horizontal partition; alphabet exclusions);

{
    // first we must segment out the data set name

    long  mark1 = source.FirstNonSpaceFollowingSpace (0,-1,1),
          mark2 = source.FindTerminator(mark1+1, "=");

    _String dsID    (source,mark1,mark2-1),
            command;

    if ( mark1==-1 || mark2==-1 || dsID.Length()==0) {
        acknError ("DataSetFilter declaration missing a valid identifier");
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


    ExtractConditions (source,mark1+1,pieces,',');
    if (!(pieces.lLength>=2UL || (pieces.lLength == 1UL && dsf->code == 6))) {
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
    ExtractConditions (source,mark1+1,pieces,',');

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
    model->parameters&&(&modelID);
    model->addAndClean (target,&pieces,0);
    return true;

}


//____________________________________________________________________________________

/*bool    _ElementaryCommand::ConstructFprintf (_String&source, _ExecutionList&target)

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
}*/

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
        p = allowedFormats.FindObject(argDesc(f));
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
bool    _ElementaryCommand::ConstructExecuteCommands (_String&source, _ExecutionList&target) {

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

    if (PeekFilePath()) {
        exc->parameters && *PeekFilePath();
    } else {
        exc->parameters && & emptyString;
    }

    if (pieces.lLength >1) {
        if (*(_String*)pieces(1) == _String("compiled")) {
            exc->simpleParameters << 1;
        } else {
          if (*(_String*)pieces(1) == _String("enclosing_namespace")) {
            exc->parameters.Delete(1);
            exc->parameters && & emptyString;
          } else {
            exc->parameters << pieces(1);
            if (pieces.lLength > 2) {
                exc->parameters << pieces(2);
            }
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

    /*else {
        _String *s0 = (_String*)pieces(0),
                 *s1 = (_String*)pieces(1);

        if (!(s0->IsValidIdentifier()&&((s1->sLength>2&&s1->getChar(0)=='"'&s1->getChar(s1->sLength-1)=='"') || s1->IsValidIdentifier()))) {
            WarnError (_String ("Both ") & *s0 & " and " & *s1 & " must be valid identifiers in call to GetInformation.");
            return     false;
        }
    }*/

    _ElementaryCommand * sp = makeNewCommand(37);
    sp->addAndClean (target, &pieces, 0);
    return true;
}

//____________________________________________________________________________________
bool    _ElementaryCommand::ConstructLF (_String&source, _ExecutionList&target)
// syntax: LikelihoodFunction id = (filter1, tree1, ..., filterN, treeN, optional compute template)
// or LikelihoodFunction3 id = (filter1, tree1, freq1, ... filterN, treeN, freqN, optional compute template)
{
    long    mark1 = source.FirstNonSpaceFollowingSpace(),
            mark2 = mark1 > 0 ? source.FindTerminator (mark1 + 1, "=") : 0;

    if ( mark1==-1 || mark2==-1 || mark1+1 > mark2  ) {
        acknError ("Likelihood function declaration missing a valid identifier");
        return false;
    }

    _String lfID (source,mark1,mark2-1);
    // now look for the opening paren

    _List pieces;
    mark2 ++;
    mark1 = source.ExtractEnclosedExpression(mark2, '(', ')', true, true);

    if ( mark1==-1 || mark2==-1 || mark1<mark2 ) {
        WarnError ("Expected: Likelihood Function ident = (tree1, datasetfilter1,...)");
        return false;
    }

    ExtractConditions (source,mark2+1,pieces,',');
   _ElementaryCommand*  dsc = new _ElementaryCommand (11);
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
  
  
    bool    isFFunction = source.beginswith (blFFunction),
            isLFunction = source.beginswith (blLFunction),
            isNameSpace = source.beginswith (blNameSpace);
  
    if (!isNameSpace) {
      if (isInFunction == _HY_FUNCTION) {
        WarnError ("Nested function declarations are not allowed");
        return false;
      }
      
    }

  

    long    mark1 = source.FirstNonSpaceIndex(isNameSpace ? blNameSpace.sLength: ((isFFunction||isLFunction)?blFFunction.sLength:blFunction.sLength),-1,1),
            mark2 = source.Find (isNameSpace ? '{' : '(', mark1, -1);


    if ( mark1==-1 || mark2==-1 || mark1>mark2-1) {
        WarnError      (_String("Function declaration missing a valid function identifier or parameter list.\n-----------\n") & source & "\n-----------\n");
      isInFunction = _HY_NO_FUNCTION;
      return false;
    }

    _String*    funcID  = new _String(source.Cut (mark1,mark2-1));

    if (!funcID->IsValidIdentifier(true)) {
      WarnError      (_String("Not a valid function/namespace identifier '") & *funcID & "'");
      isInFunction = _HY_NO_FUNCTION;
      return false;
    }
    
    *funcID = chain.AddNameSpaceToID (*funcID);

    // now look for the opening paren

    if (!isNameSpace) {
      isInFunction = _HY_FUNCTION;

      
      if ((mark1=FindBFFunctionName(*funcID)) >= 0L) {
        ReportWarning (_String("Overwritten previously defined function:'") & *funcID & '\'');
      }
      
      _List       arguments;
      _SimpleList argument_types;

      long upto = ExtractConditions (source,mark2+1,arguments,',',false);


      if (upto==source.sLength || source[upto]!='{' || source[source.sLength-1]!='}') {
          WarnError (_String("Function declaration is missing a valid function body."));
          isInFunction= _HY_NO_FUNCTION;
          return false;
      }

      _String extraNamespace;
      if (isLFunction)
          extraNamespace = _HYGenerateANameSpace();
      
      for (long k = 0UL; k < arguments.lLength; k++) {
        
          _String*   namespaced = new _String(chain.AddNameSpaceToID (*(_String*)arguments(k), & extraNamespace));
          if (namespaced->getChar(namespaced->sLength - 1L) == '&') {
            namespaced->Trim(0,namespaced->sLength-2);
            argument_types << BL_FUNCTION_ARGUMENT_REFERENCE;
          } else {
            argument_types << BL_FUNCTION_ARGUMENT_NORMAL;
          }
          arguments.Replace (k,namespaced,false);
      }
    

      _String          sfunctionBody (source, upto+1,source.Length()-2);
      _ExecutionList * functionBody;
      
      if (isLFunction) {
          _String * existing_namespace = chain.GetNameSpace();
          if (existing_namespace) {
              extraNamespace = *existing_namespace & '.' & extraNamespace;
          }
          functionBody = new _ExecutionList (sfunctionBody,&extraNamespace,true);
          if (existing_namespace) {
            functionBody->enclosingNamespace = *existing_namespace;
          }
      }
      else {
          functionBody = new _ExecutionList (sfunctionBody,chain.GetNameSpace(),true);
      }
      

      //  take care of all the return statements
      while (returnlist.lLength) {
          ((_ElementaryCommand*)(*functionBody)(returnlist(0)))->simpleParameters<<functionBody->lLength;
          returnlist.Delete(0);
      }


      if (mark1>=0) {
          batchLanguageFunctions.Replace (mark1, functionBody, false);
          batchLanguageFunctionNames.Replace (mark1, funcID, false);
          batchLanguageFunctionParameterLists.Replace (mark1, &arguments, true);
          batchLanguageFunctionParameterTypes.Replace (mark1, &argument_types, true);
        batchLanguageFunctionClassification.lData[mark1] = isLFunction ? BL_FUNCTION_LOCAL :( isFFunction? BL_FUNCTION_SKIP_UPDATE :  BL_FUNCTION_ALWAYS_UPDATE);
      } else {
          batchLanguageFunctions.AppendNewInstance(functionBody);
          batchLanguageFunctionNames.AppendNewInstance(funcID);
          batchLanguageFunctionNamesIndexed.Insert (new _String(*funcID), batchLanguageFunctions.countitems()-1, false, true);
          batchLanguageFunctionParameterLists &&(&arguments);
          batchLanguageFunctionParameterTypes &&(&argument_types);
          batchLanguageFunctionClassification <<(isLFunction ? BL_FUNCTION_LOCAL :( isFFunction? BL_FUNCTION_SKIP_UPDATE :  BL_FUNCTION_ALWAYS_UPDATE));
      }
    } else {
      if (mark2 == source.sLength || source[mark2]!='{' || source[source.sLength-1]!='}') {
        WarnError (_String("Namespace declaration is missing a body."));
        isInFunction= _HY_NO_FUNCTION;
        return false;
      }
      _String          namespace_text (source, mark2+1,source.Length()-2);
      bool             success = false;
      
      
      _ExecutionList   * namespace_payload = new _ExecutionList (namespace_text, funcID, false, &success);
      
      if (success) {
        _ElementaryCommand * nested_list = new _ElementaryCommand (HY_HBL_COMMAND_NESTED_LIST);
        nested_list->parameters.AppendNewInstance(namespace_payload);
        chain.AppendNewInstance(nested_list);
      } else {
        DeleteObject (namespace_payload);
        return false;
      }

    }


    isInFunction = _HY_NO_FUNCTION;
    return true;
}

//____________________________________________________________________________________
bool    _ElementaryCommand::ConstructReturn (_String&source, _ExecutionList&target) {
// syntax: return <statement>
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
        WarnError (_String("Could not read batch file '") & fName & "'.\nPath stack:\n\t" & GetPathStack("\n\t"));
    } else {
        _String beavis (f);

        if (beavis.beginswith("#NEXUS",false)) {
            ReadDataSetFile (f,1,nil,&fName, nil, &defaultTranslationTable, &target);
        } else {
            target.BuildList (beavis);
            target.sourceFile = fName;
        }
        fclose (f);
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
