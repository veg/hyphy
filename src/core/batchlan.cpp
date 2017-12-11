/*

 HyPhy - Hypothesis Testing Using Phylogenies.

 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
 Art FY Poon    (apoon42@uwo.ca)
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

#include <string.h>
#include <ctype.h>
#include <time.h>


#include "likefunc.h"
#include "batchlan.h"
#include "polynoml.h"
#include "scfg.h"
#include "bayesgraph.h"
#include "avllistx.h"
#include "global_object_lists.h"
#include "global_things.h"
#include "time_difference.h"
#include "global_things.h"
#include "hy_string_buffer.h"
#include "tree_iterator.h"




using namespace hy_global;
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
batchLanguageFunctions,
batchLanguageFunctionNames,
batchLanguageFunctionParameterLists,
batchLanguageFunctionParameterTypes,
compiledFormulaeParameters,
modelNames,
executionStack,
loadedLibraryPathsBackend;


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
globalPolynomialCap             ("GLOBAL_POLYNOMIAL_CAP"),
                                enforceGlobalPolynomialCap      ("ENFORCE_GLOBAL_POLYNOMIAL_CAP"),
                                dropPolynomialTerms             ("DROP_POLYNOMIAL_TERMS"),
                                maxPolyTermsPerVariable         ("MAX_POLY_TERMS_PER_VARIABLE"),
                                maxPolyExpIterates              ("MAX_POLYNOMIAL_EXP_ITERATES"),
                                polyExpPrecision                ("POLYNOMIAL_EXP_PRECISION"),
                                explicitFormMExp                ("EXPLICIT_FORM_MATRIX_EXPONENTIAL"),
                                multByFrequencies               ("MULTIPLY_BY_FREQUENCIES"),
                                defFileString                   ("DEFAULT_FILE_SAVE_NAME"),
                                VerbosityLevelString            ("VERBOSITY_LEVEL"),
                                useLastDefinedMatrix            ("USE_LAST_DEFINED_MATRIX"),
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
                                simulationFilter                ("_SIM_INTERNAL_FILTER_"),
                                prefixDS                        ("DataSet_"),
                                prefixDF                        ("Partition_"),
                                prefixLF                        ("LF_"),
                                replaceTreeStructure            ("REPLACE_TREE_STRUCTURE"),
                                matrixEvalCount                 ("MATRIX_EXPONENTIATION_COUNTS"),
                                 _hyLastExecutionError           ("LAST_HBL_EXECUTION_ERROR"),

                                kBGMData                         ("BGM_DATA_MATRIX"),

                                gdiDFAtomSize                   ("ATOM_SIZE"),
                                marginalAncestors               ("MARGINAL"),
                                doLeavesAncestors               ("DOLEAVES"),
                                 dialogPrompt,
                                hy_scanf_last_file_path,
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
blGetInformation           ("GetInformation("),
blExecuteCommands      ("ExecuteCommands("),
blExecuteAFile         ("ExecuteAFile("),
blLoadFunctionLibrary      ("LoadFunctionLibrary("),
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

void        ReportMPIError                  (int, bool);

#endif

_hy_nested_check  isInFunction = _HY_NO_FUNCTION;

hyFloat  explicitFormMatrixExponential = 0.0,
            messageLogFlag                = 1.0;


extern      _String             MATRIX_AGREEMENT,
            ANAL_COMP_FLAG;


extern      _SimpleList         freeSlots;


_AVLList    loadedLibraryPaths  (&loadedLibraryPathsBackend);

_ExecutionList
*currentExecutionList = nil;

_List const _ElementaryCommand::fscanf_allowed_formats (new _String ("Number"),
                                                        new _String ("Matrix"),
                                                        new _String ("Tree"),
                                                        new _String ("String"),
                                                        new _String ("NMatrix"),
                                                        new _String ("Raw"),
                                                        new _String ("Lines")
                                                        );

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
        HandleApplicationError (errMsg);
    }
}

#define MPI_SEND_CHUNK 0xFFFFFFL

//____________________________________________________________________________________

void    MPISendString       (_String const& theMessage, long destID, bool isError)
{

    long    messageLength = theMessage.length(),
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
        ReportMPIError(MPI_Send(theMessage.get_str()+transferCount, MPI_SEND_CHUNK, MPI_CHAR, destID, HYPHY_MPI_STRING_TAG, MPI_COMM_WORLD),true);
        transferCount += MPI_SEND_CHUNK;
    }

    if (messageLength-transferCount) {
        ReportMPIError(MPI_Send(theMessage.get_str()+transferCount, messageLength-transferCount, MPI_CHAR, destID, HYPHY_MPI_STRING_TAG, MPI_COMM_WORLD),true);
    }

    //ReportMPIError(MPI_Send(&messageLength, 1, MPI_LONG, destID, HYPHY_MPI_DONE_TAG, MPI_COMM_WORLD),true);

    _FString*    sentVal = new _FString;
    sentVal->theString = (_String*)theMessage.makeDynamic();
    _Variable *   mpiMsgVar = CheckReceptacle (&hy_env::mpi_last_sent_message, kEmptyString, false);
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
  
    // nonaggressive polling mode
  
    int message_received = 0;
    while (! message_received) {
      MPI_Iprobe (senderT, HYPHY_MPI_SIZE_TAG, MPI_COMM_WORLD, &message_received, MPI_STATUS_IGNORE);
      usleep (100);
    }

    // nonaggressive polling mode
  
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

        theMessage = new _String ((unsigned long)messageLength);

        senderT = senderID = status.MPI_SOURCE;

        while (messageLength-transferCount>MPI_SEND_CHUNK) {
            ReportMPIError(MPI_Recv(theMessage->get_str()+transferCount, MPI_SEND_CHUNK, MPI_CHAR, senderT, HYPHY_MPI_STRING_TAG, MPI_COMM_WORLD,&status),false);
            MPI_Get_count (&status,MPI_CHAR,&actualReceived);
            if (actualReceived!=MPI_SEND_CHUNK) {
                HandleApplicationError ("Failed in MPIRecvString - some data was not properly received\n");
            }
            //return    nil;
            transferCount += MPI_SEND_CHUNK;
        }

        if (messageLength-transferCount) {
            ReportMPIError(MPI_Recv(theMessage->get_str()+transferCount, messageLength-transferCount, MPI_CHAR, senderT, HYPHY_MPI_STRING_TAG, MPI_COMM_WORLD,&status),false);
            MPI_Get_count (&status,MPI_CHAR,&actualReceived);
            if (actualReceived!=messageLength-transferCount) {
                HandleApplicationError ("Failed in MPIRecvString - some data was not properly received\n");
            }
            //return    nil;
        }
        //ReportMPIError(MPI_Recv(&messageLength, 1, MPI_LONG, senderT, HYPHY_MPI_DONE_TAG, MPI_COMM_WORLD,&status),false);

        if (isError) {
            HandleApplicationError (theMessage);
        }
    }
    return theMessage;
}
#endif


//____________________________________________________________________________________

const _String GetStringFromFormula (_String const* data,_VariableContainer* theP) {

    _Formula  nameForm (*data,theP);
    _PMathObj formRes = nameForm.Compute();

    if (formRes&& formRes->ObjectClass()==STRING) {
        return ((_FString*)formRes)->get_str();
    }

    return *data;
}

//____________________________________________________________________________________

bool    numericalParameterSuccessFlag = true;

hyFloat  _ProcessNumericArgumentWithExceptions (_String& data, _VariableContainer const* theP) {
    _String   errMsg;
    _Formula  nameForm (data,theP, &errMsg);

    if (errMsg.empty()) {
        _PMathObj formRes = nameForm.Compute();
        if (formRes&& formRes->ObjectClass()==NUMBER) {
            return formRes->Value();
        } else {
            if (formRes&& formRes->ObjectClass()==STRING) {
                return _String((_String*)((_FString*)formRes)->toStr()).to_float();
            } else {
                throw (data.Enquote() & " was expected to be a numerical argument.");
            }
        }
    } else {
        throw (errMsg);
    }

    return 0.0;

}


hyFloat  ProcessNumericArgument (_String* data, _VariableContainer const* theP, _ExecutionList* currentProgram) {

    numericalParameterSuccessFlag = true;

    try {
        return _ProcessNumericArgumentWithExceptions (*data, theP);
    } catch (_String const& err) {
        if (currentProgram) {
            currentProgram->ReportAnExecutionError (err);
        } else {
            HandleApplicationError(err);
        }
    }

    numericalParameterSuccessFlag = false;
    return 0.0;
}

//____________________________________________________________________________________

_PMathObj   ProcessAnArgumentByType (_String const* expression, _VariableContainer const* theP, long objectType, _ExecutionList* currentProgram)
{
    _String   errMsg;


    _Formula  expressionProcessor (*expression, theP, currentProgram?&errMsg:nil);

    if (errMsg.nonempty() && currentProgram) {
        currentProgram->ReportAnExecutionError (errMsg);
    }
    else {
        _PMathObj expressionResult = expressionProcessor.Compute(0,theP);
        if (expressionResult && (expressionResult->ObjectClass() & objectType)) {
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
      _String result (((_FString*)getString)->get_str());
      DeleteObject(getString);
      return result;
    }

    return kEmptyString;
}

//____________________________________________________________________________________

_AssociativeList*   ProcessDictionaryArgument (_String* data, _VariableContainer* theP, _ExecutionList* currentProgram) {
  return (_AssociativeList* )ProcessAnArgumentByType (data, theP, ASSOCIATIVE_LIST, currentProgram);
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
hyBLFunctionType   GetBFFunctionType  (long idx) {
  return (hyBLFunctionType) batchLanguageFunctionClassification.Element (idx);
}

//____________________________________________________________________________________
_String const ExportBFFunction (long idx, bool recursive) {


  _StringBuffer bf (8192UL);
  if (IsBFFunctionIndexValid(idx)) {

    _String hbf_name = GetBFFunctionNameByIndex (idx);
    _ExecutionList * body = &GetBFFunctionBody(idx);

    if (body->enclosingNamespace.nonempty()) {
      bf << "namespace " << body->enclosingNamespace << " {\n";
    }

    switch (GetBFFunctionType (idx)) {
      case kBLFunctionSkipUpdate:
        bf << blFFunction;
        break;
      case kBLFunctionLocal:
        bf << blLFunction;
        break;
      default:
        bf << blFunction;
    }


    bf << hbf_name << '(';

    long argument_count = GetBFFunctionArgumentCount (idx);
    _List * argument_list = &GetBFFunctionArgumentList (idx);
    for (long argument_id = 0; argument_id < argument_count; argument_id ++) {
      if (argument_id) {
        bf << ',';
      }


      bf << body->TrimNameSpaceFromID(*(_String*)argument_list->Element (argument_id));
      if (GetBFFunctionArgumentTypes (idx).GetElement(argument_id) == kBLFunctionArgumentReference) {
        bf << '&';
      }
    }
    bf << ") {\n" << body->sourceText << "\n}";

    if (body->enclosingNamespace.nonempty()) {
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
          bf << "\n/*----- Called function '"
          << *a_name
          << "' ------*/\n"
          << ExportBFFunction (FindBFFunctionName(*a_name), false)
          << "\n\n";
        }
      }
    }
  }

  return bf;

}

//____________________________________________________________________________________
void ClearBFFunctionLists (long start_here) {
  if (start_here > 0L && start_here < batchLanguageFunctionNames.countitems()) {

    _SimpleList delete_me (batchLanguageFunctionNames.countitems()-start_here, start_here, 1L);

    batchLanguageFunctionNames.DeleteList           (delete_me);
    batchLanguageFunctions.DeleteList               (delete_me);
    batchLanguageFunctionClassification.DeleteList  (delete_me);
    batchLanguageFunctionParameterLists.DeleteList  (delete_me);
    batchLanguageFunctionParameterTypes.DeleteList  (delete_me);
  }
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
            long idx = batchLanguageFunctionNames.FindObject (&test_id);
            if (idx >= 0) {
                //s = test_id;
                return idx;
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
   return batchLanguageFunctionNames.FindObject (&s);
}


//____________________________________________________________________________________
long    FindBgmName (_String const&s) {
    return bgmNamesList.FindObject (&s);
}



//__________________________________________________________
long  AddDataSetToList (_String& theName,_DataSet* theDS) {
    theName = GenerateUniqueObjectIDByType(theName, HY_BL_DATASET);
    long k = dataSetNamesList.FindObject (&kEmptyString);
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
                      if (((_String*)likeFuncNamesList(k))->nonempty()) {
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
        
        me->UnregisterListeners();

        if (lfID<likeFuncList.lLength-1) {
            DeleteObject(likeFuncList(lfID));
            likeFuncList.lData[lfID] = 0L;
            likeFuncNamesList.Replace(lfID,new _String,false);
        } else {
            likeFuncList.Delete(lfID);
            likeFuncNamesList.Delete(lfID);
            if (lfID)
                while (((_String*)likeFuncNamesList (--lfID))->empty()) {
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
        dataSetNamesList.Replace(dsID,new _String,false);
    } else {
        dataSetList.Delete(dsID);
        dataSetNamesList.Delete(dsID);
        if (dsID)
            while (((_String*)dataSetNamesList (--dsID))->empty()) {
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
            if (k != mdID && ((_String*)modelNames(k))->nonempty()) {
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
        modelNames.Replace(mdID,new _String,false);
    } else {
        modelNames.Delete(mdID);
        modelMatrixIndices.Delete (modelMatrixIndices.lLength-1);
        modelFrequenciesIndices.Delete (modelFrequenciesIndices.lLength-1);
        modelTypeList.Delete (modelTypeList.lLength-1);

        if (mdID)
            while (((_String*)modelNames (--mdID))->empty()) {
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

BaseRef     _ExecutionList::makeDynamic (void) const {
    _ExecutionList * Res = new _ExecutionList;

    //memcpy ((char*)Res, (char*)this, sizeof (_ExecutionList));

    Res->_List::Duplicate (this);
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

void        _ExecutionList::Duplicate   (BaseRefConst source) {
    _List::Duplicate    (source);

    _ExecutionList const* s = (_ExecutionList const*)source;

    if (s->result) {
        result=(_PMathObj)s->result->makeDynamic();
    }

    errorHandlingMode  = s->errorHandlingMode;
    errorState         = s->errorState;
}


//____________________________________________________________________________________
void    _ExecutionList::ReportAnExecutionError (_String errMsg, bool doCurrentCommand, bool appendToExisting) {
    if (doCurrentCommand) {
        _ElementaryCommand *theCommand = FetchLastCommand();
        if (theCommand) {
            errMsg = errMsg & " in to " & _String ((_String*)theCommand->toStr());
        }
    }
    errorState = true;
    switch (errorHandlingMode) {
        case HY_BL_ERROR_HANDLING_SOFT:
            if (appendToExisting) {
              _FString * existing = (_FString*) FetchObjectFromVariableByType(&_hyLastExecutionError, STRING);
              if (existing) {
                errMsg = existing->get_str() & '\n' & errMsg;
              }
            }
            setParameter(_hyLastExecutionError, new _FString (errMsg, false), nil, false);

            break;
        default:
            HandleApplicationError (errMsg);
    }
}

//____________________________________________________________________________________
_String*    _ExecutionList::FetchFromStdinRedirect (void)
// grab a string from the front of the input queue
// complain if nothing is left
{
    if (!stdinRedirect) {
        HandleApplicationError ("No input buffer was given for a redirected standard input read.");
        return new _String;
    }
    long d = stdinRedirect->First();
    if (d<0) {
        HandleApplicationError ("Ran out of input in buffer during a redirected standard input read.");
        return new _String;
    }
    _String *sendBack = (_String*)stdinRedirect->GetXtra (d);
    sendBack->AddAReference();
    stdinRedirect->Delete ((*(_List*)stdinRedirect->dataList)(d),true);
    return sendBack;
}

//____________________________________________________________________________________

_String  const     _ExecutionList::GetFileName     (void)  const {
    if (sourceFile.nonempty()) {
        return sourceFile;
    } else {
        _String const *top_path = PeekFilePath();
        if (top_path)
          return *top_path;
    }
    return kEmptyString;
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

    _FString            cfp (PeekFilePath() ? *PeekFilePath () :kEmptyString),
                        *stashed = (_FString*)hy_env::EnvVariableGet(hy_env::path_to_current_bf, STRING);

    if (stashed) {
        stashed = (_FString*)stashed->makeDynamic();
    }

    hy_env::EnvVariableSet(hy_env::path_to_current_bf, &cfp, true);

    DeleteObject        (result);
    result               = nil;
    currentExecutionList = this;
    currentCommand       = 0;

    terminate_execution  = false;

    while (currentCommand<lLength) {
        if (doProfile == 1 && profileCounter) {
            long        instCounter = currentCommand;
            hyFloat  timeDiff    = 0.0;

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

        if (terminate_execution) {
            break;
        }
    }
    currentCommand = callPoints.lData[callPoints.lLength-1];
    callPoints.Delete (callPoints.lLength-1);
    currentExecutionList = stashCEL;

    if (stashed) {
        hy_env::EnvVariableSet(hy_env::path_to_current_bf, stashed, false);
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

    if (fName && !terminate_execution) {
        f = batchLanguageFunctionNames.FindObject (fName);
    }

    terminate_execution      = false;

    ClearBFFunctionLists    (g);

  return f;
}

//____________________________________________________________________________________

bool        _ExecutionList::TryToMakeSimple     (void)
{
    _SimpleList     varListAux,
                    formulaeToConvert,
                    parseCodes;
  
  
    _AVLList        varList (&varListAux);

    long            stackDepth  = 0;

    bool            status      = true;

    for (unsigned long k = 0; k<lLength && status; k++) {
        _ElementaryCommand * aStatement = (_ElementaryCommand*)(*this)(k);
        switch (aStatement->code) {
        case 0: {
            _String * formulaString = (_String*)aStatement->parameters(0);

           if ((*formulaString) (-1) != '}') {
                _Formula *f  = new _Formula,
                         *f2 = new _Formula;


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
                            f->GetIthTerm (0)->SetAVariable(mx->GetIndex());
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
                            varList.InsertNumber(fpc.assignmentRefID());
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
        varList.ReorderList();
        cli->values = new _SimpleFormulaDatum[varList.countitems()+1];
        cli->stack  = new _SimpleFormulaDatum[stackDepth+1];

        _SimpleList  avlData;
        _AVLListX    avlList (&avlData);

        for (unsigned long fi = 0; fi < formulaeToConvert.lLength; fi++) {
            ((_Formula*)formulaeToConvert(fi))->ConvertToSimple (varList);
        }

        for (unsigned long vi = 0; vi < varListAux.countitems(); vi++) {
            avlList.Insert ((BaseRef)varListAux.lData[vi], vi);
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

BaseRef  _ExecutionList::toStr (unsigned long) {
    _StringBuffer *result = new _StringBuffer (256);
    _String step ("\n\nStep ");

    _ExecutionList* stash = currentExecutionList;

    currentExecutionList = this;

    for (unsigned long i=0UL; i<countitems(); i++) {
        (*result) << &step << _String((long)i) << '.';
        result->AppendNewInstance ((_String*)GetItem(i)->toStr());
    }

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

void     _ExecutionList::SetNameSpace (_String const& nID) {
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

_String const _ExecutionList::AddNameSpaceToID (_String const & theID, _String const * extra) {
    _String name_space;

    if (extra && extra->nonempty()) {
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

_String  _ExecutionList::TrimNameSpaceFromID (_String& theID) {
    if (nameSpacePrefix) {
        _String * prefix_name = nameSpacePrefix->GetName();
        if (theID.BeginsWith(*prefix_name)) {
            return theID.Cut(prefix_name->length () +1,-1);
        }
    }
    return theID;
}

  //____________________________________________________________________________________

void  _ExecutionList::BuildChoiceList (_List * pieces, long code) {
  _ElementaryCommand * choice_list = new _ElementaryCommand (code);
  
  choice_list->parameters << pieces->GetItem(0L)
                          << pieces->GetItem(1L)
                          << pieces->GetItem(2L)
                          << pieces->GetItem(3L);
  
  if (pieces->countitems() > 5UL) { // explicit list of choices
    
    if (pieces->countitems() % 2 > 0) {
      DeleteObject(choice_list);
      throw ("Must have an even number of arguments for explicitly enumerated choice - description pairs");
    }
    
    _List * choices = new _List;
    
    for (unsigned long k = 4UL; k < pieces->countitems(); k+=2) {
        _String * selector    = new _String (*(_String*)pieces->GetItem(k)),
                * description  = new _String (*(_String*)pieces->GetItem(k+1));
        selector->StripQuotes(); description->StripQuotes();
        choices < new _List (selector,description);
    }
    
    choice_list->simpleParameters << 0L;
  } else {
    choice_list->parameters << pieces->GetItem(4L);
    choice_list->simpleParameters << 1L;
  }
  
  choice_list->addAndClean (*this);

}

  //____________________________________________________________________________________

void  _ExecutionList::BuildExecuteCommandInstruction (_List * pieces, long code) {

  _ElementaryCommand * run_source = new _ElementaryCommand (code);
  run_source->parameters<<pieces->GetItem(0);
  
  if (PeekFilePath()) {
    run_source->parameters && *PeekFilePath();
  } else {
    run_source->parameters.AppendNewInstance(new _String);
  }
  
  if (pieces->countitems() > 1UL) {
    if (*(_String*)pieces->GetItem (1UL) == _String("compiled")) {
      run_source->simpleParameters << 1;
    } else {
      if (*(_String*)pieces->GetItem (1UL) == _String("enclosing_namespace")) {
        run_source->parameters.Delete(1UL);
        run_source->parameters < new _String;
      } else {
        run_source->parameters << pieces->GetItem(1UL);
      }
    }
    
    if (pieces->countitems () > 2UL) {
      run_source->parameters << pieces->GetItem(2UL);
    }
  }
  run_source->addAndClean (*this);

}

  //____________________________________________________________________________________

void  _ExecutionList::BuildFscanf(_List * pieces, long code) {  

    static _String kFscanfRewind ("REWIND");

    long    names_vs_types_offset = 0L;
  
    _List   local_object_manager;


    _ElementaryCommand * scanf = new _ElementaryCommand (code);
    scanf->parameters << pieces->GetItem(0);

    bool                 has_rewind = *(_String*) pieces->GetItem (1) == kFscanfRewind;


      // process argument types

    local_object_manager < scanf;
    _List     argument_types;
    _String*  argument_type_spec = (_String*)pieces->GetItem(has_rewind ? 2L : 1L);
    argument_type_spec->StripQuotes();
    _ElementaryCommand::ExtractConditions(*argument_type_spec, 0, argument_types, ',');

    argument_types.ForEach ([&] (BaseRefConst t) -> void {
      long argument_type = _ElementaryCommand :: fscanf_allowed_formats.FindObject(t);
      if (argument_type == kNotFound) {
        throw (((_String*)t)->Enquote() & " is not one of the supported argument types: " & _ElementaryCommand :: fscanf_allowed_formats.Join (", "));
      }
      scanf->simpleParameters << argument_type;
    }
                            );

    for (unsigned long index = has_rewind ? 3L : 2L; index < pieces->countitems(); index ++) {
      scanf->parameters << pieces->GetItem(index);
    }

    if (scanf->parameters.countitems() + 1UL != scanf->simpleParameters.countitems()) {
      throw (_String("The numbers of parameter type descriptors (")& _String((long)scanf->simpleParameters.countitems()) &") and arguments ("
             &_String((long)(scanf->parameters.countitems() - 1UL))& ") did not match");
    }

    if (has_rewind) {
      scanf->simpleParameters << -1L;
    }

    local_object_manager.Pop();
    scanf->addAndClean (*this);
}

  //____________________________________________________________________________________


bool        _ExecutionList::BuildList   (_String& s, _SimpleList* bc, bool processed, bool empty_is_success) {
    if (terminate_execution) {
        return false;
    }

  //char const * savePointer = s.get_str();

    _SimpleList          triePath;
    _List                local_object_manager;
  
    try {

      while (s.nonempty ()) { // repeat while there is stuff left in the buffer
          _String currentLine (_ElementaryCommand::FindNextCommand (s));

          if (currentLine.get_char(0)=='}') {
              currentLine.Trim(1,kStringEnd);
          }

          if (currentLine.nonempty()) {
              continue;
          }

          triePath.Clear(false);
          long prefixTreeCode = _HY_ValidHBLExpressions.FindKey (currentLine, &triePath, true);

          _List *pieces = nil;
          _HBLCommandExtras *commandExtraInfo = nil;

          if (prefixTreeCode != kNotFound) {
              prefixTreeCode = _HY_ValidHBLExpressions.GetValue(prefixTreeCode);
              long commandExtra = _HY_HBLCommandHelper.FindLong (prefixTreeCode);
              if (commandExtra >= 0) { // pre-trim all strings as needed
                  commandExtraInfo = (_HBLCommandExtras*)_HY_HBLCommandHelper.GetXtra (commandExtra);
                  if (!commandExtraInfo->extract_conditions.empty()) {
                      local_object_manager < (pieces = new _List);
                    
                      long upto = _ElementaryCommand::ExtractConditions (currentLine, commandExtraInfo->cut_string,*pieces,commandExtraInfo->extract_condition_separator),
                           condition_index_match = commandExtraInfo->extract_conditions.Find(pieces->lLength);
                      if (condition_index_match < 0) {
                          // try to see if the command accepts a variable number of arguments (at least X)
                         if (commandExtraInfo->extract_conditions.lLength == 1 && commandExtraInfo->extract_conditions.lData[0] < 0) {
                              if (pieces->lLength < -commandExtraInfo->extract_conditions.lData[0]) {
                                   throw (_String("Incorrect number of arguments (") & (long) pieces->lLength & ") supplied: expected at least " & _String (-commandExtraInfo->extract_conditions.lData[0]) & ", while processing '"& currentLine.Cut (0, upto) & "'. ");
                              } else {
                                   throw (_String("Incorrect number of arguments (") & (long) pieces->lLength & ") supplied: expected one of " & _String ((_String*)commandExtraInfo->extract_conditions.toStr()) & ", while processing '"& currentLine.Cut (0, upto) & "'. ");
                              }
  
                         }
                      }
                    
                      if (commandExtraInfo->do_trim) {
                          currentLine.Trim (upto, kStringEnd);
                      }
                  }
              }
          }

          bool handled = true;

          switch (prefixTreeCode) {
              case HY_HBL_COMMAND_FOR:
                  _ElementaryCommand::BuildFor (currentLine, *this, pieces);
                 break;
              case HY_HBL_COMMAND_WHILE:
                  _ElementaryCommand::BuildWhile (currentLine, *this, pieces);
                  break;
              case HY_HBL_COMMAND_BREAK:
              case HY_HBL_COMMAND_CONTINUE:
                  if (bc) {
                      AppendNewInstance(new _ElementaryCommand);
                      (*bc) << ((prefixTreeCode == HY_HBL_COMMAND_BREAK) ? (countitems()-1) : (-(long)countitems()+1));
                      currentLine = kEmptyString;
                  } else {
                      throw (currentLine.Enquote() & " only makes sense in the context of a loop.");
                   }
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
              case HY_HBL_COMMAND_GET_DATA_INFO:
              case HY_HBL_COMMAND_GET_INFORMATION:
              case HY_HBL_COMMAND_REPLICATE_CONSTRAINT:
              case HY_HBL_COMMAND_MPI_SEND:
              case HY_HBL_COMMAND_MPI_RECEIVE:
              case HY_HBL_COMMAND_FIND_ROOT:
              case HY_HBL_COMMAND_ALIGN_SEQUENCES:
              case HY_HBL_COMMAND_DO_SQL:
              {
                    _ElementaryCommand::ExtractValidateAddHBLCommand (currentLine, prefixTreeCode, pieces, commandExtraInfo, *this);
                    break;
              }
              case HY_HBL_COMMAND_EXECUTE_A_FILE:
              case HY_HBL_COMMAND_EXECUTE_COMMANDS:
              case HY_HBL_COMMAND_LOAD_FUNCTION_LIBRARY: {
                  BuildExecuteCommandInstruction (pieces, prefixTreeCode);
               }
              break;
              
              case HY_HBL_COMMAND_FSCANF:
              case HY_HBL_COMMAND_SSCANF: {
                BuildFscanf (pieces, prefixTreeCode);
              }
              break;
              
              case HY_HBL_COMMAND_CHOICE_LIST: {
                BuildChoiceList(pieces, HY_HBL_COMMAND_CHOICE_LIST);
              }
              break;
              
              default :
                handled = false;

          }


          // TODO 20111212: this horrendous switch statement should be replaced with a
          // prefix tree lookup
          if (handled) {
            if (currentLine.length() > 1UL) {
              throw (currentLine.Enquote() & " contained syntax errors, possibly a missing semicolon. " );
            }
          } else {
              if (currentLine.BeginsWith (blFunction)||currentLine.BeginsWith (blFFunction)||currentLine.BeginsWith (blLFunction) || currentLine.BeginsWith (blNameSpace)) { // function declaration
                  _ElementaryCommand::ConstructFunction (currentLine, *this);
              } else if (currentLine.BeginsWithAndIsNotAnIdent (blReturnPrefix)) { // function return statement
                                                                            //StringToConsole(currentLine); NLToConsole();
                  _ElementaryCommand::ConstructReturn (currentLine, *this);
              } else if (currentLine.BeginsWith (blIf)) { // if-then-else statement
                  _ElementaryCommand::BuildIfThenElse (currentLine, *this, bc);
              } else if (currentLine.BeginsWith (blElse)) { // else clause of an if-then-else statement
                  if (lastif.countitems()) {
                      long    temp = countitems(),
                              lc   = lastif.countitems(),
                              lif  = lastif.lData[lc-1];

                      _ElementaryCommand      * stuff = new _ElementaryCommand ();
                      stuff->MakeJumpCommand  (nil,0,0,*this);
                      AppendNewInstance       (stuff);
                      currentLine.Trim        (4,-1);

                      long  index         = currentLine.length ()-1L,
                            scopeIn     = 0;

                      while (currentLine.char_at (scopeIn) =='{' && currentLine.char_at (index)=='}') {
                          scopeIn++;
                          index--;
                      }

                      if (scopeIn) {
                          currentLine.Trim (scopeIn,index);
                      }

                      BuildList (currentLine,bc,true);

                      if (lif<0 || lif>=lLength) {
                          throw ("'else' w/o an if to latch on to...");
                      }

                    
                      ((_ElementaryCommand*)((*this)(lif)))->MakeJumpCommand(nil,-1,temp+1,*this);
                      ((_ElementaryCommand*)(*this)(temp))->simpleParameters[0]=countitems();

                      while (lastif.countitems()>=lc) {
                          lastif.Delete(lastif.countitems()-1);
                      }
                  } else {
                      throw ("'else' w/o an if to latch on to...");
                  }

              } else if (currentLine.BeginsWith (blDo)) { // do {} while statement
                  _ElementaryCommand::BuildDoWhile (currentLine, *this);
              } else if (currentLine.BeginsWith (blInclude)) { // #include
                  _ElementaryCommand::ProcessInclude (currentLine, *this);
              } else if (currentLine.BeginsWith (blDataSet)) { // data set definition
                  _ElementaryCommand::ConstructDataSet (currentLine, *this);
              } else if (currentLine.BeginsWith (blDataSetFilter)) { // data set filter definition
                  _ElementaryCommand::ConstructDataSetFilter (currentLine, *this);
              } else if (currentLine.BeginsWith (blTree) || currentLine.BeginsWith (blTopology)) { // tree definition
                  _ElementaryCommand::ConstructTree (currentLine, *this);
              } else if (currentLine.BeginsWith (blLF) || currentLine.BeginsWith (blLF3)) { // LF definition
                  _ElementaryCommand::ConstructLF (currentLine, *this);
              }  else if (currentLine.BeginsWith (blCategory)) { // category variable declaration
                  _ElementaryCommand::ConstructCategory (currentLine, *this);
              } else if (currentLine.BeginsWith (blModel)) { // Model declaration
                  _ElementaryCommand::ConstructModel (currentLine, *this);
              } else if (currentLine.BeginsWith (blChoiceList)) { // choice list
                  _ElementaryCommand::ConstructChoiceList (currentLine, *this);
              } else if (currentLine.BeginsWith (blHBLProfile)) { // #profile
                  _ElementaryCommand::ConstructProfileStatement (currentLine, *this);
              } else if (currentLine.BeginsWith (blSCFG)) { // SCFG definition
                  _ElementaryCommand::ConstructSCFG (currentLine, *this);
              } else if (currentLine.BeginsWith (blBGM)) {    // Bayesian Graphical Model definition
                  _ElementaryCommand::ConstructBGM (currentLine, *this);
              }
              // plain ol' formula - parse it as such!
              else {
                  _String checker (currentLine);
                  if (_ElementaryCommand::FindNextCommand (checker).length ()==currentLine.length()) {
                      if (currentLine.length()>1)
                          while (currentLine (-1L) ==';') {
                              currentLine.Trim (0,currentLine.length()-2);
                          }
                      else {
                          continue;
                      }
                      _ElementaryCommand* oddCommand = new _ElementaryCommand(currentLine);
                      oddCommand->code = 0;
                      oddCommand->parameters&&(&currentLine);
                      AppendNewInstance (oddCommand);
                  } else {
                      while (currentLine.nonempty()) {
                          _String part (_ElementaryCommand::FindNextCommand (currentLine));
                          BuildList (part,bc,processed);
                      }
                  }
              }
           }
      }
    } catch (_String const & error) {
      if (currentExecutionList) {
        currentExecutionList->ReportAnExecutionError(error, false, true);
      } else {
        HandleApplicationError(error);
      }
    }
  //  s.sData = savePointer;
  // TODO: SLKP 20170623 why is this here? 20170704 ; for the "soft trim" situation, which we won't be using any more
    s.Clear();
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
_ElementaryCommand::~_ElementaryCommand (void) {
    if (CanFreeMe()) {
        if (code==4) {
            if (simpleParameters.lLength>2) {
                _Formula* f = (_Formula*)simpleParameters(2);
                delete (f);
            }
        } else if (code==0) {
            if (simpleParameters.lLength) {

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
BaseRef   _ElementaryCommand::makeDynamic (void) const
{
    _ElementaryCommand * nec = new _ElementaryCommand;
    nec->code = code;
    nec->Duplicate (this);
    return nec;
}

//____________________________________________________________________________________

void      _ElementaryCommand::Duplicate (BaseRefConst source)
{
    _ElementaryCommand* sec = (_ElementaryCommand*)source;
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

BaseRef   _ElementaryCommand::toStr      (unsigned long) {

    auto parameter_to_string = [&] (unsigned long i) -> _String const {
        return _String ((_String*)GetIthParameter(i)->toStr());
    };

    auto procedure = [&] (long i) -> _String const {
        return _StringBuffer (_HY_ValidHBLExpressions.RetrieveKeyByPayload(i))
                << _String ((_String*)parameters.Join (", ")) << ");";
    };

    auto assignment = [&] (long i, const _String& call) -> _String const {
        return _StringBuffer (_HY_ValidHBLExpressions.RetrieveKeyByPayload(i))
                << parameter_to_string (0)
                << " = "
                << call
                << _String ((_String*)parameters.Join (", ", 1)).Enquote('(', ')')
                << ";";
    };

    auto hash_pragma = [&] (long i) {
        return _StringBuffer (_HY_ValidHBLExpressions.RetrieveKeyByPayload(i))
        << _String ((_String*)parameters.Join (", ")) << ";";
    };

    _StringBuffer* string_form = new _StringBuffer (256);

    switch (code) {

        case HY_HBL_COMMAND_FORMULA: // formula reparser
            (*string_form) << parameter_to_string (0) << ";";
            break;

        case 4: {
            if (simpleParameters.countitems()==3 || parameters.countitems() == 1) {
                (*string_form) << "Branch under condition "
                            << parameter_to_string (0).Enquote()
                            << "\n\tto\n\t\t"
                            << _hblCommandAccessor (currentExecutionList,simpleParameters(0))
                            << "\n\telse\n\t\t"
                            << _hblCommandAccessor (currentExecutionList,simpleParameters(1));
              } else {
                (*string_form) << "Go to step " << _String (simpleParameters(0));
            }
        }
        break;

        case 5: { // data set construction
            (*string_form) << assignment (HY_HBL_COMMAND_DATA_SET, "ReadDataFile");
        }
        break;

        case 6: { // data set filter
            (*string_form) << assignment (HY_HBL_COMMAND_DATA_SET_FILTER, "CreateFilter");
        }
        break;

        case HY_HBL_COMMAND_HARVEST_FREQUENCIES:
        case HY_HBL_COMMAND_FPRINTF:
        case HY_HBL_COMMAND_OPTIMIZE: // optimize the likelihood function
        case HY_HBL_COMMAND_COVARIANCE_MATRIX:  // compute the covariance matrix
        case HY_HBL_COMMAND_EXPORT:
        case HY_HBL_COMMAND_MOLECULAR_CLOCK:
        case HY_HBL_COMMAND_CLEAR_CONSTRAINTS:
        case HY_HBL_COMMAND_SET_DIALOG_PROMPT:
        case HY_HBL_COMMAND_USE_MODEL:
        case HY_HBL_COMMAND_GET_STRING:
        case HY_HBL_COMMAND_SET_PARAMETER:
        case HY_HBL_COMMAND_DIFFERENTIATE:
        case HY_HBL_COMMAND_LFCOMPUTE:
        case HY_HBL_COMMAND_GET_URL:
        case HY_HBL_COMMAND_DELETE_OBJECT:
        case HY_HBL_COMMAND_REQUIRE_VERSION:
        case HY_HBL_COMMAND_ASSERT:
        case HY_HBL_COMMAND_FIND_ROOT:
        case HY_HBL_COMMAND_INTEGRATE:
        case HY_HBL_COMMAND_GET_DATA_INFO:
        case HY_HBL_COMMAND_CONSTRUCT_CATEGORY_MATRIX:
        case HY_HBL_COMMAND_ALIGN_SEQUENCES:
        case HY_HBL_COMMAND_REPLICATE_CONSTRAINT:
        case HY_HBL_COMMAND_MPI_RECEIVE:
        case HY_HBL_COMMAND_MPI_SEND :
        case HY_HBL_COMMAND_EXECUTE_A_FILE :
        case HY_HBL_COMMAND_EXECUTE_COMMANDS :
        case HY_HBL_COMMAND_LOAD_FUNCTION_LIBRARY :
        case HY_HBL_COMMAND_DO_SQL:
        case HY_HBL_COMMAND_CHOICE_LIST:
        case HY_HBL_COMMAND_SIMULATE_DATA_SET: {
            (*string_form) << procedure (code);
        }

        break;

        case 7: // build a tree
        case 54: { // build a tree
            (*string_form) << assignment (code == 7 ? HY_HBL_COMMAND_TREE : HY_HBL_COMMAND_TOPOLOGY, kEmptyString);

        }
        break;


        case 11: { // build the likelihood function
            (*string_form) << assignment (HY_HBL_COMMAND_LIKELIHOOD_FUNCTION, kEmptyString);
        }
        break;

 
        case 13: { // a function
            (*string_form) << "function "
                        << parameter_to_string(0)
                        << " ( "
                        << _String ((_String*)parameters.Join (", ", 1,parameters.countitems()-2L))
                        << " ) {\n" << parameter_to_string (parameters.countitems()-1L) << "\n}";
        }
        break;

        case 14: { // return statement
            (*string_form) << "return "
                        << parameter_to_string(0)
                        << ";";
        }
        break;

        case 16: { // data set merger
            (*string_form) << assignment (HY_HBL_COMMAND_DATA_SET, labs(simpleParameters(0))==1 ? "Combine" : "Concatenate");
            /*if (simpleParameters (0)<0) {
                string_form << "(deleting arguments upon completion)";
            }*/
        }
        break;


        case 20: {// category variable construction
            (*string_form) << assignment (HY_HBL_COMMAND_CATEGORY, kEmptyString);
        }
        break;

        case 24: { // select standard model
            (*string_form) << procedure (HY_HBL_COMMAND_SELECT_TEMPLATE_MODEL);
        }
        break;

        case HY_HBL_COMMAND_FSCANF: // fscanf
        case HY_HBL_COMMAND_SSCANF: { // sscanf
            (*string_form) << (code == HY_HBL_COMMAND_FSCANF ? "fscanf" : "sscanf")
                        << parameter_to_string (0)
                        << ",\"";

            long shift = 1L;
          
            for (long p = 0; p<simpleParameters.lLength; p++) {
                long theFormat = simpleParameters(p);
                if (theFormat < 0) {
                    (*string_form) << "REWIND";
                } else {
                    (*string_form) << ((_String*)_ElementaryCommand::fscanf_allowed_formats.GetItem (theFormat))->Enquote('(',')');
                }
                if (p) {
                    (*string_form) << ", ";
                }
            }
            (*string_form) << "\",";
            for (long p = 0; p<simpleParameters.lLength; p++) {
                long theFormat = simpleParameters(p);
                if (theFormat < 0) {
                    shift++;
                } else {
                    (*string_form) << parameter_to_string (p+shift);
                }
                if (p) {
                    (*string_form) << ", ";
                }
            }
            (*string_form) << ");";
        }
        break;


        case 31: { // define a model
            (*string_form) << procedure (HY_HBL_COMMAND_MODEL);
        }
        break;


        case 38: { // reconstruct ancestors
            (*string_form) << assignment (HY_HBL_COMMAND_DATA_SET, "ReconstructAncestors");

        }
        break;

        case 47: { //GetDataInfo
            (*string_form) << procedure (HY_HBL_COMMAND_STATE_COUNTER);
        }
        break;

        case 52: { //Simulate
            (*string_form) << assignment (HY_HBL_COMMAND_DATA_SET, "Simulate");
        }
        break;


        case 58: {
            (*string_form) << hash_pragma (HY_HBL_COMMAND_PROFILE);
        }
        break;

        case 61: {
            (*string_form) << assignment (HY_HBL_COMMAND_SCFG, kEmptyString);
        }
        break;

        case 64: {
            (*string_form) << assignment (HY_HBL_COMMAND_BGM, kEmptyString);
        }
        break;

        case HY_HBL_COMMAND_NESTED_LIST: {
            (*string_form) << "namespace " << parameter_to_string (0) << ";";
            break;
        }

    }
    return string_form;
}

//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase0 (_ExecutionList& chain) {
    chain.currentCommand++;

    _String * errMsg = nil;

    try {

      if (chain.cli) {
          hyFloat result = ((_Formula*)simpleParameters.lData[1])->ComputeSimple (chain.cli->stack, chain.cli->values);
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
                  simpleParameters    <<parseCode
                                      <<long (f.makeDynamic())
                                      <<long (f2.makeDynamic())
                                      <<fpc.assignmentRefID   ()
                                      <<fpc.assignmentRefType ();

                  appendCompiledFormulae (&f, &f2);

              } else {
                  ExecuteFormula(&f,&f2,parseCode,fpc.assignmentRefID(),chain.nameSpacePrefix,fpc.assignmentRefType());
                  if (terminate_execution) {
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

      if (terminate_execution) {
        errMsg = new _String ("Error computing the interpreted statement: ");
        throw 0;
      }

    } catch (int e) {
      if (errMsg) {
        HandleApplicationError (_String(errMsg) & *this);
      }
    }
}


//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase4 (_ExecutionList& chain) {
    chain.currentCommand++;

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

            if (terminate_execution && !result) {
                  _String       *s = (_String*)((_Formula*)simpleParameters(2))->toStr(kFormulaStringConversionSubstituteValues);
                  errMsg  = new _String(_String("Failed while evaluating: ") & _String((_String*)((_Formula*)simpleParameters(2))->toStr(kFormulaStringConversionNormal)) & " which expanded to  " & s);
                  throw (1);
               }

              bool conditionFalse = false;

              switch (result->ObjectClass()) {
                case NUMBER:
                    conditionFalse = result->Value()==0.0;
                    break;
                case STRING:
                    conditionFalse = ((_FString*)result)->empty();
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
          terminate_execution   = true;
          chain.currentCommand = chain.lLength;
      }
    }
    catch (int e) {
      if (expression) {
        delete expression;
      }
      if (errMsg) {
        if (e == 0) {
          HandleApplicationError (_String ("'") & *(_String*)parameters(0) & "'" & errMsg);
        } else {
          HandleApplicationError    (errMsg);
        }
        // note that errMsg will be deleted by _String (*_String) constructors
      }
    }
}

//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase5 (_ExecutionList& chain) {
    chain.currentCommand++;
    FILE*    df;
    _String  fName (*GetIthParameter(1));
    _DataSet*ds;



    if (simpleParameters.lLength == 1) {
        fName = GetStringFromFormula ((_String*)parameters(1),chain.nameSpacePrefix);
        ds = ReadDataSetFile (nil,0,&fName,nil,chain.nameSpacePrefix?chain.nameSpacePrefix->GetName():nil);
    } else {
        if (fName.Equal (&useNexusFileData)) {
            if (!lastNexusDataMatrix) {
                HandleApplicationError (useNexusFileData & " was used in ReadDataFile, and no NEXUS data matrix was available.");
                return;
            }
            ds = lastNexusDataMatrix;
        } else {
            ProcessFileName(fName, false,true,(hyPointer)chain.nameSpacePrefix);
            if (terminate_execution) {
                return;
            }
            SetStatusLine ("Loading Data");

            df = doFileOpen (fName.get_str(),"rb");
            if (df==nil) {
                // try reading this file as a string formula
                fName = GetStringFromFormula ((_String*)parameters(1),chain.nameSpacePrefix);
                ProcessFileName(fName, false,false,(hyPointer)chain.nameSpacePrefix);

                if (terminate_execution) {
                    return;
                }

                df = doFileOpen (fName.get_str(),"rb");
                if (df==nil) {
                     HandleApplicationError ((_String ("Could not find source dataset file ") & ((_String*)parameters(1))->Enquote('"')
                                & " (resolved to '" & fName & "')\nPath stack:\n\t" & GetPathStack ("\n\t")));
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
        HandleApplicationError    ("The format of the sequence file has not been recognized and may be invalid");
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
            HandleApplicationError (_String("Not a valid string matrix object passed to a _LikelihoodFunction constructor: ") & *(_String*)parameters(1));
            return;
        }
    }

    long i       = 0,
         stepper = explicitFreqs?3:2;

    for (; i<=(long)likelihoodFunctionSpec->lLength-stepper; i+=stepper) {
        _String     *dataset = (_String*)(*likelihoodFunctionSpec)(i),
                     *tree   = (_String*)(*likelihoodFunctionSpec)(i+1),
                      *freq    = explicitFreqs?(_String*)(*likelihoodFunctionSpec)(i+2):nil;

        if(GetDataFilter (AppendContainerName(*dataset,chain.nameSpacePrefix))) {
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

        if (errMsg.nonempty()) {
            DeleteObject (likelihoodFunctionSpec);
            HandleApplicationError    (errMsg);
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
            likeFuncObjectID = likeFuncNamesList.FindObject(&kEmptyString);
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

    if (tempString.nonempty()) {
        likefID = tempString;
    }

    long f  = FindLikeFuncName (likefID),
         s2 = FindSCFGName     (likefID);

    if (f==-1 && s2==-1) {
        HandleApplicationError (_String("Likelihood Function (or SCFG)")&likefID& " has not been initialized" );
        return ;
    }

    if (f>=0) {
        _DataSet  * ds = new _DataSet;

        _List     theExclusions;

        if (parameters.lLength>2) // there is a list of exclusions there
            // ';'-sep for different partitions
            // ','-sep for states in a given partition
        {
            // SLKP mod 20070622 to allow string expressions as well
            _String theExc (ProcessLiteralArgument((_String*)parameters(2),chain.nameSpacePrefix));
            if (theExc.nonempty()) {
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
                catValues = new _Matrix (1,1,false,true);
            }
        }

        if (parameters.lLength>4) {
            // a matrix to store simulated category values
            _String  matrixName (chain.AddNameSpaceToID(*(_String*)parameters(4)));
            if (!(catNameVar = CheckReceptacle(&matrixName,blSimulateDataSet,true))) {
                return;
            } else {
                catNames = new _Matrix (1,1,false,true);
            }
        }

        _String * resultingDSName = new _String (chain.AddNameSpaceToID(*(_String*)parameters(0)));

        if (!resultingDSName->IsValidIdentifier(fIDAllowCompound)) {
            errMsg = *resultingDSName & " is not a valid receptacle identifier in call to " & blSimulateDataSet;
            DeleteObject (resultingDSName);
            HandleApplicationError (errMsg);
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

void      _ElementaryCommand::ExecuteCase47 (_ExecutionList& chain) {
    chain.currentCommand++;

    _String *arg1 = GetIthParameter(0),
            *arg2 = GetIthParameter(1);

    try {
        long type = HY_BL_LIKELIHOOD_FUNCTION;
        _LikelihoodFunction const * lf = (_LikelihoodFunction const *)_HYRetrieveBLObjectByName(AppendContainerName(*arg1, chain.nameSpacePrefix), type );
        if (lf) {
            type = HY_BL_HBL_FUNCTION;
            long function_index;
            if (_HYRetrieveBLObjectByName(ProcessLiteralArgument (arg2,chain.nameSpacePrefix), type, &function_index)) {
                if (GetBFFunctionArgumentCount(function_index)!=2L) {
                    throw (arg2->Enquote() & " callback function must depend on 2 parameters ");
                } else {
                    lf->StateCounter (function_index);
                }
            } else {
                throw (arg2->Enquote() & " is not a defined user batch language function ");
            }

        } else {
            throw (arg1->Enquote() & " is not a defined likelihood function ID ");
        }

    } catch (const _String & err) {
        HandleApplicationError (err);
    }

}



//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase52 (_ExecutionList& chain) {
    chain.currentCommand++;

    long           site_count  = ProcessNumericArgument ((_String*)parameters (4),chain.nameSpacePrefix);
    _String        given_state;

    if (site_count < 1L) {
        given_state = ProcessLiteralArgument((_String*)parameters (4),chain.nameSpacePrefix);
        site_count = given_state.length();
    }

    if (site_count < 1) {
        HandleApplicationError (*(_String*)parameters (4) & " must either evaluate to a positive integer or be a non-empty string of root states");
        return;
    }

    _Variable   *  alphabet = FetchVar (LocateVarByName (AppendContainerName(*GetIthParameter(3),chain.nameSpacePrefix)), MATRIX),
                *  tree_var  = FetchVar (LocateVarByName (AppendContainerName(*GetIthParameter(1),chain.nameSpacePrefix)), TREE),
                *  freq_var  = FetchVar (LocateVarByName (AppendContainerName(*GetIthParameter(2),chain.nameSpacePrefix)), MATRIX);


    try {
        if (!alphabet) {
            throw (GetIthParameter(3)->Enquote() & " must be a defined matrix-valued variable");
        }
        if (!freq_var) {
            throw (GetIthParameter(2)->Enquote() & " must be a defined matrix-valued variable");
        }
        if (!tree_var) {
            throw (GetIthParameter(1)->Enquote() & " must be a defined tree-valued variable");
        }

        _Matrix * alphabet_matrix = (_Matrix*)alphabet->GetValue();

        if (!(alphabet_matrix->IsAStringMatrix() && alphabet_matrix->GetHDim() == 2 && alphabet_matrix->GetVDim () > 1)) {
            throw (_String("Alphabet specification variable ") & GetIthParameter(3)->Enquote() & " must be a string matrix with 2 rows and at least 2 columns");
        }

        _String base_set;

        for (unsigned long k=0UL; k < alphabet_matrix->GetVDim (); k++) {
            _FString * a_state = (_FString*)alphabet_matrix->GetFormula(0,k)->Compute();
            if (a_state) {
                if (a_state->get_str().length() == 1UL) {
                    char c = a_state->get_str().char_at(0UL);
                    if (base_set.Find(c) == -1) {
                        base_set = base_set & c;
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

        if (base_set.length() < alphabet_matrix->GetVDim ()) {
            throw ("The alphabet is mis-specified; it either has redundant characters or multi-character/non-string entries");
        }

        long unit_size = ((_FString*)alphabet_matrix->GetFormula(1,0)->Compute())->get_str().to_long();

        if (unit_size < 1L) {
            throw ("The evolutionary unit size in the alphabet matrix is mis-specified");
        }

        _Formula* exclusion_formula = alphabet_matrix->GetFormula(1,1);
        _String const* the_exclusions = &kEmptyString;

        if (exclusion_formula) {
            the_exclusions = &((_FString*)exclusion_formula->Compute())->get_str();
        }

        _TheTree * spawning_tree = (_TheTree*)tree_var;

        if (parameters.lLength>6 && (spawning_tree->CountTreeCategories()>1)) {
            throw ("Can't use spool to file option in Simulate when the tree depends on category variables.");
        }

        if (given_state.length()>1) {
        // root state
            if ((given_state.length() >= unit_size)&&(given_state.length() % unit_size == 0)) {
                site_count = given_state.length()/unit_size;
            } else {
                throw ("Root state string is either too short or has length which is not divisible by the unit size");
            }
        }

        _TranslationTable newTT (base_set);
        _DataSet * ds = new _DataSet;

        if (! newTT.IsStandardNucleotide() ) {
            ds->SetTranslationTable (&newTT);    // mod 20060113 to properly deal with non-standard alphabets
        }
        // make a dummy
        spawning_tree->AddNodeNamesToDS (ds,true,false,1);

        char    c = base_set.char_at (0);
        long    s = ds->GetNames().countitems();

        if (s<2L) {
            ds->InsertName (_String ("Root"),0L);
            s ++;
        }


        ds->AddSite(c);
        for (long u = 1L; u < s; u++) {
            ds->Write2Site(0,c);
        }
        ds->Finalize();
        ds->SetNoSpecies (s);

        unsigned long total_sites = site_count*unit_size;

        _SimpleList * the_map = & ds->GetTheMap();
        the_map->RequestSpace (total_sites);
        InitializeArray (the_map->lData, total_sites, 0L);
        the_map->lLength = total_sites;

        _DataSetFilter* new_filter = new _DataSetFilter();
        _SimpleList     h,v;

        new_filter->SetFilter     (ds,unit_size,h,v,false);
        new_filter->SetExclusions (*the_exclusions,true);
        new_filter->SetupConversion ();

        _Matrix*   root_states = nil;
        if (given_state.length()>=unit_size) {
            root_states                = new _Matrix (1,site_count,false,true);
            hyFloat*  holder         = new hyFloat [new_filter->GetDimension(false)];

            for (long cc = 0; cc < site_count; cc++) {
                unsigned long site_idex = cc*unit_size;
                _String root_char (given_state,site_idex,site_idex+unit_size-1L);
                long    root_state = new_filter->Translate2Frequencies (root_char,holder,false);
                if (root_state<0) {
                    throw (root_char & " found in the root state string at position " & _String ((long)site_idex) & " is an invalid/ambiguous state");
                    break;
                } else {
                    root_states->theData[cc] = root_state;
                }
            }
            delete [] holder;
        }


        long       filter_id = StoreDataFilter (simulationFilter, new_filter);

        spawning_tree->SetUp();
        spawning_tree->InitializeTreeFrequencies((_Matrix*)freq_var->Compute(),true);

        _String filter_specification = *GetFilterName (filter_id) & spawning_tree->GetName()->Enquote(',') & *freq_var->GetName();

        _LikelihoodFunction lf (filter_specification, nil);

        if (terminate_execution) {
            return;
        }

        bool    do_internals = parameters.countitems() > 5 ? (ProcessNumericArgument ((_String*)parameters (5),chain.nameSpacePrefix)>0.5) : nil;

        _String spool_file;

        FILE*   main_file = nil;

        if (parameters.countitems () > 6) {
            spool_file = ProcessLiteralArgument (GetIthParameter(6),chain.nameSpacePrefix);
            ProcessFileName(spool_file);
            main_file = doFileOpen (spool_file,"w");
            if (!main_file) {
                throw (_String("Failed to open ") & spool_file.Enquote() & " for writing");
            }
            if (do_internals) {
                spool_file = spool_file & ".anc";
            }
        }

        _DataSet    * sim_dataset;

        if (main_file) {
            sim_dataset = new _DataSet (main_file);
        } else {
            sim_dataset = new _DataSet (site_count);
        }

        _List exclusions;

        _String *sim_name = new _String(AppendContainerName(*GetIthParameter(0),chain.nameSpacePrefix));


        _String    rate_matrix_name       = *sim_name & ".rates";
        _Variable *category_values_id     = CheckReceptacle(&rate_matrix_name, __PRETTY_FUNCTION__);
        _Matrix*   category_values        = new _Matrix (1,1,false,true);

        _String    rate_variable_names         = *sim_name & ".rateVars";
        _Variable * category_names_id          = CheckReceptacle(&rate_variable_names, __PRETTY_FUNCTION__);
        _Matrix*    category_names       = new _Matrix (1,1,false,true);

        SetStatusLine ("Simulating Data");
        lf.Simulate (*sim_dataset, exclusions, category_values, category_names, root_states, do_internals?(main_file?&spool_file:&kEmptyString):nil);
        SetStatusLine ("Idle");

        category_values_id->SetValue(category_values, false);
        category_names_id->SetValue(category_names, false);

        StoreADataSet (sim_dataset, sim_name);
        DeleteObject (sim_name);
        DeleteDataFilter (filter_id);

        DeleteObject   (ds);
        DeleteObject (root_states);

    } catch (const _String & err) {
        HandleApplicationError(err & " in Simulate.");
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


    case 5: // data set construction

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

        hyFloat rtv = 0.0; // mod 11/19/2003
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


        _TheTree * tr = nil;

        if (treeString.get_char(0)!='(') {
            _Formula  nameForm (treeString,chain.nameSpacePrefix);
            _PMathObj formRes = nameForm.Compute();
            if (formRes) {
                if (formRes->ObjectClass () == STRING) {
                    tr = new _TheTree (treeIdent,((_FString*)formRes)->get_str(),false);
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
            HandleApplicationError("Illegal right hand side in call to Tree id = ...; it must be a string, a Newick tree spec or a topology");
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

    case 12: // data set construction by simulation

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
            HandleApplicationError (errMsg);
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
                HandleApplicationError (((_String)("Identifier ")&dsname&_String(" doesn't correspond to a valid dataset.")));
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
        _String  fName (*GetIthParameter(1));
        ProcessFileName(fName);
        if (terminate_execution) {
            return false;
        }
        FILE*   theDump = doFileOpen (fName,"rb");
        if (!theDump) {
            HandleApplicationError (((_String)("File ")&fName&_String(" couldn't be open for reading.")));
            return false;
        }

        fName = chain.AddNameSpaceToID(*(_String*)parameters(0));
        _Variable * result  = CheckReceptacle(&fName,blImport.Cut(0,blImport.length()-2),true);
        if (result) {
            _Matrix   * storage = new _Matrix (1,1,false,true);
            result->SetValue(storage,false);
            lastMatrixDeclared = result->GetIndex();
            if (!storage->ImportMatrixExp(theDump)) {
                HandleApplicationError("Matrix import failed - the file has an invalid format.");
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

    case HY_HBL_COMMAND_CONSTRUCT_CATEGORY_MATRIX: // construct the category matrix
        HandleConstructCategoryMatrix (chain);
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

    case HY_HBL_COMMAND_FSCANF: // fscanf
      HandleFscanf (chain,false);
      break;
      
    case HY_HBL_COMMAND_SSCANF: // fscanf
      HandleFscanf (chain,true);
      break;


    case HY_HBL_COMMAND_USE_MODEL:
        return HandleUseModel(chain);

    case 31:
        ExecuteCase31 (chain);
        break;

    case HY_HBL_COMMAND_CHOICE_LIST:
        HandleChoiceList (chain);
        break;

    case HY_HBL_COMMAND_GET_STRING:
        HandleGetString (chain);
        break;

    case HY_HBL_COMMAND_SET_PARAMETER:
        return HandleSetParameter(chain);

    case 38:
        ExecuteCase38 (chain, false);
        break;

    case HY_HBL_COMMAND_EXECUTE_COMMANDS:
        HandleExecuteCommandsCases(chain, false, false);
        break;

    case HY_HBL_COMMAND_EXECUTE_A_FILE:
        HandleExecuteCommandsCases(chain, true, false);
        break;

    case HY_HBL_COMMAND_LOAD_FUNCTION_LIBRARY:
        HandleExecuteCommandsCases(chain, true, true);
        break;


    case HY_HBL_COMMAND_DIFFERENTIATE:
        return HandleDifferentiate (chain);
        break;

    case HY_HBL_COMMAND_INTEGRATE:
    case HY_HBL_COMMAND_FIND_ROOT:
        return HandleFindRootOrIntegrate(chain);
        break;


    case HY_HBL_COMMAND_MPI_SEND:
        return HandleMPISend (chain);
        break;

    case HY_HBL_COMMAND_MPI_RECEIVE:
        HandleMPIReceive(chain);
        break;


    case 47:
        ExecuteCase47 (chain);
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

    case HY_HBL_COMMAND_DO_SQL:
        return HandleDoSQL(chain);

    case 54:
        ExecuteCase54 (chain);
        break;

    case HY_HBL_COMMAND_ALIGN_SEQUENCES:
        return HandleAlignSequences (chain);
        break;


    case 58:
        ExecuteCase58 (chain);
        break;

    case HY_HBL_COMMAND_DELETE_OBJECT:
        return HandleDeleteObject (chain);
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


const _String   _ElementaryCommand::FindNextCommand  (_String& input) {

    long    index     = input.length();

    if (index == 0L) {
      return kEmptyString;
    }

    bool    skipping  = false;

    enum {
        normal_text = 0,
        double_quote = 1,
        single_quote = 2

    } literal_state = normal_text;

    enum    {
        no_comment = 0,
        slash_star = 1,
        double_slash = 2
    } comment_state = no_comment;


    long    scope_depth     = 0L, // count how deep we are in {}
            matrix_depth = 0L,  // count how deep we are in {} matrix definition (different from "scope")
            parentheses_depth     = 0L,   // count how deep we are in ()
            bracket_depth   = 0L;   // count how deep we are in []

    _SimpleList is_DoWhileLoop;

    _StringBuffer result (128L);

    char    last_char = '\0';
        // a look back character


    // non printable characters at the end ?
    while (index>=0 && !isprint(input[--index])) ;
    input.Trim (0,index);

    for (index = 0L; index<input.length(); index++) {
        char c = input.char_at (index);

        if (literal_state == normal_text && c=='\t') {
            c = ' ';
        }

        // check for comments
        if (comment_state != no_comment) {
            if (comment_state == slash_star) {
                if (c=='/' && input.get_char(index-1)=='*') {
                    comment_state = no_comment;
                }
            } else if (c == '\r' || c == '\n') {
                comment_state = no_comment;
            }

            last_char  = '\0';
            continue;
        } else {
            if (literal_state == normal_text && c=='/') {
                switch (input.get_char(index+1)) {
                case '*':
                    comment_state = slash_star;
                    break;
                case '/':
                    comment_state = double_slash;
                }

                if (comment_state != no_comment) {
                    last_char  = '\0';
                    index++;
                    continue;
                }
            }
        }


        // skip spaces, except for special cases, like return and data set filters

        if (literal_state == normal_text && isspace(c)) {

          // skip/compress spaces, unless we are in a higher level HBL statement
          // where spaces can't be compressed
          // examples include
          // DataSet|DataSetFilter|return|LikelihoodFunction (something)
          // need to maintain spaces for this to work appropriately



            if (!skipping && index > 0L) {

              long trie_match = _HY_HBL_KeywordsPreserveSpaces.FindKey(input.Cut (MAX (0, index - 20), index-1).Reverse(), nil, true);
              if (trie_match != kNotFound) {
                long matched_length = _HY_HBL_KeywordsPreserveSpaces.GetValue(trie_match);
                if (matched_length == index || !(isalnum(input.get_char(index-matched_length-1)) || input.get_char(index-matched_length-1) == '_' || input.get_char(index-matched_length-1) == '.')) {
                  result << ' ';
                }
              }
            }


            skipping = true;
            continue;
        }

        if (skipping && ( isalpha(c) || c=='_') && (isalnum(last_char) || last_char =='_')) {
          // SLKP 20170704: this seems incomplete : need to check more thoroughly that this is an ident
          // this is meant to determine that we are at the beginning of a new ident-like
          // token and insert a space
            result<<' ';
        }

        skipping = false;

        result<<c;

        if (literal_state != normal_text && c == '\\') {
            // escape character \x
            result<< input.get_char(++index);
            continue;
        }

        // are we inside a string literal?

        if (c=='"') {
          if (literal_state != single_quote) {
            literal_state = literal_state == normal_text ? double_quote : normal_text;
            last_char = '\0';
            continue;
          }
        } else {
          if (c == '\'') {
              if (literal_state != double_quote) {
              literal_state = literal_state == normal_text ? single_quote : normal_text;
              last_char = '\0';
              continue;
              }
          }
        }

        if (literal_state != normal_text) {
            continue;
        }

        // maybe we are done?

        if (c==';' && scope_depth == 0L && matrix_depth == 0L && parentheses_depth == 0L && bracket_depth == 0L) {
            // SLKP 20170704 used to be parentheses_depth <= 0L && bracket_depth <= 0L
            break;
        }

        // check to see whether we are defining a matrix

        if (c=='(') {
            parentheses_depth ++;
            last_char = '\0';
            continue;
        }

        if (c==')') {
            parentheses_depth --;
            if (parentheses_depth < 0L) {
                HandleApplicationError (_String("Too many closing ')' near '") & input.Cut (MAX(0,index-32),index) & "'.");
                input.Clear();
                return kEmptyString;
            }
            last_char = '\0';
            continue;
        }

        if (c=='[') {
            bracket_depth++;
            last_char = '\0';
            continue;
        }

        if (c==']') {
            bracket_depth--;
            if (bracket_depth < 0L) {
                HandleApplicationError (_String("Too many closing ']' near '") & input.Cut (MAX(0,index-32),index) & "'.");
                input.Clear();
                return kEmptyString;
            }
            last_char = '\0';
            continue;
        }


        if (c=='{') {
            if (matrix_depth) {
                matrix_depth++;
            } else if (last_char == '=') { // a matrix def
                matrix_depth++;
            } else {
                scope_depth++;
                if (index>=2L) {
                    long t = input.FirstNonSpaceIndex (0, index-1, kStringDirectionBackward);
                    if (t>=1) {
                        if (input.get_char(t)=='o' && input.get_char(t-1)=='d') {
                            is_DoWhileLoop << scope_depth-1L;
                            //printf ("%d\n%s\n\n", isDoWhileLoop, input.Cut (t,-1).sData);
                        }
                    }
                }
            }
            last_char = '\0';
            continue;
        }

        if (c=='}') {
            if (matrix_depth) {
                matrix_depth--;
            } else {
                scope_depth--;
                if (parentheses_depth == 0L && bracket_depth == 0L) {
                    if (scope_depth >=0 && is_DoWhileLoop.lLength && is_DoWhileLoop.GetElement(-1L) == scope_depth) {
                        is_DoWhileLoop.Pop();
                    } else if (scope_depth == 0L) {
                        break;
                    }
                }

            }
            last_char = '\0';
            continue;
        }

        last_char = c;
    }


    if (scope_depth != 0L || comment_state != slash_star || literal_state != normal_text || matrix_depth != 0L || bracket_depth != 0L || parentheses_depth != 0L) {
        if (result!='}') {
            HandleApplicationError (_String("Expression appears to be incomplete/syntax error. {} scope: ") &scope_depth & ", () depth "
                       & parentheses_depth & ", matrix scope: " & matrix_depth & '.' & (literal_state == double_quote ?" In a \"\" literal. ":kEmptyString)
                       & (literal_state == single_quote?" In a '' literal. ":kEmptyString) &
                       (comment_state == slash_star ? " In a /* */ comment ":kEmptyString) & '\n' & input);
            input.Clear();
            return kEmptyString;
        } else {
            result = kEmptyString;
        }
    }

    long check_open = 0L;
    while (result.get_char(check_open)=='{') {
        check_open++;
    }

    if (check_open) {
        long index2 = result.length() - 1L;

        while (result[index2]=='}') {
            index2--;
        }

        if (result.length () - index2 - 1 < check_open) {
            HandleApplicationError ((_String)("Expression appears to be incomplete/syntax error and will be ignored:")&input);
            result.Clear ();
        } else {
            result.Trim(check_open,result.length()-1-check_open);
        }
    }

    if (index<input.length()-1) {
        input.Trim (index+1L, kStringEnd);
    } else {
        input.Clear();
    }



    return result;
}
//____________________________________________________________________________________

long _ElementaryCommand::ExtractConditions (_String const& source, long start_at, _List& receptacle, char delimiter, bool include_empty_conditions) {

    long parentheses_depth = 1L,
         // this is because extraction will work from the first character following a '(', e.g. CreateFilter([start parsing here]....)
         last_delim    = start_at,
         index             = start_at,
         curly_depth       = 0L;


    enum {
        normal_text = 0,
        single_quote = 1,
        double_quote = 2
    } quote_type;


    for (; index<source.length(); index++) {
        char c = source.char_at (index);
        if (quote_type == normal_text) {
            if (c=='(') {
                parentheses_depth++;
                continue;
            }
            if (c=='{') {
                curly_depth++;
                continue;
            }
            if (c=='}') {
                curly_depth--;
                continue;
            }
            if (c==')') {
                parentheses_depth --;
                if (!parentheses_depth == 0L) {
                    break;
                }
                continue;
            }
        }
        if (c=='"' && quote_type != single_quote) {
            if (index == start_at || source.char_at (index-1L) != '\\') {
                quote_type = quote_type == normal_text ? double_quote : normal_text;
            }
            continue;
        }
        if (c=='\'' && quote_type != double_quote) {
            if (index == start_at || source.char_at (index-1L) != '\\') {
                quote_type = quote_type == normal_text ? single_quote : normal_text;
            }
            continue;
        }
        if (c==delimiter) {
            if (parentheses_depth > 1 || quote_type != normal_text || curly_depth) {
                continue;
            }

            receptacle < new _String (source,last_delim,index-1);
            last_delim = index+1;
            continue;
        }
    }

    if (include_empty_conditions || last_delim <= index-1) {
        receptacle < new _String(source,last_delim,index-1);
    }
    return index+1L;
}

//____________________________________________________________________________________


bool       _ElementaryCommand::MakeGeneralizedLoop  (_String*p1, _String*p2, _String*p3 , bool for_or_while, _String& source, _ExecutionList&target) {

    // extract the for enclosure
    long  beginning = target.lLength,
          for_return = target.lLength,
          index;

    bool   success = true;
    bool   has_increment = false;

    _SimpleList bc;

    while (success) {

        if (p1 && p1->nonempty()) { // initialization stage
            for_return++;
            success = success && target.BuildList (*p1, nil, true); // add init step
        }

        // append condition now

        if (!success) {
            break;
        }

        if (for_or_while) {
            if (p2 && p2->nonempty()) { // condition stage
                target < new _ElementaryCommand (*p2);
            }
        }

        if (source.get_char(0)=='{') {
            source.Trim(1,kStringEnd);
        }

        if ((success = success && target.BuildList (source, &bc)) == false) { // construct the main body
            break;
        }

        if (p3 && p3->nonempty ()) { // increment stage
            success = success && target.BuildList (*p3, nil,true); // add increment step
            has_increment = true;
        }

        if (!success) {
            break;
        }

        if (for_or_while) {
            _ElementaryCommand * loopback = new _ElementaryCommand;
            success = success && loopback->MakeJumpCommand (nil,for_return,0,target);
            target < loopback;
            if (p2 && p2->nonempty()) {
                success = success && (target.GetIthCommand(for_return))->MakeJumpCommand (p2, for_return+1, target.lLength,target);
            }
        } else {
            if (p2) {
                _ElementaryCommand* loopback = new _ElementaryCommand ;
                success = success && loopback->MakeJumpCommand (p2,for_return,target.lLength+1,target);
                target < loopback;
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
                (target.GetIthCommand(loc))->MakeJumpCommand (nil, target.lLength, 0,target);
            } else { // continue
                (target.GetIthCommand(-loc))->MakeJumpCommand (nil, target.lLength-(has_increment?2:1), 0,target);
            }
        }
    }

    return true;
}

//____________________________________________________________________________________


bool       _ElementaryCommand::BuildFor (_String&source, _ExecutionList&target,  _List * pieces)

/* the for loop becomes this:

    1. initialize
    2. if (condition) then
    3. execute loop body
    4. else go to 7
    5. increment statement (if present)
    6. go to 2
    7. code following the loop
*/


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
            HandleApplicationError ("'if' header makes no sense");
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
        if (clipped.BeginsWith (blWhile)) {
            source.Trim (blDo.length(),upto);
            _List pieces;
            ExtractConditions (clipped,blWhile.length(),pieces);
            if (pieces.lLength != 1) {
                HandleApplicationError ("Malformed while clause in a do-while loop");
                return false;
            }

            if (!MakeGeneralizedLoop (nil,(_String*)pieces(0),nil,false,source,target)) {
                return false;
            }

            return true;
        }
    }
    HandleApplicationError ("Could not find a matching 'while' in the definition of a do-while loop");

    return false;
}

//____________________________________________________________________________________

bool    _ElementaryCommand::ProcessInclude      (_String&source, _ExecutionList&target)
{

    _String         fileName (source, blInclude.length (),(long)source.length () - 2L);
  
    ProcessFileName(fileName, false,false,(hyPointer)target.nameSpacePrefix);
    if (terminate_execution) {
        return false;
    }

    PushFilePath  (fileName);
    ReadBatchFile (fileName, target);
    PopFilePath   ();

    return true;
}

//____________________________________________________________________________________

_ElementaryCommand* makeNewCommand (long ccode) {
    return               new _ElementaryCommand (ccode);
}

//____________________________________________________________________________________

void _ElementaryCommand::addAndClean (_ExecutionList&target,_List* parameter_list, long start_at) {
    if (parameter_list) {
        // TODO 20170913 SLKP : check that << works as well
        for (long i = start_at; i < parameter_list->countitems(); i++) {
          parameters << parameter_list->GetItem(i);
        }
    }
    target.AppendNewInstance(this);
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
        HandleErrorWhileParsing ("DataSet declaration missing a valid identifier", source);
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
            HandleErrorWhileParsing ("DataSet declaration missing a valid filename", source);
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
            HandleErrorWhileParsing (blSimulateDataSet & "expects 1-4 parameters: likelihood function ident (needed), a list of excluded states, a matrix to store random rates in, and a matrix to store the order of random rates in (last 3 - optional).",
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
            HandleErrorWhileParsing("DataSet merging operation missing a valid list of arguments.",source);
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
            HandleErrorWhileParsing("DataSet merging operation missing a valid list of arguments.",source);
            return false;
        }

        target&&(&dsc);
        return true;

    } else {
        if (oper ==  _String("ReconstructAncestors") || oper ==  _String("SampleAncestors")) {
            _List pieces;
            ExtractConditions (source,mark1+1,pieces,',');
            if (pieces.lLength>3UL || pieces.lLength==0UL) {
                HandleErrorWhileParsing("ReconstructAncestors and SampleAncestors expects 1-4 parameters: likelihood function ident (mandatory), an matrix expression to specify the list of partition(s) to reconstruct/sample from (optional), and, for ReconstructAncestors, an optional MARGINAL flag, plus an optional DOLEAVES flag.",
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
                HandleErrorWhileParsing ("Simulate expects 4-6 parameters: tree with attached models, equilibrium frequencies, character map, number of sites|root sequence, <save internal node sequences>, <file name for direct storage>",
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
            HandleErrorWhileParsing ("Expected DataSet ident = ReadDataFile(filename); or DataSet ident = SimulateDataSet (LikelihoodFunction); or DataSet ident = Combine (list of DataSets); or DataSet ident = Concatenate (list of DataSets); or DataSet ident = ReconstructAncestors (likelihood function); or DataSet ident = SampleAncestors (likelihood function) or DataSet	  dataSetid = ReadFromString (string);",
                                   source);
        }
    }

    return false;
}
//____________________________________________________________________________________

bool    _ElementaryCommand::ConstructCategory (_String&source, _ExecutionList&target)
// category <id> = (number of int, weights, method for representation, density, cumulative, left bound, right bound);
{

    long    mark1 = source.FirstSpaceIndex  (0,kStringEnd,kStringDirectionForward),
            mark2 = source.Find             ('=', mark1, -1);

    _String catID (source,mark1+1,mark2-1);

    if (mark1==-1 || mark2==-1 || catID.empty () ) {
        HandleApplicationError("Category variable declaration missing a valid identifier");
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
    HandleApplicationError ("Expected: category <id> = (number of intervals, weights, method for representation, density, cumulative, left bound, right bound,<optional mean cumulative function>,<optional hidden markov matrix>);");
    return false;
}


//____________________________________________________________________________________
bool    _ElementaryCommand::ConstructChoiceList(_String&source, _ExecutionList&target) {
    _List args;


    ExtractConditions (source,blChoiceList.length(),args,',');
    if (args.lLength<5UL) {
        HandleApplicationError  ("ChoiceList needs at least 5 arguments");
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

bool    _ElementaryCommand::ConstructTree (_String&source, _ExecutionList&target)
// Tree   treeid = (...) or Topology = (...);
{
    long    mark1 = source.FirstSpaceIndex(0,-1,kStringDirectionForward);
    if (mark1 > 0) {
      mark1 = source.FirstNonSpaceIndex (mark1 + 1, -1);
    }


    long    mark2 = source.FindTerminator(mark1, "=");
    long    mark3 = mark2;

    if ( mark1 < 0 || mark2 < 0 || mark2 - mark1 < 1) {
        HandleApplicationError ("Tree declaration missing a valid identifier");
        return false;
    }

    _String dsID = source.Cut (mark1,mark2-1);
    // now look for the opening paren

    //(long& from, char open, char close, bool respectQuote, bool respectEscape)
    mark3 = source.ExtractEnclosedExpression (mark1, '(', ')', fExtractRespectQuote | fExtractRespectEscape);


    if (mark1 < 0 || mark3 < 0 || mark3 <= mark1) {
        mark1 = mark2+1;
        mark3 = source.FindTerminator (mark1,";")-1;
    }

    _ElementaryCommand * dsc = new _ElementaryCommand(source.BeginsWith(blTree)?7:54);

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

    long  mark1 = source.FirstNonSpaceFollowingSpace (0,-1,kStringDirectionForward),
          mark2 = source.FindTerminator(mark1+1, "=");

    _String dsID    (source,mark1,mark2-1),
            command;

    if ( mark1==-1 || mark2==-1 || dsID.empty()) {
        HandleApplicationError ("DataSetFilter declaration missing a valid identifier");
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
        HandleApplicationError ("Expected: DataSetFilter	  dataSetFilterid = CreateFilter (datasetid,unit,vertical partition,horizontal partition,alphabet exclusions); or Permute/Bootstrap (dataset/filter,<atom>,<column partition>)");
        return false;
    }


    ExtractConditions (source,mark1+1,pieces,',');
    if (!(pieces.lLength>=2UL || (pieces.lLength == 1UL && dsf->code == 6))) {
        HandleApplicationError ("Parameter(s) missing in DataSetFilter definition.");
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

    long    mark1 = source.FirstSpaceIndex(0,-1,kStringDirectionForward),
            mark2 = source.Find ('=', mark1, -1);

    _String modelID (source,mark1+1,mark2-1);

    if (mark1==-1 || mark2==-1 || !modelID.IsValidIdentifier(fIDAllowCompound)) {
        HandleApplicationError ("Model declaration missing a valid identifier.");
        return false;
    }

    // now look for the opening paren
    mark1 = source.Find ('(',mark2,-1);
    _List pieces;
    ExtractConditions (source,mark1+1,pieces,',');

    if (pieces.lLength<2) {
        HandleApplicationError ("Parameter(s) missing in Model definition. Must have a matrix and a compatible equilibrium frequencies vector.");
        return false;
    } else {
        if (pieces.lLength>3) {
            HandleApplicationError ("Too many parameters (3 max) in Model definition");
            return false;
        }
    }

    _ElementaryCommand * model = new _ElementaryCommand(31);
    model->parameters&&(&modelID);
    model->addAndClean (target,&pieces,0);
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
            HandleApplicationError("An if-then-else scoping error. Check opening and closing brackets and double else's.");
            return false;
        }
        branch1 = simpleParameters[0];
    }

    simpleParameters.Clear();
    simpleParameters<<branch1;
    simpleParameters<<branch2;
    if (source) {
        parameters && source;
    } else if (oldFla) {
        simpleParameters<<oldFla;
    }

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
        HandleApplicationError ("Likelihood function declaration missing a valid identifier");
        return false;
    }

    _String lfID (source,mark1,mark2-1);
    // now look for the opening paren

    _List pieces;
    mark2 ++;
    mark1 = source.ExtractEnclosedExpression(mark2, '(', ')', fExtractRespectQuote | fExtractRespectEscape);

    if ( mark1==-1 || mark2==-1 || mark1<mark2 ) {
        HandleApplicationError ("Expected: Likelihood Function ident = (tree1, datasetfilter1,...)");
        return false;
    }

    ExtractConditions (source,mark2+1,pieces,',');
   _ElementaryCommand*  dsc = new _ElementaryCommand (11);
    dsc->parameters&&(&lfID);

    if (source.BeginsWith(blLF3)) {
        dsc->simpleParameters << 1;
    }

    dsc->addAndClean(target,&pieces,0);
    return true;
}



//____________________________________________________________________________________
bool    _ElementaryCommand::ConstructFunction (_String&source, _ExecutionList& chain)
// syntax: function <ident> (comma separated list of parameters) {body}
{


    bool    isFFunction = source.BeginsWith (blFFunction),
            isLFunction = ! isFFunction && source.BeginsWith (blLFunction),
            isNameSpace = ! isFFunction && ! isLFunction && source.BeginsWith (blNameSpace);

    if (!isNameSpace) {
      if (isInFunction == _HY_FUNCTION) {
        HandleApplicationError ("Nested function declarations are not allowed");
        return false;
      }

    }



    long    mark1 = source.FirstNonSpaceIndex(isNameSpace ? blNameSpace.length(): ((isFFunction||isLFunction)?blFFunction.length():blFunction.length() ),-1,kStringDirectionForward),
            mark2 = source.Find (isNameSpace ? '{' : '(', mark1, -1);


    if ( mark1==-1 || mark2==-1 || mark1>mark2-1) {
        HandleApplicationError      (_String("Function declaration missing a valid function identifier or parameter list.\n-----------\n") & source & "\n-----------\n");
        isInFunction = _HY_NO_FUNCTION;
        return false;
    }

    _String*    funcID  = new _String(source.Cut (mark1,mark2-1));

    if (!funcID->IsValidIdentifier(fIDAllowCompound)) {
      HandleApplicationError      (_String("Not a valid function/namespace identifier '") & *funcID & "'");
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


      if (upto==source.length() || source[upto]!='{' || source (-1)!='}') {
          HandleApplicationError (_String("Function declaration is missing a valid function body."));
          isInFunction= _HY_NO_FUNCTION;
          return false;
      }

      _String extraNamespace;
      if (isLFunction)
          extraNamespace = _HYGenerateANameSpace();

      for (long k = 0UL; k < arguments.lLength; k++) {

          _String*   namespaced = new _String(chain.AddNameSpaceToID (*(_String*)arguments(k), & extraNamespace));
          if ((*namespaced)(-1L) == '&') {
            namespaced->Trim(0,namespaced->length() -2);
            argument_types << kBLFunctionArgumentReference;
          } else {
            argument_types << kBLFunctionArgumentNormal;
          }
          arguments.Replace (k,namespaced,false);
      }


      _String          sfunctionBody (source, upto+1,source.length ()-2);
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
        batchLanguageFunctionClassification.lData[mark1] = isLFunction ? kBLFunctionLocal :( isFFunction? kBLFunctionSkipUpdate :  kBLFunctionAlwaysUpdate);
      } else {
          batchLanguageFunctions.AppendNewInstance(functionBody);
          batchLanguageFunctionNames.AppendNewInstance(funcID);
          batchLanguageFunctionParameterLists &&(&arguments);
          batchLanguageFunctionParameterTypes &&(&argument_types);
          batchLanguageFunctionClassification <<(isLFunction ? kBLFunctionLocal :( isFFunction? kBLFunctionSkipUpdate :  kBLFunctionAlwaysUpdate));
      }
    } else {
      if (mark2 == source.length () || source[mark2]!='{' || source (-1L) !='}') {
        HandleApplicationError (_String("Namespace declaration is missing a body."));
        isInFunction= _HY_NO_FUNCTION;
        return false;
      }
      _String          namespace_text (source, mark2+1,source.length()-2);
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
    long    mark1 = source.FirstNonSpaceIndex(blReturn.length(),-1,kStringDirectionForward);

    _ElementaryCommand ret;

    ret.code = 14;

    if (mark1!=-1) {
        _String cut_s;
        if (source (-1) ==';') {
            cut_s = source.Cut (mark1,source.length () - 2L);
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

void    ReadBatchFile (_String& fName, _ExecutionList& target) {
// read/parse a file into an execution list
// THE function!!!

    ProcessFileName(fName, target.nameSpacePrefix);

    if (terminate_execution) {
        return;
    }
    /*#else
        _Variable optprec (optimizationPrecision);
        _Constant precvalue (0.01);
        FetchVar(LocateVarByName (optimizationPrecision))->SetValue(&precvalue);
    #endif*/

    FILE            *f = doFileOpen (fName.get_str (), "rb");
    SetStatusLine   ("Parsing File");
    if (!f) {
        HandleApplicationError (_String("Could not read batch file '") & fName & "'.\nPath stack:\n\t" & GetPathStack("\n\t"));
    } else {
        _String source_file (f);

        if (source_file.BeginsWith ("#NEXUS",false)) {
            ReadDataSetFile (f,1,nil,&fName, nil, &defaultTranslationTable, &target);
        } else {
            target.BuildList (source_file);
            target.sourceFile = fName;
        }
        fclose (f);
    }
}



//____________________________________________________________________________________
void    SerializeModel  (_StringBuffer & rec, long theModel, _AVLList* alreadyDone, bool completeExport)
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

    if (!alreadyDone || alreadyDone->Find ((BaseRef)tV2->GetIndex()) < 0) {
        if (alreadyDone) {
            alreadyDone->Insert ((BaseRef)tV2->GetIndex());
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

        _StringBuffer glVars (128L),
                locVars(128L);

        ExportIndVariables (glVars,locVars, &ind);
        ExportDepVariables (glVars,locVars, &dep);
        rec << glVars <<locVars;
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

    rec << "\nModel "
     << *((_String*)modelNames (theModel))
     << "=(";
    if (theExp) {
        rec << _String((_String*)(theExp->toStr(kFormulaStringConversionNormal))).Enquote();
     } else {
        rec << *tV->GetName();
    }
    rec << ',' << *tV2->GetName();
    if (theExp) {
        rec << ',' << explicitFormMExp;
    } else if (!mByF) {
        rec << ",0";
    }
    rec << ");\n";
}
