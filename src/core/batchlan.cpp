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

                                 dialogPrompt,
                                defFileNameValue;


//____________________________________________________________________________________


_String
blWhile                    ("while("),         // moved
blFunction                 ("function "),      // moved
blFFunction                ("ffunction "),     // moved
blLFunction                ("lfunction "),     // moved
blCFunction                ("cfunction "),     // moved
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

void    ReportMPIError      (int code, bool send) {
    if (code != MPI_SUCCESS) {
        _String errMsg = "MPI Error while ";
        if (send) {
            errMsg = errMsg & "sending";
        } else {
            errMsg = errMsg & "receiving";
        }

        errMsg = errMsg & _String(" with code ") & (long)code;
        HandleApplicationError (errMsg);
    }
}

#define MPI_SEND_CHUNK 0xFFFFFFL

//____________________________________________________________________________________

void    MPISendString       (_String const& theMessage, long destID, bool isError)
{

    long    messageLength = theMessage.length(),
            transferCount = 0L;

    if (isError) {
        messageLength = -messageLength;
    }

    ReportMPIError(MPI_Send(&messageLength, 1, MPI_LONG, destID, HYPHY_MPI_SIZE_TAG, MPI_COMM_WORLD),true);

    if (messageLength == 0L) {
        return;
    }

    if (isError) {
        messageLength = -messageLength;
    }

    while (messageLength-transferCount>MPI_SEND_CHUNK) {
        printf("%s",theMessage.get_str());
        ReportMPIError(MPI_Send((void*)(theMessage.get_str()+transferCount), MPI_SEND_CHUNK, MPI_CHAR, destID, HYPHY_MPI_STRING_TAG, MPI_COMM_WORLD),true);
        transferCount += MPI_SEND_CHUNK;
    }

    if (messageLength-transferCount) {
        ReportMPIError(MPI_Send((void*)(theMessage.get_str()+transferCount), messageLength-transferCount, MPI_CHAR, destID, HYPHY_MPI_STRING_TAG, MPI_COMM_WORLD),true);
    }

    //ReportMPIError(MPI_Send(&messageLength, 1, MPI_LONG, destID, HYPHY_MPI_DONE_TAG, MPI_COMM_WORLD),true);

    _FString*    sentVal = new _FString ((_String*)theMessage.makeDynamic());
    _Variable *   mpiMsgVar = CheckReceptacle (&hy_env::mpi_last_sent_message, kEmptyString, false);
    mpiMsgVar->SetValue (sentVal, false, true, NULL);
    //setParameter (mpiLastSentMsg, &sentVal);

}

//____________________________________________________________________________________
_String*    MPIRecvString       (long senderT, long& senderID) {
    _String*    theMessage = nil;
    long        messageLength = 0L,
                transferCount = 0L;

    int         actualReceived = 0;
    bool        isError       = false;

    if  (senderT<0) {
        senderT = MPI_ANY_SOURCE;
    }

    MPI_Status  status;
  
    // nonagressive polling mode
    
    //ReportWarning ("Step 1");
  
    int message_received = 0;
    while (! message_received) {
      MPI_Iprobe (senderT, HYPHY_MPI_SIZE_TAG, MPI_COMM_WORLD, &message_received, MPI_STATUS_IGNORE);
      usleep (100);
    }

    //ReportWarning ("Step 2");
    // nonagressive polling mode
  
  
    ReportMPIError(MPI_Recv(&messageLength, 1, MPI_LONG, senderT, HYPHY_MPI_SIZE_TAG, MPI_COMM_WORLD,&status),false);

    if (messageLength < 0) {
        isError = true;
        messageLength = -messageLength;
    }
    
    //ReportWarning ("Step 3");
    //printf ("MPIRecvString size tag %ld (size chunk %ld) \n",messageLength, MPI_SEND_CHUNK);

    if (!isError) {
        //MPI_Get_count (&status,MPI_CHAR,&actualReceived);

        if (messageLength==0) {
            return new _String;
        }

        theMessage = new _String ((unsigned long)messageLength);

        senderT = senderID = status.MPI_SOURCE;

        while (messageLength-transferCount>MPI_SEND_CHUNK) {
            ReportMPIError(MPI_Recv((void*)(theMessage->get_str()+transferCount), MPI_SEND_CHUNK, MPI_CHAR, senderT, HYPHY_MPI_STRING_TAG, MPI_COMM_WORLD,&status),false);
            MPI_Get_count (&status,MPI_CHAR,&actualReceived);
            if (actualReceived!=MPI_SEND_CHUNK) {
                HandleApplicationError ("Failed in MPIRecvString - some data was not properly received\n");
            }
            //return    nil;
            transferCount += MPI_SEND_CHUNK;
        }

        if (messageLength-transferCount) {
            //printf ("Clause 2 %d %d\n", messageLength-transferCount, theMessage->length());
            ReportMPIError(MPI_Recv((void*)(theMessage->get_str()+transferCount), messageLength-transferCount, MPI_CHAR, senderT, HYPHY_MPI_STRING_TAG, MPI_COMM_WORLD,&status),false);
            MPI_Get_count (&status,MPI_CHAR,&actualReceived);
            if (actualReceived!=messageLength-transferCount) {
                HandleApplicationError ("Failed in MPIRecvString - some data was not properly received\n");
            }
            //return    nil;
        }
        
        //ReportWarning ("Step 4");

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
    HBLObjectRef formRes = nameForm.Compute();

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
        HBLObjectRef formRes = nameForm.Compute();
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

/** This function returns expects that the caller will handle reference counting on the returned object;
    it will be returned with +1 reference, i.e. it needs to be deleted / managed by the caller */

HBLObjectRef   ProcessAnArgumentByType (_String const* expression, _VariableContainer const* theP, long objectType, _ExecutionList* currentProgram) {
    _String   errMsg;
    _Formula  expressionProcessor (*expression, theP, currentProgram?&errMsg:nil);

    if (errMsg.nonempty() && currentProgram) {
        currentProgram->ReportAnExecutionError (errMsg);
    }
    else {
        HBLObjectRef expressionResult = expressionProcessor.Compute(0,theP);
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
   HBLObjectRef getString = ProcessAnArgumentByType (data, theP, STRING, currentProgram);

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
  return *(_ExecutionList*)batchLanguageFunctions.GetItem (idx);
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
        if (! a_name -> Equal( hbf_name)) {
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
  if (start_here >= 0L && start_here < batchLanguageFunctionNames.countitems()) {

    _SimpleList delete_me (batchLanguageFunctionNames.countitems()-start_here, start_here, 1L);
    
    for (unsigned long k = 0UL; k < delete_me.countitems(); k++) {
      batchLanguageFunctionNamesIndexed.Delete (batchLanguageFunctionNames.GetItem (delete_me.get (k)), true);
    }

    batchLanguageFunctionNames.DeleteList           (delete_me);
    batchLanguageFunctions.DeleteList               (delete_me);
    batchLanguageFunctionClassification.DeleteList  (delete_me);
    batchLanguageFunctionParameterLists.DeleteList  (delete_me);
    batchLanguageFunctionParameterTypes.DeleteList  (delete_me);
  } else {
    if (start_here < 0) {
        batchLanguageFunctionNames.Clear           ();
        batchLanguageFunctions.Clear               ();
        batchLanguageFunctionClassification.Clear  ();
        batchLanguageFunctionParameterLists.Clear  ();
        batchLanguageFunctionParameterTypes.Clear  ();
        batchLanguageFunctionNamesIndexed.Clear    (true);
    }
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
            long idx = batchLanguageFunctionNamesIndexed.GetDataByKey(&test_id);
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
   return batchLanguageFunctionNamesIndexed.GetDataByKey (&s);
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
        dataSetList.list_data[k]     = (long)theDS;
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
                if (otherVars.BinaryFind(myVars.list_data[k])<0) {
                    wastedVars << myVars.list_data[k];
                }

            myVars.Clear();

            unsigned long component_count = me->CountObjects(kLFCountPartitions);

            for (long tree_index = 0UL; tree_index < component_count; tree_index++) {
                _TheTree* thisTree = me->GetIthTree(tree_index);
                thisTree->CompileListOfModels (myVars);
                _TreeIterator ti (thisTree, _HY_TREE_TRAVERSAL_POSTORDER);
                while (_CalcNode* tNode = ti.Next()) {
                    tNode->SetValue (new _Constant (tNode->ComputeBranchLength()),false,true,NULL);
                }
                thisTree->RemoveModel();
            }

            for (unsigned long k=0UL; k<myVars.lLength; k++)
                if (otherModels.BinaryFind (myVars.list_data[k])<0) {
                    KillModelRecord (myVars.list_data[k]);
                }

            for (unsigned long k=0UL; k<wastedVars.lLength; k++) {
                //printf ("Deleting %ld (%s)\n", wastedVars.list_data[k],  ->GetName()->getStr());
                _Variable * check_me = LocateVar(wastedVars.list_data[k]);
                if (check_me) {
                  DeleteVariable (*check_me->GetName());
                }
            }

        }
        
        me->UnregisterListeners();

        if (lfID<likeFuncList.lLength-1) {
            DeleteObject(likeFuncList(lfID));
            likeFuncList.list_data[lfID] = 0L;
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
        DeleteVariable (*LocateVar(l.list_data[k])->GetName());
    }

    l.Clear ();

    unsigned long partition_count = lf->CountObjects(kLFCountPartitions);

    for (unsigned long k=0UL; k<partition_count; k++) {
        _TheTree * thisTree = lf->GetIthTree(k);
        thisTree->CompileListOfModels (l);
        DeleteVariable (*thisTree->GetName());
    }

    for (unsigned long k=0UL; k<l.lLength; k++) {
        KillModelRecord (l.list_data[k]);
    }

    KillLFRecord (lfID);
}

//__________________________________________________________

void KillDataSetRecord (long dsID)
{
    if (dsID<dataSetList.lLength-1) {
        DeleteObject(dataSetList(dsID));
        dataSetList.list_data[dsID] = 0;
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
        if (modelTypeList.list_data[i]) {
            delete (_Formula*)(modelMatrixIndices.list_data[i]);
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


    if (modelTypeList.list_data[mdID]) {
        delete (_Formula*)(modelMatrixIndices.list_data[mdID]);
    } else {
        _Variable * modelMatrix = nil,
                    * freqMatrix  = nil;

        bool      multByFreqs   = false;



        _SimpleList  saveTheseVariablesAux;
        _AVLList     saveTheseVariables (&saveTheseVariablesAux);

        for (long k = 0; k < modelNames.lLength; k++) {
            if (k != mdID && ((_String*)modelNames(k))->nonempty()) {
                if (modelTypeList.list_data[k]) {
                    _SimpleList dependantMatrices;
                    ((_Formula*)(modelMatrixIndices.list_data[k]))->ScanFForType (dependantMatrices, MATRIX);
                    for (long k2 = 0; k2 < dependantMatrices.lLength; k2++) {
                        saveTheseVariables.Insert((BaseRef)dependantMatrices.list_data[k2]);
                    }
                } else {
                    RetrieveModelComponents(k, modelMatrix, freqMatrix, multByFreqs);

                    if (modelMatrix) {
                        saveTheseVariables.Insert((BaseRef)modelMatrix->get_index());
                    }
                    if (freqMatrix) {
                        saveTheseVariables.Insert((BaseRef)freqMatrix->get_index());
                    }
                }
            }
        }

        RetrieveModelComponents(mdID, modelMatrix, freqMatrix, multByFreqs);
        if (modelMatrix && saveTheseVariables.Find ((BaseRef)modelMatrix->get_index()) < 0) {
            DeleteVariable (*modelMatrix->GetName());
        }
        if (freqMatrix && saveTheseVariables.Find ((BaseRef)freqMatrix->get_index()) < 0) {
            DeleteVariable (*freqMatrix->GetName());
        }
    }

    if (mdID<modelNames.lLength-1) {
        modelMatrixIndices.list_data[mdID] = -1;
        modelTypeList.list_data[mdID] = 0;
        modelFrequenciesIndices.list_data[mdID] = -1;
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
    kwargs              = nil;
    kwarg_tags          = nil;
    currentKwarg        = 0;

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
void _ExecutionList::SetKWArgs(_AssociativeList* kwarg_list) {
    DeleteAndZeroObject(kwargs);
    if (kwarg_list && kwarg_list->countitems()) {
        kwargs = (_AssociativeList*)kwarg_list->makeDynamic();
    }
}

//____________________________________________________________________________________
void _ExecutionList::ClearExecutionList (void) {
    if (cli) {
        delete [] cli->values;
        delete [] cli->stack;
        delete [] cli->is_compiled;
        delete cli;
        cli = nil;
    }
    
    DeleteAndZeroObject (profileCounter);
    DeleteAndZeroObject (stdinRedirect);
    DeleteAndZeroObject (stdinRedirectAux);
    DeleteAndZeroObject (nameSpacePrefix);
    DeleteAndZeroObject (kwargs);
    DeleteAndZeroObject (kwarg_tags);
    
    ResetFormulae();
    DeleteAndZeroObject (result);
}

//____________________________________________________________________________________

_ExecutionList::~_ExecutionList (void) {
    ClearExecutionList();
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
        Res->result = (HBLObjectRef)result->makeDynamic();
    }

    return Res;
}

//____________________________________________________________________________________

void        _ExecutionList::Duplicate   (BaseRefConst source) {
    _List::Duplicate    (source);

    _ExecutionList const* s = (_ExecutionList const*)source;

    if (s->result) {
        result=(HBLObjectRef)s->result->makeDynamic();
    }

    errorHandlingMode  = s->errorHandlingMode;
    errorState         = s->errorState;
}


//____________________________________________________________________________________
void    _ExecutionList::ReportAnExecutionError (_String errMsg, bool doCurrentCommand, bool appendToExisting) {
    if (doCurrentCommand) {
        _ElementaryCommand *theCommand = FetchLastCommand();
        if (theCommand) {
            errMsg = errMsg & " in call to " & _String ((_String*)theCommand->toStr());
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
_String*    _ExecutionList::FetchFromStdinRedirect (_String const * dialog_tag, bool handle_multi_choice, bool do_echo) {
// grab a string from the front of the input queue
// complain if nothing is left
    if (! (has_stdin_redirect() || has_keyword_arguments())) {
        throw _String ("No input buffer / keyword arguments was given for a redirected standard input read.");
        return new _String;
    }
    
    /**
        if keyword arguments are present, try them first
    */
    
    HBLObjectRef user_argument = nil;
    _List        ref_manager;
    _String*     kwarg_used = nil;
    
    auto echo_argument = [&do_echo] (_String const * kwarg, _String const& value) -> void {
        if (do_echo && kwarg) {
            bool do_markdown = hy_env :: EnvVariableTrue(hy_env :: produce_markdown_output);
            if (do_markdown) {
                NLToConsole(); BufferToConsole(">"); StringToConsole(*kwarg); BufferToConsole( " –> "); BufferToConsole(value); NLToConsole();
            } else {
                StringToConsole(*kwarg); BufferToConsole( ": "); BufferToConsole(value); NLToConsole();
            }
        }
    };
    
    try {
        if (has_keyword_arguments()) {
            if (kwarg_tags && kwarg_tags->countitems() > currentKwarg) {
                _List* current_tag = (_List*)kwarg_tags->GetItem(currentKwarg++);
                
                // check to see if we need to match with the current dialog prompt
                if (current_tag -> countitems() > 3UL) {
                    _String check_against = dialog_tag ? *dialog_tag : dialogPrompt;
                    if (check_against != *(_String*)current_tag->GetItem(3)) {
                        throw (1L);
                    }
                }
                     
                if (kwargs) {
                    user_argument = kwargs->GetByKey(*(_String*)current_tag->GetItem(0));
                }

                if (user_argument) { // user argument provided
                    user_argument -> AddAReference();
                    kwarg_used = (_String*)current_tag->GetItem(0);
                    kwargs->DeleteByKey(*kwarg_used);
                    ref_manager < user_argument;
                } else { // see if there are defaults
                    if (current_tag->countitems() > 2 && ! ignore_kw_defaults) {
                        _String * default_value = (_String*)current_tag->GetItem(2);
                        kwarg_used = (_String*)current_tag->GetItem(0);
                        if (default_value) {
                            user_argument = new _FString (*(_String*)current_tag->GetItem(2));
                            ref_manager < user_argument;
                        }
                    }
                }
                
                //printf ("%ld => %s\n", currentKwarg - 1, user_argument ? _String ((_String*)user_argument->toStr()).get_str() : "''");
            }
        }
    } catch (long) {
        return FetchFromStdinRedirect (dialog_tag, handle_multi_choice, do_echo);
    }
    
    if (user_argument) {
        if (user_argument->ObjectClass() == STRING) {
            echo_argument (kwarg_used, ((_FString*)user_argument)->get_str());
            return new _String (((_FString*)user_argument)->get_str());
        } else {
            if (handle_multi_choice) {
                echo_argument (kwarg_used, _String ((_String*)user_argument->toStr()));
                user_argument->AddAReference();
                throw (user_argument);
            } else {
                throw _String ("Multi-choice keyword arguement not supported in this context");
            }
        }
    }
    
    if (has_stdin_redirect()) {
        long d = stdinRedirect->First();
        if (d<0) {
            throw _String ("Ran out of input in buffer during a redirected standard input read.");
        }

        _String *sendBack = (_String*)stdinRedirect->GetXtra (d);
        //printf ("Consumed stdin redrect %ld => %s\n", d, sendBack->get_str());
        sendBack->AddAReference();
        if (do_echo) {
            StringToConsole(*sendBack);
            NLToConsole();
        }
        stdinRedirect->Delete ((*(_List*)stdinRedirect->dataList)(d),true);
        return sendBack;
    }
    
    throw kNoKWMatch;
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

_StringBuffer const       _ExecutionList::GenerateHelpMessage(_AVLList * scanned_functions)  const {
    _StringBuffer help_message;
    
    _List ref_manager;
    bool  nested = true;
    if (!scanned_functions) {
        _List * _aux_list = new _List;
        scanned_functions = new _AVLList (_aux_list);
        ref_manager < _aux_list < scanned_functions;
        nested = false;
    }
    
    auto simplify_string = [] (_String const * s) -> const _String {
        _String sc (*s);
        if (sc.IsALiteralArgument(true)) {
            return sc;
        }
        return sc & " [computed at run time]";
    };

    ForEach ([&help_message, simplify_string, this, scanned_functions] (BaseRef command, unsigned long index) -> void {
        _ElementaryCommand * this_command = (_ElementaryCommand * )command;
        if (this_command->code == HY_HBL_COMMAND_KEYWORD_ARGUMENT) {
            _String * def_value = this_command->GetIthParameter(2L, false),
                    * applies_to = this_command->GetIthParameter(3L, false);

            if (def_value && (*def_value == kNoneToken || *def_value == kNullToken)) {
                def_value = nil;
            }
            
            help_message << simplify_string(this_command->GetIthParameter(0L)) << (def_value == nil ? (applies_to ? " [conditionally required]" : " [required]"): "") << '\n'
                         << '\t' << simplify_string(this_command->GetIthParameter(1L)) << '\n';
            
            
            if (def_value) {
                help_message << "\tdefault value: " << simplify_string(def_value) << '\n';
            }
             if (applies_to) {
                help_message << "\tapplies to: " << simplify_string(applies_to) << '\n';
            }
            help_message << '\n';
        } else {
            if (this_command->code == HY_HBL_COMMAND_FORMULA) {
                _List      hbl_functions;
                _AVLListX other_functions (&hbl_functions);
                this_command->BuildListOfDependancies(other_functions, true, *this);
                
                for (AVLListXIteratorKeyValue function_iterator : AVLListXIterator (&other_functions)) {
                    _String * function_name = (_String *)other_functions.Retrieve (function_iterator.get_index());
                    if (scanned_functions->Insert (new _String (*function_name),0,true, true) >= 0) {
                        long idx = FindBFFunctionName(*function_name);
                        if (idx >= 0) {
                            help_message << GetBFFunctionBody(idx).GenerateHelpMessage(scanned_functions);
                        }
                    }
                }
            }
        }
        
    });
    
    if (help_message.empty() && !nested) {
        help_message << "No annotated keyword arguments are available for this analysis\n";
    }
    
    return help_message;
}

//____________________________________________________________________________________

HBLObjectRef       _ExecutionList::Execute     (_ExecutionList* parent, bool ignore_CEL_kwargs) {

  //setParameter(_hyLastExecutionError, new _MathObject, nil, false);
  try{

    _ExecutionList*      stashCEL = currentExecutionList;
    callPoints << currentCommand;
    executionStack       << this;
    
    _AVLListXL * stash1 = nil;
    _List* stash2 = nil,
         * stash_kw_tags = nil;// recursion
    _AssociativeList * stash_kw = nil;
      
    currentKwarg         = (stashCEL && !ignore_CEL_kwargs) ? stashCEL->currentKwarg : 0;

    if (parent && (stdinRedirect == nil || kwargs == nil)) {
        if (!stdinRedirect) {
            stash1 = stdinRedirect = parent->stdinRedirect;
            stash2 = stdinRedirectAux = parent->stdinRedirectAux;
            if (stash1) {
                stash1->AddAReference();
                stash2->AddAReference();
            }
        }
        if (!kwargs) {
            stash_kw_tags = kwarg_tags = parent->kwarg_tags;
            stash_kw = kwargs = parent->kwargs;
            if (stash_kw_tags) stash_kw_tags->AddAReference();
            if (stash_kw) stash_kw->AddAReference();
            currentKwarg = parent->currentKwarg;
        }
    } else {
      parent = nil;
    }

    _FString            cfp (PeekFilePath() ? *PeekFilePath () :kEmptyString),
                        *stashed = (_FString*)hy_env::EnvVariableGet(hy_env::path_to_current_bf, STRING);

    if (stashed) {
        stashed = (_FString*)stashed->makeDynamic();
    }

    hy_env::EnvVariableSet(hy_env::path_to_current_bf, &cfp, true);

    DeleteAndZeroObject        (result);
    currentExecutionList = this;
    currentCommand       = 0;

    terminate_execution  = false;
      
      
    bool is_c = is_compiled();
    if (is_c) {
      //PopulateArraysForASimpleFormula (cli->varList, cli->values);
      cli->is_compiled[0] = false;
    }

    while (currentCommand<lLength) {
        
        if (is_c) {
            if ( cli->is_compiled [currentCommand+1] == false) {
               if (cli->is_compiled[0]) {
                  CopyCLIToVariables();
                  cli->is_compiled [0] = false;
               }
                
            } else {
                if (cli->is_compiled[0] == false) {
                    PopulateArraysForASimpleFormula (cli->varList, cli->values);
                    cli->is_compiled[0] = true;
                }
            }
        }
        
        if (doProfile == 1 && profileCounter) {
            long        instCounter = currentCommand;
            hyFloat  timeDiff    = 0.0;

            TimeDifference timer;
            (((_ElementaryCommand**)list_data)[currentCommand])->Execute(*this);
            timeDiff   = timer.TimeSinceStart();


          if (profileCounter) {
            // a call to _hyphy_profile_dump can set this to NULL
            profileCounter->theData[instCounter*2]   += timeDiff;
            profileCounter->theData[instCounter*2+1] += 1.0;
          }
        } else {
            (((_ElementaryCommand**)list_data)[currentCommand])->Execute(*this);
        }

        if (terminate_execution) {
            break;
        }
    }
    currentCommand = callPoints.list_data[callPoints.lLength-1];
    callPoints.Delete (callPoints.lLength-1);
    currentExecutionList = stashCEL;
      
    if (currentExecutionList && !ignore_CEL_kwargs) {
      currentExecutionList->currentKwarg = currentKwarg;
    }

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
      kwargs = nil;
      kwarg_tags = nil;
      parent->currentKwarg = currentKwarg;
      BatchDeleteObject (stash1, stash2, stash_kw, stash_kw_tags);
    }
      
    if (is_compiled(0)) {
      CopyCLIToVariables();
    }

  } catch (const _String& err) {
    HandleApplicationError(err);
  }

    return result;
}

//____________________________________________________________________________________

void        _ExecutionList::ExecuteAndClean     (long g) {
    Execute ();
    ClearBFFunctionLists    (g);
}

//____________________________________________________________________________________

bool        _ExecutionList::TryToMakeSimple     (bool partial_ok) {
    _SimpleList     varListAux,
                    formulaeToConvert,
                    parseCodes;
  
    _AVLList        varList (&varListAux);

    long            stackDepth  = 0L;
    bool            status      = true;
    unsigned long   k = 0UL;
    
    bool            *is_compiled = new bool [countitems()+1];
    InitializeArray (is_compiled, countitems()+1, false);
    
    for (; k<lLength && status; k++) {
        _ElementaryCommand * aStatement = (_ElementaryCommand*)(*this)(k);
        switch (aStatement->code) {
        case 0: { // expression
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
                            f->GetIthTerm (0)->SetAVariable(mx->get_index());
                            _Operation * last = f->GetIthTerm(assignment_length-1);
                            if (! (last->TheCode() == HY_OP_CODE_MCOORD && last->GetNoTerms() == 2)) throw 0;

                            f2->GetList() << f->GetList();
                            f->Clear();

                            _Formula *t = f2;
                            f2 = f;
                            f  = t;

                          }

                        } catch (int e) {
                          if (partial_ok) {
                              parseCodes << -1;
                              continue;
                          }
                          status = false;
                          break;
                        }
                      
                        aStatement->simpleParameters<<parseCode;
                        aStatement->simpleParameters<<(long)f;
                        aStatement->simpleParameters<<(long)f2;
                        aStatement->simpleParameters<<fpc.assignmentRefID();


                        formulaeToConvert << (long)f;
                        is_compiled [k+1] = true;


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
            
            parseCodes << -1;
            if (!partial_ok) {
                status = false;
            }
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
                    is_compiled [k+1] = true;
                } else {
                    if (partial_ok) {
                        parseCodes << -1;
                        continue;
                    }
                    status = false;
                }
            } else {
                is_compiled[k+1] = true;
            }
            break;
                
        case 14: // return statements are OK
            parseCodes << -1;
            break;
                
        default:
            status = false;
        }
        if (status == false) {
            ReportWarning (_String ("Failed to compile an execution list: offending command was\n") & _String (((_String*)aStatement->toStr())));
        }
    }

    if (status) {
        if (formulaeToConvert.nonempty()) {
            cli = new _CELInternals;
            //varList.ReorderList();
            cli->values = new _SimpleFormulaDatum[varList.countitems()+1];
            cli->stack  = new _SimpleFormulaDatum[stackDepth+1];
            cli->is_compiled = is_compiled;

            _SimpleList  avlData;
            _AVLListX    avlList (&avlData);

            for (unsigned long fi = 0; fi < formulaeToConvert.lLength; fi++) {
                ((_Formula*)formulaeToConvert(fi))->ConvertToSimple (varList);
            }

            for (unsigned long vi = 0; vi < varListAux.countitems(); vi++) {
                avlList.Insert ((BaseRef)varListAux.list_data[vi], vi);
            }

            for (unsigned long ri = 0; ri<parseCodes.lLength; ri++) {
                if (parseCodes.list_data[ri] < 0) {
                    cli->storeResults << -1;
                } else {
                    cli->storeResults << avlList.GetXtra (avlList.Find ((BaseRef) parseCodes.list_data[ri]));
                }
                //printf ("\n%ld\n",  cli->storeResults.list_data[ri]);
            }
            cli->varList.Duplicate(&varListAux);
        }
    } else {
        // clean up partially converted statements
      delete [] is_compiled;
      for (unsigned long k2 = 0UL; k2 < k; k2++) {
        GetIthCommand(k2)->simpleParameters.Pop(3UL);
      
      }
    }

    return status;
}

//____________________________________________________________________________________

void        _ExecutionList::CopyCLIToVariables(void) {
    
    cli->varList.Each ([this] (long index, unsigned long idx) -> void {
        _Variable * mv = LocateVar(index);
        if (mv->ObjectClass() == NUMBER) {
            mv->SetValue (new _Constant (this->cli->values[idx].value),false,true,NULL);
        }
    });
}

//____________________________________________________________________________________

void        _ExecutionList::ExecuteSimple       (_ExecutionList * parent) {
    //PopulateArraysForASimpleFormula (cli->varList, cli->values);
    Execute(parent);
    //CopyCLIToVariables();
}

//____________________________________________________________________________________

void        _ExecutionList::ResetFormulae       (void) {
    currentCommand = 0L;
    _SimpleList to_delete_aux;
    _AVLList to_delete (&to_delete_aux);
    while (currentCommand<lLength) {
        _ElementaryCommand* thisCommand = ((_ElementaryCommand**)list_data)[currentCommand];
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
        if (is_compiled(i+1)) {
            (*result) << "[compiled]";
        }
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
  
  if (pieces->countitems() > 5UL) { // expliit list of choices
    
    if (pieces->countitems() % 2 > 0) {
      DeleteObject(choice_list);
     throw (_String ("Must have an even number of arguments for explicitly enumerated choice - description pairs in ") & __PRETTY_FUNCTION__ & ". Had " & _String ((_String*)pieces->toStr()).Enquote());
    }
    
    _List * choices = new _List;
    
    for (unsigned long k = 4UL; k < pieces->countitems(); k+=2) {
        _String * selector    = new _String (*(_String*)pieces->GetItem(k)),
                * desription  = new _String (*(_String*)pieces->GetItem(k+1));
        selector->StripQuotes(); desription->StripQuotes();
        *choices < new _List (selector,desription);
    }
    choice_list->parameters < choices;
    choice_list->simpleParameters << 0L;
  } else {
    choice_list->parameters << pieces->GetItem(4L);
    choice_list->simpleParameters << 1L;
  }
  
  choice_list->addAndClean (*this);

}

  //____________________________________________________________________________________

void  _ExecutionList::BuildExecuteCommandInstruction (_List * pieces, long code) {

  const _String kExecuteCompiled ("compiled"),
                kExecuteEncloseingNamespace ("enclosing_namespace");
  
  _ElementaryCommand * run_source = new _ElementaryCommand (code);
  run_source->parameters<<pieces->GetItem(0);
  
  if (PeekFilePath()) {
    run_source->parameters && *PeekFilePath();
  } else {
    run_source->parameters.AppendNewInstance(new _String);
  }
  
  if (pieces->countitems() > 1UL) {
    if (*(_String*)pieces->GetItem (1UL) == kExecuteCompiled) {
      run_source->simpleParameters << 1;
    } else {
      if (*(_String*)pieces->GetItem (1UL) == kExecuteEncloseingNamespace) {
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

    static const _String kFscanfRewind ("REWIND");

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

    argument_types.ForEach ([&] (BaseRefConst t, unsigned long) -> void {
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

    if (scanf->parameters.countitems()  != scanf->simpleParameters.countitems() + 1UL) {
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

    _List                local_object_manager;
    _StringBuffer        currentLine (128UL);
  
    try {

      while (s.nonempty ()) { // repeat while there is stuff left in the buffer
          _ElementaryCommand::FindNextCommand (s,currentLine);

          if (currentLine.get_char(0)=='}') {
              currentLine.Trim(1,kStringEnd);
          }

          if (currentLine.empty()) {
              continue;
          }

          long prefixTreeCode = _HY_ValidHBLExpressions.FindKey (currentLine, nil, true);

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
                         if (commandExtraInfo->extract_conditions.countitems() == 1 && commandExtraInfo->extract_conditions.get(0) < 0) {
                              if (pieces->countitems() < -commandExtraInfo->extract_conditions.get(0)) {
                                   throw (_String("Incorrect number of arguments (") & (long) pieces->lLength & ") supplied: expected at least " & _String (-commandExtraInfo->extract_conditions.get(0)) & ", while processing '"& currentLine.Cut (0, upto) & "'. ");
                              }
                         } else {
                           throw (_String("Incorrect number of arguments (") & (long) pieces->lLength & ") supplied: expected one of " & _String ((_String*)commandExtraInfo->extract_conditions.toStr()) & ", while processing '"& currentLine.Cut (0, upto) & "'. ");
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
              case HY_HBL_COMMAND_INTEGRATE:
              case HY_HBL_COMMAND_ALIGN_SEQUENCES:
              case HY_HBL_COMMAND_CONSTRUCT_CATEGORY_MATRIX:
              case HY_HBL_COMMAND_KEYWORD_ARGUMENT:
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
              if (currentLine.BeginsWith (blFunction)||currentLine.BeginsWith (blFFunction)||currentLine.BeginsWith (blLFunction) || currentLine.BeginsWith (blNameSpace) || currentLine.BeginsWith (blCFunction)) { // function declaration
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
                              lif  = lastif.list_data[lc-1];

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
                          throw _String("'else' w/o an if to latch on to...");
                      }

                    
                      ((_ElementaryCommand*)((*this)(lif)))->MakeJumpCommand(nil,-1,temp+1,*this);
                      ((_ElementaryCommand*)(*this)(temp))->simpleParameters[0]=countitems();

                      while (lastif.countitems()>=lc) {
                          lastif.Delete(lastif.countitems()-1);
                      }
                  } else {
                      throw (_String ("'else' w/o an 'if' to latch on to..."));
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
                  _StringBuffer next_command;
                  _ElementaryCommand::FindNextCommand (checker,next_command);
                  if (next_command.length ()==currentLine.length()) {
                      if (currentLine.length()>1)
                          while (currentLine (-1L) ==';') {
                              currentLine.Trim (0,currentLine.length()-2);
                          }
                      else {
                          continue;
                      }
                      _ElementaryCommand* oddCommand = new _ElementaryCommand(currentLine);
                      oddCommand->code = 0;
                      oddCommand->parameters.AppendNewInstance (new _String (currentLine));
                      AppendNewInstance (oddCommand);
                  } else {
                      while (currentLine.nonempty()) {
                          _ElementaryCommand::FindNextCommand (currentLine,next_command);
                          BuildList (next_command,bc,processed);
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
        } else if (code == HY_HBL_COMMAND_INIT_ITERATOR) {
            if (simpleParameters.get (1) == ASSOCIATIVE_LIST && simpleParameters.countitems() > 2) {
                delete (AVLListXLIterator*)simpleParameters.get(2);
            } else {
                if ((simpleParameters.get (1) == TREE ||  simpleParameters.get (1) == TOPOLOGY) && simpleParameters.countitems() > 2) {
                    delete (node_iterator<long>*)simpleParameters.get(2);
                }
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
        
        _String command (_HY_ValidHBLExpressions.RetrieveKeyByPayload(i));
        
        if (command.EndsWith('(')) {
            return _StringBuffer (command)
            << _String ((_String*)parameters.Join (", ")) << ");";

        } else {
            return _StringBuffer (command)
                        << '(' << _String ((_String*)parameters.Join (", ")) << ");";
        }
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

        case 5: { // data set contruction
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
        case HY_HBL_COMMAND_SELECT_TEMPLATE_MODEL:
        case HY_HBL_COMMAND_KEYWORD_ARGUMENT:
        case HY_HBL_COMMAND_GET_INFORMATION:
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
            (*string_form) << (code == HY_HBL_COMMAND_FSCANF ? "fscanf(" : "sscanf(")
                        << parameter_to_string (0)
                        << ",\"";

            long shift = 1L;
          
            for (long p = 0; p<simpleParameters.lLength; p++) {
                long theFormat = simpleParameters(p);
                if (theFormat < 0) {
                    (*string_form) << "REWIND";
                } else {
                    (*string_form) << *((_String*)_ElementaryCommand::fscanf_allowed_formats.GetItem (theFormat));
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


        case 38: { // reconsruct ancestors
            (*string_form) << assignment (HY_HBL_COMMAND_DATA_SET, "ReconstuctAncestors");

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
            
        case HY_HBL_COMMAND_INIT_ITERATOR: {
            (*string_form) << "Initialize iterator on " << parameter_to_string (0) << ";";
            break;
        }

        case HY_HBL_COMMAND_ADVANCE_ITERATOR: {
            (*string_form) << "Advance iterator into ";
            (string_form -> AppendNewInstance(parameters.Join(",",0,simpleParameters.get (1)))) << ";";
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

      if (chain.is_compiled(chain.currentCommand)) {
          //if (chain.is_compiled(0) == false) {
          //    PopulateArraysForASimpleFormula (chain.cli->varList, chain.cli->values);
          //}
          hyFloat result = ((_Formula*)simpleParameters.list_data[1])->ComputeSimple (chain.cli->stack, chain.cli->values);
          long sti = chain.cli->storeResults.list_data[chain.currentCommand-1];
          if (sti>=0) {
              chain.cli->values[sti].value = result;
          }
          return;
      }

      if (!simpleParameters.lLength) { // not compiled yet
          _Formula f,
                   f2;

          _String* theFla     = (_String*)parameters(0),
                   err_msg;

          _FormulaParsingContext fpc (&err_msg, chain.nameSpacePrefix);

          long     parseCode = Parse(&f,(*theFla),fpc,&f2);
          
          //printf ("RHS = %s\n", _String ((_String*)f2.toStr(kFormulaStringConversionNormal)).get_str());

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
             errMsg = new _String (_String ("Parsing error ") & _String (err_msg.Enquote('(',')') & " while compiling the statement: "));
            throw 0;
          }
      }

      ExecuteFormula ((_Formula*)simpleParameters.list_data[1],(_Formula*)simpleParameters.list_data[2],simpleParameters.list_data[0],simpleParameters.list_data[3], chain.nameSpacePrefix, simpleParameters.list_data[4]);

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

          if (chain.is_compiled(chain.currentCommand)) {
              //if (chain.is_compiled(0) == false) {
              //    PopulateArraysForASimpleFormula (chain.cli->varList, chain.cli->values);
              //}
              if ( ((_Formula*)simpleParameters(2))->ComputeSimple(chain.cli->stack, chain.cli->values)==0.0) {
                  chain.currentCommand = simpleParameters.list_data[1];
                  return;
              }
          } else {
              HBLObjectRef result;
              if (expression) {
                  //printf ("\n*** Interpreted condition\n");
                result = expression->Compute();
              } else {
                  //printf ("\n*** Compiled condition\n");
                result = ((_Formula*)simpleParameters(2))->Compute();
              }

              // printf ("\n*** %s\n", ((_String*)result->toStr())->sData);

            if (terminate_execution && !result) {
                  _String       *s = (_String*)((_Formula*)simpleParameters(2))->toStr(kFormulaStringConversionSubstiteValues);
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
                  chain.currentCommand = simpleParameters.list_data[1];
                  return;
              }
          }
      }
      chain.currentCommand = simpleParameters.list_data[0];

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
        if (fName == useNexusFileData) {
            if (!lastNexusDataMatrix) {
                HandleApplicationError (useNexusFileData & " was used in ReadDataFile, and no NEXUS data matrix was available.");
                return;
            }
            ds = lastNexusDataMatrix;
        } else {
            ProcessFileName(fName, false,true,(hyPointer)chain.nameSpacePrefix, false, &chain, true);
            if (terminate_execution) {
                return;
            }
            SetStatusLine ("Loading Data");

            df = doFileOpen (fName.get_str(),"rb");
            if (df==nil) {
                // try reading this file as a string formula
                fName = GetStringFromFormula ((_String*)parameters(1),chain.nameSpacePrefix);
                ProcessFileName(fName, false,false,(hyPointer)chain.nameSpacePrefix, false, &chain, true);

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
                        theFreqID   = modelFrequenciesIndices.list_data[theModelID];

                        while((thisNode = ti.Next()) && !ti.IsAtRoot()) {
                            theModelID      = thisNode->GetModelIndex();
                            if (theModelID == HY_NO_MODEL) { // no model
                                done = false;
                                break;
                            }
                            if (modelFrequenciesIndices.list_data[theModelID]!=theFreqID) {
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
                likeFuncList.list_data[likeFuncObjectID] = (long)lkf;
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
              errMsg;

 
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
            // ';'-sep for different partititons
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
            catValVar->SetValue(catValues,false,true,NULL);
        }
        if (catNames) {
            catNameVar->SetValue(catNames,false,true,NULL);
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
            throw _String("The alphabet is mis-specified; it either has redundant characters or multi-character/non-string entries");
        }

        long unit_size = ((_FString*)alphabet_matrix->GetFormula(1,0)->Compute())->get_str().to_long();

        if (unit_size < 1L) {
            throw _String("The evolutionary unit size in the alphabet matrix is mis-specified");
        }

        _Formula* exclusion_formula = alphabet_matrix->GetFormula(1,1);
        _String const* the_exclusions = &kEmptyString;

        if (exclusion_formula) {
            the_exclusions = &((_FString*)exclusion_formula->Compute())->get_str();
        }

        _TheTree * spawning_tree = (_TheTree*)tree_var;

        if (parameters.lLength>6 && (spawning_tree->CountTreeCategories()>1)) {
            throw _String("Can't use spool to file option in Simulate when the tree depends on category variables.");
        }

        if (given_state.length()>1) {
        // root state
            if ((given_state.length() >= unit_size)&&(given_state.length() % unit_size == 0)) {
                site_count = given_state.length()/unit_size;
            } else {
                throw _String("Root state string is either too short or has length which is not divisible by the unit size");
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
        InitializeArray (the_map->list_data, total_sites, 0L);
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


        bool    do_internals = parameters.countitems() > 5 ? (ProcessNumericArgument ((_String*)parameters (5),chain.nameSpacePrefix)>0.5) : false;

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
        { // lf must be deleted before the referenced datafilters
            _LikelihoodFunction lf (filter_specification, nil);
            
            /*_SimpleList gl;
            lf.GetGlobalVars(gl);
            gl.Each ([] (long vi, unsigned long) -> void { StringToConsole(*LocateVar(vi)->GetName()); NLToConsole();});
            */
            
            lf.Simulate (*sim_dataset, exclusions, category_values, category_names, root_states, do_internals?(main_file?&spool_file:&kEmptyString):nil);
            SetStatusLine ("Idle");
        }
        
        
        category_values_id->SetValue(category_values, false,true,NULL);
        category_names_id->SetValue(category_names, false,true,NULL);

        
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

        
        bool replace_tree_structure = hy_env::EnvVariableTrue(replaceTreeStructure);

        _SimpleList   leftOverVars; // mod 02/03/2003
        if (varID>=0)
            if (FetchVar(varID)->ObjectClass()==TREE) {
                if (replace_tree_structure) {
                    DeleteVariable(*FetchVar(varID)->GetName());    // mod 11/19/2003
                } else {
                    DeleteTreeVariable(varID,leftOverVars);    // mod 02/03/2003
                }
            }


        _TheTree * tr = nil;

        if (treeString.get_char(0)!='(') {
            _Formula  nameForm (treeString,chain.nameSpacePrefix);
            HBLObjectRef formRes = nameForm.Compute();
            if (formRes) {
                if (formRes->ObjectClass () == STRING) {
                    tr = new _TheTree (treeIdent,((_FString*)formRes)->get_str(),false);
                } else if (formRes->ObjectClass () == TOPOLOGY) {
                    tr = new _TheTree (treeIdent,(_TreeTopology*)formRes);
                } else if (formRes->ObjectClass () == TREE) {
                    for (unsigned long i = 0; i < leftOverVars.lLength; i++) {
                        //rintf ("%s\n", LocateVar(leftOverVars.list_data[i])->GetName()->get_str());
                        DeleteVariable(leftOverVars.list_data[i], true);
                    }
                    leftOverVars.Clear();
                    tr = new _TheTree (treeIdent,(_TheTree*)formRes);
                }
            }
        } else {
            tr = new _TheTree (treeIdent,treeString,false);
        }

        if (!tr) {
            DeleteObject (tr);
            HandleApplicationError("Illegal right hand side in call to Tree id = ...; it must be a string, a Newick tree spec or a topology");
            return false;
        }

        if (leftOverVars.nonempty()) { // mod 02/03/2003 - the entire "if" block
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
                DeleteVariable (*LocateVar (holder.list_data[varID])->GetName());
            }

            tr->Clear();
 
        }

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
          
        /*if (chain.is_compiled()) {
          chain.CopyCLIToVariables();
        }*/

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

          HBLObjectRef ret_val = nil;
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
                if (dsIndex.list_data[di] != newSetID) {
                    KillDataSetRecord(dsIndex.list_data[di]);
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
            result->SetValue(storage,false,true,NULL);
            lastMatrixDeclared = result->get_index();
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

    case HY_HBL_COMMAND_KEYWORD_ARGUMENT:
        HandleKeywordArgument (chain);
        break;
          
    case HY_HBL_COMMAND_SET_PARAMETER:
        return HandleSetParameter(chain);

    case 38: // reconstruct ancestors
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

      case HY_HBL_COMMAND_GET_INFORMATION:
          return HandleGetInformation(chain);
          break;

      case HY_HBL_COMMAND_REPLICATE_CONSTRAINT:
          return HandleReplicateConstraint (chain);
          break;

    case HY_HBL_COMMAND_INTEGRATE:
    case HY_HBL_COMMAND_FIND_ROOT:
        return HandleFindRootOrIntegrate(chain, code == HY_HBL_COMMAND_INTEGRATE);
        break;


    case HY_HBL_COMMAND_MPI_SEND:
        return HandleMPISend (chain);
        break;

    case HY_HBL_COMMAND_MPI_RECEIVE:
        HandleMPIReceive(chain);
        break;


    case 47: // state counter; deprecate
        ExecuteCase47 (chain);
        break;

    case HY_HBL_COMMAND_LFCOMPUTE:
        return HandleComputeLFFunction(chain);

    case 50: // sample ancestors
        ExecuteCase38 (chain, true);
        break;

    case HY_HBL_COMMAND_GET_URL:
        return HandleGetURL (chain);

    case 52: // Simulate
        ExecuteCase52 (chain);
        break;

    case HY_HBL_COMMAND_DO_SQL:
        return HandleDoSQL(chain);

    case 54: // topology id =
        ExecuteCase54 (chain);
        break;

    case HY_HBL_COMMAND_ALIGN_SEQUENCES:
        return HandleAlignSequences (chain);
        break;


      case 58: // #profile
        ExecuteCase58 (chain);
        break;

    case HY_HBL_COMMAND_DELETE_OBJECT:
        return HandleDeleteObject (chain);
        break;

    case HY_HBL_COMMAND_REQUIRE_VERSION:
        HandleRequireVersion(chain);
        break;

    case 61: // SCFG =
        ExecuteCase61 (chain);
        break;

    case 63: // dead?
        ExecuteCase63 (chain);
        break;

    case 64: // BGM =
        ExecuteCase64 (chain);
        break;

    case HY_HBL_COMMAND_ASSERT:
        HandleAssert (chain);
        break;

    case HY_HBL_COMMAND_GET_DATA_INFO:
        HandleGetDataInfo(chain);
        break;
          
    case HY_HBL_COMMAND_INIT_ITERATOR:
        HandleInitializeIterator(chain);
        break;

    case HY_HBL_COMMAND_ADVANCE_ITERATOR:
        HandleAdvanceIterator(chain);
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


void   _ElementaryCommand::FindNextCommand  (_String& input, _StringBuffer &result) {

    long    index     = input.length();
    result.Reset();

    if (index == 0L) {
      return;
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


    char    last_char = '\0';
        // a look back character


    // non printable characters at the end ?
    while (index>0) {
      if (!isprint (input.char_at (index))) {
        index--;
      } else {
        break;
      }
    }
  
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
          // SLKP 20170704: this seems incomplete : need to check more thorougly that this is an ident
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
                result.Reset();
                return ;
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
                result.Reset();
                return ;
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


    if (scope_depth != 0L || comment_state == slash_star || literal_state != normal_text || matrix_depth != 0L || bracket_depth != 0L || parentheses_depth != 0L) {
        if (result!='}') {
            HandleApplicationError (_String("Expression appears to be incomplete/syntax error. {} scope: ") &scope_depth & ", () depth "
                       & parentheses_depth & ", matrix scope: " & matrix_depth & '.' & (literal_state == double_quote ?" In a \"\" literal. ":kEmptyString)
                       & (literal_state == single_quote?" In a '' literal. ":kEmptyString) &
                       (comment_state == slash_star ? " In a /* */ comment ":kEmptyString) & '\n' & input);
            input.Clear();
            result.Reset();
            return ;
        } else {
            result.Reset();
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

}
//____________________________________________________________________________________

long _ElementaryCommand::ExtractConditions (_String const& source, long start_at, _List& receptacle, char delimeter, bool include_empty_conditions) {

    long parentheses_depth = 1L,
         // this is because extaction will work from the first character following a '(', e.g. CreateFilter([start parsing here]....)
         last_delim    = start_at,
         index             = start_at,
         curly_depth       = 0L,
         bracket_depth     = 0L;


    enum {
        normal_text = 0,
        single_quote = 1,
        double_quote = 2
    } quote_type = normal_text;

    auto strip_last_space = [] (_String const& source, long from, long to) -> _String* {
        if (to > from) {
            if (source.char_at(to-1L) == ' ') {
                return new _String (source, from, to-2L);
            }
        }
        return new _String (source, from, to-1L);
    };

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
            if (c=='[') {
                bracket_depth++;
                continue;
            }
            if (c==']') {
                bracket_depth--;
                continue;
            }
            if (c==')') {
                parentheses_depth --;
                if (parentheses_depth == 0L) {
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
        if (c==delimeter) {
            if (parentheses_depth > 1 || quote_type != normal_text || curly_depth || bracket_depth) {
                continue;
            }

            receptacle < strip_last_space (source,last_delim,index);
            last_delim = index+1;
            continue;
        }
    }

    if (include_empty_conditions || last_delim <= index-1) {
        receptacle < strip_last_space (source,last_delim,index);
    }
    return index+1L;
}

//____________________________________________________________________________________


bool       _ElementaryCommand::MakeGeneralizedLoop  (_String*p1, _String*p2, _String*p3 , bool for_or_while, _String& source, _ExecutionList&target) {

    
    const _String kIterator ("in");
    
    // extract the for enclosure
    long  beginning = target.lLength,
          for_return = beginning;

    bool   success = true;
    bool   has_increment = false;
    bool   is_iterator = (p2 && *p2 == kIterator);

    _SimpleList bc;
    
    try {
        if (is_iterator) {
            if (! (p1 && p3)) {
                throw (kEmptyString);
            }
            for_return++;
            has_increment = true;
            _ElementaryCommand  * init    = new _ElementaryCommand (HY_HBL_COMMAND_INIT_ITERATOR);
                               
            long refer_to_init = target.countitems();
            
            init->parameters    << p3;
            target < init;
            _ElementaryCommand * advance = new _ElementaryCommand (HY_HBL_COMMAND_ADVANCE_ITERATOR);
            
            ExtractConditions (*p1,0L,advance->parameters,',',false);
            if (advance->parameters.countitems() < 1 || advance->parameters.countitems() > 3) {
                throw (_String("Must have bewteen 1 and 3 arguments for the 'for' iterator loop"));
            }
            
            advance->simpleParameters <<refer_to_init;
            advance->simpleParameters << advance->parameters.countitems();
            init->simpleParameters << target.countitems();
            init->simpleParameters << 0;
            target < advance;
            
            for_return++;
            target < new _ElementaryCommand ();

            if (source.get_char(0)=='{') {
                source.Trim(1,kStringEnd);
            }

            if (target.BuildList (source, &bc) == false) { // construct the main body
                throw (kEmptyString);
            }
            target << advance;
            
            _ElementaryCommand * loopback = new _ElementaryCommand;
            loopback->MakeJumpCommand (nil,for_return,0,target);
            target < loopback;
            
            // automatically generate
            // None != iterator_value, but we don't know
            
            _StringBuffer iterator_condition;
            iterator_condition << kEndIteration.Enquote() << "!=" << (_String*)advance->parameters.GetItem (advance->parameters.countitems() -1);
            //StringToConsole(iterator_condition);NLToConsole();
            target.GetIthCommand(for_return)->MakeJumpCommand (&iterator_condition, for_return+1, target.lLength,target);
            
        } else {
        
            if (p1 && p1->nonempty()) { // initialization stage
                for_return++;
                success = success && target.BuildList (*p1, nil, true); // add init step
            }

            // append condition now

            if (!success) {
                throw (kEmptyString);
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
                throw (kEmptyString);
            }

            if (p3 && p3->nonempty ()) { // increment stage
                success = success && target.BuildList (*p3, nil,true); // add increment step
                has_increment = true;
            }

            if (!success) {
                throw (kEmptyString);
            }
            
            if (for_or_while) {
                _ElementaryCommand * loopback = new _ElementaryCommand;
                loopback->MakeJumpCommand (nil,for_return,0,target);
                target < loopback;
                if (p2 && p2->nonempty()) {
                    target.GetIthCommand(for_return)->MakeJumpCommand (p2, for_return+1, target.lLength,target);
                }
            } else {
                if (p2) {
                    _ElementaryCommand* loopback = new _ElementaryCommand ;
                    success = success && loopback->MakeJumpCommand (p2,for_return,target.lLength+1,target);
                    target < loopback;
                }
            }
        }
    } catch (const _String& err) {
        for (long index = target.lLength - 1; index >= beginning; index--) {
            target.Delete (index);
        }
        if (err.nonempty()) {
            HandleApplicationError(err);
        }
        return false;
    }

    bc.Each ([&target, has_increment] (long loc, unsigned long) -> void {
        if (loc>0) { // break
            (target.GetIthCommand(loc))->MakeJumpCommand (nil, target.lLength, 0,target);
        } else { // continue
            (target.GetIthCommand(-loc))->MakeJumpCommand (nil, target.lLength-(has_increment?2:1), 0,target);
        }
    });
    
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

        _StringBuffer nextCommand;
        _ElementaryCommand::FindNextCommand(source,nextCommand);
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

bool    _ElementaryCommand::ProcessInclude      (_String&source, _ExecutionList&target) {
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
    
    const _String kConcat  ("Concatenate"),
                  kCombine ("Combine"),
                  kReadDataFile ("ReadDataFile"),
                  kReadFromString ("ReadFromString"),
                  kPurge ("purge"),
                  kReconstructAncestors ("ReconstructAncestors"),
                  kSampleAncestors ("SampleAncestors"),
                  kMarginalAncestors               ("MARGINAL"),
                  kDoLeavesAncestors               ("DOLEAVES"),
                  kSimulate ("Simulate");
    
    try {

        _String operation_type;
        _List   arguments;
        
        ProcessProcedureCall (source, operation_type, arguments);

        if (operation_type ==  kReadDataFile || operation_type == kReadFromString) {
            if (arguments.countitems () != 2UL) {
                throw _String ("DataSet declaration missing a valid filename/string or has extra arguments");
            }

            _ElementaryCommand * dsc = new _ElementaryCommand (5);
            
            if (operation_type == kReadFromString) {
                dsc->simpleParameters << 1;
            }
            dsc->addAndClean (target, &arguments, 0L);
        } else if (operation_type == blSimulateDataSet) {
            if ( arguments.countitems()>5UL || arguments.countitems()==1UL ) {
                throw blSimulateDataSet.Enquote() & "expects 1-4 parameters: likelihood function ident (needed), a list of excluded states, a matrix to store random rates in, and a matrix to store the order of random rates in (last 3 - optional).";
            }

            _ElementaryCommand * dsc = new _ElementaryCommand (12);
            dsc->addAndClean (target, &arguments, 0L);
        } else if ( operation_type ==  kConcat || operation_type ==  kCombine) {
            _ElementaryCommand * dsc = new _ElementaryCommand (16);
            dsc->simpleParameters<<((operation_type==kConcat)?1:2);

            if ((*(_String*)arguments.GetItem(1)) == kPurge) {
                dsc->simpleParameters[0] = - dsc->simpleParameters[0];
                arguments.Delete (1);
            }

            if (arguments.countitems() == 1UL) {
                delete (dsc);
                throw _String ("DataSet merging operation missing a valid list of arguments.");
            }
            dsc->addAndClean (target, &arguments, 0L);
            return true;

        } else {
            if (operation_type ==  kReconstructAncestors || operation_type == kSampleAncestors) {
                if (arguments.countitems()>5UL || arguments.countitems()==1L) {
                    throw  operation_type.Enquote() & " expects 1-4 parameters: likelihood function ident (mandatory), an matrix expression to specify the list of partition(s) to reconstruct/sample from (optional), and, for ReconstructAncestors, an optional MARGINAL flag, plus an optional DOLEAVES flag.";
                }
                _ElementaryCommand * dsc = new _ElementaryCommand (operation_type ==  kReconstructAncestors ? 38 : 50);
                dsc->parameters << arguments (0) << arguments (1);
               for (long optP = 2L; optP < arguments.lLength; optP++) {
                    _String * current_term = (_String*)arguments.GetItem(optP);
                    
                    if (*current_term == kMarginalAncestors) {
                        dsc->simpleParameters << -1;
                    } else if (*current_term == kDoLeavesAncestors) {
                        dsc->simpleParameters << -2;
                    } else {
                        dsc->parameters  << current_term;
                    }
                }

                dsc->addAndClean (target);
                return true;
            } else if (operation_type ==  kSimulate) {
                if ((arguments.countitems()>8)||(arguments.countitems()<5UL)) {
                    throw kSimulate.Enquote() & " expects 4-6 parameters: tree with attached models, equilibrium frequencies, character map, number of sites|root sequence, <save internal node sequences>, <file name for direct storage>";
                    
                 }

                _ElementaryCommand * dsc = new _ElementaryCommand (52);
                dsc->addAndClean (target, &arguments, 0);
                return true;
            } else {
                throw _String ("Expected DataSet ident = ReadDataFile(filename); or DataSet ident = SimulateDataSet (LikelihoodFunction); or DataSet ident = Combine (list of DataSets); or DataSet ident = Concatenate (list of DataSets); or DataSet ident = ReconstructAnscetors (likelihood function); or DataSet ident = SampleAnscetors (likelihood function) or DataSet	  dataSetid = ReadFromString (string);");
            }
        }
    } catch (const _String& err) {
        HandleErrorWhileParsing (err, source);
        return false;
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
            _String definition (source,mark1+1,mark2-1);
            _List args;
            ExtractConditions (definition,0,args,',');
            if (args.countitems()>=7UL) {
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
    //((_String*)args.list_data[1])->StripQuotes();
    cv->parameters<<args(1)
                  <<args(2)
                  <<args(3);

    if  (args.lLength>5UL) {
        _List * choices = new _List;
        for (long k = 4L; k<args.lLength-1; k+=2) {
            ((_String*)args.list_data[k])->StripQuotes();
            ((_String*)args.list_data[k+1])->StripQuotes();
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

bool    _ElementaryCommand::ConstructDataSetFilter (_String&source, _ExecutionList&target) {
// DataSetFilter      dataSetFilterid = CreateFilter (datasetid;unit;vertical partition; horizontal partition; alphabet exclusions);
    // first we must segment out the data set name
    
    const _String kCreateFilter ("CreateFilter"),
                  kPermute ("Permute"),
                  kBootstrap ("Bootstrap");
    
    _ElementaryCommand * datafilter_command = nil;

    try {

        _String operation_type;
        _List   arguments;
        
        ProcessProcedureCall (source, operation_type, arguments);

        if (operation_type == kCreateFilter) {
            datafilter_command = new _ElementaryCommand(6);
        } else if (operation_type == kPermute) {
            datafilter_command = new _ElementaryCommand(27);
        } else if (operation_type == kBootstrap) {
            datafilter_command = new _ElementaryCommand(28);
        } else {
            throw _String ("Expected: DataSetFilter	  dataSetFilterid = CreateFilter (datasetid,unit,vertical partition,horizontal partition,alphabet exclusions); or Permute/Bootstrap (dataset/filter,<atom>,<column partition>)");
        }

        if (!(arguments.countitems()>=3UL || (arguments.countitems() == 2UL && datafilter_command->code == 6))) {
            throw _String ("Parameter(s) missing in DataSetFilter definition.");
        }

        datafilter_command->addAndClean (target,&arguments);
    } catch (const _String& err) {
        HandleErrorWhileParsing (err, source);
        DeleteObject (datafilter_command);
        return false;
    }
    
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

    if (mark1==-1 || mark2==-1 || !modelID.IsValidIdentifier(fIDAllowCompound|fIDAllowFirstNumeric)) {
        HandleApplicationError ("Model declaration missing a valid identifier.");
        return false;
    }

    // now look for the opening paren
    mark1 = source.Find ('(',mark2,-1);
    _List pieces;
    ExtractConditions (source,mark1+1,pieces,',');

    if (pieces.lLength<2) {
        HandleApplicationError ("Parameter(s) missing in Model definition. Must have a matrix and a compatible eqiulibrium frequencies vector.");
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

bool      _ElementaryCommand::MakeJumpCommand       (_String* source,   long branch_true, long branch_false, _ExecutionList& parentList) {
    long existing_formula_handle = 0L;
         code        = 4L;

    if (simpleParameters.lLength==3) {
        if (source) {
            delete ((_Formula*)simpleParameters.get(2));
         } else {
            existing_formula_handle = simpleParameters.get(2);
        }
    }

    if (branch_true==-1) {
        if (simpleParameters.empty()) {
            HandleApplicationError("An if-then-else scoping error. Check opening and closing brackets and double else's.");
            return false;
        }
        branch_true = simpleParameters.get (0);
    }

    simpleParameters.Clear();
    simpleParameters<<branch_true<<branch_false;
    if (source) {
        parameters && source;
    } else if (existing_formula_handle) {
        simpleParameters<<existing_formula_handle;
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
bool    _ElementaryCommand::ConstructFunction (_String&source, _ExecutionList& chain) {
// syntax: function <ident> (comma separated list of parameters) {body}


    bool    isFFunction = source.BeginsWith (blFFunction),
            isLFunction = ! isFFunction && source.BeginsWith (blLFunction),
            isCFunction = ! isFFunction && ! isLFunction && source.BeginsWith (blCFunction),
            isNameSpace = ! isFFunction && ! isLFunction && ! isCFunction && source.BeginsWith (blNameSpace);

    _hy_nested_check save_state = isInFunction;
    
    if (!isNameSpace) {
      if (isInFunction == _HY_FUNCTION) {
        HandleApplicationError ("Nested function declarations are not allowed");
        return false;
      }
    }
    
    long mark1, mark2;
    
    if (isNameSpace) {
        mark1 = source.FirstNonSpaceIndex(blNameSpace.length(), kStringEnd, kStringDirectionForward);
        mark2 = source.Find ('{', mark1, kStringEnd);
    } else {
        _String const * prefix = &blFunction;
        
        if (isFFunction) prefix = &blFFunction; else
        if (isLFunction) prefix = &blLFunction; else
        if (isCFunction) prefix = &blCFunction;

        mark1 = source.FirstNonSpaceIndex(prefix->length(),kStringEnd,kStringDirectionForward);
        mark2 = source.Find ('(', mark1, kStringEnd);
    }
    

    if ( mark1 == kNotFound || mark2 == kNotFound || mark1 > mark2-1) {
        HandleApplicationError      (_String("Function declaration missing a valid function identifier or parameter list.\n-----------\n") & source & "\n-----------\n");
        return false;
    }

    _String*    funcID  = new _String(source,mark1,mark2-1);

    if (!funcID->IsValidIdentifier(fIDAllowCompound)) {
      HandleApplicationError      (_String("Not a valid function/namespace identifier '") & _String(funcID) & "'");
      return false;
    }

    *funcID = chain.AddNameSpaceToID (*funcID);

    // now look for the opening paren

    if (!isNameSpace) {


      _String extraNamespace;
        
      if ((mark1=FindBFFunctionName(*funcID)) >= 0L) {
          ReportWarning (_String("Overwritten previously defined function:'") & *funcID & '\'');
          if (batchLanguageFunctionClassification.get (mark1) == kBLFunctionLocal) {
              // clean up previously used namespaces
              _ExecutionList * existing = (_ExecutionList *)batchLanguageFunctions.GetItem(mark1);
              _String * nm = existing->GetNameSpace();
              if (nm) {
                  long best_index;
                  char match = variableNames.FindBest(nm, best_index);
                  //printf ("PREVIOUS LOCAL FUNCTION %s\n", nm->get_str());
                  if (best_index >= 0) {
                      _Variable * namespace_start = LocateVar (best_index);
                      if (namespace_start->GetName()->BeginsWith(*nm)) {
                          //printf ("Clearing namespace %s starting at %s for function %s \n", nm->get_str(), namespace_start->GetName()->get_str(), funcID->get_str());
                          _SimpleList locals;
                          DeleteTreeVariable (best_index, locals,nm, false);
                      }
                  }
                  if (match == 0) {
                      DeleteVariable(best_index, true, false);
                  }
              }
          }
        
      }

      _List       arguments;
      _SimpleList argument_types;

      long upto = ExtractConditions (source,mark2+1,arguments,',',false);


      if (upto==source.length() || source[upto]!='{' || source (-1)!='}') {
          HandleApplicationError (_String("Function declaration is missing a valid function body."));
          isInFunction= save_state;
          return false;
      }

      if (isLFunction && extraNamespace.empty())
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

      isInFunction = _HY_FUNCTION;
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
        
      if (isCFunction) {
          if (functionBody->TryToMakeSimple()) {
              ReportWarning(_String ("Successfully compiled code for function ") & funcID->Enquote());
          }
      }

      //  take care of all the return statements
      returnlist.Each ([functionBody] (long value, unsigned long) -> void {
        ((_ElementaryCommand*)functionBody->GetItem(value))->simpleParameters << functionBody->countitems();
      });
      returnlist.Clear();
        
      if (mark1>=0) {
          batchLanguageFunctions.Replace (mark1, functionBody, false);
          batchLanguageFunctionNames.Replace (mark1, funcID, false);
          batchLanguageFunctionParameterLists.Replace (mark1, &arguments, true);
          batchLanguageFunctionParameterTypes.Replace (mark1, &argument_types, true);
          batchLanguageFunctionClassification.list_data[mark1] = isLFunction ? kBLFunctionLocal :( isFFunction? kBLFunctionSkipUpdate :  kBLFunctionAlwaysUpdate);
      } else {
          batchLanguageFunctions.AppendNewInstance(functionBody);
          batchLanguageFunctionNamesIndexed.Insert (new _String (*funcID), batchLanguageFunctions.countitems() - 1, false, true);
          batchLanguageFunctionNames.AppendNewInstance(funcID);
          batchLanguageFunctionParameterLists &&(&arguments);
          batchLanguageFunctionParameterTypes &&(&argument_types);
          batchLanguageFunctionClassification <<(isLFunction ? kBLFunctionLocal :( isFFunction? kBLFunctionSkipUpdate :  kBLFunctionAlwaysUpdate));
      }
    } else {
      if (mark2 == source.length () || source[mark2]!='{' || source (-1L) !='}') {
        HandleApplicationError (_String("Namespace declaration is missing a body."));
        return false;
      }
      _String          namespace_text (source, mark2+1,source.length()-2);
      bool             success = false;

      isInFunction = _HY_NAMESPACE;
      _ExecutionList   * namespace_payload = new _ExecutionList (namespace_text, funcID, false, &success);
      DeleteObject (funcID);
        // 20180713 SLKP -- this was marked as deleted in one of the v2.3 branches
      if (success) {
        _ElementaryCommand * nested_list = new _ElementaryCommand (HY_HBL_COMMAND_NESTED_LIST);
        nested_list->parameters.AppendNewInstance(namespace_payload);
        chain.AppendNewInstance(nested_list);
      } else {
        DeleteObject (namespace_payload);
        return false;
      }

    }


    isInFunction = save_state;
    return true;
}

//____________________________________________________________________________________
bool    _ElementaryCommand::ConstructReturn (_String&source, _ExecutionList&target) {
// syntax: return <statement>

    long    mark1 = source.FirstNonSpaceIndex(blReturn.length(),kStringEnd,kStringDirectionForward);

    _ElementaryCommand * return_statement = new _ElementaryCommand (14);

    if (mark1 != kNotFound) { // not a trivial return statement;
        if (source (-1) ==';') {
            return_statement->parameters.AppendNewInstance (new _String (source, mark1, source.length () - 2L));
        } else {
            return_statement->parameters.AppendNewInstance (new _String (source, mark1, kStringEnd));
        }
    }

    if (isInFunction == _HY_FUNCTION) {
        returnlist<<target.countitems();
    } else {
        if (isInFunction == _HY_NAMESPACE) {
            HandleApplicationError("return statements are not allowed in namespaces");
        }
        return_statement->simpleParameters << -1;
    }
                                                            
    return_statement->addAndClean (target);
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
            ReadDataSetFile (f,1,nil,&fName, nil, &hy_default_translation_table, &target);
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

    if (modelTypeList.list_data[theModel]) {
        theExp = (_Formula*)modelMatrixIndices.list_data[theModel];
        theExp->ScanFForType(matrices, MATRIX);

        for (long mi = 0; mi < matrices.countitems(); mi++) {
            if (alreadyDone && alreadyDone->Insert ((BaseRef)matrices.list_data[mi]) < 0) {
                matrices.Delete(mi);
                mi--;
            }
        }
    } else {
        if (!alreadyDone || alreadyDone->Find ((BaseRef)modelMatrixIndices.list_data[theModel]) < 0) {
            if (alreadyDone) {
                alreadyDone->Insert ((BaseRef)modelMatrixIndices.list_data[theModel]);
            }
            matrices << modelMatrixIndices.list_data[theModel];
        }
        tV = LocateVar(modelMatrixIndices.list_data[theModel]);
    }

    long freqID = modelFrequenciesIndices.list_data[theModel];

    if (freqID>=0) {
        tV2 = LocateVar(freqID);
    } else {
        mByF = false;
        tV2 = LocateVar(-freqID-1);
    }

    if (!alreadyDone || alreadyDone->Find ((BaseRef)tV2->get_index()) < 0) {
        if (alreadyDone) {
            alreadyDone->Insert ((BaseRef)tV2->get_index());
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
            LocateVar(matrices.list_data[mi])->ScanForVariables (vlst,true);
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
            _Variable *tV = LocateVar (matrices.list_data[k]);
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
