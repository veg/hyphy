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

#include      "baseobj.h"
#include      "defines.h"
#include      "batchlan.h"
#include      "likefunc.h"
#include      "bayesgraph.h"
#include      "scfg.h"

#if defined __MAC__ || defined __WINDOZE__ || defined __HYPHY_GTK__
    #include "HYConsoleWindow.h"
    #include "HYDialogs.h"

#endif

#if defined __HYPHYQT__
#include "HYSharedMain.h"
#include "hyphy_qt_helpers.h"
#endif

_List       openFileHandlesBackend;

_AVLListX   openFileHandles     (&openFileHandlesBackend);

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleHarvestFrequencies (_ExecutionList& currentProgram) {
    currentProgram.currentCommand++;
    

    _String  freqStorageID = *(_String*)parameters(0),
             dataID        = currentProgram.AddNameSpaceToID(*(_String*)parameters(1)), 
             errMsg        ;
    
    _Variable * theReceptacle = CheckReceptacleCommandID (&AppendContainerName(freqStorageID,currentProgram.nameSpacePrefix),HY_HBL_COMMAND_HARVEST_FREQUENCIES, true, false, &currentProgram);
     if (!theReceptacle) {
        return false;
    }    
    SetStatusLine           ("Gathering Frequencies");


    long       objectType = HY_BL_DATASET|HY_BL_DATASET_FILTER;
    BaseRef    sourceObject = _HYRetrieveBLObjectByName (dataID, objectType,nil,false);
    
    long      unit      = ProcessNumericArgument((_String*)parameters(2),currentProgram.nameSpacePrefix),
              posspec   = ProcessNumericArgument((_String*)parameters(4),currentProgram.nameSpacePrefix),
              atom      = ProcessNumericArgument((_String*)parameters(3),currentProgram.nameSpacePrefix);
    
    _Matrix*            receptacle = nil;
    
    _Parameter          cghf = 1.0;
    checkParameter      (hfCountGap,cghf,1.0, currentProgram.nameSpacePrefix);
    
    if (objectType == HY_BL_DATASET) { // harvest directly from a DataSet
        _String         vSpecs,
                        hSpecs;
        
        if (parameters.lLength>5)   {
            vSpecs = *(_String*)parameters(5);
        }
        if (parameters.lLength>6)   {
            hSpecs = *(_String*)parameters(6);
        }
        
        _DataSet * dataset = (_DataSet*)sourceObject;
        _SimpleList     hL, vL;
        dataset->ProcessPartition (hSpecs,hL,false);
        dataset->ProcessPartition (vSpecs,vL,true);
        
        receptacle = dataset->HarvestFrequencies(unit,atom,posspec,hL, vL,cghf>0.5);
    } else { // harvest from a DataSetFilter
        if (objectType == HY_BL_DATASET_FILTER) {
            receptacle = ((_DataSetFilter*)sourceObject)->HarvestFrequencies(unit,atom,posspec,cghf>0.5);
        } else {
            errMsg = _String ("'") & dataID & "' is neither a DataSet nor a DataSetFilter";
        }
    }
    
    SetStatusLine           (empty);
    
    if (errMsg.sLength || receptacle == nil) {
        DeleteObject (receptacle);
        currentProgram.ReportAnExecutionError (errMsg); 
        theReceptacle->SetValue (new _MathObject, false);
        return false;
    }   
    theReceptacle->SetValue (receptacle, false);
    return true;
    
    //CheckReceptacleCommandIDAndStore (&freqStorageID,HY_HBL_COMMAND_HARVEST_FREQUENCIES,true, receptacle, false);
}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleOptimizeCovarianceMatrix (_ExecutionList& currentProgram, bool doOptimize) {
    currentProgram.currentCommand++;
    // construct the receptacle matrix
    
    _String  lfResName  (currentProgram.AddNameSpaceToID(*(_String*)parameters(0))),
             lfNameID   (currentProgram.AddNameSpaceToID(*(_String*)parameters(1)));
             
    _Variable* result = CheckReceptacleCommandID (&lfResName, doOptimize?HY_HBL_COMMAND_OPTIMIZE:HY_HBL_COMMAND_COVARIANCE_MATRIX,true);

    // Handle string variables passed as likefunc IDs?
     _String temp = ProcessLiteralArgument(&lfNameID,currentProgram.nameSpacePrefix);
    if (temp.sLength) {
        lfNameID = temp;
    }

    long       objectType = HY_BL_LIKELIHOOD_FUNCTION|HY_BL_SCFG|HY_BL_BGM;
    _LikelihoodFunction    *lkf = (_LikelihoodFunction*)_HYRetrieveBLObjectByName (lfNameID, objectType,nil,doOptimize==false);
    
    if (lkf == nil) { // will only happen if the object is a custom function 
        lkf = (_LikelihoodFunction*)checkPointer(new _CustomFunction (&lfNameID));
    }
    
    if (!doOptimize) { 
        // COVARIANCE_MATRIX 
 
        SetStatusLine (_String("Finding the cov. matrix/profile CI for ")&lfNameID);
        _String              cpl = currentProgram.AddNameSpaceToID(covarianceParameterList);
        _Variable            *  restrictVariable = FetchVar (LocateVarByName(cpl));
        _SimpleList          *  restrictor       = nil;
        
        if (objectType == HY_BL_LIKELIHOOD_FUNCTION || objectType == HY_BL_SCFG){ 
        // not a BGM
            if (restrictVariable) { // only consider some variables
                _SimpleList variableIDs;
                if (restrictVariable->ObjectClass () == ASSOCIATIVE_LIST)
                    // a list of variables stored as keys in an associative array
                {
                    checkPointer (restrictor = new _SimpleList);
                    _List*  restrictedVariables = ((_AssociativeList *)restrictVariable->GetValue())->GetKeys();
                    for (unsigned long iid = 0; iid < restrictedVariables->lLength; iid++) {
                        _String varID = currentProgram.AddNameSpaceToID(*(_String*)(*restrictedVariables)(iid));
                        variableIDs << LocateVarByName (varID);
                    }
                } else if (restrictVariable->ObjectClass () == STRING)
                    // a single variable stored in a string
                {
                    _String varID = currentProgram.AddNameSpaceToID(*((_FString*)restrictVariable->Compute())->theString);
                    variableIDs << LocateVarByName (varID);
                }
                if (variableIDs.lLength > 0) {
                    checkPointer(restrictor = new _SimpleList ());
                    for (unsigned long var_index = 0; var_index < variableIDs.lLength; var_index++) {
                        long vID = lkf->GetIndependentVars().Find(variableIDs.lData[var_index]);
                        if (vID >= 0) (*restrictor) << vID;
                    }
                    if (restrictor->lLength == 0) {
                        DeleteObject (restrictor);
                        restrictor = nil;
                    }
               }
            }   
            result->SetValue( (_Matrix*)lkf->CovarianceMatrix(restrictor),false);
            DeleteObject (restrictor);
        } else {
        // BGM
            _Matrix * optRes = (_Matrix*)lkf->CovarianceMatrix(nil);
             if (optRes) {
                result->SetValue(optRes,false);
            }
        }
    } else {
        // OPTIMIZE
        if (objectType != HY_BL_NOT_DEFINED) {
            SetStatusLine (_String("Optimizing ") & _HYHBLTypeToText (objectType) & ' ' &lfNameID);
        } else {
            SetStatusLine (_String("Optimizing user function ") &lfNameID);        
        }
        result -> SetValue(lkf->Optimize(),false);
    }

    if (objectType == HY_BL_NOT_DEFINED) {
        DeleteObject (lkf);    // delete the custom function object
    }
    
    
    SetStatusLine ("Finished with the optimization");

    return true;
}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleComputeLFFunction (_ExecutionList& currentProgram) {
        
    currentProgram.currentCommand++;

    _String *arg1       = (_String*)parameters(0),
            *arg2      = (_String*)parameters(1),
            name2Find   = AppendContainerName(*arg1,currentProgram.nameSpacePrefix);
            
    // bool isSCFG  = false;

    long                   objectType = HY_BL_LIKELIHOOD_FUNCTION|HY_BL_SCFG|HY_BL_BGM;
    _LikelihoodFunction    *lf       = (_LikelihoodFunction*)_HYRetrieveBLObjectByName (name2Find, objectType,nil, true, true);

    if (*arg2 == lfStartCompute) {
            lf->PrepareToCompute(true);
    } else if (*arg2 == lfDoneCompute) {
            lf->DoneComputing (true);
    } else {
        if (!lf->HasBeenSetup()) {
            WarnError (_String("Please call LFCompute (lf_id, ")&lfStartCompute&") before evaluating the likelihood function");
            return false;
        } else {
            _Variable* rec = CheckReceptacleCommandID (&AppendContainerName(*arg2,currentProgram.nameSpacePrefix), HY_HBL_COMMAND_LFCOMPUTE,true);
             if (!rec) {
                return false;
            }
            rec->SetValue(new _Constant (lf->Compute()),false);
        }
    }

    return true;
    
}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleSelectTemplateModel (_ExecutionList& currentProgram) {
        
    currentProgram.currentCommand++;

    SetStatusLine ("Waiting for model selection");

    _String     modelFile,
                errMsg;

    ReadModelList();

    if (((_String*)parameters(0))->Equal(&useLastModel)) {
        if (lastModelUsed.sLength) {
            PushFilePath (lastModelUsed);
        } else {
            WarnError (_String("First call to SelectTemplateModel. ") & useLastModel &" is meaningless.");
            return false;
        }
    } else {
        _String filterName (currentProgram.AddNameSpaceToID(*(_String*)parameters(0)));
        long            objectType = HY_BL_DATASET_FILTER;
        _DataSetFilter* thisDF = (_DataSetFilter*)_HYRetrieveBLObjectByName (filterName, objectType,nil,true);
                // decide what this DF is comprised of

        _String             dataType;
        long                dataDimension   = thisDF->GetDimension(),
                            unitLength      = thisDF->GetUnitLength();

        _TranslationTable*  thisTT = thisDF->GetData()->GetTT();
        
        if (unitLength==1) {
            if (thisTT->IsStandardNucleotide()) {
                dataType = "nucleotide";
            } else if (thisTT->IsStandardAA()) {
                dataType = "aminoacid";
            }
        } else {
            if (thisTT->IsStandardNucleotide())
                if (unitLength==3) {
                    dataType = "codon";
                } else if (unitLength==2) {
                    dataType = "dinucleotide";
                }
        }

        if (!dataType.sLength) {
            WarnError (_String("DataSetFilter '")&filterName&"' contains non-standard data and SelectTemplateModel is not applicable.");
            return false;
        }

        _SimpleList matchingModels;

        for (unsigned long model_index = 0; model_index < templateModelList.lLength; model_index++) {
            _List *model_components = (_List*)templateModelList(model_index);
            
            if (dataType.Equal((_String*)model_components->GetItem(3))) {
                _String * dim = (_String*)model_components->GetItem(2);
                if (*dim==_String("*")|| dataDimension == dim->toNum()) {
                    matchingModels << model_index;
                }
            }
        }

        if (!matchingModels.lLength) {
            WarnError ((_String)("DataSetFilter '")&filterName&"' could not be matched with any template models.");
            return false;
        }
        unsigned long model_id = HY_MAX_LONG_VALUE;

        if (currentProgram.stdinRedirect) {
            errMsg = currentProgram.FetchFromStdinRedirect ();
            for (model_id = 0; model_id<matchingModels.lLength; model_id++)
                if (errMsg.Equal((_String*)(*(_List*)templateModelList(matchingModels(model_id)))(0))) {
                    break;
                }

            if (model_id >= matchingModels.lLength) {
                WarnError (errMsg & " is not a valid model (with input redirect) in call to SelectTemplateModel");
                return false;
            }
        } else {
#ifdef __HEADLESS__
            WarnError ("Unhandled standard input interaction in SelectTemplateModel for headless HyPhy");
            return false;
#else
#if defined __UNIX__ && !defined __HYPHYQT__ && !defined __HYPHY_GTK__
            while (model_id == HY_MAX_LONG_VALUE) {
                printf ("\n\n               +--------------------------+\n");
                printf     ("               | Select a standard model. |\n");
                printf     ("               +--------------------------+\n\n\n");
                
                for (model_id = 0; model_id<matchingModels.lLength; model_id++) {
                    printf ("\n\t(%s):%s",((_String*)(*(_List*)templateModelList(matchingModels(model_id)))(0))->getStr(),
                            ((_String*)(*(_List*)templateModelList(matchingModels(model_id)))(1))->getStr());
                }
                printf ("\n\n Please type in the abbreviation for the model you want to use:");
                dataType.CopyDynamicString(StringFromConsole());
                dataType.UpCase();
                for (model_id = 0; model_id<matchingModels.lLength; model_id++) {
                    if (dataType.Equal((_String*)(*(_List*)templateModelList(matchingModels(model_id)))(0))) {
                        break;
                    }
                }
                if (model_id==matchingModels.lLength) {
                    model_id=HY_MAX_LONG_VALUE;
                }
            }
#endif
#if !defined __UNIX__ ||  defined __HYPHYQT__ ||  defined __HYPHY_GTK__
            _SimpleList choiceDummy (2,0,1), selDummy;
            model_id = HandleListSelection (templateModelList, choiceDummy, matchingModels, "Choose one of the standard substitution models",selDummy,1,nil);
            if (model_id==-1) {
                terminateExecution = true;
                return false;
            }
#endif
#endif
        }
        modelFile = _HYStandardDirectory (HY_HBL_DIRECTORY_TEMPLATE_MODELS) &*((_String*)(*(_List*)templateModelList(matchingModels(model_id)))(4));
        PushFilePath (modelFile, false);
    }

    _ExecutionList      stdModel;
    if (currentProgram.nameSpacePrefix) {
        stdModel.SetNameSpace (*currentProgram.nameSpacePrefix->GetName());
    }

    ReadBatchFile       (modelFile,stdModel);
    PopFilePath         ();
    lastModelUsed       = modelFile;

    stdModel.stdinRedirectAux = currentProgram.stdinRedirectAux;
    stdModel.stdinRedirect    = currentProgram.stdinRedirect;
    stdModel.Execute();
    stdModel.stdinRedirectAux = nil;
    stdModel.stdinRedirect    = nil;   

    return true;
    
}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleUseModel (_ExecutionList& currentProgram) {
    
    currentProgram.currentCommand++;
    _String namedspacedMM (currentProgram.AddNameSpaceToID(*(_String*)parameters(0)));
    long mID = FindModelName(namedspacedMM);

    if (mID<0 && !useNoModel.Equal((_String*)parameters(0))) {
        WarnError(*(_String*)parameters(0) & " does not refer to a valid defined substitution model in call to " & _HY_ValidHBLExpressions.RetrieveKeyByPayload(HY_HBL_COMMAND_USE_MODEL));
        return false;
    } else {
        lastMatrixDeclared = mID;
    }

    return true;
}

//____________________________________________________________________________________
bool      _ElementaryCommand::HandleSetParameter (_ExecutionList& currentProgram) {

    currentProgram.currentCommand++;
    /*
        first check to see if matrix parameters here are valid
    */
    
    _String *currentArgument = (_String*)parameters(0),
             nmspc           = AppendContainerName(*currentArgument,currentProgram.nameSpacePrefix),
             errMsg,
             result;

    if (currentArgument->Equal (&randomSeed)) {
        globalRandSeed = ProcessNumericArgument ((_String*)parameters(1), currentProgram.nameSpacePrefix);
        init_genrand (globalRandSeed);
        setParameter (randomSeed, ((long)globalRandSeed));
        return true;
    }

    if (currentArgument->Equal (&randomSeed)) {
        globalRandSeed = ProcessNumericArgument ((_String*)parameters(1), currentProgram.nameSpacePrefix);
        init_genrand (globalRandSeed);
        setParameter (randomSeed, ((long)globalRandSeed));
        return true;
    }
    
    if (currentArgument->Equal (&deferConstrainAssignment)) {
        bool on = ProcessNumericArgument ((_String*)parameters(1), currentProgram.nameSpacePrefix);
        if (on) {
            deferSetFormula = (_SimpleList*)checkPointer(new _SimpleList);
        } else if (deferSetFormula) {
            FinishDeferredSF ();
        }
        return true;
    }

    if (currentArgument->Equal (&_hyExecutionErrorMode)) {
        currentProgram.errorHandlingMode = ProcessNumericArgument ((_String*)parameters(1), currentProgram.nameSpacePrefix);
        return true;
    }

    if (currentArgument->Equal (&statusBarProgressValue)) {
#if !defined __UNIX__
        SetStatusLine     (empty,empty, empty,  ProcessNumericArgument ((_String*)parameters(1), currentProgram.nameSpacePrefix), HY_SL_PERCENT);
#endif
        return true;
    }

    if (currentArgument->Equal (&statusBarUpdateString)) {
        _String sbar_value = ProcessLiteralArgument ((_String*)parameters(1), currentProgram.nameSpacePrefix);

#if defined __UNIX__ 
        #if not defined __HYPHY_GTK__ && not defined __HEADLESS__
            SetStatusLineUser     (sbar_value);
        #else
            SetStatusLine (sbar_value);
        #endif 
#else
        SetStatusLine     (empty,sbar_value, empty, 0, HY_SL_TASK);
#endif
        return true;
    }


    long objectIndex,
         typeFlag    = HY_BL_ANY;
    
    BaseRef theObject      = _HYRetrieveBLObjectByName (nmspc, typeFlag, &objectIndex);
    
    switch (typeFlag)
    {
        case HY_BL_BGM: { // BGM Branch
            currentArgument = (_String*)parameters(1);

            _BayesianGraphicalModel * lkf = (_BayesianGraphicalModel *) theObject;
            // set data matrix
            if (currentArgument->Equal (&bgmData)) {
                _Matrix     * dataMx = (_Matrix *) FetchObjectFromVariableByType ( &AppendContainerName ( *(_String *) parameters (2), currentProgram.nameSpacePrefix), MATRIX, HY_HBL_COMMAND_SET_PARAMETER);
                if (dataMx) {
                    long    num_nodes = ((_BayesianGraphicalModel *)lkf)->GetNumNodes();

                    if (dataMx->GetVDim() == num_nodes) {
                        ((_BayesianGraphicalModel *)lkf)->SetDataMatrix ((_Matrix *) dataMx);
                    } else {
                        currentProgram.ReportAnExecutionError (_String("Data matrix columns (") & dataMx->GetVDim() & " ) does not match number of nodes in graph (" & num_nodes & ")");
                        return false;
                    }
                } else {
                     return false;
                }

            }

            // restore node score cache
            else if (currentArgument->Equal (&bgmScores)) {
                _AssociativeList * cacheAVL = (_AssociativeList *)FetchObjectFromVariableByType ( &AppendContainerName ( *(_String *) parameters (2), currentProgram.nameSpacePrefix), ASSOCIATIVE_LIST, HY_HBL_COMMAND_SET_PARAMETER);
                if (cacheAVL) {
                    ((_BayesianGraphicalModel *)lkf)->ImportCache (cacheAVL);
                } else {
                    return false;
                }
            }

            // set structure to user-specified adjacency matrix
            else if (currentArgument->Equal (&bgmGraph)) {
                _Matrix     * graphMx   = (_Matrix *) FetchObjectFromVariableByType ( &AppendContainerName ( *(_String *) parameters (2), currentProgram.nameSpacePrefix), MATRIX, HY_HBL_COMMAND_SET_PARAMETER);

                if (graphMx) {
                    long    num_nodes = ((_BayesianGraphicalModel *)lkf)->GetNumNodes();

                    if (graphMx->GetHDim() == num_nodes && graphMx->GetVDim() == num_nodes) {
                        ((_BayesianGraphicalModel *)lkf)->SetStructure ((_Matrix *) graphMx->makeDynamic());
                    } else {
                        currentProgram.ReportAnExecutionError("Dimension of graph does not match current graph");
                        return false;
                    }
                } else {
                    return false;
                }
            }

            // set constraint matrix
            else if (currentArgument->Equal (&bgmConstraintMx)) {
                _Matrix     * constraintMx  = (_Matrix *) FetchObjectFromVariableByType ( &AppendContainerName ( *(_String *) parameters (2), currentProgram.nameSpacePrefix), MATRIX, HY_HBL_COMMAND_SET_PARAMETER);

                if (constraintMx) {
                    long    num_nodes = ((_BayesianGraphicalModel *)lkf)->GetNumNodes();

                    if (constraintMx->GetHDim() == num_nodes && constraintMx->GetVDim() == num_nodes) {
                        ((_BayesianGraphicalModel *)lkf)->SetConstraints ((_Matrix *) constraintMx->makeDynamic());
                    } else {
                        currentProgram.ReportAnExecutionError ("Dimensions of constraint matrix do not match current graph");
                        return false;
                    }
                } else {
                    return false;
                }
            }

            // set node order
            else if (currentArgument->Equal (&bgmNodeOrder)) {
                _Matrix     * orderMx   = (_Matrix *) FetchObjectFromVariableByType ( &AppendContainerName ( *(_String *) parameters (2), currentProgram.nameSpacePrefix), MATRIX, HY_HBL_COMMAND_SET_PARAMETER);

                if (orderMx) {
                    // UNDER DEVELOPMENT  April 17, 2008 afyp
                    long    num_nodes = ((_BayesianGraphicalModel *)lkf)->GetNumNodes();

                    _SimpleList     * orderList = new _SimpleList();

                    orderList->Populate (num_nodes, 0, 0);

                    if (orderMx->GetVDim() == num_nodes) {
                        for (long i = 0; i < num_nodes; i++) {
                            orderList->lData[i] = (long) ((*orderMx) (0, i));
                        }

                        ((_BayesianGraphicalModel *)lkf)->SetNodeOrder ( (_SimpleList *) orderList->makeDynamic() );
                    } else {
                        currentProgram.ReportAnExecutionError ("Length of order vector doesn't match number of nodes in graph");
                        return false;
                    }
                } else {
                    return false;
                }
            }


            // set network parameters
            else if (currentArgument->Equal (&bgmParameters)) {
                    _AssociativeList* inAVL = (_AssociativeList*)FetchObjectFromVariableByType ( &AppendContainerName ( *(_String *) parameters (2), currentProgram.nameSpacePrefix), ASSOCIATIVE_LIST, HY_HBL_COMMAND_SET_PARAMETER);
                if (inAVL) {
                    ((_BayesianGraphicalModel *)lkf)->ImportCache (inAVL);
                } else {
                    return false;
                }
            }


            // anything else
            else {
                currentProgram.ReportAnExecutionError (*currentArgument & " is not a valid BGM parameter");
                return false;
            }
        } // end BGM
        break;
        
        case HY_BL_SCFG:
        case HY_BL_LIKELIHOOD_FUNCTION:
        {
            currentArgument = (_String*)parameters(1);
            if (typeFlag == HY_BL_SCFG && currentArgument->Equal (&scfgCorpus)) {
                ((Scfg*)theObject)->SetStringCorpus ((_String*)parameters(2));
            } else {
                _LikelihoodFunction * lkf = (_LikelihoodFunction *) theObject;
                currentArgument = (_String*)parameters(1);
                long g = ProcessNumericArgument(currentArgument,currentProgram.nameSpacePrefix);
                if (g < 0 || g >= lkf->GetIndependentVars().lLength) {
                    currentProgram.ReportAnExecutionError (*currentArgument & " (=" & g & ") is not a valid parameter index");
                    return false;
                }
                currentArgument = (_String*)parameters(2);
                lkf->SetIthIndependent (g,ProcessNumericArgument(currentArgument,currentProgram.nameSpacePrefix));
            }
        }
        break;
        // end SCFG and LF
        
        case HY_BL_DATASET:
        case HY_BL_DATASET_FILTER: {
            _DataSet * ds = nil;
            
            long f  = ProcessNumericArgument ((_String*)parameters(1),currentProgram.nameSpacePrefix);
            if (typeFlag == HY_BL_DATASET) {
                ds = (_DataSet*) theObject;
            }
            else {
                _DataSetFilter *dsf = (_DataSetFilter*)theObject;
                ds = dsf->GetData ();
                if (f >= 0 && f < dsf->theNodeMap.lLength){
                    f  = dsf->theNodeMap.lData[f];
                }
                else
                    f = -1;
            }
            
        
            _List*  dsNames = &ds->GetNames();
            
            if (f<0 || f>=dsNames->lLength) {
                currentProgram.ReportAnExecutionError (*((_String*)parameters(1)) & " (=" & f & ") is not a valid sequence index");
                return false;
            }

            dsNames->Replace(f, new _String(ProcessLiteralArgument ((_String*)parameters(2),currentProgram.nameSpacePrefix)), false);
        } // end data set and data set filter
        break; 
        // Dataset and Datasetfilter
        
        default:
            // check to see if this is a calcnode
            _CalcNode* treeNode = (_CalcNode*)FetchObjectFromVariableByType(&nmspc, TREE_NODE);
            if (treeNode) { 
                if (*((_String*)parameters(1)) == _String("MODEL")) {
                    _String modelName = AppendContainerName(*((_String*)parameters(2)),currentProgram.nameSpacePrefix);
                    long modelType = HY_BL_MODEL, modelIndex;
                    BaseRef modelObject      = _HYRetrieveBLObjectByName (modelName, modelType, &modelIndex, true);
                    if (modelObject) {
                        _VariableContainer * parentTree = treeNode->ParentTree();
                        if (!parentTree) {
                            currentProgram.ReportAnExecutionError (*((_String*)parameters(0)) & " is an orphaned tree node (the parent tree has been deleted)");
                            return false;
                        }
                        long pID, lfID = ((_TheTree*)parentTree->Compute())->IsLinkedToALF(pID);
                        if (lfID>=0){
                             currentProgram.ReportAnExecutionError ((*parentTree->GetName()) & " is linked to a likelihood function (" & *_HBLObjectNameByType (HY_BL_LIKELIHOOD_FUNCTION, lfID) &") and cannot be modified ");
                             return false;
                        }
                        
                        treeNode->ReplaceModel (modelName, parentTree);
                        break;
                    }
                    else {
                        currentProgram.ReportAnExecutionError (*((_String*)parameters(2)) & " does not appear to be a valid model name");
                        return false;
                    }
                } else {
                     currentProgram.ReportAnExecutionError (*((_String*)parameters(1)) & " is not a supported parameter type for a tree node argument");
                     return false;
                }
            }
            
            currentProgram.ReportAnExecutionError (*currentArgument & " is not a valid likelihood function/data set filter/tree topology/tree node");
            return false;

    
    } // end cases
    return true;
}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleAssert (_ExecutionList& currentProgram) {
    currentProgram.currentCommand++;

    _String assertion (*(_String*)parameters(0));

    _Formula rhs, lhs;
    _FormulaParsingContext fpc (nil, currentProgram.nameSpacePrefix);
    if (Parse (&rhs,assertion,fpc,&lhs) == HY_FORMULA_EXPRESSION) {
        _PMathObj assertionResult = rhs.Compute();
        if (assertionResult && assertionResult->ObjectClass () == NUMBER) {
            if (CheckEqual(assertionResult->Value(),0.0)) {
                _Parameter whatToDo;
                checkParameter (assertionBehavior, whatToDo, 0.0);

                _String errMsg;

                if (parameters.lLength == 1) {
                    errMsg = _String("Assertion '") & *(_String*)parameters(0) & "' failed.";
                } else {
                    errMsg = ProcessLiteralArgument((_String*)parameters(1),currentProgram.nameSpacePrefix);
                }

                if (CheckEqual (whatToDo, 1.)) {
                    StringToConsole (errMsg);
                    NLToConsole();
                    currentProgram.GoToLastInstruction ();
                } else {
                    currentProgram.ReportAnExecutionError (errMsg);
                    return false;
                }
            }
            return true;
        }
    }
    currentProgram.ReportAnExecutionError(_String("Assertion statement '") & *(_String*)parameters(0) & "' could not be computed or was not numeric.");
 
    return false;
}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleRequireVersion(_ExecutionList& currentProgram){
    currentProgram.currentCommand++;
    _String theVersion = ProcessLiteralArgument ((_String*)parameters (0),currentProgram.nameSpacePrefix);

    if (__KERNEL__VERSION__.toNum() < theVersion.toNum()) {
        currentProgram.ReportAnExecutionError (_String ("Current batch file requires at least version :")& theVersion &" of HyPhy. Please download an updated version from http://www.hyphy.org and try again.");
        return false;
    }
    return true;
}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleDeleteObject(_ExecutionList& currentProgram){
    currentProgram.currentCommand++;
    for (unsigned long objCount = 0; objCount < parameters.lLength; objCount++) {
        long       objectType = HY_BL_LIKELIHOOD_FUNCTION,
                   f = -1;
        BaseRef    sourceObject = _HYRetrieveBLObjectByName (AppendContainerName(*(_String*)parameters(objCount),currentProgram.nameSpacePrefix), objectType,&f,false);

        if  (sourceObject) {
            KillLFRecord (f,true);
        }
    }
    return true;
}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleClearConstraints(_ExecutionList& currentProgram){
    currentProgram.currentCommand++;
    for (unsigned long i = 0; i<parameters.lLength; i++) {
        _String cName (currentProgram.AddNameSpaceToID(*(_String*)parameters(i)));
        long cID = LocateVarByName (cName);
        if (cID>=0) { // variable exists
            FetchVar(cID)->ClearConstraints();
        }
    }    
    return true;
}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleMolecularClock(_ExecutionList& currentProgram){
    currentProgram.currentCommand++;
    
    _String    theBaseNode          (currentProgram.AddNameSpaceToID(*(_String*)parameters(0))),
               treeName;
    
    _Variable* theObject = FetchVar (LocateVarByName(theBaseNode));
    
    if (!theObject || (theObject->ObjectClass()!=TREE && theObject->ObjectClass()!=TREE_NODE)) {
        WarnError (_String("Not a defined tree/tree node object '") & theBaseNode & "' in call to " & _HY_ValidHBLExpressions.RetrieveKeyByPayload(HY_HBL_COMMAND_MOLECULAR_CLOCK));
        return false;
    }
    
    _TheTree *theTree = nil;
    if (theObject->ObjectClass() == TREE_NODE) {
        theTree     = (_TheTree*)((_VariableContainer*)theObject)->GetTheParent();
        if (!theTree) {
            WarnError (_String("Internal error - orphaned tree node '") & theBaseNode & "' in call to "& _HY_ValidHBLExpressions.RetrieveKeyByPayload(HY_HBL_COMMAND_MOLECULAR_CLOCK));
            return false;
            
        }
        treeName    = *theTree->GetName();
        theBaseNode = theObject->GetName()->Cut(treeName.sLength+1,-1);
    } else {
        treeName    = *theObject->GetName();
        theTree     = (_TheTree*)theObject;
        theBaseNode = empty;
    }
    
    theTree->MolecularClock(theBaseNode,parameters);
    return true;
}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleGetURL(_ExecutionList& currentProgram){
    currentProgram.currentCommand++;
    
    _String url   (ProcessLiteralArgument((_String*)parameters(1),currentProgram.nameSpacePrefix)),
            *arg1 = (_String*)parameters(0),
            *act  = parameters.lLength>2?(_String*)parameters(2):nil,
            errMsg;

    if (act==nil) {
        _Variable * rec = CheckReceptacleCommandID (&AppendContainerName(*arg1,currentProgram.nameSpacePrefix),HY_HBL_COMMAND_GET_URL, true, false, &currentProgram);

        if (!rec) {
            return false;
        }

        if (Get_a_URL(url)) {
            rec->SetValue(new _FString (url,false),false);
        } else {
            errMsg = _String ("Could not fetch '") & url & "'";
        }
    } else {
        if (act->Equal(&getURLFileFlag)) {
            _String fileName (ProcessLiteralArgument(arg1,currentProgram.nameSpacePrefix));
            fileName.ProcessFileName (true,false,(Ptr)currentProgram.nameSpacePrefix);
            if (!Get_a_URL(url, &fileName)) {
                errMsg = _String ("Could not fetch '") & url & "'";
            }
        } else {
            errMsg = "Unknown action flag";
        }
    }
    if (errMsg.sLength) {
        currentProgram.ReportAnExecutionError (errMsg); 
        return false;
    }

    return true;
}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleGetString (_ExecutionList& currentProgram){
    currentProgram.currentCommand++;
    
    _String  errMsg,
             *result = nil;

    long    f,
            sID,
            sID2 = -1;

    _Variable * theReceptacle = CheckReceptacleCommandID (&AppendContainerName(*((_String*)parameters(0)),currentProgram.nameSpacePrefix),HY_HBL_COMMAND_GET_STRING, true, false, &currentProgram);

    if (!theReceptacle) {
        return false;
    }

    sID = ProcessNumericArgument ((_String*)parameters(2),currentProgram.nameSpacePrefix);
    if (parameters.lLength>3) {
        sID2 = ProcessNumericArgument ((_String*)parameters(3),currentProgram.nameSpacePrefix);
    }

    f = _HY_GetStringGlobalTypes.Find((_String*)parameters(1));
    if (f >=0 ) {
        f = _HY_GetStringGlobalTypes.GetXtra (f);
    }

    switch (f) {

        case HY_BL_LIKELIHOOD_FUNCTION: // LikelihoodFunction
        case HY_BL_DATASET:
        case HY_BL_DATASET_FILTER:
        case HY_BL_SCFG:
        case HY_BL_BGM: {
            result = (_String*)_HBLObjectNameByType(f,sID);
            if (result) {
                result = (_String*) result->makeDynamic();
				//ReportWarning(_String((const char*)"In HandleGetString(): ") & result);
            }
            break;
        }
           
        case HY_BL_HBL_FUNCTION: // UserFunction
            result = (_String*)_HBLObjectNameByType(HY_BL_HBL_FUNCTION,sID);
            if (result) {
                _AssociativeList * resAVL = (_AssociativeList *)checkPointer(new _AssociativeList);
                resAVL->MStore ("ID", new _FString (*result), false);
                resAVL->MStore ("Arguments", new _Matrix(*(_List*)batchLanguageFunctionParameterLists(sID)), false);
                theReceptacle->SetValue (resAVL,false);
                return true;
            } 
            break;
            
        case HY_BL_TREE: { // Tree
            // 20110608 SLKP: REFACTOR into a separate function
            // I am sure this is used somewhere else (perhaps for other types)
            result = FetchMathObjectNameOfTypeByIndex (TREE, sID);
            if (result) {
                result = (_String*) result->makeDynamic();
            }
            break;
        }

        default: { // everything else...
            // decide what kind of object current argument represents
            
            _String *currentArgument = (_String*)parameters(1),
                    nmspaced       = AppendContainerName(*currentArgument,currentProgram.nameSpacePrefix);
            long    typeFlag       = HY_BL_ANY,
                    index          = -1;
                    
            BaseRef theObject      = _HYRetrieveBLObjectByName (nmspaced, typeFlag, &index);

            if (theObject) {
                switch (typeFlag) {
                case HY_BL_DATASET: {
                    _DataSet* dataSetObject = (_DataSet*)theObject;
                    if (sID>=0 && sID<dataSetObject->NoOfSpecies()) {
                        result = (_String*)(dataSetObject->GetNames())(sID)->makeDynamic();
                    } else {
                        if (sID < 0) {
                            theReceptacle->SetValue (new _Matrix (dataSetObject->GetNames()), false);
                            return true;
                        }
                    }
                    break;
                }
                case HY_BL_DATASET_FILTER: {
                    _DataSetFilter* dataSetFilterObject = (_DataSetFilter*)theObject;

                    if (sID >=0 && sID<dataSetFilterObject->NumberSpecies()) {
                        result = (_String*)(dataSetFilterObject->GetData()->GetNames())(dataSetFilterObject->theNodeMap(sID))->makeDynamic();
                    } else {
                        if (sID < 0) {
                            _List filterSeqNames,
                                  *originalNames = &dataSetFilterObject->GetData()->GetNames();
                            for (long seqID=0; seqID<dataSetFilterObject->NumberSpecies(); seqID++) {
                                filterSeqNames << (*originalNames)(dataSetFilterObject->theNodeMap(seqID));
                            }
                            theReceptacle->SetValue (new _Matrix (filterSeqNames), false);
                            return true;
                        }
                    }
                    break;
                }
                case HY_BL_BGM: {
                    //ReportWarning(_String("In HandleGetString() for case HY_BL_BGM"));
					_BayesianGraphicalModel * this_bgm      = (_BayesianGraphicalModel *) theObject;

                    switch (sID) {
                        case HY_HBL_GET_STRING_BGM_SCORE: {   // return associative list containing node score cache
                            _AssociativeList        * export_alist  = new _AssociativeList;

                            if (this_bgm -> ExportCache (export_alist)) {
                                theReceptacle -> SetValue (export_alist, false);
                                return true;
                            } else {
                                DeleteObject (export_alist);
                                errMsg = _String ("Failed to export node score cache for BGM '") & nmspaced & "'";
                            }

                            break;
                        }
                        case HY_HBL_GET_STRING_BGM_SERIALIZE: {   // return associative list with network structure and parameters
                            result = new _String (1024L, true);
                            this_bgm -> SerializeBGM (*result);
                            result->Finalize();
                            theReceptacle->SetValue (new _FString (result),false);
                            return true;
                        }
                        default: {
                            errMsg = _String ("Unrecognized index ") & sID & " for a BGM object";
                        }
                    }
                }
                case HY_BL_LIKELIHOOD_FUNCTION:
                case HY_BL_SCFG: {
                
                    _LikelihoodFunction *lf = (_LikelihoodFunction*)theObject;
                    if (sID>=0) {
                        if (sID<lf->GetIndependentVars().lLength) {
                            result = (_String*)(LocateVar(lf->GetIndependentVars().lData[sID])->GetName())->makeDynamic();
                        } else {
                            if (sID<lf->GetIndependentVars().lLength+lf->GetDependentVars().lLength) {
                                result = (_String*)(LocateVar(lf->GetDependentVars().lData[sID-lf->GetIndependentVars().lLength])->GetName())->makeDynamic();
                            }
                        }
                    } else {
                        _AssociativeList* resList = lf->CollectLFAttributes ();
                        if (typeFlag == HY_BL_SCFG) {
                            ((Scfg*)lf)->AddSCFGInfo (resList);
                        }
                        theReceptacle->SetValue (resList,false);
                        return true;
                    }
                    break;
                }
                
                case HY_BL_MODEL: {
                    if (sID>=0) {
                        // check to make see if the 
                       if (sID2 < 0) { // get the sID's parameter name
                            _SimpleList     modelP;
                            _AVLList        modelPA (&modelP);
                            ScanModelForVariables (index, modelPA, false, -1, false);
                            modelPA.ReorderList();
                            if (sID<modelP.lLength) {
                                result = (_String*)LocateVar(modelP.lData[sID])->GetName()->makeDynamic();
                            } 

                        } else { // get the formula for cell (sID, sID2)
                            if (!IsModelOfExplicitForm (index)) {
                                _Variable*      theMx = (_Variable*)theObject;
                                _Formula * cellFla = ((_Matrix*)theMx->GetValue())->GetFormula (sID,sID2);
                                if (cellFla) {
                                    result = new _String((_String*)cellFla->toStr());
                                }
                            }
                        }
                        
                    } else {
                        _Variable   * tV, * tV2;
                        bool         mByF;
                        RetrieveModelComponents (index, tV, tV2, mByF);

                        if (tV) {
                            if (sID == -1) { // branch length expression
                                result = ((_Matrix*)tV->GetValue())->BranchLengthExpression((_Matrix*)tV2->GetValue(),mByF);
                            } else
                                /*
                                    returns an AVL with keys
                                    "RATE_MATRIX" - the ID of the rate matrix
                                    "EQ_FREQS"    - the ID of eq. freq. vector
                                    "MULT_BY_FREQ" - a 0/1 flag to determine which format the matrix is in.
                                */
                            {
                                _AssociativeList * resList = new _AssociativeList;
                                resList->MStore ("RATE_MATRIX",new _FString(*tV->GetName()),false);
                                resList->MStore ("EQ_FREQS",new _FString(*tV2->GetName()),false);
                                resList->MStore ("MULT_BY_FREQ",new _Constant (mByF),false);
                                theReceptacle->SetValue (resList,false);
                                return true;
                            }
                        }
                    }
                    break;
                }
                case HY_BL_HBL_FUNCTION: {
                    _AssociativeList * resAVL = (_AssociativeList *)checkPointer(new _AssociativeList);
                    resAVL->MStore ("ID", new _FString (*_HBLObjectNameByType (HY_BL_HBL_FUNCTION, index, false)), false);
                    resAVL->MStore ("Arguments", new _Matrix(*(_List*)batchLanguageFunctionParameterLists(index)), false);
                    resAVL->MStore("Body", new _FString (((_ExecutionList*)batchLanguageFunctions(index))->sourceText,false),false);
                    theReceptacle->SetValue (resAVL,false);
                    return true;
                }
            } // end of "switch" 
        }
        else {
            if (currentArgument->Equal(&versionString)) {
                if (sID > 1.5)
#ifdef __HEADLESS__
                    result = new _String(_String ("Library version ") & __KERNEL__VERSION__);
#else
#ifdef __MAC__
                    result = new _String(_String("Macintosh ") & __KERNEL__VERSION__);
#else
#ifdef __WINDOZE__
                    result = new _String(_String("Windows ") & __KERNEL__VERSION__);
#else
                    result = new _String(_String("Source ") & __KERNEL__VERSION__);
#endif
#endif
#endif
                    else if (sID > 0.5) {
                        result = new _String(GetVersionString());
                    } else {
                        result = new _String(__KERNEL__VERSION__);
                    }
                } else if (currentArgument->Equal(&timeStamp)) {
                    result = new _String(GetTimeStamp(sID < 0.5));
                } else {
                  _Variable* theVar = FetchVar(LocateVarByName (*currentArgument));
                  if (theVar) {
                    if (theVar->IsIndependent()) {
                      result = (_String*)theVar->toStr();
                    } else {
                      if (sID == -1){
                        
                        _SimpleList vL;
                        _AVLList    vAVL (&vL);
                        theVar->ScanForVariables (vAVL, true);
                        vAVL.ReorderList();
                        _AssociativeList   * resL = (_AssociativeList *) checkPointer (new _AssociativeList);
                        _List splitVars;
                        SplitVariableIDsIntoLocalAndGlobal (vL, splitVars);
                        InsertVarIDsInList (resL, "Global", *(_SimpleList*)splitVars(0));
                        InsertVarIDsInList (resL, "Local",   *(_SimpleList*)splitVars(1));
                        
                        theReceptacle->SetValue (resL,false);
                        return true;
                      }
                      
                      else {  // formula string
                        _Matrix * formula_matrix = (sID2 >= 0 && theVar->ObjectClass() == MATRIX) ? (_Matrix*)theVar->GetValue () : nil;
                        if (formula_matrix) {
                         _Formula* cell = formula_matrix->GetFormula(sID, sID2);
                         if (cell) {
                            result = (_String*) cell->toStr();
                          }
                        } else {
                          result = (_String*)theVar->GetFormulaString ();
                        }
                      }
                    }
                  } else {
                    errMsg = _String ("'") & *currentArgument & "' is not an allowed argument type ";
                  }
                }
            }
        }
    }

    if (errMsg.sLength) {
        currentProgram.ReportAnExecutionError (errMsg); 
        DeleteObject (result);
        result = nil;    
    }
    
    if (result) {
        theReceptacle->SetValue (new _FString (result),false);
        return true;
    }

    theReceptacle->SetValue (new _MathObject(), false);
    

    return false;
}


//____________________________________________________________________________________

bool      _ElementaryCommand::HandleExport(_ExecutionList& currentProgram){

    currentProgram.currentCommand++;
    
    _String objectID (currentProgram.AddNameSpaceToID(*(_String*)parameters(1))),
            arg1 (currentProgram.AddNameSpaceToID(*(_String*)parameters(0))),
            errMsg;
    
    _Variable * theReceptacle = CheckReceptacleCommandID (&AppendContainerName(arg1,currentProgram.nameSpacePrefix),HY_HBL_COMMAND_EXPORT, true, false, &currentProgram);
    if (!theReceptacle) {
        return false;
    }    

    _FString        * outLF = new _FString (new _String (8192L,1));
    checkPointer    (outLF);
    long typeFlag = HY_BL_MODEL | HY_BL_LIKELIHOOD_FUNCTION | HY_BL_DATASET_FILTER,
             index;
        
    BaseRef objectToExport = _HYRetrieveBLObjectByName (objectID, typeFlag, &index);
    if (! objectToExport) {
        errMsg = _String ("'") & objectID & "' is not a supported type";
    } else {
        switch (typeFlag) {
            case HY_BL_LIKELIHOOD_FUNCTION: {
                ((_LikelihoodFunction*)objectToExport)->SerializeLF (*outLF->theString);
                outLF->theString->Finalize();
                break;
            }
            case HY_BL_DATASET_FILTER: {
                outLF->theString->Finalize();
                DeleteObject (outLF->theString);
                checkPointer (outLF->theString = new _String ((_String*)((_DataSetFilter*)objectToExport)->toStr()));
                break;
            }
            case HY_BL_MODEL: {
                SerializeModel (*outLF->theString,index,nil,true);
                outLF->theString->Finalize();
                break;
            }
                
        }
    }

    if (errMsg.sLength) {
        outLF->theString->Finalize();
        DeleteObject (outLF);
        currentProgram.ReportAnExecutionError (errMsg); 
        theReceptacle->SetValue (new _MathObject, false);
        return false;
    }
    
    theReceptacle->SetValue (outLF,false);

    return true;
}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleDifferentiate(_ExecutionList& currentProgram){

    currentProgram.currentCommand++;
   
    _String  arg1 (currentProgram.AddNameSpaceToID(*(_String*)parameters(0))),
             errMsg,
             expressionToParse = *(_String*)parameters(1);
             
    _Formula  *theResult = nil;

    _Variable * theReceptacle = CheckReceptacleCommandID (&AppendContainerName(arg1,currentProgram.nameSpacePrefix),HY_HBL_COMMAND_DIFFERENTIATE, true, false, &currentProgram);
    if (!theReceptacle) {
        return false;
    }    


    _Formula theExpression (expressionToParse,currentProgram.nameSpacePrefix, &errMsg);
    
    if (!theExpression.IsEmpty() && errMsg.sLength == 0) {
        long times = 1;
        if (parameters.lLength==4) {
            times = ProcessNumericArgument ((_String*)parameters(3),currentProgram.nameSpacePrefix, &currentProgram);
            if (!numericalParameterSuccessFlag) {
                return false;
            }
        }
        if (times <= 0) {
            errMsg = "The number of times to differentiate must be a non-negative integer";
        }

        theResult = theExpression.Differentiate (*(_String*)parameters(2), false);
        for (; times>1 && theResult; times--) {
            _Formula * temp = theResult->Differentiate (*(_String*)parameters(2));
            delete (theResult);
            theResult = temp;
        }
    }

    if (errMsg.sLength || theResult == nil) {
        if (theResult) { 
            delete (theResult); 
        } else {
            errMsg = _String("Differentiation of '") & *(_String*)parameters(1) & "' failed";
        }
        currentProgram.ReportAnExecutionError (errMsg); 
        theReceptacle->SetValue (new _MathObject, false);
        return false;
    }
    
    theReceptacle->SetFormula (*theResult);
    if (theResult) delete (theResult);

    return true;
}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleFprintf (_ExecutionList& currentProgram)
{
    currentProgram.currentCommand++;
    _String* targetName = (_String*)parameters(0),
             fnm;
    
    bool     doClose                 = true,
             print_to_stdout         = false;
    
    FILE*   dest = nil;
    
    try {
        bool    skipFilePathEval        = false;
        
        if (targetName->Equal(&stdoutDestination)) {
            _FString * redirect = (_FString*)FetchObjectFromVariableByType (&blFprintfRedirect, STRING);
            if (redirect && redirect->theString->sLength) {
                if (redirect->theString->Equal (&blFprintfDevNull)) {
                    return true; // "print" to /dev/null
                } else {
                    skipFilePathEval = true;
                    targetName       = redirect->theString;
                } 
            }
            else {
                print_to_stdout = true;
            }
        }
        
        checkParameter (printDigitsSpec,printDigits,0);
        
        if (!print_to_stdout) {
            fnm = *targetName;
            if (fnm.Equal(&messageLogDestination)) {
                if ((dest = globalMessageFile) == nil) {
                    return true; // requested print to MESSAGE_LOG, but it does not exist
                                 // (e.g. remote MPI nodes, or running from a read only location
                }
            } else {
                if (skipFilePathEval == false && !fnm.IsALiteralArgument()) {
                    fnm = GetStringFromFormula (&fnm,currentProgram.nameSpacePrefix);
                }
                
                if (!fnm.ProcessFileName(true,false,(Ptr)currentProgram.nameSpacePrefix, false, &currentProgram)) {
                    return false;
                }
                
                
                long k  = openFileHandles.Find (&fnm);
                
                doClose = k<0;
                
                if (!doClose) {
                    dest = (FILE*)openFileHandles.GetXtra (k);
                } else {
                    if ((dest = doFileOpen (fnm.getStr(), "a")) == nil)
                        throw  (_String  ("Could not create/open output file at path '") & fnm & "'.");
                }
            }
        }
        
        for (long i = 1; i<parameters.lLength; i++) {
            _String    *varname = (_String*)parameters(i);
            
            BaseRef    thePrintObject   =   nil;
            _Formula   f;
            
            if (varname->Equal(&clearFile)) {
                if (!print_to_stdout && dest) {
                    fclose (dest);
                    dest = doFileOpen (fnm.getStr(), "w");
                    long k = openFileHandles.Find (&fnm);
                    if (k>=0) {
                        openFileHandles.SetXtra(k, (long)dest);
                    }
                }
            } else if (varname->Equal(&keepFileOpen) && !print_to_stdout) {
                if (openFileHandles.Find (&fnm) < 0) {
                    openFileHandles.Insert (fnm.makeDynamic(), (long)dest);
                }
                doClose = false;
            } else if (varname->Equal(&closeFile) && !print_to_stdout) {
                openFileHandles.Delete (&fnm, true);
                doClose = true;
            } else if (varname->Equal(&systemVariableDump)) {
                thePrintObject=&variableNames;
            } else if (varname->Equal(&selfDump)) {
                thePrintObject=&currentProgram;
            } else {
                // check for possible string reference
                
                _String    temp    = ProcessStringArgument (varname),
                           nmspace;
                
                if (temp.sLength > 0) {
                    nmspace = AppendContainerName(temp,currentProgram.nameSpacePrefix);
                    if (nmspace.IsValidIdentifier()) {
                        thePrintObject = FetchObjectFromVariableByType (&nmspace,HY_ANY_OBJECT);
                    }
                } else {
                    nmspace = AppendContainerName(*varname,currentProgram.nameSpacePrefix);
                }
                
                
                if (thePrintObject == nil) {
                    long typeFlag = HY_BL_ANY;
                    
                    thePrintObject = _HYRetrieveBLObjectByName (nmspace, typeFlag);
                    
                    if (!thePrintObject) {
                        _String argCopy = *varname,
                                errMsg;

                        _FormulaParsingContext fpc (&errMsg, currentProgram.nameSpacePrefix);
                        if (Parse (&f,argCopy, fpc, nil) == HY_FORMULA_EXPRESSION) {
                            thePrintObject = f.Compute(0,currentProgram.nameSpacePrefix);
                        } else {
                            if (errMsg.sLength)
                                throw (errMsg);
                            else
                                throw (_String ("Argument ") & i & " is not a simple expression");
                        }
                    }
                }
            }
            
            if (thePrintObject) {
                if (!print_to_stdout) {
                    thePrintObject->toFileStr (dest);
                } else {
                    _String outS ((_String*)thePrintObject->toStr());
                    StringToConsole (outS);
                }
            }
        }
    }
    catch (_String errMsg) {
        currentProgram.ReportAnExecutionError (errMsg);
    }
    
#if !defined __UNIX__ || defined __HEADLESS__ || defined __HYPHYQT__ || defined __HYPHY_GTK__
    if (print_to_stdout) {
        yieldCPUTime();
    }
#endif
    if (dest && dest!=globalMessageFile && doClose) {
        fclose (dest);
    }
    
    return !currentProgram.IsErrorState();
}


