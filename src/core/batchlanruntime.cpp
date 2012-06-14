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

#if defined __MAC__ || defined __WINDOZE__ || defined HYPHY_GTK
    #include "HYConsoleWindow.h"
    #include "HYDialogs.h"

#endif

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleHarvestFrequencies (_ExecutionList& currentProgram) {
    currentProgram.currentCommand++;
    SetStatusLine           ("Gathering Frequencies");
    _String  freqStorageID = currentProgram.AddNameSpaceToID(*(_String*)parameters(0)),
             dataID        = currentProgram.AddNameSpaceToID(*(_String*)parameters(1));
    
    long       objectType = HY_BL_DATASET|HY_BL_DATASET_FILTER;
    BaseRef    sourceObject = _HYRetrieveBLObjectByName (dataID, objectType,nil,true);
    
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
            vSpecs = ProcessLiteralArgument((_String*)parameters(5),currentProgram.nameSpacePrefix);
        }
        if (parameters.lLength>6)   {
            hSpecs = ProcessLiteralArgument((_String*)parameters(6),currentProgram.nameSpacePrefix);
        }
        
        _DataSet * dataset = (_DataSet*)sourceObject;
        _SimpleList     hL, vL;
        dataset->ProcessPartition (hSpecs,hL,false);
        dataset->ProcessPartition (vSpecs,vL,true);
        
        receptacle = dataset->HarvestFrequencies(unit,atom,posspec,hL, vL,cghf>0.5);
    } else { // harvest from a DataSetFilter
        
        receptacle = ((_DataSetFilter*)sourceObject)->HarvestFrequencies(unit,atom,posspec,cghf>0.5);
    }
    
    SetStatusLine           (empty);
    return CheckReceptacleCommandIDAndStore (&freqStorageID,HY_HBL_COMMAND_HARVEST_FREQUENCIES,true, receptacle, false);
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
            #if not defined __AFYP_REWRITE_BGM__
                _Matrix                * optRes;
                #ifdef __AFYP_DEVELOPMENT__
                    _SimpleList *   first_order = nil;
                     optRes = (_Matrix *) lkf->CovarianceMatrix (first_order);
                        
                #else
                    optRes = (_Matrix*)lkf->CovarianceMatrix(nil);
                #endif
                    if (optRes) {
                        result->SetValue(optRes,false);
                    }
             #endif            
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
#ifdef __UNIX__
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
#ifndef __UNIX__
            _SimpleList choiceDummy (2,0,1), selDummy;
            model_id = HandleListSelection (templateModelList, choiceDummy, matchingModels, "Choose one of the standard substitution models",selDummy,1);
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

    if (currentArgument->Equal (&deferConstrainAssignment)) {
        bool on = ProcessNumericArgument ((_String*)parameters(1), currentProgram.nameSpacePrefix);
        if (on) {
            deferSetFormula = (_SimpleList*)checkPointer(new _SimpleList);
        } else if (deferSetFormula) {
            FinishDeferredSF ();
        }
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
                        ((_BayesianGraphicalModel *)lkf)->SetDataMatrix ((_Matrix *) dataMx->makeDynamic());
                    } else {
                        WarnError (_String("Data matrix columns (") & dataMx->GetVDim() & " ) does not match number of nodes in graph (" & num_nodes & ").");
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
                        WarnError("Dimension of graph does not match current graph.");
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
                        WarnError ("Dimensions of constraint matrix do not match current graph.");
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
                        WarnError("Length of order vector doesn't match number of nodes in graph.");
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
                WarnError (*currentArgument & " is not a valid BGM parameter in call to " & _HY_ValidHBLExpressions.RetrieveKeyByPayload(HY_HBL_COMMAND_SET_PARAMETER));
                return false;
            }
        } // end BGM
        break;
        
        case HY_BL_SCFG:
        case HY_BL_LIKELIHOOD_FUNCTION:
        {
            if (typeFlag == HY_BL_SCFG && currentArgument->Equal (&scfgCorpus)) {
                ((Scfg*)theObject)->SetStringCorpus ((_String*)parameters(2));
            } else {
                _LikelihoodFunction * lkf = (_LikelihoodFunction *) theObject;
                currentArgument = (_String*)parameters(1);
                long g = ProcessNumericArgument(currentArgument,currentProgram.nameSpacePrefix);
                if (g < 0 || g >= lkf->GetIndependentVars().lLength) {
                    WarnError (*currentArgument & " (=" & g & ") is not a valid parameter index in call to " & _HY_ValidHBLExpressions.RetrieveKeyByPayload(HY_HBL_COMMAND_SET_PARAMETER));
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
                WarnError (*((_String*)parameters(1)) & " (=" & f & ") is not a valid sequence index in call to " & _HY_ValidHBLExpressions.RetrieveKeyByPayload(HY_HBL_COMMAND_SET_PARAMETER));
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
                        //treeNode->SetModel ();
                    }
                    else {
                        WarnError (*((_String*)parameters(2)) & " does not appear to be a valid model name in call to "& _HY_ValidHBLExpressions.RetrieveKeyByPayload(HY_HBL_COMMAND_SET_PARAMETER));
                        return false;
                    }
                }
            }
            
            WarnError (*currentArgument & " is not a valid likelihood function/data set filter/tree topology in call to "& _HY_ValidHBLExpressions.RetrieveKeyByPayload(HY_HBL_COMMAND_SET_PARAMETER));
            return false;

    
    } // end cases
    return true;
}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleAssert (_ExecutionList& currentProgram) {
    currentProgram.currentCommand++;

    _String assertion (*(_String*)parameters(0));

    _Formula rhs, lhs;
    long     varRef;
    if (Parse (&rhs,assertion,varRef,currentProgram.nameSpacePrefix,&lhs) == HY_FORMULA_EXPRESSION) {
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
                    WarnError (errMsg);
                }
            }
            return true;
        }
    }
 
    WarnError(_String("Assertion statement '") & *(_String*)parameters(0) & "' could not be computed or was not numeric. " & _HY_ValidHBLExpressions.RetrieveKeyByPayload(HY_HBL_COMMAND_ASSERT));
    
    return false;
}

//____________________________________________________________________________________

bool      _ElementaryCommand::HandleRequireVersion(_ExecutionList& currentProgram){
    currentProgram.currentCommand++;
    _String theVersion = ProcessLiteralArgument ((_String*)parameters (0),currentProgram.nameSpacePrefix);

    if (__KERNEL__VERSION__.toNum() < theVersion.toNum()) {
        WarnError (_String ("Current batch file requires at least version :")& theVersion &" of HyPhy. Please download an updated version from http://www.hyphy.org and try again.");
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



