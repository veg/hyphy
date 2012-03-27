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

#include      "defines.h"
#include      "batchlan.h"
#include      "likefunc.h"

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
    _Matrix                * optRes;
    
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
                    for (long iid = 0; iid < restrictedVariables->lLength; iid++) {
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
                    for (long var_index = 0; var_index < variableIDs.lLength; var_index++) {
                        long vID = lkf->GetIndependentVars().Find(variableIDs.lData[var_index]);
                        if (vID >= 0) (*restrictor) << vID;
                    }
                    if (restrictor->lLength == 0) {
                        DeleteObject (restrictor);
                        restrictor = nil;
                    }
               }
            }   
            if (optRes) {
                result->SetValue( (_Matrix*)lkf->CovarianceMatrix(restrictor),false);
            }
            DeleteObject (restrictor);
        } else {
        // BGM
            #if not defined __AFYP_REWRITE_BGM__
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
        optRes = lkf->Optimize();
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
