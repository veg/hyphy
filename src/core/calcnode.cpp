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

#include <math.h>
#include <ctype.h>
#include <float.h>

#include "calcnode.h"
#include "function_templates.h"

#include "category.h"
#include "batchlan.h"
#include "tree.h"

#include "global_things.h"
#include "hbl_env.h"

using namespace hy_global;
using namespace hy_env;

//#define _UBER_VERBOSE_MX_UPDATE_DUMP
//#define _UBER_VERBOSE_MX_UPDATE_DUMP_EVAL 1

//_______________________________________________________________________________________________

_CalcNode::_CalcNode    () {
    theProbs = nil;    // default constructor, doesn't do much
    compExp = nil;
    matrixCache = nil;
    flags = 0;
    cBase = 0;
}

//_______________________________________________________________________________________________
_CalcNode::_CalcNode    (_String name, _String parms, int codeBase, _VariableContainer* theP, _AVLListXL* aCache):_VariableContainer (name, "", theP) {
    // construct a node from a string of the form
    // matrix name, <optional comma separated variable declarations, inititalizations>
    // also should be passed the pointer to a container tree
    InitializeCN (parms, codeBase, theP, aCache);
}

//_______________________________________________________________________________________________
_CalcNode::_CalcNode    (_CalcNode* sourceNode, _VariableContainer* theP):_VariableContainer (sourceNode->ContextFreeName(), "", theP) {
    _String model = *sourceNode->GetModelName();
    InitializeCN (model, 0, theP);
    if (model.nonempty()) { // copy model parameter values
        CopyMatrixParameters(sourceNode, true);
    }
}

//_______________________________________________________________________________________________
void    _CalcNode::InitializeCN     ( _String const& parms, int, _VariableContainer* theP, _AVLListXL* aCache) {
    if (theIndex < 0) return;

    cBase         = 0;
    flags         = 0;
    theProbs      = nil;
    compExp       = nil;
    matrixCache   = nil;
    
    _List          parameters;
    _ElementaryCommand::ExtractConditions(parms, 0, parameters, ',');

    InitializeVarCont (kEmptyString, *(_String*)parameters.GetItem(0), theP, aCache);
    
    // parse and instantiate all model expressions
     parameters.ForEach([=] (BaseRef expression, unsigned long) -> void {
         _Formula fg (*(_String*)expression, this);
        }, (GetModelIndex() == HY_NO_MODEL && parms.nonempty()) ? 0L : 1L);
    
    // attach the variables to the lists inside the node
    ScanAndAttachVariables();

    // check for category variables
    
    //_SimpleList local_cat_vars;
    
    ForEachLocalVariable(iVariables, [&] (long var_idx, long ref_idx, unsigned long idx) -> void {
        if (ref_idx >= 0) {
            /* TODO:
             this used to do with local category variables
             NOT TESTED and is LIKELY BROKEN
             replaced with an error message;
             */
             _Variable *test_var = (_Variable*)LocateVar(ref_idx);
             if (test_var->IsCategory()) {
                 HandleApplicationError(_String ("Local instance to category variables are not supported, ") & test_var->GetName()->Enquote());
                 return;
             }
        }
    });
    
    if (gVariables) {
        
        _SimpleList moved_to_cat;
        
        gVariables->Each ([&] (long var_idx, unsigned long idx) -> void {
            _Variable *test_v = LocateVar(var_idx);
            if (test_v->IsCategory()) {
                categoryVariables << var_idx;
                categoryIndexVars << -1;
                moved_to_cat << idx;
            }
        });
        
        if (moved_to_cat.nonempty()) {
            gVariables->DeleteList(moved_to_cat);
            if (gVariables->empty()) {
                DeleteAndZeroObject(gVariables);
            }
        }
    }
    
    variablePtrs.Replace(theIndex, this, true);
    
        /** TODO check equivalence **
         
            BaseRef temp =  variablePtrs(theIndex);
            variablePtrs[theIndex]=this->makeDynamic();
            DeleteObject(temp);
        */

}

//__________________________________________________________________________________
void    _CalcNode::SetModel (long modelID, _AVLListXL* varCache) {
   _VariableContainer::SetModel (modelID, varCache);
}

//_______________________________________________________________________________________________

long      _CalcNode::SetDependance (long var_index) {
    var_index = _VariableContainer::SetDependance (var_index);
    if (var_index >= 0L) {
        /** if the new constraint includes a category variable,
            it needs to be marked as such in this CalcNode
         */
        
 
        PopulateAndSort ([&] (_AVLList& avl) -> void {
            LocateVar (var_index)->ScanForVariables (avl,true);
        }).Each ( [&] (long var_index, unsigned long idx) -> void {
            if (LocateVar(var_index)->IsCategory() &&(categoryVariables >> var_index)) {
                categoryIndexVars<<-1;
            }
        });

        // also clear out previously computed matrix exponentials
        if (compExp) {
            DeleteAndZeroObject(compExp);
        } /*else {
            if (matrixCache) {
            }
        }*/
    }
    return var_index;
}

//_______________________________________________________________________________________________

void    _CalcNode::SetCodeBase (int codeBase) {
    if (codeBase>0) {
        if (codeBase != cBase || !theProbs) {
            if (theProbs) {
                delete [] theProbs;
            }
            theProbs = new hyFloat [codeBase];
            cBase = codeBase;
        }
        
        theProbs[0] = 1.0;
    }
}

//_______________________________________________________________________________________________
void    _CalcNode::SetCompMatrix (long categID) {
    compExp = GetCompExp (categID);
}

//_______________________________________________________________________________________________

hyFloat  _CalcNode::ProcessTreeBranchLength (_String const& branch_length) {
  hyFloat res = -1.;

  if (branch_length.nonempty()) {
    if (branch_length.char_at(0UL)==':') {
      res = branch_length.Cut(1L,-1L).to_float();
    } else {
      res = branch_length.to_float ();
    }
    res = MAX (res, 1e-10);
  }

  return res;
}

//_______________________________________________________________________________________________

void _CalcNode::Clear (void) {
    if (compExp) {
        DeleteAndZeroObject (compExp);
    }
    if (theProbs) {
        delete [] theProbs;
        theProbs = nil;
    }
    _VariableContainer::Clear();
}

//_______________________________________________________________________________________________

_CalcNode::~_CalcNode (void) {
    Clear();
}

//_______________________________________________________________________________________________

long    _CalcNode::FreeUpMemory (long) {
    long res = 0L;
    if (compExp) {
        res = compExp->GetMySize();
        DeleteAndZeroObject (compExp);
    }
    return res;
}

//__________________________________________________________________________________

void _CalcNode::RemoveModel (void) {
 
    categoryVariables.Clear();
    categoryIndexVars.Clear();
    remapMyCategories.Clear();

    Clear();
}

//__________________________________________________________________________________
_String*            _CalcNode::GetBranchSpec (void) {
    
    _StringBuffer * res = new _StringBuffer (32UL);
    *res << GetModelName();
    
    ForEachLocalVariable(iVariables, [&] (long var_index, long ref_index, unsigned long idx) -> void {
        if (idx == 0UL) {
            (*res) << (res->nonempty() ? ',' : '{');
        } else {
            (*res) << ',';
        }
        _Variable * av = LocateVar (var_index);
        if (ref_index >= 0L) {
            res->AppendAnAssignmentToBuffer(LocateVar (ref_index)->GetName(),
                                            new _String (av->Value()));
        } else {
            res->AppendAnAssignmentToBuffer(av->GetName(),
                                            new _String (av->Value()));
        }
    });
    
    ForEachLocalVariable(dVariables, [&] (long var_index, long ref_index, unsigned long) -> void {
         if (ref_index < 0L) {
             (*res) << (res->nonempty() ? ',' : '{');
             
             _Variable * av = LocateVar (var_index);
             res->AppendAnAssignmentToBuffer(av->GetName(),
                                             av->GetFormulaString(kFormulaStringConversionNormal),
                                             kAppendAnAssignmentToBufferFree | kAppendAnAssignmentToBufferAssignment);
         }
    });
    
    if (res->nonempty()) {
        (*res) << '}';
    }
    
    res->TrimSpace();
    return res;
}


//__________________________________________________________________________________

void _CalcNode::ReplaceModel (_String& modelName, _VariableContainer* parent_tree) {
    // TODO: SLKP 20171203 this is FUGLY
  RemoveModel    ();

  _TheTree * parent_tree_object = (_TheTree*)parent_tree;

  long index_in_parent = parent_tree_object->flatCLeaves._SimpleList::Find ((long)this);
  _List * container_object = nil;

  if (index_in_parent >= 0) {
    parent_tree_object->flatCLeaves.Replace (index_in_parent, nil, false);
    container_object = &parent_tree_object->flatCLeaves;
  } else {
    index_in_parent = parent_tree_object->flatTree._SimpleList::Find ((long)this);
    if (index_in_parent >= 0) {
      parent_tree_object->flatTree.Replace (index_in_parent, nil, false);
      container_object = &parent_tree_object->flatTree;
    }
  }

  long my_var_index = theIndex;

  DeleteVariable (theIndex, false); // this will clean up all the node.xx variables
  InitializeCN   (modelName, 0 , parent_tree);

  if (container_object) {
    _Variable * new_node = LocateVar(my_var_index);
    container_object->Replace (index_in_parent, new_node, false);
    new_node->AddAReference();
  }
}

//_______________________________________________________________________________________________
bool    _CalcNode::MatchSubtree (_CalcNode* mNode) {
    // TODO: SLKP 20171203 possible deprecate
    node <long>* myNode    = LocateMeInTree (),
                 * matchNode = mNode->LocateMeInTree ();
    if (myNode&&matchNode) {
        return myNode->compare_subtree(matchNode);
    }
    return      false;
}

//_______________________________________________________________________________________________

hyFloat  _CalcNode::ComputeBranchLength (void) {
    
    static const    _String kLargeMatrixBranchLengthDimension ("LARGE_MATRIX_BRANCH_LENGTH_MODIFIER_DIMENSION"),
                            kLargeMatrixBranchLengthModifier  ("LARGE_MATRIX_BRANCH_LENGTH_MODIFIER");


    if (GetModelIndex() < 0) {
        return Value();
    }

    HBLObjectRef   stencil = (HBLObjectRef)hy_env::EnvVariableGet(hy_env::branch_length_stencil, STRING | ASSOCIATIVE_LIST);

    if (stencil) {
        if (stencil->ObjectClass () == STRING) {
            if (((_FString*)stencil)->get_str() == kStringSuppliedLengths) {
                return Value();
            }
        } else {
            _AssociativeList *lookup = (_AssociativeList*)stencil;
            _String lookup_name = ContextFreeName();
            _Constant * value = (_Constant*)lookup->GetByKey (lookup_name, NUMBER);
            if (value) {
                return value->Value();
            }
        }
    }
    

    _Matrix     *freqMx,
                *theMx;

    bool        mbf;

    RetrieveModelComponents (theModel, theMx, freqMx, mbf);

    if ( ! (freqMx && theMx)) {
        return Value();
    }

    freqMx = (_Matrix*)freqMx->ComputeNumeric();

    hyFloat              result = 0.0;

    IntergrateOverAssignments (categoryVariables, true, [&] (long current_cat, hyFloat weight) -> void {
        result += fabs(ComputeModelMatrix()->ExpNumberOfSubs (freqMx, mbf))*weight;
    });

    if (freqMx->GetSize()  > hy_env::EnvVariableGetNumber(kLargeMatrixBranchLengthDimension, 20.)) {
        result /= hy_env::EnvVariableGetNumber(kLargeMatrixBranchLengthModifier, 3.);
    }
    
    return result;
}

//_______________________________________________________________________________________________

hyFloat& _CalcNode::operator[] (unsigned long i) {
    // TODO SLKP 20171203, possibly DEPRECATE?
    return theProbs [i];
}

//_______________________________________________________________________________________________

BaseRef _CalcNode::toStr (unsigned long) {
    _StringBuffer * res = new _StringBuffer (64L);
    (*res) << theName << '(' << _String(CountIndependents()) << ',' << _String(CountDependents()) << ')';
    res->TrimSpace();
    return res;
}

//__________________________________________________________________________________

void    _CalcNode::Duplicate (BaseRefConst theO) {
    _VariableContainer::Duplicate (theO);
    cBase         = 0;
    compExp       = nil;
    matrixCache   = nil;
    theProbs      = nil;
}

//_______________________________________________________________________________________________

bool        _CalcNode::HasChanged(bool) {
    
    return _VariableContainer::HasChanged() || categoryVariables.Any([&] (long cat_idx, unsigned long) -> bool {
        return LocateVar (cat_idx)->HasChanged();
    });
 }


//_______________________________________________________________________________________________

bool        _CalcNode::NeedNewCategoryExponential(long catID) const {

    if (_VariableContainer::NeedToExponentiate(catID>=0)) {
        return true;
    }

    if (catID==-1) {
        return !compExp || categoryVariables.Any([&] (long cat_idx, unsigned long) -> bool {
            return LocateVar (cat_idx)->HasChanged();
        });
    } else {
        return !GetCompExp(catID) || categoryVariables.Any([&] (long cat_idx, unsigned long i) -> bool {
            return ((_CategoryVariable*)LocateVar (cat_idx))->HaveParametersChanged(remapMyCategories.list_data[catID*(categoryVariables.countitems()+1)+i+1]);
        });
    }
    return false;
}

//_______________________________________________________________________________________________
bool        _CalcNode::RecomputeMatrix  (long categID, long totalCategs, _Matrix* storeRateMatrix, _List* queue, _SimpleList* tags, _List* bufferedOps)
{
    // assumed that NeedToExponentiate was called prior to this function

    //_Variable* curVar, *locVar;
    
    _SimpleList * var_lists [2] = {iVariables, dVariables};
    for (_SimpleList* iterable : var_lists) {
        ForEachLocalVariable(iterable,CopyModelParameterValue);
    }
  
    #ifdef _UBER_VERBOSE_MX_UPDATE_DUMP
     // if (1|| likeFuncEvalCallCount == _UBER_VERBOSE_MX_UPDATE_DUMP_LF_EVAL && gVariables) {
        for (unsigned long i=0; i<gVariables->lLength; i++) {
          _Variable* curVar = LocateVar(gVariables->GetElement(i));
          fprintf (stderr, "[_CalcNode::RecomputeMatrix] Node %s, var %s, value = %15.12g\n", GetName()->get_str(), curVar->GetName()->get_str(), curVar->Compute()->Value());
        }
      //}
    #endif

    /*
    for (unsigned long i=0; i<categoryVariables.lLength; i++) {
        if (categoryIndexVars.list_data[i]<0) {
            continue;
        }
        curVar = LocateVar (categoryIndexVars.list_data[i]);
        locVar = LocateVar (categoryVariables.list_data[i]);
        curVar->SetValue(locVar->Compute());
    } */

    if (!storeRateMatrix) {
      if (totalCategs>1) {
  #ifdef _UBER_VERBOSE_MX_UPDATE_DUMP
        fprintf (stderr, "[_CalcNode::RecomputeMatrix] Deleting category %ld for node %s at %p\n", categID, GetName()->get_str(), GetCompExp(categID));
  #endif
        if (clear_exponentials()) {
            DeleteObject(GetCompExp(categID, true));
        }
      } else {
        if (compExp) {
            if (clear_exponentials ()) {
                DeleteAndZeroObject(compExp);
            }
        }
      }
    }

    bool    isExplicitForm  = HasExplicitFormModel ();

    if (isExplicitForm && bufferedOps) {
        _Matrix * bufferedExp = (_Matrix*)GetExplicitFormModel()->Compute (0,nil, bufferedOps);
        #ifdef _UBER_VERBOSE_MX_UPDATE_DUMP
            fprintf (stderr, "[_CalcNode::RecomputeMatrix] Setting (buffered) category %ld/%ld for node %s\n", categID, totalCategs, GetName()->get_str());
         #endif
        SetCompExp ((_Matrix*)bufferedExp->makeDynamic(), totalCategs>1?categID:-1);
        return false;
    }

    unsigned long      previous_length = queue && tags ? queue->lLength: 0;
    _Matrix * myModelMatrix = GetModelMatrix(queue,tags);

    if (isExplicitForm && !myModelMatrix) { // matrix exponentiations got cached
        if (queue && queue->lLength > previous_length) {
            return true;
        } else {
            HandleApplicationError ("Internal error");
            return false;
        }
    } else {

        if (myModelMatrix->MatrixType()!=_POLYNOMIAL_TYPE) {
            _Matrix *temp = nil;
            if (isExplicitForm) {
                temp = (_Matrix*)myModelMatrix->makeDynamic();
            } else {
                temp = (_Matrix*)myModelMatrix->MultByFreqs(theModel, true);
            }
            
            // copy updated model (local) constrained parameters to their external references
            ForEachLocalVariable(dVariables, [&] (long var_idx, long ref_idx, unsigned long) -> void {
                if (ref_idx >= 0) {
                    _Variable * model_var = LocateVar (ref_idx);
                    if (!model_var -> IsIndependent()) {
                        _Variable * param_var = LocateVar (var_idx);
                        if (param_var->IsIndependent()) {
                            param_var->SetValue(param_var->Compute(),true,true,NULL);
                        }
                    }
                }
            });
            

            if (storeRateMatrix) {
                storeRateMatrix->Duplicate(temp);
                return isExplicitForm;
            }


            if (queue) {
                (*queue) << temp;
                if (tags) {
                    (*tags) << (isExplicitForm);
                }
                return isExplicitForm;
            }

            #ifdef _UBER_VERBOSE_MX_UPDATE_DUMP
                fprintf (stderr, "[_CalcNode::RecomputeMatrix] Setting category %ld/%ld for node %s\n", categID, totalCategs, GetName()->get_str());
            #endif
            SetCompExp ((_Matrix*)(isExplicitForm?temp:temp->Exponentiate(1., true)), totalCategs>1?categID:-1);

        } else {
            compExp = (_Matrix*)myModelMatrix->Evaluate(false);
        }
    }
    return false;
}

//_______________________________________________________________________________________________
void        _CalcNode::SetCompExp  (_Matrix* m, long catID, bool do_exponentiation) {
    
    _Matrix ** store_exp_here;
    
    if (catID >= 0 && matrixCache) {
        if (remapMyCategories.lLength) {
            catID = remapMyCategories.list_data[catID*(categoryVariables.lLength+1)];
        }
        store_exp_here = &(matrixCache[catID]);
    } else {
        store_exp_here = &compExp;
    }
    
    if (do_exponentiation) {
        compExp = m->Exponentiate(1., true, *store_exp_here);
        reuse_exponentials ();
    } else {
        compExp = m;
    }
    *store_exp_here = compExp;
}
//_______________________________________________________________________________________________
_Matrix*        _CalcNode::ComputeModelMatrix  (bool) {
    // assumed that NeedToExponentiate was called prior to this function
    _SimpleList * var_lists [2] = {iVariables, dVariables};
    for (_SimpleList* iterable : var_lists) {
        ForEachLocalVariable(iterable,CopyModelParameterValue);
    }

    _Matrix * modelMx = GetModelMatrix();
    if (modelMx && modelMx->ObjectClass()==MATRIX && modelMx->MatrixType()!=_POLYNOMIAL_TYPE) {
        return (_Matrix*)modelMx->ComputeNumeric();
    }

    return nil;
}

//_______________________________________________________________________________________________

_Matrix*    _CalcNode::GetCompExp       (long catID, bool doClear) const {
    if (catID==-1) {
        return compExp;
    } else {

        if (remapMyCategories.lLength) {
            catID = remapMyCategories.list_data[catID * (categoryVariables.lLength+1)];
        }

        _Matrix* ret = matrixCache?matrixCache[catID]:compExp;

        if (doClear && matrixCache) {
            if (matrixCache[catID] && matrixCache[catID]->CanFreeMe()) {
                matrixCache[catID] = nil;
            } 
        }
        return ret;
    }
}

//_______________________________________________________________________________________________

BaseRef     _CalcNode::makeDynamic(void) const {
    _CalcNode* res = new (_CalcNode);
    res->_VariableContainer::Duplicate (this);
    res->categoryVariables.Duplicate (&categoryVariables);
    res->categoryIndexVars.Duplicate (&categoryIndexVars);
    res->theValue = theValue;
    res->cBase = cBase;
    res->flags = flags;
    if (cBase) {
        res->theProbs = new hyFloat [cBase];
        CopyArray(res->theProbs, theProbs, cBase);
    } else {
        res->theProbs = nil;
    }
    res->compExp = compExp;
    if (compExp) {
        compExp->AddAReference();
    }
    return res;
}

//_______________________________________________________________________________________________

void     _CalcNode::SetupCategoryMap (_List& containerVariables, _SimpleList& classCounter, _SimpleList& multipliers) {
    // TODO 20171203, SLKP: this needs review
    
    long    totalCategories = classCounter.Element(-1),
            globalCatCount  = containerVariables.countitems(),
            localCategories = 1L,
            catCount        = categoryVariables.countitems(),
            entriesPerCat   = 1L+catCount;
    
    //for (long k = 0; k<categoryVariables.lLength;k++)
    //    printf ("%ld\n", categoryVariables(k));//, ((_Variable*)categoryVariables(k))->GetName()->sData);
    
    if (catCount == 0L) {
        remapMyCategories.Clear();
    } else {
        remapMyCategories.Populate (totalCategories*entriesPerCat,0,0);
        
        _SimpleList     remappedIDs,
                        rateMultiplers (catCount,1,0),
                        categoryValues (globalCatCount,0,0);
        
        /* find where in the global list the categories for this node are */
        for (long myCatID = 0L; myCatID < catCount; myCatID++) {
            long index = containerVariables.FindPointer(LocateVar(categoryVariables.get(myCatID)));
            if (index < 0) {
                HandleApplicationError ("Internal error in SetupCategoryMap. Please report to sergeilkp@icloud.com");
            }
            localCategories *= classCounter.get(index);
            remappedIDs << index;
        }
        
        /* generate multiplicative offsets for each category; a move by one category value
           in variable myCatID changes the compositve index by "offset"
         */
        for (long myCatID = catCount-2L; myCatID >= 0L; myCatID--) {
            rateMultiplers.list_data[myCatID] = rateMultiplers.list_data[myCatID+1]*classCounter.list_data[remappedIDs.list_data[myCatID+1]];
        }
        
        for (long currentRateCombo  = 0L; currentRateCombo < totalCategories; currentRateCombo++) {
            long copyRateCombo = currentRateCombo;
            for (long variableID = 0L; variableID < globalCatCount; variableID++) {
                categoryValues.list_data[variableID] = copyRateCombo / multipliers.list_data[variableID];
                copyRateCombo = copyRateCombo%multipliers.list_data[variableID];
                //printf ("%d %d %d %d\n", currentRateCombo, variableID, multipliers.list_data[variableID], categoryValues.list_data[variableID]);
            }
            
            long localCatID = 0L;
            
            for  (long localVariableID = 0; localVariableID<catCount; localVariableID++) {
                localCatID += rateMultiplers.list_data[localVariableID] * categoryValues.list_data[remappedIDs.list_data[localVariableID]];
            }
            
            long offset = currentRateCombo * entriesPerCat;
            remapMyCategories.list_data[offset] = localCatID;
            //printf ("[%ld] = %ld (%ld)\n", offset, localCatID, );
            
            offset++;
            for  (long localVariableID = 0; localVariableID<catCount; localVariableID++) {
                remapMyCategories[offset++] = categoryValues.list_data[remappedIDs.list_data[localVariableID]];
            }
            
        }
    }
    
    //printf ("Node remap at %s yielded %s\n", GetName()->sData, _String((_String*)remapMyCategories.toStr()).sData);
    
}

//_______________________________________________________________________________________________

node<long>* _CalcNode::LocateMeInTree (void) const {
    
    _String parentName = ParentObjectName (),
    myName     = ContextFreeName();
    
    return  ((_TreeTopology*)FetchVar(LocateVarByName(parentName)))->FindNodeByName(&myName);
    
}


//_______________________________________________________________________________________________

void _CalcNode::ConvertToSimpleMatrix (void) const {
    _Formula * mf = GetExplicitFormModel();
    if (mf) {
        mf->ConvertMatrixArgumentsToSimpleOrComplexForm (false);
    } else {
        _Matrix * mm [2] = {GetModelMatrix(), GetFreqMatrix()};
        for (_Matrix * m : mm) {
            if (m) {
                m->MakeMeSimple();
            }
        }
     }
}

//_______________________________________________________________________________________________

void _CalcNode::ConvertFromSimpleMatrix (void) {
    _Formula * mf = GetExplicitFormModel();
    if (mf) {
        mf->ConvertMatrixArgumentsToSimpleOrComplexForm (true);
    } else {
        _Matrix * mm [2] = {GetModelMatrix(), GetFreqMatrix()};
        for (_Matrix * m : mm) {
            if (m) {
                m->MakeMeGeneral();
            }
        }
    }
}

//_______________________________________________________________________________________________

_Formula*   _CalcNode::RecurseMC (long varToConstrain, node<long>* whereAmI, bool first, char rooted) {
    // TODO 20171203, SLKP: this needs review

    long descendants = whereAmI->get_num_nodes(),
    f = iVariables?iVariables->FindStepping(varToConstrain,2,1):-1,
    start = 0;
    
    if (f<0 && !first) {
        HandleApplicationError (_String ("Molecular clock constraint has failed, since variable ")
                                &LocateVar(varToConstrain)->GetName()->Enquote()
                                &" is not an independent member of the node "
                                &GetName()->Enquote()
                                );
        return nil;
    }
    
    
    if (descendants == 0) {
        if (first) {
            return nil;
        } else {
            return new _Formula (LocateVar(iVariables->get(f-1)),true);
        }
    }
    
    if (first && (!whereAmI->get_parent()) && (rooted == ROOTED_LEFT)) {
        descendants --;
    }
    if (first && (!whereAmI->get_parent()) && (rooted == ROOTED_RIGHT)) {
        start++;
    }
    
    // internal node - must do some work
    
    _Formula**  nodeConditions = new _Formula * [descendants-start];
    
    for (long k=start+1; k<=descendants; k++) {
        node<long>* downWeGo = whereAmI->go_down(k);
        if (!(nodeConditions[k-1-start] = map_node_to_calcnode(downWeGo)->RecurseMC (varToConstrain, downWeGo))) {
            for (long f2 = 0; f2 < k-start-1; f2++) {
                delete nodeConditions[f2];
            }
            
            delete[] nodeConditions;
            return nil;
        }
    }
    
    // all the conditions have been written. now check how we should resolve them
    
    long k;
    
    for (k=0; k<descendants-start; k++)
        if ((nodeConditions[k])->Length()>1) {
            break;
        }
    
    if (k==descendants-start) { // all underlying branches are "simple"
        for (long n=1; n<descendants-start; n++) {
            //printf ("Setting simple constraint at %s %s\n", LocateVar(nodeConditions[n]->GetIthTerm(0)->GetAVariable())->GetName()->get_str(),
            //        _String((_String*)nodeConditions[0]->toStr(kFormulaStringConversionNormal, nil, true)).get_str());
            LocateVar (nodeConditions[n]->GetIthTerm(0)->GetAVariable())->SetFormula (*nodeConditions[0]);
            delete (nodeConditions[n]);
            nodeConditions[n] = nil;
        }
        k = 0;
    } else {
        long l;
        for ( l=k+1; l<descendants-start; l++)
            if (nodeConditions[l]->Length()>1) {
                break;
            }
        
        if (l==descendants-start) // all but one underlying branches are "simple"
            for (long n=0; n<descendants-start; n++) {
                if (n==k) {
                    continue;
                }
                //printf ("Setting semi-simple constraint at %s : %s\n", LocateVar(nodeConditions[n]->GetIthTerm(0)->GetAVariable())->GetName()->get_str(),
                //                                                       _String((_String*)nodeConditions[k]->toStr(kFormulaStringConversionNormal, nil, true)).get_str());
                
                LocateVar (nodeConditions[n]->GetIthTerm(0)->GetAVariable())->SetFormula (*nodeConditions[k]);
                delete (nodeConditions[n]);
                nodeConditions[n] = nil;
            }
        // really bad bongos! must solve for non-additive constraint
        else
            for (long l=0; l<descendants-start; l++) {
                if (l==k) {
                    continue;
                }
                if (nodeConditions[l]->Length()==1) {
                    //printf ("Setting simple non-additive at %s %s\n", LocateVar(nodeConditions[l]->GetIthTerm(0)->GetAVariable())->GetName()->get_str(),
                    //        _String((_String*)nodeConditions[k]->toStr(kFormulaStringConversionNormal, nil, true)).get_str());
                    
                    LocateVar (nodeConditions[l]->GetIthTerm(0)->GetAVariable())->SetFormula (*nodeConditions[k]);
                } else { // solve for a non-additive constraint
                    _Variable* nonAdd = LocateVar (nodeConditions[l]->GetIthTerm(0)->GetAVariable());
                    nodeConditions[l]->GetList().Delete(0);
                    _Formula  newConstraint;
                    newConstraint.Duplicate(nodeConditions[k]);
                    for (long m=0; m<nodeConditions[l]->GetList().lLength; m++) {
                        _Operation* curOp = (_Operation*)(*nodeConditions[l]).GetList()(m);
                        if (curOp->GetNoTerms()) {
                            newConstraint.GetList().AppendNewInstance(new _Operation (HY_OP_CODE_SUB, 2L));
                        } else {
                            newConstraint.GetList()<<curOp;
                        }
                    }
                    delete (nodeConditions[l]);
                    nodeConditions[l] = nil;
                    nonAdd->SetFormula(newConstraint);
                    //printf ("Setting complex non-additive at %s %s\n", nonAdd->GetName()->get_str(),
                    //        _String((_String*)newConstraint.toStr(kFormulaStringConversionNormal, nil, true)).get_str());
                }
            }
    }
    
    
    if(!first) {
        _Formula     *result = nodeConditions[k];
        _Operation   *newVar = new _Operation;
        newVar->SetAVariable (iVariables->list_data[f-1]);
        result->GetList().AppendNewInstance( newVar);
        result->GetList().AppendNewInstance(new _Operation (HY_OP_CODE_ADD, 2L));
        
        delete [] nodeConditions;
        return result;
    }
    
    for (long k=0; k<descendants-start; k++)
        if (nodeConditions[k]) {
            delete nodeConditions[k];
        }
    
    delete [] nodeConditions;
    return nil;
}

//_______________________________________________________________________________________________

_VariableContainer*     _CalcNode::ParentTree(void) {
    _String parentTree = ParentObjectName();
    return (_VariableContainer* )FetchObjectFromVariableByType(&parentTree, TREE);
}





