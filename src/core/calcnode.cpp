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

#include "global_things.h"

using namespace hy_global;

//_______________________________________________________________________________________________

_CalcNode::_CalcNode    () {
    theProbs = nil;    // default constructor, doesn't do much
    compExp = nil;
    matrixCache = nil;
}

//_______________________________________________________________________________________________
_CalcNode::_CalcNode    (_String name, _String parms, int codeBase, _VariableContainer* theP, _AVLListXL* aCache):_VariableContainer (name, "", theP) {
    // construct a node from a string of the form
    // matrix name, <optional comma separated variable declarations, initializations>
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
    theProbs      = nil;
    compExp       = nil;
    matrixCache   = nil;
    
    _List          parameters;
    _ElementaryCommand::ExtractConditions(parms, 0, parameters, ',');


    InitializeVarCont (kEmptyString, *(_String*)parameters.GetItem(0), theP, aCache);

    
    
     parameters.ForEach([=] (BaseRef expression) -> void {
         _Formula fg (*(_String*)expression, this);
        }, (GetModelIndex() == HY_NO_MODEL && parms.nonempty()) ? 0L : 1L);
    
    // attach the variables to the lists inside the node
    ScanAndAttachVariables();

    // check for category variables
    if (iVariables) {
        for
        for (f = iVariables->lLength-2; f>=0 && iVariables->lData[f+1] >= 0; f-=2) {
            _Variable *theV = LocateVar(iVariables->lData[f+1]);
            if (theV->IsCategory()) {
                /* TODO:
                    this has to do with local category variables;
                    NOT TESTED and is LIKELY BROKEN
                */
                _CategoryVariable* theCV = (_CategoryVariable*)theV;

                _Formula           newDensity,
                                   newCumulative;

                _SimpleList        iv,
                                   iv2,
                                   dv,
                                   dv2;

                for (unsigned long k = 0; k<iVariables->lLength; k+=2) {
                    iv  << iVariables->lData[k];
                    iv2 << iVariables->lData[k+1];
                }

                if (dVariables)
                    for (unsigned long k = 0; k<dVariables->lLength; k+=2) {
                        dv  << dVariables->lData[k];
                        dv2 << dVariables->lData[k+1];
                    }


                newDensity.LocalizeFormula    (theCV->GetDensity(),   *GetName(), iv, iv2, dv,dv2);
                newCumulative.LocalizeFormula (theCV->GetCumulative(),*GetName(), iv, iv2, dv,dv2);

                _CategoryVariable newCV;
                newCV.Duplicate (theCV);
                newCV.GetDensity().Duplicate(&newDensity);
                newCV.GetCumulative().Duplicate(&newCumulative);

                theV = LocateVar(iVariables->lData[f]);
                newCV.GetName()->Duplicate (theV->GetName());
                ReplaceVar(&newCV);

                categoryVariables<<iVariables->lData[f];
                categoryIndexVars<<iVariables->lData[f+1];
                iVariables->Delete(f);
                iVariables->Delete(f);
            }
        }

        if (iVariables->lLength) {
            iVariables->TrimMemory();
        } else {
            DeleteAndZeroObject(iVariables);
        }
    }
    if (gVariables) {
        for (f = gVariables->lLength-1; f>=0; f--) {
            _Variable *theV = LocateVar(gVariables->lData[f]);
            if (theV->IsCategory()) {
                categoryVariables<<gVariables->lData[f];
                categoryIndexVars<<-1;
                gVariables->Delete(f);
            }
        }
        if (gVariables->lLength) {
            gVariables->TrimMemory();
        } else {
            delete (gVariables);
            gVariables = nil;
        }
    }

    BaseRef temp =  (variablePtrs(theIndex));
    variablePtrs[theIndex]=this->makeDynamic();
    DeleteObject(temp);

}

//__________________________________________________________________________________
void    _CalcNode::SetModel (long modelID, _AVLListXL* varCache)
{
   _VariableContainer::SetModel (modelID, varCache);
}

//_______________________________________________________________________________________________

long      _CalcNode::SetDependence (long varIndex) {
    varIndex = _VariableContainer::SetDependence (varIndex);
    if (varIndex >= 0) {
        _SimpleList checkVars;
        _AVLList    myVars (&checkVars);
        LocateVar (varIndex)->ScanForVariables (myVars,true);

        for (long k=0; k<checkVars.lLength; k++)
            if (LocateVar(checkVars.lData[k])->IsCategory() &&(categoryVariables >> checkVars.lData[k])) {
                categoryIndexVars<<-1;
            }

        // also clear out category variables
        if (compExp) {
            DeleteAndZeroObject(compExp);
        } else {
            if (matrixCache) {
            }
        }
    }
    return varIndex;
}

//_______________________________________________________________________________________________

void    _CalcNode::SetCodeBase (int codeBase)
{
    if (codeBase>0) {
        if ((codeBase != cBase)||!theProbs) {
            if (theProbs) {
                delete theProbs;
            }
            theProbs = new hyFloat [codeBase];
            cBase = codeBase;
            theProbs[0]=1.0;
        } else {
            theProbs[0]=1.0;
        }
    }
}

//_______________________________________________________________________________________________
void    _CalcNode::SetCompMatrix (long categID)
{
    compExp = GetCompExp (categID);
}

//_______________________________________________________________________________________________

hyFloat  _CalcNode::ProcessTreeBranchLength (_String const& branch_length) {
  hyFloat res = -1.;

  if (!branch_length.empty()) {
    if (branch_length.char_at(0UL)==':') {
      res = branch_length.Cut(1L,-1L).to_float();
    } else {
      res = branch_length.to_float ();
    }

    res = MAX (res, 1e-10);
    if (res < 1e-10) {
      res = 1e-10;
    }
  }

  return res;
}

//_______________________________________________________________________________________________

_CalcNode::~_CalcNode (void) {

    if (theProbs) {
        delete [] theProbs;
    }
    if (compExp && referenceNode < 0) {
        DeleteObject (compExp);
    }
}

//_______________________________________________________________________________________________

long    _CalcNode::FreeUpMemory (long)
{
    long res = 0;
    if (compExp && referenceNode < 0) {
        res = compExp->GetMySize();
        DeleteObject (compExp);
        compExp = nil;
    }
    return res;
}

//__________________________________________________________________________________

void _CalcNode::RemoveModel (void)
{

    if (compExp && referenceNode < 0) {
        DeleteAndZeroObject(compExp);
        compExp = nil;
    }

    if (matrixCache) {
    }

    categoryVariables.Clear();
    categoryIndexVars.Clear();
    remapMyCategories.Clear();

    Clear();

}

//__________________________________________________________________________________
_String*            _CalcNode::GetBranchSpec (void) {
    
    _StringBuffer * res = new _StringBuffer (32UL);
    *res << GetModelName();
    
    if (iVariables && iVariables->lLength) {
        (*res) << (res->nonempty() ? ',' : '{');
        
        
        for (unsigned long k=0UL; k < iVariables->lLength; k+=2UL) {
            if (k) {
                (*res) << ',';
            }
            
            _Variable * av = LocateVar (iVariables->lData[k]);
            if (iVariables->lData[k+1UL] >= 0L) {
                res->AppendAnAssignmentToBuffer(LocateVar (iVariables->lData[k+1UL])->GetName(),
                                                new _String (av->Value()));
            } else {
                res->AppendAnAssignmentToBuffer(av->GetName(),
                                                new _String (av->Value()));
            }
        }
    }
    
    if (dVariables && dVariables->lLength) {
        for (unsigned long k=0UL; k < dVariables->lLength; k+=2UL) {
            if (dVariables->lData[k+1UL] <= 0L) {
                (*res) << (res->nonempty() ? ',' : '{');
                
                _Variable * av = LocateVar (dVariables->lData[k]);
                res->AppendAnAssignmentToBuffer(av->GetName(),
                                                av->GetFormulaString(kFormulaStringConversionNormal),
                                                kAppendAnAssignmentToBufferFree | kAppendAnAssignmentToBufferAssignment);
                //true, false, true);
                
            }
        }
    }
    
    if (res->nonempty()) {
        (*res) << '}';
    }
    
    res->TrimSpace();
    return res;
}


//__________________________________________________________________________________

void _CalcNode::ReplaceModel (_String& modelName, _VariableContainer* parent_tree) {
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
bool    _CalcNode::MatchSubtree (_CalcNode* mNode)
{
    node <long>* myNode    = LocateMeInTree (),
                 * matchNode = mNode->LocateMeInTree ();
    if (myNode&&matchNode) {
        return myNode->compare_subtree(matchNode);
    }
    return      false;
}

//_______________________________________________________________________________________________

hyFloat  _CalcNode::ComputeBranchLength (void)
{

    if (theModel < 0) {
        return Value();
    }

    {
      _FString   *stencil = (_FString*)FetchObjectFromVariableByType (&BRANCH_LENGTH_STENCIL,STRING);

      if (stencil && stencil->get_str() == stringSuppliedLengths) {
          return Value();
      }
    }
    {
      _AssociativeList *lookup = (_AssociativeList*)FetchObjectFromVariableByType (&BRANCH_LENGTH_STENCIL,ASSOCIATIVE_LIST);
      if (lookup) {
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

    if (!freqMx || !theMx) {
        return Value();
    }

    hyFloat              weight = 1.0,
                            result = 0.0;

    long                    categoryCounter,
                            totalCategs = 1;

    _CategoryVariable* cVar = nil;

    if (categoryVariables.lLength) {
        for (categoryCounter = 0; categoryCounter<categoryVariables.lLength; categoryCounter++) {
            cVar = (_CategoryVariable*)LocateVar (categoryVariables.lData[categoryCounter]);
            cVar->Refresh();
            totalCategs *= cVar->GetNumberOfIntervals();
        }
    }

    freqMx = (_Matrix*)freqMx->ComputeNumeric();
    categoryCounter = 0;

    do {
        if (categoryVariables.lLength) {
            long c = categoryCounter;
            weight = 1.0;
            for (long k=categoryVariables.lLength-1; k>=0; k--) {
                cVar = (_CategoryVariable*)LocateVar (categoryVariables.lData[k]);
                long t = cVar->GetNumberOfIntervals();
                cVar->SetIntervalValue(c%t);
                weight*=cVar->GetIntervalWeight(c%t);
                c/=t;
            }
        }

        _Matrix*    theMx   = ComputeModelMatrix();
        hyFloat  expSubs = theMx->ExpNumberOfSubs (freqMx, mbf);

        hyFloat divisor;
        checkParameter (largeMatrixBranchLengthDimension, divisor, 20.);

        if (theMx->GetHDim()>divisor) {
            checkParameter (largeMatrixBranchLength, divisor, 3.);
            expSubs /= divisor;
        }

        categoryCounter++;
        result += fabs(expSubs)*weight;
    } while (categoryCounter<totalCategs);

    return result;
}


//_______________________________________________________________________________________________

hyFloat& _CalcNode::operator[] (unsigned long i)
{
    return theProbs [i];
}

//_______________________________________________________________________________________________

BaseRef _CalcNode::toStr (unsigned long) {
    _StringBuffer * res = new _StringBuffer (64L);
    (*res) << theName << '(';

    if (iVariables) {
        _String tempS = (long)(iVariables->lLength/2);
        (*res) << &tempS;
    } else {
        (*res) << '0';
    }
    (*res) << ',';
    if (dVariables) {
        _String tempS = (long)(dVariables->lLength/2);
        (*res) << &tempS;
    } else {
        (*res) << '0';
    }

    (*res) << ')';
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
    lastState     = -1;
    referenceNode = -1;
    slaveNodes    = 0;
}

//_______________________________________________________________________________________________

bool        _CalcNode::HasChanged(bool)
{
    if (_VariableContainer::HasChanged()) {
        return true;
    }
    for (unsigned long i = 0UL; i<categoryVariables.lLength; i++) {
        if (LocateVar (categoryVariables.lData[i])->HasChanged()) {
            return true;
        }
    }
    return false;
}

//_______________________________________________________________________________________________

long        _CalcNode::CheckForReferenceNode(void)
{
    long rN = -1,
         idx = 0;

    // check if all independents are global first
    long modIdx = GetModelIndex();

    if (modIdx != HY_NO_MODEL) {

        if (iVariables && iVariables->lLength) {
            return -1;
        }

        if (dVariables)
            for (idx = 0; idx < dVariables->lLength; idx+=2) {
                if (dVariables->lData[idx+1]>=0) {
                    bool       good = false;
                    _Variable* thisDep = LocateVar (dVariables->lData[idx]);
                    //while (thisDep->NumberOperations() == 1)
                    while (thisDep->varFormula && thisDep->varFormula->NumberOperations () == 1) {
                        _Operation* op = (_Operation*)thisDep->varFormula->GetList() (0);
                        long        isVar = op->GetIndex();
                        if (isVar >= 0) {
                            thisDep = LocateVar (isVar);
                            if (thisDep->IsIndependent()) {
                                good = true;
                                break;
                            }
                        } else {
                            break;
                        }
                    }

                    if (good) {
                        if (thisDep->IsGlobal()) {
                            continue;
                        } else {
                            _String varName = *thisDep->GetName();
                            long    dot = varName.FindBackwards ('.',0,-1);
                            if (dot > 0) {
                                varName.Trim (0,dot-1);
                                dot = LocateVarByName (varName);
                                if (dot < 0) {
                                    break;
                                }

                                if (rN == -1) {
                                    thisDep = FetchVar (dot);

                                    if (thisDep->ObjectClass () != TREE_NODE) {
                                        break;
                                    }

                                    if (((_CalcNode*)thisDep)->GetModelIndex() != modIdx) {
                                        break;
                                    }

                                    rN = thisDep->GetIndex();
                                } else {
                                    if (rN != variableNames.GetXtra(dot)) {
                                        break;
                                    }
                                }
                            } else {
                                break;
                            }
                        }
                    }
                }
            }
    }
    return rN;
}

//_______________________________________________________________________________________________

bool        _CalcNode::NeedNewCategoryExponential(long catID) const
{
    if (isInOptimize && referenceNode>=0) {
        return ((_CalcNode*)LocateVar(referenceNode))->NeedNewCategoryExponential(catID);
    }

    if (_VariableContainer::NeedToExponentiate(catID>=0)) {
        return true;
    }

    if (catID==-1) {
        if (!compExp) {
            return true;
        }

        for (unsigned long i = 0; i<categoryVariables.lLength; i++)
            if (LocateVar (categoryVariables.lData[i])->HasChanged()) {
                return true;
            }
    } else {
        if (!GetCompExp(catID)) {
            return true;
        }

        for (unsigned long i = 0; i<categoryVariables.lLength; i++)
            if (((_CategoryVariable*)LocateVar (categoryVariables.lData[i]))->HaveParametersChanged(remapMyCategories.lData[catID*(categoryVariables.lLength+1)+i+1])) {
                return true;
            }
    }
    return false;
}

//_______________________________________________________________________________________________
bool        _CalcNode::RecomputeMatrix  (long categID, long totalCategs, _Matrix* storeRateMatrix, _List* queue, _SimpleList* tags, _List* bufferedOps)
{
    // assumed that NeedToExponentiate was called prior to this function

    if (isInOptimize) {
      if (referenceNode >= 0) {
        //fprintf (stderr, "\n\n********** REFERENCE NODE ******************\n\n");
        _CalcNode* rN = (_CalcNode*)LocateVar(referenceNode);
        rN->RecomputeMatrix (categID, totalCategs, storeRateMatrix);

        if (totalCategs>1) {
          matrixCache[categID] = rN->matrixCache[categID];
          compExp = matrixCache[categID];
        } else {
          compExp = rN->compExp;
        }

        return false;
      } else {
        if (referenceNode<-1) {
          slaveNodes++;
          if (slaveNodes>1) {
            if (slaveNodes == -referenceNode) {
              slaveNodes = 0;
            }
            return false;
          }
        }
      }
    }


    _Variable* curVar, *locVar;


    if (iVariables)
        for (unsigned long i=0; i<iVariables->lLength; i+=2)
            if (iVariables->lData[i+1]>=0) {
                curVar = LocateVar (iVariables->lData[i+1]);
                if (curVar->IsIndependent()) {
                    locVar = LocateVar (iVariables->lData[i]);
                    curVar->SetValue(locVar->Compute());
                    #ifdef _UBER_VERBOSE_MX_UPDATE_DUMP
                      if (1 || likeFuncEvalCallCount == _UBER_VERBOSE_MX_UPDATE_DUMP_LF_EVAL) {
                        fprintf (stderr, "[_CalcNode::RecomputeMatrix] Node %s, var %s, value = %15.12g\n", GetName()->sData, curVar->GetName()->sData, curVar->Compute()->Value());
                      }
                    #endif
                }
            }

    if (dVariables)
        for (unsigned long i=0; i<dVariables->lLength; i+=2)
            if (dVariables->lData[i+1]>=0) {
                curVar = LocateVar (dVariables->lData[i+1]);
                if (curVar->IsIndependent()) {
                    locVar = LocateVar (dVariables->lData[i]);
                    curVar->SetValue(locVar->Compute());
                    #ifdef _UBER_VERBOSE_MX_UPDATE_DUMP
                      if (1 || likeFuncEvalCallCount == _UBER_VERBOSE_MX_UPDATE_DUMP_LF_EVAL) {
                        fprintf (stderr, "[_CalcNode::RecomputeMatrix] Node %s, var %s, value = %15.12g\n", GetName()->sData, curVar->GetName()->sData, curVar->Compute()->Value());
                      }
                    #endif
                }
            }

    #ifdef _UBER_VERBOSE_MX_UPDATE_DUMP
      if (1|| likeFuncEvalCallCount == _UBER_VERBOSE_MX_UPDATE_DUMP_LF_EVAL && gVariables) {
        for (unsigned long i=0; i<gVariables->lLength; i++) {
          _Variable* curVar = LocateVar(gVariables->GetElement(i));
          fprintf (stderr, "[_CalcNode::RecomputeMatrix] Node %s, var %s, value = %15.12g\n", GetName()->sData, curVar->GetName()->sData, curVar->Compute()->Value());
        }
      }
    #endif


    for (unsigned long i=0; i<categoryVariables.lLength; i++) {
        if (categoryIndexVars.lData[i]<0) {
            continue;
        }
        curVar = LocateVar (categoryIndexVars.lData[i]);
        locVar = LocateVar (categoryVariables.lData[i]);
        curVar->SetValue(locVar->Compute());
    }

    if (!storeRateMatrix) {
      if (totalCategs>1) {
  #ifdef _UBER_VERBOSE_MX_UPDATE_DUMP
        fprintf (stderr, "[_CalcNode::RecomputeMatrix] Deleting category %ld for node %s at %p\n", categID, GetName()->sData, GetCompExp(categID));
  #endif
        DeleteObject(GetCompExp(categID, true));

      } else {
        if (compExp) {
          DeleteObject (compExp);
          compExp = nil;
        }
      }
    }

    bool    isExplicitForm  = HasExplicitFormModel ();

    if (isExplicitForm && bufferedOps) {
        _Matrix * bufferedExp = (_Matrix*)GetExplicitFormModel()->Compute (0,nil, bufferedOps);
        #ifdef _UBER_VERBOSE_MX_UPDATE_DUMP
            fprintf (stderr, "[_CalcNode::RecomputeMatrix] Setting (buffered) category %ld/%ld for node %s\n", categID, totalCategs, GetName()->sData);
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
                temp = (_Matrix*)myModelMatrix->MultByFreqs(theModel);
            }

            if (dVariables)
                for (unsigned long i=0; i<dVariables->lLength; i+=2)
                    if (dVariables->lData[i+1]>=0) {
                        curVar = LocateVar (dVariables->lData[i+1]);
                        if (!curVar->IsIndependent()) {
                            locVar = LocateVar (dVariables->lData[i]);
                            if (locVar->IsIndependent()) {
                                locVar->SetValue (curVar->Compute());
                            }
                        }
                    }

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
                fprintf (stderr, "[_CalcNode::RecomputeMatrix] Setting category %ld/%ld for node %s\n", categID, totalCategs, GetName()->sData);
            #endif
            SetCompExp ((_Matrix*)(isExplicitForm?temp:temp->Exponentiate()), totalCategs>1?categID:-1);

        } else {
            compExp = (_Matrix*)myModelMatrix->Evaluate(false);
        }
    }
    return false;
}

//_______________________________________________________________________________________________
void        _CalcNode::SetCompExp  (_Matrix* m, long catID)
{
    compExp = m;
    if (catID >= 0 && matrixCache) {
        if (remapMyCategories.lLength) {
            catID = remapMyCategories.lData[catID*(categoryVariables.lLength+1)];
        }
        matrixCache[catID] = compExp;
    }
}
//_______________________________________________________________________________________________
_Matrix*        _CalcNode::ComputeModelMatrix  (bool)
{
    // assumed that NeedToExponentiate was called prior to this function
    _Variable   * curVar,
                *locVar;

    if (iVariables)
        for (unsigned long i=0; i<iVariables->lLength; i+=2)
            if (iVariables->lData[i+1]>=0) {
                curVar = LocateVar (iVariables->lData[i+1]);
                if (curVar->IsIndependent()) {
                    locVar = LocateVar (iVariables->lData[i]);
                    curVar->SetValue(locVar->Compute());
                }
            }

    if (dVariables)
        for (unsigned long i=0; i<dVariables->lLength; i+=2)
            if (dVariables->lData[i+1]>=0) {
                curVar = LocateVar (dVariables->lData[i+1]);
                if (curVar->IsIndependent()) {
                    locVar = LocateVar (dVariables->lData[i]);
                    curVar->SetValue(locVar->Compute());
                }
            }

    _Matrix * modelMx = GetModelMatrix();
    if (modelMx && modelMx->ObjectClass()==MATRIX && modelMx->MatrixType()!=_POLYNOMIAL_TYPE) {
        return (_Matrix*)modelMx->ComputeNumeric();
    }

    return nil;
}

//_______________________________________________________________________________________________

_Matrix*    _CalcNode::GetCompExp       (long catID, bool doClear) const
{
    if (catID==-1) {
        return compExp;
    } else {

        if (remapMyCategories.lLength) {
            catID = remapMyCategories.lData[catID * (categoryVariables.lLength+1)];
        }

        _Matrix* ret = matrixCache?matrixCache[catID]:compExp;

        if (doClear && matrixCache) {
            matrixCache[catID] = nil;
        }
        return ret;
    }
}

//_______________________________________________________________________________________________

BaseRef     _CalcNode::makeDynamic(void) const
{
    _CalcNode* res = new (_CalcNode);
    res->_VariableContainer::Duplicate (this);
    res->categoryVariables.Duplicate ((BaseRef)&categoryVariables);
    //res->randomVariables.Duplicate ((BaseRef)&randomVariables);
    res->categoryIndexVars.Duplicate ((BaseRef)&categoryIndexVars);
    //res->randomIndexVars.Duplicate ((BaseRef)&randomIndexVars);
    res->theValue = theValue;
    res->cBase = cBase;
    if (cBase) {
        res->theProbs = new hyFloat [cBase];
        memcpy (res->theProbs, theProbs, sizeof(hyFloat)*cBase);
    } else {
        res->theProbs = nil;
    }
    res->compExp = compExp;
    if (compExp) {
        compExp->AddAReference();
    }
    res->referenceNode = referenceNode;
    res->slaveNodes    = slaveNodes;
    return res;
}

//_______________________________________________________________________________________________

void     _CalcNode::SetupCategoryMap (_List& containerVariables, _SimpleList& classCounter, _SimpleList& multipliers)
{
    
    long    totalCategories = classCounter.Element(-1),
    globalCatCount  = containerVariables.lLength-1,
    localCategories = 1,
    catCount        = categoryVariables.lLength-1,
    entriesPerCat   = 2+catCount;
    
    //for (long k = 0; k<categoryVariables.lLength;k++)
    //    printf ("%ld\n", categoryVariables(k));//, ((_Variable*)categoryVariables(k))->GetName()->sData);
    
    if (catCount<0) {
        remapMyCategories.Clear();
    } else {
        
        remapMyCategories.Populate (totalCategories*entriesPerCat,0,0);
        
        _SimpleList     remappedIDs,
        rateMultiplers (categoryVariables.lLength,1,0),
        categoryValues (globalCatCount+1,0,0);
        
        for (long myCatID = 0; myCatID <= catCount; myCatID++) {
            long coordinate = containerVariables.FindPointer(LocateVar(categoryVariables.lData[myCatID]));
            if (coordinate < 0) {
                hy_global::ReportWarning ("Internal error in SetupCategoryMap. Please report to sergeilkp@icloud.com");
            }
            localCategories *= classCounter.lData[coordinate];
            remappedIDs << coordinate;
        }
        
        for (long myCatID = catCount-1; myCatID >= 0; myCatID--) {
            rateMultiplers.lData[myCatID] = rateMultiplers.lData[myCatID+1]*classCounter.lData[remappedIDs.lData[myCatID+1]];
        }
        
        for (long currentRateCombo  = 0; currentRateCombo < totalCategories; currentRateCombo++) {
            long copyRateCombo = currentRateCombo;
            for (long variableID = 0; variableID <= globalCatCount; variableID++) {
                categoryValues.lData[variableID] = copyRateCombo / multipliers.lData[variableID];
                copyRateCombo = copyRateCombo%multipliers.lData[variableID];
                //printf ("%d %d %d %d\n", currentRateCombo, variableID, multipliers.lData[variableID], categoryValues.lData[variableID]);
            }
            
            long localCatID = 0;
            
            for  (long localVariableID = 0; localVariableID<=catCount; localVariableID++) {
                localCatID += rateMultiplers.lData[localVariableID] * categoryValues.lData[remappedIDs.lData[localVariableID]];
            }
            
            long offset = currentRateCombo * entriesPerCat;
            remapMyCategories.lData[offset] = localCatID;
            //printf ("[%ld] = %ld (%ld)\n", offset, localCatID, );
            
            offset++;
            for  (long localVariableID = 0; localVariableID<=catCount; localVariableID++) {
                remapMyCategories[offset++] = categoryValues.lData[remappedIDs.lData[localVariableID]];
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
        _Matrix * mm = GetModelMatrix();
        if (mm) {
            mm->MakeMeSimple();
        }
        
        mm = GetFreqMatrix();
        if (mm) {
            mm->MakeMeSimple();
        }
    }
}

//_______________________________________________________________________________________________

void _CalcNode::ConvertFromSimpleMatrix (void) {
    _Formula * mf = GetExplicitFormModel();
    if (mf) {
        mf->ConvertMatrixArgumentsToSimpleOrComplexForm (true);
    } else {
        _Matrix * mm = GetModelMatrix();
        if (mm) {
            mm->MakeMeGeneral();
        }
        
        mm = GetFreqMatrix();
        
        if (mm) {
            mm->MakeMeGeneral();
        }
    }
}

//_______________________________________________________________________________________________

_Formula*   _CalcNode::RecurseMC (long varToConstrain, node<long>* whereAmI, bool first, char rooted) {
    long descendants = whereAmI->get_num_nodes(),
    f = iVariables?iVariables->FindStepping(varToConstrain,2,1):-1,
    k,
    l,
    start = 0;
    
    if (f<0 && !first) {
        HandleApplicationError (_String ("Molecular clock constraint has failed, since variable '")
                                &*LocateVar(varToConstrain)->GetName()
                                &"' is not an independent member of the node '"
                                & *GetName()
                                & '\''
                                );
        return nil;
    }
    
    
    if (descendants == 0) {
        if (first) {
            return nil;
        } else {
            return new _Formula (LocateVar(iVariables->lData[f-1]),true);
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
    
    for (k=start+1; k<=descendants; k++) {
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
    
    for (k=0; k<descendants-start; k++)
        if ((nodeConditions[k])->GetList().lLength>1) {
            break;
        }
    
    if (k==descendants-start) { // all underlying branches are "simple"
        for (k=1; k<descendants-start; k++) {
            LocateVar (((_Operation*)((*(nodeConditions[k])).GetList()(0)))->GetIndex())->SetFormula (*nodeConditions[0]);
            delete (nodeConditions[k]);
            nodeConditions[k] = nil;
        }
        k = 0;
    } else {
        
        for (l=k+1; l<descendants-start; l++)
            if (nodeConditions[l]->GetList().lLength>1) {
                break;
            }
        
        if (l==descendants-start) // all but one underlying branches are "simple"
            for (l=0; l<descendants-start; l++) {
                if (l==k) {
                    continue;
                }
                LocateVar (nodeConditions[l]->GetIthTerm(0)->GetIndex())->SetFormula (*nodeConditions[k]);
                delete (nodeConditions[l]);
                nodeConditions[l] = nil;
            }
        // really bad bongos! must solve for non-additive constraint
        else
            for (l=0; l<descendants-start; l++) {
                if (l==k) {
                    continue;
                }
                if (nodeConditions[l]->GetList().lLength==1) {
                    LocateVar (nodeConditions[l]->GetIthTerm(0)->GetIndex())->SetFormula (*nodeConditions[k]);
                } else { // solve for a non-additive constraint
                    _Variable* nonAdd = LocateVar (nodeConditions[l]->GetIthTerm(0)->GetIndex());
                    nodeConditions[l]->GetList().Delete(0);
                    _Formula  newConstraint;
                    newConstraint.Duplicate(nodeConditions[k]);
                    for (long m=0; m<nodeConditions[l]->GetList().lLength; m++) {
                        _Operation* curOp = (_Operation*)(*nodeConditions[l]).GetList()(m);
                        if (curOp->GetNoTerms()) {
                            newConstraint.GetList().AppendNewInstance(new _Operation ('-', 2L));
                        } else {
                            newConstraint.GetList()<<curOp;
                        }
                    }
                    delete (nodeConditions[l]);
                    nodeConditions[l] = nil;
                    nonAdd->SetFormula(newConstraint);
                }
            }
    }
    
    
    if(!first) {
        _Formula     *result = nodeConditions[k];
        
        _Operation   *newVar = new _Operation;
        
        newVar->SetAVariable (iVariables->lData[f-1]);
        
        result->GetList().AppendNewInstance( newVar);
        result->GetList().AppendNewInstance(new _Operation ('+', 2L));
        
        
        delete [] nodeConditions;
        return result;
    }
    
    for (k=0; k<descendants-start; k++)
        if (nodeConditions[k]) {
            delete nodeConditions[k];
        }
    
    delete [] nodeConditions;
    return nil;
}





