/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2009
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon              (apoon@cfenet.ubc.ca)

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

#include "defines.h"

#include "math.h"
#include "ctype.h"
#include "string.h"
#include "calcnode.h"
#include "calcnode_globals.h"
#include "thetree.h"
#include "scfg.h"
#include "legacy_parser.h"

#include "category.h"
#include "batchlan.h"
#include "likefunc.h"
#include "float.h"

#include "hy_globals.h"

extern _Parameter explicitFormMatrixExponential;
extern _String VerbosityLevelString, BRANCH_LENGTH_STENCIL;

#ifdef __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

//______________________________________________________________________________
_CalcNode::_CalcNode() {
  theProbs = nil; // default constructor, doesn't do much
  compExp = nil;
  referenceNode = -1;
}

//______________________________________________________________________________
_CalcNode::_CalcNode(_String name, _String parms, int codeBase,
                     _VariableContainer *theP, _AVLListXL *aCache)
    : _VariableContainer(name, "", theP)
      // construct a node from a string of the form
      // matrix name, <optional comma separated variable declarations,
      // inititalizations>
      // also should be passed the pointer to a container tree
      {
  InitializeCN(parms, codeBase, theP, aCache);
}

//______________________________________________________________________________
_CalcNode::_CalcNode(_CalcNode *sourceNode, _VariableContainer *theP)
    : _VariableContainer(sourceNode->ContextFreeName(), "", theP) {
  _String model = sourceNode->GetModelName();
  InitializeCN(model, 0, theP);
  if (model.sLength) { // copy model parameter values
    CopyMatrixParameters(sourceNode, true);
  }
}

//______________________________________________________________________________
void _CalcNode::InitializeCN(_String &parms, int, _VariableContainer *theP,
                             _AVLListXL *aCache) {
  cBase = 0;
  theProbs = nil;
  compExp = nil;
  referenceNode = -1;
  slaveNodes = 0;

  long f = parms.Find(','), g = -1;

  _String matrixName(parms, 0, f >= 0 ? f - 1 : -1);

  InitializeVarCont(empty, matrixName, theP, aCache);

  if (GetModelIndex() == HY_NO_MODEL && parms.Length()) {
    f = 0;
  }

  while (f != -1) {
    g = parms.Find(',', f + 1, -1);
    if (f == 0) {
      f = -1;
    }

    if (g != -1) {
      _String paramName(parms, f + 1, g - 1);
      _Formula fg(paramName, this);
    } else {
      _String paramName(parms, f + 1, g);
      _Formula fg(paramName, this);
    }
    f = g;
  }

  // attach the variables to the lists inside the node
  ScanAndAttachVariables();

  // check for category variables
  if (iVariables) {
    for (f = iVariables->lLength - 2; f >= 0 && iVariables->lData[f + 1] >= 0;
         f -= 2) {
      _Variable *theV = LocateVar(iVariables->lData[f + 1]);
      if (theV->IsCategory()) {
        /* this has to do with local category variables;
           NOT TESTED and could be BROKEN
        */
        _CategoryVariable *theCV = (_CategoryVariable *)theV;

        _Formula newDensity, newCumulative;

        _SimpleList iv, iv2, dv, dv2;

        for (unsigned long k = 0; k < iVariables->lLength; k += 2) {
          iv << iVariables->lData[k];
          iv2 << iVariables->lData[k + 1];
        }

        if (dVariables)
          for (unsigned long k = 0; k < dVariables->lLength; k += 2) {
            dv << dVariables->lData[k];
            dv2 << dVariables->lData[k + 1];
          }

        newDensity.LocalizeFormula(theCV->GetDensity(), *GetName(), iv, iv2, dv, dv2);
        newCumulative.LocalizeFormula(theCV->GetCumulative(), *GetName(), iv, iv2, dv, dv2);

        _CategoryVariable newCV;
        newCV.Duplicate(theCV);
        newCV.GetDensity().Duplicate((BaseRef) & newDensity);
        newCV.GetCumulative().Duplicate((BaseRef) & newCumulative);

        theV = LocateVar(iVariables->lData[f]);
        newCV.GetName()->Duplicate(theV->GetName());
        ReplaceVar(&newCV);

        categoryVariables << iVariables->lData[f];
        categoryIndexVars << iVariables->lData[f + 1];
        iVariables->Delete(f);
        iVariables->Delete(f);
      }
    }

    if (iVariables->lLength) {
      iVariables->TrimMemory();
    } else {
      delete (iVariables);
      iVariables = nil;
    }
  }
  if (gVariables) {
    for (f = gVariables->lLength - 1; f >= 0; f--) {
      _Variable *theV = LocateVar(gVariables->lData[f]);
      if (theV->IsCategory()) {
        categoryVariables << gVariables->lData[f];
        categoryIndexVars << -1;
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

  BaseRef temp = (variablePtrs(theIndex));
  variablePtrs[theIndex] = this->makeDynamic();
  DeleteObject(temp);

}

//______________________________________________________________________________
void _CalcNode::SetModel(long modelID, _AVLListXL *varCache) {
  _VariableContainer::SetModel(modelID, varCache);
}

//______________________________________________________________________________
long _CalcNode::SetDependance(long varIndex) {
  varIndex = _VariableContainer::SetDependance(varIndex);
  if (varIndex >= 0) {
    _SimpleList checkVars;
    _AVLList myVars(&checkVars);
    LocateVar(varIndex)->ScanForVariables(myVars, true);

    for (unsigned long k = 0; k < checkVars.lLength; k++)
      if (LocateVar(checkVars.lData[k])->IsCategory() &&
          (categoryVariables >> checkVars.lData[k])) {
        categoryIndexVars << -1;
      }

  }
  return varIndex;
}

//______________________________________________________________________________
void _CalcNode::SetCodeBase(int codeBase) {
  if (codeBase > 0) {
    if ((codeBase != cBase) || !theProbs) {
#ifndef __HYALTIVEC__
      if (theProbs) {
        delete theProbs;
      }
      theProbs = new _Parameter[codeBase];
#else
      if (theProbs) {
        vec_free(theProbs);
      }
      theProbs = (_Parameter *)VecMemAllocate(codeBase * sizeof(_Parameter));
#endif
      cBase = codeBase;
      theProbs[0] = 1.0;
    } else {
      theProbs[0] = 1.0;
    }
  }
}

//______________________________________________________________________________
void _CalcNode::SetCompMatrix(long categID) { compExp = GetCompExp(categID); }

//______________________________________________________________________________
_CalcNode::~_CalcNode(void) {

#ifndef __HYALTIVEC__
  if (theProbs) {
    delete[] theProbs;
  }
#else
  if (theProbs) {
    vec_free(theProbs);
  }
#endif
  if (compExp && referenceNode < 0) {
    DeleteObject(compExp);
  }
}

//______________________________________________________________________________
long _CalcNode::FreeUpMemory(long) {
  long res = 0;
  if (compExp && referenceNode < 0) {
    res = compExp->GetMySize();
    DeleteObject(compExp);
    compExp = nil;
  }
  return res;
}

//______________________________________________________________________________
void _CalcNode::RemoveModel(void) {
  Clear();
  if (compExp && referenceNode < 0) {
    DeleteObject(compExp);
    compExp = nil;
  }
}

//______________________________________________________________________________
void _CalcNode::ReplaceModel(_String &modelName,
                             _VariableContainer *theConext) {
  RemoveModel();
  DeleteVariable(theIndex,
                 false); // this will clean up all the node.xx variables
  InitializeCN(modelName, 0, theConext);
}

//______________________________________________________________________________
bool _CalcNode::MatchSubtree(_CalcNode *mNode) {
  node<long> *myNode = LocateMeInTree(), *matchNode = mNode->LocateMeInTree();
  if (myNode && matchNode) {
    return myNode->compare_subtree(matchNode);
  }
  return false;
}

//______________________________________________________________________________
_Parameter _CalcNode::BranchLength(void) {

  if (theModel < 0) {
    return Value();
  }

  _FString *stencil =
      (_FString *)FetchObjectFromVariableByType(&BRANCH_LENGTH_STENCIL, STRING);

  if (stencil && stencil->theString->Equal(&stringSuppliedLengths)) {
    return Value();
  }

  _Matrix *freqMx, *theMx;

  bool mbf;

  RetrieveModelComponents(theModel, theMx, freqMx, mbf);

  if (!freqMx && !theModel) {
    return Value();
  }

  _Parameter weight = 1.0, result = 0.0;

  unsigned long categoryCounter, totalCategs = 1;

  _CategoryVariable *cVar = nil;

  if (categoryVariables.lLength) {
    for (categoryCounter = 0; categoryCounter < categoryVariables.lLength; categoryCounter++) {
      cVar = (_CategoryVariable *)LocateVar(categoryVariables.lData[categoryCounter]);
      cVar->Refresh();
      totalCategs *= cVar->GetNumberOfIntervals();
    }
  }

  freqMx = dynamic_cast<_Matrix*> (freqMx->ComputeNumeric());
  categoryCounter = 0;

  do {
    if (categoryVariables.lLength) {
      long c = categoryCounter;
      weight = 1.0;
      for (long k = categoryVariables.lLength - 1; k >= 0; k--) {
        cVar = (_CategoryVariable *)LocateVar(categoryVariables.lData[k]);
        long t = cVar->GetNumberOfIntervals();
        cVar->SetIntervalValue(c % t);
        weight *= cVar->GetIntervalWeight(c % t);
        c /= t;
      }
    }

    _Matrix *theMx = ComputeModelMatrix();
    _Parameter expSubs = theMx->ExpNumberOfSubs(freqMx, mbf);

    _Parameter divisor;
    checkParameter(largeMatrixBranchLengthDimension, divisor, 20.);

    if (theMx->GetHDim() > divisor) {
      checkParameter(largeMatrixBranchLength, divisor, 3.);
      expSubs /= divisor;
    }

    categoryCounter++;
    result += fabs(expSubs) * weight;
  } while (categoryCounter < totalCategs);

  return result;
}

//______________________________________________________________________________
_Parameter &_CalcNode::operator[](unsigned long i) { return theProbs[i]; }

//______________________________________________________________________________
BaseRef _CalcNode::toStr(void) {
  _String *res = new _String(16L, true);
  checkPointer(res);
  (*res) << theName;
  (*res) << '(';

  if (iVariables) {
    _String tempS = (long)(iVariables->lLength / 2);
    (*res) << &tempS;
  } else {
    (*res) << '0';
  }
  (*res) << ',';
  if (dVariables) {
    _String tempS = (long)(dVariables->lLength / 2);
    (*res) << &tempS;
  } else {
    (*res) << '0';
  }

  (*res) << ')';
  res->Finalize();
  return res;
}

//______________________________________________________________________________
void _CalcNode::Duplicate(BaseRef theO) {
  _VariableContainer::Duplicate(theO);
  cBase = 0;
  compExp = nil;
  matrixCache = nil;
  theProbs = nil;
  lastState = -1;
  referenceNode = -1;
  slaveNodes = 0;
}

//______________________________________________________________________________
bool _CalcNode::HasChanged(void) {
  if (_VariableContainer::HasChanged()) {
    return true;
  }
  for (unsigned long i = 0; i < categoryVariables.lLength; i++)
    if (LocateVar(categoryVariables.lData[i])->HasChanged()) {
      return true;
    }
  return false;
}

//______________________________________________________________________________
long _CalcNode::CheckForReferenceNode(void) {
  long rN = -1, idx = 0;

  // check if all independents are global first
  long modIdx = GetModelIndex();

  if (modIdx != HY_NO_MODEL) {

    if (iVariables && iVariables->lLength) {
      return HY_NOT_FOUND;
    }

    if (dVariables)
      for (idx = 0; idx < dVariables->lLength; idx += 2) {
        if (dVariables->lData[idx + 1] >= 0) {
          bool good = false;
          _Variable *thisDep = LocateVar(dVariables->lData[idx]);
          //while (thisDep->NumberOperations() == 1)
          while (thisDep->varFormula && thisDep->varFormula->NumberOperations() == 1) {
            _Operation *op = thisDep->varFormula->GetIthTerm(0);
            long isVar = op->GetAVariable();
            if (isVar >= 0) {
              thisDep = LocateVar(isVar);
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
              long dot = varName.FindBackwards('.', 0, -1);
              if (dot > 0) {
                varName.Trim(0, dot - 1);
                dot = LocateVarByName(varName);
                if (dot < 0) {
                  break;
                }

                if (rN == -1) {
                  thisDep = FetchVar(dot);

                  if (thisDep->ObjectClass() != TREE_NODE) {
                    break;
                  }

                  if (((_CalcNode *)thisDep)->GetModelIndex() != modIdx) {
                    break;
                  }

                  rN = thisDep->GetAVariable();
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

//______________________________________________________________________________
bool _CalcNode::NeedToExponentiate(long catID) {
  if (isInOptimize && referenceNode >= 0) {
    return ((_CalcNode *)LocateVar(referenceNode))->NeedToExponentiate(catID);
  }

  if (_VariableContainer::NeedToExponentiate(catID >= 0)) {
    return true;
  }

  if (catID == -1) {
    if (!compExp) {
      return true;
    }

    for (unsigned long i = 0; i < categoryVariables.lLength; i++)
      if (LocateVar(categoryVariables.lData[i])->HasChanged()) {
        return true;
      }
  } else {
    if (!GetCompExp(catID)) {
      return true;
    }

    for (unsigned long i = 0; i < categoryVariables.lLength; i++)
      if (((_CategoryVariable *)LocateVar(categoryVariables.lData[i]))->HaveParametersChanged(remapMyCategories.lData[catID * (categoryVariables.lLength + 1) + i + 1])) {
        return true;
      }
  }
  return false;
}

//______________________________________________________________________________
bool _CalcNode::RecomputeMatrix(long categID, long totalCategs,
                                _Matrix *storeRateMatrix, _List *queue,
                                _SimpleList *tags, _List *bufferedOps) {
  // assumed that NeedToExponentiate was called prior to this function

  if (isInOptimize)
  {
    if (referenceNode >= 0) {
      _CalcNode *rN = (_CalcNode *)LocateVar(referenceNode);
      rN->RecomputeMatrix(categID, totalCategs, storeRateMatrix);

      if (totalCategs > 1) {
        matrixCache[categID] = rN->matrixCache[categID];
        compExp = matrixCache[categID];
      } else {
        compExp = rN->compExp;
      }
      return false;
    } else {
      if (referenceNode < -1) {
        slaveNodes++;
        if (slaveNodes > 1) {
          if (slaveNodes == -referenceNode) {
            slaveNodes = 0;
          }
          return false;
        }
      }
    }
  }

  _Variable *curVar, *locVar;

  if (iVariables)
    for (unsigned long i = 0; i < iVariables->lLength; i += 2)
      if (iVariables->lData[i + 1] >= 0) {
        curVar = LocateVar(iVariables->lData[i + 1]);
        if (curVar->IsIndependent()) {
          locVar = LocateVar(iVariables->lData[i]);
          curVar->SetValue(locVar->Compute());
        }
      }

  if (dVariables)
    for (unsigned long i = 0; i < dVariables->lLength; i += 2)
      if (dVariables->lData[i + 1] >= 0) {
        curVar = LocateVar(dVariables->lData[i + 1]);
        if (curVar->IsIndependent()) {
          locVar = LocateVar(dVariables->lData[i]);
          curVar->SetValue(locVar->Compute());
        }
      }

  for (unsigned long i = 0; i < categoryVariables.lLength; i++) {
    if (categoryIndexVars.lData[i] < 0) {
      continue;
    }
    curVar = LocateVar(categoryIndexVars.lData[i]);
    locVar = LocateVar(categoryVariables.lData[i]);
    curVar->SetValue(locVar->Compute());
  }

  if (!storeRateMatrix) {
    if (totalCategs > 1) {
      DeleteObject(GetCompExp(categID));
    } else if (compExp) {
      DeleteObject(compExp);
      compExp = nil;
    }
  }

  bool isExplicitForm = HasExplicitFormModel();

  if (isExplicitForm && bufferedOps) {
    _Matrix *bufferedExp = dynamic_cast<_Matrix*> (GetExplicitFormModel()->Compute(0, _hyDefaultExecutionContext, bufferedOps));
    SetCompExp(dynamic_cast <_Matrix*> (bufferedExp->makeDynamic()),totalCategs > 1 ? categID : -1);
    return false;
  }

  unsigned long previous_length = queue && tags ? queue->lLength : 0;
  _Matrix *myModelMatrix = GetModelMatrix(queue, tags);

  if (isExplicitForm && !myModelMatrix) { // matrix exponentiations got cached
    if (queue->lLength > previous_length) {
      return true;
    } else {
      WarnError("Internal error");
      return false;
    }
  } else {

    if (myModelMatrix->MatrixType() != _HY_MATRIX_POLYNOMIAL_TYPE) {
      _Matrix *temp = nil;
      if (isExplicitForm) {
        temp = dynamic_cast <_Matrix*>(myModelMatrix->makeDynamic());
      } else {
        temp = dynamic_cast <_Matrix*>(myModelMatrix->MultByFreqs(theModel));
      }

      if (dVariables)
        for (unsigned long i = 0; i < dVariables->lLength; i += 2)
          if (dVariables->lData[i + 1] >= 0) {
            curVar = LocateVar(dVariables->lData[i + 1]);
            if (!curVar->IsIndependent()) {
              locVar = LocateVar(dVariables->lData[i]);
              if (locVar->IsIndependent()) {
                locVar->SetValue(curVar->Compute());
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

      SetCompExp((_Matrix *)(isExplicitForm ? temp : temp->Exponentiate()),
                 totalCategs > 1 ? categID : -1);

    } else {
      compExp = dynamic_cast<_Matrix*>(myModelMatrix->Evaluate(false));
    }
  }
  return false;
}

//______________________________________________________________________________
void _CalcNode::SetCompExp(_Matrix *m, long catID) {
  compExp = m;
  if (catID >= 0 && matrixCache) {
    if (remapMyCategories.lLength) {
      catID = remapMyCategories.lData[catID * (categoryVariables.lLength + 1)];
    }
    matrixCache[catID] = compExp;
  }
}
//______________________________________________________________________________
_Matrix *_CalcNode::ComputeModelMatrix(bool) {
  // assumed that NeedToExponentiate was called prior to this function
  _Variable *curVar, *locVar;

  if (iVariables)
    for (unsigned long i = 0; i < iVariables->lLength; i += 2)
      if (iVariables->lData[i + 1] >= 0) {
        curVar = LocateVar(iVariables->lData[i + 1]);
        if (curVar->IsIndependent()) {
          locVar = LocateVar(iVariables->lData[i]);
          curVar->SetValue(locVar->Compute());
        }
      }

  if (dVariables)
    for (unsigned long i = 0; i < dVariables->lLength; i += 2)
      if (dVariables->lData[i + 1] >= 0) {
        curVar = LocateVar(dVariables->lData[i + 1]);
        if (curVar->IsIndependent()) {
          locVar = LocateVar(dVariables->lData[i]);
          curVar->SetValue(locVar->Compute());
        }
      }

  _Matrix *modelMx = GetModelMatrix();
  if (modelMx && modelMx->ObjectClass() == MATRIX &&
      modelMx->MatrixType() != _HY_MATRIX_POLYNOMIAL_TYPE) {
    return dynamic_cast<_Matrix*>(modelMx->ComputeNumeric());
  }

  return nil;
}

//______________________________________________________________________________
_Matrix *_CalcNode::GetCompExp(long catID) {
  if (catID == -1) {
    return compExp;
  } else {
    if (remapMyCategories.lLength) {
      catID = remapMyCategories.lData[catID * (categoryVariables.lLength + 1)];
    }

    return matrixCache ? matrixCache[catID] : compExp;
  }
}

//______________________________________________________________________________
_CalcNode::_CalcNode (_CalcNode& sourceNode) {
    _VariableContainer::Duplicate(&sourceNode);
    categoryVariables.Duplicate((BaseRef) & sourceNode.categoryVariables);
    //res->randomVariables.Duplicate ((BaseRef)&randomVariables);
    categoryIndexVars.Duplicate((BaseRef) & sourceNode.categoryIndexVars);
    //res->randomIndexVars.Duplicate ((BaseRef)&randomIndexVars);
    theValue = sourceNode.theValue;
    cBase = sourceNode.cBase;
    if (cBase) {
      theProbs = new _Parameter[cBase];
      memcpy(theProbs, sourceNode.theProbs, sizeof(_Parameter) * cBase);
    } else {
      theProbs = nil;
    }
    compExp = sourceNode.compExp;
    if (compExp) {
      compExp->AddAReference();
    }
    referenceNode = sourceNode.referenceNode;
    slaveNodes = sourceNode.slaveNodes;
}

//______________________________________________________________________________
BaseRef _CalcNode::makeDynamic(void) {
    return new _CalcNode (*this);
}

//______________________________________________________________________________
_PMathObj _CalcNode::Compute(void) { return this; }

#ifdef __MP__
//______________________________________________________________________________
void *MatrixUpdateFunction(void *arg) {
  ThreadMatrixTask *theTask = (ThreadMatrixTask *)arg;
  for (long k = theTask->startAt; k < theTask->endAt; k++) {
    ((_CalcNode *)(theTask->updateCN->lData[k]))
        ->RecomputeMatrix(theTask->cID, theTask->tcat);
  }
  return nil;
}
#endif

//______________________________________________________________________________
node<long> *_CalcNode::LocateMeInTree(void) {

  _String baseNode = GetName()->Cut(0, GetName()->Find('.') - 1);
  _TheTree *parentTree = (_TheTree *)FetchVar(LocateVarByName(baseNode));
  _CalcNode *curNode = parentTree->StepWiseTraversal(true);

  baseNode = GetName()->Cut(GetName()->FindBackwards('.', 0, -1), -1);
  while (curNode) {
    if (curNode->GetName()->endswith(baseNode)) {
      return &parentTree->GetCurrentNode();
    }
    curNode = parentTree->StepWiseTraversal();
  }
  return nil;
}

//______________________________________________________________________________
long _CalcNode::ConvertToSimpleMatrix(void) {
  _Formula *mf = GetExplicitFormModel();
  if (mf) {
    mf->ConvertMatrixArgumentsToSimpleOrComplexForm(false);
  } else {
    _Matrix *mm = GetModelMatrix();
    if (mm) {
      mm->MakeMeSimple();
    }

    mm = GetFreqMatrix();
    if (mm) {
      mm->MakeMeSimple();
    }
  }

  return 0;
}

//______________________________________________________________________________
void _CalcNode::ConvertFromSimpleMatrix(void) {
  _Formula *mf = GetExplicitFormModel();
  if (mf) {
    mf->ConvertMatrixArgumentsToSimpleOrComplexForm(true);
  } else {
    _Matrix *mm = GetModelMatrix();
    if (mm) {
      mm->MakeMeGeneral();
    }

    mm = GetFreqMatrix();

    if (mm) {
      mm->MakeMeGeneral();
    }
  }
}

//______________________________________________________________________________
_Formula *_CalcNode::RecurseMC(long varToConstrain, node<long> *whereAmI,
                               bool first, char rooted) {
  long descendants = whereAmI->get_num_nodes(),
       f = iVariables ? iVariables->FindStepping(varToConstrain, 2, 1) : -1, k,
       l, start = 0;

  if ((f < 0) && (!first)) { // bad bongos!
    _String errMsg("Molecular clock constraint has failed, since variable ");
    errMsg = errMsg & *LocateVar(varToConstrain)->GetName();
    errMsg = errMsg & " is not an independent member of the node ";
    errMsg = errMsg & *GetName();
    WarnError(errMsg);
    return nil;
  }

  if (descendants == 0) { // leaf node
    if (first) {
      return nil;
    } else {
      return new _Formula(LocateVar(iVariables->lData[f - 1]), true);
    }
  }

  if (first && (!whereAmI->get_parent()) && (rooted == ROOTED_LEFT)) {
    descendants--;
  }
  if (first && (!whereAmI->get_parent()) && (rooted == ROOTED_RIGHT)) {
    start++;
  }

  // internal node - must do some work

  _Formula **nodeConditions =
      (_Formula **)MemAllocate((descendants - start) * sizeof(_Formula *));

  for (k = start + 1; k <= descendants; k++) {
    node<long> *downWeGo = whereAmI->go_down(k);
    if (!(nodeConditions[k - 1 - start] = ((_CalcNode *)LocateVar(
            downWeGo->get_data()))->RecurseMC(varToConstrain, downWeGo))) {
      for (long f2 = 0; f2 < k - start - 1; f2++) {
        delete nodeConditions[f2];
      }

      free(nodeConditions);
      return nil;
    }
  }

  // all the conditions have been written. now check how we should resolve them

  for (k = 0; k < descendants - start; k++)
    if ((nodeConditions[k])->GetList().lLength > 1) {
      break;
    }

  if (k == descendants - start) { // all underlying branches are "simple"
    for (k = 1; k < descendants - start; k++) {
      LocateVar(((_Operation *)((*(nodeConditions[k])).GetList()(0)))
                    ->GetAVariable())->SetFormula(*nodeConditions[0]);
      delete (nodeConditions[k]);
      nodeConditions[k] = nil;
    }
    k = 0;
  } else {

    for (l = k + 1; l < descendants - start; l++)
      if (nodeConditions[l]->GetList().lLength > 1) {
        break;
      }

    if (l == descendants - start) // all but one underlying branches are "simple"
      for (l = 0; l < descendants - start; l++) {
        if (l == k) {
          continue;
        }
        LocateVar(((_Operation *)((*(nodeConditions[l])).GetList()(0)))->GetAVariable())->SetFormula(*nodeConditions[k]);
        delete (nodeConditions[l]);
        nodeConditions[l] = nil;
      }
    // really bad bongos! must solve for non-additive constraint
    else
      for (l = 0; l < descendants - start; l++) {
        if (l == k) {
          continue;
        }
        if (nodeConditions[l]->GetList().lLength == 1) {
          LocateVar(((_Operation *)((*(nodeConditions[l])).GetList()(0)))
                        ->GetAVariable())->SetFormula(*nodeConditions[k]);
        } else { // solve for a non-additive constraint
          _Variable *nonAdd =
              LocateVar(((_Operation *)((*(nodeConditions[l])).GetList()(0)))->GetAVariable());
          nodeConditions[l]->GetList().Delete(0);
          _Formula newConstraint;
          newConstraint.Duplicate((BaseRef) nodeConditions[k]);
          _Operation mins( HY_OP_CODE_SUB, 2L);
          for (unsigned long op_index = 0; op_index < nodeConditions[l]->GetList().lLength; op_index++) {
            _Operation *curOp = (_Operation *)(*nodeConditions[l]).GetList()(op_index);
            if (curOp->GetOpKind() == _HY_OPERATION_BUILTIN) {
              newConstraint.GetList() && &mins;
            } else {
              newConstraint.GetList() << curOp;
            }
          }
          delete (nodeConditions[l]);
          nodeConditions[l] = nil;
          nonAdd->SetFormula(newConstraint);
        }
      }
  }

  if (!first) {
    _Formula *result = nodeConditions[k];
  
    result->GetList().AppendNewInstance(new _Operation( _HY_OPERATION_VAR, 
        iVariables->lData[f - 1], _HY_OPERATION_INVALID_REFERENCE, NULL));
    result->GetList().AppendNewInstance(new _Operation( HY_OP_CODE_ADD, 2L));

    free(nodeConditions);
    return result;
  }

  for (k = 0; k < descendants - start; k++)
    if (nodeConditions[k]) {
      delete nodeConditions[k];
    }

  free(nodeConditions);
  return nil;
}

//______________________________________________________________________________
void _CalcNode::SetupCategoryMap(_List &containerVariables,
                                 _SimpleList &classCounter,
                                 _SimpleList &multipliers) {

  long totalCategories = classCounter.Element(-1),
       globalCatCount = containerVariables.lLength - 1, localCategories = 1,
       catCount = categoryVariables.lLength - 1, entriesPerCat = 2 + catCount;

  //for (long k = 0; k<containerVariables.lLength;k++)
  //  printf ("%d %s\n", k,
  // ((_Variable*)containerVariables(k))->GetName()->sData);

  if (catCount < 0) {
    remapMyCategories.Clear();
  } else {

    remapMyCategories.Populate(totalCategories * entriesPerCat, 0, 0);

    _SimpleList remappedIDs, rateMultiplers(categoryVariables.lLength, 1, 0),
        categoryValues(globalCatCount + 1, 0, 0);

    for (long myCatID = 0; myCatID <= catCount; myCatID++) {
      long coordinate = containerVariables.FindPointer(
          LocateVar(categoryVariables.lData[myCatID]));
      if (coordinate < 0) {
        WarnError("Internal error in SetupCategoryMap. Please report to "
                  "spond@ucsd.edu");
      }
      localCategories *= classCounter.lData[coordinate];
      //printf ("%d %s\n", myCatID,
      //LocateVar(categoryVariables.lData[myCatID])->GetName()->sData);
      remappedIDs << coordinate;
    }

    for (long myCatID = catCount - 1; myCatID >= 0; myCatID--) {
      rateMultiplers.lData[myCatID] =
          rateMultiplers.lData[myCatID + 1] *
          classCounter.lData[remappedIDs.lData[myCatID + 1]];
    }

    for (long currentRateCombo = 0; currentRateCombo < totalCategories;
         currentRateCombo++) {
      long copyRateCombo = currentRateCombo;
      for (long variableID = 0; variableID <= globalCatCount; variableID++) {
        categoryValues.lData[variableID] =
            copyRateCombo / multipliers.lData[variableID];
        copyRateCombo = copyRateCombo % multipliers.lData[variableID];
        //printf ("%d %d %d %d\n", currentRateCombo, variableID,
        //multipliers.lData[variableID], categoryValues.lData[variableID]);
      }

      long localCatID = 0;

      for (long localVariableID = 0; localVariableID <= catCount; localVariableID++) {
        localCatID += rateMultiplers.lData[localVariableID] *
                      categoryValues.lData[remappedIDs.lData[localVariableID]];
      }

      long offset = currentRateCombo * entriesPerCat;
      remapMyCategories.lData[offset] = localCatID;

      offset++;
      for (long localVariableID = 0; localVariableID <= catCount; localVariableID++) {
        remapMyCategories[offset++] =
            categoryValues.lData[remappedIDs.lData[localVariableID]];
      }

    }
  }

  //printf ("Node remap at %s yielded %s\n", GetName()->sData,
  //_String((_String*)remapMyCategories.toStr()).sData);

}

//______________________________________________________________________________
_VariableContainer *_CalcNode::ParentTree(void) {
  _String parentTree = ParentObjectName();
  _VariableContainer *theParent =
      (_VariableContainer *)FetchVar(LocateVarByName(parentTree));
  if (theParent && theParent->ObjectClass() == TREE) {
    return theParent;
  }
  return nil;
}

