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

#include "hy_globals.h"

#include "thetree.h"
#include "calcnode_globals.h"
#include "growingvector.h"

#include "math.h"
#include "float.h"
#include "ctype.h"
#include "string.h"
#include "calcnode.h"
#include "scfg.h"
#include "legacy_parser.h"

#include "category.h"
#include "batchlan.h"
#include "likefunc.h"

//______________________________________________________________________________
#define BAD_BRANCH_LENGTH 0.000000001

//______________________________________________________________________________
_Parameter computeChordLength(_Parameter l, _Parameter angle,
                              _Parameter *maxCoord = nil) {

  _Parameter sinV = sin(angle), cosV = cos(angle);

  if (maxCoord) {
    maxCoord[0] = MAX(maxCoord[0], cosV * l);
    maxCoord[1] = MIN(maxCoord[1], cosV * l);
    maxCoord[2] = MAX(maxCoord[2], sinV * l);
    maxCoord[3] = MIN(maxCoord[3], sinV * l);
  }

  return l / MAX(fabs(sinV), fabs(cosV));
}

_TheTree::_TheTree() {
  categoryCount = 1;
  rootIChildrenCache = nil;
  marginalLikelihoodCache = nil;
  nodeMarkers = nil;
  nodeStates = nil;
  aCache = nil;
#if USE_SCALING_TO_FIX_UNDERFLOW
  scalingForUnderflow = nil;
#endif
} // default constructor - doesn't do much

//______________________________________________________________________________
_TheTree::~_TheTree(void) {
  if (rootIChildrenCache) {
    free(rootIChildrenCache);
    rootIChildrenCache = nil;
  }
  if (marginalLikelihoodCache) {
    free(marginalLikelihoodCache);
    marginalLikelihoodCache = nil;
  }
  if (nodeMarkers) {
    free(nodeMarkers);
    nodeMarkers = nil;
  }
  if (nodeStates) {
    free(nodeStates);
    nodeMarkers = nil;
  }
  DeleteObject(aCache);

}

//______________________________________________________________________________
void _TheTree::PurgeTree(void) {
  _CalcNode *curNode = DepthWiseTraversal(TRUE), *nextNode;
  nextNode = DepthWiseTraversal();
  // this loop deletes the data & the structure
  while (nextNode) {
    DeleteVariable(*curNode->GetName());
    curNode = nextNode;
    nextNode = DepthWiseTraversal();
    delete currentNode;
  }

  DeleteObject(curNode); // error checking is implicit
}

//______________________________________________________________________________
void _TheTree::PreTreeConstructor(bool) {
  rooted = UNROOTED;
  rootIChildrenCache = nil;
  marginalLikelihoodCache = nil;
  nodeMarkers = nil;
  nodeStates = nil;
  categoryCount = 1;

  aCache = new _AVLListXL(new _SimpleList);

  convertedMatrixExpressionsL.ClearFormulasInList();
  convertedMatrixExpressions.Clear();

  iNodePrefix = "Node";
  _PMathObj iv = FetchObjectFromVariableByType(&internalNodePrefix, STRING);
  if (iv) {
    iNodePrefix = *((_FString *)iv)->theString;
  }
  checkParameter(ignoreUserINames, ignoringInternalNames, 0.0);

}

//______________________________________________________________________________
void _TheTree::PostTreeConstructor(bool dupMe) {
  _Parameter acceptRTs = 0.0;
  checkParameter(acceptRootedTrees, acceptRTs, 0.0);

  DeleteObject(aCache->dataList);
  DeleteObject(aCache);
  aCache = nil;

  convertedMatrixExpressionsL.ClearFormulasInList();
  convertedMatrixExpressions.Clear();

  while (theRoot->get_num_nodes() == 1) { 
    // dumb tree w/ an extra top level node
    // remove the root
    node<long> *node_temp = theRoot->go_down(1);
    if (!node_temp) {
      WarnError(_String("Vacuos Tree Supplied"));
      isDefiningATree = 0;
      return;
    }
    if (node_temp->get_num_nodes()) {
      _String pp = *LocateVar(theRoot->get_data())->theName;
      DeleteVariable(pp);
      delete node_temp->get_parent();
      node_temp->detach_parent();
      theRoot = node_temp;
    } else {
      break;
    }
  }

  if (theRoot->get_num_nodes() == 2) { // rooted tree - check
    if (acceptRTs < 0.1) {
      long i;
      for (i = 1; i <= 2; i++) {
        node<long> *node_temp = theRoot->go_down(i);
        if (node_temp->get_num_nodes()) { // an internal node - make it a root
          node_temp->detach_parent();
          //node_temp = theRoot->go_down((i==1)?2:1);
          _CalcNode *thecn = (_CalcNode *)LocateVar(theRoot->get_data());
          _String pp = *thecn->theName;
          DeleteVariable(pp);
          if (i == 1) {
            node_temp->add_node(*theRoot->go_down(2));
            delete theRoot;
            theRoot = node_temp;
            rooted = ROOTED_LEFT;
          } else {
            node_temp->prepend_node(*theRoot->go_down(1));
            delete theRoot;
            theRoot = node_temp;
            rooted = ROOTED_RIGHT;
          }
          thecn = (_CalcNode *)LocateVar(theRoot->get_data());
          pp = *thecn->theName;
          DeleteVariable(pp, false);
          if (i == 1) {
            ReportWarning(
                _String("Rooted tree. Removing one branch - the left root "
                        "child has been promoted to be the new root"));
          } else {
            ReportWarning(
                _String("Rooted tree. Removing one branch - the right root "
                        "child has been promoted to be the new root"));
          }
          break;
        }
      }
      if (i == 3) {
        ReportWarning((_String("One branch tree supplied - hopefully this IS "
                               "what you meant to do.")));
        node<long> *node_temp = theRoot->go_down(1);
        node_temp->detach_parent();
        _CalcNode *thecn = (_CalcNode *)LocateVar(theRoot->get_data());
        _String pp = *thecn->theName;
        DeleteVariable(pp);
        node_temp->add_node(*theRoot->go_down(2));
        delete theRoot;
        theRoot = node_temp;
        rooted = ROOTED_LEFT;
        thecn = (_CalcNode *)LocateVar(theRoot->get_data());
        pp = *thecn->theName;
        DeleteVariable(pp, false);
        ReportWarning(
            _String("Rooted tree. Removing one branch - the left root child "
                    "has been promoted to be the new root"));
        //PurgeTree();
      }
    }
  }

  if (!theRoot) {
    WarnError("Invalid tree/topology string specification.");
  } else {
    /*_CalcNode* theN = StepWiseTraversal(TRUE);
    while (theN)
    {
        theN->SetCodeBase(4);
        theN = StepWiseTraversal();
    }*/
    BaseRef temp = variablePtrs(theIndex);

    if (dupMe) {
      variablePtrs[theIndex] = this->makeDynamic();
    } else {
      variablePtrs[theIndex] = this;
    }

    DeleteObject(temp);
  }
}

//______________________________________________________________________________
_TheTree::_TheTree(_String name, _String &parms, bool dupMe)
    : _TreeTopology(&name) {
  PreTreeConstructor(dupMe);
  if (MainTreeConstructor(parms)) {
    PostTreeConstructor(dupMe);
  }
}

//______________________________________________________________________________
_TheTree::_TheTree(_String name, _TreeTopology *top) : _TreeTopology(&name) {
  PreTreeConstructor(false);
  if (top->theRoot) {
    isDefiningATree = 1;
    theRoot = top->theRoot->duplicate_tree();
    node<long> *topTraverser = DepthWiseStepTraverser(theRoot);
    while (topTraverser) {
      _Parameter nodeVal = top->compExp->theData[topTraverser->in_object];
      _String nodeVS,
          nodeName(*(_String *)top->flatTree(topTraverser->in_object)),
          nodeParams(*(_String *)top->flatCLeaves(topTraverser->in_object));

      if (!nodeName.IsValidIdentifier(true)) {
        nodeName.ConvertToAnIdent(true);
      }

      if (nodeVal != 0.0) {
        nodeVS = nodeVal;
      }

      FinalizeNode(topTraverser, 0, nodeName, nodeParams, nodeVS);
      topTraverser = DepthWiseStepTraverser((node<long> *)nil);
    }
    isDefiningATree = 0;
    PostTreeConstructor(false);
  } else {
    WarnError("Can't create an empty tree");
    return;
  }
}

//______________________________________________________________________________
_TheTree::_TheTree(_String name, _TheTree *otherTree) : _TreeTopology(&name) {
  PreTreeConstructor(false);
  if (otherTree->theRoot) {
    isDefiningATree = 1;
    theRoot = otherTree->theRoot->duplicate_tree();

    node<long> *topTraverser = DepthWiseStepTraverser(theRoot);

    while (topTraverser) {
      _CalcNode *sourceNode = (_CalcNode *)LocateVar(topTraverser->in_object),
                copiedNode(sourceNode, this);
      topTraverser->init(copiedNode.theIndex);
      topTraverser = DepthWiseStepTraverser((node<long> *)nil);
    }

    isDefiningATree = 0;
    PostTreeConstructor(false);
  } else {
    WarnError("Can't create an empty tree");
    return;
  }
}

//______________________________________________________________________________
bool _TheTree::FinalizeNode(node<long> *nodie, long number, _String &nodeName,
                            _String &nodeParameters, _String &nodeValue,
                            _String *nodeComment) {

  bool isAutoGenerated =
      (nodeName.sLength == 0 ||
       !CheckEqual(ignoringInternalNames, 0.0) && nodie->get_num_nodes() > 0);
  if (isAutoGenerated) {
    nodeName = iNodePrefix & number;
  } else {

    if (!nodeName.IsValidIdentifier(false)) {
      _String oldName(nodeName);
      nodeName.ConvertToAnIdent();
      ReportWarning(_String("Automatically renamed ") & oldName & " to " &
                    nodeName & " in order to create a valid HyPhy identifier");
    }
  }
  if (nodie == theRoot) {
    nodeParameters = empty;
    nodeValue = empty;
  } else {
    if (!nodeParameters.sLength && lastMatrixDeclared != -1) {
      nodeParameters = *(((_String **)modelNames.lData)[lastMatrixDeclared]);
    }

    if (nodeParameters.sLength) {
      ReportWarning((_String("Model ") & nodeParameters &
                     _String(" assigned to ") & nodeName));
    } else {
      ReportWarning(_String("No nodel was assigned to ") & nodeName);
    }

  }

  char saveIDT = isDefiningATree;
  isDefiningATree = isAutoGenerated ? 2 : 1;
  _CalcNode cNt(nodeName, nodeParameters, 4, this, aCache);
  isDefiningATree = saveIDT;
  nodie->init(cNt.theIndex);

  _Constant val(nodeValue.ProcessTreeBranchLength());

  if (nodeValue.Length() && takeBranchLengths) {
    if (cNt.iVariables &&
        cNt.iVariables->lLength == 2) { // can assign default values
      bool setDef = true;
      if (autoSolveBranchLengths) {
        long nodeModelID = cNt.GetModelIndex();
        if (nodeModelID != HY_NO_MODEL) {
          _Formula *expressionToSolveFor = nil;
          long alreadyConverted =
              convertedMatrixExpressions.Find((BaseRef) nodeModelID);
          if (alreadyConverted < 0) {
            _Variable *tV, *tV2;
            bool mByF;
            RetrieveModelComponents(nodeModelID, tV, tV2, mByF);
            _String *result = ((_Matrix *)tV->GetValue())
                ->BranchLengthExpression((_Matrix *)tV2->GetValue(), mByF);
            if (result->sLength) {
              expressionToSolveFor = new _Formula(*result);
              for (unsigned long cc = 0; cc < cNt.categoryVariables.lLength;
                   cc++) {
                _CategoryVariable *thisCC = (_CategoryVariable *)LocateVar(
                    cNt.categoryVariables.lData[cc]);
                thisCC->SetValue(new _Constant(thisCC->Mean()), false);

              }
            }
            DeleteObject(result);
          } else {
            expressionToSolveFor = (_Formula *)convertedMatrixExpressions.GetXtra(alreadyConverted);
          }

          if (expressionToSolveFor != nil) {
            _Variable *solveForMe = LocateVar(cNt.iVariables->lData[1]);
            _Parameter modelP = expressionToSolveFor->Brent(
                solveForMe, solveForMe->GetLowerBound(),
                solveForMe->GetUpperBound(), 1e-6, nil, val.Value());
            ReportWarning(_String("Branch parameter of ") & nodeName &
                          " set to " & modelP);
            LocateVar(cNt.iVariables->lData[0])->SetValue(new _Constant(modelP), false);
            setDef = false;
          }
        }
      }

      if (setDef) {
        LocateVar(cNt.iVariables->lData[0])->SetValue(&val);
        ReportWarning(_String("Branch parameter of ") & nodeName & " set to " &
                      nodeValue);
      }
    } else {
      ReportWarning(
          nodeName & " has " &
          _String((long)(cNt.iVariables ? cNt.iVariables->lLength / 2 : 0)) &
          " parameters - branch length not assigned");
    }
  }

  _CalcNode *nodeVar = (_CalcNode *)LocateVar(cNt.theIndex);

  nodeVar->SetValue(&val);

  nodeName = empty;
  nodeParameters = empty;
  nodeValue = empty;
  if (nodeComment && nodeComment->sLength) {
    _String commentName = *nodeVar->GetName() & "._comment";
    CheckReceptacleAndStore(&commentName, empty, false,
                            new _FString(*nodeComment));
    *nodeComment = empty;
  }

  nodeVar->categoryVariables.TrimMemory();
  nodeVar->categoryIndexVars.TrimMemory();
  nodeVar->_VariableContainer::TrimMemory();

  return true;
}

//______________________________________________________________________________
_CalcNode *_TheTree::DepthWiseTraversal(bool init) {
  DepthWiseT(init);

  if (currentNode) {
    return (_CalcNode *)(((BaseRef *)variablePtrs.lData)[currentNode->in_object]);
  }

  return nil;
}

//______________________________________________________________________________
_CalcNode *_TheTree::DepthWiseTraversalRight(bool init) {

  DepthWiseTRight(init);

  if (currentNode) {
    return (_CalcNode *)(((BaseRef *)variablePtrs.lData)[currentNode->in_object]);
  }

  return nil;
}

//______________________________________________________________________________
_CalcNode *_TheTree::LeafWiseTraversal(bool init) {
  LeafWiseT(init);
  if (currentNode) {
    return (_CalcNode *)(
        ((BaseRef *)variablePtrs.lData)[currentNode->in_object]);
  }
  return nil;
}

//______________________________________________________________________________
_CalcNode *_TheTree::StepWiseTraversal(bool init) {
  StepWiseT(init);

  if (currentNode) {
    return (_CalcNode *)(
        ((BaseRef *)variablePtrs.lData)[currentNode->in_object]);
  }

  return nil;
}

//______________________________________________________________________________
_CalcNode *_TheTree::StepWiseTraversalLevel(long &level, bool init) {
  StepWiseTLevel(level, init);

  if (currentNode) {
    return (_CalcNode *)(
        ((BaseRef *)variablePtrs.lData)[currentNode->in_object]);
  }
  return nil;
}

//______________________________________________________________________________
_CalcNode *_TheTree::DepthWiseTraversalLevel(long &level, bool init) {
  DepthWiseTLevel(level, init);

  if (currentNode) {
    return (_CalcNode *)(
        ((BaseRef *)variablePtrs.lData)[currentNode->in_object]);
  }
  return nil;
}

//______________________________________________________________________________
BaseRef _TheTree::makeDynamic(void) {
  _TheTree *res = new _TheTree;
  checkPointer(res);
  res->_CalcNode::Duplicate(this);

  res->currentNode = currentNode;
  res->rooted = rooted;
  res->categoryCount = 1;
  res->theRoot = CopyTreeStructure(theRoot, true);
  return res;
}

//______________________________________________________________________________
BaseRef _TheTree::makeDynamicCopy(_String *replacementName) {
  _TheTree *res = new _TheTree;
  checkPointer(res);

  res->rooted = rooted;
  if (theRoot) {
    _String rn = *replacementName & '.';
    res->theRoot = DuplicateTreeStructure(theRoot, &rn, true);
  } else {
    res->theRoot = nil;
  }

  res->SetIndex(variableNames.GetXtra(LocateVarByName(*replacementName)));
  res->theName = replacementName;
  res->theName->nInstances++;
  return res;
}

//______________________________________________________________________________
node<long> *_TheTree::DuplicateTreeStructure(node<long> *theNode,
                                             _String *replacementName, bool) {
  long i, j;
  node<long> *locNode = new node<long>;
  for (i = 0; i < theNode->get_num_nodes(); i++) {
    locNode->add_node(*DuplicateTreeStructure(theNode->go_down(i + 1),
                                              replacementName, false));
  }
  if (1) {
    // process this node now
    _String replacedName = *GetName() & '.', *temp;
    _CalcNode *sourceNode = (_CalcNode *)LocateVar(theNode->get_data());
    sourceNode = (_CalcNode *)sourceNode->makeDynamic();
    _String newNodeName = (LocateVar(sourceNode->GetAVariable())->GetName())
        ->Replace(replacedName, *replacementName, true);
    _Variable dummyVar(newNodeName);
    j = dummyVar.GetAVariable();
    temp = sourceNode->GetName();
    DeleteObject(temp);
    sourceNode->theName = dummyVar.GetName();
    sourceNode->theName->nInstances++;
    ReplaceVar(sourceNode);
    DeleteObject(sourceNode);
    sourceNode = (_CalcNode *)LocateVar(j);
    locNode->init(j);

#ifndef USE_POINTER_VC
    for (i = 0; i < sourceNode->independentVars.lLength; i++) {
      newNodeName = (LocateVar(sourceNode->independentVars.lData[i])->GetName())
          ->Replace(replacedName, *replacementName, true);
      _Variable dummyVar(newNodeName);
#ifndef USE_AVL_NAMES
      sourceNode->independentVars.lData[i] =
          variableReindex.lData[LocateVarByName(newNodeName)];
#else
      sourceNode->independentVars.lData[i] =
          variableNames.GetXtra(LocateVarByName(newNodeName));
#endif
    }

    // done with independents - now set the dependancies
    for (i = 0; i < sourceNode->dependentVars.lLength; i++) {
      newNodeName = (LocateVar(sourceNode->dependentVars.lData[i])->GetName())
          ->Replace(replacedName, *replacementName, true);
      _Variable dummyVar(newNodeName);
#ifndef USE_AVL_NAMES
      sourceNode->dependentVars.lData[i] =
          variableReindex.lData[LocateVarByName(newNodeName)];
#else
      sourceNode->dependentVars.lData[i] =
          variableNames.GetXtra(LocateVarByName(newNodeName));
#endif
      _String *newFormula =
          LocateVar(sourceNode->dependentVars.lData[i])->GetFormulaString();
      *newFormula = newFormula->Replace(replacedName, *replacementName, true);
      _Formula dummyF(*newFormula);
      LocateVar(sourceNode->dependentVars.lData[i])->SetFormula(dummyF);
      DeleteObject(newFormula);
    }
#else
    if (sourceNode->iVariables)
      for (i = 0; i < sourceNode->iVariables->lLength; i += 2) {
        newNodeName = (LocateVar(sourceNode->iVariables->lData[i])->GetName())
            ->Replace(replacedName, *replacementName, true);
        _Variable dummyVar(newNodeName);
#ifndef USE_AVL_NAMES
        sourceNode->iVariables->lData[i] =
            variableReindex.lData[LocateVarByName(newNodeName)];
#else
        sourceNode->iVariables->lData[i] =
            variableNames.GetXtra(LocateVarByName(newNodeName));
#endif
      }

    // done with independents - now set the dependancies

    if (sourceNode->dVariables)
      for (i = 0; i < sourceNode->dVariables->lLength; i += 2) {
        newNodeName = (LocateVar(sourceNode->dVariables->lData[i])->GetName())
            ->Replace(replacedName, *replacementName, true);
        _Variable dummyVar(newNodeName);
#ifndef USE_AVL_NAMES
        sourceNode->dVariables->lData[i] =
            variableReindex.lData[LocateVarByName(newNodeName)];
#else
        sourceNode->dVariables->lData[i] =
            variableNames.GetXtra(LocateVarByName(newNodeName));
#endif
        _String *newFormula =
            LocateVar(sourceNode->dVariables->lData[i])->GetFormulaString();
        *newFormula = newFormula->Replace(replacedName, *replacementName, true);
        _Formula dummyF(*newFormula);
        LocateVar(sourceNode->dVariables->lData[i])->SetFormula(dummyF);
        DeleteObject(newFormula);
      }

#endif
  }
  return locNode;
}

//______________________________________________________________________________
void _TheTree::GetNodeName(node<long> *n, _String &r, bool fullName) {
  if (fullName) {
    r = *((_CalcNode *)(((BaseRef *)variablePtrs.lData)[n->in_object]))
        ->GetName();
  } else {
    r = ((_CalcNode *)(((BaseRef *)variablePtrs.lData)[n->in_object]))
        ->GetName()->Cut(GetName()->sLength + 1, -1);
  }
}

//______________________________________________________________________________
BaseRef _TheTree::toStr(void) {

  _String *res = new _String((unsigned long) 128, true), num;

  _Parameter skipILabels, includeMSP;

  checkParameter(noInternalLabels, skipILabels, 0.0);
  checkParameter(includeModelSpecs, includeMSP, 0.0);

  if (IsDegenerate()) {
    _CalcNode *curNode = DepthWiseTraversal(true),
              *nextNode = DepthWiseTraversal(false);

    long l1 = GetName()->Length();

    (*res) << '(';
    num = nextNode->GetName()->Cut(l1 + 1, -1);
    (*res) << &num;
    if (includeMSP > 0.5) {
      long midx = curNode->GetModelIndex();
      if (midx != HY_NO_MODEL) {
        (*res) << '{';
        (*res) << (_String *)modelNames(midx);
        (*res) << '}';
      }
    }

    (*res) << ',';
    num = curNode->GetName()->Cut(l1 + 1, -1);
    (*res) << &num;
    if (includeMSP > 0.5) {
      long midx = nextNode->GetModelIndex();
      if (midx != HY_NO_MODEL) {
        (*res) << '{';
        (*res) << (_String *)modelNames(midx);
        (*res) << '}';
      }
    }
    (*res) << ')';
  } else {

    long level = 0, myLevel = 0, lastLevel = 0, l1 = GetName()->Length(), j;
    _CalcNode *curNode = DepthWiseTraversalLevel(myLevel, true), *nextNode;
    level = myLevel;
    bool isCTip = IsCurrentNodeATip(), isCTip2;
    nextNode = DepthWiseTraversalLevel(myLevel);
    isCTip2 = IsCurrentNodeATip();

    while (nextNode) {
      if (level > lastLevel) {
        if (lastLevel) {
          (*res) << ',';
        }
        for (j = 0; j < level - lastLevel; j++) {
          (*res) << '(';
        }
      } else if (level < lastLevel) {
        for (j = 0; j < lastLevel - level; j++) {
          (*res) << ')';
        }
      } else {
        (*res) << ',';
      }

      if ((skipILabels < 0.1) || isCTip) {
        num = curNode->GetName()->Cut(l1 + 1, -1);
        (*res) << &num;
      }
      if (includeMSP > 0.5) {
        long midx = curNode->GetModelIndex();
        if (midx != HY_NO_MODEL) {
          (*res) << '{';
          (*res) << (_String *)modelNames(midx);
          (*res) << '}';
        }
      }

      lastLevel = level;
      level = myLevel;
      curNode = nextNode;
      isCTip = isCTip2;
      nextNode = DepthWiseTraversalLevel(myLevel);
      isCTip2 = IsCurrentNodeATip();

    }
    for (j = 0; j < lastLevel - level; j++) {
      (*res) << ')';
    }
  }
  (*res) << ';';
  (*res).Finalize();
  return res;
}

//______________________________________________________________________________
void _TheTree::CompileListOfModels(_SimpleList &l) {
  _CalcNode *curNode = DepthWiseTraversal(true);
  while (curNode) {
    long modelID = curNode->GetModelIndex();
    if (modelID != HY_NO_MODEL && l.Find(modelID) == -1) {
      l << modelID;
    }
    curNode = DepthWiseTraversal(false);
  }
}

//______________________________________________________________________________
void _TheTree::SetCompMatrices(long catID) {
  _CalcNode *travNode = DepthWiseTraversal(TRUE);

  while (!IsCurrentNodeTheRoot()) {
    travNode->SetCompMatrix(catID);
    travNode = DepthWiseTraversal();
  }
}

//______________________________________________________________________________
void _TheTree::SetUp(void) {
  _CalcNode *travNode = DepthWiseTraversal(TRUE);

  if (marginalLikelihoodCache) {
    free(marginalLikelihoodCache);
    marginalLikelihoodCache = nil;
  }

  if (nodeMarkers) {
    free(nodeMarkers);
    nodeMarkers = nil;
  }

  if (nodeStates) {
    free(nodeStates);
    nodeMarkers = nil;
  }

  flatTree.Clear();
  flatNodes.Clear();
  flatLeaves.Clear();
  flatCLeaves.Clear();

  flatParents.Clear();
  _SimpleList flatINodeParents;

  while (travNode) {
    if (!IsCurrentNodeATip()) {
      flatTree << travNode;
      flatNodes << (long)(currentNode);
      travNode->lastState = -1;
      flatINodeParents << (long)(currentNode->parent);
    } else {
      flatLeaves << (long)(currentNode);
      flatCLeaves << travNode;
      flatParents << (long)(currentNode->parent);
    }
    travNode = DepthWiseTraversal();
  }

  flatParents << flatINodeParents;
  _SimpleList parentlist(flatNodes), indexer(flatNodes.lLength, 0, 1);
  SortLists(&parentlist, &indexer);
  for (long k = 0; k < flatParents.lLength; k++)
    if (flatParents.lData[k]) {
      flatParents.lData[k] =
          indexer.lData[parentlist.BinaryFind(flatParents.lData[k])];
    } else {
      flatParents.lData[k] = -1;
    }

  if (cBase > 0) {
    marginalLikelihoodCache =
        (_Parameter *)MemAllocate((flatNodes.lLength + flatLeaves.lLength) *
                                  sizeof(_Parameter) * cBase * systemCPUCount);
  }
  nodeStates = (long *)MemAllocate((flatNodes.lLength + flatLeaves.lLength) *
                                   sizeof(long) * systemCPUCount);
  nodeMarkers =
      (char *)MemAllocate(flatNodes.lLength * sizeof(char) * systemCPUCount);

  long iNodeCounter = 0, leafCounter = 0;

  travNode = DepthWiseTraversal(TRUE);

  while (travNode) {
    if (IsCurrentNodeATip()) {
      travNode->nodeIndex = leafCounter++;
    } else {
      nodeMarkers[iNodeCounter] = -1;
      for (long k = 1; k < systemCPUCount; k++) {
        nodeMarkers[iNodeCounter + k * flatNodes.lLength] = -1;
      }
      travNode->nodeIndex = flatLeaves.lLength + iNodeCounter++;
      nodeStates[travNode->nodeIndex] = -1;
      for (long m = 1; m < systemCPUCount; m++) {
        nodeStates[travNode->nodeIndex +
                   m * (flatNodes.lLength + flatLeaves.lLength)] = -1;
      }
    }
    travNode = DepthWiseTraversal();
  }

  BuildINodeDependancies();
}

//______________________________________________________________________________
bool _TheTree::AllBranchesHaveModels(long matchSize) {
  _CalcNode *travNode;
  travNode = DepthWiseTraversal(TRUE);
  if (matchSize > 0) {
    while (!IsCurrentNodeTheRoot()) {
      long mID = travNode->GetTheModelID();
      if (mID < 0) {
        return false;
      }
      travNode = DepthWiseTraversal();
    }
  } else {
    while (!IsCurrentNodeTheRoot()) {
      long mID = travNode->GetTheModelID();
      if (mID < 0) {
        return false;
      } else {
        if (travNode->GetModelMatrix()->GetHDim() != matchSize) {
          return false;
        }
      }
      travNode = DepthWiseTraversal();
    }

  }
  return true;
}

//______________________________________________________________________________
_String *_TheTree::TreeUserParams(void) {
  _String *result = new _String(16L, true);
  checkPointer(result);

  _CalcNode *travNode;
  travNode = DepthWiseTraversal(TRUE);
  while (travNode) {
    _String *nodeString = travNode->GetSaveableListOfUserParameters();
    if (nodeString->sLength) {
      *result << nodeString;
    }
    DeleteObject(nodeString);
    travNode = DepthWiseTraversal();
  }

  result->Finalize();
  return result;
}

//______________________________________________________________________________
// execute this operation with the second arg if necessary
_PMathObj _TheTree::Execute(long opCode, _PMathObj p, _PMathObj p2, _hyExecutionContext *context,
  _PMathObj ) {
  switch (opCode) {
  case HY_OP_CODE_PSTREESTRING: //PlainTreeString
    return PlainTreeString(p, p2);
    break;
  case HY_OP_CODE_TEXTREESTRING: // TEXTreeString
    return TEXTreeString(p);
    break;
  case HY_OP_CODE_TYPE: // Type
    return Type();
    break;
  }

  return _TreeTopology::Execute(opCode, p, p2, context);

}

//______________________________________________________________________________
_String _TheTree::FindMaxCommonSubTree(_TheTree *compareTo, long &sizeVar,
                                       _List *forest) {

  _List myLeaves, otherLeaves, sharedLeaves;
  _SimpleList indexer, otherIndexer, sharedLeavesIDin1, sharedLeavesIDin2;
  node<long> *myCT, *otherCT;
  _String rerootAt;
  myCT = prepTree4Comparison(myLeaves, indexer);

  otherCT = compareTo->prepTree4Comparison(otherLeaves, otherIndexer);

  sharedLeaves.Intersect(otherLeaves, myLeaves, &sharedLeavesIDin1,
                         &sharedLeavesIDin2);

  if (sharedLeaves.lLength > 1) { 
    // more than one common leaf
    // now we need to map shared leaves to a common indexing space

    _SimpleList reindexer((unsigned long) otherLeaves.lLength),
        lidx1((unsigned long) myLeaves.lLength),
        lidx2((unsigned long) otherLeaves.lLength), ldx1, ldx2;

    reindexer.lLength = (unsigned long) otherLeaves.lLength;
    lidx1.lLength = myLeaves.lLength;
    lidx2.lLength = otherLeaves.lLength;

    for (long k = 0; k < otherLeaves.lLength; k++) {
      lidx2.lData[otherIndexer.lData[k]] = k;
    }

    for (long k0 = 0; k0 < myLeaves.lLength; k0++) {
      lidx1.lData[indexer.lData[k0]] = k0;
    }

    for (long k1 = 0; k1 < reindexer.lLength; k1++) {
      reindexer.lData[k1] = -1;
    }

    for (long k2 = 0; k2 < sharedLeaves.lLength; k2++) {
      reindexer.lData[lidx2.lData[sharedLeavesIDin2.lData[k2]]] =
          lidx1.lData[sharedLeavesIDin1.lData[k2]];
    }

    // now we map actual leaf structures to their respective leaf indices
    node<long> *meNode = DepthWiseStepTraverser(myCT);

    while (meNode) {
      if (meNode->get_num_nodes() == 0) {
        ldx1 << (long) meNode;
      }
      meNode = DepthWiseStepTraverser((node<long> *)nil);
    }

    meNode = DepthWiseStepTraverser(otherCT);

    while (meNode) {
      if (meNode->get_num_nodes() == 0) {
        ldx2 << (long) meNode;
      }
      meNode = DepthWiseStepTraverser((node<long> *)nil);
    }

    // now we loop through the list of leaves and try to match them all up
    _SimpleList matchedTops, matchedSize;

    for (long k3 = 0; k3 < sharedLeaves.lLength - 1; k3++) {
      if (reindexer.lData[lidx2.lData[sharedLeavesIDin2.lData[k3]]] >= 0) {
        // leaf still available
        node<long> *ln1 = (node<long> *)ldx1
                       .lData[lidx1.lData[sharedLeavesIDin1.lData[k3]]],
                   *ln2 = (node<long> *)ldx2
                       .lData[lidx2.lData[sharedLeavesIDin2.lData[k3]]],
                   *p1 = ln1->parent, *p2 = ln2->parent;

        char cRes = 0;

        while ((internalTreeCompare(p1, p2, &reindexer, 0, myLeaves.lLength,
                                    p2->parent ? nil : ln2, compareTo) == 1) &&
               p1 && p2) {
          ln1 = p1;
          ln2 = p2;
          p1 = p1->parent;
          p2 = p2->parent;

          cRes = 1;
        }
        if (cRes) {
          _SimpleList *matchedLeaves = (_SimpleList *)ln2->in_object;

          matchedTops << (long) ln1;
          matchedSize << matchedLeaves->lLength;

          for (long k4 = 0; k4 < matchedLeaves->lLength; k4++) {
            reindexer.lData[matchedLeaves->lData[k4]] = -1;
          }
        }
      }
    }

    if (matchedSize.lLength) {
      if (forest) {
        sizeVar = 0;
        for (long k6 = 0; k6 < matchedSize.lLength; k6++) {
          long maxSz = 0;

          node<long> *meNode = DepthWiseStepTraverser(myCT),
                     *mNode = (node<long> *)matchedTops.lData[k6];

          while (meNode) {
            if (meNode->get_num_nodes()) {
              if (meNode == mNode) {
                break;
              }
              maxSz++;
            }
            meNode = DepthWiseStepTraverser((node<long> *)nil);
          }

          meNode = DepthWiseStepTraverser(theRoot);

          while (meNode) {
            if (meNode->get_num_nodes()) {
              if (maxSz == 0) {
                (*forest) << LocateVar(meNode->in_object)->GetName();
                break;
              }
              maxSz--;
            }
            meNode = DepthWiseStepTraverser((node<long> *)nil);
          }
          sizeVar += matchedSize.lData[k6];
        }
      } else {
        long maxSz = -1, maxIdx = 0;

        for (long k5 = 0; k5 < matchedSize.lLength; k5++)
          if (matchedSize.lData[k5] > maxSz) {
            maxSz = matchedSize.lData[k5];
            maxIdx = k5;
          }

        sizeVar = maxSz;

        maxSz = 0;

        node<long> *meNode = DepthWiseStepTraverser(myCT),
                   *mNode = (node<long> *)matchedTops.lData[maxIdx];

        while (meNode) {
          if (meNode->get_num_nodes()) {
            if (meNode == mNode) {
              break;
            }
            maxSz++;
          }
          meNode = DepthWiseStepTraverser((node<long> *)nil);
        }

        meNode = DepthWiseStepTraverser(theRoot);

        while (meNode) {
          if (meNode->get_num_nodes()) {
            if (maxSz == 0) {
              return *LocateVar(meNode->in_object)->GetName();
            }
            maxSz--;
          }
          meNode = DepthWiseStepTraverser((node<long> *)nil);
        }
      }
    }
  }
  return empty;
}

//______________________________________________________________________________
// compare tree topologies
_String _TheTree::CompareSubTrees(_TheTree *compareTo, node<long> *topNode) {

  _List myLeaves, otherLeaves, sharedLeaves;
  _SimpleList indexer, otherIndexer, sharedLeavesID;
  node<long> *myCT, *otherCT;
  _String rerootAt;
  myCT = prepTree4Comparison(myLeaves, indexer);
  otherCT = compareTo->prepTree4Comparison(otherLeaves, otherIndexer, topNode);

  sharedLeaves.Intersect(myLeaves, otherLeaves, &sharedLeavesID);

  // first compare the inclusion for the set of leaf labels
  if (sharedLeavesID.lLength == otherLeaves.lLength) {
    _SimpleList ilist((unsigned long) myLeaves.lLength);
    ilist.lLength = myLeaves.lLength;
    {
      //for BCC
      for (long k2 = 0; k2 < ilist.lLength; k2++) {
        ilist.lData[k2] = -1;
      }
    }
    {
      //for BCC
      for (long k2 = 0; k2 < otherIndexer.lLength; k2++) {
        ilist.lData[sharedLeavesID.lData[otherIndexer.lData[k2]]] = k2;
      }
    }
    for (long k2 = 0; k2 < indexer.lLength; k2++) {
      long lidx = ilist.lData[indexer.lData[k2]];

      if (lidx >= 0) {
        indexer.lData[k2] = lidx;
      } else {
        indexer.lData[k2] = -1;
      }
    }

    _SimpleList *reindexer = &indexer;

    // now compare explore possible subtree matchings
    // for all internal nodes of this tree except the root

    char compRes = 0;

    long tCount = 1, nc2 = topNode->get_num_nodes();

    node<long> *meNode = DepthWiseStepTraverser(myCT);
    meNode = DepthWiseStepTraverser((node<long> *)nil);

    while (meNode != myCT) {
      long nc = meNode->get_num_nodes();
      if (nc == nc2) {
        long kk;
        for (kk = 1; kk <= nc; kk++) {
          compRes = internalTreeCompare(otherCT, meNode, reindexer, 0,
                                        otherLeaves.lLength,
                                        meNode->go_down(kk), compareTo);
          if (compRes) {
            if (compRes == -1) {
              meNode = myCT;
            }
            break;
          }
        }
        if (kk > nc) {
          compRes = internalTreeCompare(otherCT, meNode, reindexer, 0,
                                        otherLeaves.lLength, nil, compareTo);
          if (compRes) {
            if (compRes == -1) {
              meNode = myCT;
            }
            break;
          }
        } else {
          break;
        }
      }
      tCount++;
      meNode = DepthWiseStepTraverser((node<long> *)nil);
    }

    if (meNode != myCT) {
      meNode = DepthWiseStepTraverser(theRoot);
      while (meNode != theRoot) {
        if (tCount == 0) {
          rerootAt = _String("Matched at the ") &
                     *LocateVar(meNode->in_object)->GetName() & '.';
          break;
        } else {
          tCount--;
        }

        meNode = DepthWiseStepTraverser((node<long> *)nil);
      }
    } else {
      long nc = myCT->get_num_nodes();

      if (nc == nc2 + 1)
        for (long kk = 1; kk <= nc; kk++) {
          compRes = internalTreeCompare(otherCT, myCT, reindexer, 0,
                                        otherLeaves.lLength, myCT->go_down(kk),
                                        compareTo);
          if (compRes == 1) {
            break;
          }
        }

      if (compRes == 1) {
        rerootAt = "Matched at the root.";
      }
    }

    if (!rerootAt.sLength) {
      rerootAt = "No match: Different topologies (matching label sets).";
    }
  } else {
    rerootAt = "No match: Unequal label sets.";
  }

  destroyCompTree(myCT);
  destroyCompTree(otherCT);

  return rerootAt;
}

//______________________________________________________________________________
void _TheTree::GetBranchLength(node<long> *n, _String &r, bool getBL) {
  if (getBL) {
    bool mbf;

    _Matrix *mm, *fv;

    RetrieveModelComponents(
        ((_CalcNode *)(((BaseRef *)variablePtrs.lData)[n->in_object]))
            ->GetModelIndex(),
        mm, fv, mbf);

    if (mm && fv) {
      r.CopyDynamicString(mm->BranchLengthExpression(fv, mbf), true);
    } else {
      r = empty;
    }

  } else {
    r = ((_CalcNode *)(((BaseRef *)variablePtrs.lData)[n->in_object]))
        ->BranchLength();
  }
}

//______________________________________________________________________________
void _TheTree::GetBranchLength(node<long> *n, _Parameter &r) {
  r = ((_CalcNode *)(((BaseRef *)variablePtrs.lData)[n->in_object]))
      ->BranchLength();
}

//______________________________________________________________________________
void _TheTree::GetBranchValue(node<long> *n, _String &r) {
  _Parameter t =
      ((_CalcNode *)(((BaseRef *)variablePtrs.lData)[n->in_object]))->Value();
  if (t != -1.) {
    r = t;
  } else {
    r = empty;
  }
}

//______________________________________________________________________________
_String *_TheTree::GetBranchSpec(node<long> *n) {
  _CalcNode *cn = (_CalcNode *)(((BaseRef *)variablePtrs.lData)[n->in_object]);

  _String *res = new _String(32L, true);
  long nModel = cn->GetTheModelID();
  if (nModel >= 0) {
    (*res) << '{';
    (*res) << LocateVar(modelMatrixIndices.lData[nModel])->GetName();
  }

  if (iVariables && iVariables->lLength) {
    if (res->sLength) {
      (*res) << ',';
    } else {
      (*res) << '{';
    }

    for (long k = 0; k < iVariables->lLength; k += 2) {
      if (k) {
        (*res) << ',';
      }
      _Variable *av = LocateVar(iVariables->lData[k]);
      if (iVariables->lData[k + 1] >= 0) {
        (*res) << LocateVar(iVariables->lData[k + 1])->GetName();
      } else {
        (*res) << av->GetName();
      }

      (*res) << '=';
      _String vval(av->Value());
      (*res) << &vval;
    }
  }

  if (dVariables && dVariables->lLength) {
    long m = 0;
    for (long k = 0; k < dVariables->lLength; k += 2) {
      if (dVariables->lData[k + 1] < 0) {
        if (m) {
          (*res) << ',';
        } else if (res->sLength) {
          (*res) << ',';
        } else {
          (*res) << '{';
        }

        m++;
        _Variable *av = LocateVar(dVariables->lData[k]);
        (*res) << av->GetName();
        (*res) << ":=";
        (*res) << '=';
        _String *vval = av->GetFormulaString();
        (*res) << vval;
        DeleteObject(vval);
      }
    }
  }

  if (res->sLength) {
    (*res) << '}';
  }

  res->Finalize();
  return res;
}

//______________________________________________________________________________
void _TheTree::GetBranchVarValue(node<long> *n, _String &r, long idx) {
  _CalcNode *travNode =
      (_CalcNode *)(((BaseRef *)variablePtrs.lData)[n->in_object]);
  long iVarValue = travNode->iVariables->FindStepping(idx, 2, 1);
  if (iVarValue > 0)
      // user variable, but not in the model
      {
    r = _String(LocateVar(travNode->iVariables->lData[iVarValue - 1])->Value());
  } else {
    _String query = _String('.') & *LocateVar(idx)->GetName();
    for (long k = 0; k < travNode->iVariables->lLength; k += 2) {
      _Variable *localVar = LocateVar(travNode->iVariables->lData[k]);
      if (localVar->GetName()->endswith(query)) {
        r = _String(localVar->Value());
        return;
      }
    }
  }
}

//______________________________________________________________________________
void _TheTree::AlignNodes(node<nodeCoord> *theNode) {
  long k = theNode->get_num_nodes();
  if (k) {
    theNode->in_object.v = (theNode->go_down(1)->in_object.v +
                            theNode->go_down(k)->in_object.v) / 2.0;
    theNode->in_object.h = 0;
    for (; k; k--) {
      _Parameter t = theNode->go_down(k)->in_object.h;
      if (t < theNode->in_object.h) {
        theNode->in_object.h = t;
      }
    }
    theNode->in_object.h -= TREE_H_SHIFT;
  } else {
    theNode->in_object.v = 0;
    theNode->in_object.h = 0;
  }
}
//______________________________________________________________________________
node<nodeCoord> *_TheTree::AlignedTipsMapping(bool first, bool respectRoot) {
  long k;
  long descendants;
  if (first) {
    treeLayoutVert = 0.0;
    descendants = theRoot->get_num_nodes();
    node<nodeCoord> *aRoot = new node<nodeCoord>; // reslove rootedness here
    aRoot->in_object.varRef = -1;
    aRoot->go_next(); // dumbass initialization for later
    if ((rooted == UNROOTED) || (!respectRoot)) {
      for (k = 1; k <= descendants; k++) {
        currentNode = theRoot->go_down(k);
        aRoot->add_node(*AlignedTipsMapping());
      }
      AlignNodes(aRoot);
      return aRoot;
    } else {
      node<nodeCoord> *aChild = new node<nodeCoord>;
      aChild->in_object.varRef = -1;
      if (rooted == ROOTED_LEFT) {
        aRoot->add_node(*aChild);
        for (k = 1; k < descendants; k++) {
          currentNode = theRoot->go_down(k);
          aChild->add_node(*AlignedTipsMapping());
        }
        currentNode = theRoot->go_down(k);
        aRoot->add_node(*AlignedTipsMapping());
      } else {
        currentNode = theRoot->go_down(1);
        aRoot->add_node(*AlignedTipsMapping());
        for (k = 2; k <= descendants; k++) {
          currentNode = theRoot->go_down(k);
          aChild->add_node(*AlignedTipsMapping());
        }
        aRoot->add_node(*aChild);
      }
      AlignNodes(aChild);
      AlignNodes(aRoot);
      return aRoot;
    }
  } else {
    node<nodeCoord> *aNode = new node<nodeCoord>;
    node<long> *saveCurrent = currentNode;
    descendants = saveCurrent->get_num_nodes();
    for (k = 1; k <= descendants; k++) {
      currentNode = saveCurrent->go_down(k);
      aNode->add_node(*AlignedTipsMapping());
    }
    if (!descendants) { // terminal node
      aNode->in_object.v = treeLayoutVert;
      aNode->in_object.h = 0;
      treeLayoutVert += TREE_V_SHIFT;
    } else {
      AlignNodes(aNode);
    }
    aNode->in_object.varRef = saveCurrent->in_object;
    currentNode = saveCurrent;
    return aNode;
  }
}

//______________________________________________________________________________
void _TheTree::ScaledBranchReMapping(node<nodeCoord> *theNode, _Parameter tw) {
  theNode->in_object.h -= tw;

  long descendants = theNode->get_num_nodes();

  for (long k = 1; k <= descendants; k++) {
    ScaledBranchReMapping(theNode->go_down(k), tw);
  }
}

//______________________________________________________________________________
_Parameter _TheTree::DetermineBranchLengthGivenScalingParameter(
    long varRef, _String &matchString, char mapMode) {

  if (mapMode == 3) {
    return 1.;
  }

  _CalcNode *travNode = (_CalcNode *)LocateVar(varRef);
  _Parameter branchLength = BAD_BRANCH_LENGTH;

  if (mapMode == 1) {
    return travNode->BranchLength();
  } else if (mapMode == 2) {
    branchLength = travNode->Value();
    if (branchLength <= 0.0) {
      branchLength = BAD_BRANCH_LENGTH;
    }
  } else {
    long j;
    if (travNode->iVariables)
      for (j = 0; j < travNode->iVariables->lLength; j += 2) {
        _Variable *curVar = LocateVar(travNode->iVariables->lData[j]);
        if (curVar->GetName()->endswith(matchString)) {
          branchLength = curVar->Compute()->Value();
          if (branchLength <= 0.0) {
            branchLength = BAD_BRANCH_LENGTH;
          } else {
            break;
          }
        }
      }

    if (((!travNode->iVariables) || j == travNode->iVariables->lLength) &&
        travNode->dVariables)
      for (j = 0; j < travNode->dVariables->lLength; j += 2) {
        _Variable *curVar = LocateVar(travNode->dVariables->lData[j]);
        if (curVar->GetName()->endswith(matchString)) {
          branchLength = curVar->Compute()->Value();
          if (branchLength <= 0.0) {
            branchLength = BAD_BRANCH_LENGTH;
          } else {
            break;
          }
        }
      }
  }

  return branchLength;
}

//______________________________________________________________________________
node<nodeCoord> *_TheTree::ScaledBranchMapping(node<nodeCoord> *theParent,
                                               _String *scalingParameter,
                                               long locDepth, long &depth,
                                               char mapMode) {

  // run a pass of aligned tip mapping then perform one more pass from the root
  // to the children
  // pre-order to remap the length of branches.
  static _Parameter treeWidth;
  bool wasRoot = !theParent;

  if (!theParent) {
    theParent = AlignedTipsMapping(true, true);
    theParent->in_object.h = 0.0;
    treeWidth = 0;
  }

  node<nodeCoord> *currentN;
  long descendants = theParent->get_num_nodes(), k = 1, j, b = -1;

  _Parameter branchLength = BAD_BRANCH_LENGTH;

  for (; k <= descendants; k++) {
    currentN = theParent->go_down(k);
    j = currentN->in_object.varRef;

    if (j >= 0) {

      branchLength = currentN->in_object.bL =
          DetermineBranchLengthGivenScalingParameter(j, *scalingParameter,
                                                     mapMode);
      branchLength += theParent->in_object.h;

      if (branchLength > treeWidth) {
        treeWidth = branchLength;
      }

      theParent->go_down(k)->in_object.h = branchLength;
      ScaledBranchMapping(theParent->go_down(k), scalingParameter, locDepth + 1,
                          depth, mapMode);

    } else {
      theParent->go_down(k)->in_object.h = 0;
      b = k;
    }

  }

  if (k == descendants + 1) {
    if (locDepth >= depth) {
      depth = locDepth + 1;
    }
  }

  if (wasRoot) {
    if (b > 0 && descendants == 2) {
      j = (b == 1) ? 2 : 1;

      //branchLength = theParent->go_down(j)->in_object.h/2;
      //treeWidth -= branchLength;
      //branchLength = theParent->go_down(j)->in_object.h;
      ScaledBranchReMapping(theParent->go_down(j), 0);
      theParent->go_down(b)->in_object.h = 0;
      ScaledBranchMapping(theParent->go_down(b), scalingParameter, locDepth,
                          depth, mapMode);
    }

    ScaledBranchReMapping(theParent, treeWidth);
    return theParent;
  }
  return nil;
}

//______________________________________________________________________________
node<nodeCoord> *_TheTree::RadialBranchMapping(
    node<long> *referenceNode, node<nodeCoord> *parentNode,
    _String *scalingParameter, _Parameter anglePerTip, long &currentTipID,
    _Parameter &maxRadius, char mapMode) {

  // label 1 stores current radial distance from the root
  // label 2 stores the angle of the line to this node
  // h and v store the Cartesian coordinates
  node<nodeCoord> *current_node = new node<nodeCoord>;

  _Parameter branchL = 0., referenceL = 0.;

  if (parentNode == nil) {
    current_node->in_object.label1 = 0.0;
    current_node->in_object.label2 = 0.0;
  } else {
    referenceL = parentNode->in_object.label1;
    branchL = DetermineBranchLengthGivenScalingParameter(
        referenceNode->in_object, *scalingParameter, mapMode);
  }

  long children = referenceNode->get_num_nodes();

  current_node->in_object.label1 = referenceL + branchL;
  if (children == 0) {
    current_node->in_object.label2 = anglePerTip * currentTipID++;
    //printf ("%d %g\n",currentTipID, current_node->in_object.label2);
  } else {
    _Parameter angleSum = 0.;
    for (long n = 1; n <= children; n++) {
      node<nodeCoord> *newChild = RadialBranchMapping(
          referenceNode->go_down(n), current_node, scalingParameter,
          anglePerTip, currentTipID, maxRadius, mapMode);
      current_node->add_node(*newChild);
      angleSum += newChild->in_object.label2;
    }
    current_node->in_object.label2 = angleSum / children;
  }

  current_node->in_object.h =
      current_node->in_object.label1 * cos(current_node->in_object.label2);
  current_node->in_object.v =
      current_node->in_object.label1 * sin(current_node->in_object.label2);
  maxRadius = MAX(maxRadius, current_node->in_object.label1);
  current_node->in_object.varRef = referenceNode->in_object;
  current_node->in_object.bL = branchL;

  return current_node;
}

//______________________________________________________________________________
void _TheTree::AssignLabelsToBranches(node<nodeCoord> *theParent,
                                      _String *scalingParameter, bool below) {

  bool wasRoot = !(theParent->parent);

  node<nodeCoord> *currentN;
  long descendants = theParent->get_num_nodes(), k = 1, j, b = -1;

  _Parameter branchLength = BAD_BRANCH_LENGTH;

  char mapMode;
  _String matchString =
      DetermineBranchLengthMappingMode(scalingParameter, mapMode);

  for (; k <= descendants; k++) {
    currentN = theParent->go_down(k);
    j = currentN->in_object.varRef;
    if (j >= 0) {

      branchLength =
          DetermineBranchLengthGivenScalingParameter(j, matchString, mapMode);
      if (below) {
        currentN->in_object.label2 = branchLength;
      } else {
        currentN->in_object.label1 = branchLength;
      }

      AssignLabelsToBranches(theParent->go_down(k), scalingParameter, below);

    } else {
      if (below) {
        currentN->in_object.label2 = 0.;
      } else {
        currentN->in_object.label1 = 0.;
      }
      b = k;
      AssignLabelsToBranches(theParent->go_down(k), scalingParameter, below);
    }

  }

  if (wasRoot) {
    if ((b > 0) && (descendants == 2)) {
      if (b == 1) {
        j = 2;
      } else {
        j = 1;
      }
      if (below) {
        theParent->in_object.label2 =
            theParent->go_down(j)->in_object.label2 / 2.;
        theParent->go_down(j)->in_object.label2 /= 2.;
      } else {
        theParent->in_object.label1 =
            theParent->go_down(j)->in_object.label1 / 2.;
        theParent->go_down(j)->in_object.label1 /= 2.;
      }
    }
  }
}

//______________________________________________________________________________
_Parameter _TheTree::PSStringWidth(_String &s) {
  _Parameter nnWidth = 0.;
  for (long cc = 0; cc < s.sLength; cc++) {
    nnWidth += _timesCharWidths[(int) s.getChar(cc)];
  }
  return nnWidth;
}

//______________________________________________________________________________
_PMathObj _TheTree::PlainTreeString(_PMathObj p, _PMathObj p2) {
  _String *res = new _String((unsigned long) 512, true);
  long treeLayout = 0;

  if (p && (p->ObjectClass() == STRING)) {
    if (p2 && (p2->ObjectClass() == MATRIX)) {
      node<nodeCoord> *newRoot, *currentNd;

      bool doEmbed = false;
      bool doSymbol = false;
      _FString *extraPS = nil, *prefixPS = nil;
      long symbolsize = 3;

      _AssociativeList *toptions =
          (_AssociativeList *)FetchObjectFromVariableByType(&treeOutputAVL,
                                                            ASSOCIATIVE_LIST);

      if (toptions) {
        _PMathObj lc = toptions->GetByKey(treeOutputLayout, NUMBER);
        if (lc) {
          treeLayout = lc->Value();
        }
        lc = toptions->GetByKey(treeOutputEmbed, NUMBER);
        if (lc) {
          doEmbed = lc->Value();
        }
        lc = toptions->GetByKey(treeOutputSymbols, NUMBER);
        if (lc) {
          doSymbol = lc->Value();
        }

        lc = toptions->GetByKey(treeOutputExtraPS, STRING);
        if (lc) {
          extraPS = (_FString *)lc->Compute();
        }

        lc = toptions->GetByKey(treeOutputPrefixPS, STRING);
        if (lc) {
          prefixPS = (_FString *)lc->Compute();
        }

        lc = toptions->GetByKey(treeOutputSymbolSize, NUMBER);
        if (lc) {
          symbolsize = lc->Value();
        }
      }

      _String *theParam = (_String *)p->toStr(), t;

      bool scaling = theParam->sLength, doLabelWidth = true;

      long tipCount = 0, fontSize = -1;

      _Matrix *dimMatrix = ((_Matrix *)p2->Compute());

      _Parameter hScale = 1.0, vScale = 1.0, labelWidth = 0.,
                 treeHeight = (*dimMatrix)(0, 1),
                 treeWidth = (*dimMatrix)(0, 0),
                 treeRotation =
                     (dimMatrix->GetVDim() > 2) ? (*dimMatrix)(0, 2) : 0.0,
                 treeRadius, treeArcAngle, totalTreeL = 0.,
                 mappedTreeHeight = 0.0, mappedTreeHeight2 = 1e+100;

      // letter size in points
      if (treeLayout == 1) {
        treeRadius = treeWidth;
        treeArcAngle = treeHeight;

        if (treeRadius <= 0.0) {
          treeRadius = 300.;
        }
        if (treeArcAngle <= 0.0 || treeArcAngle > 360.) {
          treeArcAngle = 360.0;
        }

        treeWidth = treeHeight = treeRadius * 2.;

        treeRotation = (treeRotation - floor(treeRotation / 360.) * 360.) /
                       DEGREES_PER_RADIAN;
      } else {
        if (treeHeight <= 0.0) {
          treeHeight = 792.;
        }
        if (treeWidth <= 0.0) {
          treeWidth = 576.;
        }
      }

      t = _String("%!\n%% PS file for the tree '") & *GetName() &
          "'.\n%% Generated by " & GetVersionString() & " on " & GetTimeStamp();
      (*res) << &t;
      if (treeLayout == 1) {
        (*res) << "% Radial layout\n";
        (*res) << "/righttext  {dup newpath 0 0 moveto false charpath "
                  "closepath pathbbox pop exch pop exch sub       4 -1 roll "
                  "exch sub 3 -1 roll newpath moveto show} def\n";
      }
      if (!doEmbed) {
        (*res) << "<< /PageSize [";

        (*res) << _String(treeWidth + 15);
            /*wayne added 15 to make trees sit inside the page */
        (*res) << ' ';
        (*res) << _String(treeHeight + 15);

        (*res) << "] >> setpagedevice\n";
      }

      long xtraChars = 0;

      if (toptions) {
        _Matrix *rgbColor =
            (_Matrix *)(toptions)->GetByKey(treeOutputBackground, MATRIX);
        if (rgbColor) {
          t = _String((*rgbColor)(0, 0)) & " " & _String((*rgbColor)(0, 1)) &
              " " & _String((*rgbColor)(0, 2)) & " setrgbcolor\nnewpath\n";
          (*res) << t;
          (*res) << "0 0 moveto\n";
          (*res) << "0 ";
          (*res) << _String(treeHeight);
          (*res) << " lineto\n";
          (*res) << _String(treeWidth);
          (*res) << ' ';
          (*res) << _String(treeHeight);
          (*res) << " lineto\n";
          (*res) << _String(treeWidth);
          (*res) << " 0 lineto\n";
          (*res) << "closepath\nfill\n0 0 0 setrgbcolor\n";
        }

        if (doSymbol) {
          /*add some symbol drawing postscript functions */
          (*res) << "/size ";
          (*res) << _String(symbolsize);
          (*res) << " def\n";
          (*res) << "/box { 0 -0.5 size mul rmoveto\n";
          (*res) << "1 size mul 0 size mul rlineto\n";
          (*res) << "0 size mul 1 size mul rlineto\n";
          (*res) << "-1 size mul 0 size mul rlineto\n";
          (*res) << "closepath\n";
          (*res) << "} def\n";

          (*res) << "/triangle { size size 0.5 mul rlineto 0 size -1 mul "
                    "rlineto closepath } def\n";

          (*res) << "/circle {currentpoint exch 0.5 size mul add exch 0.5 size "
                    "mul 180 540 arc\n";
          (*res) << "closepath\n";
          (*res) << "} def\n";

          (*res) << "/diamond { 0 -0.5 size mul rmoveto 0.5 size mul 0 rmoveto "
                    "45 rotate 0.707107 size mul 0 rlineto 0 size 0.707107 mul "
                    "rlineto -0.707107 size mul 0 rlineto -45 rotate  "
                    "closepath} def\n";
        }

        _Constant *fontSizeIn =
            (_Constant *)(toptions)->GetByKey(treeOutputFSPlaceH, NUMBER);
        if (fontSizeIn) {
          fontSize = fontSizeIn->Value();
        }

        fontSizeIn =
            (_Constant *)(toptions)->GetByKey(treeOutputRightMargin, NUMBER);
        if (fontSizeIn) {
          treeWidth = MAX(treeWidth / 4 + 1, treeWidth - fontSizeIn->Value());
          doLabelWidth = false;
        }

        if ((fontSizeIn = (_Constant *)(toptions)
                 ->GetByKey(treeOutputXtraMargin, NUMBER))) {
          xtraChars = fontSizeIn->Value();
        }
      }

      char mapMode = 0;
      _String scalingStringMatch;

      if (scaling) {
        scalingStringMatch =
            DetermineBranchLengthMappingMode(theParam, mapMode);
      }

      if (treeLayout == 1) {
        _PMathObj tempValue = TipCount();
        tipCount = tempValue->Value();
        DeleteObject(tempValue);
        long tipID = 0;
        hScale = 0.;
        newRoot =
            RadialBranchMapping(theRoot, nil, &scalingStringMatch,
                                treeArcAngle * pi_const / (180. * tipCount),
                                tipID, hScale, mapMode);
        totalTreeL = hScale;
      } else {
        if (scaling) {
          newRoot = ScaledBranchMapping(nil, &scalingStringMatch, 0, tipCount,
                                        mapMode);
        } else {
          newRoot = AlignedTipsMapping(true);
        }

        hScale = -treeWidth / newRoot->in_object.h;
      }

      currentNd = NodeTraverser(newRoot);

      while (currentNd) {
        if (currentNd->in_object.v > mappedTreeHeight) {
          mappedTreeHeight = currentNd->in_object.v;
        }
        if (currentNd->in_object.v < mappedTreeHeight2) {
          mappedTreeHeight2 = currentNd->in_object.v;
        }

        if (!currentNd->get_num_nodes()) {
          tipCount++;
        }

        currentNd = NodeTraverser((node<nodeCoord> *)nil);
      }

      // compute the size of the font, based on the spacing between branches.
      // 36 is max and 2 is min;

      if (fontSize < 0) {
        fontSize = treeHeight / (tipCount + 1 + (scaling > 0) * 1.5);

        if (fontSize > 36) {
          fontSize = 36;
        } else if (fontSize < 2) {
          fontSize = 2;
        }
      }

      // now recompute the width of the of the tree, including the labels
      // do not allow names to take up more that 1/2 of the width
      // try to keep full length but if the names are too wide,
      // shrink the size of the font

      currentNd = NodeTraverser(newRoot);

      _Parameter plotBounds[4];

      if (treeLayout == 1) {
        vScale = treeRadius / hScale;
        plotBounds[0] = plotBounds[1] = plotBounds[2] = plotBounds[3] = 0.;
      }

      while (currentNd) {
        if (currentNd->in_object.varRef >= 0) {
          _String nodeName(*LocateVar(currentNd->in_object.varRef)->GetName());
          nodeName.Trim(nodeName.FindBackwards('.', 0, -1) + 1, -1);

          _AssociativeList *nodeOptions = nil;
          if (toptions) {
            nodeOptions = (_AssociativeList *)toptions->GetByKey(
                nodeName, ASSOCIATIVE_LIST);
          }

          _PMathObj nodeLabel = nodeOptions ? nodeOptions->GetByKey(
                                                  treeOutputLabel, STRING)
                                            : nil,
                    nodeTLabel = nodeOptions ? nodeOptions->GetByKey(
                                                   treeOutputTLabel, STRING)
                                             : nil;

          _Parameter nnWidth = 0.0;

          if (nodeLabel) {
            nnWidth = 0.0;
          } else if (nodeTLabel) {
            nnWidth =
                1. +
                PSStringWidth(*((_FString *)nodeTLabel->Compute())->theString);
            //printf ("%g\n", nnWidth);
          } else if (currentNd->get_num_nodes() == 0) {
            nnWidth = 1. + PSStringWidth(nodeName);
          }

          nnWidth += _maxTimesCharWidth * xtraChars;
          nnWidth *= fontSize;

          if (treeLayout == 1) {
            currentNd->in_object.label2 += treeRotation;
            _Parameter chordLength = computeChordLength(
                           treeRadius, currentNd->in_object.label2, plotBounds),
                       overflow =
                           MAX(0., treeRadius +
                                       (nnWidth - chordLength) * hScale /
                                           MAX(BAD_BRANCH_LENGTH,
                                               currentNd->in_object.label1));

            //nnWidth + currentNd->in_object.label1 * vScale-chordLength);

            if (overflow > treeRadius * .5) { // shift too big
              chordLength -=
                  currentNd->in_object.label1 / (2. * hScale) * treeRadius;

              chordLength = chordLength / nnWidth;
              fontSize = fontSize * chordLength;

              if (fontSize < 2) {
                fontSize = 2;
              }

              nnWidth = treeRadius * .5;
            } else {
              nnWidth = overflow;
            }
          } else {
            if (nnWidth > treeWidth * .5) {
              fontSize = fontSize / (2. * nnWidth / treeWidth);
              if (fontSize < 2) {
                fontSize = 2;
              }
              nnWidth = treeWidth * .5;
            }
          }
          if (nnWidth > labelWidth) {
            labelWidth = nnWidth;
          }
        }
        currentNd = NodeTraverser((node<nodeCoord> *)nil);
      }

      if (!doLabelWidth) {
        labelWidth = 0;
      }

      if (scaling) {
        treeHeight -= 3.0 * fontSize;
      }

      if (treeLayout == 1) {
        hScale = (treeRadius - labelWidth) / hScale;
        vScale = 0.;
      } else {
        hScale = -(treeWidth - labelWidth) / newRoot->in_object.h;
        vScale = (treeHeight - 2 * fontSize) /
                 (mappedTreeHeight - mappedTreeHeight2);
      }
      t = _String("/Times-Roman findfont\n") & fontSize &
          " scalefont\nsetfont\n";
      (*res) << &t;

      nodeCoord dummy;

      long lw = fontSize / 6 + 1;

      (*res) << _String(lw);
      (*res) << " setlinewidth\n1 setlinecap\n";

      if (prefixPS) {
        (*res) << prefixPS->theString->Replace(treeOutputFSPlaceH,
                                               _String(fontSize), true);
      }

      if (treeLayout == 1) {
        //newRoot->in_object.h = -plotBounds[1];
        //newRoot->in_object.v = -plotBounds[3];
        //hScale                 *= 2.*treeRadius /
        //MAX(plotBounds[0]-plotBounds[1],plotBounds[2]-plotBounds[3]);
        newRoot->in_object.h = treeRadius;
        newRoot->in_object.v = treeRadius;
        vScale = 0.;

        TreePSRecurse(newRoot, (*res), hScale, vScale, ceil(treeWidth),
                      ceil(treeHeight), fontSize / 2, labelWidth - fontSize / 2,
                      toptions, 1, &treeRadius);
      } else {
        TreePSRecurse(newRoot, (*res), hScale, vScale, ceil(treeWidth),
                      ceil(treeHeight), fontSize / 2, labelWidth - fontSize / 2,
                      toptions);
      }

      if (scaling) { /* ruler */
        if (fontSize < 8) {
          fontSize = 8;
          (*res) << (_String("/Times-Roman findfont\n") & fontSize &
                     " scalefont\nsetfont\n");
        }

        if (treeWidth < 4 * fontSize) { 
          // enforce minimal ruler width
          treeWidth = 4 * fontSize;
        }

        vScale = exp(log(10.) * floor(log10(treeLayout == 1 ? totalTreeL : treeWidth / hScale))); 

        // compute the scaling factor
        if (vScale * hScale < 10.0) {
          vScale *= 10.;
        }

        _String rulerLabel(vScale, "%5.2g");

        while (vScale * hScale >
               (treeLayout == 1 ? treeRadius / 3
                                : treeWidth - newRoot->in_object.h)) {
          vScale *= 0.5;
          rulerLabel = vScale;
        }

        while (PSStringWidth(rulerLabel) * fontSize > vScale * hScale - 3) {
          vScale *= 2.0;
          rulerLabel = vScale;
        }

        long lm = newRoot->in_object.h, rm = lw + vScale * hScale;

        if (treeLayout == 1) {
          treeHeight = 2 * treeRadius - 2 * fontSize;
          lm = treeWidth - fontSize - rm;
          rm = treeWidth - fontSize - lw;
        }

        (*res) << "newpath\n";
        t = _String(lm) & ' ' & _String((long) treeHeight + 2 * lw) &
            " moveto\n";
        (*res) << &t;
        t = _String(rm) & ' ' & (long)(treeHeight + 2 * lw) & " lineto\nstroke\nnewpath\n"; // main horizontal bar
        (*res) << &t;
        t = _String(lm) & ' ' & _String((long) treeHeight) & " moveto\n";
        (*res) << &t;
        t = _String(lm) & ' ' & (long)(treeHeight + 4 * lw) & " lineto\nstroke\nnewpath\n"; // left notch
        (*res) << &t;
        t = _String(rm) & ' ' & (long)(treeHeight) & " moveto\n";
        (*res) << &t;
        t = _String(rm) & ' ' & (long)(treeHeight + 4 * lw) & " lineto\nstroke\nnewpath\n"; // right notch
        (*res) << &t;
        t = _String(lm + lw) & ' ' & _String((long) treeHeight + 3 * lw) & " moveto\n"; // main text
        (*res) << &t;
        t = _String('(') & _String(rulerLabel) & ") show\n";
        (*res) << &t;
      }

      newRoot->delete_tree();
      delete newRoot;

      if (extraPS) {
        (*res) << extraPS->theString->Replace(treeOutputFSPlaceH,
                                              _String(fontSize), true);
      }

      if (!doEmbed) {
        t = "showpage";
        (*res) << &t;
      }
      DeleteObject(theParam);
    } else {
      _String errMsg("An invalid 3rd parameter was passed to PSTreeString.");
      ReportWarning(errMsg);
    }
  } else {
    _String errMsg("An invalid 2nd parameter was passed to PSTreeString.");
    ReportWarning(errMsg);
  }

  res->Finalize();
  return new _FString(res);
}

//______________________________________________________________________________
void _TheTree::BuildINodeDependancies(void) {
  _CalcNode *travNode = DepthWiseTraversal(true);
  leftiNodes.Clear();
  long iNodeCounter = 0, leafCounter = 0;

  topLevelNodes.Clear();

  while (travNode) {
    if (IsCurrentNodeATip()) {
      leftiNodes << iNodeCounter;
      leafCounter++;
    } else {
      iNodeCounter++;
    }
    travNode = DepthWiseTraversal();
  }
}

//______________________________________________________________________________
void _TheTree::BuildTopLevelCache(void) {
  _CalcNode *travNode = DepthWiseTraversal(true);
  long iNodeCounter = 0, leafCounter = 0;

  node<long> *np, *np2;

  topLevelNodes.Clear();
  topLevelLeftL.Clear();
  topLevelRightL.Clear();

  // use cBase to count the number of leaves at or below a given node
  // use categoryIndexVars to store left and right leaves

  while (travNode) {
    if (IsCurrentNodeATip()) {
      travNode->categoryIndexVars << leafCounter;
      travNode->categoryIndexVars << leafCounter;
      leafCounter++;
      travNode->cBase = 1;
    } else {
      travNode->cBase = 0;
      for (long k = 0; k < currentNode->nodes.length; k++) {
        _CalcNode *cNode =
            (_CalcNode *)LocateVar(currentNode->nodes.data[k]->in_object);
        if (k == 0) {
          travNode->categoryIndexVars
              << cNode->categoryIndexVars[cNode->categoryIndexVars.lLength - 2];
        }
        if (k == currentNode->nodes.length - 1) {
          travNode->categoryIndexVars
              << cNode->categoryIndexVars[cNode->categoryIndexVars.lLength - 1];
        }

        travNode->cBase += cNode->cBase;
      }
      travNode->lastState = iNodeCounter;
      iNodeCounter++;
    }
    travNode = DepthWiseTraversal();
  }

  for (long level = 0; level < theRoot->nodes.length; level++) {
    np = theRoot->nodes.data[level];
    travNode = (_CalcNode *)LocateVar(np->in_object);
    if (travNode->cBase > 1) {
      topLevelNodes << travNode->lastState;
      topLevelLeftL
          << travNode
                 ->categoryIndexVars[travNode->categoryIndexVars.lLength - 2];
      topLevelRightL
          << travNode
                 ->categoryIndexVars[travNode->categoryIndexVars.lLength - 1];
      if (travNode->cBase > 4 * leafCounter / 5) {
        // one i-node hogging all the descendants
        _SimpleList sndLevel;
        for (long k = 0; k < np->nodes.length; k++) {
          np2 = np->nodes.data[k];
          if (np2->nodes.length) {
            sndLevel << (long) np2;
          }
        }
        if (sndLevel.lLength > 1) {
          topLevelLeftL.Delete(topLevelNodes.lLength - 1);
          topLevelRightL.Delete(topLevelNodes.lLength - 1);
          topLevelNodes.Delete(topLevelNodes.lLength - 1);
          for (long k = 0; k < sndLevel.lLength; k++) {
            travNode = (_CalcNode *)LocateVar(((node<long> *)sndLevel.lData[k])->in_object);

            topLevelNodes << travNode->lastState;
            topLevelLeftL << travNode->categoryIndexVars[travNode->categoryIndexVars.lLength - 2];
            topLevelRightL << travNode->categoryIndexVars[travNode->categoryIndexVars.lLength - 1];
          }
          break;
        }
      }
    }
  }

  // restore the settings of CalcNodes.
  travNode = DepthWiseTraversal(true);
  while (travNode) {
    if (!IsCurrentNodeATip()) {
      travNode->cBase = cBase;
      travNode->lastState = -1;
    }
    long k = travNode->categoryIndexVars.lLength - 2;
    travNode->categoryIndexVars.Delete(k);
    travNode->categoryIndexVars.Delete(k, true);
    travNode = DepthWiseTraversal();
  }

  if (topLevelNodes.lLength) {
    topLevelNodes << 0;
    topLevelLeftL << leafCounter;
    topLevelRightL << leafCounter - 1;
  }
}

//______________________________________________________________________________
void _TheTree::KillTopLevelCache(void) {
  topLevelNodes.Clear();
  if (rootIChildrenCache) {
    free(rootIChildrenCache);
  }
  rootIChildrenCache = nil;
}

//______________________________________________________________________________
long _TheTree::CountTreeCategories(void) {

  categoryVariables.Clear();
  {
    _AVLList cVA(&categoryVariables);
    ScanForCVariables(cVA);
    cVA.ReorderList();
  }
  categoryCount = 1;
  for (unsigned long k = 0; k < categoryVariables.lLength; k++) {
    categoryCount *=
        ((_CategoryVariable *)LocateVar(categoryVariables.lData[k]))
            ->GetNumberOfIntervals();
  }
  return categoryCount;
}

//______________________________________________________________________________
void _TheTree::AllocateResultsCache(long size) {
  if (rootIChildrenCache) {
    free(rootIChildrenCache);
  }
  rootIChildrenCache = nil;
  //topLevelNodes.Clear();
  size *= categoryCount;
  if (topLevelNodes.lLength) {
    rootIChildrenCache = (_Parameter *)MemAllocate(
        size * cBase * (topLevelNodes.lLength - 1) * sizeof(_Parameter));
  }
}

//______________________________________________________________________________
nodeCoord _TheTree::TreeTEXRecurse(node<nodeCoord> *currNode, _String &res,
                                   _Parameter hScale, _Parameter vScale,
                                   long hSize, long vSize) {

  long descendants = currNode->get_num_nodes(), vc, hc, hc1, hc2;

  _String t;

  if (descendants == 0) { // terminal node
    vc = vSize - currNode->in_object.v * vScale;
    hc = hSize + currNode->in_object.h * hScale;
    t = _String("\n\\put(") & hc & ',' & vc & "){\\circle*{2}}";
    res << &t;
    t = _String("\n\\put(") & hc + 2 & ',' & vc - 1 & "){\\makebox{\\tiny{";
    res << &t;
    t = *(LocateVar(currNode->in_object.varRef)->GetName());
    t = t.Cut(t.Find('.') + 1, -1);
    res << &t;
    res << '}';
    res << '}';
    res << '}';
  } else {
    vc = vSize - currNode->in_object.v * vScale;
    hc = hSize + currNode->in_object.h * hScale;
    for (long k = 1; k <= descendants; k++) {
      node<nodeCoord> *childNode = currNode->go_down(k);
      TreeTEXRecurse(childNode, res, hScale, vScale, hSize, vSize);

      long vcc = vSize - childNode->in_object.v * vScale,
           hcc = hSize + childNode->in_object.h * hScale;

      t = _String("\n\\put(") & hc & ',' & vcc & "){\\line(1,0){" & (hcc - hc) &
          "}}";
      res << &t;
      if (k == 1) {
        hc1 = vcc;
      } else if (k == descendants) {
        hc2 = vcc;
      }
    }
    t = _String("\n\\put(") & hc & ',' & hc2 & "){\\line(0,1){" & (hc1 - hc2) &
        "}}";
    res << &t;
    t = _String("\n\\put(") & hc & ',' & vc & "){\\circle{2}}";
    res << &t;
    if (currNode->get_parent()) {
      t = _String("\n\\put(") & hc + 2 & ',' & vc - 1 & "){\\makebox{\\tiny{";
      res << &t;
      t = *(LocateVar(currNode->in_object.varRef)->GetName());
      t = t.Cut(t.Find('.') + 1, -1);
      if (t.beginswith("Node")) {
        t = t.Cut(4, -1);
      }
      res << &t;
      res << '}';
      res << '}';
      res << '}';
    }
  }

  nodeCoord resC;
  resC.h = hc;
  resC.v = vc;
  return resC;
}

//______________________________________________________________________________
void _TheTree::TreePSRecurse(node<nodeCoord> *currNode, _String &res,
                             _Parameter hScale, _Parameter vScale, long hSize,
                             long vSize, long halfFontSize, long shift,
                             _AssociativeList *outOptions, char layout,
                             _Parameter *xtra) {

  long descendants = currNode->get_num_nodes();
  //lineW    = halfFontSize/3+1;

  _Parameter vc, hc, vcl, hcl, hc1, hc2;

  _String t, varName, colorString("0 0 0 setrgbcolor\n");

  if (currNode->parent) {
    if (currNode->in_object.varRef >= 0) {
      varName = (*LocateVar(currNode->in_object.varRef)->GetName());
    }
  } else {
    hc = GetRoot().in_object;
    if (hc >= 0) {
      varName = (*LocateVar(hc)->GetName());
    }
    if (layout == 1) {
      res << (_String(currNode->in_object.h) & ' ' &
              _String(currNode->in_object.v) & " translate\n");
    }
  }

  varName.Trim(varName.Find('.') + 1, -1);

  _AssociativeList *nodeOptions = nil;
  if (outOptions) {
    nodeOptions =
        (_AssociativeList *)outOptions->GetByKey(varName, ASSOCIATIVE_LIST);
  }

  _PMathObj nodeLabel =
                nodeOptions ? nodeOptions->GetByKey(treeOutputLabel, STRING)
                            : nil,
            nodeTLabel =
                nodeOptions ? nodeOptions->GetByKey(treeOutputTLabel, STRING)
                            : nil;

  if (layout == 1) {
    hcl = currNode->in_object.h * hScale;
    vcl = currNode->in_object.v * hScale;
  } else {
    vcl = vSize - currNode->in_object.v * vScale,
    hcl = hSize + currNode->in_object.h * hScale - shift;
  }

  if (descendants == 0 || nodeLabel) {
    // terminal node or default label
    t = empty;
    _Parameter myAngle =
        layout == 1 ? currNode->in_object.label2 * DEGREES_PER_RADIAN : 0.0;
    if (layout == 1) {
      res << (_String(myAngle) & " rotate\n");
      vc = 0;
      hcl = (vScale + currNode->in_object.bL) * hScale;
      hc = hcl + halfFontSize;
    } else {
      vc = vcl - halfFontSize;
      hc = hcl + halfFontSize;
    }

    if (!nodeTLabel) {
      res << "newpath\n";

      if (nodeLabel) {
        t = _String(hcl) & ' ' & _String(vc) & " moveto\n";

        res << &t;
        t = *((_FString *)nodeLabel->Compute())->theString;
        t = t.Replace(treeOutputNNPlaceH, varName, true);
        t = t.Replace(treeOutputFSPlaceH, _String(halfFontSize * 2), true) &
            '\n';
      }

      if (descendants == 0 && t.sLength == 0) {
        // generate the default label
        if (layout == 1 && myAngle > 90. && myAngle < 270.) {
          _Parameter xt = hc - halfFontSize / 2, yt = vc - 2 * halfFontSize / 3;
          t = _String(xt) & _String(" 0 translate 180 rotate 0 ") &
              _String(yt) & _String('(') & varName &
              ") righttext 180 rotate -" & xt & " 0 translate\n";
        } else {
          t = _String(hc - halfFontSize / 2) & ' ' &
              _String(vc - 2 * halfFontSize / 3) & " moveto\n";
          res << &t;
          t = _String('(') & varName & ") show\n";
        }
      }

      res << &t;
    }

    if (descendants == 0) {
      currNode->in_object.h = hc - halfFontSize;
    }

    if (layout == 1) {
      res << (_String(-currNode->in_object.label2 * DEGREES_PER_RADIAN) &
              " rotate\n");
    }
  }

  long minChildHC = 0x0fffffff, newV = 0;

  if (descendants) {
    vc = vSize - currNode->in_object.v * vScale;
    hc = hSize + currNode->in_object.h * hScale - shift;

    nodeCoord childCoord;
    for (long k = 1; k <= descendants; k++) {
      node<nodeCoord> *child = currNode->go_down(k);
      TreePSRecurse(child, res, hScale,
                    (layout == 1) ? vScale + currNode->in_object.bL : vScale,
                    hSize, vSize, halfFontSize, shift, outOptions, layout,
                    xtra);
      if (k == 1) {
        hc1 = layout == 1 ? child->in_object.label2 : child->in_object.v;
      }
      if (k == descendants) {
        hc2 = layout == 1 ? child->in_object.label2 : child->in_object.v;
      }
    }

    char doVLines = 3;

    for (long k = 1; k <= descendants; k++) {
      node<nodeCoord> *child = currNode->go_down(k);

      if (child->in_object.varRef >= 0) {
        t = *LocateVar(child->in_object.varRef)->GetName();
        t.Trim(t.Find('.') + 1, -1);
      } else {
        t = empty;
      }

      newV += child->in_object.v;

      _AssociativeList *childOptions = nil;

      _Parameter splitBranch = -1., lineWP = 0.0;

      _String childColor, notchColor, *blabelString = nil, childDash,
                                      linewidth1, linewidth2, linecap1,
                                      linecap2;

      _Matrix *notches = nil, *multiColor = nil;

      if (layout == 1) {
        res << (_String(child->in_object.label2 * DEGREES_PER_RADIAN) &
                " rotate\n");
      }

      if (outOptions) {
        childOptions =
            (_AssociativeList *)outOptions->GetByKey(t, ASSOCIATIVE_LIST);
        if (childOptions) {
          _PMathObj keyVal =
              childOptions->GetByKey(treeOutputThickness, NUMBER);
          if (keyVal) {
            lineWP = keyVal->Compute()->Value();
            linewidth1 =
                _String("currentlinewidth ") & lineWP & " setlinewidth\n";
            linewidth2 = "setlinewidth\n";
          }
          keyVal = childOptions->GetByKey(treeOutputLinecap, NUMBER);
          if (keyVal) {
            linecap1 = _String("currentlinecap ") & (long)
                       keyVal->Compute()->Value() & " setlinecap\n";
            linecap2 = "setlinecap\n";
          }
          keyVal = childOptions->GetByKey(treeOutputSplit, NUMBER);
          if (keyVal) {
            splitBranch = keyVal->Compute()->Value();
          }

          keyVal = childOptions->GetByKey(treeOutputNotches, MATRIX);
          if (keyVal) {
            notches = (_Matrix *)(((_Matrix *)keyVal->Compute())
                                      ->ComputeNumeric())->makeDynamic();
          }

          keyVal = childOptions->GetByKey(treeOutputColor, MATRIX);
          if (keyVal) {
            _Matrix *rgbColor = (_Matrix *)keyVal->Compute();
            if (rgbColor->GetHDim() > 1 && rgbColor->GetVDim() == 4 &&
                layout != 1) {
              multiColor = (_Matrix *)rgbColor->makeDynamic();
            } else {
              childColor = _String((*rgbColor)(0, 0)) & " " &
                           _String((*rgbColor)(0, 1)) & " " &
                           _String((*rgbColor)(0, 2)) & " setrgbcolor\n";
            }
          }

          keyVal = childOptions->GetByKey(treeOutputNotchesColor, MATRIX);
          if (keyVal) {
            _Matrix *rgbColor = (_Matrix *)keyVal->Compute();
            notchColor =
                _String((*rgbColor)(0, 0)) & " " & _String((*rgbColor)(0, 1)) &
                " " & _String((*rgbColor)(0, 2)) & " setrgbcolor\n";
          }

          keyVal = childOptions->GetByKey(treeOutputDash, MATRIX);
          if (keyVal) {
            _Matrix *dash = (_Matrix *)keyVal->Compute();
            childDash = _String('[') & _String((*dash)(0, 0)) & " " &
                        _String((*dash)(0, 1)) & "] " & _String((*dash)(0, 2)) &
                        " setdash\n";
          }
          keyVal = childOptions->GetByKey(treeOutputOLabel, STRING);
          if (keyVal) {
            blabelString = ((_FString *)keyVal->Compute())->theString;
          }
        }
      }

      if (blabelString) {
        if (layout == 1) {
          t = _String(currNode->in_object.label1 * hScale) & " 0 moveto\n";
        } else {
          t = _String(hc) & ' ' & _String(child->in_object.v) & " moveto\n";
        }
        res << &t;
        *blabelString =
            blabelString->Replace(treeOutputNNPlaceH, varName, true).Replace(
                treeOutputFSPlaceH, _String(halfFontSize * 2), true) & '\n';

        res << blabelString;
        res << '\n';
      }

      res << &childDash;
      res << &childColor;
      res << &linewidth1;
      res << &linecap1;

      if (layout == 1) {
        if (child->get_num_nodes()) { // internal node
          res << "newpath\n";
          res << (_String("0 0 ") & child->in_object.label1 * hScale & ' ' &
                  child->in_object.v & ' ' & (child->in_object.auxD) &
                  " arc \n");
          res << "stroke\n";
        }
      } else {
        if (childOptions && descendants == 2) {
          res << "newpath\n";
          t = _String(hc) & ' ' & _String(0.5 * (hc1 + hc2)) & " moveto\n";
          res << &t;
          t = _String(hc) & ' ' & _String(k == 1 ? hc1 : hc2) & " lineto\n";
          res << &t;
          res << "stroke\n";
          doVLines -= k;
        }

      }

      res << "newpath\n";
      if (layout == 1) {
        t = _String(child->in_object.label1 * hScale) & " 0 moveto\n";
        res << &t;
        t = _String(-child->in_object.bL * hScale) & " 0 rlineto\n";
      } else {

        _Parameter lineWidthInset = 0.0;

        //if (lineWP > 0.0)
        //  lineWidthInset = (lineWP-lineW)*.5;

        if (multiColor) {
          _Parameter span = child->in_object.h - hc + 2 * lineWidthInset,
                     currentX = hc - lineWidthInset;

          res << (_String(currentX) & ' ' & _String(child->in_object.v) &
                  " moveto\n");
          for (long seg = 0; seg < multiColor->GetHDim(); seg++) {
            res << (_String((*multiColor)(seg, 0)) & " " &
                    _String((*multiColor)(seg, 1)) & " " &
                    _String((*multiColor)(seg, 2)) & " setrgbcolor\n");
            _Parameter mySpan = span * (*multiColor)(seg, 3);
            res << (_String(mySpan) & " 0 rlineto\n");
            if (seg < multiColor->GetHDim() - 1) {
              currentX += mySpan;
              res << "stroke\nnewpath\n";
              res << (_String(currentX) & ' ' & _String(child->in_object.v) &
                      " moveto\n");
            }
          }
          DeleteObject(multiColor);
          t = empty;
        } else {
          t = _String(hc - lineWidthInset) & ' ' & _String(child->in_object.v) &
              " moveto\n";
          res << &t;
          t = _String(child->in_object.h) & ' ' &
              _String(child->in_object.v + lineWidthInset) & " lineto\n";
        }
        if (layout == 1) {
          minChildHC = MIN(child->in_object.label1, minChildHC);
        } else {
          minChildHC = MIN(child->in_object.h, minChildHC);
        }
      }
      res << &t;
      res << "stroke\n";
      res << &linecap2;
      res << &linewidth2;

      if (childDash.sLength) {
        res << "[] 0 setdash\n";
      }

      if (splitBranch >= 0.0 && splitBranch <= 1.0) {
        res << "newpath\n";
        _Parameter x, y;

        if (layout == 1) {
          x = (child->in_object.label1 - child->in_object.bL * splitBranch) *
              hScale;
          y = 0.;
        } else {
          x = hc + (child->in_object.h - hc) * (1. - splitBranch);
          y = child->in_object.v;
        }

        t = _String(x) & ' ' & _String(y) & " " & halfFontSize & " 0 360 arc\n";
        res << &t;
        res << "fill\n";
      }

      if (notches) {
        notches->CheckIfSparseEnough(true);
        res << notchColor;
        for (long l = 0; l < notches->GetSize(); l++) {
          _Parameter aNotch = (*notches)[l];
          if (aNotch >= 0. && aNotch <= 1.) {
            res << "newpath\n";
            _Parameter x, y;

            if (layout == 1) {
              x = (child->in_object.label1 - child->in_object.bL * aNotch) *
                  hScale;
              y = 0.;
            } else {
              x = hc + (child->in_object.h - hc) * (1. - aNotch);
              y = child->in_object.v;
            }

            t = _String(x - 0.5 * halfFontSize) & ' ' &
                _String(y - 0.5 * halfFontSize) & " moveto ";
            res << &t;
            t = _String(x + 0.5 * halfFontSize) & ' ' &
                _String(y + 0.5 * halfFontSize) & " lineto\n";
            res << &t;
            t = _String(x - 0.5 * halfFontSize) & ' ' &
                _String(y + 0.5 * halfFontSize) & " moveto ";
            res << &t;
            t = _String(x + 0.5 * halfFontSize) & ' ' &
                _String(y - 0.5 * halfFontSize) & " lineto\n";
            res << &t;
            res << "stroke\n";
          }
        }
        DeleteObject(notches);
      }

      res << &colorString;
      if (layout == 1) {
        res << (_String(-child->in_object.label2 * DEGREES_PER_RADIAN) &
                " rotate\n");
      }
    }

    if (layout == 0 && doVLines) {
      _String linewidth1, linewidth2, linecap1, linecap2;

      if (nodeOptions) {
        _PMathObj keyVal = nodeOptions->GetByKey(treeOutputThickness, NUMBER);
        if (keyVal) {
          _Parameter lineWP = keyVal->Compute()->Value();
          linewidth1 =
              _String("currentlinewidth ") & lineWP & " setlinewidth\n";
          linewidth2 = "setlinewidth\n";
        }
        keyVal = nodeOptions->GetByKey(treeOutputLinecap, NUMBER);
        if (keyVal) {
          linecap1 = _String("currentlinecap ") & (long)
                     keyVal->Compute()->Value() & " setlinecap\n";
          linecap2 = "setlinecap\n";
        }
      }

      res << &linecap1;
      res << &linewidth1;
      res << "newpath\n";

      if (doVLines == 3) {
        t = _String(hc) & ' ' & _String(hc2) & " moveto\n";
      } else {
        t = _String(hc) & ' ' & _String(0.5 * (hc1 + hc2)) & " moveto\n";
      }
      res << &t;
      if (doVLines == 3) {
        t = _String(hc) & ' ' & _String(hc1) & " lineto\n";
      } else {
        t = _String(hc) & ' ' & _String(doVLines == 1 ? hc1 : hc2) &
            " lineto\n";
      }
      res << &t;
      res << "stroke\n";
      res << &linewidth2;
      res << &linecap2;
    }

    if (layout == 0) {
      currNode->in_object.h = hc;
      currNode->in_object.v = newV / descendants;
    } else {
      currNode->in_object.auxD =
          (hc2 - currNode->in_object.label2) * DEGREES_PER_RADIAN;
      if (currNode->parent) {
        currNode->in_object.v =
            (hc1 - currNode->in_object.label2) * DEGREES_PER_RADIAN;
      }
    }
  } else {
    currNode->in_object.v = vc;
  }

  if (nodeTLabel) {
    t = *((_FString *)nodeTLabel->Compute())->theString;
    if (t.sLength) {
      _Parameter nnWidth = 1. + PSStringWidth(t), scF = 2. * halfFontSize;

      if (layout == 1) {
        res << (_String(currNode->in_object.label2 * DEGREES_PER_RADIAN) &
                " rotate\n");
      } else {
        if (descendants) {
          nnWidth = (minChildHC - hc) / nnWidth;
          vcl = currNode->in_object.v;
          //if (nnWidth < 2.*halfFontSize)
          //scF = nnWidth;
          //else
          //nnWidth = 1.;
        }
      }

      if (scF != 2. * halfFontSize) {
        res << (_String("/Times-Roman findfont\n") & scF &
                " scalefont\nsetfont\n");
      }

      res << "newpath\n";
      if (layout == 1) {
        res << (_String(currNode->in_object.label1 * hScale + halfFontSize) &
                ' ' & _String(-2 * halfFontSize / 3) & " moveto\n");
      } else {
        if (descendants) {
          hc = hcl + scF * 0.5;
          vc = vcl - scF * 0.33;

        } else {
          hc -= scF * 0.25;
          vc -= scF * 0.33;
        }
        res << (_String(hc) & ' ' & _String(vc) & " moveto\n");
      }

      res << '(';
      res << t;
      res << ") show \n";

      if (scF != 2. * halfFontSize) {
        res << (_String("/Times-Roman findfont\n") & halfFontSize * 2 &
                " scalefont\nsetfont\n");
      }

      if (layout == 1) {
        res << (_String(-currNode->in_object.label2 * DEGREES_PER_RADIAN) &
                " rotate\n");
      }
    }
  }

  if (colorString.sLength) {
    res << "0 0 0 setrgbcolor\n";
  }

  if (currNode->parent == nil && layout == 1) {
    res << (_String(-currNode->in_object.h) & ' ' &
            _String(-currNode->in_object.v) & " translate\n");
  }
}

//______________________________________________________________________________
_PMathObj _TheTree::TEXTreeString(_PMathObj p) {
  _String *res = new _String((unsigned long) 10, true);
  if (p && (p->ObjectClass() == STRING)) {
    node<nodeCoord> *newRoot;
    _String *theParam = (_String *)p->toStr(), t;

    bool scaling = theParam->sLength;

    long tipCount = 0;

    node<nodeCoord> *currentNd;

    _Parameter hScale = 1.0, vScale = 1.0, treeHeight = 0.0, treeWidth;
    if (scaling) {
      newRoot = ScaledBranchMapping(nil, theParam, 0, tipCount, 0);

      treeWidth = tipCount * WIDTH_PER_BRANCH;

      if (treeWidth < MIN_TEX_WIDTH) {
        treeWidth = MIN_TEX_WIDTH;
      } else if (treeWidth > MAX_TEX_WIDTH) {
        treeWidth = MAX_TEX_WIDTH;
      }

      hScale = -treeWidth / newRoot->in_object.h;
    } else {
      newRoot = AlignedTipsMapping(true);
      treeWidth = -newRoot->in_object.h;
      if (treeWidth < MIN_TEX_WIDTH) {
        hScale = MIN_TEX_WIDTH / treeWidth;
        treeWidth = MIN_TEX_WIDTH;
      } else if (treeWidth > MAX_TEX_WIDTH) {
        hScale = treeWidth / MAX_TEX_WIDTH;
        treeWidth = MAX_TEX_WIDTH;
      }

    }
    currentNd = newRoot;

    tipCount = newRoot->get_num_nodes();

    while (tipCount) {
      currentNd = currentNd->go_down(1);
      tipCount = currentNd->get_num_nodes();
    }

    treeHeight = currentNd->in_object.v;
    tipCount = newRoot->get_num_nodes();
    currentNd = newRoot;

    while (tipCount) {
      currentNd = currentNd->go_down(tipCount);
      tipCount = currentNd->get_num_nodes();
    }

    treeHeight = currentNd->in_object.v - treeHeight;

    tipCount = 0;

    if (treeHeight < MIN_TEX_HEIGHT) {
      vScale = (MIN_TEX_HEIGHT) / treeHeight;
      treeHeight = MIN_TEX_HEIGHT;
    } else if (treeHeight > MAX_TEX_HEIGHT) {
      vScale = treeHeight / (MAX_TEX_HEIGHT);
      treeHeight = MAX_TEX_HEIGHT;
    }

    t = "\n\\setlength{\\unitlength}{1mm}\n\\begin{picture}(";
    (*res) << &t;
    t = (long)(treeWidth + 5);
    (*res) << &t;
    (*res) << ',';
    t = (long)(treeHeight + 5);
    (*res) << &t;
    (*res) << ')';

    TreeTEXRecurse(newRoot, (*res), hScale, vScale, ceil(treeWidth),
                   ceil(treeHeight));
    newRoot->delete_tree();
    delete newRoot;

    t = "\n\\end{picture}";
    (*res) << &t;

    DeleteObject(theParam);
  } else {
    _String errMsg("An invalid 2nd parameter was passed to TEXTreeString");
  }

  res->Finalize();
  return new _FString(res);
}

//______________________________________________________________________________
void _TheTree::SetUpMatrices(long categCount) {
  _CalcNode *travNode;
  categoryCount = categCount >= 1 ? categCount : 1;
  travNode = DepthWiseTraversal(TRUE);
  while (travNode) {
    if (travNode->IsConstant()) {
      travNode->varFlags |= HY_VC_NO_CHECK;
    }
    travNode->ConvertToSimpleMatrix();
    if (categoryCount == 1) {
      travNode->matrixCache = nil;
    } else {
      checkPointer(travNode->matrixCache = (_Matrix **)MemAllocate(
          categoryCount * sizeof(_Matrix *)));
      for (long i = 0; i < categoryCount; i++) {
        travNode->matrixCache[i] = nil;
      }
    }
    travNode = DepthWiseTraversal();
  }
}

//______________________________________________________________________________
#if USE_SCALING_TO_FIX_UNDERFLOW
void _TheTree::AllocateUnderflowScalers(long sites) {
  DeleteObject(scalingForUnderflow);
  scalingForUnderflow = new _Matrix(sites, 1, false, true);
  for (long k = 0; k < sites; k++) {
    scalingForUnderflow->theData[k] = 1.0; //genrand_real2()*2.0;
  }
}

//______________________________________________________________________________
void _TheTree::DeallocateUnderflowScalers(void) {
  DeleteObject(scalingForUnderflow);
  scalingForUnderflow = nil;
}
#endif

//______________________________________________________________________________
void _TheTree::CleanUpMatrices(void) {
  _CalcNode *travNode;
  travNode = DepthWiseTraversal(TRUE);
  if (categoryCount == 1) {
    while (travNode) {

      // mod 05/03/2003 - uncomment next 5 lines
      // this breaks after ReplicateConstraint or MolecularClock is called
      // WTF?

      travNode->ConvertFromSimpleMatrix();

      if (travNode->referenceNode >= 0) {
        travNode->SetRefNode(-1);
        travNode->compExp = nil;
      } else {
        if (travNode->referenceNode < -1) {
          travNode->SetRefNode(-1);
        }
      }
      if (travNode->compExp) {
        DeleteObject(travNode->compExp);
        travNode->compExp = nil;
      }

      travNode->varFlags &= HY_VC_CLR_NO_CHECK;
      travNode = DepthWiseTraversal();
    }
  } else {
    while (travNode) {
      travNode->ConvertFromSimpleMatrix();
      if (travNode->referenceNode >= 0) {
        travNode->SetRefNode(-1);
      } else
        for (long i = 0; i < categoryCount; i++) {
          DeleteObject(travNode->matrixCache[i]);
        }

      free(travNode->matrixCache);
      travNode->matrixCache = nil;
      travNode->compExp = nil;
      travNode->varFlags &= HY_VC_CLR_NO_CHECK;
      travNode = DepthWiseTraversal();
    }
    categoryCount = 1;
  }

}

//______________________________________________________________________________
void _TheTree::RemoveModel(void) {
  _CalcNode *travNode;
  travNode = DepthWiseTraversal(TRUE);
  while (travNode) {
    travNode->RemoveModel();
    travNode = DepthWiseTraversal();
  }
  categoryCount = 1;
}

//______________________________________________________________________________
bool _TheTree::FindScalingVariables(_SimpleList &rec) {
  long i;
  rec.Clear();
  _CalcNode *travNode;
  travNode = StepWiseTraversal(TRUE);
  travNode = StepWiseTraversal();
  if (travNode) {
    if (travNode->iVariables)
      for (i = 1; i < travNode->iVariables->lLength; i += 2)
        if (travNode->iVariables->lData[i] >= 0) {
          rec << travNode->iVariables->lData[i];
        }

    if (travNode->dVariables)
      for (i = 1; i < travNode->dVariables->lLength; i += 2)
        if (travNode->dVariables->lData[i] >= 0) {
          rec << travNode->dVariables->lData[i];
        }

  }
  if (rec.lLength == 0) {
    return false;
  }
  while (travNode) {
    for (i = 0; i < rec.countitems(); i++) {
      if (((!travNode->iVariables) ||
           travNode->iVariables->FindStepping(rec.lData[i], 2, 1) < 0) &&
          ((!travNode->dVariables) ||
           travNode->dVariables->FindStepping(rec.lData[i], 2, 1) < 0)) {
        rec.Delete(i);
        if (rec.lLength == 0) {
          break;
        }
        i--;
      }
    }

    if (!((travNode->iVariables && travNode->iVariables->lLength) ||
          (travNode->dVariables && travNode->dVariables->lLength) ||
          (travNode->gVariables && travNode->gVariables->lLength))) {
      rec.Clear();
      return false;
    }

    travNode = StepWiseTraversal();
  }
  return true;
}

//______________________________________________________________________________
bool _TheTree::HaveStringBranchLengths(void) {
  _CalcNode *travNode = DepthWiseTraversal(TRUE);
  while (travNode && !IsCurrentNodeTheRoot()) {
    if (travNode->Value() < -0.9) {
      return false;
    }
    travNode = DepthWiseTraversal();
  }
  return true;
}

//______________________________________________________________________________
void _TheTree::ScanForVariables(_AVLList &l, _AVLList &l2, _AVLListX *tagger,
                                long weight) {

  unsigned long traversal_order = 0;
  _CalcNode *curNode = DepthWiseTraversal(true);
  while (curNode) {
    curNode->ScanForVariables(
        l, l2, tagger,
        weight + flatNodes.lLength + flatLeaves.lLength - traversal_order);
    curNode = DepthWiseTraversal();
    traversal_order += 1;
  }

}

//______________________________________________________________________________
void _TheTree::ScanAndAttachVariables(void) {
  _CalcNode *curNode = DepthWiseTraversal(true);
  while (curNode) {
    curNode->ScanAndAttachVariables();
    curNode = DepthWiseTraversal();
  }
}

//______________________________________________________________________________
void _TheTree::ScanForDVariables(_AVLList &l, _AVLList &l2) {
  _CalcNode *curNode = DepthWiseTraversal(true);
  while (curNode) {
    curNode->ScanForDVariables(l, l2);
    curNode = DepthWiseTraversal();
  }
}

//______________________________________________________________________________
void _TheTree::ScanForGVariables(_AVLList &li, _AVLList &ld, _AVLListX *tagger,
                                 long weight) {

  _CalcNode *curNode = DepthWiseTraversal(true);
  _SimpleList cL;
  _AVLList cLL(&cL);

  while (curNode) {

    _Formula *explicitFormMExp = curNode->GetExplicitFormModel();
    _Matrix *modelM = explicitFormMExp ? nil : curNode->GetModelMatrix();

    if (explicitFormMExp && cLL.Find((BaseRef) explicitFormMExp) < 0 ||
        modelM && cLL.Find(modelM) < 0) {
      _SimpleList temp;
      {
        _AVLList tempA(&temp);
        if (modelM) {
          modelM->ScanForVariables(tempA, true);
        } else {
          explicitFormMExp->ScanFForVariables(tempA, true, false, true, true);
        }
        tempA.ReorderList();
      }
      for (unsigned long i = 0; i < temp.lLength; i++) {
        long p = temp.lData[i];
        _Variable *v = LocateVar(p);
        if (v && v->IsGlobal()) {
          if (v->IsIndependent()) {
            li.Insert((BaseRef) p);
            if (tagger) {
              tagger->UpdateValue((BaseRef) p, weight, 0);
            }
          } else {
            ld.Insert((BaseRef) p);
          }
        }
      }
      cLL.Insert(modelM ? (BaseRef) modelM : (BaseRef) explicitFormMExp);
    }
    curNode->ScanForGVariables(li, ld);
    curNode = DepthWiseTraversal();
  }
}

//______________________________________________________________________________
void _TheTree::ScanForCVariables(_AVLList &lcat) {
  _CalcNode *curNode = DepthWiseTraversal(true);
  while (curNode) {
    for (long i = curNode->categoryVariables.lLength - 1; i >= 0; i--) {
      lcat.Insert((BaseRef) curNode->categoryVariables(i));
    }

    curNode = DepthWiseTraversal();

  }
}

//______________________________________________________________________________
bool _TheTree::HasChanged(void) {
  _CalcNode *curNode = StepWiseTraversal(true);
  while (curNode) {
    if (curNode->HasChanged()) {
      return true;
    }
    curNode = StepWiseTraversal();
  }
  return false;
}

//______________________________________________________________________________
bool _TheTree::HasChanged2(void) {
  for (long k = 0; k < categoryVariables.lLength; k++)
    if (((_CategoryVariable *)LocateVar(categoryVariables.lData[k]))
            ->HaveParametersChanged()) {
      return true;
    }
  _CalcNode *curNode = StepWiseTraversal(true);
  while (curNode) {
    if (curNode->_VariableContainer::HasChanged()) {
      return true;
    }
    curNode = StepWiseTraversal();
  }
  return false;
}

//______________________________________________________________________________
_Parameter _TheTree::Probij(long i, long j, _CalcNode *childNode) {
  // assumed that ExpMatrix has already been called for that node
  // range checking is implicit in the call of matrix (.,.) function
  if (childNode) {
    if (!childNode->GetCompExp()) {
      childNode->RecomputeMatrix();
    }
    //      if(childNode->GetCompExp())
    return (*childNode->GetCompExp())(i, j);
  }
  return 0;
}

//______________________________________________________________________________
// this will take the  matrix of frequencies and
// 1) use its dimensions to initialize tree freq holders
// 2) place global frequencies into the tree holder for later use by the
// pruning algo
// must be called before any tree pruning computations are started

void _TheTree::InitializeTreeFrequencies(_Matrix *mx, bool setDim) {
  long vecDim = mx->GetHDim() * mx->GetVDim();
  // theModel = mx;
  if (setDim) {
    SetTreeCodeBase(vecDim);
  } else
    for (long i = 0; i < vecDim; i++) {
      theProbs[i] = mx->theData[i];
    }
}

//______________________________________________________________________________
// this will take the  matrix of frequencies and
// 1) use its dimensions to initialize tree freq holders
// 2) place global frequencies into the tree holder for later use by the
// pruning algo
// must be called before any tree pruning computations are started
void _TheTree::SetTreeCodeBase(long b) {
  SetCodeBase(b);
  if (marginalLikelihoodCache) {
    free(marginalLikelihoodCache);
    marginalLikelihoodCache = nil;
  }
  if (cBase > 0)
    marginalLikelihoodCache =
        (_Parameter *)MemAllocate((flatNodes.lLength + flatLeaves.lLength) *
                                  sizeof(_Parameter) * cBase * systemCPUCount);

  _CalcNode *travNode = StepWiseTraversal(TRUE);
  while (travNode) {
    travNode->SetCodeBase(b);
    travNode = StepWiseTraversal();
  }

}

//______________________________________________________________________________
long _TheTree::IsLinkedToALF(long &pid) {
  for (long lfID = 0; lfID < likeFuncList.lLength; lfID++)
    if (likeFuncList.lData[lfID] && (pid = ((_LikelihoodFunction *)likeFuncList(lfID))->DependOnTree(*GetName())) >= 0) {
      return lfID;
    }
  return -1;
}

//______________________________________________________________________________
// set leaf value and reexp if needed
// for first entry into a datafilter
_Parameter _TheTree::ReleafTreeAndCheck(_DataSetFilter *dsf, long index,
                                        bool cache, long categID) {

#if defined __MP__
  if (systemCPUCount > 1) {
    ThreadMatrixUpdate(categID, cache);
  } else
#endif
    SerialMatrixUpdate(categID, cache);

  if (cache) {
    MatrixCacheUpdate();
  }

  if (flatLeaves.lLength == 1) {
    return ReleafTreeDegenerate(dsf, index);
  }

  if (cache) {
    return ThreadReleafTreeCache(dsf, index, -1, 0, flatLeaves.lLength - 1,
                                 categID >= 0 ? categID : 0);
  } else {
    return ReleafTree(dsf, index, -1, 0, flatLeaves.lLength - 1);
  }
}

//______________________________________________________________________________
void _TheTree::ThreadMatrixUpdate(long categID, bool cache) {
#ifdef __MP__
  _CalcNode *travNode;
  node<long> *nodeChild;
  _SimpleList *taintedNodes = new _SimpleList;

  for (long nodeCount = 0; nodeCount < flatLeaves.lLength; nodeCount++) {
    travNode = (_CalcNode *)(((BaseRef *)flatCLeaves.lData)[nodeCount]);
    bool reexpnt = travNode->NeedToExponentiate(categID);
    if (reexpnt && travNode->GetModelMatrix()) {
      (*taintedNodes) << (long) travNode;
      if (cache) {
        travNode = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[
            ((node<long> *)(flatLeaves.lData[nodeCount]))->parent->in_object]);
        travNode->cBase = -1;
      }
    } else if (categID >= 0) {
      travNode->SetCompMatrix(categID);
    }

  }

  for (long nodeCount = 0; nodeCount < flatTree.lLength; nodeCount++) {
    travNode = (_CalcNode *)(((BaseRef *)flatTree.lData)[nodeCount]);
    bool reexpnt = travNode->NeedToExponentiate(categID);
    if (reexpnt && travNode->GetModelMatrix()) {
      (*taintedNodes) << (long) travNode;
      if (cache) {
        travNode->cBase = -1;
      }
    } else if (categID >= 0) {
      travNode->SetCompMatrix(categID);
    }

    if (cache & (travNode->cBase == -1)) {
      nodeChild = ((node<long> *)(flatNodes.lData[nodeCount]))->parent;
      if (nodeChild) {
        travNode = (
            (_CalcNode *)((BaseRef *)variablePtrs.lData)[nodeChild->in_object]);
        travNode->cBase = -1;
      }
    }
  }

  if ((*taintedNodes).lLength > 1) {
    long tStep = (*taintedNodes).lLength / systemCPUCount, threadCount = 0, k,
         errCode;
    if (tStep > 0) {
      threadCount = systemCPUCount - 1;
    } else {
      tStep = 1;
      threadCount = (*taintedNodes).lLength - 1;
    }
    matrixTasks = new ThreadMatrixTask[threadCount];
    matrixThreads = new pthread_t[threadCount];

    for (k = 0; k < threadCount; k++) {
      matrixTasks[k].cID = categID;
      matrixTasks[k].tcat = categoryCount;
      matrixTasks[k].startAt = tStep * (k + 1);
      matrixTasks[k].endAt = tStep * (k + 2);
      if (k == threadCount - 1) {
        matrixTasks[k].endAt = (*taintedNodes).lLength;
      }
      matrixTasks[k].updateCN = taintedNodes;

      if (pthread_create(matrixThreads + k, NULL, MatrixUpdateFunction,
                         (void *)(matrixTasks + k))) {
        FlagError(
            "Failed to initialize a POSIX thread in ReleafTreeAndCheck()");
        exit(1);
      }
    }

    for (k = 0; k < tStep; k++) {
      ((_CalcNode *)((*taintedNodes).lData[k]))->RecomputeMatrix(categID, categoryCount);
    }

    for (k = 0; k < threadCount; k++)
      if ((errCode = pthread_join(matrixThreads[k], NULL))) {
        FlagError(_String("Failed to join POSIX threads in "
                          "ReleafTreeAndCheck(). Error Code=") & errCode);
        exit(1);
      }

    delete matrixTasks;
    delete matrixThreads;
    matrixTasks = nil;
  } else if ((*taintedNodes).lLength == 1) {
    ((_CalcNode *)((*taintedNodes).lData[0]))
        ->RecomputeMatrix(categID, categoryCount);
  }

  delete taintedNodes;
#endif
}

//______________________________________________________________________________
void _TheTree::SerialMatrixUpdate(long categID, bool cache) {
  _CalcNode *travNode;
  node<long> *nodeChild;

  for (long nodeCount = 0; nodeCount < flatLeaves.lLength; nodeCount++) {
    travNode = (_CalcNode *)(((BaseRef *)flatCLeaves.lData)[nodeCount]);
    bool reexpnt = travNode->NeedToExponentiate(categID);
    // flag all the nodes up the chain - to the root

    if (reexpnt && travNode->GetModelMatrix()) {
      travNode->RecomputeMatrix(categID, categoryCount);
      if (cache) {
        travNode = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[
            ((node<long> *)(flatLeaves.lData[nodeCount]))->parent->in_object]);
        travNode->cBase = -1;
      }
    } else if (categID >= 0) {
      travNode->SetCompMatrix(categID);
    }
  }

  {
    //for BCC
    for (long nodeCount = 0; nodeCount < flatTree.lLength; nodeCount++) {
      travNode = (_CalcNode *)(((BaseRef *)flatTree.lData)[nodeCount]);
      bool reexpnt = travNode->NeedToExponentiate(categID);
      // flag all the nodes up the chain - to the root

      if (reexpnt && travNode->GetModelMatrix()) {
        travNode->RecomputeMatrix(categID, categoryCount);
        if (cache) {
          travNode->cBase = -1;
        }
      } else if (categID >= 0) {
        travNode->SetCompMatrix(categID);
      }

      if (cache & (travNode->cBase == -1)) {
        nodeChild = ((node<long> *)(flatNodes.lData[nodeCount]))->parent;
        if (nodeChild) {
          travNode = ((_CalcNode *)(
              (BaseRef *)variablePtrs.lData)[nodeChild->in_object]);
          travNode->cBase = -1;
        }
      }
    }
  }
}

//______________________________________________________________________________
void _TheTree::MatrixCacheUpdate(void) {
  long c = 0, off = 1;
  _CalcNode *travNode;
  for (long nodeCount = 0; nodeCount < topLevelNodes.lLength - 1;
       nodeCount++, off <<= 1) {
    travNode = (_CalcNode *)(
        ((BaseRef *)flatTree.lData)[topLevelNodes.lData[nodeCount]]);
    if (travNode->cBase <= -1) {
      c |= off;
    }
  }
  topLevelNodes.lData[topLevelNodes.lLength - 1] = c;
  for (long nodeCount2 = 0; nodeCount2 < flatTree.lLength; nodeCount2++) {
    travNode = (_CalcNode *)(((BaseRef *)flatTree.lData)[nodeCount2]);
    travNode->cBase = cBase;
  }
}

//______________________________________________________________________________
_Parameter _TheTree::ReleafTreeAndCheckChar4(_DataSetFilter *dsf, long index,
                                             bool cache, long categID) {

  long nodeCount = 0, f;
  _Parameter *mCache = marginalLikelihoodCache;

  if (dsf->IsNormalFilter()) {
    char *thisState = dsf->GetColumn(index);

    for (; nodeCount < flatLeaves.lLength; nodeCount++, mCache += 4) {
      _CalcNode *travNode = (_CalcNode *)(((BaseRef *)flatCLeaves.lData)[nodeCount]);
      f = dsf->theNodeMap.lData[nodeCount];
      long *cCache = dsf->conversionCache.lData + (thisState[f] - 40) * 5;
#if USE_SCALING_TO_FIX_UNDERFLOW
      _Parameter scalingFactor = scalingForUnderflow->theData[index];
      mCache[0] = travNode->theProbs[0] = *(cCache++) * scalingFactor;
      mCache[1] = travNode->theProbs[1] = *(cCache++) * scalingFactor;
      mCache[2] = travNode->theProbs[2] = *(cCache++) * scalingFactor;
      mCache[3] = travNode->theProbs[3] = *(cCache++) * scalingFactor;
#else
      mCache[0] = travNode->theProbs[0] = *(cCache++);
      mCache[1] = travNode->theProbs[1] = *(cCache++);
      mCache[2] = travNode->theProbs[2] = *(cCache++);
      mCache[3] = travNode->theProbs[3] = *(cCache++);
#endif
      nodeStates[nodeCount] = travNode->lastState = *cCache;
    }
  } else {
    _DataSetFilterNumeric *dsfN = (_DataSetFilterNumeric *)dsf;
    for (; nodeCount < flatLeaves.lLength; nodeCount++, mCache += 4) {
      _CalcNode *travNode =
          (_CalcNode *)(((BaseRef *)flatCLeaves.lData)[nodeCount]);
      //_Parameter*        pv = dsfN->getProbabilityVector
      //(dsf->theNodeMap.lData[nodeCount],index);
      _Parameter *pv =
          (categID < 0)
              ? (dsfN->probabilityVectors.theData +
                 dsf->theNodeMap.lData[nodeCount] * dsfN->shifter + index * 4)
              : (dsfN->categoryShifter * categID +
                 dsfN->probabilityVectors.theData +
                 dsf->theNodeMap.lData[nodeCount] * dsfN->shifter + index * 4);
      mCache[0] = travNode->theProbs[0] = pv[0];
      mCache[1] = travNode->theProbs[1] = pv[1];
      mCache[2] = travNode->theProbs[2] = pv[2];
      mCache[3] = travNode->theProbs[3] = pv[3];
      nodeStates[nodeCount] = travNode->lastState = -1;
    }
  }

  if (flatLeaves.lLength == 1) {
    _CalcNode *travNode = ((_CalcNode *)(
        (BaseRef *)variablePtrs.lData)[theRoot->nodes.data[0]->in_object]);

    bool reexpnt = travNode->NeedToExponentiate(categID);

    if (reexpnt && travNode->GetModelMatrix()) {
      travNode->RecomputeMatrix(categID, categoryCount);
    } else if (categID >= 0) {
      travNode->SetCompMatrix(categID);
    }
    return ReleafTreeChar4Degenerate(dsf, index);
  }
  if (cache) {
    PruneTreeChar4Cache(categID);
    return ThreadReleafTreeChar4(dsf, index, -1, 0, flatLeaves.lLength - 1,
                                 categID < 0 ? 0 : categID);
    //return res;
  }

  return PruneTreeChar4(categID);
}
//______________________________________________________________________________
// set leaf value and re-prune w/o reexping
// for subsequent entries into a datafilter
_Parameter _TheTree::ReleafTreeDegenerate(_DataSetFilter *dsf, long index) {

  _CalcNode *
  rt = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[theRoot->in_object]),
 *tip = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[theRoot->nodes.data[0]->in_object]);

  _Parameter reslt = 0.;

  // sum over one branch in the direction from the root to the single leaf
  long state1 =
           dsf->Translate2Frequencies((*dsf)(index, 0), rt->theProbs, true),
       state2 =
           dsf->Translate2Frequencies((*dsf)(index, 1), tip->theProbs, true);

  // now perform one loop
  _Parameter *fastIdx = tip->GetCompExp()->theData;
  //*nodeProbs;
  // 4 cases are possible:
  // 1-1
  // 1-many
  // many-1
  // many-many
  if (state1 >= 0 && state2 >= 0) { // 1-1
    reslt = theProbs[state1] * fastIdx[state1 * cBase + state2];
  } else if (state1 >= 0) { // 1-many
    _Parameter tmp = 0.;

    fastIdx += state1 * cBase;

    for (long i = 0; i < cBase; i++) {
      tmp += fastIdx[i] * tip->theProbs[i];
    }

    reslt = theProbs[state1] * tmp;
  } else if (state2 >= 0) { // many to 1
    fastIdx += state2;
    for (long i = 0; i < cBase; i++, fastIdx += cBase) {
      reslt += rt->theProbs[i] * *fastIdx * theProbs[i];
    }

  } else // many to many
    for (long i = 0; i < cBase; i++) {
      _Parameter tmp = 0.0;
      for (long j = 0; j < cBase; j++, fastIdx++) {
        tmp += *fastIdx * tip->theProbs[j];
      }

      reslt += tmp * rt->theProbs[i] * theProbs[i];
    }

  return reslt <= 0.0 ? ALMOST_ZERO : reslt;

  /*_CalcNode* rt  =
  ((_CalcNode*)((BaseRef*)variablePtrs.lData)[theRoot->in_object]),
             * tip =
  ((_CalcNode*)((BaseRef*)variablePtrs.lData)[theRoot->nodes.data[0]->in_object]);

    _Parameter reslt, tmp;
    // sum over one branch in the direction from the root to the single leaf
    reslt = dsf->Translate2Frequencies ((*dsf)(index,0), rt->theProbs, true);
    tmp   = dsf->Translate2Frequencies ((*dsf)(index,1), tip->theProbs, true);
    // now perform one loop
    _Matrix* theTP = tip->GetCompExp();
    _Parameter* fastIdx = theTP->theData, *nodeProbs;
    // 4 cases are possible:
    // 1-1
    // 1-many
    // many-1
    // many-many
    if ((reslt>=0.0)&&(tmp>=0.0)) // 1-1
    {
        tmp = fastIdx[(long)reslt*cBase+(long)tmp];
        reslt = theProbs[(long)reslt]*tmp;
    }
    else
        if ((reslt>=0.0)&&(tmp<0.0)) // 1-many
        {
            tmp = 0;
            fastIdx += (long)reslt*cBase;
            nodeProbs = tip->GetProbs();
            for (long i=0; i<cBase; i++, fastIdx++,nodeProbs++)
            {
                tmp+=*fastIdx* *nodeProbs;
            }
            reslt = theProbs[(long)reslt]*tmp;
        }
        else
            if ((reslt<0.0)&&(tmp>=0.0)) // many to 1
            {
                fastIdx += (long)tmp;
                tmp = 0;
                nodeProbs  = rt->GetProbs();
                long i;
                for (i=0; i<cBase; i++, fastIdx+=cBase,nodeProbs++)
                {
                    *nodeProbs *=*fastIdx;
                }

                nodeProbs-=cBase;
                reslt = 0;
                fastIdx = GetProbs();
                for (i=0; i<cBase; i++, nodeProbs++, fastIdx++)
                {
                    reslt+=*fastIdx * *nodeProbs;
                }
            }
            else // many to many
            {
                nodeProbs  = rt->GetProbs();
                long i,j;
                for (i=0; i<cBase; i++, nodeProbs++)
                {
                    tmp = 0;
                    for (j=0; j<cBase; j++,fastIdx++)
                        tmp += *fastIdx* tip->GetProbs(j);
                    *nodeProbs*=tmp;
                }

                reslt = 0;
                fastIdx = GetProbs();
                nodeProbs-=cBase;
                for (i=0; i<cBase; i++, nodeProbs++, fastIdx++)
                {
                    reslt+=*fastIdx * *nodeProbs;
                }
            }
    if ((reslt<0.0)||(reslt==.0)) reslt = ALMOST_ZERO;
    return reslt;*/
}

//______________________________________________________________________________
// set leaf value and re-prune w/o reexping
// for subsequent entries into a datafilter
_Parameter _TheTree::ReleafTreeCharDegenerate(_DataSetFilter *dsf, long index) {

  _CalcNode *
  rt = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[theRoot->in_object]),
 *tip = ((_CalcNode *)(
     (BaseRef *)variablePtrs.lData)[theRoot->nodes.data[0]->in_object]);

  _Parameter reslt = 0.;

  // sum over one branch in the direction from the root to the single leaf
  char *thisState = dsf->GetColumn(index);

  long state1 = dsf->LookupConversion(thisState[dsf->theNodeMap.lData[0]],
                                      rt->theProbs),
       state2 = dsf->LookupConversion(thisState[dsf->theNodeMap.lData[1]],
                                      tip->theProbs);

  // now perform one loop
  _Parameter *fastIdx = tip->GetCompExp()->theData;
  //*nodeProbs;
  // 4 cases are possible:
  // 1-1
  // 1-many
  // many-1
  // many-many
  if (state1 >= 0 && state2 >= 0) { // 1-1
    reslt = theProbs[state1] * fastIdx[state1 * cBase + state2];
  } else if (state1 >= 0) { // 1-many
    _Parameter tmp = 0.;
    fastIdx += state1 * cBase;

    for (long i = 0; i < cBase; i++) {
      tmp += fastIdx[i] * tip->theProbs[i];
    }

    reslt = theProbs[state1] * tmp;
  } else if (state2 >= 0) { // many to 1
    fastIdx += state2;

    for (long i = 0; i < cBase; i++, fastIdx += cBase) {
      reslt += rt->theProbs[i] * *fastIdx * theProbs[i];
    }

  } else // many to many
    for (long i = 0; i < cBase; i++) {
      _Parameter tmp = 0.0;
      for (long j = 0; j < cBase; j++, fastIdx++) {
        tmp += *fastIdx * tip->theProbs[j];
      }

      reslt += tmp * rt->theProbs[i] * theProbs[i];
    }

  return reslt <= 0.0 ? ALMOST_ZERO : reslt;

  /*_CalcNode* rt  =
  ((_CalcNode*)((BaseRef*)variablePtrs.lData)[theRoot->in_object]),
             * tip =
  ((_CalcNode*)((BaseRef*)variablePtrs.lData)[theRoot->nodes.data[0]->in_object]);

    _Parameter reslt, tmp;
    // sum over one branch in the direction from the root to the single leaf

    char  *thisState = dsf->GetColumn(index);

    reslt = dsf->LookupConversion (thisState[dsf->theNodeMap.lData[0]],
  rt->theProbs);
    tmp   = dsf->LookupConversion (thisState[dsf->theNodeMap.lData[1]],
  tip->theProbs);
    //reslt = dsf->Translate2Frequencies ((*dsf)(index,0), rt->theProbs,
  true);
    //tmp   = dsf->Translate2Frequencies ((*dsf)(index,1), tip->theProbs,
  true);

    // now perform one loop
    _Matrix   * theTP = tip->compExp;
    _Parameter* fastIdx = theTP->theData,
              * nodeProbs;
    // 4 cases are possible:
    // 1-1
    // 1-many
    // many-1
    // many-many
    if ((reslt>=0.0)&&(tmp>=0.0)) // 1-1
    {
        tmp = fastIdx[(long)reslt*cBase+(long)tmp];
        reslt = theProbs[(long)reslt]*tmp;
    }
    else
        if ((reslt>=0.0)&&(tmp<0.0)) // 1-many
        {
            tmp = 0;
            fastIdx += (long)reslt*cBase;
            nodeProbs = tip->theProbs;
            for (long i=0; i<cBase; i++, fastIdx++,nodeProbs++)
            {
                tmp+=*fastIdx* *nodeProbs;
            }
            reslt = theProbs[(long)reslt]*tmp;
        }
        else
            if ((reslt<0.0)&&(tmp>=0.0)) // many to 1
            {
                fastIdx += (long)tmp;
                tmp = 0;
                nodeProbs  = rt->theProbs;
                long i;
                for (i=0; i<cBase; i++, fastIdx+=cBase,nodeProbs++)
                {
                    *nodeProbs *=*fastIdx;
                }

                nodeProbs-=cBase;
                reslt = 0;
                fastIdx = theProbs;
                for (i=0; i<cBase; i++, nodeProbs++, fastIdx++)
                {
                    reslt+=*fastIdx * *nodeProbs;
                }
            }
            else // many to many
            {
                nodeProbs  = rt->theProbs;
                long i,j;
                for (i=0; i<cBase; i++, nodeProbs++)
                {
                    tmp = 0;
                    for (j=0; j<cBase; j++,fastIdx++)
                        tmp += *fastIdx* tip->theProbs[j];
                    *nodeProbs*=tmp;
                }

                reslt = 0;
                fastIdx = theProbs;
                nodeProbs-=cBase;
                for (i=0; i<cBase; i++, nodeProbs++, fastIdx++)
                {
                    reslt+=*fastIdx * *nodeProbs;
                }
            }

    if ((reslt<0.0)||(reslt==.0)) reslt = ALMOST_ZERO;
    return reslt;*/
}

//______________________________________________________________________________
// set leaf value and re-prune w/o reexping
// for subsequent entries into a datafilter
_Parameter _TheTree::ReleafTree(_DataSetFilter *dsf, long index, long lastIndex,
                                long startLeaf, long endLeaf) {

  long nodeCount;
  _CalcNode *travNode, *theChildNode;
  _Parameter *fastIndex, *theProbbs, *stopper;
  node<long> *nodeChild;

  if (dsf->GetUnitLength() == 3) {
    char *c1, *c2, *c3, *o1, *o2, *o3;

    c1 = dsf->GetColumn(3 * index);
    c2 = dsf->GetColumn(3 * index + 1);
    c3 = dsf->GetColumn(3 * index + 2);

    if (lastIndex >= 0) {
      o1 = dsf->GetColumn(3 * lastIndex);
      o2 = dsf->GetColumn(3 * lastIndex + 1);
      o3 = dsf->GetColumn(3 * lastIndex + 2);
    }

    long ccount = dsf->conversionCache.lData[0], ccount2 = ccount * ccount,
         *ccodes = dsf->conversionCache.lData + 1,
         *tcodes = dsf->conversionCache.lData + 89;

    for (nodeCount = startLeaf; nodeCount <= endLeaf; nodeCount++) {
      long nMap = dsf->theNodeMap.lData[nodeCount];
      travNode = (_CalcNode *)(((BaseRef *)flatCLeaves.lData)[nodeCount]);
      char A = c1[nMap], B = c2[nMap], C = c3[nMap];

      if ((lastIndex < 0) || (A != o1[nMap]) || (B != o2[nMap]) ||
          (C != o3[nMap])) {
        A = ccodes[A - 40];
        B = ccodes[B - 40];
        C = ccodes[C - 40];

        if ((A == -1) || (B == -1) || (C == -1)) {
          travNode->lastState = dsf->Translate2Frequencies(
              (*dsf)(index, nodeCount), travNode->theProbs, true);
        } else {
          travNode->lastState = tcodes[A * ccount2 + B * ccount + C];
          if (travNode->lastState < 0) {
            travNode->lastState = dsf->Translate2Frequencies(
                (*dsf)(index, nodeCount), travNode->theProbs, true);
          } else {
            for (nMap = 0; nMap < travNode->lastState; nMap++) {
              travNode->theProbs[nMap] = 0.0;
            }
            travNode->theProbs[nMap++] = 1.0;
            for (; nMap < cBase; nMap++) {
              travNode->theProbs[nMap] = 0.0;
            }
          }
        }
        theChildNode = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[
            ((node<long> *)(flatLeaves.lData[nodeCount]))->parent->in_object]);
        theChildNode->cBase = -1;
      }
    }
  } else {
    if (dsf->GetUnitLength() == 2) {
      char *c1, *c2, *o1, *o2;

      c1 = dsf->GetColumn(2 * index);
      c2 = dsf->GetColumn(2 * index + 1);

      if (lastIndex >= 0) {
        o1 = dsf->GetColumn(2 * lastIndex);
        o2 = dsf->GetColumn(2 * lastIndex + 1);
      }

      long ccount = dsf->conversionCache.lData[0],
           *ccodes = dsf->conversionCache.lData + 1,
           *tcodes = dsf->conversionCache.lData + 89;

      for (nodeCount = startLeaf; nodeCount <= endLeaf; nodeCount++) {
        long nMap = dsf->theNodeMap.lData[nodeCount];
        travNode = (_CalcNode *)(((BaseRef *)flatCLeaves.lData)[nodeCount]);
        char A = c1[nMap], B = c2[nMap];

        if ((lastIndex < 0) || (A != o1[nMap]) || (B != o2[nMap])) {
          A = ccodes[A - 40];
          B = ccodes[B - 40];

          if ((A == -1) || (B == -1)) {
            travNode->lastState = dsf->Translate2Frequencies(
                (*dsf)(index, nodeCount), travNode->theProbs, true);
          } else {
            travNode->lastState = tcodes[A * ccount + B];
            if (travNode->lastState < 0) {
              travNode->lastState = dsf->Translate2Frequencies(
                  (*dsf)(index, nodeCount), travNode->theProbs, true);
            } else {
              for (nMap = 0; nMap < travNode->lastState; nMap++) {
                travNode->theProbs[nMap] = 0.0;
              }
              travNode->theProbs[nMap++] = 1;
              for (; nMap < cBase; nMap++) {
                travNode->theProbs[nMap] = 0.0;
              }
            }
          }
          theChildNode = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[(
              (node<long> *)(flatLeaves.lData[nodeCount]))->parent->in_object]);
          theChildNode->cBase = -1;
        }
      }
    } else {
      for (nodeCount = startLeaf; nodeCount <= endLeaf; nodeCount++) {
        travNode = (_CalcNode *)(((BaseRef *)flatCLeaves.lData)[nodeCount]);
        if ((lastIndex < 0) ||
            (!dsf->CompareTwoSites(lastIndex, index, nodeCount))) {
          travNode->lastState = dsf->Translate2Frequencies(
              (*dsf)(index, nodeCount), travNode->theProbs, true);
          theChildNode = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[(
              (node<long> *)(flatLeaves.lData[nodeCount]))->parent->in_object]);
          theChildNode->cBase = -1;
        }
      }
    }
  }

  for (nodeCount = leftiNodes.lData[startLeaf]; nodeCount < flatTree.lLength;
       nodeCount++) {
    theChildNode = (_CalcNode *)(((BaseRef *)flatTree.lData)[nodeCount]);
    // we can immediately place the changes in the parent node
    if (theChildNode->cBase < 0) {
      nodeChild = (node<long> *)flatNodes.lData[nodeCount];
      for (long i = 0; i < cBase; i++) {
        theChildNode->theProbs[i] = LIKELIHOOD_SCALER;
      }

      for (long k = 0; k < nodeChild->nodes.length; k++) {
        travNode = ((_CalcNode *)((
            BaseRef *)variablePtrs.lData)[nodeChild->nodes.data[k]->in_object]);
        fastIndex = travNode->compExp->theData;
        theProbbs = travNode->theProbs;
        long nZ = travNode->lastState;

        if (nZ >= 0) {
          _Parameter a = theProbbs[nZ];
          theProbbs = theChildNode->theProbs;
          fastIndex += nZ;
          stopper = theProbbs + cBase - cBase % 4;
          for (; theProbbs != stopper;) {
            *(theProbbs++) *= a * *fastIndex;
            fastIndex += cBase;
            *(theProbbs++) *= a * *fastIndex;
            fastIndex += cBase;
            *(theProbbs++) *= a * *fastIndex;
            fastIndex += cBase;
            *(theProbbs++) *= a * *fastIndex;
            fastIndex += cBase;
          }
          switch (cBase % 4) {
          case 1: {
            *theProbbs *= a * *fastIndex;
            break;
          }
          case 2: {
            *(theProbbs++) *= a * *fastIndex;
            fastIndex += cBase;
            *theProbbs *= a * *fastIndex;
            break;
          }
          case 3: {
            *(theProbbs++) *= a * *fastIndex;
            fastIndex += cBase;
            *(theProbbs++) *= a * *fastIndex;
            fastIndex += cBase;
            *theProbbs *= a * *fastIndex;
            break;
          }
          }
        } else {
          for (long i = 0; i < cBase; i++) {
            _Parameter tmp = 0.0;
            stopper = theProbbs + cBase - cBase % 4;
            // loop unrolled to depth 4
            for (; theProbbs != stopper;) {
              tmp += *(theProbbs++) * *(fastIndex++);
              tmp += *(theProbbs++) * *(fastIndex++);
              tmp += *(theProbbs++) * *(fastIndex++);
              tmp += *(theProbbs++) * *(fastIndex++);
            }
            switch (cBase % 4) {
            case 1: {
              tmp += *(theProbbs) * *(fastIndex++);
              break;
            }
            case 2: {
              tmp += *(theProbbs++) * *(fastIndex++);
              tmp += *(theProbbs) * *(fastIndex++);
              break;
            }
            case 3: {
              tmp += *(theProbbs++) * *(fastIndex++);
              tmp += *(theProbbs++) * *(fastIndex++);
              tmp += *(theProbbs) * *(fastIndex++);
              break;
            }
            }
            theChildNode->theProbs[i] *= tmp;
            theProbbs = travNode->theProbs;
          }
        }
      }
      theChildNode->cBase = cBase;
      nodeChild = nodeChild->parent;
      if (nodeChild) {
        ((_CalcNode *)((BaseRef *)variablePtrs.lData)[nodeChild->in_object])
            ->cBase = -1;
      }
    }
  }
  _Parameter result = 0;

  for (long i = 0; i < cBase; i++) {
    result += theProbs[i] * theChildNode->theProbs[i];
  }

  if (result <= 0.0) {
    /*_Matrix stashZeroValues (cBase, flatTree.lLength, false, true);

        _List   labels;

        for (nodeCount = 0; nodeCount<flatTree.lLength; nodeCount++)
        {
            theChildNode =
    (_CalcNode*)(((BaseRef*)flatTree.lData)[nodeCount]);
            labels <<  theChildNode->GetName();
            for (long ii=0; ii<cBase; ii++)
                stashZeroValues.Store (ii,nodeCount,
    theChildNode->theProbs[ii]);
        }

        _String cwt ("Report for site ");
        cwt = cwt & index;
        _HYChartWindow *chart = new _HYChartWindow(cwt, labels,
    stashZeroValues);
        chart->BringToFront();*/

    return ALMOST_ZERO;
  }

  return result;
}

//______________________________________________________________________________
// set leaf value and re-prune w/o reexping
// for subsequent entries into a datafilter
_Parameter _TheTree::ReleafTreeCache(_DataSetFilter *dsf, long index,
                                     long lastIndex, long startLeaf,
                                     long endLeaf, long position) {

  long nodeCount;
  _CalcNode *travNode, *theChildNode;
  _Parameter *fastIndex, *theProbbs, *stopper;
  node<long> *nodeChild;

  if (dsf->GetUnitLength() == 3) {
    char *c1, *c2, *c3, *o1, *o2, *o3;

    c1 = dsf->GetColumn(3 * index);
    c2 = dsf->GetColumn(3 * index + 1);
    c3 = dsf->GetColumn(3 * index + 2);

    if (lastIndex >= 0) {
      o1 = dsf->GetColumn(3 * lastIndex);
      o2 = dsf->GetColumn(3 * lastIndex + 1);
      o3 = dsf->GetColumn(3 * lastIndex + 2);
    }

    long ccount = dsf->conversionCache.lData[0], ccount2 = ccount * ccount,
         *ccodes = dsf->conversionCache.lData + 1,
         *tcodes = dsf->conversionCache.lData + 89;

    for (nodeCount = startLeaf; nodeCount <= endLeaf; nodeCount++) {
      long nMap = dsf->theNodeMap.lData[nodeCount];
      travNode = (_CalcNode *)(((BaseRef *)flatCLeaves.lData)[nodeCount]);
      char A = c1[nMap], B = c2[nMap], C = c3[nMap];

      if ((lastIndex < 0) || (A != o1[nMap]) || (B != o2[nMap]) ||
          (C != o3[nMap])) {
        A = ccodes[A - 40];
        B = ccodes[B - 40];
        C = ccodes[C - 40];

        if ((A == -1) || (B == -1) || (C == -1)) {
          travNode->lastState = dsf->Translate2Frequencies(
              (*dsf)(index, nodeCount), travNode->theProbs, true);
        } else {
          travNode->lastState = tcodes[A * ccount2 + B * ccount + C];
          if (travNode->lastState < 0) {
            travNode->lastState = dsf->Translate2Frequencies(
                (*dsf)(index, nodeCount), travNode->theProbs, true);
          } else {
            for (nMap = 0; nMap < travNode->lastState; nMap++) {
              travNode->theProbs[nMap] = 0.0;
            }
            travNode->theProbs[nMap++] = 1.0;
            for (; nMap < cBase; nMap++) {
              travNode->theProbs[nMap] = 0.0;
            }
          }
        }
        theChildNode = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[
            ((node<long> *)(flatLeaves.lData[nodeCount]))->parent->in_object]);
        theChildNode->cBase = -1;
      }
    }
  } else {
    if (dsf->GetUnitLength() == 2) {
      char *c1, *c2, *o1, *o2;

      c1 = dsf->GetColumn(2 * index);
      c2 = dsf->GetColumn(2 * index + 1);

      if (lastIndex >= 0) {
        o1 = dsf->GetColumn(2 * lastIndex);
        o2 = dsf->GetColumn(2 * lastIndex + 1);
      }

      long ccount = dsf->conversionCache.lData[0],
           *ccodes = dsf->conversionCache.lData + 1,
           *tcodes = dsf->conversionCache.lData + 89;

      for (nodeCount = startLeaf; nodeCount <= endLeaf; nodeCount++) {
        long nMap = dsf->theNodeMap.lData[nodeCount];
        travNode = (_CalcNode *)(((BaseRef *)flatCLeaves.lData)[nodeCount]);
        char A = c1[nMap], B = c2[nMap];

        if ((lastIndex < 0) || (A != o1[nMap]) || (B != o2[nMap])) {
          A = ccodes[A - 40];
          B = ccodes[B - 40];

          if ((A == -1) || (B == -1)) {
            travNode->lastState = dsf->Translate2Frequencies(
                (*dsf)(index, nodeCount), travNode->theProbs, true);
          } else {
            travNode->lastState = tcodes[A * ccount + B];
            if (travNode->lastState < 0) {
              travNode->lastState = dsf->Translate2Frequencies(
                  (*dsf)(index, nodeCount), travNode->theProbs, true);
            } else {
              for (nMap = 0; nMap < travNode->lastState; nMap++) {
                travNode->theProbs[nMap] = 0.0;
              }
              travNode->theProbs[nMap++] = 1;
              for (; nMap < cBase; nMap++) {
                travNode->theProbs[nMap] = 0.0;
              }
            }
          }
          theChildNode = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[(
              (node<long> *)(flatLeaves.lData[nodeCount]))->parent->in_object]);
          theChildNode->cBase = -1;
        }
      }
    } else {
      for (nodeCount = startLeaf; nodeCount <= endLeaf; nodeCount++) {
        travNode = (_CalcNode *)(((BaseRef *)flatCLeaves.lData)[nodeCount]);
        if ((lastIndex < 0) ||
            (!dsf->CompareTwoSites(lastIndex, index, nodeCount))) {
          travNode->lastState = dsf->Translate2Frequencies(
              (*dsf)(index, nodeCount), travNode->theProbs, true);
          theChildNode = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[(
              (node<long> *)(flatLeaves.lData[nodeCount]))->parent->in_object]);
          theChildNode->cBase = -1;
        }
      }
    }
  }

  long f = topLevelNodes.lLength - 1, mm = topLevelNodes.lData[f],
       startINode = leftiNodes.lData[startLeaf];

  for (long nI = 0; nI < f; nI++, mm >>= 1) {
    //if (1) {
    if (mm & 1) {
      // marked for recalculation
      long startAt = 0;

      if (nI) {
        startAt = topLevelNodes.lData[nI - 1] + 1;
      }

      if (startAt < startINode) {
        startAt = startINode;
      }

      for (nodeCount = startAt; nodeCount <= topLevelNodes.lData[nI];
           nodeCount++) {

        theChildNode = (_CalcNode *)(((BaseRef *)flatTree.lData)[nodeCount]);
        // we can immediately place the changes in the parent node
        if (theChildNode->cBase < 0) {
          nodeChild = (node<long> *)flatNodes.lData[nodeCount];
          for (long i = 0; i < cBase; i++) {
            theChildNode->theProbs[i] = LIKELIHOOD_SCALER;
          }

          for (long k = 0; k < nodeChild->nodes.length; k++) {
            travNode = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[
                nodeChild->nodes.data[k]->in_object]);
            fastIndex = travNode->compExp->theData;
            theProbbs = travNode->theProbs;
            long nZ = travNode->lastState;
            if (nZ >= 0) {
              _Parameter a = theProbbs[nZ];
              theProbbs = theChildNode->theProbs;
              fastIndex += nZ;
              stopper = theProbbs + cBase;
              for (; theProbbs != stopper; fastIndex += cBase, theProbbs++) {
                *theProbbs *= a * *fastIndex;
              }

            } else {
              /*for (long i=0; i<cBase; i++)
                            {
                                _Parameter tmp = 0.0;
                                for (long j=0; j< cBase; j++)
                                    tmp += theProbbs[j] *
              fastIndex [j];

                                theChildNode->theProbs[i]*=tmp;
                                fastIndex += cBase;
                            }*/

              long rem = cBase % 4;
              stopper = theProbbs + (cBase - rem);
              if (rem == 1) {
                for (long i = 0; i < cBase; i++) {
                  _Parameter tmp = 0.0;
                  for (; theProbbs != stopper;) {
                    tmp += *theProbbs * *fastIndex;
                    theProbbs++;
                    fastIndex++;
                    tmp += *theProbbs * *fastIndex;
                    theProbbs++;
                    fastIndex++;
                    tmp += *theProbbs * *fastIndex;
                    theProbbs++;
                    fastIndex++;
                    tmp += *theProbbs * *fastIndex;
                    theProbbs++;
                    fastIndex++;
                  }
                  tmp += *(theProbbs) * *(fastIndex++);
                  theChildNode->theProbs[i] *= tmp;
                  theProbbs = travNode->theProbs;
                }
              } else if (rem == 2) {
                for (long i = 0; i < cBase; i++) {
                  _Parameter tmp = 0.0;
                  for (; theProbbs != stopper;) {
                    tmp += *theProbbs * *fastIndex;
                    theProbbs++;
                    fastIndex++;
                    tmp += *theProbbs * *fastIndex;
                    theProbbs++;
                    fastIndex++;
                    tmp += *theProbbs * *fastIndex;
                    theProbbs++;
                    fastIndex++;
                    tmp += *theProbbs * *fastIndex;
                    theProbbs++;
                    fastIndex++;
                  }
                  tmp += *(theProbbs++) * *(fastIndex++);
                  tmp += *(theProbbs) * *(fastIndex++);
                  theChildNode->theProbs[i] *= tmp;
                  theProbbs = travNode->theProbs;
                }

              } else if (rem == 3) {
                for (long i = 0; i < cBase; i++) {
                  _Parameter tmp = 0.0;
                  for (; theProbbs != stopper;) {
                    tmp += *theProbbs * *fastIndex;
                    theProbbs++;
                    fastIndex++;
                    tmp += *theProbbs * *fastIndex;
                    theProbbs++;
                    fastIndex++;
                    tmp += *theProbbs * *fastIndex;
                    theProbbs++;
                    fastIndex++;
                    tmp += *theProbbs * *fastIndex;
                    theProbbs++;
                    fastIndex++;
                  }
                  tmp += *(theProbbs++) * *(fastIndex++);
                  tmp += *(theProbbs++) * *(fastIndex++);
                  tmp += *(theProbbs) * *(fastIndex++);
                  theChildNode->theProbs[i] *= tmp;
                  theProbbs = travNode->theProbs;
                }

              } else {
                for (long i = 0; i < cBase; i++) {
                  _Parameter tmp = 0.0;
                  for (; theProbbs != stopper;) {
                    tmp += *theProbbs * *fastIndex;
                    theProbbs++;
                    fastIndex++;
                    tmp += *theProbbs * *fastIndex;
                    theProbbs++;
                    fastIndex++;
                    tmp += *theProbbs * *fastIndex;
                    theProbbs++;
                    fastIndex++;
                    tmp += *theProbbs * *fastIndex;
                    theProbbs++;
                    fastIndex++;
                  }
                  theChildNode->theProbs[i] *= tmp;
                  theProbbs = travNode->theProbs;
                }
              }
            }
          }
          theChildNode->cBase = cBase;
          nodeChild = nodeChild->parent;
          if (nodeChild) {
            ((_CalcNode *)((BaseRef *)variablePtrs.lData)[nodeChild->in_object])
                ->cBase = -1;
          }
        }
      }

      theChildNode =
          (_CalcNode *)(((BaseRef *)flatTree.lData)[topLevelNodes.lData[nI]]);
      _Parameter *nodeProbs = theChildNode->theProbs;
      fastIndex = rootIChildrenCache +
                  position * cBase * (topLevelNodes.lLength - 1) + cBase * nI;
      for (long i = 0; i < cBase; i++) {
        *(fastIndex++) = *(nodeProbs++);
      }
    } else {
      _Parameter *nodeProbs =
          ((_CalcNode *)(((BaseRef *)flatTree.lData)[topLevelNodes.lData[nI]]))
              ->theProbs;
      fastIndex = rootIChildrenCache +
                  position * cBase * (topLevelNodes.lLength - 1) + cBase * nI;
      for (long i = 0; i < cBase; i++) {
        *(nodeProbs++) = *(fastIndex++);
      }

      nodeChild =
          ((node<long> *)(flatNodes.lData[topLevelNodes.lData[nI]]))->parent;
      ((_CalcNode *)((BaseRef *)variablePtrs.lData)[nodeChild->in_object])
          ->cBase = -1;
    }
  }

  for (nodeCount = topLevelNodes.lData[f - 1] + 1; nodeCount < flatTree.lLength;
       nodeCount++) {
    theChildNode = (_CalcNode *)(((BaseRef *)flatTree.lData)[nodeCount]);
    if (theChildNode->cBase < 0) {
      nodeChild = (node<long> *)flatNodes.lData[nodeCount];
      for (long i = 0; i < cBase; i++) {
        theChildNode->theProbs[i] = LIKELIHOOD_SCALER;
      }

      for (long k = 0; k < nodeChild->nodes.length; k++) {
        travNode = ((_CalcNode *)((
            BaseRef *)variablePtrs.lData)[nodeChild->nodes.data[k]->in_object]);
        fastIndex = travNode->compExp->theData;
        theProbbs = travNode->theProbs;
        long nZ = travNode->lastState;

        if (nZ >= 0) {
          _Parameter a = theProbbs[nZ];
          theProbbs = theChildNode->theProbs;
          fastIndex += nZ;
          stopper = theProbbs + cBase;
          for (; theProbbs != stopper; theProbbs++, fastIndex += cBase) {
            *theProbbs *= a * *fastIndex;
          }
        } else {
          /*for (long i=0; i<cBase; i++)
                    {
                        _Parameter tmp = 0.0;
                        for (long j=0;j<cBase;j++)
                            tmp += theProbbs[j] *fastIndex[j];

                        theChildNode->theProbs[i]*=tmp;
                        fastIndex += cBase;
                    }*/

          long rem = cBase % 4;
          stopper = theProbbs + (cBase - rem);
          if (rem == 1) {
            for (long i = 0; i < cBase; i++) {
              _Parameter tmp = 0.0;
              for (; theProbbs != stopper;) {
                tmp += *theProbbs * *fastIndex;
                theProbbs++;
                fastIndex++;
                tmp += *theProbbs * *fastIndex;
                theProbbs++;
                fastIndex++;
                tmp += *theProbbs * *fastIndex;
                theProbbs++;
                fastIndex++;
                tmp += *theProbbs * *fastIndex;
                theProbbs++;
                fastIndex++;
              }
              tmp += *(theProbbs) * *(fastIndex++);
              theChildNode->theProbs[i] *= tmp;
              theProbbs = travNode->theProbs;
            }
          } else if (rem == 2) {
            for (long i = 0; i < cBase; i++) {
              _Parameter tmp = 0.0;
              for (; theProbbs != stopper;) {
                tmp += *theProbbs * *fastIndex;
                theProbbs++;
                fastIndex++;
                tmp += *theProbbs * *fastIndex;
                theProbbs++;
                fastIndex++;
                tmp += *theProbbs * *fastIndex;
                theProbbs++;
                fastIndex++;
                tmp += *theProbbs * *fastIndex;
                theProbbs++;
                fastIndex++;
              }
              tmp += *(theProbbs++) * *(fastIndex++);
              tmp += *(theProbbs) * *(fastIndex++);
              theChildNode->theProbs[i] *= tmp;
              theProbbs = travNode->theProbs;
            }
          } else if (rem == 3) {
            for (long i = 0; i < cBase; i++) {
              _Parameter tmp = 0.0;
              for (; theProbbs != stopper;) {
                tmp += *theProbbs * *fastIndex;
                theProbbs++;
                fastIndex++;
                tmp += *theProbbs * *fastIndex;
                theProbbs++;
                fastIndex++;
                tmp += *theProbbs * *fastIndex;
                theProbbs++;
                fastIndex++;
                tmp += *theProbbs * *fastIndex;
                theProbbs++;
                fastIndex++;
              }
              tmp += *(theProbbs++) * *(fastIndex++);
              tmp += *(theProbbs++) * *(fastIndex++);
              tmp += *(theProbbs) * *(fastIndex++);
              theChildNode->theProbs[i] *= tmp;
              theProbbs = travNode->theProbs;
            }
          } else {
            for (long i = 0; i < cBase; i++) {
              _Parameter tmp = 0.0;
              for (; theProbbs != stopper;) {
                tmp += *theProbbs * *fastIndex;
                theProbbs++;
                fastIndex++;
                tmp += *theProbbs * *fastIndex;
                theProbbs++;
                fastIndex++;
                tmp += *theProbbs * *fastIndex;
                theProbbs++;
                fastIndex++;
                tmp += *theProbbs * *fastIndex;
                theProbbs++;
                fastIndex++;
              }
              theChildNode->theProbs[i] *= tmp;
              theProbbs = travNode->theProbs;
            }
          }
        }
      }
      theChildNode->cBase = cBase;
      nodeChild = nodeChild->parent;
      if (nodeChild) {
        ((_CalcNode *)((BaseRef *)variablePtrs.lData)[nodeChild->in_object])
            ->cBase = -1;
      }
    }
  }

  _Parameter result = 0.;

  for (long i = 0; i < cBase; i++) {
    result += theProbs[i] * theChildNode->theProbs[i];
  }

  if (result <= 0.0) {
    return ALMOST_ZERO;
  }

  return result;
}

//______________________________________________________________________________
#if USE_SCALING_TO_FIX_UNDERFLOW
_Parameter _TheTree::ThreadReleafTreeCache(_DataSetFilter *dsf, long index,
                                           long lastIndex, long startLeaf,
                                           long endLeaf, long position,
                                           long offset, long fixAttempt,
                                           _Parameter localScalingFactor)
#else
    _Parameter _TheTree::ThreadReleafTreeCache(_DataSetFilter *dsf, long index,
                                               long lastIndex, long startLeaf,
                                               long endLeaf, long position,
                                               long offset)

#endif
    // set leaf value and re-prune w/o reexping
    // for subsequent entries into a datafilter
    {

  long nodeCount,
      *ns = nodeStates + (flatLeaves.lLength + flatNodes.lLength) * offset;

  _CalcNode *travNode, *theChildNode;

  _Parameter *fastIndex, *theProbbs, *stopper,
      *mlc = marginalLikelihoodCache +
             cBase * (flatLeaves.lLength + flatNodes.lLength) * offset;

  node<long> *nodeChild;

  char *nm = nodeMarkers + (flatNodes.lLength) * offset,
       unitSize = dsf->GetUnitLength();

#if USE_SCALING_TO_FIX_UNDERFLOW
  _Parameter scalingFactor = scalingForUnderflow->theData[index];
  if (lastIndex >= 0 &&
      scalingFactor != scalingForUnderflow->theData[lastIndex]) {
    lastIndex = -1;
    startLeaf = 0;
    endLeaf = flatLeaves.lLength - 1;
  }
  bool overflowFlag = fixAttempt < 0;
  if (overflowFlag) {
    fixAttempt = -fixAttempt;
  }

#endif

  if (unitSize == 3) {

    char *c1, *c2, *c3, *o1, *o2, *o3;

    c1 = dsf->GetColumn(3 * index);
    c2 = dsf->GetColumn(3 * index + 1);
    c3 = dsf->GetColumn(3 * index + 2);

    if (lastIndex >= 0) {
      o1 = dsf->GetColumn(3 * lastIndex);
      o2 = dsf->GetColumn(3 * lastIndex + 1);
      o3 = dsf->GetColumn(3 * lastIndex + 2);
    }

    long ccount = dsf->conversionCache.lData[0], ccount2 = ccount * ccount,
         *ccodes = dsf->conversionCache.lData + 1,
         *tcodes = dsf->conversionCache.lData + 89;

    for (nodeCount = startLeaf; nodeCount <= endLeaf; nodeCount++) {
      long nMap = dsf->theNodeMap.lData[nodeCount];
      char A = c1[nMap], B = c2[nMap], C = c3[nMap];

      if (lastIndex < 0 || A != o1[nMap] || B != o2[nMap] || C != o3[nMap]) {
        A = ccodes[A - 40];
        B = ccodes[B - 40];
        C = ccodes[C - 40];

        fastIndex = mlc + cBase * nodeCount;

        if (A == -1 || B == -1 || C == -1) {
          _String codon(3, false);
          dsf->GrabSite(index, nodeCount, codon);
          ns[nodeCount] = dsf->Translate2Frequencies(codon, fastIndex, true);
#if USE_SCALING_TO_FIX_UNDERFLOW
          if (scalingFactor != 1.0)
            for (nMap = 0; nMap < cBase; nMap++) {
              fastIndex[nMap] *= scalingFactor;
            }
#endif
        } else {
          long nState = ns[nodeCount] = tcodes[A * ccount2 + B * ccount + C];
          if (nState < 0) {
            _String codon(3, false);
            dsf->GrabSite(index, nodeCount, codon);
            ns[nodeCount] = dsf->Translate2Frequencies(codon, fastIndex, true);
#if USE_SCALING_TO_FIX_UNDERFLOW
            if (scalingFactor != 1.0)
              for (nMap = 0; nMap < cBase; nMap++) {
                fastIndex[nMap] *= scalingFactor;
              }
#endif
          } else {
            for (long z = 0; z < cBase; z++) {
              fastIndex[z] = 0.0;
            }
#if USE_SCALING_TO_FIX_UNDERFLOW
            fastIndex[nState] = scalingFactor;
#else
            fastIndex[nState] = 1.0;
#endif
          }
        }
        theChildNode = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[
            ((node<long> *)(flatLeaves.lData[nodeCount]))->parent->in_object]);
        nm[theChildNode->nodeIndex - flatLeaves.lLength] = -1;
      }
    }
  } else {
    if (unitSize == 2) {
      char *c1, *c2, *o1, *o2;

      c1 = dsf->GetColumn(2 * index);
      c2 = dsf->GetColumn(2 * index + 1);

      if (lastIndex >= 0) {
        o1 = dsf->GetColumn(2 * lastIndex);
        o2 = dsf->GetColumn(2 * lastIndex + 1);
      }

      long ccount = dsf->conversionCache.lData[0],
           *ccodes = dsf->conversionCache.lData + 1,
           *tcodes = dsf->conversionCache.lData + 89;

      for (nodeCount = startLeaf; nodeCount <= endLeaf; nodeCount++) {
        long nMap = dsf->theNodeMap.lData[nodeCount];
        //travNode = (_CalcNode*)(((BaseRef*)flatCLeaves.lData)[nodeCount]);
        char A = c1[nMap], B = c2[nMap];

        if ((lastIndex < 0) || (A != o1[nMap]) || (B != o2[nMap])) {
          A = ccodes[A - 40];
          B = ccodes[B - 40];

          fastIndex = mlc + cBase * nodeCount;

          if ((A == -1) || (B == -1)) {
            _String di(2, false);
            dsf->GrabSite(index, nodeCount, di);
            ns[nodeCount] = dsf->Translate2Frequencies(di, fastIndex, true);
          } else {
            long nState = ns[nodeCount] = tcodes[A * ccount + B];

            if (nState < 0) {
              _String di(2, false);
              dsf->GrabSite(index, nodeCount, di);
              ns[nodeCount] = dsf->Translate2Frequencies(di, fastIndex, true);
            } else {
              for (nMap = 0; nMap < nState; nMap++) {
                fastIndex[nMap] = 0.0;
              }
              fastIndex[nMap++] = 1.;
              for (; nMap < cBase; nMap++) {
                fastIndex[nMap] = 0.0;
              }
            }
          }
          theChildNode = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[(
              (node<long> *)(flatLeaves.lData[nodeCount]))->parent->in_object]);
          nm[theChildNode->nodeIndex - flatLeaves.lLength] = -1;
        }
      }
    } else {
      for (nodeCount = startLeaf; nodeCount <= endLeaf; nodeCount++) {
        //travNode = (_CalcNode*)(((BaseRef*)flatCLeaves.lData)[nodeCount]);
        if ((lastIndex < 0) ||
            (!dsf->CompareTwoSites(lastIndex, index, nodeCount))) {
          _String dPoint(unitSize, false);
          dsf->GrabSite(index, nodeCount, dPoint);
          ns[nodeCount] =
              dsf->Translate2Frequencies(dPoint, mlc + cBase * nodeCount, true);
          theChildNode = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[(
              (node<long> *)(flatLeaves.lData[nodeCount]))->parent->in_object]);
          nm[theChildNode->nodeIndex - flatLeaves.lLength] = -1;
        }
      }
    }
  }

  long f = topLevelNodes.lLength - 1, mm = topLevelNodes.lData[f],
       startINode = leftiNodes.lData[startLeaf], rem = cBase % 4,
       upTo = cBase - rem;

  for (long nI = 0; nI < f; nI++, mm >>= 1) {
    //if (1) {
    if (mm & 1) {

      // marked for recalculation
      long startAt = 0;

      if (nI) {
        startAt = topLevelNodes.lData[nI - 1] + 1;
      }

      if (startAt < startINode) {
        startAt = startINode;
      }

      for (nodeCount = startAt; nodeCount <= topLevelNodes.lData[nI];
           nodeCount++) {
        if (nm[nodeCount] < 0) {
          nodeChild = (node<long> *)flatNodes.lData[nodeCount];
          _Parameter *theseProbs =
              mlc + cBase * (nodeCount + flatLeaves.lLength);
          for (long i = 0; i < cBase; i++) {
            theseProbs[i] = LIKELIHOOD_SCALER;
          }

          for (long k = 0; k < nodeChild->nodes.length; k++) {
            travNode = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[
                nodeChild->nodes.data[k]->in_object]);
            fastIndex = travNode->compExp->theData;
            _Parameter *cProbs = mlc + cBase * (travNode->nodeIndex);
            theProbbs = cProbs;
            long nZ = ns[travNode->nodeIndex];
            if (nZ >= 0) {
              _Parameter a = theProbbs[nZ];
              theProbbs = theseProbs;
              fastIndex += nZ;
              stopper = theProbbs + cBase;
              for (; theProbbs != stopper; fastIndex += cBase, theProbbs++)
#if USE_SCALING_TO_FIX_UNDERFLOW
                *theProbbs *= a * *fastIndex; // * scalingFactor;
#else
              *theProbbs *= a * *fastIndex;
#endif
            } else {
              stopper = theProbbs + (cBase - rem);
              if (rem == 1) {
                for (long i = 0; i < cBase; i++) {
                  _Parameter tmp1 = 0.0, tmp2 = 0.0, tmp3 = 0.0, tmp4 = 0.0;

                  for (long k = 0; k < upTo; k += 4) {
                    tmp1 += cProbs[k] * fastIndex[k];
                    tmp2 += cProbs[k + 1] * fastIndex[k + 1];
                    tmp3 += cProbs[k + 2] * fastIndex[k + 2];
                    tmp4 += cProbs[k + 3] * fastIndex[k + 3];
                  }

                  theseProbs[i] *= tmp1 + cProbs[upTo] * fastIndex[upTo] +
                                   tmp2 + tmp3 + tmp4;
                  fastIndex += cBase;
                }
              } else if (rem == 2) {
                for (long i = 0; i < cBase; i++) {
                  _Parameter tmp1 = 0.0, tmp2 = 0.0, tmp3 = 0.0, tmp4 = 0.0;

                  for (long k = 0; k < upTo; k += 4) {
                    tmp1 += cProbs[k] * fastIndex[k];
                    tmp2 += cProbs[k + 1] * fastIndex[k + 1];
                    tmp3 += cProbs[k + 2] * fastIndex[k + 2];
                    tmp4 += cProbs[k + 3] * fastIndex[k + 3];
                  }

                  theseProbs[i] *= tmp1 + tmp2 + tmp3 + tmp4 +
                                   +cProbs[upTo] * fastIndex[upTo] +
                                   +cProbs[upTo + 1] * fastIndex[upTo + 1];
                  fastIndex += cBase;
                }
              } else if (rem == 3)
                for (long i = 0; i < cBase; i++) {
                  _Parameter tmp1 = 0.0, tmp2 = 0.0, tmp3 = 0.0, tmp4 = 0.0;

                  for (long k = 0; k < upTo; k += 4) {
                    tmp1 += cProbs[k] * fastIndex[k];
                    tmp2 += cProbs[k + 1] * fastIndex[k + 1];
                    tmp3 += cProbs[k + 2] * fastIndex[k + 2];
                    tmp4 += cProbs[k + 3] * fastIndex[k + 3];
                  }

                  theseProbs[i] *= tmp1 + tmp2 + tmp3 + tmp4 +
                                   +cProbs[upTo] * fastIndex[upTo] +
                                   +cProbs[upTo + 1] * fastIndex[upTo + 1] +
                                   +cProbs[upTo + 2] * fastIndex[upTo + 2];

                  fastIndex += cBase;
                }
              else
                for (long i = 0; i < cBase; i++) {
                  _Parameter tmp1 = 0.0, tmp2 = 0.0, tmp3 = 0.0, tmp4 = 0.0;

                  for (long k = 0; k < upTo; k += 4) {
                    tmp1 += cProbs[k] * fastIndex[k];
                    tmp2 += cProbs[k + 1] * fastIndex[k + 1];
                    tmp3 += cProbs[k + 2] * fastIndex[k + 2];
                    tmp4 += cProbs[k + 3] * fastIndex[k + 3];
                  }

                  theseProbs[i] *= tmp1 + tmp2 + tmp3 + tmp4;
                  fastIndex += cBase;
                }
            }
#if USE_SCALING_TO_FIX_UNDERFLOW
            long nodeCheck = 0;
            for (nodeCheck = 0; nodeCheck < cBase; nodeCheck++)
              if (theseProbs[nodeCheck] > 0.0 ||
                  theseProbs[nodeCheck] == HUGE_VAL) {
                break;
              }

            if (nodeCheck == cBase || theseProbs[nodeCheck] == HUGE_VAL) {
              return doScaling(dsf, index, position, offset, fixAttempt,
                               localScalingFactor, overflowFlag,
                               nodeCheck < cBase);
            }
#endif
          }
          nm[nodeCount] = 0;
          nodeChild = nodeChild->parent;
          if (nodeChild) {
            nm[((_CalcNode *)((BaseRef *)variablePtrs.lData)[
                nodeChild->in_object])->nodeIndex - flatLeaves.lLength] = -1;
          }
        }
      }

      theChildNode =
          (_CalcNode *)(((BaseRef *)flatTree.lData)[topLevelNodes.lData[nI]]);
      _Parameter *nodeProbs = mlc + cBase * (theChildNode->nodeIndex);
      fastIndex = rootIChildrenCache +
                  position * cBase * (topLevelNodes.lLength - 1) + cBase * nI;
      for (long i = 0; i < cBase; i++) {
        *(fastIndex++) = *(nodeProbs++);
      }
    } else {
      _Parameter *nodeProbs =
          mlc +
          cBase * (((_CalcNode *)(
                      ((BaseRef *)flatTree.lData)[topLevelNodes.lData[nI]]))
                       ->nodeIndex);

      fastIndex = rootIChildrenCache +
                  position * cBase * (topLevelNodes.lLength - 1) + cBase * nI;

      for (long i = 0; i < cBase; i++) {
        *(nodeProbs++) = *(fastIndex++);
      }

      nodeChild =
          ((node<long> *)(flatNodes.lData[topLevelNodes.lData[nI]]))->parent;

      nm[((_CalcNode *)((BaseRef *)variablePtrs.lData)[nodeChild->in_object])
             ->nodeIndex - flatLeaves.lLength] = -1;
    }
  }

  for (nodeCount = topLevelNodes.lData[f - 1] + 1; nodeCount < flatTree.lLength;
       nodeCount++) {
    if (nm[nodeCount] < 0) {
      nodeChild = (node<long> *)flatNodes.lData[nodeCount];
      _Parameter *theseProbs = mlc + cBase * (nodeCount + flatLeaves.lLength);
      for (long i = 0; i < cBase; i++) {
        theseProbs[i] = LIKELIHOOD_SCALER;
      }

      for (long k = 0; k < nodeChild->nodes.length; k++) {
        travNode = ((_CalcNode *)((
            BaseRef *)variablePtrs.lData)[nodeChild->nodes.data[k]->in_object]);
        fastIndex = travNode->compExp->theData;
        _Parameter *cProbs = mlc + cBase * (travNode->nodeIndex);
        theProbbs = cProbs;
        long nZ = ns[travNode->nodeIndex];
        if (nZ >= 0) {
          _Parameter a = theProbbs[nZ];
          theProbbs = theseProbs;
          fastIndex += nZ;
          stopper = theProbbs + cBase;
          for (; theProbbs != stopper; fastIndex += cBase, theProbbs++)
#if USE_SCALING_TO_FIX_UNDERFLOW
            *theProbbs *= a * *fastIndex; // * scalingFactor;
#else
          *theProbbs *= a * *fastIndex;
#endif
        } else {
          stopper = theProbbs + (cBase - rem);
          if (rem == 1) {
            for (long i = 0; i < cBase; i++) {
              _Parameter tmp1 = 0.0, tmp2 = 0.0, tmp3 = 0.0, tmp4 = 0.0;

              for (long k = 0; k < upTo; k += 4) {
                tmp1 += cProbs[k] * fastIndex[k];
                tmp2 += cProbs[k + 1] * fastIndex[k + 1];
                tmp3 += cProbs[k + 2] * fastIndex[k + 2];
                tmp4 += cProbs[k + 3] * fastIndex[k + 3];
              }

              theseProbs[i] *=
                  tmp1 + cProbs[upTo] * fastIndex[upTo] + tmp2 + tmp3 + tmp4;
              fastIndex += cBase;
            }
          } else if (rem == 2) {
            for (long i = 0; i < cBase; i++) {
              _Parameter tmp1 = 0.0, tmp2 = 0.0, tmp3 = 0.0, tmp4 = 0.0;

              for (long k = 0; k < upTo; k += 4) {
                tmp1 += cProbs[k] * fastIndex[k];
                tmp2 += cProbs[k + 1] * fastIndex[k + 1];
                tmp3 += cProbs[k + 2] * fastIndex[k + 2];
                tmp4 += cProbs[k + 3] * fastIndex[k + 3];
              }

              theseProbs[i] *=
                  tmp1 + cProbs[upTo] * fastIndex[upTo] +
                  cProbs[upTo + 1] * fastIndex[upTo + 1] + tmp2 + tmp3 + tmp4;
              fastIndex += cBase;
            }

          } else if (rem == 3)
            for (long i = 0; i < cBase; i++) {
              _Parameter tmp1 = 0.0, tmp2 = 0.0, tmp3 = 0.0, tmp4 = 0.0;

              for (long k = 0; k < upTo; k += 4) {
                tmp1 += cProbs[k] * fastIndex[k];
                tmp2 += cProbs[k + 1] * fastIndex[k + 1];
                tmp3 += cProbs[k + 2] * fastIndex[k + 2];
                tmp4 += cProbs[k + 3] * fastIndex[k + 3];
              }

              theseProbs[i] *=
                  tmp1 + cProbs[upTo] * fastIndex[upTo] +
                  cProbs[upTo + 1] * fastIndex[upTo + 1] +
                  cProbs[upTo + 2] * fastIndex[upTo + 2] + tmp2 + tmp3 + tmp4;
              fastIndex += cBase;
            }
          else {
            for (long i = 0; i < cBase; i++) {
              _Parameter tmp1 = 0.0, tmp2 = 0.0, tmp3 = 0.0, tmp4 = 0.0;

              for (long k = 0; k < upTo; k += 4) {
                tmp1 += cProbs[k] * fastIndex[k];
                tmp2 += cProbs[k + 1] * fastIndex[k + 1];
                tmp3 += cProbs[k + 2] * fastIndex[k + 2];
                tmp4 += cProbs[k + 3] * fastIndex[k + 3];
              }

              theseProbs[i] *= tmp1 + tmp2 + tmp3 + tmp4;
              fastIndex += cBase;
            }
          }
        }
      }
      nm[nodeCount] = 0;
      nodeChild = nodeChild->parent;
      if (nodeChild) {
        nm[((_CalcNode *)((BaseRef *)variablePtrs.lData)[nodeChild->in_object])
               ->nodeIndex - flatLeaves.lLength] = -1;
      }

#if USE_SCALING_TO_FIX_UNDERFLOW
      long nodeCheck = 0;
      for (nodeCheck = 0; nodeCheck < cBase; nodeCheck++)
        if (theseProbs[nodeCheck] > 0.0 || theseProbs[nodeCheck] == HUGE_VAL) {
          break;
        }

      if (nodeCheck == cBase || theseProbs[nodeCheck] == HUGE_VAL) {
        return doScaling(dsf, index, position, offset, fixAttempt,
                         localScalingFactor, overflowFlag, nodeCheck < cBase);
      }
#endif

    }
  }

  _Parameter result = 0;
  fastIndex = mlc + cBase * (flatNodes.lLength + flatLeaves.lLength - 1);

  for (long i = 0; i < cBase; i++) {
    result += theProbs[i] * fastIndex[i];
  }

#if USE_SCALING_TO_FIX_UNDERFLOW
  if (result <= 0.0) { // underflow
    return doScaling(dsf, index, position, offset, fixAttempt,
                     localScalingFactor, overflowFlag, false);
  } else if (result == HUGE_VAL) {
    return doScaling(dsf, index, position, offset, fixAttempt,
                     localScalingFactor, overflowFlag, true);
  }
#else
  if (result <= 0.0) {
    return ALMOST_ZERO;
  }
#endif
  return result;
}

//______________________________________________________________________________
// set leaf value and re-prune w/o reexping
// for subsequent entries into a datafilter
_Parameter _TheTree::ReleafTreeChar4(_DataSetFilter *dsf, long index,
                                     long lastIndex, long startLeaf,
                                     long endLeaf, long position) {

  _CalcNode *travNode, *theChildNode;
  long nodeCount, f;
  _Parameter *fastIndex;
  node<long> *nodeChild;

  char *pastState = dsf->GetColumn(lastIndex),
       *thisState = dsf->GetColumn(index);

  for (nodeCount = startLeaf; nodeCount <= endLeaf; nodeCount++) {
    f = dsf->theNodeMap.lData[nodeCount];
    if (thisState[f] != pastState[f]) {
      travNode = (_CalcNode *)(((BaseRef *)flatCLeaves.lData)[nodeCount]);
      long *cCache = dsf->conversionCache.lData + (thisState[f] - 40) * 5;
      travNode->theProbs[0] = *(cCache++);
      travNode->theProbs[1] = *(cCache++);
      travNode->theProbs[2] = *(cCache++);
      travNode->theProbs[3] = *(cCache++);
      travNode->lastState = *cCache;
      theChildNode = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[
          ((node<long> *)(flatLeaves.lData[nodeCount]))->parent->in_object]);
      if (theChildNode->cBase > 0) {
        theChildNode->cBase = -1;
      }
    }
  }

  f = topLevelNodes.lLength - 1;
  long mm = topLevelNodes.lData[f], startINode = leftiNodes.lData[startLeaf];

  for (long nI = 0; nI < f; nI++, mm >>= 1) {
    //if (1) {
    if (mm & 1) {
      // marked for recalculation
      long startAt = 0;

      if (nI) {
        startAt = topLevelNodes.lData[nI - 1] + 1;
      }

      if (startAt < startINode) {
        startAt = startINode;
      }

      for (nodeCount = startAt; nodeCount <= topLevelNodes.lData[nI];
           nodeCount++) {
        theChildNode = (_CalcNode *)(((BaseRef *)flatTree.lData)[nodeCount]);
        if (theChildNode->cBase == -1) { // marked for update
          nodeChild = ((node<long> *)(flatNodes.lData[nodeCount]));
          theChildNode->theProbs[0] = LIKELIHOOD_SCALER;
          theChildNode->theProbs[1] = LIKELIHOOD_SCALER;
          theChildNode->theProbs[2] = LIKELIHOOD_SCALER;
          theChildNode->theProbs[3] = LIKELIHOOD_SCALER;
          theChildNode->cBase = 4;

          for (long k = 0; k < nodeChild->nodes.length; k++) {
            travNode = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[
                nodeChild->nodes.data[k]->in_object]);
            fastIndex = travNode->compExp->theData;
            if (travNode->lastState < 0) {
              _Parameter a = *(travNode->theProbs), b = travNode->theProbs[1],
                         c = travNode->theProbs[2], d = travNode->theProbs[3];

              _Parameter tmp = a * *(fastIndex++);
              tmp += b * *(fastIndex++);
              tmp += c * *(fastIndex++);
              theChildNode->theProbs[0] *= tmp + d * *(fastIndex++);

              tmp = a * *(fastIndex++);
              tmp += b * *(fastIndex++);
              tmp += c * *(fastIndex++);
              theChildNode->theProbs[1] *= tmp + d * *(fastIndex++);

              tmp = a * *(fastIndex++);
              tmp += b * *(fastIndex++);
              tmp += c * *(fastIndex++);
              theChildNode->theProbs[2] *= tmp + d * *(fastIndex++);

              tmp = a * *(fastIndex++);
              tmp += b * *(fastIndex++);
              tmp += c * *(fastIndex++);
              theChildNode->theProbs[3] *= tmp + d * *(fastIndex++);
            } else {
              fastIndex += travNode->lastState;
              theChildNode->theProbs[0] *= *fastIndex;
              theChildNode->theProbs[1] *= fastIndex[4];
              theChildNode->theProbs[2] *= fastIndex[8];
              theChildNode->theProbs[3] *= fastIndex[12];
            }
          }
          ((_CalcNode *)(
              (BaseRef *)variablePtrs.lData)[nodeChild->parent->in_object])
              ->cBase = -1;
        }
      }
      theChildNode =
          (_CalcNode *)(((BaseRef *)flatTree.lData)[topLevelNodes.lData[nI]]);
      _Parameter *nodeProbs = theChildNode->theProbs;
      fastIndex = rootIChildrenCache +
                  position * 4 * (topLevelNodes.lLength - 1) + 4 * nI;
      fastIndex[0] = nodeProbs[0];
      fastIndex[1] = nodeProbs[1];
      fastIndex[2] = nodeProbs[2];
      fastIndex[3] = nodeProbs[3];
    } else {
      _Parameter *nodeProbs =
          ((_CalcNode *)(((BaseRef *)flatTree.lData)[topLevelNodes.lData[nI]]))
              ->theProbs;

      fastIndex = rootIChildrenCache +
                  position * 4 * (topLevelNodes.lLength - 1) + 4 * nI;

      nodeProbs[0] = fastIndex[0];
      nodeProbs[1] = fastIndex[1];
      nodeProbs[2] = fastIndex[2];
      nodeProbs[3] = fastIndex[3];
      nodeChild = ((node<long> *)(flatNodes.lData[topLevelNodes.lData[nI]]))->parent;
      ((_CalcNode *)((BaseRef *)variablePtrs.lData)[nodeChild->in_object])->cBase = -1;
    }
  }

  //for (nodeCount=startINode; nodeCount<flatTree.lLength; nodeCount++) {
  for (nodeCount = topLevelNodes.lData[f - 1] + 1; nodeCount < flatTree.lLength;
       nodeCount++) {

    theChildNode = (_CalcNode *)(((BaseRef *)flatTree.lData)[nodeCount]);

    if (theChildNode->cBase == -1) { // marked for update
      nodeChild = ((node<long> *)(flatNodes.lData[nodeCount]));
      theChildNode->theProbs[0] = LIKELIHOOD_SCALER;
      theChildNode->theProbs[1] = LIKELIHOOD_SCALER;
      theChildNode->theProbs[2] = LIKELIHOOD_SCALER;
      theChildNode->theProbs[3] = LIKELIHOOD_SCALER;
      theChildNode->cBase = 4;

      for (long k = 0; k < nodeChild->nodes.length; k++) {
        travNode = ((_CalcNode *)((
            BaseRef *)variablePtrs.lData)[nodeChild->nodes.data[k]->in_object]);
        fastIndex = travNode->compExp->theData;
        if (travNode->lastState < 0) {
          _Parameter a = *(travNode->theProbs), b = travNode->theProbs[1],
                     c = travNode->theProbs[2], d = travNode->theProbs[3];

          _Parameter tmp = a * *(fastIndex++);
          tmp += b * *(fastIndex++);
          tmp += c * *(fastIndex++);
          theChildNode->theProbs[0] *= tmp + d * *(fastIndex++);

          tmp = a * *(fastIndex++);
          tmp += b * *(fastIndex++);
          tmp += c * *(fastIndex++);
          theChildNode->theProbs[1] *= tmp + d * *(fastIndex++);

          tmp = a * *(fastIndex++);
          tmp += b * *(fastIndex++);
          tmp += c * *(fastIndex++);
          theChildNode->theProbs[2] *= tmp + d * *(fastIndex++);

          tmp = a * *(fastIndex++);
          tmp += b * *(fastIndex++);
          tmp += c * *(fastIndex++);
          theChildNode->theProbs[3] *= tmp + d * *(fastIndex++);
        } else {
          fastIndex += travNode->lastState;
          theChildNode->theProbs[0] *= *fastIndex;
          theChildNode->theProbs[1] *= fastIndex[4];
          theChildNode->theProbs[2] *= fastIndex[8];
          theChildNode->theProbs[3] *= fastIndex[12];
        }
      }
      nodeChild = nodeChild->parent;
      if (nodeChild) {
        ((_CalcNode *)((BaseRef *)variablePtrs.lData)[nodeChild->in_object])
            ->cBase = -1;
      }
    }
  }

  theChildNode->cBase = 4;

  _Parameter result = theProbs[0] * theChildNode->theProbs[0] +
                      theProbs[1] * theChildNode->theProbs[1] +
                      theProbs[2] * theChildNode->theProbs[2] +
                      theProbs[3] * theChildNode->theProbs[3];

  if (result <= 0.0) {
    return ALMOST_ZERO;
  }

  return result;
}

//______________________________________________________________________________
#if USE_SCALING_TO_FIX_UNDERFLOW
_Parameter _TheTree::ThreadReleafTreeChar4(_DataSetFilter *dsf, long index,
                                           long lastIndex, long startLeaf,
                                           long endLeaf, long position,
                                           long offset, long fixAttempt,
                                           _Parameter localScalingFactor)
#else
    _Parameter _TheTree::ThreadReleafTreeChar4(_DataSetFilter *dsf, long index,
                                               long lastIndex, long startLeaf,
                                               long endLeaf, long position,
                                               long offset)
#endif
    {

  _CalcNode *travNode, *theChildNode;

  long nodeCount,
      f = topLevelNodes.lLength - 1,
      *ns = nodeStates + (flatLeaves.lLength + flatNodes.lLength) * offset;

  _Parameter *fastIndex,
      *mlc = marginalLikelihoodCache +
             cBase * (flatLeaves.lLength + flatNodes.lLength) * offset;

  node<long> *nodeChild;

  char *pastState, *thisState = dsf->GetColumn(index),
                   *nm = nodeMarkers + (flatNodes.lLength) * offset;

  if (lastIndex >= 0) {
    pastState = dsf->GetColumn(lastIndex);
  }

#if USE_SCALING_TO_FIX_UNDERFLOW
  if (lastIndex >= 0 && scalingForUnderflow->theData[index] !=
                            scalingForUnderflow->theData[lastIndex]) {
    lastIndex = -1;
    startLeaf = 0;
    endLeaf = flatLeaves.lLength - 1;
  }

  _Parameter scalingFactor = scalingForUnderflow->theData[index];
  bool overflowFlag = fixAttempt < 0;
  if (overflowFlag) {
    fixAttempt = -fixAttempt;
  }
#endif

  long mm = topLevelNodes.lData[f], startINode = leftiNodes.lData[startLeaf],
       fromL = startLeaf, toL = 0;

  // mod SLKP 20070209 to skip over constant parts of the tree

  for (long tlc = 0; tlc < topLevelNodes.lLength; tlc++, mm >>= 1) {
    if (mm & 1) {
      toL = topLevelRightL.lData[tlc];
    } else {
      toL = topLevelLeftL.lData[tlc] - 1;
    }

    if (fromL <= toL) {
      long l1 = MAX(startLeaf, fromL), l2 = MIN(endLeaf, toL);

      for (nodeCount = l1; nodeCount <= l2; nodeCount++) {
        long f2 = dsf->theNodeMap.lData[nodeCount];
        if (lastIndex == -1 || thisState[f2] != pastState[f2]
#if USE_SCALING_TO_FIX_UNDERFLOW
            || scalingFactor != scalingForUnderflow->theData[lastIndex]
#endif
            ) {
          travNode = (_CalcNode *)(((BaseRef *)flatCLeaves.lData)[nodeCount]);
          long *cCache = dsf->conversionCache.lData + (thisState[f2] - 40) * 5;

          _Parameter *fi = mlc + nodeCount * 4;

          fi[0] = cCache[0];
          fi[1] = cCache[1];
          fi[2] = cCache[2];
          fi[3] = cCache[3];

#if USE_SCALING_TO_FIX_UNDERFLOW
          if (scalingFactor != 1.0) {
            fi[0] *= scalingFactor;
            fi[1] *= scalingFactor;
            fi[2] *= scalingFactor;
            fi[3] *= scalingFactor;
          }
#endif

          ns[nodeCount] = cCache[4];

          theChildNode = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[(
              (node<long> *)(flatLeaves.lData[nodeCount]))->parent->in_object]);

          f2 = theChildNode->nodeIndex - flatLeaves.lLength;

          if (nm[f2] >= 0) {
            nm[f2] = -1;
          }
        }
      }
    }

    fromL = topLevelRightL.lData[tlc] + 1;
    if (fromL > endLeaf) {
      break;
    }
  }

  // end 20070209 mod
  mm = topLevelNodes.lData[f];

  for (long nI = 0; nI < f; nI++, mm >>= 1) {
    if (mm & 1) {
      // marked for recalculation
      long startAt = 0;

      if (nI) {
        startAt = topLevelNodes.lData[nI - 1] + 1;
      }

      if (startAt < startINode) {
        startAt = startINode;
      }

      for (nodeCount = startAt; nodeCount <= topLevelNodes.lData[nI];
           nodeCount++) {
        if (nm[nodeCount] == -1) { // marked for update
          nodeChild = ((node<long> *)(flatNodes.lData[nodeCount]));

          _Parameter *theseProbs = mlc + 4 * (nodeCount + flatLeaves.lLength);
          theseProbs[0] = LIKELIHOOD_SCALER_INT;
          theseProbs[1] = LIKELIHOOD_SCALER_INT;
          theseProbs[2] = LIKELIHOOD_SCALER_INT;
          theseProbs[3] = LIKELIHOOD_SCALER_INT;
          nm[nodeCount] = 0;

          for (long k = 0; k < nodeChild->nodes.length; k++) {
            travNode = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[
                nodeChild->nodes.data[k]->in_object]);

            _Parameter *childProbs = mlc + 4 * (travNode->nodeIndex);
            long nState = ns[travNode->nodeIndex];

            fastIndex = travNode->compExp->theData;
            if (nState < 0) {
              //marginalLFEvalsAmb ++;

              _Parameter a = *(childProbs), b = childProbs[1],
                         c = childProbs[2], d = childProbs[3];

              _Parameter tmp = a * fastIndex[0] + b * fastIndex[1] +
                               c * fastIndex[2] + d * fastIndex[3],
                         tmp2 = a * fastIndex[4] + b * fastIndex[5] +
                                c * fastIndex[6] + d * fastIndex[7],
                         tmp3 = a * fastIndex[8] + b * fastIndex[9] +
                                c * fastIndex[10] + d * fastIndex[11],
                         tmp4 = a * fastIndex[12] + b * fastIndex[13] +
                                c * fastIndex[14] + d * fastIndex[15];

              theseProbs[0] *= tmp;
              theseProbs[1] *= tmp2;
              theseProbs[2] *= tmp3;
              theseProbs[3] *= tmp4;

            } else {
              //marginalLFEvals ++;
              fastIndex += nState;
#if USE_SCALING_TO_FIX_UNDERFLOW
              theseProbs[0] *= fastIndex[0] * scalingFactor;
              theseProbs[1] *= fastIndex[4] * scalingFactor;
              theseProbs[2] *= fastIndex[8] * scalingFactor;
              theseProbs[3] *= fastIndex[12] * scalingFactor;

#else
              theseProbs[0] *= fastIndex[0];
              theseProbs[1] *= fastIndex[4];
              theseProbs[2] *= fastIndex[8];
              theseProbs[3] *= fastIndex[12];
#endif
            }
          }
          nm[((_CalcNode *)(
              (BaseRef *)variablePtrs.lData)[nodeChild->parent->in_object])
                 ->nodeIndex - flatLeaves.lLength] = -1;

#if USE_SCALING_TO_FIX_UNDERFLOW
          if (theseProbs[0] == 0.0 && theseProbs[1] == 0.0 &&
              theseProbs[2] == 0.0 && theseProbs[3] == 0.0) {
            return doChar4Scaling(dsf, index, position, offset, fixAttempt,
                                  localScalingFactor, overflowFlag, false);
          } else if (theseProbs[0] >= 1e500 || theseProbs[1] >= 1e500 ||
                     theseProbs[2] >= 1e500 || theseProbs[3] >= 1e500) {
            return doChar4Scaling(dsf, index, position, offset, fixAttempt,
                                  localScalingFactor, overflowFlag, true);
          }
#endif

        }
      }

      theChildNode =
          (_CalcNode *)(((BaseRef *)flatTree.lData)[topLevelNodes.lData[nI]]);

      _Parameter *nodeProbs = mlc + 4 * theChildNode->nodeIndex;

      fastIndex = rootIChildrenCache +
                  position * 4 * (topLevelNodes.lLength - 1) + 4 * nI;

      fastIndex[0] = nodeProbs[0];
      fastIndex[1] = nodeProbs[1];
      fastIndex[2] = nodeProbs[2];
      fastIndex[3] = nodeProbs[3];
    } else {
      _Parameter *nodeProbs =
          mlc + 4 * ((_CalcNode *)(
                        ((BaseRef *)flatTree.lData)[topLevelNodes.lData[nI]]))
                        ->nodeIndex;

      fastIndex = rootIChildrenCache +
                  position * 4 * (topLevelNodes.lLength - 1) + 4 * nI;

      nodeProbs[0] = fastIndex[0];
      nodeProbs[1] = fastIndex[1];
      nodeProbs[2] = fastIndex[2];
      nodeProbs[3] = fastIndex[3];

      nodeChild =
          ((node<long> *)(flatNodes.lData[topLevelNodes.lData[nI]]))->parent;

      nm[((_CalcNode *)((BaseRef *)variablePtrs.lData)[nodeChild->in_object])
             ->nodeIndex - flatLeaves.lLength] = -1;
    }
  }

  for (nodeCount = topLevelNodes.lData[f - 1] + 1; nodeCount < flatTree.lLength;
       nodeCount++) {

    if (nm[nodeCount] == -1) { // marked for update
      nodeChild = ((node<long> *)(flatNodes.lData[nodeCount]));

      _Parameter *theseProbs = mlc + 4 * (nodeCount + flatLeaves.lLength);
      theseProbs[0] = LIKELIHOOD_SCALER;
      theseProbs[1] = LIKELIHOOD_SCALER;
      theseProbs[2] = LIKELIHOOD_SCALER;
      theseProbs[3] = LIKELIHOOD_SCALER;
      nm[nodeCount] = 0;

      for (long k = 0; k < nodeChild->nodes.length; k++) {
        travNode = ((_CalcNode *)((
            BaseRef *)variablePtrs.lData)[nodeChild->nodes.data[k]->in_object]);

        _Parameter *childProbs = mlc + 4 * (travNode->nodeIndex);
        long nState = ns[travNode->nodeIndex];

        fastIndex = travNode->compExp->theData;
        if (nState < 0) {
          //marginalLFEvalsAmb++;

          _Parameter a = *(childProbs), b = childProbs[1], c = childProbs[2],
                     d = childProbs[3];

          _Parameter tmp = a * fastIndex[0] + b * fastIndex[1] +
                           c * fastIndex[2] + d * fastIndex[3],
                     tmp2 = a * fastIndex[4] + b * fastIndex[5] +
                            c * fastIndex[6] + d * fastIndex[7],
                     tmp3 = a * fastIndex[8] + b * fastIndex[9] +
                            c * fastIndex[10] + d * fastIndex[11],
                     tmp4 = a * fastIndex[12] + b * fastIndex[13] +
                            c * fastIndex[14] + d * fastIndex[15];

          theseProbs[0] *= tmp;
          theseProbs[1] *= tmp2;
          theseProbs[2] *= tmp3;
          theseProbs[3] *= tmp4;
        } else {
          //marginalLFEvals++;

          fastIndex += nState;
#if USE_SCALING_TO_FIX_UNDERFLOW
          theseProbs[0] *= fastIndex[0] * scalingFactor;
          theseProbs[1] *= fastIndex[4] * scalingFactor;
          theseProbs[2] *= fastIndex[8] * scalingFactor;
          theseProbs[3] *= fastIndex[12] * scalingFactor;

#else
          theseProbs[0] *= fastIndex[0];
          theseProbs[1] *= fastIndex[4];
          theseProbs[2] *= fastIndex[8];
          theseProbs[3] *= fastIndex[12];
#endif
        }
      }

#if USE_SCALING_TO_FIX_UNDERFLOW
      if (theseProbs[0] == 0.0 && theseProbs[1] == 0.0 &&
          theseProbs[2] == 0.0 && theseProbs[3] == 0.0) {
        return doChar4Scaling(dsf, index, position, offset, fixAttempt,
                              localScalingFactor, overflowFlag, false);
      } else if (theseProbs[0] >= 1e500 || theseProbs[1] >= 1e500 ||
                 theseProbs[2] >= 1e500 || theseProbs[3] >= 1e500) {
        return doChar4Scaling(dsf, index, position, offset, fixAttempt,
                              localScalingFactor, overflowFlag, true);
      }
#endif
      nodeChild = nodeChild->parent;
      if (nodeChild)
        nm[((_CalcNode *)((BaseRef *)variablePtrs.lData)[nodeChild->in_object])
               ->nodeIndex - flatLeaves.lLength] = -1;
    }
  }

  //theChildNode->cBase = 4;

  fastIndex = mlc + 4 * (flatLeaves.lLength + flatNodes.lLength - 1);

  _Parameter result = theProbs[0] * fastIndex[0] + theProbs[1] * fastIndex[1] +
                      theProbs[2] * fastIndex[2] + theProbs[3] * fastIndex[3];

#if USE_SCALING_TO_FIX_UNDERFLOW
  if (result <= 0.0) { // underflow
    return doChar4Scaling(dsf, index, position, offset, fixAttempt,
                          localScalingFactor, overflowFlag, false);
  } else if (result >= 1e500) { // overflow
    return doChar4Scaling(dsf, index, position, offset, fixAttempt,
                          localScalingFactor, overflowFlag, true);
  }

#else
  if (result <= 0.0) {
    return ALMOST_ZERO;
  }
#endif

  return result;
}

//______________________________________________________________________________
// set leaf value and re-prune w/o reexping
// for subsequent entries into a datafilter
#if USE_SCALING_TO_FIX_UNDERFLOW
_Parameter _TheTree::ReleafTreeChar4(_DataSetFilter *dsf, long index,
                                     long lastIndex, long startLeaf,
                                     long endLeaf, long fixAttempt,
                                     _Parameter localScalingFactor)
#else
    _Parameter _TheTree::ReleafTreeChar4(_DataSetFilter *dsf, long index,
                                         long lastIndex, long startLeaf,
                                         long endLeaf)
#endif

    {
  _CalcNode *travNode, *theChildNode;

  long nodeCount, f;

  _Parameter *fastIndex;
  node<long> *nodeChild;

  char *pastState = nil, *thisState = dsf->GetColumn(index);

  if (lastIndex >= 0) {
    pastState = dsf->GetColumn(lastIndex);
  }

#if USE_SCALING_TO_FIX_UNDERFLOW
  if (lastIndex >= 0 && scalingForUnderflow->theData[index] !=
                            scalingForUnderflow->theData[lastIndex]) {
    lastIndex = -1;
    startLeaf = 0;
    endLeaf = flatLeaves.lLength - 1;
  }

  _Parameter scalingFactor = scalingForUnderflow->theData[index];
  bool overflowFlag = fixAttempt < 0;

  if (overflowFlag) {
    fixAttempt = -fixAttempt;
  }
#endif

  for (nodeCount = startLeaf; nodeCount <= endLeaf; nodeCount++) {
    f = dsf->theNodeMap.lData[nodeCount];
    if (lastIndex < 0 || thisState[f] != pastState[f]
#if USE_SCALING_TO_FIX_UNDERFLOW
        || scalingFactor != scalingForUnderflow->theData[lastIndex]
#endif
        ) {
      travNode = (_CalcNode *)(((BaseRef *)flatCLeaves.lData)[nodeCount]);
      long *cCache = dsf->conversionCache.lData + (thisState[f] - 40) * 5;

#if USE_SCALING_TO_FIX_UNDERFLOW
      travNode->theProbs[0] = *(cCache++) * scalingFactor;
      travNode->theProbs[1] = *(cCache++) * scalingFactor;
      travNode->theProbs[2] = *(cCache++) * scalingFactor;
      travNode->theProbs[3] = *(cCache++) * scalingFactor;

#else
      travNode->theProbs[0] = *(cCache++);
      travNode->theProbs[1] = *(cCache++);
      travNode->theProbs[2] = *(cCache++);
      travNode->theProbs[3] = *(cCache++);
#endif

      travNode->lastState = *cCache;

      theChildNode = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[
          ((node<long> *)(flatLeaves.lData[nodeCount]))->parent->in_object]);
      if (theChildNode->cBase > 0) {
        theChildNode->cBase = -1;
      }
    }
  }

  for (nodeCount = leftiNodes.lData[startLeaf]; nodeCount < flatTree.lLength;
       nodeCount++) {
    theChildNode = (_CalcNode *)(((BaseRef *)flatTree.lData)[nodeCount]);
    if (theChildNode->cBase == -1) { // marked for update
      nodeChild = ((node<long> *)(flatNodes.lData[nodeCount]));
      theChildNode->theProbs[0] = LIKELIHOOD_SCALER;
      theChildNode->theProbs[1] = LIKELIHOOD_SCALER;
      theChildNode->theProbs[2] = LIKELIHOOD_SCALER;
      theChildNode->theProbs[3] = LIKELIHOOD_SCALER;
      theChildNode->cBase = 4;

      for (long k = 0; k < nodeChild->nodes.length; k++) {
        travNode = ((_CalcNode *)((
            BaseRef *)variablePtrs.lData)[nodeChild->nodes.data[k]->in_object]);
        fastIndex = travNode->compExp->theData;
        if (travNode->lastState < 0) {
          _Parameter a = *(travNode->theProbs), b = travNode->theProbs[1],
                     c = travNode->theProbs[2], d = travNode->theProbs[3];

          _Parameter tmp = a * *(fastIndex++);
          tmp += b * *(fastIndex++);
          tmp += c * *(fastIndex++);
          theChildNode->theProbs[0] *= tmp + d * *(fastIndex++);

          tmp = a * *(fastIndex++);
          tmp += b * *(fastIndex++);
          tmp += c * *(fastIndex++);
          theChildNode->theProbs[1] *= tmp + d * *(fastIndex++);

          tmp = a * *(fastIndex++);
          tmp += b * *(fastIndex++);
          tmp += c * *(fastIndex++);
          theChildNode->theProbs[2] *= tmp + d * *(fastIndex++);

          tmp = a * *(fastIndex++);
          tmp += b * *(fastIndex++);
          tmp += c * *(fastIndex++);
          theChildNode->theProbs[3] *= tmp + d * *(fastIndex++);
        } else {
          fastIndex += travNode->lastState;
#if USE_SCALING_TO_FIX_UNDERFLOW
          theChildNode->theProbs[0] *= *fastIndex * scalingFactor;
          theChildNode->theProbs[1] *= fastIndex[4] * scalingFactor;
          theChildNode->theProbs[2] *= fastIndex[8] * scalingFactor;
          theChildNode->theProbs[3] *= fastIndex[12] * scalingFactor;
#else
          theChildNode->theProbs[0] *= *fastIndex;
          theChildNode->theProbs[1] *= fastIndex[4];
          theChildNode->theProbs[2] *= fastIndex[8];
          theChildNode->theProbs[3] *= fastIndex[12];
#endif
        }
      }
      travNode->cBase = 4;
      nodeChild = nodeChild->parent;
      if (nodeChild) {
        ((_CalcNode *)((BaseRef *)variablePtrs.lData)[nodeChild->in_object])
            ->cBase = -1;
      }

#if USE_SCALING_TO_FIX_UNDERFLOW
      if (theChildNode->theProbs[0] == 0.0 &&
          theChildNode->theProbs[1] == 0.0 &&
          theChildNode->theProbs[2] == 0.0 &&
          theChildNode->theProbs[3] == 0.0) {
        return doChar4Scaling_nc(dsf, index, fixAttempt, localScalingFactor,
                                 overflowFlag, false);
      } else if (theChildNode->theProbs[0] == HUGE_VAL ||
                 theChildNode->theProbs[1] == HUGE_VAL ||
                 theChildNode->theProbs[2] == HUGE_VAL ||
                 theChildNode->theProbs[3] == HUGE_VAL) {
        return doChar4Scaling_nc(dsf, index, fixAttempt, localScalingFactor,
                                 overflowFlag, true);
      }
#endif

    }
  }

  theChildNode->cBase = 4;
  _Parameter result = theProbs[0] * theChildNode->theProbs[0] +
                      theProbs[1] * theChildNode->theProbs[1] +
                      theProbs[2] * theChildNode->theProbs[2] +
                      theProbs[3] * theChildNode->theProbs[3];

#if USE_SCALING_TO_FIX_UNDERFLOW
  if (result == 0.0) {
    return doChar4Scaling_nc(dsf, index, fixAttempt, localScalingFactor,
                             overflowFlag, false);
  } else if (result == HUGE_VAL) {
    return doChar4Scaling_nc(dsf, index, fixAttempt, localScalingFactor,
                             overflowFlag, true);
  }
#else
  if (result <= 0.0) {
    return ALMOST_ZERO;
  }
#endif
  return result;
}

//______________________________________________________________________________
// set leaf value and re-prune w/o reexping
// for subsequent entries into a datafilter
_Parameter _TheTree::ReleafTreeChar4Degenerate(_DataSetFilter *dsf, long index) {

  _CalcNode *
  travNode = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[theRoot->in_object]),
 *theChildNode = ((_CalcNode *)(
     (BaseRef *)variablePtrs.lData)[theRoot->nodes.data[0]->in_object]);

  char *thisState = dsf->GetColumn(index);
  long *rootState = dsf->conversionCache.lData +
                    (thisState[dsf->theNodeMap.lData[0]] - 40) * 5,
       *tipState = dsf->conversionCache.lData +
                   (thisState[dsf->theNodeMap.lData[1]] - 40) * 5,
       nodeCount = rootState[4], f = tipState[4];

  _Matrix *theTP = theChildNode->GetCompExp();

  _Parameter reslt, tmp, *fastIdx = theTP->fastIndex(), *nodeProbs;

  if ((nodeCount >= 0) && (f >= 0)) { // 1-1
    reslt = fastIdx[nodeCount * 4 + f] * theProbs[nodeCount];
  } else if (nodeCount >= 0) { // 1-many
    fastIdx += nodeCount * cBase;

    tmp = *fastIdx * *tipState;
    tmp += fastIdx[1] * tipState[1];
    tmp += fastIdx[2] * tipState[2];
    tmp += fastIdx[3] * tipState[3];

    reslt = theProbs[nodeCount] * tmp;
  } else if (f >= 0) { // many to 1
    fastIdx += f;
    nodeProbs = travNode->theProbs;
    *nodeProbs = *fastIdx * *rootState;
    nodeProbs[1] = fastIdx[4] * rootState[1];
    nodeProbs[2] = fastIdx[8] * rootState[2];
    nodeProbs[3] = fastIdx[12] * rootState[3];

    reslt = *theProbs * *nodeProbs + theProbs[1] * nodeProbs[1] +
            theProbs[2] * nodeProbs[2] + theProbs[3] * nodeProbs[3];

  } else { // many to many
    nodeProbs = travNode->theProbs;
    nodeProbs[0] =
        (fastIdx[0] * tipState[0] + fastIdx[1] * tipState[1] +
         fastIdx[2] * tipState[2] + fastIdx[3] * tipState[3]) * rootState[0];
    nodeProbs[1] =
        (fastIdx[4] * tipState[0] + fastIdx[5] * tipState[1] +
         fastIdx[6] * tipState[2] + fastIdx[7] * tipState[3]) * rootState[1];
    nodeProbs[2] =
        (fastIdx[8] * tipState[0] + fastIdx[9] * tipState[1] +
         fastIdx[10] * tipState[2] + fastIdx[11] * tipState[3]) * rootState[2];
    nodeProbs[3] =
        (fastIdx[12] * tipState[0] + fastIdx[13] * tipState[1] +
         fastIdx[14] * tipState[2] + fastIdx[15] * tipState[3]) * rootState[3];

    reslt = *theProbs * *nodeProbs + theProbs[1] * nodeProbs[1] +
            theProbs[2] * nodeProbs[2] + theProbs[3] * nodeProbs[3];
  }
  if (reslt <= 0.0) {
    reslt = ALMOST_ZERO;
  }
  //printf ("%d\t%g\n", index, reslt);
  return reslt;
}

//______________________________________________________________________________
// set leaf value and re-prune w/o reexping
// for subsequent entries into a datafilter
_Parameter _TheTree::ThreadReleafTreeCharCache(_DataSetFilter *dsf, long index,
                                               long lastIndex, long startLeaf,
                                               long endLeaf, long position,
                                               long offset) {

  _CalcNode *travNode, *theChildNode;

  long nodeCount, f,
      *ns = nodeStates + (flatLeaves.lLength + flatNodes.lLength) * offset;

  _Parameter *fastIndex, *stopper,
      *mlc = marginalLikelihoodCache +
             cBase * (flatLeaves.lLength + flatNodes.lLength) * offset;

  node<long> *nodeChild;

  char *pastState = nil, *thisState = dsf->GetColumn(index),
       *nm = nodeMarkers + flatNodes.lLength * offset;

  if (lastIndex >= 0) {
    pastState = dsf->GetColumn(lastIndex);
  }

  for (nodeCount = startLeaf; nodeCount <= endLeaf; nodeCount++) {
    f = dsf->theNodeMap.lData[nodeCount];
    if ((lastIndex == -1) || (thisState[f] != pastState[f])) {
      ns[nodeCount] =
          dsf->LookupConversion(thisState[f], mlc + cBase * nodeCount);
      theChildNode = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[
          ((node<long> *)(flatLeaves.lData[nodeCount]))->parent->in_object]);
      nm[theChildNode->nodeIndex - flatLeaves.lLength] = -1;
    }
  }

  long lm1 = topLevelNodes.lLength - 1, mm = topLevelNodes.lData[lm1],
       startINode = leftiNodes.lData[startLeaf];

  for (long nI = 0; nI < lm1; nI++, mm >>= 1) {
    if (mm & 1) {
      // marked for recalculation
      long startAt = 0;

      if (nI) {
        startAt = topLevelNodes.lData[nI - 1] + 1;
      }

      if (startAt < startINode) {
        startAt = startINode;
      }

      for (nodeCount = startAt; nodeCount <= topLevelNodes.lData[nI];
           nodeCount++) {
        if (nm[nodeCount] == -1) { // marked for update
          _Parameter *theseProbs = fastIndex =
              mlc + cBase * (nodeCount + flatLeaves.lLength);
          stopper = fastIndex + cBase;

          for (; fastIndex != stopper; *(fastIndex++) = 1.)
            ;

          nodeChild = ((node<long> *)(flatNodes.lData[nodeCount]));

          nm[nodeCount] = 0;

          for (long k = 0; k < nodeChild->nodes.length; k++) {
            travNode = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[
                nodeChild->nodes.data[k]->in_object]);

            _Parameter *cProbs = mlc + cBase * travNode->nodeIndex,
                       *stopper2 = cProbs + cBase;
            fastIndex = travNode->compExp->theData;

            f = ns[travNode->nodeIndex];
            if (f < 0) {
              for (long i = 0; i < cBase; i++) {
                _Parameter tmp = *(cProbs) * *fastIndex;
                fastIndex++;
                stopper = cProbs + 1;
                for (; stopper != stopper2; fastIndex++, stopper++) {
                  tmp += *stopper * (*fastIndex);
                }
                theseProbs[i] *= tmp;
              }
            } else {
              fastIndex += f;
              _Parameter tmp = cProbs[f];
              for (long i = 0; i < cBase; i++, fastIndex += cBase) {
                theseProbs[i] *= tmp * (*fastIndex);
              }
            }
          }
          nodeChild = nodeChild->parent;
          if (nodeChild) {
            nm[((_CalcNode *)((BaseRef *)variablePtrs.lData)[
                nodeChild->in_object])->nodeIndex - flatLeaves.lLength] = -1;
          }
        }
      }

      theChildNode =
          (_CalcNode *)(((BaseRef *)flatTree.lData)[topLevelNodes.lData[nI]]);
      _Parameter *nodeProbs = mlc + cBase * theChildNode->nodeIndex;
      fastIndex = rootIChildrenCache +
                  position * cBase * (topLevelNodes.lLength - 1) + cBase * nI;
      for (long i = 0; i < cBase; i++) {
        fastIndex[i] = nodeProbs[i];
      }
    } else {
      _Parameter *nodeProbs =
          mlc +
          cBase * ((_CalcNode *)(
                      ((BaseRef *)flatTree.lData)[topLevelNodes.lData[nI]]))
                      ->nodeIndex;
      fastIndex = rootIChildrenCache +
                  position * cBase * (topLevelNodes.lLength - 1) + cBase * nI;
      for (long i = 0; i < cBase; i++) {
        nodeProbs[i] = fastIndex[i];
      }
      nodeChild =
          ((node<long> *)(flatNodes.lData[topLevelNodes.lData[nI]]))->parent;
      nm[((_CalcNode *)((BaseRef *)variablePtrs.lData)[nodeChild->in_object])
             ->nodeIndex - flatLeaves.lLength] = -1;
    }
  }

  for (nodeCount = topLevelNodes.lData[lm1 - 1] + 1;
       nodeCount < flatTree.lLength; nodeCount++) {
    if (nm[nodeCount] == -1) { // marked for update
      _Parameter *theseProbs = fastIndex =
          mlc + cBase * (nodeCount + flatLeaves.lLength);

      for (long zeroPop = 0; zeroPop < cBase; zeroPop++) {
        fastIndex[zeroPop] = 1.;
      }

      nodeChild = ((node<long> *)(flatNodes.lData[nodeCount]));

      nm[nodeCount] = 0;

      for (long k = 0; k < nodeChild->nodes.length; k++) {
        travNode = ((_CalcNode *)((
            BaseRef *)variablePtrs.lData)[nodeChild->nodes.data[k]->in_object]);

        _Parameter *cProbs = mlc + cBase * travNode->nodeIndex,
                   *stopper2 = cProbs + cBase;

        fastIndex = travNode->compExp->theData;

        f = ns[travNode->nodeIndex];
        if (f < 0) {
          for (long i = 0; i < cBase; i++) {
            _Parameter tmp = *(cProbs) * *fastIndex;
            fastIndex++;
            stopper = cProbs + 1;
            for (; stopper != stopper2; fastIndex++, stopper++) {
              tmp += *stopper * (*fastIndex);
            }
            theseProbs[i] *= tmp;
          }
          /*if (rem==1)
                    {
                        for (long i=0; i<cBase; i++)
                        {
                            _Parameter tmp1  = 0.0,
                                       tmp2  = 0.0,
                                       tmp3  = 0.0,
                                       tmp4  = 0.0;

                            for (long k=0; k<upTo; k+=4)
                            {
                                tmp1 += cProbs[k]  *fastIndex[k];
                                tmp2 += cProbs[k+1]*fastIndex[k+1];
                                tmp3 += cProbs[k+2]*fastIndex[k+2];
                                tmp4 += cProbs[k+3]*fastIndex[k+3];
                            }

                            theseProbs[i] *= tmp1 +
          cProbs[upTo]*fastIndex[upTo] + tmp2+tmp3+tmp4;
                            fastIndex += cBase;
                        }
                    }
                    else
                        if (rem==2)
                        {
                            for (long i=0; i<cBase; i++)
                            {
                                _Parameter tmp1  = 0.0,
                                           tmp2  = 0.0,
                                           tmp3  = 0.0,
                                           tmp4  = 0.0;

                                for (long k=0; k<upTo; k+=4)
                                {
                                    tmp1 += cProbs[k]  *fastIndex[k];
                                    tmp2 +=
          cProbs[k+1]*fastIndex[k+1];
                                    tmp3 +=
          cProbs[k+2]*fastIndex[k+2];
                                    tmp4 +=
          cProbs[k+3]*fastIndex[k+3];
                                }

                                theseProbs[i] *= tmp1+tmp2+tmp3+tmp4+
                                                 +
          cProbs[upTo]*fastIndex[upTo] +
                                                 +
          cProbs[upTo+1]*fastIndex[upTo+1];
                                fastIndex += cBase;
                            }
                        }
                        else
                            if (rem==3)
                                for (long i=0; i<cBase; i++)
                                {
                                    _Parameter tmp1  = 0.0,
                                               tmp2  = 0.0,
                                               tmp3  = 0.0,
                                               tmp4  = 0.0;

                                    for (long k=0; k<upTo; k+=4)
                                    {
                                        tmp1 += cProbs[k]
          *fastIndex[k];
                                        tmp2 +=
          cProbs[k+1]*fastIndex[k+1];
                                        tmp3 +=
          cProbs[k+2]*fastIndex[k+2];
                                        tmp4 +=
          cProbs[k+3]*fastIndex[k+3];
                                    }

                                    theseProbs[i] *=
          tmp1+tmp2+tmp3+tmp4+
                                                     +
          cProbs[upTo]*fastIndex[upTo] +
                                                     +
          cProbs[upTo+1]*fastIndex[upTo+1];
                                                     +
          cProbs[upTo+2]*fastIndex[upTo+2];

                                    fastIndex += cBase;
                                }
                            else
                                for (long i=0; i<cBase; i++)
                                {
                                    _Parameter tmp1  = 0.0,
                                               tmp2  = 0.0,
                                               tmp3  = 0.0,
                                               tmp4  = 0.0;

                                    for (long k=0; k<upTo; k+=4)
                                    {
                                        tmp1 += cProbs[k]
          *fastIndex[k];
                                        tmp2 +=
          cProbs[k+1]*fastIndex[k+1];
                                        tmp3 +=
          cProbs[k+2]*fastIndex[k+2];
                                        tmp4 +=
          cProbs[k+3]*fastIndex[k+3];
                                    }

                                    theseProbs[i] *=
          tmp1+tmp2+tmp3+tmp4;
                                    fastIndex     += cBase;
                                }       */

        } else {
          fastIndex += f;
          _Parameter tmp = cProbs[f];
          for (long i = 0; i < cBase; i++, fastIndex += cBase) {
            theseProbs[i] *= tmp * (*fastIndex);
          }
        }
      }
      nodeChild = nodeChild->parent;
      if (nodeChild) {
        nm[((_CalcNode *)((BaseRef *)variablePtrs.lData)[nodeChild->in_object])
               ->nodeIndex - flatLeaves.lLength] = -1;
      }
    }
  }

  _Parameter result = 0.0;
  fastIndex = mlc + cBase * (flatLeaves.lLength + flatTree.lLength - 1);

  for (long i = 0; i < cBase; i++) {
    result += theProbs[i] * fastIndex[i];
  }

  if (result <= 0.0) {
    return ALMOST_ZERO;
  }
  return result;
}

//______________________________________________________________________________
// assign proper values to leaf conditional probability vectors
bool _TheTree::IntPopulateLeaves(_DataSetFilter *dsf, long index, long) {

  bool allGaps = true;

  for (long nodeCount = 0; nodeCount < flatLeaves.lLength; nodeCount++) {
    //if (lastIndex<0 || !dsf->CompareTwoSites(lastIndex,index,nodeCount)) {
    _CalcNode *travNode =
        (_CalcNode *)(((BaseRef *)flatCLeaves.lData)[nodeCount]);
    allGaps &= ((travNode->lastState = dsf->Translate2Frequencies(
        (*dsf)(index, nodeCount), travNode->theProbs, true)) < 0);
    if (allGaps) // check to see if all
      for (long b = 0; b < cBase; b++)
        if (travNode->theProbs[b] == 0.0) {
          allGaps = false;
          break;
        }

    _CalcNode *theChildNode = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[
        ((node<long> *)(flatLeaves.lData[nodeCount]))->parent->in_object]);
    theChildNode->cBase = -1;
  }

  return allGaps;
}

//______________________________________________________________________________
// assume current values of all parameters
// return 2 sets of vectors for each branch
//   - top-down  conditionals
//   - bottom-up conditionals
//   resultMatrix is assumed to contain
//      uniqueSites X (flatLeaves.lLength+flatTree.lLength)*cBase*2 X
// categoryCount
void _TheTree::RecoverNodeSupportStates(_DataSetFilter *dsf, long index,
                                        long lastIndex, _Matrix &resultMatrix) {

  long globalShifter = (flatLeaves.lLength + flatTree.lLength) * cBase,
       catShifer = dsf->NumberDistinctSites() * 2 * globalShifter;

  IntPopulateLeaves(dsf, index, lastIndex);

  /* pass 1; populate top-down vectors */
  /* ugly top-bottom algorithm for debuggability and compactness */

  for (long catCount = 0; catCount < categoryCount; catCount++) {
    _Parameter *currentStateVector =
                   resultMatrix.theData + 2 * globalShifter * index +
                   catShifer * catCount,
               *vecPointer = currentStateVector;

    for (long nodeCount = 0; nodeCount < flatCLeaves.lLength; nodeCount++) {
      _Parameter *leafVec =
          ((_CalcNode *)(((BaseRef *)flatCLeaves.lData)[nodeCount]))->theProbs;

      for (long cc = 0; cc < cBase; cc++, vecPointer++) {
        *vecPointer = leafVec[cc];
      }
    }

    for (long iNodeCount = 0; iNodeCount < flatTree.lLength; iNodeCount++) {

      node<long> *thisINode = (node<long> *)flatNodes.lData[iNodeCount];

      for (long cc = 0; cc < cBase; cc++, vecPointer++) {
        _Parameter tmp = 1.0;
        for (long nc = 0; nc < thisINode->nodes.length; nc++) {
          _Parameter tmp2 = 0.0;
          _CalcNode *child =
              ((_CalcNode *)((BaseRef *)variablePtrs
                                 .lData)[thisINode->nodes.data[nc]->in_object]);
          _Parameter *childSupport =
                         currentStateVector + child->nodeIndex * cBase,
                     *transMatrix =
                         child->GetCompExp(categoryCount > 1 ? catCount : (-1))
                             ->theData + cc * cBase;

          for (long cc2 = 0; cc2 < cBase; cc2++) {
            tmp2 += transMatrix[cc2] * childSupport[cc2];
          }

          tmp *= tmp2;
        }
        *vecPointer = tmp;
      }
    }
    RecoverNodeSupportStates2(&GetRoot(), currentStateVector + globalShifter,
                              currentStateVector,
                              categoryCount > 1 ? catCount : (-1));
  }

   //pass 2; populate bottom-up vectors 
   //for this we need to traverse the tree pre-order 
   //because speed is not much of a concern, use a recursive call for
   //compactness 

}

//______________________________________________________________________________
void _TheTree::RecoverNodeSupportStates2(node<long> *thisNode,
                                         _Parameter *resultVector,
                                         _Parameter *forwardVector,
                                         long catID) {
  _CalcNode *thisNodeC =
      ((_CalcNode *)((BaseRef *)variablePtrs.lData)[thisNode->in_object]);
  _Parameter *vecPointer = resultVector + thisNodeC->nodeIndex * cBase;

  if (thisNode->parent) {
    if (thisNode->parent->parent) {
      for (long cc = 0; cc < cBase; cc++, vecPointer++) {
        _Parameter tmp = 1.0;
        for (long nc = 0; nc < thisNode->parent->nodes.length; nc++) {
          _Parameter tmp2 = 0.0;
          _CalcNode *child = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[
              thisNode->parent->nodes.data[nc]->in_object]);
          bool invert = (child == thisNodeC);
          ;
          if (invert) {
            child = ((_CalcNode *)(
                (BaseRef *)variablePtrs.lData)[thisNode->parent->in_object]);
          }

          _Parameter *childSupport =
                         invert ? resultVector + cBase * child->nodeIndex
                                : forwardVector + child->nodeIndex * cBase,
                     *transMatrix =
                         child->GetCompExp(catID)->theData + cc * cBase;

          for (long cc2 = 0; cc2 < cBase; cc2++) {
            tmp2 += transMatrix[cc2] * childSupport[cc2];
          }
          tmp *= tmp2;
        }
        *vecPointer = tmp;
      }
    } else {
      for (long cc = 0; cc < cBase; cc++, vecPointer++) {
        _Parameter tmp = 1.0;
        for (long nc = 0; nc < thisNode->parent->nodes.length; nc++) {
          _Parameter tmp2 = 0.0;
          _CalcNode *child = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[
              thisNode->parent->nodes.data[nc]->in_object]);
          if (child != thisNodeC) {
            _Parameter *childSupport = forwardVector + child->nodeIndex * cBase,
                       *transMatrix =
                           child->GetCompExp(catID)->theData + cc * cBase;

            for (long cc2 = 0; cc2 < cBase; cc2++) {
              tmp2 += transMatrix[cc2] * childSupport[cc2];
            }

            tmp *= tmp2;
          }
        }
        *vecPointer = tmp;
      }
    }
  } else
    for (long cc = 0; cc < cBase; cc++) {
      vecPointer[cc] = 1.0;
    }

  for (long nc = 0; nc < thisNode->nodes.length; nc++) {
    RecoverNodeSupportStates2(thisNode->nodes.data[nc], resultVector,
                              forwardVector, catID);
  }
}

//______________________________________________________________________________
// 20090203/SLKP: clean up
_List *_TheTree::MapCBaseToCharacters(_DataSetFilter *dsf, bool normalanc) {

  _List *reply = new _List;
  checkPointer(reply);

  long unitLength = dsf->GetUnitLength();

  for (long pad = 0; pad < unitLength; pad++) {
    reply->AppendNewInstance(new _String(5, true));
  }

  _CalcNode *theChildNode = StepWiseTraversal(true);
  _String rootValue = dsf->ConvertCodeToLetters(
      dsf->CorrectCode(theChildNode->cBase), unitLength);

  for (long f = 0; f < unitLength; f++) {
    *((_String *)(*reply)(f)) << rootValue[f];
  }

  theChildNode = StepWiseTraversal(false);

  while (theChildNode) {
    if (!IsCurrentNodeATip()) {
      if (normalanc) {
        _CalcNode *travNode = (_CalcNode *)(
            (BaseRef *)variablePtrs.lData)[currentNode->parent->in_object];
        theChildNode->cBase = theChildNode->theProbs[travNode->cBase];
        theChildNode->categoryVariables
            .Delete(theChildNode->categoryVariables.lLength - 1);
      }

      _String letterValue = dsf->ConvertCodeToLetters(
          dsf->CorrectCode(theChildNode->cBase), unitLength);
      for (long i = 0; i < unitLength; i++) {
        *((_String *)(*reply)(i)) << letterValue[i];
      }
    }
    theChildNode = StepWiseTraversal(false);
  }

  for (long j = 0; j < unitLength; j++) {
    ((_String *)(*reply)(j))->Finalize();
  }

  return reply;
}

//______________________________________________________________________________
_Parameter _TheTree::ConditionalNodeLikelihood(node<long> *parentNode,
                                               node<long> *thisNode,
                                               _Parameter *scoreBelow,
                                               _Parameter *myScore,
                                               long myState, long offset) {

  if (parentNode) {
    _Parameter conditionalLikelihood = LIKELIHOOD_SCALER;

    for (long child = 0; child < thisNode->nodes.length; child++) {
      _CalcNode *childCNode =
          ((_CalcNode *)((BaseRef *)variablePtrs
                             .lData)[thisNode->nodes.data[child]->in_object]);
      conditionalLikelihood *=
          *(childCNode->compExp->theData + myState * cBase +
            childCNode->cBase) * childCNode->theValue;
    }

    myScore[myState] = conditionalLikelihood;

    return ConditionalBranchLikelihood(parentNode, thisNode, myScore,
                                       scoreBelow, -1, offset);

  } else {
    _Parameter tLik = theProbs[myState];

    for (long child = 0; child < thisNode->nodes.length; child++) {
      _CalcNode *childCNode =
          ((_CalcNode *)((BaseRef *)variablePtrs
                             .lData)[thisNode->nodes.data[child]->in_object]);
      tLik *= *(childCNode->compExp->theData + myState * cBase +
                childCNode->cBase) * childCNode->theValue;
    }

    return tLik;
  }
}

//______________________________________________________________________________
_Parameter _TheTree::ConditionalBranchLikelihood(node<long> *parentNode,
                                                 node<long> *thisNode,
                                                 _Parameter *scoreBelow,
                                                 _Parameter *myScore,
                                                 long myState, long offset) {

  for (long pstate = (myState >= 0 ? myState : 0);
       pstate < (myState >= 0 ? myState + 1 : cBase); pstate++) {
    _Parameter conditionalLikelihood = LIKELIHOOD_SCALER;
    for (long child = 0; child < parentNode->nodes.length; child++) {
      _Parameter childContribution = 0.0;

      node<long> *childNode = parentNode->nodes.data[child];
      _CalcNode *childCNode =
          ((_CalcNode *)((BaseRef *)variablePtrs.lData)[childNode->in_object]);

      _Parameter *conds,
          *transitions = childCNode->compExp->theData + pstate * cBase,
          *stopper;

      if (childNode == thisNode) { // use scoreBelow
        conds = scoreBelow;
      } else { // use existing conditionals
        if (offset >= 0) {
          conds = marginalLikelihoodCache +
                  cBase * (flatLeaves.lLength + flatNodes.lLength) * offset +
                  cBase * ((long) childCNode->theProbs[0]);
        } else {
          conds = childCNode->theProbs;
        }
      }

      long remaind = cBase % 4;

      if (remaind) {
        stopper = conds + (cBase - remaind);
        for (; conds != stopper;) {
          childContribution += *transitions * *conds;
          transitions++;
          conds++;
          childContribution += *transitions * *conds;
          transitions++;
          conds++;
          childContribution += *transitions * *conds;
          transitions++;
          conds++;
          childContribution += *transitions * *conds;
          transitions++;
          conds++;
        }

        switch (remaind) {
        case 1: {
          childContribution += *(transitions) * *(conds);
          break;
        }
        case 2: {
          childContribution += *transitions * *conds;
          transitions++;
          conds++;
          childContribution += *transitions * *conds;
          break;
        }
        case 3: {
          childContribution += *transitions * *conds;
          transitions++;
          conds++;
          childContribution += *transitions * *conds;
          transitions++;
          conds++;
          childContribution += *transitions * *conds;
          break;
        }
        }
      } else {
        stopper = conds + cBase;
        for (; conds != stopper;) {
          childContribution += *transitions * *conds;
          transitions++;
          conds++;
          childContribution += *transitions * *conds;
          transitions++;
          conds++;
          childContribution += *transitions * *conds;
          transitions++;
          conds++;
          childContribution += *transitions * *conds;
          transitions++;
          conds++;
        }
      }

      conditionalLikelihood *= childContribution;

      if (conditionalLikelihood == 0.0) {
        if (myState >= 0) {
          return 0.0;
        }
        break;
      }
    }
    myScore[pstate] = conditionalLikelihood;
  }

  if (parentNode->parent) {
    return ConditionalBranchLikelihood(parentNode->parent, parentNode, myScore,
                                       scoreBelow, -1, offset);
  }

  if (myState < 0) {
    _Parameter treeLikelihood = 0.0;

    for (long pstate = 0; pstate < cBase; pstate++) {
      treeLikelihood += theProbs[pstate] * myScore[pstate];
    }

    return treeLikelihood;
  }
  return theProbs[myState] * myScore[myState];
}

//______________________________________________________________________________
// presumes that compExp is populated and so are all the state vectors in
// internal nodes
void _TheTree::WeightedCharacterDifferences(_Parameter siteLikelihood,
                                            _Matrix *res1, _Matrix *res2,
                                            long offset) {
  if (cBase > 128) {
    WarnError("State spaces with more than 128 states are not supported in "
              "WeightedCharacterDifferences");
    return;
  }

  _Parameter vec1[128], vec2[128];

  //checkPointer (vec1);
  //checkPointer (vec2);

  //checkParameter (VerbosityLevelString,verbLevel, 0.0);
  //long          branchCount = 0;

  for (long nodeCount = 0; nodeCount < flatTree.lLength; nodeCount++) {
    node<long> *disNode = ((node<long> *)(flatNodes.lData[nodeCount]));

    for (long childCount = disNode->nodes.length - 1; childCount >= 0;
         childCount--) {
      node<long> *childNode = disNode->nodes.data[childCount];

      _CalcNode *cNode =
          ((_CalcNode *)((BaseRef *)variablePtrs.lData)[childNode->in_object]);
      _Parameter weightFactor = cNode->Value(), consistencyCheck = 0.0;

      for (long pState = 0; pState < cBase; pState++) {
        for (long cState = 0; cState < cBase; cState++) {
          for (long eraser = 0; eraser < cBase; eraser++) {
            vec1[eraser] = 0.;
            vec2[eraser] = 0.;
          }

          if (offset >= 0) {
            vec1[cState] =
                *(marginalLikelihoodCache +
                  cBase * (flatLeaves.lLength + flatNodes.lLength) * offset +
                  cBase * ((long) cNode->theProbs[0]) + cState);
          } else {
            vec1[cState] = cNode->theProbs[cState];
          }

          _Parameter cLik = (ConditionalBranchLikelihood(
              disNode, childNode, vec1, vec2, pState, offset) / siteLikelihood);
          res1->theData[pState * cBase + cState] += cLik;
          res2->theData[pState * cBase + cState] += cLik * weightFactor;
          consistencyCheck += cLik;
        }
      }

#if !defined __UNIX__ || defined __HEADLESS__ || defined __HYPHYQT__
      if ((cBase >= 20) && (offset < 1)) {
        yieldCPUTime();
        if (terminateExecution) {
          return;
        }
      }
#endif

      if (!CheckEqual(consistencyCheck, 1.)) {
        _String consistencyErra("Failed Internal Consistency Check In "
                                "WeightedCharacterDifferences at ");
        consistencyErra =
            consistencyErra & *LocateVar(disNode->in_object)->GetName() &
            " and " & *cNode->GetName() & ". Summed RLS to " & consistencyCheck;
        WarnError(consistencyErra);
      }
    }
  }

  //delete vec1;
  //delete vec2;
}

//______________________________________________________________________________
_AVLListX *_TheTree::ConstructNodeToIndexMap(bool doINodes) {
  _SimpleList *nodes = new _SimpleList,
              *whichL = doINodes ? &flatNodes : &flatLeaves;
  _AVLListX *result = new _AVLListX(nodes);

  for (unsigned long pistolero = 0; pistolero < whichL->lLength; pistolero++) {
    result->Insert((BaseRef) whichL->lData[pistolero], pistolero, false);
  }

  return result;

}

//______________________________________________________________________________
void _TheTree::MapPostOrderToInOderTraversal(_SimpleList &storeHere,
                                             bool doINodes) {
  _AVLListX *nodeMapper = ConstructNodeToIndexMap(doINodes);
  _CalcNode *traversalNode =
      doINodes ? StepWiseTraversal(true) : DepthWiseTraversal(true);

  long allNodeCount = 0;

  storeHere.Populate(doINodes ? flatTree.lLength : flatLeaves.lLength, 0, 0);

  while (traversalNode) {
    bool isTip = IsCurrentNodeATip();
    if ((isTip && !doINodes) || (!isTip && doINodes)) {
      storeHere.lData[nodeMapper->GetXtra(
          nodeMapper->Find((BaseRef)(&GetCurrentNode())))] = allNodeCount++;
    }

    traversalNode =
        doINodes ? StepWiseTraversal(false) : DepthWiseTraversal(false);
  }

  nodeMapper->DeleteAll(false);
  DeleteObject(nodeMapper);

}

//______________________________________________________________________________
void _TheTree::MarkDone(void) {
  _CalcNode *travNode = StepWiseTraversal(TRUE);

  while (travNode) {
    travNode->MarkDone();
    travNode = StepWiseTraversal();
  }
}

//______________________________________________________________________________
// notice that the limiting atomic probabilites are assumed to be stored in
// the
// model matrix of the tree itself
_Parameter _TheTree::PruneTree(long categID) {
  // do the depth first traversal to start at the leaves

  _CalcNode *curNode = DepthWiseTraversal(true);
  while (curNode) {
    bool reexpnt = curNode->NeedToExponentiate(categID);

    // flag all the nodes up the chain - to the root
    if (reexpnt && curNode->GetModelMatrix()) {
      curNode->RecomputeMatrix(categID, categoryCount);
    } else if (categID >= 0) {
      curNode->SetCompMatrix(categID);
    }

    // the probabilities below in the tree have already been computed
    long nNodes = currentNode->get_num_nodes();
    if (nNodes) {
      long j;
      for (j = 0; j < nNodes; j++) {
        _CalcNode *tC = (
            (_CalcNode *)(LocateVar(currentNode->get_node(j + 1)->get_data())));
        if (!tC->GetCompExp(categID)) {
          tC->RecomputeMatrix(categID, categoryCount);
        } else if (categID >= 0) {
          tC->SetCompMatrix(categID);
        }
      }
    }

    // the probabilities have now been set and the flag erased even if it had
    // been raised before
    curNode = DepthWiseTraversal();
  }

  return 0.0;
}

//______________________________________________________________________________
// notice that the limiting atomic probabilites are assumed to be stored in
// the
// model matrix of the tree itself
_Parameter _TheTree::PruneTreeChar4(long categID) {

  long nodeCount;
  _Parameter *fastIndex;
  node<long> *nodeChild;
  bool reexpnt;
  _CalcNode *theChildNode, *travNode;

  for (nodeCount = 0; nodeCount < flatLeaves.lLength; nodeCount++) {
    theChildNode = (_CalcNode *)(((BaseRef *)flatCLeaves.lData)[nodeCount]);
    reexpnt = theChildNode->NeedToExponentiate(categID);
    if (reexpnt && theChildNode->GetModelMatrix()) {
      theChildNode->RecomputeMatrix(categID, categoryCount);
    } else if (categID >= 0) {
      theChildNode->SetCompMatrix(categID);
    }
  }

  for (nodeCount = 0; nodeCount < flatTree.lLength; nodeCount++) {
    theChildNode = (_CalcNode *)(((BaseRef *)flatTree.lData)[nodeCount]);
    theChildNode->theProbs[0] = LIKELIHOOD_SCALER;
    theChildNode->theProbs[1] = LIKELIHOOD_SCALER;
    theChildNode->theProbs[2] = LIKELIHOOD_SCALER;
    theChildNode->theProbs[3] = LIKELIHOOD_SCALER;

    reexpnt = theChildNode->NeedToExponentiate(categID);
    if (reexpnt && theChildNode->GetModelMatrix()) {
      theChildNode->RecomputeMatrix(categID, categoryCount);
      theChildNode->cBase = -1;
    } else {
      if (categID >= 0) {
        theChildNode->SetCompMatrix(categID);
      }
    }
  }

  for (nodeCount = 0; nodeCount < flatLeaves.lLength; nodeCount++) {
    travNode = (_CalcNode *)(((BaseRef *)flatCLeaves.lData)[nodeCount]);
    theChildNode = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[
        ((node<long> *)(flatLeaves.lData[nodeCount]))->parent->in_object]);

    fastIndex = travNode->compExp->theData;
    if (travNode->lastState >= 0) {
      fastIndex += travNode->lastState;
      theChildNode->theProbs[0] *= *fastIndex;
      theChildNode->theProbs[1] *= fastIndex[4];
      theChildNode->theProbs[2] *= fastIndex[8];
      theChildNode->theProbs[3] *= fastIndex[12];
    } else {
      _Parameter a = *(travNode->theProbs), b = travNode->theProbs[1],
                 c = travNode->theProbs[2], d = travNode->theProbs[3];

      _Parameter tmp = a * *(fastIndex++);
      tmp += b * *(fastIndex++);
      tmp += c * *(fastIndex++);
      theChildNode->theProbs[0] *= tmp + d * *(fastIndex++);

      tmp = a * *(fastIndex++);
      tmp += b * *(fastIndex++);
      tmp += c * *(fastIndex++);
      theChildNode->theProbs[1] *= tmp + d * *(fastIndex++);

      tmp = a * *(fastIndex++);
      tmp += b * *(fastIndex++);
      tmp += c * *(fastIndex++);
      theChildNode->theProbs[2] *= tmp + d * *(fastIndex++);

      tmp = a * *(fastIndex++);
      tmp += b * *(fastIndex++);
      tmp += c * *(fastIndex++);
      theChildNode->theProbs[3] *= tmp + d * *(fastIndex++);
    }
  }

  for (nodeCount = 0; nodeCount < flatTree.lLength - 1; nodeCount++) {
    travNode = (_CalcNode *)(((BaseRef *)flatTree.lData)[nodeCount]);
    nodeChild = ((node<long> *)(flatNodes.lData[nodeCount]))->parent;
    theChildNode =
        ((_CalcNode *)((BaseRef *)variablePtrs.lData)[nodeChild->in_object]);

    fastIndex = travNode->compExp->theData;

    _Parameter a = *(travNode->theProbs), b = travNode->theProbs[1],
               c = travNode->theProbs[2], d = travNode->theProbs[3];

    _Parameter tmp = a * *(fastIndex++);
    tmp += b * *(fastIndex++);
    tmp += c * *(fastIndex++);
    theChildNode->theProbs[0] *= tmp + d * *(fastIndex++);

    tmp = a * *(fastIndex++);
    tmp += b * *(fastIndex++);
    tmp += c * *(fastIndex++);
    theChildNode->theProbs[1] *= tmp + d * *(fastIndex++);

    tmp = a * *(fastIndex++);
    tmp += b * *(fastIndex++);
    tmp += c * *(fastIndex++);
    theChildNode->theProbs[2] *= tmp + d * *(fastIndex++);

    tmp = a * *(fastIndex++);
    tmp += b * *(fastIndex++);
    tmp += c * *(fastIndex++);
    theChildNode->theProbs[3] *= tmp + d * *(fastIndex++);
  }

  theChildNode =
      (_CalcNode *)(((BaseRef *)variablePtrs.lData)[theRoot->in_object]);

  _Parameter result = theProbs[0] * theChildNode->theProbs[0] +
                      theProbs[1] * theChildNode->theProbs[1] +
                      theProbs[2] * theChildNode->theProbs[2] +
                      theProbs[3] * theChildNode->theProbs[3];

  if (result <= 0.0) {
    return ALMOST_ZERO;
  }
  return result;
}

//______________________________________________________________________________
// notice that the limiting atomic probabilites are assumed to be stored in
// the
// model matrix of the tree itself
_Parameter _TheTree::PruneTreeChar4Cache(long categID) {

  long nodeCount;
  _Parameter *fastIndex;
  node<long> *nodeChild;
  bool reexpnt;
  _CalcNode *theChildNode, *travNode;

  for (nodeCount = 0; nodeCount < flatLeaves.lLength; nodeCount++) {
    theChildNode = (_CalcNode *)(((BaseRef *)flatCLeaves.lData)[nodeCount]);
    reexpnt = theChildNode->NeedToExponentiate(categID);
    if (reexpnt && theChildNode->GetModelMatrix()) {
      theChildNode->RecomputeMatrix(categID, categoryCount);
      theChildNode = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[
          ((node<long> *)(flatLeaves.lData[nodeCount]))->parent->in_object]);
      theChildNode->cBase = -1;
    } else if (categID >= 0) {
      theChildNode->SetCompMatrix(categID);
    }
  }

  for (nodeCount = 0; nodeCount < flatTree.lLength; nodeCount++) {
    theChildNode = (_CalcNode *)(((BaseRef *)flatTree.lData)[nodeCount]);
    theChildNode->theProbs[0] = LIKELIHOOD_SCALER;
    theChildNode->theProbs[1] = LIKELIHOOD_SCALER;
    theChildNode->theProbs[2] = LIKELIHOOD_SCALER;
    theChildNode->theProbs[3] = LIKELIHOOD_SCALER;
    reexpnt = theChildNode->NeedToExponentiate(categID);
    if (reexpnt && theChildNode->GetModelMatrix()) {
      theChildNode->RecomputeMatrix(categID, categoryCount);
      theChildNode->cBase = -1;
    } else {
      if (categID >= 0) {
        theChildNode->SetCompMatrix(categID);
      }
    }
    if (theChildNode->cBase == -1) {
      nodeChild = ((node<long> *)(flatNodes.lData[nodeCount]))->parent;
      if (nodeChild) {
        theChildNode = (
            (_CalcNode *)((BaseRef *)variablePtrs.lData)[nodeChild->in_object]);
        theChildNode->cBase = -1;
      }
    }
  }

  for (nodeCount = 0; nodeCount < flatLeaves.lLength; nodeCount++) {
    travNode = (_CalcNode *)(((BaseRef *)flatCLeaves.lData)[nodeCount]);
    theChildNode = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[
        ((node<long> *)(flatLeaves.lData[nodeCount]))->parent->in_object]);

    fastIndex = travNode->GetCompExp()->fastIndex();
    if (travNode->lastState >= 0) {
      fastIndex += travNode->lastState;
      theChildNode->theProbs[0] *= *fastIndex;
      theChildNode->theProbs[1] *= fastIndex[4];
      theChildNode->theProbs[2] *= fastIndex[8];
      theChildNode->theProbs[3] *= fastIndex[12];
    } else {
      _Parameter a = *(travNode->theProbs), b = travNode->theProbs[1],
                 c = travNode->theProbs[2], d = travNode->theProbs[3];

      _Parameter tmp = a * *(fastIndex++);
      tmp += b * *(fastIndex++);
      tmp += c * *(fastIndex++);
      theChildNode->theProbs[0] *= tmp + d * *(fastIndex++);

      tmp = a * *(fastIndex++);
      tmp += b * *(fastIndex++);
      tmp += c * *(fastIndex++);
      theChildNode->theProbs[1] *= tmp + d * *(fastIndex++);

      tmp = a * *(fastIndex++);
      tmp += b * *(fastIndex++);
      tmp += c * *(fastIndex++);
      theChildNode->theProbs[2] *= tmp + d * *(fastIndex++);

      tmp = a * *(fastIndex++);
      tmp += b * *(fastIndex++);
      tmp += c * *(fastIndex++);
      theChildNode->theProbs[3] *= tmp + d * *(fastIndex++);
    }
  }

  for (nodeCount = 0; nodeCount < flatTree.lLength - 1; nodeCount++) {
    travNode = (_CalcNode *)(((BaseRef *)flatTree.lData)[nodeCount]);
    nodeChild = ((node<long> *)(flatNodes.lData[nodeCount]))->parent;
    theChildNode =
        ((_CalcNode *)((BaseRef *)variablePtrs.lData)[nodeChild->in_object]);

    fastIndex = travNode->GetCompExp()->fastIndex();

    _Parameter a = *(travNode->theProbs), b = travNode->theProbs[1],
               c = travNode->theProbs[2], d = travNode->theProbs[3];

    _Parameter tmp = a * *(fastIndex++);
    tmp += b * *(fastIndex++);
    tmp += c * *(fastIndex++);
    theChildNode->theProbs[0] *= tmp + d * *(fastIndex++);

    tmp = a * *(fastIndex++);
    tmp += b * *(fastIndex++);
    tmp += c * *(fastIndex++);
    theChildNode->theProbs[1] *= tmp + d * *(fastIndex++);

    tmp = a * *(fastIndex++);
    tmp += b * *(fastIndex++);
    tmp += c * *(fastIndex++);
    theChildNode->theProbs[2] *= tmp + d * *(fastIndex++);

    tmp = a * *(fastIndex++);
    tmp += b * *(fastIndex++);
    tmp += c * *(fastIndex++);
    theChildNode->theProbs[3] *= tmp + d * *(fastIndex++);
  }

  long c = 0, off = 1;
  for (nodeCount = 0; nodeCount < topLevelNodes.lLength - 1;
       nodeCount++, off <<= 1) {
    travNode = (_CalcNode *)(
        ((BaseRef *)flatTree.lData)[topLevelNodes.lData[nodeCount]]);
    if (travNode->cBase <= -1) {
      c |= off;
    }
  }
  topLevelNodes.lData[topLevelNodes.lLength - 1] = c;
  for (nodeCount = 0; nodeCount < flatTree.lLength; nodeCount++) {
    theChildNode = (_CalcNode *)(((BaseRef *)flatTree.lData)[nodeCount]);
    theChildNode->cBase = cBase;
  }

  _Parameter result = 0.0;

  theChildNode =
      (_CalcNode *)(((BaseRef *)variablePtrs.lData)[theRoot->in_object]);
  //((_CalcNode*)(LocateVar(theRoot->get_data())));

  for (nodeCount = 0; nodeCount < cBase; nodeCount++) {
    result += theProbs[nodeCount] * theChildNode->theProbs[nodeCount];
  }

  if (result <= 0.0) {
    return ALMOST_ZERO;
  }

  return result;
}

//______________________________________________________________________________
// notice that the limiting atomic probabilites are assumed to be stored in
// the
// model matrix of the tree itself
_Parameter _TheTree::PruneTreeChar(long categID) {

  long nodeCount, j;
  _Parameter *fastIndex;
  node<long> *nodeChild;
  bool reexpnt;
  _CalcNode *theChildNode, *travNode;

  for (nodeCount = 0; nodeCount < flatTree.lLength; nodeCount++) {
    theChildNode = (_CalcNode *)(((BaseRef *)flatTree.lData)[nodeCount]);
    for (j = 0; j < cBase; j++) {
      theChildNode->theProbs[j] = LIKELIHOOD_SCALER;
    }
    reexpnt = theChildNode->NeedToExponentiate(categID);
    if (reexpnt && theChildNode->GetModelMatrix()) {
      theChildNode->RecomputeMatrix(categID, categoryCount);
    } else if (categID >= 0) {
      theChildNode->SetCompMatrix(categID);
    }
  }

  for (nodeCount = 0; nodeCount < flatLeaves.lLength; nodeCount++) {
    theChildNode = (_CalcNode *)(((BaseRef *)flatCLeaves.lData)[nodeCount]);
    reexpnt = theChildNode->NeedToExponentiate(categID);
    if (reexpnt && theChildNode->GetModelMatrix()) {
      theChildNode->RecomputeMatrix(categID, categoryCount);
    } else if (categID >= 0) {
      theChildNode->SetCompMatrix(categID);
    }
  }

  for (nodeCount = 0; nodeCount < flatLeaves.lLength; nodeCount++) {
    travNode = (_CalcNode *)(((BaseRef *)flatCLeaves.lData)[nodeCount]);
    theChildNode = ((_CalcNode *)((BaseRef *)variablePtrs.lData)[
        ((node<long> *)(flatLeaves.lData[nodeCount]))->parent->in_object]);

    fastIndex = travNode->GetCompExp()->fastIndex();
    if (travNode->lastState >= 0) {
      fastIndex += travNode->lastState;
      _Parameter tmp = travNode->theProbs[travNode->lastState];
      for (long i = 0; i < cBase; i++, fastIndex += cBase) {
        theChildNode->theProbs[i] *= tmp * (*fastIndex);
      }
    } else {
      for (j = 0; j < cBase; j++) {
        _Parameter tmp = *(travNode->theProbs) * *fastIndex;
        fastIndex++;
        for (long k = 1; k < cBase; k++, fastIndex++) {
          tmp += travNode->theProbs[k] * (*fastIndex);
        }
        theChildNode->theProbs[j] *= tmp;
      }
    }
  }

  for (nodeCount = 0; nodeCount < flatTree.lLength; nodeCount++) {
    travNode = (_CalcNode *)(((BaseRef *)flatTree.lData)[nodeCount]);
    nodeChild = ((node<long> *)(flatNodes.lData[nodeCount]))->parent;
    if (nodeChild) { // not a root yet
      theChildNode =
          ((_CalcNode *)((BaseRef *)variablePtrs.lData)[nodeChild->in_object]);

      fastIndex = travNode->GetCompExp()->fastIndex();

      for (j = 0; j < cBase; j++) {
        _Parameter tmp = *(travNode->theProbs) * *fastIndex;
        fastIndex++;
        for (long k = 1; k < cBase; k++, fastIndex++) {
          tmp += travNode->theProbs[k] * (*fastIndex);
        }
        theChildNode->theProbs[j] *= tmp;
      }

    }
  }

  _Parameter result = 0;

  theChildNode =
      (_CalcNode *)(((BaseRef *)variablePtrs.lData)[theRoot->in_object]);
  //((_CalcNode*)(LocateVar(theRoot->get_data())));

  for (nodeCount = 0; nodeCount < cBase; nodeCount++) {
    result += theProbs[nodeCount] * theChildNode->theProbs[nodeCount];
  }

  if (result <= 0.0) {
    return ALMOST_ZERO;
  }
  return result;
}

//______________________________________________________________________________
long _TheTree::ComputeReleafingCostChar(_DataSetFilter *dsf, long firstIndex,
                                        long secondIndex) {

  char *pastState = dsf->GetColumn(firstIndex),
       *thisState = dsf->GetColumn(secondIndex);

  _SimpleList markedNodes(flatTree.lLength, 0, 0);

  for (long nodeID = 0; nodeID < flatLeaves.lLength; nodeID++) {
    long f = dsf->theNodeMap.lData[nodeID];
    if (thisState[f] != pastState[f]) {
      markedNodes.lData[flatParents.lData[nodeID]] = 1;
    }
  }

  long theCost = 0;

  for (long i = 0; i < flatTree.lLength; i++) {
    if (markedNodes.lData[i]) {
      long myParent = flatParents.lData[i + flatLeaves.lLength];
      if (myParent >= 0) {
        markedNodes.lData[myParent] = 1;
      }
      theCost += ((node<long> *)(flatNodes.lData[i]))->nodes.length;
    }
  }

  return theCost;

}

//______________________________________________________________________________
void _TheTree::ClearConstraints(void) {
  _CalcNode *travNode = StepWiseTraversal(true);
  while (travNode) {
    travNode->ClearConstraints();
    travNode = StepWiseTraversal();
  }
}

//______________________________________________________________________________
long _TheTree::ComputeReleafingCost(_DataSetFilter *dsf, long firstIndex,
                                    long secondIndex,
                                    _SimpleList *traversalTags,
                                    long orderIndex) {

  //static _SimpleList flatLeaves, nodeCount;
  //static _List     flatTree;
  long filterL = dsf->NumberDistinctSites();

  _SimpleList markedNodes(flatTree.lLength, 0, 0);

  for (long leafID = 0; leafID < flatLeaves.lLength; leafID++)
    if (!dsf->CompareTwoSites(firstIndex, secondIndex, leafID)) {
      markedNodes.lData[flatParents.lData[leafID]] = 1;
    }

  // now compute the cost

  long theCost = 0;

  for (long i = 0; i < flatTree.lLength; i++) {
    if (markedNodes.lData[i]) {
      long myParent = flatParents.lData[flatLeaves.lLength + i];
      if (myParent >= 0) {
        markedNodes.lData[myParent] = 1;
      }
      theCost += ((node<long> *)(flatNodes.lData[i]))->nodes.length;
    } else if (traversalTags && orderIndex) {
      long theIndex = filterL * i + orderIndex;
      traversalTags->lData[theIndex / _HY_BITMASK_WIDTH_] |=
          bitMaskArray.masks[theIndex % _HY_BITMASK_WIDTH_];
      /*printf ("%d %d %d %d %d %d\n", firstIndex, secondIndex, orderIndex, i,
      theIndex/_HY_BITMASK_WIDTH_, theIndex%_HY_BITMASK_WIDTH_);
      char b[4]; b[3] = 0;
      for (long k = 0; k < flatLeaves.lLength; k++)
      {
          dsf->GrabSite (firstIndex, dsf->theNodeMap.lData[k], b);
          printf ("%s ", b);
          dsf->GrabSite (secondIndex, dsf->theNodeMap.lData[k], b);
          printf ("%s\n", b);
      }*/
    }
  }

  return theCost;

}

//______________________________________________________________________________
void _TheTree::MarkMatches(_DataSetFilter *dsf, long firstIndex,
                           long secondIndex) {

  long n = 0, f;

  _CalcNode *travNode;

  for (n = 0; n < flatLeaves.lLength; n++) {
    travNode = (_CalcNode *)(((BaseRef *)flatCLeaves.lData)[n]);
    if (!dsf->CompareTwoSites(firstIndex, secondIndex, n)) {
      node<long> *theTreeNode = ((node<long> *)(flatLeaves.lData[n]))->parent;
      _CalcNode *cN = (
          (_CalcNode *)((BaseRef *)variablePtrs.lData)[theTreeNode->in_object]);
      cN->cBase = -1;
    }
  }
  n = 0;
  for (f = 0; f < flatTree.lLength; f++) {
    travNode = (_CalcNode *)(((BaseRef *)flatTree.lData)[f]);
    if (travNode->cBase == -1) {
      node<long> *theTreeNode = ((node<long> *)(flatNodes.lData[f]))->parent;
      if (theTreeNode) {
        _CalcNode *cN = ((_CalcNode *)(
            (BaseRef *)variablePtrs.lData)[theTreeNode->in_object]);
        cN->cBase = -1;
      }
    }
  }
  for (f = 0; f < flatTree.lLength; f++) {
    travNode = (_CalcNode *)(((BaseRef *)flatTree.lData)[f]);
    if (travNode->cBase != -1) {
      travNode->lastState = -2;
    } else {
      travNode->cBase = cBase;
    }
  }
}

//______________________________________________________________________________

long _TheTree::GetLowerBoundOnCost(_DataSetFilter *dsf) {
  long n = 0, theCost = 0;

  _CalcNode *travNode;

  for (long siteIndex = 0; siteIndex < dsf->theFrequencies.lLength;
       siteIndex++) {
    for (n = 0; n < flatTree.lLength; n++) {
      travNode = (_CalcNode *)(((BaseRef *)flatTree.lData)[n]);
      travNode->lastState = -1;
    }
    for (long matchIndex = 0; matchIndex < dsf->theFrequencies.lLength;
         matchIndex++)
      if (matchIndex != siteIndex) {
        MarkMatches(dsf, siteIndex, matchIndex);
      }
    for (n = 0; n < flatTree.lLength; n++) {
      travNode = (_CalcNode *)(((BaseRef *)flatTree.lData)[n]);
      if (travNode->lastState != -2) {
        theCost += ((node<long> *)(flatNodes.lData[n]))->nodes.length;
      }
      travNode->lastState = -1;
    }
  }
  return theCost;
}

//______________________________________________________________________________
long _TheTree::GetLowerBoundOnCostWithOrder(_DataSetFilter *dsf,
                                            _SimpleList *sl) {
  long n = 0, theCost = 0;

  _CalcNode *travNode;

  for (long siteIndex = 0; siteIndex < dsf->theFrequencies.lLength;
       siteIndex++) {
    for (n = 0; n < flatTree.lLength; n++) {
      travNode = (_CalcNode *)(((BaseRef *)flatTree.lData)[n]);
      travNode->lastState = -1;
    }
    for (long matchIndex = 0; matchIndex < siteIndex; matchIndex++)
      if (matchIndex != siteIndex) {
        MarkMatches(dsf, sl->lData[siteIndex], sl->lData[matchIndex]);
      }

    for (n = 0; n < flatTree.lLength; n++) {
      travNode = (_CalcNode *)(((BaseRef *)flatTree.lData)[n]);
      if (travNode->lastState != -2) {
        theCost += ((node<long> *)(flatNodes.lData[n]))->nodes.length;
      }
      travNode->lastState = -1;
    }
  }
  return theCost;
}

//______________________________________________________________________________
void _TheTree::DumpingOrder(_DataSetFilter *dsf, _SimpleList &receptacle) {

  _SimpleList flatLeaves, nodeCount;
  _List flatTree;
  long theCost, firstIndex, secondIndex;

  _CalcNode *travNode;

  travNode = StepWiseTraversal(TRUE);

  while (travNode) {
    travNode->GetProbs()[1] = 1;
    flatTree << travNode;
    nodeCount << currentNode->get_num_nodes();
    travNode = StepWiseTraversal();
    receptacle << receptacle.lLength;
  }

  flatLeaves.Clear();
  travNode = LeafWiseTraversal(TRUE);

  while (travNode) {
    flatLeaves << long(currentNode);
    travNode = LeafWiseTraversal();
  }

  for (firstIndex = 0, secondIndex = 1;
       secondIndex < dsf->NumberDistinctSites(); firstIndex = secondIndex++) {
    for (theCost = 0; theCost < flatLeaves.lLength; theCost++) {
      if ((*dsf)(firstIndex, theCost) != (*dsf)(secondIndex, theCost)) {
        // mark all the nodes up the ladder as "tainted"
        node<long> *theTreeNode = (node<long> *)(flatLeaves(theCost));
        while (theTreeNode) {
          ((_CalcNode *)(LocateVar(theTreeNode->get_data())))->SetSummedFlag();
          theTreeNode = theTreeNode->get_parent();
        }
      }
    }

    // now compute the cost

    theCost = 0;

    for (long i = 0; i < flatTree.lLength; i++) {
      travNode = (_CalcNode *)flatTree(i);
      if (travNode->IsSummedFlagged()) {
        travNode->RemoveSummedFlag();
        travNode->GetProbs()[1]++;
      }
    }
  }

  _SimpleList ref;
  for (theCost = 0; theCost < flatTree.lLength; theCost++) {
    ref << ((_CalcNode *)flatTree(theCost))->GetProbs()[1];
  }
  SortLists(&ref, &receptacle);

}

void _TheTree::MolecularClock(_String &baseNode, _List &varsToConstrain) {
  node<long> *topNode = nil;
  _CalcNode *curNode = StepWiseTraversal(true);
  if (!baseNode.Length()) { // called Molecular Clock on the entire tree
    topNode = &GetRoot();
    _String *childNameP;
    if (rooted ==
        ROOTED_LEFT) { // run separate constraint on the right child of the root
      childNameP = ((_CalcNode *)(((BaseRef *)variablePtrs.lData)[
          theRoot->go_down(theRoot->get_num_nodes())->in_object]))->GetName();
      _String childName = childNameP->Cut(childNameP->Find('.') + 1, -1);
      MolecularClock(childName, varsToConstrain);
    } else if (rooted == ROOTED_RIGHT) {
      childNameP = ((_CalcNode *)(
          ((BaseRef *)variablePtrs.lData)[theRoot->go_down(1)->in_object]))
          ->GetName();
      _String childName = childNameP->Cut(childNameP->Find('.') + 1, -1);
      MolecularClock(childName, varsToConstrain);
    }
  } else {
    baseNode = _String(".") & baseNode;
    while (curNode) {
      if (curNode->GetName()->endswith(baseNode)) {
        topNode = currentNode;
        break;
      }
      curNode = StepWiseTraversal();
    }
  }

  if (!topNode) {
    WarnError(_String("Molecular clock constraint has failed, since node '") &
              baseNode & "' is not a part of tree '" & *GetName() & "'");
  } else
    for (unsigned long k = 1; k < varsToConstrain.lLength; k++) {
      long varIndex = LocateVarByName(*(_String *)varsToConstrain(k));
      if (varIndex < 0) {
        WarnError(
            _String("Molecular clock constraint has failed, since variable' ") &
            *(_String *)varsToConstrain(k) & "' is undefined.");
        return;
      }
      curNode->RecurseMC(variableNames.GetXtra(varIndex), topNode, true,
                         rooted);
    }
}

void _TheTree::ScanSubtreeVars(_List &rec, char flags, _CalcNode *startAt)
    // flags = 1 - do ind
    // flags = 2 - do dep
    {
  _SimpleList scanVars;
  _VariableContainer *thisV;
  _String chop;
  long k, f;

  if (startAt) {
    thisV = startAt;
  } else {
    thisV = DepthWiseTraversal(true);
  }

  {
    _AVLList scanVarsA(&scanVars);
    if (flags & 0x01) {
      thisV->ScanForVariables(scanVarsA, scanVarsA);
    }
    if (flags & 0x02) {
      thisV->ScanForDVariables(scanVarsA, scanVarsA);
    }

    scanVarsA.ReorderList();
  }

  for (k = 0; k < scanVars.lLength; k++) {
    thisV = (_VariableContainer *)LocateVar(scanVars.lData[k]);
    f = thisV->GetName()->FindBackwards('.', 0, -1);
    if (f >= 0) {
      chop = thisV->GetName()->Cut(f + 1, -1);
      rec &&&chop;
    }
  }

  if (startAt) {
    thisV = StepWiseTraversalLevel(k, true);
    while (thisV && (thisV != startAt)) {
      thisV = StepWiseTraversalLevel(k);
    }
    if (thisV) {
      f = k;
      while ((k == f) && thisV) {
        thisV = StepWiseTraversalLevel(k);
      }
      if (thisV) {
        while ((k > f) && rec.lLength) {
          thisV->MatchParametersToList(rec, true, (flags & 0x02) != 0);
          thisV = StepWiseTraversalLevel(k);
        }
        return;
      }
    }
    rec.Clear();
  } else {
    thisV = DepthWiseTraversal();
    while (thisV && rec.lLength && (!IsCurrentNodeTheRoot())) {
      thisV->MatchParametersToList(rec, true, (flags & 0x02) != 0);
      thisV = DepthWiseTraversal();
    }
  }
}

bool _TheTree::MatchLeavesToDF(_SimpleList &tipMatches, _DataSetFilter *df,
                               bool doNumeric) {
  tipMatches.Clear();
  _CalcNode *travNode = StepWiseTraversal(true);
  _List tips;
  long j, k;

  while (travNode) {
    if (IsCurrentNodeATip()) {
      _String tipName(travNode->GetName()->Cut(
          travNode->GetName()->FindBackwards('.', 0, -1) + 1, -1));
      tips &&&tipName;
    }
    travNode = StepWiseTraversal(false);
  }

  /*for (j=0;j<tips.lLength;j++)
  {
      k = df->FindSpeciesName((_String*)tips(j));
      if (k==-1) break;
      tipMatches<<k;
  }*/

  j = df->FindSpeciesName(tips, tipMatches);

  if (doNumeric) {
    if (j != tips.lLength) {
      long sj = j;
      for (j = 0; j < tips.lLength; j++) {
        _String *thisName = (_String *)tips(j);
        k = atoi(thisName->sData);
        _String tryAgain(k);
        if ((tryAgain.Equal(thisName)) && (k <= tips.lLength)) {
          tipMatches << k;
        } else {
          break;
        }
      }
      if (j == tips.lLength) {
        if (tipMatches.Find(0) == -1) // map to indexing from 0
          for (j = 0; j < tips.lLength; j++) {
            tipMatches.lData[j]--;
          }
      } else {
        j = sj;
      }
    }
  }

  return (j == tips.lLength);
}

void _TheTree::AddNodeNamesToDS(_DataSet *ds, bool doTips, bool doInternals,
                                char dOrS) {
  if (dOrS == 2 && doTips && doInternals) {
    AddNodeNamesToDS(ds, false, true, 0);
    AddNodeNamesToDS(ds, true, false, 0);
    return;
  }

  _CalcNode *iNodeTraverser =
      dOrS ? DepthWiseTraversal(true) : StepWiseTraversal(true);

  long j = GetName()->sLength + 1;

  while (iNodeTraverser) {
    if (IsCurrentNodeATip()) {
      if (doTips) {
        ds->GetNames()
            .AppendNewInstance(new _String(*iNodeTraverser->GetName(), j, -1));
      }
    } else if (doInternals) {
      ds->GetNames()
          .AppendNewInstance(new _String(*iNodeTraverser->GetName(), j, -1));
    }

    iNodeTraverser =
        dOrS ? DepthWiseTraversal(false) : StepWiseTraversal(false);
  }
}

//______________________________________________________________________________
void _TheTree::ExponentiateMatrices(_List &expNodes, long tc, long catID) {
  _List matrixQueue, nodesToDo;

  _SimpleList isExplicitForm;
  bool hasExpForm = false;

  for (unsigned long nodeID = 0; nodeID < expNodes.lLength; nodeID++) {
    long didIncrease = matrixQueue.lLength;
    _CalcNode *thisNode = (_CalcNode *)expNodes(nodeID);
    if (thisNode->RecomputeMatrix(catID, categoryCount, nil, &matrixQueue,
                                  &isExplicitForm)) {
      hasExpForm = true;
    }
    //printf ("NodeID %d. Old length %ld, new length %ld\n", nodeID,
    //didIncrease,matrixQueue.lLength);
    if ((didIncrease = (matrixQueue.lLength - didIncrease))) {
      for (long copies = 0; copies < didIncrease; copies++) {
        nodesToDo << thisNode;
      }
    }
  }

  //printf ("%ld %d\n", nodesToDo.lLength, hasExpForm);

  unsigned long matrixID;

  _List *computedExponentials =
      hasExpForm ? new _List(matrixQueue.lLength) : nil;

#ifdef _OPENMP
  long nt = cBase < 20 ? 1 : (MIN(tc, matrixQueue.lLength / 3 + 1));
  matrixExpCount += matrixQueue.lLength;
#endif

#pragma omp parallel for default(shared) private(matrixID)                     \
    schedule(static) if (nt > 1) num_threads(nt)
  for (matrixID = 0; matrixID < matrixQueue.lLength; matrixID++) {
    if (isExplicitForm.lData[matrixID] == 0) { // normal matrix to exponentiate
      ((_CalcNode *)nodesToDo(matrixID))->SetCompExp(
          ((_Matrix *)matrixQueue(matrixID))->Exponentiate(), catID);
    } else {
      (*computedExponentials)[matrixID] =
          ((_Matrix *)matrixQueue(matrixID))->Exponentiate();
    }
  }

  if (computedExponentials) {
    _CalcNode *current_node = nil;
    _List buffered_exponentials;

    for (unsigned long mx_index = 0; mx_index < nodesToDo.lLength; mx_index++) {
      if (isExplicitForm.lData[mx_index]) {
        _CalcNode *next_node = (_CalcNode *)nodesToDo(mx_index);
        //printf ("%x %x\n", current_node, next_node);
        if (next_node != current_node) {
          if (current_node) {
            current_node->RecomputeMatrix(catID, categoryCount, nil, nil, nil,
                                          &buffered_exponentials);
          }
          current_node = next_node;
          buffered_exponentials.Clear(true);
          buffered_exponentials.AppendNewInstance(
              (*computedExponentials)(mx_index));
        } else {
          buffered_exponentials.AppendNewInstance(
              (*computedExponentials)(mx_index));
        }
      } else {
        if (current_node) {
          current_node->RecomputeMatrix(catID, categoryCount, nil, nil, nil,
                                        &buffered_exponentials);
        }
        current_node = nil;
      }
    }
    if (current_node) {
      current_node->RecomputeMatrix(catID, categoryCount, nil, nil, nil,
                                    &buffered_exponentials);
    }
    DeleteObject(computedExponentials);
  }
}

//______________________________________________________________________________
long _TheTree::DetermineNodesForUpdate(_SimpleList &updateNodes,
                                       _List *expNodes, long catID, long addOne,
                                       bool canClear) {
  nodesToUpdate.Populate(flatLeaves.lLength + flatTree.lLength - 1, 0, 0);
  _CalcNode *currentTreeNode;
  long lastNodeID = -1;

  // look for nodes with model changes and mark the path up to the root as
  // needing an update

  if (addOne >= 0) {
    nodesToUpdate.lData[addOne] = 1;
  }

  if (forceRecalculationOnTheseBranches.lLength) {
    for (unsigned long markedNode = 0;
         markedNode < forceRecalculationOnTheseBranches.lLength; markedNode++) {
      nodesToUpdate.lData[forceRecalculationOnTheseBranches.lData[markedNode]] =
          1;
    }

    if (canClear) {
      forceRecalculationOnTheseBranches.Clear();
    }
  }

  for (unsigned long nodeID = 0; nodeID < nodesToUpdate.lLength; nodeID++) {
    bool isLeaf = nodeID < flatLeaves.lLength;

    currentTreeNode =
        isLeaf ? (((_CalcNode **)flatCLeaves.lData)[nodeID])
               : (((_CalcNode **)flatTree.lData)[nodeID - flatLeaves.lLength]);

    if (currentTreeNode->NeedToExponentiate(catID)) {
      if (expNodes) {
        (*expNodes) << currentTreeNode;
        //printf ("EXP>%s\n", currentTreeNode->GetName()->sData);
        lastNodeID = nodeID;
      } else {
        currentTreeNode->RecomputeMatrix(catID, categoryCount, nil);
      }

      nodesToUpdate.lData[nodeID] = 1;
    }

    if (nodesToUpdate.lData[nodeID]) {
      nodesToUpdate.lData[flatParents.lData[nodeID] + flatLeaves.lLength] = 1;
    }
  }

  // one more pass to pick up all descendants of changed internal nodes

  for (unsigned long nodeID = 0; nodeID < nodesToUpdate.lLength; nodeID++)
    if (nodesToUpdate.lData[flatLeaves.lLength + flatParents.lData[nodeID]] &&
        nodesToUpdate.lData[nodeID] == 0) {
      nodesToUpdate.lData[nodeID] = 1;
    }

  // write out all changed nodes

  for (unsigned long nodeID = 0; nodeID < nodesToUpdate.lLength; nodeID++)
    if (nodesToUpdate.lData[nodeID]) {
      updateNodes << nodeID;
    }

  if (expNodes && expNodes->lLength == 1) {
    return lastNodeID;
  }

  return -1;
}

//______________________________________________________________________________
// this utility function will simply fill in all the conditional probability
// vectors for internal nodes,
// including those that were skipped due to column sorting optimization
// this is useful to avoid code duplication for other functions (e.g.
// ancestral sampling) that
// make use of conditional probability vectors, but would not benefit from
// subtree caching
void _TheTree::FillInConditionals(_DataSetFilter *theFilter,
                                  _Parameter *iNodeCache, _SimpleList *tcc) {
  if (!tcc) {
    return;
  }

  long alphabetDimension = theFilter->GetDimension(),
       siteCount = theFilter->NumberDistinctSites();

  for (long nodeID = 0; nodeID < flatTree.lLength; nodeID++) {
    _Parameter *conditionals =
        iNodeCache + (nodeID * siteCount) * alphabetDimension;
    long currentTCCIndex = siteCount * nodeID,
         currentTCCBit = currentTCCIndex % _HY_BITMASK_WIDTH_;

    currentTCCIndex /= _HY_BITMASK_WIDTH_;
    for (long siteID = 0; siteID < siteCount;
         siteID++, conditionals += alphabetDimension) {
      if (siteID && (tcc->lData[currentTCCIndex] &
                     bitMaskArray.masks[currentTCCBit]) > 0) {
        for (long k = 0; k < alphabetDimension; k++) {
          conditionals[k] = conditionals[k - alphabetDimension];
        }
      }
      if (++currentTCCBit == _HY_BITMASK_WIDTH_) {
        currentTCCBit = 0;
        currentTCCIndex++;
      }

    }
  }
}

//______________________________________________________________________________
#ifdef MDSOCL
_Parameter _TheTree::OCLLikelihoodEvaluator(_SimpleList &updateNodes,
                                            _DataSetFilter *theFilter,
                                            _Parameter *iNodeCache,
                                            long *lNodeFlags,
                                            _GrowingVector *lNodeResolutions,
                                            _OCLEvaluator &OCLEval) {

  _SimpleList taggedInternals(flatNodes.lLength, 0, 0);
  //printf("Launching a tree in OpenCL...\n");
  return ((_Parameter) OCLEval.launchmdsocl(
      updateNodes, flatParents, flatNodes, flatCLeaves, flatLeaves, flatTree,
      theProbs, theFilter->theFrequencies, lNodeFlags, taggedInternals,
      lNodeResolutions));
}

#endif

/*----------------------------------------------------------------------------------------------------------*/
_Parameter _TheTree::VerySimpleLikelihoodEvaluator(
    _SimpleList &updateNodes, _DataSetFilter *theFilter, _Parameter *iNodeCache,
    long *lNodeFlags, _GrowingVector *lNodeResolutions) {

  // set some useful global variables
  _SimpleList taggedInternals(flatNodes.lLength, 0, 0);
  // as we are traversing the tree, we are ensuring that each node which changes
  // forces the parent to change and preps the parent cache for computations
  // by zeroing everything out

  // flatNodes is the array of indices (post-order traveral) of _internal nodes_
  // for example ((A,C)N1,D,(B,E)N2)Root will have
  // N1 (index 0), N2 (index 1), Root (index 2) in this array

  long alphabetDimension = theFilter->GetDimension(),
       // the number of characters (and the dimension of transition matrices)
      siteCount = theFilter->NumberDistinctSites();

  // how many unique sites are there
  for (long nodeID = 0; nodeID < updateNodes.lLength; nodeID++) {
    long nodeCode = updateNodes.lData[nodeID],
         parentCode = flatParents.lData[nodeCode];

    // the INDEX of the parent node for the current node;
    // this list (a member of _TheTree is computed when the tree is created
    // for ((A,C)N1,D,(B,E)N2)Root this array will store
    // 0 (A), 0(C), 2(D), 1(B), 1(E), 2(N1), 2(N3), -1 (root)
    bool isLeaf = nodeCode < flatLeaves.lLength;

    if (!isLeaf) {
      nodeCode -= flatLeaves.lLength;
    }

    _Parameter *parentConditionals =
        iNodeCache + (parentCode * siteCount) * alphabetDimension;

    // this is a convenience pointer into the global cache
    // each node will have siteCount*alphabetDimension contiguous doubles

    // mark the parent for update and clear its conditionals if needed
    if (taggedInternals.lData[parentCode] == 0) {
      taggedInternals.lData[parentCode] = 1; // only do this once
      for (long k = 0, k3 = 0; k < siteCount; k++)
        for (long k2 = 0; k2 < alphabetDimension; k2++) {
          parentConditionals[k3++] = 1.0;
        }
    }

    // in the host code the transition matrix is retrieved from the _CalcNode
    // object
    // in the devide code there will probably be another array to grab it from
    _Parameter *tMatrix =
        (isLeaf ? ((_CalcNode *)flatCLeaves(nodeCode))
                : ((_CalcNode *)flatTree(nodeCode)))->GetCompExp(0)->theData;

    // the vector of conditional probabilities for
    // THIS node (nodeCode)
    _Parameter *childVector; 

    if (!isLeaf) {
      childVector = iNodeCache + (nodeCode * siteCount) * alphabetDimension;
    }

    // if THIS node is internal, simply look up the conditional vector in the
    // same cache
    // as the parent (just a different offset)
    for (long siteID = 0; siteID < siteCount;
         siteID++, parentConditionals += alphabetDimension) {
      _Parameter sum = 0.0;

      // the leaves do NOT have a conditional probability vector because
      // most of them are fully resolved
      // For example the codon AAG is going to map to
      // (AAA)0,(AAC)0,(AAG)1,.....,(TTT)0, so we can
      // just store this as index 2, which says to put a '1' in the 2-nd
      // index of the array (and assume the rest are 0)
      // this is stored in lNodeFlags
      // For ambiguous codons, e.g. AAR = {AAA,AAG}, the host code will
      // resolve this to 1,0,1,0,...0 at prep time (because this will never
      // change)
      // and stores it in the lNodeResolutions (say at index K) and puts
      // -K-1 into lNodeFlags, so that we now have to look it up here
      if (isLeaf) {
        long siteState = lNodeFlags[nodeCode * siteCount + siteID];

        // a single character state; sweep down the appropriate column
        // note that one of the loops (over child state) drops out, since
        // there is only one state
        if (siteState >= 0) {
          long matrixIndex = siteState;
          for (long k = 0; k < alphabetDimension;
               k++, matrixIndex += alphabetDimension) {
            parentConditionals[k] *= tMatrix[matrixIndex];
          }
          continue; // nothing else to do, move to the next site
        } else {
          childVector =
              lNodeResolutions->theData + (-siteState - 1) * alphabetDimension;
        }
        // look up the resolution for the ambugious node -- this will have to be
        // on the device as well
        // but can be in constant memory
      }

      _Parameter *matrixPointer = tMatrix;

      for (long p = 0; p < alphabetDimension; p++) {
        _Parameter accumulator = 0.0;

        for (long c = 0; c < alphabetDimension; c++) {
          accumulator += matrixPointer[c] * childVector[c];
        }

        matrixPointer += alphabetDimension;

        sum += (parentConditionals[p] *= accumulator);
      }

      // if (sum < small_number) -- handle underflow

      childVector += alphabetDimension;
      // shift the position in childvector to the next site
    }
  }

  // now just process the root and return the likelihood

  _Parameter *rootConditionals =
                 iNodeCache +
                 alphabetDimension * ((flatTree.lLength - 1) * siteCount),
             // the root is always the LAST internal node in all lists
      result = 0.0;

  for (long siteID = 0; siteID < siteCount; siteID++) {
    _Parameter accumulator = 0.;
    for (long p = 0; p < alphabetDimension; p++, rootConditionals++) {
      accumulator += *rootConditionals * theProbs[p];
    }

    //theProbs is a member variable of the tree, which basically determines
    //  what probability there is to observe a given character at the root
    //  in simple cases it is fixed for the duration of optimization, but for
    //  more
    //  complex models it may change from iteration to iteration
     
    result += log(accumulator) * theFilter->theFrequencies[siteID];
    // correct for the fact that identical alignment columns may appear more
    // than once
  }

  return result;
}

//______________________________________________________________________________
// the updateNodes flags the nodes (leaves followed by inodes in the same
// order as flatLeaves and flatNodes)
// that must be recomputed
_Parameter _TheTree::ComputeTreeBlockByBranch(
    _SimpleList &siteOrdering, _SimpleList &updateNodes, _SimpleList *tcc,
    _DataSetFilter *theFilter, _Parameter *iNodeCache, long *lNodeFlags,
    _Parameter *scalingAdjustments, _GrowingVector *lNodeResolutions,
    long &overallScaler, long siteFrom, long siteTo, long catID,
    _Parameter *storageVec, long *siteCorrectionCounts, long setBranch,
    long *setBranchTo) {

  // process the leaves first
  _SimpleList taggedInternals(flatNodes.lLength, 0, 0);
  long alphabetDimension = theFilter->GetDimension(),
       siteCount = theFilter->NumberDistinctSites(),
       alphabetDimensionmod4 = alphabetDimension - alphabetDimension % 4;

  _CalcNode *currentTreeNode;
  long localScalerChange = 0;

  if (siteTo > siteCount) {
    siteTo = siteCount;
  }

  for (unsigned long nodeID = 0; nodeID < updateNodes.lLength; nodeID++) {
    long nodeCode = updateNodes.lData[nodeID],
         parentCode = flatParents.lData[nodeCode];

    bool isLeaf = nodeCode < flatLeaves.lLength;

    if (!isLeaf) {
      nodeCode -= flatLeaves.lLength;
    }

    _Parameter *parentConditionals =
        iNodeCache + (siteFrom + parentCode * siteCount) * alphabetDimension;

    // mark the parent for update and clear its conditionals if needed
    if (taggedInternals.lData[parentCode] == 0) {
      taggedInternals.lData[parentCode] = 1;
      _Parameter _hprestrict_ *localScalingFactor =
          scalingAdjustments + parentCode * siteCount;

      bool matchSet = (parentCode == setBranch);

      if (alphabetDimension == 4) {
        long k3 = 0;
        if (matchSet)
          for (long k = siteFrom; k < siteTo; k++, k3 += 4) {
            parentConditionals[k3] = 0.;
            parentConditionals[k3 + 1] = 0.;
            parentConditionals[k3 + 2] = 0.;
            parentConditionals[k3 + 3] = 0.;
            parentConditionals[k3 + setBranchTo[siteOrdering.lData[k]]] =
                localScalingFactor[k];
          }
        else
          for (long k = siteFrom; k < siteTo; k++, k3 += 4) {
            _Parameter scaler = localScalingFactor[k];
            parentConditionals[k3] = scaler;
            parentConditionals[k3 + 1] = scaler;
            parentConditionals[k3 + 2] = scaler;
            parentConditionals[k3 + 3] = scaler;
          }
      } else {
        long k3 = 0;
        if (matchSet) {
          for (long k = siteFrom; k < siteTo; k++) {
            for (long k2 = 0; k2 < alphabetDimension; k2++) {
              parentConditionals[k3 + k2] = 0.;
            }
            parentConditionals[k3 + setBranchTo[siteOrdering.lData[k]]] =
                localScalingFactor[k];
            k3 += alphabetDimension;
          }
        } else {
          for (long k = siteFrom; k < siteTo; k++) {
            _Parameter scaler = localScalingFactor[k];
            for (long k2 = 0; k2 < alphabetDimension; k2++, k3++) {
              parentConditionals[k3] = scaler;
            }
          }
        }
      }
    }

    currentTreeNode = isLeaf ? ((_CalcNode *)flatCLeaves(nodeCode))
                             : ((_CalcNode *)flatTree(nodeCode));

    _Parameter *_hprestrict_ transitionMatrix =
        currentTreeNode->GetCompExp(catID)->theData;
    _Parameter *childVector, *lastUpdatedSite;

    if (!isLeaf) {
      childVector =
          iNodeCache + (siteFrom + nodeCode * siteCount) * alphabetDimension;
    }

    long currentTCCIndex, currentTCCBit, parentTCCIIndex, parentTCCIBit;

    if (tcc) {
      parentTCCIIndex = siteCount * parentCode + siteFrom;
      parentTCCIBit = parentTCCIIndex % _HY_BITMASK_WIDTH_;
      parentTCCIIndex = parentTCCIIndex / _HY_BITMASK_WIDTH_;
      if (!isLeaf) {
        currentTCCIndex = siteCount * nodeCode + siteFrom;
        currentTCCBit = currentTCCIndex % _HY_BITMASK_WIDTH_;
        currentTCCIndex /= _HY_BITMASK_WIDTH_;
      }
    }

    //long successiveSkips = 0;
    for (long siteID = siteFrom; siteID < siteTo;
         siteID++, parentConditionals += alphabetDimension) {
      if (tcc) {
        if (parentTCCIBit == _HY_BITMASK_WIDTH_) {
          parentTCCIBit = 0;
          parentTCCIIndex++;
        }

        if (siteID > siteFrom && (tcc->lData[parentTCCIIndex] &
                                  bitMaskArray.masks[parentTCCIBit]) > 0) {
          if (!isLeaf) {
            childVector += alphabetDimension;
            if (++currentTCCBit == _HY_BITMASK_WIDTH_) {
              currentTCCBit = 0;
              currentTCCIndex++;
            }
          }
          parentTCCIBit++;
          
          //if (likeFuncEvalCallCount == 4 || likeFuncEvalCallCount == 5) {
          //    printf ("\nSKIPPED site %ld @ eval %ld\n", siteID,
          //likeFuncEvalCallCount);
          //}

          continue;
        }
        parentTCCIBit++;
      }

      _Parameter *tMatrix = transitionMatrix, sum = 0.0;

      char didScale = 0;

      if (isLeaf) {
        long siteState;

        if (setBranch == nodeCode + flatTree.lLength) {
          siteState = setBranchTo[siteOrdering.lData[siteID]];
        } else {
          siteState =
              lNodeFlags[nodeCode * siteCount + siteOrdering.lData[siteID]];
        }
        if (siteState >= 0) {
          // a single character state; sweep down the appropriate column
          tMatrix += siteState;
          if (alphabetDimension == 4) {
            parentConditionals[0] *= tMatrix[0];
            parentConditionals[1] *= tMatrix[4];
            parentConditionals[2] *= tMatrix[8];
            parentConditionals[3] *= tMatrix[12];
          } else {
            long k = 0;
            for (; k < alphabetDimensionmod4;
                 k += 4, tMatrix += alphabetDimension + alphabetDimension +
                                    alphabetDimension + alphabetDimension) {
              parentConditionals[k] *= tMatrix[0];
              parentConditionals[k + 1] *= tMatrix[alphabetDimension];
              parentConditionals[k + 2] *=
                  tMatrix[alphabetDimension + alphabetDimension];
              parentConditionals[k + 3] *= tMatrix[
                  alphabetDimension + alphabetDimension + alphabetDimension];
            }
            for (; k < alphabetDimension; k++, tMatrix += alphabetDimension) {
              parentConditionals[k] *= *tMatrix;
            }
          }
          continue;
        } else {
          childVector =
              lNodeResolutions->theData + (-siteState - 1) * alphabetDimension;
        }
      } else {
        if (tcc) {
          if ((tcc->lData[currentTCCIndex] &
               bitMaskArray.masks[currentTCCBit]) > 0 && siteID > siteFrom) {
            for (long k = 0; k < alphabetDimension; k++) {
              childVector[k] = lastUpdatedSite[k];
            }
          }
          if (++currentTCCBit == _HY_BITMASK_WIDTH_) {
            currentTCCBit = 0;
            currentTCCIndex++;
          }
          lastUpdatedSite = childVector;
        }
      }

      if (alphabetDimension == 4) { // special case for nuc data

        _handle4x4_pruning_case(childVector, tMatrix, parentConditionals);

        // handle scaling if necessary
        // the check for sum > 0.0 is necessary for 'degenerate' log-L functions
        // (-infinity)

        sum = parentConditionals[0] + parentConditionals[1] +
              parentConditionals[2] + parentConditionals[3];
        if (sum < _lfScalingFactorThreshold && sum > 0.0) {
          _Parameter tryScale =
              scalingAdjustments[parentCode * siteCount + siteID] *
              _lfScalerUpwards;
          if (tryScale < HUGE_VAL) {
            parentConditionals[0] *= _lfScalerUpwards;
            parentConditionals[1] *= _lfScalerUpwards;
            parentConditionals[2] *= _lfScalerUpwards;
            parentConditionals[3] *= _lfScalerUpwards;
            localScalerChange +=
                theFilter->theFrequencies[siteOrdering.lData[siteID]];
            scalingAdjustments[parentCode * siteCount + siteID] = tryScale;
            didScale = 1;
          }
        } else {
          if (sum > _lfScalerUpwards) {
            scalingAdjustments[parentCode * siteCount + siteID] *=
                _lfScalingFactorThreshold;
            parentConditionals[0] *= _lfScalingFactorThreshold;
            parentConditionals[1] *= _lfScalingFactorThreshold;
            parentConditionals[2] *= _lfScalingFactorThreshold;
            parentConditionals[3] *= _lfScalingFactorThreshold;
            localScalerChange -=
                theFilter->theFrequencies[siteOrdering.lData[siteID]];
            didScale = -1;
          }
        }

        childVector += 4;
      } else {
        _Parameter sum = 0.0;

#ifndef _SLKP_SSE_VECTORIZATION_

        if (alphabetDimension > alphabetDimensionmod4) {

/*_Parameter maxParentP = _lfScalingFactorThreshold;

for (long p = 0; p < alphabetDimension; p++)
    if (parentConditionals[p] > maxParentP) {
        maxParentP = parentConditionals[p];
    }

maxParentP *= DBL_EPSILON;  */

#ifdef _SLKP_USE_SSE_INTRINSICS
          double buffer[2] __attribute__((aligned(16)));
#endif
#ifdef _SLKP_USE_AVX_INTRINSICS
          double buffer[4] __attribute__((aligned(32)));
#endif
          for (long p = 0; p < alphabetDimension; p++) {
            /*if (parentConditionals[p] < maxParentP) {
                tMatrix               += alphabetDimension;
                continue;
            }*/

            _Parameter accumulator = 0.0;

#ifdef _SLKP_USE_SSE_INTRINSICS

            __m128d buffer1, buffer2, buffer3 = _mm_setzero_pd(),
                                      buffer4 = _mm_setzero_pd(), load1, load2,
                                      load3, load4;

            if (((long int) tMatrix & 0x1111b) == 0 &&
                ((long int) childVector & 0x1111b) == 0) {
              for (long c = 0; c < alphabetDimensionmod4; c += 4) {
                load1 = _mm_load_pd(tMatrix + c);
                load2 = _mm_load_pd(tMatrix + c + 2);
                load3 = _mm_load_pd(childVector + c);
                load4 = _mm_load_pd(childVector + c + 2);
                buffer1 = _mm_mul_pd(load1, load3);
                buffer2 = _mm_mul_pd(load2, load4);
                buffer3 = _mm_add_pd(buffer1, buffer3);
                buffer4 = _mm_add_pd(buffer2, buffer4);
              }
            } else {
              for (long c = 0; c < alphabetDimensionmod4; c += 4) {
                load1 = _mm_loadu_pd(tMatrix + c);
                load2 = _mm_loadu_pd(tMatrix + c + 2);
                load3 = _mm_loadu_pd(childVector + c);
                load4 = _mm_loadu_pd(childVector + c + 2);
                buffer1 = _mm_mul_pd(load1, load3);
                buffer2 = _mm_mul_pd(load2, load4);
                buffer3 = _mm_add_pd(buffer1, buffer3);
                buffer4 = _mm_add_pd(buffer2, buffer4);
              }

            }

            buffer3 = _mm_add_pd(buffer3, buffer4);
            _mm_store_pd(buffer, buffer3);
            accumulator = buffer[0] + buffer[1];
#elif defined _SLKP_USE_AVX_INTRINSICS

            __m256d buffer1, buffer3 = _mm256_setzero_pd(), load1, load3;

            if (((long int) tMatrix & 0x11111b) == 0 &&
                ((long int) childVector & 0x11111b) == 0) {
              for (long c = 0; c < alphabetDimensionmod4; c += 4) {
                load1 = _mm256_load_pd(tMatrix + c);
                load3 = _mm256_load_pd(childVector + c);
                buffer1 = _mm256_mul_pd(load1, load3);
                buffer3 = _mm256_add_pd(buffer1, buffer3);
              }
            } else {
              for (long c = 0; c < alphabetDimensionmod4; c += 4) {
                load1 = _mm256_loadu_pd(tMatrix + c);
                load3 = _mm256_loadu_pd(childVector + c);
                buffer1 = _mm256_mul_pd(load1, load3);
                buffer3 = _mm256_add_pd(buffer1, buffer3);
              }

            }

            _mm256_store_pd(buffer, buffer3);
            accumulator = buffer[0] + buffer[1] + buffer[2] + buffer[3];
#else
            for (long c = 0; c < alphabetDimensionmod4; c += 4) {
              // 4 - unroll the loop
              _Parameter pr1 = tMatrix[c] * childVector[c],
                         pr2 = tMatrix[c + 1] * childVector[c + 1],
                         pr3 = tMatrix[c + 2] * childVector[c + 2],
                         pr4 = tMatrix[c + 3] * childVector[c + 3];
              pr1 += pr2;
              pr3 += pr4;
              accumulator += pr1 + pr3;
            }
#endif

            for (long c = alphabetDimensionmod4; c < alphabetDimension; c++) {
              accumulator += tMatrix[c] * childVector[c];
            }

            tMatrix += alphabetDimension;
            sum += (parentConditionals[p] *= accumulator);
          }
        } else
          for (long p = 0; p < alphabetDimension; p++) {
            _Parameter accumulator = 0.0;

            for (long c = 0; c < alphabetDimensionmod4;
                 c += 4) // 4 - unroll the loop
              accumulator += tMatrix[c] * childVector[c] +
                             tMatrix[c + 1] * childVector[c + 1] +
                             tMatrix[c + 2] * childVector[c + 2] +
                             tMatrix[c + 3] * childVector[c + 3];

            tMatrix += alphabetDimension;
            sum += (parentConditionals[p] *= accumulator);
          }
#else
        for (long p = 0; p < alphabetDimension; p++) {
          _Parameter accumulator = 0.0;

          for (long c = 0; c < alphabetDimension; c++) {
            accumulator += tMatrix[c] * childVector[c];
          }

          tMatrix += alphabetDimension;
          sum += (parentConditionals[p] *= accumulator);
        }

#endif

        if (sum < _lfScalingFactorThreshold && sum > 0.0) {
          _Parameter tryScale =
              scalingAdjustments[parentCode * siteCount + siteID] *
              _lfScalerUpwards;
          if (tryScale < HUGE_VAL) {
            scalingAdjustments[parentCode * siteCount + siteID] = tryScale;
            for (long c = 0; c < alphabetDimension; c++) {
              parentConditionals[c] *= _lfScalerUpwards;
            }

            localScalerChange +=
                theFilter->theFrequencies[siteOrdering.lData[siteID]];
            didScale = 1;
          }
        } else {
          if (sum > _lfScalerUpwards) {
            scalingAdjustments[parentCode * siteCount + siteID] *=
                _lfScalingFactorThreshold;
            for (long c = 0; c < alphabetDimension; c++) {
              parentConditionals[c] *= _lfScalingFactorThreshold;
            }
            localScalerChange -=
                theFilter->theFrequencies[siteOrdering.lData[siteID]];
            didScale = -1;
          }
        }
        childVector += alphabetDimension;
      }

      if (didScale) {
        if (siteCorrectionCounts) {
          siteCorrectionCounts[siteOrdering.lData[siteID]] += didScale;
        }

        //printf ("NS: site %d node %d \n", siteOrdering.lData[siteID],
        //parentCode);

        if (tcc) {
          // look ahead to see if we need to correct for downstream cached nodes
          long cparentTCCIIndex = parentTCCIIndex,
               cparentTCCIBit = parentTCCIBit;

          _Parameter scM;
          if (didScale < 0) {
            scM = _lfScalingFactorThreshold;
          } else {
            scM = _lfScalerUpwards;
          }

          for (long sid = siteID + 1; sid < siteTo; sid++, cparentTCCIBit++) {
            if (cparentTCCIBit == _HY_BITMASK_WIDTH_) {
              cparentTCCIBit = 0;
              cparentTCCIIndex++;
            }

            if ((tcc->lData[cparentTCCIIndex] &
                 bitMaskArray.masks[cparentTCCIBit]) > 0) {
              if (siteCorrectionCounts) {
                siteCorrectionCounts[siteOrdering.lData[sid]] += didScale;
              }
              scalingAdjustments[parentCode * siteCount + sid] *= scM;
              localScalerChange +=
                  didScale * theFilter->theFrequencies[siteOrdering.lData[sid]];
            } else {
              break;
            }
          }
        }
      }
    }
  }

  // assemble the entire likelihood

  _Parameter _hprestrict_ *rootConditionals =
      iNodeCache +
      alphabetDimension * (siteFrom + (flatTree.lLength - 1) * siteCount);
  _Parameter result = 0.0;

  for (long siteID = siteFrom; siteID < siteTo; siteID++) {
    _Parameter accumulator = 0.;
    if (setBranch == flatTree.lLength - 1) {
      long rootState = setBranchTo[siteOrdering.lData[siteID]];
      accumulator += rootConditionals[rootState] * theProbs[rootState];
      rootConditionals += alphabetDimension;
    } else
      for (long p = 0; p < alphabetDimension; p++, rootConditionals++) {
        accumulator += *rootConditionals * theProbs[p];
      }

    /*#pragma omp critical
    {
        if (likeFuncEvalCallCount == 0) {
            eval_buffer [siteID] = accumulator;
        } else {
            if (!CheckEqual (eval_buffer [siteID], accumulator)) {
                printf ("EVAL %ld, Site %ld/%ld, LogL %g, %g (%d)\n",
    likeFuncEvalCallCount, siteFrom, siteID, accumulator, eval_buffer
    [siteID],localScalerChange);
            }
            eval_buffer [siteID] = accumulator;
        }
    }*/

    if (storageVec) {
      storageVec[siteOrdering.lData[siteID]] = accumulator;
    } else {
      if (accumulator <= 0.0) {
        result = -A_LARGE_NUMBER;
#pragma omp critical
        {
          ReportWarning(
              _String("Site ") & (1L + siteOrdering.lData[siteID]) &
              " evaluated to a 0 probability in ComputeTreeBlockByBranch");
        }
        break;
      }
      if (theFilter->theFrequencies[siteOrdering.lData[siteID]] > 1) {
        result += log(accumulator) *
                  theFilter->theFrequencies[siteOrdering.lData[siteID]];
      } else {
        result += log(accumulator);
      }
    }
  }

  if (!storageVec && localScalerChange) {
#pragma omp atomic
    overallScaler += localScalerChange;
  }

  return result;
}

#ifdef _SLKP_DEBUG_

//______________________________________________________________________________
void echoNodeList(_SimpleList &theNodes, _SimpleList &leaves,
                  _SimpleList &iNodes) {
  for (long n = 0; n < theNodes.lLength; n++) {
    node<long> *nd = (node<long> *)(theNodes(n) < leaves.lLength
                                        ? leaves(theNodes(n))
                                        : iNodes(theNodes(n) - leaves.lLength));
    printf("%d %d %s\n", n, theNodes(n),
           LocateVar(nd->in_object)->GetName()->sData);
  }
}

#endif

//______________________________________________________________________________
void _TheTree::ComputeBranchCache(
    _SimpleList &siteOrdering, long brID, _Parameter *cache,
    _Parameter *iNodeCache, _DataSetFilter *theFilter, long *lNodeFlags,
    _Parameter *scalingAdjustments, long *siteCorrectionCounts,
    _GrowingVector *lNodeResolutions, long &overallScaler, long siteFrom,
    long siteTo, long catID, _SimpleList *tcc, _Parameter *siteRes) {

  //printf ("ComputeBranchCache\n");

  _SimpleList taggedNodes(flatLeaves.lLength + flatNodes.lLength, 0, 0),
      nodesToProcess, rootPath;

  long myParent = brID - flatLeaves.lLength,
       alphabetDimension = theFilter->GetDimension(),
       alphabetDimensionmod4 = alphabetDimension - alphabetDimension % 4,
       siteCount = theFilter->NumberDistinctSites();

  if (siteTo > siteCount) {
    siteTo = siteCount;
  }

  do {
    taggedNodes.lData[myParent + flatLeaves.lLength] = 1;
    myParent = flatParents.lData[myParent + flatLeaves.lLength];
  } while (myParent >= 0);

  for (unsigned long k = 0; k < flatLeaves.lLength + flatNodes.lLength; k++) {
    myParent = flatParents.lData[k];
    if (taggedNodes.lData[myParent + flatLeaves.lLength] == 1 &&
        taggedNodes.lData[k] == 0) {
      if (myParent != brID - flatLeaves.lLength) {
        nodesToProcess << k;
      }
    }
    if (taggedNodes.lData[k]) {
      rootPath << k;
    }
  }

  //printf ("ComputeBranchCache at branch %d; siteOdering %s\n",
  //      brID, _String((_String*)siteOrdering.toStr()).sData);

  //echoNodeList (rootPath,flatLeaves,flatNodes );
  //echoNodeList (nodesToProcess,flatLeaves,flatNodes);

  _Parameter *state = cache + alphabetDimension * siteFrom, *childVector;

  long localScalerChange = 0;

  // first populate the downward looking vector of conditionals

  if (brID < flatLeaves.lLength) { // a leaf
    for (long siteID = siteFrom; siteID < siteTo;
         siteID++, state += alphabetDimension) {
      long siteState =
          lNodeFlags[brID * siteCount + siteOrdering.lData[siteID]];
      if (siteState >= 0)
          // a single character state; sweep down the appropriate column
          {
        for (long s = 0; s < alphabetDimension; s++) {
          state[s] = 0.;
        }
        state[siteState] = 1.;
      } else {
        childVector =
            lNodeResolutions->theData + (-siteState - 1) * alphabetDimension;
        for (long s = 0; s < alphabetDimension; s++) {
          state[s] = childVector[s];
        }
      }
    }
  } else { // an internal branch
    long nodeCode = brID - flatLeaves.lLength;
    _Parameter *lastUpdated =
        iNodeCache + (nodeCode * siteCount + siteFrom) * alphabetDimension;

    long currentTCCIndex, currentTCCBit;

    if (tcc) {
      currentTCCIndex = siteCount * nodeCode + siteFrom;
      currentTCCBit = currentTCCIndex % _HY_BITMASK_WIDTH_;
      currentTCCIndex /= _HY_BITMASK_WIDTH_;
    }

    for (long siteID = siteFrom; siteID < siteTo;
         siteID++, state += alphabetDimension) {
      if (tcc) {
        if ((tcc->lData[currentTCCIndex] & bitMaskArray.masks[currentTCCBit]) ==
            0) {
          lastUpdated =
              iNodeCache + (nodeCode * siteCount + siteID) * alphabetDimension;
        }
      }

      for (long s = 0; s < alphabetDimension; s++) {
        state[s] = lastUpdated[s];
      }

      if (tcc) {
        if (++currentTCCBit == _HY_BITMASK_WIDTH_) {
          currentTCCBit = 0;
          currentTCCIndex++;
        }
      } else {
        lastUpdated += alphabetDimension;
      }
    }
  }

  state = cache + alphabetDimension * siteCount;

  taggedNodes.Populate(flatTree.lLength, 0, 0);
  rootPath.Flip();

  for (long nodeID = 0; nodeID < nodesToProcess.lLength + rootPath.lLength - 2;
       nodeID++) {
    bool notPassedRoot = nodeID < nodesToProcess.lLength;

    long nodeCode =
             notPassedRoot ? nodesToProcess.lData[nodeID]
                           : rootPath.lData[nodeID - nodesToProcess.lLength],
         parentCode =
             notPassedRoot
                 ? flatParents.lData[nodeCode]
                 : (rootPath.lData[nodeID - nodesToProcess.lLength + 1] -
                    flatLeaves.lLength);

    bool isLeaf = nodeCode < flatLeaves.lLength;

    if (!isLeaf) {
      nodeCode -= flatLeaves.lLength;
    }

    _Parameter *parentConditionals =
        iNodeCache + (siteFrom + parentCode * siteCount) * alphabetDimension;
    if (taggedNodes.lData[parentCode] == 0)
        // mark the parent for update and clear its conditionals if needed
        {
      taggedNodes.lData[parentCode] = 1;
      _Parameter _hprestrict_ *localScalingFactor =
          scalingAdjustments + parentCode * siteCount;
      if (alphabetDimension == 4) {
        long k3 = 0;
        for (long k = siteFrom; k < siteTo; k++, k3 += 4) {
          _Parameter scaler = localScalingFactor[k];
          parentConditionals[k3] = scaler;
          parentConditionals[k3 + 1] = scaler;
          parentConditionals[k3 + 2] = scaler;
          parentConditionals[k3 + 3] = scaler;
        }
      } else {
        long k3 = 0;
        for (long k = siteFrom; k < siteTo; k++) {
          _Parameter scaler = localScalingFactor[k];
          for (long k2 = 0; k2 < alphabetDimension; k2++, k3++) {
            parentConditionals[k3] = scaler;
          }
        }
      }
    }

    _CalcNode *currentTreeNode =
        isLeaf ? ((_CalcNode *)flatCLeaves(nodeCode))
               : ((_CalcNode *)flatTree(notPassedRoot ? nodeCode : parentCode));

    _Parameter *_hprestrict_ transitionMatrix =
        currentTreeNode->GetCompExp(catID)->theData;

    _Parameter *childVector, *lastUpdatedSite;

    if (!isLeaf) {
      lastUpdatedSite = childVector =
          iNodeCache + (siteFrom + nodeCode * siteCount) * alphabetDimension;
    }

    long currentTCCIndex, currentTCCBit, parentTCCIIndex, parentTCCIBit;

    if (tcc) {
      parentTCCIIndex = siteCount * parentCode + siteFrom;
      parentTCCIBit = parentTCCIIndex % _HY_BITMASK_WIDTH_;
      parentTCCIIndex = parentTCCIIndex / _HY_BITMASK_WIDTH_;
      if (!isLeaf) {
        currentTCCIndex = siteCount * nodeCode + siteFrom;
        currentTCCBit = currentTCCIndex % _HY_BITMASK_WIDTH_;
        currentTCCIndex /= _HY_BITMASK_WIDTH_;
      }
    }

    for (long siteID = siteFrom; siteID < siteTo;
         siteID++, parentConditionals += alphabetDimension) {
      _Parameter *tMatrix = transitionMatrix;

      char canScale = !notPassedRoot;

      if (isLeaf) {
        long siteState =
            lNodeFlags[nodeCode * siteCount + siteOrdering.lData[siteID]];
        if (siteState >= 0)
            // a single character state; sweep down the appropriate column
            {
          if (alphabetDimension == 4) { // special case for nuc data
            parentConditionals[0] *= tMatrix[siteState];
            parentConditionals[1] *= tMatrix[siteState + 4];
            parentConditionals[2] *= tMatrix[siteState + 8];
            parentConditionals[3] *= tMatrix[siteState + 12];
          } else {
            tMatrix += siteState;
            for (long k = 0; k < alphabetDimension;
                 k++, tMatrix += alphabetDimension) {
              parentConditionals[k] *= *tMatrix;
            }

          }
          continue;
        } else {
          childVector =
              lNodeResolutions->theData + (-siteState - 1) * alphabetDimension;
        }
        canScale = 0;
      } else {
        if (tcc && notPassedRoot) {
          if ((tcc->lData[currentTCCIndex] &
               bitMaskArray.masks[currentTCCBit]) > 0 && siteID > siteFrom)
            // the value of this conditional vector needs to be copied from a
            // previously stored site
            // subtree duplication
            for (long k = 0; k < alphabetDimension; k++) {
              childVector[k] = lastUpdatedSite[k];
            }
          else {
            lastUpdatedSite = childVector;
          }

          if (++currentTCCBit == _HY_BITMASK_WIDTH_) {
            currentTCCBit = 0;
            currentTCCIndex++;
          }
          if (++parentTCCIBit == _HY_BITMASK_WIDTH_) {
            parentTCCIBit = 0;
            parentTCCIIndex++;
          }
        }
      }

      _Parameter sum = .0;
      char didScale = 0;

      if (alphabetDimension == 4) { // special case for nuc data
        _handle4x4_pruning_case(childVector, tMatrix, parentConditionals);

        if (canScale) {
          sum = parentConditionals[0] + parentConditionals[1] +
                parentConditionals[2] + parentConditionals[3];
          if (sum < _lfScalingFactorThreshold && sum > 0.0) {
            _Parameter tryScale =
                scalingAdjustments[nodeCode * siteCount + siteID] *
                _lfScalerUpwards;
            if (tryScale < HUGE_VAL) {
              parentConditionals[0] *= _lfScalerUpwards;
              parentConditionals[1] *= _lfScalerUpwards;
              parentConditionals[2] *= _lfScalerUpwards;
              parentConditionals[3] *= _lfScalerUpwards;

              localScalerChange +=
                  theFilter->theFrequencies[siteOrdering.lData[siteID]];
              didScale = 1;
            }
          } else {
            if (sum > _lfScalerUpwards) {
              parentConditionals[0] *= _lfScalingFactorThreshold;
              parentConditionals[1] *= _lfScalingFactorThreshold;
              parentConditionals[2] *= _lfScalingFactorThreshold;
              parentConditionals[3] *= _lfScalingFactorThreshold;

              localScalerChange -=
                  theFilter->theFrequencies[siteOrdering.lData[siteID]];
              didScale = -1;
            }
          }
        }
        childVector += 4;
      } else {
#ifdef _SLKP_USE_SSE_INTRINSICS
        double buffer[2] __attribute__((aligned(16)));
#endif
#ifndef _SLKP_SSE_VECTORIZATION_
        for (long p = 0; p < alphabetDimension; p++) {
          _Parameter accumulator = 0.0;
#ifdef _SLKP_USE_SSE_INTRINSICS

          __m128d buffer1, buffer2, buffer3 = _mm_setzero_pd(),
                                    buffer4 = _mm_setzero_pd(), load1, load2,
                                    load3, load4;

          if (((long int) tMatrix & 0x1111b) == 0 &&
              ((long int) childVector & 0x1111b) == 0) {
            for (long c = 0; c < alphabetDimensionmod4; c += 4) {
              load1 = _mm_load_pd(tMatrix + c);
              load2 = _mm_load_pd(tMatrix + c + 2);
              load3 = _mm_load_pd(childVector + c);
              load4 = _mm_load_pd(childVector + c + 2);
              buffer1 = _mm_mul_pd(load1, load3);
              buffer2 = _mm_mul_pd(load2, load4);
              buffer3 = _mm_add_pd(buffer1, buffer3);
              buffer4 = _mm_add_pd(buffer2, buffer4);
            }
          } else {
            for (long c = 0; c < alphabetDimensionmod4; c += 4) {
              load1 = _mm_loadu_pd(tMatrix + c);
              load2 = _mm_loadu_pd(tMatrix + c + 2);
              load3 = _mm_loadu_pd(childVector + c);
              load4 = _mm_loadu_pd(childVector + c + 2);
              buffer1 = _mm_mul_pd(load1, load3);
              buffer2 = _mm_mul_pd(load2, load4);
              buffer3 = _mm_add_pd(buffer1, buffer3);
              buffer4 = _mm_add_pd(buffer2, buffer4);
            }

          }

          buffer3 = _mm_add_pd(buffer3, buffer4);
          _mm_store_pd(buffer, buffer3);
          accumulator = buffer[0] + buffer[1];

#else
          for (long c = 0; c < alphabetDimensionmod4;
               c += 4) // 4 - unroll the loop
              {
            _Parameter pr1 = tMatrix[c] * childVector[c],
                       pr2 = tMatrix[c + 1] * childVector[c + 1],
                       pr3 = tMatrix[c + 2] * childVector[c + 2],
                       pr4 = tMatrix[c + 3] * childVector[c + 3];
            pr1 += pr2;
            pr3 += pr4;
            accumulator += pr1 + pr3;
          }
#endif

          for (long c = alphabetDimensionmod4; c < alphabetDimension; c++) {
            accumulator += tMatrix[c] * childVector[c];
          }

          tMatrix += alphabetDimension;
          sum += (parentConditionals[p] *= accumulator);
        }
#else
        for (long p = 0; p < alphabetDimension; p++) {
          _Parameter accumulator = 0.0;

          for (long c = 0; c < alphabetDimension; c++) { // 4 - unroll the loop
            accumulator += tMatrix[c] * childVector[c];
          }

          tMatrix += alphabetDimension;
          sum += (parentConditionals[p] *= accumulator);
        }
#endif
        if (canScale) {
          if (sum < _lfScalingFactorThreshold && sum > 0.0) {
            _Parameter tryScale =
                scalingAdjustments[nodeCode * siteCount + siteID] *
                _lfScalerUpwards;
            if (tryScale < HUGE_VAL) {
              for (long c = 0; c < alphabetDimension; c++) {
                parentConditionals[c] *= _lfScalerUpwards;
              }

              localScalerChange +=
                  theFilter->theFrequencies[siteOrdering.lData[siteID]];
              didScale = 1;
            }
          } else {
            if (sum > _lfScalerUpwards) {
              for (long c = 0; c < alphabetDimension; c++) {
                parentConditionals[c] *= _lfScalingFactorThreshold;
              }

              localScalerChange -=
                  theFilter->theFrequencies[siteOrdering.lData[siteID]];
              didScale = -1;
            }
          }
        }
        childVector += alphabetDimension;
      }

      if (didScale && siteCorrectionCounts) {
        siteCorrectionCounts[siteOrdering.lData[siteID]] += didScale;
        siteRes[siteOrdering.lData[siteID]] *=
            didScale < 0 ? _lfScalingFactorThreshold : _lfScalerUpwards;
      }
    }
  }

  _Parameter _hprestrict_ *rootConditionals =
      iNodeCache + (rootPath.lData[rootPath.lLength - 2] - flatLeaves.lLength) *
                       siteCount * alphabetDimension;

  for (long ii = siteFrom * alphabetDimension; ii < alphabetDimension * siteTo;
       ii++) {
    state[ii] = rootConditionals[ii];
  }

  if (!siteCorrectionCounts && localScalerChange) {
#pragma omp atomic
    overallScaler += localScalerChange;
  }
}

//______________________________________________________________________________
_Parameter _TheTree::ComputeLLWithBranchCache(_SimpleList &siteOrdering,
                                              long brID, _Parameter *cache,
                                              _DataSetFilter *theFilter,
                                              long siteFrom, long siteTo,
                                              long catID,
                                              _Parameter *storageVec) {
  long alphabetDimension = theFilter->GetDimension(),
       //alphabetDimensionmod4  = alphabetDimension - alphabetDimension % 4,
      siteCount = theFilter->NumberDistinctSites();

  if (siteTo > siteCount) {
    siteTo = siteCount;
  }

  _Parameter _hprestrict_ *branchConditionals =
      cache + siteFrom * alphabetDimension;
  _Parameter _hprestrict_ *rootConditionals =
      branchConditionals + siteCount * alphabetDimension;
  _Parameter result = 0.0;

  //printf ("ComputeLLWithBranchCache @ %d catID = %d branchID = %d\n",
  //likeFuncEvalCallCount, catID, brID);

  _CalcNode *givenTreeNode =
      brID < flatLeaves.lLength
          ? (((_CalcNode **)flatCLeaves.lData)[brID])
          : (((_CalcNode **)flatTree.lData)[brID - flatLeaves.lLength]);

  _Parameter *_hprestrict_ transitionMatrix =
      givenTreeNode->GetCompExp(catID)->theData;

  for (long siteID = siteFrom; siteID < siteTo; siteID++) {
    _Parameter accumulator = 0.;

    if (alphabetDimension == 4) {
      accumulator = rootConditionals[0] * theProbs[0] *
                        (branchConditionals[0] * transitionMatrix[0] +
                         branchConditionals[1] * transitionMatrix[1] +
                         branchConditionals[2] * transitionMatrix[2] +
                         branchConditionals[3] * transitionMatrix[3]) +
                    rootConditionals[1] * theProbs[1] *
                        (branchConditionals[0] * transitionMatrix[4] +
                         branchConditionals[1] * transitionMatrix[5] +
                         branchConditionals[2] * transitionMatrix[6] +
                         branchConditionals[3] * transitionMatrix[7]) +
                    rootConditionals[2] * theProbs[2] *
                        (branchConditionals[0] * transitionMatrix[8] +
                         branchConditionals[1] * transitionMatrix[9] +
                         branchConditionals[2] * transitionMatrix[10] +
                         branchConditionals[3] * transitionMatrix[11]) +
                    rootConditionals[3] * theProbs[3] *
                        (branchConditionals[0] * transitionMatrix[12] +
                         branchConditionals[1] * transitionMatrix[13] +
                         branchConditionals[2] * transitionMatrix[14] +
                         branchConditionals[3] * transitionMatrix[15]);
      rootConditionals += 4;
    } else {
      long rmx = 0;
      for (long p = 0; p < alphabetDimension; p++, rootConditionals++) {
        _Parameter r2 = 0.;
        for (long c = 0; c < alphabetDimension; c++, rmx++) {
          r2 += branchConditionals[c] * transitionMatrix[rmx];
        }
        accumulator += *rootConditionals * theProbs[p] * r2;
      }

    }

    branchConditionals += alphabetDimension;

    if (storageVec) {
      storageVec[siteOrdering.lData[siteID]] = accumulator;
    } else {
      if (accumulator <= 0.0) {
        result = -A_LARGE_NUMBER;
#pragma omp critical
        {
          ReportWarning(
              _String("Site ") & (1L + siteOrdering.lData[siteID]) &
              " evaluated to a 0 probability in ComputeLLWithBranchCache");
        }
        break;
      }
      if (theFilter->theFrequencies[siteOrdering.lData[siteID]] > 1) {
        result += log(accumulator) *
                  theFilter->theFrequencies[siteOrdering.lData[siteID]];
      } else {
        result += log(accumulator);
      }
      //result += log(accumulator) * theFilter->theFrequencies
      //[siteOrdering.lData[siteID]];
    }
  }
  return result;
}

//______________________________________________________________________________
// the updateNodes flags the nodes (leaves followed by inodes in the same
// order as flatLeaves and flatNodes)
// that must be recomputed
_Parameter _TheTree::ComputeTwoSequenceLikelihood(
    _SimpleList &siteOrdering, _DataSetFilter *theFilter, long *lNodeFlags,
    _GrowingVector *lNodeResolutions, long siteFrom, long siteTo, long catID,
    _Parameter *storageVec) {

  // process the leaves first
  long alphabetDimension = theFilter->GetDimension(),
       siteCount = theFilter->NumberDistinctSites(),
       alphabetDimensionmod4 = alphabetDimension - alphabetDimension % 4;

  _CalcNode *theNode = ((_CalcNode *)flatCLeaves(0));
  _Parameter *_hprestrict_ transitionMatrix =
                               theNode->GetCompExp(catID)->theData,
                           result = 0.;

  if (siteTo > siteCount) {
    siteTo = siteCount;
  }

  for (long siteID = siteFrom; siteID < siteTo; siteID++) {
    _Parameter *tMatrix = transitionMatrix, sum = 0.;

    long siteState1 = lNodeFlags[siteOrdering.lData[siteID]],
         siteState2 = lNodeFlags[siteCount + siteOrdering.lData[siteID]];

    if (siteState1 >= 0) {
      // a single character state; sweep down the appropriate column
      if (siteState2 >= 0) { // both completely resolved;
        sum = tMatrix[siteState1 * alphabetDimension + siteState2];
      } else { // first resolved, second is not
        _Parameter *childVector =
            lNodeResolutions->theData + (-siteState2 - 1) * alphabetDimension;
        tMatrix += siteState1 * alphabetDimension;
        if (alphabetDimension == 4) { // special case for nuc data
          sum = tMatrix[0] * childVector[0] + tMatrix[1] * childVector[1] +
                tMatrix[2] * childVector[2] + tMatrix[3] * childVector[3];
        } else {
          for (long c = 0; c < alphabetDimensionmod4;
               c += 4) // 4 - unroll the loop
            sum += tMatrix[c] * childVector[c] +
                   tMatrix[c + 1] * childVector[c + 1] +
                   tMatrix[c + 2] * childVector[c + 2] +
                   tMatrix[c + 3] * childVector[c + 3];

          for (long c = alphabetDimensionmod4; c < alphabetDimension; c++) {
            sum += tMatrix[c] * childVector[c];
          }
        }
      }
      sum *= theProbs[siteState1];
    } else {
      if (siteState2 >= 0) { // second resolved, but not the first
        _Parameter *childVector =
            lNodeResolutions->theData + (-siteState1 - 1) * alphabetDimension;
        tMatrix += siteState2;
        if (alphabetDimension == 4) { // special case for nuc data
          sum = tMatrix[0] * childVector[0] * theProbs[0] +
                tMatrix[4] * childVector[1] * theProbs[1] +
                tMatrix[8] * childVector[2] * theProbs[2] +
                tMatrix[12] * childVector[3] * theProbs[3];

        } else {
          for (long c = 0; c < alphabetDimensionmod4;
               c += 4, tMatrix += 4 * alphabetDimension) // 4 - unroll the loop
            sum += tMatrix[0] * childVector[c] * theProbs[c] +
                   tMatrix[alphabetDimension] * childVector[c + 1] *
                       theProbs[c + 1] +
                   tMatrix[alphabetDimension + alphabetDimension] *
                       childVector[c + 2] * theProbs[c + 2] +
                   tMatrix[alphabetDimension + alphabetDimension +
                           alphabetDimension] * childVector[c + 3] *
                       theProbs[c + 3];

          for (long c = alphabetDimensionmod4; c < alphabetDimension;
               c++, tMatrix += alphabetDimension) {
            sum += tMatrix[0] * childVector[c] * theProbs[c];
          }
        }
      } else {
        // both unresolved
        _Parameter *childVector1 = lNodeResolutions->theData +
                                   (-siteState1 - 1) * alphabetDimension,
                   *childVector2 = lNodeResolutions->theData +
                                   (-siteState2 - 1) * alphabetDimension;

        if (alphabetDimension == 4) { // special case for nuc data
          sum =
              (tMatrix[0] * childVector2[0] + tMatrix[1] * childVector2[1] +
               tMatrix[2] * childVector2[2] + tMatrix[3] * childVector2[3]) *
                  childVector1[0] * theProbs[0] +
              (tMatrix[4] * childVector2[0] + tMatrix[5] * childVector2[1] +
               tMatrix[6] * childVector2[2] + tMatrix[7] * childVector2[3]) *
                  childVector1[1] * theProbs[1] +
              (tMatrix[8] * childVector2[0] + tMatrix[9] * childVector2[1] +
               tMatrix[10] * childVector2[2] + tMatrix[11] * childVector2[3]) *
                  childVector1[2] * theProbs[2] +
              (tMatrix[12] * childVector2[0] + tMatrix[13] * childVector2[1] +
               tMatrix[14] * childVector2[2] + tMatrix[15] * childVector2[3]) *
                  childVector1[3] * theProbs[3];

        } else {
          for (long r = 0; r < alphabetDimension; r++) { // 4 - unroll the loop
            if (childVector1[r] > 0.0) {
              _Parameter sum2 = 0.0;
              for (long c = 0; c < alphabetDimensionmod4;c += 4, tMatrix += 4) 
                // 4 - unroll the loop
                sum2 += tMatrix[0] * childVector2[c] +
                        tMatrix[1] * childVector2[c + 1] +
                        tMatrix[2] * childVector2[c + 2] +
                        tMatrix[3] * childVector2[c + 3];

              for (long c = alphabetDimensionmod4; c < alphabetDimension;c++, tMatrix++) {
                sum2 += tMatrix[0] * childVector2[c];
              }
              sum += sum2 * childVector1[r] * theProbs[r];
            } else {
              tMatrix += alphabetDimension;
            }
          }
        }
      }
    }
    if (storageVec) {
      storageVec[siteOrdering.lData[siteID]] = sum;
    } else {
      if (sum <= 0.0) {
        return -A_LARGE_NUMBER;
      } else {
        //printf ("%d: %g\n", siteID, sum);
        result +=
            log(sum) * theFilter->theFrequencies[siteOrdering.lData[siteID]];
      }
    }
  }

  return result;
}

//______________________________________________________________________________
// must be called initially with the root node
// dsf:                         the filter to sample from
// siteOrdering:                the map from cache ordering to actual
// pattern ordering
// currentNode:                 the node index to sample for
// nodeToIndex:                 an AVL that maps the address of an internal
// node pointed to by node<long> to its order in the tree postorder
// traversal
// iNodeCache:                  internal node likelihood caches
// results:                     the list that will store sampled strings
// parentStates:                sampled states for the parent of the current
// node
// expandedSiteMap:             a list of simple lists giving site indices
// for each unique column pattern in the alignment
// catAssignments:              a vector assigning a (partition specific)
// rate category to each site (nil if no rate variation)
// catCount:                    the number of rate classes
// this needs to be updated to deal with traversal caches!

void _TheTree::SampleAncestorsBySequence(
    _DataSetFilter *dsf, _SimpleList &siteOrdering, node<long> *currentNode,
    _AVLListX *nodeToIndex, _Parameter *iNodeCache, _List &result,
    _SimpleList *parentStates, _List &expandedSiteMap,
    _Parameter *catAssignments, long catCount) {

  long childrenCount = currentNode->get_num_nodes();

  if (childrenCount) {
    long siteCount = dsf->NumberDistinctSites(),
         alphabetDimension = dsf->GetDimension(),

         nodeIndex =
             nodeToIndex->GetXtra(nodeToIndex->Find((BaseRef) currentNode)),

         unitLength = dsf->GetUnitLength(),

         catBlockShifter =
             catAssignments ? (dsf->NumberDistinctSites() * GetINodeCount())
                            : 0;

    _CalcNode *currentTreeNode = ((_CalcNode *)flatTree(nodeIndex));
    _SimpleList sampledStates(dsf->GetSiteCount(), 0, 0);

    _Parameter *_hprestrict_ transitionMatrix =
        (catAssignments || !parentStates)
            ? nil
            : currentTreeNode->GetCompExp()->theData;

    _Parameter *_hprestrict_ conditionals =
        catAssignments
            ? nil
            : (iNodeCache + nodeIndex * siteCount * alphabetDimension);
    _Parameter *_hprestrict_ cache = new _Parameter[alphabetDimension];

    for (long pattern = 0; pattern < siteCount; pattern++) {
      _SimpleList *patternMap =
          (_SimpleList *)expandedSiteMap(siteOrdering.lData[pattern]);
      if (catAssignments) {
        long localCatID = catAssignments[siteOrdering.lData[pattern]];
        if (parentStates) {
          transitionMatrix = currentTreeNode->GetCompExp(localCatID)->theData;
        }

        conditionals =
            iNodeCache + localCatID * alphabetDimension * catBlockShifter +
            (pattern + nodeIndex * siteCount) * alphabetDimension;
      }

      for (long site = 0; site < patternMap->lLength; site++) {
        long siteID = patternMap->lData[site];

        _Parameter randVal = genrand_real2(), totalSum = 0.;

        _Parameter *_hprestrict_ matrixRow;

        if (parentStates == nil) {
          matrixRow = theProbs;
        } else {
          matrixRow = transitionMatrix +
                      parentStates->lData[siteID] * alphabetDimension;
        }

        for (long i = 0; i < alphabetDimension; i++) {
          totalSum += (cache[i] = matrixRow[i] * conditionals[i]);
        }

        randVal *= totalSum;
        totalSum = 0.0;
        long sampledChar = -1;
        while (totalSum < randVal) {
          sampledChar++;
          totalSum += cache[sampledChar];
        }
        sampledStates.lData[siteID] = sampledChar;
      }

      if (catAssignments == nil) {
        conditionals += alphabetDimension;
      }
    }

    delete[] cache;

    _SimpleList conversion;
    _AVLListXL conversionAVL(&conversion);

    _String *sampledSequence = new _String(siteCount * unitLength, true);
    _String letterValue(unitLength, false);
    for (long charIndexer = 0; charIndexer < sampledStates.lLength;
         charIndexer++) {
      dsf->ConvertCodeToLettersBuffered(
          dsf->CorrectCode(sampledStates.lData[charIndexer]), unitLength,
          letterValue.sData, &conversionAVL);
      (*sampledSequence) << letterValue;
    }
    sampledSequence->Finalize();
    result.AppendNewInstance(sampledSequence);
    //printf ("%d: %s\n", nodeIndex, sampledSequence->sData);

    for (long child = 1; child <= childrenCount; child++) {
      SampleAncestorsBySequence(dsf, siteOrdering, currentNode->go_down(child),
                                nodeToIndex, iNodeCache, result, &sampledStates,
                                expandedSiteMap, catAssignments, catCount);
    }
  }
}

//______________________________________________________________________________
// dsf:                         the filter to sample from
// siteOrdering:                the map from cache ordering to actual
// pattern ordering
// expandedSiteMap:             a list of simple lists giving site indices
// for each unique column pattern in the alignment
// iNodeCache:                  internal node likelihood caches
// catAssignments:              a vector assigning a (partition specific)
// rate category to each site
// catCount:                    the number of rate classes
// alsoDoLeaves:                if true, also return ML reconstruction of
// observed (or partially observed) sequences
_List *_TheTree::RecoverAncestralSequences(
    _DataSetFilter *dsf, _SimpleList &siteOrdering, _List &expandedSiteMap,
    _Parameter *iNodeCache, _Parameter *catAssignments, long catCount,
    long *lNodeFlags, _GrowingVector *lNodeResolutions, bool alsoDoLeaves) {

  long patternCount = dsf->NumberDistinctSites(),
       alphabetDimension = dsf->GetDimension(),
       unitLength = dsf->GetUnitLength(), iNodeCount = GetINodeCount(),
       leafCount = GetLeafCount(), siteCount = dsf->GetSiteCount(),
       allNodeCount = 0,
       stateCacheDim = (alsoDoLeaves ? (iNodeCount + leafCount) : (iNodeCount)),
       *stateCache =
           new long[patternCount * (iNodeCount - 1) * alphabetDimension],
       *leafBuffer = new long[(alsoDoLeaves ? leafCount * patternCount : 1) *
                              alphabetDimension];

  // a Patterns x Int-Nodes x CharStates integer table
  // with the best character assignment for node i given that its parent state
  // is j for a given site
  _Parameter *buffer = new _Parameter[alphabetDimension];

  // iNodeCache will be OVERWRITTEN with conditional pair (i,j) conditional
  // likelihoods
  checkPointer(stateCache);
  checkPointer(leafBuffer);

  _SimpleList taggedInternals(iNodeCount, 0, 0), postToIn;

  MapPostOrderToInOderTraversal(postToIn);

  // all nodes except the root
  allNodeCount = iNodeCount + leafCount - 1;

  for (long nodeID = 0; nodeID < allNodeCount; nodeID++) {
    long parentCode = flatParents.lData[nodeID], nodeCode = nodeID;

    bool isLeaf = nodeID < flatLeaves.lLength;

    if (!isLeaf) {
      nodeCode -= flatLeaves.lLength;
      AddBranchToForcedRecomputeList(nodeCode);
    }

    _Parameter *parentConditionals =
        iNodeCache + parentCode * alphabetDimension * patternCount;

    // mark the parent for update and clear its conditionals if needed
    if (taggedInternals.lData[parentCode] == 0) {
      taggedInternals.lData[parentCode] = 1;
      long k3 = 0;
      for (long k = 0; k < patternCount; k++) {
        for (long k2 = 0; k2 < alphabetDimension; k2++, k3++) {
          parentConditionals[k3] = 1.;
        }
      }
    }

    _CalcNode *currentTreeNode = isLeaf ? ((_CalcNode *)flatCLeaves(nodeCode))
                                        : ((_CalcNode *)flatTree(nodeCode));
    _Parameter *_hprestrict_ transitionMatrix =
        catAssignments ? nil : currentTreeNode->GetCompExp()->theData;

    // this will need to be toggled on a per site basis
    _Parameter *childVector;

    if (!isLeaf) {
      childVector = iNodeCache + (nodeCode * patternCount) * alphabetDimension;
    }

    for (long siteID = 0; siteID < patternCount;
         siteID++, parentConditionals += alphabetDimension) {
      if (catAssignments) {
        transitionMatrix = currentTreeNode
            ->GetCompExp(catAssignments[siteOrdering.lData[siteID]])->theData;
      }

      _Parameter _hprestrict_ *tMatrix = transitionMatrix;
      if (isLeaf) {
        long siteState =
            lNodeFlags[nodeCode * patternCount + siteOrdering.lData[siteID]];

        // a fully resolved leaf
        if (siteState >= 0) {
          tMatrix += siteState;
          for (long k = 0; k < alphabetDimension;
               k++, tMatrix += alphabetDimension) {
            parentConditionals[k] *= *tMatrix;
          }
          if (alsoDoLeaves) {
            for (long k = 0; k < alphabetDimension; k++) {
              leafBuffer[k] = siteState;
            }
            leafBuffer += alphabetDimension;
          }

          continue;
        } else {
          // an ambiguous leaf
          childVector =
              lNodeResolutions->theData + (-siteState - 1) * alphabetDimension;
        }

      }

      // now repopulate this vector as necessary -- if we are here this means
      // that the subtree below has been completely processed,
      // the i-th cell of childVector contains the likelihood of the _optimal_
      // assignment in the subtree below given that the character at the current
      // node is i.

      // hence, given parent state 'p', we optimize
      // max_i pr (p->i) childVector [i] and store it in the p cell of vector
      // childVector

      _Parameter overallMax = 0.0;

      long *stateBuffer = isLeaf ? leafBuffer : stateCache;

      // check for degeneracy

      long howManyOnes = 0;
      for (long k = 0; k < alphabetDimension; k++) {
        howManyOnes += childVector[k] == 1.;
      }

      if (howManyOnes == alphabetDimension) {
        for (long k = 0; k < alphabetDimension; k++) {
          stateBuffer[k] = -1;
        }
      } else {
        for (long p = 0; p < alphabetDimension; p++) {
          _Parameter max_lik = 0.;
          long max_idx = 0;

          for (long c = 0; c < alphabetDimension; c++) {
            _Parameter thisV = tMatrix[c] * childVector[c];
            if (thisV > max_lik) {
              max_lik = thisV;
              max_idx = c;
            }
          }

          stateBuffer[p] = max_idx;
          buffer[p] = max_lik;

          if (max_lik > overallMax) {
            overallMax = max_lik;
          }

          tMatrix += alphabetDimension;
        }

        if (overallMax > 0.0 && overallMax < _lfScalingFactorThreshold) {
          for (long k = 0; k < alphabetDimension; k++) {
            buffer[k] *= _lfScalerUpwards;
          }
        }

        // buffer[p] now contains the maximum likelihood of the tree
        // from this point forward given that parent state is p
        // and stateBuffer[p] stores the maximimizing assignment
        // for this node

        for (long k = 0; k < alphabetDimension; k++) {
          long stateValue = stateBuffer[k];
          if (stateValue >= 0) {
            parentConditionals[k] *= buffer[k];
          }
        }
      }

      if (isLeaf) {
        if (alsoDoLeaves) {
          leafBuffer += alphabetDimension;
        }
      } else {
        stateCache += alphabetDimension;
      }

      childVector += alphabetDimension;
    }
  }

  _List *result = new _List;
  for (long k = 0; k < stateCacheDim; k++) {
    result->AppendNewInstance(new _String(siteCount * unitLength, false));
  }

  _Parameter _hprestrict_ *rootConditionals =
      iNodeCache + alphabetDimension * ((iNodeCount - 1) * patternCount);
  _SimpleList parentStates(stateCacheDim, 0, 0), conversion;

  stateCache -= patternCount * (iNodeCount - 1) * alphabetDimension;
  if (alsoDoLeaves) {
    leafBuffer -= patternCount * leafCount * alphabetDimension;
  }

  _AVLListXL conversionAVL(&conversion);
  _String codeBuffer(unitLength, false);

  for (long siteID = 0; siteID < patternCount;
       siteID++, rootConditionals += alphabetDimension) {
    _Parameter max_lik = 0.;
    long max_idx = 0;

    long howManyOnes = 0;
    for (long k = 0; k < alphabetDimension; k++) {
      howManyOnes += rootConditionals[k] == 1.;
    }

    _SimpleList *patternMap =
        (_SimpleList *)expandedSiteMap(siteOrdering.lData[siteID]);

    if (howManyOnes != alphabetDimension) {
      for (long c = 0; c < alphabetDimension; c++) {
        _Parameter thisV = theProbs[c] * rootConditionals[c];
        if (thisV > max_lik) {
          max_lik = thisV;
          max_idx = c;
        }
      }

      parentStates.lData[iNodeCount - 1] = max_idx;
      for (long nodeID = iNodeCount - 2; nodeID >= 0; nodeID--) {
        long parentState =
            parentStates.lData[flatParents.lData[nodeID + flatLeaves.lLength]];
        if (parentState == -1) {
          parentStates.lData[nodeID] = -1;
        } else {
          parentStates.lData[nodeID] =
              stateCache[(patternCount * nodeID + siteID) * alphabetDimension +
                         parentState];
        }
      }
      if (alsoDoLeaves)
        for (long nodeID = 0; nodeID < leafCount; nodeID++) {
          long parentState = parentStates.lData[flatParents.lData[nodeID]];
          if (parentState == -1) {
            parentStates.lData[nodeID + iNodeCount] = -1;
          } else {
            parentStates.lData[nodeID + iNodeCount] =
                leafBuffer[(patternCount * nodeID + siteID) *
                               alphabetDimension + parentState];
          }
        }
    } else {
      parentStates.Populate(stateCacheDim, -1, 0);
    }

    for (long nodeID = 0; nodeID < stateCacheDim; nodeID++) {
      dsf->ConvertCodeToLettersBuffered(
          dsf->CorrectCode(parentStates.lData[nodeID]), unitLength,
          codeBuffer.sData, &conversionAVL);
      _String *sequence = (_String *)(*result)(
          nodeID < iNodeCount ? postToIn.lData[nodeID] : nodeID);

      for (long site = 0; site < patternMap->lLength; site++) {
        char *storeHere =
            sequence->sData + patternMap->lData[site] * unitLength;
        for (long charS = 0; charS < unitLength; charS++) {
          storeHere[charS] = codeBuffer.sData[charS];
        }
      }
    }
  }

  delete[] stateCache;
  delete[] leafBuffer;
  delete[] buffer;

  return result;
}

//______________________________________________________________________________
void _TheTree::SetupCategoryMapsForNodes(_List &containerVariables,
                                         _SimpleList &classCounter,
                                         _SimpleList &multipliers) {
  _CalcNode *iterator = DepthWiseTraversal(true);
  while (iterator) {
    iterator->SetupCategoryMap(containerVariables, classCounter, multipliers);
    iterator = DepthWiseTraversal(false);
  }
}

//______________________________________________________________________________
_Parameter _TheTree::Process3TaxonNumericFilter(_DataSetFilterNumeric *dsf,
                                                long catID) {

  _Parameter *l0 = dsf->probabilityVectors.theData +
                   dsf->categoryShifter * catID +
                   dsf->theNodeMap.lData[0] * dsf->shifter,
             *l1 = dsf->probabilityVectors.theData +
                   dsf->categoryShifter * catID +
                   dsf->theNodeMap.lData[1] * dsf->shifter,
             *l2 = dsf->probabilityVectors.theData +
                   dsf->categoryShifter * catID +
                   dsf->theNodeMap.lData[2] * dsf->shifter,
             *matrix0 =
                 ((_CalcNode *)(LocateVar(theRoot->nodes.data[0]->in_object)))
                     ->GetCompExp(catID)->theData,
             *matrix1 =
                 ((_CalcNode *)(LocateVar(theRoot->nodes.data[1]->in_object)))
                     ->GetCompExp(catID)->theData,
             *matrix2 =
                 ((_CalcNode *)(LocateVar(theRoot->nodes.data[2]->in_object)))
                     ->GetCompExp(catID)->theData,
             overallResult = 0.;

  long patternCount = dsf->NumberDistinctSites();

  _Parameter currentAccumulator = 1.;

  for (long patternIndex = 0; patternIndex < patternCount;
       patternIndex++, l0 += 4, l1 += 4, l2 += 4) {
    _Parameter rp0 = l0[0] * matrix0[0] + l0[1] * matrix0[1] +
                     l0[2] * matrix0[2] + l0[3] * matrix0[3];
    _Parameter rp1 = l0[0] * matrix0[4] + l0[1] * matrix0[5] +
                     l0[2] * matrix0[6] + l0[3] * matrix0[7];
    _Parameter rp2 = l0[0] * matrix0[8] + l0[1] * matrix0[9] +
                     l0[2] * matrix0[10] + l0[3] * matrix0[11];
    _Parameter rp3 = l0[0] * matrix0[12] + l0[1] * matrix0[13] +
                     l0[2] * matrix0[14] + l0[3] * matrix0[15];

    rp0 *= l1[0] * matrix1[0] + l1[1] * matrix1[1] + l1[2] * matrix1[2] +
           l1[3] * matrix1[3];
    rp1 *= l1[0] * matrix1[4] + l1[1] * matrix1[5] + l1[2] * matrix1[6] +
           l1[3] * matrix1[7];
    rp2 *= l1[0] * matrix1[8] + l1[1] * matrix1[9] + l1[2] * matrix1[10] +
           l1[3] * matrix1[11];
    rp3 *= l1[0] * matrix1[12] + l1[1] * matrix1[13] + l1[2] * matrix1[14] +
           l1[3] * matrix1[15];

    rp0 *= l2[0] * matrix2[0] + l2[1] * matrix2[1] + l2[2] * matrix2[2] +
           l2[3] * matrix2[3];
    rp1 *= l2[0] * matrix2[4] + l2[1] * matrix2[5] + l2[2] * matrix2[6] +
           l2[3] * matrix2[7];
    rp2 *= l2[0] * matrix2[8] + l2[1] * matrix2[9] + l2[2] * matrix2[10] +
           l2[3] * matrix2[11];
    rp3 *= l2[0] * matrix2[12] + l2[1] * matrix2[13] + l2[2] * matrix2[14] +
           l2[3] * matrix2[15];

    _Parameter result = theProbs[0] * rp0 + theProbs[1] * rp1 +
                        theProbs[2] * rp2 + theProbs[3] * rp3;

    if (result <= 0.0) {
      return -A_LARGE_NUMBER;
    }

    long patternFreq = dsf->theFrequencies[patternIndex];
    for (long freqIterator = 0; freqIterator < patternFreq; freqIterator++) {
      _Parameter tryMultiplication = currentAccumulator * result;
      if (tryMultiplication > 1.e-300) {
        currentAccumulator = tryMultiplication;
      } else {
        overallResult += myLog(currentAccumulator);
        currentAccumulator = result;
      }
    }
  }
  return overallResult + myLog(currentAccumulator);
}
