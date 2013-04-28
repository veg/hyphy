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
#include "treetopology.h"
#include "thetree.h"
#include "calcnode_globals.h"

#include "math.h"
#include "ctype.h"
#include "string.h"
#include "calcnode.h"
#include "scfg.h"
#include "parser.h"

#include "category.h"
#include "batchlan.h"
#include "likefunc.h"
#include "float.h"

long *nonZeroNodes = nil, nonZeroNodesDim = 0;

#ifdef __MP__
pthread_t *matrixThreads = nil;
ThreadMatrixTask *matrixTasks = nil;
pthread_mutex_t matrixMutex = PTHREAD_MUTEX_INITIALIZER;
#endif

char takeBranchLengths = 0, autoSolveBranchLengths = 0;

_Parameter ignoringInternalNames = 0.0;

_Parameter treeLayoutVert, scalingLogConstant = log(1.e300);

_SimpleList convertedMatrixExpressionsL;
_AVLListX convertedMatrixExpressions (&convertedMatrixExpressionsL);

_String expectedNumberOfSubs = "EXPECTED_NUMBER_OF_SUBSTITUTIONS",
        stringSuppliedLengths = "STRING_SUPPLIED_LENGTHS",
        noInternalLabels = "NO_INTERNAL_LABELS",
        includeModelSpecs = "INCLUDE_MODEL_SPECS",
        acceptRootedTrees = "ACCEPT_ROOTED_TREES",
        acceptBranchLengths = "ACCEPT_BRANCH_LENGTHS",
        autoConvertBL = "AUTOMATICALLY_CONVERT_BRANCH_LENGTHS",
        splitNodeNames = "SPLIT_NODE_NAMES",
        internalNodePrefix = "INTERNAL_NODE_PREFIX",
        ignoreUserINames = "IGNORE_INTERNAL_NODE_LABELS", cotNode = "COT_NODE",
        cotSplit = "COT_SPLIT", cotBranchLength = "COT_BRANCH_LENGTH",
        cotDistance = "COT_DISTANCE", cotCDF = "COT_CDF",
        cotSamples = "COT_SAMPLES", cotSampler = "COT_SAMPLER",
        cotToNode = "COT_TO_NODE", treeOutputAVL = "TREE_OUTPUT_OPTIONS",
        treeOutputBackground = "TREE_OUTPUT_BACKGROUND",
        treeOutputRightMargin = "TREE_OUTPUT_RIGHT_MARGIN",
        treeOutputEmbed = "TREE_OUTPUT_EMBED",
        treeOutputXtraMargin = "TREE_OUTPUT_XTRA_MARGIN",
        treeOutputSplit = "TREE_OUTPUT_BRANCH_SPLIT",
        treeOutputNotchesColor = "TREE_OUTPUT_BRANCH_NOTCHES_COLOR",
        treeOutputNotches = "TREE_OUTPUT_BRANCH_NOTCHES",
        treeOutputLabel = "TREE_OUTPUT_BRANCH_LABEL",
        treeOutputTLabel = "TREE_OUTPUT_BRANCH_TLABEL",
        treeOutputColor = "TREE_OUTPUT_BRANCH_COLOR",
        treeOutputThickness = "TREE_OUTPUT_BRANCH_THICKNESS",
        treeOutputLinecap = "TREE_OUTPUT_BRANCH_LINECAP",
        treeOutputDash = "TREE_OUTPUT_BRANCH_DASH",
        treeOutputOLabel = "TREE_OUTPUT_OVER_BRANCH",
        treeOutputSymbols = "TREE_OUTPUT_SYMBOLS",
        treeOutputSymbolSize = "TREE_OUTPUT_SYMBOL_SIZE",
        treeOutputExtraPS = "TREE_OUTPUT_EXTRA_POSTSCRIPT",
        treeOutputPrefixPS = "TREE_OUTPUT_PREFIX_POSTSCRIPT",
        treeOutputLayout = "TREE_OUTPUT_LAYOUT",
        treeOutputNNPlaceH = "__NODE_NAME__",
        treeOutputFSPlaceH = "__FONT_SIZE__",
        largeMatrixBranchLengthDimension =
            "LARGE_MATRIX_BRANCH_LENGTH_MODIFIER_DIMENSION",
        largeMatrixBranchLength = "LARGE_MATRIX_BRANCH_LENGTH_MODIFIER",
        newNodeGraftName = "NAME", newNodeGraftWhere = "WHERE",
        newNodeGraftParent = "PARENT", eqWithReroot = "Equal with reroot at ",
        eqWithoutReroot = "Equal without rerooting", iNodePrefix;

_Parameter _timesCharWidths[256] =
    { // Hardcoded relative widths of all 255 characters in the Times font, for
      // the use of PSTreeString
      0, 0.721569, 0.721569, 0.721569, 0.721569, 0.721569, 0.721569, 0.721569,
      0, 0.25098, 0.721569, 0.721569, 0.721569, 0, 0.721569, 0.721569, 0.721569,
      0.721569, 0.721569, 0.721569, 0.721569, 0.721569, 0.721569, 0.721569,
      0.721569, 0.721569, 0.721569, 0.721569, 0.721569, 0, 0.721569, 0.721569,
      0.25098, 0.333333, 0.407843, 0.501961, 0.501961, 0.831373, 0.776471,
      0.180392, 0.333333, 0.333333, 0.501961, 0.564706, 0.25098, 0.333333,
      0.25098, 0.278431, 0.501961, 0.501961, 0.501961, 0.501961, 0.501961,
      0.501961, 0.501961, 0.501961, 0.501961, 0.501961, 0.278431, 0.278431,
      0.564706, 0.564706, 0.564706, 0.443137, 0.921569, 0.721569, 0.666667,
      0.666667, 0.721569, 0.611765, 0.556863, 0.721569, 0.721569, 0.333333,
      0.388235, 0.721569, 0.611765, 0.890196, 0.721569, 0.721569, 0.556863,
      0.721569, 0.666667, 0.556863, 0.611765, 0.721569, 0.721569, 0.945098,
      0.721569, 0.721569, 0.611765, 0.333333, 0.278431, 0.333333, 0.470588,
      0.501961, 0.333333, 0.443137, 0.501961, 0.443137, 0.501961, 0.443137,
      0.333333, 0.501961, 0.501961, 0.278431, 0.278431, 0.501961, 0.278431,
      0.776471, 0.501961, 0.501961, 0.501961, 0.501961, 0.333333, 0.388235,
      0.278431, 0.501961, 0.501961, 0.721569, 0.501961, 0.501961, 0.443137,
      0.478431, 0.2, 0.478431, 0.541176, 0.721569, 0.721569, 0.721569, 0.666667,
      0.611765, 0.721569, 0.721569, 0.721569, 0.443137, 0.443137, 0.443137,
      0.443137, 0.443137, 0.443137, 0.443137, 0.443137, 0.443137, 0.443137,
      0.443137, 0.278431, 0.278431, 0.278431, 0.278431, 0.501961, 0.501961,
      0.501961, 0.501961, 0.501961, 0.501961, 0.501961, 0.501961, 0.501961,
      0.501961, 0.501961, 0.4, 0.501961, 0.501961, 0.501961, 0.34902, 0.454902,
      0.501961, 0.760784, 0.760784, 0.980392, 0.333333, 0.333333, 0.54902,
      0.890196, 0.721569, 0.713725, 0.54902, 0.54902, 0.54902, 0.501961,
      0.576471, 0.494118, 0.713725, 0.823529, 0.54902, 0.27451, 0.27451,
      0.309804, 0.768627, 0.666667, 0.501961, 0.443137, 0.333333, 0.564706,
      0.54902, 0.501961, 0.54902, 0.611765, 0.501961, 0.501961, 1, 0.25098,
      0.721569, 0.721569, 0.721569, 0.890196, 0.721569, 0.501961, 1, 0.443137,
      0.443137, 0.333333, 0.333333, 0.54902, 0.494118, 0.501961, 0.721569,
      0.168627, 0.745098, 0.333333, 0.333333, 0.556863, 0.556863, 0.501961,
      0.25098, 0.333333, 0.443137, 1, 0.721569, 0.611765, 0.721569, 0.611765,
      0.611765, 0.333333, 0.333333, 0.333333, 0.333333, 0.721569, 0.721569,
      0.788235, 0.721569, 0.721569, 0.721569, 0.721569, 0.278431, 0.333333,
      0.333333, 0.333333, 0.333333, 0.333333, 0.333333, 0.333333, 0.333333,
      0.333333, 0.333333 },
           _maxTimesCharWidth = 0.980392;

//______________________________________________________________________________
_TreeTopology::_TreeTopology() { rooted = UNROOTED; }

//______________________________________________________________________________
_TreeTopology::~_TreeTopology(void) {
  if (theRoot) {
    theRoot->delete_tree();
    delete (theRoot);
    theRoot = nil;
  }
  if (compExp) {
    DeleteObject(compExp);
    compExp = nil;
  }
}

//______________________________________________________________________________
void _TreeTopology::PreTreeConstructor(bool) {
  rooted = UNROOTED;
  compExp = (_Matrix *)checkPointer(new _GrowingVector);

  iNodePrefix = "Node";
  _PMathObj iv = FetchObjectFromVariableByType(&internalNodePrefix, STRING);
  if (iv) {
    iNodePrefix = *((_FString *)iv)->theString;
  }
  checkParameter(ignoreUserINames, ignoringInternalNames, 0.0);
}

//______________________________________________________________________________
node<long> *_TreeTopology::CopyTreeStructure(node<long> *theNode, bool) {
  long i;
  node<long> *locNode = new node<long>;
  for (i = 0; i < theNode->get_num_nodes(); i++) {
    locNode->add_node(*CopyTreeStructure(theNode->go_down(i + 1), false));
  }

  locNode->init(theNode->in_object);
  return locNode;
}

//______________________________________________________________________________
bool _TreeTopology::IsCurrentNodeATip(void) {
  return currentNode ? currentNode->get_num_nodes() == 0 : false;
}

//______________________________________________________________________________
bool _TreeTopology::IsCurrentNodeTheRoot(void) {
  return currentNode->parent == nil;
}

//______________________________________________________________________________
bool _TreeTopology::IsDegenerate(void) { return theRoot->get_num_nodes() == 1; }

//______________________________________________________________________________
void _TreeTopology::GetNodeName(node<long> *n, _String &r, bool fullName) {
  if (fullName) {
    r = *GetName() & '.' & *(_String *)(flatTree(n->in_object));
  } else {
    r = *(_String *)(flatTree(n->in_object));
  }
}

//______________________________________________________________________________
void _TreeTopology::SetLeafName(long res, _String *newName) {
  long count = 0;

  LeafWiseT(true);
  while (currentNode) {
    if (res == count) {
      flatTree.Replace(currentNode->in_object, newName, true);
      break;
    }
    count++;
    LeafWiseT(false);
  }
}

//______________________________________________________________________________
BaseRef _TreeTopology::toStr(void) {
  _String *res = new _String((unsigned long) 128, true), num;

  _Parameter skipILabels, includeMSP;

  checkParameter(noInternalLabels, skipILabels, 0.0);
  checkParameter(includeModelSpecs, includeMSP, 0.0);

  if (IsDegenerate()) {

    DepthWiseT(true);

    (*res) << '(';

    GetNodeName(theRoot, num);

    (*res) << &num;
    if (includeMSP > 0.5) {
      _String *mSpec = (_String *)flatCLeaves(theRoot->in_object);
      if (mSpec->sLength) {
        (*res) << '{';
        (*res) << mSpec;
        (*res) << '}';
      }
    }
    (*res) << ',';
    GetNodeName(currentNode, num);
    (*res) << &num;
    if (includeMSP > 0.5) {
      _String *mSpec = (_String *)flatCLeaves(currentNode->in_object);
      if (mSpec->sLength) {
        (*res) << '{';
        (*res) << mSpec;
        (*res) << '}';
      }
    }
    (*res) << ')';
  } else {

    long level = 0, myLevel = 0, lastLevel = 0, j;

    DepthWiseTLevel(myLevel, true);

    node<long> *curNode = currentNode, *nextNode;

    level = myLevel;

    bool isCTip = IsCurrentNodeATip(), isCTip2;

    DepthWiseTLevel(myLevel);
    nextNode = currentNode;

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
        GetNodeName(curNode, num);
        (*res) << &num;
      }

      if (includeMSP > 0.5) {
        _String *mSpec = (_String *)flatCLeaves(curNode->in_object);
        if (mSpec->sLength) {
          (*res) << '{';
          (*res) << mSpec;
          (*res) << '}';
        }
      }

      lastLevel = level;
      level = myLevel;
      curNode = nextNode;
      isCTip = isCTip2;
      DepthWiseTLevel(myLevel);
      nextNode = currentNode;
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
void _TreeTopology::toFileStr(FILE *f) {
  _String *s = (_String *)toStr();
  fprintf(f, "%s", s->sData);
  DeleteObject(s);
}

//______________________________________________________________________________
// execute this operation with the second arg if necessary
_PMathObj _TreeTopology::Execute(long opCode, _PMathObj p, _PMathObj p2, _hyExecutionContext *context) {

  switch (opCode) {
  case HY_OP_CODE_IDIV: { // Split ($) - 2nd argument
    if (p->ObjectClass() != NUMBER) {
      context->ReportError(
          "Invalid (not a number) 2nd argument is call to $ for trees.");
      return new _MathObject;
    }
    _Constant *cc = (_Constant *)TipCount();
    long size = cc->Value() / p->Value();

    if ((size <= 4) || (size > cc->Value() / 2)) {
      context->ReportError("Poor choice of the 2nd numeric agrument in to $ "
                           "for tree. Either the resulting cluster size is too "
                           "big(>half of the tree), or too small (<4)!");
      return new _MathObject;
    }

    long checkSize = 1, tol = 0;

    while (tol < size - 2) {
      _List *resL = SplitTreeIntoClusters(size, tol);

      checkSize = cc->Value();

      if (resL->lLength) {
        _Matrix *mRes = new _Matrix(resL->lLength, 2, false, true);
        checkPointer(mRes);

        for (long k = 0; k < resL->lLength; k++) {
          _List *thisList = (_List *)(*resL)(k);
          long nL = ((_Constant *)(*thisList)(1))->Value();
          mRes->Store(k, 0, nL);
          mRes->Store(k, 1, thisList->lLength - 2);
          checkSize -= nL;
        }

        if (checkSize == 0) {
          DeleteObject(cc);
          _Matrix selMatrix(1, resL->lLength, false, true);
          _List sortedList;
          for (long k = 0; k < resL->lLength; k++) {
            _List *thisList = (_List *)(*resL)(k);
            sortedList << (_String *)(*thisList)(0);
            _FString *choiceString = new _FString(*(_String *)(*thisList)(0));
            _Formula sf(choiceString);
            selMatrix.MStore(0, k, sf);
          }
          sortedList.Sort();
          for (long m = 0; m < sortedList.lLength; m++) {
            _FString *choiceString = new _FString(*(_String *)sortedList(m));
            _Formula sf(choiceString);
            selMatrix.MStore(0, m, sf);
          }

          CheckReceptacle(&splitNodeNames, empty, false)->SetValue(&selMatrix);
          DeleteObject(resL);
          return mRes;
        }

        DeleteObject(mRes);
      }

      DeleteObject(resL);
      tol++;
    }

    DeleteObject(cc);
    return new _Matrix(1, 1, false, true);
  } break;

  case HY_OP_CODE_MUL: // compute the strict consensus between T1 and T2
    if (p)
      return SplitsIdentity(p);
    break;

  case HY_OP_CODE_ADD: // +
    if (!p) {
      return Sum();
    }
    AddANode(p);
    return new _Constant(0.0);
    break;

  case HY_OP_CODE_LEQ: { // MatchPattern (<=)
    if ((p->ObjectClass() != TREE) && (p->ObjectClass() != TOPOLOGY)) {
      context->ReportError("Invalid (not a tree/topology) 2nd argument is call "
                           "to <= for trees/topologies.");
      return new _MathObject;
    }
    _String res(((_TreeTopology *)p)->MatchTreePattern(this));
    return new _Constant(!res.beginswith("Unequal"));
    break;
  }
  case HY_OP_CODE_EQ: // ==
    return Compare(p);
    break;
  case HY_OP_CODE_ABS: // Abs
    return FlatRepresentation();
    break;
  case HY_OP_CODE_BRANCHCOUNT: //BranchCount
    return BranchCount();
    break;
  case HY_OP_CODE_BRANCHLENGTH: //BranchLength
    return BranchLength(p);
    break;
  case HY_OP_CODE_BRANCHNAME: //BranchName
    return BranchName(p);
    break;
  case HY_OP_CODE_FORMAT: { // Format
    currentNode = theRoot;
    _String *tStr = new _String((unsigned long) 1024, true);
    SubTreeString(*tStr, p->Compute()->Value() > 0.1,
                  p2->Compute()->Value() > 0.1 ? -3 : -1, nil);
    tStr->Finalize();
    return new _FString(tStr);
  }
  case HY_OP_CODE_MACCESS: // MAccess
    return BranchName(p, true, p2);
    break;
  case HY_OP_CODE_MIN: // COT (Min)
    return FindCOT(p);
    break;
  case HY_OP_CODE_REROOTTREE: // RerootTree
    return RerootTree(p);
    break;
  case HY_OP_CODE_TEXTREESTRING: // TEXTreeString
                                 //return TEXTreeString(p);
    break;
  case HY_OP_CODE_TIPCOUNT: // TipCount
    return TipCount();
    break;
  case HY_OP_CODE_TIPNAME: // TipName
    return TipName(p);
    break;
  case HY_OP_CODE_TYPE: // Type
    return Type();
    break;
  case HY_OP_CODE_POWER: //^
    if (p)
      return AVLRepresentation(p);
  }

  WarnNotDefined(this, opCode, context);
  return nil;

}

//______________________________________________________________________________
void _TreeTopology::FindCOTHelper(node<long> *aNode, long parentIndex,
                                  _Matrix &distances, _Matrix &rootDistances,
                                  _Matrix &branchLengths, _List &childLists,
                                  _AVLListX &addressToIndexMap2, _Parameter d) {

  long myIndex =
           addressToIndexMap2.GetXtra(addressToIndexMap2.Find((BaseRef) aNode)),
       leafCount = distances.GetVDim();

  //bool        isRoot      = parentIndex < 0;

  _SimpleList *childLeaves = (_SimpleList *)childLists(myIndex);

  _Matrix *lookup = parentIndex >= 0 ? &distances : &rootDistances;

  if (parentIndex < 0) {
    parentIndex = 0;
  }

  long ci2 = 0;

  _Parameter myLength = branchLengths.theData[myIndex];

  for (long ci = 0; ci < leafCount; ci++) {
    if (ci == childLeaves->lData[ci2]) {
      if (ci2 < childLeaves->lLength - 1) {
        ci2++;
      }
    } else {
      distances.Store(myIndex, ci, (*lookup)(parentIndex, ci) + myLength);
    }
  }

  for (long ci3 = aNode->get_num_nodes(); ci3; ci3--) {
    FindCOTHelper(aNode->go_down(ci3), myIndex, distances, rootDistances,
                  branchLengths, childLists, addressToIndexMap2, myLength);
  }

}

//______________________________________________________________________________
void _TreeTopology::FindCOTHelper2(node<long> *aNode, _Matrix &branchSpans,
                                   _Matrix &branchLengths,
                                   _AVLListX &addressToIndexMap2,
                                   node<long> *referrer, _Parameter d) {
  long myIndex =
      aNode->parent
          ? addressToIndexMap2.GetXtra(addressToIndexMap2.Find((BaseRef) aNode))
          : -1;
  _Parameter myLength = myIndex >= 0 ? branchLengths.theData[myIndex] : 0.0;

  for (long ci = aNode->get_num_nodes(); ci; ci--) {
    node<long> *daChild = aNode->go_down(ci);
    if (daChild != referrer) {
      FindCOTHelper2(daChild, branchSpans, branchLengths, addressToIndexMap2,
                     nil, d + myLength);
    }
  }
  if (aNode->parent) { // not the root
    if (d >= 0.0) {
      branchSpans.Store(myIndex, 0, d);
    } else {
      branchSpans.Store(myIndex, 0, 0.);
    }

    branchSpans.Store(myIndex, 1, d + myLength);
    if (referrer) {
      FindCOTHelper2(aNode->parent, branchSpans, branchLengths,
                     addressToIndexMap2, aNode, d + myLength);
    }
  }
}

//______________________________________________________________________________
_AssociativeList *_TreeTopology::FindCOT(_PMathObj p)
    // Find the Center of the Tree (COT) location
    // using an L_p metric (L_2 works well)
    {
  _Parameter power = p->Compute()->Value(), totalTreeLength = 0.0;

  _AssociativeList *resList = new _AssociativeList;

  if (power <= 0.) {
    WarnError(_String("Invalid power argument in call to COT finder (Min on "
                      "trees). Must be positive, had :") & power);
    return resList;
  }

  _SimpleList avlSL, avlSL2, listOfNodes;

  _List avlSL3, childLists;

  _AVLListX addressToIndexMap(&avlSL), // leaves only
      addressToIndexMap2(&avlSL2),     // complete index
      lengthToIndexMap(&avlSL3);

  long leafCount = 0, branchCount = 0, tIndex = 0;

  DepthWiseT(true);

  while (currentNode->parent) {
    if (IsCurrentNodeATip()) {
      addressToIndexMap.Insert((BaseRef) currentNode, leafCount++);
    } else {
      branchCount++;
    }

    addressToIndexMap2.Insert((BaseRef) currentNode, tIndex++);
    DepthWiseT(false);
  }

  // allocate the matrix of path lengths with hardwired (traversal order)
  // indices
  // also allocate a list of sorted lists to store children nodes
  // and a map of (longed) node addresses to post order traversal indices

  _Matrix distances(branchCount + leafCount, leafCount, false, true),
      rootDistances(1, leafCount, false, true),
      branchLengths(1, branchCount + leafCount, false, true),
      branchSpans(branchCount + leafCount + 1, 2, false, true);

  _String nodeName;

  // pass 1: fill up the nodes up to the root (i.e. below any internal node)

  DepthWiseT(true);
  tIndex = 0;

  while (currentNode->parent) {

    _Parameter myLength;
    GetBranchLength(currentNode, myLength);
    lengthToIndexMap.Insert(new _String(totalTreeLength), tIndex, false, true);
    totalTreeLength += myLength;

    branchLengths.Store(0, tIndex++, myLength);
    listOfNodes << (long) currentNode;
    _SimpleList *childIndices = (_SimpleList *)checkPointer(new _SimpleList);

    if (IsCurrentNodeATip()) {
      (*childIndices) << addressToIndexMap.GetXtra(
                             addressToIndexMap.Find((BaseRef) currentNode));
      //leafNames<<LocateVar(currentNode->in_object)->GetName();
    } else {
      long myIndex = addressToIndexMap2.GetXtra(
          addressToIndexMap2.Find((BaseRef) currentNode));

      _SimpleList mappedLeaves(leafCount, 0, 0);
      for (long ci = currentNode->get_num_nodes(); ci; ci--) {
        long childIndex = addressToIndexMap2.GetXtra(
            addressToIndexMap2.Find((BaseRef) currentNode->go_down(ci)));
        _SimpleList *childLeaves = (_SimpleList *)childLists(childIndex);

        myLength = branchLengths.theData[childIndex];

        for (long ci2 = 0; ci2 < childLeaves->lLength; ci2++) {
          long ttIndex = childLeaves->lData[ci2];
          mappedLeaves.lData[ttIndex] = 1;
          distances.Store(myIndex, ttIndex,
                          distances(childIndex, ttIndex) + myLength);
        }
      }

      for (long ci2 = 0; ci2 < leafCount; ci2++)
        if (mappedLeaves.lData[ci2]) {
          (*childIndices) << ci2;
        }
    }

    childLists.AppendNewInstance(childIndices);
    DepthWiseT(false);
  }

  // pass 2: fill the root vector

  //nodeName = "COT_DM1";
  //setParameter (nodeName, &distances);

  for (long ci = theRoot->get_num_nodes(); ci; ci--) {
    long childIndex = addressToIndexMap2.GetXtra(
        addressToIndexMap2.Find((BaseRef) theRoot->go_down(ci)));
    _SimpleList *childLeaves = (_SimpleList *)childLists(childIndex);
    _Parameter myLength = branchLengths.theData[childIndex];
    for (long ci2 = 0; ci2 < childLeaves->lLength; ci2++) {
      tIndex = childLeaves->lData[ci2];
      rootDistances.Store(0, tIndex, distances(childIndex, tIndex) + myLength);
      //printf ("root->%s = %g\n", ((_String*)leafNames(tIndex))->sData,
      //distances (childIndex, tIndex) + myLength);
    }
  }

  // pass 3: fill in the "other site" branch lengths

  for (long ci3 = theRoot->get_num_nodes(); ci3; ci3--) {
    FindCOTHelper(theRoot->go_down(ci3), -1, distances, rootDistances,
                  branchLengths, childLists, addressToIndexMap2, 0);
  }

  //nodeName = "COT_DM2";
  //setParameter (nodeName, &distances);

  // now traverse the tree and look for the maximum

  tIndex = 0; // stores the index of current min

  _Parameter currentMin = 1e100, currentBranchSplit = 0;

  for (long ci = distances.GetHDim() - 1; ci >= 0; ci--) {
    _Parameter T = branchLengths.theData[ci];
    _SimpleList *childLeaves = (_SimpleList *)childLists(ci);
    long ci2 = 0;

    if (CheckEqual(power, 2.0)) {
      _Parameter sumbT = 0., sumbT2 = 0., suma = 0., suma2 = 0.;

      for (long ci3 = 0; ci3 < leafCount; ci3++) {
        _Parameter tt = distances(ci, ci3);

        /*
        printf ("%s->%s = %g\n", ttt2.sData, ((_String*)leafNames(ci3))->sData,
        tt);
        */

        if (ci3 == childLeaves->lData[ci2]) {
          if (ci2 < childLeaves->lLength - 1) {
            ci2++;
          }

          suma += tt;
          suma2 += tt * tt;
        } else {
          sumbT += tt;
          sumbT2 += tt * tt;
        }
      }

      _Parameter tt = (sumbT - suma) / leafCount; /*(sumbT-suma)/leafCount*/
      ;
      if (tt < 0.0) {
        tt = 0.;
      } else if (tt > T) {
        tt = T;
      }

      sumbT = tt * tt * leafCount + 2 * tt * (suma - sumbT) + suma2 + sumbT2;

      if (sumbT < currentMin) {
        tIndex = ci;
        currentBranchSplit = tt;
        currentMin = sumbT;
      }
    } else {
      _Parameter step = T > 0.0 ? T * 0.0001 : 0.1, currentT = 0.;

      while (currentT < T) {
        _Parameter dTT = 0.0;

        ci2 = 0;

        for (long ci3 = 0; ci3 < leafCount; ci3++) {
          _Parameter tt = distances(ci, ci3);
          if (ci3 == childLeaves->lData[ci2]) {
            if (ci2 < childLeaves->lLength - 1) {
              ci2++;
            }

            dTT += pow(tt + currentT, power);
          } else {
            dTT += pow(tt - currentT, power);
          }
        }

        if (dTT < currentMin) {
          tIndex = ci;
          currentBranchSplit = currentT;
          currentMin = dTT;
        }
        currentT += step;
      }
    }
  }

  node<long> *cotBranch = (node<long> *)listOfNodes.lData[tIndex];
  GetNodeName(cotBranch, nodeName);

  resList->MStore(cotNode, new _FString(nodeName, false), false);
  resList->MStore(cotSplit, new _Constant(currentBranchSplit), false);
  resList->MStore(cotDistance, new _Constant(currentMin), false);
  resList->MStore(cotBranchLength, new _Constant(branchLengths.theData[tIndex]),
                  false);

  //  compute the distribution of lengths away from COT

  FindCOTHelper2(cotBranch, branchSpans, branchLengths, addressToIndexMap2, nil,
                 currentBranchSplit - branchLengths.theData[tIndex]);
  if (cotBranch->parent) {
    _Parameter adjuster = branchLengths.theData[tIndex] - currentBranchSplit;
    branchSpans.Store(branchCount + leafCount, 1, adjuster);
    node<long> *cotParent = cotBranch->parent;
    if (cotParent->parent) {
      adjuster -= branchLengths.theData[addressToIndexMap2.GetXtra(
          addressToIndexMap2.Find((BaseRef) cotParent))];
    }
    FindCOTHelper2(cotParent, branchSpans, branchLengths, addressToIndexMap2,
                   cotBranch, adjuster);
  }

  _List timeSplits;
  _AVLListX timeSplitsAVL(&timeSplits);

  for (long sc = 0; sc <= branchCount + leafCount; sc++) {
    char buffer[256];
    snprintf(buffer, 256, "%.15f", branchSpans(sc, 0));
    nodeName = buffer;
    branchSpans.Store(sc, 0, nodeName.toNum());
    timeSplitsAVL.Insert(nodeName.makeDynamic(), 0, false, true);
    snprintf(buffer, 256, "%.15f", branchSpans(sc, 1));
    nodeName = buffer;
    branchSpans.Store(sc, 1, nodeName.toNum());
    timeSplitsAVL.Insert(nodeName.makeDynamic(), 0, false, true);
  }

  _Matrix cotCDFPoints(timeSplitsAVL.countitems(), 3, false, true);

  _AssociativeList *ctl =
      (_AssociativeList *)checkPointer(new _AssociativeList());

  _SimpleList tcache;

  long iv, k = timeSplitsAVL.Traverser(tcache, iv, timeSplitsAVL.GetRoot());

  for (long pc = 0; k >= 0; k = timeSplitsAVL.Traverser(tcache, iv), pc++) {
    timeSplitsAVL.SetXtra(k, pc);
    cotCDFPoints.Store(
        pc, 0, ((_String *)(*((_List *)timeSplitsAVL.dataList))(k))->toNum());
  }

  for (long mxc = 0; mxc <= branchCount + leafCount; mxc++) {
    _Parameter T0 = branchSpans(mxc, 0), T1 = branchSpans(mxc, 1);

    if (mxc < branchCount + leafCount) {
      GetNodeName((node<long> *)listOfNodes(mxc), nodeName);
      ctl->MStore(nodeName, new _Constant(T1), false);
    }

    char buffer[256];
    snprintf(buffer, 256, "%.15f", T0);
    nodeName = buffer;
    tcache.Clear();
    long startingPos = timeSplitsAVL.Find(&nodeName, tcache),
         k = timeSplitsAVL.Next(startingPos, tcache);

    for (long pc = timeSplitsAVL.GetXtra(k); k >= 0;
         k = timeSplitsAVL.Next(k, tcache), pc++) {
      _Parameter ub =
          ((_String *)(*((_List *)timeSplitsAVL.dataList))(k))->toNum();
      if (ub < T1 || CheckEqual(T1, ub)) {
        cotCDFPoints.Store(pc, 1, cotCDFPoints(pc, 1) + 1);
      } else {
        break;
      }
    }
  }

  for (long mxc = 1; mxc < cotCDFPoints.GetHDim(); mxc++) {
    cotCDFPoints.Store(mxc, 1,
                       cotCDFPoints(mxc, 1) *
                           (cotCDFPoints(mxc, 0) - cotCDFPoints(mxc - 1, 0)) /
                           totalTreeLength);
    cotCDFPoints.Store(mxc, 2, cotCDFPoints(mxc, 1) + cotCDFPoints(mxc - 1, 2));
  }

  timeSplitsAVL.Clear(true);

  resList->MStore(cotCDF, &cotCDFPoints, true);
  resList->MStore(cotToNode, ctl, false);

  //  sample  random branch placement
  if (totalTreeLength > 0.0) {
    _Parameter sampler;
    checkParameter(cotSamples, sampler, 0.0);

    tIndex = sampler;
    if (tIndex >= 1) {
      _Matrix sampledDs(tIndex, 1, false, true);

      for (long its = 0; its < tIndex; its++) {
        _Parameter tSample = (genrand_real2() * totalTreeLength);
        nodeName = tSample;
        long branchIndex = 0;
        if (lengthToIndexMap.FindBest(&nodeName, branchIndex) <= 0) {
          branchIndex--;
        }

        _Parameter T = branchLengths.theData[branchIndex];
        _SimpleList *childLeaves = (_SimpleList *)childLists(branchIndex);

        _Parameter dTT = 0.0;
        tSample -=
            ((_String *)(*(_List *)lengthToIndexMap.dataList)(branchIndex))
                ->toNum();

        long ci2 = 0;
        for (long ci3 = 0; ci3 < leafCount; ci3++) {
          _Parameter tt = distances(branchIndex, ci3);
          if (ci3 == childLeaves->lData[ci2]) {
            if (ci2 < childLeaves->lLength - 1) {
              ci2++;
            }

            dTT += pow(tt + tSample, power);
          } else {
            dTT += pow(tt + T - tSample, power);
          }
        }

        sampledDs.Store(its, 0, dTT);
      }
      setParameter(cotSampler, &sampledDs);
    }

  }
  return resList;
}

//______________________________________________________________________________
// compare tree topologies
_FString *_TreeTopology::Compare(_PMathObj p) {
  _FString *res = new _FString;

  long objClass = p->ObjectClass();

  if (objClass == TREE || objClass == TOPOLOGY) {
    _String cmp = CompareTrees((_TreeTopology *)p);
    if (cmp.startswith(eqWithReroot)) {
      (*res->theString) = cmp.Cut(
          eqWithReroot.sLength + ((_TreeTopology *)p)->GetName()->sLength + 1,
          cmp.sLength - 2);
    } else if (cmp.startswith(eqWithoutReroot)) {
      (*res->theString) = _String(' ');
    }
  }
  return res;
}

//______________________________________________________________________________
// compares whether the nodes create the same subpartition
char _TreeTopology::internalNodeCompare(node<long> *n1, node<long> *n2,
                                        _SimpleList &subTreeMap,
                                        _SimpleList *reindexer, bool cangoup,
                                        long totalSize, node<long> *n22,
                                        _TreeTopology *tree2, bool isPattern) {

  // count the number of children
  long nc1 = n1->get_num_nodes(),
       nc2 = n2->get_num_nodes() + (cangoup && n2->parent) -
             (n22 && (!n2->parent));

  if ((nc1 == nc2) || (isPattern && (nc2 < nc1))) {

    // now see if the descendants in all directions can be matched
    // first prepare the list of subtrees in the 2nd tree
    _List nodeMap;

    _SimpleList *complement = nil;
    _SimpleList stSizes;

    if ((cangoup || n22) && n2->parent) {
      complement = new _SimpleList((unsigned long) totalSize);
      checkPointer(complement);
      complement->lLength = totalSize;

      for (long k3 = 0; k3 < totalSize; k3++) {
        complement->lData[k3] = 1;
      }
    }

    nc2 = n2->get_num_nodes();

    long skippedChild = -1, complementCorrection = 0;

    for (long k = 1; k <= nc2; k++) {
      node<long> *meChild = n2->go_down(k);
      _SimpleList *childLeaves = (_SimpleList *)meChild->in_object;
      if (meChild == n22) {
        if (complement)
          if (reindexer)
            for (long k2 = 0; k2 < childLeaves->lLength; k2++) {
              long lidx = reindexer->lData[childLeaves->lData[k2]];
              if (lidx >= 0) {
                complement->lData[lidx] = 0;
              } else {
                complementCorrection++;
              }
            }
          else
            for (long k2 = 0; k2 < childLeaves->lLength; k2++) {
              complement->lData[childLeaves->lData[k2]] = 0;
            }

        skippedChild = k - 1;
      } else {
        _SimpleList *subTreeLeaves = new _SimpleList((unsigned long) totalSize);
        checkPointer(subTreeLeaves);
        subTreeLeaves->lLength = totalSize;

        //subTreeMap << ((_SimpleList*)n1->go_down(k)->in_object)->lLength;

        if (complement) {
          if (reindexer) {
            for (long k2 = 0; k2 < childLeaves->lLength; k2++) {
              long lidx = reindexer->lData[childLeaves->lData[k2]];
              if (lidx < 0) {
                delete complement;
                delete subTreeLeaves;
                return 0;
              }
              subTreeLeaves->lData[lidx] = 1;
              complement->lData[lidx] = 0;
            }
          } else {
            for (long k2 = 0; k2 < childLeaves->lLength; k2++) {
              subTreeLeaves->lData[childLeaves->lData[k2]] = 1;
              complement->lData[childLeaves->lData[k2]] = 0;
            }
          }
        } else if (reindexer) {
          for (long k2 = 0; k2 < childLeaves->lLength; k2++) {
            long lidx = reindexer->lData[childLeaves->lData[k2]];
            if (lidx < 0) {
              delete subTreeLeaves;
              return 0;
            }
            subTreeLeaves->lData[lidx] = 1;
          }
        } else {
          for (long k2 = 0; k2 < childLeaves->lLength; k2++) {
            subTreeLeaves->lData[childLeaves->lData[k2]] = 1;
          }
        }
        stSizes << childLeaves->lLength;
        nodeMap << subTreeLeaves;
        DeleteObject(subTreeLeaves);
      }
    }

    if (complement) {
      nodeMap << complement;
      DeleteObject(complement);
      stSizes << totalSize;

      for (long k4 = 0; k4 < stSizes.lLength - 1; k4++) {
        stSizes.lData[stSizes.lLength - 1] -= stSizes.lData[k4];
      }

      if (n22) {
        stSizes.lData[stSizes.lLength - 1] -=
            ((_SimpleList *)n22->in_object)->lLength - complementCorrection;
      }
    }

    // now go through the subtrees of the first tree and match up (if we can)
    // with the 2nd 1
    _SimpleList unmatchedPatterns;

    for (long k5 = 1; k5 <= nc1; k5++) {
      _SimpleList *childNodes = (_SimpleList *)n1->go_down(k5)->in_object;
      long k6;
      for (k6 = 0; k6 < stSizes.lLength; k6++)
        // potential subtree match
        if (stSizes.lData[k6] == childNodes->lLength) {
          _SimpleList *potMap = (_SimpleList *)nodeMap(k6);
          long k7;
          for (k7 = 0; k7 < childNodes->lLength; k7++)
            if (potMap->lData[childNodes->lData[k7]] == 0) {
              break;
            }

          if (k7 == childNodes->lLength) {
            stSizes.lData[k6] = -1;

            if (complement && (k6 == stSizes.lLength - 1)) {
              subTreeMap << -1;
            } else {
              if (n22 && (k6 >= skippedChild)) {
                subTreeMap << k6 + 1;
              } else {
                subTreeMap << k6;
              }
            }
            break;
          }
        }

      if (k6 == stSizes.lLength)
        if (isPattern) {
          unmatchedPatterns << k5;
          subTreeMap << -2;
        } else {
          return 0;
        }
    }

    if (isPattern && unmatchedPatterns.lLength) {
      _SimpleList rematchedPatterns;
      for (long k7 = 0; k7 < unmatchedPatterns.lLength; k7++) {
        rematchedPatterns << -1;
      }

      for (long k8 = 0; k8 < stSizes.lLength; k8++) {
        if (stSizes.lData[k8] > 0) {
          _SimpleList *potMap = (_SimpleList *)nodeMap(k8);
          for (long k9 = 0; k9 < unmatchedPatterns.lLength; k9++) {
            if (rematchedPatterns.lData[k9] < 0) {
              _SimpleList *childNodes = (_SimpleList *)n1
                  ->go_down(unmatchedPatterns.lData[k9])->in_object;
              long k10 = 0;
              for (; k10 < childNodes->lLength; k10++)
                if (potMap->lData[childNodes->lData[k10]] == 0) {
                  break;
                }
              if (k10 == childNodes->lLength) {
                rematchedPatterns.lData[k9] = k8;
                if (complement && (k8 == stSizes.lLength - 1)) {
                  subTreeMap.lData[unmatchedPatterns.lData[k9] - 1] =
                      -(0xfffffff);
                } else {
                  if (n22 && (k8 >= skippedChild)) {
                    subTreeMap.lData[unmatchedPatterns.lData[k9] - 1] = -k8 - 3;
                  } else {
                    subTreeMap.lData[unmatchedPatterns.lData[k9] - 1] = -k8 - 2;
                  }
                }
              }
            }
          }
        } else {
          stSizes.lData[k8] = 0;
        }
      }

      if (rematchedPatterns.Find(-1) >= 0) {
        return 0;
      }

      for (long k11 = 0; k11 < rematchedPatterns.lLength; k11++) {
        stSizes.lData[rematchedPatterns.lData[k11]] -=
            ((_SimpleList *)n1->go_down(unmatchedPatterns.lData[k11])
                 ->in_object)->lLength;
      }

      for (long k12 = 0; k12 < stSizes.lLength; k12++)
        if (stSizes.lData[k12]) {
          return 0;
        }
    }

    /*_String nm1, nm2;
    tree2->GetNodeName (n1,nm1);
    GetNodeName (n2,nm2);
    printf ("Node match %s, %s\n",  nm1.sData, nm2.sData, "\n");*/

    return 1;
  }

  return 0;
}

//long   itcCount = 0;
//______________________________________________________________________________
// compare tree topologies
char _TreeTopology::internalTreeCompare(node<long> *n1, node<long> *n2,
                                        _SimpleList *reindexer, char compMode,
                                        long totalSize, node<long> *n22,
                                        _TreeTopology *tree2, bool isPattern) {
  if (n1->get_num_nodes() == 0) {
    return 1;
  } else {
    _SimpleList mapper;
    if (internalNodeCompare(n1, n2, mapper, reindexer, compMode, totalSize, n22,
                            tree2, isPattern)) {
      long nc1 = n1->get_num_nodes();

      _SimpleList furtherMatchedPatterns;
      _List patternList;

      for (long k1 = 0; k1 < nc1; k1++) {
        if (mapper.lData[k1] >= 0) {
          if (internalTreeCompare(n1->go_down(k1 + 1),
                                  n2->go_down(mapper.lData[k1] + 1), reindexer,
                                  0, totalSize, nil, tree2, isPattern) < 1) {
            return -1;
          }
        } else if (mapper.lData[k1] == -1) {
          if (internalTreeCompare(n1->go_down(k1 + 1), n2->parent, reindexer, 0,
                                  totalSize, n2, tree2, isPattern) < 1) {
            return -1;
          }
        } else {
          long idx = -mapper.lData[k1] - 2, k;

          _SimpleList *patched;
          if ((k = furtherMatchedPatterns.Find(idx)) < 0) {
            k = furtherMatchedPatterns.lLength;
            furtherMatchedPatterns << idx;
            patched = new _SimpleList;
            checkPointer(patched);
            patternList << patched;
            DeleteObject(patched);
          }

          patched = (_SimpleList *)patternList(k);
          (*patched) << k1;
        }
      }
      for (long k2 = 0; k2 < furtherMatchedPatterns.lLength; k2++) {
        node<long> *dummy = new node<long>;
        checkPointer(dummy);
        dummy->parent = n1->parent;
        _SimpleList *children = (_SimpleList *)patternList(k2),
                    *newLeaves = new _SimpleList;

        checkPointer(newLeaves);
        dummy->in_object = (long) newLeaves;

        for (long k3 = 0; k3 < children->lLength; k3++) {
          node<long> *aChild = n1->go_down(children->lData[k3] + 1);
          dummy->nodes.add(aChild);
          _SimpleList t;
          t.Union(*newLeaves, *(_SimpleList *)aChild->in_object);
          newLeaves->Clear();
          newLeaves->Duplicate(&t);
        }

        long k4 = furtherMatchedPatterns.lData[k2];
        char res = 1;

        if (newLeaves->lLength > 1) {
          if (k4 < n2->get_num_nodes()) {
            res = internalTreeCompare(dummy, n2->go_down(k4 + 1), reindexer, 0,
                                      totalSize, nil, tree2, isPattern);
          } else {
            res = internalTreeCompare(dummy, n2->parent, reindexer, 0,
                                      totalSize, n2, tree2, isPattern);
          }
        }

        DeleteObject(newLeaves);
        delete (dummy);

        if (res < 1) {
          return -1;
        }

      }
    } else {
      return 0;
    }
  }
  return 1;
}

//______________________________________________________________________________
void _TreeTopology::EdgeCount(long &leaves, long &internals) {
  leaves = 0;
  internals = 0;
  DepthWiseT(true);
  while (currentNode) {
    if (IsCurrentNodeATip()) {
      leaves++;
    } else {
      internals++;
    }

    DepthWiseT(false);
  }

}

//______________________________________________________________________________
_PMathObj _TreeTopology::TipCount(void) {
  long leaves, ints;
  EdgeCount(leaves, ints);
  return new _Constant(leaves);
}

//______________________________________________________________________________
_PMathObj _TreeTopology::BranchCount(void) {
  long leaves, ints;
  EdgeCount(leaves, ints);
  return new _Constant(ints - 1);
}

//______________________________________________________________________________
_PMathObj _TreeTopology::FlatRepresentation(void) {
  _SimpleList flatTree;

  node<long> *tNode = DepthWiseStepTraverser(theRoot);

  long count = 0;

  while (tNode) {
    flatTree << tNode->in_object;
    tNode->in_object = count++;
    tNode = DepthWiseStepTraverser((node<long> *)nil);
  }

  _Matrix *res = new _Matrix(1, count, false, true);
  checkPointer(res);

  tNode = DepthWiseStepTraverser(theRoot);
  count = 0;

  while (tNode) {
    if (tNode->parent) {
      res->theData[count] = tNode->parent->in_object;
    } else {
      res->theData[count] = -1;
    }

    tNode->in_object = flatTree.lData[count++];
    tNode = DepthWiseStepTraverser((node<long> *)nil);
  }

  return res;
}

//______________________________________________________________________________
_PMathObj _TreeTopology::AVLRepresentation(_PMathObj layoutOption) {

  if (layoutOption->ObjectClass() == NUMBER) {
    bool preOrder = layoutOption->Compute()->Value() > 0.5;

    _AssociativeList *masterList =
        (_AssociativeList *)checkPointer(new _AssociativeList());
    _FString nameHolder;
    //             arrayKey;

    _Constant lengthHolder;

    long rootIndex = 0, nodeLevel = 0;

    _SimpleList nodeList;
    _AVLListX nodeIndexList(&nodeList);

    node<long> *tNode =
        preOrder ? StepWiseTraverser(theRoot) : DepthWiseStepTraverser(theRoot);

    while (tNode) {
      nodeIndexList.Insert((BaseObj *)tNode, nodeIndexList.countitems() + 1);

      if (tNode->parent == nil) {
        rootIndex = nodeIndexList.countitems();
      }

      tNode = preOrder ? StepWiseTraverser((node<long> *)nil)
                       : DepthWiseStepTraverser((node<long> *)nil);
    }

    tNode = preOrder ? StepWiseTraverserLevel(nodeLevel, theRoot)
                     : DepthWiseStepTraverserLevel(nodeLevel, theRoot);

    while (tNode) {
      _AssociativeList *nodeList =
          (_AssociativeList *)checkPointer(new _AssociativeList());
      GetNodeName(tNode, *nameHolder.theString);
      nodeList->MStore("Name", &nameHolder, true);
      GetBranchLength(tNode, lengthHolder.theValue);
      nodeList->MStore("Length", &lengthHolder, true);
      lengthHolder.theValue = nodeLevel;
      nodeList->MStore("Depth", new _Constant(nodeLevel), false);
      if (tNode->parent) {
        nodeList->MStore("Parent",
                         new _Constant(nodeIndexList.GetXtra(
                             nodeIndexList.Find((BaseObj *)tNode->parent))),
                         false);
      }

      long nCount = tNode->get_num_nodes();
      if (nCount) {
        _AssociativeList *childList = new _AssociativeList();
        checkPointer(childList);
        for (long k = 1; k <= nCount; k = k + 1) {
          childList->MStore(
              _String((long)(k - 1)),
              new _Constant(nodeIndexList.GetXtra(
                  nodeIndexList.Find((BaseObj *)tNode->go_down(k)))),
              false);
        }
        nodeList->MStore("Children", childList, false);
      }
      masterList->MStore(_String((long) nodeIndexList.GetXtra(
                             nodeIndexList.Find((BaseObj *)tNode))),
                         nodeList, false);
      tNode =
          preOrder ? StepWiseTraverserLevel(nodeLevel, (node<long> *)nil)
                   : DepthWiseStepTraverserLevel(nodeLevel, (node<long> *)nil);
    }

    _AssociativeList *headerList = new _AssociativeList();
    checkPointer(headerList);

    headerList->MStore("Name", new _FString(*GetName()), false);
    headerList->MStore("Root", new _Constant(rootIndex), false);
    masterList->MStore("0", headerList, false);

    return masterList;
  }
  return new _Constant(0.0);
}

//______________________________________________________________________________
_PMathObj _TreeTopology::TipName(_PMathObj p) {
  _String resString;

  if (p && p->ObjectClass() == NUMBER) {
    long res = p->Value(), count = 0;

    _List *allLeaves = nil;

    if (res < 0) {
      allLeaves = (_List *)checkPointer(new _List);
    }

    LeafWiseT(true);

    while (currentNode) {
      if (res < 0) {
        GetNodeName(currentNode, resString);
        (*allLeaves) && &resString;
      } else if (res == count) {
        //resString = travNode->GetName()->Cut(travNode->GetName()->Find
        //('.')+1,-1);
        GetNodeName(currentNode, resString);
        break;
      }
      count++;
      LeafWiseT(false);
    }

    if (res < 0) {
      _Matrix *res = new _Matrix(*allLeaves);
      DeleteObject(allLeaves);
      return res;
    }
  }

  return new _FString(resString, false);
}

//______________________________________________________________________________
_PMathObj _TreeTopology::BranchLength(_PMathObj p) {
  _Parameter resValue = HY_INVALID_RETURN_VALUE;

  if (p) {
    if (p->ObjectClass() == NUMBER) {
      long res = p->Value(), count = 0;

      if (res < 0)
          // get ALL branch lengths
          {
        EdgeCount(count, res);
        _Matrix *branchLengths =
            (_Matrix *)checkPointer(new _Matrix(1, count + res, false, true));

        count = 0;
        DepthWiseT(true);
        while (!IsCurrentNodeTheRoot()) {
          GetBranchLength(currentNode, branchLengths->theData[count++]);
          DepthWiseT(false);
        }
        return branchLengths;
      } else
          // get a branch length
          {
        DepthWiseT(true);
        while (currentNode && res != count) {
          count++;
          DepthWiseT(false);
        }
        if (currentNode && !IsCurrentNodeTheRoot()) {
          GetBranchLength(currentNode, resValue);
        }
      }
    } else {
      if (p->ObjectClass() == STRING) {
        _List *twoIDs = ((_FString *)p->Compute())->theString->Tokenize(";");
        if (twoIDs->lLength == 2 || twoIDs->lLength == 1) {
          _String *node1 = (_String *)(*twoIDs)(0),
                  *node2 = twoIDs->lLength > 1 ? (_String *)(*twoIDs)(1) : nil;

          node<long> *n1 = nil, *n2 = nil;

          long l1 = 0, l2 = 0, l = 0;

          DepthWiseTLevel(l, true);

          _String cBranchName;

          while (currentNode && (!n1 || !n2)) {
            GetNodeName(currentNode, cBranchName);
            if (cBranchName.Equal(node1)) {
              n1 = currentNode;
              l1 = l;
            } else if (node2 && cBranchName.Equal(node2)) {
              n2 = currentNode;
              l2 = l;
            }

            DepthWiseTLevel(l, false);
          }

          if (n1 && n2) {
            _Parameter p1 = 0, p2 = 0, p;

            while (l1 < l2) {
              GetBranchLength(n2, p);
              p2 += p;
              n2 = n2->parent;
              l2--;
            }

            while (l2 < l1) {
              GetBranchLength(n1, p);
              p1 += p;
              n1 = n1->parent;
              l1--;
            }

            while (n1 != n2) {
              GetBranchLength(n1, p);
              p1 += p;
              GetBranchLength(n2, p);
              p2 += p;
              n2 = n2->parent;
              n1 = n1->parent;
            }

            resValue = p1 + p2;
          } else if (n1)
            if (node2) {
              if (node1->Equal(node2)) {
                resValue = 0.0;
              } else if (node2->Equal(&expectedNumberOfSubs)) {
                _String bl;
                GetBranchLength(n1, bl, true);
                if (bl.sLength) {
                  DeleteObject(twoIDs);
                  return new _FString(bl);
                }
              }
            } else {
              GetBranchLength(n1, resValue);
            }
        }
        DeleteObject(twoIDs);
      }
    }
  }

  if (isnan(resValue)) {
    return new _MathObject();
  }

  return new _Constant(resValue);

}

//______________________________________________________________________________
_PMathObj _TreeTopology::BranchName(_PMathObj p, bool subtree, _PMathObj p2) {
  _String resString;

  if (p) {
    if (p->ObjectClass() == NUMBER) {
      long res = p->Value(), count = -1;

      if (res >= 0) {
        DepthWiseT(true);
        while (currentNode) {
          if (!IsCurrentNodeATip()) {
            count++;
          }

          if (res == count) {
            if (subtree) {
              _String st(128L, true);
              char mapMode = -1;
              if (p2) {
                _String *t = (_String *)p2->Compute()->toStr();
                DetermineBranchLengthMappingMode(t, mapMode);
                DeleteObject(t);
                switch (mapMode) {
                case 3:
                  mapMode = -1;
                  break;
                case 1:
                  mapMode = -3;
                  break;
                case 2:
                  mapMode = -2;
                  break;
                }
              }
              SubTreeString(st, true, mapMode);
              st.Finalize();
              resString = st;
            } else
                //resString = travNode->GetName()->Cut(travNode->GetName()->Find
                //('.')+1,-1);
                {
              GetNodeName(currentNode, resString);
            }
            break;
          }
          DepthWiseT(false);
          if (IsCurrentNodeTheRoot()) {
            return new _MathObject();
          }
        }
      } else {
        count = 0;
        DepthWiseT(true);
        while (currentNode) {
          count++;
          DepthWiseT(false);
        }

        _Matrix *branchLengths = new _Matrix(1, count, false, true);
        checkPointer(branchLengths);
        branchLengths->Convert2Formulas();

        count = 0;
        DepthWiseT(true);

        //long  cutAt = GetName()->sLength+1;
        while (currentNode) {
          _String bs; //(travNode->GetName()->Cut(cutAt,-1));
          GetNodeName(currentNode, bs);
          _FString *bName = new _FString(bs);
          checkPointer(bName);
          _Formula bNamef(bName, false);

          branchLengths->StoreFormula(0, count++, bNamef);
          DepthWiseT(false);
        }
        return branchLengths;
      }
    } else {
      if (p->ObjectClass() == STRING) {
        _List *twoIDs = ((_FString *)p->Compute())->theString->Tokenize(";");

        if (twoIDs->lLength == 2 || twoIDs->lLength == 1) {

          _String *node1 = (_String *)(*twoIDs)(0),
                  *node2 = twoIDs->lLength > 1 ? (_String *)(*twoIDs)(1) : nil;

          if (twoIDs->lLength == 1) {
            _AssociativeList *resList =
                (_AssociativeList *)checkPointer(new _AssociativeList);

            long level = 0, masterLevel = 0;
            StepWiseTLevel(level, true);
            _String givenNodeName;

            while (currentNode) {
              GetNodeName(currentNode, givenNodeName);
              if (givenNodeName.Equal(node1)) {
                masterLevel = level;
                resList->MStore(givenNodeName,
                                new _Constant(currentNode->get_num_nodes()));
                do {
                  GetNodeName(currentNode, givenNodeName);
                  resList->MStore(
                      givenNodeName,
                      new _Constant(1 + (currentNode->get_num_nodes() > 0)));
                  StepWiseTLevel(level, false);
                } while (currentNode && level > masterLevel);
                break;
              }
              StepWiseTLevel(level, false);
            }
            if (resList->avl.countitems() == 0) {
              // fail
              DeleteObject(resList);
              return new _MathObject;
            }
            return resList;
          } else {
            node<long> *n1 = nil, *n2 = nil;

            long l1 = 0, l2 = 0, l = 0;

            DepthWiseTLevel(l, true);

            _String cBranchName;

            while (currentNode && (!n1 || !n2)) {
              GetNodeName(currentNode, cBranchName);
              if (cBranchName.Equal(node1)) {
                n1 = currentNode;
                l1 = l;
              } else if (node2 && cBranchName.Equal(node2)) {
                n2 = currentNode;
                l2 = l;
              }
              DepthWiseTLevel(l, false);
            }

            if (n1 && n2) {
              _List prefix, suffix;

              while (l1 < l2) {
                GetNodeName(n2, cBranchName);
                suffix.AppendNewInstance(cBranchName.makeDynamic());
                n2 = n2->parent;
                l2--;
              }

              while (l2 < l1) {
                GetNodeName(n1, cBranchName);
                prefix.AppendNewInstance(cBranchName.makeDynamic());
                n1 = n1->parent;
                l1--;
              }

              while (n1 != n2) {
                GetNodeName(n2, cBranchName);
                suffix.AppendNewInstance(cBranchName.makeDynamic());
                GetNodeName(n1, cBranchName);
                prefix.AppendNewInstance(cBranchName.makeDynamic());
                n2 = n2->parent;
                n1 = n1->parent;
              }

              suffix.Flip();
              prefix << suffix;
              return new _Matrix(prefix);
            } else if (n1) {
              return new _Matrix();
            }
            return new _MathObject();
          }
        }
      }
    }
  }
  return new _FString(resString);
}

//______________________________________________________________________________
void _TreeTopology::RerootTreeInternalTraverser(long originator,
                                                bool passedRoot, _String &res,
                                                long blOption, bool firstTime) {
  if (passedRoot) {
    SubTreeString(res);
  } else {
    // move to parent now
    node<long> *myParent = currentNode->get_parent(), *saveCurrent;
    _String t;
    if (myParent->get_parent()) { // not root yet
      res << '(';
      saveCurrent = currentNode;
      currentNode = myParent;
      RerootTreeInternalTraverser(currentNode->get_child_num(), false, res,
                                  blOption);
      for (long i = 1; i <= myParent->get_num_nodes(); i++) {
        if (i == originator) {
          continue;
        }
        currentNode = myParent->go_down(i);
        res << ',';
        SubTreeString(res, false, blOption);
      }
      res << ')';
      currentNode = saveCurrent;
      if (!firstTime) {
        GetNodeName(currentNode /*myParent*/, t);
        if (!t.startswith(iNodePrefix)) {
          res << &t;
        }
      }
      PasteBranchLength(currentNode, res, blOption);
    } else { 
      // passing old root
      // create a new root with >=2 children nodes - this, and one more
      // containing all other children (>=2 of them)
      long count = 0, rootChildren = theRoot->get_num_nodes();

      if (rootChildren > 2) {
        res << '(';
      }

      node<long> *stashOriginator;

      for (long k = 1; k <= theRoot->get_num_nodes(); k++) {
        currentNode = theRoot->go_down(k);
        if (k == originator) {
          stashOriginator = currentNode;
          continue;
        }
        if (count) {
          res << ',';
        }
        count++;
        SubTreeString(res, false, blOption);
      }
      if (rootChildren > 2) {
        res << ')';
      }

      /*if (stashOriginator->get_num_nodes())
      {
          GetNodeName (stashOriginator, t);
          if (!t.startswith(iNodePrefix))
              res<<&t;
      }*/
      PasteBranchLength(stashOriginator, res, blOption);
    }
  }
}

//______________________________________________________________________________
void _TreeTopology::SubTreeString(_String &res, bool allNames,
                                  long branchLengths, _AVLListXL *subs) {

  long level = 0, lastLevel = 0, j;
  _String t;
  node<long> *saveCurrent = currentNode;
  currentNode = DepthWiseStepTraverserLevel(level, currentNode);

  while (currentNode) {
    if (level > lastLevel) {
      if (lastLevel) {
        res << ',';
      }
      for (j = 0; j < level - lastLevel; j++) {
        res << '(';
      }
    } else if (level < lastLevel) {
      for (j = 0; j < lastLevel - level; j++) {
        res << ')';
      }
    } else if (lastLevel) {
      res << ',';
    }

    GetNodeName(currentNode, t);
    if (subs) {
      long mapIdx = subs->Find(&t);
      if (mapIdx >= 0) {
        t = *(_String *)subs->GetXtra(mapIdx);
      }
    }

    lastLevel = level;
    if (!IsCurrentNodeTheRoot()) {
      if (allNames || (!t.startswith(iNodePrefix))) {
        res << &t;
      }
      PasteBranchLength(currentNode, res, branchLengths);

    }
    currentNode = DepthWiseStepTraverserLevel(level, (node<long> *)nil);
  }

  currentNode = saveCurrent;
}

//______________________________________________________________________________
void _TreeTopology::PasteBranchLength(node<long> *currentNode, _String &res,
                                      long branchLengths, _Parameter factor) {

  if (branchLengths != -1) {
    _String t;
    if (branchLengths == -2) {
      GetBranchValue(currentNode, t);
    } else if (branchLengths == -3) {
      GetBranchLength(currentNode, t);
    } else {
      GetBranchVarValue(currentNode, t, branchLengths);
    }

    if (t.sLength) {
      t = t.toNum() * factor;
      res << ':';
      res << &t;
    }
  }

}

//______________________________________________________________________________
void _TreeTopology::GetBranchLength(node<long> *n, _String &r, bool getBL) {
  if (getBL) {
    r = empty;
  } else {
    r = compExp->theData[n->in_object];
  }
}

//______________________________________________________________________________
void _TreeTopology::GetBranchLength(node<long> *n, _Parameter &r) {
  r = compExp->theData[n->in_object];
}

//______________________________________________________________________________
void _TreeTopology::GetBranchValue(node<long> *, _String &r) { r = empty; }

//______________________________________________________________________________
void _TreeTopology::GetBranchVarValue(node<long> *, _String &r, long) {
  r = empty;
}

//______________________________________________________________________________
_PMathObj _TreeTopology::RerootTree(_PMathObj p) {

  _String *res = new _String((unsigned long) 256, true);

  iNodePrefix = "Node";
  _PMathObj iv = FetchObjectFromVariableByType(&internalNodePrefix, STRING);

  if (iv) {
    iNodePrefix = *((_FString *)iv)->theString;
  }

  if (p && p->ObjectClass() == STRING) {
    if (rooted == UNROOTED) {
      ReportWarning("Reroot was called with an unrooted tree. Rerooting was "
                    "still performed.");
    }

    _String *tNodeN = (_String *)p->toStr();

    node<long> *rerootNode = FindNodeByName(tNodeN);

    if (rerootNode) { // good node name, can reroot
      if (rerootNode->parent) {
        node<long> *saveCN = rerootNode;
        (*res) << '('; // opening tree (
        RerootTreeInternalTraverser(rerootNode->get_child_num(), false, *res,
                                    -2, true);
        (*res) << ',';
        currentNode = saveCN;
        SubTreeString(*res, 0, -2);
        (*res) << ')';
      } else {
        SubTreeString(*res, 0, -2);
      }
    }
    DeleteObject(tNodeN);
  } else {
    _String errMsg("Reroot Tree was passed an invalid branch argument.");
    WarnError(errMsg);
  }

  res->Finalize();
  return new _FString(res);
}

//______________________________________________________________________________
_String _TreeTopology::DetermineBranchLengthMappingMode(_String *param,
                                                        char &mapMode) {
  mapMode = 3;
  if (param) {
    if (param->Equal(&expectedNumberOfSubs)) {
      mapMode = 1;
    } else if (param->Equal(&stringSuppliedLengths)) {
      mapMode = 2;
    } else {
      mapMode = 0;
      return _String('.') & *param;
    }

  }
  return empty;
}

//______________________________________________________________________________
void _TreeTopology::PostTreeConstructor(bool dupMe) {
  BaseRef temp = variablePtrs(theIndex);

  if (dupMe) {
    variablePtrs[theIndex] = this->makeDynamic();
  } else {
    variablePtrs[theIndex] = this;
  }

  DeleteObject(temp);
}

//______________________________________________________________________________
_TreeTopology::_TreeTopology(_TheTree *top)
    : _CalcNode(*top->GetName(), empty) {
  PreTreeConstructor(false);
  if (top->theRoot) {
    isDefiningATree = 1;
    theRoot = top->theRoot->duplicate_tree();
    node<long> *topTraverser = DepthWiseStepTraverser(theRoot);
    while (topTraverser && topTraverser->parent) {
      _String nodeVS, *nodeSpec, nodeName;

      top->GetBranchValue(topTraverser, nodeVS);
      top->GetNodeName(topTraverser, nodeName);
      nodeSpec = top->GetBranchSpec(topTraverser);

      FinalizeNode(topTraverser, 0, nodeName, *nodeSpec, nodeVS);
      DeleteObject(nodeSpec);
      topTraverser = DepthWiseStepTraverser((node<long> *)nil);
    }
    isDefiningATree = 0;
  } else {
    WarnError("Can't create an empty tree");
    return;
  }
}

//______________________________________________________________________________
// builds a tree from a string
_TreeTopology::_TreeTopology(_String name, _String &parms, bool dupMe)
    : _CalcNode(name, empty) {
  PreTreeConstructor(dupMe);
  if (MainTreeConstructor(parms, false)) {
    PostTreeConstructor(dupMe);
  } else {
    DeleteObject(compExp);
    compExp = nil;
  }
}

//______________________________________________________________________________
_TreeTopology::_TreeTopology(_String *name) : _CalcNode(*name, empty) {}

//______________________________________________________________________________
bool _TreeTopology::MainTreeConstructor(_String &parms, bool checkNames) {
  long i, nodeCount = 0, lastNode;

  _Parameter checkABL;

  checkParameter(acceptBranchLengths, checkABL, 1.0);
  takeBranchLengths = !CheckEqual(checkABL, 0.0);

  checkParameter(autoConvertBL, checkABL, 0.0);
  autoSolveBranchLengths = CheckEqual(checkABL, 1.0);

  _SimpleList nodeStack, nodeNumbers;

  _String nodeName, nodeParameters, nodeValue, nodeComment;

  char lastChar = 0;
  bool isInLiteral = false;

  node<long> *currentNode = theRoot = nil, *newNode = nil, *parentNode = nil;

  isDefiningATree = 1;

  for (i = 0; i < parms.sLength; i++) {
    switch (parms[i]) {
    case '(': { // creating a new internal node one level down
                // a new node
      newNode = new node<long>;
      checkPointer(newNode);
      if (lastChar == '(' || lastChar == ',') {
        currentNode->add_node(*newNode);
      } else {
        if (theRoot) {
          parentNode = currentNode->get_parent();
          if (!parentNode) {
            WarnError((_String("'(' is out of context: ...") &
                       parms.Cut(i > 31 ? i - 32 : 0, i) & "?" &
                       parms.Cut(i + 1, parms.sLength - i > 32 ? i + 32 : -1)));
            isDefiningATree = 0;
            return false;
          } else {
            parentNode->add_node(*newNode);
          }
        } else {
          theRoot = newNode;
        }
        currentNode = newNode;
        nodeStack << (long) currentNode;
        nodeNumbers << nodeCount;
        newNode = new node<long>;
        checkPointer(newNode);
        currentNode->add_node(*newNode);
        nodeCount++;
      }
      currentNode = newNode;
      nodeStack << (long) currentNode;
      nodeNumbers << nodeCount;
      nodeCount++;
      break;
    }

    case ',':
    case ')':
        { // creating a new node on the same level and finishes updating the
          // list of parameters
      lastNode = nodeStack.lLength - 1;
      if (lastNode < 0) {
        WarnError((_String(parms[i]) & _String(" is out of context:") &
                   parms.Cut(i > 31 ? i - 32 : 0, i) & "?" &
                   parms.Cut(i + 1, parms.sLength - i > 32 ? i + 32 : -1)));
        //PurgeTree();
        isDefiningATree = 0;
        return false;
      }
      parentNode = (node<long> *)nodeStack(lastNode);
      FinalizeNode(parentNode, nodeNumbers(lastNode), nodeName, nodeParameters,
                   nodeValue, &nodeComment);
      nodeStack.Delete(lastNode, false);
      nodeNumbers.Delete(lastNode, false);

      if (parms[i] == ',') { // also create a new node on the same level
        checkPointer(newNode = new node<long>);
        if (!(parentNode = parentNode->get_parent())) {
          WarnError((_String("',' is out of context:") &
                     parms.Cut(i > 31 ? i - 32 : 0, i) & "?" &
                     parms.Cut(i + 1, parms.sLength - i > 32 ? i + 32 : -1)));
          isDefiningATree = 0;
          return false;
        }
        currentNode = newNode;
        parentNode->add_node(*currentNode);
        nodeStack << (long) currentNode;
        nodeNumbers << nodeCount;
        nodeCount++;
      }
      break;
    }

    case '{': { // parameter list definition
      lastNode = parms.Find("}", i + 1, -1);
      if (lastNode < 0) {
        WarnError((_String("'{' has no matching '}':") &
                   parms.Cut(i > 31 ? i - 32 : 0, i) & "?" &
                   parms.Cut(i + 1, parms.sLength - i > 32 ? i + 32 : -1)));
        isDefiningATree = 0;
        return false;
      }
      nodeParameters = parms.Cut(i + 1, lastNode - 1);
      i = lastNode;
      break;
    }

    case '[': { // hackish Newick annotation
      lastNode = parms.Find("]", i + 1, -1);
      if (lastNode < 0) {
        WarnError((_String("'[' has no matching ']':") &
                   parms.Cut(i > 31 ? i - 32 : 0, i) & "?" &
                   parms.Cut(i + 1, parms.sLength - i > 32 ? i + 32 : -1)));
        isDefiningATree = 0;
        return false;
      }
      nodeComment = parms.Cut(i + 1, lastNode - 1);
      i = lastNode;
      break;
    }

    case ':': { // tree branch definition
      lastNode = i + 1;
      char c = parms[lastNode];

      while (isspace(c))
        if (lastNode < parms.sLength) {
          c = parms[++lastNode];
        } else {
          break;
        }

      if (lastNode < parms.sLength)
        while ((c <= '9' && c >= '0') || c == '.' || c == '-' || c == '+' ||
               c == 'e' || c == 'E') {
          if (lastNode < parms.sLength) {
            c = parms[++lastNode];
          } else {
            break;
          }
        }

      nodeValue = parms.Cut(i, lastNode - 1);
      i = lastNode - 1;
      break;
    }

    default: { // node name
      lastNode = i;

      char c = parms[lastNode];

      if (c == ';') {
        break;
      }

      if (isspace(c)) {
        continue;
      }

      if (c == '\'') {
        c = parms[lastNode++];
        isInLiteral = true;
        i++;
      }

      if (!isInLiteral && !(isalnum(c) || c == '_')) {
        WarnError(
            (_String("Node names should begin with a letter, a number, or an "
                     "underscore. Had:") & parms.Cut(i > 31 ? i - 32 : 0, i) &
             "?" & parms.Cut(i + 1, parms.sLength - i > 32 ? i + 32 : -1)));
        isDefiningATree = 0;
        return false;
      }

      if (checkNames)
        while (isalnum(c) || c == '_')
          if (lastNode < parms.sLength) {
            lastNode++;
            c = parms[lastNode];
          } else {
            break;
          }
      else
        while (isInLiteral || !(c == ',' || c == ':' || c == ')' || c == '(' ||
                                c == '{' || c == '}' || isspace(c)))
          if (lastNode < parms.sLength) {
            lastNode++;
            c = parms[lastNode];
            if (c == '\'') {
              if (isInLiteral) {
                break;
              } else {
                WarnError(
                    (_String("Unxpected \'. Had:") &
                     parms.Cut(lastNode > 31 ? lastNode - 32 : 0, lastNode) &
                     "?" & parms.Cut(lastNode + 1, parms.sLength - lastNode > 32
                                                       ? lastNode + 32
                                                       : -1)));
                isDefiningATree = 0;
                return false;
              }
            }
          } else {
            break;
          }

      if (isInLiteral) {
        if (c != '\'') {
          WarnError(
              (_String("Unterminated \'. Had:") &
               parms.Cut(lastNode > 31 ? lastNode - 32 : 0, lastNode) & "?" &
               parms.Cut(lastNode + 1,
                         parms.sLength - lastNode > 32 ? lastNode + 32 : -1)));
          isDefiningATree = 0;
          return false;
        }
        nodeName = parms.Cut(i, lastNode - 1);
        i = lastNode;
        isInLiteral = false;
      } else {
        nodeName = parms.Cut(i, lastNode - 1);
        i = lastNode - 1;
      }
      break;
    }
    }

    if ((lastChar = parms[i]) == ';') {
      break;
    }
  }

  lastNode = nodeStack.lLength - 1;

  while (lastNode >= 0) {
    parentNode = (node<long> *)nodeStack(lastNode);
    FinalizeNode(parentNode, nodeNumbers(lastNode), nodeName, nodeParameters,
                 nodeValue, &nodeComment);
    lastNode--;
  }

  if (!theRoot) {
    isDefiningATree = 0;
    WarnError("Can't create empty trees.");
    return false;
  }

  isDefiningATree = 0;
  return true;

}

//______________________________________________________________________________
bool _TreeTopology::FinalizeNode(node<long> *nodie, long number,
                                 _String &nodeName, _String &nodeParameters,
                                 _String &nodeValue, _String *nodeComment) {
  if (!nodeName.sLength ||
      !CheckEqual(ignoringInternalNames, 0.0) && nodie->get_num_nodes() > 0) {
    nodeName = iNodePrefix & number;
  }

  if (nodie == theRoot) {
    nodeParameters = "";
    nodeValue = "";
  }

  nodie->in_object = flatTree.lLength;
  flatTree &&&nodeName;
  flatCLeaves &&&nodeParameters;

  ((_GrowingVector *)compExp)->Store(nodeValue.ProcessTreeBranchLength());

  nodeName = empty;
  nodeParameters = empty;
  nodeValue = empty;
  if (nodeComment)
    *nodeComment = empty;

  return true;
}

//______________________________________________________________________________
node<long> *_TreeTopology::FindNodeByName(_String *match) {
  DepthWiseT(true);

  _String nn;
  while (currentNode) {
    GetNodeName(currentNode, nn);
    if (match->Equal(&nn)) {
      return currentNode;
    }
    DepthWiseT();
  }

  return nil;

}

//______________________________________________________________________________
void _TreeTopology::AddANode(_PMathObj newNode) {
  if (newNode->ObjectClass() == ASSOCIATIVE_LIST) {
    _AssociativeList *newNodeSpec = (_AssociativeList *)newNode;
    _FString *newName =
                 (_FString *)newNodeSpec->GetByKey(newNodeGraftName, STRING),
             *newLocation =
                 (_FString *)newNodeSpec->GetByKey(newNodeGraftWhere, STRING),
             *newParent =
                 (_FString *)newNodeSpec->GetByKey(newNodeGraftParent, STRING);

    if (!newName) {
      WarnError(_String("Missing/invalid mandatory argument (\"") &
                newNodeGraftName & "\") in call to _TreeTopology::AddANode");
      return;
    }
    if (!newLocation) {
      WarnError(_String("Missing/invalid mandatory argument (\"") &
                newNodeGraftWhere & "\") in call to _TreeTopology::AddANode");
      return;
    }
    if (!newParent) {
      WarnError(_String("Missing/invalid mandatory argument (\"") &
                newNodeGraftParent & "\") in call to _TreeTopology::AddANode");
      return;
    }

    node<long> *graftAt = FindNodeByName(newLocation->theString);
    if (!graftAt || graftAt->get_parent() == nil) {
      WarnError("Attachment node must be an exiting non-root node in call to "
                "_TreeTopology::AddANode");
      return;
    }

    node<long> *newp = (node<long> *)checkPointer(new node<long>),
               *curp = graftAt->get_parent();

    newp->set_parent(*curp);
    newp->add_node(*graftAt);
    curp->replace_node(graftAt, newp);

    if (!newName->IsEmpty()) {
      node<long> *newt = (node<long> *)checkPointer(new node<long>);
      newp->add_node(*newt);
      FinalizeNode(newt, 0, *newName->theString, empty, empty);
    }

    FinalizeNode(newp, 0, *newParent->theString, empty, empty);

  } else {
    WarnError("An invalid argument (not an associative array) supplied to "
              "_TreeTopology::AddANode");
  }

}

//______________________________________________________________________________
void _TreeTopology::DepthWiseT(bool init, _HYTopologyTraversalFunction *handler,
                               Ptr extra) {
  if (init) {
    currentNode = DepthWiseStepTraverser(theRoot);
  } else {
    currentNode = (DepthWiseStepTraverser((node<long> *)nil));
  }

  if (handler)
    if (!(*handler)(currentNode, extra)) {
      currentNode = nil;
    }

}

//______________________________________________________________________________
void _TreeTopology::DepthWiseTRight(bool init) {
  if (init) {
    currentNode = DepthWiseStepTraverserRight(theRoot);
  } else {
    currentNode = (DepthWiseStepTraverserRight((node<long> *)nil));
  }
}

//______________________________________________________________________________
void _TreeTopology::LeafWiseT(bool init) {
  if (init) {
    currentNode = DepthWiseStepTraverser(theRoot);
  } else {
    currentNode = DepthWiseStepTraverser((node<long> *)nil);
  }

  while (currentNode && currentNode->get_num_nodes()) {
    currentNode = DepthWiseStepTraverser((node<long> *)nil);
  }
}

//______________________________________________________________________________
void _TreeTopology::StepWiseT(bool init, _HYTopologyTraversalFunction *handler,
                              Ptr extra) {
  if (init) {
    currentNode = StepWiseTraverser(theRoot);
  } else {
    currentNode = (StepWiseTraverser((node<long> *)nil));
  }

  if (handler)
    if (!(*handler)(currentNode, extra)) {
      currentNode = nil;
    }
}

//______________________________________________________________________________
void _TreeTopology::StepWiseTLevel(long &level, bool init) {
  if (init) {
    currentNode = StepWiseTraverserLevel(level, theRoot);
  } else {
    currentNode = (StepWiseTraverserLevel(level, (node<long> *)nil));
  }

}

//______________________________________________________________________________
void _TreeTopology::DepthWiseTLevel(long &level, bool init) {
  if (init) {
    currentNode = DepthWiseStepTraverserLevel(level, theRoot);
  } else {
    currentNode = (DepthWiseStepTraverserLevel(level, (node<long> *)nil));
  }
}

//______________________________________________________________________________
BaseRef _TreeTopology::makeDynamic(void) {
  _TreeTopology *res = new _TreeTopology;
  checkPointer(res);
  res->_CalcNode::Duplicate(this);

  res->flatTree.Duplicate(&flatTree);
  res->flatCLeaves.Duplicate(&flatCLeaves);
  if (compExp) {
    res->compExp = (_Matrix *)compExp->makeDynamic();
  } else {
    res->compExp = nil;
  }

  res->currentNode = currentNode;
  res->rooted = rooted;
  res->theRoot = CopyTreeStructure(theRoot, true);
  return res;
}

//______________________________________________________________________________
_String _TreeTopology::CompareTrees(_TreeTopology *compareTo) {
  _List myLeaves, otherLeaves;

  _SimpleList indexer, otherIndexer;

  node<long> *myCT, *otherCT;

  _String rerootAt;

  myCT = prepTree4Comparison(myLeaves, indexer);
  otherCT = compareTo->prepTree4Comparison(otherLeaves, otherIndexer);

  // first compare the set of leaf labels

  if (myLeaves.Equal(otherLeaves)) {
    _SimpleList *reindexer = nil;

    if (!indexer.Equal(otherIndexer)) {
      _SimpleList ilist((unsigned long) myLeaves.lLength);
      ilist.lLength = myLeaves.lLength;

      for (long k2 = 0; k2 < indexer.lLength; k2++) {
        ilist.lData[indexer.lData[k2]] = k2;
      }

      for (long k3 = 0; k3 < otherIndexer.lLength; k3++) {
        otherIndexer.lData[k3] = ilist.lData[otherIndexer.lData[k3]];
      }

      reindexer = &otherIndexer;
    }

    char compRes;

    if ((compRes = internalTreeCompare(myCT, otherCT, reindexer, 1,
                                       myLeaves.lLength, nil, compareTo)) > 0) {
      rerootAt = eqWithoutReroot;
    } else {
      long tCount = 0;

      node<long> *meNode = DepthWiseStepTraverser(otherCT);

      while (meNode != otherCT) {
        if (meNode->get_num_nodes()) {
          compRes = internalTreeCompare(myCT, meNode, reindexer, 1,
                                        myLeaves.lLength, nil, compareTo);
          if (compRes > 0) {
            break;
          } else if (compRes) {
            meNode = otherCT;
            break;
          }
        }

        tCount++;
        meNode = DepthWiseStepTraverser((node<long> *)nil);
      }

      if (meNode != otherCT) {
        meNode = DepthWiseStepTraverser(compareTo->theRoot);
        while (meNode != theRoot) {
          if (tCount == 1) {
            rerootAt =
                eqWithReroot & *LocateVar(meNode->in_object)->GetName() & '.';
            break;
          } else {
            tCount--;
          }

          meNode = DepthWiseStepTraverser((node<long> *)nil);
        }
      }
    }
    if (!rerootAt.sLength) {
      rerootAt = "Unequal topologies (matching label sets).";
    }
  } else {
    rerootAt = "Unequal label sets.";
  }

  destroyCompTree(myCT);
  destroyCompTree(otherCT);

  return rerootAt;
}

//______________________________________________________________________________
node<long> *_TreeTopology::prepTree4Comparison(_List &leafNames,
                                               _SimpleList &mapping,
                                               node<long> *topNode) {
  node<long> *res = topNode ? topNode->duplicate_tree()
                            : theRoot->duplicate_tree(),
             *meNode;

  checkPointer(res);

  meNode = DepthWiseStepTraverser(res);

  _SimpleList indexer;

  while (meNode) {
    long numChildren = meNode->get_num_nodes();

    _SimpleList *descendants = new _SimpleList;
    checkPointer(descendants);
    if (numChildren) {
      for (long k = 1; k <= numChildren; k++) {
        node<long> *aChild = meNode->go_down(k);
        (*descendants) << *(_SimpleList *)aChild->in_object;
      }
    } else {
      (*descendants) << leafNames.lLength;
      indexer << leafNames.lLength;
      _String *dd = (_String *)checkPointer(
          new _String); // (*LocateVar
                        // (meNode->in_object)->GetName(),treeNameLength,-1);
      GetNodeName(meNode, *dd);
      leafNames.AppendNewInstance(dd);
    }

    meNode->in_object = (long) descendants;
    meNode = DepthWiseStepTraverser((node<long> *)nil);
  }

  // now sort leaf names and build the indexer

  mapping.Clear();
  mapping.Duplicate(&indexer);

  SortLists(&leafNames, &indexer);
  SortLists(&indexer, &mapping);

  return res;
}

//______________________________________________________________________________
void _TreeTopology::destroyCompTree(node<long> *compRoot) {
  long nc = compRoot->get_num_nodes();
  for (int i = 1; i <= nc; i++) {
    destroyCompTree(compRoot->go_down(i));
  }

  DeleteObject((BaseRef) compRoot->in_object);
  delete (compRoot);
}

// returns the list of nodes included in the splits BELOW trav1
_List *_TreeTopology::SplitTreeIntoClustersInt(node<long> *, _List *,
                                               _AVLListX &, long, long) {
  /*_List * dependancies = nil;

    for (long k = trav1->get_num_nodes(); k; k--)
    {
        _List * me = SplitTreeIntoClustersInt(trav1->go_down(k), res, counts,
  size, tol);
        if (me)
        {
            if (!dependancies)
                checkPointer (dependancies = new _List);

            (*dependancies) << (*me);
            DeleteObject (me);
        }
    }

    long distinctElements = counts.GetXtra
  (counts.Find((BaseRef)trav1->in_object)),
         dLength = dependancies?dependancies->lLength:0;

    if
  (((trav1->parent==nil)||(distinctElements+dLength>=size-tol))&&(distinctElements+dLength<=size+tol))
    {
        _List   entry;
        _String nName;
        GetNodeName (trav1, nName);
        entry && & nName;
        entry && & _Constant (distinctElements);
        if (dependancies)
            entry << (*dependancies);
        (*res) && & entry;

        trav1 = trav1->parent;
        while (trav1)
        {
            long aidx = counts.Find((BaseRef)trav1->in_object);
            counts.SetXtra(aidx,counts.GetXtra (aidx) - distinctElements);
            trav1 = trav1->parent;
        }

        if (dependancies)
            dependancies->Clear();
        else
            checkPointer (dependancies = new _List);

        (*dependancies) && & nName;
    }

    return dependancies;*/

  // decide if child or i-node

  /*if (trav1->get_num_nodes())
    {

    }
    else
    {
        node <long>* pn = trav1->parent;

        while (pn)
        {
            long pidx = counts.Find((BaseRef)pn->in_object);

        }
    }*/
  // TBI
  return new _List;
}

//______________________________________________________________________________
// assume fixed rooting for now
// returns the list of list speccing pointers to calcnodes rooting the
// splits
// the member list contains 2 or more entries for each node:
// itself, the number of (independent)leaves encompassed by that node,
// and an optional references to the node whose
// results it depends on.
// assumes that size-tol>=2
_List *_TreeTopology::SplitTreeIntoClusters(unsigned long size,
                                            unsigned long tol) {
  _SimpleList counts;
  _AVLListX cavl(&counts);

  DepthWiseT(true);

  while (currentNode) {
    long nC = currentNode->get_num_nodes();
    if (nC) {
      long c = 0;
      for (long k = 1; k <= nC; k++) {
        c += counts.lData[currentNode->go_down(k)->in_object];
      }
      cavl.Insert((BaseRef) currentNode->in_object, c);
    } else {
      cavl.Insert((BaseRef) currentNode->in_object, 1);
    }

    DepthWiseT(false);
  }

  _List *result = new _List;
  checkPointer(result);

  DeleteObject(SplitTreeIntoClustersInt(theRoot, result, cavl, size, tol));

  return result;
}

//______________________________________________________________________________
// the pattern is this
// compare to is the potential match
_String _TreeTopology::MatchTreePattern(_TreeTopology *compareTo) {

  _List myLeaves, otherLeaves, overlappedLeaves;
  _SimpleList indexer, otherIndexer;
  node<long> *myCT, *otherCT;
  _String rerootAt;

  myCT = prepTree4Comparison(myLeaves, indexer);
  otherCT = compareTo->prepTree4Comparison(otherLeaves, otherIndexer);

  _SimpleList matchedLeaves;

  overlappedLeaves.Intersect(myLeaves, otherLeaves, &matchedLeaves);

  if ((myLeaves.lLength >= otherLeaves.lLength) &&
      (overlappedLeaves.lLength == otherLeaves.lLength)) {
    if (myLeaves.lLength > otherLeaves.lLength) {
      // map leaves back to ordering space

      // trim the pattern tree
      _SimpleList allowedLeaves((unsigned long) myLeaves.lLength),
          recordTransfer((unsigned long) myLeaves.lLength),
          invertedMap((unsigned long) myLeaves.lLength);

      allowedLeaves.lLength = myLeaves.lLength;
      recordTransfer.lLength = myLeaves.lLength;
      invertedMap.lLength = myLeaves.lLength;

      for (long giantsSuck = 0; giantsSuck < myLeaves.lLength; giantsSuck++) {
        invertedMap.lData[indexer.lData[giantsSuck]] = giantsSuck;
      }

      for (long padresSuck = 0; padresSuck < matchedLeaves.lLength;
           padresSuck++) {
        allowedLeaves
            .lData[invertedMap.lData[matchedLeaves.lData[padresSuck]]] = 1;
      }

      long slider = 0;

      for (long dodgersSuck = 0; dodgersSuck < recordTransfer.lLength;
           dodgersSuck++)
        if (allowedLeaves.lData[dodgersSuck]) {
          recordTransfer[dodgersSuck] = slider++;
        }

      // pass 1 - delete all the superfluous leaves
      node<long> *padresStillSuck = DepthWiseStepTraverser(myCT);
      while (padresStillSuck != myCT) {
        if (padresStillSuck->get_num_nodes() == 0) {
          _SimpleList *descendants =
              ((_SimpleList *)padresStillSuck->in_object);

          if ((descendants && (allowedLeaves.lData[descendants->lData[0]] ==
                               0)) || (!descendants)) {
            // mark this node for deletion
            if (descendants) {
              DeleteObject(descendants);
            }

            node<long> *sacLamb = padresStillSuck;
            padresStillSuck = DepthWiseStepTraverser((node<long> *)nil);

            if (sacLamb->parent->get_num_nodes() == 1) {
              DeleteObject((BaseRef) sacLamb->parent->in_object);
              sacLamb->parent->in_object = nil;
            }

            sacLamb->parent->detach_child(sacLamb->get_child_num());
            delete (sacLamb);
            continue;
          }
        }
        padresStillSuck = DepthWiseStepTraverser((node<long> *)nil);
      }

      // pass 2 - prune internal nodes with exactly one child

      // >O
      padresStillSuck = DepthWiseStepTraverser(myCT);

      while (padresStillSuck) {
        long nn = padresStillSuck->get_num_nodes();
        if (nn == 1) {
          node<long> *loneChild = padresStillSuck->go_down(1);
          DeleteObject((_SimpleList *)padresStillSuck->in_object);
          padresStillSuck->in_object = loneChild->in_object;
          padresStillSuck->detach_child(1);
          delete (loneChild);
        } else {
          if (nn > 1) {
            _SimpleList *myDescs = (_SimpleList *)padresStillSuck->in_object;
            myDescs->Clear();

            for (long cc = 1; cc <= nn; cc++) {
              _SimpleList temp;

              temp.Union(*myDescs, *(_SimpleList *)padresStillSuck->go_down(cc)
                                       ->in_object);

              myDescs->Clear();
              myDescs->Duplicate(&temp);
            }

            /*for (long cc2 = 0; cc2 < myDescs->lLength; cc2++)
                myDescs->lData[cc2] =
                recordTransfer.lData[myDescs->lData[cc2]];*/
          } else {
            ((_SimpleList *)padresStillSuck->in_object)
                ->lData[0] = recordTransfer
                .lData[((_SimpleList *)padresStillSuck->in_object)->lData[0]];
          }

        }
        padresStillSuck = DepthWiseStepTraverser((node<long> *)nil);
      }

      _List newLeaves;
      recordTransfer.Clear();
      invertedMap.Clear();

      for (long lc = 0; lc < allowedLeaves.lLength; lc++)
        if (allowedLeaves.lData[lc]) {
          recordTransfer << newLeaves.lLength;
          invertedMap << newLeaves.lLength;
          newLeaves << myLeaves(indexer.lData[lc]);
        }

      SortLists(&newLeaves, &recordTransfer);
      SortLists(&recordTransfer, &invertedMap);
      indexer.Clear();
      indexer.Duplicate(&invertedMap);
      myLeaves.Clear();
      myLeaves.Duplicate(&newLeaves);

      // finally check whether the root still has 3 or more children

      if (myCT->get_num_nodes() == 2) {
        node<long> *promoteMe = nil;

        if (myCT->go_down(1)->get_num_nodes()) {
          promoteMe = myCT->go_down(1);
        } else {
          promoteMe = myCT->go_down(2);
        }

        long nn = promoteMe->get_num_nodes();
        if (nn) {
          for (long cc = 1; cc <= nn; cc++) {
            myCT->add_node(*promoteMe->go_down(cc));
          }

          myCT->detach_child(promoteMe->get_child_num());
          DeleteObject((BaseRef) promoteMe->in_object);
          delete promoteMe;
        } else {
          WarnError("Internal tree pattern error in MatchTreePattern");
          return "Unequal: Error";
        }
      }
    }

    _SimpleList *reindexer = nil;

    if (!indexer.Equal(otherIndexer)) {
      _SimpleList ilist((unsigned long) myLeaves.lLength);
      ilist.lLength = myLeaves.lLength;

      for (long k2 = 0; k2 < indexer.lLength; k2++) {
        ilist.lData[indexer.lData[k2]] = k2;
      }

      for (long k3 = 0; k3 < otherIndexer.lLength; k3++) {
        otherIndexer.lData[k3] = ilist.lData[otherIndexer.lData[k3]];
      }

      reindexer = &otherIndexer;
    }

    char compRes;

    if ((compRes =
             internalTreeCompare(myCT, otherCT, reindexer, 1, myLeaves.lLength,
                                 nil, compareTo, true)) > 0) {
      rerootAt = "Equal w/o rerooting.";
    } else {
      long tCount = 0;

      node<long> *meNode = DepthWiseStepTraverser(otherCT);

      while (meNode != otherCT) {
        if (meNode->get_num_nodes()) {
          compRes = internalTreeCompare(myCT, meNode, reindexer, 1,
                                        myLeaves.lLength, nil, compareTo, true);
          if (compRes > 0) {
            break;
          } else if (compRes) {
            meNode = otherCT;
            break;
          }
        }

        tCount++;
        meNode = DepthWiseStepTraverser((node<long> *)nil);
      }

      if (meNode != otherCT) {
        meNode = DepthWiseStepTraverser(compareTo->theRoot);
        while (meNode != theRoot) {
          if (tCount == 1) {
            rerootAt = _String("Equal with reroot at ") &
                       *LocateVar(meNode->in_object)->GetName() & '.';
            break;
          } else {
            tCount--;
          }

          meNode = DepthWiseStepTraverser((node<long> *)nil);
        }
      }
    }
    if (!rerootAt.sLength) {
      rerootAt = "Unequal topologies (matching label sets).";
    }
  } else {
    rerootAt = "Unequal label sets.";
  }

  destroyCompTree(myCT);
  destroyCompTree(otherCT);
  return rerootAt;
}


//______________________________________________________________________________
// Tree cluster comparison functions
//______________________________________________________________________________
void _TreeTopology::ComputeClusterTable(_SimpleList &result,
                                        _SimpleList &pswRepresentation) {

  long leafCount = pswRepresentation.Element(-2), leafCode = 0, L, R;

  result.Clear();
  result.Populate(3 * leafCount, -1, 0);

  for (long k = 0; k < pswRepresentation.lLength - 2; k += 2) {
    if (pswRepresentation.lData[k] < leafCount) { // is a leaf
      R = leafCode++;
    } else {
      long row;
      L = pswRepresentation.lData[k - 2 * pswRepresentation.lData[k + 1]];
      if (k == pswRepresentation.lLength - 4 /* root */
          || pswRepresentation.lData[k + 3] == 0) {
        row = R;
      } else {
        row = L;
      }

      result.lData[row * 3] = L;
      result.lData[row * 3 + 1] = R;
    }
  }
}

//______________________________________________________________________________
_String *_TreeTopology::ConvertFromPSW(_AVLListX &nodeMap,
                                       _SimpleList &pswRepresentation) {
  _String *result = new _String(128L, true);
  if (pswRepresentation.lLength > 4) {
    long leafCount = pswRepresentation.Element(-2);
    // traverse backwards
    bool lastLeaf = false;
    _SimpleList bounds;

    for (long k = pswRepresentation.lLength - 4; k >= 0; k -= 2) {
      if (lastLeaf) {
        (*result) << ',';
      }
      if (pswRepresentation.lData[k] >= leafCount) { //
        (*result) << ')';
        lastLeaf = false;
        bounds << k - 2 * pswRepresentation.lData[k + 1];
      } else {
        _String nodeName(
            *(_String *)nodeMap.Retrieve(pswRepresentation.lData[k]));
        nodeName.Flip();
        (*result) << nodeName;
        while (bounds.Element(-1) == k && bounds.lLength) {
          (*result) << '(';
          bounds.Pop();
        }
        lastLeaf = true;
      }
    }
  }
  result->Finalize();
  result->Flip();
  return result;
}

//______________________________________________________________________________
bool _TreeTopology::ConvertToPSW(_AVLListX &nodeMap, _List *inames,
                                 _SimpleList &pswRepresentation,
                                 bool reference) {
  if (reference == false) {
    nodeMap.Clear();
  }

  pswRepresentation.Clear();

  long leafIndex = 0, iNodeCount = -1, level = 0;

  _String nodeName;

  DepthWiseTLevel(level, &GetRoot());
  _SimpleList levelBuffer;

  while (currentNode) {
    GetNodeName(currentNode, nodeName, false);

    while (levelBuffer.countitems() <= level) {
      levelBuffer << 0;
    }

    if (IsCurrentNodeATip()) {

      pswRepresentation << leafIndex;
      pswRepresentation << 0;
      if (reference) {
        long remapped = nodeMap.Find(&nodeName);
        if (remapped < 0) {
          return false;
        } else {
          remapped = nodeMap.GetXtra(remapped);
          if (remapped >= 0) {
            pswRepresentation << remapped;
          } else {
            return false;
          }
        }

        leafIndex++;
      } else {
        nodeMap.Insert(nodeName.makeDynamic(), leafIndex++, false);
      }
    } else {
      pswRepresentation << iNodeCount;
      pswRepresentation << levelBuffer.lData[level];
      if (reference) {
        pswRepresentation << 0;
      } else {
        (*inames) && &nodeName;
      }

      iNodeCount--;
    }
    if (level) {
      levelBuffer.lData[level - 1] += levelBuffer.lData[level] + 1;
    }
    levelBuffer.lData[level] = 0;
    DepthWiseTLevel(level, nil);
  }

  for (long k = 0; k < pswRepresentation.lLength; k += (reference ? 3 : 2))
    if (pswRepresentation.lData[k] < 0) {
      pswRepresentation.lData[k] = leafIndex - pswRepresentation.lData[k] - 1;
    }

  pswRepresentation << leafIndex;
  pswRepresentation << (-iNodeCount - 1);

  return true;
}

//______________________________________________________________________________
// compare tree topologies
_AssociativeList *_TreeTopology::SplitsIdentity(_PMathObj p) {

  _Matrix *result = (_Matrix *)checkPointer(new _Matrix(2, 1, false, true)),
          *result2 = nil;

  _FString *treeR = new _FString();

  _Constant *bc = (_Constant *)BranchCount();
  result->theData[0] = bc->Value();
  result->theData[1] = -1;

  if (p && (p->ObjectClass() == TOPOLOGY || p->ObjectClass() == TREE)) {
    _List avlSupport, iNames;

    _AVLListX nameMap(&avlSupport);

    _SimpleList psw, psw2, clusters, inodeList;

    ConvertToPSW(nameMap, &iNames, psw);
    ComputeClusterTable(clusters, psw);
    if (((_TreeTopology *)p)->ConvertToPSW(nameMap, nil, psw2, true)) {
      _SimpleList workSpace;
      long leafCount = psw.Element(-2);

      for (unsigned long k = 0; k < psw2.lLength - 3; k += 3) {

        if (psw2.lData[k] < leafCount) {
          workSpace << 1;
          workSpace << 1;
          workSpace << psw2.lData[k + 2];
          workSpace << psw2.lData[k + 2];
        } else {
          _SimpleList quad;
          quad << leafCount + 1;
          quad << 0;
          quad << 0;
          quad << 1;

          long w = psw2.lData[k + 1];
          while (w > 0) {
            _SimpleList quad2;
            quad2 << workSpace.Pop();
            quad2 << workSpace.Pop();
            quad2 << workSpace.Pop();
            quad2 << workSpace.Pop();
            w -= quad2.lData[3];
            quad.lData[0] = MIN(quad2.lData[0], quad.lData[0]);
            quad.lData[1] = MAX(quad2.lData[1], quad.lData[1]);
            quad.lData[2] += quad2.lData[2];
            quad.lData[3] += quad2.lData[3];
          }

          if (quad.lData[2] == quad.lData[1] - quad.lData[0] + 1) {
            if (clusters.lData[3 * quad.lData[0]] == quad.lData[0] &&
                clusters.lData[3 * quad.lData[0] + 1] == quad.lData[1]) {
              clusters.lData[3 * quad.lData[0] + 2] = 1;
            } else if (clusters.lData[3 * quad.lData[1]] == quad.lData[0] &&
                       clusters.lData[3 * quad.lData[1] + 1] == quad.lData[1]) {
              clusters.lData[3 * quad.lData[1] + 2] = 1;
            }
          }
          quad.Flip();
          workSpace << quad;
        }
      }

      psw2.Clear();
      long matchCount = 0, iNodeCount = 0;

      long L, R;

      _SimpleList leafSpans(leafCount, 0, 0), iNodesTouched;

      for (unsigned long k = 0; k < psw.lLength - 2; k += 2) {
        if (psw.lData[k] < leafCount) {
          R = psw.lData[k];
          psw2 << R;
          psw2 << 0;
          leafSpans.lData[R] = (psw2.lLength >> 1);
        } else {
          long ll = k - 2 * psw.lData[k + 1];
          L = psw.lData[ll];
          if ((clusters.lData[3 * L] == L && clusters.lData[3 * L + 1] == R &&
               clusters.lData[3 * L + 2] > 0) ||
              (clusters.lData[3 * R] == L && clusters.lData[3 * R + 1] == R &&
               clusters.lData[3 * R + 2] > 0)) {
            L = (psw2.lLength >> 1) - leafSpans.lData[L] + 1;
            psw2 << leafCount + iNodeCount++;
            psw2 << L;

            iNodesTouched << psw.lData[k];
          }
        }
      }

      for (unsigned long k = 0; k < psw2.lLength; k += 2)
        if (psw2.lData[k] < leafCount) {
          psw2.lData[k + 1] = 0;
        } else {
          matchCount++;
        }

      psw2 << leafCount;
      psw2 << iNodeCount;

      result->theData[0] = psw.Element(-1);
      result->theData[1] = matchCount;

      *treeR->theString = ConvertFromPSW(nameMap, psw2);

      _List sharedNames;
      for (long k = 0; k < iNodesTouched.lLength; k++) {
        sharedNames << iNames(iNodesTouched(k) - leafCount);
      }

      result2 = new _Matrix(sharedNames);
    }

  }

  DeleteObject(bc);

  _AssociativeList *resultList = new _AssociativeList;
  resultList->MStore("CLUSTERS", result, false);
  if (result2) {
    resultList->MStore("NODES", result2, false);
  }
  resultList->MStore("CONSENSUS", treeR, false);
  return resultList;
}
