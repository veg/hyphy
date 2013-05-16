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

#include <tr1/tuple>
#include <iostream>
#include "gtest/gtest.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

//Generate necessary includes from the respective implementation file
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

namespace {

// The fixture for testing class Foo.
class _TreeTopologyTest : public ::testing::Test {

protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  _TreeTopologyTest() {
    // You can do set-up work for each test here.
    // Create objects of every type needed. Performance doesn't matter.

    FILEtest = fopen ("./tests/gtests/res/HIV_gp120.nex" , "r");

    longtest = new long();
    unsignedtest = new unsigned();
    _TreeTopologytest = new _TreeTopology();
    _hyExecutionContexttest = new _hyExecutionContext(_VariableContainertest, _Stringtest);
    _Listtest = new _List();
    _HYTopologyTraversalFunctiontest = new _HYTopologyTraversalFunction();
    node<long>test = new node<long>();
    _Stringtest = new _String(FILEtest);
    _Matrixtest = new _Matrix();
    _AVLListXtest = new _AVLListX(_SimpleListtest);
    _SimpleListtest = new _SimpleList();
    booltest = new bool();
    _Parametertest = new _Parameter();
    _PMathObjtest = new _PMathObj();
  }

  virtual ~_TreeTopologyTest() {
    // You can do clean-up work that doesn't throw exceptions here.
  }

  // If the constructor and destructor are not enough for setting up
  // and cleaning up each test, you can define the following methods:

  virtual void SetUp() {
    // Code here will be called immediately after the constructor (right
    // before each test).
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test (right
    // before the destructor).
    delete longtest;
    delete unsignedtest;
    delete _TreeTopologytest;
    delete _hyExecutionContexttest;
    delete _Listtest;
    delete _HYTopologyTraversalFunctiontest;
    delete node<long>test;
    delete _Stringtest;
    delete _Matrixtest;
    delete _AVLListXtest;
    delete _SimpleListtest;
    delete booltest;
    delete _Parametertest;
    delete _PMathObjtest;
    fclose (FILEtest);
  }

  FILE* FILEtest;
  long* longtest;
  unsigned* unsignedtest;
  _TreeTopology* _TreeTopologytest;
  _hyExecutionContext* _hyExecutionContexttest;
  _List* _Listtest;
  _HYTopologyTraversalFunction* _HYTopologyTraversalFunctiontest;
  node<long>* node<long>test;
  _String* _Stringtest;
  _Matrix* _Matrixtest;
  _AVLListX* _AVLListXtest;
  _SimpleList* _SimpleListtest;
  bool* booltest;
  _Parameter* _Parametertest;
  _PMathObj* _PMathObjtest;
};


TEST_F(_TreeTopologyTest, AVLRepresentationTest) {

  _PMathObj result_PMathObj = _TreeTopologytest->AVLRepresentation(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(_TreeTopologyTest, AddANodeTest) {

  _TreeTopologytest->AddANode(*_PMathObjtest);
  //EXPECT_EQ (_TreeTopologytest, 0);

}


TEST_F(_TreeTopologyTest, BranchCountTest) {

  _PMathObj result_PMathObj = _TreeTopologytest->BranchCount();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(_TreeTopologyTest, BranchLengthTest) {

  _PMathObj result_PMathObj = _TreeTopologytest->BranchLength(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(_TreeTopologyTest, BranchNameTest) {

  _PMathObj result_PMathObj = _TreeTopologytest->BranchName(*_PMathObjtest, *booltest, *_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(_TreeTopologyTest, CompareTest) {

  _FString* result_FString = _TreeTopologytest->Compare(*_PMathObjtest);
  //EXPECT_EQ (result_FString*, 0);

}


TEST_F(_TreeTopologyTest, CompareTreesTest) {

  _String result_String = _TreeTopologytest->CompareTrees(_TreeTopologytest);
  //EXPECT_EQ (result_String, 0);

}


TEST_F(_TreeTopologyTest, ComputeClusterTableTest) {

  _TreeTopologytest->ComputeClusterTable(*_SimpleListtest);
  //EXPECT_EQ (_TreeTopologytest, 0);

}


TEST_F(_TreeTopologyTest, ConvertFromPSWTest) {

  _String* result_String = _TreeTopologytest->ConvertFromPSW(*_AVLListXtest);
  //EXPECT_EQ (result_String*, 0);

}


TEST_F(_TreeTopologyTest, ConvertToPSWTest) {

  bool resultbool = _TreeTopologytest->ConvertToPSW(*_AVLListXtest, _Listtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_TreeTopologyTest, CopyTreeStructureTest) {

  node<long>* resultnode<long> = _TreeTopologytest->CopyTreeStructure(node<long>test, *booltest);
  //EXPECT_EQ (resultnode<long>*, 0);

}


TEST_F(_TreeTopologyTest, DepthWiseTTest) {

  _TreeTopologytest->DepthWiseT(*booltest, _HYTopologyTraversalFunctiontest);
  //EXPECT_EQ (_TreeTopologytest, 0);

}


TEST_F(_TreeTopologyTest, DepthWiseTLevelTest) {

  _TreeTopologytest->DepthWiseTLevel(*longtest, *booltest);
  //EXPECT_EQ (_TreeTopologytest, 0);

}


TEST_F(_TreeTopologyTest, DepthWiseTRightTest) {

  _TreeTopologytest->DepthWiseTRight(*booltest);
  //EXPECT_EQ (_TreeTopologytest, 0);

}


TEST_F(_TreeTopologyTest, DetermineBranchLengthMappingModeTest) {

  _String result_String = _TreeTopologytest->DetermineBranchLengthMappingMode(_Stringtest);
  //EXPECT_EQ (result_String, 0);

}


TEST_F(_TreeTopologyTest, EdgeCountTest) {

  _TreeTopologytest->EdgeCount(*longtest, *longtest);
  //EXPECT_EQ (_TreeTopologytest, 0);

}


TEST_F(_TreeTopologyTest, ExecuteTest) {

  _PMathObj result_PMathObj = _TreeTopologytest->Execute(*longtest, *_PMathObjtest, *_PMathObjtest, _hyExecutionContexttest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(_TreeTopologyTest, FinalizeNodeTest) {

  bool resultbool = _TreeTopologytest->FinalizeNode(node<long>test, *longtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_TreeTopologyTest, FindCOTTest) {

  _AssociativeList* result_AssociativeList = _TreeTopologytest->FindCOT(*_PMathObjtest);
  //EXPECT_EQ (result_AssociativeList*, 0);

}


TEST_F(_TreeTopologyTest, FindCOTHelperTest) {

  _TreeTopologytest->FindCOTHelper(node<long>test, *longtest);
  //EXPECT_EQ (_TreeTopologytest, 0);

}


TEST_F(_TreeTopologyTest, FindCOTHelper2Test) {

  _TreeTopologytest->FindCOTHelper2(node<long>test, *_Matrixtest);
  //EXPECT_EQ (_TreeTopologytest, 0);

}


TEST_F(_TreeTopologyTest, FindNodeByNameTest) {

  node<long>* resultnode<long> = _TreeTopologytest->FindNodeByName(_Stringtest);
  //EXPECT_EQ (resultnode<long>*, 0);

}


TEST_F(_TreeTopologyTest, FlatRepresentationTest) {

  _PMathObj result_PMathObj = _TreeTopologytest->FlatRepresentation();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(_TreeTopologyTest, GetBranchLengthTest) {

  _TreeTopologytest->GetBranchLength(node<long>test, *_Parametertest);
  //EXPECT_EQ (_TreeTopologytest, 0);

}


TEST_F(_TreeTopologyTest, GetBranchLength1Test) {

  _TreeTopologytest->GetBranchLength(node<long>test, *_Stringtest, *booltest);
  //EXPECT_EQ (_TreeTopologytest, 0);

}


TEST_F(_TreeTopologyTest, GetBranchValueTest) {

  _TreeTopologytest->GetBranchValue(node<long>test, *_Stringtest);
  //EXPECT_EQ (_TreeTopologytest, 0);

}


TEST_F(_TreeTopologyTest, GetBranchVarValueTest) {

  _TreeTopologytest->GetBranchVarValue(node<long>test, *_Stringtest, *longtest);
  //EXPECT_EQ (_TreeTopologytest, 0);

}


TEST_F(_TreeTopologyTest, GetNodeNameTest) {

  _TreeTopologytest->GetNodeName(node<long>test, *_Stringtest, *booltest);
  //EXPECT_EQ (_TreeTopologytest, 0);

}


TEST_F(_TreeTopologyTest, IsCurrentNodeATipTest) {

  bool resultbool = _TreeTopologytest->IsCurrentNodeATip();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_TreeTopologyTest, IsCurrentNodeTheRootTest) {

  bool resultbool = _TreeTopologytest->IsCurrentNodeTheRoot();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_TreeTopologyTest, IsDegenerateTest) {

  bool resultbool = _TreeTopologytest->IsDegenerate();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_TreeTopologyTest, LeafWiseTTest) {

  _TreeTopologytest->LeafWiseT(*booltest);
  //EXPECT_EQ (_TreeTopologytest, 0);

}


TEST_F(_TreeTopologyTest, MainTreeConstructorTest) {

  bool resultbool = _TreeTopologytest->MainTreeConstructor(*_Stringtest, *booltest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_TreeTopologyTest, MatchTreePatternTest) {

  _String result_String = _TreeTopologytest->MatchTreePattern(_TreeTopologytest);
  //EXPECT_EQ (result_String, 0);

}


TEST_F(_TreeTopologyTest, PasteBranchLengthTest) {

  _TreeTopologytest->PasteBranchLength(node<long>test, *_Stringtest);
  //EXPECT_EQ (_TreeTopologytest, 0);

}


TEST_F(_TreeTopologyTest, PostTreeConstructorTest) {

  _TreeTopologytest->PostTreeConstructor(*booltest);
  //EXPECT_EQ (_TreeTopologytest, 0);

}


TEST_F(_TreeTopologyTest, PreTreeConstructorTest) {

  _TreeTopologytest->PreTreeConstructor(*booltest);
  //EXPECT_EQ (_TreeTopologytest, 0);

}


TEST_F(_TreeTopologyTest, RerootTreeTest) {

  _PMathObj result_PMathObj = _TreeTopologytest->RerootTree(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(_TreeTopologyTest, RerootTreeInternalTraverserTest) {

  _TreeTopologytest->RerootTreeInternalTraverser(*longtest);
  //EXPECT_EQ (_TreeTopologytest, 0);

}


TEST_F(_TreeTopologyTest, SetLeafNameTest) {

  _TreeTopologytest->SetLeafName(*longtest, _Stringtest);
  //EXPECT_EQ (_TreeTopologytest, 0);

}


TEST_F(_TreeTopologyTest, SplitTreeIntoClustersTest) {

  _List* result_List = _TreeTopologytest->SplitTreeIntoClusters(*longtest);
  //EXPECT_EQ (result_List*, 0);

}


TEST_F(_TreeTopologyTest, SplitTreeIntoClustersIntTest) {

  _List* result_List = _TreeTopologytest->SplitTreeIntoClustersInt(node<long>test, _Listtest);
  //EXPECT_EQ (result_List*, 0);

}


TEST_F(_TreeTopologyTest, SplitsIdentityTest) {

  _AssociativeList* result_AssociativeList = _TreeTopologytest->SplitsIdentity(*_PMathObjtest);
  //EXPECT_EQ (result_AssociativeList*, 0);

}


TEST_F(_TreeTopologyTest, StepWiseTTest) {

  _TreeTopologytest->StepWiseT(*booltest, _HYTopologyTraversalFunctiontest);
  //EXPECT_EQ (_TreeTopologytest, 0);

}


TEST_F(_TreeTopologyTest, StepWiseTLevelTest) {

  _TreeTopologytest->StepWiseTLevel(*longtest, *booltest);
  //EXPECT_EQ (_TreeTopologytest, 0);

}


TEST_F(_TreeTopologyTest, SubTreeStringTest) {

  _TreeTopologytest->SubTreeString(*_Stringtest, *booltest);
  //EXPECT_EQ (_TreeTopologytest, 0);

}


TEST_F(_TreeTopologyTest, TipCountTest) {

  _PMathObj result_PMathObj = _TreeTopologytest->TipCount();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(_TreeTopologyTest, TipNameTest) {

  _PMathObj result_PMathObj = _TreeTopologytest->TipName(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(_TreeTopologyTest, destroyCompTreeTest) {

  _TreeTopologytest->destroyCompTree(node<long>test);
  //EXPECT_EQ (_TreeTopologytest, 0);

}


TEST_F(_TreeTopologyTest, internalNodeCompareTest) {

  char resultchar = _TreeTopologytest->internalNodeCompare(node<long>test, node<long>test);
  //EXPECT_EQ (resultchar, 0);

}


TEST_F(_TreeTopologyTest, internalTreeCompareTest) {

  char resultchar = _TreeTopologytest->internalTreeCompare(node<long>test, node<long>test);
  //EXPECT_EQ (resultchar, 0);

}


TEST_F(_TreeTopologyTest, makeDynamicTest) {

  BaseRef resultBaseRef = _TreeTopologytest->makeDynamic();
  //EXPECT_EQ (resultBaseRef, 0);

}


TEST_F(_TreeTopologyTest, prepTree4ComparisonTest) {

  node<long>* resultnode<long> = _TreeTopologytest->prepTree4Comparison(*_Listtest);
  //EXPECT_EQ (resultnode<long>*, 0);

}


TEST_F(_TreeTopologyTest, toFileStrTest) {

  _TreeTopologytest->toFileStr(FILEtest);
  //EXPECT_EQ (_TreeTopologytest, 0);

}


TEST_F(_TreeTopologyTest, toStrTest) {

  BaseRef resultBaseRef = _TreeTopologytest->toStr();
  //EXPECT_EQ (resultBaseRef, 0);

}


}
