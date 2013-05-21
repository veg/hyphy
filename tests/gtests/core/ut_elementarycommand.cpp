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
#include "alignment.h"
#include "batchlan.h"
#include "batchlan_globals.h"
#include "baseobj.h"
#include "bayesgraph.h"
#include "elementarycommand.h"
#include "executionlist.h"
#include "hy_globals.h"
#include "HYNetInterface.h"
#include "thetree.h"
#include "datasetfilter.h"
#include "likefunc.h"
#include "ctype.h"
#include "scfg.h"
#include "customfunction.h"
#include "thetree.h"
#include "sqlite3.h"

namespace {

// The fixture for testing class Foo.
class _ElementaryCommandTest : public ::testing::Test {

protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  _ElementaryCommandTest() {
    // You can do set-up work for each test here.
    // Create objects of every type needed. Performance doesn't matter.

    FILEtest = fopen ("./tests/gtests/res/HIV_gp120.nex" , "r");

    _ElementaryCommandtest = new _ElementaryCommand(*_Stringtest);
    longtest = new long();
    _ExecutionListtest = new _ExecutionList();
    _Listtest = new _List();
    _Stringtest = new _String(FILEtest);
    BaseReftest = new BaseRef();
    booltest = new bool();
  }

  virtual ~_ElementaryCommandTest() {
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
    delete _ElementaryCommandtest;
    delete longtest;
    delete _ExecutionListtest;
    delete _Listtest;
    delete _Stringtest;
    delete BaseReftest;
    delete booltest;
    fclose (FILEtest);
  }

  FILE* FILEtest;
  _ElementaryCommand* _ElementaryCommandtest;
  long* longtest;
  _ExecutionList* _ExecutionListtest;
  _List* _Listtest;
  _String* _Stringtest;
  BaseRef* BaseReftest;
  bool* booltest;
};


TEST_F(_ElementaryCommandTest, BuildDoWhileTest) {

  bool resultbool = _ElementaryCommandtest->BuildDoWhile(*_Stringtest, *_ExecutionListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, BuildForTest) {

  //bool resultbool = _ElementaryCommandtest->BuildFor(*_Stringtest, *_ExecutionListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, BuildIfThenElseTest) {

  //bool resultbool = _ElementaryCommandtest->BuildIfThenElse(*_Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, BuildWhileTest) {

  //bool resultbool = _ElementaryCommandtest->BuildWhile(*_Stringtest, *_ExecutionListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ConstructAlignSequencesTest) {

  //bool resultbool = _ElementaryCommandtest->ConstructAlignSequences(*_Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ConstructBGMTest) {

  bool resultbool = _ElementaryCommandtest->ConstructBGM(*_Stringtest, *_ExecutionListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ConstructCategoryTest) {

  //bool resultbool = _ElementaryCommandtest->ConstructCategory(*_Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ConstructCategoryMatrixTest) {

  //bool resultbool = _ElementaryCommandtest->ConstructCategoryMatrix(*_Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ConstructChoiceListTest) {

  //bool resultbool = _ElementaryCommandtest->ConstructChoiceList(*_Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ConstructDataSetTest) {

  //bool resultbool = _ElementaryCommandtest->ConstructDataSet(*_Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ConstructDataSetFilterTest) {

  //bool resultbool = _ElementaryCommandtest->ConstructDataSetFilter(*_Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ConstructDoSQLTest) {

  bool resultbool = _ElementaryCommandtest->ConstructDoSQL(*_Stringtest, *_ExecutionListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ConstructExecuteCommandsTest) {

  //bool resultbool = _ElementaryCommandtest->ConstructExecuteCommands(*_Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ConstructFindRootTest) {

  //bool resultbool = _ElementaryCommandtest->ConstructFindRoot(*_Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ConstructFscanfTest) {

  //bool resultbool = _ElementaryCommandtest->ConstructFscanf(*_Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ConstructFunctionTest) {

  //bool resultbool = _ElementaryCommandtest->ConstructFunction(*_Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ConstructGetDataInfoTest) {

  //bool resultbool = _ElementaryCommandtest->ConstructGetDataInfo(*_Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ConstructGetInformationTest) {

  //bool resultbool = _ElementaryCommandtest->ConstructGetInformation(*_Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ConstructGetNeutralNullTest) {

  //bool resultbool = _ElementaryCommandtest->ConstructGetNeutralNull(*_Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ConstructLFTest) {

  bool resultbool = _ElementaryCommandtest->ConstructLF(*_Stringtest, *_ExecutionListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ConstructMPIReceiveTest) {

  //bool resultbool = _ElementaryCommandtest->ConstructMPIReceive(*_Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ConstructMPISendTest) {

  //bool resultbool = _ElementaryCommandtest->ConstructMPISend(*_Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ConstructModelTest) {

  bool resultbool = _ElementaryCommandtest->ConstructModel(*_Stringtest, *_ExecutionListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ConstructNNTest) {

  bool resultbool = _ElementaryCommandtest->ConstructNN(*_Stringtest, *_ExecutionListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ConstructOpenDataPanelTest) {

  //bool resultbool = _ElementaryCommandtest->ConstructOpenDataPanel(*_Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ConstructOpenWindowTest) {

  //bool resultbool = _ElementaryCommandtest->ConstructOpenWindow(*_Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ConstructProfileStatementTest) {

  //bool resultbool = _ElementaryCommandtest->ConstructProfileStatement(*_Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ConstructReplicateConstraintTest) {

  //bool resultbool = _ElementaryCommandtest->ConstructReplicateConstraint(*_Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ConstructReturnTest) {

  //bool resultbool = _ElementaryCommandtest->ConstructReturn(*_Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ConstructSCFGTest) {

  bool resultbool = _ElementaryCommandtest->ConstructSCFG(*_Stringtest, *_ExecutionListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ConstructSpawnLFTest) {

  //bool resultbool = _ElementaryCommandtest->ConstructSpawnLF(*_Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ConstructStateCounterTest) {

  //bool resultbool = _ElementaryCommandtest->ConstructStateCounter(*_Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ConstructTreeTest) {

  bool resultbool = _ElementaryCommandtest->ConstructTree(*_Stringtest, *_ExecutionListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, DuplicateTest) {

  _ElementaryCommandtest->Duplicate(*BaseReftest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteTest) {

  //bool resultbool = _ElementaryCommandtest->Execute();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase0Test) {

  _ElementaryCommandtest->ExecuteCase0(*_ExecutionListtest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase11Test) {

  //_ElementaryCommandtest->ExecuteCase11();
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase12Test) {

  _ElementaryCommandtest->ExecuteCase12(*_ExecutionListtest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase21Test) {

  _ElementaryCommandtest->ExecuteCase21(*_ExecutionListtest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase25Test) {

  _ElementaryCommandtest->ExecuteCase25(*_ExecutionListtest, *booltest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase26Test) {

  _ElementaryCommandtest->ExecuteCase26(*_ExecutionListtest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase31Test) {

  _ElementaryCommandtest->ExecuteCase31(*_ExecutionListtest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase32Test) {

  _ElementaryCommandtest->ExecuteCase32(*_ExecutionListtest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase36Test) {

  _ElementaryCommandtest->ExecuteCase36(*_ExecutionListtest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase37Test) {

  _ElementaryCommandtest->ExecuteCase37(*_ExecutionListtest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase38Test) {

  _ElementaryCommandtest->ExecuteCase38(*_ExecutionListtest, *booltest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase39Test) {

  _ElementaryCommandtest->ExecuteCase39(*_ExecutionListtest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase4Test) {

  _ElementaryCommandtest->ExecuteCase4(*_ExecutionListtest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase40Test) {

  _ElementaryCommandtest->ExecuteCase40(*_ExecutionListtest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase41Test) {

  _ElementaryCommandtest->ExecuteCase41(*_ExecutionListtest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase43Test) {

  _ElementaryCommandtest->ExecuteCase43(*_ExecutionListtest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase44Test) {

  _ElementaryCommandtest->ExecuteCase44(*_ExecutionListtest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase45Test) {

  _ElementaryCommandtest->ExecuteCase45(*_ExecutionListtest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase46Test) {

  _ElementaryCommandtest->ExecuteCase46(*_ExecutionListtest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase47Test) {

  _ElementaryCommandtest->ExecuteCase47(*_ExecutionListtest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase5Test) {

  _ElementaryCommandtest->ExecuteCase5(*_ExecutionListtest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase52Test) {

  _ElementaryCommandtest->ExecuteCase52(*_ExecutionListtest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase53Test) {

  _ElementaryCommandtest->ExecuteCase53(*_ExecutionListtest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase54Test) {

  _ElementaryCommandtest->ExecuteCase54(*_ExecutionListtest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase55Test) {

  _ElementaryCommandtest->ExecuteCase55(*_ExecutionListtest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase57Test) {

  _ElementaryCommandtest->ExecuteCase57(*_ExecutionListtest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase58Test) {

  _ElementaryCommandtest->ExecuteCase58(*_ExecutionListtest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase61Test) {

  _ElementaryCommandtest->ExecuteCase61(*_ExecutionListtest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase63Test) {

  _ElementaryCommandtest->ExecuteCase63(*_ExecutionListtest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteCase64Test) {

  _ElementaryCommandtest->ExecuteCase64(*_ExecutionListtest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExecuteDataFilterCasesTest) {

  _ElementaryCommandtest->ExecuteDataFilterCases(*_ExecutionListtest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, ExtractConditionsTest) {

  //long resultlong = _ElementaryCommandtest->ExtractConditions(*_Stringtest, *longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_ElementaryCommandTest, ExtractValidateAddHBLCommandTest) {

  //bool resultbool = _ElementaryCommandtest->ExtractValidateAddHBLCommand();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, FindNextCommandTest) {

  _String result_String = _ElementaryCommandtest->FindNextCommand(*_Stringtest, *booltest);
  //EXPECT_EQ (result_String, 0);

}


TEST_F(_ElementaryCommandTest, HandleAssertTest) {

  bool resultbool = _ElementaryCommandtest->HandleAssert(*_ExecutionListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, HandleClearConstraintsTest) {

  bool resultbool = _ElementaryCommandtest->HandleClearConstraints(*_ExecutionListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, HandleComputeLFFunctionTest) {

  bool resultbool = _ElementaryCommandtest->HandleComputeLFFunction(*_ExecutionListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, HandleDeleteObjectTest) {

  bool resultbool = _ElementaryCommandtest->HandleDeleteObject(*_ExecutionListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, HandleDifferentiateTest) {

  bool resultbool = _ElementaryCommandtest->HandleDifferentiate(*_ExecutionListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, HandleExportTest) {

  bool resultbool = _ElementaryCommandtest->HandleExport(*_ExecutionListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, HandleFprintfTest) {

  bool resultbool = _ElementaryCommandtest->HandleFprintf(*_ExecutionListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, HandleGetStringTest) {

  bool resultbool = _ElementaryCommandtest->HandleGetString(*_ExecutionListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, HandleGetURLTest) {

  bool resultbool = _ElementaryCommandtest->HandleGetURL(*_ExecutionListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, HandleMolecularClockTest) {

  bool resultbool = _ElementaryCommandtest->HandleMolecularClock(*_ExecutionListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, HandleOptimizeCovarianceMatrixTest) {

  //bool resultbool = _ElementaryCommandtest->HandleOptimizeCovarianceMatrix();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, HandleRequireVersionTest) {

  bool resultbool = _ElementaryCommandtest->HandleRequireVersion(*_ExecutionListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, HandleSetParameterTest) {

  bool resultbool = _ElementaryCommandtest->HandleSetParameter(*_ExecutionListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, HandleUseModelTest) {

  bool resultbool = _ElementaryCommandtest->HandleUseModel(*_ExecutionListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, MakeGeneralizedLoopTest) {

  //bool resultbool = _ElementaryCommandtest->MakeGeneralizedLoop(_Stringtest, _Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, MakeJumpCommandTest) {

  //bool resultbool = _ElementaryCommandtest->MakeJumpCommand(_Stringtest, *longtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, ProcessIncludeTest) {

  //bool resultbool = _ElementaryCommandtest->ProcessInclude(*_Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ElementaryCommandTest, addAndCleanTest) {

  //_ElementaryCommandtest->addAndClean(*_ExecutionListtest, _Listtest);
  //EXPECT_EQ (_ElementaryCommandtest, 0);

}


TEST_F(_ElementaryCommandTest, makeDynamicTest) {

  BaseRef resultBaseRef = _ElementaryCommandtest->makeDynamic();
  //EXPECT_EQ (resultBaseRef, 0);

}


TEST_F(_ElementaryCommandTest, toStrTest) {

  BaseRef resultBaseRef = _ElementaryCommandtest->toStr();
  //EXPECT_EQ (resultBaseRef, 0);

}


}
