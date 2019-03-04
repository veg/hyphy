ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "sscanf";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Not in docs but appears to take three arguments:
  //   1. An input stream
  //   2. The type of input stream ("Number", "Matrix", "Tree", "String", "NMatrix", "Raw", "Lines")
  //   3. A variable identifier to save the input stream
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);
  list1 = {"key1": "val1"};
  matrix1 = {{1,2}{3,4}};

  sscanf('3', "Number", three);
  assert(three == 3, "Failed to sscanf an integer");
  sscanf('3.01', "Number", threePointOhOne);
  assert(threePointOhOne == 3.01, "Failed to sscanf a float");
  
  sscanf('blahblah', "Raw", rawHolder);
  sscanf("blahblah", "String", stringHolder);
  assert(rawHolder == stringHolder, "Failed to sscanf a raw string as a string");
  
  sscanf("{{1,2}{3,4}}", "Matrix", matrixHolder);
  assert(matrixHolder == {{1,2}{3,4}}, "Failed to sscanf a matrix");

  sscanf("((1,2),(3,4),5)", "Tree", treeHolder);
  assert(BranchCount(treeHolder) == 2, "Failed to sscanf a tree (the resulting tree didn't have the expected number of branches)");

  // TODO Unsure of how NMatrix differs from Matrix.
  sscanf("{{1,2}{3,4}}", "NMatrix", matrixHolder);
  assert(matrixHolder == {{1,2}{3,4}}, "Failed to sscanf a matrix");

  // The "Lines" option isn't tested because: (1) It's unclear what the functionality should be and (2) It is only used in one file, "help/Commands/query.bf"

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  assert (runCommandWithSoftErrors ("sscanf('3', 'notValid', three);", "is not one of the supported argument types: Number, Matrix, Tree, String, NMatrix, Raw, Lines"), "Failed error checking for passing an invalid argument type (second argument) to sscanf");
  
  testResult = 1;

  return testResult;
}
