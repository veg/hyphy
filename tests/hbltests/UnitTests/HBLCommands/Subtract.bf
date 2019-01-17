ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "-";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Subtract numbers 
  assert(2-1 == 1, "Failed to subtract two integers");
  assert(3.33333-1.00001 == 2.33332, "Failed to subtract two floats");
  assert(-100-1 == -101, "Failed to subtract an int from a negative int");
  // Subtract matrices
  assert({{1,2,3}} - {{3,2,1}} == {{-2,0,2}}, "Failed to subtract two 1d matrices");
  assert({{1,2}{3,4}} - {{1,2}{1,2}} == {{0,0}{2,2}}, "Failed to subtract two 2d matrices");
  // Subtract two associative arrays (second elemet is removed from first element)
  list1 = {"key1": "value1", "key2": "value2", "key3": "value3"};
  list2 = {"key2": "value2"};
  expectedSubtractedList = {"key1": "value1", "key3": "value3"};
  list1 - list2;


  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  Topology T2 = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);
  oneDMatrix = {{1,2}};
  twoDMatrix = {{3,4}{5,6}};

  assert (runCommandWithSoftErrors ('T2-T2', "An invalid argument"), "Failed error checking for trying to add two topologies (you can only add a node, i.e. associative array, to a topology)");
  assert (runCommandWithSoftErrors ('TT-TT', "An invalid argument"), "Failed error checking for trying to add trees");
  assert (runCommandWithSoftErrors ('oneDMatrix - twoDMatrix', "Incompatible dimensions when trying to add or subtract matrices"), "Failed error checking for adding matrices of different dimensions");
  assert (runCommandWithSoftErrors ('"string1" - "string2"', "Operation '-' is not implemented/defined for a String"), "x");
  

  testResult = 1;

  return testResult;
}
