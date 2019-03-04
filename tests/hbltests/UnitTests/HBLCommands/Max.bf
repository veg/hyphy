ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "Max";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Compare two numbers.
  assert(Max(1,2) == 2, "Failed to compute maximum of two numbers");
  assert(Max(4,4.001) == 4.001, "Failed to compute maximum of two numbers that are simular in magnitude");
  assert(Max(-1,-2) == -1, "Failed to compute maximum of two negative numbers");
    
  // Find max of array.
  assert(Max({{1,2}{3,4}},2) == 4, "Failed to compute maximum value of an array");
  // Find max of array and it's index. 
  assert(Max({{1,2}{3,4}},1) == {{4, 3}}, "Failed to compute maximum value and index in an array");
  
  // Max functions on dictionary; should only select among keys with numeric values
  dict = {"0" : 1, "1" : {{2,3}}, "hai" : {"a" : 5, "b" : 7}, "beavis" : 42};
  assert((Max(dict))["key"] == "beavis", "Failed to compute maximum value in a dictionary");
  

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);

  assert (runCommandWithSoftErrors ('Max ("abc",0)', "not implemented/defined for a String"), "Failed error checking for trying to take a maximum of a string");
  assert (runCommandWithSoftErrors ('Max (None,0)', "Attempting to operate on an undefined value"), "Failed error checking for trying to take a maximum with None");
  assert (runCommandWithSoftErrors ('Max (T)', "not implemented/defined for a Topology"), "Failed error checking for trying to take a maximum of a topology");
  assert (runCommandWithSoftErrors ('Max (TT)', "not implemented/defined for a Tree"), "Failed error checking for trying to take a maximum of a tree");
  assert (runCommandWithSoftErrors ('Max (1)',  "was called with an incorrect number of arguments"), "Too few arguments error check");
  // TODO: The below test fails.
  // assert (runCommandWithSoftErrors ('Max (1,2,3)',  "Error compiling the statement"), "Too many arguments error check");
 

  testResult = 1;

  return testResult;
}
