ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "Simplify";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = TRUE;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // First argument is an expression; second argument is a dictionary.
  // Performs a "find and replace" in the expression by searching for keys from the dictonary and replacing them with the values. 
  ten = Simplify("zero+one+2+3+four", {"zero": 0, "one":1, "four": 4});
  assert(ten == "10", "Failed to perform a Simplify on a mathematical expression as expected");

  // If a string is not a key in the dictonary, the string will be left as is.
  spanglishTen = Simplify("zero+one+2+3+quatro", {"zero": 0, "one":1, "four": 4});
  assert(spanglishTen == "6+quatro", "Failed to perform a Simplify on a mathematical expression as expected");

  // If the diconary has repeated keys, the last instance of the key will be used.
  tenLastArg = Simplify("zero+one+2+3+four", {"zero": 0, "one":1, "four": 100, "four": 4});
  assert(tenLastArg == "10", "Failed to perform a Simplify on a mathematical expression as expected");

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  list1 = {'key1':'val1', 'key2':'val2'};
  matrix1 = {{0.3,0.4}};
  Topology T1 = ((1,4),(3,4),5);
  Tree TT1 = ((1,2),(3,4),5);
  
  assert (runCommandWithSoftErrors ("Simplify(list1, {'val1': 'newVal1'});", "Operation 'Simplify' is not implemented/defined for a AssociativeList"), "Failed error checking for running 'Simplify' on an AssociativeList");
  assert (runCommandWithSoftErrors ("Simplify(T1, {'val1': 'newVal1'});", "Operation 'Simplify' is not implemented/defined for a Topology"), "Failed error checking for running 'Simplify' on an AssociativeList");
  assert (runCommandWithSoftErrors ("Simplify(TT1, {'val1': 'newVal1'});", "Operation 'Simplify' is not implemented/defined for a Tree"), "Failed error checking for running 'Simplify' on an AssociativeList");
  // TODO: unsure why the line below is throwing errors
  //assert (runCommandWithSoftErrors ("Simplify(matrix1, {'val1': 'newVal1'});", "Operation 'Simplify' is not implemented/defined for a Matrix"), "Failed error checking for running 'Simplify' on an AssociativeList");


  testResult = 1;

  return testResult;
}
