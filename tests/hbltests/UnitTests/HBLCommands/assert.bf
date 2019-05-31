ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "assert";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = TRUE;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Assert is used for unit testing. The expression must evaluate to a Numeric, and it will fail if it returns 0 (a numerical false). Non-numeric expressions will cause an execution terminating error.
  //    If the assert passes, it will return nothing. The batch file will gleefully move on and execute subsequent commands.
  //    If the assert fails, it will greet the respective batch file executor with a generic (if error_statement is missing) or a custom error message. It will state what is in the error_statement. The executor will then either choose to move on, or halt the script from processing.

  // Passes with any non zero numeric value
  assert(1, "assert didn't work as expected, it shouldn't fail with value of 1");
  assert(0.01, "assert didn't work as expected, it shouldn't fail with value of 0.01");
  assert(100, "assert didn't work as expected, it shouldn't fail with value of 100");
  assert(-1, "assert didn't work as expected, it shouldn't fail with value of -1");

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  list1 = {'key1':'val1', 'key2':'val2'};
  matrix1 = {{0.3,0.4}};
  Topology T1 = ((1,4),(3,4),5);
  Tree TT1 = ((1,2),(3,4),5);
  
  assert (runCommandWithSoftErrors ('assert(list1);', "did not evaluate to a Number in call to assert"), "Failed error checking for assert returning a list");
  assert (runCommandWithSoftErrors ('assert(matrix1);', "did not evaluate to a Number in call to assert"), "Failed error checking for assert returning a matrix");
  assert (runCommandWithSoftErrors ('assert(T1);', "did not evaluate to a Number in call to assert"), "Failed error checking for assert returning a topology");
  assert (runCommandWithSoftErrors ('assert(TT1);', "did not evaluate to a Number in call to assert"), "Failed error checking for assert returning a tree");
  
  testResult = 1;

  return testResult;
}
