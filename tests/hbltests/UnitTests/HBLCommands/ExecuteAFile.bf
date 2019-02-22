ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "ExecuteAFile";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = TRUE;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Executes a set of HBL commands from within another file

  ExecuteAFile("./../../data/ExecuteAFileSetXTo5.bf");
  assert(X==5, "Failed to set a variable from a file executed with ExecuteAFile");

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  list1 = {'key1':'val1', 'key2':'val2'};
  assert (runCommandWithSoftErrors ("ExecuteAFile(list1);", "did not evaluate to a String in call to ExecuteAFile"), "Failed error checking for trying to perfor ExecuteAFile with a non-string object");
  assert (runCommandWithSoftErrors ("ExecuteAFile('');", "Empty/missing source code string in call to ExecuteAFile"), "Failed error checking for passing an empty string into ExecuteCommands");
  
  testResult = 1;

  return testResult;
}
