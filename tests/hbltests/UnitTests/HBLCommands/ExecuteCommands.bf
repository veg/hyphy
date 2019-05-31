ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "ExecuteCommands";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = TRUE;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Executes a set of HBL commands in string format

  x = 'testString';
  ExecuteCommands("y = 'testString';");
  assert(x==y, "failed to set a variable inside ExecuteCommands");

  a = 'testStringTwo';
  ExecuteCommands("tempStringOne = 'testString'; tempStringTwo = 'Two'; b = tempStringOne+tempStringTwo;");
  assert(a==b, "failed to execute three commands in series inside ExecuteCommands");

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  list1 = {'key1':'val1', 'key2':'val2'};
  assert (runCommandWithSoftErrors ("ExecuteCommands(list1);", "did not evaluate to a String in call to ExecuteCommands"), "Failed error checking for trying to perfor ExecuteCommands with a non-string object");
  assert (runCommandWithSoftErrors ("ExecuteCommands('');", "Empty/missing source code string in call to ExecuteCommands"), "Failed error checking for passing an empty string into ExecuteCommands");
  
  testResult = 1;

  return testResult;
}
