ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "do";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = TRUE;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // do...while is a control flow tool. It executes the code block within the loop at least once, and then continues while the while condition is satisfied.
  
  i = 0;
  do
  {
      i=i+1;
  }  
  while (i<=100);

  assert(i == 101, "`do` behaved unexpectedly");

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  assert (runCommandWithSoftErrors ('do{};;', "Could not find a matching 'while' in the definition of a do-while loop"), "Failed error checking for trying to execute `do` without a `while` definition.");
  
  testResult = 1;

  return testResult;
}
