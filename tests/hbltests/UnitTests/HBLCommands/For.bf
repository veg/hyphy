ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "For";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = TRUE;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // The general for loop construct. Very much like the C for, initial_statement consists of one statement (typically assignment), 
  // condition is a logical statement (compound in general), and increment is the statement executed at the end of each loop iteration.

  x = 0;
  y = 0;

  for (i=0; i<10; i=i+1) {
    x = x+1;
    for (j=o; j<10; j=j+1) {
      y = y+1;
    }
  }

  assert(x == 10, "A for loop behaved unexpectedly");
  assert(y == 100, "A nested for loop behaved unexpectedly");

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  assert (runCommandWithSoftErrors ('for(i=0; i<10) {i=i+1;};', "Incorrect number of arguments"), "Failed error checking for too few arguments in a for loop.");
  
  testResult = 1;

  return testResult;
}
