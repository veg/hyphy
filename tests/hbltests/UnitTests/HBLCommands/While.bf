ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "While";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = TRUE;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // The general while loop construct. Very much like the C while. 
  // It executes the code block within the loop while the `while` condition is satisfied.

  x = 0;
  y = 0;

  i = 0;
  while (i<10) {
    x = x+1;
    j=0;
    while (j<10) {
      y = y+1;
      j=j+1;
    }
    i=i+1;
  }

  assert(x == 10, "A while loop behaved unexpectedly");
  assert(y == 100, "A nested while loop behaved unexpectedly");

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  assert (runCommandWithSoftErrors ('while(i=0; i<10) {i=i+1;};', "Incorrect number of arguments"), "Failed error checking for too many arguments in a while loop.");
  
  testResult = 1;

  return testResult;
}
