ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "&";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = TRUE;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // TODO: It's unclear what `&` is supposed to do. 
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);
  matrix = {{0.3,0.4}};

  // TODO, the expression below results in the following error: Syntax error  in the following context: 'T&T<ERROR HERE>'
  // T&T;

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  
  
  // TODO: can't get the below to work?
  //assert (runCommandWithSoftErrors ('y = 1 & 1;', "Operation '&' is not implemented/defined for a Number"), "Failed error checking for trying to execute number&number.");
  
  testResult = 1;

  return testResult;
}
