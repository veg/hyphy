ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "Call";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = TRUE;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // TODO: Not in docs. I'm not sure what the expected functionality should be.
  // Used in `/libv3/models/protein/REF.bf like this:
  //__rp = Call (__rate_variation[terms.rate_variation.distribution], __rate_variation[terms.rate_variation.options], namespace);

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  list1 = {'key1':'val1', 'key2':'val2'};
  matrix1 = {{0.3,0.4}};
  Topology T1 = ((1,4),(3,4),5);
  Tree TT1 = ((1,2),(3,4),5);
  
  assert (runCommandWithSoftErrors ('Call(1);', "Operation 'Call' is not implemented/defined for a Number"), "Failed error checking for executing `Call` on a number");
  assert (runCommandWithSoftErrors ('Call(list1);', "Operation 'Call' is not implemented/defined for a AssociativeList"), "Failed error checking for executing `Call` on an associative list");
  assert (runCommandWithSoftErrors ('Call(T1);', "Operation 'Call' is not implemented/defined for a Topology"), "Failed error checking for executing `Call` on a topology");
  assert (runCommandWithSoftErrors ('Call(TT1);', "Operation 'Call' is not implemented/defined for a Tree"), "Failed error checking for executing `Call` on a tree");
  
  // TODO: When `Call(matrix1);` is run on it's own it has the error listed below... when run in the runCommandWithSoftErrors it has this error, `'Call(matrix1)' evaluated with errors. in call to Call(matrix1);`.
  //assert (runCommandWithSoftErrors ('Call(matrix1);', "Operation 'Call' is not implemented/defined for a Matrix"), "Failed error checking for executing `Call` on a matrix");
  
  testResult = 1;

  return testResult;
}
