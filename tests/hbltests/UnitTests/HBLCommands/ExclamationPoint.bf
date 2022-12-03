ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "!";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // For numbers `!` should return 1 if the number is zero or return 0 otherwise.
  assert(!0 == 1, "Failed to return 1 when evaluating !0");
  assert(!1 == 0, "Failed to return 0 when evaluating !1");
  assert(!0.0 == 1, "Failed to return 1 when evaluatinng !0.0"); 
  assert(!0.0001 == 0, "Failed to return 0 when evaluationg !0.0001");
  negativeOne = -1;
  assert(!negativeOne == 0, "Failed to return 0 when evaluating !-1");
  // For strings `!` should return 1 if the file exists or return 0 if it doesn't.
  validPath = PATH_TO_CURRENT_BF + "TestTools.ibf";
  invalidPath = "/thisFileDoesNotExist";
  assert(!validPath == 1, "Failed to return 1 when evaluating `!` on a valid path string (`validPath`)");
  assert(!invalidPath == 0, "Failed to return 0 when evaluation `!` on an invalid path string (/thisFileDoesNotExist)");
  // The array functionality isn't documented but intuitively it seems like it sould return an array with `!elm` for each element in the original array.
  fprintf (stdout, "!{{1,1}{0,0}}: ", !{{1,1}{0,0}}, "\n");
  // TODO: the test below fails.
  // assert(!{{1,1}{0,0}} == {{0,0}{1,1}}, "Failed to evaluate ! on a matrix intuitively... The desired functionality is unknown");


  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);

  assert (runCommandWithSoftErrors ('!None', "Attempting to operate on an undefined value"), "Failed error checking for trying to evaluate `!None`");
  assert (runCommandWithSoftErrors ('!T', "not implemented/defined for a Topology"), "Failed error checking for trying to take a maximum of a topology");
  assert (runCommandWithSoftErrors ('!TT', "not implemented/defined for a Tree"), "Failed error checking for trying to take a maximum of a tree");
 

  testResult = 1;

  return testResult;
}
