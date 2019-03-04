ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "CGammaDist";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Calculate the incomplete gamma function (CGammaDist(a,b,c) --> iGamma(b,a*c)).
  assert((Abs(CGammaDist(4,20,5) - 0.529743)) < 1e-6, "Failed to return the expected value of CGammaDist(4,20,5)");
  
  // Comparing to none (CGammaDist treats none as a zero and therefore always returns 1 if either argument is none)
  assert(CGammaDist(none,5,1) == CGammaDist(0,5,1), "Failed to return treat 'none' as '0'");
  assert(CGammaDist(5,none,1) == CGammaDist(5,0,1), "Failed to return treat 'none' as '0'");
  assert(CGammaDist(5,5,none) == CGammaDist(5,5,0), "Failed to return treat 'none' as '0'");

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);
  list1 = {"key1": "val1"};
  matrix1 = {{1,2}{3,4}};

  assert (runCommandWithSoftErrors ('CGammaDist(T,1,1)', "Operation 'CGammaDist' is not implemented/defined for a Topology"), "Failed error checking for trying to take CGammaDist of topologies");
  assert (runCommandWithSoftErrors ('CGammaDist(TT,1,1)', "Operation 'CGammaDist' is not implemented/defined for a Tree"), "Failed error checking for trying to take CGammaDist of trees");
  assert (runCommandWithSoftErrors ('CGammaDist(list1,1,1)', "Operation 'CGammaDist' is not implemented/defined for a AssociativeList"), "Failed error checking for trying to take CGammaDist of lists");

  // TODO: Not sure why the assert below doesn't work.
  //assert (runCommandWithSoftErrors ('CGammaDist(matrix1,1,1)', "Operation 'CGammaDist' is not implemented/defined for a Matrix"), "Failed error checking for trying to take CGammaDist of matrices");
  

  testResult = 1;

  return testResult;
}
