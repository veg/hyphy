ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "LUSolve";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // LUSolve is a matrix built-in function whose purpose is to find x in the formula: Ax = b.
  // If parameters are LUSolve(Matrix lu, Matrix b), it will solve for x in: Ax = b. lu is LUDecompose(A).

  lu = LUDecompose({{1.2,5.6,7}{3,4,12}{12,3.23,8}});
  b = {{1}{2}{3}};
  simpleExample = LUSolve(lu,b);
  expectedSimpleExample = {{0.167947882078}{-0.0227434688339}{0.132260852425}};
  for (i = 0; i < 3; i += 1) {
    assert(Abs(simpleExample[i] - expectedSimpleExample[i]) < 1e-6, "Failed to correctly LUSolve a 3x 3 matrix");
  }


  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  list1 = {"key1": "val1"};
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);

  assert (runCommandWithSoftErrors ('LUSolve(list1, list1)', "Operation 'LUSolve' is not implemented/defined for a AssociativeList"), "Failed error checking for trying to take LUSolve of associative lists");
  assert (runCommandWithSoftErrors ('LUSolve(T, T)', "Operation 'LUSolve' is not implemented/defined for a Topology"), "Failed error checking for trying to take LUSolve of Topologies");
  assert (runCommandWithSoftErrors ('LUSolve(TT, TT)', "Operation 'LUSolve' is not implemented/defined for a Tree"), "Failed error checking for trying to take LUSolve of Trees");
  assert (runCommandWithSoftErrors ('LUSolve(1, 1)', "Operation 'LUSolve' is not implemented/defined for a Number"), "Failed error checking for trying to take LUSolve of numbers");

  testResult = 1;

  return testResult;
}

