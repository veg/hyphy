ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "LUDecompose";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Returns the LU decomposition of a matrix using Crout's algorithm with partial pivoting.
  // The return object is an nx(n+1) matrix which contains the LU decomposition followed by a vector of row interchanges.

  simpleExample = LUDecompose({{1.2,5.6,7}{3, 4,12}{12,3.23,8}});
  expectedSimpleExample = {{12, 3.23,8,2}{0.1,5.277,6.2,2}{0.25,0.604983892363,6.24909986735,2}};
  for (i = 0; i < 3; i += 1) {
    for (j = 0; j < 4; j += 1) {
      assert(Abs(simpleExample[i][j] - expectedSimpleExample[i][j]) < 1e-6, "Failed to correctly LUDecompose a 3x 3 matrix");

    }
  }


  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  list1 = {"key1": "val1"};
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);

 
  assert (runCommandWithSoftErrors ('LUDecompose (1)', "Operation 'LUDecompose' is not implemented/defined for a Number"), "Failed error checking for trying to take LUDecompose of a number");
  assert (runCommandWithSoftErrors ('LUDecompose(list1)', "Operation 'LUDecompose' is not implemented/defined for a AssociativeList"), "Failed error checking for trying to take LUDecompose of associative list");
  assert (runCommandWithSoftErrors ('LUDecompose(T)', "Operation 'LUDecompose' is not implemented/defined for a Topology"), "Failed error checking for trying to take LUDecompose of Topology");
  assert (runCommandWithSoftErrors ('LUDecompose(TT)', "Operation 'LUDecompose' is not implemented/defined for a Tree"), "Failed error checking for trying to take LUDecompose of Tree");

  testResult = 1;

  return testResult;
}

