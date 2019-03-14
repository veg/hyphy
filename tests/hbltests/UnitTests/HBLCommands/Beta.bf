ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "Beta";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = FALSE;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Calculate the beta function of two numbers.
  // TODO: (this was passing before ~1/20/2019) assert(Beta(1,1) == 1, "Failed to return 1 for beta(1,1)");
  
  for (k = 0; k < 100; k+=1) {
     x = Random (1,100);
     y = Random (1,100);
     assert(Beta(x,y) == Beta(y,x), "Beta(x,y) = Beta (y,x) check failed"); 
  }
  

  // beta(2,2) acording to wolfram-alpha is 0.166...
  // beta (2,2) = Gamma (2) * Gamma (2) / Gamma (4) = 1! * 1! / 3! = 1/6

  beta_2_2 = 1/6;
  assert (Abs(Beta(2,2)-beta_2_2) < 1e-6, "Failed to return the expected value of beta(2,2) within float percision");

  // Comparing to none (Beta treats none as a zero and therefore always returns 1 if either argument is none)
  // TODO: (the three tests below were passing before ~1/20/2019)
  

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  Topology T2 = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);
  list1 = {"key1": "val1"};
  list2 = {"key2": "val2"};
  matrix1 = {{1,2}{3,4}};
  matrix2 = {{5,6}{7,8}};

  assert (runCommandWithSoftErrors ('Beta(T2,T2)', "Operation 'Beta' is not implemented/defined for a Topology"), "Failed error checking for trying to take Beta of topologies");
  assert (runCommandWithSoftErrors ('Beta(TT,TT)', "Operation 'Beta' is not implemented/defined for a Tree"), "Failed error checking for trying to take Beta of trees");
  assert (runCommandWithSoftErrors ('Beta(list1,list2)', "Operation 'Beta' is not implemented/defined for a AssociativeList"), "Failed error checking for trying to take Beta of lists");
  assert (runCommandWithSoftErrors ('Beta(matrix1, matrix2)', "Operation 'Beta' is not implemented/defined for a Matrix"), "Failed error checking for trying to take Beta of matrices");
  
  testResult = TRUE;

  return testResult;
}
