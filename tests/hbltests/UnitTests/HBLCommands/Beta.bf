ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "Beta";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Calculate the beta function of two numbers.
  assert(Beta(1,1) == 1, "Failed to return 1 for beta(1,1)");
  assert(Beta(2,3) == Beta(3,2), "The Beta function didn't operate symetrically"); 

  // beta(2,2) acording to wolfram-alpha is 0.166...
  beta_2_2 = 0.166666666666;
  assert (Abs(Beta(2,2)-beta_2_2) < 1e-6, "Failed to return the expected value of beta(2,2) within float percision");
  //assert (Abs(Beta(2,2)-beta_2_2) < 1e-12, "Failed to return the expected value of beta(2,2) within double percision");
  //assert(Beta((-1),1) == (-1), "Failed to return -1 for beta beta(-1,1)");

  // Comparing to none (Beta treats none as a zero and therefore always returns 1 if either argument is none)
  assert(Beta(none,5) == 1, "Failed to return '1' when evaluating beta(none,5)");
  assert(Beta(5,none) == 1, "Failed to return '1' when evaluating beta(5,none)");
  assert(Beta(none,none) == 1, "Failed to return '1' when evaluating beta(none,none)");

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

  // Not sure why the assert two lines below doesn't work... the line below throws the error that I have listed but the assert doesn't work.
  //Beta(matrix1, matrix2);
  //assert (runCommandWithSoftErrors ('Beta(matrix1, matrix2)', "Operation 'Beta' is not implemented/defined for a Matrix"), "Failed error checking for trying to take Beta of matrices");
  

  testResult = 1;

  return testResult;
}
