ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "Branchcount";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = FALSE;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Count the number of internal branches in a topology or tree. Uses the ACCEPT_ROOTED_TREES environmental variable.
  Topology T = ((a,b),(c,d));
  Topology largeT = ((1,2),(3,4),5,(10,(15,3)));
  Tree TT = ((a,b),(c,d));

  ACCEPT_ROOTED_TREES = 0;
  unrooted = BranchCount(TT);
  assert(unrooted == 1, "Failed to calculate the number of internal branches in an unrooted tree " );
  assert(BranchCount(T) == BranchCount(TT), "Branch count of matching tree and topology are the same");
  assert(BranchCount(largeT) == 4, "Branch count of larger tree is correct");

  ACCEPT_ROOTED_TREES = 1;
  Tree TT = ((a,b),(c,d));
  rooted = BranchCount(TT);
  assert(rooted == 2, "Failed to calculate the number of internal branches in a rooted tree");
  assert(BranchCount(largeT) == 4, "Branch count of larger tree is correct");


  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  list1 = {"key1": "val1"};
  matrix1 = {{1,2}{3,4}};

  assert (runCommandWithSoftErrors ('BranchCount(none)', "Operation 'BranchCount' is not implemented/defined for a Number"), "Failed error checking for running BranchCount on none");
  assert (runCommandWithSoftErrors ('BranchCount(10)', "Operation 'BranchCount' is not implemented/defined for a Number"), "Failed error checking for running BranchCount on a number");
  assert (runCommandWithSoftErrors ('BranchCount(list1)', "Operation 'BranchCount' is not implemented/defined for a AssociativeList"), "Failed error checking for running BranchCount on an associative list");
  assert (runCommandWithSoftErrors ('BranchCount(matrix1)', "Operation 'BranchCount' is not implemented/defined for a Matrix"), "Failed error checking for running BranchCount on a matrix");


  testResult = TRUE;

  return testResult;
}
