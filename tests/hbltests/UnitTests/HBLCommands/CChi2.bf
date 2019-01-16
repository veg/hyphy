ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "CChi2";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Computes Fisher's exact test. If number is greater than 1, then the p-value will be computed with Cochran-Armitage test for trend.
  twoCategories = CChi2({{1,2}{3,0}},5);
  assert((Abs(tw0Categories) - 4) < 1e-6, "Failed to return the expected value for a two category chi squared test");

  // TODO: see below
  fprintf (stdout, 'there doesnt seem to be a difference between a number greater than or less than one (i.e. it looks like its always computing Cochran-Armitage).\n');
  threeCategoriesMatrixTrend = {{1,2,3}{10,21,29}};
  threeCategoriesMatrixNoTrend = {{1,3,2}{29,21,10}};
  threeCategoriesChiTrendCochran = CChi2(threeCategoriesMatrixTrend, 2);
  threeCategoriesChiTrendNonCochran = CChi2(threeCategoriesMatrixTrend, 1);
  threeCategoriesChiNoTrendCochran = CChi2(threeCategoriesMatrixNoTrend, 3);
  threeCategoriesChiNoTrendNonCochran = CChi2(threeCategoriesMatrixNoTrend, 0);
  fprintf (stdout, 'threeCategoriesChiTrendCochran: ', threeCategoriesChiTrendCochran, '\n');
  fprintf (stdout, 'threeCategoriesChiTrendNonCochran: ', threeCategoriesChiTrendNonCochran, '\n');
  fprintf (stdout, 'threeCategoriesChiNoTrendCochran: ', threeCategoriesChiNoTrendCochran, '\n');
  fprintf (stdout, 'threeCategoriesChiNoTrendNonCochran: ', threeCategoriesChiNoTrendNonCochran, '\n');
  // Confirm the matrix shape doesn't matter
  threeCategoriesFlipped = {{1,10}{10,21},{3,29}};
  fprintf (stdout, 'CChi2(threeCategoriesFlipped, 1): ', CChi2(threeCategoriesFlipped, 1), '\n');
  fprintf (stdout, 'CChi2(threeCategoriesFlipped, 2): ', CChi2(threeCategoriesFlipped, 2), '\n');

  //Computes Chi-square probability function. chsq is chi-square and df is the number of degrees of freedom.
  numericChiProbOneDeg = CChi2(1.44,1);
  assert((Abs(numericChiProbOneDeg) - 0.769861) < 1e-6, "Failed to return the expected value for a chi squared probability function with one degree of freedom");
  numericChiProbTwoDeg = CChi2(1.44,2);
  assert((Abs(numericChiProbTwoDeg) - 0.5132477) < 1e-6, "Failed to return the expected value for a chi squared probability function with two degrees of freedom");
  assert(CChi2(1.44, 0) == 0, "Failed to return the expected value for a chi squared probability function with zero degrees of freedom");
  assert(CChi2(1.44, none) == 0, "Failed to return the expected value for a chi squared probability function with zero (none) degrees of freedom");

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  list1 = {"key1": "val1"};
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);

  
  assert (runCommandWithSoftErrors ('CChi2(list1, 1)', "Operation 'CChi2' is not implemented/defined for a AssociativeList"), "Failed error checking for trying to take Chi squared of associative list");
  assert (runCommandWithSoftErrors ('CChi2(TT,1)', "Operation 'CChi2' is not implemented/defined for a Tree"), "Failed error checking for trying to take Chi squared of Tree");
  assert (runCommandWithSoftErrors ('CChi2(T,1)', "Operation 'CChi2' is not implemented/defined for a Topology"), "Failed error checking for trying to take Chi squared of Topology");

  testResult = 1;

  return testResult;
}
