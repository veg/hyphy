ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "Erf";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Calculate the error function of a number.
  erf1 = Erf(0.75);
  expectedErf1 = 0.7111556336535151315989378345914107773742059540965372; // From wolfram alpha
  assert(Abs(erf1 - expectedErf1) < 1e-6, "Failed to compute error function of 0.75 within 6 digits");
  // TODO: minor issue below.
  //assert(Abs(erf1 - expectedErf1) < 1e-11, "Failed to compute error function of 0.75 within 11 digits");
  erf2 = Erf(-0.1);
  expectedErf2 = -0.112462916018284892203275071743968383221696299159702; // From wolfram alpha
  assert(Abs(erf2 - expectedErf2) < 1e-6, "Failed to compute error function of -0.1 within 6 digits");
  assert(Abs(erf2 - expectedErf2) < 1e-11, "Failed to compute error function of -0.1 within 11 digits");
  erf3 = Erf(5);
  expectedErf3 = 0.9999999999984625402055719651498; // From wolfram alpha
  assert(Abs(erf3 - expectedErf3) < 1e-6, "Failed to compute error function of 5 within 6 digits");
  assert(Abs(erf3 - expectedErf3) < 1e-11, "Failed to compute error function of 5 within 11 digits");
  assert(Erf(0) == 0, "Failed to compute error function of zero");
  assert(Erf(none) == 0, "Failed to compute error function of zero (none)");
  

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);
  list1 = {"key1": "value1", "key2": "value2"};

  assert (runCommandWithSoftErrors ('Erf (T)', "not implemented/defined for a Topology"), "Failed error checking for trying to take a Erfonential of a topology");
  assert (runCommandWithSoftErrors ('Erf (TT)', "not implemented/defined for a Tree"), "Failed error checking for trying to take a Erfonential of a tree");
  assert (runCommandWithSoftErrors ('Erf (list1)', "not implemented/defined for a AssociativeList"), "Failed error checking for trying to take a Erfonential of a AssociativeList");
  assert (runCommandWithSoftErrors ('Erf ("string")', "not implemented/defined for a String"), "Failed error checking for trying to take a Erfonential of a stirng");


  testResult = 1;

  return testResult;
}
