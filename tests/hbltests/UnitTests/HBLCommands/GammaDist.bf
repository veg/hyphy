ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "GammaDist";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Computes a point from the Gamma Distribution. GammaDist(x,k,Theta) k=shape, Theta=scale

  // TODO: doesn't match expected results.
  simpleExample = GammaDist(2,1,2);
  expectedSimpleExampleFromDocs = 0.366313;
  // Couldn't find a standard way to calculate gammadist but I'm using these sources:
  //  https://en.wikipedia.org/wiki/Gamma_distribution
  //  https://www.wolframalpha.com/input/?i=gamma+distribution(1,2)
  //  https://getcalc.com/statistics-gamma-distribution-calculator.htm?shape=1&scale=2&x=2#workout
  expectedSimpleExampleFromOtherSources = 0.18;
  

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  list1 = {"key1": "val1"};
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);

  
  // TODO: The three lines below should be correct but I'm getting 'Unconsumed values on the stack'.
  assert (runCommandWithSoftErrors ('GammaDist(list1)', "GammaDist' is not implemented/defined for a AssociativeList"), "Failed error checking for trying to take GammaDist of associative list");
  assert (runCommandWithSoftErrors ('GammaDist(T)', "Operation 'GammaDist' is not implemented/defined for a Topology"), "Failed error checking for trying to take GammaDist of Topology");
  assert (runCommandWithSoftErrors ('GammaDist(TT,1)', "Operation 'GammaDist' is not implemented/defined for a Tree"), "Failed error checking for trying to take GammaDist of Tree");

  testResult = 1;

  return testResult;
}

