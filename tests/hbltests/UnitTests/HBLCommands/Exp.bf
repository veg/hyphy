ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "Exp";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  assert(Exp(1) == 2.718281828459045, "Failed to compute exponential of a number");
  assert(Exp(0.0001) == 1.000100005000167, "Failed to compute exponential of a small number (0.0001)");
  assert(Exp(100) == 2.688117141816136e+43, "Failed to compute exponential of a large number (100)");

  // Find exponent of matrix.
  assert(Abs (Exp({{1,2}{2,1}}) - {{10.22670818217949, 9.858828741008054}{9.858828741008054, 10.22670818217949}}) < 1e-14, "Failed to compute exponential value of an array");

  // Exp function on string; should return the length of the Lempel Ziv Production History.
  assert(Exp("1001111011000010") == 6, "Failed to compute exponential (Lempel Ziv Production History) of a string");
  

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);

  assert (runCommandWithSoftErrors ('Exp (None)', "Attempting to operate on an undefined value"), "Failed error checking for trying to take a exponential with None");
  assert (runCommandWithSoftErrors ('Exp (T)', "not implemented/defined for a Topology"), "Failed error checking for trying to take a exponential of a topology");
  assert (runCommandWithSoftErrors ('Exp (TT)', "not implemented/defined for a Tree"), "Failed error checking for trying to take a exponential of a tree");

  // TODO: DOESNT PASS THE TOO MANY ARGUMENTS ERROR CHECK.
  //assert (runCommandWithSoftErrors ('Exp (3,1)',  "Unconsumed values on the stack"), "Failed too many arguments error check");
  //assert (runCommandWithSoftErrors ('Exp ()',  "evaluated with errors"), "Failed too few arguments");
  //assert (runCommandWithSoftErrors ('Exp ({{1,2,3}{3,2,1}})',  "not defined for non-square matrices"), "Failed too few arguments");
 

  testResult = 1;

  return testResult;
}
