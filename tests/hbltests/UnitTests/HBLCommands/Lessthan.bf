ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "<";
}


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;


  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Check if one integer is less than another
  assert (1<2 == 1, "Failed to return true when evaluating a less than expression (1<2)");
  assert (-1<1 == 1, "Failed to return true when evaluating a less than expresion (-1<1)");
  assert (1<1 == 0, "Failed to return false when evaluationg a less than exprssion (1<1)");
  assert(100<1 == 0, "Failed to return false when evaluating a less than expression (100<1)");
  // Matrices (the docs state: Returns the Path Log Likelihood assuming that matrix A is a 3xK matrix, where each column is of the form A: integer in [0,N-1], B: integer in [0,N-1], T: real >= 0, and matrix B is an NxN RATE matrix for a Markov chain.)
  matrix1 = {
             {0,1}
             {1,2}
             {5,6}
            };
  allOnesRateMatrix = {4,4}["0.1"]; // all equal rate matrix
  allOnesRateMatrix [0][0] = -0.3;
  allOnesRateMatrix [1][1] = -0.3;
  allOnesRateMatrix [2][2] = -0.3;
  allOnesRateMatrix [3][3] = -0.3;
  x = matrix1<allOnesRateMatrix;
  assert(x == -3.013102130631045, "Failed to return the path log likelihood");
  // Strings: Lexicographic (Generalized alphabetic order) comparison between strings
  assert("Battlestar Galactica"<"Bears" == 1, "Failed to return true for two strings in alphabetical order");
  assert("zoro" < "Aladin" == 0, "Failed to return false for two strings not in alphabetical order");
  // Comparing to none
  assert(1<none == 0, "Failed to return false for an int less than none");
  assert(none<1 == 1, "Failed to return true for none less than an int");


  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  Topology T2 = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);
  oneDMatrix = {{1,2}};
  twoDMatrix = {{3,4}{5,6}};
  list1 = {"key1": "val1"};
  list2 = {"key2": "val2"};

  assert (runCommandWithSoftErrors ('T2<T2', "Operation '<' is not implemented/defined for a Topology"), "Failed error checking for trying to compare topologies (<)");
  assert (runCommandWithSoftErrors ('TT<TT', "Operation '<' is not implemented/defined for a Tree"), "Failed error checking for trying to compare trees (<)");
  assert (runCommandWithSoftErrors ('list1<list2', "Operation '<' is not implemented/defined for a AssociativeList"), "Failed error checking for trying to compare associateive lists (<)");
  // TODO: Looks like an uncaught error for the line below: "libc++abi.dylib: terminating with uncaught exception of type char const*"
  //assert (runCommandWithSoftErrors ('oneDMatrix < twoDMatrix', "Incompatible dimensions when trying to add or subtract matrices"), "Failed error checking for addig matrices of different dimensions");


  testResult = 1;

  return testResult;
}
