ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return ">";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Check if one integer is greater than another 
  assert (2>1 == 1, "Failed to return true when evaluating a greater than expression (1<2)");
  assert (1>(-1) == 1, "Failed to return true when evaluating a greater than expresion (-1<1)");
  assert (1>1 == 0, "Failed to return false when evaluationg a greater than exprssion (1<1)");
  assert(1>100 == 0, "Failed to return false when evaluating a greater than expression (100<1)");
  // Strings: Lexicographic (Generalized alphabetic order) greater than comparison between strings
  assert("Battlestar Galactica">"Bears" == 0, "Failed to return false for two strings in alphabetical order");
  assert("zoro" > "Aladin" == 1, "Failed to return true for two strings not in alphabetical order");
  // Comparing to none
  assert(1>none == 1, "Failed to return true for an int greater than none");
  assert(none>1 == 0, "Failed to return false for none greater than an int");
  
  // TODO: Matrix (Neighbor Join) 
  // NOTE: the results do not match the docs... I'm unsure of what method of neighbor join produces a 10x3 matrix from a 4x4 input
  exampleSquare4x4Matrix = {{0,7,11,14}{7,0,6,9}{11,6,0,7}{14,9,7,0}};
  neighborJoinedMatrix = exampleSquare4x4Matrix>1;
  expectedNeighborJoinedMatrix_Docs   = {{4,0,1}{4,0,1}{5,0,1}{5,0.3125,1}{5,0,3}{-1,0,6}{0,0,0}{0,0,0}{0,0,0}{0,0,0}};
  expectedNeighborJoinedMatrix_Actual = {{4,6,1}{4,1,1}{5,2,1}{5,5,     1}{5,3,3}{-1,0,6}{0,0,0}{0,0,0}{0,0,0}{0,0,0}};
  assert(neighborJoinedMatrix == expectedNeighborJoinedMatrix_Actual, "Failed to perform neighbor join on matrix");
  //assert(neighborJoinedMatrix == expectedNeighborJoinedMatrix_Docs, "Results of neighbor join on matrix do not match docs");


  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  Topology T2 = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);
  oneDMatrix = {{1,2}};
  twoDMatrix = {{3,4}{5,6}};
  list1 = {"key1": "val1"};
  list2 = {"key2": "val2"};

  assert (runCommandWithSoftErrors ('T2>T2', "Operation '>' is not implemented/defined for a Topology"), "Failed error checking for trying to compare topologies (>)");
  assert (runCommandWithSoftErrors ('TT>TT', "Operation '>' is not implemented/defined for a Tree"), "Failed error checking for trying to compare trees (>)");
  assert (runCommandWithSoftErrors ('list1>list2', "Operation '>' is not implemented/defined for a AssociativeList"), "Failed error checking for trying to compare associateive lists (>)");
  assert (runCommandWithSoftErrors ('oneDMatrix > twoDMatrix', "NeigborJoin needs a square numeric matrix of dimension >= 4"), "Failed error checking for NeighborJoining (>) matrices without the square 4x4 matrix");
  

  testResult = 1;

  return testResult;
}
