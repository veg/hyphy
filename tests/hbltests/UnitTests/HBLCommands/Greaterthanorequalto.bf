ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return ">=";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Check if one integer is greater than or equal to another 
  assert (2>=1 == 1, "Failed to return true when evaluating a greater than or equal to expression (2>=1)");
  assert (1>=(-1) == 1, "Failed to return true when evaluating a greater than or equal to expresion (1>=-1)");
  assert (1>=1 == 1, "Failed to return false when evaluationg a greater than or equal to exprssion (1>=1)");
  assert(1>=100 == 0, "Failed to return false when evaluating a greater than expression (1>=100)");
  // Strings: Lexicographic (Generalized alphabetic order) greater than or equal comparison between strings
  assert("Bears">="Battlestar Galactica" == 1, "Failed to return true for two strings in reverse alphabetical order");
  assert("Aladin" >= "Zoro" == 0, "Failed to return false for two strings in alphabetical order");
  assert("abcde" >= "abcde" == 1, "Failed to return true for two identical strings");
  // Comparing to none
  assert(none>=1 == 0, "Failed to return false for none greater than or equal to int");
  assert(1>=none == 1, "Failed to return true for int greater than or equal to an none");
  assert(none>=-1 == 1, "test");
  assert(none>=none == 1, "Failed to return true for none less than or equal to none");

  // Matrices (Creates a tree from the parent matrix passed in)
  matrix1 = {{0,7,11}{7,0,6}{11,6,0}};
  constructedTree = matrix1 >= 3;
  expectedTree = {{4, 0, 0, 0, -1}{0, 0, 0, 0, 0} {0, 0, 0, 0, 0} {0, 0, 0, 0, 0}{0, 0, 0, 11, 0}{1, 1, 0, 0, 0} {0, 0, 0, 0, 0} {0, 0, 0, 0, 0}};
  assert(constructedTree == expectedTree, "Faild to construct a tree from a parent matrix");


  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  Topology T2 = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);
  list1 = {"key1": "val1"};
  list2 = {"key2": "val2"};

  assert (runCommandWithSoftErrors ('T2>=T2', "Operation '>=' is not implemented/defined for a Topology"), "Failed error checking for trying to compare topologies (>=)");
  assert (runCommandWithSoftErrors ('TT>=TT', "Operation '>=' is not implemented/defined for a Tree"), "Failed error checking for trying to compare trees (>=)");
  assert (runCommandWithSoftErrors ('list1>=list2', "Operation '>=' is not implemented/defined for a AssociativeList"), "Failed error checking for trying to compare associateive lists (>=)");
  

  testResult = 1;

  return testResult;
}
