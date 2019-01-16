ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "<=";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Check if one integer is less than or equal to another 
  assert (1<=2 == 1, "Failed to return true when evaluating a less than expression (1<=2)");
  assert (-1<=1 == 1, "Failed to return true when evaluating a less than expresion (-1<=1)");
  assert (1<=1 == 1, "Failed to return false when evaluationg a less than exprssion (1<=1)");
  assert(100<=1 == 0, "Failed to return false when evaluating a less than expression (100<=1)");
  // Strings: Lexicographic (Generalized alphabetic order) less than or equal comparison between strings
  assert("Battlestar Galactica"<="Bears" == 1, "Failed to return true for two strings in alphabetical order");
  assert("zoro" <= "Aladin" == 0, "Failed to return false for two strings not in alphabetical order");
  assert("abcde" <= "abcde" == 1, "Failed to return true for two identical strings");
  // Comparing to none
  assert(1<=none == 0, "Failed to return false for an int less than or equal to none");
  assert(none<=1 == 1, "Failed to return true for none less than or equal to an int");
  assert(none<=none == 1, "Failed to return true for none less than or equal to none");

  // TODO: Matrices (seems to implement k-means clustering based on the error message; undocumented; potentially unimplemented)
  matrix1 = {{1,2}{3,4}{5,6}};
  // Error message indicating k-means
  //y = matrix1<=3;


  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  Topology T2 = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);
  list1 = {"key1": "val1"};
  list2 = {"key2": "val2"};

  // TODO: Both lines below result in, "HYPHYMP(70363,0x7fff96910380) malloc: *** error for object 0x7fad8de03fb0: pointer being freed was not allocated"
  //assert (runCommandWithSoftErrors ('T2<=T2', "Operation '<' is not implemented/defined for a Topology"), "Failed error checking for trying to compare topologies (<)");
  //assert (runCommandWithSoftErrors ('TT<=TT', "Operation '<' is not implemented/defined for a Tree"), "Failed error checking for trying to compare trees (<)");

  assert (runCommandWithSoftErrors ('list1<=list2', "Operation '<=' is not implemented/defined for a AssociativeList"), "Failed error checking for trying to compare associateive lists (<=)");
  

  testResult = 1;

  return testResult;
}
