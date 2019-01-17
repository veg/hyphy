ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "/";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Divide numbers 
  assert(6/3 == 2, "Failed to divide two integers that are evenly divisable (no remainder)");
  assert(18/4 == 4.5, "Failed to divide two integers that are not evenly divisable (the floating point value should be returned)"); 
  assert(18.654321/4.56789 == 4.083793830411853, "Failed to divide two floats");
  assert(-20/5 == -4, "Failed to divide a negative number by a positive number; expected the result to be negative");
  assert(20/(-5) == -4, "Failed to divide a positive number by a negative number; expected the result to be negative");
  assert(-20/(-5) ==4, "Failed to divide two negative numbers; expected the result to be positive");
  // Check if two strings are equal (the second argument can use wild card characger "*" as a substitue for any character)
  assert("string1"/"string1" == 1, "Failed to compare two strings that were the same");
  assert("string1"/"String1" == 0, "Failed to compare two strings taht were the same except for case");
  assert("string1"/"string2" == 0, "Failed to comare two string that were different");
  assert("string1"/"string*" == 1, "Failed to compare two strings with the second string having a character replaced with *");
  assert("string1"/"*trin**" == 1, "Failed to compare two strings with the second string having multiple characters replaced with * (each replaced character has it's own *)");
  assert("string1"/"*tri*" == 1, "Failed to compare two strings with the second string having a multipe characters replaced with * (more than one consecutive characters was repalced with just one *)");
  // Divide Matrices (element wise division) (NOT IN DOCS)
  assert({{4,2}{10,8}} / {{1,1}{5,2}} =={{4,2}{2,4}}, "Failed to perform element wise division on matrices");
  // TODO: Divide two associative lists (I'm not sure what this is supposed to do) (NOT IN DOCS)
  list1 = {"key": "val"};
  list2 = {"otherKey": "otherVal"};
  x = list1/list2;
  //fprintf (stdout, "x : ", x , "\n");

  
  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  Topology T2 = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);
  oneDMatrix = {{1,2}};
  twoDMatrix = {{3,4}{5,6}};


  assert (runCommandWithSoftErrors ('T2/T2', "Operation '/' is not implemented/defined for a Topology"), "Failed error checking for trying to add two topologies (you can only add a node, i.e. associative array, to a topology)");
  assert (runCommandWithSoftErrors ('TT/TT', "Operation '/' is not implemented/defined for a Tree"), "Failed error checking for trying to add trees");
  assert (runCommandWithSoftErrors ('oneDMatrix/twoDMatrix', "Element-wise multiplication/division requires matrixes of the same dimension"), "Failed error checking for dividing matrices of differnet dimensions");

  testResult = 1;

  return testResult;
}
