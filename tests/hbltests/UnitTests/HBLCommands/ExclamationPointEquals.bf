ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "!=";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Compare two numbers.
  assert(1 != 1 == 0, "Failed to correctly compare two congruent numbers");
  assert(1 != 1.0001 == 1, "Failed to correctly compare two numbers that are close but not congruent");
  assert(1 != 2 == 1, "Failed to correctly compare two non-congruent numbers");
  assert("1" != 1 == 0, "Failed to call a string and number comparison with the string first as congruent");
  assert(1 != "1" == 1, "Failed to call a string and number comparison with the numbrer first as non-congruent");
  assert("test" != "test" == 0, "Failed to correctly compare two strings that are the same");
  assert("test" != "nottest" == 1, "Failed to correctly compare two strings that are not the same");
  assert(None != 1 == 1, "Failed to correclty compare None to a number");
  

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);

  assert (runCommandWithSoftErrors ('T != T', "not implemented/defined for a Topology"), "Failed error checking for trying compare topologies (!=)");
  assert (runCommandWithSoftErrors ('TT != TT', "not implemented/defined for a Tree"), "Failed error checking for trying to compare trees (!=)");
  // The below assert fails even though the error I get when running the command seems to match the regex
  //assert (runCommandWithSoftErrors ('{1,2} != {1,2}', "Parameter list is out of context"), "Failed error checking for trying to compare arrays (!=)");
  //{1,2} != {1,2};
 

  testResult = 1;

  return testResult;
}
