ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "||"
}		


function runTest () {
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // BASIC FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  assert((0||0) == 0, "False or False should be False");
  assert((1||0) == 1, "True or False should be True");
  assert((0||1) == 1, "False or True should be True");
  assert((1||1) == 1, "True or True should be True");

  imprecise_expression = 1 >= 2 || 0;
  proper_operator_precedence = (1 >= 2) || 0;
  assert(imprecise_expression == proper_operator_precedence, "Improper handling of operator precedence");

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);
  matrix = {{0.3,0.4}};

  assert (runCommandWithSoftErrors ('T||T', "is not implemented/defined for a Topology"), "Failed error checking for trying to compare topologies (||)");
  assert (runCommandWithSoftErrors ('TT||TT', "is not implemented/defined for a Tree"), "Failed error checking for trying to compare trees (||)");
  assert (runCommandWithSoftErrors ('4||"String"', "where 'X' is not a number"), "Failed error checking for trying to compare number||string");
  // TODO: Fails error checking for matrices. Not sure why.
  //assert (runCommandWithSoftErrors ('matrix||matrix', "is not implemented/defined for a Matrix"), "Failed error checking for trying to compare matrix||matrix");

  testResult = 1;

  return testResult;
}
