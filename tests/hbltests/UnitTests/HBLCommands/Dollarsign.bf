ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "$";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Integer division.
  assert(1$1 == 1, "Failed to correctly perform integer division on 1$1");
  assert(90000$4 == 22500, "Failed to corretly perform integer divion on 90000$4");
  assert(9.99$3.99 == 3, "Failed to correctly perform integer division on 9.99$3.99");
  assert(8.99$3 == 2, "Failed to correctly perform integer division on 8.99$3");
  assert(4$9000 == 0, "Failed to correctly perform integer division on 4$9000");
  assert(0$1 == 0, "Failed to correctly perform integer division on 0$1");
  assert(1$0 == 0, "Failed to correctly perform integer division on 1$0");
  assert(4$none == 0, "Failed to return 0 when performing integer division with a denominator of none");
  // Matrix element-wise multiplication.
  assert({{1,2,3}}${{2,3,4}} == {{2,6,12}}, "Failed to correctly perform element-wise multiplication");
  assert({{1,2,3}{2,2,2}}${{2,3,4}{2,4,6}} == {{2,6,12}{4,8,12}}, "Failed to correctly perform element-wise multiplication on a 2d matrix");
  // Reg-ex matching.
  assert("ATATA"$"ATAT" == {{0}{3}}, "Failed to correctly match regex to string containing the regex");
  assert("ATATA"$"CCC" == {{-1}{-1}}, "Failed to correctly match regex to string that doesn't contain the regex");
  assert("ATATA"$"(TA)T" =={{1}{3}{1}{2}}, "Failed to correctly match regex to string when string contained parenthetical expression");
  

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);

  assert (runCommandWithSoftErrors ('T$T', "2nd argument is call to"), "Failed error checking for trying integer divide/compare topologies ($)");
  assert (runCommandWithSoftErrors ('TT$TT', "2nd argument is call to"), "Failed error checking for trying integer divide/compare trees ($)");
  assert (runCommandWithSoftErrors ('"Test"$4', "2nd argument in call to string"), "Failed error checking for trying integer divide/compare string$number");
 

  testResult = 1;

  return testResult;
}
