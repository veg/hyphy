ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "%";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Remainder from integer division (modulo).
  assert(29%5 == 4, "Failed to correctly perform modulo on 8.99%3");
  assert(1%1 == 0, "Failed to correctly perform modulo on 1%1");
  assert(90000%4 == 0, "Failed to corretly perform modulo on 90000%4");
  assert(9.99%3.99 == 0, "Failed to correctly perform modulo on 9.99%3.99");
  assert(0%1 == 0, "Failed to correctly perform modulo on 0%1");
  assert(0%1000 == 0, "Failed to correctly perform modulo on 0%1000");
  assert(1%0 == 1, "Failed to correctly perform modulo on 1%0");
  assert(1000%0 == 1000, "Failed to correctly perform modulo on 1000%0");
  assert(10%none == 10, "Failed to correctly perform modulo on 10%none");
  assert("string"%4 == 0, "Failed to correctly perform modulo on string%number");
  // Sort matrices.
  unsortedMatrix = {{1,2,3}{100,101,102}{0.1,0.8,0.2}{0.0001,0.00002,1000}};
  sortedMatrix_lastColumn = unsortedMatrix%2;
  assert(sortedMatrix_lastColumn == {{0.1, 0.8, 0.2}{1, 2, 3}{100, 101, 102}{0.0001, 2e-05, 1000}}, "Failed to sort a matrix on the last column");
  sortedMatrix_secondColumn = sortedMatrix_lastColumn%1;
  assert(sortedMatrix_secondColumn == {{0.0001, 2e-05, 1000}{0.1, 0.8, 0.2}{1, 2, 3}{100, 101, 102}}, "Failed to sort a matrix on the second column");
  // Compare strings in a case-insensitive way.
  assert("HyPhy"%"HYPHY" == 1, "Failed to correctly compare two strings that are spelled the same but have differnet case");
  assert("HyPhy"%"HyPhy" == 1, "Failed to correctly compare two strings that are spelled the same and are the same case");
  assert("HyPhy"%"AAAAAA" == 0, "Failed to correctly compare two strings that are completly different");
  assert("HyPhy"%"HyPhy2" == 0, "Failed to correctly compare two strings that are only slightly different");
  

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);

  assert (runCommandWithSoftErrors ('T%T', "is not implemented/defined for a Topology"), "Failed error checking for trying modulo/compare topologies (%)");
  assert (runCommandWithSoftErrors ('TT%TT', "is not implemented/defined for a Tree"), "Failed error checking for trying modulo/compare trees (%)");
  assert (runCommandWithSoftErrors ('4$"String"', "where 'X' is not a number"), "Failed error checking for trying module/compare number%string");
 

  testResult = 1;

  return testResult;
}
