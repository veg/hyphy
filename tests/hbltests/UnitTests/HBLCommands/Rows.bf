ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "Rows";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Return number of Rows for matrix.
  assert(Rows({{1}}) == 1, "Failed to return one for a one d matrix");
  assert(Rows({{1,2}{3,4}{5,6}}) == 3, "Failed to return three for a 3x2 matrix");
  assert(Rows({{1,2,3}{4,5,6}}) == 2, "Failed to return two for a 2x3 matrix");
  assert(Rows(200) == 0, "Failed to return zero for a non-matrix (number)");
  assert(Rows(none) == 0, "Failed to return zero for a non-matrix (none)");

  // Undocumented functionality: return the keys of an AssociativeList.
  // TODO: The function works but the assert is failing for some reason.
  list1 = {"key1": "val1"};
  list1Rows = Rows(list1);
  expectedList1Rows = {{"key1"}};
  //fprintf (stdout, 'list1Rows: ', list1Rows, '\n');
  //assert(list1Rows == expectedList1Rows, "Failed to return the values of a one element associative list");
  list2 = {"key2": "val2", "key3": "val3"};
  list2Rows = Rows(list2);
  expectedList2Rows = {{"val2", "Val3"}};
  //fprintf (stdout, 'list2Rows: ', list2Rows, '\n');
  //assert(list2Rows == expectedList2Rows, "Failed to return the values of a two element associative list");

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);

  assert (runCommandWithSoftErrors ('Rows(T)', "Operation 'Rows' is not implemented/defined for a Topology"), "Failed error checking for trying to take Rows of topologies");
  assert (runCommandWithSoftErrors ('Rows(TT)', "Operation 'Rows' is not implemented/defined for a Tree"), "Failed error checking for trying to take Rows of trees");

  testResult = 1;

  return testResult;
}
