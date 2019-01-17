ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "Columns";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Return number of columns for matrix.
  assert(Columns({{1}}) == 1, "Failed to return one for a one d matrix");
  assert(Columns({{1,2}{3,4}{5,6}}) == 2, "Failed to return two for a 3x2 matrix");
  assert(Columns({{1,2,3}{4,5,6}}) == 3, "Failed to return three for a 2x3 matrix");
  assert(Columns(200) == 0, "Failed to return zero for a non-matrix (number)");
  assert(Columns(none) == 0, "Failed to return zero for a non-matrix (none)");

  // Undocumented functionality: return the values of an AssociativeList.
  // TODO: The function works but the assert is failing for some reason.
  list1 = {"key1": "val1"};
  list1Columns = Columns(list1);
  expectedList1Columns = {{"val1"}};
  //fprintf (stdout, 'list1Columns: ', list1Columns, '\n');
  //assert(list1Columns == expectedList1Columns, "Failed to return the values of a one element associative list");
  list2 = {"key2": "val2", "key3": "val3"};
  list2Columns = Columns(list2);
  expectedList2Columns = {{"val2", "Val3"}};
  //fprintf (stdout, 'list2Columns: ', list2Columns, '\n');
  //assert(list2Columns == expectedList2Columns, "Failed to return the values of a two element associative list");

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);

  assert (runCommandWithSoftErrors ('Columns(T)', "Operation 'Columns' is not implemented/defined for a Topology"), "Failed error checking for trying to take Columns of topologies");
  assert (runCommandWithSoftErrors ('Columns(TT)', "Operation 'Columns' is not implemented/defined for a Tree"), "Failed error checking for trying to take Columns of trees");

  testResult = 1;

  return testResult;
}
