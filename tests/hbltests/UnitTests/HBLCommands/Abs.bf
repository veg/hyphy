ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();

function getTestName () {
	return "Absolute value";
}

function runTest () {
  testResult = FALSE;
  assert(Abs({"key":"value", "key2":"value2"}) == 2, "Does not return number of items in an associative array");
  assert(Abs({{1,2,3}}) == 3.741657386773941, "Does not return L2 norm of a vector");
  assert(Abs({{1,2}{3,-4}}) == 4, "Does not return largest absolute value among all elements of a matrix");
  assert(Abs(-7.5) == 7.5, "Does not return absolute value of a number");
  assert(Abs("HyPhy") == 5, "Does not return length of a string");
  Topology T = ((a,b)N1,c,d);
  parentPostOrderTraversalIndices = Abs(T);
  assert(parentPostOrderTraversalIndices == {{ 2, 2, 5, 5, 5, -1 }}, "Does not return tree parent indices for a test case");

  assert(runCommandWithSoftErrors ("Abs(None)", "Attempting to operate on an undefined value"), "Failed error checking Abs (None)");
  assert(runCommandWithSoftErrors ("Abs(1,2)", "Unconsumed values on the stack"), "Failed error checking Abs (1,2)");


  testResult = TRUE;
  return testResult;
}
