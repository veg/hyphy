ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();

function getTestName ()
{
	return "Absolute value";
}

function runTest ()
{
  testResult = 0;
  assert(Abs({"key":"value", "key2":"value2"}) == 2, "Return number of items in an associative array");
  assert(Abs({{1,2,3}}) == 3.741657386773941, "Return L2 norm of a vector");
  assert(Abs({{1,2}{3,-4}}) == 4, "Return largest absolute value among all elements of a matrix");
  assert(Abs(-7.5) == 7.5, "Return absolute value of a number");
  assert(Abs("HyPhy") == 5, "Return length of a string");
  Topology T = ((a,b)N1,c,d);
  parentPostOrderTraversalIndices = Abs(T);
  assert(parentPostOrderTraversalIndices[0] == 2);
  assert(parentPostOrderTraversalIndices[1] == 2);
  assert(parentPostOrderTraversalIndices[2] == 5);
  assert(parentPostOrderTraversalIndices[3] == 5);
  assert(parentPostOrderTraversalIndices[4] == 5);
  assert(parentPostOrderTraversalIndices[5] == -1);
	testResult = 1;
	return testResult;
}
