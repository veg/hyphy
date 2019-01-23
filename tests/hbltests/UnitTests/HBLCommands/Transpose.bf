ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();

function getTestName ()
{
	return "Transpose";
}

function runTest()
{
  testResult = 0;
  // Vector
  x = {{1, 2, 3}};
  xT = {{1}{2}{3}};
  assert(Transpose(x)==xT, "Does not agree with a vector that was manually transposed");
  assert(x==Transpose(xT), "Does not agree with a vector that was manually transposed");

  // Rectangular matrix
  X = {{1, 2, 3}{4, 5, 6}};
  XT = {{1, 4}{2, 5}{3, 6}};
  assert(Transpose(X)==XT, "Does not agree with a rectangular matrix that was manually transposed");
  assert(X==Transpose(XT), "Does not agree with a rectangular matrix that was manually transposed");

  // Square matrix
  Y = {{10, 20}{30, 40}};
  YT = {{10, 30}{20, 40}};
  assert(Transpose(Y)==YT, "Does not agree with a square matrix that was manually transposed");
  assert(Y==Transpose(YT), "Does not agree with a square matrix that was manually transposed");

  // Invalid types
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);

  assert(runCommandWithSoftErrors('Transpose("")', "Operation 'Transpose' is not implemented/defined for a String"), "Transpose String argument is invalid");
  assert(runCommandWithSoftErrors('Transpose({})', "Operation 'Transpose' is not implemented/defined for a AssociativeList"), "Transpose Associative list argument is invalid");
  assert(runCommandWithSoftErrors('Transpose(T)', "Operation 'Transpose' is not implemented/defined for a Topology"), "Transpose Topology argument is invalid");
  assert(runCommandWithSoftErrors('Transpose(TT)', "Operation 'Transpose' is not implemented/defined for a Tree"), "Transpose Tree argument is invalid");

  testResult = 1;
	return testResult;
}
