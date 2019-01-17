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
  assert(Transpose(x)==xT);
  assert(x==Transpose(xT));

  // Rectangular matrix
  X = {{1, 2, 3}{4, 5, 6}};
  XT = {{1, 4}{2, 5}{3, 6}};
  assert(Transpose(X)==XT);
  assert(X==Transpose(XT));

  // Square matrix
  Y = {{10, 20}{30, 40}};
  YT = {{10, 30}{20, 40}};
  assert(Transpose(Y)==YT);
  assert(Y==Transpose(YT));

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
