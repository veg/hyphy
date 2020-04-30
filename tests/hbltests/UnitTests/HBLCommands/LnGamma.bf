ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();

function getTestName ()
{
	return "Logarithm of Gamma function";
}

function runTest()
{
  testResult = 0;

  // x = np.linspace(1, 10, 10)
  x = {{ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 }};
  // gammaln(x)
  lngamma_x = {{ 0.0, 0.0, 0.693147, 1.791759, 3.178054, 4.787492, 6.579251, 8.525161, 10.604603, 12.801827 }};
  assert(Abs(x["LnGamma(_MATRIX_ELEMENT_VALUE_)"] - lngamma_x) < 1e-6, "Does not agree with existing numerical computing frameworks");

  //Invalid types
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);

  assert(runCommandWithSoftErrors('LnGamma("")', "Operation 'LnGamma' is not implemented/defined for a String"), "LnGamma String argument is invalid");
  assert(runCommandWithSoftErrors('LnGamma({})', "Operation 'LnGamma' is not implemented/defined for a AssociativeList"), "LnGamma Associative list argument is invalid");
  assert(runCommandWithSoftErrors('LnGamma({{"a","b"}})', "Operation 'LnGamma' is not implemented/defined for a Matrix"), "LnGamma Matrix argument is invalid");
  assert(runCommandWithSoftErrors('LnGamma(T)', "Operation 'LnGamma' is not implemented/defined for a Topology"), "LnGamma Topology argument is invalid");
  assert(runCommandWithSoftErrors('LnGamma(TT)', "Operation 'LnGamma' is not implemented/defined for a Tree"), "LnGamma Tree argument is invalid");

  testResult = 1;
	return testResult;
}
