ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();

function getTestName ()
{
	return "Square root";
}

function runTest()
{
  testResult = 0;
  // Critical values
  assert(Sqrt(0) == 0, "Square root of zero is not 0");
  assert(Sqrt(1) == 1, "Square root of 1 is not 1");

  x = {{ // np.linspace(0, 10, 25)
    0.0, 0.4166666666666667, 0.8333333333333334, 1.25, 1.6666666666666667,
    2.0833333333333335, 2.5, 2.916666666666667, 3.3333333333333335, 3.75,
    4.166666666666667, 4.583333333333334, 5.0, 5.416666666666667, 5.833333333333334,
    6.25, 6.666666666666667, 7.083333333333334, 7.5, 7.916666666666667,
    8.333333333333334, 8.75, 9.166666666666668, 9.583333333333334, 10.0
  }};
  sqrt_x = {{
    0.0, 0.6454972243679028, 0.9128709291752769, 1.118033988749895, 1.2909944487358056,
    1.4433756729740645, 1.5811388300841898, 1.7078251276599332, 1.8257418583505538, 1.9364916731037085,
    2.041241452319315, 2.1408720964441885, 2.23606797749979, 2.327373340628157, 2.41522945769824,
    2.5, 2.581988897471611, 2.6614532371118855, 2.7386127875258306, 2.8136571693556887,
    2.886751345948129, 2.958039891549808, 3.0276503540974917, 3.095695936834452, 3.1622776601683795
  }}; 

  assert(Abs(x["Sqrt(_MATRIX_ELEMENT_VALUE_)"] - sqrt_x) < 1e-8, "Does not agree with existing numerical computing frameworks");

  //Invalid types
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);

  assert(runCommandWithSoftErrors('Sqrt("")', "Operation 'Sqrt' is not implemented/defined for a String"), "Sqrt String argument is invalid");
  assert(runCommandWithSoftErrors('Sqrt({})', "Operation 'Sqrt' is not implemented/defined for a AssociativeList"), "Sqrt Associative list argument is invalid");
  assert(runCommandWithSoftErrors('Sqrt({{"a","b"}})', "Operation 'Sqrt' is not implemented/defined for a Matrix"), "Sqrt Matrix argument is invalid");
  assert(runCommandWithSoftErrors('Sqrt(T)', "Operation 'Sqrt' is not implemented/defined for a Topology"), "Sqrt Topology argument is invalid");
  assert(runCommandWithSoftErrors('Sqrt(TT)', "Operation 'Sqrt' is not implemented/defined for a Tree"), "Sqrt Tree argument is invalid");

	testResult = 1;
	return testResult;
}
