ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();

function getTestName ()
{
	return "Sine";
}

function runTest()
{
  testResult = 0;
  // Critical values
  assert(Sin(0) == 0, "Sine of zero is not 0");
  pi = 4*Arctan(1);
  assert(Abs(Sin(pi/2) - 1) < 1e-8, "Sine of pi over 2 is not 1");
  assert(Abs(Sin(pi) - 0) < 1e-8, "Sine of pi is not 0");
  assert(Abs(Sin(3*pi/2) - (-1)) < 1e-8, "Sine of 3*pi/2 is not -1");

  x = {{ // np.linspace(-5, 5, 25)
    -5.0, -4.583333333333333, -4.166666666666667, -3.75, -3.333333333333333,
    -2.9166666666666665, -2.5, -2.083333333333333, -1.6666666666666665, -1.25,
    -0.833333333333333, -0.4166666666666661, 0.0, 0.41666666666666696, 0.8333333333333339,
    1.25, 1.666666666666667, 2.083333333333334, 2.5, 2.916666666666667,
    3.333333333333334, 3.75, 4.166666666666668, 4.583333333333334, 5.0
  }};
  sine_x = {{
    0.958924274663, 0.991683871943, 0.854752607239, 0.571561318742, 0.190567962875,
    -0.22303421401, -0.598472144104, -0.871503191441, -0.995407957752, -0.948984619356,
    -0.740176853196, -0.404714563561, 0.0, 0.404714563561, 0.740176853196,
    0.948984619356, 0.995407957752, 0.871503191441, 0.598472144104, 0.22303421401,
    -0.190567962875, -0.571561318742, -0.854752607239, -0.991683871943, -0.958924274663
  }}; 

  assert(Abs(x["Sin(_MATRIX_ELEMENT_VALUE_)"] - sine_x) < 1e-8, "Does not agree with existing numerical computing frameworks");

  //Invalid types
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);

  assert(runCommandWithSoftErrors('Sin("")', "Operation 'Sin' is not implemented/defined for a String"), "Sin String argument is invalid");
  assert(runCommandWithSoftErrors('Sin({})', "Operation 'Sin' is not implemented/defined for a AssociativeList"), "Sin Associative list argument is invalid");
  assert(runCommandWithSoftErrors('Sin({{"a","b"}})', "Operation 'Sin' is not implemented/defined for a Matrix"), "Sin Matrix argument is invalid");
  assert(runCommandWithSoftErrors('Sin(T)', "Operation 'Sin' is not implemented/defined for a Topology"), "Sin Topology argument is invalid");
  assert(runCommandWithSoftErrors('Sin(TT)', "Operation 'Sin' is not implemented/defined for a Tree"), "Sin Tree argument is invalid");

	testResult = 1;
	return testResult;
}
