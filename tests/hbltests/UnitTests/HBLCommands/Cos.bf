ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();

function getTestName ()
{
	return "Cosine";
}

function runTest()
{
  testResult = 0;
  // Critical values
  assert(Cos(0) == 1, "Cosine of zero is not 1");
  pi = 4*Arctan(1);
  assert(Abs(Cos(pi/2)) < 1e-8, "Cosine of pi over 2 is not 0");
  assert(Abs(Cos(pi) - (-1)) < 1e-8, "Cosine of pi is not -1");
  assert(Abs(Cos(3*pi/2)) < 1e-8, "Cosine of 3*pi/2 is not 0");

  x = {{ // np.linspace(-5, 5, 25)
    -5.0, -4.583333333333333, -4.166666666666667, -3.75, -3.333333333333333,
    -2.9166666666666665, -2.5, -2.083333333333333, -1.6666666666666665, -1.25,
    -0.833333333333333, -0.4166666666666661, 0.0, 0.41666666666666696, 0.8333333333333339,
    1.25, 1.666666666666667, 2.083333333333334, 2.5, 2.916666666666667,
    3.333333333333334, 3.75, 4.166666666666668, 4.583333333333334, 5.0
  }};
  cosine_x = {{
    0.283662185463, -0.128697700557, -0.519035625385, -0.82055935734, -0.981674004711,
    -0.974810617187, -0.801143615547, -0.490389831978, -0.0957235480144, 0.315322362395,
    0.672412244083, 0.914443066594, 1.0, 0.914443066594, 0.672412244083,
    0.315322362395, -0.0957235480144, -0.490389831978, -0.801143615547, -0.974810617187,
    -0.981674004711, -0.82055935734, -0.519035625385, -0.128697700557, 0.283662185463
  }}; 

  assert(Abs(x["Cos(_MATRIX_ELEMENT_VALUE_)"] - cosine_x) < 1e-8, "Does not agree with existing numerical computing frameworks");

  //Invalid types
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);

  assert(runCommandWithSoftErrors('Cos("")', "Operation 'Cos' is not implemented/defined for a String"), "Cos String argument is invalid");
  assert(runCommandWithSoftErrors('Cos({})', "Operation 'Cos' is not implemented/defined for a AssociativeList"), "Cos Associative list argument is invalid");
  assert(runCommandWithSoftErrors('Cos({{"a","b"}})', "Operation 'Cos' is not implemented/defined for a Matrix"), "Cos Matrix argument is invalid");
  assert(runCommandWithSoftErrors('Cos(T)', "Operation 'Cos' is not implemented/defined for a Topology"), "Cos Topology argument is invalid");
  assert(runCommandWithSoftErrors('Cos(TT)', "Operation 'Cos' is not implemented/defined for a Tree"), "Cos Tree argument is invalid");

	testResult = 1;
	return testResult;
}
