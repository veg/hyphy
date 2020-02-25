ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();

function getTestName ()
{
	return "Tangent";
}

function runTest()
{
  testResult = 0;
  // Critical values
  assert(Tan(0) == 0, "Tangent of zero is not 0");
  pi = 4*Arctan(1);
  assert(Abs(Tan(pi) - 0) < 1e-8, "Tangent of pi is not 0");

  x = {{ // np.linspace(-5, 5, 25)
    -5.0, -4.583333333333333, -4.166666666666667, -3.75, -3.333333333333333,
    -2.9166666666666665, -2.5, -2.083333333333333, -1.6666666666666665, -1.25,
    -0.833333333333333, -0.4166666666666661, 0.0, 0.41666666666666696, 0.8333333333333339,
    1.25, 1.666666666666667, 2.083333333333334, 2.5, 2.916666666666667,
    3.333333333333334, 3.75, 4.166666666666668, 4.583333333333334, 5.0
  }};
  tangent_x = {{
    3.380515006246585, -7.705529062699944, -1.6468091310792805, -0.6965508511114602, -0.1941255059835998,
    0.22879748135394737, 0.7470222972386602, 1.7771640735810772, 10.39877834033354, -3.0095696738628313,
    -1.100778368789801, -0.44258038400206584, 0.0, 0.44258038400206695, 1.1007783687898027,
    3.0095696738628313, -10.398778340333491, -1.7771640735810736, -0.7470222972386602, -0.2287974813539469,
    0.19412550598360073, 0.6965508511114602, 1.6468091310792838, 7.705529062699998, -3.380515006246585
  }}; 

  assert(Abs(x["Tan(_MATRIX_ELEMENT_VALUE_)"] - tangent_x) < 1e-8, "Does not agree with existing numerical computing frameworks");

  //Invalid types
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);

  assert(runCommandWithSoftErrors('Tan("")', "Operation 'Tan' is not implemented/defined for a String"), "Tan String argument is invalid");
  assert(runCommandWithSoftErrors('Tan({})', "Operation 'Tan' is not implemented/defined for a AssociativeList"), "Tan Associative list argument is invalid");
  assert(runCommandWithSoftErrors('Tan({{"a","b"}})', "Operation 'Tan' is not implemented/defined for a Matrix"), "Tan Matrix argument is invalid");
  assert(runCommandWithSoftErrors('Tan(T)', "Operation 'Tan' is not implemented/defined for a Topology"), "Tan Topology argument is invalid");
  assert(runCommandWithSoftErrors('Tan(TT)', "Operation 'Tan' is not implemented/defined for a Tree"), "Tan Tree argument is invalid");

	testResult = 1;
	return testResult;
}
