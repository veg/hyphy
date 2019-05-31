ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();

function getTestName ()
{
	return "Inverse tangent";
}

function runTest()
{
  testResult = 0;
  // Critical values
  assert(4*Arctan(1) == 3.141592653589793, "Arctan of pi over 4 is not 1");
  assert(Arctan(0) == 0, "Arctan of 0 is not 0");

  // Comparison with NumPy
  x = {{ // np.linspace(-5, 5, 25)
    -5.0, -4.583333333333333, -4.166666666666667, -3.75, -3.333333333333333,
    -2.9166666666666665, -2.5, -2.083333333333333, -1.6666666666666665, -1.25,
    -0.833333333333333, -0.4166666666666661, 0.0, 0.41666666666666696, 0.8333333333333339,
    1.25, 1.666666666666667, 2.083333333333334, 2.5, 2.916666666666667,
    3.333333333333334, 3.75, 4.166666666666668, 4.583333333333334, 5.0
  }};
  arctan_x = {{
    -1.373400766945016, -1.3559809263932379, -1.3352513460740334, -1.3101939350475558, -1.2793395323170293,
    -1.240498971965643, -1.1902899496825317, -1.1232763516377267, -1.0303768265243125, -0.8960553845713439,
    -0.694738276196703, -0.394791119699761, 0.0, 0.39479111969976177, 0.6947382761967036,
    0.8960553845713439, 1.0303768265243125, 1.1232763516377269, 1.1902899496825317, 1.240498971965643,
    1.2793395323170296, 1.3101939350475558, 1.3352513460740334, 1.355980926393238, 1.373400766945016
  }};
    
  assert(Abs(x["Arctan(_MATRIX_ELEMENT_VALUE_)"] - arctan_x) < 1e-8, "Does not agree with existing numerical computing frameworks");
  
  //  Invalid types
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);

  assert(runCommandWithSoftErrors('Arctan("")', "Operation 'Arctan' is not implemented/defined for a String"), "Arctan String argument is invalid");
  assert(runCommandWithSoftErrors('Arctan({})', "Operation 'Arctan' is not implemented/defined for a AssociativeList"), "Arctan Associative list argument is invalid");
  assert(runCommandWithSoftErrors('Arctan({{"a","b"}})', "Operation 'Arctan' is not implemented/defined for a Matrix"), "Arctan Matrix argument is invalid");
  assert(runCommandWithSoftErrors('Arctan(T)', "Operation 'Arctan' is not implemented/defined for a Topology"), "Arctan Topology argument is invalid");
  assert(runCommandWithSoftErrors('Arctan(TT)', "Operation 'Arctan' is not implemented/defined for a Tree"), "Arctan Tree argument is invalid");
	testResult = 1;
	return testResult;
}
