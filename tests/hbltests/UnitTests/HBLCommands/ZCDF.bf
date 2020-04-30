ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();

function getTestName ()
{
	return "Standard normal cumulative distribution function";
}

function runTest()
{
  testResult = 0;
  // Critical values
  assert(ZCDF(0) == .5, "ZCDF of zero is not .5");

  x = {{ // np.linspace(-5, 5, 25)
    -5.0, -4.583333333333333, -4.166666666666667, -3.75, -3.333333333333333,
    -2.9166666666666665, -2.5, -2.083333333333333, -1.6666666666666665, -1.25,
    -0.833333333333333, -0.4166666666666661, 0.0, 0.41666666666666696, 0.8333333333333339,
    1.25, 1.666666666666667, 2.083333333333334, 2.5, 2.916666666666667,
    3.333333333333334, 3.75, 4.166666666666668, 4.583333333333334, 5.0
  }};
  zcdf_x = {{
    2.8665157187919333e-07, 2.2881082869581217e-06, 1.5454296882295967e-05, 8.841728520080377e-05, 0.0004290603331968372,
    0.0017689682391110525, 0.006209665325776132, 0.01861042518988636, 0.0477903522728147, 0.10564977366685535,
    0.20232838096364308, 0.3384611195106898, 0.5, 0.6615388804893104, 0.7976716190363572,
    0.8943502263331446, 0.9522096477271853, 0.9813895748101137, 0.9937903346742238, 0.9982310317608889,
  0.9995709396668032, 0.9999115827147992, 0.9999845457031177, 0.999997711891713, 0.9999997133484281
  }}; 

  assert(Abs(x["ZCDF(_MATRIX_ELEMENT_VALUE_)"] - zcdf_x) < 1e-8, "Does not agree with existing numerical computing frameworks");

  //Invalid types
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);

  assert(runCommandWithSoftErrors('ZCDF("")', "Operation 'ZCDF' is not implemented/defined for a String"), "ZCDF String argument is invalid");
  assert(runCommandWithSoftErrors('ZCDF({})', "Operation 'ZCDF' is not implemented/defined for a AssociativeList"), "ZCDF Associative list argument is invalid");
  assert(runCommandWithSoftErrors('ZCDF({{"a","b"}})', "Operation 'ZCDF' is not implemented/defined for a Matrix"), "ZCDF Matrix argument is invalid");
  assert(runCommandWithSoftErrors('ZCDF(T)', "Operation 'ZCDF' is not implemented/defined for a Topology"), "ZCDF Topology argument is invalid");
  assert(runCommandWithSoftErrors('ZCDF(TT)', "Operation 'ZCDF' is not implemented/defined for a Tree"), "ZCDF Tree argument is invalid");

	testResult = 1;
	return testResult;
}
