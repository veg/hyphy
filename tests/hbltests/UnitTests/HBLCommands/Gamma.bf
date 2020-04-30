ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();

function getTestName ()
{
	return "Gamma function";
}

function runTest()
{
  testResult = 0;

  x = {{ //np.concatenate([np.linspace(-10.5, -.5, 10), np.linspace(1, 10, 15)])
    -10.5, -9.38888888888889, -8.277777777777779, -7.166666666666666, -6.055555555555555,
    -4.944444444444445, -3.833333333333333, -2.7222222222222214, -1.6111111111111107, -0.5,
    1.0, 1.6428571428571428, 2.2857142857142856, 2.928571428571429, 3.5714285714285716,
    4.214285714285714, 4.857142857142858, 5.5, 6.142857142857143, 6.7857142857142865,
    7.428571428571429, 8.071428571428573, 8.714285714285715, 9.357142857142858, 10.0
  }};
  gamma_x = {{
    -2.640121820547716e-07, 3.807937180089331e-06, -5.587096541725678e-05, 0.0008893140184485769, -0.022639090719436133,
    -0.16570696708613508, 0.33545433798219826, -0.9580642143945058, 2.319024629802521, -3.5449077018110318,
    1.0, 0.8990555517429472, 1.156817798360793, 1.874329490574376, 3.598822676556153,
    7.9035848077182616, 19.39812956539834, 52.34277778455352, 153.40129668244913, 483.71930429998673,
    1628.9243800135312, 5822.437791883003, 21978.98383450019, 87247.96540737787, 362880.0
  }}; 
  abs_gamma_x = gamma_x["Abs(_MATRIX_ELEMENT_VALUE_)"];
  assert(Abs((x["Gamma(_MATRIX_ELEMENT_VALUE_)"] - gamma_x)/abs_gamma_x) < 1e-8, "Does not agree with existing numerical computing frameworks");

  //Invalid types
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);

  assert(runCommandWithSoftErrors('Gamma("")', "Operation 'Gamma' is not implemented/defined for a String"), "Gamma String argument is invalid");
  assert(runCommandWithSoftErrors('Gamma({})', "Operation 'Gamma' is not implemented/defined for a AssociativeList"), "Gamma Associative list argument is invalid");
  assert(runCommandWithSoftErrors('Gamma({{"a","b"}})', "Operation 'Gamma' is not implemented/defined for a Matrix"), "Gamma Matrix argument is invalid");
  assert(runCommandWithSoftErrors('Gamma(T)', "Operation 'Gamma' is not implemented/defined for a Topology"), "Gamma Topology argument is invalid");
  assert(runCommandWithSoftErrors('Gamma(TT)', "Operation 'Gamma' is not implemented/defined for a Tree"), "Gamma Tree argument is invalid");

	testResult = 1;
	return testResult;
}
