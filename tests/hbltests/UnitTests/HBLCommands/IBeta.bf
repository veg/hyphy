ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();

function getTestName ()
{
	return "Incomplete Beta function";
}

function runTest()
{
  testResult = 0;
  tol = 1e-6;
  // np.linspace(0, 1, 11)
  x = {{ 0. , 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1. }};

  // betainc(.5, .5, x)
  B = {{ 0.0, 0.204833, 0.295167, 0.36901, 0.435906, 0.5, 0.564094, 0.63099, 0.704833, 0.795167, 1.0 }};
  for(i=0; i<11; i=i+1) {
    assert(Abs(B[i] - IBeta(x[i], .5, .5)) < tol, "Does not agree with existing numerical computing frameworks");
  }
  // betainc(1, .5, x)
  B = {{ 0.0, 0.051317, 0.105573, 0.16334, 0.225403, 0.292893, 0.367544, 0.452277, 0.552786, 0.683772, 1.0 }};
  for(i=0; i<11; i=i+1) {
    assert(Abs(B[i] - IBeta(x[i], 1, .5)) < tol, "Does not agree with existing numerical computing frameworks");
  }
  // betainc(2, .5, x)
  B = {{ 0.0, 0.003883, 0.01613, 0.037841, 0.070484, 0.116117, 0.177808, 0.260575, 0.373901, 0.54147, 1.0 }};
  for(i=0; i<11; i=i+1) {
    assert(Abs(B[i] - IBeta(x[i], 2, .5)) < tol, "Does not agree with existing numerical computing frameworks");
  }
  // betainc(.5, 1, x)
  B = {{ 0.0, 0.316228, 0.447214, 0.547723, 0.632456, 0.707107, 0.774597, 0.83666, 0.894427, 0.948683, 1.0 }};
  for(i=0; i<11; i=i+1) {
    assert(Abs(B[i] - IBeta(x[i], .5, 1)) < tol, "Does not agree with existing numerical computing frameworks");
  }
  // betainc(1, 1, x)
  B = {{ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 }};
  for(i=0; i<11; i=i+1) {
    assert(Abs(B[i] - IBeta(x[i], 1, 1)) < tol, "Does not agree with existing numerical computing frameworks");
  }
  // betainc(2, 1, x)
  B = {{ 0.0, 0.01, 0.04, 0.09, 0.16, 0.25, 0.36, 0.49, 0.64, 0.81, 1.0 }};
  for(i=0; i<11; i=i+1) {
    assert(Abs(B[i] - IBeta(x[i], 2, 1)) < tol, "Does not agree with existing numerical computing frameworks");
  }
  // betainc(.5, 2, x)
  B = {{ 0.0, 0.45853, 0.626099, 0.739425, 0.822192, 0.883883, 0.929516, 0.962159, 0.98387, 0.996117, 1.0 }};
  for(i=0; i<11; i=i+1) {
    assert(Abs(B[i] - IBeta(x[i], .5, 2)) < tol, "Does not agree with existing numerical computing frameworks");
  }
  // betainc(1, 2, x)
  B = {{ 0.0, 0.19, 0.36, 0.51, 0.64, 0.75, 0.84, 0.91, 0.96, 0.99, 1.0 }};
  for(i=0; i<11; i=i+1) {
    assert(Abs(B[i] - IBeta(x[i], 1, 2)) < tol, "Does not agree with existing numerical computing frameworks");
  }
  // betainc(2, 2, x)
  B = {{ 0.0, 0.028, 0.104, 0.216, 0.352, 0.5, 0.648, 0.784, 0.896, 0.972, 1.0 }};
  for(i=0; i<11; i=i+1) {
    assert(Abs(B[i] - IBeta(x[i], 2, 2)) < tol, "Does not agree with existing numerical computing frameworks");
  }

  //Invalid types
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);

  assert(runCommandWithSoftErrors('IBeta("", 1, 1)', "Operation 'IBeta' is not implemented/defined for a String"), "IBeta String argument is invalid");
  assert(runCommandWithSoftErrors('IBeta({}, 1, 1)', "Operation 'IBeta' is not implemented/defined for a AssociativeList"), "IBeta Associative list argument is invalid");
  assert(runCommandWithSoftErrors('IBeta({{"a","b"}}, 1, 1)', "Operation 'IBeta' is not implemented/defined for a Matrix"), "IBeta Matrix argument is invalid");
  assert(runCommandWithSoftErrors('IBeta(T, 1, 1)', "Operation 'IBeta' is not implemented/defined for a Topology"), "IBeta Topology argument is invalid");
  assert(runCommandWithSoftErrors('IBeta(TT, 1, 1)', "Operation 'IBeta' is not implemented/defined for a Tree"), "IBeta Tree argument is invalid");

  testResult = 1;
	return testResult;
}
