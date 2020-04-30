ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();

function getTestName ()
{
	return "Incomplete Gamma function";
}

function runTest()
{
  testResult = 0;

  // np.arange(.5, 5, .5)
  s = {{ 0.5, 1. , 1.5, 2. , 2.5, 3. , 3.5, 4. , 4.5 }};
  // np.linspace(0, 10, 11)
  x = {{ 0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9., 10. }};
  // gammainc(s[i], x[j])
  G = {
    { 0.0, 0.842701, 0.9545, 0.985694, 0.995322, 0.998435, 0.999468, 0.999817, 0.999937, 0.999978, 0.999992 },
    { 0.0, 0.632121, 0.864665, 0.950213, 0.981684, 0.993262, 0.997521, 0.999088, 0.999665, 0.999877, 0.999955 },
    { 0.0, 0.427593, 0.738536, 0.88839, 0.953988, 0.981434, 0.992617, 0.997095, 0.998866, 0.99956, 0.99983 },
    { 0.0, 0.264241, 0.593994, 0.800852, 0.908422, 0.959572, 0.982649, 0.992705, 0.996981, 0.998766, 0.999501 },
    { 0.0, 0.150855, 0.450584, 0.693781, 0.843764, 0.924765, 0.965212, 0.984391, 0.993156, 0.997054, 0.99875 },
    { 0.0, 0.080301, 0.323324, 0.57681, 0.761897, 0.875348, 0.938031, 0.970364, 0.986246, 0.993768, 0.997231 },
    { 0.0, 0.04016, 0.220223, 0.460251, 0.667406, 0.811427, 0.899441, 0.948819, 0.974884, 0.98803, 0.99443 },
    { 0.0, 0.018988, 0.142877, 0.352768, 0.56653, 0.734974, 0.848796, 0.918235, 0.95762, 0.978774, 0.989664 },
    { 0.0, 0.008532, 0.088587, 0.260082, 0.465854, 0.649515, 0.786691, 0.877675, 0.933118, 0.964826, 0.982088 }
  };

  for(i=0; i<9; i=i+1) {
    for(j=0; j<11; j=j+1) {
      assert(Abs(G[i][j] - IGamma(s[i], x[j])) < 1e-6, "Does not agree with existing numerical computing frameworks");
    }
  }

  //Invalid types
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);

  assert(runCommandWithSoftErrors('IGamma("", 1)', "Operation 'IGamma' is not implemented/defined for a String"), "IGamma String argument is invalid");
  assert(runCommandWithSoftErrors('IGamma({}, 1)', "Operation 'IGamma' is not implemented/defined for a AssociativeList"), "IGamma Associative list argument is invalid");
  assert(runCommandWithSoftErrors('IGamma({{"a","b"}}, 1)', "Operation 'IGamma' is not implemented/defined for a Matrix"), "IGamma Matrix argument is invalid");
  assert(runCommandWithSoftErrors('IGamma(T, 1)', "Operation 'IGamma' is not implemented/defined for a Topology"), "IGamma Topology argument is invalid");
  assert(runCommandWithSoftErrors('IGamma(TT, 1)', "Operation 'IGamma' is not implemented/defined for a Tree"), "IGamma Tree argument is invalid");

  testResult = 1;
	return testResult;
}
