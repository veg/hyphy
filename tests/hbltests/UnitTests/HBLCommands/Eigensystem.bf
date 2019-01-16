ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "Eigensystem";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Compute the eigenvectors and eigenvalues of a matrix.
  matrix1 = {{19,3}{-2,26}};
  eigenSystem1 = Eigensystem(matrix1);
  eigenVector1 = eigenSystem["0"];
  eigenValues1 = eigenSystem["1"];
  expectedEigenVector1 = {{20}{25}};
  expectedEigenValues1 = {{0}{0}};

  // TODO: see below.
  fprintf(stdout, 'The eigen values and eigen vectors are printing as <HyPhy Base Object> which I dont know how to assert against\n');
  fprintf (stdout, 'eigenSystem1: ', eigenSystem1, '\n');
  fprintf (stdout, 'eigenValues1: ', eigenValues1, '\n');
  fprintf (stdout, 'eigenValues1: ', eigenValues1, '\n');

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  Topology T = ((1,2),(3,4),5);
  Tree TT = ((1,2),(3,4),5);
  list1 = {"key1": "value1", "key2": "value1"};

  //Eigensystem(none);
  assert (runCommandWithSoftErrors ('Eigensystem(none)', "Operation 'Eigensystem' is not implemented/defined for a Number"), "Failed error checking for trying to take a Eigensystem of None");
  assert (runCommandWithSoftErrors ('Eigensystem(3)', "Operation 'Eigensystem' is not implemented/defined for a Number"), "Failed error checking for trying to take a Eigensystem of a number");
  assert (runCommandWithSoftErrors ('Eigensystem(T)', "not implemented/defined for a Topology"), "Failed error checking for trying to take a Eigensystem of a topology");
  assert (runCommandWithSoftErrors ('Eigensystem(TT)', "not implemented/defined for a Tree"), "Failed error checking for trying to take a Eigensystem of a tree");
  assert (runCommandWithSoftErrors ('Eigensystem(list1)', "not implemented/defined for a AssociativeList"), "Failed error checking for trying to take a Eigensystem of an associative list");

  testResult = 1;

  return testResult;
}
