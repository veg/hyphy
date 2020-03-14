ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "Random";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // TODO: These test do not attempt to test the "randomness" of the pseudo random number generator, only its functionality.
  exampleMatrix = {{1,2,3}{4,5,6}{7,8,9}};

  // Generate random number between two values. Random(Number min,Number max). Returns a random float out to 4 decimal places between the min and the max.
  assert(Random(1,2) < 2, "Random number generated greater than max");
  assert(Random(1,2) > 1, "Random number generated less than min");
  assert(Random(-2,-1) < -1, "Negative random number generated greater than max");
  assert(Random(0,100) != Random(0,100), "Two randomly generated numbers were equal (not impossible but extremly rare)");

  // Suffle a matrix. Random(Matrix m, Number n). If n>1, then allow duplicates(permute) when reshuffling.
  sumOfFirstElementsTenTimesAfterShuffling = 0;
  sumOfMaxInMatrix = 0;
  for (x=0; x<100; x+=1) {
    shuffledExampleMatrix = Random(exampleMatrix, 0);
    sumOfFirstElementsTenTimesAfterShuffling += shuffledExampleMatrix[0];
    shuffledWithReplacementExampleMatrix = Random(exampleMatrix, 2);
    sumOfMaxInMatrix += Max(shuffledWithReplacementExampleMatrix, 2);
  }
  assert(sumOfFirstElementsTenTimesAfterShuffling != 100, "Failed to shuffle a matrix without resampling (shuffled 100 times and tested if the first element was the same each time)");
  assert(sumOfMaxInMatrix != 1000, "Failed to shuffle a matrix with resampling (shuffled 100 times and tested that the max changed at least once)");

  // Latin hypercube sampling if str is "LHS". Random(Matrix m, String str).
  // TODO: confirm that the observed behaviour is desired.
  // By observation, the resulting column index stays the same but the row index is randomly suffled.
  sumOfFirstColumn = 0;
  sumOfFirstRow = 0;
  for (j=0; j<100; j+=1) {
    LHSedExampleMatrix = Random(exampleMatrix, "LHS");
    sumOfFirstColumn += (LHSedExampleMatrix[0][0] + LHSedExampleMatrix[0][1] + LHSedExampleMatrix[0][2]);
    sumOfFirstRow += (LHSedExampleMatrix[0][0] + LHSedExampleMatrix[1][0] + LHSedExampleMatrix[2][0]);
  }
  assert(sumOfFirstColumn != 600, "Failed to run Random(matrix, LHS) as expected; the set of values in the first row stayed the same over 100 iterations.");
  assert(sumOfFirstRow == 1200, "Failed to run Random(matrix, LHS) as expected; the set of values in the first column changed");

  // A matrix m and Associative List of the form {"PDF":PDF_SORT,"ARG0":ARG0,"ARG1":ARG1...,"ARGN":ARGN}
  // PDF Type	Explanation: Dirichlet:	Dirichlet distribution; Gaussian:	multivariate Gaussian distribution; Wishart:	Wishart distribution; InverseWishart:	Inverse-Wishart distribution; Multinomial	Multinomial distribution

  // TODO: testing the multivariate gaussian or dirichlet sampling returns the following error:
  // HYPHYMP(53058,0x7fffa3bcb380) malloc: *** error for object 0x7ffeed534218: pointer being freed was not allocated
  // *** set a breakpoint in malloc_error_break to debug

  // The code below returns the errors
  mean1 = {{1,1}};
  cov1 = {{1,0}{0,1}};
  a1 = {"PDF":"Gaussian","ARG0":cov1};
  z = Random(mean1,a1);

  mean2 = {{1,1,1}};
  cov2 = {{1,0}{0,1}};
  a2 = {"PDF":"Dirichlet","ARG0":cov2};
  b = Random(mean2,a2);
  assert(+b == 1, "Dirichet deviate did not sum up to 1");
   

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------

  exampleString = "abcd";
  Tree exampleTree = ((1,2),(3,4),5);
  Topology exampleTopology = ((1,2),(3,4),5);
  exampleList = {"key1": "val1", "key2": "val2"};

  assert (runCommandWithSoftErrors ('Random(exampleList, exampleList)', "Unsupported argument type"), "Failed error checking for trying to run Random(AssociativeList, AssociativeList)");
  assert (runCommandWithSoftErrors ('Random(exampleTree, exampleTree)', "Unsupported agrument type 'Tree'"), "Failed error checking for trying to run Random(Tree, Tree)");
  assert (runCommandWithSoftErrors ('Random(exampleTopology, exampleTopology)', "Unsupported agrument type 'Topology'"), "Failed error checking for trying to run Random(Topology, Topology)");
  assert (runCommandWithSoftErrors ('Random(exampleString, exampleString)', "Operation 'Random' is not implemented/defined for a String"), "Failed error checking for trying to run Random(string, string)");

  testResult = 1;

  return testResult;
}

