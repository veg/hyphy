ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "fscanf_fprintf";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = TRUE;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------

  // Generate a unique file Path.
  tempFilePathString = './../../data/tempFileTesting-fprintf_fscanf' + Random(0,1);
  tempFilePathMatrix = './../../data/tempFileTesting-fprintf_fscanf_matrix' + Random(0,1);
  // Write to the file.
  fprintf(tempFilePathString, "this string is written to and read from a file");
  fprintf(tempFilePathMatrix, matrix1);
  // Read from the file.
  fscanf(tempFilePathString, String, testString);
  fscanf(tempFilePathString, Matrix, testMatrix);
  // Make sure it worked as expected.
  expectedTestString = 'this string is written to and read from a file';
  assert(testString == expectedTestString, "fscanf and fprintf String from and to files did not work as expected");
  assert(testMatrix == matrix1, "fscanf and fprintf Matrix from and to files did not work as expected");

  // Same for Matrix.
  matrix1 = {{0.3,0.4}};
  tempFilePathMatrix = './../../data/tempFileTesting-fprintf_fscanf_matrix' + Random(0,1);
  fprintf(tempFilePathMatrix, matrix1);
  fscanf(tempFilePathMatrix, Matrix, testMatrix);
  assert(testMatrix == matrix1, "fscanf and fprintf Matrix from and to files did not work as expected");

  // Same for Tree.
  /* TODO: results in seg fault.
  Tree TT1 = ((1,2),(3,4),5);
  tempFilePathTree = './../../data/tempFileTesting-fprintf_fscanf_tree' + Random(0,1);
  fprintf(tempFilePathTree, TT1);
  fscanf(tempFilePathTree, Tree, testTree);
  assert(testTree == TT1, "fscanf and fprintf Tree from and to files did not work as expected");
  */

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  assert (runCommandWithSoftErrors ("fscanf('./thisFileDoesNotExist', String, x);", "could not be opened for reading by fscanf"), "Failed error checking for trying to read from a file that does not exist");
  

  
  testResult = 1;

  return testResult;
}
