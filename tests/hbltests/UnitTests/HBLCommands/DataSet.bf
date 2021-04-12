ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "DataSet";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = TRUE;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Test that DataSet can load the supported filetypes without errors (correctness of the loaded files is tested in the DataSetFilter.bf test)  
  DataSet cd2nex = ReadDataFile (PATH_TO_CURRENT_BF + '/../../data/CD2.nex');
  DataSet 2fas = ReadDataFile (PATH_TO_CURRENT_BF  + '/../../data/2.fas');
  DataSet cd2Phylip = ReadDataFile(PATH_TO_CURRENT_BF  + '/../../data/CD2.phylip');
  
  

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  list1 = {'key1':'val1', 'key2':'val2'};
  matrix1 = {{0.3,0.4}};
  Topology T1 = ((1,4),(3,4),5);
  Tree TT1 = ((1,2),(3,4),5);

  assert (runCommandWithSoftErrors ('DataSet list_ds = ReadFromString(list1);', "The format of the sequence file has not been recognized and may be invalid"), "Failed error checking for trying to create a data set with ReadFromString(list)");
  assert (runCommandWithSoftErrors ('DataSet list_ds = ReadFromString(matrix1);', "The format of the sequence file has not been recognized and may be invalid"), "Failed error checking for trying to create a data set with ReadFromString(matrix)");
  assert (runCommandWithSoftErrors ('DataSet list_ds = ReadFromString(T1);', "The format of the sequence file has not been recognized and may be invalid"), "Failed error checking for trying to create a data set with ReadFromString(topology)");
  assert (runCommandWithSoftErrors ('DataSet list_ds = ReadFromString(TT1);', "The format of the sequence file has not been recognized and may be invalid"), "Failed error checking for trying to create a data set with ReadFromString(tree)");
  

  assert (runCommandWithSoftErrors ("DataSet thisIsntAValidFilePath = ReadDataFile('./ThisFileDoesNotExist.txt');", "Could not find source dataset file"), "Failed error checking for trying to create a data set with ReadDataFile with an invalid path");

  assert (runCommandWithSoftErrors ("DataSet thisFileIsntInAValidFormat = ReadDataFile(PATH_TO_CURRENT_BF + '/assert.bf');", "The format of the sequence file has not been recognized and may be invalid"), "Failed error checking for trying to create a data set with ReadDataFile with a file in an invalid format");
  assert (runCommandWithSoftErrors ("DataSet newickFile = ReadDataFile (PATH_TO_CURRENT_BF + '/../../data/CD2.newick');", "The format of the sequence file has not been recognized and may be invalid"), "Failed error checking for trying to create a data set with ReadDataFile with a file file containing a newick string but no sequences");
  
  testResult = 1;

  return testResult;
}
