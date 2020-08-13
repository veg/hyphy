ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "DataSetFilter";
}


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;


  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // A Data Filter specifies a part or parts of a data set. HyPhy provides powerful tools to select sites and sequences from a data set to analyze.
  // Syntax: DataSetFilter dataSetFilterid = CreateFilter (datasetid;unit;vertical partition; horizontal partition; alphabet exclusions);

  DataSet cd2nex = ReadDataFile (PATH_TO_CURRENT_BF + '/../../data/CD2.nex');

  DataSetFilter unChanged = CreateFilter (cd2nex,1,"","","");
  DataSetFilter removedAAA = CreateFilter (cd2nex,3,"","","AAA");
  DataSetFilter onlyFiveThroughTen = CreateFilter (cd2nex,1,"4-9");
  DataSetFilter onlyLiveStock = CreateFilter (cd2nex,1,"","4-6");

  // Test that the filters worked by comparing the frequencies of A
  HarvestFrequencies (freqsUnChanged, unChanged, 1, 1, 1);
  HarvestFrequencies (freqsRemovedAAA, removedAAA, 1, 1, 1);
  HarvestFrequencies (freqsOnlyFiveThroughTen, onlyFiveThroughTen, 1, 1, 1);
  HarvestFrequencies (freqsOnlyLiveStock, onlyLiveStock, 1, 1, 1);

  assert(Abs((freqsUnChanged[0] - freqsRemovedAAA[0]) - 0.073822) < 0.0001 , "Failed to remove AAA codons with DataSetFilter");
  assert(Abs((freqsUnChanged[0] - freqsOnlyFiveThroughTen[0]) + 0.0145053) < 0.001, "Failed to filter based on site index with DataSetFilter");
  assert(Abs((freqsUnChanged[0] - freqsOnlyLiveStock[0]) - 0.00130718) < 0.001 , "Failed to filter based on sequence index with DataSetFilter");

  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  list1 = {'key1':'val1', 'key2':'val2'};

  assert (runCommandWithSoftErrors ("DataSetFilter test1 = CreateFilter (list1,1,'','4-6');", "Expected a non-empty string matrix as the FILTER_NAMES argument in call to CreateFilter"), "Failed error checking for trying to perform a dataset filter on a non-file object");

  //TODO{low-priority}: The line below results in "Floating point exception: 8"
  //assert (runCommandWithSoftErrors ("DataSetFilter test2 = CreateFilter (cd2nex,list1,'','4-6');", "was expected to be a numerical argument"), "Failed error checking for trying to perform a dataset filter on a non-file object");

  testResult = 1;

  return testResult;
}
