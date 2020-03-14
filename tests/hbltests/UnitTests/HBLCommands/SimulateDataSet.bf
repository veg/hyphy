ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "SiulateDataSet";
}	


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;

  // NOTE: This just tests the functionality of the SimulateDataSet function, not the underlying logic or the quality of the generated datasets.

  // Simple example likelihood function for use later.
  DataSet         nucleotideSequences = ReadDataFile (PATH_TO_CURRENT_BF + "../../data/CD2_reduced.fna");
  DataSetFilter   filteredData = CreateFilter (nucleotideSequences,1);
  HarvestFrequencies (observedFreqs, filteredData, 1, 1, 1);
  F81RateMatrix = 
        {{*,mu,mu,mu}
         {mu,*,mu,mu}
         {mu,mu,*,mu}
         {mu,mu,mu,*}};
  Model   F81 = (F81RateMatrix, observedFreqs);
  Tree    givenTree = DATAFILE_TREE;
  LikelihoodFunction  LF = (filteredData, givenTree);
  Optimize (paramValues, LF);
  
  // Simulate the dataset twice
  DataSet simulatedDataSet1 = SimulateDataSet(LF);
  DataSet simulatedDataSet2 = SimulateDataSet(LF);

  // Create dataset filters for each and get the sequences 
  DataSetFilter filteredDataSet1 = CreateFilter (simulatedDataSet1,1);
  DataSetFilter filteredDataSet2 = CreateFilter (simulatedDataSet2,1);
  GetInformation(simulatedSequences1, filteredDataSet1);
  GetInformation(simulatedSequences2, filteredDataSet2);

  // Test to confirm that the simulated data was the same dimmensions (sites and sequences) as the data used to creaet the likelihood function
  expectedNumberOfSequences = 10;
  assert(Columns(simulatedSequences1) == expectedNumberOfSequences, "Failed to simulate a dataset that had the same number of sequences (tips) as the original dataset on which the likelihood function was fit");
  expectedNumberOfSites = 18;
  numberOfSites = 0;
  firstSequence = simulatedSequences1[0];
  for (i = 0; i < expectedNumberOfSites+1; i+=1) {
    if (firstSequence[i]) {
      numberOfSites+=1;
    }
  }
  assert(numberOfSites == expectedNumberOfSites, "Failed to simulate a dataset that had the same number of sites as the orignal dataset on which the likelihood function was fit");

  // Test to confirm that the two simulated datasets are different
  numberOfMatchingSequences = 0;
  for (i = 0; i < expectedNumberOfSequences; i+=1) {
    if (simulatedSequences1 == simulatedDataSet2) {
      numberOfMatchingSequences+=1;
    }
  }
  assert(numberOfMatchingSequences < 8, "Two simulated datasets, based on the same likelihood function had 80% of the timp sequences matching; this shouldn't happen");

  //DataSet notFromLF = SimulateDataSet(1);
  assert (runCommandWithSoftErrors ('DataSet notFromLF = SimulateDataSet(1);', 'has not been initialized'), "Failed error checking for trying to simulate a dataset with an object other than an initialized likelihood function");
  assert (runCommandWithSoftErrors ('DataSet notFromLF = SimulateDataSet();', 'has not been initialized'), "Failed error checking for trying to simulate a dataset with no argument");

  testResult = 1;

  return testResult;
}
