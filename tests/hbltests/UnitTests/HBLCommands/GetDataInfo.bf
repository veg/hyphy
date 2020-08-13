function getTestName () {
	return "GetDataInfo";
}

function getTestedFunctions () {
	return {{"_ExecutionList::ExecuteCase46","_DataSetFilter::FindUniqueSequences"}};
}

function runTest () {
	/* define some auxiliary objects here */

	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult  	   = 0;

	DataSet 			fiveSeqs = ReadDataFile (PATH_TO_CURRENT_BF + "../../data/5.fas");	
	DataSetFilter		nucF	 = CreateFilter (fiveSeqs,1);

	GetDataInfo 		(seqInfo, nucF, -1);

	assert (seqInfo["UNIQUE_SEQUENCES"] == 5, "Expected 5 unique sequences with strict filtering");
	GetDataInfo 		(seqInfo, nucF, -2);
	assert (seqInfo["UNIQUE_SEQUENCES"] == 4, "Expected 4 unique sequences with strict+gap filtering");
	GetDataInfo 		(seqInfo, nucF, -3);
	assert (seqInfo["UNIQUE_SEQUENCES"] == 3, "Expected 3 unique sequences with superset filtering");
	GetDataInfo 		(seqInfo, nucF, -4);
	assert (seqInfo["UNIQUE_SEQUENCES"] == 2, "Expected 2 unique sequences with partial match filtering");

	DataSetFilter		dinucF	 = CreateFilter (fiveSeqs,2);

	GetDataInfo 		(seqInfo, dinucF, -2);
	assert (seqInfo["UNIQUE_SEQUENCES"] == 5, "Expected 5 unique sequences with strict+gap filtering (dinuc)");

	testResult = 1;
	return testResult;
}

/* execution stub */

fprintf    (stdout, "[Running COVERAGE TEST '", getTestName(), "']\n");
result  =  runTest();

if (result) {
	fprintf (stdout, "[TEST PASSED]\n");
} else {
	fprintf (stdout, "[TEST FAILED]\n");
}
