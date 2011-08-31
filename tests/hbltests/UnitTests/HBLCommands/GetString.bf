function getTestName ()
{
	return "GetString";
}		

function getTestedFunctions ()
{
	return {{"_ExecutionList::ExecuteCase33"}};
}	

function runTest ()
{
	/* define some auxiliary objects here */

	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult  	   = 0;
	
	//-----------------------------------------------------------------------------------------------------------------
	// VERSION STRINGS
	//-----------------------------------------------------------------------------------------------------------------

	GetString (versionString, HYPHY_VERSION, 0);
	assert ((versionString$"^[0-9]+\\.[0-9a-zA-Z]+$")[0]==0, "The short version string must be of the form major.minor[beta]. Had " + versionString);
	
	GetString (versionString, HYPHY_VERSION, 1);
	assert ((versionString$"^HYPHY\\ [0-9]+\\.[0-9a-zA-Z]+.+\\ for .+$")[0]==0, "The full version string must be of the form major.minor[beta][MP] for platform description. Had " + versionString);

	GetString (versionString, HYPHY_VERSION, 2);
	assert ((versionString$"^.+\\ [0-9]+\\.[0-9a-zA-Z]+$")[0]==0, "The intermediate version string must be of the form build type major.minor[beta]. Had " + versionString);

	//-----------------------------------------------------------------------------------------------------------------
	// TIME STAMPS
	//-----------------------------------------------------------------------------------------------------------------

	GetString (timeStamp, TIME_STAMP, 0);
	assert ((timeStamp$"^[0-9]{4}/[0-9][0-9]?/[0-9][0-9]?\\ [0-9][0-9]?\\:[0-9][0-9]?$")[0]==0, "The GMT version of the time stamp must be of the form YYYY/M/D H:M. Had " + timeStamp);

	GetString (timeStamp, TIME_STAMP, 1);
	assert (Type (timeStamp) == "String", "The local version of the time stamp must be a string. Had " + Type (timeStamp));

	//-----------------------------------------------------------------------------------------------------------------
	// DATA SET 
	//-----------------------------------------------------------------------------------------------------------------

	DataSet 			_testDataSet = ReadDataFile ("../../data/mtDNA.fas");

	GetString (sequenceNames, _testDataSet, -1);
	assert (Type(sequenceNames) == "Matrix" && Type (sequenceNames[0]) == "String" && Rows(sequenceNames) == 1 && Columns (sequenceNames) == _testDataSet.species, "Retrieve all sequence names from a DataSet");
	
	GetString (sequenceName, _testDataSet, 2);
	assert (sequenceNames[2] == sequenceName, "Retrieve a sequence name from DataSet by index");

	GetString (sequenceName, _testDataSet, 1024);
	assert (Type (sequenceName) == "Unknown", "Retrieve an invalid sequence index from a DataSet");

	//-----------------------------------------------------------------------------------------------------------------
	// DATA SET FILTER
	//-----------------------------------------------------------------------------------------------------------------

	DataSetFilter 		_testDataSetFilter = CreateFilter (_testDataSet,1,"",""+(_testDataSet.species-1)+",1,0");

	GetString (filterNames, _testDataSetFilter, -1);
	assert (Type(filterNames) == "Matrix" && Type (filterNames[0]) == "String" && Rows(filterNames) == 1 && Columns (filterNames) == _testDataSetFilter.species, "Retrieve all sequence names from a DataSetFilter");

	GetString (filterSeqName, _testDataSetFilter, 0);
	assert (sequenceNames[_testDataSet.species-1] == filterSeqName, "Retrieve a sequence name from DataSetFilter by index");

	GetString (filterSeqName, _testDataSet, 1024);
	assert (Type (filterSeqName) == "Unknown", "Retrieve an invalid sequence index from a DataSet");

	//-----------------------------------------------------------------------------------------------------------------
	// VARIABLES
	//-----------------------------------------------------------------------------------------------------------------

	testVariable = Log(2);
	
	GetString (variableInfo, testVariable, 0);
	assert (Type (variableInfo) == "String" && variableInfo / "0.693*", "Retrieve the string representation of the value stored in a variable");
	
	global		  Z = 0;
	testVariable := Z * (2+Y);
	
	GetString (variableInfo, testVariable, 0);
	assert (Type (variableInfo) == "String" && variableInfo == "Z*(2+Y)", "Retrieve the expression for the constraint");

	GetString (variableInfo, testVariable, -1);
	
	assert (Type (variableInfo) == "AssociativeList" && (variableInfo["Global"])[0] == "Z" && (variableInfo["Local"])[0] == "Y", "Retrieve the variables invovled in the constraint");

	//Topology T 			   = ((a:0.1,b:0.2):0.4,c:0.15,d:0.33);

	
	
	testResult = 1;
		
	return testResult;
}

/* execution stub */

fprintf    (stdout, "[Running COVERAGE TEST '", getTestName(), "']\n");
result  =  runTest();

if (result)
{
	fprintf (stdout, "[TEST PASSED]\n");
}
else
{
	fprintf (stdout, "[TEST FAILED]\n");
}