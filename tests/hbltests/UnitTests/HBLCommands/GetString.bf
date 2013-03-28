ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();

function getTestName ()
{
	return "GetString";
}		

function getTestedFunctions ()
{
	return {{"_ExecutionList::HandleGetString"}};
}	

function factorial (N) {
    f = 1; 
    for (k = 2; k <= N; k+=1) {
        f = f * k;
    }
    return f;
}

function runTest ()
{
	/* define some auxiliary objects here */
	
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult  	   = 0;

	//-----------------------------------------------------------------------------------------------------------------
	// ERROR HANDLING
	//-----------------------------------------------------------------------------------------------------------------
	
    assert (runCommandWithSoftErrors ("GetString (a+b, HYPHY_VERSION, 0)", " is not a valid variable identifier in call to GetString"), "Failed error checking for an invalid receptacle");
    assert (runCommandWithSoftErrors ("GetString (data, this_object_better_not_exist, 0)", " is not an allowed argument type"), "Failed error checking for an invalid argument");

	
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

	//-----------------------------------------------------------------------------------------------------------------
	// HBL User Function
	//-----------------------------------------------------------------------------------------------------------------

	GetString (functionInfo, factorial, -1);
	assert (Type (functionInfo) == "AssociativeList", "Retrieve information about a user function");
	ExecuteCommands ("function " + functionInfo["ID"] + "(" + Join(",", functionInfo["Arguments"]) + ") {" + functionInfo["Body"] + "} fact5 = factorial (5);");
	assert (fact5 == 120, "HBL function component retrieval 5! != " + fact5);
	
	//-----------------------------------------------------------------------------------------------------------------
	// SUBSTITUTION MODELS
	//-----------------------------------------------------------------------------------------------------------------

    global R = 1;
    binaryModel = {{*,rate1}{R*rate2,*}};
    freqs       = {{0.5,0.5}};
    Model         testModel1 = (binaryModel, freqs, 0);
    Model         testModel2 =("Exp(binaryModel)",freqs,EXPLICIT_FORM_MATRIX_EXPONENTIAL);

    GetString (modelInfo1, testModel1,       1);
    GetString (modelInfo2, testModel1,       2);
    GetString (modelBL, testModel1,         -1);
    GetString (modelBits, testModel1,       -2);
    GetString (modelInfo1_0, testModel1,    1,0);

 
	assert (Type (modelInfo2) == "Unknown" && modelInfo1 == "rate2" && modelBL == "0.5*rate1+0.5*R*rate2" && modelInfo1_0 == "R*rate2" &&
	        Type (modelBits) == "AssociativeList" && modelBits["EQ_FREQS"] == "freqs", "Retrieve information about a substitution model");

    GetString (modelInfo1, testModel2,       0);
    GetString (modelBL, testModel2,         -1);
    GetString (modelBits, testModel2,       -2);
    GetString (modelInfo1_0, testModel2,    1,0);

	assert ( modelInfo1 == "rate1" && Type (modelInfo2) == "Unknown" && Type (modelBL) == "Unknown" && Type (modelInfo1_0) == "Unknown" , "Retrieve information about a substitution model in explicit form");

	//-----------------------------------------------------------------------------------------------------------------
	// Likelihood Functions
	//-----------------------------------------------------------------------------------------------------------------

    ExecuteAFile (PATH_TO_CURRENT_BF  + "res" + DIRECTORY_SEPARATOR + "test_likefunc.nex");
    ExecuteAFile (PATH_TO_CURRENT_BF  + "res" + DIRECTORY_SEPARATOR + "test_likefunc2.nex");

    GetString (likelihoodFunctionInfoArgument, lf, 10);
    GetString (likelihoodFunctionInfoArgumentDep, lf, 18);
    GetString (likelihoodFunctionInfoArgumentFail, lf, 65536);
	assert (Type (likelihoodFunctionInfoArgumentFail) == "Unknown" && likelihoodFunctionInfoArgument == "givenTree.B_US_83_RF_ACC_M17451.synRate" &&
	              likelihoodFunctionInfoArgumentDep == "GT", "Retrieve information about likelihood function parameters");

    GetString (likelihoodFunctionInfoArgumentTotal, lf, -1);

 	//-----------------------------------------------------------------------------------------------------------------
	// Global Retrieval by Type
	//-----------------------------------------------------------------------------------------------------------------


    GetString (objectID0, LikelihoodFunction, 0);
    GetString (objectID1, DataSet, 1);
    GetString (objectID2, DataSetFilter, 2000);
    GetString (objectInfo, UserFunction, 0);
    GetString (objectInfoInvalid, UserFunction, 5000);
    GetString (treeInfo, Tree, 1);
    GetString (treeInfoInvalid, Tree, 8192); 
    DeleteObject (lf);
    GetString (objectID0_2, LikelihoodFunction, 0);

	assert (Type (objectID2) == "Unknown" && objectID0 == "lf" && objectID1 == "ds" && Type (objectInfo) == "AssociativeList" && Type (objectInfoInvalid) == "Unknown" &&
	        objectID0_2 == "IntermediateCodon_AA_LF", "Retrieve identifiers of HBL objects by index");


	//-----------------------------------------------------------------------------------------------------------------
	// SCFG
	//-----------------------------------------------------------------------------------------------------------------
    
    ExecuteAFile (PATH_TO_CURRENT_BF + "res" + DIRECTORY_SEPARATOR + "SCFG" + DIRECTORY_SEPARATOR + "scfgG6c.bf", {"0": "small.txt"});
    GetString (scfgID, SCFG, 0);
    assert (scfgID == "G6", "Retrieve an SCFG identifier");
    ExecuteCommands ("GetString (scfgInfoTotal, `scfgID`, -1)");
    assert ((scfgInfoTotal["TERMINALS"])[2] == "<", "Retrieve SCFG components");

	testResult = 1;
		
	return testResult;
}
