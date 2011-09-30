/* test preamble */

	_testDescription 		= " a series of model fits on two HIV-1 RT sequences (nucleotide, protein and codon)";
	expectedLL = {{-947.794202692758,-947.794927641445,-662.04841954136,-930.569844169857}};
	ExecuteAFile 			("../Shared/TestInstrumentation.bf");
	startTestTimer 			(_testDescription);

/* end test preamble */


modelLL    = {1,4};

runTimer = Time(1);

fprintf (stdout, "\nRunning a series of model fits on two HIV-1 RT sequences\n");

inputOptions = {};
inputOptions["0"] = PATH_TO_CURRENT_BF+".." +DIRECTORY_SEPARATOR+"data"+DIRECTORY_SEPARATOR+"2.fas";
inputOptions["1"] = "HKY85";
inputOptions["2"] = "Global";
inputOptions["6"] = PATH_TO_CURRENT_BF+".."+DIRECTORY_SEPARATOR+"data"+DIRECTORY_SEPARATOR+"2.tree";

ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "AnalyzeNucProtData.bf", inputOptions);

modelLL [0] = res[1][0];

inputOptions["2"] = "Global w/variation";
inputOptions["3"] = "General Discrete";
inputOptions["4"] = "2";
inputOptions["7"] = "Don't Display";

ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "AnalyzeNucProtData.bf", inputOptions);

modelLL [1] = res[1][0];

inputOptions["0"] = PATH_TO_CURRENT_BF+".."+DIRECTORY_SEPARATOR+"data"+DIRECTORY_SEPARATOR+"2.prot";
inputOptions["1"] = "HIVBETWEEN";
inputOptions["2"] = "Rate variation";

ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "AnalyzeNucProtData.bf", inputOptions);

modelLL [2] = res[1][0];


inputOptions["1"] = PATH_TO_CURRENT_BF+".."+DIRECTORY_SEPARATOR+"data"+DIRECTORY_SEPARATOR+"2.fas";
inputOptions["0"] = "Universal";
inputOptions["2"] = "MG94";
inputOptions["3"] = "Local";
inputOptions["4"] = inputOptions["6"];
inputOptions["5"] = inputOptions["7"];

ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "AnalyzeCodonData.bf", inputOptions);

modelLL [3] = res[1][0];


/* test epilogue */
	timeMatrix = endTestTimer 				  (_testDescription);
	if (logTestResult (Abs (expectedLL-modelLL) < 2*OPTIMIZATION_PRECISION))
	{
		return timeMatrix;
	}
	return 0;
/* end test epilogue */