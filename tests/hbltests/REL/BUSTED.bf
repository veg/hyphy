/* test preamble */

	_testDescription 		= "Whole tree BUSTED on the alignment of 10 CD2 sequences";
	_expectedLL = 			{{-3431.21988,-3438.42464}};
	ExecuteAFile 			("../Shared/TestInstrumentation.bf");
	startTestTimer 			(_testDescription);

/* end test preamble */

VERBOSITY_LEVEL			   = 1;

runTimer = Time (1);

inputOptions = {
"00": "Universal",
"01": PATH_TO_CURRENT_BF + ".." + DIRECTORY_SEPARATOR + "data" + DIRECTORY_SEPARATOR + "CD2.nex",
"02": "Internal",
"03": ""
};


ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "SelectionAnalyses" + DIRECTORY_SEPARATOR + "BUSTED.bf", inputOptions);



/* test epilogue */
	timeMatrix = endTestTimer 				  (_testDescription);
	if (logTestResult (Abs (((_BUSTED_json["fits"])["Unconstrained model"])["log-likelihood"] - _expectedLL[0]) < 2*OPTIMIZATION_PRECISION)) {
	    if (logTestResult (Abs (((_BUSTED_json["fits"])["Constrained model"])["log-likelihood"] - _expectedLL[1]) < 2*OPTIMIZATION_PRECISION)) {
		    return timeMatrix;
		}
	}
	return 0;
/* end test epilogue */



