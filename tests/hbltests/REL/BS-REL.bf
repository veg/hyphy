/* test preamble */

	_testDescription 		= "Branch-site REL on the alignment of 10 CD2 sequences";
	_expectedLL = 			{{-3450.39628591,-3410.18927}};
	ExecuteAFile 			("../Shared/TestInstrumentation.bf");
	startTestTimer 			(_testDescription);
	
/* end test preamble */

VERBOSITY_LEVEL			   = 1;
USE_ADAPTIVE_VARIABLE_STEP = 1;
OPTIMIZATION_METHOD        = 4;

runTimer = Time (1);

inputOptions = {};
inputOptions ["00"] = "Universal";
inputOptions ["01"] = PATH_TO_CURRENT_BF + ".." + DIRECTORY_SEPARATOR + "data" + DIRECTORY_SEPARATOR + "CD2.nex"; 
inputOptions ["02"] = "y";
inputOptions ["03"] = PATH_TO_CURRENT_BF + ".." + DIRECTORY_SEPARATOR + "Results" + DIRECTORY_SEPARATOR + "CD2.bsrel.out";



//ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "BranchSiteREL.bf", inputOptions);
ExecuteAFile ("../../../res/TemplateBatchFiles/BranchSiteREL.bf", inputOptions);

fittedLL = {1,2};
fittedLL[0] = res_base[1][0];
fittedLL[1] = res_three_LF[1][0];



/* test epilogue */
	timeMatrix = endTestTimer 				  (_testDescription);
	if (logTestResult (Abs (fittedLL - _expectedLL) < 2*OPTIMIZATION_PRECISION))
	{
		return timeMatrix;
	}
	return 0;
/* end test epilogue */



