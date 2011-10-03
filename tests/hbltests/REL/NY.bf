/* test preamble */

	_testDescription 		= " M0,M1,M2 PAML style fit out an HIV env alignment with 13 sequences";
	_expectedLL = 			{{-1137.68873768577,-1114.64174822308,-1106.44533632047}};
	ExecuteAFile 			("../Shared/TestInstrumentation.bf");
	startTestTimer 			(_testDescription);
	
/* end test preamble */

VERBOSITY_LEVEL			   = 1;
USE_ADAPTIVE_VARIABLE_STEP = 1;
OPTIMIZATION_METHOD        = 4;

fprintf (stdout, "\nRunning an \n");
runTimer = Time (1);

inputOptions = {};
inputOptions ["00"] = "Universal";
inputOptions ["01"] = PATH_TO_CURRENT_BF + ".." + DIRECTORY_SEPARATOR + "data" + DIRECTORY_SEPARATOR + "HIVenvSweden.seq"; 
inputOptions ["02"] = "y";
inputOptions ["03"] = "Run Custom";
inputOptions ["04"] = "Single Rate";
inputOptions ["05"] = "Neutral";
inputOptions ["06"] = "Selection";
inputOptions ["07"] = "";
inputOptions ["08"] = "Default";
inputOptions ["09"] = "GY94 3x4";
inputOptions ["10"] = "0.9";
inputOptions ["11"] = PATH_TO_CURRENT_BF + ".." + DIRECTORY_SEPARATOR + "Results" + DIRECTORY_SEPARATOR + "HIVSweden.out";

ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "NielsenYang.bf", inputOptions);


fittedLL = {1,3};
fittedLL[0] = modelLL[-1];
fittedLL[1] = modelLL[0];
fittedLL[2] = modelLL[1];

/* test epilogue */
	timeMatrix = endTestTimer 				  (_testDescription);
	if (logTestResult (Abs (fittedLL - _expectedLL) < 2*OPTIMIZATION_PRECISION))
	{
		return timeMatrix;
	}
	return 0;
/* end test epilogue */



