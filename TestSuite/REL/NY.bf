fprintf (stdout, "\nRunning an M0,M1,M2 PAML style fit out an HIV alignment with 13 sequences\n");
runTimer = Time (1);

inputOptions = {};
inputOptions ["00"] = "Universal";
inputOptions ["01"] = PATH_TO_CURRENT_BF + "HIVenvSweden.seq"; 
inputOptions ["02"] = "y";
inputOptions ["03"] = "Run Custom";
inputOptions ["04"] = "Single Rate";
inputOptions ["05"] = "Neutral";
inputOptions ["06"] = "Selection";
inputOptions ["07"] = "";
inputOptions ["08"] = "Default";
inputOptions ["09"] = "GY94 3x4";
inputOptions ["10"] = "0.9";
inputOptions ["11"] = PATH_TO_CURRENT_BF + "Results" + DIRECTORY_SEPARATOR + "HIVSweden.out";

ExecuteAFile (HYPHY_BASE_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "NielsenYang.bf", inputOptions);

expectedLL = {{-1137.68873768577,-1114.64174822308,-1106.44533632047}};

fprintf (stdout, "\nThe analysis took ", Time (1)-runTimer, " seconds.\n", diffLL , " difference between obtained and expected likelihood\n\nAbs(expected LL - obtained LL)\n");

for (k = 0; k < 3; k = k+1)
{
	fprintf (stdout, "Model ", k+1, " : ", Abs(expectedLL[k] - modelLL[k-1]), "\n");
}


