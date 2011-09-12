bivariateFitOptions = {};

SetDialogPrompt ("Datafile:");
DataSet ds = ReadDataFile (PROMPT_FOR_FILE);
_bivariateFilePath = LAST_FILE_PATH;

ExecuteAFile (HYPHY_LIB_DIRECTORY + 
				"TemplateBatchFiles" + 
				DIRECTORY_SEPARATOR  + 
				"TemplateModels" + 
				DIRECTORY_SEPARATOR  + 
				"chooseGeneticCode.def");

_genCodeID = _geneticCodeOptionMatrix [modelType][0];


ChoiceList (branchLengths,"Branch Lengths",1,SKIP_NONE,
            "Codon Model","Jointly optimize rate parameters and branch lengths (slow and thorough)",
            "Nucleotide Model","Estimate branch lengths once, using an appropriate nucleotide model (quick and dirty)."
            );

assert (branchLengths >= 0);

fprintf (stdout, "Model string:");
fscanf  (stdin,  "String",mstring);

ExecuteAFile ("TemplateModels/MGwAA.ibf");

currentLFSpool = _bivariateFilePath + ".fit.1";

bivariateFitOptions ["00"] = "New run";
if (branchLengths == 0)
{
    bivariateFitOptions ["01"] = "Codon Model";
}
else
{
    bivariateFitOptions ["01"] = "Nucleotide Model";
}
bivariateFitOptions ["02"] = "1";
bivariateFitOptions ["03"] = _genCodeID;
bivariateFitOptions ["04"] = _bivariateFilePath;
if (Abs(DATAFILE_TREE))
{
	bivariateFitOptions ["06"] = "y";
}
else
{
	bivariateFitOptions ["06"] = "";
}

bivariateFitOptions ["07"] = mstring;
bivariateFitOptions ["08"] = "Default";
bivariateFitOptions ["09"] = "1";
bivariateFitOptions ["10"] = "Unconstrained";
bivariateFitOptions ["11"] = currentLFSpool;

ExecuteAFile (  "dNdSBivariateRateAnalysis.bf",
				bivariateFitOptions
			 );
			 
currentCAIC   = bivariateReturnAVL["cAIC"];
bestCAICsoFar = 1e100;
currentRateCount = 1;
bestClassCount = 1;

gateauxOptions = {};
gateauxOptions ["02"] = "-1";


while (currentCAIC < bestCAICsoFar)
{
	bestClassCount = currentRateCount;
	gateauxOptions ["01"] = currentLFSpool;
	fprintf (stdout, "[FINISHED FITTING A MODEL WITH ", currentRateCount, " RATES]\n");
	fprintf (stdout, "[CURRENT c-AIC = ", currentCAIC, ". BEST c-AIC SO FAR = ", bestCAICsoFar, "]\n");
	bestCAICsoFar = currentCAIC;

	DeleteObject (lf);

	ExecuteAFile (  HYPHY_LIB_DIRECTORY + 
				"TemplateBatchFiles" + 
				DIRECTORY_SEPARATOR  + 
				"GateauxMR.bf",
				gateauxOptions
			 );

	improved = bivariateReturnAVL ["DIFF"];
	if (improved == 0)
	{
		fprintf (stdout, "[NO GATEAUX IMPROVEMENT. TERMINATING...]\n");
		break;
	}
	lastLFSpool		 = currentLFSpool;
	currentRateCount = currentRateCount + 1;
	currentLFSpool = currentLFSpool + "." + currentRateCount;

	bivariateFitOptions ["00"] = "Continued run";
	bivariateFitOptions ["01"] = currentLFSpool;
	bivariateFitOptions ["02"] = "Unconstrained";

	ExecuteAFile (  HYPHY_LIB_DIRECTORY + 
				"TemplateBatchFiles" + 
				DIRECTORY_SEPARATOR  + 
				"dNdSBivariateRateAnalysis.bf",
				bivariateFitOptions
			 	);
			 
	currentCAIC   = bivariateReturnAVL["cAIC"];
	
	
}

fprintf (stdout, "[BEST-FITTING MODEL HAS ", bestClassCount, " RATES]\n");
fprintf (stdout, "[USE CodonBivariateRateProcessor.bf ON ", lastLFSpool, " TO PROCESS THE RESULTS]\n");





