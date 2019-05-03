/* 

THIS HyPhy BATCH LANGUAGE FILE IS PROVIDED AS AN EXAMPLE/TOOL FOR 
THE 2nd EDITION OF 'The Phylogenetic Handbook'

Written by 
Sergei L Kosakovsky Pond (spond@ucsd.edu)
Art FY Poon				 (apoon@ucsd.edu)
Simon DW Frost			 (sdfrost@ucsd.edu)

*/


//RequireVersion ("2.220140505");
VERBOSITY_LEVEL = -1;

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

fprintf (stdout, "\n\nThis example analysis fits the local Muse-Gaut 94 model with the \nuser selected nucleotide correction to a small codon dataset",
				 "\nestimates omega for each branch and reports",
				 "\nconfidence intervals and sampling properties for each",
				 "\nbased on asymptotic normality, profile likelihood and",
				 "\nsampling importance resampling\n");

ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Utility"+DIRECTORY_SEPARATOR+"DescriptiveStatistics.bf");
ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"chooseGeneticCode.def");

ModelMatrixDimension = 64;
for (k=0; k<64; k=k+1)
{
	if (_Genetic_Code[k] == 10)
	{
		ModelMatrixDimension = ModelMatrixDimension -1;
	}
}

ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"2RatesAnalyses"+DIRECTORY_SEPARATOR+"MG94xREV.mdl");

ChoiceList (freqs, "Alignment", 1, SKIP_NONE, "Flu H5N1 HA", 		"Use the example Influenza H5N1 heamagglutinin alignment (5 sequences)",
											  "Drospophila adh", 	"Use the example Drosophila ADH alignment (6 sequences).",
											  "Custom", 			"Load your own coding alignment.");
													 
if (freqs < 0)
{
	return 0;
}
													 
if (freqs == 2)
{
	SetDialogPrompt     ("Choose a nucleotide alignment");
	DataSet ds        = ReadDataFile (PROMPT_FOR_FILE);
}
else
{
	if (freqs == 1)
	{
		DataSet ds = ReadDataFile (PATH_TO_CURRENT_BF + '/datasets/Drosophilia_adh.nex');
	}
	else
	{
		DataSet ds = ReadDataFile (PATH_TO_CURRENT_BF + '/datasets/H5N1_HA_5.nex');
	}
}

			  			  			  			  
DataSetFilter	  	filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);

SKIP_MODEL_PARAMETER_LIST = 1;
done 					  = 0;

while (!done)
{
	fprintf (stdout,"\nPlease enter a 6 character model designation (e.g:010010 defines HKY85):");
	fscanf  (stdin,"String", modelDesc);
	if (Abs(modelDesc)==6)
	{	
		done = 1;
	}
}	
modelType 				  = 0;
ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"MG94custom.mdl");
SKIP_MODEL_PARAMETER_LIST = 0;

ExecuteAFile 		(HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"queryTree.bf");

brNames	= BranchName (givenTree,-1);
cpm = {};
for (k=0; k < Columns (brNames)-1; k=k+1)
{
	ExecuteCommands ("givenTree."+brNames[k]+".nonSynRate:=givenTree."+brNames[k]+".omega*givenTree."+brNames[k]+".synRate;");
	cpm["givenTree."+brNames[k]+".omega"] = 1;
}

LikelihoodFunction  theLnLik = (filteredData, givenTree);

fprintf 					   (stdout, "Fitting the local model to the data...\n");
Optimize 					   (res, theLnLik);
fprintf						   (stdout, theLnLik,"\n\n");

cpm2 = {};
for (k=0; k < Columns (brNames)-1; k=k+1)
{
	pName = "givenTree."+brNames[k]+".omega";
	ExecuteCommands ("pVal = "+pName+";");
	if (pVal < 1e-10)
	{
		fprintf (stdout, "Excluding ", pName, " because the value (",pVal,") is on the boundary of the parameter space...\n");
	}
	else
	{
		cpm2 [pName] = 1;
	}
}

fprintf 					   (stdout, "\nEstimating the variance-covariance matrix\n");
COVARIANCE_PARAMETER 			= cpm2;
COVARIANCE_PRECISION 			= 2;
CovarianceMatrix 				(cmx, theLnLik);

pNames = Rows (cpm2);
for (k=0; k < Columns(pNames); k=k+1)
{
	pName = pNames[k];
	ExecuteCommands ("pVal = "+pName+";");
	fprintf (stdout, pName, "\n\tMLE = ", Format(pVal,10,3), ", Sampling Variance = ", Format(cmx[k][k],10,3), ", 95% CI: (",Format(Max(0,pVal-1.96*Sqrt(cmx[k][k])),10,3),",",Format(pVal+1.96*Sqrt(cmx[k][k]),10,3),"),\n");
}


COVARIANCE_PARAMETER = cpm;
SAMPLE_N 			 = 1000;
SAMPLE_M 			 = 100;
LF_NAME  			 = "theLnLik";
VERBOSITY_LEVEL 	 = 0;

ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Samplers"+DIRECTORY_SEPARATOR+"lhc-ErrorEst.bf");

pNames = Rows (cpm2);

fprintf (stdout, "\nConfidence Intervals using profile likelihood and sampling - resampling\n");

pNames = Rows (cpm);
for (k=0; k < Columns(pNames); k=k+1)
{
	pName = pNames[k];
	fprintf (stdout, pName, "\n\tMLE = ", Format (covMx[k][1],10,3), " 95% CI: (",Format(covMx[k][0],10,3),",",Format(covMx[k][2],10,3),")");
	dsc = GatherDescriptiveStats (sampledPoints[-1][k+1]);
	fprintf (stdout, "\n\tSampled mean = ", Format (dsc["Mean"],10,3), ", variance ", Format (dsc["Variance"],10,3)," and 95% CI: (",Format(dsc["2.5%"],10,3),",",Format(dsc["97.5%"],10,3),")\n");
}