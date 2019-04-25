/* 

THIS HyPhy BATCH LANGUAGE FILE IS PROVIDED AS AN EXAMPLE/TOOL FOR 
THE 2nd EDITION OF 'The Phylogenetic Handbook'

Written by 
Sergei L Kosakovsky Pond (spond@ucsd.edu)
Art FY Poon				 (apoon@ucsd.edu)
Simon DW Frost			 (sdfrost@ucsd.edu)

*/


RequireVersion ("0.9920060830");
VERBOSITY_LEVEL = -1;

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

fprintf (stdout, "\n\nThis example analysis fits the local & global Muse-Gaut 94 model with the \nuser selected nucleotide correction to a small codon dataset",
				 "\nestimates omega for each branch and reports",
				 "\nconfidence intervals and sampling properties for each.",
				 "\nThe analysis also performs a LRT for branch-by-branch rate variation.\n");
				 
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
		GetURL			(dataString, "http://www.hyphy.org/phylohandbook/data/Drosophila_adh.nex");	
	}
	else
	{
		GetURL			(dataString, "http://www.hyphy.org/phylohandbook/data/H5N1_HA_5.nex");
	}	
	DataSet ds		 = ReadFromString (dataString);
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

brNames				= BranchName (givenTree,-1);
COVARIANCE_PARAMETER 				= {};
global global_OMEGA = 1;
for (k=0; k < Columns (brNames)-1; k=k+1)
{
	ExecuteCommands ("givenTree."+brNames[k]+".nonSynRate:=givenTree."+brNames[k]+".omega*givenTree."+brNames[k]+".synRate;");
	COVARIANCE_PARAMETER["givenTree."+brNames[k]+".omega"] = 1;
}

LikelihoodFunction  theLnLik = (filteredData, givenTree);

for (k=0; k < Columns (brNames)-1; k=k+1)
{
	ExecuteCommands ("givenTree."+brNames[k]+".omega:=global_OMEGA;");
}

fprintf 					   (stdout, "Fitting the global model to the data...\n");
Optimize 					   (res_global, theLnLik);
fprintf						   (stdout, theLnLik,"\n\n");

for (k=0; k < Columns (brNames)-1; k=k+1)
{
	ExecuteCommands ("givenTree."+brNames[k]+".omega=global_OMEGA;");
}

fprintf 					   (stdout, "Fitting the local model to the data...\n");
Optimize 					   (res_local, theLnLik);
fprintf						   (stdout, theLnLik,"\n\n");

LR = 2(res_local[1][0]-res_global[1][0]);
DF = res_local[1][1]-res_global[1][1];

fprintf (stdout, "\nLRT for variable omega across the tree\n\tLR = ", LR, "\n\tConstaints = ", DF, "\n\tp-value = ", 1-CChi2(LR,DF),"\n\n");

COVARIANCE_PRECISION = 0.95;
CovarianceMatrix (covMx, theLnLik);

VERBOSITY_LEVEL = 0;

for (k=0; k < Columns (brNames)-1; k=k+1)
{
	fprintf (stdout, "Branch :", brNames[k], "\n\tomega MLE = ", Format (covMx[k][1],6,3), "\n\t95% CI = (",Format (covMx[k][0],6,3), ",", Format (covMx[k][2],6,3), ")\n");
}