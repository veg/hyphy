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

internalNodes = BranchCount(givenTree);
leafNodes     = TipCount(givenTree);
choiceMatrix  = {internalNodes+leafNodes,2};

for (bc=0; bc<internalNodes; bc=bc+1)
{
	choiceMatrix[bc][0] = BranchName(givenTree,bc);
	choiceMatrix[bc][1] = "Internal Branch Rooting " + givenTree[bc];
}

for (bc=0; bc<leafNodes; bc=bc+1)
{
	choiceMatrix[bc+internalNodes][0] = TipName(givenTree,bc);
	choiceMatrix[bc+internalNodes][1] = "Terminal branch endin in " + choiceMatrix[bc+internalNodes][0];
}

ChoiceList  (bOption,"Choose the branch(es) to test",0,NO_SKIP,choiceMatrix);

if (bOption[0] < 0)
{
	return -1;
}

selectedBranches = {};
for (k=0; k<Columns(bOption); k=k+1)
{
	selectedBranches[choiceMatrix[bOption[k]][0]] = 1;
	fprintf (stdout, "Selected branch ", choiceMatrix[bOption[k]][0], " as foreground\n");
}
selectedBranchesKey = Rows (selectedBranches);

brNames								= BranchName (givenTree,-1);

global global_OMEGA 				= 1;

ClearConstraints (givenTree);

ReplicateConstraint 		   ("this1.?.nonSynRate:=this2.?.synRate*global_OMEGA", givenTree, givenTree);
LikelihoodFunction  theLnLik = (filteredData, givenTree);

fprintf 					   (stdout, "Fitting the global model to the data...\n");
Optimize 					   (res_global, theLnLik);
fprintf						   (stdout, "Estimated tree-wide omega to be ", global_OMEGA, "\n");

for (k=0; k < Columns (brNames)-1; k=k+1)
{
	if (selectedBranches[brNames[k]] == 1)
	{
		ExecuteCommands ("ClearConstraints (givenTree."+brNames[k]+");");
	}
}

fprintf 					   (stdout, "\nFitting the a priori branch model (with uniform background) to the data\n");
Optimize 					   (res_apriori, theLnLik);

fprintf						   (stdout, "Estimated background omega to be ", global_OMEGA, "\n");
for (k=0; k < Columns (selectedBranchesKey); k=k+1)
{
	reportBranchOmega (selectedBranchesKey[k]);
}

LR = 2*(res_apriori[1][0]-res_global[1][0]);
DF = res_apriori[1][1]-res_global[1][1];

fprintf (stdout, "LRT test (foreground is different from background)\n\tH_0 LogL   : ", Format (res_global[1][0],8,3),
						 "\n\tH_A LogL   : ", Format (res_apriori[1][0],8,3),
				         "\n\tLR         : ", Format (LR,8,3), 
				         "\n\tConstraints: ", Format (DF,8,0), 
				         "\n\tp-value    : ", Format(1-CChi2(LR,DF),8,3), "\n");


ClearConstraints			   (givenTree);
for (k=0; k < Columns (selectedBranchesKey); k=k+1)
{
	ExecuteCommands ("givenTree."+selectedBranchesKey[k]+".nonSynRate:=givenTree."+selectedBranchesKey[k]+".omega*givenTree."+selectedBranchesKey[k]+".synRate");
	ExecuteCommands ("givenTree."+selectedBranchesKey[k]+".omega:<1;");
	
}

fprintf 					   (stdout, "\nFitting the a priori branch model (non-uniform background) to the data\n");
LikelihoodFunction  theLnLik = (filteredData, givenTree);
Optimize 					   (res_apriori2, theLnLik);


ClearConstraints			   (givenTree);
LikelihoodFunction  theLnLik = (filteredData, givenTree);
fprintf 					   (stdout, "Fitting the local model to the data\n");
Optimize 					   (res_local, theLnLik);

LR = 2*(res_local[1][0]-res_apriori2[1][0]);
DF = Columns (selectedBranchesKey);

for (k=0; k < Columns (selectedBranchesKey); k=k+1)
{
	reportBranchOmega (selectedBranchesKey[k]);
}

if (DF == 1)
{
	fprintf (stdout, "LRT test (positive selection on foreground branches)\n\tH_0 LogL   : ", Format (res_apriori2[1][0],8,3),
							 "\n\tH_A LogL   : ", Format (res_local[1][0],8,3),
							 "\n\tLR         : ", Format (LR,8,3), 
							 "\n\tConstraints: ", Format (DF,8,0), " (one sided)",
							 "\n\tp-value    : ", Format(0.5(1-CChi2(LR,DF)),8,3), "\n");
}
else
{
	fprintf (stdout, "LRT test\n\tH_0 LogL   : ", Format (res_apriori2[1][0],8,3),
							 "\n\tH_A LogL   : ", Format (res_local[1][0],8,3),
							 "\n\tLR         : ", Format (LR,8,3), 
							 "\n\tConstraints: ", Format (DF,8,0), " (one sided)",
							 "\n\tp-value    : ", Format(1-CChi2(LR,DF),8,3), " (conservative test)\n");
}

VERBOSITY_LEVEL = 0;

/*---------------------------------------------------------------------------------------------------------------------------------------------*/
function reportBranchOmega (myBrName)
{
	ExecuteCommands ("v1=givenTree."+myBrName+".nonSynRate;v2=givenTree."+myBrName+".synRate");
	if (v2>0)
	{
		fprintf						   (stdout, "Estimated branch omega for ", myBrName, " to be ", v1/v2, "\n");
	}
	else
	{
		if (v1>0)
		{
			fprintf						   (stdout, "Estimated branch omega for ", myBrName, " to be infinite (alpha=0)\n");	
		}
		else
		if (v1>0)
		{
			fprintf						   (stdout, "Estimated branch omega for ", myBrName, " to be indeterminate (alpha and beta =0)\n");	
		}
	}
	return 0;
}