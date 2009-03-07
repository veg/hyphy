#include "relratioBootstrap.bf";

RESTORE_GLOBALS = 1;

function RestoreGlobalValues (lfIndex)
{
	if (lfIndex==0)
	{
		for (i=0;i<SAVE_GLOBALS;i=i+1)
		{
			SetParameter (lf,i,globalSpoolMatrix[i]);
		}
	}
	if (lfIndex==1)
	{
		for (i=0;i<SAVE_GLOBALS;i=i+1)
		{
			SetParameter (lf2,i,globalSpoolMatrix1[i]);
		}
	}
	if (lfIndex==2)
	{
		for (i=0;i<SAVE_GLOBALS2;i=i+1)
		{
			SetParameter (lfConstrained,i,globalSpoolMatrix2[i]);
		}
	}
	return 0;
}

fprintf(stdout,"\n ---- RUNNING RELATIVE RATIO ANALYSIS ---- \n");

ChoiceList (dataType,"Data type",1,SKIP_NONE,"Nucleotide/Protein","Nucleotide or amino-acid (protein).",
				     "Codon","Codon (several available genetic codes).");

if (dataType<0) 
{
	return;
}
if (dataType)
{
	NICETY_LEVEL = 3;
	#include "TemplateModels/chooseGeneticCode.def";
}

SetDialogPrompt ("Choose the 1st data file:");

DataSet ds = ReadDataFile (PROMPT_FOR_FILE);

SetDialogPrompt ("Choose the 2nd data file:");

treeRead = 0;

if (IS_TREE_PRESENT_IN_DATA)
{
	treeString = DATAFILE_TREE;
	treeRead = 1;
}

DataSet ds2 = ReadDataFile (PROMPT_FOR_FILE);

while (ds2.species<ds1.species)
{
	SetDialogPrompt ("Too few taxa in 2nd file. Choose again:");
	DataSet ds2 = ReadDataFile (PROMPT_FOR_FILE);
}

fprintf (stdout,"The following data was read:\n",ds,"\n", ds2,"\n");

if (dataType)
{
	DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);
}
else
{
	DataSetFilter filteredData = CreateFilter (ds,1);
}

SelectTemplateModel(filteredData);

if (treeRead)
{
	fprintf (stdout, "\n\nA tree was found in the data file:\n",treeString,"\n\nWould you like to use it:(Y/N)?");
	fscanf (stdin, "String", response);
	if ((response=="n")||(response=="N"))
	{
		treeRead = 0;
	}
	fprintf (stdout, "\n\n");
}

if (!treeRead)
{
	SetDialogPrompt ("Please select a tree file for the data:");

	fscanf (PROMPT_FOR_FILE, "String", treeString);
	
}

global RelRatio = 1.0;

relationString = ":=RelRatio*";

#include "selectModelParameters.bf";

Tree tr1  = treeString;
Tree tr1c = treeString;

LikelihoodFunction lf = (filteredData,tr1);

Optimize (res,lf);

SAVE_GLOBALS = res[1][2];

if (SAVE_GLOBALS)
{
	globalSpoolMatrix = {1,SAVE_GLOBALS};
	for (i=0;i<SAVE_GLOBALS;i=i+1)
	{
		globalSpoolMatrix[i]=res[0][i];
	}
}

separator = "*-----------------------------------------------------------*";

fprintf (stdout, "\n", separator, "\nFULL MODEL RESULTS:\nDataset 1:",lf);

fullModelLik = res[1][0];

fullVars = 2*res[1][1];

if (dataType)
{
	DataSetFilter filteredData2 = CreateFilter (ds2,3,"",speciesIndex<ds.species,GeneticCodeExclusions);
}
else
{
	DataSetFilter filteredData2 = CreateFilter (ds2,1,"",speciesIndex<ds.species);
}

HarvestFrequencies (vectorOfFrequencies2,filteredData2,1,1,0);

if (FREQUENCY_SENSITIVE)
{
	modelMatrix2 = 0;
	if (USE_POSITION_SPECIFIC_FREQS)
	{
		HarvestFrequencies (vectorOfFrequencies2,filteredData2,3,1,1);
	}
	MULTIPLY_BY_FREQS = PopulateModelMatrix ("modelMatrix2",vectorOfFrequencies2);
	if (dataType)
	{
		CodonEFV2 = BuildCodonFrequencies (vectorOfFrequencies2);
		Model model2 = (modelMatrix2,CodonEFV2,MULTIPLY_BY_FREQS);
	}
	else
	{
		Model model2 = (modelMatrix2,vectorOfFrequencies2,MULTIPLY_BY_FREQS);
	}	
}

Tree tr2  = treeString;
Tree tr2c = treeString;

LikelihoodFunction lf2 = (filteredData2,tr2);

Optimize (res2,lf2);

if (SAVE_GLOBALS)
{
	globalSpoolMatrix1 = {1,SAVE_GLOBALS};
	for (i=0;i<SAVE_GLOBALS;i=i+1)
	{
		globalSpoolMatrix1[i]=res2[0][i];
	}
}

fullModelLik = fullModelLik+res2[1][0];

fullModelParam = 2*res[1][1];

fprintf (stdout, "\n\nDataset 2:",lf2, "\n\nJoint ln-likelihood=",fullModelLik);


/* now specify the constraint */

if (FREQUENCY_SENSITIVE)
{
	modelMatrix1 = 0;
	MULTIPLY_BY_FREQS = PopulateModelMatrix ("modelMatrix1",vectorOfFrequencies);
	if (dataType)
	{
		CodonEFV = BuildCodonFrequencies (vectorOfFrequencies);
		Model model1 = (modelMatrix1,CodonEFV,MULTIPLY_BY_FREQS);
	}
	else
	{
		Model model1 = (modelMatrix1,vectorOfFrequencies,MULTIPLY_BY_FREQS);
	}	
}

ReplicateConstraint (constraintString,tr2c,tr1c);

/*
mxTreeSpec = {5,1};
mxTreeSpec [0] = "tr2c";
mxTreeSpec [1] = "8240";
mxTreeSpec [2] = "10,40,-10,-175,1";
mxTreeSpec [3] = "";
mxTreeSpec [4] = "";
OpenWindow (TREEWINDOW, mxTreeSpec);							
internalNodes = BranchCount(tr2c);
leafNodes	  = TipCount (tr2c);
choiceMatrix = {leafNodes+internalNodes,2};
for (bc=0; bc<leafNodes; bc=bc+1)
{
	choiceMatrix[bc][0] = TipName(tr2c,bc);
	choiceMatrix[bc][1] = "Taxon " + choiceMatrix[bc][0];
}
for (bc=0; bc<internalNodes; bc=bc+1)
{
	choiceMatrix[bc+leafNodes][0] = BranchName(tr2c,bc);
	choiceMatrix[bc+leafNodes][1] = "Internal Branch Rooting " + tr2c[bc];
}
ChoiceList  (stOption,"Do not constrain:",0,NO_SKIP,choiceMatrix);

if (stOption[0] >= 0)
{
	for (k=0; k<Columns (stOption) * Rows (stOption); k=k+1)
	{
		ExecuteCommands ("ClearConstraints (tr2c."+choiceMatrix[(stOption[k])][0]+");");
	}
}
*/

LikelihoodFunction lfConstrained = (filteredData2,tr2c,filteredData,tr1c);

Optimize (res1,lfConstrained);

SAVE_GLOBALS2 = res1[1][2];

globalSpoolMatrix2 = {1,SAVE_GLOBALS2};

for (i=0;i<SAVE_GLOBALS2;i=i+1)
{
	globalSpoolMatrix2[i]=res1[0][i];
}

fprintf (stdout, "\n", separator,"\n\nCONSTRAINED MODEL RESULTS:\n",lfConstrained);

lnLikDiff = 2(fullModelLik-res1[1][0]);

degFDiff = fullModelParam-res1[1][1];

fprintf (stdout, "\n", separator,"\n\n-2(Ln Likelihood Ratio)=",lnLikDiff,"\n","Constrained parameters:",Format(degFDiff,0,0));

fprintf (stdout, "\nP-Value:",1-CChi2(lnLikDiff,degFDiff));

