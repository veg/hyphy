VERBOSITY_LEVEL = -1;
ALLOW_SEQUENCE_MISMATCH = 1;

#include "relrateBootstrap.bf";

fprintf(stdout,"\n ---- RUNNING RELATIVE RATE ANALYSIS ---- \n");

ChoiceList (dataType,"Data type",1,SKIP_NONE,"Nucleotide/Protein","Nucleotide or amino-acid (protein).",
				     "Codon","Codon (several available genetic codes).");

if (dataType<0) 
{
	return;
}
if (dataType)
{
	NICETY_LEVEL = 3;
	SetDialogPrompt ("Please choose a codon data file:");
	#include "TemplateModels/chooseGeneticCode.def";
}
else
{
	SetDialogPrompt ("Please choose a nucleotide or amino-acid data file:");
}

fileCount = 0;
speciesCount = 0;

if ((OPTIMIZATION_PRECISION==0)||(OPTIMIZATION_PRECISION>0.001))
{
	OPTIMIZATION_PRECISION = 0.001;
}

while (speciesCount < 3)
{
	if (fileCount>=1)
	{
		SetDialogPrompt ("Please choose an additional data file:");
		DataSet ds2 = ReadDataFile (PROMPT_FOR_FILE);
		speciesCount = speciesCount + ds2.species;
		DataSet ds = Combine (ds,ds2);
	}
	else
	{
		DataSet ds = ReadDataFile (PROMPT_FOR_FILE);
		speciesCount = speciesCount + ds.species;
	}
	
	fileCount = fileCount+1;
}

firstSpec = 0;
secondSpec = 1;
thirdSpec = 2;
fprintf (stdout,"The following data were read:\n",ds,"\n");

if (speciesCount>3) 
{
	speciesNames = {1,3};
	ChoiceList (twoSpec, "Choose 2 taxa to share rates:",2,SKIP_NONE,ds);
	if (twoSpec[0]<0)
	{
		return ;
	}
	firstSpec = twoSpec[0];
	secondSpec = twoSpec[1];
	ChoiceList (thirdSpec, "Choose the outgroup:",1,twoSpec,ds);
	
}
else
	{		
		ChoiceList (thirdSpec,"Choose the outgroup:",1,SKIP_NONE,ds);
		
		if (thirdSpec == 1)
		{
			firstSpec = 0;
			secondSpec = 2;
		}
		else
		{
			if (thirdSpec == 2)
			{
				firstSpec = 0;
				secondSpec = 1;
			}
			else
			{
				firstSpec = 1;
				secondSpec = 2;
				thirdSpec = 0;
			}

		}
			
	}

		
if (dataType)
{
	DataSetFilter filteredData = CreateFilter (ds,3,"",(speciesIndex==firstSpec)||(speciesIndex==secondSpec)||(speciesIndex==thirdSpec),GeneticCodeExclusions);
}
else
{
	DataSetFilter filteredData = CreateFilter (ds,1,"",(speciesIndex==firstSpec)||(speciesIndex==secondSpec)||(speciesIndex==thirdSpec));
}

SelectTemplateModel(filteredData);

if ((firstSpec<thirdSpec)&&(secondSpec<thirdSpec))
{
	treeString = "(Ingroup1,Ingroup2,OutGroup)";
}
else
{
	if ((firstSpec>thirdSpec)&&(secondSpec>thirdSpec))
	{
		treeString = "(OutGroup,Ingroup1,Ingroup2)";
	}
	else
	{
		treeString = "(Ingroup1,OutGroup,Ingroup2)";
	}
}

relationString = ":=";

#include "selectModelParameters.bf";

Tree threeTaxaTree = treeString;

LikelihoodFunction lf = (filteredData,threeTaxaTree);

GetString (specName,ds,thirdSpec);
fprintf (stdout,"\nOutgroup = ",specName);
GetString (specName,ds,firstSpec);
fprintf (stdout,"\nIngroup1 = ",specName);
GetString (specName,ds,secondSpec);
fprintf (stdout,"\nIngroup2 = ",specName);


Optimize (res,lf);

fprintf (stdout, "\n_________________RESULTS_________________\n\nFull Model Results:\n",lf);

Tree constrained3TaxaTree = treeString;

/* now specify the constraint */

LikelihoodFunction lfConstrained = (filteredData,constrained3TaxaTree);

ReplicateConstraint (constraintString,constrained3TaxaTree.Ingroup1, constrained3TaxaTree.Ingroup2);

Optimize (res1,lfConstrained);

fprintf (stdout, "\n\nConstrained Model Results:\n",lfConstrained);

lnLikDiff = -2(res1[1][0]-res[1][0]);

degFDiff = res[1][1]-res1[1][1];

fprintf (stdout, "\n\n-2(Ln Likelihood Difference)=",lnLikDiff,"\n","Difference in number of parameters:",Format(degFDiff,0,0));

fprintf (stdout, "\nP-Value:",1-CChi2(lnLikDiff,degFDiff));
