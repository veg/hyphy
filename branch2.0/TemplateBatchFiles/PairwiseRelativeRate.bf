ALLOW_SEQUENCE_MISMATCH = 1;
VERBOSITY_LEVEL = -1;

fprintf(stdout,"\n ---- RUNNING PAIRWISE RELATIVE RATE ANALYSIS ---- \n");

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

speciesCount = 0;

if ((OPTIMIZATION_PRECISION==0)||(OPTIMIZATION_PRECISION>0.001))
{
	OPTIMIZATION_PRECISION = 0.001;
}

while (speciesCount < 3)
{
	SetDialogPrompt ("Choose a data file with 3 or more sequences");
	DataSet ds = ReadDataFile (PROMPT_FOR_FILE);
	speciesCount = ds.species;
}

ChoiceList (thirdSpec,"Choose the outgroup:",1,SKIP_NONE,ds);

if (thirdSpec<0)
{
	return;
}

GetString (outgroupName,ds,thirdSpec);

if (dataType)
{
	DataSetFilter filteredData = CreateFilter (ds,3);
}
else
{
	DataSetFilter filteredData = CreateFilter (ds,1);
}

SelectTemplateModel(filteredData);

relationString = ":=";

#include "selectModelParameters.bf";

SetDialogPrompt ("Save full results to:");

modelParameterCount = Rows("LAST_MODEL_PARAMETER_LIST");

fprintf (PROMPT_FOR_FILE,CLEAR_FILE);

MESSAGE_LOGGING = 0;

maxLength = 0;


for (firstSpec = 0; firstSpec < speciesCount; firstSpec = firstSpec+1)
{
	GetString (firstName,ds,firstSpec);

	tLength = Abs(outgroupName)+Abs(firstName);
	
	for (secondSpec = firstSpec+1; secondSpec < speciesCount; secondSpec = secondSpec+1)
	{
		GetString (secondName,ds,secondSpec);
		tLength2 = tLength + Abs (secondName) + 30;
		if (tLength2>maxLength)
		{
			maxLength = tLength2;
		}
	}
}

maxLength = maxLength+6;
fprintf(stdout,"\n\n In the summary table below, branch lengths denote");

if (parameter2Constrain)
{
	fprintf (stdout," the value of parameter ", funnyString, " for that branch.");
}
else
{
	fprintf (stdout," the expected number of substituions per that branch.");
}

fprintf (stdout,"\n\n (*)   corresponds to the .05 significance level\n (**)  corresponds to the .01 significance level\n (***) corresponds to the .001 significance level.\n\n");

fprintf (stdout,"\n\n  Taxa Triplet");

for (tLength = 15; tLength<maxLength; tLength=tLength+1)
{
	fprintf (stdout," ");
}

fprintf (stdout,"      LRT        P-Value \n");

separatorString = "-";

for (tLength = 0; tLength<=maxLength; tLength=tLength+1)
{
	separatorString = separatorString + "-";
}

separatorString = separatorString + "+";

for (tLength = 0; tLength<12; tLength=tLength+1)
{
	separatorString = separatorString + "-";
}

separatorString = separatorString + "+";

for (tLength = 0; tLength<12; tLength=tLength+1)
{
	separatorString = separatorString + "-";
}



tLength2 = 0;

for (firstSpec = 0; firstSpec < speciesCount; firstSpec = firstSpec+1)
{

	if (firstSpec == thirdSpec)
	{
		continue;
	}
	
	fprintf (stdout,separatorString);

	GetString (firstName,ds,firstSpec);
	
	for (secondSpec = firstSpec+1; secondSpec < speciesCount; secondSpec = secondSpec+1)
	{
		if (secondSpec == thirdSpec)
		{
			continue;
		}
		
		if (dataType)
		{
			DataSetFilter filteredData = CreateFilter (ds,3,"",(speciesIndex==firstSpec)||(speciesIndex==secondSpec)||(speciesIndex==thirdSpec),GeneticCodeExclusions);
		}
		else
		{
			DataSetFilter filteredData = CreateFilter (ds,1,"",(speciesIndex==firstSpec)||(speciesIndex==secondSpec)||(speciesIndex==thirdSpec));
		}

		HarvestFrequencies (rrEFV,filteredData,1,1,1);
		
		if (FREQUENCY_SENSITIVE)
		{
			rrModelMatrix = 0;
			if (USE_POSITION_SPECIFIC_FREQS)
			{
				HarvestFrequencies (vectorOfFrequencies,filteredData,3,1,1);
			}
			MULTIPLY_BY_FREQS = PopulateModelMatrix ("rrModelMatrix",rrEFV);
			if (dataType)
			{
				rrCodonEFV = BuildCodonFrequencies (rrEFV);
				Model rrModel = (rrModelMatrix,rrCodonEFV,MULTIPLY_BY_FREQS);	
			}
			else
			{
				Model rrModel = (rrModelMatrix,rrEFV,MULTIPLY_BY_FREQS);
			}	
		}
		
		if ((firstSpec<thirdSpec)&&(secondSpec<thirdSpec))
		{
			treeString = "(FirstSpecies,SecondSpecies,OutGroup)";
			i1 = 0;
			i2 = 1;
			i3 = 2;
		}
		else
		{
			if ((firstSpec>thirdSpec)&&(secondSpec>thirdSpec))
			{
				treeString = "(OutGroup,FirstSpecies,SecondSpecies)";
				i1 = 1;
				i2 = 2;
				i3 = 0;
			}
			else
			{
				treeString = "(FirstSpecies,OutGroup,SecondSpecies)";
				i1 = 0;
				i2 = 2;
				i3 = 1;
			}
		}

		Tree threeTaxaTree = treeString;

		LikelihoodFunction lf = (filteredData,threeTaxaTree);

		Optimize (res,lf);

		Tree constrained3TaxaTree = treeString;

		/* now specify the constraint */

		LikelihoodFunction lfConstrained = (filteredData,constrained3TaxaTree);

		ReplicateConstraint (constraintString,constrained3TaxaTree.FirstSpecies, constrained3TaxaTree.SecondSpecies);

		Optimize (res1,lfConstrained);

		lnLikDiff = -2(res1[1][0]-res[1][0]);

		degFDiff = res[1][1]-res1[1][1];
		
		if (lnLikDiff>0)
		{
			pValue = 1.0-CChi2(lnLikDiff,degFDiff);
		}
		else
		{
			pValue = 1;
			fprintf (MESSAGE_LOG,"\nA negative LRT statistic encoutered. You may want to increase the optimization precision settings to resolve numerical apporximation errors");
		}

		GetString (secondName,ds,secondSpec);
		
		if (parameter2Constrain==0)
		{
			b_out = BranchLength(threeTaxaTree,i3);
			b_1 = BranchLength(threeTaxaTree,i1);
			b_2 = BranchLength(threeTaxaTree,i2);
		}
		else
		{
			b_out = res[0][modelParameterCount*i3+parameter2Constrain-1];
			b_1 = res[0][modelParameterCount*i1+parameter2Constrain-1];
			b_2 = res[0][modelParameterCount*i2+parameter2Constrain-1];
		}
		
		printTreeString = "("+outgroupName+":"+b_out+",("
		+firstName+":"+b_1+","+secondName+":"+b_2+"))";

		fprintf (stdout, "\n", printTreeString);
		
		for (tLength=Abs(printTreeString);tLength<maxLength;tLength=tLength+1)
		{
			fprintf (stdout," ");
		}
		
		fprintf (stdout, "  | ", Format(lnLikDiff,9,4),"  | ", pValue );
		
		if (pValue<0.05)
		{
			if (pValue<0.01)
			{
				if (pValue<0.001)
				{
				    fprintf (stdout," (***) ");
				}
				else
				{
			     	fprintf (stdout," (**) ");
			    }
			}
			else
			{
			     fprintf (stdout," (*) ");			
			}
		}
		if (!tLength2)
		{
			tLength2 = 1;
			dataDimension = Columns(res);
			fprintf (LAST_FILE_PATH,",,,Unconstrained Optimizations");
			for (i=0;i<dataDimension;i=i+1)
			{
				fprintf	  (LAST_FILE_PATH,", ");
			}
			fprintf (LAST_FILE_PATH,"Constrained Optimizations");
			for (i=0;i<dataDimension;i=i+1)
			{
				fprintf	  (LAST_FILE_PATH,", ");
			}
			fprintf (LAST_FILE_PATH,"\nOutgroup,Taxon 1,Taxon 2,Ln-likelihood");
			for (i=0;i<dataDimension;i=i+1)
			{
				GetString (argName,lf,i);
				fprintf	  (LAST_FILE_PATH,",",argName);
			}
			fprintf (LAST_FILE_PATH,",Ln-likelihood");
			for (i=0;i<dataDimension;i=i+1)
			{
				GetString (argName,lfConstrained,i);
				fprintf	  (LAST_FILE_PATH,",",argName);
			}	
			fprintf (LAST_FILE_PATH,",LR,P-Value");
		}	
		fprintf (LAST_FILE_PATH,"\n",outgroupName,",",firstName,",",secondName,",",res[1][0]);
		for (i=0;i<dataDimension;i=i+1)
		{
			fprintf (LAST_FILE_PATH,",",res[0][i]);
		}
		fprintf (LAST_FILE_PATH,",",res1[1][0]);
		for (i=0;i<dataDimension;i=i+1)
		{
			fprintf (LAST_FILE_PATH,",",res1[0][i]);
		}
		fprintf (LAST_FILE_PATH,",",lnLikDiff ,",",pValue);
	}
	fprintf (stdout,"\n");
}

MESSAGE_LOGGING = 1;
