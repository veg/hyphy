fprintf(stdout,"\n ---- RUNNING PAIRWISE RELATIVE RATIO ANALYSIS ---- \n");

ChoiceList (dataType,"Data type",1,SKIP_NONE,"Nucleotide/Protein","Nucleotide or amino-acid (protein).",
				     "Codon","Codon (several available genetic codes).");

if (dataType<0) 
{
	return;
}

if (dataType)
{
	NICETY_LEVEL = 3;
	ChoiceList (codeTables,"Genetic Code",1,SKIP_NONE,"Same Code","All files use the same genetic code. The index file contains file names only. ",
				     	"Different Codes","Each file name line in the index file is followed by a line with a HYPHY genetic table index. See Help for details.");
	if (codeTables<0)
	{
		return;
	}
	if (codeTables>0)
	{
		skipCodeSelectionStep = 1;
	}
	#include "TemplateModels/chooseGeneticCode.def";
}

#include "readIndexFile.bf";

fileCount = Rows(stringMatrix);

if (fileCount<=1)
{
	fprintf (stdout,"\n\n\n****** EMPTY INDEX FILE - NOTHING TO DO! **********\n\n\n");
	return;
}

/* read the first file */

DataSet ds = ReadDataFile (stringMatrix[0]);

referenceSpecCount = ds.species;

treeRead = 0;

if (IS_TREE_PRESENT_IN_DATA)
{
	treeString = DATAFILE_TREE;
	treeRead = 1;
}

fprintf(stdout,"\n File 1:",stringMatrix[0]," - ", Format (ds.species,0,0)," sequences, ", Format (ds.sites,0,0)," sites.");

for (counter = 1; counter<fileCount; counter = counter+1)
{
	DataSet ds2 = ReadDataFile (stringMatrix[counter]);
	if (ds2.species!=referenceSpecCount)
	{
		fprintf (stdout,"\n\n\n****** SEQUENCE COUNT MISMATCH **********\n");
		fprintf (stdout,"\n\n File ",stringMatrix[counter]," had ",Format(ds2.species,0,0)," sequences, whereas ",Format(referenceSpecCount,0,0)," were expected.\n\n\n");
		return;
	}
	fprintf(stdout,"\n File ",Format(counter+1,0,0),":",stringMatrix[counter]," - ", Format (ds2.species,0,0)," sequences, ", Format (ds2.sites,0,0)," sites.");
	if (!treeRead)
	{
		if (IS_TREE_PRESENT_IN_DATA)
		{
			treeString = DATAFILE_TREE;
			treeRead = 1;
		}
	}
}

fprintf (stdout,"\n\n** All data files were read successfully **\n\n");

if (dataType)
{
	if (codeTables)
	{
		dummy = ApplyGeneticCodeTable (codeTableMatrix[0]);
		ModelMatrixDimension = 0;
	}
	DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);
}
else
{
	DataSetFilter filteredData = CreateFilter (ds,1);
}


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

SelectTemplateModel(filteredData);

global  RelRatio;

relationString = ":=RelRatio*";

#include "selectModelParameters.bf";

SetDialogPrompt ("Save full results to:");

modelParameterCount = Rows("LAST_MODEL_PARAMETER_LIST");

fprintf (PROMPT_FOR_FILE,CLEAR_FILE);

tabulatedFileName = LAST_FILE_PATH;


singleFileResults = {fileCount,1};

fprintf (stdout,"\n\n***** RUNNING SINGLE FILE ANALYSES *****\n\n");

fullParameterCount = 0;

/*MESSAGE_LOGGING = 0;*/

OPTIMIZATION_PRECISION = OPTIMIZATION_PRECISION/10;

timer = Time(0);

for (counter = 1; counter<= fileCount; counter = counter+1)
{
	HarvestFrequencies (vectorOfFrequencies,filteredData,1,1,1);
	/*fprintf (stdout,"\n",GeneticCodeExclusions,"\n",_Genetic_Code,"\n");*/
	if (FREQUENCY_SENSITIVE)
	{
		modelMatrix = 0;
		if (USE_POSITION_SPECIFIC_FREQS)
		{
			HarvestFrequencies (vectorOfFrequencies,filteredData,3,1,1);
		}
		MULTIPLY_BY_FREQS = PopulateModelMatrix ("modelMatrix",vectorOfFrequencies);
		if (dataType)
		{
			CodonEFV = BuildCodonFrequencies (vectorOfFrequencies);
			Model RRmodel = (modelMatrix,CodonEFV,MULTIPLY_BY_FREQS);
		}
		else
		{
			Model RRmodel = (modelMatrix,vectorOfFrequencies,MULTIPLY_BY_FREQS);
		}	
	}	
	Tree firstFileTree = treeString;
	
	LikelihoodFunction lf = (filteredData,firstFileTree);
	Optimize (res,lf);

	if (counter<fileCount)
	{
		DataSet ds = ReadDataFile (stringMatrix[counter]);
		if (dataType)
		{
			if (codeTables)
			{
				dummy = ApplyGeneticCodeTable (codeTableMatrix[counter]);
				ModelMatrixDimension = 0;
			}
			DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);
		}
		else
		{
			DataSetFilter filteredData = CreateFilter (ds,1);
		}	
	}
	else
	{
		fullParameterCount = res[1][1];
	}
	singleFileResults [counter-1] = res[1][0];
	fprintf (stdout,"\nFile ",stringMatrix[counter-1]," : ln-likelihood = ",res[1][0]);
	if (counter==1)
	{
		if (codeTables)
		{
			fprintf (tabulatedFileName,"\n\nSingle File Results\n\nFile Number,File Path,Code Table,Ln-Likelihood");
		}		
		else
		{
			fprintf (tabulatedFileName,"\n\nSingle File Results\n\nFile Number,File Path,Ln-Likelihood");
		}
		dataDimension = Columns(res);
		for (counter3 = 0; counter3 < dataDimension; counter3=counter3+1)
		{
			GetString (thisLine,lf,counter3);
			fprintf (tabulatedFileName,",",thisLine);
		}
	}
	
	if (codeTables)
	{
		fprintf (tabulatedFileName,"\n",Format(counter,0,0),",",stringMatrix[counter-1],",",Format(codeTableMatrix[counter-1],0,0),",",res[1][0]);	
	}
	else
	{
		fprintf (tabulatedFileName,"\n",Format(counter,0,0),",",stringMatrix[counter-1],",",res[1][0]);	
	}
	for (counter3 = 0; counter3 < dataDimension; counter3=counter3+1)
	{
		fprintf (tabulatedFileName,",",res[0][counter3]);
	}
}

/*OPTIMIZATION_PRECISION = OPTIMIZATION_PRECISION*10;*/

timer = (Time(0)-timer)/fileCount;

counter = fileCount*(fileCount+1)/2;

fprintf (stdout,"\n\n***** RUNNING PAIRWISE CONSTRAINED ANALYSES (",Format(counter,0,0), " PAIR) *****\n");

timer = counter*timer*2.5;

fprintf (stdout,"\n***** ESTIMATED TIME TO COMPLETE ALL ANALYSES: ",Format(timer$3600,0,0)," hours, ",Format((timer%3600)$60+1,0,0), " minutes. *****\n\n");

fprintf(stdout,"\n\n In the summary table below");

fprintf (stdout,"\n\n (*)   corresponds to the .05 significance level\n (**)  corresponds to the .01 significance level\n (***) corresponds to the .001 significance level.\n\n");

separator = "+--------+--------+--------------+--------------+";
fprintf (stdout,separator,"\n| File 1 | File 2 |      LRT     |    P-Value   |\n",separator);
for (counter = 0; counter< fileCount; counter = counter+1)
{
	DataSet ds = ReadDataFile (stringMatrix[counter]);
	if (dataType)
	{
		if (codeTables)
		{
			dummy = ApplyGeneticCodeTable (codeTableMatrix[counter]);
			ModelMatrixDimension = 0;
		}
		DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);
	}
	else
	{
		DataSetFilter filteredData = CreateFilter (ds,1);
	}		
	HarvestFrequencies (vectorOfFrequencies,filteredData,1,1,1);
	if (FREQUENCY_SENSITIVE)
	{
		modelMatrix = 0;
		if (USE_POSITION_SPECIFIC_FREQS)
		{
			HarvestFrequencies (vectorOfFrequencies,filteredData,3,1,1);
		}
		MULTIPLY_BY_FREQS = PopulateModelMatrix ("modelMatrix",vectorOfFrequencies);
		if (dataType)
		{
			CodonEFV = BuildCodonFrequencies (vectorOfFrequencies);
			Model RRmodel = (modelMatrix,CodonEFV,MULTIPLY_BY_FREQS);
		}
		else
		{
			Model RRmodel = (modelMatrix,vectorOfFrequencies,MULTIPLY_BY_FREQS);
		}	
	}	
	Tree firstFileTree = treeString;
	for (counter2 = counter+1; counter2< fileCount; counter2 = counter2+1)
	{
		fprintf (stdout,"\n| ",Format (counter+1,6,0)," | ", Format (counter2+1,6,0));
		DataSet ds2 = ReadDataFile (stringMatrix[counter2]);
		if (dataType)
		{
			if (codeTables)
			{
				dummy = ApplyGeneticCodeTable (codeTableMatrix[counter2]);
				ModelMatrixDimension = 0;
			}
			DataSetFilter filteredData2 = CreateFilter (ds2,3,"","",GeneticCodeExclusions);
		}
		else
		{
			DataSetFilter filteredData2 = CreateFilter (ds2,1);
		}		
		HarvestFrequencies (vectorOfFrequencies2,filteredData2,1,1,1);
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
				Model RRmodel2 = (modelMatrix2,CodonEFV2,MULTIPLY_BY_FREQS);
			}
			else
			{
				Model RRmodel2 = (modelMatrix2,vectorOfFrequencies2,MULTIPLY_BY_FREQS);
			}	
		}	
		
		Tree secondFileTree = treeString;	
			
		ReplicateConstraint (constraintString,secondFileTree,firstFileTree);

		LikelihoodFunction lfConstrained = (filteredData2,secondFileTree,filteredData,firstFileTree);
		
		Optimize (res1,lfConstrained);
			
		LRT = 2*(singleFileResults[counter]+singleFileResults[counter2]-res1[1][0]);

		degFDiff = 2*fullParameterCount-res1[1][1];
		
		if (LRT>0)
		{
			pValue = 1.0-CChi2(LRT,degFDiff);
		}
		else
		{
			pValue = 1;
			fprintf (MESSAGE_LOG,"\nA negative LRT statistic encoutered. You may want to increase the optimization precision settings to resolve numerical apporximation errors");
		}
		
		fprintf (stdout," | ",Format (LRT,12,5)," | ", Format (pValue,12,8)," |");	
		
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
		
		if ((counter==0)&&(counter2==1))
		{
			fprintf (tabulatedFileName,"\n\nRelative Ratio Tests\n\nFile 1,File 2,Ln-Likelihood,LRT,pValue");
			dataDimension = Columns(res1);
			for (counter3 = 0; counter3 < dataDimension; counter3=counter3+1)
			{
				GetString (thisLine,lfConstrained,counter3);
				fprintf (tabulatedFileName,",",thisLine);
			}
		}
		
		fprintf (tabulatedFileName,"\n",Format(counter+1,0,0),",",Format(counter2+1,0,0),",",res1[1][0],",",LRT,",",pValue);
		for (counter3 = 0; counter3 < dataDimension; counter3=counter3+1)
		{
			fprintf (tabulatedFileName,",",res1[0][counter3]);
		}
	}
}

fprintf (stdout,"\n",separator,"\n");

MESSAGE_LOGGING = 1;
