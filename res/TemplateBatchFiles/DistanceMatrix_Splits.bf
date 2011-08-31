ChoiceList (distanceChoice, "Distance Computation",1,SKIP_NONE,
			"Distance formulae","Use one of the predefined distance measures based on data comparisons. Fast.",
			"Full likelihood","Estimate distances using pairwise MLE. More choices but slow.");
			
if (distanceChoice < 0)
{
	return 0;
}


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


DataSet ds = ReadDataFile (PROMPT_FOR_FILE);

fprintf (stdout, "Read:\n", ds, "\n");

window_size = 0;
while (window_size < 1 || window_size > ds.sites)
{
	fprintf (stdout, "Window Size: ");
	fscanf  (stdin,  "Number", window_size);
}

for (windowPosition = 0; windowPosition < ds.sites; windowPosition=windowPosition+window_size)
{
	if (dataType)
	{
		DataSetFilter filteredData = CreateFilter (ds,3,siteIndex>=windowPosition && siteIndex<windowPosition+window_size,"",GeneticCodeExclusions);
	}
	else
	{
		DataSetFilter filteredData = CreateFilter (ds,1,siteIndex>=windowPosition && siteIndex<windowPosition+window_size,"");
	}

	if (distanceChoice > 0)
	{
		SelectTemplateModel(filteredData);
	}

	if (windowPosition==0)
	{
		ChoiceList (distanceFormat, "Matrix Format",1,SKIP_NONE,
					"PHYLIP Style","Symmetric PHYLIP style matrix.",
					"HyPhy Matrix","Use HyPhy matrix syntax (suitable for input to other HyPhy analyses) and is written to disk incrementally");
					
		if (distanceFormat < 0)
		{
			return 0;
		}
	}

	distanceMatrix = {ds.species,ds.species};
	if (windowPosition == 0)
	{
		if (distanceFormat == 0)
		{
			SetDialogPrompt ("Save distance matrix (symmetric PHYLIP style):");
			fprintf (PROMPT_FOR_FILE,CLEAR_FILE);
		}
		else
		{
			SetDialogPrompt ("Save distance matrix (HyPhy style):");
			fprintf (PROMPT_FOR_FILE,CLEAR_FILE);
			BASE_PATH = LAST_FILE_PATH;
		}
	}
	LAST_FILE_PATH = BASE_PATH + "." + windowPosition;
	if (distanceFormat == 0)
	{
		SetDialogPrompt ("Save distance matrix (symmetric PHYLIP style):");
		fprintf (LAST_FILE_PATH,CLEAR_FILE,ds.species,"\n");
	}
	else
	{
		SetDialogPrompt ("Save distance matrix (HyPhy style):");
		fprintf (LAST_FILE_PATH,CLEAR_FILE);
	}
		

	if (distanceChoice)
	{
		fprintf (stdout,"\nHYPHY is computing pairwise maximum likelihood distance estimates. A total of ", Format(ds.species*(ds.species-1)/2,0,0),
					    " estimations will be performed.\n");

		MESSAGE_LOGGING = 0;

		for (i = 0; i<ds.species-1; i=i+1)
		{
			for (j = 0; j<=i; j = j+1)
			{
				if (dataType)
				{
					DataSetFilter twoSpecFilter = CreateFilter (filteredData,3,"",(speciesIndex==i+1)||(speciesIndex==j),GeneticCodeExclusions);
				}
				else
				{
					DataSetFilter twoSpecFilter = CreateFilter (filteredData,1,"",(speciesIndex==i+1)||(speciesIndex==j));
				}
				if (FREQUENCY_SENSITIVE)
				{
					if (USE_POSITION_SPECIFIC_FREQS)
					{
						HarvestFrequencies (vectorOfFrequencies,filteredData,3,1,1);
					}
					else
					{
						HarvestFrequencies (vectorOfFrequencies,twoSpecFilter,1,1,0);
					}
				}
				if (FREQUENCY_SENSITIVE)
				{
					MULTIPLY_BY_FREQS = PopulateModelMatrix ("modelMatrix",vectorOfFrequencies);
					if (dataType)
					{
						codonFrequencies = BuildCodonFrequencies (vectorOfFrequencies);
						Model pairModel = (modelMatrix, codonFrequencies, MULTIPLY_BY_FREQS);
					}
					else
					{
						Model pairModel = (modelMatrix, vectorOfFrequencies, MULTIPLY_BY_FREQS);
					}
				}		
				else
				{
					if (i+j==0)
					{
						MULTIPLY_BY_FREQS = PopulateModelMatrix ("modelMatrix",equalFreqs);
						Model pairModel = (modelMatrix, equalFreqs, MULTIPLY_BY_FREQS);
					}
				}
				
				Tree Inferred_Tree = (1,2);
				LikelihoodFunction lf = (twoSpecFilter,Inferred_Tree);
				Optimize (res,lf);
				k = BranchLength (Inferred_Tree,0);
				distanceMatrix[j][i+1] = k;
				distanceMatrix[i+1][j] = k;
			}
		}
	}
	else
	{
		tdc = 0;
		tdp = 0;
		ldp = 0;
		togo = ds.species*(ds.species-1)/2;


		if (windowPosition == 0)
		{
			DISTANCE_PROMPTS = 1;
			#include "chooseDistanceFormula.def";
			dummy = InitializeDistances (0);
		}
		
		fprintf (stdout,"\nHYPHY is computing pairwise distance estimates. A total of ", Format(togo,0,0),
					    " estimations will be performed.\n");
		
		if (distanceFormat == 0)
		{
			for (i = 0; i<ds.species-1; i=i+1)
			{
				for (j = 0; j<=i; j = j+1)
				{
					k = ComputeDistanceFormula (i+1,j);
					distanceMatrix[j][i+1] = k;
					distanceMatrix[i+1][j] = k;
				}
				tdc = tdc+i+1;
				tdp = (tdc/togo * 100)$1;
				if (tdp>ldp)
				{
					ldp = tdp;
					fprintf (stdout, ldp, "% done\n");
				}
			}
		}
		else
		{
			fprintf (LAST_FILE_PATH, CLEAR_FILE, "{\n");
			for (i = 0; i<ds.species; i=i+1)
			{
				outRow = "{0";
				for (j = 0; j<i; j = j+1)
				{
					outRow = outRow+",0";
				}
				for (j = i+1; j<ds.species; j = j+1)
				{
					k = ComputeDistanceFormula (i,j);
					outRow = outRow+","+k;
				}
				fprintf (LAST_FILE_PATH, outRow, "}\n");
				tdc = tdc+(ds.species-i-1);
				tdp = (tdc/togo * 100)$1;
				if (tdp>ldp)
				{
					ldp = tdp;
					fprintf (stdout, ldp, "% done\n");
				}
			}	
			fprintf (LAST_FILE_PATH, "}\n");
		}

		DISTANCE_PROMPTS = 0;
	}
		

	if (distanceFormat == 0)
	{
		spacer = "           ";

		for (i = 0; i<ds.species; i=i+1)
		{
			GetString	(seqString,ds,i);
			j = Abs(seqString);
			if (j>=10)
			{
				seqString = seqString[0][9];
			}
			else
			{
				seqString = seqString+spacer[0][10-j-1];
			}
			for (j = 0; j<ds.species; j = j+1)
			{
				seqString = seqString + "\t" + (distanceMatrix[i][j]);
			}
			fprintf (LAST_FILE_PATH, seqString, "\n");
		}
	}
	else
	{
		if (distanceChoice)
		{
			fprintf (LAST_FILE_PATH,CLEAR_FILE,distanceMatrix);
		}
	}
}
		
