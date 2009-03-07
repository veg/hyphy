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

if (dataType)
{
	DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);
}
else
{
	DataSetFilter filteredData = CreateFilter (ds,1,"","");
}

if (distanceChoice > 0)
{
	SelectTemplateModel(filteredData);
}

ChoiceList (distanceFormat, "Matrix Format",1,SKIP_NONE,
			"PHYLIP Style","Symmetric PHYLIP style matrix.",
			"HyPhy Matrix","Use HyPhy matrix syntax (suitable for input to other HyPhy analyses) and is written to disk incrementally",
			"TAB","Use tab separated format, suitable for importing into a statistical analysis package and written to disk incrementally");
			
if (distanceFormat < 0)
{
	return 0;
}

if (distanceFormat == 0)
{
	SetDialogPrompt ("Save distance matrix (symmetric PHYLIP style):");
	fprintf (PROMPT_FOR_FILE,CLEAR_FILE,KEEP_OPEN,ds.species);
	distanceMatrix = {ds.species,ds.species};
}
else
{
	if (distanceFormat == 1)
	{
		SetDialogPrompt ("Save distance matrix (HyPhy style):");
		if (distanceChoice==1)
		{
			distanceMatrix = {ds.species,ds.species};
		}
		fprintf (PROMPT_FOR_FILE,CLEAR_FILE,KEEP_OPEN);
	}
	else
	{
		SetDialogPrompt ("Save distance matrix (tab format):");
		fprintf (PROMPT_FOR_FILE,CLEAR_FILE,KEEP_OPEN);
		distanceMatrix = {ds.species,ds.species};
	}
}
	

if (distanceChoice)
{
	fprintf 	 (stdout,"\nHYPHY is computing pairwise maximum likelihood distance estimates. A total of ", Format(ds.species*(ds.species-1)/2,0,0),
				    	 " estimations will be performed.\n");
	_pddVF = 1;
	ExecuteAFile ("pairwiseDistanceEstimator.ibf");
}
else
{
	DISTANCE_PROMPTS = 1;
	
	#include "chooseDistanceFormula.def";
	
	dummy = InitializeDistances (0);
	
				    
	tdc = 0;
	tdp = 0;
	ldp = 0;
	togo = ds.species*(ds.species-1)/2;

	fprintf (stdout,"\nHYPHY is computing pairwise distance estimates. A total of ", Format(togo,0,0),
				    " estimations will be performed.\n");
	
	if (distanceFormat != 1)
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
			outRow = "";
			outRow * 256;
			outRow * "{0";
			for (j = 0; j<i; j = j+1)
			{
				outRow * ",0";
			}
			for (j = i+1; j<ds.species; j = j+1)
			{
				k = ComputeDistanceFormula (i,j);
				outRow * (","+k);
			}
			outRow * 0;
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
	

if (distanceFormat != 1)
{
	if (distanceFormat == 2)
	{
		GetString	(seqString,ds,0);
		fprintf (LAST_FILE_PATH, CLEAR_FILE, seqString);
		for (i = 1; i<ds.species; i=i+1)
		{
			GetString	(seqString,ds,i);
			fprintf (LAST_FILE_PATH, "\t", seqString);
		}
	}
	else
	{
		spacer = "           ";
	}
	
	for (i = 0; i<ds.species; i=i+1)
	{
		GetString	(seqString,ds,i);
		if (distanceFormat == 0)
		{
			j = Abs(seqString);
			if (j>=10)
			{
				seqString = seqString[0][9];
			}
			else
			{
				seqString = seqString+spacer[0][10-j-1];
			}
		}
		fprintf (LAST_FILE_PATH, "\n", seqString);
		
		for (j = 0; j<ds.species; j = j+1)
		{
			fprintf (LAST_FILE_PATH,"\t",distanceMatrix[i][j]);
		}
	}
}
else
{
	if (distanceChoice)
	{
		fprintf (LAST_FILE_PATH,CLEAR_FILE,distanceMatrix);
	}
}
	
fprintf (LAST_FILE_PATH, CLOSE_FILE);