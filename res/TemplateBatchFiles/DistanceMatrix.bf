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


DataSet ds   = ReadDataFile (PROMPT_FOR_FILE);

mpiPrefix    = ""; mpiPrefix * 128;
mpiPrefix	 * ("DataSet ds   = ReadDataFile (\"" + LAST_FILE_PATH + "\");\n");
 

if (dataType)
{
	DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);
	mpiPrefix	 * ("DataSetFilter filteredData = CreateFilter (ds,3,\"\",\"\",\""+GeneticCodeExclusions+"\");");
}
else
{
	DataSetFilter filteredData = CreateFilter (ds,1,"","");
	mpiPrefix	 * ("DataSetFilter filteredData = CreateFilter (ds,1);\n");
}

if (distanceChoice > 0)
{
	SelectTemplateModel(filteredData);
}

ChoiceList (distanceFormat, "Matrix Format",1,SKIP_NONE,
			"PHYLIP Style","Symmetric PHYLIP style matrix.",
			"HyPhy Matrix","Use HyPhy matrix syntax (suitable for input to other HyPhy analyses) and is written to disk incrementally",
			"TAB","Use tab separated format, suitable for importing into a statistical analysis package and written to disk incrementally",
			"Pairwise","Write a tab format in the form: seq1 seq2 distance");
			
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
	mpiPrefix * 0;
	fprintf 	 (stdout,"\nHYPHY is computing pairwise maximum likelihood distance estimates. A total of ", Format(ds.species*(ds.species-1)/2,0,0),
				    	 " estimations will be performed.\n");
	_pddVF = 1;
	
	pairwiseAlign = 0;
	
	GetDataInfo (charInfo, filteredData, "CHARACTERS");
	ChoiceList (pairwiseAlign, "Peform pairwise alignment prior to comparison",1,SKIP_NONE,
				"No","Assume that the sequences are already aligned.",
				"Yes","Align each pair of sequences before computing pairwise distances");

	if (pairwiseAlign < 0)
	{
		return 0;
	}
	
			
	
	if (pairwiseAlign)
	{
		if (Columns (charInfo) == 20)
		{
			_skipPredefsSeqAlignShared = 1;
			ExecuteAFile ("SeqAlignShared.ibf");
		}
		else
		{
            ChoiceList (codonAlign, "Use codon-based alignment",1,SKIP_NONE,
                                        "No","Align nucleotide sequences directly.",
                                        "Yes","Align nucleotide sequences in the codon-space fixing frameshifts if necessary");
                                        
            if (codonAlign < 0)
            {
                return 0;
            }
            if (codonAlign == 0)
            {
                alignOptions = {};
                scoreMatrix = {
                {5,-4,-4,-4}
                {-4,5,-4,-4}
                {-4,-4,5,-4}
                {-4,-4,-4,5}
                };
                alignOptions ["SEQ_ALIGN_SCORE_MATRIX"] = 	scoreMatrix;
                alignOptions ["SEQ_ALIGN_GAP_OPEN"]		= 	10;
                alignOptions ["SEQ_ALIGN_GAP_OPEN2"]	= 	5;
                alignOptions ["SEQ_ALIGN_GAP_EXTEND"]	= 	1;
                alignOptions ["SEQ_ALIGN_GAP_EXTEND2"]	= 	1;
                alignOptions ["SEQ_ALIGN_AFFINE"]		=   1;
                alignOptions ["SEQ_ALIGN_CHARACTER_MAP"]=   "ACGT";
            }
            else
            {
                _skipPredefsSeqAlignShared = 1;
				LoadFunctionLibrary ("SeqAlignmentCodonShared");
			}
		}
	}

    LoadFunctionLibrary ("pairwiseDistanceEstimator");	
}
else
{
	DISTANCE_PROMPTS = 1;
	
	#include "chooseDistanceFormula.def";
	
	InitializeDistances (0);
	
	mpiPrefix * ("ExecuteAFile(\""+_distanceFilePath+"\"); InitializeDistances (0);\n");
	mpiPrefix * 0;
				    
	tdc = 0;
	tdp = 0;
	ldp = 0;
	togo = ds.species*(ds.species-1)/2;

	fprintf (stdout,"\nHYPHY is computing pairwise distance estimates. A total of ", Format(togo,0,0),
				    " estimations will be performed.\n");
	
	distanceMatrix = {ds.species,ds.species};
	
	if (MPI_NODE_COUNT > 1)
	{
		MPI_NODE_INFO  = {MPI_NODE_COUNT-1,2};
		perNodeFile    = togo / (MPI_NODE_COUNT-1);
		
		
		accumulator = 0;
		lastRowSent = 0;
		totalSent	= 0;
		
		for (i = 0; i < ds.species-1; i = i + 1)
		{
			accumulator = accumulator + (ds.species-1-i);
			
			if (accumulator >= perNodeFile || i == ds.species-2)
			{
				fprintf (stdout, "[SENT MATRIX ROWS ", lastRowSent, "-", i, " to MPI node ", totalSent + 1,"]\n");
				MPI_NODE_INFO[totalSent][0] = lastRowSent;
				MPI_NODE_INFO[totalSent][1] = i;
				toSend	= mpiPrefix + "vecSize = " + accumulator + ";" +
												  "fromRow = " + lastRowSent + ";" + 
												  "toRow   = " + i + ";" +
												  "ExecuteAFile (HYPHY_LIB_DIRECTORY + \"TemplateBatchFiles\" + DIRECTORY_SEPARATOR + \"pairwiseDistanceEstimatorCounter.ibf\"); return result;";
												  
				MPISend (totalSent+1, toSend);
				
				accumulator = 0;
				lastRowSent	= i+1;
				totalSent   = totalSent + 1;
			}
		}
		
		for (i = 0; i < totalSent; i = i+1)
		{
			MPIReceive		(-1,fromNode,theVector);
			fromNode		= fromNode - 1;
			ExecuteCommands ("partialVector = " + theVector);
			fromRow			= MPI_NODE_INFO[fromNode][0];
			toRow			= MPI_NODE_INFO[fromNode][1];
			fprintf (stdout, "[GOT MATRIX ROWS ", fromRow, "-", toRow, " from MPI node ", fromNode,"]\n");
			j = 0;
			for (r = fromRow; r <= toRow; r = r + 1)
			{
				for (c = r+1; c < ds.species; c = c + 1)
				{
					distanceMatrix [r][c] = partialVector[j];
					distanceMatrix [c][r] = partialVector[j];
					j = j+1;
				}
			}
		}
		
	}
	else
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

	DISTANCE_PROMPTS = 0;
}
	

if (distanceFormat != 1 && distanceFormat != 3)
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
	if (distanceFormat == 1)
	{
		fprintf (LAST_FILE_PATH,CLEAR_FILE,distanceMatrix);
	}
	else
	{
		GetString (seqString,ds,-1);
		fprintf (LAST_FILE_PATH,CLEAR_FILE,"Sequence1\tSequence2\tDistance");
		for (i = 0; i<ds.species; i=i+1)
		{
			for (j = i+1; j<ds.species; j = j+1)
			{
				fprintf (LAST_FILE_PATH,"\n",seqString[i],"\t",seqString[j],"\t",distanceMatrix[i][j]);
			}
		}		
	}
}
	
fprintf (LAST_FILE_PATH, CLOSE_FILE);
