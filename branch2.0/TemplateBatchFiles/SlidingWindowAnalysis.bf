fprintf(stdout,"\n ---- RUNNING SLIDING DATA WINDOW ANALYSIS ---- \n");

ChoiceList (dataType,"Data type",1,SKIP_NONE,"Nucleotide/Protein","Nucleotide or amino-acid (protein).",
				     "Codon","Codon (several available genetic codes).");

if (dataType<0) 
{
	return;
}

if (dataType)
{
	SetDialogPrompt ("Please choose a codon data file:");
	#include "TemplateModels/chooseGeneticCode.def";
}
else
{
	SetDialogPrompt ("Please choose a nucleotide or amino-acid data file:");
}

DataSet ds = ReadDataFile (PROMPT_FOR_FILE);

fprintf (stdout,"\n\n**** Read file ", LAST_FILE_PATH," ****\n**** ", Format(ds.sites,0,0)," sites (", 
				Format(ds.unique_sites,0,0)," of them are distinct), and ", Format(ds.species,0,0), " species. ****"); 

if (IS_TREE_PRESENT_IN_DATA)
{
	treeString = DATAFILE_TREE;
	fprintf (stdout, "\n\nA tree was found in the data file:\n",treeString,"\n\nWould you like to use it:(Y/N)?");
	fscanf (stdin, "String", response);
	if ((response=="n")||(response=="N"))
	{
		IS_TREE_PRESENT_IN_DATA = 0;
	}
	fprintf (stdout, "\n\n");
}

if (!IS_TREE_PRESENT_IN_DATA)
{
	SetDialogPrompt ("Please select a tree file for the data:");
	fscanf (PROMPT_FOR_FILE, "String", treeString);
}

if (dataType)
{
	DataSetFilter filteredData = CreateFilter (ds,3);
}
else
{
	DataSetFilter filteredData = CreateFilter (ds,1);
}

SetDialogPrompt		("Save results to:");
fprintf 			(PROMPT_FOR_FILE,CLEAR_FILE);
SelectTemplateModel (filteredData);

windowWidth:>0;
windowShift:>0;

windowWidth = 0;
windowShift = 0;

if (dataType)
{
	while (windowWidth<3 || windowWidth> filteredData.sites*3-3 || windowWidth%3)
	{
		fprintf (stdout,"\nEnter the width (in bp, divisble by 3) of the data window [3-",filteredData.sites*3-3,"]:");
		fscanf  (stdin,"Number",windowWidth);
	}
	firstTime = filteredData.sites*3-windowWidth;
	while (windowShift<3 || windowShift> firstTime || windowShift%3)
	{
		fprintf (stdout,"\nEnter the stride (in bp, divisble by 3) at each step [3-",firstTime,"]:");
		fscanf  (stdin,"Number",windowShift);
	}
}
else
{
	while (windowWidth<1 || windowWidth> filteredData.sites-1)
	{
		fprintf (stdout,"\nEnter the width (in bp) of the data window [1-",filteredData.sites-1,"]:");
		fscanf  (stdin,"Number",windowWidth);
	}
	firstTime = filteredData.sites-windowWidth;
	while (windowShift<1 || windowShift> firstTime)
	{
		fprintf (stdout,"\nEnter the stride (in bp) at each step [1-",firstTime,"]:");
		fscanf  (stdin,"Number",windowShift);
	}
}

fprintf (stdout,"\n");

MESSAGE_LOGGING = 0;
VERBOSITY_LEVEL = -1;

firstTime = 1;

if (MPI_NODE_COUNT>1)
{
	MPINodeState  = {MPI_NODE_COUNT-1,1};
	MPINodeBounds = {MPI_NODE_COUNT-1,2};
	OPTIMIZE_SUMMATION_ORDER = 0;
}

for (currentStart = 0; currentStart < ds.sites; currentStart = currentStart + windowShift)
{
	currentEnd = currentStart + windowWidth;
	if (currentEnd>=ds.sites)
	{
		currentEnd = ds.sites - 1;
	}
	
	if (dataType)
	{
		DataSetFilter filteredData = CreateFilter (ds,3,(siteIndex>=currentStart)&&(siteIndex<=currentEnd),"",GeneticCodeExclusions);
	}
	else
	{
		DataSetFilter filteredData = CreateFilter (ds,1,(siteIndex>=currentStart)&&(siteIndex<=currentEnd));
	}

	HarvestFrequencies (swEFV,filteredData,1,1,1);
	
	if (FREQUENCY_SENSITIVE)
	{
		swModelMatrix = 0;
		if (USE_POSITION_SPECIFIC_FREQS)
		{
			HarvestFrequencies (swEFV,filteredData,3,1,1);
		}
		MULTIPLY_BY_FREQS = PopulateModelMatrix ("swModelMatrix",swEFV);
		if (dataType)
		{
			swCodonEFV    = BuildCodonFrequencies (swEFV);
			Model swModel = (swModelMatrix,swCodonEFV,MULTIPLY_BY_FREQS);	
		}
		else
		{
			Model swModel = (swModelMatrix,swEFV,MULTIPLY_BY_FREQS);
		}	
	}
	
	Tree dataTree 		  = treeString;
	LikelihoodFunction lf = (filteredData,dataTree);

	if (MPI_NODE_COUNT>1)
	{
		for (mpiNode = 0; mpiNode < MPI_NODE_COUNT-1; mpiNode = mpiNode+1)
		{
			if (MPINodeState[mpiNode]==0)
			{
				break;	
			}
		}
		
		if (mpiNode==MPI_NODE_COUNT-1)
		/* all nodes busy */
		{
			mpiNode = ReceiveJobs (1);
		}
		else
		{
			MPISend (mpiNode+1,lf);
			MPINodeState[mpiNode] = 1;
			MPINodeBounds[mpiNode][0] = currentStart;
			MPINodeBounds[mpiNode][1] = currentEnd;
		}
	}
	else
	{
		Optimize (lf_MLES,lf);
		segmentStart = currentStart;
		segmentEnd   = currentEnd;
		dummy = ReceiveJobs (0);
	}
}

if (MPI_NODE_COUNT>1)
{
	while (1)
	{
		for (nodeCounter = 0; nodeCounter < MPI_NODE_COUNT-1; nodeCounter = nodeCounter+1)
		{
			if (MPINodeState[nodeCounter]==1)
			{
				fromNode = ReceiveJobs (0);
				break;	
			}
		}
		if (nodeCounter == MPI_NODE_COUNT-1)
		{
			break;
		}
	}	
	OPTIMIZE_SUMMATION_ORDER = 1;
	
	for (nodeCounter = 0; nodeCounter*windowShift < ds.sites; nodeCounter=nodeCounter+1)
	{
		fprintf (LAST_FILE_PATH,"\n",Format(nodeCounter*windowShift,0,0)," , ",
									 Format(Min((1+nodeCounter)*windowShift,ds.sites-1),0,0),",",MPI_Results_Cache[dataDimension][nodeCounter]);	
		
		for (i=0;i<dataDimension;i=i+1)
		{
			fprintf (LAST_FILE_PATH,",",MPI_Results_Cache[i][nodeCounter]);
		}			
	}
}


function ReceiveJobs (sendOrNot)
{
	if (MPI_NODE_COUNT>1)
	{
		MPIReceive (-1, fromNode, result_String);
		
		segmentStart = MPINodeBounds[fromNode-1][0];
		segmentEnd   = MPINodeBounds[fromNode-1][1];
		
		if (sendOrNot)
		{
			MPISend (fromNode,lf);
			MPINodeBounds[fromNode-1][0] = currentStart;
			MPINodeBounds[fromNode-1][1] = currentEnd;
		}
		else
		{
			MPINodeState[fromNode-1] = 0;
			MPINodeBounds[fromNode-1][0] = -1;
			MPINodeBounds[fromNode-1][1] = -1;
		}
		ExecuteCommands (result_String);
	}
	
	if (firstTime)
	{
		firstTime = 0;
		
		dataDimension = Columns(lf_MLES);
		
		fprintf (LAST_FILE_PATH,"\nStarting bp,Ending bp,Ln-likelihood");
		for (i=0;i<dataDimension;i=i+1)
		{
			GetString (argName,lf,i);
			fprintf	  (LAST_FILE_PATH,",",argName);
		}
		
		if (MPI_NODE_COUNT>1)
		{
			MPI_Results_Cache = {dataDimension+1,ds.sites$windowShift+1};
		}
	}	
	
	fprintf (stdout,"Finshed Segment ",Format(segmentStart,0,0)," - ",Format(segmentEnd,0,0),"\n");
	if (MPI_NODE_COUNT<2)
	{
		fprintf (LAST_FILE_PATH,"\n",Format(segmentStart,0,0)," , ",Format(segmentEnd,0,0),",",lf_MLES[1][0]);
		
		for (i=0;i<dataDimension;i=i+1)
		{
			fprintf (LAST_FILE_PATH,",",lf_MLES[0][i]);
		}			
	}
	else
	{
		cacheIndex = segmentStart $ windowShift;
		for (i=0;i<dataDimension;i=i+1)
		{
			MPI_Results_Cache[i][cacheIndex] = lf_MLES[0][i];
		}			
		MPI_Results_Cache[i][cacheIndex] = lf_MLES[1][0];
	}
	return fromNode-1;
}


MESSAGE_LOGGING = 1;
