IS_BOOTSTRAP_AVAILABLE 		= 1;
IS_NPBOOTSTRAP_AVAILABLE    = 1;

function PadString (padLength,padChar)
{
	for (padCounter=0;padCounter<padLength;padCounter=padCounter+1)
	{
		fprintf (stdout,padChar);
	}
	return padLength;
}

/* mod 11/10/2005 to spool branch lengths as well */
#include 	"dSdNTreeTools.ibf";
ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Utility"+DIRECTORY_SEPARATOR+"CodonTools.bf");
/* end mod */

/*----------------------------------------------------------------------------------------------------------------------------*/

function BootStrapFunction (bsIterates, tabulatedFileName, parametricOrNot)
{
	VERBOSITY_LEVEL				= -1;
	/* assume that data set is filteredData, tree string is treeString */
	SAVE_GLOBALS = res[1][2];
	if (SAVE_GLOBALS)
	{
		globalSpoolMatrix = {1,SAVE_GLOBALS};
		for (bsCounter = 0; bsCounter< SAVE_GLOBALS; bsCounter = bsCounter+1)
		{
			globalSpoolMatrix[bsCounter] = res[0][bsCounter];
		}
	}
		
	if (categoriesUsed&&parametricOrNot)
	{
		ChoiceList (CATEGORY_SIMULATION_METHOD,"Sampling Options",1,SKIP_NONE,
					"Discrete","Random rate variability will be sampled from the discretized disrtibution.",
					"Continuous","Random rate variability will be sampled from the continuous disrtibution (if possible).",
		);
		if (CATEGORY_SIMULATION_METHOD<0)
		{
			return 0;
		}
		else
		{
			CATEGORY_SIMULATION_METHOD = CATEGORY_SIMULATION_METHOD+1;		
		}
	}	
	
	/* mod 11/10/2005 */
	GetString (lf_summary,lf,-1);
	tree_ids = lf_summary ["Trees"];
	generatingTree = tree_ids[0];
	ExecuteCommands(
	"BranchNamesVector = BranchName ("+generatingTree+", -1);"
	);
	branchTotalCount = Columns (BranchNamesVector)-1; /* do not count the 'root' node */
	/* end mod */
	
	fprintf 		(tabulatedFileName,KEEP_OPEN,"Iteration,Ln-likelihood");
	dataDimension = Columns(res);
	dataMatrix 	  = {2,dataDimension+1};

	/* map likelihood function parameters to indices */
	_variableMap  = {};
	for (bsCounter=0; bsCounter < dataDimension; bsCounter = bsCounter + 1)
	{
		GetString (_i, lf, bsCounter);
		_i = _i^{{"givenTree\\.",""}};
		_variableMap[_i] = Abs (_variableMap);
	}
	
	nameWidth = 11;
	if (Abs(_Genetic_Code)) /* have coding data */
	{
		_baseAVL = _computeSNSSites ("filteredData", _Genetic_Code, vectorOfFrequencies, 0);		
		fprintf	  (tabulatedFileName,",Syn_sites,NS_sites,Total_sites");
	}
	for (bsCounter=0;bsCounter<dataDimension;bsCounter=bsCounter+1)
	{
		GetString (argName,lf,bsCounter);
		fprintf	  (tabulatedFileName,",",argName);
		temp = Abs(argName);
		if (temp>nameWidth)
		{
			nameWidth = temp;
		}
	}
	/* mod 11/10/2005 */
	if (Abs(_Genetic_Code)) /* have coding data */
	{
		scalingStencils = ComputeScalingStencils (0); /* defined in dSdNTreeTools.ibf */
		for (bsCounter=0;bsCounter<branchTotalCount;bsCounter=bsCounter+1)
		{
			fprintf	  (tabulatedFileName,",TotalLength(",BranchNamesVector[bsCounter],"),",
										 "SynLength(",BranchNamesVector[bsCounter],"),",
										 "NonSynLength(",BranchNamesVector[bsCounter],")");
		}	
		BaseBranchLengthsVector = ReturnVectorsOfCodonLengths (scalingStencils,generatingTree);
	}
	else
	{
		for (bsCounter=0;bsCounter<branchTotalCount;bsCounter=bsCounter+1)
		{
			fprintf	  (tabulatedFileName,",BranchLength(",BranchNamesVector[bsCounter],")");
		}
		ExecuteCommands(
		"BaseBranchLengthsVector = BranchLength ("+generatingTree+", -1);"
		);		
	}
	
	/* end mod */
	
	if (MPI_NODE_COUNT>1)
	{
		MPINodeState = {MPI_NODE_COUNT-1,1};
		finishedJobs = 0;
		OPTIMIZE_SUMMATION_ORDER = 0;
	}
	
	for (bsCounter=1;bsCounter<=bsIterates;bsCounter=bsCounter+1)
	{
		if (MPI_NODE_COUNT<=1)
		{
			fprintf (stdout,"\nIteration ",Format(bsCounter,0,0),"/",Format(bsIterates,0,0)," ");
			fprintf (tabulatedFileName,"\n",Format(bsCounter,0,0));		
		}
		
		if (parametricOrNot)
		{
			DataSet simulatedDataSet = SimulateDataSet (lf);
			if (Abs(_Genetic_Code))
			{
				DataSetFilter simulatedDataFilter = CreateFilter (simulatedDataSet,3,"","",GeneticCodeExclusions);
			}
			else
			{
				DataSetFilter simulatedDataFilter = CreateFilter (simulatedDataSet,1);
			}
		}	
		else
		{
			if (Abs(_Genetic_Code))
			{
				DataSetFilter simulatedDataFilter = Bootstrap (filteredData,3);
			}
			else
			{
				DataSetFilter simulatedDataFilter = Bootstrap (filteredData,1);
			}
		}	
		
		HarvestFrequencies (simulatedEFV,simulatedDataFilter,1,1,1);
		
		if (FREQUENCY_SENSITIVE)
		{
			simulatedModelMatrix = 0;
			if (USE_POSITION_SPECIFIC_FREQS)
			{
				HarvestFrequencies (simulatedEFV,simulatedDataFilter,3,1,1);
			}
			MULTIPLY_BY_FREQS = PopulateModelMatrix ("simulatedModelMatrix",simulatedEFV);
			if (Abs(_Genetic_Code))
			{
				simulatedCodonEFV 	 = BuildCodonFrequencies (simulatedEFV);
				Model simulatedModel = (simulatedModelMatrix,simulatedCodonEFV,MULTIPLY_BY_FREQS);	
			}
			else
			{
				Model simulatedModel = (simulatedModelMatrix,simulatedEFV,MULTIPLY_BY_FREQS);
			}	
		}
		
		Tree simulatedTree = treeString;
		LikelihoodFunction simulatedLF = (simulatedDataFilter,simulatedTree);
		
		if (MPI_NODE_COUNT>1)
		{
			for (nodeCounter = 0; nodeCounter < MPI_NODE_COUNT-1; nodeCounter = nodeCounter+1)
			{
				if (MPINodeState[nodeCounter]==0)
				{
					break;	
				}
			}
			
			if (nodeCounter==MPI_NODE_COUNT-1)
			{
				nodeCounter = ProcessIterate (1);
			}
			else
			{
				MPISend (nodeCounter+1,simulatedLF);
				MPINodeState[nodeCounter] = 1;
			}		
		}
		else
		{
			Optimize (simulatedResults,simulatedLF);
			dummy = ProcessIterate (0);
		}
		
		if (SAVE_GLOBALS)
		{
			for (i = 0; i< SAVE_GLOBALS; i = i+1)
			{
				SetParameter(lf,i,globalSpoolMatrix[i]);
			}
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
					fromNode = ProcessIterate (0);
					break;	
				}
			}
			if (nodeCounter == MPI_NODE_COUNT-1)
			{
				break;
			}
		}		
		OPTIMIZE_SUMMATION_ORDER = 1;
	}

	fprintf (tabulatedFileName,"\nMLE,",res[1][0]);
	
	temp = bsIterates-1;
	
	fprintf (stdout,"\n\n\t\tBOOTSTRAPPING SUMMARY\n\n");
	
	tableWidth = nameWidth+49;
	fprintf (stdout,"+");
	PadString (tableWidth-2,"-");
	fprintf (stdout,"+\n| Parameter");
	PadString (nameWidth-9," ");
	fprintf (stdout,"|    MLE      |     Mean     |    Variance    |");
	fprintf (stdout,"\n+");
	PadString (tableWidth-2,"-");
	fprintf (stdout,"+");

	if (Abs(_Genetic_Code)) /* have coding data */
	{
		fprintf (tabulatedFileName,",",_baseAVL["SSites"],",",_baseAVL["NSSites"],",",_baseAVL["Sites"]);
	}

	for (i=0;i<=dataDimension;i=i+1)
	{
		if (i<dataDimension)
		{
			GetString (argName,lf,i);
			fprintf (tabulatedFileName,",",res[0][i]);
			mle = res[0][i];
		}
		else
		{
			mle = res[1][0];
			argName = "Ln-Lklhood";
		}
		fprintf (stdout,"\n| ",argName);
		PadString(-Abs(argName)+nameWidth," ");
		fprintf (stdout,"|",Format(mle,12,6)," | ",Format(dataMatrix[0][i]/bsIterates,12,6)," | ",
						    Format((dataMatrix[1][i]-dataMatrix[0][i]*dataMatrix[0][i]/bsIterates)/temp,14,7)," |");
	}
	/* mod 11/10/2005; print MLE branch lengths */
	blString = "";
	blString * 128;
	if (Abs(_Genetic_Code)) /* have coding data */
	{
		tv = BaseBranchLengthsVector["Total"];
		sv = BaseBranchLengthsVector["Syn"];
		nv = BaseBranchLengthsVector["NonSyn"];
		
		for (bsCounter=0;bsCounter<branchTotalCount;bsCounter=bsCounter+1)
		{
			blString * (","+tv[bsCounter]+","+sv[bsCounter]+","+nv[bsCounter]);
		}
	}
	else
	{
		for (bsCounter=0;bsCounter<branchTotalCount;bsCounter=bsCounter+1)
		{
			blString * (","+BaseBranchLengthsVector[bsCounter]);
		}
	}
	blString * 0;
	fprintf (tabulatedFileName, blString,CLOSE_FILE);
	/* end */
	fprintf (stdout,"\n+");
	dumb = PadString (tableWidth-2,"-");
	fprintf (stdout,"+\n\n");
	
	return 0;
}

/*----------------------------------------------------------------------------------*/

function ProcessIterate (sendMore)
{
	if (MPI_NODE_COUNT>1)
	{
		finishedJobs = finishedJobs + 1;
		MPIReceive (-1, fromNode, result_String);
		if (sendMore)
		{
			MPISend (fromNode,simulatedLF);
		}
		else
		{
			MPINodeState[fromNode-1] = 0;
		}
		ExecuteCommands (result_String);
		simulatedResults = simulatedLF_MLES;
		fprintf (stdout,"\nIteration ",Format(finishedJobs,0,0),"/",Format(bsIterates,0,0)," ");
		fprintf (tabulatedFileName,"\n",Format(finishedJobs,0,0));		
	}
	
	
	/* mod 11/10/2005; print MLE branch lengths */
	blString = "";
	blString * 128;
	blString * (","+simulatedResults[1][0]);
	
	tMap = {dataDimension,1};
	for (i=0;i<dataDimension;i=i+1)
	{
		GetString (_i, simulatedLF, i);
		_i = _variableMap [_i^{{"simulatedTree\\.",""}}];
		temp 				= 	simulatedResults[0][i];
		dataMatrix[0][_i]	=	dataMatrix[0][_i]	+	temp;
		dataMatrix[1][_i]	=	dataMatrix[1][_i]	+	temp*temp;
		tMap [_i] = temp;
	}
	
	if (Abs(_Genetic_Code)) /* have coding data */
	{
		if (FREQUENCY_SENSITIVE)
		{
			_snsAVL	   			 = _computeSNSSites ("simulatedDataFilter", _Genetic_Code, simulatedCodonEFV, 1);
		}
		else
		{
			_snsAVL				 = _baseAVL;
		}
		blString * (","+_snsAVL["SSites"]+","+_snsAVL["NSSites"]+","+_snsAVL["Sites"]);
	}

	for (i=0;i<dataDimension;i=i+1)
	{
		blString * (","+tMap [i]);
	}

	if (Abs(_Genetic_Code)) /* have coding data */
	{
		BranchLengthsVector = ReturnVectorsOfCodonLengths (scalingStencils,"simulatedTree");
		tv = BranchLengthsVector["Total"];
		sv = BranchLengthsVector["Syn"];
		nv = BranchLengthsVector["NonSyn"];
		
		for (i=0;i<branchTotalCount;i=i+1)
		{
			blString * (","+tv[i]+","+sv[i]+","+nv[i]);
		}
		
	}
	else
	{
		BranchLengthsVector = BranchLength (simulatedTree,-1);
		for (i=0;i<branchTotalCount;i=i+1)
		{
			blString * (","+BranchLengthsVector[i]);
		}
	}
	blString * 0;
	fprintf (tabulatedFileName, blString);
	/* end */

	temp = simulatedResults[1][0];
	dataMatrix[0][dataDimension]=dataMatrix[0][dataDimension]+temp;
	dataMatrix[1][dataDimension]=dataMatrix[1][dataDimension]+temp*temp;
	
	return fromNode-1;
}
