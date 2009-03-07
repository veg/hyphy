IS_BOOTSTRAP_AVAILABLE = 1;
IS_NPBOOTSTRAP_AVAILABLE = 1;

function PadString (padLength,padChar)
{
	for (padCounter=0;padCounter<padLength;padCounter=padCounter+1)
	{
		fprintf (stdout,padChar);
	}
	return padLength;
}

function BootStrapFunction (bsIterates, tabulatedFileName, parametricOrNot)
{
	VERBOSITY_LEVEL = -1;
	/* assume that data set is filteredData, tree string is treeString */
	LRCounter = 0;
	LRsum = 0;
	LRsquaresum = 0;
	fprintf (tabulatedFileName,"Iteration,Likelihood Ratio");
	
	fprintf (stdout,"\nIteration      LRT       Current P-Value\n");

	if (MPI_NODE_COUNT>1)
	{
		MPINodeState   = {MPI_NODE_COUNT-1,3};
		MPIResultCache = {bsIterates, 2};
		for (bsCounter = 0; bsCounter<bsIterates; bsCounter=bsCounter+1)
		{
			MPIResultCache[bsCounter][0] = "";
			MPIResultCache[bsCounter][1] = "";
		}
		
		iteratesScanFrom  = 0;
	}

	iteratesDoneCount = 1;
	
	for (bsCounter=1;bsCounter<=bsIterates;bsCounter=bsCounter+1)
	{
		if (parametricOrNot)
		{
			if (SAVE_GLOBALS)
			{
				for (i = 0; i< SAVE_GLOBALS; i = i+1)
				{
					SetParameter(lfConstrained,i,globalSpoolMatrix[i]);
				}
			}
			DataSet simulatedDataSet = SimulateDataSet (lfConstrained);
			if (NICETY_LEVEL == 3)
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
			if (NICETY_LEVEL == 3)
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
			if (NICETY_LEVEL==3)
			{
				simulatedCodonEFV = BuildCodonFrequencies (simulatedEFV);
				Model simulatedModel = (simulatedModelMatrix,simulatedCodonEFV,MULTIPLY_BY_FREQS);	
				Tree simulatedTree = treeString;
			}
			else
			{
				Model simulatedModel = (simulatedModelMatrix,simulatedEFV,MULTIPLY_BY_FREQS);
				Tree simulatedTree = treeString;
			}	

		}
		else
		{
			Tree simulatedTree = treeString;
		}
		
		LikelihoodFunction simulatedLF = (simulatedDataFilter,simulatedTree);
		if (MPI_NODE_COUNT>1)
		{
			dummy = SendIterate (0);
			ExecuteCommands ("MolecularClock (simulatedTree,"+parameter2ConstrainString+");");
			dummy = SendIterate (1);
			dummy = ManageMPIReturns (0);
		}
		else
		{
			Optimize (simulatedResults,simulatedLF);
			ExecuteCommands ("MolecularClock (simulatedTree,"+parameter2ConstrainString+");");
			Optimize (simulatedResults2,simulatedLF);
			fromNode = HandleAnIterate (0);	
		}
				
	}
	
	if (MPI_NODE_COUNT>1)
	{
		while (1)
		{
			for (nodeCounter = 0; nodeCounter < MPI_NODE_COUNT-1; nodeCounter = nodeCounter+1)
			{
				if (MPINodeState[nodeCounter][0]==1)
				{
					fromNode = ReceiveJobs (0,0);
					break;	
				}
			}
			if (nodeCounter == MPI_NODE_COUNT-1)
			{
				break;
			}
		}	
		OPTIMIZE_SUMMATION_ORDER = 1;
		dummy = ManageMPIReturns (0);
	}	
	
	fprintf (tabulatedFileName,"\nMLE,",lnLikDiff);
	fprintf (stdout,"\n\n\t\tBOOTSTRAPPING SUMMARY\n\n");
	
	fprintf (stdout,"\n\nLikelihood Ratio Statistics:\nMEAN     = ",Format(LRsum/bsIterates,12,7),
		    "\nVARIANCE = ",Format((LRsquaresum-LRsum*LRsum/bsIterates)/(bsIterates-1),12,7),
			"\nProportion larger that the original likelihood ratio=",1-LRCounter/bsIterates,"\n");
			
	return 0;
}

function HandleAnIterate (dummy)
{
	fprintf (stdout,"\n",Format(iteratesDoneCount,8,0));
	fprintf (tabulatedFileName,"\n",Format(iteratesDoneCount,0,0));
	SimLR = 2*(simulatedResults[1][0]-simulatedResults2[1][0]);
	if (SimLR<0)
	{
		fprintf (MESSAGE_LOG,"\nA negative LRT statistic encoutered. You may want to increase the optimization precision settings to resolve numerical apporximation errors");
	}
	if (SimLR<lnLikDiff)
	{
		LRCounter = LRCounter+1;
	}
	LRsum = LRsum+SimLR;
	LRsquaresum = LRsquaresum+SimLR*SimLR;
	fprintf (tabulatedFileName,",",SimLR);
	fprintf (stdout, "  ",Format(SimLR,12,7),
					  "    ",Format(1-LRCounter/iteratesDoneCount,12,7));
					  
	iteratesDoneCount = iteratesDoneCount + 1;

	return 0;
}

function ReceiveJobs (sendOn,hyp)
{
	MPIReceive (-1, fromNode, result_String);
	jobIterateCount = MPINodeState[fromNode-1][1];	
	jobHypothesis   = MPINodeState[fromNode-1][2];
	
	if (sendOn)
	{	
		MPISend (fromNode,simulatedLF);
		MPINodeState[fromNode-1][1] = bsCounter-1;	
		MPINodeState[fromNode-1][2] = hyp;	
	}
	else
	{
		MPINodeState[fromNode-1][0] = 0;		
	}
	
	MPIResultCache[jobIterateCount][jobHypothesis] = result_String;
	result_String = "";
	
	return 0;
}

function SendIterate (dummy)
{
	for (mpiNode = 0; mpiNode < MPI_NODE_COUNT-1; mpiNode = mpiNode+1)
	{
		if (MPINodeState[mpiNode][0]==0)
		{
			break;	
		}
	}
	
	if (mpiNode==MPI_NODE_COUNT-1)
	/* all nodes busy */
	{
		mpiNode = ReceiveJobs (1,dummy);
	}
	else
	{
		MPISend (mpiNode+1,simulatedLF);
		MPINodeState[mpiNode][0] = 1;
		MPINodeState[mpiNode][1] = bsCounter-1;
		MPINodeState[mpiNode][2] = dummy;
	}
	return 0;
}

function ManageMPIReturns (dummy)
{
	incrementLowerBound = 1;
	for (jobIterateCount = iteratesScanFrom; jobIterateCount < bsCounter-1; jobIterateCount = jobIterateCount+1)
	{
		if (Abs(MPIResultCache[jobIterateCount][0])*Abs(MPIResultCache[jobIterateCount][1]))
		{
			if (incrementLowerBound)
			{
				iteratesScanFrom = iteratesScanFrom+1;
			}
			ExecuteCommands (MPIResultCache[jobIterateCount][0]);
			simulatedResults = simulatedLF_MLES;
			ExecuteCommands (MPIResultCache[jobIterateCount][1]);
			simulatedResults2 = simulatedLF_MLES;
			MPIResultCache[jobIterateCount][0] = "";
			MPIResultCache[jobIterateCount][1] = "";
			dummy = HandleAnIterate (0);
		}
		else
		{
			incrementLowerBound = 0;
		}
	}
	return 0;
}

