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
	/* assume that data set is filteredData, tree string is treeString */
	SAVE_GLOBALS = res[1][2];
	LRcounter = 0;
	if (SAVE_GLOBALS)
	{
		globalSpoolMatrix = {1,SAVE_GLOBALS};
		for (bsCounter = 0; bsCounter< SAVE_GLOBALS; bsCounter = bsCounter+1)
		{
			globalSpoolMatrix[bsCounter] = res1[0][bsCounter];
		}
	}
	dataDimension = Columns(res);
	dataMatrix = {2,2*dataDimension+3};
	nameWidth = 11;
	fprintf (tabulatedFileName,",Unconstrained Optimizations,");
	for (bsCounter=0;bsCounter<dataDimension;bsCounter=bsCounter+1)
	{
		fprintf	  (tabulatedFileName,", ");
	}
	fprintf (tabulatedFileName,"Constrained Optimizations");
	for (bsCounter=0;bsCounter<dataDimension;bsCounter=bsCounter+1)
	{
		fprintf	  (tabulatedFileName,", ");
	}
	fprintf (tabulatedFileName,"\nIteration,Ln-likelihood");
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
	fprintf (tabulatedFileName,",Ln-likelihood");
	for (bsCounter=0;bsCounter<dataDimension;bsCounter=bsCounter+1)
	{
		GetString (argName,lfConstrained,bsCounter);
		fprintf	  (tabulatedFileName,",",argName);
		temp = Abs(argName);
		if (temp>nameWidth)
		{
			nameWidth = temp;
		}
	}	
	fprintf (tabulatedFileName,",LR");
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
			dummy = SendIterate (0);
			ReplicateConstraint (constraintString,simulatedTree.Ingroup1, simulatedTree.Ingroup2);
			dummy = SendIterate (1);
			dummy = ManageMPIReturns (0);
		}
		else
		{
			Optimize (simulatedResults,simulatedLF);
			ReplicateConstraint (constraintString,simulatedTree.Ingroup1, simulatedTree.Ingroup2);
			Optimize (simulatedResults2,simulatedLF);
			fromNode = HandleAnIterate (0);	
		}

		if (SAVE_GLOBALS)
		{
			for (i = 0; i< SAVE_GLOBALS; i = i+1)
			{
				SetParameter(lfConstrained,i,globalSpoolMatrix[i]);
			}
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

	fprintf (tabulatedFileName,"\nMLE,",res[1][0]);
	
	temp = bsIterates-1;
	
	fprintf (stdout,"\n\n\t\tBOOTSTRAPPING SUMMARY\n\n\t\tUNCONSTRAINED OPTIMIZATION\n\n");
	
	tableWidth = nameWidth+49;
	fprintf (stdout,"+");
	dumb = PadString (tableWidth-2,"-");
	fprintf (stdout,"+\n| Parameter");
	dumb = PadString (nameWidth-9," ");
	fprintf (stdout,"|    MLE      |     Mean     |    Variance    |");
	fprintf (stdout,"\n+");
	dumb = PadString (tableWidth-2,"-");
	fprintf (stdout,"+");
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
		dumb = PadString(-Abs(argName)+nameWidth," ");
		fprintf (stdout,"|",Format(mle,12,6)," | ",Format(dataMatrix[0][i]/bsIterates,12,6)," | ",
						    Format((dataMatrix[1][i]-dataMatrix[0][i]*dataMatrix[0][i]/bsIterates)/temp,14,7)," |");
	}
	fprintf (stdout,"\n+");
	dumb = PadString (tableWidth-2,"-");
	fprintf (stdout,"+\n\n");

	fprintf (stdout,"\n\n\t\tRATE CONSTRAINED OPTIMIZATION\n\n");
	
	tableWidth = nameWidth+49;
	fprintf (stdout,"+");
	dumb = PadString (tableWidth-2,"-");
	fprintf (stdout,"+\n| Parameter");
	dumb = PadString (nameWidth-9," ");
	fprintf (stdout,"|    MLE      |     Mean     |    Variance    |");
	fprintf (stdout,"\n+");
	dumb = PadString (tableWidth-2,"-");
	fprintf (stdout,"+");
	fprintf (tabulatedFileName,",",res1[1][0]);
	for (i=0;i<=dataDimension;i=i+1)
	{
		if (i<dataDimension)
		{
			GetString (argName,lfConstrained,i);
			fprintf (tabulatedFileName,",",res1[0][i]);
			mle = res1[0][i];
		}
		else
		{
			mle = res1[1][0];
			argName = "Ln-Lklhood";
		}
		fprintf (stdout,"\n| ",argName);
		dumb = PadString(-Abs(argName)+nameWidth," ");
		fprintf (stdout,"|",Format(mle,12,6)," | ",Format(dataMatrix[0][dataDimension+i+1]/bsIterates,12,6)," | ",
						    Format((dataMatrix[1][dataDimension+i+1]-dataMatrix[0][dataDimension+i+1]*dataMatrix[0][dataDimension+i+1]/bsIterates)/temp,14,7)," |");
	}
	fprintf (stdout,"\n+");
	dumb = PadString (tableWidth-2,"-");
	fprintf (stdout,"+\n\n");
	
	fprintf (tabulatedFileName,",",lnLikDiff);
	
	temp=2*dataDimension+2;
	fprintf (stdout,"\n\nLikelihood Ratio Statistics:\nMEAN=",dataMatrix[0][temp]/bsIterates,
			"\nVARIANCE=",(dataMatrix[1][temp]-dataMatrix[0][temp]*dataMatrix[0][temp]/bsIterates)/(bsIterates-1),
			"\nProportion larger that the original likelihood ratio=",1-LRcounter/bsIterates,"\n");

	return 0;
}

function HandleAnIterate (dummy)
{
	fprintf (stdout,"\n",Format(iteratesDoneCount,8,0));
	fprintf (tabulatedFileName,"\n",Format(iteratesDoneCount,0,0));
	fprintf (tabulatedFileName,",",simulatedResults[1][0]);
	for (i=0;i<dataDimension;i=i+1)
	{
		temp = simulatedResults[0][i];
		dataMatrix[0][i]=dataMatrix[0][i]+temp;
		dataMatrix[1][i]=dataMatrix[1][i]+temp*temp;
		fprintf (tabulatedFileName,",",temp);
	}
	temp = simulatedResults[1][0];
	dataMatrix[0][i]=dataMatrix[0][i]+temp;
	dataMatrix[1][i]=dataMatrix[1][i]+temp*temp;
	fprintf (tabulatedFileName,",",simulatedResults2[1][0]);
	for (i=0;i<dataDimension;i=i+1)
	{
		temp = simulatedResults2[0][i];
		dataMatrix[0][dataDimension+1+i]=dataMatrix[0][dataDimension+1+i]+temp;
		dataMatrix[1][dataDimension+1+i]=dataMatrix[1][dataDimension+1+i]+temp*temp;
		fprintf (tabulatedFileName,",",temp);
	}
	temp = simulatedResults2[1][0];
	dataMatrix[0][dataDimension+1+i]=dataMatrix[0][dataDimension+1+i]+temp;
	dataMatrix[1][dataDimension+1+i]=dataMatrix[1][dataDimension+1+i]+temp*temp;
	
	/* Likelihood Ratio here */
	simLR = 2*(simulatedResults[1][0]-simulatedResults2[1][0]);
	if (simLR<0)
	{
		fprintf (MESSAGE_LOG,"\nA negative LRT statistic encoutered. You may want to increase the optimization precision settings to resolve numerical apporximation errors");
	}
	if (simLR<lnLikDiff)
	{
		LRcounter=LRcounter+1;
	}
	fprintf (tabulatedFileName,",",simLR);
	temp = 2*dataDimension+2;
	dataMatrix[0][temp]=dataMatrix[0][temp]+simLR;
	dataMatrix[1][temp]=dataMatrix[1][temp]+simLR*simLR;
	
	fprintf (stdout, "  ",Format(simLR,12,7),
					  "    ",Format(1-LRcounter/iteratesDoneCount,12,7));
					  
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
