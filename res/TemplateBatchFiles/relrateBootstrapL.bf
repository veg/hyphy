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
	SAVE_GLOBALS = res1[1][2];
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
	dataMatrix = {2,2*dataDimension+4};
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
	for (bsCounter=1;bsCounter<=bsIterates;bsCounter=bsCounter+1)
	{
		fprintf (stdout,"\n",Format(bsCounter,8,0));
		fprintf (tabulatedFileName,"\n",Format(bsCounter,0,0));
		if (parametricOrNot)
		{
			DataSet simulatedDataSet = SimulateDataSet (lfConstrained);
			if (NICETY_LEVEL == 3)
			{
				DataSetFilter simulatedDataFilter = CreateFilter (simulatedDataSet,3,"","",GeneticCodeExclusions);
			}
		}	
		else
		{
			if (NICETY_LEVEL == 3)
			{
				DataSetFilter simulatedDataFilter = Bootstrap (filteredData,3);
			}
		}	

		if (modelType)
		{
			HarvestFrequencies (simulatedEFV,simulatedDataFilter,3,1,1);
		}
		else
		{
			HarvestFrequencies (simulatedEFV,simulatedDataFilter,1,1,1);
		}

		simulatedCodonEFV = BuildCodonFrequencies (simulatedEFV);
		MULTIPLY_BY_FREQS = PopulateModelMatrix ("simulatedModelMatrix",simulatedEFV);
		Model simulatedModel = (simulatedModelMatrix,simulatedCodonEFV,MULTIPLY_BY_FREQS);

		Tree simulatedTree = treeString;
		LikelihoodFunction simulatedLF = (simulatedDataFilter,simulatedTree);
		Optimize (simulatedResults,simulatedLF);
		
		simulatedTree.Ingroup1.nonSynRate := R*simulatedTree.Ingroup1.synRate;
		simulatedTree.Ingroup2.nonSynRate := R*simulatedTree.Ingroup2.synRate;

		Optimize (simulatedResults2,simulatedLF);
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
		for (i=0;i<=dataDimension;i=i+1)
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
		temp = 2*dataDimension+3;
		dataMatrix[0][temp]=dataMatrix[0][temp]+simLR;
		dataMatrix[1][temp]=dataMatrix[1][temp]+simLR*simLR;
		
		if (SAVE_GLOBALS)
		{
			for (i = 0; i< SAVE_GLOBALS; i = i+1)
			{
				SetParameter(lfConstrained,i,globalSpoolMatrix[i]);
			}
		}
		
		fprintf (stdout, "  ",Format(simLR,12,7),
						  "    ",Format(1-LRcounter/bsCounter,12,7));
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
	for (i=0;i<=dataDimension+1;i=i+1)
	{
		if (i<=dataDimension)
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
