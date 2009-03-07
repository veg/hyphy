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
	LRCounter = 0;
	LRsum = 0;
	LRsquaresum = 0;
	fprintf (tabulatedFileName,"Iteration,Likelihood Ratio");
	
	fprintf (stdout,"\nIteration      LRT       Current P-Value\n");
	for (bsCounter=1;bsCounter<=bsIterates;bsCounter=bsCounter+1)
	{
		fprintf (stdout,"\n",Format(bsCounter,8,0));
		fprintf (tabulatedFileName,"\n",Format(bsCounter,0,0));
		if (parametricOrNot)
		{
			if (SAVE_GLOBALS2)
			{
				for (i = 0; i< SAVE_GLOBALS2; i = i+1)
				{
					SetParameter(lfConstrained,i,globalSpoolMatrix2[i]);
				}
			}
			DataSet simulatedDataSet = SimulateDataSet (lfConstrained);

			if (NICETY_LEVEL == 3)
			{
				DataSetFilter simulatedDataFilter  = CreateFilter (simulatedDataSet,3,siteIndex<ds.sites,"",GeneticCodeExclusions);
				DataSetFilter simulatedDataFilter2 = CreateFilter (simulatedDataSet,3,siteIndex>=ds.sites,"",GeneticCodeExclusions);
			}
			else
			{
				DataSetFilter simulatedDataFilter  = CreateFilter (simulatedDataSet,1,siteIndex<ds.sites);
				DataSetFilter simulatedDataFilter2 = CreateFilter (simulatedDataSet,1,siteIndex>=ds.sites);
			}
		}	
		else
		{
			if (NICETY_LEVEL == 3)
			{
				DataSetFilter simulatedDataFilter  = Bootstrap (filteredData,3);
				DataSetFilter simulatedDataFilter2 = Bootstrap (filteredData2,3);
			}
			else
			{
				DataSetFilter simulatedDataFilter  = Bootstrap (filteredData,1);
				DataSetFilter simulatedDataFilter2 = Bootstrap (filteredData2,1);
			}
		}	
		HarvestFrequencies (simulatedEFV,simulatedDataFilter,1,1,1);
		HarvestFrequencies (simulatedEFV2,simulatedDataFilter2,1,1,1);

		if (FREQUENCY_SENSITIVE)
		{
			simulatedModelMatrix = 0;
			simulatedModelMatrix2 = 0;
			if (USE_POSITION_SPECIFIC_FREQS)
			{
				HarvestFrequencies (simulatedEFV,simulatedDataFilter,3,1,1);
				HarvestFrequencies (simulatedEFV2,simulatedDataFilter2,3,1,1);
			}
			MULTIPLY_BY_FREQS = PopulateModelMatrix ("simulatedModelMatrix",simulatedEFV);
			MULTIPLY_BY_FREQS2 = PopulateModelMatrix ("simulatedModelMatrix2",simulatedEFV2);
			if (NICETY_LEVEL==3)
			{
				simulatedCodonEFV = BuildCodonFrequencies (simulatedEFV);
				simulatedCodonEFV2 = BuildCodonFrequencies (simulatedEFV2);
				Model simulatedModel = (simulatedModelMatrix,simulatedCodonEFV,MULTIPLY_BY_FREQS);	
				Tree simulatedTree = treeString;
				Model simulatedModel2 = (simulatedModelMatrix2,simulatedCodonEFV2,MULTIPLY_BY_FREQS);	
				Tree simulatedTree2 = treeString;
			}
			else
			{
				Model simulatedModel = (simulatedModelMatrix,simulatedEFV,MULTIPLY_BY_FREQS);
				Tree simulatedTree = treeString;
				Model simulatedModel2 = (simulatedModelMatrix2,simulatedEFV2,MULTIPLY_BY_FREQS);
				Tree simulatedTree2 = treeString;
			}	

		}
		else
		{
			Tree simulatedTree = treeString;
			Tree simulatedTree2 = treeString;
		}
		
		LikelihoodFunction simulatedLF = (simulatedDataFilter,simulatedTree);
		Optimize (simulatedResults,simulatedLF);
		LikelihoodFunction simulatedLF2 = (simulatedDataFilter2,simulatedTree2);
		Optimize (simulatedResults1,simulatedLF2);
		ReplicateConstraint (constraintString,simulatedTree2, simulatedTree);
		LikelihoodFunction lfSimConstrained = (simulatedDataFilter2,simulatedTree2,simulatedDataFilter,simulatedTree);
		Optimize (simulatedResults2,lfSimConstrained);
		SimLR = 2*(simulatedResults[1][0]+simulatedResults1[1][0]-simulatedResults2[1][0]);
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
						  "    ",Format(1-LRCounter/bsCounter,12,7));
				
	}
	fprintf (tabulatedFileName,"\nMLE,",lnLikDiff);
	fprintf (stdout,"\n\n\t\tBOOTSTRAPPING SUMMARY\n\n");
	
	fprintf (stdout,"\n\nLikelihood Ratio Statistics:\nMEAN=",LRsum/bsIterates,
			"\nVARIANCE=",(LRsquaresum-LRsum*LRsum/bsIterates)/(bsIterates-1),
			"\nProportion larger that the original likelihood ratio=",1-LRCounter/bsIterates,"\n");
			
	return 0;
}

