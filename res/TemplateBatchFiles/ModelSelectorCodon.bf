fprintf (stdout, "\nStart with this many rate classes:");
fscanf  (stdin,  "Number", startWithRateClasses);

rateClassesCount = startWithRateClasses; 
autoStepFlag	 = 1;
maxRateClasses   = 75;
modelSpawnPrefix = "";
modelSpawnPrefix * 128;
modelSpawnSuffix = "";
modelSpawnSuffix * 128;
whichIC			 = 0;

orderString				= "FLIMVSPTAYHQNKDECWRG";

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

produceOffspring		= MPI_NODE_COUNT-1;
populationSize  		= 2*produceOffspring;
incestDistance  		= 0;
generationCount		  	= 5000;
maxSampleTries			= populationSize*10;
mutationThreshhold		= 0.0001;
mutationProb			= 0.15;
mutationProbDecrease	= 0.95;
annealingPhase			= 100;
SHORT_MPI_RETURN		= 1;
totalSampleCounter		= 0;
localMutationRate		= 0.05;
localMutationInterval	= 20;

stoppingCriterion		= 50;
sampleCount				= 0;
familyControlSize		= produceOffspring$6;

predef   				= {};

MasterList				= {};
verboseFlag				= 0;
stateVectorDimension    = 190;
branchUpdateThreshold   = 50;

predef   = {};
/* charge  */
predef[0] = {{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0}};
/* polarity  */
predef[1] = {{ 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1}};
/* polairty and size */
predef[2] = {{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
predef   = {};

modelComplexityPenalty = 0;

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function StirlingNumberOf2ndKind (N,k)
{
	SN  :< 1e300;
	kF	:< 1e300;
	CF  :< 1e300;
	
	SN  = 0;
	kF  = 1;
	m1t = 1;
	CF	= 1;
	
	for (i=0; i<k; i=i+1)
	{
		SN = SN + m1t*CF*(k-i)^N;
		m1t = -m1t;
		CF	= CF * (k-i)/(i+1);
		kF = kF*(i+1);
	}
	
	return SN/kF;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function returnIC (logL, df, sampleCount)
{
	if (whichIC)
	{
		return -2*logL + 2 * Log(sampleCount) * df + modelComplexityPenalty;
	}
	return -2*logL + 2*df * (sampleCount/(sampleCount-df-1));
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function ExportAMatrix (fileName, rateMatrix, doClear, doOrderString)
{
	outString = "";
	outString * 256;
	
	if (doOrderString)
	{
		outString * orderString;
	}
	
	outString * "\n{";
	for (h=0; h<21; h=h+1)
	{
		if (h!=10)
		{
			outString * "\n{";
			for (v=0; v<h; v=v+1)
			{
				if (v!=10)
				{
					if (v)
					{
						outString * (","+rateMatrix[h][v]);
					}
					else
					{
						outString * (""+rateMatrix[h][v]);				
					}
				}
			}
			if (v)
			{
				outString * ",*";
			}
			else
			{
				outString * "*";				
			}
			for (v=h+1; v<21; v=v+1)
			{
				if (v!=10)
				{
					outString * (","+rateMatrix[h][v]);				
				}
			}
			outString * "}";
		}
	}
	outString * "\n}";
	outString * 0;
	if (doClear)
	{
		fprintf (fileName, CLEAR_FILE,outString);
	}
	else
	{
		fprintf (fileName,"\n",outString);
	}
	return 0;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function MatrixToString (rateMatrix)
{
	outString = "";
	outString * 256;
	for (h=0; h<stateVectorDimension; h=h+1)
	{
		outString * (","+rateMatrix[h]);				
	}
	outString * 0;
	return outString;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function CleanUpMPI (dummy)
{
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
	}
	return 0;
}


/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function ReceiveJobs (sendOrNot, ji)
{
	if (MPI_NODE_COUNT>1)
	{
		MPIReceive (-1, fromNode, result_String);
		mji = MPINodeState[fromNode-1][1];
		mdf	= MPINodeState[fromNode-1][2];
		
		if (sendOrNot)
		{
			MPISend (fromNode,msg2Send);
			
			/*
			fprintf ("lastm.in", CLEAR_FILE, msg2Send);
			*/
			
			MPINodeState[fromNode-1][1] = ji;			
			MPINodeState[fromNode-1][2] = modelDF;			
		}
		else
		{
			MPINodeState[fromNode-1][0] = 0;
			MPINodeState[fromNode-1][1] = -1;		
		}
		ExecuteCommands (result_String);
		/*
		fprintf ("last.in", CLEAR_FILE, result_String);
		*/
		myDF	  = lf_MLES[1][1] + baseParams - 1;
		myAIC 	  = -returnIC (lf_MLES[1][0],myDF,sampleCount);
		myLFScore = lf_MLES[1][0];
		ji = mji;
	}
	else
	{
		myDF	  = modelDF + baseParams;
		myAIC 	  = -returnIC (lf_MLES[1][0],myDF,sampleCount);
		myLFScore = res[1][0];
	}
	
	if (resultProcessingContext==0)
	{
		sortedScores[ji][0] = myAIC;
	}
	else
	{
		intermediateProbs[ji][0] = myAIC;	
	}
	if (ji>=0)
	{
		if (resultProcessingContext==0)
		{
			ExportAMatrix (detailsOutputFile,StringToMatrix(currentPopulation[ji]),0,0);
			sortedBP = MatrixToString (currentPopulation[ji]);
			MasterList [sortedBP] = myAIC;
		}
		else
		{
			ExportAMatrix (detailsOutputFile,StringToMatrix(children[ji-populationSize]),0,0);	
			sortedBP = MatrixToString (children[ji-populationSize]);
			MasterList [sortedBP] = myAIC;
		}
		
		fprintf (detailsOutputFile,"\n{{",Format(myLFScore,25,10),",",myDF,",",Format(-myAIC,25,10));
		if (MPI_NODE_COUNT > 1)
		{
			for (spoolC = 0; spoolC < lf_MLES[1][1]; spoolC = spoolC + 1)
			{
				aVal = lf_MLE_VALUES["NSR"+spoolC];
				fprintf (detailsOutputFile,",",aVal);
			}
		}		
		else
		{
			for (spoolC = 0; spoolC < lf_MLES[1][1]; spoolC = spoolC + 1)
			{
				execCommand = "fprintf(detailsOutputFile,\",\",NSR"+spoolC+");";
				ExecuteCommands (execCommand);
			}
		}
		fprintf (detailsOutputFile,"}}\n");

		totalSampleCounter = totalSampleCounter + 1;
	}
	return fromNode-1;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function RunASample (modelDF, jobIndex)
{
	sampleString = MatrixToString (cString);
	myAIC = MasterList[sampleString];
	
	if (myAIC < (-0.1))
	{
		if (resultProcessingContext==0)
		{
			sortedScores[jobIndex][0] = myAIC;
		}
		else
		{
			intermediateProbs[jobIndex][0] = myAIC;	
		}			
		return 0;
	}
		
	msg2Send = "";
	msg2Send * 128;
	msg2Send * ("aaRateMultipliers = " + aaRateMultipliers + ";\n");
	msg2Send * ("mgPlainFactor = \"" + mgPlainFactor + "\";\n");
	msg2Send * ("presetBranchParameters = {}\n;");
	
	for (treeCount = 1; treeCount <= dataPartsRead; treeCount = treeCount + 1)
	{
		ExecuteCommands ("tcn = TipCount (codon_tree_"+treeCount+")-1; bcn = BranchCount (codon_tree_"+treeCount+")-1;");
		msg2Send * ("presetBranchParameters["+treeCount+"] = {}\n;");
		for (brCount = tcn; brCount >= 0; brCount = brCount - 1)
		{
			ExecuteCommands("tn = TipName (codon_tree_"+treeCount+", brCount);");
			ExecuteCommands ("ssr = codon_tree_"+treeCount+"." + tn + ".synRate;");
			msg2Send * ("(presetBranchParameters["+treeCount+"])[\"" + tn + "\"] = " + ssr + ";\n");
		}
		for (brCount = bcn; brCount >= 0; brCount = brCount - 1)
		{
			ExecuteCommands("tn = BranchName (codon_tree_"+treeCount+", brCount);");
			ExecuteCommands ("ssr = codon_tree_"+treeCount+"." + tn + ".synRate;");
			msg2Send * ("(presetBranchParameters["+treeCount+"])[\"" + tn + "\"] = " + ssr + ";\n");
		}
	}
	
	msg2Send * modelFreqSpec;
	msg2Send * modelSpawnPrefix;
	msg2Send * ("AC:="+AC+";\n");
	msg2Send * ("AT:="+AT+";\n");
	msg2Send * ("CG:="+CG+";\n");
	msg2Send * ("CT:="+CT+";\n");
	msg2Send * ("GT:="+GT+";\n");
	msg2Send * modelSpawnSuffix;
	msg2Send * 0;
	
	/* fprintf ("last.sent", CLEAR_FILE, msg2Send); */
	
	if ((MPI_NODE_COUNT>1) && (jobIndex>=0))
	{
		for (mpiNode = 0; mpiNode < MPI_NODE_COUNT-1; mpiNode = mpiNode+1)
		{
			if (MPINodeState[mpiNode][0]==0)
			{
				break;	
			}
		}
		if (mpiNode==MPI_NODE_COUNT-1)
		{
			mpiNode = ReceiveJobs (1,jobIndex);
		}
		else
		{
			MPISend (mpiNode+1,msg2Send);
			MPINodeState[mpiNode][0] = 1;
			MPINodeState[mpiNode][1] = jobIndex;
			MPINodeState[mpiNode][2] = modelDF;
		}
	}
	else
	{
		Optimize (res,lf);
		if (jobIndex>=0)
		{
			mpiNode = ReceiveJobs (1, jobIndex);
		}
		else
		{
			myAIC 	  = -returnIC (res[1][0],modelDF,sampleCount);
 		}
	}
	return 0;
}






/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function SpawnRandomString (clsCnt)
{
	rModel = {stateVectorDimension,1};
	for (h=0; h<stateVectorDimension; h=h+1)
	{
		rModel[h] = Random(0,clsCnt)$1;
	}
	
	return MakeStringCanonical(rModel,clsCnt);
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function MakeStringCanonical (randomModel, classCount)
{
	compressedString = {classCount,1};
	for (h=0; h<stateVectorDimension; h=h+1)
	{
		v = randomModel[h];
		compressedString[v] = 1;
	}
	compressedString[0] = 0;
	for (h=1; h<classCount; h=h+1)
	{
		compressedString[h] = compressedString[h]+compressedString[h-1];
	}
	for (h=0; h<stateVectorDimension; h=h+1)
	{
		v = randomModel[h];
		v = compressedString[v];
		randomModel[h] = v;
	}
	v = compressedString[classCount-1]+1;
	if (v>1)
	{
		sortedOrder = {v,2};
		for (h=0; h<v; h=h+1)
		{
			sortedOrder[h][0] = -1;
		}
		cc = 0;
		for (h=0; h<stateVectorDimension; h=h+1)
		{
			hshift = randomModel[h];
			if (sortedOrder[hshift][0] < 0)
			{
				sortedOrder[hshift][0] = cc;
				sortedOrder[hshift][1] = hshift;
				cc = cc+1;
			}
		}
		sortedOrder = sortedOrder%1;
		for (h=0; h<stateVectorDimension; h=h+1)
		{
			v = randomModel[h];
			randomModel[h] = sortedOrder[v][0];
		}
		
	}
	return randomModel;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function StringToMatrix (stringSpec)
{
	matrixSpec = {21,21};
	cc = 0;
	
	for (h=0; h<21; h=h+1)
	{
		for (v=h+1; v<21; v=v+1)
		{
			m_idx   = {{h,v}};
			map_idx = mapUpTo2121[m_idx];
			if (map_idx>0)
			{
				map_idx = 0+stringSpec[map_idx - 1];
				matrixSpec[h][v] = map_idx;
				matrixSpec[v][h] = map_idx;
			}
 		}
	}
	return matrixSpec;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function IsChildViable (putativeChild)
{
	sampleString 	= MatrixToString (putativeChild);
	myAIC 			= MasterList[sampleString];
	testChild 		= putativeChild;
	mutPassCount 	= 1;

	partStringsMain  = {};
	partStringsChild = {};
	
	for (_lc = 0; _lc < populationSize && myAIC > (-0.1); _lc = _lc + 1)
	{
		if (Rows (currentPopulation[_lc]))
		{
			if (partStringsMain[_lc] == 0)
			{
				partStringsMain[_lc] = MatrixToString (currentPopulation[_lc]);
			}
			myAIC = -(sampleString == partStringsMain[_lc]);
		}
	} 
	for (_lc = 0; _lc < Abs(children) && myAIC > (-0.1); _lc = _lc + 1)
	{
		if (Rows (children[_lc]))
		{	
			if (partStringsChild[_lc] == 0)
			{
				partStringsChild [_lc] = MatrixToString (children[_lc]);
			}
			myAIC = -(sampleString == partStringsChild [_lc]);
		}
	} 

	while (myAIC < (-0.1))
	{
		if (verboseFlag)
		{
			fprintf (stdout,"Adjusting the child to avoid a duplicate. Pass ", mutPassCount, "\n");
		}
		
		mutPassCount = mutPassCount + 1;
		
		sampleString = Min(Random(0,stateVectorDimension)$1,stateVectorDimension-1);
		myAIC 		 = testChild[sampleString];
		
		newValue 	 = Random (0,rateClassesCount-0.0000001)$1;
		
		while (newValue == myAIC)
		{
			newValue = Random (0,rateClassesCount-0.0000001)$1;
		}
		
		testChild [sampleString] = newValue;
		sampleString 			 = MatrixToString (testChild);
		myAIC 					 = MasterList[sampleString];

		if (myAIC > (-0.1))
		{
			for (_lc = 0; _lc < populationSize && myAIC > (-0.1); _lc = _lc + 1)
			{
				if (Rows (currentPopulation[_lc]))
				{
					myAIC = -(sampleString == partStringsMain[_lc]);
				}
			} 
			for (_lc = 0; _lc < Abs(children) && myAIC > (-0.1); _lc = _lc + 1)
			{
				if (Rows (children[_lc]))
				{	
					myAIC = -(sampleString == partStringsChild [_lc]);
				}
			} 
		
		}
	}
	return testChild;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function UpdateBL (dummy)
{
	bestInd							= sortedScores[populationSize-1][1];
	aaRateMultipliers 				= StringToMatrix(currentPopulation[bestInd]);
	MULTIPLY_BY_FREQS 				= PopulateModelMatrix ("MG94MLEFREQS", paramFreqs,1);
	Model MG94MLEFREQS_MODEL	    = ( MG94MLEFREQS, paramCodonFreqs, 0);

	AC=AC;
	AT=AT;
	CG=CG;
	CT=CT;
	GT=GT;
	
	populateTrees ("codon_tree", dataPartsRead);
	for (k=1; k<=dataPartsRead; k=k+1)
	{
		ExecuteCommands("ReplicateConstraint (\"this1.?.synRate:=3*this2.?.t\", codon_tree_"+k+", nuc_tree_"+k+");ClearConstraints(codon_tree_"+k+");");
	}
	ExecuteCommands(constructLF("lf","filteredData","codon_tree",dataPartsRead));
	if (verboseFlag)
	{
		VERBOSITY_LEVEL 			= 1;
	}
	Optimize 						(res,lf);
	VERBOSITY_LEVEL 				= 0;
	
	AC:=AC__;
	AT:=AT__;
	CG:=CG__;
	CT:=CT__;
	GT:=GT__;
	
	for (k = 0; k < 4; k = k+1)
	{
		for (p = 0; p < 3; p=p+1)
		{
			ExecuteCommands ("positionFrequencies[k][p] = " + paramFreqs[k][p]);
		}
	}
		
	vectorOfFrequencies                     = BuildCodonFrequencies (positionFrequencies);
	modelFreqSpec							= "observedFreq = " + positionFrequencies + ";\nvectorOfFrequencies = " + vectorOfFrequencies + ";\n";
	global									 mgPlainFactor = FindBranchLengthExpression(0,"MG94MLEFREQS");
	
	myLFScore								 = res[1][0];
	myDF									 = res[1][1];	
	myAIC 	  								 = -returnIC(res[1][0],myDF,sampleCount);
	
	fprintf									 (stdout, "\nUpdated BLs\nDF=", myDF,"\nDelta IC = ", Format(myAIC-sortedScores[populationSize-1][0],20,5), "\n", lf, "\n", positionFrequencies, "\n", mgPlainFactor, "\n");
	sortedScores[populationSize-1][0] 		 = myAIC;
	MasterList[sampleString] 				 = myAIC;
			
	ExportAMatrix (detailsOutputFile,StringToMatrix(currentPopulation[bestInd]),0,0);	
	fprintf (detailsOutputFile,"\n{{",Format(myLFScore,25,10),",",myDF,",",Format(-myAIC,25,10));

	for (spoolC = 0; spoolC <= myDF-baseParams; spoolC = spoolC + 1)
	{
		execCommand = "fprintf(detailsOutputFile,\",\",NSR"+spoolC+");";
		ExecuteCommands (execCommand);
	}
	fprintf (detailsOutputFile,"}}\n");	

	return 0;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

ExecuteAFile ("_MSCOneStep.ibf");

fprintf (stdout, "\n", stateVectorDimension, " one step reachable states\n");

ExecuteAFile 	("_MFReader_.ibf");

dataPartsRead = fileCount;

/* 
data sets will be in ds_ID
nuc filters 				- nucData_ID
codon filters				- filteredData_ID
trees 						- treeStrings[ID]
   
ID ranges from 1 to dataPartsRead
*/

/* print some temp files - one per part */

if (Abs (_cmsShared) == 0)
{
	fscanf			("_CMS_Aux.ibf", "Raw", _cmsShared);
}
ExecuteCommands   (_cmsShared);
modelSpawnPrefix * _cmsShared;
modelSpawnPrefix * ("GeneticCodeExclusions = \"" + GeneticCodeExclusions + "\";\n");
modelSpawnPrefix * ("_Genetic_Code = " + _Genetic_Code + ";\n");


baseDataPath = LAST_FILE_PATH;
for (k=1; k<=dataPartsRead; k=k+1)
{
	h				 				= baseDataPath + "." + k;
	DATAFILE_TREE 					= treeStrings[k];
	IS_TREE_PRESENT_IN_DATA 		= 1;
	
	ExecuteCommands ("fprintf (h, CLEAR_FILE, filteredData_" + k+ ");");
	modelSpawnPrefix * ("DataSet 	ds_"+k+"  = ReadDataFile (\"" + h + "\");\n");
	modelSpawnPrefix * ("DataSetFilter	filteredData_"+k+" = CreateFilter(ds_"+k+",3,\"\",\"\",GeneticCodeExclusions);\n");
}

modelSpawnPrefix * 0;


fprintf						  (stdout, "\nModel string for nucleotide biases:");
fscanf						  (stdin,"String",modelDesc);

ExecuteAFile("Utility/GrabBag.bf");
icSampleSizeMultiplier = prompt_for_a_value ("Select the sample size multiplier (x the number of characters)", 1, 0.1, 2.0, 0);


sampleCount 				 = totalCharCount*icSampleSizeMultiplier;
fprintf						  (stdout, "Using ", sampleCount, " as a sample count.\n");


/*---------------------------------------------------------------------------------------------------------------------------------------------*/


global AC = 1;
global AT = 1;
global CG = 1;
global CT = 1;
global GT = 1;

catCounterAL = {stateVectorDimension,1};

MGCustomRateBiasTerms = {{"AC*","","AT*","CG*","CT*","GT*"}};	

		
paramCount	   = 0;
_nucBiasTerms  = {4,4};

_nucBiasTerms[0][0] = "";


if (modelDesc[0]==modelDesc[1])
{
	MGCustomRateBiasTerms[0] = MGCustomRateBiasTerms[1];
}

_nucBiasTerms[1][0] = MGCustomRateBiasTerms[0];
_nucBiasTerms[0][1] = MGCustomRateBiasTerms[0];
_nucBiasTerms[2][0] = MGCustomRateBiasTerms[1];
_nucBiasTerms[0][2] = MGCustomRateBiasTerms[1];

h = 0;
v = 3;

for (customLoopCounter2=2; customLoopCounter2<6; customLoopCounter2=customLoopCounter2+1)
{
	for (customLoopCounter=0; customLoopCounter<customLoopCounter2; customLoopCounter=customLoopCounter+1)
	{
		if (modelDesc[customLoopCounter]==modelDesc[customLoopCounter2])
		{
			_nucBiasTerms[h][v] = MGCustomRateBiasTerms[customLoopCounter];
			_nucBiasTerms[v][h] = MGCustomRateBiasTerms[customLoopCounter];
			break;
		}
	}
	if (customLoopCounter == customLoopCounter2)
	{
		_nucBiasTerms[h][v] = MGCustomRateBiasTerms[customLoopCounter2];
		_nucBiasTerms[v][h] = MGCustomRateBiasTerms[customLoopCounter2];
	}
	
	v = v+1;
	if (v==4)
	{
		h=h+1;
		v=h+1;
	}
}

_nucRateMatrix = {4,4};

modelSpawnPrefix * ("_nucBiasTerms = {4,4};_nucBiasTerms[0][0] = \"\";\n");

for (customLoopCounter = 0; customLoopCounter < 4; customLoopCounter = customLoopCounter + 1)
{
	for (customLoopCounter2 = customLoopCounter+1; customLoopCounter2 < 4; customLoopCounter2 = customLoopCounter2 + 1)
	{
		ExecuteCommands ("_nucRateMatrix[" + customLoopCounter  + "][" + customLoopCounter2 + "]:=" + _nucBiasTerms[customLoopCounter][customLoopCounter2] + "t;");
		ExecuteCommands ("_nucRateMatrix[" + customLoopCounter2 + "][" + customLoopCounter  + "]:=" + _nucBiasTerms[customLoopCounter2][customLoopCounter] + "t;");
		modelSpawnPrefix * ("_nucBiasTerms[" + customLoopCounter + "][" + customLoopCounter2 + "] = \"" + _nucBiasTerms[customLoopCounter][customLoopCounter2] + "\";");
		modelSpawnPrefix * ("_nucBiasTerms[" + customLoopCounter2 + "][" + customLoopCounter + "] = \"" + _nucBiasTerms[customLoopCounter2][customLoopCounter] + "\";");
	}
}


fprintf (stdout, "\nFitting a nucleotide model to approximate branch lengths...\n");
Model 	nucModel = (_nucRateMatrix,overallFrequencies,1);

populateTrees ("nuc_tree", dataPartsRead);
ExecuteCommands(constructLF ("lf","nucData","nuc_tree",dataPartsRead));

Optimize (nuc_res,lf);

fprintf (stdout, "\n", lf, "\n");

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

SetDialogPrompt ("Save the best model to:");
fprintf (PROMPT_FOR_FILE,CLEAR_FILE);
modelFile = LAST_FILE_PATH;

detailsOutputFile = modelFile+".details";
fprintf (detailsOutputFile,CLEAR_FILE,KEEP_OPEN,orderString,"\n");

aaRateMultipliers = {21,21};


icSampleSizeMultiplier = 1.0;

ChoiceList (whichIC, "Which IC to use for scoring?",1,SKIP_NONE,
			"c-AIC","Small sample AIC.",
			"BIC", "BIC.");

if (whichIC<0)
{
	return 0;
}


ChoiceList (matingChoice, "Mating Weights",1,SKIP_NONE,
			"Akaike Weights","Individuals are chosen to reproduce with based on their Akaike weights.",
			"Equiprobable", "All individuals are equally likely to mate.",
			"Rank","Mating probabilities are proportional to the c-AIC rank of the individual",
			"Random","One of the three mating schemes in picked at random for every generation");

if (matingChoice<0)
{
	return 0;
}

ChoiceList (seedChoice, "Seeds",1,SKIP_NONE,
			"Random","The entire starting population is seeded at random.",
			"Load", "Load some predefined partitions",
			"Completely Random","The entire starting population is seeded at random, skipping even hardwired presets."
);

if (seedChoice<0)
{
	return 0;
}

if (seedChoice == 2)
{
	predef = {};
}

ChoiceList (reprFlag, "Mating Representation",1,SKIP_NONE,
			"General","Use n-ary representation (n=number of rate classes).",
			"Binary", "Use binary representation");

if (reprFlag<0)
{
	return 0;
}

if (seedChoice == 1)
{
	SetDialogPrompt ("File with partition specifications:");
	END_OF_FILE = 0;
	k = Abs(predef);
	fscanf (PROMPT_FOR_FILE,"NMatrix",seeds);
	predef[Abs(predef)] = seeds;
	while (!END_OF_FILE)
	{
		fscanf (LAST_FILE_PATH,"NMatrix",seeds);
		if (Abs(seeds))
		{
			predef[Abs(predef)] = seeds;
		}
	}
	fprintf (stdout, "Loaded ", Abs(predef) - k, " model definitions\n");
}




observedFreq					= positionFrequencies;
ConstructParametricFrequencies	 (positionFrequencies);
PopulateModelMatrix				 ("MG94MLEFREQS",paramFreqs,1);
BuildParametricCodonFrequencies	  (paramFreqs,"paramCodonFreqs");
 
Model MG94MLEFREQS_MODEL	= ( MG94MLEFREQS, paramCodonFreqs, 0);

populateTrees			("codon_tree", dataPartsRead);
 
for (k=1; k<=dataPartsRead; k=k+1)
{
	ExecuteCommands("ReplicateConstraint (\"this1.?.synRate:=3*this2.?.t\", codon_tree_"+k+", nuc_tree_"+k+");ClearConstraints(codon_tree_"+k+");");
}

ExecuteCommands(constructLF ("lf","filteredData","codon_tree",dataPartsRead));

currentPopulation  = {};
USE_LAST_RESULTS=1;
if (verboseFlag)
{
	VERBOSITY_LEVEL 				= 1;
}
Optimize (res,lf);

fprintf  (stdout, lf);

for (k = 0; k < 4; k = k+1)
{
	for (p = 0; p < 3; p=p+1)
	{
		ExecuteCommands ("positionFrequencies[k][p] = " + paramFreqs[k][p]);
	}
}
	
vectorOfFrequencies                     = BuildCodonFrequencies (positionFrequencies);
mgPlainFactor							= FindBranchLengthExpression(0,"MG94MLEFREQS");
fprintf									(stdout, positionFrequencies, "\n", mgPlainFactor,"\n");

/*
MULTIPLY_BY_FREQS                       = PopulateModelMatrix ("MG94plain", positionFrequencies,0);
*/

modelFreqSpec							= "observedFreq = " + positionFrequencies + ";\nvectorOfFrequencies = " + vectorOfFrequencies + ";\n";
modelSpawnPrefix						* ("global AC = 1; global AT = 1; global CG = 1; global CT = 1; global GT = 1;\n");
modelSpawnPrefix						* 0;
modelSpawnSuffix						* ("MULTIPLY_BY_FREQS = PopulateModelMatrix (\"MG94custom\", observedFreq,0);\n");
modelSpawnSuffix						* ("Model MG94customModel = (MG94custom,vectorOfFrequencies,0);\n");

for (k=1; k<=dataPartsRead; k=k+1)
{
	modelSpawnSuffix *	("Tree codon_tree_"+k+" = " + treeStrings[k] + ";\n");
	modelSpawnSuffix *	("replicateBranchLengths ("+k+");");
}

modelSpawnSuffix *	("\nExecuteCommands(constructLF(\"lf\",\"filteredData\",\"codon_tree\","+dataPartsRead+"));\nOptimize (lf_MLES,lf);\nreturn makeReturnValue(0);\n");
modelSpawnSuffix *  0;


VERBOSITY_LEVEL = 0;
USE_LAST_RESULTS=0;


baseParams 		   = res[1][1];
fprintf 			(stdout, "\nBase model has ", baseParams, " parameters\n", lf);
crapAIC			   = returnIC (res[1][0], baseParams, sampleCount);

sortedScores	   = {populationSize,2};
sortedScores[0][0] = -crapAIC;
sortedScores[0][1] = 0;
currentPopulation [0] = {stateVectorDimension,1};

ExportAMatrix (detailsOutputFile,StringToMatrix(currentPopulation[0]),0,0);
fprintf		  (detailsOutputFile,"\n{{",Format(res[1][0],25,10),",",baseParams,",",Format(crapAIC,25,10));
ExecuteCommands ("fprintf(detailsOutputFile,\",\",NSR0);");
fprintf (detailsOutputFile,"}}\n");


AC:=AC__;
AT:=AT__;
CG:=CG__;
CT:=CT__;
GT:=GT__;




/* INCLUDE THE CHC CORE */

if (reprFlag)
{
	ibfPath = HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "GA_CHC_Binary.ibf";
}
else
{
	ibfPath = HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "GA_CHC.ibf";
}

byBPImprovement		  = {};
byBPImprovement[0]    = crapAIC;
bestIndividualOverall = sortedScores[0][0];
currentBEST_IC		  = crapAIC;


for (currentBPC = startWithRateClasses; currentBPC < maxRateClasses; currentBPC = currentBPC + 1)
{
	rateClassesCount 		= currentBPC;
	fprintf (stdout, "\n\nStarting GA with ", rateClassesCount, " rate classes\n");
	addOnLine = " with " + rateClassesCount + " rate classes.";
	
	/*
	modelComplexityPenalty = Log(StirlingNumberOf2ndKind (stateVectorDimension, currentBPC));
	*/
	
	fprintf (stdout, "[INFO: MODEL COMPLEXITY PENALTY TERM SET TO ", modelComplexityPenalty, "]\n");
	
	if (currentBPC > startWithRateClasses)
	{
		children = {};
		for (individual=populationSize $ 2; individual<populationSize-1; individual=individual+1)
		{
			cString    = currentPopulation[individual];
			toggleRate = Random (0,currentBPC-1.00001)$1;
			for (bitP = 0; bitP < stateVectorDimension; bitP = bitP + 1)
			{
				if (cString [bitP] == toggleRate)
				{
					if (Random(0,1)<0.5)
					{
						cString [bitP] = currentBPC-1;
					}
				}
			}
			cString = MakeStringCanonical(cString,currentBPC);
			aaRateMultipliers = StringToMatrix(cString);
			sortedScores[individual][1] = individual;
			currentPopulation[individual] = cString;
			RunASample (compressedString[rateClassesCount-1],individual);
		}
	}	
	
	mutationProb = 0.15;
	ExecuteCommands ("#include \"" + ibfPath + "\";");
	
	kf						 = -sortedScores[populationSize-1][0];
	
	if (currentBEST_IC > kf)
	{
		currentBEST_IC = kf;
	}
	else
	{
		break;
	}
}


/* try a consensus matrix */

ExportAMatrix (modelFile,StringToMatrix(currentPopulation[populationSize-1]),1,1);
fpath = modelFile+".bestAIC";

fprintf (fpath, CLEAR_FILE, -sortedScores[populationSize-1][0]);
fprintf (detailsOutputFile,CLOSE_FILE);


