RequireVersion("0.9920061128");

if (MPI_NODE_COUNT <= 1)
{
	fprintf (stdout, "[ERROR] This analysis requires an MPI environment to run\n");
	return 0;
}

RequireVersion			(".9920060901");
partCount				= 2;

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
shortestAllowedSegment  = 0;

stoppingCriterion		= 50;
sampleCount				= 0;
familyControlSize		= produceOffspring$6;

verboseFlag				= 0;
rateClassesCount		= 2;
MESSAGE_LOGGING			= 0;
cAICPenaltyPerSite		= 100;
adjustAICScores			= 1;
maxTimeAllowed			= 3600*1000;
startAtBreakpoint		= 1;

/* ________________________________________________________________________________________________*/

global AC 				= 1;
global AT 				= 1;
global CG 				= 1;
global CT 				= 1;
global GT 				= 1;



/* ________________________________________________________________________________________________*/

MasterList						= {};
REPLACE_TREE_STRUCTURE  		= 1;
bppMap							= {};
SHORT_MPI_RETURN				= 1;
totalBitSize					= 0;
LIKELIHOOD_FUNCTION_OUTPUT 		= 0;
FILE_SEPARATOR			   		= "__FILE_SEPARATOR__";

ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Utility"+DIRECTORY_SEPARATOR+"NJ.bf");


/* ________________________________________________________________________________________________*/

function ComputeDistanceFormula (s1,s2)
{
	GetDataInfo (siteDifferenceCount, filteredData, s1, s2, DIST);
	
	totalSitesCompared = Transpose(summingVector)*(siteDifferenceCount*summingVector);
	totalSitesCompared = totalSitesCompared[0];
	
	if (_useK2P)
	{
		_dTransitionCounts 	 =    siteDifferenceCount[0][2]+siteDifferenceCount[2][0]  /* A-G and G-A */
								 +siteDifferenceCount[1][3]+siteDifferenceCount[3][1]; /* C-T and T-C */
							
		_dTransversionCounts = (siteDifferenceCount[0][0]+siteDifferenceCount[1][1]+siteDifferenceCount[2][2]+siteDifferenceCount[3][3])+_dTransitionCounts;
		
		_dTransitionCounts	 = _dTransitionCounts/totalSitesCompared;
		_dTransversionCounts = 1-_dTransversionCounts/totalSitesCompared;
		
		_d1C = 1-2*_dTransitionCounts-_dTransversionCounts;
		_d2C = 1-2*_dTransversionCounts;
		
		if (_d1C>0 && _d2C>0)
		{
			return -(0.5*Log(_d1C)+.25*Log(_d2C));	
		}
	}
	else
	{
		_dAGCounts 	 =    siteDifferenceCount[0][2]+siteDifferenceCount[2][0]  /* A-G and G-A */;
		_dCTCounts	 = 	  siteDifferenceCount[1][3]+siteDifferenceCount[3][1]; /* C-T and T-C */
							
		_dTransversionCounts = (siteDifferenceCount[0][0]+siteDifferenceCount[1][1]+siteDifferenceCount[2][2]+
								siteDifferenceCount[3][3])+_dAGCounts+_dCTCounts;
		
		_dAGCounts	 = _dAGCounts/totalSitesCompared;
		_dCTCounts	 = _dCTCounts/totalSitesCompared;
		
		_dTransversionCounts = 1-_dTransversionCounts/totalSitesCompared;
		
		_d1C = 1-_dAGCounts/_d_TN_K1-0.5*_dTransversionCounts/_d_fR;
		_d2C = 1-_dCTCounts/_d_TN_K2-0.5*_dTransversionCounts/_d_fY;
		_d3C = 1-0.5*_dTransversionCounts/_d_fY/_d_fR;
		
		if ((_d1C>0)&&(_d2C>0)&&(_d3C>0))
		{
			return -_d_TN_K1*Log(_d1C)-_d_TN_K2*Log(_d2C)-_d_TN_K3*Log(_d3C);
		}
	}
	
	return 1000;
}


/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function StringToMatrix (zz)
{
	return zz;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function ExportAMatrix (fileName, rateMatrix,dummy)
{
	if (Abs(rateMatrix))
	{
		sortedBP 	= ConvertToPart (rateMatrix);
	}
	else
	{
		v = 0;
	}
	theAVL	= {};
	
	theAVL ["BP"]    = {};
	theAVL ["Trees"] = {};
	
	if (dataType)
	{
		bpF  = 0;
		bpF2 = 0;
	}	
	else
	{
		bpF 	= -1;
		bpF2	= -1;
	}
	
	USE_DISTANCES = 0;

	lfDef = "";
	lfDef * 128;
	lfDef  * "LikelihoodFunction multiPart  = (";


	partSpecs = {};
	
	for (h=0; h<v; h=h+1)
	{
		bpF2 = sortedBP[h];
		bpF2 = bppMap[bpF2];
		if (dataType)
		{
			filterString = ""+(3*bpF)+"-"+(3*bpF2-1);
		}
		else
		{
			filterString = ""+(bpF+1)+"-"+bpF2;		
			(theAVL ["BP"])[h]   = ""+(bpF+2)+"-"+(bpF2+1);

		}
		DataSetFilter filteredData = CreateFilter(ds,1,filterString);
		InferTreeTopology (0);
		treeString=TreeMatrix2TreeString(0);

		ExecuteCommands ("DataSetFilter part_" + h + " = CreateFilter (ds, 1, \"" + filterString + "\");");
		ExecuteCommands ("Tree tree_" 		   + h + " = " + treeString + ";");
				
		if (h)
		{
			lfDef * ",";
		}
		lfDef  * ("part_"+h+",tree_"+h);
		partSpecs [h] = filterString;
		bpF = bpF2;
	}
	
	if ((bpF2<ds.sites && (dataType == 0)) || (3*bpF2<ds.sites && dataType))
	{
		if (dataType)
		{
			filterString = ""+(3*bpF2)+"-"+(ds.sites-1);
		}
		else
		{
			filterString = ""+(bpF2+1)+"-"+(ds.sites-1);
			(theAVL ["BP"])[h]   = ""+(bpF+2)+"-"+ds.sites;
		}
		DataSetFilter filteredData = CreateFilter(ds,1,filterString);
		partSpecs [h] = filterString;
		InferTreeTopology 		  (0);
		treeString				   = TreeMatrix2TreeString(0);
		ExecuteCommands 			("DataSetFilter part_" + h + " = CreateFilter (ds, 1, \"" + filterString + "\");");
		ExecuteCommands 			("Tree tree_" + h + " = " + treeString + ";");
		if (h)
		{
			lfDef * ",";
		}
		lfDef  * ("part_"+h+",tree_"+h);
	}

	lfDef  * ");";
	lfDef  * 0;
	ExecuteCommands (lfDef);
	Optimize (res, multiPart);

	ConstraintString = "";
	ConstraintString * 8192;
	ConstraintString * "\n";
	
	for (h=0; h<Abs (partSpecs); h=h+1)
	{
		ConstraintString * ("\n" + partSpecs [h] + "\n");
		ExecuteCommands    ("filterString = Format (tree_" + h + ",0,1);");
		ConstraintString * filterString;
		(theAVL ["Trees"]) [h]    = filterString;
	}

	ConstraintString * 0;
	fprintf (fileName,ConstraintString);

	return theAVL;
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

function adjustAICScore (theLL,bpMatrix)
{
	daAICScore = 2*(baseParams*(baseSites/(baseSites-baseParams-1)) - theLL) ;
	lastBpos   = 0;
	
	for (_pid = 0; _pid < Rows(bpMatrix); _pid = _pid+1)
	{
		thisSpan = bppMap[bpMatrix[_pid]] - lastBpos+1;
		lastBpos = bppMap[bpMatrix[_pid]];
		if (thisSpan > perPartParameterCount + 1)
		{
			daAICScore = daAICScore + 2*(perPartParameterCount*(thisSpan/(thisSpan-perPartParameterCount-1)));
		}
		else
		{
			daAICScore = daAICScore + 2*perPartParameterCount*cAICPenaltyPerSite;
		}
	}
	
	thisSpan = baseSites-lastBpos;
	if (thisSpan > perPartParameterCount + 1)
	{
		daAICScore = daAICScore + 2*(perPartParameterCount*(thisSpan/(thisSpan-perPartParameterCount-1)));
	}
	else
	{
		daAICScore = daAICScore + 2*perPartParameterCount*cAICPenaltyPerSite;
	}
	
	return -daAICScore;
}


/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function ReceiveJobs (sendOrNot, ji)
{
	if (MPI_NODE_COUNT>1)
	{
		MPIReceive (-1, fromNode, result_String);
		mji = MPINodeState[fromNode-1][1];
		
		if (sendOrNot)
		{
			MPISend 	(fromNode,mpiStringToSend);
			MPINodeState[fromNode-1][1] = ji;			
		}
		else
		{
			MPINodeState[fromNode-1][0] = 0;
			MPINodeState[fromNode-1][1] = -1;		
		}
		ExecuteCommands (result_String);
		myDF  = lf_MLES[1][1]+baseParams;
		myAIC = 2*(lf_MLES[1][0]-myDF*(baseSites/(baseSites-myDF-1)));
		myLL  = lf_MLES[1][0];
		ji 	  = mji;
	}
	else
	{
		myDF  = res[1][1]+baseParams;
		myAIC = 2*(res[1][0]-myDF*(baseSites/(baseSites-myDF-1)));
	}
	
	sortedBP = {{-1}};
	
	if (resultProcessingContext==0)
	{
		sortedScores[ji][0] = myAIC;
		if (ji>=0)
		{
			jobPrint = ConvertToPartString (currentPopulation[ji]);
			sortedBP = ConvertToPart 	   (currentPopulation[ji]);
			if (adjustAICScores)
			{
				myAIC	 = adjustAICScore (myLL, sortedBP);
			}
			v 		 = Rows (sortedBP);
			sortedScores[ji][0] = myAIC;
		}
		if (verboseFlag)
		{
			fprintf (stdout, "Individual ",ji," AIC-c = ",-myAIC," ");
		}
	}
	else
	{
		intermediateProbs[ji][0] = myAIC;	
		if (ji>=0)
		{
			jobPrint = ConvertToPartString (children[ji-populationSize]);
			sortedBP = ConvertToPart 	   (children[ji-populationSize]);
			if (adjustAICScores)
			{
				myAIC	 = adjustAICScore (myLL, sortedBP);
			}
			v = Rows (sortedBP);
			intermediateProbs[ji][0] = myAIC;	
		}
		if (verboseFlag)
		{
			fprintf (stdout, "Offspring ",ji," AIC-c = ",-myAIC," ");
		}
	}
	
	if (sortedBP[0]>=0)
	{
		if (dataType)
		{
			bpF = 0;
		}
		else
		{
			bpF = -1;
		}
		filterString = "";
		for (h=0; h<v; h=h+1)
		{
			bpF2 = sortedBP[h];
			bpF2 = bppMap[bpF2];
			if (dataType)
			{
				filterString = filterString+" "+(3*bpF)+"-"+(3*bpF2-1);
			}
			else
			{
				filterString = filterString+" "+(bpF+1)+"-"+bpF2;			
			}
			bpF = bpF2;
		}
		
		if ((bpF2<ds.sites && (dataType == 0)) || (bpF2<3*ds.sites && dataType))
		{
			if (dataType)
			{
				filterString = filterString+" "+(3*bpF2)+"-"+(ds.sites-1);
			}
			else
			{
				filterString = filterString+" "+(bpF2+1)+"-"+(ds.sites-1);			
			}
		}

		if (verboseFlag)
		{
			fprintf (stdout, " ", filterString, "\n");
		}

		MasterList [jobPrint] = myAIC;
	}
	return fromNode-1;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

compressedString = {{1,1}};

function MakeStringCanonical (someString, dummy)
{
	return someString;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function ConvertToPartString (pString)
{
	sortedBP = ConvertToPart (pString);
	if (dataType)
	{
		bpF = 0;	
	}
	else
	{
		bpF = -1;
	}
	
	minPartLength 	  = 1e100;
	
	_ConstraintString = "";
	_ConstraintString * 256;
	

	for (h=0; h<v; h=h+1)
	{
		bpF2 = bppMap[sortedBP[h]];
	
		if (h)
		{
			_ConstraintString * ",";
		}
		if (dataType)
		{
			_ConstraintString * (""+(3*bpF)+"-"+(3*bpF2-1));		
			curSegLength = 3*(bpF2-bpF);
		}		
		else
		{
			_ConstraintString * (""+(bpF+1)+"-"+bpF2);		
			curSegLength = bpF2-bpF;
		}
		bpF = bpF2;
		
		if (curSegLength < minPartLength && curSegLength>0)
		{
			minPartLength = curSegLength;
		}
	}
	
	if (bpF2<ds.sites && dataType == 0 || bpF2<3*ds.sites && dataType)
	{
		if (dataType)
		{
			_ConstraintString * (","+(3*bpF2)+"-"+(ds.sites-1));
			curSegLength = ds.sites-3*bpF2;
		}
		else
		{
			_ConstraintString * (","+(bpF2+1)+"-"+(ds.sites-1));		
			curSegLength = ds.sites-bpF2;
		}
		if (curSegLength < minPartLength && curSegLength>0)
		{
			minPartLength = curSegLength;
		}
	}

	_ConstraintString * 0;
	return _ConstraintString;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function ConvertToPart (pString)
{
	partitionHits    = 		{};
	h 				 = 		0; 
	v 				 = 		bppSize;
	
	for (mpiNode=0; mpiNode<partCount; mpiNode=mpiNode+1)
	{
		aBP    = 0;
		bpF	   = 2^(bppSize-1);
		
		for (; h<v; h=h+1)
		{
			aBP = aBP + bpF*(0+pString[h]);
			bpF = bpF/2;
		}
		
		if (aBP>=Abs(bppMap)-1)
		{
			aBP = aBP - 2^(bppSize-1);
		}
		
		v = v + bppSize;
		partitionHits[aBP] = 1;
	}

	meKeys	 = Rows(partitionHits);
	v 		 = Columns(meKeys);
	sortedBP = {v,1};
	
	for (h=0; h<v; h=h+1)
	{
		sortedBP [h] = 0+meKeys[h];
	}
	
	sortedBP = sortedBP % 0;
	return 	   sortedBP;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function RunASample (dummy,jobIndex)
{	
	myAIC	 = MasterList[ConvertToPartString (cString)];

	if (myAIC<0)
	{		
		if (resultProcessingContext==0)
		{
			sortedScores[jobIndex][0] = myAIC;
			if (verboseFlag)
			{
				fprintf (stdout, "Individual ",jobIndex," AIC-c = ",-myAIC, "\n");
			}
		}
		else
		{
			intermediateProbs[jobIndex][0] = myAIC;	
			if (verboseFlag)
			{
				fprintf (stdout, "Offspring ",jobIndex," AIC-c = ",-myAIC,"\n");
			}
		}	
		return 0;
	}

	if (dataType)
	{
		bpF = 0;	
	}
	else
	{
		bpF = -1;
	}
	
	mpiStringToSend 			 = "";
	mpiStringToSend				* 8192;
	
	ConstraintString			 = "";
	LikelihoodFunctionString 	 = "";
	ConstraintString 			* 8192;
	LikelihoodFunctionString 	* 256;
	LikelihoodFunctionString	* "LikelihoodFunction lf=(";
	

	for (h=0; h<v; h=h+1)
	{
		bpF2 = bppMap[sortedBP[h]];
		if (dataType)
		{
			filterString = ""+(3*bpF)+"-"+(3*bpF2-1);
		}
		else
		{
			filterString = ""+(bpF+1)+"-"+bpF2;		
		}
		ConstraintString * ("DataSetFilter filteredData = CreateFilter(ds,1,\""+filterString+"\");");
		if (dataType)
		{
			ConstraintString * ("DataSetFilter filteredData"+h+" = CreateFilter (ds,3,\""+filterString+"\",\"\",GeneticCodeExclusions);");		
		}
		else
		{
			ConstraintString * ("DataSetFilter filteredData"+h+" = CreateFilter (ds,1,\""+filterString+"\");");
		}
		ConstraintString * ("InferTreeTopology (0);treeString=TreeMatrix2TreeString(0);Tree givenTree"+h+" = treeString;");
		if (h)
		{
			LikelihoodFunctionString * ",";
		}
		LikelihoodFunctionString * ("filteredData" + h + ",givenTree"+h);
		bpF = bpF2;
	}
	
	if ((bpF2<ds.sites && (dataType == 0)) || (bpF2<3*ds.sites && dataType))
	{
		if (dataType)
		{
			filterString = ""+(3*bpF2)+"-"+(ds.sites-1);
		}
		else
		{
			filterString = ""+(bpF2+1)+"-"+(ds.sites-1);		
		}

		ConstraintString * ("DataSetFilter filteredData = CreateFilter(ds,1,\""+filterString+"\");");
		if (dataType)
		{
			ConstraintString * ("DataSetFilter filteredData"+h+" = CreateFilter (ds,3,\""+filterString+"\",\"\",GeneticCodeExclusions);");		
		}
		else
		{
			ConstraintString * ("DataSetFilter filteredData"+h+" = CreateFilter (ds,1,\""+filterString+"\");");
		}
		ConstraintString * ("InferTreeTopology (0);treeString=TreeMatrix2TreeString(0);Tree givenTree"+h+" = treeString;");

		LikelihoodFunctionString * (",filteredData" + h + ",givenTree"+h);
	}


	LikelihoodFunctionString * (");");
	ConstraintString 		 * 0;
	LikelihoodFunctionString * 0;	
	
	mpiStringToSend * _mpiPrefixString;
	mpiStringToSend * ConstraintString;
	mpiStringToSend * LikelihoodFunctionString;
	
	
	if ((MPI_NODE_COUNT>1) && (jobIndex>=0))
	{
		mpiStringToSend * ("Optimize(res,lf);return \"lf_MLES=\"+res+\";\";");
		mpiStringToSend * 0;
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
			MPISend (mpiNode+1,mpiStringToSend);
			MPINodeState[mpiNode][0] = 1;
			MPINodeState[mpiNode][1] = jobIndex;
		}
	}
	else
	{
		mpiStringToSend * 0;
		ExecuteCommands (mpiStringToSend);
		Optimize (res,lf);
		if (jobIndex>=0)
		{
			mpiNode = ReceiveJobs (1, jobIndex);
		}
		else
		{
			myAIC = 2*(res[1][0]-res[1][1]-baseParams);
		}
	}
	return 0;	
}




/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function SpawnRandomString (clsCnt)
{
	rModel = {totalBitSize,1};
	for (h=0; h<totalBitSize; h=h+1)
	{
		rModel[h] = Random(0,2)$1;
	}
	return rModel;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function IsChildViable (putativeChild)
{
	sampleString = 	ConvertToPartString (putativeChild);
	myAIC		 = MasterList[sampleString];
	testChild 	 = putativeChild;
	mutPassCount = 1;
	
	for (_lc = 0; _lc < populationSize && myAIC > (-0.1); _lc = _lc + 1)
	{
		if (Rows (currentPopulation[_lc]))
		{
			myAIC = -(sampleString == ConvertToPartString (currentPopulation[_lc]));
		}
	} 
	for (_lc = 0; _lc < Abs(children) && myAIC > (-0.1); _lc = _lc + 1)
	{
		if (Rows (children[_lc]))
		{
			myAIC = -(sampleString == ConvertToPartString (children[_lc]));
		}
	} 
	
	while ((myAIC<(-0.1) || minPartLength<shortestAllowedSegment)&& mutPassCount < 25)
	{
		if (verboseFlag > 1)
		{
			fprintf (stdout,"Adjusting the child to avoid a duplicate. Min(fragment) = ",minPartLength,".  Pass ", mutPassCount, "\n");
		}
		
		mutPassCount 	= mutPassCount + 1;
		sampleString 	= Min(Random(0,stateVectorDimension)$1,stateVectorDimension-1);
		myAIC 			= testChild[sampleString];
		newValue 		= Random (0,rateClassesCount-0.0000001)$1;
		
		while (newValue == myAIC)
		{
			newValue 	= Random (0,rateClassesCount-0.0000001)$1;
		}
		
		testChild [sampleString]	 = newValue;
		sampleString 				 = ConvertToPartString (testChild);
		myAIC 						 = MasterList		   [sampleString];
		for (_lc = 0; _lc < populationSize && myAIC > (-0.1); _lc = _lc + 1)
		{
			if (Rows (currentPopulation[_lc]))
			{
				myAIC = -(sampleString == ConvertToPartString (currentPopulation[_lc]));
			}
		} 
		for (_lc = 0; _lc < Abs(children) && myAIC > (-0.1); _lc = _lc + 1)
		{
			if (Rows (children[_lc]))
			{
				myAIC = -(sampleString == ConvertToPartString (children[_lc]));
			}
		} 
	}
	return testChild;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function UpdateBL (dummy)
{
	return 0;
}


/*---------------------------------------------------------------------------------------------------------------------------------------------*/


dataType = 0;

fprintf	    (stdout, "Initialized GARD on ", MPI_NODE_COUNT, " MPI nodes.\nPopulation size is ", populationSize, " models\n");
SetDialogPrompt 	("Nucleotide file to screen:");
DataSet 	ds    = ReadDataFile (PROMPT_FOR_FILE);

baseFilePath 	  = LAST_FILE_PATH;

done = 0;

while (!done)
{
	fprintf (stdout,"\nPlease enter a 6 character model designation (e.g:010010 defines HKY85):");
	fscanf  (stdin,"String", modelDesc);
	if (Abs(modelDesc)==6)
	{	
		done = 1;
	}
}

ChoiceList (rvChoice, "Rate variation options", 1, SKIP_NONE,
					  "None", 					"Homogeneous rates across sites (fastest)",
					  "General Discrete", 		"General discrete distribution on N-bins",
					  "Beta-Gamma", 			"Adaptively discretized gamma N-bins. Try is GDD doesn't converge well, or if you suspect a lot of rate classes and do not want to overparameterize the model");
					  
if (rvChoice<0)
{
	return 0;
}

if (rvChoice)
{
	rateClasses = 0;
	while (rateClasses < 2 || rateClasses > 32)
	{
		fprintf (stdout, "How many distribution bins [2-32]?:");
		fscanf  (stdin,  "Number", rateClasses);
		rateClasses = rateClasses $ 1;
	}
	fprintf (stdout, "\nUsing ", rateClasses, " distribution bins\n");
}

DEFAULT_FILE_SAVE_NAME 	 	= baseFilePath + ".html";

SetDialogPrompt			   ("Save results to:");
fprintf 				   (PROMPT_FOR_FILE, CLEAR_FILE, KEEP_OPEN);
outputFilePath 	  		  = LAST_FILE_PATH;



DataSetFilter filteredData  = CreateFilter (ds,1);

InferTreeTopology (0);
treeString = TreeMatrix2TreeString(0);

/* find "informative sites" */

if (dataType==0)
{
	for (h=0; h<filteredData.sites; h=h+1)
	{
		filterString = Format(h,20,0);
		DataSetFilter siteFilter = CreateFilter (filteredData,1,filterString);

		HarvestFrequencies (f1, siteFilter, 1, 1, 0);
		m1 = 0;
		for (mpiNode=0; (mpiNode < 4) && (m1<=1) ; mpiNode=mpiNode+1)
		{
			if (f1[mpiNode]>0)
			{
				m1=m1+1;
			}
		}	
		if (m1>1)
		{
			bppMap[Abs(bppMap)] = h;
		}
	}
}


bppSize = (Log(Abs(bppMap))/Log(2)+1)$1;
fprintf (stdout, "There are ",Abs(bppMap)," potential breakpoints. Bit size of the sample is ", bppSize,".\n");


partCount = 2;
h = Abs(bppMap);

if (h <= partCount)
{
	fprintf (outputFilePath,   "ERROR: There are too few potential break points to support ", partCount-1, " recombination events.\n");
	return 0;
}


maxBPP 	   = h-1;
ModelTitle = ""+modelDesc[0];
			
rateBiasTerms = {{"AC","1","AT","CG","CT","GT"}};
paramCount	  = 0;

modelConstraintString = "";

for (customLoopCounter2=1; customLoopCounter2<6; customLoopCounter2=customLoopCounter2+1)
{
	for (customLoopCounter=0; customLoopCounter<customLoopCounter2; customLoopCounter=customLoopCounter+1)
	{
		if (modelDesc[customLoopCounter2]==modelDesc[customLoopCounter])
		{
			ModelTitle  = ModelTitle+modelDesc[customLoopCounter2];	
			if (rateBiasTerms[customLoopCounter2] == "1")
			{
				modelConstraintString = modelConstraintString + rateBiasTerms[customLoopCounter]+":="+rateBiasTerms[customLoopCounter2]+";";
			}
			else
			{
				modelConstraintString = modelConstraintString + rateBiasTerms[customLoopCounter2]+":="+rateBiasTerms[customLoopCounter]+";";			
			}
			break;
		}
	}
	if (customLoopCounter==customLoopCounter2)
	{
		ModelTitle = ModelTitle+modelDesc[customLoopCounter2];	
	}
}	

if (Abs(modelConstraintString))
{
	ExecuteCommands (modelConstraintString);
}


if (rvChoice)
{
	rateClasses = rateClasses $ 1;
	if (rateClasses < 2)
	{
		rateClasses = 2;
	}
	else
	{
		if (rateClasses > 8)
		{
			rateClasses = 8;
		}
	}	
	
	if (ds.species < 5)
	{
		rateClasses = Min (3,rateClasses);
	}
	else
	{
		if (ds.species < 10)
		{
			rateClasses = Min (4,rateClasses);	
		}
		else
		{
			if (ds.species < 30)
			{
				rateClasses = Min (5,rateClasses);	
			}
		}
	}
	
	if (rvChoice == 1)
	{
		gdDefString = "";
		gdDefString * 1024;
		for (mi=1; mi<rateClasses; mi=mi+1)
		{
			gdDefString*("global PS_"+mi+" = 1/"+((rateClasses+1)-mi)+";\nPS_"+mi+":<1;\n");
		}
		
		gdDefString*("\n\nglobal RS_1 = .3;\nRS_1:<1;RS_1:>0.000000001;\n");

		for (mi=3; mi<=rateClasses; mi=mi+1)
		{
			gdDefString*("global RS_"+mi+" = 1.5;"+"\nRS_"+mi+":>1;RS_"+mi+":<100000;\n");
		} 

		rateStrMx    = {rateClasses,1};
		rateStrMx[0] = "RS_1";
		rateStrMx[1] = "1";

		for (mi=3; mi<=rateClasses; mi=mi+1)
		{
			rateStrMx[mi-1] = rateStrMx[mi-2]+"*RS_"+mi;
		} 	

		freqStrMx    = {rateClasses,1};
		freqStrMx[0] = "PS_1";

		for (mi=1; mi<rateClasses-1; mi=mi+1)
		{
			freqStrMx[mi] = "";
			for (mi2=1;mi2<=mi;mi2=mi2+1)
			{
				freqStrMx[mi] = freqStrMx[mi]+"(1-PS_"+mi2+")";		
			}
			freqStrMx[mi] = freqStrMx[mi]+"PS_"+(mi+1);	
		}	

		freqStrMx[mi] = "";
		for (mi2=1;mi2<mi;mi2=mi2+1)
		{
			freqStrMx[mi] = freqStrMx[mi]+"(1-PS_"+mi2+")";		
		}
		freqStrMx[mi] = freqStrMx[mi]+"(1-PS_"+mi+")";	


		gdDefString*("\n\nglobal c_scale:="+rateStrMx[0]+"*"+freqStrMx[0]);

		for (mi=1; mi<rateClasses; mi=mi+1)
		{
			gdDefString*("+"+rateStrMx[mi]+"*"+freqStrMx[mi]);
		}

		gdDefString*(";c_scale :< 1e100; \ncategFreqMatrix={{"+freqStrMx[0]);

		for (mi=1; mi<rateClasses; mi=mi+1)
		{
			gdDefString*(","+freqStrMx[mi]);
		}

		gdDefString*("}};\ncategRateMatrix={{"+rateStrMx[0]+"/c_scale");

		for (mi=1; mi<rateClasses; mi=mi+1)
		{
			gdDefString*(","+rateStrMx[mi]+"/c_scale");
		}

		gdDefString*("}};\n\ncategory c  = ("+rateClasses+", categFreqMatrix , MEAN, ,categRateMatrix, 0, 1e25);\n\n");
		gdDefString*0;
		ExecuteCommands (gdDefString);	
	}
	else
	{
		global betaP = 1;
		global betaQ = 1;
		betaP:>0.05;betaP:<85;
		betaQ:>0.05;betaQ:<85;
		category pc = (rateClasses-1, EQUAL, MEAN, 
						_x_^(betaP-1)*(1-_x_)^(betaQ-1)/Beta(betaP,betaQ), /* density */
						IBeta(_x_,betaP,betaQ), /*CDF*/
						0, 				   /*left bound*/
						1, 			   /*right bound*/
						IBeta(_x_,betaP+1,betaQ)*betaP/(betaP+betaQ)
					   );
		
		global alpha = .5;
		alpha:>0.01;alpha:<100;
		category c = (rateClasses, pc, MEAN, 
						GammaDist(_x_,alpha,alpha), 
						CGammaDist(_x_,alpha,alpha), 
						0 , 
				  	    1e25,
				  	    CGammaDist(_x_,alpha+1,alpha)
				  	 );
			
	}
	NucleotideMatrix	 = {{*,c*AC*t,c*t,c*AT*t}
							{c*AC*t,*,c*CG*t,c*CT*t}
							{c*t,c*CG*t,*,c*GT*t}
							{c*AT*t,c*CT*t,c*GT*t,*}};
}
else
{
	NucleotideMatrix	 = {{*,AC*t,t,AT*t}{AC*t,*,CG*t,CT*t}{t,CG*t,*,GT*t}{AT*t,CT*t,GT*t,*}};
}

HarvestFrequencies 		  (nucEFV, filteredData, 1, 1, 1);
Model nucModel   		= (NucleotideMatrix, nucEFV, 1);
Tree givenTree 			= treeString;
branchNames	  		    = BranchName (givenTree,-1);
LikelihoodFunction lf   = (filteredData,givenTree);


/* check parameter counts */

GetString 				(varList, lf, -1);
baseParams 		   		= Columns(varList["Global Independent"])+3;
perPartParameterCount	= Columns(varList["Local Independent"]);
baseSites		   		= filteredData.sites;

if (baseParams + (partCount+1) * perPartParameterCount >= baseSites - 1)
{
	fprintf (stdout,   "ERROR: Too few sites for c-AIC inference.\n");
	return 0;
}

if (baseParams + 2 * perPartParameterCount >= baseSites - 1)
{
	fprintf (stdout,   "ERROR: Too few sites for reliable c-AIC inference.\n");
	return 0;
}


fprintf (stdout, 		"\nFitting a baseline nucleotide model...\n");
Optimize 		   (res,lf);
baseLL			   = res[1][0];

currentPopulation  = {};
sortedScores	   = {populationSize,2};
/*baseParams 		   = res[1][2];*/

myDF 			   = baseParams+perPartParameterCount;

sortedScores[0][0] = 2*(res[1][0]-myDF*(baseSites/(baseSites-myDF-1)));
sortedScores[0][1] = 0;

fprintf (stdout, CLEAR_FILE, "Done with single partition analysis. ",
				 "Log(L) = ",lf,
				 ", c-AIC  = ",-sortedScores[0][0], "\n");


if (baseParams>3)
{

	ConstraintString = "";
	ConstraintString*256;
	for (h=0; h<baseParams-3; h=h+1)
	{
		GetString (v,lf,h);
		ConstraintString * (v+":="+v+"__;\n");
	}
	ConstraintString*0;
	ExecuteCommands (ConstraintString);
}

_mpiPrefixString = "";
_mpiPrefixString * 256;

_mpiPrefixString * ("ExecuteAFile (HYPHY_LIB_DIRECTORY+\"TemplateBatchFiles\"+DIRECTORY_SEPARATOR+\"Utility\"+DIRECTORY_SEPARATOR+\"NJ.bf\");");
_mpiPrefixString * ("DataSet 	ds    = ReadDataFile (\""+baseFilePath+"\");");
Export (modelString, USE_LAST_MODEL);
_mpiPrefixString * modelString;
_mpiPrefixString * 0;

crapAIC 		 = -sortedScores[0][0];
startTimer		 = Time (1);

if (MPI_NODE_COUNT>1)
{
	MPINodeState 			= {MPI_NODE_COUNT-1,3};
	MPINodeBreakpoints		= {};
}

currentBEST_IC 			 = crapAIC;
ibfPath 				 = "GARD_GA_CHC.ibf";
current_BPP_done_counter = 0;

lastImprovedBPC = 0;

byBPImprovement = {};
byBPSplits		= {};

DataSetFilter allData = CreateFilter (ds, 1);

byBPImprovement[0]    = crapAIC;
bestIndividualOverall = 0;

fprintf (stdout, "Starting the GA...\n");

for (currentBPC = startAtBreakpoint; currentBPC < maxBPP; currentBPC = currentBPC + 1)
{
	partCount 				= currentBPC;
	totalBitSize 			= bppSize * partCount;
	stateVectorDimension 	= totalBitSize;
	
	if (currentBPC == startAtBreakpoint)
	{
		currentPopulation [0] = {totalBitSize,1};
	}
	else
	{
		for (k=0; k<populationSize; k=k+1)
		{
			oldVector 			 = currentPopulation[k];
			newVector			 = {totalBitSize,1};
			currentPopulation[k] = newVector ["oldVector[_MATRIX_ELEMENT_ROW_%(totalBitSize-bppSize)]"];
		}
		children = {};
	}	
	
	totalModelCounter 		 = 1;
	kf 						 = 1;

	for (k=1; k <= partCount; k=k+1)
	{
		totalModelCounter = totalModelCounter * (Abs(bppMap)-k+1);
		kf 				  = kf * k;
	} 
	totalModelCounter = totalModelCounter / kf;

	ExecuteAFile			 (ibfPath);	 
	
	current_BPP_done_counter = Abs (MasterList);
	kf						 = -sortedScores[populationSize-1][0];
	
	if (currentBEST_IC > kf || currentBPC == 1)
	{
		byBPSplits		[currentBPC] = ConvertToPart 	(currentPopulation[populationSize-1]);
		byBPImprovement [currentBPC] = currentBEST_IC-kf;
		if (currentBEST_IC > kf)
		{
			lastImprovedBPC       = currentBPC;
			bestIndividualOverall = currentPopulation[populationSize-1];
		}
		currentBEST_IC = Min(kf,currentBEST_IC);
	}
	else
	{
		break;
	}
}


fprintf (outputFilePath, "<!DOCTYPE HTML PUBLIC '-//W3C//DTD HTML 4.01 Transitional//EN'><html><head><META http-equiv='Content-Style-Type' content='text/css'><meta http-equiv='Content-Type' content='text/html; charset=ISO-8859-1'><title>GARD Results</title><LINK REL=STYLESHEET TYPE='text/css' HREF='http://www.datamonkey.org/GARD/styles/styles.css'></head><body>");
fprintf (outputFilePath, "<DIV CLASS = 'RepClassSM'>Processed alignment <b>", baseFilePath,"</b> with ", ds.species, " sequences and ", ds.sites, " sites.",
					      "<p> Processor time taken: ", Time(1)-startTimer, " seconds.</DIV>");
					      
if (timeLeft < 0)
{
	fprintf (outputFilePath, "<DIV CLASS = 'ErrorTag'> This analysis was stopped before convergence, because the CPU time limit per job was reached. The reported results may, therefore be incomplete, as there may be more breakpoints, for example. ",
							 "Consider reducing the size of the alignment, either by removing sequences or selecting shorted sequence regions, or choosing a simpler model (no rate variation etc). </DIV>\n");
}

fprintf (outputFilePath, "<DIV class = 'RepClassSM'><u>Results</u><p>", reportImprovementHTML (0), "</DIV>");

						 

fprintf (outputFilePath, "<DIV CLASS='RepClassSM'><span class = 'font-variant:small-caps;'>Substitution process parameters</span><p><b>Nucleotide Model (",modelDesc,
		") Fit Results</b><p>Log likelihood : ",outputSpan (150,baseLL),"<br>",
		"c-AIC : ", outputSpan(150, crapAIC),"<br>");

if (rvChoice)
{
	if (rvChoice == 1)
	{
		fprintf (outputFilePath, "<p>Using general discrete distribution of rates across sites");
	}	
	else
	{
		fprintf (outputFilePath, "<p>Using the beta-gamma distribution of rates across sites");	
	}
	GetInformation (cI,c);
	for (k = 0; k < Columns (cI); k=k+1)
	{
		fprintf (outputFilePath, "<br> Rate : ", Format (cI[0][k],7,3), outputSpan (150, "Weight : " + Format (cI[1][k],7,3)));
	}
}
else
{
	fprintf (outputFilePath, "<p>Homogeneous rate distribution across sites");
}

fprintf (outputFilePath, "<p><u>Nucleotide substitution rate matrix</u><p>",
		"<TABLE><TR CLASS = 'TRReportT'><TH>&nbsp;</TH><TH>A</TH><TH>C</TH><TH>G</TH><TH>T</TH></TR>",
		"<TR CLASS = 'TRReport1'><TH>A</TH><TD>&#42;</TD><TD>",AC,"</TD><TD>1</TD><TD>",AT,"</TD></TR>",		
		"<TR CLASS = 'TRReport2'><TH>C</TH><TD>&#45;</TD><TD>&#42</TD><TD>",CG,"</TD><TD>",CT,"</TD></TR>",		
		"<TR CLASS = 'TRReport1'><TH>G</TH><TD>&#45;</TD><TD>&#45</TD><TD>&#42</TD><TD>",GT,"</TD></TR>",		
		"<TR CLASS = 'TRReport2'><TH>T</TH><TD>&#45;</TD><TD>&#45</TD><TD>&#45</TD><TD>&#42</TD></TR></TABLE>",		
		"</DIV></html>");


fprintf  (outputFilePath, "</body></html>", CLOSE_FILE);
fprintf	  (stdout, "Performing the final optimization...\n");

rawOutFile = outputFilePath + "_splits";
fprintf (rawOutFile,CLEAR_FILE);



if (Rows(bestIndividualOverall))
{
	totalBitSize = Rows (bestIndividualOverall)*Columns(bestIndividualOverall);
	partCount    = totalBitSize/bppSize;
	outAVL 		 = ExportAMatrix (rawOutFile,bestIndividualOverall,0);
}
else
{
	outAVL 		= ExportAMatrix (rawOutFile,0,0);
}


rawOutFile = outputFilePath + "_ga_details";
fprintf (rawOutFile,CLEAR_FILE, KEEP_OPEN);

masterKeys = Rows(MasterList);
for (h=Rows(masterKeys)*Columns(masterKeys)-1; h>=0; h=h-1)
{
	aKey = masterKeys[h];
	fprintf (rawOutFile,(-MasterList[aKey]),"\n",aKey,"\n");
}
fprintf (rawOutFile,CLOSE_FILE);

IS_TREE_PRESENT_IN_DATA = 0;
DATA_FILE_PRINT_FORMAT	= 6;

tempFile 			  = outputFilePath + "_finalout";

fprintf 	 		 (tempFile,CLEAR_FILE,KEEP_OPEN,allData,"\nBEGIN TREES;\n");
bpL					 = outAVL ["BP"];
bpT					 = outAVL ["Trees"];
for 				 (rvChoice = 0; rvChoice < Abs(bpT); rvChoice = rvChoice + 1)
{
	fprintf (tempFile, "\tTREE part_", rvChoice+1, " = ", bpT[rvChoice], ";\n");
}
fprintf 	 		 (tempFile,"\nEND;\n\nBEGIN ASSUMPTIONS;\n");
for 				 (rvChoice = 0; rvChoice < Abs(bpL); rvChoice = rvChoice + 1)
{
	fprintf (tempFile, "\tCHARSET span_", rvChoice+1, " = ", bpL[rvChoice], ";\n");
}
fprintf				 (tempFile,"\nEND;\n",CLOSE_FILE);


/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function	reportImprovementHTML (_IC)
{
	if (lastImprovedBPC)
	{
		return "<DIV class = 'RepClassSM'><b>GARD found evidence of "+lastImprovedBPC+" breakpoints</b><p>"+spoolAICTable(0) + "</DIV>\n";
								
	}
	return "<DIV class = 'RepClassSM'>GARD found no evidence of recombination</DIV>";
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function	outputSpan (_offset, _text)
{
	return "<span style = 'position: absolute; left: " + (140+_offset) + "px'>" + _text + "</span>";
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function 	spoolAICTable (dummy)
{
	colorList 	= {{"red","black","blue","green","purple","orange"}};
	fcolorList 	= {{"white","white","white","white","white","black"}};

	htmlAICTable = "";
	htmlAICTable * 128;


	htmlAICTable * "<html><body><div style = 'width: 580px; border: black solid 1px; '>";
	htmlAICTable * "<table style = 'width: 100%;font-size: 10px;text-align:left;'><tr><th>BPs</th><th>c-AIC</th><th>&Delta; c-AIC</th><th width = '70%'>Segments</th></tr>";


	currentAIC = byBPImprovement [0];

	for (_partCount = 0; _partCount <Abs (byBPImprovement); _partCount = _partCount + 1)
	{
		if (_partCount)
		{
			currentAIC = currentAIC - byBPImprovement [_partCount];
			ci 		   = byBPImprovement [_partCount];
			bpLocs2    = byBPSplits		 [_partCount];
			pxPerSpan  = 406/allData.sites;
			sp		   = "<table style = 'padding: 0px; spacing: 0px;'><tr>";
			
			prv 	   = Rows (bpLocs2);
			bpLocs	   = {1,prv+1};
			for (k = 0; k < prv; k=k+1)
			{
				bpLocs[k]	 = bppMap[bpLocs2[k]];
			}
			
			bpLocs[prv]  = allData.sites-1;
			prv 	     = 1;
			
			for (k=0; k<Columns (bpLocs); k=k+1)
			{
				sp = sp + "<td style = 'width:"+
					 pxPerSpan*(bpLocs[k]-prv+1)$1+
					 "px; background-color: "+
					 colorList[k%Columns(colorList)]+
					 "; color: "+
					 fcolorList[k%Columns(colorList)]+
					 "; text-align: right; font-size: 10px;'>";
					 
				if (k<Columns (bpLocs)-1)
				{
					sp = sp + (bpLocs[k] + 1);	
				}
				else
				{
					sp = sp + "&nbsp";	
				}
				sp = sp + "</td>";
				prv = bpLocs[k];
			}	
			sp = sp + "</tr></table>";
		}
		else
		{
			ci 		   = "";
			sp		   = "<table><tr><td style = 'font-size:10px;width: 406px;background-color:"+colorList[0]+"; color:"+fcolorList[0]+"'>1-"+allData.sites+"</td></tr></table>";
		}
		htmlAICTable * ("\n<tr><td>"+ _partCount+ 
							  "</td><td><div style = 'width: "+100*currentAIC/byBPImprovement [0]$1+"%; background-color: purple; color: white;'>"+currentAIC+ 
							  "</div></td><td>"+ ci+ 
							  "</td><td>"+ sp+
							  "</td></tr>");
		
	}

	htmlAICTable * "\n</table></div></body></html>";
	htmlAICTable * 0;
	return htmlAICTable;
}
