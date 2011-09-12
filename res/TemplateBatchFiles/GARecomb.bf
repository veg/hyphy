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
shortestAllowedSegment  = 20;

stoppingCriterion		= 100;
sampleCount				= 0;
familyControlSize		= produceOffspring$6;

verboseFlag				= 0;
stateVectorDimension    = 120;
rateClassesCount		= 2;

/* ________________________________________________________________________________________________*/

MasterList				= {};
REPLACE_TREE_STRUCTURE  = 1;
bppMap					= {};
SHORT_MPI_RETURN		= 1;
totalBitSize			= 0;

/* ________________________________________________________________________________________________*/


function InitializeDistances (dummy)
{
	HarvestFrequencies (_dNucFreq,filteredData,1,1,0);
	_d_fR = _dNucFreq[0]+_dNucFreq[2];
	_d_fY = _dNucFreq[1]+_dNucFreq[3];
	
	if (_dNucFreq[0] == 0 || _dNucFreq[1] == 0 || _dNucFreq[2] == 0 || _dNucFreq[3] == 0)
	{
		_useK2P = 1;
	}
	else
	{
		_d_TN_K1 = 2*_dNucFreq[0]*_dNucFreq[2]/_d_fR;
		_d_TN_K2 = 2*_dNucFreq[1]*_dNucFreq[3]/_d_fY;
		_d_TN_K3 = 2*(_d_fR*_d_fY-_dNucFreq[0]*_dNucFreq[2]*_d_fY/_d_fR-_dNucFreq[1]*_dNucFreq[3]*_d_fR/_d_fY);
		_useK2P = 0;
	}
	
	
	summingVector = {{1}{1}{1}{1}};

	return 0;
}

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


/* ________________________________________________________________________________________________*/

function TreeMatrix2TreeString (doLengths)
{
	treeString = "";
	p = 0;
	k = 0;
	m = treeNodes[0][1];
	n = treeNodes[0][0];
	d = treeString*(Rows(treeNodes)*25);

	while (m)
	{	
		if (m>p)
		{
			if (p)
			{
				d = treeString*",";
			}
			for (j=p;j<m;j=j+1)
			{
				d = treeString*"(";
			}
		}
		else
		{
			if (m<p)
			{
				for (j=m;j<p;j=j+1)
				{
					d = treeString*")";
				}
			}	
			else
			{
				d = treeString*",";
			}	
		}
		if (n<filteredData.species)
		{
			GetString (nodeName, filteredData, n);
			d = treeString*nodeName;
		}
		if (doLengths>.5)
		{
			nodeName = ":"+treeNodes[k][2];
			d = treeString*nodeName;
		}
		k=k+1;
		p=m;
		n=treeNodes[k][0];
		m=treeNodes[k][1];
	}

	for (j=m;j<p;j=j+1)
	{
		d = treeString*")";
	}
	
	d=treeString*0;
	return treeString;
}

/* ________________________________________________________________________________________________*/



function InferTreeTopology(verbFlag)
{
	distanceMatrix = {filteredData.species,filteredData.species};
	dummy = InitializeDistances (0);
		
	for (i = 0; i<filteredData.species; i=i+1)
	{
		for (j = i+1; j<filteredData.species; j = j+1)
		{
			distanceMatrix[i][j] = ComputeDistanceFormula (i,j);
		}
	}

	MESSAGE_LOGGING 		 	= 1;
	cladesMade 					= 1;
	

	if (filteredData.species == 2)
	{
		d1 = distanceMatrix[0][1]/2;
		treeNodes = {{0,1,d1__},
					 {1,1,d1__},
					 {2,0,0}};
					 
		cladesInfo = {{2,0}};
	}
	else
	{
		if (filteredData.species == 3)
		{
			d1 = (distanceMatrix[0][1]+distanceMatrix[0][2]-distanceMatrix[1][2])/2;
			d2 = (distanceMatrix[0][1]-distanceMatrix[0][2]+distanceMatrix[1][2])/2;
			d3 = (distanceMatrix[1][2]+distanceMatrix[0][2]-distanceMatrix[0][1])/2;
			treeNodes = {{0,1,d1__},
						 {1,1,d2__},
						 {2,1,d3__}
						 {3,0,0}};
						 
			cladesInfo = {{3,0}};		
		}
		else
		{	
			njm = (distanceMatrix > methodIndex)>=filteredData.species;
				
			treeNodes 		= {2*(filteredData.species+1),3};
			cladesInfo	    = {filteredData.species-1,2};
			
			for (i=Rows(treeNodes)-1; i>=0; i=i-1)
			{
				treeNodes[i][0] = njm[i][0];
				treeNodes[i][1] = njm[i][1];
				treeNodes[i][2] = njm[i][2];
			}

			for (i=Rows(cladesInfo)-1; i>=0; i=i-1)
			{
				cladesInfo[i][0] = njm[i][3];
				cladesInfo[i][1] = njm[i][4];
			}
			
			njm = 0;
		}
	}
	distanceMatrix = 0;
	
	return 1.0;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function StringToMatrix (zz)
{
	return zz;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function ExportAMatrix (fileName, rateMatrix,dummy)
{
	if (dummy == 1)
	{
		return 0;
	}
	
	sortedBP = ConvertToPart (rateMatrix);
	
	if (dataType)
	{
		bpF = 0;
	}	
	else
	{
		bpF = -1;
	}
	ConstraintString = "";
	ConstraintString * 8192;
	
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
		}
		DataSetFilter filteredData = CreateFilter(ds,1,filterString);
		InferTreeTopology (0);
		treeString=TreeMatrix2TreeString(0);
		ConstraintString * (filterString+"\n"+treeString+"\n\n");
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
		}
		DataSetFilter filteredData = CreateFilter(ds,1,filterString);
		InferTreeTopology (0);
		treeString=TreeMatrix2TreeString(0);
		ConstraintString * (filterString+"\n"+treeString+"\n\n");
	}

	ConstraintString * 0;
	
	fprintf (fileName, CLEAR_FILE,ConstraintString);
	return 0;
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
		
		if (sendOrNot)
		{
			MPISend (fromNode,lf);
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
		ji = mji;
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
			sortedBP = ConvertToPart (currentPopulation[ji]);
			v = Rows (sortedBP);
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
			sortedBP = ConvertToPart (children[ji-populationSize]);
			v = Rows (sortedBP);
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

function ConvertToPart (pString)
{
	partitionHits    = {};
	h = 0; 
	v = bppSize;
	
	for (mpiNode=0; mpiNode<partCount; mpiNode=mpiNode+1)
	{
		aBP    = 0;
		bpF	   = 2^(bppSize-1);
		
		for (; h<v; h=h+1)
		{
			aBP = aBP + bpF*(0+pString[h]);
			bpF = bpF/2;
		}
		
		if (aBP>Abs(bppMap))
		{
			aBP = aBP - 2^(bppSize-1);
		}
		
		v = v + bppSize;
		
		if (partitionHits[aBP] == 0)
		{
			partitionHits[aBP] = 1;
		}
	}

	meKeys	 = Rows(partitionHits);
	
	v = Columns(meKeys);
	
	sortedBP = {v,1};
	
	for (h=0; h<v; h=h+1)
	{
		sortedBP [h] = 0+meKeys[h];
	}
	sortedBP = sortedBP % 0;
	return sortedBP;
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
		bpF2 = sortedBP[h];
		bpF2 = bppMap[bpF2];
		if ((h>0) + sortedBP[h])
		{
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
		}
		else
		{
			bpF2 = 0;
			curSegLength = 0;
		}
		if (curSegLength < minPartLength && curSegLength>0)
		{
			minPartLength = curSegLength;
		}
	}
	
	if ((bpF2<ds.sites && (dataType == 0)) || (bpF2<3*ds.sites && dataType))
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

function RunASample (dummy,jobIndex)
{	
	sortedBP = ConvertToPart (cString);
	v = Rows (sortedBP);

	if (dataType)
	{
		bpF = 0;	
	}
	else
	{
		bpF = -1;
	}
	ConstraintString = "";
	LikelihoodFunctionString = "";
	ConstraintString * 8192;
	LikelihoodFunctionString * 256;
	LikelihoodFunctionString* "LikelihoodFunction lf=(";
	

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

	ConstraintString * 0;
	LikelihoodFunctionString * 0;
	
	partMap = ConvertToPartString (cString);
	myAIC = MasterList[partMap];
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
	}
	else
	{
		ExecuteCommands (ConstraintString);
		ExecuteCommands (LikelihoodFunctionString);
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
				MPISend (mpiNode+1,lf);
				MPINodeState[mpiNode][0] = 1;
				MPINodeState[mpiNode][1] = jobIndex;
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
				myAIC = 2*(res[1][0]-res[1][1]-baseParams);
			}
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

	myAIC = MasterList[sampleString];
	testChild = putativeChild;
	mutPassCount = 1;
	while ((myAIC<(-0.1) || minPartLength<shortestAllowedSegment)&& mutPassCount < 25)
	{
		if (verboseFlag > 1)
		{
			fprintf (stdout,"Adjusting the child to avoid a duplicate. Min(fragment) = ",minPartLength,".  Pass ", mutPassCount, "\n");
		}
		
		mutPassCount = mutPassCount + 1;
		
		sampleString = Min(Random(0,stateVectorDimension)$1,stateVectorDimension-1);
		myAIC = testChild[sampleString];
		
		newValue = Random (0,rateClassesCount-0.0000001)$1;
		
		while (newValue == myAIC)
		{
			newValue = Random (0,rateClassesCount-0.0000001)$1;
		}
		
		testChild [sampleString] = newValue;
		sampleString = 	ConvertToPartString (testChild);
		myAIC = MasterList[sampleString];
	}
	return testChild;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function UpdateBL (dummy)
{
	return 0;
}


/*---------------------------------------------------------------------------------------------------------------------------------------------*/


ChoiceList (dataType,"Data type",1,SKIP_NONE,"Nucleotide","Nucleotide data.",
				     "Codon","Codon (several available genetic codes).");

if (dataType<0) 
{
	return;
}

if (dataType==0)
{
	SetDialogPrompt 	  		  ("Locate a nucleotide data  file:");
}
else
{
	SetDialogPrompt 	  		  ("Locate a codon data file:");
}

DataSet 		ds  		 = ReadDataFile (PROMPT_FOR_FILE);
fprintf (stdout, "\nRead ", ds.species, " sequences and ", ds.sites, " sites.\n");
DataSetFilter	filteredData = CreateFilter(ds,1);

InferTreeTopology (0);
treeString = TreeMatrix2TreeString(0);

/* find "informative sites" */

if (dataType==0)
{
	for (h=0; h<filteredData.sites; h=h+1)
	{
		filterString = "" + h;
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
else
{
	incFileName = HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"chooseGeneticCode.def";
	ExecuteCommands  ("#include \""+incFileName+"\";");
	
	DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);
	for (h=0; h<filteredData.sites; h=h+1)
	{
		filterString = "" + 3*h + "-" + (3*h+2);
		DataSetFilter siteFilter = CreateFilter (filteredData,3,filterString,"",GeneticCodeExclusions);

		HarvestFrequencies (f1, siteFilter, 3, 3, 1);
		m1 = 0;
		for (mpiNode=0; (mpiNode < 64) && (m1<=1); mpiNode=mpiNode+1)
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
fprintf (stdout, "\nThere are ",Abs(bppMap)," potential breakpoints. Bit size of the sample is ", bppSize,".\n");

fprintf (stdout, "How many breakpoints?");
fscanf (stdin, "Number", partCount);

h = Abs(bppMap);

if (h <= partCount)
{
	fprintf (stdout, "\nThere are too few potential break points to support ", partCount-1, " recombination events.\n");
	return 0;
}

totalModelCounter = 1;
kf = 1;

for (k=1; k <= partCount; k=k+1)
{
	totalModelCounter = totalModelCounter * (h-k);
	kf = kf * k;
} 
totalModelCounter = totalModelCounter / kf;
fprintf (stdout, "\nThere are ", Format(totalModelCounter,20,0), " total possible models.\n");

totalBitSize = bppSize * partCount;
stateVectorDimension = totalBitSize;

SelectTemplateModel (filteredData);

SetDialogPrompt ("Save the best model to:");
fprintf (PROMPT_FOR_FILE,CLEAR_FILE);
modelFile = LAST_FILE_PATH;

Tree givenTree = treeString;

branchNames	  = BranchName (givenTree,-1);
LikelihoodFunction lf = (filteredData,givenTree);

Optimize 				(res,lf);

/* check parameter counts */


currentPopulation  = {};
sortedScores	   = {populationSize,2};
baseSites		   = filteredData.sites;
baseParams 		   = res[1][2];

sortedScores[0][0] = 2*(res[1][0]-res[1][1]*(baseSites/(baseSites-res[1][1]-1)));
sortedScores[0][1] = 0;


currentPopulation [0] = {totalBitSize,1};
fprintf (stdout, "\n______________SINGLE RATE______________\n",lf,"\nAIC=",-sortedScores[0][0],"\n");

perPartParameterCount = (res[1][1]-res[1][2]);
if (baseParams + (partCount+1) * perPartParameterCount >= baseSites - 1)
{
	fprintf (stdout, "\nERROR: Too few sites/sequences for c-AIC inference\n");	
	nullAICFile = modelFile + ".NULL_AIC";
	fprintf (nullAICFile,CLEAR_FILE,-sortedScores[0][0]);

	nullAICFile = modelFile + ".BEST_AIC";
	fprintf (nullAICFile,CLEAR_FILE,-sortedScores[0][0]);
	return 0;
}
	
if (baseParams>0)
{

	ConstraintString = "";
	ConstraintString*256;
	for (h=0; h<baseParams; h=h+1)
	{
		GetString (v,lf,h);
		ConstraintString * (v+":="+v+"__;\n");
	}
	ConstraintString*0;
	ExecuteCommands (ConstraintString);
}

ibfPath = HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "GA_CHC.ibf";
ExecuteCommands ("#include \"" + ibfPath + "\";");


ExportAMatrix (modelFile,currentPopulation[populationSize-1],0);

nullAICFile = modelFile + ".NULL_AIC";
fprintf (nullAICFile,CLEAR_FILE,crapAIC);

nullAICFile = modelFile + ".BEST_AIC";
fprintf (nullAICFile,CLEAR_FILE,lastBestAIC);

modelFile = modelFile + ".samples";

masterKeys = Rows(MasterList);
outString = "";
outString * 65536;
for (h=Rows(masterKeys)*Columns(masterKeys)-1; h>=0; h=h-1)
{
	aKey = masterKeys[h];
	outString * (""+(-MasterList[aKey])+"\n");
	outString * aKey;
	outString * "\n";
}
outString * 0;
fprintf (modelFile,CLEAR_FILE,outString);

fprintf (stdout, "\n\nGA search explored ", Abs(MasterList), "/", Format(totalModelCounter,20,0), "(", Format (Abs(MasterList)/totalModelCounter*100,10,8),"%) models.\n");
