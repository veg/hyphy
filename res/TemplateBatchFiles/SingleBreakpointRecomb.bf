partCount				= 2;

/* ________________________________________________________________________________________________*/

REPLACE_TREE_STRUCTURE  = 1;
bppMap					= {};
tree1AVL				= {};
tree2AVL				= {};
SHORT_MPI_RETURN		= 1;
VERBOSITY_LEVEL			= -1;

/*--------------------------------------------------------------------------------------*/

function testLRT (vec1, vec2)
{
	size1 = Columns(vec1);
	
	sumVec1 = {size1,1};
	jvec	= {2,size1};
	resMx1	= {itCount,1};
	resMx2	= {itCount,1};
	
	for (k=0; k<size1; k=k+1)
	{
		sumVec1 [k]	   = 1;
		jvec	[0][k] = Log(vec1[k]);
		jvec	[1][k] = Log(vec2[k]);
	}
	
	
	for (k=0; k<itCount; k=k+1)
	{
		resampled = Random(jvec,1);
		resampled = resampled*sumVec1;
		resMx1[k] = resampled[0];
		resMx2[k] = resampled[1];
	}
	
	resMx1 = (resMx1-resMx2)*2;
	resMx1 = resMx1 % 0;
	for (k=0; k<Rows(resMx1); k=k+1)
	{
		if (resMx1[k]>0)
		{
			break;
		}
	}	
	return k/itCount;
}

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
					fromNode = ReceiveJobs (0,0,0);
					break;	
				}
			}
			if (nodeCounter == MPI_NODE_COUNT-1)
			{
				break;
			}
		}			
		if (khOption)
		{
			for (bpi=0; bpi<Abs(bppMap); bpi=bpi+1)
			{
				bpv = bppMap[bpi];
				if (ResamplesDone[bpv] == 0)
				{
					if (Abs (MatrixList1[bpi]) && Abs (MatrixList2[bpi]))
					{
						runKHResampler (bpi);
					}
				}
			}		
		}
	}
	return 0;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function runKHResampler (ji)
{
	psize1 = bppMap[ji]+1;
	psize2 = siteCount - psize1;
	v1_1 = {1,psize1};
	v1_2 = {1,psize1};
	v2_1 = {1,psize2};
	v2_2 = {1,psize2};
	s1 	 = MatrixList1[ji];
	s2 	 = MatrixList2[ji];
	L_L  = {2,1};
	for (i = 0; i < psize1; i=i+1)
	{
		v1_1[i] = s1[i];
		v1_2[i] = s2[i];
		L_L[0] = L_L[0]+Log(s1[i]/s2[i]);
	}
	for (i = 0; i < psize2; i=i+1)
	{
		i2 = i+psize1;
		v2_1[i] = s1[i2];
		v2_2[i] = s2[i2];
		L_L[1] = L_L[1]+Log(s1[i2]/s2[i2]);
	}
	MasterList[ji-1][2] = testLRT (v1_1,v1_2);
	MasterList[ji-1][3] = testLRT (v2_1,v2_2);
	ResamplesDone[ji] = 1;
	
	fprintf (stdout, "\nKH::::Breakpoint at position ", psize1-1, ". Part 1 p=",  MasterList[ji-1][2], " (delta LnL = ", L_L[0],
														"). Part 2 p=", MasterList[ji-1][3], " (delta LnL = ", L_L[1],")\n");
	fprintf (khFilePath, psize1-1, ",", MasterList[ji-1][2], ",", L_L[0], ",", MasterList[ji-1][3], ",", L_L[1], "\n");
	
	return 0;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function ReceiveJobs (sendOrNot, ji, jobKind)
{
	if (MPI_NODE_COUNT>1)
	{
		MPIReceive (-1, fromNode, result_String);
		mji = MPINodeState[fromNode-1][1];
		mjk = MPINodeState[fromNode-1][2];
		
		if (sendOrNot)
		{
			if (jobKind < 2)
			{
				MPISend (fromNode,lf);
			}
			else
			{
				MPISend (fromNode,lf2);
			}
			MPINodeState[fromNode-1][1] = ji;			
			MPINodeState[fromNode-1][2] = jobKind;			
		}
		else
		{
			MPINodeState[fromNode-1][0] = 0;
			MPINodeState[fromNode-1][1] = -1;		
			MPINodeState[fromNode-1][2] = 0;		
		}
		
		ExecuteCommands (result_String);
		ji = mji;
		jobKind = mjk; 
		if (jobKind < 2)
		{
			myDF 	= lf_MLES[1][1]+baseParams;
			myAICc  = -2*(lf_MLES[1][0]-myDF*(baseSites/(baseSites-myDF-1)));
			myAIC   = -2*(lf_MLES[1][0]-myDF);
			myBIC	= -2*(lf_MLES[1][0]-myDF*Log(baseSites));
		}
	}
	else
	{
		if (jobKind < 2)
		{
			myDF 	= res[1][1]+baseParams;
			myAICc  = -2*(res[1][0]-myDF*(baseSites/(baseSites-myDF-1)));
			myAIC   = -2*(res[1][0]-myDF);
			myBIC	= -2*(res[1][0]-myDF*Log(baseSites));
		}
	}
	
	if (jobKind == 2)
	{
		ConstructCategoryMatrix (siteLikelihoods, lf2, COMPLETE);		
		MatrixList2 [ji] = siteLikelihoods;
	}
	else
	{
		fprintf (stdout, "\nBreakpoint at position ", Format(bppMap[ji],6,0), ". dAIC =  ", Format(nullAIC-myAIC,10,2), " dAICc = ",Format(nullAICc-myAICc,10,2)," dBIC = ",Format(nullBIC-myBIC,10,2), );	
		fprintf (resFilePath, bppMap[ji], ",", myAIC, ",", nullAIC-myAIC, "\n");
	
		MasterList [ji-1][0] = ji;
		MasterList [ji-1][1] = myAIC;
		MasterList [ji-1][4] = myAICc;
		MasterList [ji-1][5] = myBIC;
		if (jobKind == 1)
		{
			ConstructCategoryMatrix (siteLikelihoods, lf, COMPLETE);		
			MatrixList1 [ji] = siteLikelihoods;	
		}	
	}
	
	if (jobKind)
	{
		for (bpi=0; bpi<Abs(bppMap); bpi=bpi+1)
		{
			if (ResamplesDone[bpi] == 0)
			{
				if (Abs (MatrixList1[bpi]) > 0 && Abs (MatrixList2[bpi]) > 0)
				{
					runKHResampler (bpi);
				}
			}
		}
	}
	return fromNode-1;
}



/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function RunASample (jobIndex)
{	

	ConstraintString = "";
	LikelihoodFunctionString = "";
	ConstraintString * 8192;
	
	bpF2 = bppMap[jobIndex];
	if (dataType)
	{
		filterString = "0-"+(3*bpF2-1);
	}
	else
	{
		filterString = "0-"+bpF2;		
	}
	
	ConstraintString * ("DataSetFilter filteredData = CreateFilter(ds,1,\""+filterString+"\");");
	if (dataType)
	{
		ConstraintString * ("DataSetFilter filteredData0 = CreateFilter (ds,3,\""+filterString+"\",\"\",GeneticCodeExclusions);");		
	}
	else
	{
		ConstraintString * ("DataSetFilter filteredData0 = CreateFilter (ds,1,\""+filterString+"\");");
	}
	
	ConstraintString * ("InferTreeTopology (0);treeString=TreeMatrix2TreeString(0);Tree givenTree0 = treeString;");
					
	
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
		ConstraintString * ("DataSetFilter filteredData1 = CreateFilter (ds,3,\""+filterString+"\",\"\",GeneticCodeExclusions);");		
	}
	else
	{
		ConstraintString * ("DataSetFilter filteredData1 = CreateFilter (ds,1,\""+filterString+"\");");
	}
	ConstraintString * ("InferTreeTopology (0);treeString=TreeMatrix2TreeString(0);Tree givenTree1 = treeString;");
	ConstraintString * 0;

	ExecuteCommands (ConstraintString);
	
	fprintf (treeFilePath1,bpF2, "\t", Format (givenTree0,0,0), "\n");
	tree1AVL [jobIndex] = Format(givenTree0,0,0);
	fprintf (treeFilePath2,bpF2, "\t", Format (givenTree1,0,0), "\n");
	tree2AVL [jobIndex] = Format(givenTree1,0,0);

	if ((MPI_NODE_COUNT>1) && (jobIndex>=0))
	{
		OPTIMIZE_SUMMATION_ORDER = 0;
	}
	LikelihoodFunction lf  = (filteredData0,givenTree0,filteredData1,givenTree1);
	
	if (khOption > 1)
	{
		Tree j_tree_0 = baseTreeString;
		Tree j_tree_1 = baseTreeString;
	}
	
		
	if ((MPI_NODE_COUNT>1) && (jobIndex>=0))
	{
		OPTIMIZE_SUMMATION_ORDER = 1;
		for (mpiNode = 0; mpiNode < MPI_NODE_COUNT-1; mpiNode = mpiNode+1)
		{
			if (MPINodeState[mpiNode][0]==0)
			{
				break;	
			}
		}
		if (mpiNode==MPI_NODE_COUNT-1)
		{
			mpiNode = ReceiveJobs (1,jobIndex,khOption>0);
		}
		else
		{
			MPISend (mpiNode+1,lf);
			MPINodeState[mpiNode][0] = 1;
			MPINodeState[mpiNode][1] = jobIndex;
			MPINodeState[mpiNode][2] = khOption>0;
		}
		
		if (khOption)
		{
			if (khOption == 1)
			{
				LikelihoodFunction lf2 = (filteredData0,givenTree1,filteredData1,givenTree0);
			}
			else
			{
				LikelihoodFunction lf2 = (filteredData0,j_tree_0,filteredData1,j_tree_1);			
			}
			
			for (mpiNode = 0; mpiNode < MPI_NODE_COUNT-1; mpiNode = mpiNode+1)
			{
				if (MPINodeState[mpiNode][0]==0)
				{
					break;	
				}
			}
			if (mpiNode==MPI_NODE_COUNT-1)
			{
				mpiNode = ReceiveJobs (1,jobIndex,2);
			}
			else
			{
				MPISend (mpiNode+1,lf2);
				MPINodeState[mpiNode][0] = 1;
				MPINodeState[mpiNode][1] = jobIndex;
				MPINodeState[mpiNode][2] = 2;
			}		
		}
	}
	else
	{
		Optimize (res,lf);
		
		if (jobIndex>=0)
		{
			mpiNode = ReceiveJobs (1, jobIndex, khOption>0);
			if (khOption)
			{
				if (khOption == 1)
				{
					LikelihoodFunction lf2 = (filteredData0,givenTree1,filteredData1,givenTree0);
				}
				else
				{
					LikelihoodFunction lf2 = (filteredData0,j_tree_0,filteredData1,j_tree_1);			
				}
				Optimize (res2,lf2);
				mpiNode = ReceiveJobs (1, jobIndex, 2);
			}
		}
		else
		{
			myAIC = 2*(res[1][0]-res[1][1]-baseParams);
		}
	}
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

DataSetFilter	filteredData = CreateFilter(ds,1);


siteCount = filteredData.sites;

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

if (Abs(bppMap) <= partCount)
{
	fprintf (stdout, "\nThere are too few potential break points to support ", partCount-1, " recombination events.\n");
	return 0;
}

bppSize = (Log(Abs(bppMap))/Log(2)+1)$1;
fprintf (stdout, "\nThere are ",Abs(bppMap)," potential breakpoints.\n");

ChoiceList (khOption,"KH Testing",1,SKIP_NONE,"Skip","Use only AIC to measure goodness of fit.",
				     "Run 1","Verify conflicting phylogenetic signal with KH resampling, swapping trees between partitions for the test",
				     "Run 2","Verify conflicting phylogenetic signal with KH resampling, using the joint tree as the null.");


if (khOption < 0)
{
	return 0;
}

if (khOption > 0)
{
	itCount = 0;
	while (itCount < 1)
	{
		fprintf (stdout, "How many KH samples should be drawn per breakpoint ?");
		fscanf (stdin, "Number", itCount);
	}
	SHORT_MPI_RETURN = 0;	
}

SelectTemplateModel (filteredData);

SetDialogPrompt ("Save results to:");

fprintf (PROMPT_FOR_FILE,CLEAR_FILE);

resFilePath = LAST_FILE_PATH;
treeFilePath1 = resFilePath + ".trees1";
treeFilePath2 = resFilePath + ".trees2";

fprintf (treeFilePath1,CLEAR_FILE);
fprintf (treeFilePath2,CLEAR_FILE);
if (khOption)
{
	khFilePath = resFilePath + ".kh";
	fprintf (khFilePath, CLEAR_FILE, "Breakpoint, Part 1 KH p-value, Part 1 Delta LogL, Part 2 KH p-value, Part 2 Delta LogL\n");
}

baseTreeString = treeString;
Tree givenTree = treeString;

branchNames	  = BranchName (givenTree,-1);
LikelihoodFunction lf = (filteredData,givenTree);
Optimize (res,lf);
currentPopulation  = {};

baseSites  = filteredData.sites;
baseParams = res[1][2] + 3; /* for frequencies */

myDF	 = res[1][1];
nullAIC  = -2*(res[1][0]-myDF);
nullAICc = -2*(res[1][0]-myDF*(baseSites/(baseSites-myDF-1)));
nullBIC	 = -2*(res[1][0]-myDF*Log(baseSites));

MasterList    = {Abs(bppMap)-1,6};
MatrixList1   = {};
MatrixList2   = {};
ResamplesDone = {};

fprintf (stdout, "\n1). Single partition analysis\n",lf,
				 "\n AIC  = ",nullAIC,  "\n",
				 "\nc-AIC = ",nullAICc, "\n",
				 "\n BIC  = ",nullBIC,  "\n",
				 "2). Looking for a breakpoint...\n");

baseParams = res[1][2];
	
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
if (MPI_NODE_COUNT>1)
{
	MPINodeState = {MPI_NODE_COUNT-1,3};
}


for (individual=1; individual<Abs(bppMap); individual=individual+1)
{
	RunASample (individual);
}

CleanUpMPI (0);

/* process the results */

sortedScores = MasterList % 1;

fprintf (stdout, "\nAIC");

if (sortedScores[0][1] < nullAIC)
{
	individual = sortedScores[0][0];
	fprintf (stdout, "\n\nBest supported breakpoint is located at position ", bppMap[individual], "\nAIC = ", sortedScores[0][1], " : an improvement of ", 
				     nullAIC - sortedScores[0][1], " AIC points\n");
	bpF2 = bppMap[individual];
	compactResultPath = resFilePath + "_AIC.splits";
	fprintf (compactResultPath, CLEAR_FILE, "\n\n0-", bpF2, "\n", tree1AVL[individual], "\n", bpF2+1,"-", ds.sites-1, "\n", tree2AVL[individual]);
}
else
{
	fprintf (stdout, "\n\nThere seems to be NO recombination in this alignment\n\n");
}

sortedScores = MasterList % 4;

fprintf (stdout, "\nAIC-c");

if (sortedScores[0][4] < nullAICc)
{
	individual = sortedScores[0][0];
	fprintf (stdout, "\n\nBest supported breakpoint is located at position ", bppMap[individual], "\nAIC = ", sortedScores[0][4], " : an improvement of ", 
				     nullAICc - sortedScores[0][4], " AIC points\n");
	bpF2 = bppMap[individual];
	compactResultPath = resFilePath + "_cAIC.splits";
	fprintf (compactResultPath, CLEAR_FILE, "\n\n0-", bpF2, "\n", tree1AVL[individual], "\n", bpF2+1,"-", ds.sites-1, "\n", tree2AVL[individual]);
}
else
{
	fprintf (stdout, "\n\nThere seems to be NO recombination in this alignment\n\n");
}

sortedScores = MasterList % 5;

fprintf (stdout, "\nBIC");

if (sortedScores[0][5] < nullBIC)
{
	individual = sortedScores[0][0];
	fprintf (stdout, "\n\nBest supported breakpoint is located at position ", bppMap[individual], "\nAIC = ", sortedScores[0][5], " : an improvement of ", 
				     nullBIC - sortedScores[0][5], " AIC points\n");
	bpF2 = bppMap[individual];
	compactResultPath = resFilePath + "_BIC.splits";
	fprintf (compactResultPath, CLEAR_FILE, "\n\n0-", bpF2, "\n", tree1AVL[individual], "\n", bpF2+1,"-", ds.sites-1, "\n", tree2AVL[individual]);
}
else
{
	fprintf (stdout, "\n\nThere seems to be NO recombination in this alignment\n\n");
}	
	
return 0;



ChoiceList (dataType,"Verify by LRT",1,SKIP_NONE,"Sure","Verify the statistical significance of the breakpoint using the likelihood ratio test; this may be slow",
				     "Skip","Skip the LRT");
				     
if (dataType != 0)
{
	return 0;
}

itCount = 0;
while (itCount < 1)
{
	fprintf (stdout, "How many samples should be used to tabulate the LRT statistic ?");
	fscanf (stdin, "Number", itCount);
}

bestTree1 = tree1AVL[individual];
bestTree2 = tree2AVL[individual];

if (dataType)
{
	DataSetFilter filter1 = CreateFilter (ds,3,siteIndex<3*bpF2,,GeneticCodeExclusions);
	DataSetFilter filter2 = CreateFilter (ds,3,siteIndex>=3*bpF2,,GeneticCodeExclusions);
}
else
{
	DataSetFilter filter1 = CreateFilter (ds,1,siteIndex<=bpF2);
	DataSetFilter filter2 = CreateFilter (ds,1,siteIndex>bpF2);
}

SelectTemplateModel (filteredData);

lrtFilePath = resFilePath + ".lrt";
fprintf (lrtFilePath, CLEAR_FILE);

fprintf (stdout, "\nFitting the null hypothesis (same tree, different models)\n");

if (FREQUENCY_SENSITIVE)
{
	HarvestFrequencies (partFreqs1,filter1,1+2*USE_POSITION_SPECIFIC_FREQS,1,1);
	HarvestFrequencies (partFreqs2,filter2,1+2*USE_POSITION_SPECIFIC_FREQS,1,1);

	MULTIPLY_BY_FREQS = PopulateModelMatrix ("partMatrix1",partFreqs1);
	MULTIPLY_BY_FREQS = PopulateModelMatrix ("partMatrix2",partFreqs2);
	
	if (dataType)
	{
		codonEFV1 = BuildCodonFrequencies (partFreqs1);
		codonEFV2 = BuildCodonFrequencies (partFreqs2);
		Model partModel1 = (partMatrix1,codonEFV1,MULTIPLY_BY_FREQS);	
		Model partModel2 = (partMatrix2,codonEFV2,MULTIPLY_BY_FREQS);	
	}
	else
	{
		Model partModel1 = (partMatrix1,partFreqs1,MULTIPLY_BY_FREQS);	
		Model partModel2 = (partMatrix2,partFreqs2,MULTIPLY_BY_FREQS);	
	}	
}
else
{
	fprintf (stdout, "ERROR: Models which do not incorporate empirical base frequencies have not yet been included in this analysis\n\n");
	return 0;
}

UseModel (partModel1);
Tree	 part1TreeNull = baseTreeString;
LikelihoodFunction lf_null_1 = (filter1, part1TreeNull);
Optimize (res_null_1, lf_null_1);
fprintf (stdout, "\nPart 1 fit:\n",lf_null_1);

UseModel (partModel2);
Tree	 part2TreeNull = baseTreeString;
LikelihoodFunction lf_null_2 = (filter2, part2TreeNull);
Optimize (res_null_2, lf_null_2);
fprintf (stdout, "\n\nPart 2 fit:\n",lf_null_2);

fprintf (stdout, "\n\nFitting the alternative hypothesis (different tree, different models)\n");

UseModel (partModel1);
Tree	 part1TreeAlt = bestTree1;
LikelihoodFunction lf = (filter1, part1TreeAlt);
Optimize (res_alt_1, lf);
fprintf (stdout, "\nPart 1 fit:\n",lf);

UseModel (partModel2);
Tree	 part2TreeAlt = bestTree2;
LikelihoodFunction lf = (filter2, part2TreeAlt);
Optimize (res_alt_2, lf);
fprintf (stdout, "\n\nPart 2 fit:\n",lf);

obsLRT = 2(res_alt_1[1][0]+res_alt_2[1][0]-res_null_1[1][0]-res_null_2[1][0]);
fprintf (stdout, "\n\nObserved LRT:", obsLRT, "\n\n"); 

fprintf (stdout, "\nStarting the simulations...\n"
				 "\nSimulation     LRT Running p-value\n");

fprintf (lrtFilePath, "Simulation, LRT\n");
fprintf (lrtFilePath, "Observed", obsLRT, "\n");

SimLRTs = {itCount,5};

SHORT_MPI_RETURN = 1;
finishedIterates = 0;
runningPValue	 = 0;

LIKELIHOOD_FUNCTION_OUTPUT = 6;

for (itCounter=0; itCounter<itCount; itCounter=itCounter+1)
{
	/* restore global variables of lf_null_1 */
	for (k2 = 0; k2 < res_null_1[1][2]; k2 = k2+1)
	{
		SetParameter (lf_null_1,k2,res_null_1[0][k2]);
	}
	
	/* simulate null data in partition 1 */
	DataSet simmed1 = SimulateDataSet (lf_null_1); 
	/*fprintf (stdout, "Simulation ", itCounter,"/1\n",lf_null_1);*/

	/* restore global variables of lf_null_2 */
	for (k2 = 0; k2 < res_null_2[1][2]; k2 = k2+1)
	{
		SetParameter (lf_null_2,k2,res_null_2[0][k2]);
	}
	/*fprintf (stdout, "Simulation ", itCounter,"/2\n",lf_null_2);*/
	
	/* simulate null data in partition 2 */
	DataSet simmed2 = SimulateDataSet		(lf_null_2); 
	
	DataSetFilter filteredData = CreateFilter (simmed1,1);
	InferTreeTopology (0);
	simTree1 = TreeMatrix2TreeString(0);
	
	seqNameMap = {};
	GetInformation (filterStrings, filteredData);
	/* merge the alignments first, by appropriately mapping sequence names */
	for (k2 = 0; k2 < filteredData.species; k2=k2+1)
	{
		GetString (seqName, filteredData, k2);
		seqNameMap [seqName] = filterStrings[k2];
	}

	DataSetFilter filteredData = CreateFilter (simmed2,1);
	InferTreeTopology (0);
	simTree2 = TreeMatrix2TreeString(0);
	
	GetInformation (filterStrings, filteredData);
	/* merge the alignments first, by appropriately mapping sequence names */
	for (k2 = 0; k2 < filteredData.species; k2=k2+1)
	{
		GetString (seqName, filteredData, k2);
		seqNameMap [seqName] = seqNameMap [seqName]+filterStrings[k2];
	}
	
	assembledString = "";
	assembledString * 8192;
	
	seqNames = Rows (seqNameMap);
	
	for (k2 = 0; k2 < Abs (seqNameMap); k2 = k2+1)
	{
		seqName = seqNames[k2];
		assembledString * (">"+seqName+"\n");
		assembledString * seqNameMap[seqName];
		assembledString * "\n";
	}
	assembledString * 0;
	DataSet jointD = ReadFromString (assembledString);
	DataSetFilter filteredData = CreateFilter (jointD,1);
	InferTreeTopology (0);
	simTreeJ = TreeMatrix2TreeString(0);
	
	/* make a joint tree as well */

	if (dataType)
	{
		DataSetFilter simFilter1 = CreateFilter (simmed1,3,"","",GeneticCodeExclusions);
		DataSetFilter simFilter2 = CreateFilter (simmed2,3,"","",GeneticCodeExclusions);
	}
	else
	{
		DataSetFilter simFilter1 = CreateFilter (simmed1,1);
		DataSetFilter simFilter2 = CreateFilter (simmed2,1);
	}

	if (FREQUENCY_SENSITIVE)
	{
		HarvestFrequencies (simFreqs1,simFilter1,1+2*USE_POSITION_SPECIFIC_FREQS,1,1);
		HarvestFrequencies (simFreqs2,simFilter2,1+2*USE_POSITION_SPECIFIC_FREQS,1,1);

		MULTIPLY_BY_FREQS = PopulateModelMatrix ("simMatrix1",simFreqs1);
		MULTIPLY_BY_FREQS = PopulateModelMatrix ("simMatrix2",simFreqs2);
		
		if (dataType)
		{
			simCodonEFV1 = BuildCodonFrequencies (simFreqs1);
			simCodonEFV2 = BuildCodonFrequencies (simFreqs2);
			Model simModel1 = (simMatrix1,simCodonEFV1,MULTIPLY_BY_FREQS);	
			Model simModel2 = (simMatrix2,simCodonEFV2,MULTIPLY_BY_FREQS);	
		}
		else
		{
			Model simModel1 = (simMatrix1,simFreqs1,MULTIPLY_BY_FREQS);	
			Model simModel2 = (simMatrix2,simFreqs2,MULTIPLY_BY_FREQS);	
		}	
	}
	RunASample_LRT (itCounter,0);
	RunASample_LRT (itCounter,1);
	RunASample_LRT (itCounter,2);
	RunASample_LRT (itCounter,3);
}

CleanUpMPI_LRT (0);

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function RunASample_LRT (iterateIndex, lfKind)
{	

	if (lfKind == 0) /* null 1st part */
	{
		UseModel (simModel1);
		Tree	 daTree = simTreeJ;
		LikelihoodFunction lf = (simFilter1, daTree);
	}
	if (lfKind == 1) /* null 2nd part */
	{
		UseModel (simModel2);
		Tree	 daTree = simTreeJ;
		LikelihoodFunction lf = (simFilter2, daTree);
	}	
	if (lfKind == 2) /* alt 1st part */
	{
		UseModel (simModel1);
		Tree	 daTree = simTree1;
		LikelihoodFunction lf = (simFilter1, daTree);
	}
	if (lfKind == 3) /* alt 2nd part */
	{
		UseModel (simModel2);
		Tree	 daTree = simTree2;
		LikelihoodFunction lf = (simFilter2, daTree);
	}
	
		
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
			mpiNode = ReceiveJobs_LRT (1,iterateIndex,lfKind);
		}
		else
		{
			MPISend (mpiNode+1,lf);
			MPINodeState[mpiNode][0] = 1;
			MPINodeState[mpiNode][1] = iterateIndex;
			MPINodeState[mpiNode][2] = lfKind;
		}
	}
	else
	{
		Optimize (lf_MLES,lf);
		/*fName = lrtFilePath + "." + iterateIndex + "." + lfKind;
		fprintf (fName, CLEAR_FILE,lf);*/
		mpiNode = ReceiveJobs_LRT (0,iterateIndex,lfKind);
	}
	return 0;	
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function CleanUpMPI_LRT (dummy)
{
	if (MPI_NODE_COUNT>1)
	{
		while (1)
		{
			for (nodeCounter = 0; nodeCounter < MPI_NODE_COUNT-1; nodeCounter = nodeCounter+1)
			{
				if (MPINodeState[nodeCounter][0]==1)
				{
					fromNode = ReceiveJobs_LRT (0,0,0);
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

function ReceiveJobs_LRT (sendOrNot, ji, jobKind)
{
	if (MPI_NODE_COUNT>1)
	{
		MPIReceive (-1, fromNode, result_String);
		mji = MPINodeState[fromNode-1][1];
		mjk = MPINodeState[fromNode-1][2];
		
		if (sendOrNot)
		{
			MPISend (fromNode,lf);
			/*fTemp = resFilePath + "." + ji + "." + jobKind;
			fprintf (fTemp, CLEAR_FILE, MPI_LAST_SENT_MSG);*/
			MPINodeState[fromNode-1][1] = ji;			
			MPINodeState[fromNode-1][2] = jobKind;			
		}
		else
		{
			MPINodeState[fromNode-1][0] = 0;
			MPINodeState[fromNode-1][1] = -1;		
			MPINodeState[fromNode-1][2] = 0;		
		}
		
		ExecuteCommands (result_String);
		
		ji = mji;
		jobKind = mjk; 
	}

	SimLRTs[ji][jobKind] = lf_MLES[1][0];
	
	if (SimLRTs[ji][0]<0 && SimLRTs[ji][1]<0 && SimLRTs[ji][2]<0 && SimLRTs[ji][3]<0)
	{
		finishedIterates = finishedIterates + 1;
		/* finished an iterate */	
		SimLRTs[ji][4] = 2*(SimLRTs[ji][2]+SimLRTs[ji][3]-SimLRTs[ji][0]-SimLRTs[ji][1]);
		if (SimLRTs[ji][4] > obsLRT)
		{
			runningPValue = runningPValue + 1;
		}
		fprintf (stdout, Format (finishedIterates,10,0), " ",
						 Format (SimLRTs[ji][4], 7, 3), " ",
						 Format (runningPValue/finishedIterates, 15, 6), "\n");
						 
		fprintf (lrtFilePath,  Format (finishedIterates,10,0), ",",Format (SimLRTs[ji][4], 7, 3),"\n");
	}
	
	return fromNode-1;
}
