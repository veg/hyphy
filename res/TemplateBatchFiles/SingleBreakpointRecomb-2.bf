partCount				= 2;

/* ________________________________________________________________________________________________*/

REPLACE_TREE_STRUCTURE  = 1;
SHORT_MPI_RETURN		= 1;

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

	_d_TN_K1 = 2*_dNucFreq[0]*_dNucFreq[2]/_d_fR;
	_d_TN_K2 = 2*_dNucFreq[1]*_dNucFreq[3]/_d_fY;
	_d_TN_K3 = 2*(_d_fR*_d_fY-_dNucFreq[0]*_dNucFreq[2]*_d_fY/_d_fR-_dNucFreq[1]*_dNucFreq[3]*_d_fR/_d_fY);
	
	return 0;
}

/* ________________________________________________________________________________________________*/

function ComputeDistanceFormula (s1,s2)
{
	GetDataInfo (siteDifferenceCount, filteredData, s1, s2, DIST);
	
	_dAGCounts 	 =    siteDifferenceCount[0][2]+siteDifferenceCount[2][0]  /* A-G and G-A */;
	_dCTCounts	 = 	  siteDifferenceCount[1][3]+siteDifferenceCount[3][1]; /* C-T and T-C */
						
	_dTransversionCounts = (siteDifferenceCount[0][0]+siteDifferenceCount[1][1]+siteDifferenceCount[2][2]+
							siteDifferenceCount[3][3])+_dAGCounts+_dCTCounts;
	
	_dAGCounts	 = _dAGCounts/filteredData.sites;
	_dCTCounts	 = _dCTCounts/filteredData.sites;
	
	_dTransversionCounts = 1-_dTransversionCounts/filteredData.sites;
	
	_d1C = 1-_dAGCounts/_d_TN_K1-0.5*_dTransversionCounts/_d_fR;
	_d2C = 1-_dCTCounts/_d_TN_K2-0.5*_dTransversionCounts/_d_fY;
	_d3C = 1-0.5*_dTransversionCounts/_d_fY/_d_fR;
	
	if ((_d1C>0)&&(_d2C>0)&&(_d3C>0))
	{
		return -_d_TN_K1*Log(_d1C)-_d_TN_K2*Log(_d2C)-_d_TN_K3*Log(_d3C);
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

function ReceiveJobs (sendOrNot, ji, jobKind)
{
	if (MPI_NODE_COUNT>1)
	{
		MPIReceive (-1, fromNode, result_String);
		mji = MPINodeState[fromNode-1][1];
		mjk = MPINodeState[fromNode-1][2];
		
		if (sendOrNot)
		{
			MPISend (fromNode,lf);
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
		myDF 	= lf_MLES[1][1]+baseParams;
		myAICc  = -2*(lf_MLES[1][0]-myDF*(baseSites/(baseSites-myDF-1)));
		myAIC   = -2*(lf_MLES[1][0]-myDF);
		myBIC	= -2*(lf_MLES[1][0]-myDF*Log(baseSites));
		ji = mji;
		jobKind = mjk; 
	}
	else
	{
		myDF 	= res[1][1]+baseParams;
		myAICc  = -2*(res[1][0]-myDF*(baseSites/(baseSites-myDF-1)));
		myAIC   = -2*(res[1][0]-myDF);
		myBIC	= -2*(res[1][0]-myDF*Log(baseSites));
	}
	
	if (jobKind == 2)
	{
		ConstructCategoryMatrix (siteLikelihoods, lf2, COMPLETE);		
		MatrixList2 [ji] = siteLikelihoods;
	}
	else
	{
		fprintf (stdout, "\nBreakpoint at position ", Format(bppMap[ji],6,0), ". dAIC =  ", Format(nullAIC-myAIC,10,2), " dAICc = ",Format(nullAICc-myAICc,10,2)," dBIC = ",Format(nullBIC-myBIC,10,2), );	
	
		MasterList [ji-1][0] = ji;
		MasterList [ji-1][1] = myAIC;
		MasterList [ji-1][2] = myAICc;
		MasterList [ji-1][3] = myBIC;
		if (jobKind == 1)
		{
			ConstructCategoryMatrix (siteLikelihoods, lf, COMPLETE);		
			runKHResampler (ji,siteLikelihoods);
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
	fprintf (treeFilePath2,bpF2, "\t", Format (givenTree1,0,0), "\n");

	LikelihoodFunction lf  = (filteredData0,givenTree0,filteredData1,givenTree1);
	
		
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
			mpiNode = ReceiveJobs (1,jobIndex,0);
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
			mpiNode = ReceiveJobs (1, jobIndex, khOption);
			if (khOption)
			{
				LikelihoodFunction lf2 = (filteredData0,givenTree1,filteredData1,givenTree0);
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


/*ChoiceList (dataType,"Data type",1,SKIP_NONE,"Nucleotide","Nucleotide data.",
				     "Codon","Codon (several available genetic codes).");

if (dataType<0) 
{
	return;
}*/

dataType = 0;

if (dataType==0)
{
	SetDialogPrompt 	  		  ("Locate a multi-alignment nucleotide data  file:");
}
else
{
	SetDialogPrompt 	  		  ("Locate a multi-alignment codon data file:");
}

fscanf  (PROMPT_FOR_FILE,"Lines",inFile);

alignmentStride = 11;
alignmentCount	= Columns (inFile)/alignmentStride;
fprintf (stdout, "Read ", alignmentCount, " alignments.\n");

alStr    = 0;
khOption = 0;

resultMatrix = {alignmentCount, 6};

for (alignmentIndex = 0; alignmentIndex < alignmentCount; alignmentIndex = alignmentIndex + 1)
{
	fprintf (stdout, "Working on alignment ", alignmentIndex+1, "\n");
	
	bppMap					= {};
	
	alignmentString = "";
	alignmentString * 1024;
	for (inStr=0; inStr < alignmentStride; inStr = inStr+1)
	{
		alignmentString * (inFile[alStr] + "\n");
		alStr = alStr + 1;
	}
	alignmentString * 0;
	
	DataSet 		ds  		 = ReadFromString (alignmentString);
	DataSetFilter	filteredData = CreateFilter(ds,1);

	siteCount = filteredData.sites;

	InferTreeTopology (0);
	treeString = TreeMatrix2TreeString(0);

	baseSites  = filteredData.sites;

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

	/*ChoiceList (khOption,"KH Testing",1,SKIP_NONE,"Skip","Use only AIC to measure goodness of fit.",
					     "Run","Verify conflicting phylogenetic signal with KH resampling.");


	if (khOption < 0)
	{
		return 0;
	}

	if (khOption == 1)
	{
		itCount = 0;
		while (itCount < 1)
		{
			fprintf (stdout, "How many KH samples should be drawn per breakpoint ?");
			fscanf (stdin, "Number", itCount);
		}
		SHORT_MPI_RETURN = 0;	
	}*/

	if (alignmentIndex == 0)
	{
		SelectTemplateModel (filteredData);	
	}
	else
	{
		HarvestFrequencies (vectorOfFrequencies,filteredData,1,1,1);
		mbf = PopulateModelMatrix ("myMatrix",vectorOfFrequencies);
		Model myModel = (myMatrix,vectorOfFrequencies,mbf);
	}

	/*SetDialogPrompt ("Save results to:");*/

	/*fprintf (PROMPT_FOR_FILE,CLEAR_FILE);

	resFilePath = LAST_FILE_PATH;
	treeFilePath1 = resFilePath + ".trees1";
	treeFilePath2 = resFilePath + ".trees2";

	fprintf (treeFilePath1,CLEAR_FILE);
	fprintf (treeFilePath2,CLEAR_FILE);
	if (khOption)
	{
		khFilePath = resFilePath + ".kh";
		fprintf (khFilePath, CLEAR_FILE, "Breakpoint, Part 1 KH p-value, Part 2 KH p-value\n");
	}*/

	Tree givenTree = treeString;

	branchNames	  = BranchName (givenTree,-1);
	LikelihoodFunction lf = (filteredData,givenTree);
	Optimize (res,lf);
	currentPopulation  = {};
	baseParams = res[1][2] + 3; /* for frequencies */
	
	myDF	 = res[1][1];
	nullAIC  = -2*(res[1][0]-myDF);
	nullAICc = -2*(res[1][0]-myDF*(baseSites/(baseSites-myDF-1)));
	nullBIC	 = -2*(res[1][0]-myDF*Log(baseSites));
	
	MasterList    = {Abs(bppMap)-1,4};
	MatrixList1   = {};
	MatrixList2   = {};
	ResamplesDone = {};

	fprintf (stdout, "\n1). Single partition analysis\n",lf,"\nAIC=",nullAIC,"\n\n2). Looking for a breakpoint...\n");

		
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
		fprintf (stdout, "\n\nBest supported breakpoing is located at position ", bppMap[individual], "\nAIC = ", sortedScores[0][1], " : an imrovement of ", 
					     nullAIC - sortedScores[0][1], " AIC points\n");
		resultMatrix [alignmentIndex][0] = bppMap[individual];
		resultMatrix [alignmentIndex][1] = nullAIC - sortedScores[0][1];
	}
	else
	{
		fprintf (stdout, "\n\nThere seems to be NO recombination in this alignment\n\n");
		resultMatrix [alignmentIndex][0] = -1;
		resultMatrix [alignmentIndex][1] = 0;
	}

	sortedScores = MasterList % 2;
	
	fprintf (stdout, "\nAIC-c");
	
	if (sortedScores[0][2] < nullAICc)
	{
		individual = sortedScores[0][0];
		fprintf (stdout, "\n\nBest supported breakpoing is located at position ", bppMap[individual], "\nAIC = ", sortedScores[0][2], " : an imrovement of ", 
					     nullAICc - sortedScores[0][2], " AIC points\n");
		resultMatrix [alignmentIndex][2] = bppMap[individual];
		resultMatrix [alignmentIndex][3] = nullAICc - sortedScores[0][1];
	}
	else
	{
		fprintf (stdout, "\n\nThere seems to be NO recombination in this alignment\n\n");
		resultMatrix [alignmentIndex][2] = -1;
		resultMatrix [alignmentIndex][3] = 0;
	}
	
	sortedScores = MasterList % 3;
	
	fprintf (stdout, "\nBIC");
	
	if (sortedScores[0][3] < nullBIC)
	{
		individual = sortedScores[0][0];
		fprintf (stdout, "\n\nBest supported breakpoing is located at position ", bppMap[individual], "\nAIC = ", sortedScores[0][3], " : an imrovement of ", 
					     nullBIC - sortedScores[0][3], " AIC points\n");
		resultMatrix [alignmentIndex][4] = bppMap[individual];
		resultMatrix [alignmentIndex][5] = nullBIC - sortedScores[0][3];
	}
	else
	{
		fprintf (stdout, "\n\nThere seems to be NO recombination in this alignment\n\n");
		resultMatrix [alignmentIndex][4] = -1;
		resultMatrix [alignmentIndex][5] = 0;
	}	
}

columnHeaders = {{"AIC site","dAIC","AICc site","dAICc","BIC site","dBIC",";1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50;51;52;53;54;55;56;57;58;59;60;61;62;63;64;65;66;67;68;69;70;71;72;73;74;75;76;77;78;79;80;81;82;83;84;85;86;87;88;89;90;91;92;93;94;95;96;97;98;99;100"}};
OpenWindow (CHARTWINDOW,{{"Recombination Run"}
		{"columnHeaders"}
		{"resultMatrix"}
		{"Scatterplot"}
		{"Index"}
		{"AIC site;AICc site;BIC site"}
		{""}
		{""}
		{""}
		{"0"}
		{""}
		{"0;0"}
		{"10;1.309;0.785398"}
		{"Times:12:0;Times:10:0;Times:12:2"}
		{"0;0;13816530;16777215;0;0;6579300;11842740;13158600;14474460;0;3947580;16777215;15670812;6845928;16771158;2984993;9199669;7018159;1460610;16748822;11184810;14173291"}
		{"16,0,0"}
		},
		"816;639;70;70");
		
SetDialogPrompt ("Save results to:");
fprintf (PROMPT_FOR_FILE,CLEAR_FILE,resultMatrix);
