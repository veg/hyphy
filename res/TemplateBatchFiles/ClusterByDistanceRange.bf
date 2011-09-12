ExecuteAFile ("Utility/GrabBag.bf");
ExecuteAFile ("SequentialAddition.ibf");


/*___________________________________________________________________________________________________________*/

function FindMaxDistance (rIndex, fIndex)
{
	max 	= 0;
	oldIdx 	= 0;
	for (k=0; k<rIndex; k=k+1)
	{
		t = distanceMatrix[rIndex][k];
		if (t>max)
		{
			oldIdx = min;
			max = t;
			min = k;
		}
	}
	for (k=rIndex+1; k<mDim; k=k+1)
	{
		t = distanceMatrix[rIndex][k];
		if (t>max)
		{
			oldIdx = min;
			max = t;
			min = k;
		}
	}
	if (min == fIndex)
	{
		return oldIdx;
	}
	return min;
}
/*___________________________________________________________________________________________________________*/

function PadString (padLength,padChar)
{
	for (padCounter=0;padCounter<padLength;padCounter=padCounter+1)
	{
		fprintf (stdout,padChar);
	}
	return padLength;
}

/*___________________________________________________________________________________________________________*/

function	PrintASCIITable (dataMatrix, titleMatrix)
{
	columnWidths = {1,Columns(titleMatrix)};
	for (counter1=0; counter1<Columns(titleMatrix); counter1 = counter1+1)
	{
		counter2 = Abs (titleMatrix[0][counter1])+2;
		if (counter2<12)
		{
			counter2 = 12;
		}
		columnWidths[0][counter1] = counter2;
	}
	fprintf (stdout, "\n");
	for (counter2 = 0; counter2 < Columns (titleMatrix); counter2 = counter2+1)
	{
		fprintf (stdout,"+-");
		dummy = PadString (columnWidths[0][counter2],"-");
		fprintf (stdout,"-");
	}
	fprintf (stdout,"+\n| ");
	
	for (counter1=0; counter1<Columns(titleMatrix); counter1 = counter1+1)
	{
		fprintf (stdout, titleMatrix[counter1]);
		dummy = PadString (columnWidths[0][counter1]-Abs(titleMatrix[counter1])," ");
		fprintf (stdout, " | ");
	}
	
	fprintf (stdout, "\n");
	
	for (counter1=-1; counter1<Rows(dataMatrix); counter1 = counter1 + 1)
	{
		if (counter1>=0)
		{
			fprintf (stdout,"| ");
			fprintf (stdout,Format(counter1+1,columnWidths[0][0],0));
			for (counter2 = 1; counter2 < Columns (titleMatrix); counter2 = counter2+1)
			{
				fprintf (stdout," | ");
				fprintf (stdout,Format(dataMatrix[counter1][counter2-1],columnWidths[0][counter2],-1));
			}
			fprintf (stdout," ");
			fprintf (stdout, "|\n");
		}
		for (counter2 = 0; counter2 < Columns (titleMatrix); counter2 = counter2+1)
		{
			fprintf (stdout,"+-");
			dummy = PadString (columnWidths[0][counter2],"-");
			fprintf (stdout,"-");
		}
		fprintf (stdout, "+\n");
	}
	
	return 1;
}

/*___________________________________________________________________________________________________________*/

ChoiceList (dataType,"Data type",1,SKIP_NONE,"Nucleotide/Protein","Nucleotide or amino-acid (protein).",
				     "Codon","Codon (several available genetic codes).");

if (dataType<0) 
{
	return;
}
if (dataType)
{
	SetDialogPrompt ("Please choose a codon data file:");
	ExecuteAFile("TemplateModels/chooseGeneticCode.def");
}
else
{
	SetDialogPrompt ("Please choose a nucleotide or amino-acid data file:");
}


DataSet ds = ReadDataFile (PROMPT_FOR_FILE);

if (dataType)
{
	DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);
}
else
{
	DataSetFilter filteredData = CreateFilter (ds,1,"","");
}

SetDialogPrompt 	("Load a distance matrix:");
fscanf 				(PROMPT_FOR_FILE,"NMatrix",distanceMatrix);


/* compute some matrix statistics */

sum  = 0;
sum2 = 0;
max	 = 0;
min  = 1e100;

mDim = Rows (distanceMatrix);

if (mDim != ds.species)
{
	fprintf (stdout, "[ERROR] Read ", ds.species, " but the distance matrix had an incompatible dimension: ", mDim, "\n");
	return 0;
}

dDim = mDim*(mDim-1)/2;

ChoiceList (doHist,"Histogram",1,SKIP_NONE,"Generate","Compute and report distribution statistics of pairwise distances.",
				     "Skip","Skip this step");
	
if (doHist<0)
{
	return 0;
}

if (doHist == 0)
{
	vec  = {dDim,1};

	k	  = 0;
	for (r=0; r< mDim; r=r+1)
	{
		for (c=r+1; c<mDim; c=c+1)
		{
			t                    = Max(distanceMatrix[r][c],0);
			distanceMatrix[c][r] = t;
			distanceMatrix[r][c] = t;
			vec[k]				 = t;
			k					 = k+1;
			sum				     = sum+t;
		}
	}

	vec		 = vec%0;
	min		 = vec[0];
	max		 = vec[dDim-1];
	sum      = sum/dDim;


	IQR = vec[3*dDim$4] - vec[dDim$4];

	h    = 2*IQR/dDim^(1/3);
	b25	 = vec[25*dDim$1000];
	t975 = vec[975*dDim$1000];
	span = t975 - b25;
	bins = (span/h+0.5)$1;
	h	 = span/bins;


	histogramTable		 = {bins,4};
	histogramTable[0][0] = min;
	histogramTable[0][1] = b25+h;

	for (c=1; c<bins;c=c+1)
	{
		histogramTable[c][0] = histogramTable[c-1][1];
		histogramTable[c][1] = histogramTable[c-1][1] + h;
	}

	histogramTable[c-1][1] = max;
	bm1				       = bins-1;

	for (r=0; r< mDim; r=r+1)
	{
		for (c=r+1; c<mDim; c=c+1)
		{
			t     = distanceMatrix[r][c];
			sum2  = sum2+(t-sum)^2;
			t	  = Max(0,Min(((t-b25)/h$1),bm1));
			
			
			histogramTable[t][2] = histogramTable[t][2] + 1;
		}
	}

	for (c=0; c<25;c=c+1)
	{
		histogramTable[c][3] = histogramTable[c][2]/dDim;
	}

	fprintf (stdout, "\nDistance matrix info:\n",
					 "Mean        : ", sum, "\n",
					 "Variance    : ", sum2/(dDim-1), "\n",
					 "IQR         : ", IQR, "\n",
					 "Minimum     : ", min, "\n",
					 "Maximum     : ", max, "\n");

	labels = {{"Bin","From","To", "Count", "Proportion"}};

	PrintASCIITable (histogramTable, labels);
}
else
{
	/*for (r=0; r< mDim; r=r+1)
	{
		for (c=r+1; c<mDim; c=c+1)
		{
			t = Max(distanceMatrix[r][c], 0);
			distanceMatrix[c][r] = t;
			distanceMatrix[r][c] = t;
		}
	}*/	
	
	max = Max(distanceMatrix,0);

}

clumpingL = -1;


ChoiceList (clusteringType,"Clustering algorithm",1,SKIP_NONE,
						   "At least one","A sequence in a cluster is sufficiently close to AT LEAST ONE of the other sequences in the cluster.",
						   "All",		  "A sequence in a cluster is sufficiently close to ALL OTHER sequences in the cluster.");
						   
if (clusteringType < 0)
{
		return 0;
}	
		
while (clumpingL < 0)
{
	fprintf (stdout, "Please enter the lower distance bound (>=0, or -1 to plot clusters vs cutoff):");
	fscanf  (stdin,"Number",clumpingL);
	if (clumpingL == (-1))
	{
		break;
	}
}

if (clumpingL>=0)
{
	do
	{
		SetDialogPrompt ("Set comma separated results to:");
		DEFAULT_FILE_SAVE_NAME = "Clustering.csv";	
		fprintf (PROMPT_FOR_FILE,CLEAR_FILE);	
		outFile = LAST_FILE_PATH;
		outString = "";
		outString * 8192;
		clumpingU	 = prompt_for_a_value ("Please enter the upper distance bound", clumpingL + 0.01, clumpingL, max, 0);
		clusterCount = doClustering (clumpingL, clumpingU, clusteringType);
		
		
		fprintf (stdout, "\nFound ", clusterCount, " clusters on total ", totalEdgeCount," edges\n");
		
		sortedCluster = {mDim,2};
		degreeDistro  = {};
		ones		  = {mDim,1}["1"];
		maxDegree	  = 0;
		
		for (k=0; k<mDim; k=k+1)
		{
			sortedCluster[k][0] = visited[k];
			sortedCluster[k][1] = k;
			thisRow 			= ((gatedDistances[k][-1])*ones)[0] + 1;
			degreeDistro[thisRow] = degreeDistro[thisRow] + 1;
		}	
		
		diffDegrees = Abs(degreeDistro);
		logLog 		= {diffDegrees,2};
		keys		= Rows (degreeDistro);
		
		for (k=0; k<diffDegrees; k=k+1)
		{
			logLog[k][0] = 0+keys[k];
			logLog[k][1] = degreeDistro[keys[k]];
		}
	
		
		logLog = logLog % 0;
		logLog[diffDegrees-1][1] = logLog[diffDegrees-1][1] / mDim;
		for (k=diffDegrees-2; k>=0; k=k-1)
		{
			logLog[k][1] = logLog[k][1]/mDim+logLog[k+1][1];
		}
		
		logLog = Log(logLog);
		
columnHeaders = {{"Log[Degree]","Log[Prob]"}};
OpenWindow (CHARTWINDOW,{{"Log-Log degree plot"}
		{"columnHeaders"}
		{"logLog"}
		{"Scatterplot"}
		{"Log[Degree]"}
		{"Log[Prob]"}
		{"Log(Degree)"}
		{""}
		{"Log(Probability)"}
		{"0"}
		{""}
		{"-1;-1"}
		{"10;1.309;0.785398"}
		{"Times:12:0;Times:10:0;Times:12:2"}
		{"0;0;16777215;8355711;0;0;6579300;11842740;13158600;14474460;0;3947580;16777215;15670812;6845928;16771158;2984993;9199669;7018159;1460610;16748822;11184810;14173291"}
		{"16,0,0"}
		},
		"480;665;70;70");
		
		sortedCluster = sortedCluster % 0;
		
		stringLabel = "";
		graphVizSt  = ""; graphVizSt * 128;
		
		edgesPrinted = 0;
		nodesMade	 = {};
		
		outString  * ("SequenceID,ClusterID_"+clumpingL+"_"+clumpingU);
		
		lastClusterID     = 1;
		clusterSpan	      = 0;
		clusterSizes      = {};
		clusterMembership = {};
		definedNode = {};
		
		for (k=0; k<=mDim; k=k+1)
		{
			if (k < mDim)
			{
				clusterID		= sortedCluster[k][0];
				visited [k]		= clusterID;
				GetString		(seqName, ds,sortedCluster[k][1]);
				stringLabel		= stringLabel + ";" + seqName;
				outString	    * ("\n"+seqName+","+clusterID);
			}
			else
			{
				clusterID 		= lastClusterID+1;
			}	
			
			if (lastClusterID != clusterID)
			{
				clusterSize						= k-clusterSpan;
				clusterSizes[clusterSize]		= clusterSizes[clusterSize] + 1;
				clusterMembership[clusterSize]  = clusterMembership[clusterSize] + clusterSize;
				
				graphVizSt * ("\n\n/* cluster " + lastClusterID + "*/\n\n");
				
				
				for (n1 = 0; n1 < clusterSize; n1=n1+1)
				{
					id1  = sortedCluster[clusterSpan + n1][1];
					GetString		(seqName, ds,id1); 
										
					for (n2 = n1 + 1; n2 < clusterSize; n2 = n2+1)
					{
						id2 = sortedCluster[clusterSpan + n2][1];
						d = distanceMatrix[id1][id2];
						/*fprintf (stdout, id1, ":", id2, ":", d, "\n");*/
						if (d <= clumpingU && d>= clumpingL)
						{
							if (Abs(definedNode[seqName]) == 0)
							{
								graphVizSt * ("\"" + seqName + "\" [label = \"\"];\n");
								definedNode[seqName]  = 1;
							}
							GetString		(seqName2, ds,id2); 
							if (Abs(definedNode[seqName2]) == 0)
							{
								graphVizSt * ("\"" + seqName2 + "\" [label = \"\"];\n");
								definedNode[seqName2]  = 1;
							}
							graphVizSt * ("\"" + seqName + "\" -- \"" + seqName2 + "\";\n");
							edgesPrinted = edgesPrinted + 1;
							nodesMade [seqName] = 1;
							nodesMade [seqName2] = 1;
						}
					}
				}
				clusterSpan = k;
				lastClusterID = clusterID;
			}
		}
		
		outString * 0;
		fprintf (outFile, outString);
		fprintf (stdout, "Graphviz edges/nodes = ", edgesPrinted, "/", Abs(nodesMade), "\n");
		
		fprintf (stdout, "\nCluster size distribution\n");
		_printAnAVLNumericTotal (clusterSizes, ".");

		fprintf (stdout, "\nCluster membership distribution\n");
		_printAnAVLNumericTotal (clusterMembership, ".");
		
		graphVizSt * 0;
		
		SetDialogPrompt ("Save GraphViz file to:");
		fprintf			(PROMPT_FOR_FILE, CLEAR_FILE, "graph G{\n", graphVizSt, "\n};");
		
		fprintf (stdout, "Continue with another threshold (y/n)?");
		fscanf  (stdin,"String", shouldCont);
		
	}
	while (shouldCont[0] != "n" && shouldCont[0] != "N");
	
	columnHeaders = {{"Cluster ID",  stringLabel}};

	OpenWindow (CHARTWINDOW,{{"Cluster Allocation"}
		{"columnHeaders"}
		{"visited"}
		{"None"}
		{"Cutoff"}
		{"Clusters"}
		{"Cutoff"}
		{""}
		{"Cluster Count"}
		{"0"}
		{""}
		{"-1;-1"}
		{"10;1.309;0.785398"}
		{"Times:14:0;Times:10:0;Times:14:1"}
		{"0;0;16777215;8421504;0;0;6579300;11842740;13158600;14474460;0;3947580;16777215;0;6845928;16771158;2984993;9199669;7018159;1460610;16748822;11184810;14173291"}
		{"16,0,0"}
		},
		"SCREEN_WIDTH-100;SCREEN_HEIGHT-100;50;50");
}
else
{
	SetDialogPrompt ("Set comma separated results to:");

	DEFAULT_FILE_SAVE_NAME = "Clustering.csv";

	fprintf (PROMPT_FOR_FILE,CLEAR_FILE);

	outFile = LAST_FILE_PATH;
	outString = "";
	outString * 8192;

	maxSteps = 500;
	
	distanceMx = {maxSteps,2} ["(_MATRIX_ELEMENT_COLUMN_==0)*(_MATRIX_ELEMENT_ROW_)*max/maxSteps"];
	
	/*size 	 = mDim*(mDim-1)/2;
	sizeStep = size$maxSteps;
	
	if (sizeStep<1)
	{
		sizeStep = 1;
	}

	size = Min(maxSteps,size);
	
	distanceMx = {size,2};
	
	k  = 0;
	k2 = 0;
	for (r=0; r< mDim; r=r+1)
	{
		for (c=r+1; c<mDim; c=c+1)
		{
			k2 = k2+1;
			if (k2%sizeStep == 0)
			{
				distanceMx [k][0] = distanceMatrix[r][c];
				k = k+1;
			}
		}
	}*/
	
	lastTime = Time (0);
	outString * "Cutoff, Cluster Count";

	distanceMx = distanceMx%0;
	distanceMx[0][1] = doClustering (0, distanceMx[0][0], clusteringType);
	outString * ("\n"+distanceMx[0][0]+","+distanceMx[0][1]);
	for (rz = 1; rz < Rows(distanceMx); rz = rz+1)
	{
		if (Time(0)-lastTime>5)
		{
			fprintf (stdout, "Step ", rz, "/", Rows(distanceMx),"\n");
			lastTime = Time(0);
		}

		if (distanceMx[rz][0]>distanceMx[rz-1][0])
		{
			distanceMx[rz][1] = doClustering (0, distanceMx[rz][0], clusteringType);
		}
		else
		{
			distanceMx[rz][1] = distanceMx[rz-1][1];
		}
		outString * ("\n"+distanceMx[rz][0]+","+distanceMx[rz][1]);
	}
	
	outString * 0;
	fprintf (outFile, outString);
	columnHeaders = {{"Cutoff","Clusters"}};
	OpenWindow (CHARTWINDOW,{{"Cluster Count Plot"}
		{"columnHeaders"}
		{"distanceMx"}
		{"Line Plot"}
		{"Cutoff"}
		{"Clusters"}
		{"Cutoff"}
		{""}
		{"Cluster Count"}
		{"0"}
		{""}
		{"-1;-1"}
		{"10;1.309;0.785398"}
		{"Times:14:0;Times:10:0;Times:14:1"}
		{"0;0;16777215;8421504;0;0;6579300;11842740;13158600;14474460;0;3947580;16777215;0;6845928;16771158;2984993;9199669;7018159;1460610;16748822;11184810;14173291"}
		{"16,0,0"}
		},
		"SCREEN_WIDTH-100;SCREEN_HEIGHT-100;50;50");
}



/*___________________________________________________________________________________________________________*/

function doClustering (from, to, clusteringAlg)
{
	visited 	   = {mDim,1};
	gatedDistances = distanceMatrix["_MATRIX_ELEMENT_VALUE_ >= from && _MATRIX_ELEMENT_VALUE_ <= to && _MATRIX_ELEMENT_ROW_!=_MATRIX_ELEMENT_COLUMN_"];

	clusterCount = 1;
	done		 = 0;
	vertexCount  = 0;
	
	totalEdgeCount = (+gatedDistances) $ 2;
	
	if (clusteringAlg == 0)	
	{
		for (currentVertex=0; currentVertex < mDim; currentVertex=currentVertex+1)
		{
			if (visited[currentVertex] == 0)
			{
				visited[currentVertex] = clusterCount;
				for (k=currentVertex+1; k<mDim; k=k+1)
				{
					if (visited[k] == 0)
					{
						if (gatedDistances[currentVertex][k])
						{
							visited[k] = clusterCount;
							doAVertex (k, 0);
						}
					}
				}
				clusterCount += 1;
			}
		}
	}
	else
	{
		for (currentVertex=0; currentVertex < mDim; currentVertex +=1 )
		{
			//fprintf (stdout, visited, "\n");
			for (k=0; k<currentVertex; k += 1 )
			{
				if (gatedDistances[currentVertex][k])
				{
					tryThisCluster = visited[k];
					//fprintf (stdout, "\n", currentVertex, " testing: ", tryThisCluster, "\n");
					for (k2 = 0; k2 < currentVertex; k2+=1)
					{
						if (visited[k2] == tryThisCluster)
						{
							if (gatedDistances[currentVertex][k2] == 0)
							{
								//fprintf (stdout, currentVertex, " failed: ", k2, "/", distanceMatrix[currentVertex][k2], "\n");
								k2 = -1;
								break;
							}
						}
					}
					if (k2 > 0)
					{
						//fprintf (stdout, currentVertex, " passed: ", tryThisCluster, "\n");
						visited[currentVertex] = tryThisCluster;
						k = -1;
						break;
					}
				}
			}
			if (k >= 0)
			{
				//fprintf (stdout, currentVertex, " assigned to a new cluster: ", clusterCount, "\n");
				visited[currentVertex] = clusterCount;
				clusterCount += 1;			
			}
		}
		
	}
	return clusterCount-1;
}


/*___________________________________________________________________________________________________________*/

function 	doAVertex (vertexID, currentNeighbor)
{
	for (;currentNeighbor < mDim; currentNeighbor=currentNeighbor+1)
	{
		if (visited[currentNeighbor]==0)
		{
			if (gatedDistances[vertexID][currentNeighbor])
			{
				visited[currentNeighbor] = clusterCount;
				ExecuteCommands("doAVertex ("+currentNeighbor+", 0);");
			}
		}
	}
	return 0;
}
