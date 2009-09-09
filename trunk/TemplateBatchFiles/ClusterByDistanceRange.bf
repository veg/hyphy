#include "SequentialAddition.ibf";

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
	#include "TemplateModels/chooseGeneticCode.def";
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
fscanf 				(PROMPT_FOR_FILE,"Matrix",distanceMatrix);

/* compute some matrix statistics */

sum  = 0;
sum2 = 0;
max	 = 0;
min  = 1e100;

mDim = Rows (distanceMatrix);
dDim = mDim*(mDim-1)/2;

for (r=0; r< mDim; r=r+1)
{
	for (c=r+1; c<mDim; c=c+1)
	{
		t = distanceMatrix[r][c];
		if (t < 0)
		{
			distanceMatrix[r][c] = 0;
			distanceMatrix[c][r] = 0;
			t = 0;
		}	
		else
		{
			distanceMatrix[c][r] = distanceMatrix[r][c];
		}
		
		min = Min(min,t);
		max = Max(max,t);
		sum  = sum  + t;
	}
}

sum = sum/dDim;
r   = (max-min)/25;

histogramTable = {25,4};
histogramTable[0][1] = r;

for (c=1; c<25;c=c+1)
{
	histogramTable[c][0] = histogramTable[c-1][1];
	histogramTable[c][1] = histogramTable[c-1][1] + r;
}

histogramTable[c-1][1] = max;

for (r=0; r< mDim; r=r+1)
{
	for (c=r+1; c<mDim; c=c+1)
	{
		t = distanceMatrix[r][c];
		sum2  = sum2+(t-sum)*(t-sum);
		t = t-min;
		k = 0;
		while (histogramTable[k][1]<t)
		{
			k=k+1;
		}
		histogramTable[k][2] = histogramTable[k][2] + 1;
	}
}

for (c=0; c<25;c=c+1)
{
	histogramTable[c][3] = histogramTable[c][2]/dDim;
}

fprintf (stdout, "\nDistance matrix info:\n",
				 "Mean        : ", sum, "\n",
				 "Variance    : ", sum2/(dDim-1), "\n",
				 "Minimum     : ", min, "\n",
				 "Maximum     : ", max, "\n");

labels = {{"Bin","From","To", "Count", "Proportion"}};

dummy = PrintASCIITable (histogramTable, labels);

clumpingL = -1;

SetDialogPrompt ("Set comma separated results to:");

DEFAULT_FILE_SAVE_NAME = "Clustering.csv";

fprintf (PROMPT_FOR_FILE,CLEAR_FILE);

outFile = LAST_FILE_PATH;
outString = "";
outString * 8192;

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
	clumpingU = -1;

	while (clumpingU <= clumpingL)
	{
		fprintf (stdout, "Please enter the upper distance bound (>=",clumpingL,"):");
		fscanf  (stdin,"Number",clumpingU);
	}

	clusterCount = doClustering (clumpingL, clumpingU);
	fprintf (stdout, "\nFound ", clusterCount, " clusters\n");
	
	sortedCluster = {mDim,2};
	
	for (k=0; k<mDim; k=k+1)
	{
		sortedCluster[k][0] = visited[k];
		sortedCluster[k][1] = k;
	}	
	
	sortedCluster = sortedCluster % 0;
	
	stringLabel = "";
	outString * "Sequence, Cluster ID";
	
	for (k=0; k<mDim; k=k+1)
	{
		visited [k] = sortedCluster[k][0];
		GetString (seqName, ds,sortedCluster[k][1]);
		stringLabel = stringLabel + ";" + seqName;
		outString * ("\n"+seqName+","+visited[k]);
	}
	
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
	distanceMx[0][1] = doClustering (0, distanceMx[0][0]);
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
			distanceMx[rz][1] = doClustering (0, distanceMx[rz][0]);
		}
		else
		{
			distanceMx[rz][1] = distanceMx[rz-1][1];
		}
		outString * ("\n"+distanceMx[rz][0]+","+distanceMx[rz][1]);
	}
	
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

outString * 0;

fprintf (outFile, outString);


/*___________________________________________________________________________________________________________*/

function doClustering (from, to)
{
	visited 	   = {mDim,1};
	gatedDistances = distanceMatrix["_MATRIX_ELEMENT_VALUE_ >= from && _MATRIX_ELEMENT_VALUE_ <= to && _MATRIX_ELEMENT_ROW_!=_MATRIX_ELEMENT_COLUMN_"];

	clusterCount = 1;
	done		 = 0;
	vertexCount  = 0;

	for (currentVertex=0; currentVertex < mDim; currentVertex=currentVertex+1)
	{
		if (visited[currentVertex] == 0)
		{
			visited[currentVertex] = clusterCount;
			for (k=currentVertex+1; k<mDim; k=k+1)
			{
				if (visited[k] == 0 && gatedDistances[currentVertex][k])
				{
					visited[k] = clusterCount;
					doAVertex (k, 0);
				}
			}
			clusterCount = clusterCount + 1;
		}
	}
	return clusterCount-1;
}


/*___________________________________________________________________________________________________________*/

function 	doAVertex (vertexID, currentNeighbor)
{
	for (;currentNeighbor < mDim; currentNeighbor=currentNeighbor+1)
	{
		if (visited[currentNeighbor]==0 && gatedDistances[vertexID][currentNeighbor])
		{
			visited[currentNeighbor] = clusterCount;
			if (clusterCount == 20)
			{
				GetString (seqName1, ds, vertexID);
				GetString (seqName2, ds, currentNeighbor);
				fprintf (stdout, seqName1, "=>", seqName2, "\n");
			}
			ExecuteCommands("doAVertex ("+currentNeighbor+", 0);");
		}
	}
	return 0;
}