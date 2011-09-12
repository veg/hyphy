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

SelectTemplateModel	(filteredData);

SetDialogPrompt 	("Load a distance matrix:");

fscanf (PROMPT_FOR_FILE,"Matrix",distanceMatrix);

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
		sum  = sum  + t;
		if (t<min)
		{
			min = t;
		}
		if (t>max)
		{
			max = t;
		}
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

clumping = -1;

while (clumping <= 0)
{
	fprintf (stdout, "Please enter the clustering threshold (>=0):");
	fscanf  (stdin,"Number",clumping);
}

stillInPlay = {mDim,1};

for (r=0; r<mDim; r=r+1)
{
	stillInPlay[r] = 1;
}	

reductionAmount = 0;

SetDialogPrompt ("Save smaller data set to:");

fprintf (PROMPT_FOR_FILE,CLEAR_FILE);
compactedSeqFile = LAST_FILE_PATH;

GetInformation (dataCharacters, filteredData);

for (rowCount=0; rowCount<mDim; rowCount=rowCount+1)
{
	if (stillInPlay[rowCount])
	{
		putativeClump     = {mDim,1};
		putativeClumpSize = 1;
		putativeClump[0]  = rowCount;
		for (c=0; c<mDim; c=c+1)
		{
			if (stillInPlay[c] && (c!=rowCount))
			{
				t = distanceMatrix[rowCount][c];
				if (t <= clumping)
				{
					putativeClump[putativeClumpSize] = c;
					putativeClumpSize = putativeClumpSize+1;
				}
			}
		}
		
		if (putativeClumpSize>1)
		{
			/* check pairwise distances */
			for (c=1; c<putativeClumpSize; c=c+1)
			{
				min = putativeClump[c];
				if (min>=0)
				{
					for (t=c+1;t<putativeClumpSize;t=t+1)
					{
						k = putativeClump[t];
						if (k>=0)
						{
							max = distanceMatrix[k][min];
							if (max>clumping)
							{
								putativeClump[t] = -1;
							}
						}
					}
				}
			}
			
			GetString (sName,ds,putativeClump[0]);
			fprintf (stdout, "\n",sName);
			stillInPlay[rowCount] 	= 0;
			actualClumpSize = 1;
			actualClump 	= {putativeClumpSize,1};
			for (c=1; c<putativeClumpSize; c=c+1)
			{
				min = putativeClump[c];
				if (min>=0)
				{
					GetString (sName,ds,min);
					fprintf (stdout, ",", sName);
					stillInPlay[min] = 0;
					reductionAmount = reductionAmount+1;
					actualClump[actualClumpSize] = min;
					actualClumpSize = actualClumpSize+1;
				}
			}
			fprintf (stdout, "\n");
			nSName = "";
			if (actualClumpSize == 2)
			/* add 2 sequences for the outgroup */
			{
				s1 	 = actualClump[0];
				s2 	 = actualClump[1];
				
				GetString (sName, ds, s1);
				aTreeString = "(("+sName+",";
				nSName = sName;
				GetString (sName, ds, s2);
				nSName = nSName+"_"+sName;
				aTreeString = aTreeString+sName+"),";
				o1  = FindMaxDistance (s1,-1);
				GetString (sName, ds, o1);
				aTreeString = aTreeString+sName+",";
				o2  = FindMaxDistance (s2,o1);
				GetString (sName, ds, o2);
				aTreeString = aTreeString+sName+")";
				pString = ""+s1+(","+s2)+(","+o1)+(","+o2);
				getIndex = 1;
			}
			else
			{
				if (actualClumpSize == 3)
				{
					s1 	 = actualClump[0];
					s2 	 = actualClump[1];
					o1 	 = actualClump[2];
					
					GetString (sName, ds, s1);
					nSName = sName;
					aTreeString = "("+sName+",";
					GetString (sName, ds, s2);
					nSName = nSName+"_"+sName;
					aTreeString = aTreeString+sName+",";
					GetString (sName, ds, o1);
					nSName = nSName+"_"+sName;
					aTreeString = aTreeString+sName+")";
					pString = ""+s1+(","+s2)+(","+o1);
					getIndex = 0;				
				}
				else
				{
					pString = ""+actualClump[0];
					GetString (nSName, ds, actualClump[0]);
					for (k=1; k<actualClumpSize; k=k+1)
					{
						GetString (sName, ds, actualClump[k]);
						nSName = nSName+"_"+sName;
						pString = pString+(","+actualClump[k]);
					}
					if (dataType)
					{
						DataSetFilter filteredData = CreateFilter (ds,3,"",pString,GeneticCodeExclusions);
					}
					else
					{
						DataSetFilter filteredData = CreateFilter (ds,1,"",pString);
					}
					randomOption 	   = 0;
					haveTreeConstraint = 0;
					doNNIOption		   = 1;
					nniPeriod		   = 1;
					methodIndex		   = 0;
					first3Taxa		   = {{0,1,2}};
					k		   		   = InferTreeTopology (0.0);
					aTreeString		   = bestTree;
					getIndex 		   = 0;				
				}
			}
			
			/*fprintf (stdout, "\n", aTreeString, "\n", pString);*/
			Tree aTree = aTreeString;
			if (dataType)
			{
				DataSetFilter filteredData = CreateFilter (ds,3,"",pString,GeneticCodeExclusions);
			}
			else
			{
				DataSetFilter filteredData = CreateFilter (ds,1,"",pString);
			}
			LikelihoodFunction lfClump = (filteredData,aTree);
			DataSet 		   aDS 	   = ReconstructAncestors (lfClump);
			if (dataType)
			{
				DataSetFilter filteredData = CreateFilter (aDS,3,"","",GeneticCodeExclusions);
			}
			else
			{
				DataSetFilter filteredData = CreateFilter (aDS,1,"","");
			}
			GetInformation 	  (aDataCharacters, filteredData);
			GetString (sName,aDS,getIndex);
			fprintf (compactedSeqFile,">",nSName,"\n",aDataCharacters[getIndex],"\n");
		}
		else
		{
			GetString (sName,ds,rowCount);
			fprintf (compactedSeqFile,">",sName,"\n",dataCharacters[rowCount],"\n");
		}
	}
}

fprintf (stdout, "\nTotal sequences pruned: ", reductionAmount,"\n");

distanceMatrix = 0;
dataCharacters = 0;
