#include "binomial.ibf";

/*___________________________________________________________________________________________________________*/

function PadString (padLength,padChar)
{
	padString = "";
	d=padString*padLength;
	for (padCounter=0;padCounter<padLength;padCounter=padCounter+1)
	{
		d=padString*padChar;
	}
	d=padString*0;
	fprintf (stdout,padString);
	return padString;
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

nucCharacters="ACGT";
function codonString (ccode)
{
	return nucCharacters[ccode$16]+nucCharacters[(ccode%16)$4]+nucCharacters[ccode%4];
}

/*___________________________________________________________________________________________________________*/

function	PrintTableToFile (dataMatrix, titleMatrix, promptOrNot)
{
	SetDialogPrompt ("Export tab separated data to:");
	
	if (promptOrNot)
	{
		fprintf (PROMPT_FOR_FILE, CLEAR_FILE);
	}
	
	fprintf (LAST_FILE_PATH, titleMatrix[0][0]);

	for (counter1=1; counter1<Columns(titleMatrix); counter1 = counter1+1)
	{
		fprintf (LAST_FILE_PATH, "\t", titleMatrix[0][counter1]);
	}
	
	for (counter1=0; counter1<Rows(dataMatrix); counter1 = counter1 + 1)
	{
		fprintf (LAST_FILE_PATH,"\n",dataMatrix[counter1][0]);
		for (counter2 = 1; counter2 < Columns (dataMatrix); counter2 = counter2+1)
		{
			fprintf (LAST_FILE_PATH,"\t",dataMatrix[counter1][counter2]);
		}
	}
	
	return 1;
}

/*----------------------------------------------------------------------------*/

if (pipeThroughFlag == 0)
{
	ChoiceList (ambChoice, "Treatment of Ambiguities",1,SKIP_NONE,
				"Averaged","All possible resolutions are considered and averaged.",
				"Resolved","The most frequent (for that site) resolution is chosen.");

	if (ambChoice<0)
	{
		return;
	}

	if (useCustomCountingBias)
	{
		incFileName = "Distances/CodonToolsMain.def";
	}
	else
	{
		incFileName = "Distances/CodonTools.def";
	}

	ExecuteCommands  ("#include \""+incFileName+"\";");
}



DataSet dsA					= ReconstructAncestors (lf);
DataSetFilter filteredDataA = CreateFilter (dsA,3,"","",GeneticCodeExclusions);

HarvestFrequencies			  (observedCEFV,filteredData,3,3,0);
seqToBranchMap = {stateCharCount,1};

DataSet		   dsJoint = Combine(dsA,ds);
DataSetFilter filteredDataJ = CreateFilter (dsJoint,3,"","",GeneticCodeExclusions);

hShift = 0;

for (k=0; k<64; k=k+1)
{
	if (_Genetic_Code[k]==10)
	{
		hShift = hShift+1;
	}
	else
	{
		seqToBranchMap[k-hShift] = observedCEFV[k];
	}
}

observedCEFV = seqToBranchMap;


branchNames = BranchName (givenTree,-1);
h = Columns (branchNames);

seqToBranchMap 	= {h, 2};
/* maps sequence names to branch order in column 1 
   and the other way around in column 2 */

for (k=0; k<filteredData.species; k=k+1)
{
	GetString (seqName, filteredData, k);
	seqToBranchMap[k][0] = -1;
	for (v=0; v<h; v=v+1)
	{
		if (branchNames[v] % seqName)
		{
			seqToBranchMap[k][0] = v;
			seqToBranchMap[v][1] = k;
			break;
		}
	}
}

seqToBranchMap[filteredData.species][0] = h-1;
seqToBranchMap[h-1][1] = filteredData.species;

internalBranchNames = {};

for (k=1; k<filteredDataA.species; k=k+1)
{
	GetString (seqName, filteredDataA, k);
	internalBranchNames[seqName] = 1;
	seqToBranchMap[filteredData.species+k][0] = -1;
	for (v=0; v<h; v=v+1)
	{
		if (branchNames[v] % seqName)
		{
			seqToBranchMap[k+filteredData.species][0] = v;
			seqToBranchMap[v][1] = k+filteredData.species;
			break;
		}
	}
}

/* total tree length */

totalTreeLength = 0;
branchLengths   = BranchLength(givenTree,-1);
terminalBL 	    = 0;
internalBL		= 0;

for (v=Columns(branchLengths)-1; v>=0; v=v-1)
{
	totalTreeLength = totalTreeLength + branchLengths[v];
	bName = branchNames[v];
	if (internalBranchNames[bName])
	{
		internalBL = internalBL + branchLengths[v];
	}
}

terminalBL = totalTreeLength - internalBL;

/* get codon matrix */

codonInfo  = {filteredData.species, filteredData.unique_sites};
codonInfo2 = {filteredDataA.species, filteredDataA.unique_sites};

GetDataInfo    (dupInfo, filteredData);
GetDataInfo	   (dupInfoA, filteredDataA);

matrixTrick  = {1,stateCharCount};
matrixTrick2 = {1,stateCharCount};

for (h=Columns(matrixTrick)-1; h>=0; h=h-1)
{
	matrixTrick  [h] = h;
	matrixTrick2 [h] = 1;
}

for (v=0; v<filteredData.unique_sites;v=v+1)
{
	for (h=0; h<filteredData.species;h=h+1)
	{
		GetDataInfo (siteInfo, filteredData, h, v);
		_SITE_ES_COUNT = matrixTrick2 * siteInfo; 
		if (_SITE_ES_COUNT[0] == 1)
		{
			siteInfo = matrixTrick * siteInfo;
			codonInfo[h][v] = siteInfo[0];
		}
		else
		{
			codonInfo[h][v] = -1;
		}
	}
}

for (v=0; v<filteredDataA.unique_sites;v=v+1)
{
	for (h=0; h<filteredDataA.species;h=h+1)
	{
		GetDataInfo (siteInfo, filteredDataA, h, v);
		siteInfo = matrixTrick * siteInfo;
		codonInfo2[h][v] = siteInfo[0];
	}
}

_SITE_RESULTS = {4,filteredData.sites};
flatTreeRep	  = Abs (givenTree);

resultMatrix  = {filteredData.sites,12};
resultMatrixInternal = {filteredData.sites,12};

for (v=0; v<filteredData.sites;v=v+1)
{
	_SITE_ES_COUNT = {stateCharCount,stateCharCount};
	_SITE_EN_COUNT = {stateCharCount,stateCharCount};
	_SITE_OS_COUNT = {stateCharCount,stateCharCount};
	_SITE_ON_COUNT = {stateCharCount,stateCharCount};
	
	_SITE_ES_COUNTInternals = {stateCharCount,stateCharCount};
	_SITE_EN_COUNTInternals = {stateCharCount,stateCharCount};
	_SITE_OS_COUNTInternals = {stateCharCount,stateCharCount};
	_SITE_ON_COUNTInternals = {stateCharCount,stateCharCount};

	/* do internal nodes first */
	
	k = filteredData.species+1;
	
	/* first sequence is always the root */
	c1 = dupInfoA[v];
	if (codonInfo2[1][c1] >= stateCharCount)
	{
		continue;
	}
	for (h=1; h<filteredDataA.species; h=h+1)
	{
		p1 = seqToBranchMap[k][0];
		p2 = flatTreeRep[p1];
		p2 = seqToBranchMap[p2][1]-filteredData.species;
		
		branchFactor = branchLengths[p1]/internalBL;
		
		cd1 = codonInfo2[h] [c1];
		cd2 = codonInfo2[p2][c1];
		
		if (cd1 >= 0 && cd2 >= 0)
		{
            _SITE_OS_COUNTInternals[cd1][cd2] = _SITE_OS_COUNTInternals[cd1][cd2] + 1;		
            _SITE_ON_COUNTInternals[cd1][cd2] = _SITE_ON_COUNTInternals[cd1][cd2] + 1;		
            _SITE_ES_COUNTInternals[cd1][cd2] = _SITE_ES_COUNTInternals[cd1][cd2] + branchFactor;		
            _SITE_EN_COUNTInternals[cd1][cd2] = _SITE_EN_COUNTInternals[cd1][cd2] + branchFactor;		
		}
        		
		k=k+1;
	}
	
	/* now do the leaves */
	
	observedCEFV = {{0}};
	
	c2 = dupInfo[v];
	for (h=0; h<filteredData.species; h=h+1)
	{
		p1 = seqToBranchMap[h][0];
		p2 = flatTreeRep[p1];
		p2 = seqToBranchMap[p2][1]-filteredData.species;
		
		cd2 = codonInfo2[p2][c1];
		cd1 = codonInfo[h] [c2];
		
		branchFactor = branchLengths[p1]/terminalBL;

		if (cd1>=0)
		/* no ambiguities */
		{
			_SITE_OS_COUNT[cd1][cd2] = _SITE_OS_COUNT[cd1][cd2] + 1;		
			_SITE_ON_COUNT[cd1][cd2] = _SITE_ON_COUNT[cd1][cd2] + 1;		
			_SITE_ES_COUNT[cd1][cd2] = _SITE_ES_COUNT[cd1][cd2] + branchFactor;		
			_SITE_EN_COUNT[cd1][cd2] = _SITE_EN_COUNT[cd1][cd2] + branchFactor;		

		}	
		else
		/* ambiguities here */
		{
			GetDataInfo    (ambInfo, filteredData, h, c2);	
			if (Rows(observedCEFV) == 1)
			{
				siteFilter = ""+(v*3)+"-"+(v*3+2);
				DataSetFilter filteredDataSite = CreateFilter (dsJoint,3,siteFilter,"",GeneticCodeExclusions);
				HarvestFrequencies			  (observedCEFV,filteredDataSite,3,3,0);
				tempMx = {stateCharCount,1};

				hShift = 0;

				for (k=0; k<64; k=k+1)
				{
					if (_Genetic_Code[k]==10)
					{
						hShift = hShift+1;
					}
					else
					{
						tempMx[k-hShift] = observedCEFV[k];
					}
				}	
				observedCEFV = tempMx;		
			}
			
			weightFactor = matrixTrick2*ambInfo;
			if (weightFactor[0]<stateCharCount)
			{
				ambInfo  	 = ambInfo$observedCEFV;
				
				if (ambChoice)
				{
					weightFactor = 0;
					tempMx = -1;
					for (k=0; k<stateCharCount; k=k+1)
					{
						if (ambInfo[k]>weightFactor)
						{
							weightFactor = ambInfo[k];
							tempMx = k;
						}
					}
					if (tempMx>=0)
					{
						_SITE_OS_COUNT[tempMx][cd2] = _SITE_OS_COUNT[tempMx][cd2] + 1;		
						_SITE_ON_COUNT[tempMx][cd2] = _SITE_ON_COUNT[tempMx][cd2] + 1 ;		
						_SITE_ES_COUNT[tempMx][cd2] = _SITE_ES_COUNT[tempMx][cd2] + branchFactor;		
						_SITE_EN_COUNT[tempMx][cd2] = _SITE_EN_COUNT[tempMx][cd2] + branchFactor;		
					}
				}
				else
				{
					weightFactor = matrixTrick2*ambInfo;
					weightFactor = weightFactor[0];
					if (weightFactor)
					{
						ambInfo		 = ambInfo * (1/weightFactor);
						for (k=0; k<stateCharCount; k=k+1)
						{
							weightFactor = ambInfo[k];
							if (weightFactor>0)
							{
								_SITE_OS_COUNT[k][cd2] = _SITE_OS_COUNT[k][cd2] + weightFactor;		
								_SITE_ON_COUNT[k][cd2] = _SITE_ON_COUNT[k][cd2] + weightFactor;		
								_SITE_ES_COUNT[k][cd2] = _SITE_ES_COUNT[k][cd2] + weightFactor*branchFactor;		
								_SITE_EN_COUNT[k][cd2] = _SITE_EN_COUNT[k][cd2] + weightFactor*branchFactor;
							}
						}
					}
				}
			}
		}
	}
	
	_SITE_OS_COUNT = matrixTrick2*(_OBSERVED_S_$_SITE_OS_COUNT)*Transpose(matrixTrick2);
	_SITE_ON_COUNT = matrixTrick2*(_OBSERVED_NS_$_SITE_ON_COUNT)*Transpose(matrixTrick2);
	_SITE_ES_COUNT = matrixTrick2*(_PAIRWISE_S_$_SITE_ES_COUNT)*Transpose(matrixTrick2);
	_SITE_EN_COUNT = matrixTrick2*(_PAIRWISE_NS_$_SITE_EN_COUNT)*Transpose(matrixTrick2);
	
	resultMatrix[v][0] = _SITE_OS_COUNT[0];
	resultMatrix[v][1] = _SITE_ON_COUNT[0];
	resultMatrix[v][2] = _SITE_ES_COUNT[0];
	resultMatrix[v][3] = _SITE_EN_COUNT[0];
	
	weight = _SITE_EN_COUNT[0]+_SITE_ES_COUNT[0];
	
	p = _SITE_ES_COUNT[0]/weight;
	
	resultMatrix[v][5] = p;
	
	p2 = resultMatrix[v][0]+resultMatrix[v][1];
	
	resultMatrix[v][7] = resultMatrix[v][1]/resultMatrix[v][3];
	
	if (resultMatrix[v][2])
	{
		resultMatrix[v][6]  = resultMatrix[v][0]/resultMatrix[v][2];
		resultMatrix[v][8]  = resultMatrix[v][7]-resultMatrix[v][6];	
		resultMatrix[v][11] = resultMatrix[v][8]/totalTreeLength;	
	}
	
	if (p2)
	{
		resultMatrix[v][4]  = resultMatrix[v][0]/p2;		
		resultMatrix[v][9]  = extendedBinTail (p2,p,resultMatrix[v][0]);
		if (resultMatrix[v][0]>=1)
		{
			resultMatrix[v][10] = 1-extendedBinTail(p2,p,resultMatrix[v][0]-1);
		}
		else
		{
			resultMatrix[v][10] = 1-extendedBinTail (p2,p,0);
		}
	}		

	_SITE_OS_COUNTInternals = matrixTrick2*(_OBSERVED_S_$_SITE_OS_COUNTInternals)*Transpose(matrixTrick2);
	_SITE_ON_COUNTInternals = matrixTrick2*(_OBSERVED_NS_$_SITE_ON_COUNTInternals)*Transpose(matrixTrick2);
	_SITE_ES_COUNTInternals = matrixTrick2*(_PAIRWISE_S_$_SITE_ES_COUNTInternals)*Transpose(matrixTrick2);
	_SITE_EN_COUNTInternals = matrixTrick2*(_PAIRWISE_NS_$_SITE_EN_COUNTInternals)*Transpose(matrixTrick2);
	
	resultMatrixInternal[v][0] = _SITE_OS_COUNTInternals[0];
	resultMatrixInternal[v][1] = _SITE_ON_COUNTInternals[0];
	resultMatrixInternal[v][2] = _SITE_ES_COUNTInternals[0];
	resultMatrixInternal[v][3] = _SITE_EN_COUNTInternals[0];
	
	weight = _SITE_EN_COUNTInternals[0]+_SITE_ES_COUNTInternals[0];
	
	p = _SITE_ES_COUNTInternals[0]/weight;
	
	resultMatrixInternal[v][5] = p;
	
	p2 = resultMatrixInternal[v][0]+resultMatrixInternal[v][1];
	
	resultMatrixInternal[v][7] = resultMatrixInternal[v][1]/resultMatrixInternal[v][3];
	
	if (resultMatrixInternal[v][2])
	{
		resultMatrixInternal[v][6]  = resultMatrixInternal[v][0]/resultMatrixInternal[v][2];
		resultMatrixInternal[v][8]  = resultMatrixInternal[v][7]-resultMatrixInternal[v][6];	
		resultMatrixInternal[v][11] = resultMatrixInternal[v][8]/totalTreeLength;	
	}
	
	if (p2)
	{
		resultMatrixInternal[v][4]  = resultMatrixInternal[v][0]/p2;		
		resultMatrixInternal[v][9]  = extendedBinTail (p2,p,resultMatrixInternal[v][0]);
		if (resultMatrixInternal[v][0]>=1)
		{
			resultMatrixInternal[v][10] = 1-extendedBinTail(p2,p,resultMatrixInternal[v][0]-1);
		}
		else
		{
			resultMatrixInternal[v][10] = 1-extendedBinTail (p2,p,0);
		}
	}		
}

labelMatrix =     {1,12};
labelMatrix[0] = "Observed S Changes";
labelMatrix[1] = "Observed NS Changes";
labelMatrix[2] = "E[S Sites]";
labelMatrix[3] = "E[NS Sites]";
labelMatrix[4] = "Observed S. Prop.";
labelMatrix[5] = "P{S}";
labelMatrix[6] = "dS";
labelMatrix[7] = "dN";
labelMatrix[8] = "dN-dS";
labelMatrix[9]  = "P{S leq. observed}";
labelMatrix[10] = "P{S geq. observed}";
labelMatrix[11] = "Scaled dN-dS";

sigLevel = -1;

while ((sigLevel<=0)||(sigLevel>=1))
{
	fprintf (stdout, "\nSignificance level for a site to be classified as positively/negatively selected?");
	fscanf  (stdin, "Number", sigLevel);
}

posSelected = 0;
negSelected = 0;

p = Rows(resultMatrix);

for (p2=0; p2<p; p2=p2+1)
{
	v = resultMatrix [p2][8];
	if (v>0)
	{
		if (resultMatrix [p2][9] < sigLevel)
		{
			posSelected = posSelected+1;
		}
	}
	else
	{
		if (v<0)
		{
			if (resultMatrix [p2][10] < sigLevel)
			{
				negSelected = negSelected+1;
			}
		}
	}
}

selLabelMatrix = {{"Index","Site Index","dN-dS","p-value"}};

if (posSelected)
{
	psMatrix = {posSelected, 3};
	h = 0;
	for (p2=0; p2<p; p2=p2+1)
	{
		v = resultMatrix [p2][8];
		if (v>0)
		{
			if (resultMatrix [p2][9] < sigLevel)
			{
				psMatrix[h][0] = p2+1;
				psMatrix[h][1] = v;
				psMatrix[h][2] = resultMatrix [p2][9];
				h = h+1;
			}
		}
	}
	
	fprintf (stdout,"\n******* FOUND ", posSelected, " POSITIVELY SELECTED SITES ALONG TERMINAL BRANCHES ********\n\n");
	dummy = PrintASCIITable  (psMatrix, selLabelMatrix);
}
else
{
	fprintf (stdout,"\n******* FOUND NO POSITIVELY SELECTED SITES ALONG TERMINAL BRANCHES ********\n\n");
}

if (negSelected)
{
	psMatrix = {negSelected, 3};
	h = 0;
	for (p2=0; p2<p; p2=p2+1)
	{
		v = resultMatrix [p2][8];
		if (v<0)
		{
			if (resultMatrix [p2][10] < sigLevel)
			{
				psMatrix[h][0] = p2+1;
				psMatrix[h][1] = v;
				psMatrix[h][2] = resultMatrix [p2][10];
				h = h+1;
			}
		}
	}
	
	fprintf (stdout,"\n******* FOUND ", negSelected, " NEGATIVELY SELECTED SITES ALONG TERMINAL BRANCHES ********\n\n");
	dummy = PrintASCIITable  (psMatrix, selLabelMatrix);
}
else
{
	fprintf (stdout,"\n******* FOUND NO NEGATIVELY SELECTED SITES ALONG TERMINAL BRANCHES ********\n\n");
}

ChoiceList (outputChoice, "Output Options",1,SKIP_NONE,
			"ASCII Table",	 "Output is printed to the console as an ASCII table.",
			"Export to File","Output is spooled to a tab separated file.",
			"Chart","A HYPHY chart window is displayed (Mac/Doze Only).");
			
if (outputChoice<0)
{
	return;
}

if (outputChoice==0)
{
	dummy = PrintASCIITable  (resultMatrix, labelMatrix);
}
else
{
	if (outputChoice == 1)
	{
		SetDialogPrompt ("Save result table to:");
		dummy = PrintTableToFile  (resultMatrix, labelMatrix, !pipeThroughFlag);
	}
	else
	{
		OpenWindow (CHARTWINDOW,{{"SLAC Results For Terminal Branches"}
								   {"labelMatrix"},
								   {"resultMatrix"},
								   {"Contrast Bars"},
								   {"Index"},
								   {labelMatrix[6]+";"+labelMatrix[7]},
								   {"Site Index"},
								   {"dN"},
								   {"dS"},
								   {"0"}},
								   "SCREEN_WIDTH-60;SCREEN_HEIGHT-50;30;50");
	}
}

posSelected = 0;
negSelected = 0;

p = Rows(resultMatrixInternal);

for (p2=0; p2<p; p2=p2+1)
{
	v = resultMatrixInternal [p2][8];
	if (v>0)
	{
		if (resultMatrixInternal [p2][9] < sigLevel)
		{
			posSelected = posSelected+1;
		}
	}
	else
	{
		if (v<0)
		{
			if (resultMatrixInternal [p2][10] < sigLevel)
			{
				negSelected = negSelected+1;
			}
		}
	}
}

selLabelMatrix = {{"Index","Site Index","dN-dS","p-value"}};

if (posSelected)
{
	psMatrix = {posSelected, 3};
	h = 0;
	for (p2=0; p2<p; p2=p2+1)
	{
		v = resultMatrixInternal [p2][8];
		if (v>0)
		{
			if (resultMatrixInternal [p2][9] < sigLevel)
			{
				psMatrix[h][0] = p2+1;
				psMatrix[h][1] = v;
				psMatrix[h][2] = resultMatrixInternal [p2][9];
				h = h+1;
			}
		}
	}
	
	fprintf (stdout,"\n******* FOUND ", posSelected, " POSITIVELY SELECTED SITES ALONG INTERNAL BRANCHES ********\n\n");
	dummy = PrintASCIITable  (psMatrix, selLabelMatrix);
}
else
{
	fprintf (stdout,"\n******* FOUND NO POSITIVELY SELECTED SITES ALONG INTERNAL BRANCHES ********\n\n");
}

if (negSelected)
{
	psMatrix = {negSelected, 3};
	h = 0;
	for (p2=0; p2<p; p2=p2+1)
	{
		v = resultMatrixInternal [p2][8];
		if (v<0)
		{
			if (resultMatrixInternal [p2][10] < sigLevel)
			{
				psMatrix[h][0] = p2+1;
				psMatrix[h][1] = v;
				psMatrix[h][2] = resultMatrixInternal [p2][10];
				h = h+1;
			}
		}
	}
	
	fprintf (stdout,"\n******* FOUND ", negSelected, " NEGATIVELY SELECTED SITES ALONG INTERNAL BRANCHES ********\n\n");
	dummy = PrintASCIITable  (psMatrix, selLabelMatrix);
}
else
{
	fprintf (stdout,"\n******* FOUND NO NEGATIVELY SELECTED SITES ALONG INTERNAL BRANCHES ********\n\n");
}

if (outputChoice==0)
{
	dummy = PrintASCIITable  (resultMatrixInternal, labelMatrix);
}
else
{
	if (outputChoice == 1)
	{
		SetDialogPrompt ("Save results for internal branches to:");
		dummy = PrintTableToFile  (resultMatrixInternal, labelMatrix, !pipeThroughFlag);
	}
	else
	{
		OpenWindow (CHARTWINDOW,{{"SLAC Results For Internal Branches"}
								   {"labelMatrix"},
								   {"resultMatrixInternal"},
								   {"Contrast Bars"},
								   {"Index"},
								   {labelMatrix[6]+";"+labelMatrix[7]},
								   {"Site Index"},
								   {"dN"},
								   {"dS"},
								   {"0"}},
								   "SCREEN_WIDTH-60;SCREEN_HEIGHT-50;30;50");
	}
}
