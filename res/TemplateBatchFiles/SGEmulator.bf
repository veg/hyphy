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

function mergeSimulatedNulls (nullAVL1, nullAVL2)
{
	nullAVLJoint = {};
	
	ccode = Abs(nullAVL1);
	
	for (counter1 = 0; counter1 < ccode; counter1 = counter1 + 1)
	{
		av1 		= nullAVL1[counter1];
		av2 		= nullAVL2[counter1];
		avJ 		= {};
		doneKeys 	= {};
		keys1		= Rows(av1);
		keys2 		= Rows(av2);
		keyC		= Columns(keys1);
		
		for (counter2 = 0; counter2 < keyC; counter2 = counter2 + 1)
		{
			key1Value = keys1[counter2];
			mx1 = av1[key1Value];
			mx2 = av2[key1Value];
			if (Abs(mx2))
			{
				/* merge the two matrices */
				
				d = Rows(mx1);
				mx3 = {d,2};
				mx3 [0][0] = mx1[0][0] + mx2[0][0];
				w1 = mx1[0][0] / mx3 [0][0];
				w2 = mx2[0][0] / mx3 [0][0];
				
				for (counter3 = 1; counter3 < d; counter3 = counter3 + 1)
				{
					mx3[counter3][0] = mx1[counter3][0];
					mx3[counter3][1] = mx1[counter3][1]*w1 + mx2[counter3][1]*w2;
				}
				avJ [key1Value] = mx3;
			}
			else
			{
				avJ [key1Value] = mx1;
			}
			doneKeys[key1Value] = 1;
		}		

		keyC		= Columns(keys2);
		for (counter2 = 0; counter2 < keyC; counter2 = counter2 + 1)
		{
			key1Value = keys2[counter2];
			if (doneKeys[key1Value] < 0.5)
			{
				avJ [key1Value] = av2[key1Value];
			}
		}
		
		nullAVLJoint[counter1] = avJ;
	}
	
	return nullAVLJoint;
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
	
	bufferString = "";
	bufferString * 8192;
	bufferString * titleMatrix[0][0];

	for (counter1=1; counter1<Columns(titleMatrix); counter1 = counter1+1)
	{
		bufferString * ("\t"+titleMatrix[0][counter1]);
	}
	
	for (counter1=0; counter1<Rows(dataMatrix); counter1 = counter1 + 1)
	{
		bufferString * ("\n"+dataMatrix[counter1][0]);
		for (counter2 = 1; counter2 < Columns (dataMatrix); counter2 = counter2+1)
		{
		    bufferString * ("\t"+dataMatrix[counter1][counter2]);
		}
	}
	
	bufferString * 0;
	fprintf (LAST_FILE_PATH,bufferString);
	bufferString = 0;
	
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
	
	ChoiceList	(distribChoice, "Test Statistic",1,SKIP_NODE,
				 "Approximate", "Use the approximate extended binomial distribution (fast)",
				 "Simulated Null", "Dynamically generate the null distribution from the data (very slow).");
	
	if (distribChoice < 0)
	{
		return;
	}
}


if (distribChoice > 0)
{
	STORE_ROOT_SUPPORT = 1;
}

DataSet dsA					= ReconstructAncestors (lf);

if (distribChoice > 0)
{
	STORE_ROOT_SUPPORT = 0;
	SUPPORT_MATRIX_LIST = Transpose (SUPPORT_MATRIX_LIST[0]);
}
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

for (k=1; k<filteredDataA.species; k=k+1)
{
	GetString (seqName, filteredDataA, k);
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

/*for (k=0; k<Rows(seqToBranchMap); k=k+1)
{
	v = seqToBranchMap[k][0];
	fprintf (stdout, "\n", branchNames[v], "=>");
	if (k<filteredData.species)
	{
		GetString (seqName, filteredData, k);
	}
	else
	{
		GetString (seqName, filteredDataA, k-filteredData.species);
	}
	fprintf (stdout, seqName);	
}

fprintf (stdout, "\n\n");

for (k=0; k<Rows(seqToBranchMap); k=k+1)
{
	fprintf (stdout, "\n", branchNames[k], "=>");
	v = seqToBranchMap[k][1];
	if (v<filteredData.species)
	{
		GetString (seqName, filteredData, v);
	}
	else
	{
		GetString (seqName, filteredDataA, v-filteredData.species);
	}
	fprintf (stdout, seqName);	
}

fprintf (stdout, "\n\n");

return 0;*/

/* total tree length */

totalTreeLength = 0;
branchLengths   = BranchLength(givenTree,-1);

for (v=Columns(branchLengths)-1; v>=0; v=v-1)
{
	totalTreeLength = totalTreeLength + branchLengths[v];
}

/* get codon matrix */

codonInfo  = {filteredData.species, filteredData.unique_sites};
codonInfo2 = {filteredDataA.species, filteredDataA.unique_sites};

GetDataInfo    (dupInfo, filteredData);
GetDataInfo	   (dupInfoA, filteredDataA);

matrixTrick  = {1,stateCharCount};
matrixTrick2 = {1,stateCharCount};
matrixTrick  = matrixTrick["_MATRIX_ELEMENT_COLUMN_"];
matrixTrick2 = matrixTrick2["1"];

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
		if ((matrixTrick2*siteInfo)[0] > 1)
		{
			codonInfo2[h][v] = -1;
		}
		else
		{
			codonInfo2[h][v] = (matrixTrick * siteInfo)[0];
		}
	}
}

_SITE_RESULTS = {4,filteredData.sites};
flatTreeRep	  = Abs (givenTree);

resultMatrix = {filteredData.sites,12};
perBranchMatrix = {Columns(branchNames),2};



if (distribChoice)
{
	dNdS = 	1;
	branchLengthsNeutral   = BranchLength(codonTree,-1);
	branchNamesList		   = BranchName (codonTree,-1);
	
	ClearConstraints (codonTree);
	/*fprintf (stdout, "\nTree before: ", Format (codonTree,1,1),"\n");*/
	
	for (v=Columns(branchLengthsNeutral)-1; v>=0; v=v-1)
	{
		ExecuteCommands ("codonTree."+branchNamesList[v]+".synRate = codonTree."+branchNamesList[v]+".synRate*branchLengths[v]/branchLengthsNeutral[v];"); 
	}	
	
	/*fprintf (stdout, "\nTree after: ", Format (codonTree,1,1),"\n");*/

	fprintf (stdout, "Simulating the null distribution for average branch lengths...\n");
	
	GetNeutralNull (simulatedNullAVL1, lf, _OBSERVED_S_, _OBSERVED_NS_, 33333);

	for (v=Columns(branchLengthsNeutral)-1; v>=0; v=v-1)
	{
		ExecuteCommands ("codonTree."+branchNamesList[v]+".synRate = codonTree."+branchNamesList[v]+".synRate*0.125"); 
	}	

	fprintf (stdout, "Simulating the null distribution for short branch lengths...\n");
	GetNeutralNull (simulatedNullAVL2, lf, _OBSERVED_S_, _OBSERVED_NS_, 33333);
	
	simulatedNullAVL = mergeSimulatedNulls(simulatedNullAVL1,simulatedNullAVL2); 
	
	for (v=Columns(branchLengthsNeutral)-1; v>=0; v=v-1)
	{
		ExecuteCommands ("codonTree."+branchNamesList[v]+".synRate = codonTree."+branchNamesList[v]+".synRate*32"); 
	}	

	fprintf (stdout, "Simulating the null distribution for long branch lengths...\n");
	GetNeutralNull (simulatedNullAVL3, lf, _OBSERVED_S_, _OBSERVED_NS_, 33333);
	simulatedNullAVL = mergeSimulatedNulls(simulatedNullAVL,simulatedNullAVL3); 
}

for (v=0; v<filteredData.sites;v=v+1)
{
	_SITE_ES_COUNT = {stateCharCount,stateCharCount};
	_SITE_EN_COUNT = {stateCharCount,stateCharCount};
	_SITE_OS_COUNT = {stateCharCount,stateCharCount};
	_SITE_ON_COUNT = {stateCharCount,stateCharCount};
	
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
		
		branchFactor = branchLengths[p1]/totalTreeLength;
		
		cd1 = codonInfo2[h] [c1];
		cd2 = codonInfo2[p2][c1];
		
		if (cd1 >= 0 && cd2 >= 0)
		{
			_SITE_OS_COUNT[cd1][cd2] = _SITE_OS_COUNT[cd1][cd2] + 1;		
			_SITE_ON_COUNT[cd1][cd2] = _SITE_ON_COUNT[cd1][cd2] + 1;		
			_SITE_ES_COUNT[cd1][cd2] = _SITE_ES_COUNT[cd1][cd2] + branchFactor;		
			_SITE_EN_COUNT[cd1][cd2] = _SITE_EN_COUNT[cd1][cd2] + branchFactor;		
			
			if (_PAIRWISE_S_[cd1][cd2])
			{
				perBranchMatrix[p1][0] = perBranchMatrix[p1][0]+_OBSERVED_S_[cd1][cd2]/_PAIRWISE_S_[cd1][cd2];
			}
			perBranchMatrix[p1][1] = perBranchMatrix[p1][1]+_OBSERVED_NS_[cd1][cd2]/_PAIRWISE_NS_[cd1][cd2];
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
		
		branchFactor = branchLengths[p1]/totalTreeLength;

		if (cd1>=0)
		/* no ambiguities */
		{
			/*p2 = flatTreeRep[p1];
			fprintf (stdout, "Change from ", codonString (cd1), " to ", codonString (cd2), " along ", branchNames[p1], " and ", branchNames[p2], "(",branchFactor,")\n");
			*/
			_SITE_OS_COUNT[cd1][cd2] = _SITE_OS_COUNT[cd1][cd2] + 1;		
			_SITE_ON_COUNT[cd1][cd2] = _SITE_ON_COUNT[cd1][cd2] + 1;		
			_SITE_ES_COUNT[cd1][cd2] = _SITE_ES_COUNT[cd1][cd2] + branchFactor;		
			_SITE_EN_COUNT[cd1][cd2] = _SITE_EN_COUNT[cd1][cd2] + branchFactor;		

			if (_PAIRWISE_S_[cd1][cd2])
			{
				perBranchMatrix[p1][0] = perBranchMatrix[p1][0]+_OBSERVED_S_[cd1][cd2]/_PAIRWISE_S_[cd1][cd2];
			}
			perBranchMatrix[p1][1] = perBranchMatrix[p1][1]+_OBSERVED_NS_[cd1][cd2]/_PAIRWISE_NS_[cd1][cd2];
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
						if (_PAIRWISE_S_[tempMx][cd2])
						{
							perBranchMatrix[p1][0] = perBranchMatrix[p1][0]+_OBSERVED_S_[tempMx][cd2]/_PAIRWISE_S_[tempMx][cd2];
						}
						perBranchMatrix[p1][1] = perBranchMatrix[p1][1]+_OBSERVED_NS_[tempMx][cd2]/_PAIRWISE_NS_[tempMx][cd2];
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
										
								if (_PAIRWISE_S_[k][cd2])
								{
									perBranchMatrix[p1][0] = perBranchMatrix[p1][0]+weightFactor*_OBSERVED_S_[k][cd2]/_PAIRWISE_S_[k][cd2];
								}
								perBranchMatrix[p1][1] = perBranchMatrix[p1][1]+weightFactor*_OBSERVED_NS_[k][cd2]/_PAIRWISE_NS_[k][cd2];
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
		if (distribChoice)
		{
			p3 = p2$1;
			p4 = ((resultMatrix[v][0]*6)+0.5)$1;
			
			waP          = 0;
			waN			 = 0;

			if (p3 != p2)
			{
				/* use averaging */
				fprintf (stdout, "Interpolating at site ", v, ": ", p3, " ", p3+1, " ", p2, "\n");
				
				w1 = p2-p3;
				w2 = 1-w1;
				
				for (kk2 = 0; kk2 < stateCharCount; kk2 = kk2+1)
				{			
					ambInfo = SUPPORT_MATRIX_LIST[kk2][v];
					if (ambInfo>0.000001)
					{
						subMatrixAVL = simulatedNullAVL[kk2];
						subMatrixAVL1 = subMatrixAVL[p3];
						subMatrixAVL2 = subMatrixAVL[p3+1];
						if (Abs(subMatrixAVL1) && Abs(subMatrixAVL2))
						{
							pv1 = subMatrixAVL1[p4+1][1];
							pv2 = subMatrixAVL2[p4+1][1];
							if (subMatrixAVL1[0][0] >= 100 && subMatrixAVL2[0][0] >= 100)
							{
								waP = waP + (pv1*w1+pv2*w2)*ambInfo;
								if (p4 == 0)
								{
									waN = waN + ambInfo;
								}
								else
								{
									waN = waN + (1 + w1*(pv1 - 2*subMatrixAVL1[p4][1]) + w2*(pv2 - 2*subMatrixAVL2[p4][1])) * ambInfo;
								}
							}
							else
							{
								fprintf (stdout, "Bailing at site ", v, " codon ", kk2, " sim count ", subMatrixAVL[0][0], "\n");
								break;
							}
						}
						else
						{
							fprintf (stdout, "Bailing at site ", v, " codon ", kk2 , " sim count 0 \n");
							break;					
						}
					}
				}				
			}
			else
			{
				for (kk2 = 0; kk2 < stateCharCount; kk2 = kk2+1)
				{			
					ambInfo = SUPPORT_MATRIX_LIST[kk2][v];
					if (ambInfo>0.000001)
					{
						subMatrixAVL = simulatedNullAVL[kk2];
						subMatrixAVL = subMatrixAVL[p3];
						if (Abs(subMatrixAVL))
						{
							if (subMatrixAVL[0][0] >= 100)
							{
								waP = waP + subMatrixAVL[p4+1][1]*ambInfo;
								if (p4 == 0)
								{
									waN = waN + ambInfo;
								}
								else
								{
									waN = waN + (1 + subMatrixAVL[p4+1][1] - 2*subMatrixAVL[p4][1]) * ambInfo;
								}
							}
							else
							{
								fprintf (stdout, "Bailing at site ", v, " codon ", kk2, " sim count ", subMatrixAVL[0][0], "\n");
								break;
							}
						}
						else
						{
							fprintf (stdout, "Bailing at site ", v, " codon ", kk2 , " sim count 0 \n");
							break;					
						}
					}
				}						
			}



			if (kk2 == stateCharCount)
			{
				resultMatrix[v][9]  = waP;
				resultMatrix[v][10] = waN;
				continue;
			}
		}
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
}

perBranchMatrix = perBranchMatrix * (1/filteredData.sites);
if (pipeThroughFlag == 0)
{

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
	
	fprintf (stdout, "\nTotal codon tree length (subs/nuc/unit time): ", Format (totalTreeLength,10,5), "\n");

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
		
		fprintf (stdout,"\n******* FOUND ", posSelected, " POSITIVELY SELECTED SITES ********\n\n");
		dummy = PrintASCIITable  (psMatrix, selLabelMatrix);
	}
	else
	{
		fprintf (stdout,"\n******* FOUND NO POSITIVELY SELECTED SITES ********\n\n");
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
		
		fprintf (stdout,"\n******* FOUND ", negSelected, " NEGATIVELY SELECTED SITES ********\n\n");
		dummy = PrintASCIITable  (psMatrix, selLabelMatrix);
	}
	else
	{
		fprintf (stdout,"\n******* FOUND NO NEGATIVELY SELECTED SITES ********\n\n");
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
			SetDialogPrompt ("Save result table to");
			dummy = PrintTableToFile  (resultMatrix, labelMatrix, !pipeThroughFlag);
		}
		else
		{
			OpenWindow (CHARTWINDOW,{{"Rates by sites"}
									{"labelMatrix"}
									{"resultMatrix"}
									{"Contrast Bars"}
									{"Index"}
									{labelMatrix[6]+";"+labelMatrix[7]}
									{"Site Index"}
									{"dN"}
									{"dS"}
									{"0"}
									{""}
									{"-1;-1"}
									{"10;1.309;0.785398"}
									{"Times:12:0;Times:10:0;Times:12:2"}
									{"0;0;16777215;1644825;0;0;6579300;11842740;13158600;14474460;0;3947580;16777215;5000268;6845928;16771158;2984993;9199669;7018159;1460610;16748822;11184810;14173291"}
									{"16,0,0"}
									},
									"(SCREEN_WIDTH-60)/2;(SCREEN_HEIGHT-50)/2;50;50");
		}
	}
	
	defString = "";
	defString * 8192;
	for (h=Rows(perBranchMatrix)-2; h>=0; h=h-1)
	{
		defString * ("_branchScaledTree."+branchNames[h]+".dS = "+perBranchMatrix[h][0]);
		defString * (";_branchScaledTree."+branchNames[h]+".dN = "+perBranchMatrix[h][1]+";");
		if (perBranchMatrix[h][1]+perBranchMatrix[h][0])
		{
			defString * ("_branchScaledTree."+branchNames[h]+".dNdS = "+Min(perBranchMatrix[h][1]/perBranchMatrix[h][0],5)+";");
		}
		else
		{
			defString * ("_branchScaledTree."+branchNames[h]+".dNdS = 0;");
		}
	}
	defString * 0;
	ExecuteCommands (defString);
	UseModel (USE_NO_MODEL);
	REPLACE_TREE_STRUCTURE = 1;
	Tree _branchScaledTree = ""+givenTree;
	REPLACE_TREE_STRUCTURE = 0;
	OpenWindow (TREEWINDOW,{{"_branchScaledTree"}{"8211"}{""}{"dNdS"}},"SCREEN_WIDTH/2-30;SCREEN_HEIGHT-120;100;100");
}
