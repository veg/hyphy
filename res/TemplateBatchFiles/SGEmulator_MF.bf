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

ChoiceList (ambChoice, "Treatment of Ambiguities",1,SKIP_NONE,
			"Averaged","All possible resolutions are considered and averaged.",
			"Resolved","The most frequent (for that site) resolution is chosen.");

if (ambChoice<0)
{
	return;
}

if (useCustomCountingBias)
{
	ExecuteAFile("Distances/CodonToolsMain.def");
}
else
{
	ExecuteAFile("Distances/CodonTools.def");
}

observedCEFV 		  = {64,1};
for (fileID = 1; fileID <= fileCount; fileID = fileID + 1)
{
	ExecuteCommands 	  		("HarvestFrequencies (tp, filteredData_"+fileID+",3,3,0);cfs = filteredData_"+fileID+".sites;");
	observedCEFV 				= observedCEFV 		 + tp*(cfs/totalCodonCount);
}


seqToBranchMap 				  			= {stateCharCount,1};
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

vOffset 	 = 0;
resultMatrix = {totalCodonCount,12};

treeLengthArray = {};

for (fileID = 1; fileID <= fileCount; fileID = fileID + 1)
{
	ExecuteCommands ("LikelihoodFunction tempLF = (filteredData_"+fileID+",codonTree_"+fileID+");");
	DataSet 		dsA		 				= ReconstructAncestors (tempLF);
	ExecuteCommands ("DataSet		   	dsJoint 				= Combine(dsA,ds_"+fileID+");");
	ExecuteCommands ("DataSetFilter		filteredData 			= CreateFilter (ds_"+fileID+",3,\"\",\"\",GeneticCodeExclusions);");
	
	DataSetFilter 	filteredDataA = CreateFilter (dsA,3,"","",GeneticCodeExclusions);
	DataSetFilter  filteredDataJ  = CreateFilter (dsJoint,3,"","",GeneticCodeExclusions);

	ExecuteCommands ("branchNames = BranchName (codonTree_"+fileID+",-1);");
	h 			= Columns (branchNames);

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


	/* total tree length */

	totalTreeLength = 0;
	
	ExecuteCommands ("branchLengths   = BranchLength(nucTree_"+fileID+",-1);");
	for (v=Columns(branchLengths)-1; v>=0; v=v-1)
	{
		totalTreeLength = totalTreeLength + branchLengths[v];
	}
	
	treeLengthArray [fileID] = totalTreeLength;

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
	ExecuteCommands ("flatTreeRep	  = Abs (nucTree_"+fileID+");");

	perBranchMatrix = {Columns(branchNames),2};

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
		
		shiftedV = v+vOffset;
		
		resultMatrix[shiftedV][0] = _SITE_OS_COUNT[0];
		resultMatrix[shiftedV][1] = _SITE_ON_COUNT[0];
		resultMatrix[shiftedV][2] = _SITE_ES_COUNT[0];
		resultMatrix[shiftedV][3] = _SITE_EN_COUNT[0];
		
		weight = _SITE_EN_COUNT[0]+_SITE_ES_COUNT[0];
		
		p = _SITE_ES_COUNT[0]/weight;
		
		resultMatrix[shiftedV][5] = p;
		
		p2 = resultMatrix[shiftedV][0]+resultMatrix[shiftedV][1];
		
		resultMatrix[shiftedV][7] = resultMatrix[shiftedV][1]/resultMatrix[shiftedV][3];
		
		if (resultMatrix[shiftedV][2])
		{
			resultMatrix[shiftedV][6]  = resultMatrix[shiftedV][0]/resultMatrix[shiftedV][2];
			resultMatrix[shiftedV][8]  = resultMatrix[shiftedV][7]-resultMatrix[shiftedV][6];	
			resultMatrix[shiftedV][11] = resultMatrix[shiftedV][8]/totalTreeLength;	
		}
		
		if (p2)
		{
			resultMatrix[shiftedV][4]  = resultMatrix[shiftedV][0]/p2;		
			if (distribChoice)
			{
				p3 = p2$1;
				p4 = ((resultMatrix[shiftedV][0]*6)+0.5)$1;
				
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
					resultMatrix[shiftedV][9]  = waP;
					resultMatrix[shiftedV][10] = waN;
					continue;
				}
			}
			resultMatrix[shiftedV][9]  = extendedBinTail (p2,p,resultMatrix[shiftedV][0]);
			if (resultMatrix[shiftedV][0]>=1)
			{
				resultMatrix[shiftedV][10] = 1-extendedBinTail(p2,p,resultMatrix[shiftedV][0]-1);
			}
			else
			{
				resultMatrix[shiftedV][10] = 1-extendedBinTail (p2,p,0);
			}
		}		
	}
	vOffset = vOffset + filteredData.sites;
}

perBranchMatrix = perBranchMatrix * (1/totalCodonCount);
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

for (fileID = 1; fileID <= fileCount; fileID = fileID + 1)
{
	fprintf (stdout, "\nTotal length (subs/nuc/unit time) for tree ", fileID, ": ", Format (treeLengthArray[fileID],10,5));
}

fprintf (stdout, "\n");

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
