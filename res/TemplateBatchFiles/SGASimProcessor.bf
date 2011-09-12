#include "binomial.ibf";

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
	incFileName = "Distances/CodonToolsMain.def";
}
else
{
	incFileName = "Distances/CodonTools.def";
}

ExecuteCommands  ("#include \""+incFileName+"\";");

SetDialogPrompt ("Please specify the first simulated file:");

DataSet 	  dsA 		   = ReadDataFile (PROMPT_FOR_FILE);
FILE_BASE 				   = LAST_FILE_PATH;

iterates = 0;

while (iterates<1)
{
	fprintf (stdout, "\nHow many data files should I process?");
	fscanf  (stdin, "Number", iterates);
}

it 	  = 2;
itidx = 0;

tailEnds = -1;

while ((tailEnds<=0) || (tailEnds >=.5))
{
	fprintf (stdout, "\nDistribution tails to tabulate (from 0 to 0.5)?");
	fscanf  (stdin, "Number", tailEnds);
}

tailEnds = tailEnds * 100;

ESResultMatrix = {iterates,filteredData.sites};
ENResultMatrix = {iterates,filteredData.sites};
OSResultMatrix = {iterates,filteredData.sites};
ONResultMatrix = {iterates,filteredData.sites};
dNmdS		   = {iterates,filteredData.sites};
pValuesPS	   = {iterates,filteredData.sites};
pValuesNS	   = {iterates,filteredData.sites};

while (1)
{
	DataSetFilter filteredDataA = CreateFilter (dsA,3,"","",GeneticCodeExclusions);

	HarvestFrequencies			  (observedCEFV,filteredData,3,3,0);
	seqToBranchMap 				= {stateCharCount,1};

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

	resultMatrix = {filteredData.sites,4};

	for (v=0; v<filteredData.sites;v=v+1)
	/*for (v=0; v<1;v=v+1)*/
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
						
			_SITE_OS_COUNT[cd1][cd2] = _SITE_OS_COUNT[cd1][cd2] + 1;		
			_SITE_ON_COUNT[cd1][cd2] = _SITE_ON_COUNT[cd1][cd2] + 1;		
			_SITE_ES_COUNT[cd1][cd2] = _SITE_ES_COUNT[cd1][cd2] + branchFactor;		
			_SITE_EN_COUNT[cd1][cd2] = _SITE_EN_COUNT[cd1][cd2] + branchFactor;		
			
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
		
		_SITE_OS_COUNT = matrixTrick2*(_OBSERVED_S_$_SITE_OS_COUNT)*Transpose(matrixTrick2);
		_SITE_ON_COUNT = matrixTrick2*(_OBSERVED_NS_$_SITE_ON_COUNT)*Transpose(matrixTrick2);
		_SITE_ES_COUNT = matrixTrick2*(_PAIRWISE_S_$_SITE_ES_COUNT)*Transpose(matrixTrick2);
		_SITE_EN_COUNT = matrixTrick2*(_PAIRWISE_NS_$_SITE_EN_COUNT)*Transpose(matrixTrick2);
		
		OSResultMatrix[itidx][v] = _SITE_OS_COUNT[0];
		ONResultMatrix[itidx][v] = _SITE_ON_COUNT[0];
		ESResultMatrix[itidx][v] = _SITE_ES_COUNT[0];
		ENResultMatrix[itidx][v] = _SITE_EN_COUNT[0];
		
		weight = _SITE_EN_COUNT[0]+_SITE_ES_COUNT[0];
		
		p = _SITE_ES_COUNT[0]/weight;
		p2 = OSResultMatrix[itidx][v]+ONResultMatrix[itidx][v];
		ds = OSResultMatrix[itidx][v]/ESResultMatrix[itidx][v];
		dn = ONResultMatrix[itidx][v]/ENResultMatrix[itidx][v];
		
		if (ESResultMatrix[itidx][v])
		{
			dNmdS[itidx][v] = dn-ds;	
		}
		
		if (p2)
		{
			pv  = extendedBinTail (p2,p,OSResultMatrix[itidx][v]);
			if (OSResultMatrix[itidx][v]>=1)
			{
				pv2 = 1-pv+(pv-extendedBinTail(p2,p,OSResultMatrix[itidx][v]-1));
			}
			else
			{
				pv2 = 1-pv+(pv-extendedBinTail (p2,p,0));
			}
			pValuesPS[itidx][v] = pv;
			pValuesNS[itidx][v] = pv2;
		}				
	}
	if (it <= iterates)
	{		
		fprintf (stdout, "\nIterate ", it, "/", iterates);
		fName = FILE_BASE+"_"+it;
		DataSet 	  dsA 		   = ReadDataFile (fName);
		it = it+1;
		itidx = itidx + 1;
	}
	else
	{
		break;
	}
}

fprintf (stdout, "\n");

/* generate medians and tailEnd% and (100-tailEnd)% profiles */

statisticalOverview = {filteredData.sites,12};

tailEnd1 = " "+tailEnds+"%";
tailEnd2 = " "+(100-tailEnds)+"%";

sovLabels = {{"","ES Median","","","EN Median","","","OS Median","","","ON Median",""}};

sovLabels[0]  = "ES"+tailEnd1;
sovLabels[2]  = "ES"+tailEnd2;
sovLabels[3]  = "EN"+tailEnd1;
sovLabels[5]  = "EN"+tailEnd2;
sovLabels[6]  = "OS"+tailEnd1;
sovLabels[8]  = "OS"+tailEnd2;
sovLabels[9]  = "ON"+tailEnd1;
sovLabels[11] = "ON"+tailEnd2;

ColumnToSort = {iterates,1};

lowBar  = (iterates*tailEnds*0.01)$1;
highBar = (iterates*(100-tailEnds)*0.01)$1;

/*if ((iterates*9)%10==0)
{
	highBar = highBar-1;
}*/

median	= iterates$2;
doAv	= (iterates%2==0);

for (v=0;v<filteredData.sites;v=v+1)
{
	for (h=0;h<iterates;h=h+1)
	{
		ColumnToSort[h]=ESResultMatrix[h][v];
	}
	ColumnToSort = ColumnToSort%0;
	statisticalOverview[v][0] = ColumnToSort[lowBar];
	statisticalOverview[v][2] = ColumnToSort[highBar];
	if (doAv)
	{
		statisticalOverview[v][1] = (ColumnToSort[median]+ColumnToSort[median+1])/2;	
	}
	else
	{
		statisticalOverview[v][1] = ColumnToSort[median];
	}
}

for (v=0;v<filteredData.sites;v=v+1)
{
	for (h=0;h<iterates;h=h+1)
	{
		ColumnToSort[h]=ENResultMatrix[h][v];
	}
	/*dummy = doTheSort (iterates);*/
	ColumnToSort = ColumnToSort%0;
	statisticalOverview[v][3] = ColumnToSort[lowBar];
	statisticalOverview[v][5] = ColumnToSort[highBar];
	if (doAv)
	{
		statisticalOverview[v][4] = (ColumnToSort[median]+ColumnToSort[median+1])/2;	
	}
	else
	{
		statisticalOverview[v][4] = ColumnToSort[median];
	}
}

for (v=0;v<filteredData.sites;v=v+1)
{
	for (h=0;h<iterates;h=h+1)
	{
		ColumnToSort[h]=OSResultMatrix[h][v];
	}
	ColumnToSort = ColumnToSort%0;
	statisticalOverview[v][6] = ColumnToSort[lowBar];
	statisticalOverview[v][8] = ColumnToSort[highBar];
	if (doAv)
	{
		statisticalOverview[v][7] = (ColumnToSort[median]+ColumnToSort[median+1])/2;	
	}
	else
	{
		statisticalOverview[v][7] = ColumnToSort[median];
	}
}

for (v=0;v<filteredData.sites;v=v+1)
{
	for (h=0;h<iterates;h=h+1)
	{
		ColumnToSort[h]=ONResultMatrix[h][v];
	}
	ColumnToSort = ColumnToSort%0;
	statisticalOverview[v][9] = ColumnToSort[lowBar];
	statisticalOverview[v][11] = ColumnToSort[highBar];
	if (doAv)
	{
		statisticalOverview[v][10] = (ColumnToSort[median]+ColumnToSort[median+1])/2;	
	}
	else
	{
		statisticalOverview[v][10] = ColumnToSort[median];
	}
}

selLabelMatrix = {1,filteredData.sites};

for (h=1; h<=filteredData.sites;h=h+1)
{
	selLabelMatrix[h-1] = "Site " + h;
}

OpenWindow (CHARTWINDOW,{{"Simulated ES"}
						   {"selLabelMatrix"},
						   {"ESResultMatrix"},
						   {"None"},
						   {"Index"},
						   {""},
						   {""},
						   {""},
						   {""},
						   {"0"}},
						   "(SCREEN_WIDTH-60)/2;(SCREEN_HEIGHT-60)/2;15;30");

OpenWindow (CHARTWINDOW,{{"Simulated EN"}
						   {"selLabelMatrix"},
						   {"ENResultMatrix"},
						   {"None"},
						   {"Index"},
						   {""},
						   {""},
						   {""},
						   {""},
						   {"0"}},
						   "(SCREEN_WIDTH-60)/2;(SCREEN_HEIGHT-60)/2;30+(SCREEN_WIDTH-60)/2;30");

OpenWindow (CHARTWINDOW,{{"Simulated OS"}
						   {"selLabelMatrix"},
						   {"OSResultMatrix"},
						   {"None"},
						   {"Index"},
						   {""},
						   {""},
						   {""},
						   {""},
						   {"0"}},
						   "(SCREEN_WIDTH-60)/2;(SCREEN_HEIGHT-60)/2;15;45+(SCREEN_HEIGHT-60)/2");

OpenWindow (CHARTWINDOW,{{"Simulated ON"}
						   {"selLabelMatrix"},
						   {"ONResultMatrix"},
						   {"None"},
						   {"Index"},
						   {""},
						   {""},
						   {""},
						   {""},
						   {"0"}},
						   "(SCREEN_WIDTH-60)/2;(SCREEN_HEIGHT-60)/2;30+(SCREEN_WIDTH-60)/2;45+(SCREEN_HEIGHT-60)/2");

OpenWindow (CHARTWINDOW,{{"Summary"}
						   {"sovLabels"},
						   {"statisticalOverview"},
						   {"Line Plot"},
						   {"Index"},
						   {""},
						   {""},
						   {""},
						   {""},
						   {"0"}},
						   "(SCREEN_WIDTH-200);(SCREEN_HEIGHT-200);100;100");

