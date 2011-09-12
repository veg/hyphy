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

ExecuteAFile (incFileName);
	
DataSet dsA					= ReconstructAncestors (lf);
DataSetFilter filteredDataA = CreateFilter (dsA,3,"","",GeneticCodeExclusions);

HarvestFrequencies			  (observedCEFV,filteredData,3,3,0);
seqToBranchMap 				  = {stateCharCount,1};

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

_SITE_RESULTS = {Rows(seqToBranchMap)-1,filteredData.sites};
flatTreeRep	  = Abs (givenTree);


for (v=0; v<filteredData.sites;v=v+1)
{	
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
		
		cd1 = codonInfo2[h] [c1];
		cd2 = codonInfo2[p2][c1];
		
		if (_OBSERVED_NS_[cd1][cd2] > 0.5)
		{
			_SITE_RESULTS[p1][v] = 1; 
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
		
		observedCEFV = {{0}};
		
		if (cd1>=0)
		/* no ambiguities */
		{
			if (_OBSERVED_NS_[cd1][cd2] > 0.5)
			{
				_SITE_RESULTS[p1][v] = 1; 
			}		
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
						if (_OBSERVED_NS_[tempMx][cd2] > 0.5)
						{
							_SITE_RESULTS[p1][v] = 1; 
						}		
					}
				}
				else
				{
					weightFactor = matrixTrick2*ambInfo;
					weightFactor = weightFactor[0];
					if (weightFactor)
					{
						ambInfo		 = ambInfo * (1/weightFactor);
						wk			 = 0;
						for (k=0; k<stateCharCount; k=k+1)
						{
							wk = wk + _OBSERVED_NS_[k][cd2]*ambInfo[k];
						}
						if (wk > 0.5)
						{
							_SITE_RESULTS[p1][v] = 1; 
						}		
					}
				}
			}
		}
	}	
}
