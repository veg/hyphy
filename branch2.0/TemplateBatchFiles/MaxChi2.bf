R_2 = {{1,1}};
C_2 = {{1}{1}};

SetDialogPrompt ("Please specify a nucleotide or amino-acid data file:");

DataSet ds = ReadDataFile (PROMPT_FOR_FILE);

fprintf (stdout, "Read ", ds.species, " sequences and ", ds.sites, " sites.");

/* find sites which are not constant */

DataSetFilter filteredData = CreateFilter (ds,1,"/^(.)\1+$/");

vsc = ds.sites-filteredData.sites;
fprintf (stdout, " Found ", vsc, " variable sites");

/* make a list of variable sites */

k = 0;
i = 0; 
variableMap 	= {vsc,1};
map 			= {};
toVariableMap   = {};

for (k2 = 0; k2 < ds.sites; k2=k2+1)
{
	if (filteredData.site_map[i] == k2)
	{
		if (i < filteredData.sites - 1)
		{
			i=i+1;
		}
	}
	else
	{
		map[k2] = 1;
		variableMap[k] = k2;
		toVariableMap[k2] = k;
		k = k+1;
	}
}

pairDiffMap = {};

DataSetFilter filteredData = CreateFilter (ds,1,map[siteIndex]);

for (k=0; k<filteredData.species; k=k+1)
{
	for (k2=k+1; k2<filteredData.species; k2=k2+1)
	{
		DataSetFilter seqFilter    = CreateFilter (filteredData,1,"",speciesIndex == k || speciesIndex == k2);
		DataSetFilter pairedFilter = CreateFilter (seqFilter,1,"/^(.)\1+$/");
		mismatchSites = {1,filteredData.sites};	

		i  = 0;
		for (k3 = 0; k3 < filteredData.sites; k3=k3+1)
		{
			k4 = filteredData.site_map[k3];
			if (pairedFilter.site_map[i] == k4)
			{
				if (i < pairedFilter.sites - 1) 
				{
					i=i+1;
				}
			}
			else
			{
				k4 = toVariableMap[k4];
				mismatchSites[k4] = 1;
			}
		}
		
		pairDiffMap[""+k+","+k2] = mismatchSites;
	}
}

ChoiceList	(useFisherExact,"Test to use",1,SKIP_NONE,
					   "Chi square",      "Compute association bias using the chi-squared approximation (faster).",
					   "Fisher's exact",  "Compute association bias using Fisher's exact test approximation (slower)");

if (useFisherExact<0)
{
	return 0;
}

minWindowSize  = promptForValue ("Minimum partition size ", 1, filteredData.sites-1,filteredData.sites$3,1);
iterates	   = promptForValue ("Boostrap sample size ", 1, 100000000 ,1000,1);
reportP	   	   = promptForValue ("p-value threshold ", 0, 1 ,0.05,0);

SetDialogPrompt ("Save results to:");
fprintf 		(PROMPT_FOR_FILE, CLEAR_FILE, "Sequence 1,Sequence 2,Breakpoint,Variable_Left,Variable_Right,Test p-value,Resampled p-value");
saveTo	= LAST_FILE_PATH;

windowCount	   = filteredData.sites-2*minWindowSize+1;

unitVector	= {filteredData.sites,1};
stencil  	= {filteredData.sites,1};

pairTotal = filteredData.species*(filteredData.species-1)/2;
completeReport = {pairTotal*windowCount,7};

for (bppos = 0; bppos < filteredData.sites; bppos = bppos + 1)
{
	unitVector[bppos] = 1;
}

rowIndex  = 0;
lastSpool = 0;

timer = Time(0);
pairCount = 0;

for (k=0; k<filteredData.species; k=k+1)
{
	GetString (seq1,filteredData,k);
	for (k2=k+1; k2<filteredData.species; k2=k2+1)
	{
		GetString (seq2,filteredData,k2);
		if (pairCount)
		{
			fprintf (stdout, "\nPair (", k+1, ",", k2+1, "). Estimated Remaining Time (seconds): ", (Time(0)-timer)*(pairTotal-pairCount)/pairCount, "\n");
		}
		else
		{
			fprintf (stdout, "\nPair (", k+1, ",", k2+1, ")\n");		
		}
		pairCount = pairCount + 1;
		pairwiseDiffs = pairDiffMap[""+k+","+k2];
		td = pairwiseDiffs*unitVector;
		td = td[0];
		stencil  	= {filteredData.sites,1};
		for (bppos = 0; bppos < minWindowSize; bppos = bppos + 1)
		{
			stencil[bppos] = 1;
		}

		for (bppos = minWindowSize-1; bppos < filteredData.sites-minWindowSize; bppos = bppos+1)
		{
			stencil[bppos] = 1;
			completeReport[rowIndex][0] = k+1;
			completeReport[rowIndex][1] = k2+1;
			completeReport[rowIndex][2] = variableMap[bppos];
			pdc = pairwiseDiffs*stencil;
			pdc = pdc[0];
			completeReport[rowIndex][3] = pdc;
			completeReport[rowIndex][4] = td-pdc;
			twoByTwo = {2,2};
			twoByTwo[0][0] = pdc;
			twoByTwo[0][1] = td-pdc;
			twoByTwo[1][0] = bppos+1;
			twoByTwo[1][1] = filteredData.sites-twoByTwo[1][0];
			
			myP = getP (twoByTwo);
			
			completeReport[rowIndex][5] = myP;
			simP = 0;
			for (it = 0; it < iterates; it=it+1)
			{
				randomized = Random(pairwiseDiffs,0);
				pdc = randomized*stencil;
				pdc = pdc[0];
				twoByTwo[0][0] = pdc;
				twoByTwo[0][1] = td-pdc;
				simP2 = getP (twoByTwo);
				if (simP2 <= myP)
				{
					simP = simP+1;
				}
			}
			simP = simP / iterates;
			completeReport[rowIndex][6] = simP;
			if (simP <= reportP)
			{
				fprintf (stdout, "Pair (", k+1, ",", k2+1, ") position ", completeReport[rowIndex][2], " p-value:", simP, " (*)\n");
			}
			rowIndex = rowIndex+1;
		}
		outString = "";
		outString * 512;
		for (bppos = lastSpool; bppos < rowIndex; bppos = bppos+1)
		{
			outString * ("\n"+seq1+","+seq2+","+completeReport[bppos][2]+","+completeReport[bppos][3]+","+completeReport[bppos][4]+","
						+completeReport[bppos][5]+","+completeReport[bppos][6]);
		}
		lastSpool = rowIndex;
		outString * 0;
		fprintf (saveTo, outString);
	}
}

/*------------------------------------------------------------------------------------------------------------*/

function getP (aMatrix)
{
	if (useFisherExact)
	{
		return CChi2(aMatrix,0);
	}
	else
	{
		R_S = aMatrix * C_2;
		C_S = R_2*aMatrix;
		T = C_S[0]+C_S[1];
		
		n_ij = R_S[0]*C_S[0]/T;
		chi2Sum = (aMatrix[0][0]-n_ij)^2 / n_ij;
		n_ij = R_S[0]*C_S[1]/T;
		chi2Sum = chi2Sum + (aMatrix[0][1]-n_ij)^2 / n_ij;
		n_ij = R_S[1]*C_S[0]/T;
		chi2Sum = chi2Sum + (aMatrix[1][0]-n_ij)^2 / n_ij;
		n_ij = R_S[1]*C_S[1]/T;
		chi2Sum = chi2Sum + (aMatrix[1][1]-n_ij)^2 / n_ij;
		
		return 1-CChi2(chi2Sum,1);
	}
}

/*------------------------------------------------------------------------------------------------------------*/

function promptForValue (pString, minV, maxV, defV,roundMe)
{
	while (1)
	{
		inString = "";
		fprintf (stdout, "\n", pString, "[", minV, "-", maxV,"] (default ", defV, "):");
		fscanf (stdin, "String", inString); 
		if (Abs(inString))
		{
			unitVector = (0+inString);
			if (roundMe)
			{
				unitVector = unitVector $ 1;
			}
			if (unitVector >= minV && unitVector <= maxV)
			{
				return unitVector;
			}
		}
		else
		{
			return defV;
		}
	}

}
