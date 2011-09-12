ExecuteAFile  ( HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"chooseGeneticCode.def");

/*___________________________________________________________________________________________________________*/

_S_NS_POSITIONS_ = {3,64};

hShift = 0;

for (h=0; h<64; h=h+1)
{
	myAA = _Genetic_Code[h];
	
	if (myAA != 10) /* not a stop codon */
	{
		norm_factor  = 0.0;
		sSites	 	 = 0.0;
		nsSites 	 = 0.0;
		
		/* first position change */
		/* actual first position */
		
		p1 = h$16; /* 0->A, 1->C, 2->G, 3->T */
		p2 = h%16; /* remainder - i.e. positions 2 and 3*/
				
		for (n1 = 0; n1 < 4; n1=n1+1)
		{
			if (n1 != p1) /* a change */
			{
				/* new codon */
				nc = n1*16 + p2;
				newAA = _Genetic_Code[nc];
				
				if (newAA!=10) /* not a stop codon */
				{
					if (newAA == myAA) /* syn. change */
					{
						sSites = sSites + 1;
					}
					else
					{
						nsSites = nsSites + 1;
					}
					norm_factor = norm_factor + 1;
				}
			}
		}
		
		if (norm_factor)
		{
			_S_NS_POSITIONS_[0][h-hShift] = _S_NS_POSITIONS_[0][h-hShift] + sSites;
			_S_NS_POSITIONS_[1][h-hShift] = _S_NS_POSITIONS_[1][h-hShift] + nsSites;
			_S_NS_POSITIONS_[2][h-hShift] = _S_NS_POSITIONS_[2][h-hShift] + norm_factor;
		}
		
		norm_factor  = 0.0;
		sSites	 	 = 0.0;
		nsSites 	 = 0.0;

		/* second position change */
		/* actual second position */
		
		p1 = (h%16)$4;
		p2 = (h$16)*16+h%4; /* remainder - i.e. positions 1 and 3*/
		
		for (n1 = 0; n1 < 4; n1=n1+1)
		{
			if (n1 != p1) /* a change */
			{
				/* new codon */
				nc = n1*4 + p2;
				newAA = _Genetic_Code[nc];
				
				if (newAA!=10) /* not a stop codon */
				{
					if (newAA == myAA) /* syn. change */
					{
						sSites = sSites + 1;
					}
					else
					{
						nsSites = nsSites + 1;
					}
					norm_factor = norm_factor + 1;
				}
			}
		}

		/* 3rd position change */
		/* actual 3rd position */
		
		if (norm_factor)
		{
			_S_NS_POSITIONS_[0][h-hShift] = _S_NS_POSITIONS_[0][h-hShift] + sSites;
			_S_NS_POSITIONS_[1][h-hShift] = _S_NS_POSITIONS_[1][h-hShift] + nsSites;
			_S_NS_POSITIONS_[2][h-hShift] = _S_NS_POSITIONS_[2][h-hShift] + norm_factor;
		}
		
		norm_factor  = 0.0;
		sSites	 	 = 0.0;
		nsSites 	 = 0.0;

		p1 = h%4;
		p2 = (h$4)*4; /* remainder - i.e. positions 1 and 2*/
		
		for (n1 = 0; n1 < 4; n1=n1+1)
		{
			if (n1 != p1) /* a change */
			{
				/* new codon */
				nc = n1 + p2;
				newAA = _Genetic_Code[nc];
				
				if (newAA!=10) /* not a stop codon */
				{
					if (newAA == myAA) /* syn. change */
					{
						sSites = sSites + 1;
					}
					else
					{
						nsSites = nsSites + 1;
					}
					norm_factor = norm_factor + 1;
				}
			}
		}
		
		if (norm_factor)
		{
			_S_NS_POSITIONS_[0][h-hShift] = _S_NS_POSITIONS_[0][h-hShift] + sSites;
			_S_NS_POSITIONS_[1][h-hShift] = _S_NS_POSITIONS_[1][h-hShift] + nsSites;
			_S_NS_POSITIONS_[2][h-hShift] = _S_NS_POSITIONS_[2][h-hShift] + norm_factor;
		}
	}
	else
	{
		hShift = hShift+1;
	}
}

senseCodons = 64-hShift;

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


/*----------------------------------------------------------------------------*/


SetDialogPrompt				("Select a data file:");

DataSet		  ds		   = ReadDataFile	(PROMPT_FOR_FILE);
DataSetFilter filteredData = CreateFilter   (ds,3,"","",GeneticCodeExclusions);
HarvestFrequencies			  (observedCEFV,filteredData,3,3,0);

hShift = 0;

stateCharCount = 0;

for (k=0; k<64; k=k+1)
{
	if (_Genetic_Code[k]!=10)
	{
		stateCharCount = stateCharCount + 1;
	}
}

hShift = 0;

seqToBranchMap = {stateCharCount,1};
codeRemap	   = {stateCharCount,1};

for (k=0; k<64; k=k+1)
{
	if (_Genetic_Code[k]==10)
	{
		hShift = hShift+1;
	}
	else
	{
		seqToBranchMap[k-hShift] = observedCEFV[k];
		codeRemap	  [k-hShift] = k;
	}
}

observedCEFV = seqToBranchMap;


matrixTrick  = {1,stateCharCount};
matrixTrick2 = {1,stateCharCount};

for (h=Columns(matrixTrick)-1; h>=0; h=h-1)
{
	matrixTrick  [h] = h;
	matrixTrick2 [h] = 1;
}

ChoiceList (allowed, "Allowed ambiguities",1,SKIP_NONE,
			"2",	 		"Two way ambiguities only",
			"2,3", 		 	"Two and three way ambiguities only.",
			"All",	 		"All ambiguities");

if (allowed < 0)
{
	return;
}

allowedChars = "";

if (allowed == 0)
{
	allowedChars = "ACGTUKMRSWY";
}
else
{
	if (allowed == 1)
	{
		allowedChars = "ACGTUKMRSWYBDHV";
	}
}

allowedStates = {};

for (v=0; v<Abs(allowedChars); v=v+1)
{
	h = allowedChars[v];
	allowedStates[h] = 1;
}

SetDialogPrompt				("Specify output file:");

fprintf (PROMPT_FOR_FILE,CLEAR_FILE,"Sequence ,Site, Synonymous Ambigs, Non-Synonymous Ambigs, Ambiguous Codon,Syn neighbors, Non-syn neighbors, Sense neighbors\n");

codonInfo  = {filteredData.species, filteredData.unique_sites};

disallowed  = 0;
siteLookups = {filteredData.unique_sites, senseCodons};

for (v=0; v<filteredData.unique_sites;v=v+1)
{
	for (h=0; h<filteredData.species;h=h+1)
	{
		GetDataInfo (siteInfo, filteredData, h, v);
		_SITE_ES_COUNT = matrixTrick2 * siteInfo; 
		if (_SITE_ES_COUNT[0] == 1)
		{
			siteInfo 		   = matrixTrick * siteInfo;
			k2 				   = siteInfo[0];
			codonInfo[h][v]    = k2;
			siteLookups[v][k2] = 1;
		}
		else
		{
			codonInfo[h][v] = -1;
		}
	}
}

GetDataInfo    (dupInfo, filteredData);
GetInformation (sequenceData, filteredData);

_SITE_RESULTS = {4,filteredData.sites};

ambigList 		  = {};
ambigTabulator	  = {};

fprintf (stdout, "Sequence\tSyn Ambigs\tNon-syn ambigs\tTotal ambigs\n");


for (h=0; h<filteredData.species; h=h+1)
{		
	outString  = "";
	outString2 = "";
	dummy = outString*1024;
	GetString (seqName,filteredData,h);
	
	fprintf (stdout, "\n", seqName,"\t");
	
	seqString = sequenceData[h];
	tsCount    = 0;
	tnsCount   = 0;
	
	for (v=0; v<filteredData.sites;v=v+1)
	{
		c2  = dupInfo[v];
		cd1 = codonInfo [h]  [c2];
		
		if (cd1>=0)
		/* no ambiguities */
		{
			outString*(seqName+","+(v+1)+",0,0,No ambiguities,"+_S_NS_POSITIONS_[0][cd1]+","+_S_NS_POSITIONS_[1][cd1]+","+_S_NS_POSITIONS_[2][cd1]+"\n");
		}	
		else
		{
						
			aCodon = seqString[v*3][v*3+2];
			ambigTabulator[aCodon] = ambigTabulator[aCodon]+1;
			if (Abs(ambigList[aCodon]) == 0)
			{
				sCount 	= aCodon[0];
				nsCount	= aCodon[1];
				cd1		= aCodon[2];
				
				sCount  = allowedStates[sCount];
				nsCount = allowedStates[nsCount];
				cd1		= allowedStates[cd1];
				
				if ((allowed>1)||(sCount&&nsCount&&cd1))
				{
					nsCount = 0;
					sCount	= 0;
					sSites  = 0;
					nsSites = 0;
					
					
					GetDataInfo (ambInfo , filteredData, h, c2);
					ambInfo  	 = ambInfo$observedCEFV;
					weightFactor = (matrixTrick2*ambInfo)[0];
					if (weightFactor)
					{
						ambInfo		 = ambInfo * (1/weightFactor);
						resolutions  = {};
						for (k=0; k<stateCharCount; k=k+1)
						{
							if (ambInfo[k]>0)
							{
								resolutions[Abs(resolutions)] = k;
							}
						}
						
						cd1 = Abs(resolutions);
						
						if (cd1)
						{
							weightFactor = 0;
							for (k=0; k<cd1;k=k+1)
							{
								f1 = resolutions[k];
								sSites  = sSites  + _S_NS_POSITIONS_[0][f1];
								nsSites = nsSites + _S_NS_POSITIONS_[1][f1];
								c2 = codeRemap [f1];
								c2 = _Genetic_Code[c2];
								
								for (k2=k+1;k2<cd1; k2=k2+1)
								{
									f2 = resolutions[k2];
									
									dummy = ambInfo[f1]*ambInfo[f2];
									weightFactor = weightFactor + dummy;
									
									f2 = codeRemap[f2];
																			
									if (c2 == _Genetic_Code[f2])
									{
										sCount = sCount + dummy;
									}
									else
									{
										nsCount = nsCount + dummy;
									}
								}
							}
						}
						
						sCount = sCount/weightFactor;
						nsCount = nsCount/weightFactor;
						sSites  = sSites/cd1;
						nsSites  = nsSites/cd1;
					}
				}
				else
				{
					nsCount = 0;
					sCount	= 0;
					sSites  = 0;
					nsSites = 0;
					disallowed = disallowed + 1;
				}
				
				cd1 = {5,1};
				cd1[0] = sCount;
				cd1[1] = nsCount;
				cd1[2] = sSites;
				cd1[3] = nsSites;
				cd1[4] = 1;
				ambigList[aCodon] = cd1;
			}
			else
			{
				cd1 = ambigList[aCodon];
				sCount  = cd1[0];
				nsCount = cd1[1];
				sSites  = cd1[2];
				nsSites = cd1[3];
			}
			
			outString*(seqName+","+(v+1)+","+sCount+","+nsCount+","+aCodon+","+sSites+","+nsSites+","+(sSites+nsSites)+"\n");
			tsCount  = tsCount  + sCount;
			tnsCount = tnsCount + nsCount;
		}
	}
	
	outString*0;
	fprintf (LAST_FILE_PATH, outString);
	fprintf (stdout, tsCount, "\t", tnsCount, "\t", tsCount + tnsCount);
}

fprintf (stdout, "\n\n");

fprintf (stdout, "\nA total of ", Abs(ambigList), " distinct ambiguous patterns were found\n", disallowed, " were disallowed");

cd1 = Rows(ambigTabulator);
k 	= Columns (cd1);
cd2 = {k,2};

for (h=0; h<k; h=h+1)
{
	k2  = cd1[h];
	cd2 [h][0] = h;
	cd2 [h][1] = ambigTabulator[k2];
}

cd2 = cd2%1;

fprintf (stdout, "\n\nAmbiguities:\n\n");

for (h=k-1; h>=0; h=h-1)
{
	v  = cd2[h][0];
	k2 = cd1[v];
	
	fprintf (stdout, k2);
	k2 = ambigList [k2];
	if (k2[2]+k2[3] == 0)
	{
		fprintf (stdout, " Disallowed "); 
	}
	else
	{
		fprintf (stdout, "            "); 
	}
	fprintf (stdout, Format(cd2[h][1],10,0), " instances\n");
}
