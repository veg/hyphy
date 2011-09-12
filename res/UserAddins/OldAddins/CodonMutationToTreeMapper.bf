MESSAGE_LOGGING = 0;

branchColors = {};

branchColors [0] = {{255,0,0}}*(1/255); 		/* red */
branchColors [3] = {{64,0,128}}*(1/255); 		/* eggplant */
branchColors [2] = {{0,0,0}};   	/* light gray */
branchColors [1] = {{0,128,0}}*(1/255);   		/* clover */

codonTo3 = {};

codonOffset = 0;

codonTo3[0] = "F";
codonTo3[1] = "L";
codonTo3[2] = "I";
codonTo3[3] = "M";
codonTo3[4] = "V";
codonTo3[5] = "S";
codonTo3[6] = "P";
codonTo3[7] = "T";
codonTo3[8] = "A";
codonTo3[9] = "Y";
codonTo3[10] = "Stop";
codonTo3[11] = "H";
codonTo3[12] = "Q";
codonTo3[13] = "N";
codonTo3[14] = "K";
codonTo3[15] = "D";
codonTo3[16] = "E";
codonTo3[17] = "C";
codonTo3[18] = "W";
codonTo3[19] = "R";
codonTo3[20] = "G";

nucCharacters = "ACGT";

function codeToLetters (codonCode)
{
	return nucCharacters[codonCode$16]+nucCharacters[(codonCode%16)$4]+nucCharacters[codonCode%4];	
}

function setupMapToTree (dummy)
{
	skipCodeSelectionStep = 1;
	dummy = HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"chooseGeneticCode.def";
	ExecuteCommands ("#include \""+dummy+"\";");
	ApplyGeneticCodeTable (0);

	DataSet			ancestralSeqs  = ReconstructAncestors (lf);
	DataSetFilter	filteredDataA  = CreateFilter(ancestralSeqs,3,"","",GeneticCodeExclusions);

	HarvestFrequencies			  (observedCEFV,filteredData,3,3,0);

	stateCharCount = 64;

	for (h=0; h<64; h=h+1)
	{
		if (_Genetic_Code[h] == 10)
		{
			stateCharCount = stateCharCount -1;
		}
	}
	ambChoice = 1;

	seqToBranchMap = {stateCharCount,1};

	hShift = 0;

	_GC_ = {stateCharCount,1};
	correctCode = {stateCharCount,1};

	for (k=0; k<64; k=k+1)
	{
		if (_Genetic_Code[k]==10)
		{
			hShift = hShift+1;
		}
		else
		{
			seqToBranchMap[k-hShift] = observedCEFV[k];
			_GC_[k-hShift] = _Genetic_Code[k];
			correctCode[k-hShift] = k;
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

	GetDataInfo    (dupInfo, filteredData);
	GetDataInfo	   (dupInfoA, filteredDataA);

	matrixTrick  = {1,stateCharCount};
	matrixTrick2 = {1,stateCharCount};

	for (h=Columns(matrixTrick)-1; h>=0; h=h-1)
	{
		matrixTrick  [h] = h;
		matrixTrick2 [h] = 1;
	}

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

	flatTreeRep	  = Abs (givenTree);
	GetInformation (seqStrings, filteredData);
	return 0;
}


function mapSiteToTree (site, mapToAA)
{
	TREE_OUTPUT_OPTIONS = {};
	k = filteredData.species+1;
	
	countSubstitutionsInt  = {};
	countSubstitutionsLeaf = {};
	
	/* first sequence is always the root */
	c1 = dupInfoA[site];
	for (h=1; h<filteredDataA.species; h=h+1)
	{
		p1  = seqToBranchMap[k][0];
		pid = flatTreeRep[p1];
		parentName = branchNames [pid];
		if (Abs(TREE_OUTPUT_OPTIONS[parentName])==0)
		{
			parentName = 0;
		}
		p2  = seqToBranchMap[pid][1]-filteredData.species;
		
		cd1 = codonInfo2[h] [c1];
		cd2 = codonInfo2[p2][c1];
		
		nodeSpec  = {};
		if (mapToAA)
		{
			cd3 = _GC_[cd1];
			cd3 = codonTo3 [cd3];
		}		
		else
		{
			cd3 = codeToLetters(correctCode[cd1]);
		}
		
		if (cd1 == cd2)
		{
			nodeSpec ["TREE_OUTPUT_BRANCH_COLOR"] = branchColors[2];
			nodeSpec ["TREE_OUTPUT_BRANCH_LABEL"] = ""+cd3;
		}
		else
		{
			nodeSpec ["TREE_OUTPUT_BRANCH_LABEL"] = "__FONT_SIZE__ 2 idiv\n__FONT_SIZE__ 3 idiv\nneg\nrmoveto\n("+cd3+") show";
			
			if (Abs(parentName))
			{
				pinfoavl  = TREE_OUTPUT_OPTIONS[parentName];
				pbname    = pinfoavl ["TREE_OUTPUT_BRANCH_LABEL"];
				if (Abs (pbname) < 3)
				{
					pinfoavl ["TREE_OUTPUT_BRANCH_LABEL"] = "__FONT_SIZE__ 2 idiv\n__FONT_SIZE__ 3 idiv\nneg\nrmoveto\n("+pbname+") show";
					TREE_OUTPUT_OPTIONS[parentName] = pinfoavl;
				}
			}
			
			if (_GC_[cd1] == _GC_[cd2])
			{
				nodeSpec ["TREE_OUTPUT_BRANCH_COLOR"] = branchColors[1];
			}
			else
			{
				nodeSpec ["TREE_OUTPUT_BRANCH_COLOR"] = branchColors[0];	
				nodeSpec ["TREE_OUTPUT_BRANCH_DASH"] = {{1,4,1}};		
			}
			aa1 = _GC_[cd1];
			aa1 = codonTo3[aa1];
			aa2 = _GC_[cd2];
			aa2 = codonTo3[aa2];
			cd1 = codeToLetters(correctCode[cd1]);
			cd2 = codeToLetters(correctCode[cd2]);
			aa1 = cd2+"("+aa2+") to " + cd1 + "(" + aa1 + ")";
			countSubstitutionsInt[aa1] = countSubstitutionsInt[aa1] + 1;			
		}
		
		cd1 = branchNames[p1];
		TREE_OUTPUT_OPTIONS[cd1] = nodeSpec;
		
		k=k+1;
	}
	
	/* now do the leaves */
	
	c2 = dupInfo[site];
	for (h=0; h<filteredData.species; h=h+1)
	{
		p1 = seqToBranchMap[h][0];
		pid = flatTreeRep[p1];
		parentName = branchNames [pid];
		if (Abs(TREE_OUTPUT_OPTIONS[parentName])==0)
		{
			parentName = 0;
		}
		p2 = seqToBranchMap[pid][1]-filteredData.species;
		
		cd2 = codonInfo2[p2][c1];
		cd1 = codonInfo[h] [c2];
		

		if (cd1>=0)
		/* no ambiguities */
		{
			nodeSpec  = {};
			if (mapToAA)
			{
				cd3 = _GC_[cd1];
				cd3 = codonTo3 [cd3];
			}		
			else
			{
				cd3 = codeToLetters(correctCode[cd1]);
			}
			
			
			if (cd1 == cd2)
			{
				nodeSpec ["TREE_OUTPUT_BRANCH_COLOR"] = branchColors[2];
				nodeSpec ["TREE_OUTPUT_BRANCH_LABEL"] = "";
			}
			else
			{
				if (Abs(parentName))
				{
					pinfoavl  = TREE_OUTPUT_OPTIONS[parentName];
					pbname    = pinfoavl ["TREE_OUTPUT_BRANCH_LABEL"];
					if (Abs (pbname) < 3)
					{
						pinfoavl ["TREE_OUTPUT_BRANCH_LABEL"] = "__FONT_SIZE__ 2 idiv\n__FONT_SIZE__ 3 idiv\nneg\nrmoveto\n("+pbname+") show";
						TREE_OUTPUT_OPTIONS[parentName] = pinfoavl;
					}
				}
				nodeSpec ["TREE_OUTPUT_BRANCH_LABEL"] = "__FONT_SIZE__ 2 idiv\n__FONT_SIZE__ 3 idiv\nneg\nrmoveto\n("+cd3+") show";
				if (_GC_[cd1] == _GC_[cd2])
				{
					nodeSpec ["TREE_OUTPUT_BRANCH_COLOR"] = branchColors[1];
				}
				else
				{
					nodeSpec ["TREE_OUTPUT_BRANCH_COLOR"] = branchColors[0];			
					nodeSpec ["TREE_OUTPUT_BRANCH_DASH"] = {{1,4,1}};		
				}
				aa1 = _GC_[cd1];
				aa1 = codonTo3[aa1];
				aa2 = _GC_[cd2];
				aa2 = codonTo3[aa2];
				cd1 = codeToLetters(correctCode[cd1]);
				cd2 = codeToLetters(correctCode[cd2]);
				aa1 = cd2+"("+aa2+") to " + cd1 + "(" + aa1 + ")";
				countSubstitutionsLeaf[aa1] = countSubstitutionsLeaf[aa1] + 1;
			}
			
			cd1 = branchNames[p1];
			TREE_OUTPUT_OPTIONS[cd1] = nodeSpec;
		}	
		else
		/* ambiguities here */
		{
			seqString = seqStrings[h];
			GetDataInfo    (ambInfo, filteredData, h, c2);	
			nodeSpec  = {};
			haveSyn = 0;
			haveNS  = 0;
			
			if (Abs(parentName))
			{
				pinfoavl  = TREE_OUTPUT_OPTIONS[parentName];
				pbname    = pinfoavl ["TREE_OUTPUT_BRANCH_LABEL"];
				if (Abs (pbname) < 3)
				{
					pinfoavl ["TREE_OUTPUT_BRANCH_LABEL"] = "__FONT_SIZE__ 2 idiv\n__FONT_SIZE__ 3 idiv\nneg\nrmoveto\n("+pbname+") show";
					TREE_OUTPUT_OPTIONS[parentName] = pinfoavl;
				}
			}
			
			
			ambBName = "";
			alreadyDone = {};
			for (k=0; k<stateCharCount; k=k+1)
			{
				if (ambInfo[k])
				{
					if (k == cd2)
					{
						haveSyn = 1;
					} 
					else
					{
						haveNS = 1;
					}
					aa1 = _GC_[k];
					aa1 = codonTo3[aa1];
					if (alreadyDone[aa1] == 0)
					{	
						alreadyDone[aa1] = 1;
						if (Abs(ambBName))
						{
							ambBName = ambBName+"/"+aa1;
						}	
						else
						{	
							ambBName = aa1;
						}
					}
				}
			}	

			if (mapToAA)
			{
				nodeSpec ["TREE_OUTPUT_BRANCH_LABEL"] = "__FONT_SIZE__ 2 idiv\n__FONT_SIZE__ 3 idiv\nneg\nrmoveto\n("+ambBName+") show";
			}
			else
			{
				nodeSpec ["TREE_OUTPUT_BRANCH_LABEL"] = "__FONT_SIZE__ 2 idiv\n__FONT_SIZE__ 3 idiv\nneg\nrmoveto\n("+seqString[3*site][3*site+2]+" ) show";
			}

			if (haveSyn || haveNS)
			{
				aa2 = _GC_[cd2];
				aa2 = codonTo3[aa2];
				if (Abs(alreadyDone)>1)
				{
					nodeSpec ["TREE_OUTPUT_BRANCH_COLOR"] = branchColors[3];
					aa1 = "AMB";
					nodeSpec ["TREE_OUTPUT_BRANCH_DASH"] = {{1,4,1}};		
				}
				else
				{
					nodeSpec ["TREE_OUTPUT_BRANCH_COLOR"] = branchColors[1];
					aa1 = aa2;
				}
				cd1 = seqString[3*site][3*site+2];
				cd2 = codeToLetters(correctCode[cd2]);
				aa1 = cd2+"("+aa2+") to " + cd1 + "(" + aa1 + ")";
				countSubstitutionsLeaf[aa1] = countSubstitutionsLeaf[aa1] + 1;
			}
			else
			{
				nodeSpec ["TREE_OUTPUT_BRANCH_COLOR"] = branchColors[0];			
			}
			cd1 = branchNames[p1];
			TREE_OUTPUT_OPTIONS[cd1] = nodeSpec;
		}
	}
		
	for (aa2 = 0; aa2 < Columns (branchNames); aa2=aa2+1)
	{
		parentName = branchNames[aa2];
		pinfoavl   = TREE_OUTPUT_OPTIONS[parentName];
		if (Abs(pinfoavl))
		{
			pbname    = pinfoavl ["TREE_OUTPUT_BRANCH_LABEL"];
			if (Abs (pbname) < 3)
			{
				pinfoavl ["TREE_OUTPUT_BRANCH_LABEL"] = "";
				TREE_OUTPUT_OPTIONS[parentName] = pinfoavl;
			}
		}
	}	
		
	psString = PSTreeString (givenTree,"EXPECTED_NUMBER_OF_SUBSTITUTIONS",{{300,600}});
	repMx = {{"setfont"}{"setfont\n3 setlinewidth\n1 setlinecap"}};
	psString = psString ^ repMx;
	repMx = {{"[0-9] scalefont"}{"11 scalefont"}};
	psString = psString ^ repMx;

	outResult = {};
	outResult ["Tree string"] 	  = psString;
	outResult ["Leaf subs"]   	  = countSubstitutionsLeaf;
	outResult ["Internal subs"]   = countSubstitutionsInt;
	
	return outResult;
}

setupMapToTree(0);

while (1)
{
	fprintf (stdout, "0-based site to map (-1 to stop:):");
	fscanf  (stdin, "Number", siteToMap);
	if (siteToMap < 0)
	{
		break;
	}
	z = mapSiteToTree (siteToMap$1, 1);
	SetDialogPrompt ("Write tree to:");
	fprintf (PROMPT_FOR_FILE,CLEAR_FILE,z["Tree string"]);
	fprintf (stdout, "\nLeaf: ", z["Leaf subs"], "\nInternal: ", z["Internal subs"], "\n");
}
