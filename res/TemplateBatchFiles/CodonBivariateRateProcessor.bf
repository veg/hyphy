/*---------------------------------------------------------------------------------------------------------------------------------------------------*/

function echoCatVar (ps,values)
{
	DD 			= Rows(ps);
	EE 			= 0.0;
	sampleVar 	= 0.0;
	
	for (k=0; k<DD; k=k+1)
	{
		EE = ps[k]*values[k]+EE;
		sampleVar = sampleVar+ps[k]*values[k]*values[k];
	}
		
	sampleVar = sampleVar-EE*EE;
	
	fprintf  (stdout,  "Mean     = ",EE, 
					 "\nVariance = ",sampleVar,
					 "\nCOV      = ", Sqrt(sampleVar)/EE,"\n");
					 
	for (k=0; k<DD; k=k+1)
	{
		fprintf (stdout,"\nRate[",Format(k,0,0),"]=",Format(values[k],12,8), " (weight=", 
						  Format(ps[k],9,7),")");
	}
	return EE;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------------*/

function echoCovariance (ps,values1,values2)
{
	DD 			= Rows(ps);
	EE 			= 0.0;
	EE2			= 0.0;
	sampleVar 	= 0.0;
	sampleVar2 	= 0.0;
	
	for (k=0; k<DD; k=k+1)
	{
		EE  = ps[k]*values1[k]+EE;
		EE2 = ps[k]*values2[k]+EE2;
		sampleVar = sampleVar+ps[k]*values1[k]*values1[k];
		sampleVar2 = sampleVar2+ps[k]*values2[k]*values2[k];
	}
		
	sampleVar = sampleVar-EE*EE;
	sampleVar2 = sampleVar2-EE2*EE2;
	
	cov  = 0;
	cov2 = 0;
	for (k=0; k<DD; k=k+1)
	{
		cov  = cov  + ps[k]*(values1[k]-EE);
		cov2 = cov2 + ps[k]*(values2[k]-EE2);
	}

	return cov*cov2/Sqrt(sampleVar*sampleVar2);
}

/*---------------------------------------------------------------------------------------------------------------------------------------------------*/

ExecuteAFile(HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"chooseGeneticCode.def");
ExecuteAFile(HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Utility"+DIRECTORY_SEPARATOR+"PS_Plotters.bf");

SetDialogPrompt		("Choose a model fit:");
ExecuteAFile		(PROMPT_FOR_FILE);
basePath		    = LAST_FILE_PATH;
GetInformation		(vars,"^P_[0-9]+");
rateCount			= Columns (vars)+1;

GetString			(lfInfo, lf, -1);
fileCount			= Columns(lfInfo["Trees"])/rateCount;

fprintf (stdout, "\nLoaded a fit on ", fileCount, " data sets with ", rateCount, " rates\n");

ps 			 = {rateCount,1};
rateInfo	 = {rateCount,4};
rateInfo	 = rateInfo["1"];

for (mi = 0; mi < rateCount-1; mi=mi+1)
{
	ExecuteCommands ("ps["+mi+"]="+"P_"+(mi+1)+";");
}

for (mi = 0; mi < rateCount; mi=mi+1)
{
	ExecuteCommands ("rateInfo[mi][1]="+"S_"+mi+"/c_scale;");
	ExecuteCommands ("rateInfo[mi][2]="+"NS_"+mi+"/c_scale;");
	rateInfo[mi][3] = rateInfo[mi][2]-rateInfo[mi][1];
}


for (mi=0; mi<rateCount-1; mi=mi+1)
{
	for (mi2 = 0; mi2 < mi; mi2=mi2+1)
	{
		rateInfo[mi][0] = rateInfo[mi][0] * (1-ps[mi2]);
	}
	rateInfo[mi][0] = rateInfo[mi][0] * ps[mi];
}


for (mi2 = 0; mi2 < mi; mi2=mi2+1)
{
	rateInfo[mi][0] = rateInfo[mi][0] * (1-ps[mi2]);
}

fprintf (stdout, "\n\n1).Synonymous rates:\n\n");
echoCatVar (rateInfo[-1][0],rateInfo[-1][1]);
fprintf (stdout, "\n\n2).Non-synonymous rates:\n\n");
echoCatVar (rateInfo[-1][0],rateInfo[-1][2]);

columnHeaders = {{"P","dS","dN","dN-dS"}};
OpenWindow (CHARTWINDOW,{{"Rates"}
		{"columnHeaders"}
		{"rateInfo"}
		{"Scatterplot"}
		{"dS"}
		{"dN"}
		{"dS"}
		{""}
		{"dN"}
		{"0"}
		{""}
		{"-1;-1"}
		{"10;1.309;0.785398"}
		{"Times:12:0;Times:10:0;Times:12:2"}
		{"0;0;16777215;0;0;0;6579300;11842740;13158600;14474460;0;3947580;16777215;16711680;6845928;16771158;2984993;9199669;7018159;1460610;16748822;11184810;14173291"}
		{"16,0,0"}
		},
		"749;663;70;70");
		
		

		
ConstructCategoryMatrix (cm, lf, COMPLETE);

site_count 		= Columns (cm)/rateCount;
posteriorProbs  = {site_count, rateCount};
columnHeaders   = {1,rateCount};

for (rate_enumerator = 0; rate_enumerator < rateCount; rate_enumerator = rate_enumerator + 1)
{
	columnHeaders[rate_enumerator] = "Rate " + (rate_enumerator+1);
}

for (site_enumerator = 0; site_enumerator < site_count; site_enumerator = site_enumerator + 1)
{
	sum = 0; 
	
	smallestScaler = 1e100;
	
	for (rate_enumerator = 0; rate_enumerator < rateCount; rate_enumerator = rate_enumerator + 1)
	{
		smallestScaler = Min(smallestScaler,cm.site_scalers[rate_enumerator*site_count+site_enumerator]);
	}
	
	for (rate_enumerator = 0; rate_enumerator < rateCount; rate_enumerator = rate_enumerator + 1)
	{
		v = cm[rate_enumerator*site_count+site_enumerator] * rateInfo[rate_enumerator][0] * Exp(cm.log_scale_multiplier*(smallestScaler-cm.site_scalers[rate_enumerator*site_count+site_enumerator]));
		posteriorProbs[site_enumerator] [rate_enumerator]= v;
		sum = sum + v;
	}

	for (rate_enumerator = 0; rate_enumerator < rateCount; rate_enumerator = rate_enumerator + 1)
	{
		posteriorProbs[site_enumerator] [rate_enumerator]= posteriorProbs[site_enumerator] [rate_enumerator]/sum;
	}
}

distrInfo = "";
distrInfo * 128;
distrInfo * "dN_minus_dS";
for (k = 0; k <= mi; k = k+1)
{
	distrInfo * (":" + rateInfo[k][0]);
}
for (k = 0; k <= mi; k = k+1)
{
	distrInfo * (":" + rateInfo[k][3]);
}
distrInfo * 0;

{"SmallCodon_part_Categ:0.25:0.25:0.25:0.25:0.00116849:0.0498505:0.447572:3.50141"}

OpenWindow (DISTRIBUTIONWINDOW,{{"Posteriors"}
		{"columnHeaders"}
		{"posteriorProbs"}
		{"None"}
		{""}
		{""}
		{""}
		{""}
		{""}
		{"0"}
		{""}
		{"-1;-1"}
		{"10;1.309;0.785398"}
		{"Times:12:0;Times:10:0;Times:12:2"}
		{"0;0;16777215;0;0;0;6579300;11842740;13158600;14474460;0;3947580;16777215;16711680;6845928;16771158;2984993;9199669;7018159;1460610;16748822;11184810;14173291"}
		{"16,0,0"}
		{distrInfo}		
		},
		"749;663;200;200");


nonStopCount = 0;

for (h1 = 0; h1<64; h1=h1+1)
{
	if (_Genetic_Code[h1]!=10)
	{
		nonStopCount = nonStopCount + 1;
	}
}

/* make syn and non-syn template matrices */

ExecuteAFile(HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Distances"+DIRECTORY_SEPARATOR+"CodonTools.def");

for (fileID = 1; fileID <= fileCount; fileID = fileID+1)
{
	fprintf (stdout, "\n\nWorking on file ", fileID, "\n\n");
	ExecuteCommands ("codonCount=filteredData_"+fileID+".sites;");

	sSites  = 0;
	nsSites = 0;

	synM    = {nonStopCount,nonStopCount};
	nonSynM = {nonStopCount,nonStopCount};

	vertOnes = {nonStopCount,1};
	horOnes  = {1,nonStopCount};

	for (h1 = 0; h1<nonStopCount; h1=h1+1)
	{
		vertOnes [h1] = 1;
		horOnes  [h1] = 1;
	}

	hShift = 0;
	for (h1 = 0; h1 < 64; h1=h1+1)
	{
		gc1 = _Genetic_Code[h1];
		if (gc1 == 10)
		{
			hShift = hShift+1;
		}
		else
		{
			sSites = sSites   + codonCount * _S_NS_POSITIONS_[0][h1] * vectorOfFrequencies[h1-hShift];
			nsSites = nsSites + codonCount * _S_NS_POSITIONS_[1][h1] * vectorOfFrequencies[h1-hShift];
			
			vShift = hShift;
			for (v1 = h1+1; v1 < 64; v1=v1+1)
			{
				gc2 = _Genetic_Code[v1];
				if (gc2 == 10)
				{
					vShift = vShift + 1;
				}
				else
				{
					if (gc1 == gc2)
					{
						synM [h1-hShift][v1-vShift] = vectorOfFrequencies[h1-hShift];
						synM [v1-vShift][h1-hShift] = vectorOfFrequencies[v1-vShift];
					}
					else
					{
						nonSynM [h1-hShift][v1-vShift] = vectorOfFrequencies[h1-hShift];
						nonSynM [v1-vShift][h1-hShift] = vectorOfFrequencies[v1-vShift];
					}
				}
			}
		}
	}


	ExecuteAFile(HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TreeTools.ibf");

	fprintf (stdout, "\nTotal nucleotide sites :", codonCount*3,
					 "\nSynonymous  sites      :", sSites, 
					 "\nNonsynonymous  sites   :", nsSites, "\n");
					 
	sSites  = codonCount/sSites;
	nsSites = codonCount/nsSites;

	ExecuteCommands ("branchNames = BranchName (tree_"+fileID+"_0,-1);");
	T 			= Columns    (branchNames);

	synSubsAVL = {};
	dSAVL	   = {};
	nsSubsAVL  = {};
	dNAVL	   = {};

	for (treeCounter = 0; treeCounter < rateCount; treeCounter = treeCounter + 1)
	{
		for (h1=0; h1 < T-1; h1=h1+1)
		{
			abn = branchNames[h1];
			ExecuteCommands("GetInformation (aRateMx, tree_"+fileID+"_"+treeCounter+"."+abn+");");
			synSubs  = (horOnes*(aRateMx$synM))*vertOnes;
			nsynSubs = (horOnes*(aRateMx$nonSynM))*vertOnes;
			synSubs = synSubs[0]/3;
			nsynSubs = nsynSubs[0]/3;
			
			synSubsAVL[abn] = synSubsAVL[abn] + synSubs*rateInfo[treeCounter][0];
			nsSubsAVL [abn] = nsSubsAVL [abn] + nsynSubs*rateInfo[treeCounter][0];
			dSAVL[abn]	    = dSAVL[abn] + synSubs *sSites*rateInfo[treeCounter][0];
			dNAVL[abn]	    = dNAVL[abn] + nsynSubs*nsSites*rateInfo[treeCounter][0];
		}
	}

	ExecuteCommands ("treeAVL = tree_"+fileID+"_0^0");

	synTreeString 		= PostOrderAVL2StringDistances (treeAVL, synSubsAVL); 
	nonSynTreeString	= PostOrderAVL2StringDistances (treeAVL, nsSubsAVL);
	dSTreeString 		= PostOrderAVL2StringDistances (treeAVL, dSAVL); 
	dNTreeString	    = PostOrderAVL2StringDistances (treeAVL, dNAVL);


	/*fprintf (stdout, "\nE[Syn subs/nucleotide site] tree: \n\t",     synTreeString, 	   "\n");
	fprintf (stdout, "\nE[Non-syn subs/nucleotide site] tree: \n\t", nonSynTreeString, "\n");
	fprintf (stdout, "\ndS tree: \n\t", dSTreeString, "\n");
	fprintf (stdout, "\ndN tree: \n\t", dNTreeString, "\n");*/

	UseModel (USE_NO_MODEL);

	ExecuteCommands ("Tree 	synSubsTree_"+fileID+"= synTreeString;Tree	nonsynSubsTree_"+fileID+" 	= nonSynTreeString;Tree 	dSTree_"+fileID+" 			= dSTreeString;Tree	dNTree_"+fileID+" = dNTreeString;");

	ProcessATree ("synSubsTree_"+fileID);
	ProcessATree ("nonsynSubsTree_"+fileID);
	ProcessATree ("dSTree_"+fileID);
	ProcessATree ("dNTree_"+fileID);
	
	mxTreeSpec  = {5,1};

	mxTreeSpec [0] = "nonsynSubsTree_"+fileID;
	mxTreeSpec [3] = "";
	mxTreeSpec [4] = "Inferred_Tree."+nodeName;
	mxTreeSpec [1] = "8211";
	mxTreeSpec [2] = "";


	OpenWindow (TREEWINDOW, mxTreeSpec,"(SCREEN_WIDTH-50)/2;(SCREEN_HEIGHT-50)/2;10;40");

	mxTreeSpec [0] = "synSubsTree_"+fileID;
	OpenWindow (TREEWINDOW, mxTreeSpec,"(SCREEN_WIDTH-50)/2;(SCREEN_HEIGHT-50)/2;30+(SCREEN_WIDTH-30)/2;40");

	mxTreeSpec [0] = "dSTree_"+fileID;
	OpenWindow (TREEWINDOW, mxTreeSpec,"(SCREEN_WIDTH-50)/2;(SCREEN_HEIGHT-50)/2;10;45+(SCREEN_HEIGHT-50)/2");

	mxTreeSpec [0] = "dNTree_"+fileID;
	OpenWindow (TREEWINDOW, mxTreeSpec,"(SCREEN_WIDTH-50)/2;(SCREEN_HEIGHT-50)/2;30+(SCREEN_WIDTH-30)/2;45+(SCREEN_HEIGHT-50)/2");
}

logRates = {rateCount,3};
epsilon  = Exp(-4);
for (mi = 0; mi < rateCount; mi=mi+1)
{
	logRates[mi][0] = Min(Log(rateInfo[mi][1]+epsilon),4);
	logRates[mi][1] = Min(Log(rateInfo[mi][2]+epsilon),4);
	logRates[mi][2] = rateInfo[mi][0];
}
		
psCode = ScaledDensityPlot ("logRates", {{-4,4}{-4,4}}, "Courier", {{400,400,14,40}},
								"_dNdSDensityPlot",
								{{"","log (alpha)", "log (beta)"}}, 1, 1);

psFile = basePath + ".ps";
fprintf (psFile, CLEAR_FILE, psCode);


/*---------------------------------------------------------*/

function ProcessATree (treeName)
{
	ExecuteCommands ("treeAVL2 = "+treeName + " ^ 0;leafCount=TipCount("+treeName+");"); 
	
	multFactors = {};
	for (k=1; k<Abs(treeAVL2); k=k+1)
	{
		aNode = treeAVL2[k];
		aNodeName = aNode["Name"];
		parentIndex = aNode["Parent"];
		k2 = Abs(aNode["Children"]);
		if (k2)
		{
			currentDepth = aNode["Below"];
			multFactors[aNodeName] = currentDepth;		
			if (parentIndex > 0)
			{
				pInfo = treeAVL2[parentIndex];
				pInfo ["Below"] = pInfo ["Below"] + currentDepth;
				treeAVL2[parentIndex] = pInfo;
			}
		}
		else
		{
			multFactors[aNodeName] = 1;
			pInfo = treeAVL2[parentIndex];
			pInfo ["Below"] = pInfo ["Below"] + 1;
			treeAVL2[parentIndex] = pInfo;
		}
		
	}

	pKeys 			= Rows(multFactors);

	for (k=0; k<Columns(pKeys); k=k+1)
	{
		aNodeName = pKeys[k];
		multFactors[aNodeName] = multFactors[aNodeName] * (leafCount-multFactors[aNodeName]);
	}

	divInfo 		=	 computeTotalDivergence (treeName);
	pInfo 			= 	2*divInfo[0]/leafCount/(leafCount-1);
	currentDepth	= 	divInfo[1]/(Abs(treeAVL2)-2);
	
	fprintf (stdout, "Mean pairwise divergence for ",treeName, " is ", pInfo, 	   "\n");
	fprintf (stdout, "Mean branch length for ",      treeName, " is ", currentDepth, "\n");
	return 0;
}
	
/*---------------------------------------------------------*/

function computeTotalDivergence (treeID)
{
	ExecuteCommands ("bNames = BranchName   ("+treeID+",-1);");
	ExecuteCommands ("bLen   = BranchLength ("+treeID+",-1);");
	
	sum  = 0;
	sum2 = 0;
	
	for (k=0; k<Columns(bNames); k=k+1)
	{
		aNodeName = bNames[k];
		sum  = sum + bLen[k]*multFactors[aNodeName];
		sum2 = sum2 + bLen[k];
	}	
	return {{sum,sum2}};
}
