/****************************************************

ModelSelectorCodonProcessor.bf

Sergei L Kosakovsky Pond (sergeilkp@mac.com)
May 2007-Aug 2009

This HyPhy batch file is to be used to process
the output for GA codon model selection runs.

****************************************************/

_BF_RESULTS = {};


ExecuteAFile  ("ModelSelectorCodon.ibf");



/****************************************************/

autoSave = 1; /* do not prompt the user for file names */

ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "Utility" + DIRECTORY_SEPARATOR + "DBTools.ibf");
/* need this include to generate output tables file paths */

ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "Utility" + DIRECTORY_SEPARATOR + "GrabBag.bf");
/* need this include to process file paths */

ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "Utility" + DIRECTORY_SEPARATOR + "DescriptiveStatistics.bf");

ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "_MSCOneStep.ibf");
/* need this include to determine all 1-to-1 substitutions */

ExecuteAFile(HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Utility"+DIRECTORY_SEPARATOR+"PostScript.bf");
	
lcapFile = HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"EmpiricalAA" + DIRECTORY_SEPARATOR+"LCAP";
fscanf (lcapFile,"NMatrix,NMatrix,NMatrix,NMatrix,NMatrix", lcap1,lcap2,lcap3,lcap4,lcap5);


/*SetDialogPrompt   ("Locate the database file to create/add to:");
resultsDatabase = _openCacheDB ("");

if (_TableExists (resultsDatabase,"DATASET"))
{*/

	SetDialogPrompt ("Locate the .details file produced by the GA:");
	fscanf 			(PROMPT_FOR_FILE, "Raw", detailedFile);

	pathParts		= splitFilePath(LAST_FILE_PATH);
	dir_prefix		= pathParts["DIRECTORY"];
	file_name		= pathParts["FILENAME"];


	/* splitFilePath is defined in GrabBag.bf */

	sscanf (detailedFile, "String", aminoacidOrdering); 
	/* read the 0-19 to amino-acid code indexing string */

	modelCount 		= 0;  /* how many models have been read */
	credibleModels	= 0;  /* how many models are credible   */

	classMatrices 	= {}; /* stores class allocation matrices for each model  */
	rateMatrices  	= {}; /* matrices storing c-AIC and rate estimates */
	AICScores		= {}; /* the c-AIC score for every model */
	rateClasses		= {}; /* how many rate classes are there for each model */
	compressedRates = {};

	topNCount		= 10;
	bestScore		= 1e100;
	bestModelID		= 0;
	bestRates		= 0;
	byRateClass		= {}; 
	evidenceRCut	= 0.01; /* smallest evidence ratio against the best model to be 
							   included in the credible set */

	/*-----------------------------------------------------------------------------*/

	fprintf (stdout, "[PHASE 1.] Parsing the .details file\n");

	/* first pass to determine the best AIC score
	   and to identify the criterion for inclusion in the credible set of models
	*/
	
	byRateBestIC = {};

	ScoreProfile = {};
	
	while (!END_OF_FILE)
	{
		sscanf (detailedFile, "NMatrix,NMatrix",classMatrix, rateInfo);
		
		/* check for the end-of-file condition */
		rateCount = Columns(rateInfo);
		if (rateCount == 0)
		{
			break;
		}
		
		modelAIC				   = rateInfo[2];
		modelRates				   = Columns(rateInfo)-3;
		if (byRateBestIC[modelRates] == 0)
		{
			byRateBestIC[modelRates] = rateInfo[0];
		}
		else
		{
			byRateBestIC[modelRates] = Max(rateInfo[0],byRateBestIC[modelRates]);
		}
		byRateClass [modelRates]   = byRateClass[modelRates] + 1;
		
		if (modelAIC < bestScore)
		{
			bestScore   = modelAIC;	
			bestModelID = modelCount;
			bestRates	= modelRates;
		}
		
		ScoreProfile[modelCount]   = modelAIC;
		modelCount 				   = modelCount + 1;
		
	}
	
	ScoreProfileM = avlToMatrix ("ScoreProfile"); ScoreProfile = 0;

	ratesRead   = avlKeysToMatrix (byRateClass);
	scoreCutoff = bestScore-2*Log(evidenceRCut);
	fprintf (stdout, "\tRead ", modelCount, " models\n",
					 "\tBest score = ", Format(bestScore,15,4), " achieved with ", bestRates, " rates\n",
					 "\tScore cutoff of ", scoreCutoff, " to be included in the credible set at ", evidenceRCut, " level \n",
					 "\t", Abs (byRateClass), " different rate counts\n"
					 );
					 
	
	_BF_RESULTS ["TOTAL MODEL COUNT"]				= modelCount;
	_BF_RESULTS ["BEST MODEL SCORE"]				= Format(bestScore,15,4);
	_BF_RESULTS ["BEST MODEL RATE COUNT"]			= bestRates;
	_BF_RESULTS ["SCORE CUTOFF"]					= evidenceRCut;
					 				 
	for (k=0; k<Abs(byRateClass); k=k+1)
	{
		fprintf (stdout, "\t", Format(byRateClass[ratesRead[k]],8,0), " models with ", Format(ratesRead[k],4,0), " rate classes\n");
	}
	
	sscanf (detailedFile, REWIND, "String", k); 

	/*-----------------------------------------------------------------------------*/

	fprintf (stdout, "[PHASE 2.] Building the set of credible models\n");

	modelCount = 0;
	byRateClass		= {}; 

	while (!END_OF_FILE)
	{
		sscanf (detailedFile, "NMatrix,NMatrix",classMatrix, rateInfo);
		
		/* check for the end-of-file condition */
		rateCount = Columns(rateInfo);
		if (rateCount == 0)
		{
			break;
		}
		modelAIC				   		= rateInfo[2];
		if (modelAIC <= scoreCutoff) /* model is in the credible set */
		{
			if (modelAIC == bestScore)
			{
				bestModelID = modelCount;
				bestNRates  = rateInfo;
			}
			
			AICScores[modelCount] 		= modelAIC;
			qMatrix						= classMatrix;
			
			modelRates				   = Columns(rateInfo)-3;
			byRateClass [modelRates]   = byRateClass[modelRates] + 1;
			compressedM 			   = {stateVectorDimension,1};
			/* stateVectorDimension is the number of valid one-step substitution rates */

			k3 = 0;
			for (k=0; k<20; k=k+1)
			{
				classMatrix[k][k] = -1;
				qMatrix    [k][k] = -1;
				
				for (k2=k+1; k2<20; k2=k2+1)
				{
					if (isOneStepSub[k][k2])
					{
						qMatrix[k][k2]  = rateInfo[3+qMatrix[k][k2]];
						compressedM[k3] = qMatrix[k][k2];
						k3 = k3+1;
					}
					else
					{
						qMatrix[k][k2]     = -1;
						classMatrix[k][k2] = -1;
						classMatrix[k2][k] = -1;
					}
					qMatrix[k2][k] = qMatrix[k][k2];
				}
			}

			compressedRates[modelCount]		    	= compressedM;
			
			if (modelAIC == bestScore)
			{
				bestModelByRate     = {};
				bestCompressedMR    = {stateVectorDimension,1};
				
				k3 = 0;
				for (k=0; k<20; k=k+1)
				{
					for (k2=k+1; k2<20; k2=k2+1)
					{
						if (isOneStepSub[k][k2])
						{
							bestCompressedMR[k3] = classMatrix[k][k2];
							k3 = k3+1;
						}
						qMatrix[k2][k] = qMatrix[k][k2];
					}
				}
					
				for (k=0; k<stateVectorDimension; k=k+1)
				{
					bestModelByRate [compressedM[k]] = bestModelByRate [compressedM[k]] + 1;
				}
				
			}

			classMatrices[modelCount]   			= classMatrix;
			rateMatrices[modelCount]				= qMatrix;
			modelCount 				    			= modelCount + 1;
		}
	}

	ratesRead   = avlKeysToMatrix (byRateClass);
	fprintf (stdout, "\tFound ", modelCount, " credible models\n",
					 "\t", Abs (byRateClass), " different rate counts\n");
					 

	_BF_RESULTS ["CREDIBLE MODELS"]					= modelCount;
	_BF_RESULTS ["ONE STEP SUBSTITUTIONS"]			= stateVectorDimension;
	
	
		
	for (k=0; k<Abs(byRateClass); k=k+1)
	{
		fprintf (stdout, "\t", Format(byRateClass[ratesRead[k]],8,0), " models with ", Format(ratesRead[k],4,0), " rate classes\n");
		fprintf (stdout, "\t Best model score in this class: ", byRateBestIC[ratesRead[k]], "\n");
	}

	fprintf (stdout, "\tThe best model has the following rate class assignments\n");
	for (k=0; k<Abs(bestModelByRate); k=k+1)
	{
		fprintf (stdout, "\tRate ", Format(bestNRates[k+3],8,4), " with ", Format(bestModelByRate[bestNRates[k+3]],4,0), " substitutions\n");
	}

	bestClassMatrix 	  = classMatrices[bestModelID];
	bestRateMatrix		  = rateMatrices [bestModelID];

	DEFAULT_FILE_SAVE_NAME = file_name + ".best_matrix";
	ExportAMatrix ("bestClassMatrix","","Save the best-scoring matrix to:");

	/*-----------------------------------------------------------------------------*/

	fprintf (stdout, "[PHASE 3.] Computing a model averaged numerical matrix\n");

	scalingFactorSum = 1;
	modelAveragedM	 = {20,20};
	akaikeWeights	 = {modelCount,1};

	for (k=0; k<modelCount; k=k+1)
	{
		akaikeWeight = Exp((bestScore-AICScores[k])*0.5);
		modelAveragedM = modelAveragedM + rateMatrices[k] * akaikeWeight;
		akaikeWeights [k] = akaikeWeight;
	}

	modelAveragedM = modelAveragedM * (1/(-modelAveragedM[0][0]));
	akaikeWeights  = akaikeWeights*(1/({1,modelCount}["1"]*akaikeWeights)[0]);


	_BF_RESULTS ["BEST MATRIX"]						= {};
	
	//aminoacidOrdering;
	k3 = 0;
	
	for (k=0; k<20; k=k+1)
	{
		for (k2=k+1; k2<20; k2=k2+1)
		{
			if (isOneStepSub[k][k2])
			{
				thisRecord = {"FROM": aminoacidOrdering[k], "TO": aminoacidOrdering[k2], 
							  "CLASS": bestCompressedMR[k3],
							  "NUMERIC": bestNRates [bestCompressedMR[k3]+3],
							  "AVERAGED": modelAveragedM [k][k2]};
				_BF_RESULTS ["BEST MATRIX"] + thisRecord;			
				k3 = k3+1;
			}
			qMatrix[k2][k] = qMatrix[k][k2];
		}
	}
				
	DEFAULT_FILE_SAVE_NAME = file_name + ".ma_matrix";
	ExportAMatrix ("modelAveragedM","","Save the numerical model-averaged rate matrix to:");

	DEFAULT_FILE_SAVE_NAME = file_name + ".ma_matrix.ps";
	SerDialogPrompt ("Save the PS plot for the model-averaged rate matrixto:");
	if (autoSave)
	{
		fileName = dir_prefix + DEFAULT_FILE_SAVE_NAME;
		fprintf (fileName,CLEAR_FILE);	
	}
	else
	{
		fprintf (PROMPT_FOR_FILE,CLEAR_FILE);
		fileName = LAST_FILE_PATH;
	}
	
	fprintf (fileName, rateMatrixToPS (aminoacidOrdering, modelAveragedM, "model-averaged model"));


	/*-----------------------------------------------------------------------------*/

	polarity = {};
	charge   = {};
	stanfel	 = {};
	
	storeProfile ("charge","RHK",1);
	storeProfile ("charge","DE",-1);
	storeProfile ("charge","ANCQGILMFPSTWYV",0);
	
	storeProfile ("polarity","RNDCQEGHKSTY",1);
	storeProfile ("polarity","AILMFPWV",0);
	
	storeProfile ("stanfel","ACGILMPSTV",1);
	storeProfile ("stanfel","DENQ",2);
	storeProfile ("stanfel","FWY",3);
	storeProfile ("stanfel","HKR",4);

	fprintf (stdout, "[PHASE 4.] Computing rates vs structure assignments\n");

	DEFAULT_FILE_SAVE_NAME = file_name + "_reliability.csv";
	SerDialogPrompt ("Save the top-N model .csv to:");
	if (autoSave)
	{
		fileName = dir_prefix + DEFAULT_FILE_SAVE_NAME;
		fprintf (fileName,CLEAR_FILE,KEEP_OPEN,"AA1,AA2,Structural,Averaged,StanfelChange,PolarityChange,ChargeChange,ChemicalComposition,Polarity,Volume,IsoelectricPoint,Hydropathy");	
	}
	else
	{
		fprintf (PROMPT_FOR_FILE,CLEAR_FILE,KEEP_OPEN,"AA1,AA2,Structural,Averaged,StanfelChange,PolarityChange,ChargeChange,ChemicalComposition,Polarity,Volume,IsoelectricPoint,Hydropathy");
		fileName = LAST_FILE_PATH;
	}
	for (aa1=0;aa1<20;aa1=aa1+1)
	{
		for (aa2=aa1+1;aa2<20;aa2=aa2+1)
		{
			if (isOneStepSub[aa1][aa2])
			{
				stats = GatherDescriptiveStats (rateInfo);
				aal1 = aminoacidOrdering[aa1];
				aal2 = aminoacidOrdering[aa2];
				fprintf (fileName,"\n",aal1,",",
											 aal2,",",
											 bestRateMatrix[aa1][aa2],",",
											 modelAveragedM[aa1][aa2],",",
											 stanfel[aal1] == stanfel[aal2],",",
											 polarity[aal1] == polarity[aal2],",",
											 charge[aal1] == charge[aal2],",",										 
											 Abs(lcap1[aa1][aa2]),",",
											 Abs(lcap2[aa1][aa2]),",",
											 Abs(lcap3[aa1][aa2]),",",
											 Abs(lcap4[aa1][aa2]),",",
											 Abs(lcap5[aa1][aa2])
											 );
			}
		}
	}

	fprintf (fileName,CLOSE_FILE);

	/*fprintf (stdout, "[PHASE 4.] Computing variability in rate assignments\n");

	DEFAULT_FILE_SAVE_NAME = file_name + "_reliability.csv";
	SerDialogPrompt ("Save the top-N model .csv to:");
	if (autoSave)
	{
		fileName = dir_prefix + DEFAULT_FILE_SAVE_NAME;
		fprintf (fileName,CLEAR_FILE,KEEP_OPEN,"AA1,AA2,Mean,Median,2.5%,97.5%,Variance,COV,Unreliable\n");	
	}
	else
	{
		fprintf (PROMPT_FOR_FILE,CLEAR_FILE,KEEP_OPEN,"AA1,AA2,Mean,Median,2.5%,97.5%,Variance,COV,Unreliable\n");
		fileName = LAST_FILE_PATH;
	}
	for (aa1=0;aa1<20;aa1=aa1+1)
	{
		for (aa2=aa1+1;aa2<20;aa2=aa2+1)
		{
			if (isOneStepSub[aa1][aa2])
			{
				rateInfo = {modelCount,1};
				for (k=0; k<modelCount; k=k+1)
				{
					rateInfo[k] = (rateMatrices[k])[aa1][aa2];
				}
				stats = GatherDescriptiveStats (rateInfo);
				fprintf (fileName,"\n",aminoacidOrdering[aa1],",",
											 aminoacidOrdering[aa2],",",
											 stats["Mean"],",",
											 stats["Median"],",",
											 stats["2.5%"],",",
											 stats["97.5%"],",",
											 stats["Variance"],",",
											 stats["COV"],",",
											 stats["COV"] > 0.1
											 );
			}
		}
	}

	fprintf (fileName,CLOSE_FILE);*/


	/*-----------------------------------------------------------------------------*/

	fprintf (stdout, "[PHASE 5.] Computing a 4-category [0-0.1,0.1-0.5,0.5-1,1+] approximation\n");

	binByType = {stateVectorDimension,4};

	for (k=0; k<modelCount; k=k+1)
	{
		compressedM = compressedRates[k];
		myAW		= akaikeWeights  [k];
		for (aa1=0;aa1<stateVectorDimension;aa1=aa1+1)
		{
			if (compressedM[aa1] < 0.1)
			{
				binByType [aa1][0] = binByType [aa1][0] + myAW;
			}
			else{
				if (compressedM[aa1] < 0.5)
				{
					binByType [aa1][1] = binByType [aa1][1] + myAW;
				}
				else{
					if (compressedM[aa1] < 1.00)
					{
						binByType [aa1][2] = binByType [aa1][2] + myAW;
					}
					else
					{
						binByType [aa1][3] = binByType [aa1][3] + myAW;
					}
				}
			}
		}
	}

	globbedMatrix = {20,20};
	shortToLong   = {stateVectorDimension,1};

	k = 0;
	for (aa1=0;aa1<20;aa1=aa1+1)
	{
		for (aa2=aa1+1;aa2<20;aa2=aa2+1)
		{
			if (isOneStepSub[aa1][aa2])
			{
				maxV = 0; maxI = 0;
				for (k2=0; k2<4; k2=k2+1){if (binByType[k][k2]>maxV) {maxV = binByType[k][k2]; maxI = k2;}} 
				globbedMatrix[aa1][aa2] = maxI; 
				globbedMatrix[aa2][aa1] = maxI; 
				shortToLong [k] = aminoacidOrdering[aa1]+aminoacidOrdering[aa2];
				k=k+1;
			}
		}
	}

	DEFAULT_FILE_SAVE_NAME = file_name + ".4_matrix";
	ExportAMatrix ("globbedMatrix","","Save the 4-class matrix to:");

	/*-----------------------------------------------------------------------------*/

	fprintf (stdout, "[PHASE 6.] Computing the list of rates derived from top ",topNCount," models\n");

	topNList = ({modelCount,2}["_MATRIX_ELEMENT_ROW_*(_MATRIX_ELEMENT_COLUMN_)+AICScores[_MATRIX_ELEMENT_ROW_]*(_MATRIX_ELEMENT_COLUMN_==0)"])%0;

	DEFAULT_FILE_SAVE_NAME = file_name + "_topModels.csv";
	SetDialogPrompt ("Save the top-N model .csv to:");
	if (autoSave)
	{
		fileName = dir_prefix + DEFAULT_FILE_SAVE_NAME;
		fprintf (fileName,CLEAR_FILE,KEEP_OPEN,"AA1,AA2,Rate,Model_Rank\n");
	}
	else
	{
		fprintf (PROMPT_FOR_FILE,CLEAR_FILE,KEEP_OPEN,"AA1,AA2,Rate,Model_Rank\n");
		fileName = LAST_FILE_PATH;
	}
	for (aa1=0;aa1<20;aa1=aa1+1)
	{
		for (aa2=aa1+1;aa2<20;aa2=aa2+1)
		{
			if (isOneStepSub[aa1][aa2])
			{
				for (k=0; k<topNCount; k=k+1)
				{
					fprintf (fileName,"\n",aminoacidOrdering[aa1],",",aminoacidOrdering[aa2],",",(rateMatrices[topNList[k][1]])[aa1][aa2],",",k+1);
				}
			}
		}
	}
	fprintf (fileName,CLOSE_FILE);
	
	DEFAULT_FILE_SAVE_NAME = file_name + "_init.models";
	SetDialogPrompt ("Save the top-N models in a format suitable for seeding another GA run:");
	if (autoSave)
	{
		fileName = dir_prefix + DEFAULT_FILE_SAVE_NAME;
		fprintf (fileName,CLEAR_FILE,KEEP_OPEN);
	}
	else
	{
		fprintf (PROMPT_FOR_FILE,CLEAR_FILE,KEEP_OPEN);	
		fileName = LAST_FILE_PATH;
	}
	for (k3=0; k3<topNCount; k3=k3+1)
	{
		fprintf (fileName,ConvertMatrixToStateVector(classMatrices[topNList[k3][1]]),"\n");
	}
	fprintf (fileName,CLOSE_FILE);


	/*-----------------------------------------------------------------------------*/

	fprintf (stdout, "[PHASE 7.] Finding stable structures in the set of credible models\n");

	consensusStructure = {stateVectorDimension,stateVectorDimension}; 
						/* only one step-substitutions */

	bestMatrix = compressedRates[bestModelID];

	consensusStructure = consensusStructure["bestMatrix[_MATRIX_ELEMENT_ROW_]==bestMatrix[_MATRIX_ELEMENT_COLUMN_]"]*akaikeWeights  [bestModelID];


	for (k=0; k<modelCount; k=k+1)
	{
		/*matchCount 	   = ({1,stateVectorDimension}["1"]*(consensusStructure*{stateVectorDimension,1}["1"]))[0];
		matchCount 	   = (matchCount-stateVectorDimension)/2;
		fprintf (stdout, k, ":", matchCount, "\n");*/
		
		if (k!=bestModelID)
		{
			akaikeWeight       = Exp((bestScore-AICScores[k])*0.5);
			bestMatrix         = compressedRates[k];
			consensusStructure = consensusStructure + (consensusStructure["bestMatrix[_MATRIX_ELEMENT_ROW_]==bestMatrix[_MATRIX_ELEMENT_COLUMN_]"])*akaikeWeights  [k];
		}
	}
	
	clusterSupport = {stateVectorDimension,1};
	
	for (k = 0; k < stateVectorDimension; k=k+1)
	{
		myRate = bestCompressedMR[k];
		for (k2 = 0; k2 < stateVectorDimension; k2 = k2 + 1)
		{
			if (k2 != k && myRate == bestCompressedMR[k2])
			{
				clusterSupport[k] = clusterSupport[k] + consensusStructure[k][k2];
			}
		}
		bmr = (bestModelByRate[bestNRates[myRate+3]]);
		if (bmr>1)
		{
			clusterSupport[k] = clusterSupport[k]/(bestModelByRate[bestNRates[myRate+3]]-1);
		}
	}

	

	graphs = {};
	aaByNode ={};
	for (k = 0; k < bestRates; k = k+1)
	{
		graphs[k] = ""; graphs[k] * 128; graphs[k] * ("subgraph cluster_G"+k+"{label = \"Rate = " + bestNRates[k+3]+ "\" fontsize = \"24\" penwidth=\"0\"; \n\n");
		aaByNode[k] = {};
	}

	/*for (k=0; k<20; k=k+1)
	{
		graphString * ("\nsubgraph cluster_"+k+"{\n\""+aminoacidOrdering[k]+"\"[ shape = \"record\" label = \" <r0> ");
		for (k2 = 1; k2 < bestRates; k2=k2+1)
		{
			graphString * ("| <r" + k2 + ">");
		}
		graphString * ( "\"]; label = \""+aminoacidOrdering[k]+"\";} \n");	
	}
	
	for (k=0; k<20; k=k+1)
	{	
		for (k2=k+1; k2<20; k2 = k2+1)
		{
			ratePlug = bestClassMatrix[k][k2];
			if (ratePlug >=0)
			{
				graphString * ( aminoacidOrdering[k] + ": r" + ratePlug  +" -- " + aminoacidOrdering[k2] + ": r"+ratePlug + ";\n");
			}
		}
	}*/
	
	
	k3 = 0;
	for (k=0; k<20; k=k+1)
	{	
		for (k2=k+1; k2<20; k2 = k2+1)
		{
			ratePlug = bestClassMatrix[k][k2];
			if (ratePlug >=0)
			{
				if ((aaByNode[ratePlug])[k] == 0)
				{
					(aaByNode[ratePlug])[k] = 1; 
					graphs[ratePlug] * ("" + aminoacidOrdering[k] + ratePlug + 
						"[label = \"" + 
						aminoacidOrdering[k] + 
						"\"" +
						residueStyle(stanfel[aminoacidOrdering[k]],polarity[aminoacidOrdering[k]],charge[aminoacidOrdering[k]]) +
						"];");
				}
				if ((aaByNode[ratePlug])[k2] == 0)
				{
					(aaByNode[ratePlug])[k2] = 1; 
					graphs[ratePlug] * ("" + aminoacidOrdering[k2] + ratePlug + 
						"[label = \"" + 
						aminoacidOrdering[k2] + 
						"\"" +
						residueStyle(stanfel[aminoacidOrdering[k2]],polarity[aminoacidOrdering[k2]],charge[aminoacidOrdering[k2]]) +
						"];");
				}
				if (clusterSupport[k3] < 0.5)
				{
					style = "dotted";
				}
				else
				{
					if (clusterSupport[k3] < .9)
					{
						style = "dashed";
					}				
					else
					{
						style = "solid";
					}
				}
				
				graphs[ratePlug] * ( aminoacidOrdering[k]  + ratePlug + " -- " + aminoacidOrdering[k2] + ratePlug+ " [style = \""+style+"\" label = \""+Format(modelAveragedM[k][k2],5,2)+"\"];\n");
				
				k3 = k3+1;
			}
		}
	}
	
	sortBestRates = {bestRates,2}["_MATRIX_ELEMENT_ROW_"];
	for (k=0; k<bestRates; k=k+1)
	{
		sortBestRates[k][1]=bestNRates[k+3];
	}
	sortBestRates = sortBestRates % 1;
	
	for (k = 0; k < bestRates; k = k+1)
	{
		graphs[k] * "}\n\n"; graphs[k] * 0; 
	}

	SetDialogPrompt ("Save the GraphViz .dot file to:");
	DEFAULT_FILE_SAVE_NAME = file_name + ".dot";
	if (autoSave)
	{
		fileName = dir_prefix + DEFAULT_FILE_SAVE_NAME;
		fprintf (fileName,CLEAR_FILE,KEEP_OPEN);
	}
	else
	{
		fprintf (PROMPT_FOR_FILE,CLEAR_FILE,KEEP_OPEN);	
		fileName = LAST_FILE_PATH;
	}
	
	fprintf (fileName, "graph G{remincross = \"true\" rankdir=\"LR\" size = \"12,12\";\nnode [fontsize=14,width=".4", height=".4", margin=0];graph[fontsize=14];");
	
	for (k = 0; k < bestRates; k = k+1)
	{
		fprintf (fileName, graphs[sortBestRates[k][0]],"\n"); 
	}
	
	fprintf (fileName, "\n}", CLOSE_FILE);


/*}

else
{
	fprintf (stdout, "[ERROR] The database file must contain the 'DATASET' table storing alignment records"); 
}

_closeCacheDB (resultsDatabase);
*/

/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/

function ExportAMatrix				(theMatrix&, fileName, theprompt)
{
	if (Abs(fileName) == 0)
	{
		if (autoSave)
		{
			fileName = dir_prefix+DEFAULT_FILE_SAVE_NAME;
			fprintf (fileName, CLEAR_FILE);
		}
		else
		{
			SetDialogPrompt (theprompt);
			fprintf (PROMPT_FOR_FILE, CLEAR_FILE, KEEP_OPEN);
			fileName = LAST_FILE_PATH;
		}
	}
	else
	{
		fprintf (fileName, CLEAR_FILE, KEEP_OPEN);
	}
	
	fprintf (fileName, aminoacidOrdering, "\n", theMatrix, CLOSE_FILE);
	return 0;
	
}

function ConvertMatrixToStateVector (theMatrix)
{
	vector = {1,stateVectorDimension};
	for (k = 0; k < 20; k=k+1)
	{
		for (k2 = k+1; k2<20; k2=k2+1)
		{
			if (matrixToVectorMap[k][k2]>=0)
			{
				vector[matrixToVectorMap[k][k2]] = theMatrix[k][k2];
			}
		}
	}
	return vector;
}

function 	storeProfile (recp&, string, value)
{
	for (h=0; h<Abs(string); h=h+1)
	{
		recp[string[h]] = value;
	}
	return 0;
}

function 	residueStyle (s,p,c)
{
	if (s == 1)
	{
		 color = "#B03060";
		 labelcolor = "#FFFFFF";
	}
	if (s == 2)
	{
		 color = "#00A86B";
		 labelcolor = "#000000";
	}
	if (s == 3)
	{
		 color = "#FF8C00";
		 labelcolor = "#000000";
	}
	if (s == 4)
	{
		 color = "#4B0082";
		 labelcolor = "#FFFFFF";
	}
	
	if (p == 0)
	{
		if (c == 0)
		{
			shape = "rect";
		}
		else
		{
			if (c == 1)
			{
				shape = "trapezium";
			}
			else
			{
				shape = "invtrapezium";
			}
		}
	}
	else
	{
		if (c == 0)
		{
			shape = "diamond";
		}
		else
		{
			if (c == 1)
			{
				shape = "triangle";
			}
			else
			{
				shape = "invtriangle";
			}
		}	
	}
	
	return " style=\"filled\" color=\"" + color + "\" fontcolor=\"" + labelcolor + "\" shape=\"" + shape + "\"";
}


function rateMatrixToPS (chars,rateMatrix, title)
{
	
	codon_order = "FLIMVSPTAYHQNKDECWRG";
	codon_idx   = {};
	
	for (k=0; k<=20; k=k+1)
	{
		codon_idx[codon_order[k]] = k;
	}
	
	psFigure = "";
	psFigure * 8192;

	psFigure * _HYPSPageHeader (445,500, "Protein Rate Matrix Plot");
	psFigure * "\n";
	psFigure * _HYPSSetFont ("Times-Roman", 12);
	psFigure * "\n";
	psFigure * _HYPSTextCommands(0);

	psFigure * "/box {\n0 20 rlineto \n20 0 rlineto \n0 -20 rlineto \nclosepath } def\n 0 30 translate\n";
	
	offset = 24;
	
	maxVal = Max(rateMatrix,0);
	
	
	
	charsToCodonOrder 	= {};
	
	for (k=0; k<=20; k=k+1)
	{
		charsToCodonOrder[k] = codon_idx[chars[k]];
	}
	
	reordering = {};
	reordering["A"] = 0;
	reordering["C"] = 1;
	reordering["G"] = 2;
	reordering["I"] = 3;
	reordering["L"] = 4;
	reordering["M"] = 5;
	reordering["P"] = 6;
	reordering["S"] = 7;
	reordering["T"] = 8;
	reordering["V"] = 9;
	reordering["D"] = 10;
	reordering["E"] = 11;
	reordering["N"] = 12;
	reordering["Q"] = 13;
	reordering["F"] = 14;
	reordering["W"] = 15;
	reordering["Y"] = 16;
	reordering["H"] = 17;
	reordering["K"] = 18;
	reordering["R"] = 19;
	
	rk = Rows (reordering);
	
	for (h=0; h<20; h=h+1)
	{
		h2 = reordering[chars[h]];
		label = rk[h];
		psFigure * ("0 setgray\n10 "+(offset+(19-h)*20+6)+" \n("+label+") centertext\n");
		psFigure * ("0 setgray\n"+(offset+410)+" "+(offset+(19-h)*20+6)+" ("+label+") centertext\n");
		psFigure * ("0 setgray\n"+(offset+h*20+10)+" 10 ("+label+") centertext\n");
		psFigure * ("0 setgray\n"+(offset+h*20+10)+" "+(405+offset)+" ("+label+") centertext\n");
		for (v=0; v<20; v=v+1)
		{
			if (h!=v)
			{
				v2 = reordering[chars[v]];
				label2 = chars[v];
				psFigure * ("newpath\n"+(offset+v2*20)+" "+(offset+(19-h2)*20)+" moveto\n");
				greyColor = 1-rateMatrix[h][v]/maxVal;
				psFigure * (""+greyColor+" setgray\nbox fill\n");
				if (isOneStepSub[charsToCodonOrder[h]][charsToCodonOrder[v]])
				{
					if (greyColor>0.5)
					{
						psFigure * ("0 setgray\n");
					}
					else
					{
						psFigure * ("1 setgray\n");				
					}
					psFigure * ("\nnewpath\n"+(offset+v2*20+10)+" "+(offset+(19-h2)*20+10)+" 5 0 360 arc\n\nstroke\nclosepath\n");
				}
			}
		}
	}
	
	psFigure * ("\n"+offset+" "+offset+" translate\n0 setgray\nnewpath\n0 0 moveto\n0 400 lineto\n400 400 lineto\n400 0 lineto\n0 0 lineto\nstroke\nclosepath");
	
	psFigure * ("\n\nnewpath\n0 200 moveto\n200 200 lineto\n200 400 lineto\nstroke");
	psFigure * ("\n\nnewpath\n200 120 moveto\n200 200 lineto\n280 200 lineto\n280 120 lineto\nclosepath stroke");
	psFigure * ("\n\nnewpath\n280 120 moveto\n340 120 lineto\n340 60 lineto\n280 60 lineto\nclosepath stroke");
	psFigure * ("\n\nnewpath\n340 60 moveto\n400 60 lineto\n400 0 lineto\n340 0 lineto\nclosepath stroke");
	psFigure * ("\n0 -30 translate\n 225.5 460 (Rate matrix plot for " + title + ") centertext\n 0 -20 translate\n");
	tableText = {3,1};
	tableText [0] = "Shading indicates relative substitution rates (black = max, white = 0)";
	tableText [1] = "Circles show residue pairs that can be exchanged with one nucleotide substitution";
	tableText [2] = "Amino-acids are grouped into 4 Stanfel classification clusters";
	psFigure * _HYPSTextTable (400,30,12,tableText,tableText["0"]);
	psFigure * ("\nshowpage");
	psFigure * 0;
	return psFigure;
}


/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/
