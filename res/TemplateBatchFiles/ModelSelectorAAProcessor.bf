/****************************************************

ModelSelectorAAProcessor.bf

Sergei L Kosakovsky Pond (sergeilkp@mac.com)
October 2007

This HyPhy batch file is to be used to process
the output for GA codon model selection runs.

****************************************************/


ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "Utility" + DIRECTORY_SEPARATOR + "GrabBag.bf");
/* need this include to process file paths */

ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "Utility" + DIRECTORY_SEPARATOR + "DescriptiveStatistics.bf");


SetDialogPrompt ("Locate the .details file produced by the GA:");
fscanf 			(PROMPT_FOR_FILE, "Raw", detailedFile);

pathParts  	= splitFilePath(LAST_FILE_PATH);
dir_prefix	= pathParts["DIRECTORY"];
file_name   = pathParts["FILENAME"];


/* splitFilePath is defined in GrabBag.bf */

/*
sscanf (detailedFile, "String", aminoacidOrdering); 
*/

aminoacidOrdering = "ACDEFGHIKLMNPQRSTVWY";
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
						   
stateVectorDimension = 190;

/*-----------------------------------------------------------------------------*/

fprintf (stdout, "[PHASE 1.] Parsing the .details file\n");

/* first pass to determine the best AIC score
   and to identify the criterion for inclusion in the credible set of models
*/

while (!END_OF_FILE)
{
	sscanf (detailedFile, "Number,NMatrix,NMatrix",modelAIC,classMatrix, rateInfo);
	
	/* check for the end-of-file condition */
	modelRates = Rows(rateInfo);
	if (modelRates == 0)
	{
		break;
	}
	
	byRateClass [modelRates]   = byRateClass[modelRates] + 1;
	
	if (modelAIC < bestScore)
	{
		bestScore   = modelAIC;	
		bestModelID = modelCount;
		bestRates	= modelRates;
	}
	
	modelCount 				   = modelCount + 1;
	
}

ratesRead   = avlKeysToMatrix (byRateClass);
scoreCutoff = bestScore-2*Log(evidenceRCut);
fprintf (stdout, "\tRead ", modelCount, " models\n",
				 "\tBest score = ", bestScore," achieved with ", bestRates, " rates\n",
				 "\tScore cutoff of ", scoreCutoff, " to be included in the credible set at ", evidenceRCut, " level \n",
				 "\t", Abs (byRateClass), " different rate counts\n"
				 );
				 
				 
for (k=0; k<Abs(byRateClass); k=k+1)
{
	fprintf (stdout, "\t", Format(byRateClass[ratesRead[k]],8,0), " models with ", Format(ratesRead[k],4,0), " rate classes\n");
}

/* hack to reset END_OF_FILE */

sscanf (aminoacidOrdering,"String",k);

/*-----------------------------------------------------------------------------*/

fprintf (stdout, "[PHASE 2.] Building the set of credible models\n");

checkAIC = {modelCount,1};
counter2 = 0;

modelCount = 0;
byRateClass		= {}; 

while (!END_OF_FILE)
{
	sscanf (detailedFile, "Number,NMatrix,NMatrix",modelAIC, classMatrix, rateInfo);
	
	/* check for the end-of-file condition */
	modelRates = Rows(rateInfo);
	if (modelRates == 0)
	{
		break;
	}
	checkAIC[counter2] = modelAIC;
	
	counter2 = counter2 + 1;
	if (modelAIC <= scoreCutoff) /* model is in the credible set */
	{
		if (modelAIC == bestScore)
		{
			bestModelID = modelCount;
			bestNRates  = rateInfo;
		}
		
		AICScores[modelCount] 		= modelAIC;
		qMatrix						= classMatrix;
		
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
				qMatrix[k][k2]  = rateInfo[qMatrix[k][k2]];
				compressedM[k3] = qMatrix[k][k2];
				k3 = k3+1;
				qMatrix[k2][k] = qMatrix[k][k2];
			}
		}

		compressedRates[modelCount]		    	= compressedM;
		
		if (modelAIC == bestScore)
		{
			bestModelByRate = {};
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
				 
topNCount = Min (topNCount, modelCount);
				 
for (k=0; k<Abs(byRateClass); k=k+1)
{
	fprintf (stdout, "\t", Format(byRateClass[ratesRead[k]],8,0), " models with ", Format(ratesRead[k],4,0), " rate classes\n");
}

fprintf (stdout, "\tThe best model has the following rate class assignments\n");
for (k=0; k<Abs(bestModelByRate); k=k+1)
{
	fprintf (stdout, "\tRate ", Format(bestNRates[k],8,4), " with ", Format(bestModelByRate[bestNRates[k]],4,0), " substitutions\n");
}

bestMatrix 	  = classMatrices[bestModelID];

DEFAULT_FILE_SAVE_NAME = file_name + ".best_matrix";
ExportAMatrix ("bestMatrix","","Save the best-scoring matrix to:");

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

DEFAULT_FILE_SAVE_NAME = file_name + ".ma_matrix";
ExportAMatrix ("modelAveragedM","","Save the numerical model-averaged rate matrix to:");

/*-----------------------------------------------------------------------------*/

fprintf (stdout, "[PHASE 4.] Computing variability in rate assignments\n");

DEFAULT_FILE_SAVE_NAME = file_name + "_reliability.csv";
SerDialogPrompt ("Save the top-N model .csv to:");
fprintf (PROMPT_FOR_FILE,CLEAR_FILE,KEEP_OPEN,"AA1,AA2,Mean,Median,2.5%,97.5%,Variance,COV,Unreliable\n");

for (aa1=0;aa1<20;aa1=aa1+1)
{
	for (aa2=aa1+1;aa2<20;aa2=aa2+1)
	{
		rateInfo = {modelCount,1};
		for (k=0; k<modelCount; k=k+1)
		{
			rateInfo[k] = (rateMatrices[k])[aa1][aa2];
		}
		stats = GatherDescriptiveStats (rateInfo);
		fprintf (LAST_FILE_PATH,"\n",aminoacidOrdering[aa1],",",
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

fprintf (LAST_FILE_PATH,CLOSE_FILE);


/*-----------------------------------------------------------------------------*/

fprintf (stdout, "[PHASE 5.] Computing a 4-category approximation\n");

binByType = {stateVectorDimension,4};

for (k=0; k<modelCount; k=k+1)
{
	compressedM = compressedRates[k];
	myAW		= akaikeWeights  [k];
	for (aa1=0;aa1<stateVectorDimension;aa1=aa1+1)
	{
		if (compressedM[aa1] <0.1)
		{
			binByType [aa1][0] = binByType [aa1][0] + myAW;
		}
		else{
			if (compressedM[aa1] <0.5)
			{
				binByType [aa1][1] = binByType [aa1][1] + myAW;
			}
			else{
				if (compressedM[aa1] <1.25)
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

k = 0;
for (aa1=0;aa1<20;aa1=aa1+1)
{
	for (aa2=aa1+1;aa2<20;aa2=aa2+1)
	{
		maxV = 0; maxI = 0;
		for (k2=0; k2<4; k2=k2+1){if (binByType[k][k2]>maxV) {maxV = binByType[k][k2]; maxI = k2;}} 
		globbedMatrix[aa1][aa2] = maxI; 
		globbedMatrix[aa2][aa1] = maxI; 
		k=k+1;
	}
}

DEFAULT_FILE_SAVE_NAME = file_name + ".4_matrix";
ExportAMatrix ("globbedMatrix","","Save the 4-class matrix to:");

/*-----------------------------------------------------------------------------*/

fprintf (stdout, "[PHASE 6.] Computing the list of rates derived from top ",topNCount," models\n");

topNList = ({modelCount,2}["_MATRIX_ELEMENT_ROW_*(_MATRIX_ELEMENT_COLUMN_)+AICScores[_MATRIX_ELEMENT_ROW_]*(_MATRIX_ELEMENT_COLUMN_==0)"])%0;

DEFAULT_FILE_SAVE_NAME = file_name + "_topModels.csv";
SerDialogPrompt ("Save the top-N model .csv to:");
fprintf (PROMPT_FOR_FILE,CLEAR_FILE,KEEP_OPEN,"AA1,AA2,Rate,Model_Rank\n");

for (aa1=0;aa1<20;aa1=aa1+1)
{
	for (aa2=aa1+1;aa2<20;aa2=aa2+1)
	{
		for (k=0; k<topNCount; k=k+1)
		{
			fprintf (LAST_FILE_PATH,"\n",aminoacidOrdering[aa1],",",aminoacidOrdering[aa2],",",(rateMatrices[topNList[k][1]])[aa1][aa2],",",k+1);
		}
	}
}
fprintf (LAST_FILE_PATH,CLOSE_FILE);

return 0;


/*-----------------------------------------------------------------------------*/

fprintf (stdout, "[PHASE 7.] Finding stable structures in the set of credible models\n");

consensusStructure = {stateVectorDimension,stateVectorDimension}; 
					/* only one step-substitutions */

bestMatrix = compressedRates[bestModelID];

consensusStructure = consensusStructure["bestMatrix[_MATRIX_ELEMENT_ROW_]==bestMatrix[_MATRIX_ELEMENT_COLUMN_]"];


for (k=0; k<modelCount; k=k+1)
{
	/*matchCount 	   = ({1,stateVectorDimension}["1"]*(consensusStructure*{stateVectorDimension,1}["1"]))[0];
	matchCount 	   = (matchCount-stateVectorDimension)/2;
	fprintf (stdout, k, ":", matchCount, "\n");*/
	if (k!=bestModelID)
	{
		bestMatrix = compressedRates[k];
		consensusStructure = consensusStructure$consensusStructure["bestMatrix[_MATRIX_ELEMENT_ROW_]==bestMatrix[_MATRIX_ELEMENT_COLUMN_]"];
	}
}

matchCount 	   = ({1,stateVectorDimension}["1"]*(consensusStructure*{stateVectorDimension,1}["1"]))[0];
matchCount 	   = (matchCount-stateVectorDimension)/2;


/*-----------------------------------------------------------------------------*/

function ExportAMatrix (theMatrix&, fileName, theprompt)
{
	if (Abs(fileName) == 0)
	{
		SetDialogPrompt (theprompt);
		fprintf (PROMPT_FOR_FILE, CLEAR_FILE, KEEP_OPEN);
		fileName = LAST_FILE_PATH;
	}
	else
	{
		fprintf (fileName, CLEAR_FILE, KEEP_OPEN);
	}
	
	fprintf (fileName, aminoacidOrdering, "\n", theMatrix, CLOSE_FILE);
	return 0;
	
}
