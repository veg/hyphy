startWithRateClasses	= 2;
maxRateClasses			= 10;
rateClassesCount 		= startWithRateClasses; 
autoStepFlag	 		= 1;

modelSpawnPrefix = "";
modelSpawnPrefix * 128;
modelSpawnSuffix = "";
modelSpawnSuffix * 128;

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

produceOffspring		= MPI_NODE_COUNT-1;
populationSize  		= 2*produceOffspring;
incestDistance  		= 0;
generationCount		  	= 5000;
maxSampleTries			= populationSize*10;
mutationThreshhold		= 0.0001;
mutationProb			= 0.15;
mutationProbDecrease	= 0.95;
annealingPhase			= 100;
SHORT_MPI_RETURN		= 1;
totalSampleCounter		= 0;
localMutationRate		= 0.05;
localMutationInterval	= 20;

stoppingCriterion		= 30;
sampleCount				= 0;
familyControlSize		= produceOffspring$6;

predef   				= {};
MasterList				= {};
verboseFlag				= 0;
stateVectorDimension    = 190;

predef   				= {};



/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function ExportAMatrix2 (fileName, rateMatrix, aic_score, rateAVL)
{
	fprintf (fileName, aic_score, "\n{");
	for (kr=0; kr<Rows (rateMatrix); kr=kr+1)
	{
		fprintf (fileName, "{", rateAVL["NSR"+rateMatrix[kr]], "}\n");
	}
	fprintf (fileName, "}\n");
	return 0;
}
/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function ExportAMatrix (fileName, rateMatrix, doClear)
{
	outString = "";
	outString * 256;
	outString * "branchAssignments={};\n";
	for (nodeCounter = 0; nodeCounter < stateVectorDimension; nodeCounter = nodeCounter+1)
	{
		outString * ("branchAssignments[\""+branchNames[nodeCounter]+"\"] = " + rateMatrix[nodeCounter] + ";\n");
	}
	outString * 0;
	if (doClear)
	{
		fprintf (fileName, CLEAR_FILE,outString);
	}
	else
	{
		fprintf (fileName,"\n",outString);
	}
	return 0;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function MatrixToString (rateMatrix)
{
	outString = "";
	outString * 256;
	for (h=0; h<stateVectorDimension; h=h+1)
	{
		outString * (","+rateMatrix[h]);				
	}
	outString * 0;
	return outString;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function CleanUpMPI (dummy)
{
	if (MPI_NODE_COUNT>1)
	{
		while (1)
		{
			for (nodeCounter = 0; nodeCounter < MPI_NODE_COUNT-1; nodeCounter = nodeCounter+1)
			{
				if (MPINodeState[nodeCounter][0]==1)
				{
					fromNode = ReceiveJobs (0,0);
					break;	
				}
			}
			if (nodeCounter == MPI_NODE_COUNT-1)
			{
				break;
			}
		}			
	}
	return 0;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function ReceiveJobs (sendOrNot, ji)
{
	if (MPI_NODE_COUNT>1)
	{
		MPIReceive (-1, fromNode, result_String);
		mji = MPINodeState[fromNode-1][1];
		mdf	= MPINodeState[fromNode-1][2];
		
		if (sendOrNot)
		{
			MPISend (fromNode,msg2Send);
			MPINodeState[fromNode-1][1] = ji;			
			MPINodeState[fromNode-1][2] = modelDF;			
		}
		else
		{
			MPINodeState[fromNode-1][0] = 0;
			MPINodeState[fromNode-1][1] = -1;		
		}
		ExecuteCommands (result_String);
		myDF	  = lf_MLES[1][1] + baseParams - 1;
		myAIC 	  = 2*(lf_MLES[1][0]-myDF*sampleCount/(sampleCount-myDF-1));
		myLFScore = lf_MLES[1][0];
		/*fprintf (stdout, myLFScore, " ", lf_MLES[1][1], " ", myAIC, "\n");*/
		ji = mji;
	}
	else
	{
		myDF	  = modelDF + baseParams;
		myAIC 	  = 2*(lf_MLES[1][0]-myDF*sampleCount/(sampleCount-myDF-1));
		myLFScore = res[1][0];
	}
	
	if (resultProcessingContext==0)
	{
		sortedScores[ji][0] = myAIC;
	}
	else
	{
		intermediateProbs[ji][0] = myAIC;	
	}
	if (ji>=0)
	{
		if (resultProcessingContext==0)
		{
			ExportAMatrix2 (detailedResultPath,StringToMatrix(currentPopulation[ji]),myAIC,lf_MLE_VALUES);
			sortedBP = MatrixToString (currentPopulation[ji]);
			MasterList [sortedBP] = myAIC;
		}
		else
		{
			ExportAMatrix2 (detailedResultPath,StringToMatrix(children[ji-populationSize]),myAIC,lf_MLE_VALUES);	
			sortedBP = MatrixToString (children[ji-populationSize]);
			MasterList [sortedBP] = myAIC;
		}
		
		/*fprintf (fpath,"\n{{",myLFScore,",",-myAIC);
		if (MPI_NODE_COUNT > 1)
		{
			for (spoolC = 0; spoolC < lf_MLES[1][1]; spoolC = spoolC + 1)
			{
				aVal = lf_MLE_VALUES["NSR"+spoolC];
				fprintf (fpath,",",aVal);
			}
		}		
		else
		{
			for (spoolC = 0; spoolC < lf_MLES[1][1]; spoolC = spoolC + 1)
			{
				execCommand = "fprintf(fpath,\",\",NSR"+spoolC+");";
				ExecuteCommands (execCommand);
			}
		}
		fprintf (fpath,"\n}}\n");*/

		totalSampleCounter = totalSampleCounter + 1;
	}
	return fromNode-1;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function RunASample (modelDF, jobIndex)
{
	sampleString = MatrixToString (cString);
	myAIC = MasterList[sampleString];
	
	if (myAIC < (-0.1))
	{
		if (resultProcessingContext==0)
		{
			sortedScores[jobIndex][0] = myAIC;
		}
		else
		{
			intermediateProbs[jobIndex][0] = myAIC;	
		}			
		return 0;
	}
		
	msg2Send = "";
	msg2Send * 128;
	msg2Send * ("presetBranchParameters = "+currentBLEstimates+";bClasses="+cString+";");
	
	msg2Send * modelSpawnPrefix;
	msg2Send * ("AC:="+AC+";\n");
	msg2Send * ("AT:="+AT+";\n");
	msg2Send * ("CG:="+CG+";\n");
	msg2Send * ("CT:="+CT+";\n");
	msg2Send * ("GT:="+GT+";\n");
	msg2Send * modelSpawnSuffix;
	msg2Send * 0;
	
	
	if ((MPI_NODE_COUNT>1) && (jobIndex>=0))
	{
		for (mpiNode = 0; mpiNode < MPI_NODE_COUNT-1; mpiNode = mpiNode+1)
		{
			if (MPINodeState[mpiNode][0]==0)
			{
				break;	
			}
		}
		if (mpiNode==MPI_NODE_COUNT-1)
		{
			mpiNode = ReceiveJobs (1,jobIndex);
		}
		else
		{
			MPISend (mpiNode+1,msg2Send);
			MPINodeState[mpiNode][0] = 1;
			MPINodeState[mpiNode][1] = jobIndex;
			MPINodeState[mpiNode][2] = modelDF;
		}
	}
	else
	{
		Optimize (res,lf);
		if (jobIndex>=0)
		{
			mpiNode = ReceiveJobs (1, jobIndex);
		}
		else
		{
			myAIC 	  = 2*(res[1][0]-(baseParams+modelDF)*sampleCount/(sampleCount-baseParams-modelDF-1));
		}
	}
	return 0;
}


/*---------------------------------------------------------------------------------------------------------------------------------------------*/


function BuildCodonFrequencies (obsF)
{
	PIStop = 1.0;
	result = {ModelMatrixDimension,1};
	hshift = 0;

	for (h=0; h<64; h=h+1)
	{
		first = h$16;
		second = h%16$4;
		third = h%4;
		if (_Genetic_Code[h]==10) 
		{
			hshift = hshift+1;
			PIStop = PIStop-obsF[first][0]*obsF[second][1]*obsF[third][2];
			continue; 
		}
		result[h-hshift][0]=obsF[first][0]*obsF[second][1]*obsF[third][2];
	}
	return result*(1.0/PIStop);
}



/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function SpawnRandomString (clsCnt)
{
	rModel = {stateVectorDimension,1};
	for (h=0; h<stateVectorDimension; h=h+1)
	{
		rModel[h] = Random(0,clsCnt)$1;
	}
	
	return MakeStringCanonical(rModel,clsCnt);
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function MakeStringCanonical (randomModel, classCount)
{
	compressedString = {classCount,1};
	for (h=0; h<stateVectorDimension; h=h+1)
	{
		v = randomModel[h];
		compressedString[v] = 1;
	}
	compressedString[0] = 0;
	for (h=1; h<classCount; h=h+1)
	{
		compressedString[h] = compressedString[h]+compressedString[h-1];
	}
	for (h=0; h<stateVectorDimension; h=h+1)
	{
		v = randomModel[h];
		v = compressedString[v];
		randomModel[h] = v;
	}
	v = compressedString[classCount-1]+1;
	if (v>1)
	{
		sortedOrder = {v,2};
		for (h=0; h<v; h=h+1)
		{
			sortedOrder[h][0] = -1;
		}
		cc = 0;
		for (h=0; h<stateVectorDimension; h=h+1)
		{
			hshift = randomModel[h];
			if (sortedOrder[hshift][0] < 0)
			{
				sortedOrder[hshift][0] = cc;
				sortedOrder[hshift][1] = hshift;
				cc = cc+1;
			}
		}
		sortedOrder = sortedOrder%1;
		for (h=0; h<stateVectorDimension; h=h+1)
		{
			v = randomModel[h];
			randomModel[h] = sortedOrder[v][0];
		}
		
	}
	return randomModel;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function StringToMatrix (stringSpec)
{
	return stringSpec;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function IsChildViable (putativeChild)
{
	sampleString = MatrixToString (putativeChild);
	myAIC = MasterList[sampleString];
	testChild = putativeChild;
	mutPassCount = 1;
	while (myAIC < (-0.1))
	{
		if (verboseFlag)
		{
			fprintf (stdout,"Adjusting the child to avoid a duplicate. Pass ", mutPassCount, "\n");
		}
		
		mutPassCount = mutPassCount + 1;
		
		sampleString = Min(Random(0,stateVectorDimension)$1,stateVectorDimension-1);
		myAIC = testChild[sampleString];
		
		newValue = Random (0,rateClassesCount-0.0000001)$1;
		
		while (newValue == myAIC)
		{
			newValue = Random (0,rateClassesCount-0.0000001)$1;
		}
		
		testChild [sampleString] = newValue;
		sampleString = MatrixToString (testChild);
		myAIC = MasterList[sampleString];
	}
	return testChild;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function UpdateBL (dummy)
{
	bestInd							= sortedScores[populationSize-1][1];
	dNdSBranchClasses 				= StringToMatrix(currentPopulation[bestInd]);
	MG94custom = 0;
	MULTIPLY_BY_FREQS 				= PopulateModelMatrix ("MG94custom", observedFreq);
	Model MG94customModel 			= (MG94custom,vectorOfFrequencies,0);

	AC=AC;
	AT=AT;
	CG=CG;
	CT=CT;
	GT=GT;
	
	Tree givenTree				    = treeString;
	ClearConstraints				(givenTree);
	ReplicateConstraint 			("this1.?.synRate:=3*this2.?.t", givenTree, nucTree); 
	ClearConstraints				(givenTree);
	
	bNames  	 = BranchName (givenTree,-1);
	
	definedOrNot = {};
	for (brCount = 0; brCount < stateVectorDimension; brCount = brCount + 1)
	{
		cbc	= dNdSBranchClasses[brCount];
		if (definedOrNot[cbc] == 0)
		{
			definedOrNot[cbc] = 0;
			ExecuteCommands ("global NSR"+cbc+"=1;");
		}
		tn = bNames[brCount];
		ExecuteCommands ("givenTree." + tn + ".nonSynRate:= NSR" + cbc + "*givenTree." + tn +".synRate;");
	}
	
	LikelihoodFunction 				lf = (filteredData,givenTree);
	
	if (filteredData.species > 10 && filteredData.unique_sites/(MPI_NODE_COUNT) > 20)
	{
		AUTO_PARALLELIZE_OPTIMIZE		= 1;
	}
	if (verboseFlag)
	{
		VERBOSITY_LEVEL 			= 1;
	}
	Optimize 						(res,lf);
	AUTO_PARALLELIZE_OPTIMIZE		= 0;	
	VERBOSITY_LEVEL 				= 0;
	
	AC:=AC__;
	AT:=AT__;
	CG:=CG__;
	CT:=CT__;
	GT:=GT__;

	currentBLEstimates 						= BranchLength (givenTree,-1)*3;
	myDF									= res[1][1] + 9;	
	myAIC 	  								 = 2*(res[1][0]-myDF*sampleCount/(sampleCount-myDF-1));
	fprintf									 (stdout, "\nUpdated BLs\nAICs: ", -myAIC, "\t",  -sortedScores[populationSize-1][0], "\n", lf, "\n");
	sortedScores[populationSize-1][0] 		 = myAIC;
	MasterList[sampleString] 				 = myAIC;
	
	return 0;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

incFileName = HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"chooseGeneticCode.def";
ExecuteCommands  ("#include \""+incFileName+"\";");


SetDialogPrompt 	  		  ("Locate a codon file:");
DataSet 		ds  		 = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter	filteredData = CreateFilter(ds,3,"","",GeneticCodeExclusions);

fscanf			("_BMS_Aux.ibf", "Raw", sampleCount);
ExecuteCommands (sampleCount);
modelSpawnPrefix * sampleCount;
modelSpawnPrefix * ("DataSet 	ds  = ReadDataFile (\"" + LAST_FILE_PATH + "\");\n");
modelSpawnPrefix * ("GeneticCodeExclusions = \"" + GeneticCodeExclusions + "\";\n");
modelSpawnPrefix * ("_Genetic_Code = " + _Genetic_Code + ";\n");
modelSpawnPrefix * ("DataSetFilter	filteredData = CreateFilter(ds,3,\"\",\"\",GeneticCodeExclusions);\n");
modelSpawnPrefix * 0;


fprintf						  (stdout, "Model string for nucleotide biases:");
fscanf						  (stdin,"String",modelDesc);

sampleCount 				 = filteredData.sites;


/*---------------------------------------------------------------------------------------------------------------------------------------------*/


global AC = 1;
global AT = 1;
global CG = 1;
global CT = 1;
global GT = 1;

catCounterAL = {stateVectorDimension,1};

MGCustomRateBiasTerms = {{"AC*","","AT*","CG*","CT*","GT*"}};	

		
paramCount	   = 0;
_nucBiasTerms  = {4,4};

_nucBiasTerms[0][0] = "";


if (modelDesc[0]==modelDesc[1])
{
	MGCustomRateBiasTerms[0] = MGCustomRateBiasTerms[1];
}

_nucBiasTerms[1][0] = MGCustomRateBiasTerms[0];
_nucBiasTerms[0][1] = MGCustomRateBiasTerms[0];
_nucBiasTerms[2][0] = MGCustomRateBiasTerms[1];
_nucBiasTerms[0][2] = MGCustomRateBiasTerms[1];

h = 0;
v = 3;

for (customLoopCounter2=2; customLoopCounter2<6; customLoopCounter2=customLoopCounter2+1)
{
	for (customLoopCounter=0; customLoopCounter<customLoopCounter2; customLoopCounter=customLoopCounter+1)
	{
		if (modelDesc[customLoopCounter]==modelDesc[customLoopCounter2])
		{
			_nucBiasTerms[h][v] = MGCustomRateBiasTerms[customLoopCounter];
			_nucBiasTerms[v][h] = MGCustomRateBiasTerms[customLoopCounter];
			break;
		}
	}
	if (customLoopCounter == customLoopCounter2)
	{
		_nucBiasTerms[h][v] = MGCustomRateBiasTerms[customLoopCounter2];
		_nucBiasTerms[v][h] = MGCustomRateBiasTerms[customLoopCounter2];
	}
	
	v = v+1;
	if (v==4)
	{
		h=h+1;
		v=h+1;
	}
}

_nucRateMatrix = {4,4};

modelSpawnPrefix * ("_nucBiasTerms = {4,4};_nucBiasTerms[0][0] = \"\";\n");

for (customLoopCounter = 0; customLoopCounter < 4; customLoopCounter = customLoopCounter + 1)
{
	for (customLoopCounter2 = customLoopCounter+1; customLoopCounter2 < 4; customLoopCounter2 = customLoopCounter2 + 1)
	{
		ExecuteCommands ("_nucRateMatrix[" + customLoopCounter  + "][" + customLoopCounter2 + "]:=" + _nucBiasTerms[customLoopCounter][customLoopCounter2] + "t;");
		ExecuteCommands ("_nucRateMatrix[" + customLoopCounter2 + "][" + customLoopCounter  + "]:=" + _nucBiasTerms[customLoopCounter2][customLoopCounter] + "t;");
		modelSpawnPrefix * ("_nucBiasTerms[" + customLoopCounter + "][" + customLoopCounter2 + "] = \"" + _nucBiasTerms[customLoopCounter][customLoopCounter2] + "\";");
		modelSpawnPrefix * ("_nucBiasTerms[" + customLoopCounter2 + "][" + customLoopCounter + "] = \"" + _nucBiasTerms[customLoopCounter2][customLoopCounter] + "\";");
	}
}


DataSetFilter	nucFilter = CreateFilter(ds,1);
HarvestFrequencies (nucFreq,nucFilter,1,1,1);
Model 	nucModel = (_nucRateMatrix,nucFreq,1);

incFileName = HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"queryTree.bf";
ExecuteAFile (incFileName);

Tree 	nucTree = treeString;

branchNames 		 = BranchName (nucTree,-1);
stateVectorDimension = Columns (branchNames)-1;

fprintf (stdout, "\nA total of ", stateVectorDimension, " branches\n");
fprintf (stdout, "\nFitting a nucleotide model to approximate branch lengths...\n");

LikelihoodFunction lf = (nucFilter,nucTree);
Optimize (nuc_res,lf);

currentBLEstimates   = BranchLength (nucTree,-1) * 3;

fprintf (stdout, "\n", lf, "\n");

HarvestFrequencies (observedFreq,filteredData,3,1,1);

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

SetDialogPrompt ("Save the best model to:");
fprintf (PROMPT_FOR_FILE,CLEAR_FILE);
modelFile = LAST_FILE_PATH;

detailedResultPath = modelFile+".details";
fprintf 		  (detailedResultPath,CLEAR_FILE,KEEP_OPEN,treeString,"\n");

dNdSBranchClasses = {stateVectorDimension,1};

MG94plain = 0;
MULTIPLY_BY_FREQS 				= PopulateModelMatrix ("MG94plain", observedFreq);
vectorOfFrequencies 			= BuildCodonFrequencies (observedFreq);
modelSpawnPrefix 				* ("observedFreq = " + observedFreq + ";\n");
modelSpawnPrefix 				* ("vectorOfFrequencies = " + vectorOfFrequencies + ";\n");
modelSpawnPrefix				* ("global AC = 1; global AT = 1; global CG = 1; global CT = 1; global GT = 1;\n");
modelSpawnPrefix				* 0;



Model MG94plainmodel	 		= (MG94plain,vectorOfFrequencies,0);


/*treeString	   =   DATAFILE_TREE;*/
Tree givenTree =   treeString;

modelSpawnSuffix *	("MULTIPLY_BY_FREQS = PopulateModelMatrix (\"MG94custom\", observedFreq);\n");
modelSpawnSuffix *  ("Model MG94customModel = (MG94custom,vectorOfFrequencies,0);\n");
modelSpawnSuffix *	("Tree givenTree2 = " + treeString + ";\n");
modelSpawnSuffix *	("replicateBranchLengths (0);\nLikelihoodFunction lf = (filteredData,givenTree2);\n");
modelSpawnSuffix *	("Optimize (lf_MLES,lf);\nreturn makeReturnValue(0);\n");
modelSpawnSuffix *  0;

global 		   codonBranchScaler = 1;
ReplicateConstraint ("this1.?.synRate:=3*this2.?.t", givenTree, nucTree); 
ClearConstraints(givenTree);
LikelihoodFunction lf = (filteredData,givenTree);

ChoiceList (matingChoice, "Mating Weights",1,SKIP_NONE,
			"Akaike Weights","Individuals are chosen to reproduce with based on their Akaike weights.",
			"Equiprobable", "All individuals are equally likely to mate.",
			"Rank","Mating probabilities are proportional to the c-AIC rank of the individual",
			"Random","One of the three mating schemes in picked at random for every generation");

if (matingChoice<0)
{
	return 0;
}

ChoiceList (reprFlag, "Mating Representation",1,SKIP_NONE,
			"General","Use n-ary representation (n=number of rate classes).",
			"Binary", "Use binary representation");

if (reprFlag<0)
{
	return 0;
}


currentPopulation  = {};
USE_LAST_RESULTS=1;
if (verboseFlag)
{
	VERBOSITY_LEVEL 				= 1;
}
if (filteredData.species > 10 && filteredData.unique_sites/(MPI_NODE_COUNT) > 20)
{
	AUTO_PARALLELIZE_OPTIMIZE		= 1;
}

global NS0 = 1;
ReplicateConstraint ("this1.?.nonSynRate:=NS0*this2.?.synRate", givenTree, givenTree);

Optimize (res,lf);
AUTO_PARALLELIZE_OPTIMIZE = 0;

currentBLEstimates   = BranchLength (givenTree,-1)*3;

VERBOSITY_LEVEL  = 0;
USE_LAST_RESULTS = 0;



baseParams 		   = res[1][1]+9;
fprintf (stdout, "\nBase model has ", baseParams, " parameters\n", lf);
crapAIC			   = -2*(res[1][0]-baseParams*sampleCount/(sampleCount-baseParams-1));
sortedScores	   = {populationSize,2};
sortedScores[0][0] = -crapAIC;
sortedScores[0][1] = 0;
currentPopulation [0] = {stateVectorDimension,1};

AC:=AC__;
AT:=AT__;
CG:=CG__;
CT:=CT__;
GT:=GT__;


/* INCLUDE THE CHC CORE */

ibfPath = HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "GA_CHC.ibf";


byBPImprovement		  = {};
byBPImprovement[0]    = crapAIC;
bestIndividualOverall = sortedScores[0][0];
currentBEST_IC		  = crapAIC;

maxRateClasses		  = Min(maxRateClasses,stateVectorDimension);

for (currentBPC = startWithRateClasses; currentBPC < maxRateClasses; currentBPC = currentBPC + 1)
{
	rateClassesCount 		= currentBPC;
	fprintf (stdout, "\n\nStarting GA with ", rateClassesCount, " rate classes\n");
	addOnLine = " with " + rateClassesCount + " rate classes.";
	
	if (currentBPC > startWithRateClasses)
	{
		children = {};
		individual = 0;
		
		for (individual=populationSize $ 2; individual<populationSize-1; individual=individual+1)
		{
			cString    = currentPopulation[populationSize-1];
			toggleRate = Random (0,stateVectorDimension-0.00001)$1;
			cString [toggleRate] = currentBPC-1;
			cString 						= MakeStringCanonical(cString,rateClassesCount);
			cString							= IsChildViable (cString);
			sortedScores[individual][1] 	= individual;
			currentPopulation[individual] 	= cString;
			RunASample (compressedString[rateClassesCount-1],individual);
		}
	}
	
	
	mutationProb = 0.15;
	ExecuteAFile (ibfPath);
	
	kf						 = -sortedScores[populationSize-1][0];
	
	if (currentBEST_IC > kf)
	{
		currentBEST_IC = kf;
	}
	else
	{
		break;
	}
}


/* try a consensus matrix */

UpdateBL (0);
ExportAMatrix (modelFile,StringToMatrix(currentPopulation[populationSize-1]),1);
fprintf (detailedResultPath, CLOSE_FILE);
