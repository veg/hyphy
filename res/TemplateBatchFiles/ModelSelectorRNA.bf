fprintf (stdout, "\nHow many rate classes?");
fscanf  (stdin,  "Number", rateClassesCount);

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

stoppingCriterion		= 100;
sampleCount			= 0;
familyControlSize		= produceOffspring$6;

predef   = {};

MasterList				= {};
verboseFlag				= 1;
stateVectorDimension    = 120;

predef [0] = {{1,1,1,1,0,0,0,1,0,0,0,1,0,0,0,1,1,0,1,0,0,0,1,0,0,0,1,0,0,1,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,1,1,1,1,1,0,0,0,1,0,0,0,1,1,0,1,0,0,0,1,0,0,1,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,1,1,1,1,1,0,0,0,1,1,0,1,0,0,1,0,0,1,0,0,0,0,1,1,1,1,1,1,1}};

if (rateClassesCount>3)
{
	predef [1] = {{2,2,1,2,0,0,0,2,0,0,0,1,0,0,0,2,1,0,2,0,0,0,1,0,0,0,2,0,0,1,0,0,1,0,0,0,2,0,0,0,1,0,0,0,0,1,0,0,0,3,0,0,0,1,2,1,2,2,0,0,0,1,0,0,0,1,2,0,1,0,0,0,2,0,0,1,0,0,1,0,0,0,3,0,0,0,0,1,0,0,0,2,1,2,1,1,0,0,0,1,3,0,1,0,0,1,0,0,1,0,0,0,0,1,1,3,1,1,2,1}};
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function ExportAMatrix (fileName, rateMatrix,doClear)
{
	outString = "";
	outString * 256;
	outString * "\n{";
	for (h=0; h<16; h=h+1)
	{
		outString * "\n{";
		for (v=0; v<h; v=v+1)
		{
			if (v)
			{
				outString * (","+rateMatrix[h][v]);
			}
			else
			{
				outString * (""+rateMatrix[h][v]);				
			}
		}
		if (v)
		{
			outString * ",*";
		}
		else
		{
			outString * "*";				
		}
		for (v=h+1; v<16; v=v+1)
		{
			outString * (","+rateMatrix[h][v]);				
		}
		outString * "}";
	}
	outString * "\n}";
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
	for (h=0; h<120; h=h+1)
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
			MPISend (fromNode,lf);
			MPINodeState[fromNode-1][1] = ji;			
			MPINodeState[fromNode-1][2] = modelDF;			
		}
		else
		{
			MPINodeState[fromNode-1][0] = 0;
			MPINodeState[fromNode-1][1] = -1;		
		}
		ExecuteCommands (result_String);
		myDF	  = mdf;
		myAIC 	  = 2*(lf_MLES[1][0]-(baseParams+mdf)*sampleCount/(sampleCount-baseParams-mdf-1));
		myLFScore = lf_MLES[1][0];
		ji = mji;
	}
	else
	{
		myDF	  = modelDF;
		myAIC 	  = 2*(res[1][0]-(baseParams+myDF)*sampleCount/(sampleCount-baseParams-myDF-1));
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
		fpath = modelFile+".details";
		if (resultProcessingContext==0)
		{
			ExportAMatrix (fpath,StringToMatrix(currentPopulation[ji]),0);
			sortedBP = MatrixToString (currentPopulation[ji]);
			MasterList [sortedBP] = myAIC;
		}
		else
		{
			ExportAMatrix (fpath,StringToMatrix(children[ji-populationSize]),0);	
			sortedBP = MatrixToString (children[ji-populationSize]);
			MasterList [sortedBP] = myAIC;
		}
		fprintf (fpath,"\n{{",myLFScore,",",-myAIC);
		if (MPI_NODE_COUNT > 1)
		{
			/*fprintf (stdout, lf_MLE_VALUES);*/
			for (spoolC = 0; spoolC <= myDF; spoolC = spoolC + 1)
			{
				aVal = lf_MLE_VALUES["R_"+spoolC];
				fprintf (fpath,",",aVal);
			}
		}		
		else
		{
			for (spoolC = 0; spoolC <= myDF; spoolC = spoolC + 1)
			{
				execCommand = "fprintf(fpath,\",\",R_"+spoolC+");";
				ExecuteCommands (execCommand);
			}
		}
		fprintf (fpath,"\n}}\n");
		fpath = modelFile+".AIC";
		if (totalSampleCounter)
		{
			fprintf (fpath,"\n",-myAIC);
		}
		else
		{
			fprintf (fpath,CLEAR_FILE,"\n",-myAIC);		
		}
		totalSampleCounter = totalSampleCounter + 1;
	}
	return fromNode-1;
}

/*outCounter = 0;*/
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
	}
	else
	{
		AAModelMatrix 					= 0;
		MULTIPLY_BY_FREQS 				= PopulateModelMatrix ("AAModelMatrix", vectorOfFrequencies);
		Model AAModel	     			= (AAModelMatrix,vectorOfFrequencies,1);
		Tree givenTree2 = treeString;
		if (optimizeBLFlag < 1)
		{
			ReplicateConstraint("this1.?.t:=Null_Scaler__*this2.?.t__/branchLengthScaler",givenTree2,givenTree);
		}
		LikelihoodFunction lf = (filteredData,givenTree2);
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
				MPISend (mpiNode+1,lf);
				/*fName = "spool."+outCounter;
				fprintf (fName,CLEAR_FILE,MPI_LAST_SENT_MSG);
				outCounter = outCounter + 1;*/
				MPINodeState[mpiNode][0] = 1;
				MPINodeState[mpiNode][1] = jobIndex;
				MPINodeState[mpiNode][2] = modelDF;
			}
		}
		else
		{
			/*LIKELIHOOD_FUNCTION_OUTPUT = 6;
			fprintf ("spool",CLEAR_FILE,lf);*/
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
	}
	return 0;	
}



/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function PopulateModelMatrix (ModelMatrixName&, EFV)
{
	ModelMatrixName = {16,16}; 

	hshift = 0;
	
	modelDefString = "";
	modelDefString*16384;
	globalDefString = "";
	globalDefString*128;
	
	catCounterAL = {};
	
	/* build a branch length constraint */
	
	global meanOneConstraint = 1;
	meanOneConstraintString = "";
	
	for (aa1 = 0; aa1 < 15; aa1 = aa1+1)
	{
		for (aa2 = aa1+1; aa2 < 16; aa2 = aa2+1)
		{
			bt = aaRateMultipliers[aa1][aa2];
			if (catCounterAL[bt] == 0)
			{
				globalDefString*("\nglobal R_"+bt+"=1;");
				meanOneConstraintString = meanOneConstraintString + "+R_"+bt;
			}
			catCounterAL[bt] = catCounterAL[bt]+2*EFV[aa2]*EFV[aa1];
			bt = "R_"+bt+"*meanOneConstraint*";
			modelDefString*("ModelMatrixName["+aa1+"]["+aa2+"] := "+bt+"t;\nModelMatrixName["+aa2+"]["+aa1+"] := "+bt+"t;\n");
		}	
    }		
	modelDefString*0;
	/*bt = aaRateMultipliers[0][3];
	globalDefString*("\nR_"+bt+":=1;\n");*/
	globalDefString*0;
	ExecuteCommands (globalDefString);
	ExecuteCommands (modelDefString);
	ExecuteCommands ("global meanOneConstraint:=1/("+meanOneConstraintString[1][Abs(meanOneConstraintString)-1]+");");
	
	modelDefString = "";
	modelDefString*128;
	modelDefString*("global branchLengthScaler := "+catCounterAL[0]+"*R_0");
	
	for (aa1 = 1; aa1 < Abs(catCounterAL); aa1 = aa1+1)
	{
		modelDefString * ("+" + catCounterAL[aa1]+"*R_" + aa1);
	}
	modelDefString*0;
	ExecuteCommands (modelDefString+";");
	return 0;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function SpawnRandomString (clsCnt)
{
	rModel = {120,1};
	for (h=0; h<120; h=h+1)
	{
		rModel[h] = Random(0,clsCnt)$1;
	}
	
	return MakeStringCanonical(rModel,clsCnt);
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function MakeStringCanonical (randomModel, classCount)
{
	compressedString = {classCount,1};
	for (h=0; h<120; h=h+1)
	{
		v = randomModel[h];
		compressedString[v] = 1;
	}
	compressedString[0] = 0;
	for (h=1; h<classCount; h=h+1)
	{
		compressedString[h] = compressedString[h]+compressedString[h-1];
	}
	for (h=0; h<120; h=h+1)
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
		for (h=0; h<120; h=h+1)
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
		for (h=0; h<120; h=h+1)
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
	matrixSpec = {16,16};
	cc = 0;
	for (h=0; h<16; h=h+1)
	{
		for (v=h+1; v<16; v=v+1)
		{
			matrixSpec[h][v] = stringSpec[cc];
			matrixSpec[v][h] = stringSpec[cc];
			cc = cc+1;
		}
	}
	return matrixSpec;
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
	aaRateMultipliers 				= StringToMatrix(currentPopulation[bestInd]);
	AAModelMatrix 					= 0;
	MULTIPLY_BY_FREQS 				= PopulateModelMatrix ("AAModelMatrix", vectorOfFrequencies);
	Model AAModel	     			= (AAModelMatrix,vectorOfFrequencies,1);
	Tree givenTree				    = treeString;
	LikelihoodFunction 				lf = (filteredData,givenTree);
	
	fprintf 						(progressFile,"\nUpdating branch lengths...");
	Optimize (res,lf);
	global Null_Scaler 				= 	branchLengthScaler;
	
	myAIC 	  = 2*(res[1][0]-(res[1][1])*sampleCount/(sampleCount-res[1][1]-1));
	fprintf							(stdout, "\nUpdated BLs\nAICs: ", -myAIC, "\t",  -sortedScores[populationSize-1][0], "\n", lf, "\n");
	sortedScores[populationSize-1][0] 		 = myAIC;
	MasterList[sampleString] 				 = myAIC;
	
	return 0;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/



SetDialogPrompt 	  		  ("Locate a dinucleotide file:");
DataSet 		ds  		 = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter	filteredData = CreateFilter(ds,2);

sampleCount 				 = filteredData.sites;

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

HarvestFrequencies (vectorOfFrequencies,filteredData,2,2,0);

SetDialogPrompt ("Save the best model to:");
fprintf (PROMPT_FOR_FILE,CLEAR_FILE);
modelFile = LAST_FILE_PATH;

fpath = modelFile+".details";
fprintf (fpath,CLEAR_FILE);

aaRateMultipliers = {16,16};

aa_plain = 0;
MULTIPLY_BY_FREQS 				= PopulateModelMatrix ("aa_plain", vectorOfFrequencies);
Model AAPlainModel		 		= (aa_plain,vectorOfFrequencies,1);

incFileName = HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"queryTree.bf";
ExecuteCommands  ("#include \""+incFileName+"\";");
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

ChoiceList (seedChoice, "Seeds",1,SKIP_NONE,
			"Random","The entire starting population is seeded at random.",
			"Load", "Load some predefined partitions",
			"Completely Random","The entire starting population is seeded at random, skipping even hardwired presets."
);

if (seedChoice<0)
{
	return 0;
}

if (seedChoice == 2)
{
	predef = {};
}

ChoiceList (reprFlag, "Mating Representation",1,SKIP_NONE,
			"General","Use n-ary representation (n=number of rate classes).",
			"Binary", "Use binary representation");

if (reprFlag<0)
{
	return 0;
}

if (seedChoice == 1)
{
	SetDialogPrompt ("File with partition specifications:");
	END_OF_FILE = 0;
	fscanf (PROMPT_FOR_FILE,"NMatrix",seeds);
	predef[Abs(predef)] = seeds;
	while (!END_OF_FILE)
	{
		fscanf (LAST_FILE_PATH,"NMatrix",seeds);
		if (Abs(seeds))
		{
			predef[Abs(predef)] = seeds;
		}
	}
}

currentPopulation  = {};
Optimize (res,lf);
baseParams = res[1][1];
sortedScores	   = {populationSize,2};
sortedScores[0][0] = 2*(res[1][0]-baseParams*sampleCount/(sampleCount-baseParams-1));
sortedScores[0][1] = 0;
currentPopulation [0] = {stateVectorDimension,1};

global Null_Scaler = branchLengthScaler;


/* INCLUDE THE CHC CORE */

if (reprFlag)
{
	ibfPath = HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "GA_CHC_Binary.ibf";
}
else
{
	ibfPath = HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "GA_CHC.ibf";
}
ExecuteCommands ("#include \"" + ibfPath + "\";");

/* try a consensus matrix */

ExportAMatrix (modelFile,StringToMatrix(currentPopulation[populationSize-1]),1);
fpath = modelFile+".bestAIC";

fprintf (fpath, CLEAR_FILE, -sortedScores[populationSize-1][0]);
