RequireVersion ("0.9920070101");
VERBOSITY_LEVEL = -1;

ExecuteAFile(HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "TreeTools.ibf");
ExecuteAFile(HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "Utility" + DIRECTORY_SEPARATOR + "GrabBag.bf");
#include "GADatedClock.ibf";
#include "_tipDater.ibf";

/*________________________________________________________________________________________________*/

produceOffspring		= MPI_NODE_COUNT-1;
if (produceOffspring <= 0)
{
	produceOffspring = 15;
}

populationSize  		= 2*produceOffspring;
incestDistance  		= 0;
generationCount		  	= 5000;
maxSampleTries			= populationSize*10;
mutationThreshhold		= 0.0001;
mutationProb			= 0.15;
mutationProbDecrease	= 0.95;
annealingPhase			= 100;
totalSampleCounter		= 0;
localMutationRate		= 0.05;
localMutationInterval	= 15;

stoppingCriterion		= 50;
sampleCount				= 0;
familyControlSize		= produceOffspring$6;

verboseFlag				= 0;
stateVectorDimension    = 120;
rateClassesCount		= 2;
startWithRateClasses	= 2;
autoStepFlag			= 1;

/* ________________________________________________________________________________________________*/

MasterList				= {};
REPLACE_TREE_STRUCTURE  = 1;
SHORT_MPI_RETURN		= 1;

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function StringToMatrix (zz)
{
	return zz;
}


/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function ExportAMatrix (fileName, rateMatrix,dummy, dummy2)
{
	if (dummy == 1)
	{
		return 0;
	}
	
	fprintf (fileName, CLEAR_FILE,ConvertToPartString(rateMatrix));
	return 0;
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

/*-----------------------------------------------------------*/

function generateBLVector (treeNameID)
{
	ExecuteCommands ("treePostOrderAVL = "+treeNameID+"^0;");
	_blVector = {};
		
	for (nodeIndex = 1; nodeIndex < nodeCount; nodeIndex = nodeIndex+1)
	{
		nodeInfo 	= treePostOrderAVL[nodeIndex];
		nodeNameS	= nodeInfo["Name"];
		nodeParent = nodeInfo["Parent"];
		if (nodeParent)
		{
			parentNodeID = nodeParent;
			nodeParent = treePostOrderAVL[nodeParent];
		}
		if (Abs(nodeParent))
		{
			pName    = nodeParent["Name"];
			if (Abs(nodeInfo["Children"]))
			{
				ExecuteCommands ("_thisBL = " + treeNameID+"_"+nodeNameS+"_BL;tipDateAVL[\""+nodeNameS+"\"]="+ treeNameID+"_"+nodeNameS+"_T;");
			}
			else
			{
				ExecuteCommands ("_thisBL = "+ tipDateAVL[nodeNameS] +" -" + treeNameID+"_"+pName+"_T;");			
			}
			_blVector [nodeNameS] = _thisBL;
		}
	}
	
	
	return _blVector;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function ReceiveJobs (sendOrNot, ji)
{
	if (MPI_NODE_COUNT>1)
	{
		MPIReceive (-1, fromNode, result_String);
		mji = MPINodeState[fromNode-1][1];
		
		if (sendOrNot)
		{
			MPISend (fromNode,mpiMessageToSend);
			MPINodeState[fromNode-1][1] = ji;			
		}
		else
		{
			MPINodeState[fromNode-1][0] = 0;
			MPINodeState[fromNode-1][1] = -1;		
		}
		ExecuteCommands (result_String);
		myAIC 			= 2*(_hyphyAssociativeArray["LogL"]-_hyphyAssociativeArray["DF"]-baseParams);
		tsv				= _hyphyAssociativeArray["tsv"];
		rsv				= _hyphyAssociativeArray["rsv"];
		mtd				= _hyphyAssociativeArray["tree"];
		ji 				= mji;
	}
	else
	{
		myAIC 			= 2*(res[1][0]-res[1][1]-baseParams);
		tsv 			= generateTimeStampVector("clockTree");
		if (ji>=0)
		{
			if (resultProcessingContext==0)
			{
				rsv = generateRatesVector ("clockTree",Max(currentPopulation[ji],0));
			}
			else
			{
				rsv = generateRatesVector ("clockTree",Max(children[ji-populationSize],0));
			}
			mtd = Format (clockTree,1,1);
		}
	}

	if (resultProcessingContext==0)
	{
		sortedScores[ji][0] = myAIC;
		if (ji>=0)
		{
			jobPrint = ConvertToPartString (currentPopulation[ji]);
		}
		if (verboseFlag)
		{
			fprintf (stdout, "Individual ",ji," AIC = ",-myAIC," ");
		}
	}
	else
	{
		intermediateProbs[ji][0] = myAIC;	
		if (ji>=0)
		{
			jobPrint = ConvertToPartString (children[ji-populationSize]);
		}
		if (verboseFlag)
		{
			fprintf (stdout, "Offspring ",ji," AIC = ",-myAIC," ");
		}
	}
	
	if (Abs(jobPrint))
	{
		MasterList [jobPrint] = myAIC;
		if ( (-myAIC) < bestModelIC)
		{
			bestModelIC = -myAIC;
			bestTreeString = mtd;
		}
		
		if (verboseFlag)
		{
			fprintf (stdout, " ", jobPrint, "\n");
		}
		fprintf (detailedFile, "\n{{", jobPrint, "}}\n", -myAIC,"\n", tsv, "\n",rsv);
	}
	
	return fromNode-1;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function ConvertToPartString (modelSpec)
{
	outString = "";
	outString * 128;
	outString * (""+modelSpec[0]);
	
	for (zzz = 1; zzz < stateVectorDimension; zzz = zzz + 1)
	{
		outString * (","+modelSpec[zzz]);
	}
	
	outString * 0;
	return outString;
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

function RunASample (dummy,jobIndex)
{	
	myAIC = MasterList[cString];
	if (myAIC<0)
	{		
		if (resultProcessingContext==0)
		{
			sortedScores[jobIndex][0] = myAIC;
			if (verboseFlag)
			{
				fprintf (stdout, "Individual ",jobIndex," AIC = ",-myAIC, "\n");
			}
		}
		else
		{
			intermediateProbs[jobIndex][0] = myAIC;	
			if (verboseFlag)
			{
				fprintf (stdout, "Offspring ",jobIndex," AIC = ",-myAIC,"\n");
			}
		}	
	}
	else
	{
		ExecuteCommands ("treePostOrderAVL = "+treeNameID+"^0;");
		nodeCount 	=  Abs(treePostOrderAVL);

		if ((MPI_NODE_COUNT>1) && (jobIndex>=0))
		{
			mpiMessageToSend = "\ncString="+cString+ 
							   ";\nExecuteCommands 	(generateDatedTipConstraints  (\"clockTree\",parameter2ConstrainString,tipDateAVL,cString,singleModelValues));"+
							   "LikelihoodFunction lf			 = (filteredData, clockTree);Optimize (res,lf); return makeReturnAVL(0);";

			mpiMessageToSend = _mpiPrefixString + mpiMessageToSend;
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
				MPISend (mpiNode+1,mpiMessageToSend);
				/*fprintf ("last.sent", CLEAR_FILE, MPI_LAST_SENT_MSG);*/
				MPINodeState[mpiNode][0] = 1;
				MPINodeState[mpiNode][1] = jobIndex;
			}
		}
		else
		{
			ClearConstraints 				   (clockTree);
			ExecuteCommands 				   (generateDatedTipConstraints  ("clockTree",parameter2ConstrainString,tipDateAVL,cString,singleModelValues));
			LikelihoodFunction lf			 = (filteredData, clockTree);
			
			Optimize (res,lf);
			
			if (jobIndex>=0)
			{
				mpiNode = ReceiveJobs (1, jobIndex);
			}
			else
			{
				myAIC = 2*(res[1][0]-res[1][1]-baseParams);
			}
		}
	}
	return 0;	
}




/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function SpawnRandomString (clsCnt)
{
	rModel = {stateVectorDimension,1};
	rModel = rModel ["Random(0,rateClassesCount)$1"];
	return 	IsChildViable(MakeStringCanonical(rModel,rateClassesCount));
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function IsChildViable (putativeChild)
{
	myAIC 		 = MasterList[ConvertToPartString(putativeChild)];
	testChild 	 = putativeChild;
	mutPassCount = 1;
	
	while (myAIC < (-0.1) && mutPassCount < 25)
	{
		if (verboseFlag)
		{
			fprintf (stdout,"Adjusting the child to avoid a duplicate. Pass ", mutPassCount, "\n");
		}
		
		mutPassCount = mutPassCount + 1;
		
		sampleString = Min(Random(0,stateVectorDimension)$1,stateVectorDimension-1);
		myAIC 		 = testChild[sampleString];
		
		newValue = Random (0,rateClassesCount-0.0000001)$1;
		
		while (newValue == myAIC)
		{
			newValue = Random (0,rateClassesCount-0.0000001)$1;
		}
		
		testChild [sampleString] = newValue;
		myAIC 					 = MasterList[ConvertToPartString(testChild)];
	}
	return testChild;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function UpdateBL (dummy)
{
	return 0;
}


/*---------------------------------------------------------------------------------------------------------------------------------------------*/


/*
ChoiceList (dataType,"Data type",1,SKIP_NONE,"Nucleotide","Nucleotide data.",
				     "Codon","Codon (several available genetic codes).");
*/

dataType = 0;
if (dataType<0) 
{
	return;
}

if (dataType==0)
{
	SetDialogPrompt 	  		  ("Locate a nucleotide data  file:");
}
else
{
	SetDialogPrompt 	  		  ("Locate a codon data file:");
}

DataSet 		ds  		 = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter	filteredData = CreateFilter(ds,1);

_mpiPrefixString = "";
_mpiPrefixString * 256;

_mpiPrefixString * ("ExecuteAFile(HYPHY_LIB_DIRECTORY + \"TemplateBatchFiles\" + DIRECTORY_SEPARATOR + \"TreeTools.ibf\");\nExecuteAFile(HYPHY_LIB_DIRECTORY + \"TemplateBatchFiles\" + DIRECTORY_SEPARATOR + \"GADatedClock.ibf\");");
_mpiPrefixString * ("DataSet 	ds    = ReadDataFile (\""+LAST_FILE_PATH+"\");");
_mpiPrefixString * ("DataSetFilter	filteredData = CreateFilter(ds,1);ACCEPT_ROOTED_TREES=1;");

SelectTemplateModel 	(filteredData);
ExecuteAFile  			(HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"queryTree.bf");

ChoiceList			   (df, "Date Format", 1, SKIP_NONE,
							"TipDate","Taxon names include sampling dates, e.g. ((s1_85,s2_90),t_83)",
							"Branch Lengths", "Branch lengths are the sampling dates, e.g. ((s1:85,s2:90),t:83)");

if (df<0)
{
	return 0;
}			
else
{				
	if (df == 0)
	{							
		tipDateAVL 			   = getTipDatesFromNames1 ("givenTree");
	}
	else
	{
		tipDateAVL 			   = getTipDatesFromNames2 ("givenTree");	
	}
	_mpiPrefixString * ("tipDateAVL="+tipDateAVL+";");
}

fprintf (stdout, "What units are the dates measured in (e.g. months. This is only used for reporting the results.)?");
fscanf  (stdin,"String", dateUnit);


if (Abs(tipDateAVL) == 0)
{
	fprintf (stdout, "\nERROR: Input tree must contain date information in the chosen format for all sequences.\n",
					  "Branch lengths for internal nodes are ignored.\n");
	return 0;
}
else
{
	fprintf (stdout, "Read the following dates: \n");
	seqNames = Rows (tipDateAVL);
	for (sid = 0; sid < Columns (seqNames); sid = sid + 1)
	{
		fprintf (stdout, seqNames[sid], ":\t", retAVL[seqNames[sid]]/maxV,"\t",dateUnit,"\n");
	}
}

parameter2Constrain = 0;

ACCEPT_ROOTED_TREES = 1;
Tree clockTree = treeString;


if (Rows("LAST_MODEL_PARAMETER_LIST")>1)
{
	ChoiceList (parameter2Constrain, "Parameter(s) to constrain:",1,SKIP_NONE,LAST_MODEL_PARAMETER_LIST);

	if (parameter2Constrain<0)
	{
		return;
	}
	if (parameter2Constrain==0)
	{
		parameter2ConstrainString = "";
		for (parameter2Constrain=Rows("LAST_MODEL_PARAMETER_LIST")-1; parameter2Constrain; parameter2Constrain = parameter2Constrain-1)
		{
			GetString (funnyString,LAST_MODEL_PARAMETER_LIST,parameter2Constrain);
			parameter2ConstrainString = parameter2ConstrainString + funnyString + ",";
		}
		GetString (funnyString,LAST_MODEL_PARAMETER_LIST,0);
		parameter2ConstrainString = parameter2ConstrainString + funnyString;
	}
	else
	{
		GetString (parameter2ConstrainString,LAST_MODEL_PARAMETER_LIST,parameter2Constrain-1);
	}
}
else
{
	GetString (parameter2ConstrainString,LAST_MODEL_PARAMETER_LIST,0);
}

_mpiPrefixString * ("parameter2ConstrainString=\"" + parameter2ConstrainString + "\";maxV="+maxV+";minV="+minV+";\nUSE_DISTANCES=0;\nUSE_LAST_RESULTS = 1;\nMAXIMUM_ITERATIONS_PER_VARIABLE=1000;OPTIMIZATION_PRECISION=0.01;");

SetDialogPrompt ("Save the best model to:");
fprintf 		(PROMPT_FOR_FILE,CLEAR_FILE);
modelFile 		= LAST_FILE_PATH;

Tree givenTree 				= treeString;
stateVectorDimension 		= Columns(BranchLength (clockTree,-1));


USE_LAST_RESULTS 				= 1;
LikelihoodFunction lf			= (filteredData,givenTree);
Optimize						(res_free,lf);


baseParams = res_free[1][2];
	
if (baseParams>0)
{
	ConstraintString = "";
	ConstraintString*256;
	for (h=0; h<baseParams; h=h+1)
	{
		GetString (v,lf,h);
		ConstraintString * (v+":="+v+"__;\n");
	}
	ConstraintString*0;
	ExecuteCommands (ConstraintString);
}

Export (modelString, USE_LAST_MODEL);
_mpiPrefixString * modelString;
_mpiPrefixString * ("Tree clockTree = " + treeString + ";");


IC_BOUND = 2*(res_free[1][1]-res_free[1][0]);
fprintf (stdout, "\n______________FREE RATES______________\n",lf,"\nBEST AIC=",IC_BOUND,"\n");


Tree clockTree 		= treeString;
cString 			= {stateVectorDimension,1};

global				clockTree_scaler_0;

ExecuteCommands 	(generateDatedTipConstraints  ("clockTree",parameter2ConstrainString,tipDateAVL,cString,0));


fprintf (stdout, "Fitting the single rate model...\n");

LikelihoodFunction lf 		= (filteredData,clockTree);

/*
LIKELIHOOD_FUNCTION_OUTPUT							  = 7;
fprintf ("/Users/sergei/Desktop/dump.lf", CLEAR_FILE, lf);
*/

Optimize 					(res,lf);

rateClassesCount = 1;
bestRES			 = res;
bestSnap		 = generateLFSnapshot("clockTree", cString);

MAXIMUM_ITERATIONS_PER_VARIABLE = 10000;
OPTIMIZATION_PRECISION			= 0.001;

for (k=0; k<3; k=k+1)
{
	fprintf 					(stdout, "Pass ", rateClassesCount,"(",k,")", " log-L = ", res[1][0], "\n");
	ExecuteCommands 			(generateDatedTipConstraints  ("clockTree",parameter2ConstrainString,tipDateAVL,cString,0));
	Optimize 					(res,lf);
	if (bestRES[1][0] < res[1][0] - OPTIMIZATION_PRECISION)
	{
		bestRES 		 = res;
		k				 = -1;
		bestSnap		 = generateLFSnapshot("clockTree", cString);
	}
	rateClassesCount = rateClassesCount + 1;
} 
VERBOSITY_LEVEL = 0;
ExecuteCommands (bestSnap);
res = bestRES;

MAXIMUM_ITERATIONS_PER_VARIABLE = 1000;
OPTIMIZATION_PRECISION			= 0.01;
USE_LAST_RESULTS 				= 1;


currentPopulation  = {};
sortedScores	   = {populationSize,2};
sortedScores[0][0] = 2*(res[1][0]-res[1][1]);
currentBEST_IC	   = -sortedScores[0][0];
sortedScores[0][1] = 0;

currentPopulation [0] = cString;
fprintf (stdout, "\n______________SINGLE RATE______________\n",lf,"\nAIC=",currentBEST_IC,"\n");
maxRateClasses	      = stateVectorDimension$2;

fprintf (stdout, "Maximum of ", maxRateClass, " rate classes will be considered");

/* store initial guesses from the single rate model to feed into the multiple rate model  */

singleModelValues 		  = {};
singleModelValues["Rate"] = clockTree_scaler_0;

for (nodeIndex = 1; nodeIndex < nodeCount; nodeIndex = nodeIndex+1)
{
	nodeInfo 	= treePreOrderAVL[nodeIndex];
	if (Abs(nodeInfo["Children"]))
	{
		nodeNameS	= nodeInfo["Name"];
		ExecuteCommands ("singleModelValues[\""+nodeNameS+"\"]= clockTree_"+nodeNameS+"_BL;");
	}
}	

nodeNameS			= (treePreOrderAVL[1])["Name"];
ExecuteCommands 	  ("singleModelValues[\""+nodeNameS+"\"]= clockTree_"+nodeNameS+"_T;");
bestTreeString 		= Format (clockTree,1,1);

detailedFile 	   = modelFile + ".samples";
fprintf				(detailedFile, CLEAR_FILE, KEEP_OPEN, Format(givenTree,1,1), "\n", bestTreeString, "\n", dateUnit);

tc            = TipCount(givenTree);
tipDateMatrix = {tc,1};
for (k=0; k<tc; k=k+1)
{
	tipDateMatrix [k] = tipDateAVL[TipName(givenTree,k)]/maxV;
}

fprintf (detailedFile, "\n", tipDateMatrix);


_mpiPrefixString * ("singleModelValues=" + singleModelValues + ";");
_mpiPrefixString * 0;

bestModelIC = currentBEST_IC;

for (rateClassesCount = startWithRateClasses; rateClassesCount < maxRateClasses; rateClassesCount = rateClassesCount + 1)
{
	fprintf (stdout, "\n\nStarting GA with ", rateClassesCount, " rate classes\n");
	addOnLine = " with " + rateClassesCount + " rate classes.";
	
	resultProcessingContext = 0;
	
	if (rateClassesCount > startWithRateClasses)
	{
		children = {};
		for (individual=populationSize $ 2; individual<populationSize-1; individual=individual+1)
		{
			cString    = currentPopulation[individual];
			toggleRate = Random (0,rateClassesCount-1.00001)$1;
			for (bitP = 0; bitP < stateVectorDimension; bitP = bitP + 1)
			{
				if (cString [bitP] == toggleRate)
				{
					if (Random(0,1)<0.5)
					{
						cString [bitP] = rateClassesCount-1;
					}
				}
			}
			cString = IsChildViable(MakeStringCanonical(cString,rateClassesCount));
			sortedScores[individual][1] 		= individual;
			currentPopulation[individual] 		= cString;
			RunASample (compressedString[rateClassesCount-1],individual);
		}
	}	
	
	ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "GA_CHC.ibf");
	kf				= -sortedScores[populationSize-1][0];
	
	if (currentBEST_IC > kf)
	{
		currentBEST_IC = kf;
		/*if (IC_BOUND > currentBEST_IC)
		better than the free rates model 
		{
			break;
		}*/
	}
	else
	{
		break;
	}
}

ExportAMatrix (modelFile,currentPopulation[populationSize-1],0, 0);
fprintf 	  (detailedFile, CLOSE_FILE);
fscanf		  (detailedFile, "Raw", MasterList);
fprintf		  (detailedFile, CLEAR_FILE, bestTreeString, "\n", MasterList);
