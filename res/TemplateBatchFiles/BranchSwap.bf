A_DISTANCE_METHOD 	   = 1;
SHORT_MPI_RETURN 	   = 1;

/*---------------------------------------------------------------*/

function cutTheStringUp		  (delimiterChar)
{
	while (treeStringIndex<treeStringLength)
	{
		aTreeChar = aTreeString[treeStringIndex];
		if ((aTreeChar == delimiterChar)&&(parenthesesDepth==0))
		{
			break;
		}
		else
		{
			if (aTreeChar == ")")
			{
				parenthesesDepth = parenthesesDepth-1;
			}
			else
			{
				if (aTreeChar == "(")
				{
					parenthesesDepth = parenthesesDepth+1;
				}
			}
		}
		treeStringIndex = treeStringIndex+1;
	}
	return 0;
}

/*---------------------------------------------------------------*/

function getTheClustersToSwap (aTreeString)
{
	clusterOne   = "";
	clusterTwo   = "";
	clusterThree = "";
	clusterFour  = "";
	
	treeStringLength = Abs (aTreeString);
	aTreeString		 = aTreeString [2][treeStringLength-3];
	parenthesesDepth = 0;
	treeStringIndex  = 0;
	treeStringLength = treeStringLength - 4;
	startCuttingAt   = 0;
	cutTheStringUp (",");
	clusterOne = aTreeString[startCuttingAt][treeStringIndex-1];
	if (Abs(clusterOne)==0)
	{
		return 1;
	}
	treeStringIndex = treeStringIndex+1;
	startCuttingAt   = treeStringIndex;
	cutTheStringUp (")");
	clusterTwo = aTreeString[startCuttingAt][treeStringIndex-1];
	if (Abs(clusterTwo)==0)
	{
		return 2;
	}
	treeStringIndex = treeStringIndex + 1;
	cutTheStringUp (",");
	treeStringIndex = treeStringIndex+2;
	startCuttingAt   = treeStringIndex;
	cutTheStringUp (",");
	clusterThree = aTreeString[startCuttingAt][treeStringIndex-1];
	if (Abs(clusterThree)==0)
	{
		return 3;
	}
	startCuttingAt = treeStringIndex+1;
	cutTheStringUp (")");
	clusterFour = aTreeString[startCuttingAt][treeStringIndex-1];
	return 4;
}

/*---------------------------------------------------------------*/

function ReceiveJobs (sendOrNot)
{
	if (MPI_NODE_COUNT>1)
	{
		MPIReceive (-1, fromNode, result_String);
		branchSwap_2 = MPINodeTree[fromNode-1][0];
		anIntBranch2 = MPINodeTree[fromNode-1][1];
		if (sendOrNot)
		{
			MPISend (fromNode,theLF);
			MPINodeTree[fromNode-1][0] = branchSwap;
			MPINodeTree[fromNode-1][1] = anIntBranch;
		}
		else
		{
			MPINodeState[fromNode-1]    = 0;
			MPINodeTree[fromNode-1][0]  = "";
			MPINodeTree[fromNode-1][1]  = "";
		}
		
		anIntBranch = anIntBranch2;
		branchSwap	= branchSwap_2;
		ExecuteCommands (result_String);
		
		res = theLF_MLES;
	}
	
	fprintf (stdout, totalRearrangementsTried, "). Swap at:", anIntBranch, ".Log-L = ", Format(res[1][0],10,5)," (", Format(res[1][0]-originalValueToBeat,10,5), ")\n");		
	
	dummy = CheckForImprovement(0);
	
	return fromNode-1;
}

/*---------------------------------------------------------------*/

function CheckForImprovement (dummy)
{
	if (res[1][0]>valueToBeat+0.5*OPTIMIZATION_PRECISION)
	{
		if (globalOption == 2)
		{
			ClearConstraints (testTree);
			ExecuteCommands (setString);
			Optimize (res, theLF);
		}
		totalSupportWeight = totalSupportWeight*Exp(valueToBeat-res[1][0])+1;
				
		valueToBeat     = res[1][0];
		currentBestTree = branchSwap;
		doSwapping 		= 1; 
		savedMLEs		= res;
		
		if (!swapType)
		{
			ibc = intBranchCount;
		}
	}
	else
	{
		totalSupportWeight = totalSupportWeight + Exp(res[1][0]-valueToBeat);
	}
	return 0;
}

/*---------------------------------------------------------------*/

function StartSwappingJob (dummy)
{
	USE_LAST_RESULTS = 1;
	
	branchSwap = RerootTree (branchSwap,0);
	
	Tree 	   testTree = branchSwap;		
	
	/* map parameter values from one tree to the next */
	
	GetInformation (params1, "^givenTree\.[^\.]+\..+$");
	GetInformation (params2, "^testTree\.[^\.]+\..+$");
	
	initString = "";
	initString * 8192;
	existingData = {};
	
	for (mpiNode=Columns(params1)-1; mpiNode >= 0; mpiNode = mpiNode-1)
	{
		eVar = params1[mpiNode];
		eVar = eVar[10][Abs(eVar)-1];
		existingData [eVar] = 1;
	}
	
	if (globalOption == 2)
	{
		for (mpiNode=Columns(params2)-1; mpiNode >= 0; mpiNode = mpiNode-1)
		{
			eVar = params2[mpiNode];
			eVar = eVar[9][Abs(eVar)-1];
			eIdx = existingData[eVar];
			if (eIdx)
			{
				initString * ("testTree."+eVar+":=givenTree."+eVar+"__;\n");
			}
		}
	}
	else
	{
		for (mpiNode=Columns(params2)-1; mpiNode >= 0; mpiNode = mpiNode-1)
		{
			eVar = params2[mpiNode];
			eVar = eVar[9][Abs(eVar)-1];
			eIdx = existingData[eVar];
			if (eIdx)
			{
				initString * ("testTree."+eVar+"=givenTree."+eVar+";\n");
			}
		}
		
	}
	
	
	initString * 0;
	if (Abs(initString))
	{
		ExecuteCommands (initString);
	}
	
	/* done with mapping */
	
	LikelihoodFunction theLF = (filteredData, testTree);
	if (MPI_NODE_COUNT>1)
	{
		for (mpiNode = 0; mpiNode < MPI_NODE_COUNT-1; mpiNode = mpiNode+1)
		{
			if (MPINodeState[mpiNode]==0)
			{
				break;	
			}
		}
		
		if (mpiNode==MPI_NODE_COUNT-1)
		/* all nodes busy */
		{
			mpiNode = ReceiveJobs (1);
		}
		else
		{
			MPISend (mpiNode+1,theLF);
			MPINodeState[mpiNode] = 1;
			MPINodeTree[mpiNode][0] = branchSwap;
			MPINodeTree[mpiNode][1] = anIntBranch;
		}
	}
	else
	{
		Optimize (res, theLF);
		mpiNode = ReceiveJobs (0);
	}
	
	totalRearrangementsTried = totalRearrangementsTried + 1;
	return 0;
}

/*---------------------------------------------------------------*/

ChoiceList (swapType,"Swapping Strategy",1,SKIP_NONE,"Greedy","As soon as branch swapping yields a better tree, the swapping procedure is restarted on that tree.",
				     "Patient","After ALL branch swaps have been attempted on a given tree, the procedure is restarted on the tree with the best likelihood improvement over the original (if there is one, of course).");

if (swapType<0) 
{
	return;
}


ChoiceList (dataType,"Data type",1,SKIP_NONE,"Nucleotide/Protein","Nucleotide or amino-acid (protein).",
				     "Codon","Codon (several available genetic codes).");

if (dataType<0) 
{                   
	return;
}

#include "heuristicMethodNPBootstrap.bf";

MESSAGE_LOGGING = 0;

if (dataType)
{
	NICETY_LEVEL = 3;
	SetDialogPrompt ("Please choose a codon data file:");
	#include "TemplateModels/chooseGeneticCode.def";
}
else
{
	SetDialogPrompt ("Please choose a nucleotide or amino-acid data file:");
}

totalRearrangementsTried = 0;
totalSupportWeight		 = 1;

DataSet ds = ReadDataFile (PROMPT_FOR_FILE);

if (dataType)
{
	DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);
}
else
{
	DataSetFilter filteredData = CreateFilter (ds,1,"","");
}

SelectTemplateModel(filteredData);

_DO_TREE_REBALANCE_ = 1;

_KEEP_I_LABELS_ = 1;
#include "queryTree.bf";
_KEEP_I_LABELS_ = 0;

ChoiceList (globalOption,"Global parameters",1,SKIP_NONE,"Re-estimate","Re-estimate all global parameters for every rearrangement (slower).",
				     "Do once","Use the estimates obtained from the original tree for all rearrangements (faster, but less precise).",
				     "Quick and Dirty","Only re-estimate the branch involved in the swap, fix other parameter estimates");
				     
if (globalOption < 0)
{
	return 0;
}


SetDialogPrompt ("If a better tree is found, save it to:");
fprintf (PROMPT_FOR_FILE,CLEAR_FILE);
treeFileSave = LAST_FILE_PATH;


LikelihoodFunction theLF = (filteredData, givenTree);
Optimize (res, theLF);

if (globalOption > 0)
{
	setString 	= "";
	resetString	= "";
	if (res[1][2]>0)
	{
		setString * 128;
		resetString *128;
	
		for (k=0; k<res[1][2]; k=k+1)
		{
			GetString (eVar, theLF, k);
			setString * (eVar + ":="+eVar+"__;\n");
			resetString * (eVar + "="+eVar+";\n");
		}
	
		setString * 0;
		resetString *0;
		
		ExecuteCommands (setString);
	}
}

doSwapping 			= 1;
valueToBeat  		= res [1][0];
originalValueToBeat = valueToBeat;

svl = VERBOSITY_LEVEL;
VERBOSITY_LEVEL = -1;

fprintf	 (stdout, "\n\n****** ORIGINAL TREE AND LIKELIHOOD ******\n\n", 
		    	  theLF,"\n\n****** RUNNING BRANCH SWAPPING ANALYSIS ******\n\n");
		    	  
phaseCounter 	= 1;
currentBestTree = treeString;

Tree originalTree = treeString;

if (MPI_NODE_COUNT>1)
{
	MPINodeState = {MPI_NODE_COUNT-1,1};
	MPINodeTree  = {MPI_NODE_COUNT-1,2};
	MPINodeTree[0]  = "";
}

while (doSwapping)
{
	doSwapping    = 0;
	intBranchCount  = BranchCount (givenTree);

	fprintf	 (stdout, "\n> PHASE ",phaseCounter,"\n\n"); 

	for (ibc = 0; ibc < intBranchCount; ibc=ibc+1)
	{
		anIntBranch   = BranchName (givenTree,ibc);
		aRerootedTree = RerootTree (givenTree, anIntBranch);
			
		dummy = getTheClustersToSwap (aRerootedTree);
		branchSwap = "(("+clusterOne+","+clusterThree+"),"+clusterTwo+","+clusterFour+")";
		dummy = StartSwappingJob(0);
		if (ibc == intBranchCount)
		{
			break;
		}
		branchSwap = "(("+clusterOne+","+clusterFour+"),"+clusterTwo+","+clusterThree+")";
		dummy = StartSwappingJob(0);
	}		
	if (MPI_NODE_COUNT>1)
	{
		while (1)
		{
			for (nodeCounter = 0; nodeCounter < MPI_NODE_COUNT-1; nodeCounter = nodeCounter+1)
			{
				if (MPINodeState[nodeCounter]==1)
				{
					fromNode = ReceiveJobs (0);
					break;	
				}
			}
			if (nodeCounter == MPI_NODE_COUNT-1)
			{
				break;
			}
		}	
	}	
	if (doSwapping)
	{
		Tree givenTree = currentBestTree;
		fprintf (stdout, "\n\n**** AFTER PHASE ",phaseCounter,"\n\n>Best Tree:\n",currentBestTree,"\n>Log-likelihood = ",valueToBeat," (Improvement of ", valueToBeat-originalValueToBeat,")\n\n"); 
		phaseCounter = phaseCounter+1;
	}	
}

if (globalOption > 0)
{
	if (Abs(resetString)>0)
	{
		ExecuteCommands (resetString);
	}
}

fprintf (stdout, "\n\n***** ", totalRearrangementsTried, " TREE REARRANGEMENTS WERE EXPOLORED *****\n\n",
				 "\n\n***** RELATIVE TREE SUPPORT (WITH UN-INFORMATIVE PRIOR): ", 1/totalSupportWeight, " ******\n\n");

if (phaseCounter>1)
{
	fprintf (stdout, "\n\n***** BRANCH SWAPPING FOUND A BETTER TREE! *****\n\n");
	Tree			   newBestTree = currentBestTree;
	NO_INTERNAL_LABELS = 1;
	fprintf (treeFileSave,newBestTree);	
	fprintf (stdout, newBestTree,"\n");
	NO_INTERNAL_LABELS = 0;
	fprintf (stdout, "\n\nA likelihood improvement of ",valueToBeat-originalValueToBeat,"\n\n");
}
else
{
	fprintf (stdout, "\n\n***** BRANCH SWAPPING FAILED TO YIELD A BETTER TREE *****\n\n");
	NO_INTERNAL_LABELS = 1;
	fprintf (treeFileSave,givenTree);	
	NO_INTERNAL_LABELS = 0;
}


VERBOSITY_LEVEL = svl;
