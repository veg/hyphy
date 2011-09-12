function PadString (padLength,padChar)
{
	for (padCounter=0;padCounter<padLength;padCounter=padCounter+1)
	{
		fprintf (stdout,padChar);
	}
	return padLength;
}

fprintf(stdout,"\n ---- RUNNING LOCAL MOLECULAR CLOCKS ANALYSIS ---- \n");

ChoiceList (dataType,"Data type",1,SKIP_NONE,"Nucleotide/Protein","Nucleotide or amino-acid (protein).",
				     "Codon","Codon (several available genetic codes).");

if (dataType<0) 
{
	return;
}
if (dataType)
{
	NICETY_LEVEL = 3;
	#include "TemplateModels/chooseGeneticCode.def";
}

SetDialogPrompt ("Choose the data file:");

DataSet ds = ReadDataFile (PROMPT_FOR_FILE);

fprintf (stdout,"The following data was read:\n",ds,"\n");

if (dataType)
{
	DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);
}
else
{
	DataSetFilter filteredData = CreateFilter (ds,1);
}

SelectTemplateModel(filteredData);

#include "queryTree.bf";

global RelRatio;

RelRatio = 1.0;

relationString = ":=RelRatio*";

parameter2Constrain = 0;

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

/*internalNodes = BranchCount(givenTree);

choiceMatrix = {internalNodes,2};

for (bc=0; bc<internalNodes; bc=bc+1)
{
	choiceMatrix[bc][0] = BranchName(givenTree,bc);
	choiceMatrix[bc][1] = "Internal Branch Rooting " + givenTree[bc];
}

ChoiceList  (bOption,"Choose root of subtree",1,NO_SKIP,choiceMatrix);

if (bOption < 0)
{
	return -1;
}*/


LikelihoodFunction lf = (filteredData,givenTree);

Optimize (res,lf);

fprintf (stdout, "\n\nRESULTS WITHOUT THE CLOCK:\n",lf);

fullModelLik = res[1][0];

fullVars = res[1][1];

/* now specify the constraints */

fprintf (stdout, "\n\n\nRESULTS WITH LOCAL CLOCKS\n\n(*)   - Clock rejected at significance level 0.05\n(**)  - Clock rejected at significance level 0.01\n(***) - Clock rejected at significance level 0.001\n\n");

Tree clockTree = treeString;

intBranchCount  = BranchCount (clockTree);

nodeLength = 10;

for (nodeCounter = 0; nodeCounter < intBranchCount; nodeCounter = nodeCounter+1)
{
	anIntBranch   = BranchName (givenTree,nodeCounter);
	anIntBranch   = Abs (anIntBranch)+2;
	if (anIntBranch>nodeLength)
	{
		nodeLength = anIntBranch;
	}
}

separator = "+----------";

for (nodeCounter = 10; nodeCounter < nodeLength; nodeCounter = nodeCounter+1)
{
	separator = separator + "-";
}

separator = separator + "+--------------+-------------+------------+\n";

fprintf (stdout, separator, "| Clock At ");
dummy = PadString (nodeLength-10," ");
fprintf (stdout, "| LR Statistic | Constraints |  P-Value   |\n",separator);

if (MPI_NODE_COUNT>1)
{
	MPINodeState = {MPI_NODE_COUNT-1,1};
	MPINodeRoot  = {MPI_NODE_COUNT-1,1};
	MPINodeRoot[0]  = "";
	OPTIMIZE_SUMMATION_ORDER = 0;
}

for (nodeCounter = 0; nodeCounter < intBranchCount; nodeCounter = nodeCounter+1)
{
	anIntBranch   = BranchName (givenTree,nodeCounter);
	
	ExecuteCommands ("MolecularClock (clockTree."+anIntBranch+","+parameter2ConstrainString+");");

	LikelihoodFunction lfConstrained = (filteredData, clockTree);

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
			MPISend (mpiNode+1,lfConstrained);
			MPINodeState[mpiNode] = 1;
			MPINodeRoot[mpiNode] = anIntBranch;
		}
	}
	else
	{
		Optimize (lfConstrained_MLES,lfConstrained);
		dummy = ReceiveJobs (0);
	}
	
	ClearConstraints (clockTree);
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
	OPTIMIZE_SUMMATION_ORDER = 1;
}

function ReceiveJobs (sendOrNot)
{
	if (MPI_NODE_COUNT>1)
	{
		MPIReceive (-1, fromNode, result_String);
		anIntBranch2 = MPINodeRoot[fromNode-1];
		if (sendOrNot)
		{
			MPISend (fromNode,lfConstrained);
			MPINodeRoot[fromNode-1] = anIntBranch;
		}
		else
		{
			MPINodeState[fromNode-1] = 0;
			MPINodeRoot[fromNode-1]  = "";
		}
		
		anIntBranch = anIntBranch2;
		
		ExecuteCommands (result_String);
	}
	
	lnLikDiff = 2(fullModelLik-lfConstrained_MLES[1][0]);

	degFDiff = fullVars - lfConstrained_MLES[1][1];
	
	fprintf (stdout, "|", anIntBranch);

	dummy = PadString (nodeLength-Abs(anIntBranch)," ");
	
	pValue = 1-CChi2(lnLikDiff,degFDiff);

	fprintf (stdout, "| ", Format (lnLikDiff,12,6), " | ", Format (degFDiff,11,0), " | ", Format (pValue,10,8), " |");

	if (pValue<0.001)
	{
		fprintf (stdout, " (***)");
	}
	else
	{
		if (pValue<0.01)
		{
			fprintf (stdout, " (**)");
		}
		else
		{
			if (pValue<0.05)
			{
				fprintf (stdout, " (*)");
			}
		}
	}

	fprintf (stdout, "\n",separator);
			
	return fromNode-1;
}
