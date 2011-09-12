SHORT_MPI_RETURN = 1;

/*------------------------------------------------------------------------*/

function ReportSiteST (siteI, siteM)
{
	fullSites[siteI][0] = doneSites[siteM][0];
	fullSites[siteI][1] = doneSites[siteM][1];
	fullSites[siteI][2] = doneSites[siteM][2];
	fullSites[siteI][3] = doneSites[siteM][0]/doneSites[siteM][2];
	fullSites[siteI][4] = doneSites[siteM][1]/doneSites[siteM][2];
	
	fullSites[siteI][5] = doneSites[siteM][3];
	fullSites[siteI][6] = doneSites[siteM][4];
	fullSites[siteI][7] = doneSites[siteM][3]/doneSites[siteM][4];
	fullSites[siteI][8] = doneSites[siteM][5];
	fullSites[siteI][9] = doneSites[siteM][6];


	fprintf (stdout, "Site ", Format(siteI+1,4,0)," dN Tree= ", Format(fullSites[siteI][0],6,3),
					 " dN Clade ", Format(fullSites[siteI][1],6,3),
					 " dS = ", Format(fullSites[siteI][2],6,3),
					 " Joint dN ",Format(fullSites[siteI][5],6,3),
					 " Joint dS ",Format(fullSites[siteI][6],6,3),
					 " LRT ", Format(fullSites[siteI][8],6,3),
					 " p-value = ",Format(fullSites[siteI][9],7,4));		
					 
	if (fullSites[siteI][9]<pValue)
	{
		fprintf (stdout, " *");
	}
	fprintf (stdout, "\n");
	return 0;
}

/*------------------------------------------------------------------------*/

function ReceiveJobsST (sendOrNot, nullAlt)
{
	MPIReceive (-1, fromNode, result_String);
	
	siteIndex = MPINodeState[fromNode-1][1];
	siteNA	  = MPINodeState[fromNode-1][2];
	
	if (sendOrNot)
	{
		MPISend (fromNode,siteLikelihood);
		MPINodeState[fromNode-1][1] = siteCount;			
		MPINodeState[fromNode-1][2] = nullAlt;			
	}
	else
	{
		MPINodeState[fromNode-1][0] = 0;
		MPINodeState[fromNode-1][1] = -1;		
	}
	
	siteMap = dupInfo[siteIndex];
	
	/*
	echoFile = "echo."+siteIndex+"."+siteNA;
	fprintf (echoFile,CLEAR_FILE,result_String);
	*/
	
	ExecuteCommands (result_String);
	
	if (siteNA)
	{
		doneSites[siteMap][0] = siteLikelihood_MLE_VALUES["nFactor_Tree"];
		doneSites[siteMap][1] = siteLikelihood_MLE_VALUES ["nFactor_Clade"];
		doneSites[siteMap][2] = siteLikelihood_MLE_VALUES ["sFactor"];
		doneSites[siteMap][5] = doneSites[siteMap][5]+2*siteLikelihood_MLES[1][0];
	}
	else
	{
		doneSites[siteMap][5] = doneSites[siteMap][5]-2*siteLikelihood_MLES[1][0];	
		doneSites[siteMap][4] = siteLikelihood_MLE_VALUES["sFactor"];
		doneSites[siteMap][3] = siteLikelihood_MLE_VALUES["nFactor_Tree"];
	}

	if (doneSites[siteMap][6] == 0)
	{
		doneSites[siteMap][6] = -1;
	}
	else
	{
		if (doneSites[siteMap][6] == (-1))
		{
			doneSites[siteMap][6] = 1-CChi2(doneSites[siteMap][5],1);						
			dummy = ReportSiteST (siteIndex, siteMap);
		}
	}
	
	return fromNode-1;
}

/*------------------------------------------------------------------------*/

stdr = _DO_TREE_REBALANCE_;
_DO_TREE_REBALANCE_ = 0;

#include "qndhelper1.ibf";

_DO_TREE_REBALANCE_ = stdr;


fprintf (stdout,"\n\n______________READ THE FOLLOWING DATA______________\n",ds,
				"\n\nPhase 1:Nucleotide Model (",ModelTitle,") Model Fit\n\n");


if (nrChoice == 0)
{
	LikelihoodFunction nucLF = (nucData,givenTree);
	Optimize (res,nucLF);
	stashLOF = LIKELIHOOD_FUNCTION_OUTPUT ;
	LIKELIHOOD_FUNCTION_OUTPUT  = 6;
	fprintf (NUC_FILE_PATH,nucLF);
	LIKELIHOOD_FUNCTION_OUTPUT  = stashLOF;
}

mxTreeSpec = {5,1};
mxTreeSpec [0] = "givenTree";
mxTreeSpec [1] = "8240";
mxTreeSpec [2] = "10,40,-10,-175,1";
mxTreeSpec [3] = "";
mxTreeSpec [4] = "";
OpenWindow (TREEWINDOW, mxTreeSpec);


#include "qndhelper3.ibf";

global			sFactor 	  = 1;
global			nFactor_Tree  = 1;
global			nFactor_Clade = 1;

doneSites    = {filteredData.unique_sites,7};
fullSites    = {filteredData.sites,10};
Tree		   siteTree = treeString;



/* display the tree */

ChoiceList  (stOption,"Branch Subset",1,NO_SKIP,"Subtree","Select a subtree rooted at an internal node",
												"Internal vs Tips", "Internal nodes vs leaves of the tree");
												
if (stOptions < 0)
{
	return 0;
}
												
ReplicateConstraint ("this1.?.synRate   :=sFactor*this2.?.t__",siteTree,givenTree);

if (stOption == 0)
{
	internalNodes = BranchCount(siteTree);

	choiceMatrix = {internalNodes,2};

	for (bc=0; bc<internalNodes; bc=bc+1)
	{
		choiceMatrix[bc][0] = BranchName(siteTree,bc);
		choiceMatrix[bc][1] = "Internal Branch Rooting " + siteTree[bc];
	}

	ChoiceList  (bOption,"Choose root of subtree",1,NO_SKIP,choiceMatrix);

	if (bOption < 0)
	{
		return -1;
	}

	cladeConstraintCommand = "ReplicateConstraint(\"this1.?.nonSynRate:=nFactor_Clade*this2.?.t__\",siteTree."+BranchName(siteTree,bOption)+",givenTree."+BranchName(givenTree,bOption)+");";																					      
}
else
{
	cladeConstraintCommand = "";
	cladeConstraintCommand * 8192;
	internalNodes = BranchCount(siteTree);
	for (bc=0; bc<internalNodes; bc=bc+1)
	{
		ibn = BranchName(siteTree,bc);
		cladeConstraintCommand * ("siteTree."+ibn+".nonSynRate:=nFactor_Clade*givenTree."+ibn+".t__;\n");
	}
	cladeConstraintCommand * 0;
}

pValue = 0;
while ((pValue<=0)||(pValue>=1))
{	
	fprintf (stdout, "\nSignificance level for Likelihood Ratio Tests (between 0 and 1)?");
	fscanf  (stdin,"Number", pValue);
}
GetDataInfo    (dupInfo, filteredData);
fprintf (stdout, "\n");



ExecuteCommands (cladeConstraintCommand);
ReplicateConstraint ("this1.?.nonSynRate:=nFactor_Tree*this2.?.t__",siteTree,givenTree);


labels = {{"dN Tree","dN Clade", "dS","dN/dS Tree", "dN/dS Clade",  "dN Tree = dN Clade", "Joint dS", "Joint dN/dS", "LRT","p-value"}};

if (MPI_NODE_COUNT<=1)
{
	for (siteCount = 0; siteCount < filteredData.sites; siteCount = siteCount+1)
	{
		siteMap = dupInfo[siteCount];
		if (doneSites[siteMap][0] == 0)
		{
			filterString = "";
			filterString = filterString + (siteCount*3) + "-" + (siteCount*3+2);
			DataSetFilter siteFilter = CreateFilter (ds,3,filterString,"",GeneticCodeExclusions);
			HarvestFrequencies (f1, siteFilter, 3, 3, 0);
			m1 = 0;
			for (mpiNode=0; mpiNode < 64; mpiNode=mpiNode+1)
			{
				if (f1[mpiNode]>0)
				{
					m1=m1+1;
				}
			}	
			
			if (m1>1)
			{
				LikelihoodFunction siteLikelihood = (siteFilter, siteTree);
				sFactor = 1;
				nFactor_Clade	= 1;
				nFactor_Tree	= 1;
				Optimize (site_res, siteLikelihood);
				doneSites[siteMap][0] = nFactor_Tree;
				doneSites[siteMap][1] = nFactor_Clade;
				doneSites[siteMap][2] = sFactor;
				nFactor_Clade:=nFactor_Tree;
				Optimize (site_resN, siteLikelihood);
				doneSites[siteMap][4] = sFactor;
				doneSites[siteMap][3] = nFactor_Tree;
				doneSites[siteMap][5] = 2*(site_res[1][0]-site_resN[1][0]);
				doneSites[siteMap][6] = 1-CChi2(doneSites[siteMap][5],1);		
			}
			else
			{
				doneSites[siteMap][0] = 1;
				doneSites[siteMap][6] = 1;
			}				
		}
		dummy = ReportSiteST (siteCount, siteMap);		
	}	
}
else
{
	MPINodeState = {MPI_NODE_COUNT-1,3};
	alreadyDone = {filteredData.unique_sites,1};
	for (siteCount = 0; siteCount < filteredData.sites; siteCount = siteCount+1)
	{
		siteMap = dupInfo[siteCount];
		if (alreadyDone[siteMap] == 0)
		{
			filterString = "";
			filterString = filterString + (siteCount*3) + "-" + (siteCount*3+2);
			DataSetFilter siteFilter = CreateFilter (filteredData,3,filterString,"",GeneticCodeExclusions);
			
			HarvestFrequencies (f1, siteFilter, 3, 3, 0);
			m1 = 0;
			for (mpiNode=0; mpiNode < 64; mpiNode=mpiNode+1)
			{
				if (f1[mpiNode]>0)
				{
					m1=m1+1;
				}
			}	
			
			alreadyDone[siteMap] = 1;				
			
			if (m1>1)
			{
				LikelihoodFunction siteLikelihood = (siteFilter, siteTree);				
				sFactor = 1;
				nFactor_Clade	= 1;
				nFactor_Tree	= 1;
				for (mpiNode = 0; mpiNode < MPI_NODE_COUNT-1; mpiNode = mpiNode+1)
				{
					if (MPINodeState[mpiNode][0]==0)
					{
						break;	
					}
				}
				
				if (mpiNode==MPI_NODE_COUNT-1)
				/* all nodes busy */
				{
					mpiNode = ReceiveJobsST (1,1);
				}
				else
				{
					MPISend (mpiNode+1,siteLikelihood);
					MPINodeState[mpiNode][0] = 1;
					MPINodeState[mpiNode][1] = siteCount;
					MPINodeState[mpiNode][2] = 1;
				}
				
				nFactor_Clade:=nFactor_Tree;

				for (mpiNode = 0; mpiNode < MPI_NODE_COUNT-1; mpiNode = mpiNode+1)
				{
					if (MPINodeState[mpiNode][0]==0)
					{
						break;	
					}
				}
				
				DataSetFilter siteFilter = CreateFilter (filteredData,3,filterString,"",GeneticCodeExclusions);
				LikelihoodFunction siteLikelihood = (siteFilter, siteTree);				
				if (mpiNode==MPI_NODE_COUNT-1)
				/* all nodes busy */
				{
					mpiNode = ReceiveJobsST (1,0);
				}
				else
				{
					MPISend (mpiNode+1,siteLikelihood);
					MPINodeState[mpiNode][0] = 1;
					MPINodeState[mpiNode][1] = siteCount;
					MPINodeState[mpiNode][2] = 0;
				}
			}
			else
			{
				doneSites[siteMap][6] = 1;
				ReportSiteST (siteCount, siteMap);				 
			}				
		}
	}					
	while (1)
	{
		for (nodeCounter = 0; nodeCounter < MPI_NODE_COUNT-1; nodeCounter = nodeCounter+1)
		{
			if (MPINodeState[nodeCounter][0]==1)
			{
				fromNode = ReceiveJobsST (0,0);
				break;	
			}
		}
		if (nodeCounter == MPI_NODE_COUNT-1)
		{
			break;
		}
	}					
	fprintf (stdout, "\n\n\n");
	for (siteCount = 0; siteCount < filteredData.sites; siteCount = siteCount+1)
	{
		siteMap = dupInfo[siteCount];
		dummy = ReportSiteST (siteCount, siteMap);				 
	}
}

OpenWindow (CHARTWINDOW,{{"ALS Results"}
					   {"labels"},
					   {"fullSites"},
					   {"Bar Chart"},
					   {"Index"},
					   {labels[0]},
					   {"Site Index"},
					   {""},
					   {labels[0]},
					   {"0"}},
					   "SCREEN_WIDTH-60;SCREEN_HEIGHT-50;30;50");
					   
SetDialogPrompt ("Save site-by-site LRT results to:");

siteCount = Columns (fullSites);
fprintf (PROMPT_FOR_FILE,CLEAR_FILE,labels[0]);
for (nodeCounter=1; nodeCounter<siteCount; nodeCounter=nodeCounter+1)
{
	fprintf (LAST_FILE_PATH,",",labels[nodeCounter]);
}

for (nodeCounter=0; nodeCounter < Rows (fullSites); nodeCounter = nodeCounter+1)
{
	fprintf (LAST_FILE_PATH,"\n",fullSites[nodeCounter][0]);
	for (mpiNode=1; mpiNode<siteCount; mpiNode=mpiNode+1)
	{
		fprintf (LAST_FILE_PATH,",",fullSites[nodeCounter][mpiNode]);
	}
}
