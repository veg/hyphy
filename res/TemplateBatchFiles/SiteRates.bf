/*------------------------------------------------------------------------*/

function ReportSite (siteI, siteM)
{
	fullSites[siteI][0] = doneSites[siteM][0];
	fullSites[siteI][1] = doneSites[siteM][1];


	fprintf (stdout, "Site ", Format(siteI+1,4,0),
					 " Rate = ", Format(fullSites[siteI][0],7,4),
					 " Log(L) ", Format(fullSites[siteI][1],7,4),"\n");		
					 
	return 0;
}

/*------------------------------------------------------------------------*/

function ReceiveJobs (sendOrNot)
{
	MPIReceive (-1, fromNode, result_String);
	
	siteIndex = MPINodeState[fromNode-1][1];
	
	if (sendOrNot)
	{
		MPISend (fromNode,siteLikelihood);
		MPINodeState[fromNode-1][1] = siteCount;			
	}
	else
	{
		MPINodeState[fromNode-1][0] = 0;
		MPINodeState[fromNode-1][1] = -1;		
	}
	
	siteMap = dupInfo[siteIndex];
	
	ExecuteCommands (result_String);
	
	doneSites[siteMap][0] = siteRate;
	doneSites[siteMap][1] = siteLikelihood_MLES[1][0];

	ReportSite (siteIndex, siteMap);
	
	return fromNode-1;
}

/*------------------------------------------------------------------------*/

SetDialogPrompt ("Please specify a nucleotide or amino-acid data file:");

DataSet ds = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter filteredData = CreateFilter (ds,1);

fprintf (stdout,"\n______________READ THE FOLLOWING DATA______________\n",ds);

categoriesUsed = 0;

SelectTemplateModel(filteredData);

if (categoriesUsed)
{
	fprintf (stdout,"\nPlease use only models without explicit rate variation components for this analysis\n");
}
else
{
	_DO_TREE_REBALANCE_ = 1;
	#include "queryTree.bf";
	LikelihoodFunction lf = (filteredData,givenTree);
	Optimize (resNull,lf);
	fprintf (stdout, "\n\nGlobal fit",lf);
	
	fprintf (stdout, "\n\nSite by site fits:\n\n");

	global siteRate = 1;
	siteRate:<100;
	
	doneSites    = {filteredData.unique_sites,2};
	fullSites    = {filteredData.sites,2};
	Tree		   siteTree = treeString;
	
	GetDataInfo    (dupInfo, filteredData);
	alreadyDone	= {};
	
	ReplicateConstraint ("this1.?.?:=siteRate*this2.?.?__",siteTree,givenTree);
	
	ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles" + DIRECTORY_SEPARATOR + "Utility" + DIRECTORY_SEPARATOR + "GrabBag.bf");
	fixGlobalParameters  ("lf");
		
	labels = {{"Rate","Log[L]"}};
	
	if (MPI_NODE_COUNT<=1)
	{
		VERBOSITY_LEVEL = -1;
		for (siteCount = 0; siteCount < filteredData.sites; siteCount = siteCount+1)
		{
			siteMap = dupInfo[siteCount];
			if (alreadyDone[siteMap] == 0)
			{
				filterString = "" + siteCount;
				DataSetFilter siteFilter = CreateFilter (ds,1,filterString);
				LikelihoodFunction siteLikelihood = (siteFilter, siteTree);
				siteRate = 1;
				Optimize (site_res, siteLikelihood);
				alreadyDone[siteMap]  = 1;
				doneSites[siteMap][0] = siteRate;
				doneSites[siteMap][1] = site_res[1][0];
			}
			dummy = ReportSite (siteCount, siteMap);				 
		}	
		VERBOSITY_LEVEL = 0;
	}
	else
	{
		MPINodeState = {MPI_NODE_COUNT-1,2};
		for (siteCount = 0; siteCount < filteredData.sites; siteCount = siteCount+1)
		{
			siteMap = dupInfo[siteCount];
			if (alreadyDone[siteMap] == 0)
			{
				filterString = "" + siteCount;
				DataSetFilter siteFilter = CreateFilter (ds,1,filterString);
				LikelihoodFunction siteLikelihood = (siteFilter, siteTree);				
				alreadyDone[siteMap] = 1;				
				siteRate = 1;
				
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
					mpiNode = ReceiveJobs (1);
				}
				else
				{
					MPISend (mpiNode+1,siteLikelihood);
					MPINodeState[mpiNode][0] = 1;
					MPINodeState[mpiNode][1] = siteCount;
				}
			}
		}					
		while (1)
		{
			for (nodeCounter = 0; nodeCounter < MPI_NODE_COUNT-1; nodeCounter = nodeCounter+1)
			{
				if (MPINodeState[nodeCounter][0]==1)
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
		fprintf (stdout, "\n\n");
		for (siteCount = 0; siteCount < filteredData.sites; siteCount = siteCount+1)
		{
			siteMap = dupInfo[siteCount];
			dummy = ReportSite (siteCount, siteMap);				 
		}
	}

	nodeCounter = 0;
	likelihoodBound = 0;

	for (siteCount = 0; siteCount < filteredData.sites; siteCount = siteCount+1)
	{
		nodeCounter 		= nodeCounter     + fullSites[siteCount][0];
		likelihoodBound		= likelihoodBound + fullSites[siteCount][1];
	}

	nodeCounter = nodeCounter/filteredData.sites; 

	for (siteCount = 0; siteCount < filteredData.sites; siteCount = siteCount+1)
	{
		fullSites[siteCount][0] = fullSites[siteCount][0]/nodeCounter;
	}

	OpenWindow (CHARTWINDOW,{{"Data Rates"}
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
							   

	fprintf (stdout, "\n\nLikelihood lower bound = ", resNull[1][0]," AIC=", 2*(resNull[1][1]-resNull[1][0]), "\n");
	fprintf (stdout, "\n\nApproximate likelihood upper bound = ", likelihoodBound," AIC=", 2*(resNull[1][1]+filteredData.unique_sites-likelihoodBound), "\n");


	SetDialogPrompt ("Save rate results to:");
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

	for (siteCount = 0; siteCount < resNull[1][2]; siteCount = siteCount+1)
	{
		GetString (globalVarName,lf,siteCount);
		ExecuteCommands (globalVarName+"="+globalVarName+";");
	}
	
}
