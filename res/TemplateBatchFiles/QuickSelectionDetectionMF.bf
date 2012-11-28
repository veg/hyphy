RequireVersion ("2.12");

ExecuteAFile ("qndhelper1_mf.ibf");
LoadFunctionLibrary ("GrabBag.bf");

REPLACE_TREE_STRUCTURE = 1;


ChoiceList  (rOptions,"dN/dS bias parameter options",1,NO_SKIP,
			 "Neutral","dN/dS=1",
			 "User","Custom dN/dS value",
			 "Estimate", "Estimate from data with branch corrections(slower).",
			 "Estimate + CI", "Estimate dN/dS and a 95% confidence interval from data (slower).",
			 "Estimate dN/dS only", "Estimate from data without branch corrections."
);
	

if (rOptions < 0) {
	return;
}

if (rOptions == 1) {
    RR = prompt_for_a_value ("Enter the value for dNdS",0.25,0,10000,0);
	dNdS = RR;
}

ChoiceList  (cOptions,"Which method?",1,NO_SKIP,
			 "SLAC","Single likelihood ancestor counting",
			 "FEL","Fixed effects likelihood",
			 "MEME","Mixed effects model of evolution");


if (cOptions < 0) {
	return;
}

ExecuteAFile("qndhelper2_mf.ibf");

fprintf (stdout, "\n\nPhase 4: Selection Analysis\n\n");


pooledFreqs = {4,1};

for (k=0; k<4; k+=1) {
	pooledFreqs[k] = (+positionFrequencies[k][-1])/3;
}

_EFV_MATRIX0_ = {{1,AC__*pooledFreqs[1],pooledFreqs[2],AT__*pooledFreqs[3]}
				{AC__*pooledFreqs[0],1,CG__*pooledFreqs[2],CT__*pooledFreqs[3]}
				{pooledFreqs[0],CG__*pooledFreqs[3],1,GT__*pooledFreqs[3]}
				{AT__*pooledFreqs[0],CT__*pooledFreqs[3],GT__*pooledFreqs[2],1}};

_EFV_MATRIX1_ = _EFV_MATRIX0_;
_EFV_MATRIX2_ = _EFV_MATRIX0_;				
useCustomCountingBias = 1;

if (cOptions == 0) {		
	ExecuteAFile("SGEmulator_MF.bf");
}
else { // FEL 
    pValue = prompt_for_a_value ("Significance level for Likelihood Ratio Tests",0.1,0,1,0)
					
	SHORT_MPI_RETURN = 1;
	FEL_RUN_TIMER    = Time(1);
	FEL_MASTER_TIMER = Time(0);
			
	ChoiceList  (brOptions,"Branch Options",1,NO_SKIP,
						 "All","Test for non-neutral evolution on all branches",
						 "Internal Only","Test for non-neutral evolution on internal branches only");
					 
	if (brOptions < 0) {
		return 0;
	}

	ExecuteCommands ("#include\"qndhelper3.ibf\";");				
	global			sFactor = 1;
	global			nFactor = 1;
	
	if (brOptions > 0)
	{
		doneSites    = {totalUniqueSites,8};
		fullSites    = {totalCodonCount ,8};						
		labels = {{"dN","dS","dN/dS","dS=dN","LRT","p-value","Full Log(L)","dN_other"}};
	}
	else
	{
		doneSites    = {totalUniqueSites,7};
		fullSites    = {totalCodonCount ,7};
		labels = {{"dN","dS","dN/dS","dS=dN","LRT","p-value","Full Log(L)"}};
	}
	
	if (MPI_NODE_COUNT>1)
	{
		MPINodeState = {MPI_NODE_COUNT-1,5};
	}
	
	vOffset  = 0;
	vuOffset = 0;
	
	alreadyDone = {totalUniqueSites,1};
	
	for (fileID = 1; fileID <= fileCount; fileID = fileID+1)
	{
		ClearConstraints (siteTree);
		Tree		   	  siteTree = treeStrings[fileID];
		
		ExecuteCommands ("GetDataInfo  (dupInfo, filteredData_"+fileID+");");			
		ExecuteCommands ("thisFilterSize  = filteredData_"+fileID+".sites;");			
		ExecuteCommands ("thisFilterSizeU = filteredData_"+fileID+".unique_sites;");			
		
		fprintf (stdout, "\nWorking on file ", fileID, " with ", thisFilterSize, " codons (", thisFilterSizeU, " unique patterns)\n");
		
		ExecuteCommands ("ReplicateConstraint(\"this1.?.synRate:=sFactor*this2.?.synRate__\",siteTree,codonTree_"+fileID+");");
		
		if (brOptions == 1)
		{
			global nFactorOther = 1;
			ExecuteCommands ("ReplicateConstraint(\"this1.Node?.nonSynRate:=nFactor*this2.Node?.synRate__\",siteTree,codonTree_"+fileID+");");
			ExecuteCommands ("ReplicateConstraint(\"this1.?.nonSynRate:=nFactorOther*this2.?.synRate__\",siteTree,codonTree_"+fileID+");");
		}
		else
		{
			ExecuteCommands ("ReplicateConstraint(\"this1.?.nonSynRate:=nFactor*this2.?.synRate__\",siteTree,codonTree_"+fileID+");");
		}
		
		lfSpawnDone = 0;
					
		if (MPI_NODE_COUNT<=1)
		{
			for (siteCount = 0; siteCount < thisFilterSize; siteCount = siteCount+1)
			{
				siteMap = dupInfo[siteCount];
				if (alreadyDone[siteMap+vuOffset] == 0)
				{
					alreadyDone[siteMap+vuOffset] = 1;					
					filterString 		 = "";
					filterString 		 = filterString + (siteCount*3) + "-" + (siteCount*3+2);
					ExecuteCommands ("DataSetFilter siteFilter = CreateFilter (ds_"+fileID+",3,filterString,\"\",GeneticCodeExclusions);");
					
					HarvestFrequencies (f1, siteFilter, 3, 3, 0);
					m1 = 0;
					for (mpiNode=0; mpiNode < 64; mpiNode=mpiNode+1)
					{
						if (f1[mpiNode]>0)
						{
							m1=m1+1;
						}
					}
					
					siteMap = siteMap + vuOffset;
					if (m1>1)
					{
						if (lfSpawnDone == 0)
						{
							LikelihoodFunction siteLikelihood = (siteFilter, siteTree);
							lfSpawnDone = 1;
						}
						
						sFactor = 1;
						nFactor	= 1;
						Optimize (site_res, siteLikelihood);
						doneSites[siteMap][0] = nFactor;
						doneSites[siteMap][1] = sFactor;
						m1 = sFactor+nFactor;
						sFactor 	= m1/2;
						nFactor:=sFactor;
						Optimize (site_resN, siteLikelihood);
						
						doneSites[siteMap][2] = sFactor;
						doneSites[siteMap][3] = 2*(site_res[1][0]-site_resN[1][0]);
						doneSites[siteMap][4] = 1-CChi2(doneSites[siteMap][3],1);															
						doneSites[siteMap][5] = site_res[1][0];															
						doneSites[siteMap][6] = doneSites[siteMap][0]/doneSites[siteMap][1];	
						if (brOptions > 0)
						{
							doneSites[siteMap][7] = nFactorOther;
						}													
					}
					else
					{	
						doneSites[siteMap][0] = 0;
						doneSites[siteMap][1] = 0;
						doneSites[siteMap][2] = 0;
						doneSites[siteMap][3] = 0;
						doneSites[siteMap][4] = 1;															
						doneSites[siteMap][5] = 0;
						if (brOptions > 0)
						{
							doneSites[siteMap][7] = 0;
						}
					}	
				}
				else
				{
					siteMap = siteMap + vuOffset;
				}
				dummy = ReportSite2 (siteCount+vOffset, siteMap);				 
			}	
		}
		else
		{
			
			for (siteCount = 0; siteCount < thisFilterSize; siteCount = siteCount+1)
			{
				siteMap = dupInfo[siteCount];
				if (alreadyDone[siteMap+vuOffset] == 0)
				{
					filterString = "";
					filterString = filterString + (siteCount*3) + "-" + (siteCount*3+2);
					ExecuteCommands ("DataSetFilter siteFilter = CreateFilter (ds_"+fileID+",3,filterString,\"\",GeneticCodeExclusions);");

					HarvestFrequencies (f1, siteFilter, 3, 3, 0);
					m1 = 0;
					for (mpiNode=0; mpiNode < 64; mpiNode=mpiNode+1)
					{
						if (f1[mpiNode]>0)
						{
							m1=m1+1;
						}
					}
					
					siteMap = siteMap + vuOffset;
					alreadyDone[siteMap] = 1;				
					if (m1>1)
					{
						sFactor = 1;
						nFactor	= 1;
						
						if (lfSpawnDone == 0)
						{
							LikelihoodFunction siteLikelihood = (siteFilter, siteTree);	
							LIKELIHOOD_FUNCTION_OUTPUT = 7;
							lfSpawnDone = 1;
						}
									
						for (mpiNode = 0; mpiNode < MPI_NODE_COUNT-1; mpiNode = mpiNode+1)
						{
							if (MPINodeState[mpiNode][0]==0)
							{
								break;	
							}
						}
						
						if (mpiNode==MPI_NODE_COUNT-1)
						{
							mpiNode = ReceiveJobs2 (1,1,siteCount+vOffset,siteMap);
						}
						else
						{
							MPISend (mpiNode+1,siteLikelihood);
							MPINodeState[mpiNode][0] = 1;
							MPINodeState[mpiNode][1] = siteCount+vOffset;
							MPINodeState[mpiNode][2] = 1;
							MPINodeState[mpiNode][3] = siteMap;
							MPINodeState[mpiNode][4] = MPINodeState[mpiNode][4] + 1;
							
						}
						
						m1 = sFactor+nFactor;
						sFactor 	= m1/2;
						nFactor:=sFactor;

						for (mpiNode = 0; mpiNode < MPI_NODE_COUNT-1; mpiNode = mpiNode+1)
						{
							if (MPINodeState[mpiNode][0]==0)
							{
								break;	
							}
						}
						
						if (mpiNode==MPI_NODE_COUNT-1)
						{
							mpiNode = ReceiveJobs2 (1,0,siteCount+vOffset,siteMap);
						}
						else
						{
							MPISend (mpiNode+1,siteLikelihood);
							MPINodeState[mpiNode][0] = 1;
							MPINodeState[mpiNode][1] = siteCount+vOffset;
							MPINodeState[mpiNode][3] = siteMap;
							MPINodeState[mpiNode][2] = 0;
						}
					}
					else
					{
	                    doneSites[siteMap][0] = 0;
	                    doneSites[siteMap][1] = 0;
	                    doneSites[siteMap][2] = 0;
	                    doneSites[siteMap][3] = 0;
	                    doneSites[siteMap][4] = 1;									
					}
				}
			}					
		}
		vOffset 	= vOffset  + thisFilterSize;
		vuOffset 	= vuOffset + thisFilterSizeU;
	}
		
	if (MPI_NODE_COUNT>1)
	{
		while (1)
		{
			for (nodeCounter = 0; nodeCounter < MPI_NODE_COUNT-1; nodeCounter = nodeCounter+1)
			{
				if (MPINodeState[nodeCounter][0]==1)
				{
					fromNode = ReceiveJobs2 (0,0,0,0);
					break;	
				}
			}
			if (nodeCounter == MPI_NODE_COUNT-1)
			{
				break;
			}
		}					
		fprintf (stdout, "\n\n\n");
		vOffset  = 0;
		vuOffset = 0;
		
		alreadyDone = {totalUniqueSites,1};
		
		for (fileID = 1; fileID <= fileCount; fileID = fileID+1)
		{
			ClearConstraints (siteTree);
			ExecuteCommands ("GetDataInfo  (dupInfo, filteredData_"+fileID+");");			
			ExecuteCommands ("thisFilterSize  = filteredData_"+fileID+".sites;");			
			ExecuteCommands ("thisFilterSizeU = filteredData_"+fileID+".unique_sites;");			
			for (siteCount = 0; siteCount < thisFilterSize; siteCount = siteCount+1)
			{
				siteMap = dupInfo[siteCount];
				dummy = ReportSite2 (siteCount+vOffset, siteMap+vuOffset);				 
			}
			vOffset 	= vOffset  + thisFilterSize;
			vuOffset 	= vuOffset + thisFilterSizeU;
		}
	}

	OpenWindow (CHARTWINDOW,{{"FEL Results"}
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
							   
	SHORT_MPI_RETURN = 0;
							   
	fprintf (stdout, "\n\nWall-clock run-time (seconds): ", Time(1) - FEL_RUN_TIMER, 
					 "\nProcess time on the master node: ", Time(0) - FEL_MASTER_TIMER,
					 "\nMPI nodes used: ", MPI_NODE_COUNT, "\n",
					 "\nWall-clock time per site (seconds): ", (Time(1) - FEL_RUN_TIMER)/totalUniqueSites, "\n");
					 
			
	SetDialogPrompt ("Save site-by-site LRT results to:");
	
	siteCount = Columns (fullSites);
	fprintf (PROMPT_FOR_FILE,CLEAR_FILE,labels[0]);
	
	outString = "";
	outString * 8192;
	
	for (nodeCounter=1; nodeCounter<siteCount; nodeCounter=nodeCounter+1)
	{
		outString * ("," + labels[nodeCounter]);
	}
	
	for (nodeCounter=0; nodeCounter < Rows (fullSites); nodeCounter = nodeCounter+1)
	{
		outString * ("\n"+fullSites[nodeCounter][0]);
		for (mpiNode=1; mpiNode<siteCount; mpiNode=mpiNode+1)
		{
			outString * (","+fullSites[nodeCounter][mpiNode]);
		}
	}
	outString * 0;
	fprintf (LAST_FILE_PATH, outString);
	outString = 0;
}

useCustomCountingBias = 0;
