/* ____________________________________________*/

function TreeMatrix2TreeString (levelIndex)
{
	treeString = "";
	p = 0;
	k = 0;
	m = treeNodes[0][levelIndex+1];
	n = treeNodes[0][levelIndex];

	while (m)
	{	
		if (m>p)
		{
			if (p)
			{
				treeString = treeString+",";
			}
			for (j=p;j<m;j=j+1)
			{
				treeString = treeString+"(";
			}
		}
		else
		{
			if (m<p)
			{
				for (j=m;j<p;j=j+1)
				{
					treeString = treeString+")";
				}
			}	
			else
			{
				treeString = treeString+",";
			}	
		}
		if (n<ds.species)
		{
			GetString (nodeName, ds, n);
			treeString = treeString+nodeName;
		}
		k=k+1;
		p=m;
		n=treeNodes[k][levelIndex];
		m=treeNodes[k][levelIndex+1];
	}

	for (j=m;j<p;j=j+1)
	{
		treeString = treeString+")";
	}
	
	return treeString;
}


/* ____________________________________________*/

function  _PrepareForTreeSearch (treesToBeSearched)
{
	bestTreesStash    = {10,2};
	globalTreeCounter = 0;
	treeStatistics    = {treesToBeSearched, 1};
	for (ii=0; ii<10; ii=ii+1)
	{
		bestTreesStash [ii][1] = -1e100;
		bestTreesStash [ii][0] = "";
	}
	return 1;
}

/* ____________________________________________*/

function  _AddTreeToResults		(currentTreeString, currentLFValue)
{
	treeStatistics [globalTreeCounter][0] = currentLFValue;
	globalTreeCounter = globalTreeCounter+1;
	
	for (ii = 0; ii<10; ii=ii+1)
	{
		if (currentLFValue>bestTreesStash[ii][1])
		{
			break;
		}
	}
	if (ii<10)
	{
		for (ii2 = 8; ii2>=ii; ii2=ii2-1)
		{
			bestTreesStash [ii2+1][1] = bestTreesStash[ii2][1];
			bestTreesStash [ii2+1][0] = bestTreesStash[ii2][0];
		}
		bestTreesStash [ii][0] = currentTreeString;
		bestTreesStash [ii][1] = currentLFValue;
	}
	return 1;
}

/*---------------------------------------------------------------*/

function ReceiveJobs (sendOrNot)
{
	if (MPI_NODE_COUNT>1)
	{
		MPIReceive (-1, fromNode, result_String);
		
		saveTree = MPINodeTree[fromNode-1];
		saveIdx	 = MPINodeInfo[fromNode-1][1];
		
		if (sendOrNot)
		{
			MPISend 	(fromNode,current_MPIJOB);
			MPINodeTree[fromNode-1]    = thisTree;
			MPINodeInfo[fromNode-1][1] = treeCounter;
		}
		else
		{
			MPINodeInfo[fromNode-1][0]  = 0;
			MPINodeInfo[fromNode-1][1]  = -1;
			MPINodeTree[fromNode-1]     = "";
		}
		
		thisTree    = saveTree;
		ExecuteCommands ("currentLL = " + result_String + ";");
				
		saveCounter = treeCounter;
		treeCounter	= saveIdx;
	}
	else
	{
		currentLL = res[1][0];
	}
	
	dummy = _AddTreeToResults (thisTree, currentLL);
	if (currentLL>bestValue)
	{
		bestValue = currentLL;
		bestTree = thisTree;
	}
	fprintf (stdout, " ==> logLhd = ", currentLL);
	
	if (MPI_NODE_COUNT>1)
	{
		treeCounter = saveCounter;
	}
	
	return fromNode-1;
}

/* ____________________________________________*/

function  _ReportTreeStatistics		(currentLFValue)
{
	ii = 0;
	fprintf (stdout, "\n\n**************************\n",
					     "*     TREE REPORT	       *\n",
					     "**************************\n\n");
					     
	fprintf (stdout, "\n#### BEST TREES #####\n\n");
					     
	for (ii=0; ii<10; ii = ii+1)
	{
		if (bestTreesStash[ii][1]==(-1e100))
		{
			break;
		}
		fprintf (stdout, ii+1, ").");
		
		if (ii>0)
		{
			fprintf (stdout, " Worse by: ", bestTreesStash[ii][1]-currentLFValue);
		}
		fprintf (stdout,"\n",  bestTreesStash[ii][0], "\nLog-likelihood = ", bestTreesStash[ii][1], "\n\n");
	}
	
	fprintf (stdout, "\n#### STATISTICS #####\n\n");
	
	bestTreesStash [0][0] = 0.1;
	bestTreesStash [1][0] = 0.5;
	bestTreesStash [2][0] = 1;
	bestTreesStash [3][0] = 5;
	bestTreesStash [4][0] = 10;
	bestTreesStash [5][0] = 50;
	bestTreesStash [6][0] = 100;
	bestTreesStash [7][0] = 1000;
	bestTreesStash [8][0] = 10000;
	bestTreesStash [9][0] = 1e100;

	bestTreesStash [0][1] = 0;
	bestTreesStash [1][1] = 0;
	bestTreesStash [2][1] = 0;
	bestTreesStash [3][1] = 0;
	bestTreesStash [4][1] = 0;
	bestTreesStash [5][1] = 0;
	bestTreesStash [6][1] = 0;
	bestTreesStash [7][1] = 0;
	bestTreesStash [8][1] = 0;
	bestTreesStash [9][1] = 0;
	
	probabilityOfTheData = 0;
	
	for (i=0; i<globalTreeCounter; i=i+1)
	{
		diff = currentLFValue-treeStatistics[i];
		j = 0;
		while (diff>bestTreesStash[j][0])
		{
			j=j+1;
		}
		bestTreesStash [j][1] = bestTreesStash [j][1] + 1;
		probabilityOfTheData = probabilityOfTheData+Exp(-diff);
	}
	
	bestTreesStash [0][1] = bestTreesStash [0][1]-1;
	
	ii = "+---------------+---------------+---------------+---------------+\n";
	fprintf (stdout, "\n\n", ii, 
							    "| From Best +   |  To Best +    |   Tree Count  |  % of total	  |\n",
							 ii);
	for (i=0; i<10; i=i+1)
	{	
		if (i)
		{
			fprintf (stdout, "| " , Format (bestTreesStash [i-1][0],13,1));
		}
		else
		{
			fprintf (stdout, "|             0");
		}
		if (i<9)
		{
			fprintf (stdout, " | " , Format (bestTreesStash [i][0],13,1));
		}
		else
		{
			fprintf (stdout, " |      Infinity");
		}		
		fprintf (stdout, " | ", Format (bestTreesStash [i][1],13,0), " | ", Format (100*bestTreesStash [i][1]/globalTreeCounter,13,8), " |\n",ii); 
	}
	
	fprintf (stdout, "\n\nPosterior probability of the best tree (with uninformative prior) = ",1./probabilityOfTheData,"\n\n");
	
	fprintf (stdout, "\n\n***********Save full tree statistics to a file (y/n)?");

	fscanf  (stdin, "String", resp);

	if ((resp!="n")&&(resp!="N"))
	{
		SetDialogPrompt ("Write tree stats string to:");
		fprintf (PROMPT_FOR_FILE,CLEAR_FILE,treeStatistics);
	}
	treeStatistics = 0;

	return 1;
}

/* ____________________________________________*/

MESSAGE_LOGGING = 0;
VERBOSITY_LEVEL = -1;

SetDialogPrompt ("Please choose a nucleotide or amino-acid data file:");

DataSet ds 				   = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter filteredData = CreateFilter (ds,1,"","");
SelectTemplateModel			(filteredData);

if (MPI_NODE_COUNT > 1)
{
	mpiJobPrefix = "";
	mpiJobPrefix * 128;
	mpiJobPrefix * ("DataSet ds = ReadDataFile (\""+LAST_FILE_PATH+"\");DataSetFilter filteredData = CreateFilter (ds,1);");
	Export		 (modelExp, USE_LAST_MODEL);
	mpiJobPrefix * modelExp;
	mpiJobPrefix * 0;
	mpiJobSuffix = "LikelihoodFunction lf = (filteredData,givenTree);Optimize(res,lf);return res[1][0];";
}

treeNodes 		= {2*(ds.species+1),2*(ds.species-2)};
cladesInfo 		= {ds.species,2*(ds.species-2)};

branchIndex		= {ds.species-3,1};
currentLevel 	= 0;

done 			= false;

i = 2*ds.species-5;
j = 1;
while (i>1)
{
	j = j*i;
	i = i-2;
}

dummy = _PrepareForTreeSearch (j);

treeNodes[0][0]=0;
treeNodes[0][1]=1;
treeNodes[1][0]=1;
treeNodes[1][1]=1;
treeNodes[2][0]=2;
treeNodes[2][1]=1;
treeNodes[3][0]=ds.species;
treeNodes[3][1]=0;
cladesInfo[0][0]=0;
cladesInfo[0][1]=4;

bestTree ="";
bestValue=-1e20;

done = 0;
treeCounter = 0;

if (MPI_NODE_COUNT>1)
{
	MPINodeInfo  = {MPI_NODE_COUNT-1,2};
	MPINodeTree  = {MPI_NODE_COUNT-1,1};
	MPINodeTree[0]  = "";
}

while (!done)
{
	if (branchIndex[currentLevel]<2*currentLevel+3)
	{
		i = 0;
		shift = 0;
		j = 2*currentLevel;
		k = j+2;
		m = j+1;
		while (treeNodes[i][m])
		{
			/*copy tree from prev level to this level */
			if (i==branchIndex[currentLevel])
			/*insert new branch*/
			{
				shift = 2;
				if (treeNodes[i][j]<ds.species)
				/* simple branch */
				{
					treeNodes[i][k]=treeNodes[i][j];
					treeNodes[i][k+1]=treeNodes[i][m]+1;
					treeNodes[i+1][k]=currentLevel+3;
					treeNodes[i+1][k+1]=treeNodes[i][m]+1;
					treeNodes[i+2][k]=currentLevel+ds.species+1;
					treeNodes[i+2][k+1]=treeNodes[i][m];
					cladesInfo[currentLevel+1][k] = i;
					cladesInfo[currentLevel+1][k+1] = 3;					
				}
				else
				{
					/* update node depths for the entire clade now*/
					l = treeNodes[i][j]-ds.species;
					s = cladesInfo[l][j];
					for (p=s+cladesInfo[l][m]-1; p>=s; p=p-1)
					{
						treeNodes[i][k]=treeNodes[i][j];
						treeNodes[i][k+1]=treeNodes[i][m]+1;						
						i=i-1;
					}
					i=i+cladesInfo[l][m];
					/* new clade record */
					cladesInfo[currentLevel+1][k] = cladesInfo[l][j];
					cladesInfo[currentLevel+1][k+1] = cladesInfo[l][m]+2;
					/* now we need to insert two more nodes */
					treeNodes[i+1][k]=currentLevel+3;
					treeNodes[i+1][k+1]=treeNodes[i][m]+1;
					treeNodes[i+2][k]=currentLevel+ds.species+1;
					treeNodes[i+2][k+1]=treeNodes[i][m];
				}
				for (p=0; p<=currentLevel; p=p+1)
				{
					if (cladesInfo[p][j]>i)
					{
						cladesInfo[p][k] = cladesInfo[p][j]+2;
					}
					else
					{
						cladesInfo[p][k] = cladesInfo[p][j];
					}
					
					if ((cladesInfo[p][j]<=i)&&((cladesInfo[p][j]+cladesInfo[p][m])>i+1))
					{
						cladesInfo[p][k+1] = cladesInfo[p][m]+2;
					}
					else
					{
						cladesInfo[p][k+1] = cladesInfo[p][m];
					}
				}				
			}
			else
			{
				treeNodes[i+shift][k]=treeNodes[i][j];
				treeNodes[i+shift][k+1]=treeNodes[i][m];
			}
			i = i+1;
		}
		treeNodes[i+2][k]=treeNodes[i][j];
		treeNodes[i+2][k+1]=treeNodes[i][j+1];
		if (currentLevel<ds.species-4)
		{
			currentLevel = currentLevel+1;
		}
		else
		{
			thisTree = TreeMatrix2TreeString (2*(currentLevel+1));
			branchIndex[currentLevel]=branchIndex[currentLevel]+1;
			fprintf (stdout, "\nTree#",Format(treeCounter,0,0)," ", thisTree);
			current_MPIJOB = mpiJobPrefix + "Tree givenTree = " + thisTree + ";" + mpiJobSuffix;
			if (MPI_NODE_COUNT>1)
			{
				for (mpiNode = 0; mpiNode < MPI_NODE_COUNT-1; mpiNode = mpiNode+1)
				{
					if (MPINodeInfo[mpiNode][0]==0)
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
					MPINodeInfo[mpiNode][0] = 1;
					MPINodeInfo[mpiNode][1] = treeCounter;
					MPINodeTree[mpiNode] = thisTree;
					MPISend (mpiNode+1,current_MPIJOB);
				}
			}
			else
			{
				Tree    treeVar = thisTree;
				LikelihoodFunction lf = (filteredData, treeVar);
				Optimize (res, lf);
				ReceiveJobs (0);
			}
			treeCounter = treeCounter+1;
		}
	}
	else
	{
		branchIndex[currentLevel]=0;
		if (currentLevel==0)
		{
			done = 1;
		}
		else
		{
			currentLevel = currentLevel-1;
			branchIndex[currentLevel]=branchIndex[currentLevel]+1;
		}
	}
}

if (MPI_NODE_COUNT>1)
{
	while (1)
	{
		for (nodeCounter = 0; nodeCounter < MPI_NODE_COUNT-1; nodeCounter = nodeCounter+1)
		{
			if (MPINodeInfo[nodeCounter][0]==1)
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


VERBOSITY_LEVEL = 0;

fprintf (stdout,"\n\n --------------------- RESULTS --------------------- \n\n");
fprintf (stdout,"\n\n BestTree =", bestTree);

dummy = _ReportTreeStatistics (bestValue);

Tree	tr = bestTree;
LikelihoodFunction lf = (filteredData, tr);
Optimize (res,lf);
fprintf (stdout, "\n",lf, "\n\n***********Save this tree to a file (y/n)?");

fscanf  (stdin, "String", resp);

if ((resp!="n")&&(resp!="N"))
{
	SetDialogPrompt ("Write tree string to:");
	fprintf (PROMPT_FOR_FILE,CLEAR_FILE,bestTree,";");
}

