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
	if (Abs(statFileName)>0)
	{
		fprintf (statFileName, currentLFValue, ";", currentTreeString, "\n");
	}
	
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
	
	for (i=0; i<globalTreeCounter; i=i+1)
	{
		diff = currentLFValue-treeStatistics[i];
		j = 0;
		while (diff>bestTreesStash[j][0])
		{
			j=j+1;
		}
		bestTreesStash [j][1] = bestTreesStash [j][1] + 1;
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
	
	treeStatistics = 0;

	return 1;
}

/* ____________________________________________*/

MESSAGE_LOGGING = 0;

SetDialogPrompt ("Please choose a nucleotide or amino-acid data file:");
DataSet ds = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter filteredData = CreateFilter (ds,1,"","");

SetDialogPrompt ("Please choose a tree constraint:");
fscanf			(PROMPT_FOR_FILE, "String", Tree_Constraint_String);

Tree			Tree_Constraint = Tree_Constraint_String;

SelectTemplateModel(filteredData);
treeNodes = {2*(ds.species+1),2*(ds.species-2)};
cladesInfo = {ds.species,2*(ds.species-2)};
branchIndex= {ds.species-3,1};
currentLevel = 0;

done = false;

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

fprintf (stdout, "\n\n***********Save full tree statistics to a file (y/n)?");
fscanf  (stdin, "String", resp);

statFileName = "";

if ((resp!="n")&&(resp!="N"))
{
	SetDialogPrompt ("Write tree stats strings to:");
	fprintf (PROMPT_FOR_FILE,CLEAR_FILE);
	statFileName = LAST_FILE_PATH;
}

fprintf (stdout, "\n\n*********** RUNNING CONSTRAINED TREE SEARCH ***********\nConstraint:", Tree_Constraint_String, "\n\n");


treeCounter = 0;

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
			Tree    treeVar = thisTree;
			if (treeVar <= Tree_Constraint)
			{
				fprintf (stdout, "\nTree#",Format(treeCounter,0,0)," ", thisTree);
				LikelihoodFunction lf = (filteredData, treeVar);
				Optimize (res,lf);
				dummy = _AddTreeToResults (thisTree, res[1][0]);
				if (res[1][0]>bestValue)
				{
					bestValue = res[1][0];
					bestTree = thisTree;
				}
				fprintf (stdout, " ==> logLhd = ", res[1][0]);
				treeCounter = treeCounter+1;
			}
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

if (globalTreeCounter)
{
	fprintf (stdout,"\n\n --------------------- RESULTS --------------------- \n\n");
	fprintf (stdout," Total tree count =", Format(treeCounter,0,0));
	fprintf (stdout,"\n\n BestTree =", bestTree);
	fprintf (stdout,"\n\nTree Constraint:\n",Tree_Constraint,"\n");
 

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
}
else
{
	fprintf (stdout, "\n\n***********************************\n",
					     "* NO TREES MATCHED THE CONSTRAINT	*\n",
					     "***********************************\n\n");
}


