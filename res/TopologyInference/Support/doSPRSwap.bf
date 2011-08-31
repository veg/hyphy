l = 0;
while (cladesInfo[l][1])
{
	l = l+1;
}

if (verbFlag)
{
	fprintf (stdout, "\n\n\tPerforming subtree pruning and regrafting on the best tree so far.\n\n");
}
nniDidBetter = 0;
for (l2 = 0; l2 < l; l2=l2+1)
{
	dummy = splitCladesOnIBranch (_NUMBER_OF_SEQUENCES,_NUMBER_OF_SEQUENCES+l2);
	l3s = 0;
	l3f = 2;
	if ((splitTreeInfo[0][1]==1)&&(splitTreeInfo[0][3]==1))
	{
		l3s = 2;
		l3f = 4;
	}
	
	for (l3=l3s; l3<l3f; l3=l3+1)
	{
		l4 = 0;
		while (splitTreeInfo[l4+1][2*l3+1])
		{
			/* check for degeneracies 
			if ((l3==1)&&(splitTreeInfo[0][1]==1)&&(splitTreeInfo[0][3]==1))
			{
				break;
			}*/
			if (l3<2)
			{
				dummy = doSPRSwap (2,3,l3,l4,_NUMBER_OF_SEQUENCES);
			}
			else
			{
				dummy = doSPRSwap (0,1,l3,l4,_NUMBER_OF_SEQUENCES);								
			}
			
			thisTree = TreeMatrix2TreeString (2,0);
			if (verbFlag)
			{
				fprintf (stdout, "\n\tSPR tree ",thisTree);
			}
			Tree    Inferred_Tree = thisTree;
			if (starDecomposition)
			{
				SpawnLikelihoodFunction ("_INF_LF_", "Inferred_Tree", INFERENCE_DATA_WINDOW, _all_sequence_matrix);
			}		
			else
			{
				SpawnLikelihoodFunction ("_INF_LF_", "Inferred_Tree", INFERENCE_DATA_WINDOW, _subset_to_filter);
			}

			if (globalParametersOption && Abs(_cString))
			{	
				ExecuteCommands (_cString);	
			}
			Optimize (res,_INF_LF_);
			
			if (verbFlag)
			{
				fprintf (stdout, " => ", res[1][0]);
			}
			if (res[1][0]>bestValue+OPTIMIZATION_PRECISION)
			{	
				nniDidBetter = 1;
				bestValue = res[1][0];
				bestTree = thisTree;
				for (l4=2*taxonCounter-2;l4>=0;l4=l4-1)
				{
					bestTreeNodes[l4][0]= treeNodes[l4][2];
					bestTreeNodes[l4][1]= treeNodes[l4][3];
					treeNodes[l4][0]    = treeNodes[l4][2];
					treeNodes[l4][1]    = treeNodes[l4][3];
				}
				for (l4=0; l4<taxonCounter-2; l4=l4+1)
				{
					bestCladesInfo[l4][0]= cladesInfo[l4][2];
					bestCladesInfo[l4][1]= cladesInfo[l4][3];
					cladesInfo[l4][0]    = cladesInfo[l4][2];
					cladesInfo[l4][1]    = cladesInfo[l4][3];
				}
				l = 0;
				while (cladesInfo[l][1])
				{
					l = l+1;
				}				
				l2 = -1;
				l3 = l3f;
				break;
			}
			l4 = l4+1;
		}
	}
}
if (nniDidBetter)
{
	for (l=Rows(bestTreeNodes)-1;l>=0;l=l-1)
	{
		treeNodes[l][0]=bestTreeNodes[l][0];
		treeNodes[l][1]=bestTreeNodes[l][1];
	}
	for (l=0; l<Rows(bestCladesInfo); l=l+1)
	{
		cladesInfo[l][0]=bestCladesInfo[l][0];
		cladesInfo[l][1]=bestCladesInfo[l][1];
	}
	if (verbFlag)
	{
		fprintf (stdout,"\n\n\t Improved SPR ", Format(taxonCounter,0,0)," taxa tree is ", bestTree, " with log likelihood of ", bestValue);
	}
}
