l = 0;
while (cladesInfo[l][1])
{
	l = l+1;
}

if (verbFlag)
{
	fprintf (stdout, "\n\n\tPerforming NNI branch swapping on the best tree so far.\n\tAt least ",l*2," trees will be examined.\n\n");
}

nniDidBetter = 0;

for (l2 = 0; l2 < l; l2=l2+1)
{
	dummy = splitCladesOnIBranch (_NUMBER_OF_SEQUENCES,_NUMBER_OF_SEQUENCES+l2);
	for (l3=2; l3<4; l3=l3+1)
	{
		dummy = doNNISwap (l3,_NUMBER_OF_SEQUENCES);
		thisTree = TreeMatrix2TreeString (2,0);
		if (verbFlag)
		{
			fprintf (stdout, "\n\tNNI tree ",thisTree);
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
			break;
		}
	}
}
if (nniDidBetter)
{
	for (l=2*taxonCounter-2;l>=0;l=l-1)
	{
		treeNodes[l][0]=bestTreeNodes[l][0];
		treeNodes[l][1]=bestTreeNodes[l][1];
	}
	for (l=0; l<taxonCounter-2; l=l+1)
	{
		cladesInfo[l][0]=bestCladesInfo[l][0];
		cladesInfo[l][1]=bestCladesInfo[l][1];
	}
	if (verbFlag)
	{
		fprintf (stdout,"\n\n\t Improved NNI ", Format(taxonCounter,0,0)," taxa tree is ", bestTree, " with log likelihood of ", bestValue);
	}
}
