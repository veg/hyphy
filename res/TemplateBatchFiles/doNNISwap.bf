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

while(1)
{
	for (l2 = 0; l2 < l; l2=l2+1)
	{
		dummy = splitCladesOnIBranch (filteredData.species,filteredData.species+l2);
		for (l3=2; l3<4; l3=l3+1)
		{
			dummy = doNNISwap (l3,filteredData.species);
			thisTree = TreeMatrix2TreeString (2,0);
			if (verbFlag )
			{	
				if (MPI_NODE_COUNT<2)
				{
					fprintf (stdout, "\n\tNNI tree ",thisTree);
				}
				else
				{
					fprintf (stdout, "\n\tNNI swap ",l3-1, " at node ", l2, "\n");			
				}
			}
			Tree    Inferred_Tree = thisTree;
			if (haveTreeConstraint)
			{
				if ((Inferred_Tree<=_topologyPattern) == 0)
				{
					if (verbFlag)
					{
						fprintf (stdout, " => Rejected by the topology constraint.");
					}
					continue;
				}
				else
				{
					if (verbFlag)
					{
						fprintf (stdout, " => Accepted by the topology constraint.");
					}
				}
			}
			if (starDecomposition)
			{
				LikelihoodFunction lf = (filteredData, Inferred_Tree);
				Optimize (res,lf);
				if (verbFlag)
				{
					fprintf (stdout, " => ", res[1][0]);
				}
				if (res[1][0]>bestValue)
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
					if (verbFlag)
					{
						OpenWindow (TREEWINDOW, {{"Inferred_Tree"}});
						fprintf (stdout, "\n\tRestarting NNI on the better tree.");
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
			else
			{
				dummy = SendOffTreeJob (1);
			}
		}
	}

	if ((starDecomposition==0)&&(MPI_NODE_COUNT>1))
	{
		inBestValue = bestValue;
		while (1)
		{
			for (nodeCounter = 0; nodeCounter < MPI_NODE_COUNT-1; nodeCounter = nodeCounter+1)
			{
				if (MPINodeState[nodeCounter]==1)
				{
					fromNode = ReceiveJobs (0,1);
					break;	
				}
			}
			if (nodeCounter == MPI_NODE_COUNT-1)
			{
				break;
			}
		}	
		if (bestValue>inBestValue)
		{
			l = 0;
			while (cladesInfo[l][1])
			{
				l = l+1;
			}				
			l2 = -1;
			continue;
		}
	}	
	break;
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
