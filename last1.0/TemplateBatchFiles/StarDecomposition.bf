MESSAGE_LOGGING = 0;
twHasBeenOpened = 0;

#include "heuristicMethodNPBootstrap.bf";
#include "branchSwappingFunctions.bf";

ChoiceList (doNNIOption,"Branch Swapping",1,SKIP_NONE,
			"No Swapping","No branch swapping is performed.",
			"Global NNI","Nearest neighbor interchange is performed after ALL the sequences have been added. Order (sequences)^1 additional trees are examined.",
			"Global SPR","Subtree pruning and regrafting is performed after ALL the sequences have been added. Order (sequences)^2 additional trees are examined.");

if (doNNIOption<0)
{
	return;
}
ChoiceList (dataType,"Data type",1,SKIP_NONE,"Nucleotide/Protein","Nucleotide or amino-acid (protein).",
				     "Codon","Codon (several available genetic codes).");

if (dataType<0) 
{
	return;
}

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


DataSet ds = ReadDataFile (PROMPT_FOR_FILE);

if (dataType)
{
	DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);
}
else
{
	DataSetFilter filteredData = CreateFilter (ds,1,"","");
}

haveTreeConstraint = 0;

function InferTreeTopology (verbFlag)
{
	treeNodes = {2*(ds.species+1),4};
	/* stores previous best tree in columns 0 and 1, and current working tree in columns 2,3 */
	bestTreeNodes = {2*(ds.species+1),2};
	/* holder for current best tree */

	cladesInfo = {ds.species,4};
	/* stores previous best tree clade info in columns 0 and 1, and current working tree clade info in columns 2,3 */
	bestCladesInfo = {ds.species,2};
	/* holder for current best tree */

	eligibleNodesCount = ds.species;
	eligibleNodes = {ds.species,1};

	done = false;

	for (k=0; k<ds.species; k=k+1)
	{
		treeNodes[k][0] = k;
		treeNodes[k][1] = 1;
		eligibleNodes[k] = k;
	}

	treeNodes[ds.species][0] = ds.species;
	treeNodes[ds.species][1] = 0;

	cladesInfo[0][0]=0;
	cladesInfo[0][1]=ds.species+1;

	bestTree="";
	bestValue=-1e20;

	currentTaxon = 0;
	taxonCounter = 3;
	cladesCount = 1;

	if (verbFlag)
	{
		fprintf (stdout, "**** Running star decomposition on ", LAST_FILE_PATH, " ****");
	}
	while (eligibleNodesCount>3)
	{
		bestNode1 = 0;
		bestNode2 = 0;
		/* stores the nodes which will be grouped at this stage */
		
		if (verbFlag)
		{
			fprintf (stdout, "\n\n\nSTEP ",Format(ds.species-eligibleNodesCount+1,0,0),"/",Format(ds.species-3,0,0),"). ");
		}
		bestValue=-1e20;
		if (eligibleNodesCount>4)
		{
			loop1 = 0;
			if (verbFlag)
			{
				fprintf (stdout, Format((eligibleNodesCount-1)*eligibleNodesCount/2,0,0), " trees for this step.\n");
			}
		}
		else
		{
			loop1 = 1;
			if (verbFlag)
			{
				fprintf (stdout, Format(3,0,0), " trees for this step.\n");
			}
		}
		for (; loop1 < eligibleNodesCount-1; loop1 = loop1+1)
		{
			for (loop2 = loop1+1; loop2 < eligibleNodesCount; loop2 = loop2+1)
			{
				n = eligibleNodes[loop1];
				m = eligibleNodes[loop2];
				
				/* 3 cases to consider */
				/* first - see what happens when both nodes grouped are leaves */
				if ((n<ds.species)&&(m<ds.species))
				{
					/* a few things to do here:
					   copy all nodes from the old tree to the new tree excluding the two being added
					   add a new clade at the end of the tree
					   update clade lengths of the original tree
					*/
					i=0;
					shift = 0;
					/* copying nodes other than n and m */
					while (treeNodes[i][1])
					{
						p = treeNodes[i][0];
						if ((p!=m)&&(p!=n))
						{
							treeNodes[i-shift][2]=treeNodes[i][0];
							treeNodes[i-shift][3]=treeNodes[i][1];
						}
						else
						{
							shift = shift+1;
						}
						i = i+1;
					}
					
					i = i-2;
					treeNodes[i][2]=m;
					treeNodes[i][3]=2;
					i = i+1;
					treeNodes[i][2]=n;
					treeNodes[i][3]=2;
					i = i+1;
					treeNodes[i][2]=cladesCount+ds.species;
					treeNodes[i][3]=1;
					i = i+1;
					treeNodes[i][2]=treeNodes[i-1][0];
					treeNodes[i][3]=treeNodes[i-1][1];
					
					cladesInfo[cladesCount][2] = i-3;
					cladesInfo[cladesCount][3] = 3;
					
					/* copy all clades - all existing clades will have their tree positions shifted by 2 to the left, except
					   for the "global clade", i.e the entire tree.*/
					   
					cladesInfo[0][2] = cladesInfo[0][0];
					cladesInfo[0][3] = cladesInfo[0][1]+1;
					
					i = 1;
					while (cladesInfo[i][1])
					{
						cladesInfo[i][2] = cladesInfo[i][0]-2;
						cladesInfo[i][3] = cladesInfo[i][1];
						i = i+1;
					}				
				}
				else
				{
					if ((n<ds.species)||(m<ds.species))
					{
						/* one node is a leaf, and the other one is "proper" clade */
						/* arrage so that n is the leaf */
						if (m<n)
						{
							i=m;
							m=n;
							n=i;
						}
						/* again - a few things to do 
						   copy the tree except for the leaf, which is added at the end of the existing clade
						   shift clade info if necessary
						*/
						m = m-ds.species;
						i=0;
						shift = 0;
						j = 0;
						s = cladesInfo[m][0]+cladesInfo[m][1];
						t_count = cladesInfo[m][0];
						/* copying nodes other than n, adding leaf at the end of the clade*/
						while (1)
						{
							p = treeNodes[i][0];
							if (p==n)
							{
								shift = -1;
							}
							else
							{	
								if (i==s)
								/* insert new node and new clade*/
								{
									treeNodes[i-1][2] = n;
									treeNodes[i-1][3] = 2;
									treeNodes[i][2] = cladesCount+ds.species;
									treeNodes[i][3] = 1;
									cladesInfo[cladesCount][2] = cladesInfo[m][0]-1;
									cladesInfo[cladesCount][3] = cladesInfo[m][1]+2;
									shift = 1;
									j = 0;
								}
								else 
								{
									if (i==t_count)
									{
										j=1;
										/* send nodes of clade m one level deeper */
									}
								}
								treeNodes[i+shift][2]=treeNodes[i][0];
								treeNodes[i+shift][3]=treeNodes[i][1]+j;							
							}
							if (treeNodes[i][1]==0) 
							{
								break;
							}
							i = i+1;
						}
						
						treeNodes[i+1][2]=treeNodes[i][0];
						treeNodes[i+1][3]=treeNodes[i][1];							
						
						/* last thing to do:
							copy/update clade info matrix => clades starting before newly added clade will
							get shifted 1 to the left (their starting location)
							clades starting after the enlarged will be shifted to the right 1.
						*/
						s = cladesInfo[m][0]+cladesInfo[m][1];
						t_count = cladesInfo[m][0];
						for (k=1; k<cladesCount; k=k+1)
						{
							if (cladesInfo[k][0]<=t_count)
							{
								cladesInfo[k][2] = cladesInfo[k][0]-1;
							}
							else
							{
								cladesInfo[k][2] = cladesInfo[k][0]+1;
							}
							cladesInfo[k][3] = cladesInfo[k][1];
						}
						
						cladesInfo[0][2]=0;
						cladesInfo[0][3]=cladesInfo[0][1]+1;
						
					}
					else
					/* both nodes are "proper" clades 
					   Arrange that node n is to the left of node m. */
					{
						/* start up the merging procedure */
						/* copy all nodes of the tree up to the beginning of the first clade */
						m = m-ds.species;
						n = n-ds.species;
						if (cladesInfo[m][0]<cladesInfo[n][0])
						{
							i=m;
							m=n;
							n=i;
						}
						for (k=cladesInfo[n][0]-1;k>=0; k=k-1)
						{
							treeNodes[k][2] = treeNodes[k][0];
							treeNodes[k][3] = treeNodes[k][1];
						}
						/* copy the nodes of the first clade sending them one level deeper */

						for (k=cladesInfo[n][0]+cladesInfo[n][1]-1;k>=cladesInfo[n][0]; k=k-1)
						{
							treeNodes[k][2] = treeNodes[k][0];
							treeNodes[k][3] = treeNodes[k][1]+1;
						}
						/* copy the nodes of the 2nd clade sending them one level deeper */
						p = cladesInfo[m][0];
						s = cladesInfo[n][0]+cladesInfo[n][1]+cladesInfo[m][1];
							
						for (k=cladesInfo[n][0]+cladesInfo[n][1];k<s; k=k+1)
						{
							treeNodes[k][2] = treeNodes[p][0];
							treeNodes[k][3] = treeNodes[p][1]+1;
							p = p+1;
						}
						
						/* add new node */
						
						treeNodes[s][2]=ds.species+cladesCount;
						treeNodes[s][3]=1;
						
						k=s+1;
						
						/* copy the nodes between the first and second clades */
						
						p = cladesInfo[m][0];
						for (j=cladesInfo[n][0]+cladesInfo[n][1]; j<p; j=j+1)
						{
							treeNodes[k][2]=treeNodes[j][0];
							treeNodes[k][3]=treeNodes[j][1];
							k=k+1;
						}
											
						/* finally, copy all the nodes after the end of the 2nd clade */
						for (j=cladesInfo[m][0]+cladesInfo[m][1]; treeNodes[j][1]; j=j+1)
						{
							treeNodes[k][2]=treeNodes[j][0];
							treeNodes[k][3]=treeNodes[j][1];
							k=k+1;
						}
						
						k=cladesInfo[0][1];
						/*copy root*/
						treeNodes[k][2]=treeNodes[k-1][0];
						treeNodes[k][3]=treeNodes[k-1][1];
						
						/* create new clade info */
						cladesInfo[cladesCount][2] = cladesInfo[n][0];
						cladesInfo[cladesCount][3] = cladesInfo[n][1]+cladesInfo[m][1]+1;
						
						/* lastly - update info for other clades */
						s = cladesInfo[n][0];
						p = cladesInfo[n][1]+s;
						t_count = cladesInfo[m][0];
						j = cladesInfo[m][1]+1;
						for (k=0; k<cladesCount; k=k+1)
						{
							cladesInfo[k][3] = cladesInfo[k][1];
							l = cladesInfo[k][0];
							if (l<=s)
							{
								cladesInfo[k][2] = l;
							}
							else
							{
								if ((l>=p)&&(l<t_count))
								{
									cladesInfo[k][2] = l+j;
								}
								else
								{
									if (l>=t_count+j-1)
									{
										cladesInfo[k][2] = l+1;
									}
									else
									{
										if (l>=t_count)
										{
											cladesInfo[k][2] = l-t_count+p;
										}
										else
										{
											cladesInfo[k][2] = l;
										}
									}
								}
							}
						}
						cladesInfo[m][2] = p;
						cladesInfo[0][3] = cladesInfo[0][1]+1;
					}
				}
				
				thisTree = TreeMatrix2TreeString (2,0);
				if (verbFlag)
				{
					fprintf (stdout, "\n\tTree ",thisTree);
				}				
				Tree    Inferred_Tree = thisTree;
				LikelihoodFunction lf = (filteredData, Inferred_Tree);
			    if (categoriesUsed) {alpha = 0.5;}
				Optimize (res,lf);
				if (res[1][0]>bestValue)
				{
					bestValue = res[1][0];
					bestTree = thisTree;
					bestNode1 = eligibleNodes[loop1];
					bestNode2 = eligibleNodes[loop2];
					
					for (l=0;treeNodes[l][3];l=l+1)
					{
						bestTreeNodes[l][0]=treeNodes[l][2];
						bestTreeNodes[l][1]=treeNodes[l][3];
					}
					
					bestTreeNodes[l][0]=treeNodes[l][2];
					bestTreeNodes[l][1]=treeNodes[l][3];
									
					for (l=0; l<=cladesCount; l=l+1)
					{
						bestCladesInfo[l][0]=cladesInfo[l][2];
						bestCladesInfo[l][1]=cladesInfo[l][3];
					}
					if (verbFlag)
					{
						mxTreeSpec = {4,1};
						mxTreeSpec [0] = "Inferred_Tree";
						mxTreeSpec [3] = "";
						if (twHasBeenOpened)
						{
							mxTreeSpec [1] = "";
							mxTreeSpec [2] = "";
						}
						else
						{
							mxTreeSpec [1] = "8211";
							mxTreeSpec [2] = "10,40,-10,-175,1";
							twHasBeenOpened = 1;
						}
						OpenWindow (TREEWINDOW, mxTreeSpec );
					}

				}
				if (verbFlag)
				{
					fprintf (stdout, " =>", res[1][0]);
				}
			}

		}
		

		/* update eligible nodes lists */
		shift = 0;
		for (l=0;l<eligibleNodesCount;l=l+1)
		{
			if ((eligibleNodes[l]==bestNode1)||(eligibleNodes[l]==bestNode2))
			{
				shift = shift+1;
			}
			else
			{
				eligibleNodes[l-shift]=eligibleNodes[l];
			}
		}
		eligibleNodesCount = eligibleNodesCount-1;
		eligibleNodes[eligibleNodesCount-1]=cladesCount+ds.species;
		
		for (l=0;bestTreeNodes[l][1];l=l+1)
		{
			treeNodes[l][0]=bestTreeNodes[l][0];
			treeNodes[l][1]=bestTreeNodes[l][1];
		}

		treeNodes[l][0]=bestTreeNodes[l][0];
		treeNodes[l][1]=bestTreeNodes[l][1];
		
		for (l=0; l<=cladesCount; l=l+1)
		{
			cladesInfo[l][0]=bestCladesInfo[l][0];
			cladesInfo[l][1]=bestCladesInfo[l][1];
		}
		
		if (verbFlag)
		{
			if (bestNode1<ds.species)
			{
				GetString (s1,ds,bestNode1);
			}
			else
			{
				s1 = "Node"+Format(bestNode1,0,0);
			}
			if (bestNode2<ds.species)
			{
				GetString (s2,ds,bestNode2);
			}
			else
			{
				s2 = "Node"+Format(bestNode2,0,0);
			}
			fprintf (stdout,"\n\n\t ----> Best tree obtained by grouping ",s1," and ", s2, " is\n", bestTree, " with log likelihood of ", bestValue);
		}
		cladesCount = cladesCount+1;
	}
	
	if (doNNIOption)
	{
		taxonCounter = ds.species;
		/* strip the root off for compatibility */
		for (l=Rows(treeNodes)-1; l>0; l=l-1)
		{
			if (treeNodes[l][0]>=ds.species)
			{
				if (treeNodes[l][0] == ds.species)
				{
					treeNodes[l][0] = 0;
				}
				else
				{
					treeNodes[l][0] = treeNodes[l][0]-1;
				}
			}
		}
		for (l=0; l<Rows(cladesInfo)-1; l=l+1)
		{
			cladesInfo[l][0] = cladesInfo[l+1][0];
			cladesInfo[l][1] = cladesInfo[l+1][1];
		}
		cladesInfo[l][0] = 0;
		cladesInfo[l][1] = 0;
		
		starDecomposition = 1;
		if (doNNIOption==1)
		{
			#include "doNNISwap.bf";
		}
		else
		{
			#include "doSPRSwap.bf";
		}
	}
	return 1;
}

l = InferTreeTopology (1.0);

fprintf (stdout,"\n\n --------------------- RESULTS --------------------- \n\n");
fprintf (stdout,"BestTree =", bestTree);

Tree	Inferred_Tree = bestTree;
LikelihoodFunction lf = (filteredData, Inferred_Tree);
if (categoriesUsed) {alpha = 0.5;}
Optimize (res,lf);

if (verbFlag)
{
	OpenWindow (TREEWINDOW, {{"Inferred_Tree"}} );
}

fprintf (stdout, "\n",lf, "\n\n***********Save this tree to a file (y/n)?");
fscanf  (stdin, "String", resp);

if ((resp!="n")&&(resp!="N"))
{
	SetDialogPrompt ("Write tree string to:");
	fprintf (PROMPT_FOR_FILE,CLEAR_FILE,bestTree,";");
}

saveTreeNodes = {2*(ds.species+1),3};
for (i=2*ds.species+1;i>=0;i=i-1)
{
	saveTreeNodes[i][0] = bestTreeNodes[i][0];
	saveTreeNodes[i][1] = bestTreeNodes[i][1];
}
