ChoiceList (distanceChoice, "Distance Computation",1,SKIP_NONE,
			"Distance formulae","Use one of the predefined distance measures based on data comparisons. Fast.",
			"Full likelihood","Estimate distances using pairwise MLE. More choices but slow.");
			
if (distanceChoice < 0)
	return 0;

ChoiceList (dataType,"Data type",1,SKIP_NONE,"Nucleotide/Protein","Nucleotide or amino-acid (protein).",
				     "Codon","Codon (several available genetic codes).");
				     
A_DISTANCE_METHOD = 1;

#include "distanceRMethodNPBootstrap.bf";

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

if (distanceChoice > 0)
{
	SelectTemplateModel(filteredData);
}

ChoiceList (methodIndex,"Clustering Method",1,SKIP_NONE,
			"UPGMA","Unweighted Pair Group Method using Arithmetic Averages. (new distance)=1/2(sum of old distances).",
			"WPGMA","Weighted Pair Group Method using Arithmetic Averages. (new distance)=((Taxa in clade 1)dist1+(Taxa in clade 2)dist2)/(Taxa in clade 1+Taxa in clade 2).",
			"Complete Linkage","Complete Linkage Method. (new distance)=max(old distances).",
			"Single Linkage","Single Linkage Method. (new distance)=min(old distances).");
			
if (methodIndex<0)
{
	return;
}
function InferTreeTopology(verbFlag)
{
	columnMin 		= {ds.species, 1};
	indexMin  		= {ds.species, 1};
	distanceMatrix  = {ds.species,ds.species};
	useColumn 		= {ds.species, 1};
	columnIndex 	= {ds.species, 1};
	parentInfo		= {2*(ds.species+1),3};
	netDivergence 	= {ds.species, 1};
	cladesInfo 		= {2*(ds.species+1),1};
	taxaCount		= {ds.species-1, 1};

	/* first column has the node indices in pre-order */
	/* second column contains the depths of the nodes (root has depth 0) */

	/* first column has the clades size (in nodes)
	   second column has the clades position in the tree "stack" */
	   
	/* used to count how many taxa are in a given clade */
	   
	treeLength = 0;

	MESSAGE_LOGGING = 0;

	if (distanceChoice)
	{
		if (verbFlag)
		{
			fprintf (stdout,"\nHYPHY is computing pairwise maximum likelihood distance estimates. A total of ", Format(ds.species*(ds.species-1)/2,0,0),
						   " estimations will be performed.\n");
		}
		_pddVF = 1;
		ExecuteAFile ("pairwiseDistanceEstimator.ibf");

		for (i = 0; i<ds.species; i=i+1)
		{
			min = 1e10;
			minIndex = -1;
			for (j = 0; j<i; j = j+1)
			{
				k = distanceMatrix[j][i];			
				if (k<min)
				{
					min 	 = k;
					minIndex = j;
				}
			}
			columnMin   [i] = min;
			indexMin    [i] = minIndex;
			useColumn   [i] = 1;
			columnIndex [i] = i;
			cladesInfo  [i] = 1;
		}
	}
	else
	{
		#include "chooseDistanceFormula.def";
		
		dummy = InitializeDistances (0);
		
		for (i = 0; i<ds.species; i=i+1)
		{
			min = 1e10;
			minIndex = -1;
			for (j = 0; j<i; j = j+1)
			{
				k = ComputeDistanceFormula (i,j);
				distanceMatrix[j][i] = k;			
				if (k<min)
				{
					min 	 = k;
					minIndex = j;
				}
			}
			columnMin	[i] = min;
			indexMin 	[i] = minIndex;
			useColumn	[i] = 1;
			columnIndex	[i] = i;
			cladesInfo  [i] = 1;
		}
	}

	MESSAGE_LOGGING = 1;

	cladesMade = 1;

	while (cladesMade < ds.species)
	{
		/*fprintf (stdout,"\nStage ",cladesMade,"\n"); */
		min 		= 1e9;
		minIndex 	= -1;
		for (i=0;i<ds.species;i=i+1)
		{
			if (useColumn[i])
			{
				if (columnMin[i]<min)
				{
					min 	 = columnMin[i];
					minIndex = i;
				}
			}
		}

		d 		 = min;
		j 		 = minIndex;
		minIndex = indexMin[minIndex];
		
		/* update tree description */
		
		m = columnIndex[j];
		n = columnIndex[minIndex];

		k = ds.species+cladesMade-1;		

		parentInfo[n][0] = k; /* parent */
		parentInfo[m][0] = k; /* parent */
		cladesInfo[k]    = cladesInfo[n]+cladesInfo[m]+1;
		
		p = parentInfo[n][2];
		l = parentInfo[m][2];
		
		d2 = (d+p-l)/2;
		
		if (d2>d)
		{
			d2 = d;
		}
		if (d2<0)
		{
			d2 = 0;
		}
		
		/*fprintf (stdout, "Cluster ", n, ":", d-d2, " and ", m, ":", d2, "\n");*/
				
		parentInfo[m][1] = d2;
		parentInfo[n][1] = d-d2;
				
		parentInfo[k][2] = parentInfo[n][1]+parentInfo[n][2];
		
		k = cladesMade - 1;
		
		if (n<ds.species)
		{
			taxaCount [k] = taxaCount [k]+1;
		}
		else
		{
			taxaCount [k] = taxaCount [k]+taxaCount[n-ds.species];		
		}
		
		if (m<ds.species)
		{
			taxaCount [k] = taxaCount [k]+1;
		}
		else
		{
			taxaCount [k] = taxaCount [k]+taxaCount[m-ds.species];		
		}
		

		useColumn[j] = 0;

		if (methodIndex == 0)
		{
			/*UPGMA*/
			for (k=0;k<minIndex; k=k+1)
			{
				if (useColumn[k])
				{
					distanceMatrix[k][minIndex] = (distanceMatrix[k][minIndex]+distanceMatrix[k][j])/2;
				}
			}
						
			for (k=minIndex+1;k<j;k=k+1)
			{
				if (useColumn[k])
				{
					distanceMatrix[minIndex][k] = (distanceMatrix[minIndex][k] + distanceMatrix[k][j])/2;
				}
			}
			
			for (k=j+1;k<ds.species;k=k+1)
			{
				if (useColumn[k])
				{
					distanceMatrix[minIndex][k] = (distanceMatrix[minIndex][k] + distanceMatrix[j][k])/2;
				}
			}	
		}
		else
		{
			if (methodIndex == 1)
			{
				/*WPGMA*/
				m = columnIndex[j];
				n = columnIndex[minIndex];
				if (m>=ds.species)
				{
					p = taxaCount[m-ds.species];
				}
				else
				{
					p = 1;
				}
				if (n>=ds.species)
				{
					l = taxaCount[n-ds.species];
				}
				else
				{
					l = 1;
				}
				m = l+p;
				for (k=0;k<minIndex; k=k+1)
				{
					if (useColumn[k])
					{
						distanceMatrix[k][minIndex] = (distanceMatrix[k][minIndex]*l+distanceMatrix[k][j]*p)/m;
					}
				}
				
				for (k=minIndex+1;k<j;k=k+1)
				{
					if (useColumn[k])
					{
						distanceMatrix[minIndex][k] = (distanceMatrix[minIndex][k]*l+distanceMatrix[k][j]*p)/m;
					}
				}
				for (k=j+1;k<ds.species;k=k+1)
				{
					if (useColumn[k])
					{
						distanceMatrix[minIndex][k] = (distanceMatrix[minIndex][k]*l+distanceMatrix[j][k]*p)/m;
					}
				}	
			}	
			else
			{
				if (methodIndex == 2)
				{
					/*complete linkage*/
					for (k=0;k<minIndex; k=k+1)
					{
						if (useColumn[k])
						{
							if (distanceMatrix[k][j]>distanceMatrix[k][minIndex])
							{
								distanceMatrix[k][minIndex] = distanceMatrix[k][j];
							}
						}
					}

					for (k=minIndex+1;k<j;k=k+1)
					{
						if (useColumn[k])
						{
							if (distanceMatrix[minIndex][k]<distanceMatrix[k][j])
							{
								distanceMatrix[minIndex][k]=distanceMatrix[k][j];
							}
						}
					}
					
					for (k=j+1;k<ds.species;k=k+1)
					{
						if (useColumn[k])
						{
							if (distanceMatrix[minIndex][k]<distanceMatrix[j][k])
							{
								distanceMatrix[minIndex][k]=distanceMatrix[j][k];
							}
						}
					}	
				}
				else
				{
					/* single linkage */
					for (k=0;k<minIndex; k=k+1)
					{
						if (useColumn[k])
						{
							if (distanceMatrix[k][j]<distanceMatrix[k][minIndex])
							{
								distanceMatrix[k][minIndex] = distanceMatrix[k][j];
							}
						}
					}
					for (k=minIndex+1;k<j;k=k+1)
					{
						if (useColumn[k])
						{
							if (distanceMatrix[minIndex][k]>distanceMatrix[k][j])
							{
								distanceMatrix[minIndex][k]=distanceMatrix[k][j];
							}
						}
					}
					for (k=j+1;k<ds.species;k=k+1)
					{
						if (useColumn[k])
						{
							if (distanceMatrix[minIndex][k]>distanceMatrix[j][k])
							{
								distanceMatrix[minIndex][k]=distanceMatrix[j][k];
							}
						}
					}				
				}	
			}
		}
		
		columnIndex[minIndex] = ds.species+cladesMade-1;
		
		for (i=0; i<j; i=i+1)
		{
			distanceMatrix[i][j]=0;
		}

		for (i=j+1; i<ds.species; i=i+1)
		{
			distanceMatrix[j][i]=0;
		}
		
		cladesMade = cladesMade+1;
		
		for (i = 1; i<ds.species; i=i+1)
		{
			if (useColumn[i])
			{
				min = distanceMatrix[0][i];
				minIndex = 0;
				
				for (j = 0; j<i; j = j+1)
				{
					if (useColumn[j])
					{
						if (distanceMatrix[j][i]<min)
						{
							min = distanceMatrix[j][i];
							minIndex = j;
						}
					}
				}
				columnMin[i] = min;
				indexMin [i] = minIndex;
			}
		}
	}
	
	taxaCount		= 0;
	
	columnIndex     = {2*(ds.species+1),1};
	distanceMatrix  = {ds.species-1,2};

	treeNodes 		= {2*(ds.species+1),3};
	useColumn		= 0;
	columnMin		= 0;
	indexMin		= 0;

	   
	parentInfo	[ds.species+cladesMade-2][0] = -1; 

	for (j=0; j<ds.species-1; j=j+1)
	{
		distanceMatrix[j][1] = -1;
	}

	for (j=0; j<ds.species; j=j+1)
	{
		i = placeParentNodes (j);
	}

	for (j=0; j<ds.species-1; j=j+1)
	{
		distanceMatrix[j][0] = cladesInfo[ds.species+j];
	}

	distanceMatrix[cladesMade-2][1] = 0;
	treeNodes[useColumn][0] = cladesMade-2+ds.species;
	treeNodes[useColumn][1] = 0;

	cladesInfo 	   = distanceMatrix;  

	distanceMatrix = 0;
	parentInfo 	   = 0;
	columnIndex	   = 0;
	
	return 1;
}

DISTANCE_PROMPTS  = 1;
p = InferTreeTopology(1.0);
DISTANCE_PROMPTS  = 0;


/* now with the treeNodes matrix ready we can convert it into a Newick string */

treeString = TreeMatrix2TreeString (1);

fprintf (stdout,"\n\n --------------------- INFERRED TREE --------------------- \n\n", treeString);

fprintf (stdout, "\n\n***********Save this tree to a file (y/n)?");

fscanf  (stdin, "String", resp);

if ((resp!="n")&&(resp!="N"))
{
	SetDialogPrompt ("Write tree string to:");
	fprintf (PROMPT_FOR_FILE,CLEAR_FILE,treeString,";");
}

UseModel (USE_NO_MODEL);

Tree Inferred_Tree = treeString; 

bestTreeNodes = treeNodes;
/* now with the treeNodes matrix ready we can convert it into a Newick string */

/*treeString = "";
p = 0;
k = 0;
m = treeNodes[0][1];
n = treeNodes[0][0];

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
	treeString = treeString+":"+treeNodes[k][2];
	k=k+1;
	p=m;
	n=treeNodes[k][0];
	m=treeNodes[k][1];
}

for (j=m;j<p;j=j+1)
{
	treeString = treeString+")";
}

fprintf (stdout,"\n\n --------------------- INFERRED TREE --------------------- \n\n", treeString);

fprintf (stdout, "\n\n***********Save this tree to a file (y/n)?");

fscanf  (stdin, "String", resp);

if ((resp!="n")&&(resp!="N"))
{
	SetDialogPrompt ("Write tree string to:");
	fprintf (PROMPT_FOR_FILE,CLEAR_FILE,treeString,";");
}*/
