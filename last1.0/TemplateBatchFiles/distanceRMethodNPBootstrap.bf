
IS_NPBOOTSTRAP_AVAILABLE = 1;

function BootStrapFunction (bsIterates, tabulatedFileName, parametricOrNot)
{
	ChoiceList (bsTreeFormat,"Bootstrap Summary Tree String Format",1,SKIP_NONE,
				"Proportions","Branch lengths of internal nodes reflect the proportion of times the clade starting at that node was contained in replicate trees. Tip branches have length 1.0.",
				"Counts","Branch lengths of internal nodes reflect the number of times the clade starting at that node was contained in replicate trees. Tip branches have length equal to the total number of replicates.");
	if (bsTreeFormat<0)
	{
		return 0.0;
	}
	DataSetFilter filteredData2 = CreateFilter (ds,1,"","");
	nm = ExtractCladeInfo (treeNodes);
	for (bsCounter=1;bsCounter<=bsIterates;bsCounter=bsCounter+1)
	{
		fprintf (stdout,"\nIteration ",Format(bsCounter,0,0),"/",Format(bsIterates,0,0)," ");

		DataSetFilter filteredData = Bootstrap (filteredData2,1);
		p = InferTreeTopology (0.0);
		nm2 = ExtractCladeInfo (treeNodes);
		p = MatchClades (nm2,nm);
		treeString = TreeMatrix2TreeString (0);
		fprintf (tabulatedFileName, treeString,";\n");

	}
	treeNodes = bestTreeNodes;
	i = 0;
	j = 0;
	n = treeNodes[0][1];
	if (bsTreeFormat<.5)
	{
		while (n)
		{
			if (treeNodes[i][0]>=ds.species)
			{
				treeNodes[i][2] = nm[0][j]/bsIterates;
				j=j+1;
			}
			else
			{
				treeNodes[i][2] = 1.0;
			}
			n = treeNodes[i][1];
			i = i+1;
		}
	}
	else
	{
		while (n)
		{
			if (treeNodes[i][0]>=ds.species)
			{
				treeNodes[i][2] = nm[0][j];
				j=j+1;
			}
			else
			{
				treeNodes[i][2] = bsIterates;
			}
			n = treeNodes[i][1];
			i = i+1;
		}
	}

	treeString = TreeMatrix2TreeString (1);
	UseModel(USE_NO_MODEL);
	Tree Bootstrap_Tree = treeString;
	
	if (bsTreeFormat<.5)
	{
		fprintf (stdout,"\n\n************* Bootstrapping summary tree (with ",Format(bsIterates,4,0)," iterates). \nInternal branch 'lengths' indicate the proportion of replicate trees which contained the clade starting at that node.\n\n", treeString,"\n");		
		fprintf (tabulatedFileName,"\n\nBootstrapping summary tree (with ",Format(bsIterates,4,0)," iterates).  \nInternal branch 'lengths' indicate the proportion of replicate trees which contained the clade starting at that node.\n\n", treeString);		
	}
	else
	{
		fprintf (stdout,"\n\n************* Bootstrapping summary tree (with ",Format(bsIterates,4,0)," iterates). \nInternal branch 'lengths' indicate the number of times the clade starting at that node was contained in replicate trees.\n\n", treeString,"\n");		
		fprintf (tabulatedFileName,"\n\nBootstrapping summary tree (with ",Format(bsIterates,4,0)," iterates).  \nInternal branch 'lengths' indicate the number of times the clade starting at that node was contained in replicate trees.\n\n", treeString);		
	}
	treeNodes = bestTreeNodes;
	treeString = TreeMatrix2TreeString (1);
	Tree Inferred_Tree = treeString;
	return 0;
}

/* ____________________________________________*/

function  ExtractCladeInfo (treeNodeMatrix)
{
	i = 0;
	n = treeNodeMatrix[0][1]; 
	s = 0;
	matrixOfClades = {ds.species+1,ds.species-1};
	while (n)
	{
		if (treeNodeMatrix[i][0]>=ds.species)
		{
			j=i-1;
			m = 1;
			while (treeNodeMatrix[j][1]>n)
			{
				if (treeNodeMatrix[j][0]<ds.species)
				{
					m = m+1;
					matrixOfClades[m][s]=treeNodeMatrix[j][0];
				}
				j = j-1;
				if (j<0)
				{
					break;
				}
			}
			matrixOfClades[1][s]=m-1;
			s = s+1;
		}
		i = i+1;
		n = treeNodeMatrix[i][1];
	}
	return matrixOfClades;
}

/* ____________________________________________*/
/* see if clades in matrix2 exist in matrix1 and if so, increment the counters
   for relevant clades.
*/

function   MatchClades(treeNodeMatrix1,treeNodeMatrix2)
{
	n = Columns(treeNodeMatrix2);
	for (i=0;i<n;i=i+1)
	{
		m = treeNodeMatrix1[1][i];
		if (m==0)
		{
			break;
		}
		for (j=0;j<n;j=j+1)
		{
			if (treeNodeMatrix2[1][j]==m)
			/*might work*/
			{
				for (s=m+1; s>1 ; s=s-1)
				/* for each node in the clade */
				{
					k = treeNodeMatrix1[s][i];
					for (p=m+1; p>1; p=p-1)
					{
						if (treeNodeMatrix2[p][j]==k)
						{
							break;
						}
					}
					if (p==1)
					{
						break;
					}
				}
				if (s==1)
				{
					/*fprintf (stdout, "\n  ",j);*/
					treeNodeMatrix2[0][j] = treeNodeMatrix2[0][j]+1;
					break;
				}
			}
		}
		if ((j==n)&&(m>2))
		/* clade not found, try complement */
		{
			t = ds.species-m;
			if (m==0)
			{
				break;
			}
			for (j=0;j<n;j=j+1)
			{
				if (treeNodeMatrix2[1][j]==t)
				/*might work*/
				{
					for (s=m+1; s>1 ; s=s-1)
					/* for each node in the clade */
					{
						k = treeNodeMatrix1[s][i];
						for (p=t+1; p>1; p=p-1)
						{
							if (treeNodeMatrix2[p][j]==k)
							{
								break;
							}
						}
						if (p>1)
						{
							break;
						}
					}
					if (s==1)
					{
						/*fprintf (stdout, "\nC ",j,"(",m,")");*/
						treeNodeMatrix2[0][j] = treeNodeMatrix2[0][j]+1;
						break;
					}
				}
			}		
		}
	}
	return 1.0;
}

/* ____________________________________________*/

function TreeMatrix2TreeString (doLengths)
{
	treeString = "";
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
		if (doLengths>.5)
		{
			treeString = treeString+":"+treeNodes[k][2];
		}
		k=k+1;
		p=m;
		n=treeNodes[k][0];
		m=treeNodes[k][1];
	}

	for (j=m;j<p;j=j+1)
	{
		treeString = treeString+")";
	}
	
	return treeString;
}

/* ____________________________________________*/

function	placeParentNodes (nodeID)
{
	nodeDepth = 0;
	
	saveNodeID = nodeID;
	parentID   = parentInfo[nodeID][0];
	
	saveNodeID = nodeID;
	parentID   = parentInfo[nodeID][0];
	
	layoutOffset = useColumn;

	while (parentID>=0)
	{
		n = distanceMatrix[parentID-ds.species][1];
		if (n >= 0)
		/* it has been laid out */
		{	
			layoutOffset = n+distanceMatrix[parentID-ds.species][0];
			break;
		}
		parentID  = parentInfo[parentID][0];
	}

	parentID   = parentInfo[nodeID][0];

	while (parentID>=0)
	{
		n = parentID-ds.species;
		m = cladesInfo[nodeID];
		if (distanceMatrix[n][1] < 0)
		/* this node hasn't been laid out yet */
		{
			if (parentInfo[parentID][0]>=0)
			{
				distanceMatrix[n][1] = layoutOffset; /* where the layout for the clade begins */
				distanceMatrix[n][0] = m; /* offset for that layout */
			}
						
			m = layoutOffset + m - 1;
			
			treeNodes[m][0] = nodeID;
			treeNodes[m][2] = parentInfo[nodeID][1];
			
			columnIndex[nodeID] = m;
		}
		else
		/* it has been laid out */
		{	
			m = distanceMatrix[n][0]+distanceMatrix[n][1]+ m - 1;
			
			treeNodes[m][0] = nodeID;
			treeNodes[m][2] = parentInfo[nodeID][1];
			
			distanceMatrix[n][0] = m + cladesInfo[nodeID];
			columnIndex[nodeID] = m;
			nodeDepth = nodeDepth+1;
			
			break;
		}
		nodeDepth = nodeDepth+1;
		nodeID    = parentID;
		parentID  = parentInfo[nodeID][0];
	}
	
	/* update levels of nodes */
	
	if (parentID<0)
	{
		nodeID   = saveNodeID;
		parentID = parentInfo[nodeID][0];
		
		while (parentID>=0)
		{
			m = columnIndex[nodeID];
			treeNodes[m][1] = nodeDepth;
			nodeDepth = nodeDepth-1;
			saveNodeID = nodeID;
			nodeID     = parentID;
			parentID   = parentInfo[nodeID][0];
		}
		
		useColumn = useColumn+cladesInfo[saveNodeID];
	}
	else
	{
		m = columnIndex[parentID]; 
		
		n = treeNodes[m][1];/* depth of the parent */
		
		nodeID   = saveNodeID;

		while (nodeDepth >= 0)
		{
			m = columnIndex[nodeID];
			treeNodes[m][1] = nodeDepth+n;
			
			nodeDepth = nodeDepth-1;
			nodeID  = parentInfo[nodeID][0];
		}
	}
	return 0;
}

