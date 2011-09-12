/* NNI Branch Swapping Functions */

function	copyTreeStructure	 (terminationLevel,levelModifier)
{
	while (1)
	{
		scbj = treeNodes[scbi][0];
		splitTreeInfo  [scbk][scbp] = treeNodes[scbi][1]+levelModifier;
		
		if (scbj<spCount)
		{
			splitTreeInfo [scbk][scbn] = scbj;
		}
		else
		{
			splitTreeInfo   [scbk][scbn] = spCount+scbm;
			splitCladesInfo [scbm][scbn] = cladesInfo [scbj-spCount][0]-scbo;
			splitCladesInfo [scbm][scbp] = cladesInfo [scbj-spCount][1];
			scbm = scbm+1;
		}
		if (treeNodes[scbi][1]==terminationLevel)
		{
			break;
		}
		scbk = scbk+1;		
		scbi = scbi+1;
	}
	return 0;
}

function	splitCladesOnIBranch (spCount, iNodeIndex)
{
	splitTreeInfo      = {Rows(treeNodes),8};
	splitCladesInfo    = {Rows(cladesInfo),8};
	/* do the two easy ones 1st */
	iNodeIndex = iNodeIndex - spCount;
	treeInfoSplitIndex = cladesInfo[iNodeIndex][1]+cladesInfo[iNodeIndex][0]-1;
	treeInfoSplitDepth = treeNodes[treeInfoSplitIndex][1];
	scbk 		  = 0;
	scbm 		  = 0;
	scbn		  = 0;
	scbp		  = 1;
	scbo		  = cladesInfo[iNodeIndex][0];
	newINodeIndex = 0;
	for (scbi=cladesInfo[iNodeIndex][0];scbi<treeInfoSplitIndex;scbi=scbi+1)
	{
		scbj = treeNodes[scbi][0];
		splitTreeInfo  [scbk][scbp] = treeNodes[scbi][1]-treeInfoSplitDepth;
		
		if (scbj<spCount)
		{
			splitTreeInfo [scbk][scbn] = scbj;
		}
		else
		{
			splitTreeInfo   [scbk][scbn] = spCount+scbm;
			splitCladesInfo [scbm][scbn] = cladesInfo [scbj-spCount][0]-scbo;
			splitCladesInfo [scbm][scbp] = cladesInfo [scbj-spCount][1];
			scbm = scbm+1;
		}
		if ((splitTreeInfo[scbk][1]==1)&&(scbn==0))
		{
			scbo = scbi+1;
			scbn = 2;
			scbp = 3;
			scbk = 0;
			scbm = 0;
			continue;
		}
		scbk = scbk+1;
	}
	/* build the third clade by taking the 2nd child of the parent */
	scbn = 4;
	scbp = 5;
	scbm = 0;
	scbk = 0;
	scbi = scbi+1;
	if (treeNodes[scbi][1] < treeInfoSplitDepth)
	/* was the right child */
	{
		scbi = cladesInfo[iNodeIndex][0]-1;
		scbj = treeNodes[scbi][0];
		if (scbj>=spCount)
		{
			scbi = cladesInfo[scbj-spCount][0];
		}
	}
	scbo = scbi;
	scbn = copyTreeStructure (treeInfoSplitDepth,1-treeInfoSplitDepth);
	scbn = 6;
	scbp = 7;
	if (treeInfoSplitDepth == 1)
	/* root level split */
	{
		scbi = scbi+1;
		if ((treeNodes[scbi][1] == 0)||(scbi==cladesInfo[iNodeIndex][0]))
		{
			scbi = 0;
		}
		scbo = scbi;
		scbm = 0;
		scbk = 0;
		scbj = copyTreeStructure (1,0);
	}
	else
	/* not a root level split */
	{
		scbj = treeInfoSplitDepth-1;
		pathToTheRoot = {scbj,1};
		pathToTheRoot [0] = L1ParentIndex+1;
		/* find a level 1 parent of the split point */
		for (L1ParentIndex = treeInfoSplitIndex;;L1ParentIndex=L1ParentIndex+1)
		{
			if (treeNodes[L1ParentIndex][1] == scbj)
			{
				pathToTheRoot[scbj-1] = L1ParentIndex;
				if (scbj == 1)
				{
					break;
				} 
				scbj = scbj - 1;
			}
		}
		/* handle the top level separately */
		scbk = 0;
		scbm = 0;
		
		scbi = 0;
		scbo = 0;
		scbq = pathToTheRoot[0];
		scbq = treeNodes[scbq][0]-spCount;
		scbq = cladesInfo[scbq][0];
		for (scbr = 0; scbr < 3; scbr = scbr+1)
		{
			if (scbi == scbq)
			{
				scbo = pathToTheRoot[0]-scbi+1;
				scbi = pathToTheRoot[0];
			}
			else
			{
				dummy = copyTreeStructure (1,treeInfoSplitDepth-1);
				scbk = scbk+1;
			}
			scbi = scbi+1;
		}
		
		/* add a new internal node to join the other two */
		splitTreeInfo   [scbk][scbn] = spCount+scbm;
		splitTreeInfo   [scbk][scbp] = treeInfoSplitDepth-1;
		splitCladesInfo [scbm][scbn] = 0;
		scbk = scbk+1;
		splitCladesInfo [scbm][scbp] = scbk;
		scbm = scbm+1;

		for (scbq = 1; scbq<Rows(pathToTheRoot); scbq = scbq+1)
		{
			scbj = pathToTheRoot [scbq];
			scbo = pathToTheRoot [scbq-1];
			if (scbo==scbj+1)
			/* came from the right child */
			{
				scbi = treeNodes  [scbo][0]-spCount;
				scbi = cladesInfo[scbi][0];
			}
			/* came from the left child */
			else
			{
				scbi = scbj+1;
			}
			scbs = scbi;
			scbo = scbi-scbk;
			dummy = copyTreeStructure (scbq+1,treeInfoSplitDepth-2*scbq-1);			
			scbk = scbk+1;
			splitTreeInfo   [scbk][scbn] = spCount+scbm;
			splitTreeInfo   [scbk][scbp] = Rows(pathToTheRoot)-scbq;
			splitCladesInfo [scbm][scbn] = 0;
			scbk = scbk+1;
			splitCladesInfo [scbm][scbp] = scbk;
			scbm = scbm+1;
		} 		
	}
	return 0;
}

function	doNNISwap (b2group,spCount)
{
	if (b2group==3)
	{
		branchOrderMatrix = {{0,3,1,2}};	
	}
	else
	{
		branchOrderMatrix = {{0,2,1,3}};
	}
	
	scbo = 0;
	scbn = 0;
	scbp = 1;
	scbj = 0;
	scbk = 0;
	
	for (scbi=0; scbi<4; scbi = scbi+1)
	{
		scbt = branchOrderMatrix[scbi];
		scbn = scbt*2;
		scbp = scbn+1;
		scbq = 0;
		while (splitTreeInfo[scbq][scbp])
		{
			if (splitTreeInfo[scbq][scbn]<spCount)
			{
				treeNodes [scbj][2] = splitTreeInfo[scbq][scbn];
			}
			else
			{
				scbm = splitTreeInfo[scbq][scbn]-spCount;
				treeNodes	[scbj][2] = spCount+scbk;
				cladesInfo	[scbk][2] = splitCladesInfo[scbm][scbn] + scbo;
				cladesInfo  [scbk][3] = splitCladesInfo[scbm][scbp];
				scbk=scbk+1;
			}
			if (scbi<2)
			{
				treeNodes[scbj][3] = splitTreeInfo[scbq][scbp]+1;
			}
			else
			{
				treeNodes[scbj][3] = splitTreeInfo[scbq][scbp];
			}
			scbj = scbj+1;
			scbq = scbq+1;
		}
		scbo = scbo + scbq;
		if (scbi==1)
		{
			scbo = scbo+1;
			treeNodes	[scbj][2] = spCount+scbk;
			treeNodes	[scbj][3] = 1;
			cladesInfo	[scbk][2] = 0;
			cladesInfo  [scbk][3] = scbo;
			scbj = scbj+1;
			scbk = scbk+1;			
		}
	}
	treeNodes[scbj][2] = 0;
	treeNodes[scbj][3] = 0;
	return 0;
}

function	doSubtreeCopy (startAt, upTo, depthCorrection)
{
	if (upTo>=(-1))
	{
		for   (scbi=startAt; scbi<=upTo; scbi=scbi+1)
		{
			if (splitTreeInfo[scbi][scbn]<spCount)
			{
				treeNodes [scbj][2] = splitTreeInfo[scbi][scbn];
			}
			else
			{
				scbm = splitTreeInfo[scbi][scbn]-spCount;
				treeNodes	[scbj][2] = spCount+scbk;
				cladesInfo	[scbk][2] = splitCladesInfo[scbm][scbn] + scbo;
				cladesInfo  [scbk][3] = splitCladesInfo[scbm][scbp];
				scbk=scbk+1;
			}
			treeNodes[scbj][3] = splitTreeInfo[scbi][scbp]+depthCorrection;
			scbj = scbj+1;
		}
		return 0;
	}
	else
	{
		scbi = startAt;
		while (splitTreeInfo[scbi][scbp])
		{
			if (splitTreeInfo[scbi][scbn]<spCount)
			{
				treeNodes [scbj][2] = splitTreeInfo[scbi][scbn];
			}
			else
			{
				scbm = splitTreeInfo[scbi][scbn]-spCount;
				treeNodes	[scbj][2] = spCount+scbk;
				cladesInfo	[scbk][2] = splitCladesInfo[scbm][scbn] + scbo;
				cladesInfo  [scbk][3] = splitCladesInfo[scbm][scbp];
				scbk=scbk+1;
			}
			treeNodes[scbj][3] = splitTreeInfo[scbi][scbp]+depthCorrection;
			scbj = scbj+1;		
			scbi = scbi+1;
		}
		return scbi-startAt;
	}
	
}


function	doSPRSwap (bFgroup1, bFgroup2, b2group,where2group,spCount)
{	
	scbn = b2group*2;
	scbp = scbn+1;
		
	scbo = 0;
	scbj = 0;
	scbk = 0;
	scbd = splitTreeInfo[where2group][scbp];
	
	if (splitTreeInfo[where2group][scbn]<spCount)
	/* grafting onto a leaf */
	{
		dummy = doSubtreeCopy (0,where2group-1,0);
		scbr  = where2group;
		dummy = doSubtreeCopy (where2group,where2group,1);
		scbn  = bFgroup1*2;
		scbp  = scbn+1;
		scbo  = where2group+1;
		scbq  = scbj;
		dummy = doSubtreeCopy (0,-2,scbd+1);
		scbo  = scbo+dummy;
		scbn  = bFgroup2*2;
		scbp  = scbn+1;
		dummy = doSubtreeCopy (0,-2,scbd+1);
		scbo  = scbo+dummy;	
	}
	else
	/* grafting onto an internal branch */
	{
		scbq  = splitTreeInfo  [where2group][scbn]-spCount;
		scbq  = splitCladesInfo[scbq][scbn];
		dummy = doSubtreeCopy (0,scbq-1,0);
		scbr  = scbj;
		dummy = doSubtreeCopy (scbq,where2group,1);
		scbo  = scbj;
		scbn  = bFgroup1*2;
		scbp  = scbn+1;
		scbo  = where2group+1;
		scbq  = scbj;
		dummy = doSubtreeCopy (0,-2,scbd+1);
		scbo  = scbo+dummy;
		scbn  = bFgroup2*2;
		scbp  = scbn+1;
		dummy = doSubtreeCopy (0,-2,scbd+1);
		scbo  = scbo+dummy;		
	}
	
	scbv = scbo-scbq+2;
	
	/* create 2 new internal nodes */
		
	for (scbi=1; scbi>=0; scbi = scbi-1)
	{
		scbo = scbo+1;
		treeNodes	[scbj][2] = spCount+scbk;
		treeNodes	[scbj][3] = scbd+scbi;
		cladesInfo	[scbk][2] = scbq;
		cladesInfo  [scbk][3] = scbj-scbq+1;
		scbj = scbj+1;
		scbk = scbk+1;		
		scbq = scbr;
	}	
	
	
	scbn  = b2group*2;
	scbp  = scbn+1;	
	scbo  = scbo-where2group-1;
	dummy = doSubtreeCopy (where2group+1,-2,0);
	
	scbr = scbd-1;
	/* update the number of nodes in affected clades */
	for (scbi=where2group+1; scbi<scbj; scbi=scbi+1)
	{
		if ((treeNodes[scbi][3]==scbr)&&(treeNodes[scbi][2]>spCount))
		{
			scbq = treeNodes[scbi][2]-spCount;
			cladesInfo[scbq][3] = cladesInfo[scbq][3]+scbv;
			cladesInfo[scbq][2] = cladesInfo[scbq][2]-scbo;
			scbr = scbr-1;
		}
	}

	scbj  = scbj-1;
	scbk  = scbk-1;
	scbo  = scbo+where2group+dummy;
	for (scbi=0; scbi<scbj; scbi=scbi+1)
	{
		treeNodes[scbi][3] = treeNodes[scbi][3]-1;
	}
	
	
	for (scbi=0; scbi<4; scbi=scbi+1)
	{
		if ((scbi==bFgroup1)||(bFgroup2==scbi)||(b2group==scbi))
		{
			continue;
		}
		break;
	}
	
	scbn  = scbi*2;
	scbp  = scbn+1;
	dummy = doSubtreeCopy (0,-2,0);

	treeNodes[scbj][2] = 0;
	treeNodes[scbj][3] = 0;
	cladesInfo[scbk][2] = 0;
	cladesInfo[scbk][3] = 0;
	return 0;
}
