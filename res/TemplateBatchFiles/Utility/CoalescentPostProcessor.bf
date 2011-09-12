VERBOSITY_LEVEL = -1;
incFileName = HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "Utility" + DIRECTORY_SEPARATOR + "NJ.bf";
ExecuteCommands ("#include \""+incFileName + "\";");

/* ________________________________________________________________________________________________*/


function processAPartition (dataFilePath, partitionMatrix)
{
	if (Abs (dataFilePath))
	{
		DataSet 	ds = ReadDataFile (dataFilePath);
	}
	readPCount 	   = Columns (partitionMatrix);
	/*fprintf (stdout, "\nLoaded ", readPCount , " partitions\n");*/
	
	outAVL = {};
	
	partitionStrings = {};
	treeStrings		= {};

	lastP = 0;
	for (h=0; h<readPCount; h=h+1)
	{
		pString 					= ""+lastP+"-"+partitionMatrix[h];
		lastP 						= partitionMatrix[h]+1;
		partitionStrings [h] 		= pString;
		DataSetFilter filteredData 	= CreateFilter (ds,1,pString);
		treeStrings[h] 				= InferTreeTopology(filteredData);
	}
	
	pString = ""+lastP+"-"+(ds.sites-1);
	partitionStrings [h] = pString;
	DataSetFilter filteredData = CreateFilter (ds,1,pString);
	treeStrings[h] = InferTreeTopology(filteredData);

	resp = 4;

	global betaP = 1;
	global betaQ = 1;
	betaP:>0.05;betaP:<85;
	betaQ:>0.05;betaQ:<85;

	category pc = (resp-1, EQUAL, MEAN, 
					_x_^(betaP-1)*(1-_x_)^(betaQ-1)/Beta(betaP,betaQ), /* density */
					IBeta(_x_,betaP,betaQ), /*CDF*/
					0, 				   /*left bound*/
					1, 			   /*right bound*/
					IBeta(_x_,betaP+1,betaQ)*betaP/(betaP+betaQ)
				   );


	global alpha = .5;
	alpha:>0.01;alpha:<100;
	category c = (resp, pc, MEAN, 
					GammaDist(_x_,alpha,alpha), 
					CGammaDist(_x_,alpha,alpha), 
					0 , 
			  	    1e25,
			  	    CGammaDist(_x_,alpha+1,alpha)
			  	 );

	global		 AC = 1;
	global		 AT = 1;
	global		 CG = 1;
	global		 CT = 1;
	global		 GT = 1;

	GTR_Matrix = {{*,AC*t*c,t*c,AT*t*c}
				 {AC*t*c,*,CG*t*c,CT*t*c}
				 {t*c,CG*t*c,*,GT*t*c}
				 {AT*t*c,CT*t*c,GT*t*c,*}};

	DataSetFilter	filteredData = CreateFilter (ds,1);

	bppMap = {};
	for (h=0; h<filteredData.sites; h=h+1)
	{
		filterString = "" + h;
		DataSetFilter siteFilter = CreateFilter (filteredData,1,filterString);

		HarvestFrequencies (f1, siteFilter, 1, 1, 0);
		m1 = 0;
		for (mpiNode=0; (mpiNode < 4) && (m1<=1) ; mpiNode=mpiNode+1)
		{
			if (f1[mpiNode]>0)
			{
				m1=m1+1;
			}
		}	
		if (m1>1)
		{
			bppMap[Abs(bppMap)] = h;
		}
	}

	/*fprintf (stdout, "\nSequences :", filteredData.species,
					 "\nSites     :", filteredData.sites,
					 "\nVariable  :", Abs(bppMap), "\n"); */
					 
	outAVL ["VS"] = Abs(bppMap);

	HarvestFrequencies (baseFreqs,ds,1,1,1);

	/*fprintf (stdout, "\n\nf(A) = ", baseFreqs[0],
					 "\nf(C) = ", baseFreqs[1],
					 "\nf(G) = ", baseFreqs[2],
					 "\nf(T) = ", baseFreqs[3],"\n");*/
					
	Model GTR_Model =  (GTR_Matrix, baseFreqs, 1);

	treeString 	= InferTreeTopology   (0);
	Tree givenTree = treeString;

	LikelihoodFunction singlePart = (filteredData, givenTree);
	Optimize (res, singlePart);

	baseParams = res[1][2];

	if (baseParams>0)
	{

		ConstraintString = "";
		ConstraintString*256;
		for (h=0; h<baseParams; h=h+1)
		{
			GetString (v,singlePart,h);
			ConstraintString * (v+":="+v+"__;\n");
		}
		ConstraintString*0;
		ExecuteCommands (ConstraintString);
	}

	nullCAIC = -2(res[1][0]-res[1][1]*filteredData.sites/(filteredData.sites-res[1][1]-1));

	/*fprintf  (stdout, "\n\nc-AIC = ", nullCAIC, "\n", singlePart);*/
	
	outAVL["AC"] = AC;
	outAVL["AT"] = AT;
	outAVL["CG"] = CG;
	outAVL["CT"] = CT;
	outAVL["GT"] = GT;
	outAVL["alpha"] = alpha;

	m1 = computeMeanDivergence ("givenTree");

	/*fprintf  (stdout, "\nMean divergence : ", m1*100, "%\n");*/
	outAVL ["Divergence"] = m1*100;

	USE_DISTANCES = 0;

	if (readPCount == 0)
	{
		return outAVL;
	}
	
	lfDef = "";
	lfDef * 128;
	lfDef  * "LikelihoodFunction multiPart  = (";

	for (pccounter = 0; pccounter <= readPCount; pccounter = pccounter + 1)
	{
		ExecuteCommands ("DataSetFilter part_" + pccounter + " = CreateFilter (ds, 1, \"" + partitionStrings[pccounter] + "\");");
		ExecuteCommands ("Tree tree_" + pccounter + " = " + treeStrings[pccounter] + ";");
		if (pccounter)
		{
			lfDef * ",";
		}
		lfDef  * ("part_"+pccounter+",tree_"+pccounter);
	}
	lfDef  * ");";
	lfDef  * 0;
	ExecuteCommands (lfDef);
	Optimize (res2, multiPart);

	/*fprintf (stdout, multiPart,"\n",res,"\n", res2);*/
	
	myDF = baseParams + res2[1][1];

	fullCAIC = -2(res2[1][0]-myDF*filteredData.sites/(filteredData.sites-myDF-1));
	/*fprintf  (stdout, "\n\nc-AIC = ", fullCAIC, "\nDelta AIC = ", nullCAIC-fullCAIC,"\n\n");*/
	outAVL ["DAIC"]  = nullCAIC-fullCAIC;
	outAVL ["DLogL"] = res2[1][0]-res[1][0];
	
	if (outAVL ["DAIC"] < 0)
	{
		return outAVL;
	}
	
	fragMatrix = {2,readPCount+1};

	for (pccounter = 0; pccounter <=  readPCount; pccounter = pccounter + 1)
	{
		ExecuteCommands ("siteCount = part_" + pccounter + ".sites;");
		/*fprintf (stdout, "Partition ", pccounter+1, " : ", siteCount, " sites\n");*/
		fragMatrix [0][pccounter] = siteCount;
	}

	pairWiseSplitMatch = {readPCount+1,readPCount+1};

	for (pccounter = 0; pccounter <=  readPCount; pccounter = pccounter + 1)
	{
		for (pc2 = pccounter+1; pc2 <=  readPCount; pc2 = pc2 + 1)
		{
			ExecuteCommands ("siteCount = compareTreesForSplits(\"tree_" +pccounter + "\",\"tree_"+pc2+"\");");
			pairWiseSplitMatch[pccounter][pc2] = siteCount[2]/Min(siteCount[0],siteCount[1]);
		}
	}

	/* now do SH testing */

	conditionals 	 = {};
	likelihoodScores = {readPCount+1,readPCount+1};
	pairwiseP		 = {readPCount+1,readPCount+1};

	treeSplitMatches 	 = 0;
	khIterations		 = 1000;

	for (pccounter = 0; pccounter <=  readPCount; pccounter = pccounter + 1)
	{
		for (pc2 = 0; pc2 <=  readPCount; pc2 = pc2 + 1)
		{
			if (Abs(pc2-pccounter) <= 1)
			{
				DataSetFilter 		aPart = CreateFilter (ds,1,partitionStrings[pccounter]);
				Tree		  		aTree = treeStrings[pc2];
				LikelihoodFunction	aLF	= (aPart,aTree);
				
				/*fprintf 		(stdout, "\n\nFitting tree ", pc2+1, " to partition ", pccounter+1, "\n");*/
				Optimize 		(aRes, aLF);	
				/*fprintf			(stdout, aLF);*/
				GetInformation  (cI,c);
				cI 				= cI[1][-1];
				ConstructCategoryMatrix (conds, aLF);
				conds 			= Log(cI*conds);
				conditionals [""+pccounter+","+pc2] = conds;
				likelihoodScores[pccounter][pc2] = aRes[1][0];
				treeSplitMatches = treeSplitMatches + pairWiseSplitMatch[pccounter][pc2];
			}
		}
		
		/*fprintf (stdout, "\nKH Testing partition ", pccounter+1, "\n");*/
		partTreeConds = conditionals[""+pccounter+","+pccounter];
		
		for (pc2 = 0; pc2 <=  readPCount; pc2 = pc2 + 1)
		{
			if (Abs(pc2-pccounter) == 1)
			{
				otherPartTree = conditionals[""+pccounter+","+pc2];
				baseLRT = 2*(likelihoodScores[pccounter][pccounter]-likelihoodScores[pccounter][pc2]);
				/*fprintf (stdout, "Tree ", pc2+1, " base LRT = ", baseLRT, ". p-value = ");*/
				textMx = testLRT(partTreeConds,otherPartTree,khIterations) % 0;
				for (kk=0; kk<khIterations; kk=kk+1)
				{	
					if (textMx[kk]>=0)
					{
						break;
					}
				}
				pval = Max(1,kk)/khIterations;
				/*
				fprintf (stdout, pval , "\n");
				*/
				
				pairwiseP[pccounter][pc2] = pval;
			}
		}
	}

	treeSplitMatches = treeSplitMatches*2/readPCount/(readPCount+1);
	totalComparisons = readPCount*2;
	threshP			 = 0.01/totalComparisons;

	significantCouplings = 0;

	for (pccounter = 1; pccounter <=  readPCount; pccounter = pccounter + 1)
	{
		if (pairwiseP[pccounter][pccounter-1] <= threshP && pairwiseP[pccounter-1][pccounter] <= threshP)
		{
			significantCouplings = significantCouplings + 1;
		}
		fragMatrix [1][pccounter] = totalComparisons*Max (pairwiseP[pccounter][pccounter-1],pairwiseP[pccounter-1][pccounter]);
	}

	/*fprintf (stdout, "A total of ", significantCouplings, "/", readPCount, " significant couplings\n");
	fprintf (stdout, "Mean splits identify: ",treeSplitMatches, "\n" );*/
	outAVL["SIG_CPL"] = significantCouplings;
	outAVL["SPLITS"] = treeSplitMatches;
	outAVL["Fragments"] = fragMatrix;

	return outAVL;
}

/* ________________________________________________________________________________________________*/

function computeMeanDivergence (treeID)
{
	ExecuteCommands ("treeAVL2="+treeID+"^0;leafCount=TipCount("+treeID+")");
	multFactors = {};
	
	for (k=1; k<Abs(treeAVL2); k=k+1)
	{
		aNode = treeAVL2[k];
		aNodeName = aNode["Name"];
		parentIndex = aNode["Parent"];
		k2 = Abs(aNode["Children"]);
		if (k2)
		{
			currentDepth = aNode["Below"];
			multFactors[aNodeName] = currentDepth;		
			if (parentIndex > 0)
			{
				pInfo = treeAVL2[parentIndex];
				pInfo ["Below"] = pInfo ["Below"] + currentDepth;
				treeAVL2[parentIndex] = pInfo;
			}
		}
		else
		{
			multFactors[aNodeName] = 1;
			pInfo = treeAVL2[parentIndex];
			pInfo ["Below"] = pInfo ["Below"] + 1;
			treeAVL2[parentIndex] = pInfo;
		}
	}

	pKeys 			= Rows(multFactors);

	for (k=0; k<Columns(pKeys); k=k+1)
	{
		aNodeName = pKeys[k];
		multFactors[aNodeName] = multFactors[aNodeName] * (leafCount-multFactors[aNodeName]);
	}
	
	ExecuteCommands ("bNames = BranchName   ("+treeID+",-1);");
	ExecuteCommands ("bLen   = BranchLength ("+treeID+",-1);");
	
	sum = 0;
	bc  = Columns(bNames)-1;
	
	for (k=0; k<bc; k=k+1)
	{
		aNodeName = bNames[k];
		sum = sum + bLen[k]*multFactors[aNodeName];
	}	
	return sum/(leafCount*(leafCount-1)*0.5);
}

/* ________________________________________________________________________________________________*/

function compareTreesForSplits (tName, tName2)
{
	/* obtain an AVL data structure of the tree, post-order layout */

	ExecuteCommands   ("_astavl_="+tName+"^1;");
	_tree_size_ = Abs (_astavl_);


	nodeMapAVL = {};
	
	for (_a_node = 2; _a_node < _tree_size_; _a_node = _a_node + 1)
	{
		_node_info = _astavl_[_a_node];
		myDegree = Abs(_node_info["Children"]);
		
		if (myDegree == 0)
		{
			nodeName = _node_info["Name"];
			nodeMapAVL [nodeName] = Abs(nodeMapAVL)+1;
		}
	}
	
	split1 = getSplits(0);
	ExecuteCommands   ("_astavl_="+tName2+"^1;");
	_tree_size_ = Abs (_astavl_);
	split2 = getSplits(0);
	
	s1c = Abs(split1);
	s2c = Abs(split2);
		
	matches = {};
	match1  = {};
	match2  = {};
	for (r1 = 0; r1 < s1c; r1 = r1 + 1)
	{
		if (match1[r1] == 0)
		{
			for (r2 = 0; r2 < s2c; r2 = r2 + 1)
			{
				if (match2[r2] == 0)
				{
					splitCheck = compareSplits (split1[r1],split2[r2]);
					if (splitCheck)
					{
						mr = {{r1,r2}};
						matches [Abs(matches)] = mr;
						match1[r1] = 1;
						match2[r2] = 1;
					}
				}
			}
		}
	}
	
	
	return {{s1c,s2c,Abs(matches)}};
}

/* ________________________________________________________________________________________________*/

function getSplits (dummy)
{
	treeSplits = {};
	for (_a_node = 2; _a_node < _tree_size_; _a_node = _a_node + 1)
	{
		_node_info = _astavl_[_a_node];
		myDegree = Abs(_node_info["Children"]);
		myDepth  = _node_info["Depth"];
		
			
		if (myDegree && _node_info["Length"]>0)
		{
			nodeName = _node_info["Name"];
			aSplit = {Abs(nodeMapAVL),1};
			for (_b_node = _a_node + 1; _b_node < _tree_size_; _b_node = _b_node + 1)
			{
				_bnode_info = _astavl_[_b_node];
					if (_bnode_info["Depth"] <= myDepth)
				{
					break;
				}
				if (Abs(_bnode_info["Children"])==0)
				{
					_bnode_info = _bnode_info["Name"];
					_bnode_info = nodeMapAVL[_bnode_info];
					if (_bnode_info == 0)
					{
						fprintf (stdout, "Error: extraneous node name\n");
						return 0;
					}
					aSplit [_bnode_info-1] = 1;
				}
			}
			treeSplits [Abs(treeSplits)] = aSplit;
		}
	}
	return treeSplits;
}

/* ________________________________________________________________________________________________*/

function compareSplits (s1, s2)
{
	sl            = Rows(s1);
	positiveCheck = (s1[0] == s2[0]);
	for (k=1; k<sl; k=k+1)
	{
		if (s1[k] != s2[k])
		{
			if (positiveCheck)
			{
				return 0;
			}
		}
		else
		{
			if (!positiveCheck)
			{
				return 0;
			}
		}
	}
	return -1+2*positiveCheck;	
}

/*--------------------------------------------------------------------------------------*/

function testLRT (vec1, vec2, itCount)
{
	size1 = Columns(vec1);
	
	sumVec1 = {size1,1};
	jvec	= {2,size1};
	resMx1	= {itCount,1};
	resMx2	= {itCount,1};
	
	for (k=0; k<size1; k=k+1)
	{
		sumVec1 [k]	   = 1;
		jvec	[0][k] = vec1[k];
		jvec	[1][k] = vec2[k];
	}
	
	
	for (k=0; k<itCount; k=k+1)
	{
		resampled = Random(jvec,1);
		resampled = resampled*sumVec1;
		resMx1[k] = resampled[0];
		resMx2[k] = resampled[1];
	}
	
	return (resMx1-resMx2)*2;
}
