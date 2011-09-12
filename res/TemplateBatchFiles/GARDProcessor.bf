VERBOSITY_LEVEL = -1;

SetDialogPrompt ("Please load a nucleotide data file:");
DataSet 	ds = ReadDataFile (PROMPT_FOR_FILE);
baselineDSPath = LAST_FILE_PATH;

SetDialogPrompt ("Please load a GA partition analysis result file:");
fscanf  (PROMPT_FOR_FILE,REWIND,"Lines",partLines);

readPCount = Columns (partLines)/2-1;
fprintf (stdout, "\nLoaded ", readPCount , " partitions\n");

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

/*AT:=AC;
CG:=AC;
CT:=1;
GT:=AC;*/

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

fprintf (stdout, "\nSequences :", filteredData.species,
				 "\nSites     :", filteredData.sites,
				 "\nVariable  :", Abs(bppMap), "\n"); 

HarvestFrequencies (baseFreqs,ds,1,1,1);

fprintf (stdout, "\n\nf(A) = ", baseFreqs[0],
				 "\nf(C) = ", baseFreqs[1],
				 "\nf(G) = ", baseFreqs[2],
				 "\nf(T) = ", baseFreqs[3],"\n");
				
Model GTR_Model =  (GTR_Matrix, baseFreqs, 1);

InferTreeTopology   (0);
treeString 	= TreeMatrix2TreeString(0);
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

fprintf  (stdout, "\n\nc-AIC = ", nullCAIC, "\n", singlePart);

m1 = computeMeanDivergence ("givenTree");

fprintf  (stdout, "\nMean divergence : ", m1*100, "%\n");

fprintf  (stdout, "\n\nFitting a single-tree, multiple partition model\n");

USE_DISTANCES = 0;

lfDef = "";
lfDef * 128;
lfDef  * "LikelihoodFunction multiPart  = (";

for (pccounter = 0; pccounter < readPCount; pccounter = pccounter + 1)
{
	ExecuteCommands ("DataSetFilter part_" + pccounter + " = CreateFilter (ds, 1, \"" + partLines[2+pccounter*2] + "\");");
	ExecuteCommands ("Tree tree_" + pccounter + " = treeString;");
	if (pccounter)
	{
		lfDef * ",";
		ExecuteCommands ("global S_"+pccounter+"=1;\nReplicateConstraint (\"this1.?.?:=S_"+pccounter+"*this2.?.?\",tree_"+pccounter+",tree_0);");
	}
	lfDef  * ("part_"+pccounter+",tree_"+pccounter);
}
lfDef  * ");";
lfDef  * 0;
ExecuteCommands (lfDef);
Optimize (res2, multiPart);
myDF2 = baseParams + res2[1][1];
intCAIC = -2(res2[1][0]-myDF2*filteredData.sites/(filteredData.sites-myDF2-1));
fprintf (stdout, multiPart);

USE_DISTANCES = 0;

lfDef = "";
lfDef * 128;
lfDef  * "LikelihoodFunction multiPart  = (";

fprintf  (stdout, "\n\nFitting a mutilple tree, multiple partition model\n");

for (pccounter = 0; pccounter < readPCount; pccounter = pccounter + 1)
{
	ExecuteCommands ("DataSetFilter part_" + pccounter + " = CreateFilter (ds, 1, \"" + partLines[2+pccounter*2] + "\");");
	ExecuteCommands ("Tree tree_" + pccounter + " = " + partLines[pccounter*2+3] + ";");
	if (pccounter)
	{
		lfDef * ",";
	}
	lfDef  * ("part_"+pccounter+",tree_"+pccounter);
}
lfDef  * ");";
lfDef  * 0;
ExecuteCommands (lfDef);
Optimize (res, multiPart);

fprintf (stdout, multiPart);

lfout = baselineDSPath + "_multi.fit";
LIKELIHOOD_FUNCTION_OUTPUT = 7;
fprintf (lfout, CLEAR_FILE, multiPart);
LIKELIHOOD_FUNCTION_OUTPUT = 2;

myDF = baseParams + res[1][1];

fullCAIC = -2(res[1][0]-myDF*filteredData.sites/(filteredData.sites-myDF-1));

fprintf  (stdout, "\n\nVersus the single partition model: c-AIC = ", fullCAIC, "\nDelta AIC = ", nullCAIC-fullCAIC,"\n\n");
fprintf  (stdout, "\n\nVersus the single tree/multiple partition model: Delta AIC = ", intCAIC-fullCAIC,"\n\n");

bpLocations = {readPCount, 1}; 
for (pccounter = 0; pccounter <  readPCount; pccounter = pccounter + 1)
{
	if (pccounter)
	{
		bpLocations[pccounter] = siteCount + bpLocations[pccounter-1];
	}
	ExecuteCommands ("siteCount = part_" + pccounter + ".sites;");
	fprintf (stdout, "Partition ", pccounter+1, " : ", siteCount, " sites\n");
}

pairWiseSplitMatch = {readPCount,readPCount};

for (pccounter = 0; pccounter <  readPCount; pccounter = pccounter + 1)
{
	for (pc2 = pccounter+1; pc2 <  readPCount; pc2 = pc2 + 1)
	{
		ExecuteCommands ("siteCount = compareTreesForSplits(\"tree_" +pccounter + "\",\"tree_"+pc2+"\");");
		pairWiseSplitMatch[pccounter][pc2] = siteCount[2]/Min(siteCount[0],siteCount[1]);
	}
}

/* now do SH testing */

conditionals 	 = {};
likelihoodScores = {readPCount,readPCount};
pairwiseP		 = {readPCount,readPCount};

treeSplitMatches 	 = 0;
khIterations		 = 10000;

if (OPTIMIZATION_PRECISION == 0)
{
	cutThreshold		= 0.001;
}
else
{
	cutThreshold		= 2 * OPTIMIZATION_PRECISION;
}


for (pccounter = 0; pccounter <  readPCount; pccounter = pccounter + 1)
{
	for (pc2 = 0; pc2 <  readPCount; pc2 = pc2 + 1)
	{
		if (Abs(pc2-pccounter) <= 1)
		{
			DataSetFilter 		aPart = CreateFilter (ds,1,partLines[pccounter*2+2]);
			Tree		  		aTree = partLines[pc2*2+3];
			LikelihoodFunction	aLF	= (aPart,aTree);
			
			fprintf				(stdout, "\n\nFitting tree ", pc2+1, " to partition ", pccounter+1, "\n");
			Optimize			(aRes, aLF);	
			fprintf				(stdout, aLF);
			LIKELIHOOD_FUNCTION_OUTPUT = 2;
			GetInformation		(cI,c);
			cI					= cI[1][-1];
			ConstructCategoryMatrix (conds, aLF);
			conds				= Log(cI*conds);
			conditionals		[pccounter*readPCount+pc2] = conds;
			likelihoodScores	[pccounter][pc2] = aRes[1][0];
			treeSplitMatches    = treeSplitMatches + pairWiseSplitMatch[pccounter][pc2];
		}
	}
	
	fprintf (stdout, "\nKH Testing partition ", pccounter+1, "\n");
	partTreeConds = conditionals[pccounter*readPCount+pccounter];
	
	for (pc2 = 0; pc2 <  readPCount; pc2 = pc2 + 1)
	{
		if (Abs(pc2-pccounter) == 1)
		{
			otherPartTree = conditionals[pccounter*readPCount+pc2];
			baseLRT = 2*(likelihoodScores[pccounter][pccounter]-likelihoodScores[pccounter][pc2]);
			fprintf (stdout, "Tree ", pc2+1, " base LRT = ", baseLRT, ". p-value = ");
			textMx = testLRT(partTreeConds,otherPartTree,khIterations) % 0;
			for (kk=khIterations-1; kk>=0; kk=kk-1)
			{	
				if (textMx[kk] < baseLRT)
				{
					break;
				}
			}
			pval = Max(1,(khIterations-1-kk))/khIterations;
			fprintf (stdout, pval , "\n");
			pairwiseP[pccounter][pc2] = pval;
		}
	}
}

treeSplitMatches = treeSplitMatches*2/readPCount/(readPCount-1);
totalComparisons = (readPCount-1)*2;
threshP			 = 0.01/totalComparisons;

bypvalue         = {{0.01, 0}{0.05, 0}{0.1, 0}} * (1/totalComparisons);


fprintf (stdout, "\nBreakpoint | LHS Raw p | LHS adjusted p | RHS Raw p | RHS adjusted p \n");

for (pccounter = 1; pccounter <  readPCount; pccounter = pccounter + 1)
{
	lhs     = pairwiseP[pccounter][pccounter-1];
	rhs     = pairwiseP[pccounter-1][pccounter];
						 
	for (k = 0; k < Rows (bypvalue); k = k + 1)
	{
		threshP = bypvalue[k][0];
		if (lhs <= threshP && rhs <= threshP)
		{
			bypvalue[k][1] = bypvalue[k][1] + 1;
		}
	}
	fprintf (stdout, Format(bpLocations[pccounter],10,0), " | ",
					 Format(lhs,9,5), " | ",
					 Format(Min(1,lhs*totalComparisons),14,5), " | ",
					 Format(rhs,9,5), " | ",
					 Format(Min(1,rhs*totalComparisons),14,5), "\n");
}

fprintf (stdout, "\n");
for (k = 0; k < Rows (bypvalue); k = k + 1)
{
	fprintf (stdout, "At p = ", bypvalue[k][0]*totalComparisons, " there are ", bypvalue[k][1], " significant breakpoints\n");
}
fprintf (stdout, "\nMean splits identify: ",treeSplitMatches, "\n" );







/* AUX STUFF */

/*lfDef = "";
lfDef * 128;
lfDef  * "LikelihoodFunction multiPart2  = (";

for (pccounter = 0; pccounter < readPCount; pccounter = pccounter + 1)
{
	ExecuteCommands ("DataSetFilter part2_" + pccounter + " = CreateFilter (ds, 1, \"" + partLines[pccounter*2+2] + "\");");
	ExecuteCommands ("Tree tree2_" + pccounter + " = " + treeString + ";");
	if (pccounter)
	{
		lfDef * ",";
	}
	lfDef  * ("part2_"+pccounter+",tree2_"+pccounter);
}

lfDef  * ");";
lfDef  * 0;
ExecuteCommands (lfDef);
Optimize (res, multiPart2);

fprintf (stdout, multiPart2);

full2CAIC = -2(res[1][0]-res[1][1]*filteredData.sites/(filteredData.sites-res[1][1]-1));
fprintf  (stdout, "\n\nc-AIC = ", full2CAIC, "\nDelta Alt-AIC = ", full2CAIC-fullCAIC,"\n\n");*/


/* ________________________________________________________________________________________________*/


function InitializeDistances (dummy)
{
	HarvestFrequencies (_dNucFreq,filteredData,1,1,0);
	_d_fR = _dNucFreq[0]+_dNucFreq[2];
	_d_fY = _dNucFreq[1]+_dNucFreq[3];
	
	if (_dNucFreq[0] == 0 || _dNucFreq[1] == 0 || _dNucFreq[2] == 0 || _dNucFreq[3] == 0)
	{
		_useK2P = 1;
	}
	else
	{
		_d_TN_K1 = 2*_dNucFreq[0]*_dNucFreq[2]/_d_fR;
		_d_TN_K2 = 2*_dNucFreq[1]*_dNucFreq[3]/_d_fY;
		_d_TN_K3 = 2*(_d_fR*_d_fY-_dNucFreq[0]*_dNucFreq[2]*_d_fY/_d_fR-_dNucFreq[1]*_dNucFreq[3]*_d_fR/_d_fY);
		_useK2P = 0;
	}
	
	
	summingVector = {{1}{1}{1}{1}};

	return 0;
}

/* ________________________________________________________________________________________________*/

function ComputeDistanceFormula (s1,s2)
{
	GetDataInfo (siteDifferenceCount, filteredData, s1, s2, DIST);
	
	totalSitesCompared = Transpose(summingVector)*(siteDifferenceCount*summingVector);
	totalSitesCompared = totalSitesCompared[0];
	
	if (_useK2P)
	{
		_dTransitionCounts 	 =    siteDifferenceCount[0][2]+siteDifferenceCount[2][0]  /* A-G and G-A */
								 +siteDifferenceCount[1][3]+siteDifferenceCount[3][1]; /* C-T and T-C */
							
		_dTransversionCounts = (siteDifferenceCount[0][0]+siteDifferenceCount[1][1]+siteDifferenceCount[2][2]+siteDifferenceCount[3][3])+_dTransitionCounts;
		
		_dTransitionCounts	 = _dTransitionCounts/totalSitesCompared;
		_dTransversionCounts = 1-_dTransversionCounts/totalSitesCompared;
		
		_d1C = 1-2*_dTransitionCounts-_dTransversionCounts;
		_d2C = 1-2*_dTransversionCounts;
		
		if (_d1C>0 && _d2C>0)
		{
			return -(0.5*Log(_d1C)+.25*Log(_d2C));	
		}
	}
	else
	{
		_dAGCounts 	 =    siteDifferenceCount[0][2]+siteDifferenceCount[2][0]  /* A-G and G-A */;
		_dCTCounts	 = 	  siteDifferenceCount[1][3]+siteDifferenceCount[3][1]; /* C-T and T-C */
							
		_dTransversionCounts = (siteDifferenceCount[0][0]+siteDifferenceCount[1][1]+siteDifferenceCount[2][2]+
								siteDifferenceCount[3][3])+_dAGCounts+_dCTCounts;
		
		_dAGCounts	 = _dAGCounts/totalSitesCompared;
		_dCTCounts	 = _dCTCounts/totalSitesCompared;
		
		_dTransversionCounts = 1-_dTransversionCounts/totalSitesCompared;
		
		_d1C = 1-_dAGCounts/_d_TN_K1-0.5*_dTransversionCounts/_d_fR;
		_d2C = 1-_dCTCounts/_d_TN_K2-0.5*_dTransversionCounts/_d_fY;
		_d3C = 1-0.5*_dTransversionCounts/_d_fY/_d_fR;
		
		if ((_d1C>0)&&(_d2C>0)&&(_d3C>0))
		{
			return -_d_TN_K1*Log(_d1C)-_d_TN_K2*Log(_d2C)-_d_TN_K3*Log(_d3C);
		}
	}
	
	return 1000;
}

/* ________________________________________________________________________________________________*/

function TreeMatrix2TreeString (doLengths)
{
	treeString = "";
	p = 0;
	k = 0;
	m = treeNodes[0][1];
	n = treeNodes[0][0];
	d = treeString*(Rows(treeNodes)*25);

	while (m)
	{	
		if (m>p)
		{
			if (p)
			{
				d = treeString*",";
			}
			for (j=p;j<m;j=j+1)
			{
				d = treeString*"(";
			}
		}
		else
		{
			if (m<p)
			{
				for (j=m;j<p;j=j+1)
				{
					d = treeString*")";
				}
			}	
			else
			{
				d = treeString*",";
			}	
		}
		if (n<filteredData.species)
		{
			GetString (nodeName, filteredData, n);
			d = treeString*nodeName;
		}
		if (doLengths>.5)
		{
			nodeName = ":"+treeNodes[k][2];
			d = treeString*nodeName;
		}
		k=k+1;
		p=m;
		n=treeNodes[k][0];
		m=treeNodes[k][1];
	}

	for (j=m;j<p;j=j+1)
	{
		d = treeString*")";
	}
	
	d=treeString*0;
	return treeString;
}

/* ________________________________________________________________________________________________*/

function InferTreeTopology(verbFlag)
{
	distanceMatrix = {filteredData.species,filteredData.species};
	dummy = InitializeDistances (0);
		
	for (i = 0; i<filteredData.species; i=i+1)
	{
		for (j = i+1; j<filteredData.species; j = j+1)
		{
			distanceMatrix[i][j] = ComputeDistanceFormula (i,j);
		}
	}

	MESSAGE_LOGGING 		 	= 1;
	cladesMade 					= 1;
	

	if (filteredData.species == 2)
	{
		d1 = distanceMatrix[0][1]/2;
		treeNodes = {{0,1,d1__},
					 {1,1,d1__},
					 {2,0,0}};
					 
		cladesInfo = {{2,0}};
	}
	else
	{
		if (filteredData.species == 3)
		{
			d1 = (distanceMatrix[0][1]+distanceMatrix[0][2]-distanceMatrix[1][2])/2;
			d2 = (distanceMatrix[0][1]-distanceMatrix[0][2]+distanceMatrix[1][2])/2;
			d3 = (distanceMatrix[1][2]+distanceMatrix[0][2]-distanceMatrix[0][1])/2;
			treeNodes = {{0,1,d1__},
						 {1,1,d2__},
						 {2,1,d3__}
						 {3,0,0}};
						 
			cladesInfo = {{3,0}};		
		}
		else
		{	
			njm = (distanceMatrix > methodIndex)>=filteredData.species;
				
			treeNodes 		= {2*(filteredData.species+1),3};
			cladesInfo	    = {filteredData.species-1,2};
			
			for (i=Rows(treeNodes)-1; i>=0; i=i-1)
			{
				treeNodes[i][0] = njm[i][0];
				treeNodes[i][1] = njm[i][1];
				treeNodes[i][2] = njm[i][2];
			}

			for (i=Rows(cladesInfo)-1; i>=0; i=i-1)
			{
				cladesInfo[i][0] = njm[i][3];
				cladesInfo[i][1] = njm[i][4];
			}
			
			njm = 0;
		}
	}
	distanceMatrix = 0;
	
	return 1.0;
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
	
	jvec	= {2,size1};
	resMx1	= {itCount,1};
	resMx2	= {itCount,1};
	
	for (k=0; k<size1; k=k+1)
	{
		jvec	[0][k] = vec1[k];
		jvec	[1][k] = vec2[k];
	}
	
	
	for (k=0; k<itCount; k=k+1)
	{
		resampled = Random(jvec,1);
		resMx1[k] = +(resampled[0][-1]);
		resMx2[k] = +(resampled[1][-1]);
        SetParameter (STATUS_BAR_STATUS_STRING, "Drawing resampled likelihoods for segments " + pccounter + " and " + pc2 + "("+k+"/"+itCound+" done)",0);
	}
	
	resMx = (resMx1-resMx2)*2;
    return resMx + (+resMx) * (-1/itCount);
}
