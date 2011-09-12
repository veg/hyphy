ChoiceList (testType,"Test for",1,SKIP_NONE,"Mean branch length","Compare mean branch lengths between two or more non-nested clades.",
				     "Mean pairwise divergence","Compare mean within-clade pairwise divergence between two or more non-nested clades.");

if (testType < 0)
{
	return 0;
}

if (testType)
{
	echoString = "mean within-clade divergence";
}
else
{
	echoString = "mean within-clade branch length (or component branch length)";
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
	incFileName = HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"chooseGeneticCode.def";
	ExecuteCommands  ("#include \""+incFileName+"\";");
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
	DataSetFilter filteredData = CreateFilter (ds,1);
}

fprintf (stdout,"\n______________READ THE FOLLOWING DATA______________\n",ds);


_DO_TREE_REBALANCE_ = 0;

SelectTemplateModel(filteredData);

if (Rows("LAST_MODEL_PARAMETER_LIST")>1)
{
	ChoiceList (parameter2Constrain, "Parameter(s) to constrain:",1,SKIP_NONE,LAST_MODEL_PARAMETER_LIST);

	if (parameter2Constrain<0)
	{
		return 0;
	}
	if (parameter2Constrain==0)
	{
		fprintf (stdout, "ERROR: Multiple parameter constraints are not supported by this analysis; sorry!\n");
		return 0;
	}
}
else
{
	parameter2Constrain = 1;
}

GetString (funnyString,LAST_MODEL_PARAMETER_LIST,parameter2Constrain-1);

_DO_TREE_REBALANCE_ = 0;
incFileName = HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"queryTree.bf";
ExecuteCommands  ("#include \""+incFileName+"\";");

treeAVL  = givenTree ^ 1; /* pre-order traversal */

if (testType)
{
	treeAVL2 = givenTree ^ 0; /* post-order traversal */

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
	treeAVL2 		= 0;
}

branchNameList 	= BranchName (givenTree,-1);
mapToBL		    = {};

for (k=0; k<Columns(branchNameList); k=k+1)
{
	aNode = branchNameList[k];
	mapToBL [aNode] = k;
}

cladeNodeLists 	= {};
cladeIDLists    = {};
cladeSumVar		= {};
nameToAVLIdx	= {};
cladeLeafCount  = {};

for (k=1; k<Abs(treeAVL); k=k+1)
{
	aNode 	  = treeAVL[k];
	aNodeName = aNode["Name"]; /*upcase*/
	nameToAVLIdx [aNodeName] = k;
	match = (aNodeName&&1)$"^CLADE[0-9]+$";
	if (match[0]>=0 && Abs (aNode["Children"])>0)
	/* valid clade */
	{
		currentDepth = aNode["Depth"];		
		nodeList = {};
		
		leafCount = 0;
		
		for (k2=k+1; k2<Abs(treeAVL); k2=k2+1)
		{
			aNode2 = treeAVL[k2];
			if (aNode2["Depth"] <= currentDepth)
			{
				break;
			}
			if (Abs(aNode2["Children"]) == 0)
			{
				leafCount = leafCount + 1;
			}
			minCladeIndex = aNode2["Name"];
			nodeList[Abs(nodeList)] = minCladeIndex;
			nameToAVLIdx [minCladeIndex] = k2;
		}
				
		k = k2-1;
		cladesFound = Abs(nodeList);
		cladeNodeLists[aNodeName] = nodeList;
		indexMatrix = {cladesFound,1};
		for (k2 = 0; k2 < cladesFound; k2=k2+1)
		{
			aNode2 = nodeList[k2];
			indexMatrix [k2] = mapToBL[aNode2];
		}
		
		cladeIDLists[aNodeName] = indexMatrix;
		cladeLeafCount[aNodeName] = leafCount;
	}
}


cladesFound = Abs (cladeNodeLists);
if (cladesFound < 2)
{ 
	fprintf (stdout, "ERROR: Couldn't find 2 or more (non-nested) clades with roots labeled by CladeN, where N is a number. Found ", cladesFound, " labeled clades\n");
	return 0;
}

cladeRoots = Rows (cladeNodeLists);

minCladeSize  = 1e100;
minCladeIndex = 0;

mFactors  = {};
bmFactors = {};

fprintf (stdout, "\n\nFound ", cladesFound, " labeled clades\n");
for (k=0; k<cladesFound; k=k+1)
{
	cladeName = cladeRoots[k];
	cladeSize =  Abs(cladeNodeLists[cladeName]);
	leafCount =  Abs(cladeLeafCount[cladeName]);
	
	fprintf (stdout, "Clade rooted at ", cladeName, " with ", cladeSize , " branches and ", leafCount, " leaves.\n");
	nodeList = cladeNodeLists[cladeName];
	
	if (testType)
	{
		if (minCladeSize > leafCount)
		{
			minCladeSize = leafCount;
			minCladeIndex = k;
		}
	}
	else
	{
		if (minCladeSize > cladeSize)
		{
			minCladeSize = cladeSize;
			minCladeIndex = k;
		}
	}
	
	bmFactor = {Columns(branchNameList),1};
	cString = "";
	cString * 256;
	for (k2 = 0; k2<Abs(nodeList); k2=k2+1)
	{
		aNodeName = nodeList[k2];
		aNode = nameToAVLIdx[aNodeName];
		aNode = treeAVL[aNode];
		if (testType)
		{
			cCount = Abs(aNode["Children"]);
			if (cCount)
			{	
				cCount = multFactors[aNodeName];
				mFactors[aNodeName] = (leafCount-cCount)*cCount
			}
			else
			{
				mFactors[aNodeName] = leafCount-1;
			}
		}
		else
		{
			mFactors[aNodeName] = 1;
		}
		k3 = mapToBL[aNodeName];
		bmFactor [k3] = mFactors[aNodeName];
		if (k2)
		{
			cString * "+";
		}
		cString * ("TREE_PLACE_HOLDER."+aNodeName + "." + funnyString + "*" + mFactors[aNodeName]);
	}
	cString * 0;
	cladeSumVar[cladeName] = cString;
	bmFactors[cladeName] = bmFactor;
	
}


fprintf (stdout, "1). Running an unconstrained fit\n");

Tree aTree = treeString;

USE_LAST_RESULTS = 0;

LikelihoodFunction lf = (filteredData,aTree);

Optimize (res_free,lf);
fprintf  (stdout, lf, "\n");

blVector = BranchLength (aTree,-1);
ReportMeans (blVector);

USE_LAST_RESULTS = 1;

fprintf (stdout, "2). Running a completely constrained fit\n");

Tree mirrorTree = treeString;


ReplicateConstraint ("this1.?.?:=this2.?.?",mirrorTree,aTree);

ClearConstraints (mirrorTree);

for (k=0; k<cladesFound; k=k+1)
{
	if (k==minCladeIndex)
	{
		smallCladeName = cladeRoots[minCladeIndex];
		cladeConstraint = cladeSumVar[smallCladeName] ^ {{"TREE_PLACE_HOLDER","aTree"}};
		ExecuteCommands ("global SUM_"+smallCladeName+":="+cladeConstraint+";");
	}
	else
	{
		cladeName = cladeRoots[k];
		cladeConstraint = cladeSumVar[cladeName] ^ {{"TREE_PLACE_HOLDER","mirrorTree"}};
		ExecuteCommands ("global SUM_"+cladeName+":="+cladeConstraint+";");	
	}
}


for (k=0; k<cladesFound; k=k+1)
{
	if (k!=minCladeIndex)
	{
		cladeName = cladeRoots[k];
		if (testType)
		{
			cladeSize =  cladeLeafCount[cladeName];
			ExecuteCommands ("ReplicateConstraint (\"this1.?."+funnyString+":= SUM_"+smallCladeName+"*"+(cladeSize*(cladeSize-1)/minCladeSize/(minCladeSize-1))+"/SUM_"+cladeName+"*this2.?."+funnyString+"\",aTree."+cladeName+",mirrorTree."+cladeName+");");
		}
		else
		{
			cladeSize =  Abs(cladeNodeLists[cladeName]);
			ExecuteCommands ("ReplicateConstraint (\"this1.?."+funnyString+":= SUM_"+smallCladeName+"*"+(cladeSize/minCladeSize)+"/SUM_"+cladeName+"*this2.?."+funnyString+"\",aTree."+cladeName+",mirrorTree."+cladeName+");");		
		}
		ExecuteCommands ("aTree."+cladeName+"."+funnyString+"=0.1;");		
	}
}

LikelihoodFunction lf = (filteredData,aTree);
Optimize (res_contsr,lf);

fprintf  (stdout, "\n\nLRT for difference in mean divergences\n",
				  "\nLR = ", 2*(res_free[1][0]-res_contsr[1][0]), 
				  "\np-value = ", 1-CChi2(2*(res_free[1][0]-res_contsr[1][0]),cladesFound-1),"\n");
				  
fprintf  (stdout, lf, "\n");

blVector = BranchLength (aTree,-1);
ReportMeans (blVector);

comparisonIndex = 3;

ClearConstraints (aTree);


if (cladesFound > 2)
{
	for (c1 = 0; c1<cladesFound-1; c1 = c1+1)
	{
		smallCladeName = cladeRoots[c1];
		if (testType)
		{
			minCladeSize =  cladeLeafCount[smallCladeName];
		}
		else
		{
			minCladeSize =  Abs(cladeNodeLists[smallCladeName]);
		}

		cladeConstraint = cladeSumVar[smallCladeName] ^ {{"TREE_PLACE_HOLDER","aTree"}};
		ExecuteCommands ("global SUM_"+smallCladeName+":="+cladeConstraint+";");
		
		for (c2 = c1+1; c2 < cladesFound; c2=c2+1)
		{
			cladeName = cladeRoots[c2];

			cladeConstraint = cladeSumVar[cladeName] ^ {{"TREE_PLACE_HOLDER","mirrorTree"}};
			ExecuteCommands ("global SUM_"+cladeName+":="+cladeConstraint+";");	

			if (testType)
			{
				cladeSize =  cladeLeafCount[cladeName];
				ExecuteCommands ("ReplicateConstraint (\"this1.?."+funnyString+":= SUM_"+smallCladeName+"*"+(cladeSize*(cladeSize-1)/minCladeSize/(minCladeSize-1))+"/SUM_"+cladeName+"*this2.?."+funnyString+"\",aTree."+cladeName+",mirrorTree."+cladeName+");");
			}
			else
			{
				cladeSize =  Abs(cladeNodeLists[cladeName]);
				ExecuteCommands ("ReplicateConstraint (\"this1.?."+funnyString+":= SUM_"+smallCladeName+"*"+(cladeSize/minCladeSize)+"/SUM_"+cladeName+"*this2.?."+funnyString+"\",aTree."+cladeName+",mirrorTree."+cladeName+");");
			}
			
			ExecuteCommands ("aTree."+cladeName+".t=0.1;");		
			fprintf (stdout, comparisonIndex, "). Running a constrained fit on clades ", smallCladeName, " and ", cladeName, "\n");
		
			LikelihoodFunction lf = (filteredData,aTree);
			Optimize (res_2c,lf);

fprintf  (stdout, "\n\nLRT for difference in mean divergences\n",
				  "\nLR = ", 2*(res_free[1][0]-res_2c[1][0]), 
				  "\np-value = ", 1-CChi2(2*(res_free[1][0]-res_2c[1][0]),1),"\n");

			fprintf  (stdout, lf, "\n");

			blVector = BranchLength (aTree,-1);
			ReportMeans (blVector);
			comparisonIndex = comparisonIndex + 1;
			
			ClearConstraints (aTree);
		}
		
	}
}

USE_LAST_RESULTS = 0;

/*--------------------------------------------------------------------------------------------------------------*/

function ReportMeans (branches)
{
	for (k=0; k<cladesFound; k=k+1)
	{
		cladeName = cladeRoots[k];
		sum = branches*bmFactors[cladeName];
		
		if (testType)
		{
			leafCount = cladeLeafCount[cladeName];
			sum = sum[0]*2/leafCount/(leafCount-1);
		}
		else
		{
			cladeSize =  Abs(cladeNodeLists[cladeName]);
			sum = sum[0]/cladeSize;		
		}
				
		fprintf (stdout, "Clade rooted at ", cladeName, " has ", echoString, " of ", sum, "\n");
		
	}
	return 0;
}
