/*--------------------------------------------------- */

function factorial (x)
{
	res = 1;
	for (x2 = 2; x2 <= x; x2=x2+1)
	{
		res = res * x2;
	}
	return res;
}

/*--------------------------------------------------- */

function computeMigrationEvents (assignmentVector)
{

	assignmentMatrices = {};
	nodeState 		   = {};

	for (specIndex = 0; specIndex < leafCount; specIndex = specIndex + 1)
	{
		if (assignmentVector[specIndex])
		{
			aMx = {{1,1}{1,0}};
		}
		else
		{
			aMx = {{0,0}{0,1}};	
		}
		specName = TipName (givenTree, specIndex);
		assignmentMatrices [specName] = aMx;
	}

	for (node = 1; node < treeSize; node = node + 1)
	{
		nodeInfo 		= treeAVL[node];
		nodeChildren	= nodeInfo ["Children"];
		cCount			= Abs(nodeChildren);
		
		if (cCount)
		{
			localMatrices = {};
			
			nodeName = nodeInfo["Name"];
			
			/*fprintf (stdout, "\n@----------\nNode ", nodeName, "\n");*/
			for (s1 = 0; s1<cCount; s1=s1+1)
			{
				childIndex = nodeChildren[s1];
				childIndex = treeAVL	 [childIndex];
				childName  = childIndex  ["Name"];
				localMatrices[s1] = assignmentMatrices[childName];
				/*fprintf (stdout, "\nChild ", childName, localMatrices[s1], "\n");*/
			}
			
			twoWay = {2,1};
			for (s2 = 0; s2 < 2; s2 = s2+1)
			{
				lc = 0;
				for (s3 = 0; s3<cCount; s3=s3+1)
				{
					aMx = localMatrices[s3];
					lc = lc + aMx[1][s2];
				}
				twoWay[s2] = lc;
			}
			
			/*fprintf (stdout, ">> TWO-WAY ", twoWay, "\n");*/

			if (nodeInfo["Parent"])
			{
				aMx = {2,2};
				if (twoWay[1]+1<twoWay[0])
				{
					aMx[0][0]=1;
					aMx[1][0]=twoWay[1]+1;
				}
				else
				{
					aMx[0][0]=0;
					aMx[1][0]=twoWay[0];		
				}
				
				if (twoWay[0]+1<twoWay[1])
				{
					aMx[0][1]=0;
					aMx[1][1]=twoWay[0]+1;
				}
				else
				{
					aMx[0][1]=1;
					aMx[1][1]=twoWay[1];		
				}
				
				/*fprintf (stdout, "NODE MATRIX", aMx, "\n");*/
				assignmentMatrices [nodeName] = aMx;
			}
			else
			{
				totalCost = Min(twoWay[0], twoWay[1]);
				nodeState [nodeName] = (twoWay[1]<twoWay[0]);
			}
		}
	}
	return totalCost;
}

/*--------------------------------------------------- */

function spoolBranchClass (classTag, classAVL)
{
	k2 = Abs(classAVL);
	if (k2)
	{
		classSpool = "";
		classSpool * 256;
		classSpool * ("_branchClasses[\"" + classTag + "\"]={{");
		classSpool * ("\""+classAVL[0]+"\"");
		for (k=1; k<k2;k=k+1)
		{
			classSpool * (",\""+classAVL[k]+"\"");
		}
		classSpool * ("}};\n");
		classSpool * 0;
		fprintf (LAST_FILE_PATH, classSpool);
	}
	return k2;
}
/*--------------------------------------------------- */


SetDialogPrompt ("Load a Newick tree file");
ACCEPT_ROOTED_TREES = 1;
fscanf (PROMPT_FOR_FILE,"Tree",givenTree);
ACCEPT_ROOTED_TREES = 0;
leafCount = TipCount (givenTree);

fprintf (stdout, "Read tree: ", Format (givenTree,0,0),"\n");
goOn = 1;

while (goOn)
{
	leafAllocs = {1,leafCount};

	fprintf (stdout,"\nEnter a reg exp to separate the sequences into two clades:");
	fscanf  (stdin,"String",theRegExp);

	cladeA = {};
	cladeB = {};

	st1 = "";
	st2 = "";


	for (specIndex = 0; specIndex < leafCount; specIndex = specIndex + 1)
	{
		specName = TipName (givenTree, specIndex);
		specMatch = specName $ theRegExp;
		
		if (specMatch[0]>=0)
		{
			cladeA [specName] = 1;
			if (Abs(st1))
			{
				st1 = st1 + "," + specName;
			}
			else
			{
				st1 = specName;
			}
			
		}
		else
		{
			leafAllocs[specIndex] = 1;
			cladeB [specName] = 1;
			
			if (Abs(st2))
			{
				st2 = st2 + "," + specName;
			}
			else
			{
				st2 = specName;
			}
		}
	}
	
	clASize = Abs (cladeA);
	clBSize = Abs (cladeB);
	
	if (clASize == 0 || clBSize == 0 || clASize + clBSize < ds.species)
	{
		fprintf (stdout, "\nERROR: invalid sequence partitionings - one of the clades is empty or there were duplicate sequence names\n");
		return 0;
	}

	fprintf (stdout, "\nSet 1 (TYPE 1) includes ", clASize," sequences:\n");
	cladeKeys = Rows (cladeA);
	for (specIndex = 0; specIndex < clASize; specIndex = specIndex + 1)
	{
		fprintf (stdout, "\t", cladeKeys[specIndex],"\n");
	}

	fprintf (stdout, "\nSet 2 (TYPE 2) includes ", clBSize," sequences:\n");
	cladeKeys = Rows (cladeB);
	for (specIndex = 0; specIndex < clBSize; specIndex = specIndex + 1)
	{
		fprintf (stdout, "\t", cladeKeys[specIndex],"\n");
	}
	
	fprintf (stdout, "\nIs this partitioning correct (y/n)");
	fscanf (stdin, "String", goOn);
	goOn = (goOn[0] == "n" || goOn[0] == "N");
}
f_1 = clASize/leafCount;
f_2 = clBSize/leafCount;

fprintf (stdout, "Please enter a descriptive name for TYPE 1 sequences:");
fscanf	(stdin, "String", class1Name);

fprintf (stdout, "Please enter a descriptive name for TYPE 2 sequences:");
fscanf	(stdin, "String", class2Name);

fprintf (stdout, "\nProportion of ",class1Name," sequences: ", f_1,
				 "\nProportion of ",class2Name," sequences: ", f_2,"\n");




treeAVL 	= givenTree^0;
treeAVL2 	= givenTree^1;
treeSize    = Abs (treeAVL);

observedEvents = computeMigrationEvents (leafAllocs);

fprintf (stdout, "\nInferred ", observedEvents, " migration events\n");

transitionNodes12 = {};
transitionNodes21 = {};
type1Nodes		  = {};
type2Nodes		  = {};

for (node = 2; node < treeSize; node = node + 1)
{
	nodeInfo 		= treeAVL2[node];
	nodeChildren	= nodeInfo ["Children"];
	parentNode = nodeInfo["Parent"];
	parentNode = treeAVL2[parentNode];
	parentNode = parentNode["Name"];
	parentNode = nodeState[parentNode];
	nodeName   = nodeInfo["Name"];
	myState    = assignmentMatrices[nodeName];
	myState    = myState[0][parentNode];
	if (parentNode != myState)
	{
		if (parentNode == 0)
		{
			transitionNodes12[Abs(transitionNodes12)] = nodeName;
		}
		else
		{
			transitionNodes21[Abs(transitionNodes21)] = nodeName;		
		}
	}
	else
	{
		if (parentNode == 0)
		{
			type1Nodes[Abs(type1Nodes)] = nodeName;
		}
		else
		{
			type2Nodes[Abs(type2Nodes)] = nodeName;		
		}
	}
	nodeState [nodeName] = myState;
}

if (Abs(transitionNodes12)+Abs(transitionNodes21))
{
	fprintf (stdout, "\nThe following branches have migration events:\n\n");
	if (Abs(transitionNodes12))
	{
		fprintf (stdout, "\n",class1Name," --> ",class2Name,"\n");
		for (node = 0; node < Abs(transitionNodes12); node = node + 1)
		{
			fprintf (stdout, transitionNodes12[node], "\n");
		}
	}
	
	if (Abs(transitionNodes21))
	{
		fprintf (stdout, "\n",class2Name," --> ",class1Name,"\n");
		for (node = 0; node < Abs(transitionNodes21); node = node + 1)
		{
			fprintf (stdout, transitionNodes21[node], "\n");
		}
	}
}

SetDialogPrompt ("Write a tree with branch partitions to:");
fprintf (PROMPT_FOR_FILE,CLEAR_FILE,"Tree theTree=", Format (givenTree,1,0),";\n_branchClasses={};\n");

spoolBranchClass(class1Name,type1Nodes);
spoolBranchClass(class2Name,type2Nodes);
spoolBranchClass(class1Name+"->"+class2Name,transitionNodes12);
spoolBranchClass(class2Name+"->"+class1Name,transitionNodes21);


ChoiceList (resample,"Permutation Test",1,SKIP_NONE,
					 "Skip","Do not perform a permutation test.",
				     "Most certainly","Randomly allocate sequences into classes and tabulate the distribution of migration events.");
				 

if (resample < 1)
{
	return 0;
}

sampleCount = 0;
totalPossible = factorial(leafCount)/factorial(clASize)/factorial (clBSize);
while (sampleCount < 1)
{
	fprintf (stdout, "\nHow many permutations would you like to perform (total of ", totalPossible, " unique permutations available)?");
	fscanf (stdin,"Number",sampleCount);
}

histogram = {};

step = sampleCount$20;
for (sampleID = 0; sampleID < sampleCount; sampleID = sampleID + 1)
{
	aSample = Random(leafAllocs,0);
	k = computeMigrationEvents(aSample);
	histogram[k] = histogram[k] + 1;

	if ((sampleID+1)%step == 0)
	{
		fprintf (stdout, Format ((sampleID+1)*100/sampleCount, 6, 2), "% done\n");
	}
}

k2 = Abs (histogram);
countMatrix = {k2,3};
s1 = Rows (histogram);
for (sampleID = 0; sampleID < k2; sampleID = sampleID+1)
{
	s2 = 0+s1[sampleID];
	countMatrix[sampleID][0] = s2;
	countMatrix[sampleID][1] = histogram[s2];
}

countMatrix = countMatrix % 0;
pVal = 0;

for (sampleID = 0; sampleID < k2; sampleID = sampleID+1)
{
	if (countMatrix[sampleID][0] <= observedEvents)
	{
		pVal = pVal+countMatrix[sampleID][1];
	}
	if (sampleID)
	{
		countMatrix[sampleID][2] = countMatrix[sampleID][1]/sampleCount + countMatrix[sampleID-1][2];
	}
	else
	{
		countMatrix[sampleID][2] = countMatrix[sampleID][1]/sampleCount;	
	}
}

fprintf (stdout, "\n\nProb{as many or fewer migration events by chance} = ", pVal/sampleCount,"\n\n");

labels = {{"Events","Count","Cumulative Weight"}};

OpenWindow (CHARTWINDOW,{{"Simulated distribution of migration events"}
		{"labels"}
		{"countMatrix"}
		{"Bar Chart"}
		{"Events"}
		{"Count"}
		{"Number of Events"}
		{""}
		{"Counts"}
		{""}
		{""}
		{"-1;-1"}
		{"10;1.309;0.785398"}
		{"Times:14:0;Times:12:0;Times:14:1"}
		{"16777215;16777215;16512;11776947;0;16777215;16711680;11842740;13158600;14474460;0;3947580;79;16744448;16777215;2984993;9199669;7018159;1460610;16748822;11184810;14173291;14173291"}
		{"16"}
		},
		"SCREEN_WIDTH-100;SCREEN_HEIGHT-100;50;50");
