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
		aMx = {2,kindCount};
		s2  = assignmentVector[specIndex];
		for (k=0; k<kindCount; k=k+1)
		{
			aMx[0][k] = s2;
			aMx[1][k] = 1-(k==s2);	
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
			
			/*
			fprintf (stdout, "\n@----------\nNode ", nodeName, "\n");
			*/
			
			for (s1 = 0; s1<cCount; s1=s1+1)
			{
				childIndex = nodeChildren[s1];
				childIndex = treeAVL	 [childIndex];
				childName  = childIndex  ["Name"];
				localMatrices[s1] = assignmentMatrices[childName];
				/*
				fprintf (stdout, "\nChild ", childName, localMatrices[s1], "\n");
				*/
			}
			
			twoWay = {kindCount,1};
			
			for (s2 = 0; s2 < kindCount; s2 = s2+1)
			{
				lc = 0;
				for (s3 = 0; s3<cCount; s3=s3+1)
				{
					aMx = localMatrices[s3];
					lc  = lc + aMx[1][s2];
				}
				twoWay[s2] = lc;
			}
			
			/*
			fprintf (stdout, ">> TWO-WAY ", twoWay, "\n");
			*/
			
			if (nodeInfo["Parent"])
			{
				aMx = {2,kindCount};
				
				for (s2 = 0; s2 < kindCount; s2 = s2+1)
				{
					minV = 1e100;
					minI = 0;
					
					for (s3 = 0; s3 < kindCount; s3 = s3+1)
					{
						aCost = (s3!=s2) + twoWay[s3];
						if (minV > aCost)
						{
							minV  = aCost;
							minI  = s3;
						}	
					}
					
					aMx[0][s2] = minI;
					aMx[1][s2] = minV;
				}
				
				/*
				fprintf (stdout, "NODE MATRIX", aMx, "\n");
				*/
				assignmentMatrices [nodeName] = aMx;
			}
			else
			{
				totalCost = 1e100;
				nodeState [nodeName] = 0;
				
				
				for (s2 = 0; s2 < kindCount; s2 = s2+1)
				{
					if (twoWay[s2] < totalCost)
					{
						totalCost = twoWay[s2];
						nodeState [nodeName] = s2;
					}	
				}
			}
		}
	}
	return totalCost;
}

/*--------------------------------------------------- */

function spoolBranchClass (classTag, classAVL, doPrint)
{
	k2 = Abs(classAVL);
	if (k2)
	{
		classSpool = "";
		classSpool * 256;
		classSpool * ("_branchClasses[\"" + classTag + "\"]={{");
		classSpool * ("\""+classAVL[0]+"\"");
		if (doPrint)
		{
			fprintf (stdout, classAVL[0], "\n");
		}
		for (k=1; k<k2;k=k+1)
		{
			classSpool * (",\""+classAVL[k]+"\"");
			if (doPrint)
			{
				fprintf (stdout, classAVL[k], "\n");
			}
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

kindCount = 0;
while (kindCount < 2)
{
	fprintf (stdout, "How many sequence types: (>=2):");
	fscanf  (stdin, "Number", kindCount);
	kindCount = kindCount $ 1;
}

goOn  = 1;


while (goOn)
{
	defClades 	 = 0;
	clades  	 = {};
	strings 	 = {};
	matchedSoFar = {leafCount,1};
	leafAllocs 	 = {1,leafCount};

	while (defClades < kindCount)
	{

		st     = "";
		aClade = {};
		
		if (defClades < kindCount-1)
		{
			fprintf (stdout,"\nEnter a reg exp used to define clade ",defClades+1,":");
			fscanf  (stdin,"String",theRegExp);

			for (specIndex = 0; specIndex < leafCount; specIndex = specIndex + 1)
			{
				specName = TipName (givenTree, specIndex);
				specMatch = specName $ theRegExp;
				
				if (specMatch[0]>=0 && matchedSoFar[specIndex] == 0)
				{
					aClade [specName] = 1;
					if (Abs(st))
					{
						st = st + "," + specName;
					}
					else
					{
						st = specName;
					}
					
					matchedSoFar [specIndex] = 1;
					leafAllocs   [specIndex] = defClades;
				}
			}
		}
		else
		{
			for (specIndex = 0; specIndex < leafCount; specIndex = specIndex + 1)
			{
				if (matchedSoFar[specIndex] == 0)
				{
					specName = TipName (givenTree, specIndex);
					aClade [specName] = 1;
					if (Abs(st))
					{
						st = st + "," + specName;
					}
					else
					{
						st = specName;
					}
					matchedSoFar [specIndex] = 1;
					leafAllocs   [specIndex] = kindCount-1;
				}
			}		
		}
		
		if (Abs(aClade) == 0)
		{
			fprintf (stdout, "ERROR: an empty clade for reg-exp ", goOn, "\n");
			defClades = kindCount;
			break;
		}
		else
		{
			fprintf (stdout, "Matched: ",st,"\n");	
		}
		strings[Abs(strings)] = st;
		clades [Abs(clades) ] = aClade;
		defClades = defClades + 1;
	}

	if (Abs(clades) == kindCount)
	{
		for (k=0; k<kindCount; k=k+1)
		{
			aClade = clades[k];
			clASize = Abs(aClade);
			fprintf (stdout, "\nSet ",k+1," (TYPE ",k+1,") includes ", clASize," sequences:\n");
			cladeKeys = Rows (aClade);
			for (specIndex = 0; specIndex < clASize; specIndex = specIndex + 1)
			{
				fprintf (stdout, "\t", cladeKeys[specIndex],"\n");
			}
		}
		
		fprintf (stdout, "\nIs this partitioning correct (y/n)");
		fscanf (stdin, "String", goOn);
		goOn = (goOn[0] == "n" || goOn[0] == "N");
	}
	else
	{
		goOn = 1;
	}
}

descriptives = {};
for (k=0; k<kindCount; k=k+1)
{
	fprintf (stdout, "Please enter a descriptive name for TYPE ",k+1," sequences:");
	fscanf	(stdin, "String", className);
	descriptives [k] = className;

	fprintf (stdout, "\nProportion of ",className," sequences: ", Abs(clades[k])/leafCount, "\n");
}

treeAVL 	= givenTree^0;
treeAVL2 	= givenTree^1;
treeSize    = Abs (treeAVL);

observedEvents = computeMigrationEvents (leafAllocs);

fprintf (stdout, "\nInferred ", observedEvents, " migration events\n");

SetDialogPrompt ("Write a tree with branch partitions to:");
fprintf (PROMPT_FOR_FILE,CLEAR_FILE,"Tree theTree=", Format (givenTree,1,0),";\n_branchClasses={};\n");

transitions = {};
sameState   = {};

for (node = 2; node < treeSize; node = node + 1)
{
	nodeInfo 		= treeAVL2[node];
	nodeChildren	= nodeInfo ["Children"];
	parentNode 		= nodeInfo["Parent"];
	parentNode 		= treeAVL2[parentNode];
	parentNode 		= parentNode["Name"];
	parentNode 		= nodeState[parentNode];
	nodeName   		= nodeInfo["Name"];
	myState    		= assignmentMatrices[nodeName];
	myState    		= myState[0][parentNode];
	
	if (parentNode != myState)
	{
		aKey = descriptives[parentNode] + " --> " + descriptives[myState];
		aList = transitions[aKey];
		if (Abs(aList) == 0)
		{
			aList = {};
		}
		aList[Abs(aList)] = nodeName;
		transitions[aKey] = aList;
	}
	else
	{
		aKey = descriptives[parentNode];
		aList = sameState[aKey];
		if (Abs(aList) == 0)
		{
			aList = {};
		}
		aList[Abs(aList)] = nodeName;
		sameState[aKey] = aList;
	}
	nodeState [nodeName] = myState;
}

returnAVL = {"MIGRATIONS": Abs(transitions)};

if (Abs(transitions))
{
	fprintf (stdout, "\nThe following branches have migration events:\n\n");
	
	keys = Rows (transitions);
	
	for (cc = 0; cc < Columns (keys); cc = cc + 1)
	{
		aKey = keys[cc];
		fprintf (stdout, "\n", aKey, ":\n");
		aList = transitions[aKey];
		spoolBranchClass(aKey,aList,1);
	}
}

if (Abs(sameState))
{
	keys = Rows (sameState);
	
	for (cc = 0; cc < Columns (keys); cc = cc + 1)
	{
		aKey = keys[cc];
		aList = sameState[aKey];
		spoolBranchClass(aKey,aList,0);
	}
}


ChoiceList (resample,"Permutation Test",1,SKIP_NONE,
					 "Skip","Do not perform a permutation test.",
				     "Most certainly","Randomly allocate sequences into classes and tabulate the distribution of migration events.");
				 

if (resample)
{
	sampleCount = 0;
	totalPossible = factorial(leafCount);
	
	for (k=0; k<kindCount; k=k+1)
	{
		totalPossible = totalPossible/factorial(Abs(clades[k]));
	}
	
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
			
	returnAVL["p"] = pVal/sampleCount;
}

return returnAVL;