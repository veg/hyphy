RequireVersion ("0.9920060601");
VERBOSITY_LEVEL = -1;

/*************************************************************************************************/

function compute1mf 	  (aNum)
{
	maxV = Max (maxV,aNum);
	sumV = sumV + aNum;
	return aNum;
}

/*************************************************************************************************/

function computeSimmondsD (leafAssignmentVector, mapVec, catCount)
{
	_d           = 0;
	_currentLeaf = 0;
	for (_nodeIdx = 1; _nodeIdx < Abs (treeAVL); _nodeIdx = _nodeIdx + 1)
	{
		_currentNode 	= treeAVL[_nodeIdx];
		_childrenCount  = Abs (_currentNode["Children"]);
		_parentNode		= _currentNode ["Parent"];
		if (_parentNode)
		{
			if (Columns((treeAVL[_parentNode])["Counts"])==0)
			{
				(treeAVL[_parentNode])["Counts"] = {1,catCount};
			}
			if (_childrenCount) /* internal node */
			{
				(treeAVL[_parentNode])["N"] = (treeAVL[_parentNode])["N"]+_currentNode["N"];
				_currentCounts = _currentNode["Counts"];
				(treeAVL[_parentNode])["Counts"] = (treeAVL[_parentNode])["Counts"] + 
												  _currentCounts;
												  
				maxV = 0;
				sumV = 0;
				_currentCounts ["compute1mf(_MATRIX_ELEMENT_VALUE_)"];
				/*fprintf (stdout, 1-maxV/sumV, "\n");*/
				_d = _d + (1-maxV/sumV)*2^(1-_currentNode["N"]);
			
			}
			else /* terminal node */
			{
				(treeAVL[_parentNode])["N"] = (treeAVL[_parentNode])["N"]+1;
				((treeAVL[_parentNode])["Counts"])[leafAssignmentVector[mapVec[_currentLeaf]]] =
				((treeAVL[_parentNode])["Counts"])[leafAssignmentVector[mapVec[_currentLeaf]]] + 1; 
				_currentLeaf = _currentLeaf + 1;
			}
		}
	}
	
	return _d;
}

/*************************************************************************************************/

SetDialogPrompt ("Load a nucleotide sequence file:");

DataSet ds = ReadDataFile (PROMPT_FOR_FILE);
leafCount  = ds.species;

fprintf (stdout, "Read file: ", ds,"\n");

DataSetFilter filteredData = CreateFilter (ds,1);

goOn = 1;

kindCount = 0;
while (kindCount < 2)
{
	fprintf (stdout, "How many sequence types: (>=2):");
	fscanf  (stdin, "Number", kindCount);
	kindCount = kindCount $ 1;
}

goOn  = 1;

seqNamesInFile = {};
choiceMatrix   = {leafCount,2};

for (specIndex = 0; specIndex < leafCount; specIndex = specIndex + 1)
{
	GetString (specName, ds, specIndex);
	seqNamesInFile[specName] = specIndex;
	choiceMatrix[specIndex][0] = specName;
	choiceMatrix[specIndex][1] = "Use "+specName+" as the outgroup";
}

ChoiceList (rerootAt, "Outgroup", 1, SKIP_NONE, choiceMatrix);

if (rerootAt < 0)
{
	return 0;
}

choiceMatrix = choiceMatrix[rerootAt][0];

fprintf (stdout, "Using ", choiceMatrix, " as outgroup\n");

ACCEPT_ROOTED_TREES = 1;

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
				GetString (specName, ds, specIndex);
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
					GetString (specName, ds, specIndex);
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

/*descriptives = {};*/
for (k=0; k<kindCount; k=k+1)
{
	/*fprintf (stdout, "Please enter a descriptive name for TYPE ",k+1," sequences:");
	fscanf	(stdin, "String", className);
	descriptives [k] = className;*/

	fprintf (stdout, "\nProportion of sequences in group ",k,": ", Abs(clades[k])/leafCount, "\n");
}

k = HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "Utility" + DIRECTORY_SEPARATOR + "NJ.bf";
ExecuteCommands ("#include \""+k+"\";");
ts	 = InferTreeTopology (0);
Tree givenTree = ts;
ts	 = RerootTree (givenTree, choiceMatrix);
Tree givenTree = ts;
/*fprintf (stdout, givenTree, "\n");*/


fprintf		(stdout, "How many relabelings per sample (default 10):?");
fscanf		(stdin, "String",s);
shuffleIt	= 0 + s;
if (shuffleIt < 1)
{
	shuffleIt = 10;
}

fprintf		(stdout, "How many tree bootstrap samples (default 100):?");
fscanf		(stdin, "String",s);
treeIt	= 0 + s;
if (treeIt < 1)
{
	treeIt = 100;
}

fprintf		(stdout, "Proportion of reshufflings less associated than the sample needed for significance (default 2/3)?");
fscanf		(stdin, "String",s);
propSig 	= 0 + s;
if (propSig <= 0 || propSig>=1)
{
	propSig = 2/3;
}

fprintf (stdout, "Using ", treeIt, " tree bootstraps and ", shuffleIt, " relabelings per sample with significance called at ", propSig, "\n");

treeAVL 	= givenTree^0;
treeAVL2 	= givenTree^1;

baseD		= runATreeSample (0);
fprintf   	(stdout, "\nBaseline d = ", baseD[0], "\n");

totalRes    = {treeIt,3};

totalRes[0][0] = baseD[0];
totalRes[0][1] = baseD[1];
totalRes[0][2] = baseD[2];

fprintf		(stdout, "Running tree simulations...\n");
meanO 		= baseD[0];
meanS		= baseD[1];
sigB		= (propSig < baseD[2]);

for (it = 0; it < treeIt-1; it = it + 1)
{
	DataSetFilter filteredData = Bootstrap (ds,1);
	ts	 = InferTreeTopology (0);
	Tree givenTree = ts;
	ts	 = RerootTree (givenTree, choiceMatrix);
	Tree givenTree = ts;
	treeAVL 	= givenTree^0;
	treeAVL2 	= givenTree^1;
	simD		= runATreeSample (0);
	totalRes[it+1][0] = simD[0];
	totalRes[it+1][1] = simD[1];
	totalRes[it+1][2] = simD[2];
	meanO = meanO + simD[0];
	meanS = meanS + simD[1];
	sigB  = sigB + (propSig < simD[2]);
}

fprintf (stdout, "\n\nAssociation Index: ", meanO/meanS, "\nBootstrap significance :" , sigB, "/", treeIt, "\n");

columnHeaders = {{"Observed","Mean Control","Proportion Control d > Observed d"}};

OpenWindow (CHARTWINDOW,{{"Simmonds AI Summary"}
		{"columnHeaders"}
		{"totalRes"}
		{"Scatterplot"}
		{"Observed"}
		{"Mean Control"}
		{"Observed d"}
		{""}
		{"Control d"}
		{"0"}
		{"_x_"}
		{"-1;-1"}
		{"10;1.309;0.785398"}
		{"Times:12:0;Times:10:0;Times:12:2"}
		{"0;0;16777215;0;0;0;6579300;11842740;13158600;14474460;0;3947580;16777215;15670812;6845928;16771158;2984993;9199669;7018159;1460610;16748822;11184810;14173291"}
		{"16,0,0"}
		},
		"925;642;70;70");
		
ACCEPT_ROOTED_TREES = 0;

/*************************************************************************************************/

function runATreeSample (dummy)
{
	mapVec	    = {1,leafCount}["_MATRIX_ELEMENT_COLUMN_"];
	myLeafAlloc = {1,leafCount};
	for (_k = 0; _k < leafCount; _k = _k + 1)
	{
		leafName = TipName (givenTree,_k);
		myLeafAlloc [_k] = leafAllocs [seqNamesInFile[leafName]];
	}
	baseD 		= computeSimmondsD (myLeafAlloc, mapVec, kindCount);
	gte 		= 0;
	meanRat		= 0;
	for (_k=0; _k<shuffleIt; _k = _k + 1)
	{
		treeAVL 	= givenTree^0;
		rsD			= computeSimmondsD (myLeafAlloc, Random(mapVec,0), kindCount);
		if (rsD > baseD) 
		{
			gte = gte +1;
		}
		meanRat		= meanRat + rsD;
	}	
	outMx = {{baseD, meanRat/shuffleIt, gte/shuffleIt}};
	return outMx;
}
