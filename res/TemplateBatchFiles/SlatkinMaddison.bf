/*--------------------------------------------------- */

function factorial (x)
{
	res = 1;
	for (x2 = 2; x2 <= x; x2 += 1) {
		res = res * x2;
	}
	return res;
}

/*--------------------------------------------------- */

function computeMigrationEvents (assignmentVector) {
    /* 
        SLKP 20180125
        updated this for cases when internal node names are not unique, e.g., when
        bootstrap values are supplied
    */
	
	
	leaf_count = 0;

	for (node = 1; node < treeSize; node += 1) {
		nodeChildren	= (treeAVL[node]) ["Children"];
		cCount			= Abs(nodeChildren);

        /*
             for each node populate a 2xN (N = number of character states)
             element (1,N) stores the cost of the subtree starting at this node, assuming parent has state N
             element (0,N) stores the assignment at this node, which realizes the score
         
        */

		if (cCount == 0) { // leaf
		
		    s = assignmentVector[leaf_count];
		    (treeAVL[node])["_PARSIMONY_"] = {2,kindCount} ["(_MATRIX_ELEMENT_ROW_==0)*s+(_MATRIX_ELEMENT_ROW_==1)*(1-(_MATRIX_ELEMENT_COLUMN_==s))"];
		    leaf_count += 1;
		    
		} else { // internal node

			subtreeCost = {kindCount,1};
			
			// total cost of setting this node to state s2 

			for (s2 = 0; s2 < kindCount; s2 += 1) {
				lc = 0;
				for (s3 = 0; s3<cCount; s3 += 1) {
				    lc  += ((treeAVL[nodeChildren[s3]])["_PARSIMONY_"])[1][s2];
				}
				subtreeCost[s2] = lc;
			}


			if ((treeAVL[node])["Parent"]) {
			    // not the root
			    
				aMx = {2,kindCount};

				for (s2 = 0; s2 < kindCount; s2 += 1) {
                    score_matrix = Min ({1,kindCount} ["(_MATRIX_ELEMENT_COLUMN_!=s2)+subtreeCost[_MATRIX_ELEMENT_COLUMN_]"], 1);
					aMx[0][s2] = score_matrix[1];
					aMx[1][s2] = score_matrix[0];
				}
				
				(treeAVL[node])["_PARSIMONY_"] = aMx;
				
				
			} else {
				subtreeCost = Min (subtreeCost, 1);
				(treeAVL[node])["_PARSIMONY_"] = subtreeCost[1];
			}
		}
	}
	
	
	return subtreeCost[0];
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
fscanf (PROMPT_FOR_FILE,"String",treeString);

Topology givenTree = treeString;

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


while (goOn) {
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

			for (specIndex = 0; specIndex < leafCount; specIndex = specIndex + 1) {
				specName = TipName (givenTree, specIndex);
				specMatch = specName $ theRegExp;

				if (specMatch[0]>=0 && matchedSoFar[specIndex] == 0) {
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
					leafAllocs   [specIndex] = Random (0,1) > 0.05;//defClades;
				}
			}
		}
		else {
			for (specIndex = 0; specIndex < leafCount; specIndex = specIndex + 1)
			{
				if (matchedSoFar[specIndex] == 0) {
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
					leafAllocs   [specIndex] = Random (0,1) > 0.05; //kindCount-1;
				}
			}
		}

		if (Abs(aClade) == 0) {
			fprintf (stdout, "ERROR: an empty clade for reg-exp ", goOn, "\n");
			defClades = kindCount;
			break;
		}
		else
		{
			fprintf (stdout, "Matched: ",st,"\n");
		}
		strings + st;
		clades + aClade;
		defClades += 1;
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
for (k=0; k<kindCount; k=k+1) {
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

for (node = treeSize - 2; node > 0 ; node = node - 1) {
    
    
    parentState = (treeAVL[((treeAVL[node])["Parent"])])["_PARSIMONY_"];
    myState     = ((treeAVL[node])["_PARSIMONY_"])[0][parentState];
    
    (treeAVL[node])["_PARSIMONY_"] = myState;
    
	if (parentState != myState) {
		aKey = descriptives[parentState] + " --> " + descriptives[myState];
		if (Abs(transitions[aKey]) == 0) {
			transitions[aKey] = {};
		}
		transitions[aKey] + ((treeAVL[node])["Name"]);
	}
	else {
		aKey = descriptives[parentNode];
		if (Abs(sameState[aKey]) == 0) {
			sameState[aKey] = {};
		}
		sameState[aKey] + ((treeAVL[node])["Name"]);
	}
}

returnAVL = {"MIGRATIONS": Abs(transitions)};

if (Abs(transitions))
{
	fprintf (stdout, "\nThe following branches have migration events:\n\n");

	keys = Rows (transitions);

	for (cc = 0; cc < Columns (keys); cc += 1) {
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

	for (k=0; k<kindCount; k += 1) {
		totalPossible = totalPossible/factorial(Abs(clades[k]));
	}

	while (sampleCount < 1)
	{
		fprintf (stdout, "\nHow many permutations would you like to perform (total of ", totalPossible, " unique permutations available)?");
		fscanf (stdin,"Number",sampleCount);
	}

	histogram = {};

	step = sampleCount$20;
	for (sampleID = 0; sampleID < sampleCount; sampleID += 1) {
		k = computeMigrationEvents(Random(leafAllocs,0));
		fprintf (stdout, "Replicate ", (sampleID+1), " has ", k, " migration events\n");
		histogram[k] = histogram[k] + 1;
	}

	k2 = Abs (histogram);
	countMatrix = {k2,3};
	s1 = Rows (histogram);
	for (sampleID = 0; sampleID < k2; sampleID = sampleID+1)
	{
		s2 = +s1[sampleID];
		countMatrix[sampleID][0] = s2;
		countMatrix[sampleID][1] = histogram[s2];
	}

	countMatrix = countMatrix % 0;
	pVal = 0;

	for (sampleID = 0; sampleID < k2; sampleID = sampleID+1)
	{
		if (countMatrix[sampleID][0] <= observedEvents) {
			pVal = pVal+countMatrix[sampleID][1];
		}
		if (sampleID) {
			countMatrix[sampleID][2] = countMatrix[sampleID][1]/sampleCount + countMatrix[sampleID-1][2];
		}
		else {
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
