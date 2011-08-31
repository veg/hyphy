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

function computeGCC (matrix1&, matrix2&)
{
	dim = Rows (matrix1);
	offDiag = ({dim,dim})["_MATRIX_ELEMENT_ROW_!=_MATRIX_ELEMENT_COLUMN_"];
	rowM	= {1,dim}["1"];
	colM	= {dim,1}["1"];
	
	sum1   = (rowM*(matrix1$offDiag)*colM)[0]/(dim*(dim-1));
	sum2   = (rowM*(matrix2$offDiag)*colM)[0]/(dim*(dim-1));
	
	tmx1   = (matrix1 - offDiag*sum1);
	tmx2   = (matrix2 - offDiag*sum2);
	
	numerator 	= (rowM*(tmx1$tmx2)*colM)[0];
	denominator	= Sqrt((rowM*(tmx1$tmx1)*colM)[0] * (rowM*(tmx2$tmx2)*colM)[0]);
	return numerator/denominator;
}

/*--------------------------------------------------- */

function computeTreeDistanceMetrics (treeID)
{
	ExecuteCommands ("_treeAVL = "+ treeID + " ^ 0;_leafCount=TipCount("+treeID+");"); 
	
	branchDistances 	= {_leafCount,_leafCount};
	geneticDistances	= {_leafCount,_leafCount};
	
	_leafIndexer		= 0;
	
	for (_k=1; _k<Abs(_treeAVL)-1; _k=_k+1)
	{
		aNode 	    = _treeAVL[_k];
		aNodeName   = aNode["Name"];
		parentIndex = aNode["Parent"];
		bL			= aNode["Length"]; 

		if (Rows((_treeAVL[parentIndex])["Split"]) == 0)
		{
			(_treeAVL[parentIndex])["Split"] = {_leafCount,1};
		} 

		if (Abs(aNode["Children"]) == 0)
		{
			for (_k2 = 0; _k2 < _leafCount; _k2 = _k2 + 1)
			{
				if (_leafIndexer != _k2)
				{
					branchDistances [_k2][_leafIndexer]  = branchDistances [_k2][_leafIndexer] + 1;
					branchDistances [_leafIndexer][_k2]  = branchDistances [_k2][_leafIndexer];
					geneticDistances [_k2][_leafIndexer] = geneticDistances [_k2][_leafIndexer] + bL;
					geneticDistances [_leafIndexer][_k2] = geneticDistances [_k2][_leafIndexer];
				}
			}
			((_treeAVL[parentIndex])["Split"])[_leafIndexer] = 1;
			_leafIndexer = _leafIndexer + 1;
		}
		else
		{
			_inL 	= {};
			_outL	= {};
			
			for (_k2 = 0; _k2 < _leafCount; _k2 = _k2 + 1)
			{
				if ((aNode["Split"])[_k2])
				{
					_inL[Abs(_inL)] = _k2;
				}
				else
				{
					_outL[Abs(_outL)] = _k2;
				}
			}
				
			for (_k2 = Abs(_outL)-1; _k2 >=0 ; _k2 = _k2 - 1)
			{
				_idx1 = _outL[_k2];
				for (_k3 = Abs(_inL)-1; _k3 >=0 ; _k3 = _k3 - 1)
				{
					_idx2 = _inL[_k3];
					branchDistances [_idx1][_idx2]   = branchDistances [_idx1][_idx2] + 1;
					branchDistances [_idx2][_idx1]   = branchDistances [_idx2][_idx1] + 1;
					geneticDistances [_idx1][_idx2]  = geneticDistances [_idx1][_idx2] + bL;
					geneticDistances [_idx2][_idx1]  = geneticDistances [_idx2][_idx1] + bL;
				}
			}
			
			(_treeAVL[parentIndex])["Split"] = (_treeAVL[parentIndex])["Split"] + aNode["Split"];
		}
		
	}

	_returnValue = {};
	_returnValue ["GD"] = geneticDistances;
	_returnValue ["BD"] = branchDistances;
	
	return 		 _returnValue;
	
}



/*--------------------------------------------------- */

fprintf (stdout, "\n+-------------------------------------+\n",
				   "| r and r_b compartmentalization test |\n",
				   "| implementing the procedure given by |\n",
				   "|       D.E. Critchlow et al. in      |\n",
				   "| Mathematical and Computer Modeling  |\n",
				   "|          32 (2000) 69-81            |\n",
				   "+-------------------------------------+\n");
				   
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

goOn  		 = 1;

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

splitsMatrix   	   = {leafCount,leafCount};

for (k=0; k<leafCount; k=k+1)
{
	for (k2=0; k2<leafCount; k2=k2+1)
	{
		if (leafAllocs[k] != leafAllocs[k2])
		{
			splitsMatrix[k][k2] = 1;
			splitsMatrix[k2][k] = 1;
		}
	}
}

distanceMatrix 	   = computeTreeDistanceMetrics ("givenTree");
bdmx			   = distanceMatrix["BD"];
gdmx			   = distanceMatrix["GD"];

ccb				   = computeGCC ("bdmx","splitsMatrix");
ccg				   = computeGCC ("gdmx","splitsMatrix");

fprintf			   (stdout, "\nCorrelation coefficients:\n",
						    "\n\tBranch counts (r_b) : ", ccb, 
						    "\n\tPath lengths  (r)   : ", ccg, "\n");

ChoiceList (resample,"Permutation Test",1,SKIP_NONE,
					 "Skip","Do not perform a permutation test.",
				     "Most certainly","Randomly allocate sequences into classes and tabulate the distribution of migration events.");
				 

if (resample < 1)
{
	return 0;
}

sampleCount = 0;
totalPossible = factorial(leafCount);

for (k=0; k<kindCount; k=k+1)
{
	totalPossible = totalPossible/factorial(Abs(clades[k]));
}

while (sampleCount < 1)
{
	fprintf (stdout, "\nHow many permutations would you like to perform?");
	fscanf (stdin,"Number",sampleCount);
}


distribution_R  = {sampleCount,1};
distribution_RB = {sampleCount,1};
pvR 			= 0;
pvRB			= 0;

step = sampleCount$10;

for (sampleID = 0; sampleID < sampleCount; sampleID = sampleID + 1)
{
	aSample 			   = Random(leafAllocs,0);
	splitsMatrixR   	   = {leafCount,leafCount};

	for (k=0; k<leafCount; k=k+1)
	{
		for (k2=0; k2<leafCount; k2=k2+1)
		{
			if (aSample[k] != aSample[k2])
			{
				splitsMatrixR[k][k2] = 1;
				splitsMatrixR[k2][k] = 1;
			}
		}
	}

	ccbR				   = computeGCC ("bdmx","splitsMatrixR");
	ccgR				   = computeGCC ("gdmx","splitsMatrixR");

	distribution_R  [sampleID] = ccgR;
	distribution_RB [sampleID] = ccgB;
	
	pvR  = pvR + (ccgR>=ccg);
	pvRB = pvRB + (ccbR>=ccb);
	
	if ((sampleID+1)%step == 0)
	{
		fprintf (stdout, Format ((sampleID+1)*100/sampleCount, 6, 2), "% done\n");
	}
}

fprintf (stdout, "\n\nProb{r_b random >= r_b observed} < ", (pvR+1)/(1+sampleCount),"");
fprintf (stdout, "\nProb{r random   >= r observed  } < ", (pvRB+1)/(1+sampleCount),"\n\n");
