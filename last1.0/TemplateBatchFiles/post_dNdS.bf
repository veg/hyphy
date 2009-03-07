function valueGrab (varName&)
{
	return varName;
}

T = TipCount (givenTree);

maxLength = 0;

for (h = 0; h<T; h=h+1)
{
	bName = Abs(TipName (givenTree,h));
	if (bName>maxLength)
	{
		maxLength = bName;
	}
}

T = BranchCount (givenTree);

for (h = 0; h<T; h=h+1)
{
	bName = Abs(BranchName (givenTree,h));
	if (bName>maxLength)
	{
		maxLength = bName;
	}
}

fprintf (stdout,"\nLeaf branches: \n\n");

T = TipCount (givenTree);

for (h = 0; h<T; h=h+1)
{
	bName = TipName (givenTree,h);
	fprintf (stdout,bName);
	for (v = Abs(bName); v<maxLength; v=v+1)
	{
		fprintf (stdout," ");
	}
	fprintf (stdout,": ",GetBranchDNDS (bName),"\n");
}

T = BranchCount (givenTree);

if (T>0)
{
	fprintf (stdout,"\nInternal Branches:\n\n");
	for (h = 0; h<T; h=h+1)
	{
		bName = BranchName (givenTree,h);
		fprintf (stdout,bName);
		for (v = Abs(bName); v<maxLength; v=v+1)
		{
			fprintf (stdout," ");
		}
		fprintf (stdout,": ",GetBranchDNDS (bName),"\n");
	}
}

branchNames = BranchName (givenTree,-1);
T = Columns (branchNames);

ExecuteCommands(
"GetInformation (aRateMx, givenTree."+branchNames[0]+");");

/* make syn and non-syn template matrices */

nonStopCount = Columns (aRateMx);

synM    = {nonStopCount,nonStopCount};
nonSynM = {nonStopCount,nonStopCount};

vertOnes = {nonStopCount,1};
horOnes  = {1,nonStopCount};

for (h1 = 0; h1<nonStopCount; h1=h1+1)
{
	vertOnes [h1] = 1;
	horOnes  [h1] = 1;
}

hShift = 0;
for (h1 = 0; h1 < 64; h1=h1+1)
{
	gc1 = _Genetic_Code[h1];
	if (gc1 == 10)
	{
		hShift = hShift+1;
	}
	else
	{
		vShift = hShift;
		for (v1 = h1+1; v1 < 64; v1=v1+1)
		{
			gc2 = _Genetic_Code[v1];
			if (gc2 == 10)
			{
				vShift = vShift + 1;
			}
			else
			{
				if (gc1 == gc2)
				{
					synM [h1-hShift][v1-vShift] = vectorOfFrequencies[h1-hShift];
					synM [v1-vShift][h1-hShift] = vectorOfFrequencies[v1-vShift];
				}
				else
				{
					nonSynM [h1-hShift][v1-vShift] = vectorOfFrequencies[h1-hShift];
					nonSynM [v1-vShift][h1-hShift] = vectorOfFrequencies[v1-vShift];
				}
			}
		}
	}
}

synSubsAVL = {};
nsSubsAVL  = {};

#include "TreeTools.ibf";

for (h1=0; h1 < T-1; h1=h1+1)
{
	abn = branchNames[h1];
	ExecuteCommands("GetInformation (aRateMx, givenTree."+abn+");");
	synSubs  = (horOnes*(aRateMx$synM))*vertOnes;
	nsynSubs = (horOnes*(aRateMx$nonSynM))*vertOnes;
	synSubs = synSubs[0]/3;
	nsynSubs = nsynSubs[0]/3;
	
	synSubsAVL[abn] = synSubs;
	nsSubsAVL [abn] = nsynSubs;
}

treeAVL = givenTree ^ 0;

synTreeString 		= PostOrderAVL2StringDistances (treeAVL, synSubsAVL); 
nonSynTreeString	= PostOrderAVL2StringDistances (treeAVL, nsSubsAVL);

fprintf (stdout, "\nE[Syn subs/nucleotide site] tree: \n\t",    synTreeString, 	   "\n");
fprintf (stdout, "\nE[Non-syn subs/nucleotide site] tree: \n\t", nonSynTreeString, "\n");

