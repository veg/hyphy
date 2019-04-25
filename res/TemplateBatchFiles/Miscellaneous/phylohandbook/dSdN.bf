/* 

THIS HyPhy BATCH LANGUAGE FILE IS PROVIDED AS AN EXAMPLE/TOOL FOR 
THE 2nd EDITION OF 'The Phylogenetic Handbook'

Written by 
Sergei L Kosakovsky Pond (spond@ucsd.edu)
Art FY Poon				 (apoon@ucsd.edu)
Simon DW Frost			 (sdfrost@ucsd.edu)

*/

RequireVersion ("0.9920060830");
VERBOSITY_LEVEL = -1;

/*---------------------------------------------------------------------------------------------------------------------------------------------*/


fprintf (stdout, "\n\nThis example analysis fits the local MG94xREV model \n",
				 "(separate alpha and beta on each branch) to a codon dataset,\n",
				 "and computes several distance measures between sequences:\n",
				 " - expected number of substitutions per codon\n",
				 " - expected number of syn. substitutions per codon\n",
				 " - expected number of nonsyn. substitutions per codon\n",
				 " - dS\n",
				 " - dN\n",
				 "and provides the example contrasting the estimation of dN and dS\n",
				 "between a pair of sequences based on the phylogenetic tree and\n",
				 "direct pairwise estimation.\n"); 

ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Utility"+DIRECTORY_SEPARATOR+"DescriptiveStatistics.bf");
ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"chooseGeneticCode.def");
ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"dSdNTreeTools.ibf");

ChoiceList (freqs, "Alignment", 1, SKIP_NONE, "Flu H5N1 HA", 		"Use the example Influenza H5N1 hemagglutinin alignment (5 sequences)",
											  "Drospophila adh", 	"Use the example Drosophila ADH alignment (6 sequences).",
											  "Custom", 			"Load your own alignment and tree.");
													 
if (freqs < 0)
{
	return 0;
}
													 
if (freqs == 2)
{
	SetDialogPrompt     ("Choose a nucleotide alignment");
	DataSet ds        = ReadDataFile (PROMPT_FOR_FILE);
}
else
{
	if (freqs == 1)
	{
		DataSet ds = ReadDataFile (PATH_TO_CURRENT_BF + '/datasets/Drosophilia_adh.nex');
	}
	else
	{
		DataSet ds = ReadDataFile (PATH_TO_CURRENT_BF + '/datasets/H5N1_HA_5.nex');
	}	
}

DataSetFilter	  	filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);

SKIP_MODEL_PARAMETER_LIST = 1;
modelDesc 				  = "012345";
modelType 				  = 0;
ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"MG94custom.mdl");
SKIP_MODEL_PARAMETER_LIST = 0;

ExecuteAFile 		(HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"queryTree.bf");

LikelihoodFunction  theLnLik = (filteredData, givenTree);

fprintf (stdout, "Fitting local MG94xREV to observed data... \n");
Optimize (res, theLnLik);

nonStopCount = ModelMatrixDimension;
ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Distances"+DIRECTORY_SEPARATOR+"CodonTools.def");

sSites  = 0;
nsSites = 0;

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
		sSites = sSites   + filteredData.sites * _S_NS_POSITIONS_[0][h1] * vectorOfFrequencies[h1-hShift];
		nsSites = nsSites + filteredData.sites * _S_NS_POSITIONS_[1][h1] * vectorOfFrequencies[h1-hShift];
		
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
dSAVL	   = {};
nsSubsAVL  = {};
dNAVL	   = {};

ExecuteAFile 		(HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TreeTools.ibf");

fprintf (stdout, "\nTotal nucleotide sites :", filteredData.sites*3,
				 "\nSynonymous  sites      :", sSites, 
				 "\nNonsynonymous  sites   :", nsSites, "\n");
				 
sSites  = 3*filteredData.sites/sSites;
nsSites = 3*filteredData.sites/nsSites;
branchNames = BranchName (givenTree,-1);
total = Columns(branchNames)-1;

branchReport = {total,6};
labelString  = "";

totalLength  = {4,1};

for (h1=0; h1 < total; h1=h1+1)
{
	abn = branchNames[h1];
	ExecuteCommands("GetInformation (aRateMx, givenTree."+abn+");");
	synSubs  = (horOnes*(aRateMx$synM))*vertOnes;
	nsynSubs = (horOnes*(aRateMx$nonSynM))*vertOnes;
	synSubs = synSubs[0]/3;
	nsynSubs = nsynSubs[0]/3;
	
	synSubsAVL[abn] = synSubs;
	nsSubsAVL [abn] = nsynSubs;
	dSAVL[abn]	    = synSubs *sSites;
	dNAVL[abn]	    = nsynSubs*nsSites;
	
	ExecuteCommands("localOmega =  givenTree."+abn+".nonSynRate/givenTree."+abn+".synRate;");
	
	branchReport[h1][0] = synSubs;
	branchReport[h1][1] = nsynSubs;
	branchReport[h1][2] = dSAVL[abn];
	branchReport[h1][3] = dNAVL[abn];
	branchReport[h1][4] = dNAVL[abn]/dSAVL[abn];
	branchReport[h1][5] = localOmega;
	
	for (h2 = 0; h2<4; h2=h2+1)
	{
		totalLength[h2] = totalLength[h2] + branchReport[h1][h2];
	}
	labelString = labelString+";"+abn;
}

fprintf (stdout, "\nTotal tree lengths:",
				 "\n\tE[syn]        : ", totalLength[0], 
				 "\n\tE[non-syn]    : ", totalLength[1], 
				 "\n\tdS            : ", totalLength[2], 
				 "\n\tdN            : ", totalLength[3], "\n"); 
				 

treeAVL = givenTree ^ 0;

synTreeString 		= PostOrderAVL2StringDistances (treeAVL, synSubsAVL); 
nonSynTreeString	= PostOrderAVL2StringDistances (treeAVL, nsSubsAVL);
dSTreeString 		= PostOrderAVL2StringDistances (treeAVL, dSAVL); 
dNTreeString	    = PostOrderAVL2StringDistances (treeAVL, dNAVL);


fprintf (stdout, "\n\nE[Syn subs/nucleotide site] tree: \n\t",     synTreeString, 	   "\n");
fprintf (stdout, "\nE[Non-syn subs/nucleotide site] tree: \n\t", nonSynTreeString, "\n");
fprintf (stdout, "\ndS tree: \n\t", dSTreeString, "\n");
fprintf (stdout, "\ndN tree: \n\t", dNTreeString, "\n");

UseModel (USE_NO_MODEL);

Tree 	synSubsTree 	= synTreeString;
Tree	nonsynSubsTree 	= nonSynTreeString;
Tree 	dSTree 	= dSTreeString;
Tree	dNTree 	= dNTreeString;

mxTreeSpec  = {5,1};

mxTreeSpec [0] = "nonsynSubsTree";
mxTreeSpec [3] = "";
mxTreeSpec [4] = "Inferred_Tree."+nodeName;
mxTreeSpec [1] = "8211";
mxTreeSpec [2] = "";



mxTreeSpec [0] = "synSubsTree";

mxTreeSpec [0] = "dSTree";

mxTreeSpec [0] = "dNTree";

/* start: constrain global parameters for pairwise estimates */

GetString 			  (likelihoodInfo, theLnLik, -1);   
globalVariableList 	= likelihoodInfo["Global Independent"];
globalVariableCount	= Columns (globalVariableList);

if (globalVariableCount)
{
	stashed_GV = {globalVariableCount,1};
	for (vc = 0; vc < globalVariableCount; vc = vc + 1)
	{
		ExecuteCommands (globalVariableList[vc]+":="+globalVariableList[vc]+"__;");
	}
}

/* end: constrain global parameters for pairwise estimates */

taxa = {{0,0}};
theStencil = ComputeScalingStencils (0);
tc = filteredData.species*(filteredData.species-1)/2;

fprintf (stdout, "\nEstimating pairwise dN/dS: total of ", tc, " analyses\n");

branchReport2 = {tc,4};
labelString  = "";
tc = 0;

for (tax1 = 0; tax1 < filteredData.species; tax1=tax1+1)
{
	for (tax2 = tax1+1; tax2 < filteredData.species; tax2=tax2+1)
	{
		GetString (tx1, filteredData, tax1);
		GetString (tx2, filteredData, tax2);
		BRANCH_LENGTH_STENCIL = theStencil["Syn"];
		dS = BranchLength(givenTree, tx1+";"+tx2);
		BRANCH_LENGTH_STENCIL = theStencil["NonSyn"];
		dN = BranchLength(givenTree, tx1+";"+tx2);
		BRANCH_LENGTH_STENCIL = 0;
		labelString = labelString + ";" + tx1 + "," + tx2;
		branchReport2 [tc][0] = dS*sSites;
		branchReport2 [tc][1] = dN*nsSites; 
						 
		DataSetFilter twoSeqFilter = CreateFilter (ds,3,"",speciesIndex == tax1 || speciesIndex == tax2,GeneticCodeExclusions);
		twoTreeString = "("+tx1+","+tx2+")";
		UseModel (MG94customModel);
		Tree t2 = twoTreeString;
		LikelihoodFunction  twoLnLik = (twoSeqFilter, t2);
		Optimize (res2, twoLnLik);
		v1 = ReturnVectorsOfCodonLengths (theStencil,"t2");
		branchReport2 [tc][2] = (scaledVectors["Syn"])[0]*sSites;
		branchReport2 [tc][3] = (scaledVectors["NonSyn"])[0]*nsSites; 
		tc = tc + 1;
	}
}

VERBOSITY_LEVEL = 0;