ChoiceList  (response,"Output options",1,NO_SKIP,
			 "Display","Display trees on the console (and graphically if GUI is present).",
			 "Save to file","Save tree strings to a file (syn followed by non-syn)");
			 
if (response<0)
{
	return 0;
}

branchNames = BranchName (givenTree,-1);
T 			= Columns    (branchNames);




synSubsAVL = {};
dSAVL	   = {};
nsSubsAVL  = {};
dNAVL	   = {};

ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Utility"+DIRECTORY_SEPARATOR+"CodonTools.bf");
_snsAVL	   = _computeSNSSites ("filteredData", _Genetic_Code, vectorOfFrequencies, 0);

ExecuteAFile ("TreeTools.ibf");

nonStopSites = sSites+nsSites;				 

fprintf (stdout, "\nTotal nucleotide sites :", filteredData.sites*3,
				 "\nTotal sense codon sites:", nonStopSites,
				 "\nSynonymous  sites      :", sSites, 
				 "\nNonsynonymous  sites   :", nsSites, "\n");

sSites  = nonStopSites/sSites;
nsSites = nonStopSites/nsSites;
				 

for (h1=0; h1 < T-1; h1=h1+1)
{
	abn = branchNames[h1];
	ExecuteCommands("GetInformation (aRateMx, givenTree."+abn+");");
	synSubs  = ((horOnes*(aRateMx$synM))*vertOnes)[0]/3;
	nsynSubs = ((horOnes*(aRateMx$nonSynM))*vertOnes)[0]/3;
	
	synSubsAVL[abn] = synSubs;
	nsSubsAVL [abn] = nsynSubs;
	dSAVL[abn]	    = synSubs *sSites;
	dNAVL[abn]	    = nsynSubs*nsSites;
}

treeAVL = givenTree ^ 0;

synTreeString 		= PostOrderAVL2StringDistances (treeAVL, synSubsAVL); 
nonSynTreeString	= PostOrderAVL2StringDistances (treeAVL, nsSubsAVL);
dSTreeString 		= PostOrderAVL2StringDistances (treeAVL, dSAVL); 
dNTreeString	    = PostOrderAVL2StringDistances (treeAVL, dNAVL);


if (response == 0)
{
	fprintf (stdout, "\nE[Syn subs/nucleotide site] tree: \n\t",     synTreeString, 	   "\n");
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


	OpenWindow (TREEWINDOW, mxTreeSpec,"(SCREEN_WIDTH-50)/2;(SCREEN_HEIGHT-50)/2;10;40");

	mxTreeSpec [0] = "synSubsTree";
	OpenWindow (TREEWINDOW, mxTreeSpec,"(SCREEN_WIDTH-50)/2;(SCREEN_HEIGHT-50)/2;30+(SCREEN_WIDTH-30)/2;40");

	mxTreeSpec [0] = "dSTree";
	OpenWindow (TREEWINDOW, mxTreeSpec,"(SCREEN_WIDTH-50)/2;(SCREEN_HEIGHT-50)/2;10;45+(SCREEN_HEIGHT-50)/2");

	mxTreeSpec [0] = "dNTree";
	OpenWindow (TREEWINDOW, mxTreeSpec,"(SCREEN_WIDTH-50)/2;(SCREEN_HEIGHT-50)/2;30+(SCREEN_WIDTH-30)/2;45+(SCREEN_HEIGHT-50)/2");
}
else
{
	SetDialogPrompt ("Save trees to file:");
	fprintf (PROMPT_FOR_FILE, CLEAR_FILE, synTreeString, ";\n", nonSynTreeString, ";\n",dSTreeString, ";\n", dNTreeString, ";\n");
}
