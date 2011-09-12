treeCount = Rows ("Tree");
lf_Count  = Rows ("LikelihoodFunction");

if (treeCount > 0 && lf_Count > 0)
{
	treeChoices = {treeCount,2};
	for (k=0; k<treeCount; k=k+1)
	{
		GetString (treeID, Tree, k);
		treeChoices [k][0] = treeID;
		treeChoices [k][1] = "Tree " + treeID;
	}
	
	chosenTree = 0;
	if (treeCount > 1)
	{
		ChoiceList (chosenTree, "Which tree?", 1, SKIP_NONE, treeChoices);
		if (chosenTree < 0)
		{
			return;
		}
	}
	
	_treeID = treeChoices[chosenTree][0];

	for (lf_ID = 0; lf_ID < lf_Count; lf_ID = lf_ID + 1)
	{
		GetString (treeID, LikelihoodFunction,lf_ID);
		ExecuteCommands ("GetString(lfInfo,"+treeID+",-1);");
		lfTrees = lfInfo["Trees"];
		for (k = 0; k<Columns(lfTrees); k=k+1)
		{
			if (lfTrees[k] == treeChoices[chosenTree][0])
			{
				break;
			}
		}
		if (k < Columns(lfTrees))
		{
			fprintf (stdout, "\nTree ",Columns(lfTrees)," is a part of likelihood function ", treeID, "\n");
			dfName = (lfInfo["Datafilters"])[k];
			break;
		}
		else
		{
			return 0;
		}
		
	}
		
	ExecuteAFile(HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"chooseGeneticCode.def");

	ChoiceList  (response,"Output options",1,NO_SKIP,
				 "Display","Display trees on the console (and graphically if GUI is present).",
				 "Save to file","Save tree strings to a file (syn followed by non-syn)");
				 
	if (response<0)
	{
		return 0;
	}

	ExecuteCommands ("branchNames = BranchName ("+_treeID+",-1);");
	T = Columns (branchNames);

	ExecuteCommands(
	"GetInformation (aRateMx, "+_treeID+"."+branchNames[0]+");");

	/* make syn and non-syn template matrices */

	nonStopCount = Columns (aRateMx);
	ExecuteAFile(HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Distances"+DIRECTORY_SEPARATOR+"CodonTools.def"");

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

	#include "TreeTools.ibf";

	fprintf (stdout, "\nTotal nucleotide sites :", filteredData.sites*3,
					 "\nSynonymous  sites      :", sSites, 
					 "\nNonsynonymous  sites   :", nsSites, "\n");
					 
	sSites  = filteredData.sites/sSites;
	nsSites = filteredData.sites/nsSites;

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
		fprintf (PROMPT_FOR_FILE, CLEAR_FILE, synTreeString, "\;\n", nonSynTreeString, "\;\n",dSTreeString, "\;\n", dNTreeString, "\;\n");
	}
}
