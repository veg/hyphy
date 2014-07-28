LoadFunctionLibrary ("TreeTools");

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

lfunction extractBranchLengthsFromTreeAsDict (tree_id) {
    _treeLengthDict = {};
    _bls = BranchLength (^tree_id, -1);
    _bns = BranchName  (^tree_id, -1);
    
    for (k = 0; k < Columns (_bns) - 1; k += 1) {
        _treeLengthDict[_bns[k]] = _bls[k];
    }
    
    return _treeLengthDict;
}


/*---------------------------------------------------------------------------------------------------------------------------------------------*/

lfunction getNucRevBranchLengthsAndParameters (dataset_id, tree_id) {

   DataSetFilter nucs = CreateFilter          (^dataset_id, 1);
   HarvestFrequencies   (Freqs, 	nucs, 1,1,1);

   global 	nAC = 1;
   global 	nAT = 1;
   global 	nCG = 1;
   global 	nCT = 1;
   global 	nGT = 1;
    
   revRateMatrix =  {{*,	nAC*t,t,	nAT*t}
                      {	nAC*t,*,	nCG*t,	nCT*t}
                      {t,	nCG*t,*,	nGT*t}
                      {	nAT*t,	nCT*t,	nGT*t,*}};
                    
   Model 	revQ = (revRateMatrix, 	Freqs);

   ExecuteCommands           ("Tree 	tree = " + Eval("Format (^tree_id,1,1)"));
   LikelihoodFunction 	LF = (nucs,	tree);
   Optimize                  (res,LF);
   res = {"AC": nAC, "AT": nAT, "CG": nCG, "CT": nCT, "GT": nGT};
   res ["lengths"] = extractBranchLengthsFromTreeAsDict (&tree);
   DeleteObject              (LF);
   
   return res;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function killZeroBranchesBasedOnNucFit (dataset_id, string) {

   //fprintf (stdout, "CHECK 1\n");

	

   DataSetFilter nucs = CreateFilter          (^dataset_id, 1);

   //fprintf (stdout, "CHECK 2\n");

   HarvestFrequencies   (Freqs, 	nucs, 1,1,1);

   //fprintf (stdout, "CHECK 3\n");

   global 	nAC = 1;
   global 	nAT = 1;
   global 	nCG = 1;
   global 	nCT = 1;
   global 	nGT = 1;
    
   revRateMatrix =  {{*,	nAC*t,t,	nAT*t}
                      {	nAC*t,*,	nCG*t,	nCT*t}
                      {t,	nCG*t,*,	nGT*t}
                      {	nAT*t,	nCT*t,	nGT*t,*}};
                    
   Model 	revQ = (revRateMatrix, 	Freqs);

   //fprintf (stdout, "CHECK 4\n");

   ExecuteCommands           ("Tree 	tree = " + string);
   LikelihoodFunction 	LF = (nucs,	tree);
   Optimize                  (res,LF);
   
   
   treeAVL = tree ^ 0;
   
   res = {"AC": nAC, "AT": nAT, "CG": nCG, "CT": nCT, "GT": nGT};
   res [ "_collapsed_tree" ]  =  KillInternalZeroBranchLengths (treeAVL);
   
   
   //DeleteObject              (LF);
   
   return res;
}

