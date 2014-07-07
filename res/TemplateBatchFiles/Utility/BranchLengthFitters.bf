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

   global 	AC = 1;
   global 	AT = 1;
   global 	CG = 1;
   global 	CT = 1;
   global 	GT = 1;
    
   revRateMatrix =  {{*,	AC*t,t,	AT*t}
                      {	AC*t,*,	CG*t,	CT*t}
                      {t,	CG*t,*,	GT*t}
                      {	AT*t,	CT*t,	GT*t,*}};
                    
   Model 	revQ = (revRateMatrix, 	Freqs);

   ExecuteCommands           ("Tree 	tree = " + Eval("Format (^tree_id,1,1)"));
   LikelihoodFunction 	LF = (nucs,	tree);
   Optimize                  (res,LF);
   DeleteObject              (LF);
   
   return {"AC": AC, "AT": AT, "CG": CG, "CT": CT, "GT": GT, "lengths": extractBranchLengthsFromTreeAsDict (&tree)};
}
