/*---------------------------------------------------------------------------------------------------------------------------------------------*/

lfunction extractBranchLengthsFromTreeAsDict (tree_id) {
    _treeLengthDict = {};
    _bls = BranchLength (^tree_id, -1);
    _bns = BranchName  (^tree_id, -1);
    
    return _bls;
}


/*---------------------------------------------------------------------------------------------------------------------------------------------*/

lfunction getNucRevBranchLengthsAndParameters (datafilter_id, tree_id) {

   DataSetFilter nucs = CreateFilter          (^datafilter_id, 1);
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
   fprintf                   (stdout, LF, "\n", extractBranchLengthsFromTreeAsDict (&tree));
   
   return 0;
}