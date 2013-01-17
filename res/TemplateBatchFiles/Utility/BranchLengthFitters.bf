function getNucRevBranchLengthsAndParameters (datafilter_id, tree_id) {

   ExecuteCommands ("DataSetFilter _getNucRevBranchLengthsAndParameters.nucs = CreateFilter (`datafilter_id`, 1)");
   HarvestFrequencies(_getNucRevBranchLengthsAndParameters.Freqs, _getNucRevBranchLengthsAndParameters.nucs, 1,1,1);

   global _getNucRevBranchLengthsAndParameters.AC = 1;
   global _getNucRevBranchLengthsAndParameters.AT = 1;
   global _getNucRevBranchLengthsAndParameters.CG = 1;
   global _getNucRevBranchLengthsAndParameters.CT = 1;
   global _getNucRevBranchLengthsAndParameters.GT = 1;
    
   _getNucRevBranchLengthsAndParameters.revRateMatrix =  {{*,_getNucRevBranchLengthsAndParameters.AC*t,t,_getNucRevBranchLengthsAndParameters.AT*t}
                      {_getNucRevBranchLengthsAndParameters.AC*t,*,_getNucRevBranchLengthsAndParameters.CG*t,_getNucRevBranchLengthsAndParameters.CT*t}
                      {t,_getNucRevBranchLengthsAndParameters.CG*t,*,_getNucRevBranchLengthsAndParameters.GT*t}
                      {_getNucRevBranchLengthsAndParameters.AT*t,_getNucRevBranchLengthsAndParameters.CT*t,_getNucRevBranchLengthsAndParameters.GT*t,*}};
                    
   Model _getNucRevBranchLengthsAndParameters.revQ = (_getNucRevBranchLengthsAndParameters.revRateMatrix, _getNucRevBranchLengthsAndParameters.Freqs);

   ExecuteCommands ("Tree _getNucRevBranchLengthsAndParameters.Tree = " + Eval("Format (`tree_id`,1,1)"));
   LikelihoodFunction _getNucRevBranchLengthsAndParameters.LF = (_getNucRevBranchLengthsAndParameters.nucs,_getNucRevBranchLengthsAndParameters.Tree);
   Optimize (_getNucRevBranchLengthsAndParameters.res, _getNucRevBranchLengthsAndParameters.LF);
   
  
   
   return 0;
}