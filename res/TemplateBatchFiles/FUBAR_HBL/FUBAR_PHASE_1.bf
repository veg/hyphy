fprintf (stdout, "[FUBAR PHASE 1] Optimizing relative branch lengths under the nucleotide REV model\n");

vectorOfFrequencies = overallFrequencies;
// 'overallFrequencies' is supplied by _MFReader_.ibf

SKIP_HARVEST_FREQ = 1;
LoadFunctionLibrary ("GRM", {"00":"Global"});

populateTrees ("nuc_tree", fileCount);
ExecuteCommands(constructLF ("nucLF", "nucData", "nuc_tree", fileCount));
Optimize (nuc_res, nucLF);

fprintf (stdout, "[FUBAR PHASE 1 FINISHED] log(L) = ", nuc_res[1][0], "\n");
for (k = 1; k <= fileCount; k += 1)
{
    fprintf (stdout, "\tLength of tree ", k, " (substitutions/site) = ", +Eval("BranchLength(nuc_tree_" + k + ",-1)"),"\n");
}

LF_NEXUS_EXPORT_EXTRA = "positionalFrequencies = " + positionFrequencies + ";";
SetDialogPrompt ("Write nucleotide model fit to");
LIKELIHOOD_FUNCTION_OUTPUT = 7;
fprintf (PROMPT_FOR_FILE, CLEAR_FILE, nucLF);
DeleteObject (nucLF);
