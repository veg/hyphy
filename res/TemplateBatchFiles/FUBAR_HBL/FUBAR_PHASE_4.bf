fscanf              (stdin, "String", nuc_fit_file);
fscanf              (stdin, "String", grid_file);
fscanf              (stdin, "String", sample_base_file);
fscanf              (stdin, "Number", _chainCount);
fscanf              (stdin, "String", results_file);

// for PHASE 5
fscanf (stdin, "String", sim_fit_file);
fscanf (stdin, "String", sim_grid_info);
fscanf (stdin, "String", codon_fit_file);

ExecuteAFile        (PATH_TO_CURRENT_BF + "FUBAR_tools.ibf");
LoadFunctionLibrary ("GrabBag");
LoadFunctionLibrary ("WriteDelimitedFiles");

ExecuteAFile        (nuc_fit_file);
sequences = nucData_1.species;
GetInformation (treeCount,"^nuc_tree_[0-9]+$");
fileCount       = Columns (treeCount);
treeLengths     = {fileCount,1};

for (fileID = 1; fileID <= fileCount; fileID += 1)
{
	treeLengths [fileID-1] = + Eval("BranchLength(nuc_tree_"+fileID+",-1)");
}

fscanf (grid_file, REWIND, "NMatrix,Raw", grid, site_probs);
site_probs = Eval (site_probs);
sites   = Columns (site_probs["conditionals"]);


readMCMCSamples (sample_base_file,_chainCount);


notPositiveSelection = {points,1} ["grid[_MATRIX_ELEMENT_ROW_][0]>=grid[_MATRIX_ELEMENT_ROW_][1]"];
nonPositiveCount     = +notPositiveSelection;

priorMean            = {1, points};
sampleFromThisDistro = {nonPositiveCount,2};

tabulateGridResults (points, sites, samples, _chainCount);

from = 0;
for (_point = 0; _point < points; _point += 1) {
    priorMean [_point] = (+jointSamples[-1][_point])/samples;
    if (notPositiveSelection [_point]) {
        sampleFromThisDistro [from][0] = _point;
        sampleFromThisDistro [from][1] = priorMean [_point];
        from += 1;
    }
}

priorNN = +(sampleFromThisDistro [-1][1]);
fubar_results = reportSiteResults   (sites, 0, priorNN, _fubar_do_simulations);
fubarRowCount     = Rows (fubar_results);

if (_fubar_do_simulations) {
    simPatterns = Random(sampleFromThisDistro, {"PDF":"Multinomial","ARG0":1000});
    fprintf (sample_base_file, CLEAR_FILE, _chainCount, "\n", jointLogL, "\n", jointSamples);
    
    fprintf (stdout, "\n[FUBAR PHASE 5] Performing 1,000 simulations under the data-derived composite null model to derive False Discovery Rate (FDR) estimates\n");
    ExecuteAFile (PATH_TO_CURRENT_BF + "FUBAR_PHASE_5.bf");
    
    // sort posteriorsUnderNN on the posterior prob of pos.sel. column
    
    posteriorsNNPP    = (posteriorsUnderNN[-1][3]) % 0;
    sortedFubarP      = ({Rows (fubar_results), 2} ["_MATRIX_ELEMENT_ROW_*(_MATRIX_ELEMENT_COLUMN_==0)+fubar_results[_MATRIX_ELEMENT_ROW_][3]*(_MATRIX_ELEMENT_COLUMN_==1)"])%1;
    
    currentNNRows     = Rows (posteriorsUnderNN);
    
    currentNNIndex    = currentNNRows - 1;
    
    
    for (currentFubarIndex =  fubarRowCount - 1; currentFubarIndex >= 0; currentFubarIndex += -1) {
        currentFubarPosteriorP = sortedFubarP[currentFubarIndex][1];
        
        while (currentNNIndex > 0 && posteriorsNNPP[currentNNIndex] > currentFubarPosteriorP) {
            currentNNIndex = currentNNIndex - 1;
        }
        
        FDR = Min (1,(currentNNRows-currentNNIndex)/currentNNRows * priorNN / ((fubarRowCount-currentFubarIndex) / fubarRowCount));
        
        if (currentNNIndex == 0) {
            break;
        }
        fubar_results[sortedFubarP[currentFubarIndex][0]][7] = FDR;
        
    }
    
    for (; currentFubarIndex >= 0; currentFubarIndex += -1) {
        FDR = Min (1, priorNN / ((fubarRowCount-currentFubarIndex) / fubarRowCount));
        fubar_results[sortedFubarP[currentFubarIndex][0]][7] = FDR;
    }
}

site_counter = {};
for (currentFubarIndex = 0; currentFubarIndex < fubarRowCount; currentFubarIndex += 1) {
    site_counter + (currentFubarIndex+1);
}

if (_fubar_do_simulations) {
    WriteSeparatedTable (results_file, {{"Codon","alpha","beta","beta-alpha","Prob[alpha<beta]", "Prob[alpha>beta]", "BayesFactor","PSRF", "Neff", "FDR"}}, fubar_results, site_counter, ",");
} else {
    WriteSeparatedTable (results_file, {{"Codon","alpha","beta","beta-alpha","Prob[alpha<beta]", "Prob[alpha>beta]", "BayesFactor","PSRF", "Neff"}}, fubar_results, site_counter, ",");
}
 