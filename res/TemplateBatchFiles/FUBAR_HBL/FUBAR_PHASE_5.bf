ExecuteAFile (codon_fit_file);

codonCharactersAgrument = {{"A","C","G","T"}
			  			   {"3",GeneticCodeExclusions,"",""}};

nn_sample_count = + (simPatterns[-1][1]);
byGridPoint     = {nn_sample_count, 3};

GetString         (lfInfo, codonLF, -1);
_treeCount      = Columns (lfInfo["Trees"]);
_filterNames    = lfInfo ["Datafilters"];

baseFrequencies = Eval ((lfInfo["Base frequencies"])[0]);
GetString (targetSpeciesOrdering, existingFit.ds_1,-1);

sequenceCount = existingFit.ds_1.species;
indexer       = 0;


accumulatedSequences = {};

for (k = 0; k < _treeCount; k+=1) {
    accumulatedSequences [k] = {};
    for (k2 = 0; k2 < sequenceCount; k2+=1) {
        (accumulatedSequences [k])[k2] = "";
        (accumulatedSequences [k])[k2] *128;
    }
}

for (k = 0; k < Rows(simPatterns); k += 1) {
    samples_for_this_value = simPatterns[k][1];
    alpha_beta  = grid[simPatterns[k][0]][-1];
    for (k2 = 0; k2 < samples_for_this_value; k2+=1) {
        byGridPoint[indexer][0] = alpha_beta[0];
        byGridPoint[indexer][1] = alpha_beta[1];
        byGridPoint[indexer][2] = 1+Random(0,_treeCount-0.000000001)$1;
        indexer += 1;
    }
}

byGridPoint  = byGridPoint % {{0,1,2}};

//fprintf (sim_grid_info, CLEAR_FILE, byGridPoint);

extraFunctions        = exportFunctionDefinition("simulateFromTreeGivenRates") + exportFunctionDefinition("mapSets");
LF_NEXUS_EXPORT_EXTRA = "\n;PRESERVE_SLAVE_NODE_STATE=1;MPI_NEXUS_FILE_RETURN=Rows(UserFunction);\nbaseFrequencies= " + baseFrequencies + ";\ncodonCharactersAgrument="+codonCharactersAgrument+"\n;targetSpeciesOrdering="+targetSpeciesOrdering+";\n"; 
                        

MPI_NODE_STATUS = {MPI_NODE_COUNT,1};

Export(lfExport,codonLF);

for (k = 1; k < MPI_NODE_COUNT; k+=1) {
    MPISend (k, lfExport);
}

for (k = 1; k < MPI_NODE_COUNT; k+=1) {
    MPIReceive (-1, fromNode, result);
 }

//------------------------------------------------------------------------------

currentIndex = 0;
t0 = Time(1);
for (k = 1; k < nn_sample_count; k+=1)
{
    if (Abs(byGridPoint [k][-1] - byGridPoint[k-1][-1]) > 1e-10) {
        treeID         = byGridPoint[k-1][2];
        dispatchABlock (byGridPoint[k-1][2],byGridPoint[k-1][0], byGridPoint[k-1][1],  k-currentIndex);
        currentIndex = k;
        SetParameter (STATUS_BAR_STATUS_STRING, "Simulating neg/neutral sites "+ (k) + "/" + nn_sample_count + " " + _formatTimeString(Time(1)-t0),0);
    }
}
   
dispatchABlock (byGridPoint[k-1][2],byGridPoint[k-1][0], byGridPoint[k-1][1],  k-currentIndex);

still_busy = +(MPI_NODE_STATUS["_MATRIX_ELEMENT_VALUE_>0"]);
for (k = 0; k < still_busy; k += 1) {
    processABlock (0);
}


for (k = 0; k < _treeCount; k+=1) {
    newData = ""; newData * 128;
    for (k2 = 0; k2 < sequenceCount; k2+=1) {
         (accumulatedSequences [k])[k2] *0;
         newData * (">" + targetSpeciesOrdering[k2] + "\n" + (accumulatedSequences [k])[k2] + "\n");
    }
    newData * 0;
    ExecuteCommands ("DataSet simulated_data_" + (k+1)  +  " = ReadFromString (newData)");
    ExecuteCommands ("GetDataInfo (info, " + _filterNames[k] +",\"PARAMETERS\");");
    ExecuteCommands ("DataSetFilter " + _filterNames[k] + " = CreateFilter (simulated_data_" + (k+1) + "," + info["ATOM_SIZE"] + ",,\"" + info["SEQUENCES_STRING"] + "\",\"" + info["EXCLUSIONS"] + "\");");
}

LIKELIHOOD_FUNCTION_OUTPUT = 7;
LF_NEXUS_EXPORT_EXTRA = "";
fprintf (sim_fit_file, CLEAR_FILE, codonLF);
DeleteObject (codonLF);
ExecuteAFile (sim_fit_file);
site_probs = computeLFOnGrid ("codonLF", grid,1);

fprintf (sim_grid_info, CLEAR_FILE, byGridPoint, "\n", site_probs);

tabulateGridResults (points, nn_sample_count, samples, _chainCount);
posteriorsUnderNN = reportSiteResults (nn_sample_count, 0, 0, 1);

return posteriorsUnderNN;

//------------------------------------------------------------------------------

function   dispatchABlock (treeID, alphaV, betaV, howMany) {
    if (MPI_NODE_COUNT <= 1) {
        simulatedChunk = simulateFromTreeGivenRates (treeID,alphaV,betaV,howMany, "baseFrequencies", "codonCharactersAgrument", targetSpeciesOrdering);
        return processABlock (treeID);
    }

    //fprintf(stdout, MPI_NODE_STATUS, "\n");
    
    for (node_id = 1; node_id < MPI_NODE_COUNT; node_id += 1) {
        if (MPI_NODE_STATUS[node_id] == 0) {
            break;
        }
    }
    //fprintf (stdout, "Node id:", node_id,"\n");
    
    if (node_id == MPI_NODE_COUNT) {
        node_id = processABlock(0);
    }
    
    mpi_command = extraFunctions + "\nreturn simulateFromTreeGivenRates ("+treeID+","+alphaV+","+betaV+","+howMany+",\"baseFrequencies\", \"codonCharactersAgrument\", targetSpeciesOrdering);";
    //fprintf (stdout, "[DEBUG:] To node ", node_id, "=>", mpi_command, "\n");
    MPISend (node_id, mpi_command);
    //fprintf (stdout, MPI_LAST_SENT_MSG);
    MPI_NODE_STATUS[node_id] = treeID;
    return 0;
}

//------------------------------------------------------------------------------

function   processABlock (treeID) {
   if (MPI_NODE_COUNT > 1) {
        MPIReceive (-1,fromNode,result);
        simulatedChunk = Eval (result);
        //fprintf (stdout, "Result:",result,"\n");
        treeID = MPI_NODE_STATUS[fromNode];
        MPI_NODE_STATUS[fromNode] = 0;
   }    
   
   for (k2 = 0; k2 < sequenceCount; k2+=1) {
        ((accumulatedSequences [treeID-1])[k2]) * simulatedChunk[k2];
   }           
   
   return fromNode;
}

//------------------------------------------------------------------------------

function simulateFromTreeGivenRates (treeID, alphaValue, betaValue, codons, frequencies, characters, speciesOrdering) {
    alpha = alphaValue;
    beta  = betaValue;
    ExecuteCommands ("DataSet theData = Simulate (codon_tree_"+treeID+",`frequencies`,`characters`,codons,0);");
    GetString (similatedOrdering, theData,-1);
    reordering = mapSets (speciesOrdering, similatedOrdering);
    DataSetFilter reporterFilter = CreateFilter (theData,1,"",Join(",",reordering));
    GetInformation (allSeqs, reporterFilter);
    GetString (similatedOrdering, reporterFilter,-1);
    return allSeqs;
}
