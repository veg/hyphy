RequireVersion ("2.1320130313");

VERBOSITY_LEVEL				= 0;
maximum_number_of_omegas   = 10;
skipCodeSelectionStep 		= 0;

LoadFunctionLibrary("chooseGeneticCode");
LoadFunctionLibrary("GrabBag");
LoadFunctionLibrary("dSdNTreeTools");
LoadFunctionLibrary("CF3x4");
LoadFunctionLibrary("BranchSiteTemplate");
LoadFunctionLibrary("TreeTools");



DataSet 			ds 				= ReadDataFile(PROMPT_FOR_FILE);
DataSetFilter 		dsf 			= CreateFilter(ds,3,"","",GeneticCodeExclusions);
HarvestFrequencies	(nuc3, dsf, 3, 1, 1);


nucCF						= CF3x4	(nuc3, GeneticCodeExclusions);

for (matrix_id = 1; matrix_id <= maximum_number_of_omegas; matrix_id += 1) {
    PopulateModelMatrix			  ("MGMatrix" + matrix_id,  nucCF, "t", "omega" + matrix_id, "");
}

PopulateModelMatrix			  ("MGMatrixLocal",  nucCF, "syn", "", "nonsyn");
codon3x4					= BuildCodonFrequencies (nucCF);

Model		MGL				= (MGMatrixLocal, codon3x4, 0);

LoadFunctionLibrary("queryTree");


totalBranchCount			 = BranchCount(givenTree) + TipCount (givenTree);
bNames						 = BranchName (givenTree, -1);

//******* SELECT WHICH BRANCHES TO TEST ********//

selectTheseForTesting = {totalBranchCount + 2, 2};
selectTheseForTesting [0][0] = "None"; selectTheseForTesting [0][1] = "Just fit the branch-site REL model";
selectTheseForTesting [1][0] = "All";  selectTheseForTesting [1][1] = "Test all branches";

for (k = 0; k < totalBranchCount; k = k+1) {
    selectTheseForTesting [k+2][0] = bNames[k];
    selectTheseForTesting [k+2][1] = "Test branch '" + bNames[k] + "'";
    
}

ChoiceList  (whichBranchesToTest,"Which branches to test?",0,NO_SKIP,selectTheseForTesting);

if (whichBranchesToTest[0] < 0) {
    return;
}

selectedBranches = {};
for (k = 0; k < Columns (whichBranchesToTest); k += 1) {
    if (whichBranchesToTest [k] == 0) {
        selectedBranches = {}; break;
    }
    if (whichBranchesToTest [k] == 1) {
        for (k = 0; k <  totalBranchCount; k += 1) {
            selectedBranches[k] = 1;
        }
        break;
    }
    selectedBranches [whichBranchesToTest [k] - 2] = 1;
}

fprintf (stdout, "Selected the following branches for testing");
selectedBranches["printSelectedBranches"][""];

function printSelectedBranches (key, value) {
    fprintf (stdout, "\n\t", bNames[0+key]);
    return 0;
}

fprintf (stdout, "\n");
 
//******* BASE MODEL FITTING ********//
 

SetDialogPrompt ("Save analysis results to");

lrt_column = 4;
p_uncorrected_column = 5;
p_corrected_column = 6;
branch_length_column = 7;

fprintf (PROMPT_FOR_FILE, CLEAR_FILE, KEEP_OPEN,"Branch,Mean_dNdS,RateClasses,OmegaOver1,WtOmegaOver1,LRT,p,p_Holm,BranchLength");
csvFilePath = LAST_FILE_PATH;

fprintf 					  (stdout, "[PHASE 0] Fitting the local MG94 (no site-to-site variation) to obtain initial parameter estimates\n");

LikelihoodFunction	base_LF	 = (dsf, givenTree);

Optimize					  (res_base,base_LF);

lfOut	= csvFilePath + ".mglocal.fit";
LIKELIHOOD_FUNCTION_OUTPUT = 7;
fprintf (lfOut, CLEAR_FILE, base_LF);
LIKELIHOOD_FUNCTION_OUTPUT = 2;

localLL						 = res_base[1][0];
localParams					 = res_base[1][1] + 9;

LoadFunctionLibrary			 ("DescriptiveStatistics");

//GetInformation	   		   (varNames, "givenTree\\..*\\.omega1");
//localOmegaValues			 = {totalBranchCount,1}["Eval(varNames[_MATRIX_ELEMENT_ROW_])"];

pValueByBranch				  = {totalBranchCount,8};

for (k = 0; k < totalBranchCount; k = k+1) {
	srate  = Eval ("givenTree." + bNames[k] + ".syn");
	nsrate = Eval ("givenTree." + bNames[k] + ".nonsyn");
	if (srate > 0) {
		pValueByBranch [k][0] = Min (10, nsrate/srate);
	}
	else {
		pValueByBranch [k][0] = 10;
	}	
}

omegaStats					 = GatherDescriptiveStats (pValueByBranch[-1][0]);
current_parameter_count = localParams;
current_log_L           = localLL;
sample_size             = dsf.sites * dsf.species;
current_IC             = getIC (current_log_L, current_parameter_count, sample_size);


fprintf						 (stdout, "\nLog L = ", localLL, " with ", localParams, " degrees of freedom. IC = ", current_IC, "\n");

PrintDescriptiveStats		 ("Branch omega values", omegaStats);


for (matrix_id = 2; matrix_id <= maximum_number_of_omegas; matrix_id += 1) {
    freq_weights = generateFrequencyParameterization ("Paux", matrix_id);
    for (k = 0; k < matrix_id; k += 1) {
        freq_weights[k] =  "Exp(MGMatrix" + (k+1) + ")*" + freq_weights[k];
    }
    ExecuteCommands ("Model 		MG" + matrix_id + "=(\"" + Join("+",freq_weights) + "\",codon3x4,EXPLICIT_FORM_MATRIX_EXPONENTIAL);");
}


ASSUME_REVERSIBLE_MODELS	  = 1;
USE_LAST_RESULTS              = 1;
OPTIMIZATION_METHOD           = 0;
Tree stepupTree               = givenTree;



fprintf 					  (stdout, "[PHASE 1] Fitting Branch Site REL models to one branch at a time\n");

mg94bls   = BranchLength (givenTree,-1);
sortedBLs = {totalBranchCount, 2}["mg94bls[_MATRIX_ELEMENT_ROW_]*(_MATRIX_ELEMENT_COLUMN_==0)+_MATRIX_ELEMENT_ROW_*(_MATRIX_ELEMENT_COLUMN_==1)"];
sortedBLs = sortedBLs%0;

rate_classes_by_branch = {};

for (k = 0; k < totalBranchCount; k+=1) {
    reordered_index = sortedBLs[totalBranchCount-k-1][1];
    //reordered_index = sortedBLs[k][1];
    local_branch_name = bNames[reordered_index];
    accepted_rates_count = 1;
    
    branch_name_to_test_base = "stepupTree.`local_branch_name`";
    ExecuteCommands ("SetParameter (`branch_name_to_test_base`, MODEL, MG2);");
    
    Tree         mixtureTree = stepupTree;

    for (l = 0; l < totalBranchCount; l+=1) {
        if (l != reordered_index) {
            _constrainVariablesAndDescendants ("mixtureTree." + bNames[l]);
        }
    }
    
    branch_name_to_test = "mixtureTree." + bNames[reordered_index];
    setBSRELBounds (branch_name_to_test, 2);
    copyomegas (branch_name_to_test, Eval ("givenTree.`local_branch_name`.syn"), Eval ("givenTree.`local_branch_name`.nonsyn"));
    
    
    LikelihoodFunction stepupLF = (dsf, mixtureTree);
    fixGlobalParameters           ("stepupLF");
    Optimize                      (res, stepupLF);
    test_IC = getIC (res[1][0], current_parameter_count + 2, sample_size);
    fprintf 					  (stdout, "\n[PHASE 1] Branch ", local_branch_name, " log(L) = ", res[1][0], ", IC = ", test_IC, "\n\t2 rate clases\n\t");
    printNodeDesc ("mixtureTree.`local_branch_name`", 2);
    
    while (test_IC < current_IC) {
        accepted_rates_count += 1;
        current_parameter_count += 2;
        current_IC = test_IC;
        saved_MLEs = saveNodeMLES (branch_name_to_test, accepted_rates_count);
        //fprintf (stdout, saved_MLEs, "\n");
        ExecuteCommands ("SetParameter (`branch_name_to_test_base`, MODEL, MG" + (accepted_rates_count+1) +");");
        
        Tree         mixtureTree = stepupTree;
        for (l = 0; l < totalBranchCount; l+=1) {
            if (l != reordered_index) {
                _constrainVariablesAndDescendants ("mixtureTree." + bNames[l]);
            }
        }
        setBSRELBounds (branch_name_to_test, accepted_rates_count+1);
        initNodeMLESPlus1 (branch_name_to_test, saved_MLEs);
        LikelihoodFunction stepupLF = (dsf, mixtureTree);
        fixGlobalParameters           ("stepupLF");
        Optimize                      (res, stepupLF);
        test_IC = getIC (res[1][0], current_parameter_count + 2, sample_size);
        fprintf 					  (stdout, "\n[PHASE 1] Branch ", local_branch_name, " log(L) = ", res[1][0], ", IC = ", test_IC, "\n\t", accepted_rates_count+1, " rate clases\n\t");
        printNodeDesc ("mixtureTree.`local_branch_name`", accepted_rates_count + 1);
    }
    
    if (accepted_rates_count > 1) {
        ExecuteCommands ("SetParameter (`branch_name_to_test_base`, MODEL, MG" + (accepted_rates_count) +");");
        Tree         mixtureTree = stepupTree;
        restoreNodeMLE (branch_name_to_test, saved_MLEs);
    }

    Tree stepupTree = mixtureTree;
    
    if (accepted_rates_count == 1) {
        ExecuteCommands ("SetParameter (`branch_name_to_test_base`, MODEL, MGL);");
        ExecuteCommands ("`branch_name_to_test_base`.syn = givenTree.`local_branch_name`.syn");
        ExecuteCommands ("`branch_name_to_test_base`.nonsyn = givenTree.`local_branch_name`.nonsyn");
    }

    rate_classes_by_branch [reordered_index] = accepted_rates_count;
}


Tree mixtureTree = stepupTree;
ClearConstraints (mixtureTree);

fprintf (stdout, "\n\n[INFERRED MODEL COMPLEXITY]");

for (k = 0; k < totalBranchCount; k+=1) {
    branch_name_to_test = "mixtureTree." + bNames[k];
    pValueByBranch [k][1] = rate_classes_by_branch[k];
    if (rate_classes_by_branch [k] != 1) {
        setBSRELBounds (branch_name_to_test, rate_classes_by_branch[k]);
    }
    fprintf (stdout, "\n\t", branch_name_to_test, " has ", rate_classes_by_branch [k], " site rate classes");
}
fprintf (stdout, "\n");

LikelihoodFunction three_LF   = (dsf,mixtureTree);
unconstrainGlobalParameters ("three_LF");

VERBOSITY_LEVEL     = 0;
USE_LAST_RESULTS    = 1;
OPTIMIZATION_METHOD = 0;



fprintf 					  (stdout, "\n\n[PHASE 2] Fitting the full LOCAL alternative model (no constraints)\n");
Optimize					  (res_three_LF,three_LF);
fprintf						  (stdout, "\nLog L = ", res_three_LF[1][0], " with ", res_three_LF[1][1] + 9, " degrees of freedom, IC = ", getIC (res_three_LF[1][0], res_three_LF[1][1], sample_size), "\n");

lfOut	= csvFilePath + ".fit";
LIKELIHOOD_FUNCTION_OUTPUT = 7;
fprintf (lfOut, CLEAR_FILE, three_LF);
LIKELIHOOD_FUNCTION_OUTPUT = 2;

bsrel_bls = extractBranchLengthsFromBSREL ("mixtureTree");

tavl         = mixtureTree ^ 0;
renderString = PostOrderAVL2StringDistances (tavl, bsrel_bls);

fprintf (stdout, renderString, "\n");

UseModel (USE_NO_MODEL);
Tree T = renderString;

for	(k = 0; k < totalBranchCount; k += 1) {
    pValueByBranch[k][branch_length_column] = bsrel_bls[bNames[k]];
    

    ref = "mixtureTree."+bNames[k];
    
    rate_classes = rate_classes_by_branch[k];
    node_omegas = getNodeOmegaDistribution (ref, rate_classes);
    
    fprintf (stdout, "\n");
    printNodeDesc (ref, rate_classes_by_branch[k]);
    if (node_omegas[rate_classes-1][0] > 1) {
        pValueByBranch[k][2] = node_omegas[rate_classes-1][0];
        pValueByBranch[k][3] = node_omegas[rate_classes-1][1];
    }
    
        
    if (node_omegas[rate_classes-1][0] > 1 && node_omegas[rate_classes-1][1] > 1e-6 && selectedBranches[k])
    {
        fprintf (stdout, "...Testing for selection at this branch\n");
        _stashLF = saveLF ("three_LF");
        
        if (rate_classes > 1) {
            testing_parameter =   ref + ".omega" + rate_classes;
            *testing_parameter :< 1;
        } else {
            global _mg94omega = 1;
            _mg94omega :< 1;
            *(ref + ".nonsyn") := *(ref + ".syn") * _mg94omega;
        }
        
        Optimize					  (res_three_current,three_LF);
        
        fprintf (stdout, "\n");
        printNodeDesc (ref, rate_classes_by_branch[k]);
        pValueByBranch[k][lrt_column]			  = 2*(res_three_LF[1][0] - res_three_current[1][0]);				 
        pValueByBranch[k][p_uncorrected_column]			  = (1-CChi2 (pValueByBranch[k][lrt_column],1))*.5;
        fprintf (stdout, "p-value = ", pValueByBranch[k][p_uncorrected_column], "\n");
        
        if (rate_classes > 1) {
            *testing_parameter :< 10000;
        }
        
        if (pValueByBranch[k][lrt_column] < (-0.5)) {
            fprintf 					  (stdout, "[PHASE 2/REPEAT] Detected a convergence problem; refitting the LOCAL alternative model with new starting values\n");
            lfOut	= csvFilePath + ".fit";
            Optimize					  (res_three_LF,three_LF);
            LIKELIHOOD_FUNCTION_OUTPUT = 7;
            fprintf (lfOut, CLEAR_FILE, three_LF);
            LIKELIHOOD_FUNCTION_OUTPUT = 2;
            _stashLF = saveLF ("three_LF");
            k = 0;
        }
        else {
            _stashLF ["restoreLF"][""];
        }
    }
    else
    {
        pValueByBranch[k][p_uncorrected_column] = 1.0;
    }
}

OPTIMIZATION_METHOD = 4;


pValueSorter = {totalBranchCount,2};
pValueSorter = pValueSorter["_MATRIX_ELEMENT_ROW_*(_MATRIX_ELEMENT_COLUMN_==0)+pValueByBranch[_MATRIX_ELEMENT_ROW_][p_uncorrected_column]*(_MATRIX_ELEMENT_COLUMN_==1)"];
pValueSorter = pValueSorter % 1;
pValueSorter = pValueSorter["_MATRIX_ELEMENT_VALUE_*(_MATRIX_ELEMENT_COLUMN_==0)+_MATRIX_ELEMENT_VALUE_*(totalBranchCount-_MATRIX_ELEMENT_ROW_)*(_MATRIX_ELEMENT_COLUMN_==1)"];

fprintf (stdout,"\n\nSummary of branches under episodic selection (", Abs(selectedBranches)," were tested) :\n");
hasBranchesUnderSelection = 0;

pthreshold = 0.05;

for		(k = 0; k < totalBranchCount; k = k+1)
{
    pValueByBranch[pValueSorter[k][0]][p_corrected_column] = Min (1,pValueSorter[k][1]);
    if (pValueSorter[k][1] <= pthreshold)
    {
        fprintf (stdout, "\t", bNames[pValueSorter[k][0]], " p = ", pValueByBranch[pValueSorter[k][0]][p_corrected_column], "\n");
        hasBranchesUnderSelection += 1;
    }
}


if (hasBranchesUnderSelection == 0) {
    fprintf (stdout, "\tNo branches were found to be under selection at p <= ", pthreshold, "\n");
}


for		(k = 0; k < totalBranchCount; k = k+1) {
    fprintf (csvFilePath, "\n", bNames[k], ",", Join(",",pValueByBranch[k][-1]));
}


fprintf (csvFilePath, CLOSE_FILE);

//---- TREE RENDERING -----

LoadFunctionLibrary ("ReadDelimitedFiles");
treeString = Format (givenTree, 1,1);
resultTable=ReadCSVTable (csvFilePath, 2);

rows	= Rows (resultTable[2]);
columns = Columns(resultTable[1]);


TREE_OUTPUT_OPTIONS = {"__FONT_SIZE__":14};

tavl = T^0;
for (k = 1; k < Abs (tavl)-1; k+=1) {
	TREE_OUTPUT_OPTIONS[(tavl[k])["Name"]] = {};
}

for (k = 1; k < Abs (tavl)-1; k+=1) {
	thisP = (resultTable[1])[k-1][p_corrected_column];
	
    ref = "mixtureTree."+bNames[k-1];
    rate_classes = rate_classes_by_branch[k-1];
    myRateDistro = getNodeOmegaDistribution (ref, rate_classes);
    
    parentName = (tavl[((tavl[k])["Parent"])])["Name"];
		
	myRateDistro = myRateDistro % 0;
	
	colors = {rate_classes,4};
	for (cc = 0; cc < rate_classes; cc += 1) {
	    color1 = makeDNDSColor (myRateDistro[cc][0]);
	    colors [cc][0] = color1[0];
	    colors [cc][1] = color1[1];
	    colors [cc][2] = color1[2];
	    colors [cc][3] = myRateDistro[cc][1];
	}
	
	(TREE_OUTPUT_OPTIONS[(tavl[k])["Name"]])["TREE_OUTPUT_BRANCH_COLOR"] 		= colors;
	(TREE_OUTPUT_OPTIONS[(tavl[k])["Name"]])["TREE_OUTPUT_BRANCH_LINECAP"]  = 	thisP > 0.05;
	
	if (thisP <= 0.05)
	{
		(TREE_OUTPUT_OPTIONS[(tavl[k])["Name"]])["TREE_OUTPUT_BRANCH_THICKNESS"] = 	14;

	}
	if (Abs((tavl[k])["Children"]) > 0)
	{
		(TREE_OUTPUT_OPTIONS[(tavl[k])["Name"]])["TREE_OUTPUT_BRANCH_TLABEL"] = (tavl[k])["Name"]; 
	}
}

height = TipCount (T) * 36;
psTree = PSTreeString (T,"STRING_SUPPLIED_LENGTHS",{{400,height}});

treePath = csvFilePath + ".ps";

fprintf (treePath, CLEAR_FILE, psTree);

return pValueByBranch;

//------------------------------------------------------------------------------------------------------------------------
function makeDNDSColor (omega) {
	if (omega < 1) {
		return {{0.5*omega__,0.5*omega__,1-0.5*omega__}};
	}
	omega = Min (omega,5);
	return {{0.5+0.125*(omega__-1),0.5-(omega__-1)*0.125,0.5-(omega__-1)*0.125}};
}

//------------------------------------------------------------------------------------------------------------------------
function saveLF (ID) {
	ExecuteCommands ("GetString(_lfInfo,"+ID+",-1)");
	_stashLF = {};
	for (_k = 0; _k < Columns (_lfInfo["Global Independent"]); _k+=1)
	{
		_stashLF [(_lfInfo["Global Independent"])[_k]] = Eval ((_lfInfo["Global Independent"])[_k]);
	}
	for (_k = 0; _k < Columns (_lfInfo["Local Independent"]); _k+=1)
	{
		_stashLF [(_lfInfo["Local Independent"])[_k]] = Eval ((_lfInfo["Local Independent"])[_k]);
	}
	
	return _stashLF;
}

//------------------------------------------------------------------------------------------------------------------------

function restoreLF (key, value) {
	ExecuteCommands (key + " = " + value);
	return 0;
}


//------------------------------------------------------------------------------------------------------------------------

function gEval (expr) {
    return Eval (expr);
}


//------------------------------------------------------------------------------------------------------------------------

lfunction printNodeDesc (ref, rate_classes) {
    if (rate_classes > 1) {
        wts = generateFrequencyParameterization ("Paux", rate_classes);
     
        fprintf (stdout, "Node: ", ref, 
                     "\n\tLength parameter = ", ^(ref+".t"));
                 
        for (k = 1; k < rate_classes; k+=1) {
            *("Paux" + k) = ^(ref+".Paux" + k);
        }

        for (k = 0; k < rate_classes; k+=1) {
            fprintf (stdout, "\n\tClass ", k+1,": omega = ", ^(ref+".omega" + (k+1)), ", weight = ", Eval(wts[k]));
        }
                  
        fprintf (stdout, "\n"); 
    } else {
        fprintf (stdout, "Node: ", ref, 
                     "\n\tLength parameter = ", ^(ref+".syn"));
        
        fprintf (stdout, "\n\tClass ", rate_classes ,": omega = ", _standardizeRatio(^(ref+".nonsyn"),^(ref+".syn")), ", weight = 1.\n");
    }    
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------

function getIC (logl,params,samples) {
    return -2*logl + 2*samples/(samples-params-1)*params;
    //return -2*logl + Log (samples)*params;
}

//------------------------------------------------------------------------------------------------------------------------

function copyomegas (node, syn, nonsyn) {
    local_omega = =  Min (10, nonsyn/syn);
    Eval("`node`.omega2 = 1");
    if (local_omega < 1) {
        Eval("`node`.Paux1 = 1");
        Eval("`node`.omega1 =  local_omega");
        Eval("`node`.omega2 = 10");
    } else {
        Eval("`node`.Paux1 = 0");
        Eval("`node`.omega2 =  local_omega");
        Eval("`node`.omega1 = 0");    
    }
    Eval("`node`.t = syn");
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------

lfunction setBSRELBounds (node, count) {
    for (k = 1; k < count; k += 1) {
        ^("`node`.Paux" + k ) :< 1;   
        ^("`node`.omega" + k ) :<1;   
    }
    return 0;
}


//------------------------------------------------------------------------------------------------------------------------

lfunction getNodeOmegaDistribution (ref, rate_classes) {
    values = {rate_classes, 2};
    if (rate_classes > 1) {
        wts = generateFrequencyParameterization ("Paux", rate_classes);
        for (k = 0; k < rate_classes; k+=1) {
            values[k][0] = ^(ref+".omega" + (k+1));
            values[k][1] = gEval(wts[k]);
        }
    } else {
        values [0][0] =   ^("`ref`.nonsyn")/^("`ref`.syn");
        values [0][1] =   1;        
    }
    return values;
}

//------------------------------------------------------------------------------------------------------------------------

lfunction saveNodeMLES (node, count) {
    values = {count, 2};
    for (k = 0; k < count - 1; k += 1) {
        values [k][0] = ^("`node`.omega" + (k+1));
        values [k][1] = ^("`node`.Paux" + (k+1));        
    }
    values [count-1][0] = ^("`node`.omega" + count);
    values [count-1][1] = ^("`node`.t");
    return values;
}

//------------------------------------------------------------------------------------------------------------------------

lfunction initNodeMLESPlus1 (node, values) {
    count = Rows (values);
    for (k = 0; k < count - 1; k += 1) {
        ^("`node`.omega" + (k+1)) = values[k][0];
        ^("`node`.Paux"  + (k+1)) = values[k][1];        
    }
    if (values[count-1][0] > 1) {
        ^("`node`.omega" + (count+1)) = values[count-1][0];
        ^("`node`.Paux"  + (count))   = 0.01;        
        ^("`node`.omega" + (count))   = 0.5;
    } else {
        ^("`node`.omega" + (count+1)) = 2;
        ^("`node`.Paux"  + (count))   = 0.99;        
        ^("`node`.omega" + (count))   = values[count-1][0];
    }
    ^("`node`.t") = values[count-1][1];  
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------

lfunction restoreNodeMLE (node, values) {
    count = Rows (values);
    for (k = 0; k < count - 1; k += 1) {
         ^("`node`.omega" + (k+1)) = values[k][0];
         ^("`node`.Paux"  + (k+1)) = values[k][1];        
    }
    ^("`node`.omega" + count) = values[count-1][0];  
    ^("`node`.t") = values[count-1][1];  
    return 0;
}


//------------------------------------------------------------------------------------------------------------------------

lfunction generateFrequencyParameterization (prefix, count) {
    result = {1, count};
    if (count > 1) {
        result[0] = "`prefix`1";
        buffer = "(1-`prefix`1)";
        for (i = 1; i < count - 1; i+=1) {
            result[i] = buffer + "*" + "`prefix`" + (i+1);
            buffer   += "*(1-`prefix`" + (i+1) + ")";
        } 
        result[count-1] = buffer;
    } else {
        result [0] = "1";
    }
    return result;
}
