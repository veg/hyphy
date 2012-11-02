skipCodeSelectionStep 		= 0;


LoadFunctionLibrary("chooseGeneticCode");
LoadFunctionLibrary("GrabBag");
LoadFunctionLibrary("dSdNTreeTools");
LoadFunctionLibrary("CF3x4");
LoadFunctionLibrary("BranchSiteTemplate");
LoadFunctionLibrary("DescriptiveStatistics");
LoadFunctionLibrary("TreeTools");


DataSet 			ds 				= ReadDataFile(PROMPT_FOR_FILE);
DataSetFilter 		dsf 			= CreateFilter(ds,3,"","",GeneticCodeExclusions);

HarvestFrequencies	(nuc3, dsf, 3, 1, 1);
nucCF			  = CF3x4	(nuc3, GeneticCodeExclusions);

rate_class_count = prompt_for_a_value ("How many rate classes", 3, 2, 8, 1);
fprintf (stdout, "\nUsing ", rate_class_count, " rate classes\n");

PopulateModelMatrix			  ("MGMatrixLocal",  nucCF, "syn", "", "nonsyn");
codon3x4					= BuildCodonFrequencies (nucCF);
Model		MGL				= (MGMatrixLocal, codon3x4, 0);

_DO_TREE_REBALANCE_ = 0;
LoadFunctionLibrary			  ("queryTree");

SetDialogPrompt ("Save analysis results to");
fprintf (PROMPT_FOR_FILE, CLEAR_FILE, KEEP_OPEN,"Model\tLogL\tNP\tBIC\tTree String");
csvFilePath = LAST_FILE_PATH;


// ---- PHASE 0 : LOCAL (no omega variation model) FIT. -----

fprintf 					  (stdout, "[BS-REL PHASE 0] Fitting the local MG94 (no site-to-site variation) to obtain initial parameter estimates\n");


LikelihoodFunction	base_LF	 = (dsf, givenTree);

/*
ExecuteAFile (csvFilePath + ".mglocal.fit");
res_base = {2,2};
*/

Optimize					  (res_base,base_LF);
writeTheLF (".mglocal.fit", "base_LF");
baseParameters                = 9;
localLL						 = res_base[1][0];
localParams					 = res_base[1][1] + baseParameters;

totalBranchCount			 = BranchCount(givenTree) + TipCount (givenTree);

pValueByBranch				  = {totalBranchCount,10};
bNames						  = BranchName (givenTree, -1);

for (k = 0; k < totalBranchCount; k = k+1) {
	srate  = Eval ("givenTree." + bNames[k] + ".syn");
	nsrate = Eval ("givenTree." + bNames[k] + ".nonsyn");
	if (srate > 0) {
		pValueByBranch [k][0] = Min (10, nsrate/srate);
	} else {
		pValueByBranch [k][0] = 10;
	}	
}

sample_size                    = dsf.sites * dsf.species;
omegaStats					 = GatherDescriptiveStats (pValueByBranch[-1][0]);
fprintf						 (stdout, "Local omega model (no rate variation) ");


localBIC                    = BIC (localLL, localParams, sample_size);
PrintDescriptiveStats		 ("Branch omega values", omegaStats);

for (rcc = 1; rcc <= rate_class_count; rcc += 1) {
    PopulateModelMatrix			  ("MGMatrix"+rcc,  nucCF, "t", "omega"+rcc, "");
    ExecuteCommands               ("global omegaG" + rcc + " = " + (0.3*rcc - 0.2));
    PopulateModelMatrix			  ("MGMatrixG"+rcc,  nucCF, "t", "omegaG"+rcc, "");
}
    
// ---- PHASE 1 : FULL GLOBAL MODEL FIT ---- //

generate_bs_rel_gdd_freqs (rate_class_count, "global_freqsG", "model_mixingG", "PauxG", "MGMatrixG",1);
generate_bs_rel_gdd_freqs (rate_class_count, "global_freqs",  "model_mixing",  "Paux", "MGMatrix",0);

sortedLocalDNDS              = (pValueByBranch[-1][0])%0;

for (rcc = 1; rcc <= rate_class_count; rcc += 1) {
    ExecuteCommands ("omegaG" + rcc + "=" +  sortedLocalDNDS[Max(0,totalBranchCount*rcc$rate_class_count-1)]);
}


ExecuteCommands ("Model 		MGG	=(\"`model_mixingG`\",codon3x4,EXPLICIT_FORM_MATRIX_EXPONENTIAL);");
Tree						   mixtureTreeG = treeString;

ExecuteCommands ("Model 		MGLocal=(\"`model_mixing`\",codon3x4,EXPLICIT_FORM_MATRIX_EXPONENTIAL);");
Tree						   mixtureTree = treeString;

ReplicateConstraint 		  ("this1.?.t:=this2.?.syn",mixtureTreeG,givenTree);
ReplicateConstraint 		  ("this1.?.t:=this2.?.syn",mixtureTree,givenTree);
ClearConstraints			  (mixtureTree);
ClearConstraints			  (mixtureTreeG);


ASSUME_REVERSIBLE_MODELS	  = 1;
OPTIMIZATION_METHOD           = 4;

BIC_scores                    = {};
BranchLengthEstimates         = {};

mixTreeAVL                    = mixtureTreeG ^ 0;

LikelihoodFunction three_LFG   = (dsf,mixtureTreeG);
fprintf 					  (stdout, "[BS-REL PHASE 1] Fitting a GLOBAL branch-site matrix mixture\n");


blG = getBranchLengthExpression (1);
removeT = "/(" + (blG ^{{"\\*t",""}}) + ")";

fprintf (stdout, base_LF);

blLocal = {};
for (k = 0; k < totalBranchCount; k += 1) {
    blLocal [bNames[k]] = Eval("BranchLength(givenTree,\""+bNames[k]+"\")")*3;
    ExecuteCommands ("mixtureTreeG." + bNames[k] + ".t := " + blLocal [bNames[k]] + removeT );    
}


USE_LAST_RESULTS			= 1;
SHORT_MPI_RETURN            = 1;
VERBOSITY_LEVEL				= 0;

OPTIMIZATION_PRECISION      = 0.1;

Optimize					  (res_three_LF_global,three_LFG);

fprintf (stdout, "\n\nApproximation phase 1 (use MG-local branch lengths): ", res_three_LF_global[1][0], "\n");
fprintf     (stdout,"Tree branch lengths:\n",PostOrderAVL2StringDistances(mixTreeAVL,getBranchLengths ("mixtureTreeG", 1)));

OPTIMIZATION_PRECISION      = 0.01;

for (k = 0; k < totalBranchCount; k += 1) {
    ExecuteCommands ("mixtureTreeG." + bNames[k] + ".t = mixtureTreeG." + bNames[k] + ".t");    
    bl = blLocal [bNames[k]] ;
    ExecuteCommands ("FindRoot (z,`blG`-"+bl+",t,0,10000);");
    ExecuteCommands ("mixtureTreeG." + bNames[k] + ".t :< " + 10*z);    

}

Optimize					  (res_three_LF_global,three_LFG);
fprintf (stdout, "\n\nApproximation phase 2 (maximum branch lengths are limited to 10x those of the null model): ", res_three_LF_global[1][0], "\n");
fprintf     (stdout,"Tree branch lengths:\n",PostOrderAVL2StringDistances(mixTreeAVL,getBranchLengths ("mixtureTreeG", 1)));

OPTIMIZATION_PRECISION      = 0.001;

for (k = 0; k < totalBranchCount; k += 1) {
    ExecuteCommands ("mixtureTreeG." + bNames[k] + ".t :< 10000");    
}

Optimize					  (res_three_LF_global,three_LFG);

//ExecuteAFile (csvFilePath + ".relglobal.fit");


blStashByName = {};
for (k = 0; k < totalBranchCount; k += 1) {
    bl = blLocal [bNames[k]] ;
    ExecuteCommands ("FindRoot (z,`blG`-"+bl+",t,0,10000);");
    blStashByName[bNames[k]] = z;
    //fprintf (stdout, "\n", bNames[k], " : ", z, " (", bl/3, ")\n");
}

writeTheLF (".relglobal.fit", "three_LFG");


fprintf (stdout, "Global model fit:");

BIC_scores ["global model"] = BIC (res_three_LF_global[1][0], res_three_LF_global[1][1]+baseParameters, sample_size);
BranchLengthEstimates ["global model"] = getBranchLengths ("mixtureTreeG", 1);

fprintf     (csvFilePath, "\nGlobal\t", res_three_LF_global[1][0],"\t", res_three_LF_global[1][1]+baseParameters, "\t", BIC_scores ["global model"], "\t", PostOrderAVL2StringDistances(mixTreeAVL,BranchLengthEstimates ["global model"]));

fprintf (stdout, "\nInferred this global omega distribution: ");
reportOmegaDistro ("");
fprintf     (stdout,"Tree branch lengths:\n",PostOrderAVL2StringDistances(mixTreeAVL,BranchLengthEstimates ["global model"]));
globalState = saveLF ("three_LFG");


for (k = 0; k < totalBranchCount; k += 1) {
    constrainABranch (bNames[k], 0);
}

// ---- PHASE 2 : LOCAL MODEL FITS ---- //

fprintf (stdout, "\n[BS-REL PHASE 2] Fitting a GLOBAL branch-site matrix mixture with a SINGLE unconstrained branch\n");

OPTIMIZATION_METHOD = 0;

LikelihoodFunction three_LF   = (dsf,mixtureTree);

branchValues = {};

if (MPI_NODE_COUNT > 1) {
    MPI_NODE_STATE = {MPI_NODE_COUNT-1,1};
    MPI_NODE_STATE[0] = "";
}


for (k = 0; k < totalBranchCount; k+=1) {
    thisBranchName = bNames[k];
    fprintf (stdout, "\n[BS-REL PHASE 2. Branch '", thisBranchName, "']\n");
    globalState["restoreLF"][""];
    if (k) {
        if (MPI_NODE_COUNT > 1) {
            for (nodeID = 0; nodeID < totalBranchCount; nodeID += 1) {
                constrainABranch (bNames[nodeID],0);
            }
        } else {
            constrainABranch (bNames[k-1],0);    
        }
    }
    unConstrainABranch (thisBranchName);
    
    if (MPI_NODE_COUNT > 1) {
        for (nodeID = 0; nodeID < MPI_NODE_COUNT-1; nodeID += 1) {
            if (Abs(MPI_NODE_STATE[nodeID]) == 0) {
                MPISend (nodeID+1, three_LF);
                fprintf (stdout, "\n[SENT TO NODE ", nodeID, "]\n");
                MPI_NODE_STATE [nodeID] = thisBranchName;
                break;
            }   
        }
        if (nodeID == MPI_NODE_COUNT-1) {
            processABranch (thisBranchName, 1);
        }
    } else {
        OPTIMIZATION_PRECISION = 0.05;
        maxBL = blStashByName[thisBranchName]*10;
        t = maxBL;
        fprintf (stdout, thisBranchName, ": maximum t value: ", maxBL, ". Maximum achieved branch length: ", Eval (blG), "\n");
        constrainABranch (thisBranchName,maxBL);
        Optimize (localBranchRes, three_LF);
        fprintf (stdout, "Approximation phase 1 (constrained branch length): ", localBranchRes[1][0], "\n");
        constrainABranch (thisBranchName,-1);
        OPTIMIZATION_PRECISION = 0.001;
        Optimize (localBranchRes, three_LF);
        processABranch (thisBranchName,0);
    }
}


leftOver = 0;
for (nodeID = 0; nodeID < MPI_NODE_COUNT-1; nodeID += 1) {
    leftOver += Abs(MPI_NODE_STATE[nodeID])>0;
}

//fprintf (stdout, MPI_NODE_STATE, "\n", leftOver, "\n");


for (nodeID = 0; nodeID < leftOver; nodeID += 1) {
    processABranch ("", 0);
}

branchValues ["restoreLF"][""];

fprintf (stdout, "\n[BS-REL PHASE 3] Fitting a LOCAL branch-site matrix mixture model\n");
Optimize (res_local, three_LF);

writeTheLF (".local.fit", "three_LF");
 
fprintf (stdout, "Local model fit:");   
BIC_scores ["local model"] = BIC (res_local[1][0], res_local[1][1]+baseParameters, sample_size);
BranchLengthEstimates ["local model"] = getBranchLengths ("mixtureTree", 0);
fprintf     (stdout, "\nTree: ", PostOrderAVL2StringDistances(mixTreeAVL,BranchLengthEstimates ["local model"]), "\n");
fprintf     (csvFilePath, "\nLocal\t", res_local[1][0],"\t", res_local[1][1]+baseParameters, "\t", BIC_scores ["local model"], "\t", PostOrderAVL2StringDistances(mixTreeAVL,BranchLengthEstimates ["local model"]), CLOSE_FILE);

min_BIC = 1e100;

function compute_min_BIC (key, value) {
    min_BIC = Min (min_BIC, value);
    return 0;
}

function akaike_weights_pass1 (key, value) {
    akaike_scores [key] = Exp (0.5*(min_BIC-value));
    return 0;
}


function akaike_weights_pass2 (key, value) {
    BIC_scores [key] = value/norm_factor;
    return 0;
}

fprintf (stdout, "\n[BS-REL PHASE 4] Generating a model averaged branch lengths\n");
akaike_scores = {};
BIC_scores ["compute_min_BIC"][""];
BIC_scores ["akaike_weights_pass1"][""];
norm_factor = +akaike_scores;
akaike_scores ["akaike_weights_pass2"][""];

reweighted_branch_lengths = {};

model_names = Rows(BranchLengthEstimates);
for (k = 0; k < Abs (BranchLengthEstimates); k+=1) {
    key = model_names [k];
    bl_avl = BranchLengthEstimates[key];
    bnames = Rows (bl_avl);
    
    for (b = 0; b < Abs (bl_avl); b += 1) {
        reweighted_branch_lengths[bnames[b]] += BIC_scores[key]*bl_avl[bnames[b]];
    }
}

fprintf (stdout,  PostOrderAVL2StringDistances(mixTreeAVL,reweighted_branch_lengths), "\n");

//------------------------------------------------------------------------------------------------------------------------
function processABranch (thisBranchName, doSend) {
    if (MPI_NODE_COUNT > 1) {
        MPIReceive (-1,fromNode,resStr);
	    prevBranch = MPI_NODE_STATE[fromNode-1];
        fprintf (stdout, "\n[RECEIVED " + prevBranch +" FROM NODE ", (fromNode-1), "]\n");
        if (doSend) {
            MPISend (fromNode, three_LF);
            MPI_NODE_STATE [fromNode-1] = thisBranchName;
            fprintf (stdout, "\n[SENT " + thisBranchName +" TO NODE ", (fromNode-1), "]\n");
        } else {
            MPI_NODE_STATE [fromNode-1] = "";
        }
        thisBranchName = prevBranch;
        ExecuteCommands (resStr);
        three_LF_MLE_VALUES ["restoreLF"][""];
        localBranchRes = three_LF_MLES;
    }
    writeTheLF (".`thisBranchName`.fit", "three_LF");
    fprintf (stdout, "\nModel fit:");   
    BIC_scores [thisBranchName] = BIC (localBranchRes[1][0], localBranchRes[1][1]+baseParameters, sample_size);
    BranchLengthEstimates [thisBranchName] = getBranchLengths ("mixtureTree", 0);
    fprintf     (csvFilePath, "\nGlobal+",thisBranchName,"\t", localBranchRes[1][0],"\t", localBranchRes[1][1]+baseParameters, "\t", BIC_scores [thisBranchName], "\t", PostOrderAVL2StringDistances(mixTreeAVL,BranchLengthEstimates [thisBranchName]));
    fprintf (stdout, "\nLocal branch omega distribution: ");
    reportOmegaDistro(thisBranchName);
    fprintf (stdout, "Global omega distribution on the rest of the branches: ");
    reportOmegaDistro ("");
    stashBranchValues (thisBranchName, "branchValues");
    pv = 1-CChi2 (2*(localBranchRes[1][0]-res_three_LF_global[1][0]),5);
    fprintf (stdout, "\nLRT p-value for branch deviation from the global pattern = ", pv, "\n");
    fprintf     (stdout,"\n", PostOrderAVL2StringDistances(mixTreeAVL,BranchLengthEstimates [thisBranchName]));
    return 0;

}

//------------------------------------------------------------------------------------------------------------------------
function writeTheLF (fileNameExt,lfID) {
    lfOut	= csvFilePath + fileNameExt;//".local.fit";    
    LIKELIHOOD_FUNCTION_OUTPUT = 7;
    ExecuteCommands ("fprintf (lfOut, CLEAR_FILE, `lfID`);");
    fprintf (stdout, "[WROTE MODEL FIT TO ", lfOut, "]\n");
    LIKELIHOOD_FUNCTION_OUTPUT = 2;
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------
function constrainABranch (branch_name, z) {
    if (z > 0) {
        ExecuteCommands ("mixtureTree." + branch_name + ".t :< " + z);
        return 0;
    } else {
        if (z < 0) {
            ExecuteCommands ("mixtureTree." + branch_name + ".t :< 10000");
            return 0;
        }
    }
    ExecuteCommands ("mixtureTree." + branch_name + ".t = mixtureTreeG." + branch_name + ".t");
    for (rcc = 1; rcc <= rate_class_count; rcc += 1) {
        ExecuteCommands ("mixtureTree." + branch_name + ".omega" + rcc + ":= omegaG" + rcc);
    }        

    for (rcc = 1; rcc < rate_class_count; rcc += 1) {
        ExecuteCommands ("mixtureTree." + branch_name + ".Paux" + rcc + ":= PauxG" + rcc);
    }        
        
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------
function stashBranchValues (branch_name, storage&) {

    locals = {"0" : "t"};
    
    for (rcc = 1; rcc <= rate_class_count; rcc += 1) {
        locals + ("omega" + rcc);
        if (rcc < rate_class_count) {
            locals + ("Paux" + rcc);
        }
    }
    for (_varID = 0; _varID < Abs (locals); _varID += 1) {
        storage ["mixtureTree." + branch_name + "." + locals[_varID]] = Eval ("mixtureTree." + branch_name + "." + locals[_varID]);
    }
    
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------
function unConstrainABranch (branch_name) {

    ExecuteCommands ("mixtureTree." + branch_name + ".t = mixtureTreeG." + branch_name + ".t");
    for (rcc = 1; rcc <= rate_class_count; rcc += 1) {
        ExecuteCommands ("mixtureTree." + branch_name + ".omega" + rcc + "= omegaG" + rcc);
    }        

    for (rcc = 1; rcc < rate_class_count; rcc += 1) {
        ExecuteCommands ("mixtureTree." + branch_name + ".Paux" + rcc + "= PauxG" + rcc);
    }        

    return 0;
}

//------------------------------------------------------------------------------------------------------------------------
function BIC (logL, params, sample_size) {
    BICvalue = -2*logL + 2*params*(sample_size-params-1);
    fprintf (stdout, " log(L) = ", logL, " with ", params, " parameters, yielding BIC = ", BICvalue, " assuming the sample size of ", sample_size, ".");
    return BICvalue;
}

//------------------------------------------------------------------------------------------------------------------------
function reportOmegaDistro (branchName) {
    if (Abs(branchName)==0) {
        for (rcc = 1; rcc <= rate_class_count; rcc += 1) {
            fprintf (stdout, "\n\t omega = ", Eval("omegaG"+rcc), ", p = ", Eval(global_freqsG[rcc-1]));
        }
    } else {
        for (rcc = 1; rcc < rate_class_count; rcc += 1) {
            ExecuteCommands ("Paux" + rcc + "=mixtureTree.`branchName`.Paux" +rcc);
        }        

        for (rcc = 1; rcc <= rate_class_count; rcc += 1) {
            fprintf (stdout, "\n\t omega = ", Eval("mixtureTree.`branchName`.omega"+rcc), ", p = ", Eval(global_freqs[rcc-1]));
        }    
    }
    fprintf (stdout, "\n");
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------

function getBranchLengthExpression (modelType) {
    bl_expression = ""; bl_expression * 128;
    
    if (modelType == 0) {
        modelPrefix = "MGMatrix";
        freqMatrix  = global_freqs;
    } else {
        modelPrefix = "MGMatrixG";
        freqMatrix  = global_freqsG;
    }
    for (rcc = 1; rcc <= rate_class_count; rcc += 1) {
        if (rcc > 1) {
            bl_expression * "+";
        }
        ExecuteCommands ("Model _temp = (`modelPrefix`" + rcc + ", codon3x4, 0);");
        GetString (modelBL, _temp, -1);
        bl_expression * (freqMatrix[rcc-1]+"*(`modelBL`)");
    }   
    bl_expression * 0;
    return bl_expression;
}

//------------------------------------------------------------------------------------------------------------------------

function getBranchLengths (treeID, modelType) {
    blByName        = {};
    bnames_list     = Eval ("BranchName (`treeID`,-1)");

    bl_expression = getBranchLengthExpression(modelType);
 
    if (modelType == 0) {
        locals = {"0" : "t"};
        
        for (rcc = 1; rcc <= rate_class_count; rcc += 1) {
            locals + ("omega" + rcc);
            if (rcc < rate_class_count) {
                locals + ("Paux" + rcc);
            }
        }
      
        for (_bID = 0; _bID < Columns (bnames_list) - 1; _bID += 1) {
            branch_name = bnames_list[_bID];
         
             for (_varID = 0; _varID < Abs (locals); _varID += 1) {
                ExecuteCommands (locals[_varID] + " = Eval (\"`treeID`.`branch_name`." + locals[_varID] + "\")");
            }
            blByName [branch_name] = Eval (bl_expression)/3;
        }
    } else {
        for (_bID = 0; _bID < Columns (bnames_list) - 1; _bID += 1) {
            branch_name = bnames_list[_bID];
            t = Eval ("`treeID`.`branch_name`.t");        
            blByName [branch_name] = Eval (bl_expression)/3;
        }    
    }
    
    return blByName;
}


//------------------------------------------------------------------------------------------------------------------------
function saveLF (ID)
{
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
function restoreLF (key, value)
{
	ExecuteCommands (key + " = " + value);
	return 0;
}


//------------------------------------------------------------------------------------------------------------------------
function generate_bs_rel_gdd_freqs (numberOfRates, freqs&, model_mixing&, probPrefix, modelPrefix, is_global)
{

	freqs    	 = {numberOfRates,1};
	model_mixing	 = ""; model_mixing * 128; 
	
	if (numberOfRates == 1) {
		freqs[0] = "1";
	}
	else {
	    global_prefix = "";
	    
	    if (is_global) 
	        { global_prefix = "global "; }
	    
        for (mi=1; mi<numberOfRates; mi += 1) {
            ExecuteCommands (global_prefix+probPrefix+mi+":<1;"+probPrefix+mi+" = 1/" + (numberOfRates-mi+1));
            ExecuteCommands (global_prefix+probPrefix+mi+":>0;");
        }
		
		freqs[0] 	 = ""+probPrefix+"1";
		for (mi=1; mi<numberOfRates-1; mi+=1) {
			freqs[mi] = "";
			for (mi2=1;mi2<=mi;mi2+=1) {
				freqs[mi] = freqs[mi] + "(1-"+probPrefix+mi2+")*";		
			}
			freqs[mi] = freqs[mi] + probPrefix+(mi+1);	
		}	
	
		freqs[mi] = "";
		for (mi2=1;mi2<mi;mi2+=1)
		{
			freqs[mi] = freqs[mi] + "(1-"+probPrefix+mi2+")*";		
		}
		freqs[mi] = freqs[mi] + "(1-"+probPrefix+mi+")";	
	}
	
	model_mixing * ("Exp(`modelPrefix`1)*"+freqs[0]);
	for (mi = 1; mi < numberOfRates; mi=mi+1) {
		model_mixing * ("+Exp(`modelPrefix`"+(mi+1)+")*" + freqs[mi]);
	}
	model_mixing * 0;
	return 0;
}
