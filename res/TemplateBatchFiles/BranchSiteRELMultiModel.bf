skipCodeSelectionStep 		= 0;
VERBOSITY_LEVEL				= 0;


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



omega1 = 0.2;
omega2 = 0.5;
omega3 = 1.0;

PopulateModelMatrix			  ("MGMatrix1",  nucCF, "t", "omega1", "");
PopulateModelMatrix			  ("MGMatrix2",  nucCF, "t", "omega2", "");
PopulateModelMatrix			  ("MGMatrix3",  nucCF, "t", "omega3", "");

global	                        omegaG1 = 0.2; omegaG1 :< 1;
global	                        omegaG2 = 0.5; omegaG2 :< 1;
global	                        omegaG3 = 2.0;

PopulateModelMatrix			  ("MGMatrix1G",  nucCF, "t", "omegaG1", "");
PopulateModelMatrix			  ("MGMatrix2G",  nucCF, "t", "omegaG2", "");
PopulateModelMatrix			  ("MGMatrix3G",  nucCF, "t", "omegaG3", "");

PopulateModelMatrix			  ("MGMatrixLocal",  nucCF, "syn", "", "nonsyn");

codon3x4					= BuildCodonFrequencies (nucCF);
Model		MGL				= (MGMatrixLocal, codon3x4, 0);

LoadFunctionLibrary			  ("queryTree");

SetDialogPrompt ("Save analysis results to");
fprintf (PROMPT_FOR_FILE, CLEAR_FILE, KEEP_OPEN,"Model\tLogL\tNP\tBIC\tTree String");
csvFilePath = LAST_FILE_PATH;


// ---- PHASE 0 : LOCAL (no omega variation model) FIT. ----- //

fprintf 					  (stdout, "[BS-REL PHASE 0] Fitting the local MG94 (no site-to-site variation) to obtain initial parameter estimates\n");


LikelihoodFunction	base_LF	 = (dsf, givenTree);
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

// ---- PHASE 1 : FULL GLOBAL MODEL FIT ---- //

omegaG1						 = omegaStats["2.5%"];
omegaG2						 = omegaStats["Median"];
omegaG3						 = omegaStats["97.5%"];

Paux1 						 = 0.3;
Paux1 						 :< 1;
Paux2 						 = 0.4;
Paux2 						 :< 1;

global Paux1G 				  = 0.3;
global Paux2G 				  = 0.4;


Model 		MGG		=		  ("Exp(MGMatrix1G)*Paux1G+Exp(MGMatrix2G)*(1-Paux1G)*Paux2G+Exp(MGMatrix3G)*(1-Paux1G)*(1-Paux2G)",codon3x4,EXPLICIT_FORM_MATRIX_EXPONENTIAL);
Tree						   mixtureTreeG = treeString;

Model 		MG1		=		  ("Exp(MGMatrix1)*Paux1+Exp(MGMatrix2)*(1-Paux1)*Paux2+Exp(MGMatrix3)*(1-Paux1)*(1-Paux2)",codon3x4,EXPLICIT_FORM_MATRIX_EXPONENTIAL);
Tree						   mixtureTree = treeString;



ReplicateConstraint 		  ("this1.?.t:=this2.?.syn",mixtureTreeG,givenTree);
ReplicateConstraint 		  ("this1.?.t:=this2.?.syn",mixtureTree,givenTree);

ClearConstraints			  (mixtureTree);
ClearConstraints			  (mixtureTreeG);

omegaG1						 :< 1;
omegaG2						 :< 1;
Paux1G 						 :< 1;
Paux2G 						 :< 1;

ASSUME_REVERSIBLE_MODELS	  = 1;
OPTIMIZATION_METHOD           = 0;

BIC_scores                    = {};
BranchLengthEstimates         = {};

mixTreeAVL                    = mixtureTreeG ^ 0;

LikelihoodFunction three_LFG   = (dsf,mixtureTreeG);
fprintf 					  (stdout, "[BS-REL PHASE 1] Fitting a GLOBAL branch-site matrix mixture\n");

for (k = 0; k < totalBranchCount; k += 1) {
    if (k == 0) {
        expr            = Eval("BranchLength(givenTree,\""+bNames[0]+";EXPECTED_NUMBER_OF_SUBSTITUTIONS\")");
        syn             = 1; nonsyn = 0;
        synM            = Eval(expr);
        syn             = 0; nonsyn = 1;
        nonsynM         = Eval(expr);
    }
    
 	srate  = Eval ("givenTree." + bNames[k] + ".syn");
	nsrate = Eval ("givenTree." + bNames[k] + ".nonsyn");
    bl = Eval("BranchLength(givenTree,\""+bNames[k]+"\")")*3;
    
    if (srate > 0) {
        baseOmega = nsrate/srate;
    } else {
        baseOmega = 10000;
    }
        
    bl = bl / (synM + nonsynM * baseOmega);  
    ExecuteCommands ("mixtureTreeG." + bNames[k] + ".t = bl");
}

USE_LAST_RESULTS			  = 1;

SHORT_MPI_RETURN = 1;

Optimize					  (res_three_LF_global,three_LFG);
baseParameters                += 5; // to count nucleotide rates

writeTheLF (".relglobal.fit", "three_LFG");

fprintf (stdout, "Global model fit:");

BIC_scores ["global model"] = BIC (res_three_LF_global[1][0], res_three_LF_global[1][1]+baseParameters, sample_size);
BranchLengthEstimates ["global model"] = getBranchLengths ("mixtureTreeG", 1);

fprintf     (csvFilePath, "\nGlobal\t", res_three_LF_global[1][0],"\t", res_three_LF_global[1][1]+baseParameters, "\t", BIC_scores ["global model"], "\t", PostOrderAVL2StringDistances(mixTreeAVL,BranchLengthEstimates ["global model"]));

fprintf (stdout, "\nInferred this global omega distribution: ");
reportOmegaDistro (omegaG1,omegaG2,omegaG3,Paux1G,Paux2G);
globalState = saveLF ("three_LFG");

for (k = 0; k < totalBranchCount; k += 1) {
    constrainABranch (bNames[k]);
}

// ---- PHASE 2 : LOCAL MODEL FITS ---- //

fprintf (stdout, "\n[BS-REL PHASE 2] Fitting a GLOBAL branch-site matrix mixture with a SINGLE unconstrained branch\n");


LikelihoodFunction three_LF   = (dsf,mixtureTree);

branchValues = {};

if (MPI_NODE_COUNT > 1) {
    MPI_NODE_STATE = {MPI_NODE_COUNT-1,1};
    MPI_NODE_STATE[0] = "";
}

for (k = 0; k < totalBranchCount; k+=1) {
    fprintf (stdout, "\n[BS-REL PHASE 2. Branch '", bNames[k], "']\n");
    thisBranchName = bNames[k];
    globalState["restoreLF"][""];
    if (k) {
        if (MPI_NODE_COUNT > 1) {
            for (nodeID = 0; nodeID < totalBranchCount; nodeID += 1) {
                constrainABranch (bNames[nodeID]);
            }
        } else {
            constrainABranch (bNames[k-1]);    
        }
    }
    unConstrainABranch (bNames[k]);
    
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

fprintf (stdout, "\n[BS-REL PHASE 3] Fitting a LOCAL branch-site matrix mixture with a SINGLE unconstrained branch\n");
Optimize (res_local, three_LF);

writeTheLF (".local.fit", "three_LF");
 
fprintf (stdout, "Local model fit:");   
BIC_scores ["local model"] = BIC (res_local[1][0], res_local[1][1]+baseParameters, sample_size);
BranchLengthEstimates ["local model"] = getBranchLengths ("mixtureTree", 0);
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
        fprintf (stdout, resStr, "\n");
        three_LF_MLE_VALUES ["restoreLF"][""];
        localBranchRes = three_LF_MLES;
    }
    writeTheLF (".`thisBranchName`.fit", "three_LF");
    fprintf (stdout, "\nModel fit:");   
    BIC_scores [thisBranchName] = BIC (localBranchRes[1][0], localBranchRes[1][1]+baseParameters, sample_size);
    BranchLengthEstimates [thisBranchName] = getBranchLengths ("mixtureTree", 0);
    fprintf     (csvFilePath, "\nGlobal+",thisBranchName,"\t", localBranchRes[1][0],"\t", localBranchRes[1][1]+baseParameters, "\t", BIC_scores [thisBranchName], "\t", PostOrderAVL2StringDistances(mixTreeAVL,BranchLengthEstimates [thisBranchName]));
    fprintf (stdout, "\nLocal branch omega distribution: ");
    ExecuteCommands ("reportOmegaDistro(mixtureTree.`thisBranchName`.omega1,mixtureTree.`thisBranchName`.omega2,mixtureTree.`thisBranchName`.omega3, mixtureTree.`thisBranchName`.Paux1, mixtureTree.`thisBranchName`.Paux2)");
    fprintf (stdout, "Global omega distribution on the rest of the branches: ");
    reportOmegaDistro (omegaG1,omegaG2,omegaG3,Paux1G,Paux2G);
    stashBranchValues (thisBranchName, "branchValues");
    pv = 1-CChi2 (2*(localBranchRes[1][0]-res_three_LF_global[1][0]),5);
    fprintf (stdout, "\nLRT p-value for branch deviation from the global pattern = ", pv, "\n");
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
function constrainABranch (branch_name) {
    ExecuteCommands ("mixtureTree." + branch_name + ".t = mixtureTreeG." + branch_name + ".t");
    ExecuteCommands ("mixtureTree." + branch_name + ".omega1 := omegaG1");
    ExecuteCommands ("mixtureTree." + branch_name + ".omega2 := omegaG2");
    ExecuteCommands ("mixtureTree." + branch_name + ".omega3 := omegaG3");
    ExecuteCommands ("mixtureTree." + branch_name + ".Paux1  := Paux1G");
    ExecuteCommands ("mixtureTree." + branch_name + ".Paux2  := Paux2G");
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------
function stashBranchValues (branch_name, storage&) {

    locals = {{"t","omega1","omega2","omega3","Paux1", "Paux2"}};
    for (_varID = 0; _varID < Columns (locals); _varID += 1) {
        storage ["mixtureTree." + branch_name + "." + locals[_varID]] = Eval ("mixtureTree." + branch_name + "." + locals[_varID]);
    }
    
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------
function unConstrainABranch (branch_name) {
    ExecuteCommands ("mixtureTree." + branch_name + ".t = mixtureTreeG." + branch_name + ".t");
    ExecuteCommands ("mixtureTree." + branch_name + ".omega1 = omegaG1");
    ExecuteCommands ("mixtureTree." + branch_name + ".omega1 :< 1;");
    ExecuteCommands ("mixtureTree." + branch_name + ".omega2 = omegaG2");
    ExecuteCommands ("mixtureTree." + branch_name + ".omega2 :< 1;");
    ExecuteCommands ("mixtureTree." + branch_name + ".omega3 = omegaG3");
    ExecuteCommands ("mixtureTree." + branch_name + ".Paux1  = Paux1G");
    ExecuteCommands ("mixtureTree." + branch_name + ".Paux2  = Paux2G");
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------
function BIC (logL, params, sample_size) {
    BICvalue = -2*logL + Log(sample_size) * params;
    fprintf (stdout, " log(L) = ", logL, " with ", params, " parameters, yielding BIC = ", BICvalue, " assuming the sample size of ", sample_size, ".");
    return BICvalue;
}

//------------------------------------------------------------------------------------------------------------------------
function reportOmegaDistro (w1,w2,w3,p1,p2) {
    fprintf (stdout, "\n\t omega = ", w1, ", p = ", p1,
                     "\n\t omega = ", w2, ", p = ", (1-p1)*p2,
                     "\n\t omega = ", w3, ", p = ", (1-p1)*(1-p2),"\n");
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------

function getBranchLengths (treeID, modelType) {
    blByName        = {};
    bnames_list     = Eval ("BranchName (`treeID`,-1)");
    
    
    if (modelType == 0) {
        Model _temp = (MGMatrix1, codon3x4, 0);
        GetString (model1BL, _temp, -1);
        Model _temp = (MGMatrix2, codon3x4, 0);
        GetString (model2BL, _temp, -1);
        Model _temp = (MGMatrix3, codon3x4, 0);
        GetString (model3BL, _temp, -1);
        bl_expression = "Paux1*(`model1BL`)+(1-Paux1)*Paux2*(`model2BL`)+(1-Paux1)*(1-Paux2)*(`model3BL`)";
        locals = {{"t","omega1","omega2","omega3","Paux1", "Paux2"}};
      
        for (_bID = 0; _bID < Columns (bnames_list) - 1; _bID += 1) {
            branch_name = bnames_list[_bID];
         
             for (_varID = 0; _varID < Columns (locals); _varID += 1) {
                ExecuteCommands (locals[_varID] + " = Eval (\"`treeID`.`branch_name`." + locals[_varID] + "\")");
            }
            blByName [branch_name] = Eval (bl_expression)/3;
        }
    } else {
        Model _temp = (MGMatrix1G, codon3x4, 0);
        GetString (model1BL, _temp, -1);
        Model _temp = (MGMatrix2G, codon3x4, 0);
        GetString (model2BL, _temp, -1);
        Model _temp = (MGMatrix3G, codon3x4, 0);
        GetString (model3BL, _temp, -1);
        bl_expression = "Paux1G*(`model1BL`)+(1-Paux1G)*Paux2G*(`model2BL`)+(1-Paux1G)*(1-Paux2G)*(`model3BL`)";
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

function restoreLF (key, value)
{
	ExecuteCommands (key + " = " + value);
	return 0;
}

