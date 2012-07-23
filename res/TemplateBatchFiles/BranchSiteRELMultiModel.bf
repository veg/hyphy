skipCodeSelectionStep 		= 0;
LoadFunctionLibrary("chooseGeneticCode");

LoadFunctionLibrary("GrabBag");
LoadFunctionLibrary("dSdNTreeTools");
LoadFunctionLibrary("CF3x4");
LoadFunctionLibrary("BranchSiteTemplate");



DataSet 			ds 				= ReadDataFile(PROMPT_FOR_FILE);
DataSetFilter 		dsf 			= CreateFilter(ds,3,"","",GeneticCodeExclusions);

HarvestFrequencies	(nuc3, dsf, 3, 1, 1);

omega1 = 0.2;
omega2 = 0.5;
omega3 = 1.0;

nucCF						= CF3x4	(nuc3, GeneticCodeExclusions);

PopulateModelMatrix			  ("MGMatrix1",  nucCF, "t", "omega1", "");
PopulateModelMatrix			  ("MGMatrix2",  nucCF, "t", "omega2", "");
PopulateModelMatrix			  ("MGMatrix3",  nucCF, "t", "omega3", "");

global	omegaG1 = 0.2;
omegaG1 :< 1;
global	omegaG2 = 0.5;
omegaG2 :< 1;
global	omegaG3 = 2.0;
omegaG3 :> 1;

PopulateModelMatrix			  ("MGMatrix1G",  nucCF, "t", "omegaG1", "");
PopulateModelMatrix			  ("MGMatrix2G",  nucCF, "t", "omegaG2", "");
PopulateModelMatrix			  ("MGMatrix3G",  nucCF, "t", "omegaG3", "");



PopulateModelMatrix			  ("MGMatrixLocal",  nucCF, "syn", "", "nonsyn");

codon3x4					= BuildCodonFrequencies (nucCF);
Model		MGL				= (MGMatrixLocal, codon3x4, 0);

LoadFunctionLibrary			  ("queryTree");

SetDialogPrompt ("Save analysis results to");
fprintf (PROMPT_FOR_FILE, CLEAR_FILE, KEEP_OPEN,"Branch,Mean_dNdS,Omega1,P1,Omega2,P2,Omega3,P3,LRT,p,p_Holm");
csvFilePath = LAST_FILE_PATH;

fprintf 					  (stdout, "[BS-REL PHASE 0] Fitting the local MG94 (no site-to-site variation) to obtain initial parameter estimates\n");
VERBOSITY_LEVEL				 = 0;

LikelihoodFunction	base_LF	 = (dsf, givenTree);
Optimize					  (res_base,base_LF);

AC := AC__;
AT := AT__;
CG := CG__;
CT := CT__;
GT := GT__;

lfOut	= csvFilePath + ".mglocal.fit";
LIKELIHOOD_FUNCTION_OUTPUT = 7;
fprintf (lfOut, CLEAR_FILE, base_LF);
LIKELIHOOD_FUNCTION_OUTPUT = 2;

baseParameters                = 9;
localLL						 = res_base[1][0];
localParams					 = res_base[1][1] + baseParameters;

LoadFunctionLibrary			 ("DescriptiveStatistics");

totalBranchCount			 = BranchCount(givenTree) + TipCount (givenTree);

pValueByBranch				  = {totalBranchCount,10};
bNames						  = BranchName (givenTree, -1);

for (k = 0; k < totalBranchCount; k = k+1)
{
	srate  = Eval ("givenTree." + bNames[k] + ".syn");
	nsrate = Eval ("givenTree." + bNames[k] + ".nonsyn");
	if (srate > 0)
	{
		pValueByBranch [k][0] = Min (10, nsrate/srate);
	}
	else
	{
		pValueByBranch [k][0] = 10;
	}	
}

omegaStats					 = GatherDescriptiveStats (pValueByBranch[-1][0]);

fprintf						 (stdout, "\nLog L = ", localLL, " with ", localParams, " degrees of freedom\n");

PrintDescriptiveStats		 ("Branch omega values", omegaStats);

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


LikelihoodFunction three_LFG   = (dsf,mixtureTreeG);
fprintf 					  (stdout, "[BS-REL PHASE 1] Fitting a GLOBAL branch-site matrix mixture\n");

for (k = 0; k < totalBranchCount; k = k+1)
{
    if (k == 0)
    {
        expr            = Eval("BranchLength(givenTree,\""+bNames[0]+";EXPECTED_NUMBER_OF_SUBSTITUTIONS\")");
        syn             = 1; nonsyn = 0;
        synM            = Eval(expr);
        syn             = 0; nonsyn = 1;
        nonsynM         = Eval(expr);
    }
    
 	srate  = Eval ("givenTree." + bNames[k] + ".syn");
	nsrate = Eval ("givenTree." + bNames[k] + ".nonsyn");
    bl = Eval("BranchLength(givenTree,\""+bNames[k]+"\")")*3;
    
    if (srate > 0)
    {
        baseOmega = nsrate/srate;
    }
    else
    {
        baseOmega = 10000;
    }
        
    bl = bl / (synM + nonsynM * baseOmega);  
    ExecuteCommands ("mixtureTreeG." + bNames[k] + ".t = bl");
}

USE_LAST_RESULTS			  = 1;
VERBOSITY_LEVEL               = 1;

Optimize					  (res_three_LF_global,three_LFG);

baseParameters                += 5;
sample_size                    = dsf.sites * dsf.species;



lfOut	= csvFilePath + ".relglobal.fit";

LIKELIHOOD_FUNCTION_OUTPUT = 7;
fprintf (lfOut, CLEAR_FILE, three_LF);
LIKELIHOOD_FUNCTION_OUTPUT = 2;

fprintf (stdout, "Global model fit:");

BIC_scores ["global"] = BIC (res_three_LF_global[1][0], res_three_LF_global[1][1]+baseParameters, sample_size);
fprintf (stdout, "\nInferred this global omega distribution: ");
reportOmegaDistro (omegaG1,omegaG2,omegaG3,Paux1G,Paux2G);
globalState = saveLF ("three_LFG");

for (k = 0; k < totalBranchCount; k = k+1)
{
    constrainABranch (bnames[k-1]);
}

fprintf (stdout, "\n[BS-REL PHASE 2] Fitting a GLOBAL branch-site matrix mixture with a SINGLE unconstrained branch\n");


LikelihoodFunction three_LF   = (dsf,mixtureTree);

for (k = 0; k < totalBranchCount; k = k+1) {
    fprintf (stdout, "[BS-REL PHASE 2. Branch '", bnames[k], "']\n");
    globalState["restoreLF"][""];
    if (k) {
        constrainABranch (bnames[k-1]);
    }
    unConstrainABranch (bnames[k]);
    Optimize (localBranchRes, three_LF);
}

//------------------------------------------------------------------------------------------------------------------------
function constrainABranch (branch_name) {
    ExecuteCommands ("mixtureTree." + branch_name + ".t = mixtureTreeG." + bNames[k] + ".t");
    ExecuteCommands ("mixtureTree." + branch_name + ".omega1 := omegaG1");
    ExecuteCommands ("mixtureTree." + branch_name + ".omega2 := omegaG2");
    ExecuteCommands ("mixtureTree." + branch_name + ".omega3 := omegaG3");
    ExecuteCommands ("mixtureTree." + branch_name + ".Paux1  := Paux1G");
    ExecuteCommands ("mixtureTree." + branch_name + ".Paux2  := Paux2G");
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------
function unConstrainABranch (branch_name) {
    ExecuteCommands ("mixtureTree." + branch_name + ".t = mixtureTreeG." + bNames[k] + ".t");
    ExecuteCommands ("mixtureTree." + branch_name + ".omega1 = omegaG1");
    ExecuteCommands ("mixtureTree." + branch_name + ".omega2 = omegaG2");
    ExecuteCommands ("mixtureTree." + branch_name + ".omega3 = omegaG3");
    ExecuteCommands ("mixtureTree." + branch_name + ".Paux1  = Paux1G");
    ExecuteCommands ("mixtureTree." + branch_name + ".Paux2  = Paux2G");
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------
function BIC (logL, params, sample_size) {
    BICvalue = -2*logL + Log(sample_size) * params;
    fprintf (stdout, " log(L) = ", logL, " with ", params, " parameters, yielding BIC = ", BICvalue);
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

