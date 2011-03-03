skipCodeSelectionStep 		= 0;
LoadFunctionLibrary("chooseGeneticCode");

LoadFunctionLibrary("GrabBag");
LoadFunctionLibrary("dSdNTreeTools");
LoadFunctionLibrary("CF3x4");

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function BuildCodonFrequencies (obsF)
{
	PIStop = 1.0;
	result = {ModelMatrixDimension,1};
	hshift = 0;

	for (h=0; h<64; h=h+1)
	{
		first = h$16;
		second = h%16$4;
		third = h%4;
		if (_Genetic_Code[h]==10) 
		{
			hshift = hshift+1;
			PIStop = PIStop-obsF[first][0]*obsF[second][1]*obsF[third][2];
			continue; 
		}
		result[h-hshift][0]=obsF[first][0]*obsF[second][1]*obsF[third][2];
	}
	return result*(1.0/PIStop);
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

global AC = 1;
global AT = 1;
global CG = 1;
global CT = 1;
global GT = 1;

MGCustomRateBiasTerms = {{"AC*","","AT*","CG*","CT*","GT*"}};	

_nucBiasTerms = {4,4};
_nucBiasTerms[0][0] = "";

hv = 0;

for (h=0; h<4; h=h+1)
{
	for (v=h+1; v<4; v=v+1)
	{
		_nucBiasTerms[h][v] = MGCustomRateBiasTerms[hv];
		_nucBiasTerms[v][h] = MGCustomRateBiasTerms[hv];
		hv = hv + 1;	
	}
}

h=0;
v=0;

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function PopulateModelMatrix (ModelMatrixName&, EFV, synrateP, globalP, nonsynRateP)
{
	if (!ModelMatrixDimension)
	{
		ModelMatrixDimension = 64;
		for (h = 0 ;h<64; h=h+1)
		{
			if (_Genetic_Code[h]==10)
			{
				ModelMatrixDimension = ModelMatrixDimension-1;
			}
		}
	}
	
	ModelMatrixName = {ModelMatrixDimension,ModelMatrixDimension}; 

	hshift = 0;
	
	modelDefString = "";
	modelDefString*16384;
	
	catCounterAL = {};
	
	for (h=0; h<64; h=h+1)
	{
		if (_Genetic_Code[h]==10) 
		{
			hshift = hshift+1;
			continue; 
		}
		vshift = hshift;
		for (v = h+1; v<64; v=v+1)
		{
			diff = v-h;
			if (_Genetic_Code[v]==10) 
			{
				vshift = vshift+1;
				continue; 
			}
			nucPosInCodon = 2;
			if ((h$4==v$4)||((diff%4==0)&&(h$16==v$16))||(diff%16==0))
			{
				if (h$4==v$4)
				{
					transition = v%4;
					transition2= h%4;
				}
				else
				{
					if(diff%16==0)
					{
						transition = v$16;
						transition2= h$16;
						nucPosInCodon = 0;
					}
					else
					{
						transition = v%16$4;
						transition2= h%16$4;
						nucPosInCodon = 1;
					}
				}
				hs  = Format(h-hshift,0,0);
				vs  = Format(v-vshift,0,0);
				ts  = Format(transition,0,0);
				ts2 = Format(transition2,0,0);
				ps  = Format(nucPosInCodon,0,0);
				aa1 = _Genetic_Code[0][h];
				aa2 = _Genetic_Code[0][v];
				if (aa1==aa2) 
				{
					modelDefString*("ModelMatrixName["+hs+"]["+vs+"] := "+_nucBiasTerms[transition][transition2]+""+synrateP+"*EFV__["+ts+"]["+ps+"];\n"+
									"ModelMatrixName["+vs+"]["+hs+"] := "+_nucBiasTerms[transition][transition2]+""+synrateP+"*EFV__["+ts2+"]["+ps+"];\n");
				}
				else
				{
					if (Abs(globalP))
					{
						modelDefString*("ModelMatrixName["+hs+"]["+vs+"] := "+_nucBiasTerms[transition][transition2]+""+synrateP+"*"+globalP+"*EFV__["+ts+"]["+ps+"];\n"+
										"ModelMatrixName["+vs+"]["+hs+"] := "+_nucBiasTerms[transition][transition2]+""+synrateP+"*"+globalP+"*EFV__["+ts2+"]["+ps+"];\n");	
					
					}
					else
					{
						modelDefString*("ModelMatrixName["+hs+"]["+vs+"] := "+_nucBiasTerms[transition][transition2]+""+nonsynRateP+"*EFV__["+ts+"]["+ps+"];\n"+
										"ModelMatrixName["+vs+"]["+hs+"] := "+_nucBiasTerms[transition][transition2]+""+nonsynRateP+"*EFV__["+ts2+"]["+ps+"];\n");	
					}
				}
			}
	    }
    }		
	modelDefString*0;
	ExecuteCommands (modelDefString);
	return 0;
}


DataSet 			ds 				= ReadDataFile(PROMPT_FOR_FILE);
DataSetFilter 		dsf 			= CreateFilter(ds,3,"","",GeneticCodeExclusions);

HarvestFrequencies	(nuc3, dsf, 3, 1, 1);

omega1 = 0.2;
omega2 = 0.5;
omega3 = 1.0;

nucCF						= CF3x4	(nuc3, GeneticCodeExclusions);

PopulateModelMatrix			  ("MGMatrix1",  nucCF, "t1", "omega1", "");
PopulateModelMatrix			  ("MGMatrix2",  nucCF, "t2", "omega2", "");
PopulateModelMatrix			  ("MGMatrix3",  nucCF, "t3", "omega3", "");

global	omegaG1 = 0.2;
omegaG1 :< 1;
global	omegaG2 = 0.5;
omegaG2 :< 1;
global	omegaG3 = 2.0;
omegaG3 :> 1;

PopulateModelMatrix			  ("MGMatrix1G",  nucCF, "t1", "omegaG1", "");
PopulateModelMatrix			  ("MGMatrix2G",  nucCF, "t2", "omegaG2", "");
PopulateModelMatrix			  ("MGMatrix3G",  nucCF, "t3", "omegaG3", "");



PopulateModelMatrix			  ("MGMatrixLocal",  nucCF, "syn", "", "nonsyn");

codon3x4					= BuildCodonFrequencies (nucCF);
Model		MGL				= (MGMatrixLocal, codon3x4, 0);

LoadFunctionLibrary			  ("queryTree");

SetDialogPrompt ("Save analysis results to");
fprintf (PROMPT_FOR_FILE, CLEAR_FILE, KEEP_OPEN,"Branch,Mean_dNdS,Omega1,P1,Omega2,P2,Omega3,P3,LRT,p,p_Holm");
csvFilePath = LAST_FILE_PATH;

fprintf 					  (stdout, "[PHASE 0] Fitting the local MG94 (no site-to-site variation) to obtain initial parameter estimates\n");
VERBOSITY_LEVEL				 = 0;

LikelihoodFunction	base_LF	 = (dsf, givenTree);
Optimize					  (res_base,base_LF);

lfOut	= csvFilePath + ".mglocal.fit";
LIKELIHOOD_FUNCTION_OUTPUT = 7;
fprintf (lfOut, CLEAR_FILE, base_LF);
LIKELIHOOD_FUNCTION_OUTPUT = 2;

localLL						 = res_base[1][0];
localParams					 = res_base[1][1] + 9;

LoadFunctionLibrary			 ("DescriptiveStatistics");

totalBranchCount			 = BranchCount(givenTree) + TipCount (givenTree);

//GetInformation	   		   (varNames, "givenTree\\..*\\.omega1");
//localOmegaValues			 = {totalBranchCount,1}["Eval(varNames[_MATRIX_ELEMENT_ROW_])"];

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



ReplicateConstraint 		  ("this1.?.t1:=this2.?.syn",mixtureTreeG,givenTree);
ReplicateConstraint 		  ("this1.?.t1:=this2.?.syn",mixtureTree,givenTree);

ClearConstraints			  (mixtureTree);
ClearConstraints			  (mixtureTreeG);

omega1G						 :< 1;
omega2G						 :< 1;
Paux1G 						 :< 1;
Paux2G 						 :< 1;

ReplicateConstraint 		  ("this1.?.t2:=this2.?.t1",mixtureTree,mixtureTree);
ReplicateConstraint 		  ("this1.?.t3:=this2.?.t1",mixtureTree,mixtureTree);
ReplicateConstraint 		  ("this1.?.t2:=this2.?.t1",mixtureTreeG,mixtureTreeG);
ReplicateConstraint 		  ("this1.?.t3:=this2.?.t1",mixtureTreeG,mixtureTreeG);

ASSUME_REVERSIBLE_MODELS	  = 1;

/*LikelihoodFunction three_LF   = (dsf,mixtureTreeG);


fprintf 					  (stdout, "[PHASE 1] Fitting a GLOBAL branch-site matrix mixture\n");

Optimize					  (res_three_LF_global,three_LF);
fprintf						  (stdout,"\n",three_LF);

lfOut	= csvFilePath + ".relglobal.fit";
LIKELIHOOD_FUNCTION_OUTPUT = 7;
fprintf (lfOut, CLEAR_FILE, three_LF);
LIKELIHOOD_FUNCTION_OUTPUT = 2;

global Paux1G=0.6417071370308534;
global Paux2G=1;
global omegaG3=1.064891502388613;
global omegaG2=0.9203043460462096;
global omegaG1=0.3863674819501042;


wg1 = Max(0.05,Paux1G);
wg2 = Max(0.05,(1-Paux1G)*Paux2G);
wg3 = Max(0.05,(1-Paux1G)*(1-Paux2G));

sum = wg1+wg2+wg3;

wg1 = wg1/sum;
wg2 = wg2/sum;
wg3 = wg3/sum;

Paux1G = wg1;
Paux2G = wg2/(1-Paux1G);*/

USE_LAST_RESULTS			  = 1;

LikelihoodFunction three_LF   = (dsf,mixtureTree);


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
    
    ExecuteCommands ("mixtureTree." + bNames[k] + ".t1 = bl");
    ExecuteCommands ("mixtureTree." + bNames[k] + ".omega1 :< 1;");
	ExecuteCommands ("mixtureTree." + bNames[k] + ".omega2 :< 1;");
    if (baseOmega > 1)
    {
        ExecuteCommands ("mixtureTree." + bNames[k] + ".omega1 = 0.1;");
        ExecuteCommands ("mixtureTree." + bNames[k] + ".omega2 = 1;");
        ExecuteCommands ("mixtureTree." + bNames[k] + ".omega3 = baseOmega;");

        ExecuteCommands ("mixtureTree." + bNames[k] + ".Paux1 = 0.01;");
        ExecuteCommands ("mixtureTree." + bNames[k] + ".Paux2 = 0.01;");
    }
    else
    {
        ExecuteCommands ("mixtureTree." + bNames[k] + ".omega1 = baseOmega;");
        ExecuteCommands ("mixtureTree." + bNames[k] + ".omega2 = 1;");
        ExecuteCommands ("mixtureTree." + bNames[k] + ".omega3 = 2;");

        ExecuteCommands ("mixtureTree." + bNames[k] + ".Paux1 = 0.98;");
        ExecuteCommands ("mixtureTree." + bNames[k] + ".Paux2 = 0.5;");    
    }
}


//VERBOSITY_LEVEL = 10;
OPTIMIZATION_METHOD = 0;

fprintf 					  (stdout, "[PHASE 2] Fitting the full LOCAL alternative model (no constraints)\n");
Optimize					  (res_three_LF,three_LF);
fprintf						  (stdout,"\n",three_LF);


lfOut	= csvFilePath + ".fit";
LIKELIHOOD_FUNCTION_OUTPUT = 7;
fprintf (lfOut, CLEAR_FILE, three_LF);
LIKELIHOOD_FUNCTION_OUTPUT = 2;



for							  (k = 0; k < totalBranchCount; k = k+1)
{
	ref = "mixtureTree."+bNames[k];
	
	thisOmega3 = Eval (ref+".omega3");
	wt3        = Eval ("(1-"+ref+".Paux1)*(1-"+ref+".Paux2)");

	pValueByBranch [k][1] = Eval (ref+".omega1");
	pValueByBranch [k][2] = Eval (ref+".Paux1");
	pValueByBranch [k][3] = Eval (ref+".omega2");
	pValueByBranch [k][4] = Eval ("(1-"+ref+".Paux1)*"+ref+".Paux2");
	pValueByBranch [k][5] = thisOmega3;
	pValueByBranch [k][6] = wt3;
	
	fprintf (stdout, "\nNode: ", ref, 
					 "\n\tClass 1: omega = ", Eval (ref+".omega1"), " weight = ", Eval (ref+".Paux1"),
					 "\n\tClass 2: omega = ", Eval (ref+".omega2"), " weight = ", Eval ("(1-"+ref+".Paux1)*"+ref+".Paux2"),
					 "\n\tClass 3: omega = ", thisOmega3, " weight = ", wt3 , "\n"
					 );
		
	if (thisOmega3 > 1 && wt3 > 1e-6)
	{
		fprintf (stdout, "...Testing for selection at this branch\n");
		_stashLF = saveLF ("three_LF");

		ExecuteCommands ("mixtureTree." + bNames[k] + ".omega3 := 1");
		Optimize					  (res_three_current,three_LF);
		
		fprintf (stdout, "\nNode: ", ref, 
						 "\n\tClass 1: omega = ", Eval (ref+".omega1"), " weight = ", Eval (ref+".Paux1"),
						 "\n\tClass 2: omega = ", Eval (ref+".omega2"), " weight = ", Eval ("(1-"+ref+".Paux1)*"+ref+".Paux2"),
						 "\n\tClass 3: omega = ", Eval (ref+".omega3"), " weight = ", Eval ("(1-"+ref+".Paux1)*(1-"+ref+".Paux2)"), "\n"
						 );
		pValueByBranch[k][7]			  = 2*(res_three_LF[1][0] - res_three_current[1][0]);				 
		pValueByBranch[k][8]			  = (1-CChi2 (pValueByBranch[k][7],1))*.5;
		fprintf (stdout, "\np-value = ", pValueByBranch[k][8],"\n\n", three_LF, "\n");
		
		ExecuteCommands ("mixtureTree." + bNames[k] + ".omega3 :< 1e26");
		
		if (pValueByBranch[k][7] < (-0.5))
		{
			fprintf 					  (stdout, "[PHASE 2/REPEAT] Detected a convergence problem; refitting the LOCAL alternative model with new starting values\n");
			lfOut	= csvFilePath + ".fit";
			Optimize					  (res_three_LF,three_LF);
			LIKELIHOOD_FUNCTION_OUTPUT = 7;
			fprintf (lfOut, CLEAR_FILE, three_LF);
			LIKELIHOOD_FUNCTION_OUTPUT = 2;
			_stashLF = saveLF ("three_LF");
			k = 0;
		}
		else
		{
			_stashLF ["restoreLF"][""];
		}
	}
	else
	{
		pValueByBranch[k][8] = 1.0;
	}
}

OPTIMIZATION_METHOD = 4;


pValueSorter = {totalBranchCount,2};
pValueSorter = pValueSorter["_MATRIX_ELEMENT_ROW_*(_MATRIX_ELEMENT_COLUMN_==0)+pValueByBranch[_MATRIX_ELEMENT_ROW_][8]*(_MATRIX_ELEMENT_COLUMN_==1)"];
pValueSorter = pValueSorter % 1;
pValueSorter = pValueSorter["_MATRIX_ELEMENT_VALUE_*(_MATRIX_ELEMENT_COLUMN_==0)+_MATRIX_ELEMENT_VALUE_*(totalBranchCount-_MATRIX_ELEMENT_ROW_)*(_MATRIX_ELEMENT_COLUMN_==1)"];

fprintf (stdout,"\n\nSummary of branches under episodic selection:\n");
hasBranchesUnderSelection = 0;

pthreshold = 0.05;

for		(k = 0; k < totalBranchCount; k = k+1)
{
	pValueByBranch[pValueSorter[k][0]][9] = Min (1,pValueSorter[k][1]);
	if (pValueSorter[k][1] <= pthreshold)
	{
		fprintf (stdout, "\t", bNames[pValueSorter[k][0]], " p = ", pValueByBranch[pValueSorter[k][0]][9], "\n");
		hasBranchesUnderSelection += 1;
	}
}


if (hasBranchesUnderSelection == 0)
{
	fprintf (stdout, "\tNo branches found to be under selection at p <= ", threshold, "\n");
}

for		(k = 0; k < totalBranchCount; k = k+1)
{
	fprintf (csvFilePath, "\n", bNames[k], ",", Join(",",pValueByBranch[k][-1]));
}
fprintf (csvFilePath, CLOSE_FILE);

//---- TREE RENDERING -----

LoadFunctionLibrary ("ReadDelimitedFiles");
treeString = Format (givenTree, 1,1);
resultTable=ReadCSVTable (csvFilePath, 2);

rows	= Rows (resultTable[2]);
columns = Columns(resultTable[1]);

UseModel (USE_NO_MODEL);
Tree T = treeString;

TREE_OUTPUT_OPTIONS = {"__FONT_SIZE__":14};

tavl = T^0;
for (k = 1; k < Abs (tavl)-1; k+=1)
{
	TREE_OUTPUT_OPTIONS[(tavl[k])["Name"]] = {};
}

for (k = 1; k < Abs (tavl)-1; k+=1)
{
	thisP = (resultTable[1])[k-1][9];
	
	parentName = (tavl[((tavl[k])["Parent"])])["Name"];
	
	myRateDistro = {3,2};
	myRateDistro [0][0] = (resultTable[1])[k-1][1];
	myRateDistro [0][1] = (resultTable[1])[k-1][2];
	myRateDistro [1][0] = (resultTable[1])[k-1][3];
	myRateDistro [1][1] = (resultTable[1])[k-1][4];
	myRateDistro [2][0] = (resultTable[1])[k-1][5];
	myRateDistro [2][1] = (resultTable[1])[k-1][6];
	
	myRateDistro = myRateDistro % 0;
	
	color1 = makeDNDSColor (myRateDistro[0][0]);
	color2 = makeDNDSColor (myRateDistro[1][0]);
	color3 = makeDNDSColor (myRateDistro[2][0]);
	
	(TREE_OUTPUT_OPTIONS[(tavl[k])["Name"]])["TREE_OUTPUT_BRANCH_COLOR"] 		= {{color1__[0],color1__[1],color1__[2],myRateDistro__[0][1]}
																				   {color2__[0],color2__[1],color2__[2],myRateDistro__[1][1]}
																				   {color3__[0],color3__[1],color3__[2],myRateDistro__[2][1]}};
	(TREE_OUTPUT_OPTIONS[(tavl[k])["Name"]])["TREE_OUTPUT_BRANCH_LINECAP"]  = 	0;
	
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
function makeDNDSColor (omega)
{
	if (omega < 1)
	{
		return {{0.5*omega__,0.5*omega__,1-0.5*omega__}};
	}
	omega = Min (omega,5);
	return {{0.5+0.125*(omega__-1),0.5-(omega__-1)*0.125,0.5-(omega__-1)*0.125}};
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

