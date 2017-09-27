ModelMatrixDimension = 0;
_DO_TREE_REBALANCE_ = 0;

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

function PopulateModelMatrix (ModelMatrixName&, EFV)
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
	
	if (modelType > 0)
	{
		catCounterAL = {};
	}
	
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
				hs = Format(h-hshift,0,0);
				vs = Format(v-vshift,0,0);
				ts = Format(transition,0,0);
				ts2= Format(transition2,0,0);
				ps = Format(nucPosInCodon,0,0);
				aa1 = _Genetic_Code[0][h];
				aa2 = _Genetic_Code[0][v];
				if (aa1==aa2) 
				{
					modelDefString*("ModelMatrixName["+hs+"]["+vs+"] := "+_nucBiasTerms[transition][transition2]+"synRate*EFV__["+ts+"]["+ps+"];\n"+
													 "ModelMatrixName["+vs+"]["+hs+"] := "+_nucBiasTerms[transition][transition2]+"synRate*EFV__["+ts2+"]["+ps+"];\n");
				}
				else
				{
					bt = "nsClass"+aaRateMultipliers[aa1][aa2];
					modelDefString*("ModelMatrixName["+hs+"]["+vs+"] := "+_nucBiasTerms2[transition][transition2]+bt+"*synRate*EFV__["+ts+"]["+ps+"];\n"+
														 "ModelMatrixName["+vs+"]["+hs+"] := "+_nucBiasTerms2[transition][transition2]+bt+"*synRate*EFV__["+ts2+"]["+ps+"];\n");						
				}
			}
	    }
    }		
	modelDefString*0;
	ExecuteCommands (modelDefString);
	return 0;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

incFileName = HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"chooseGeneticCode.def";
ExecuteCommands  ("#include \""+incFileName+"\";");

SetDialogPrompt ("Select a codon data file");
DataSet ds 	= ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);


/*---------------------------------------------------------------------------------------------------------------------------------------------*/

ChoiceList	(rateType,"Site-to-site rate variation model",1,SKIP_NONE,
					   "None", "No site-by-site rate variation.",
					   "dN only",  "Non-synonymous rate varies across sites via a gamma/beta distribution, synonymous rate is fixed at 1.",
					   "Complete",  "Both synonymous and non-synonymous rate varies across sites via two independent gamma/beta distributions.");

ClearConstraints (NC);
ClearConstraints (SC);

if (rateType < 0)
{
	return;
}

if (rateType == 0)
{
	MGCustomRateBiasTermsS = {{"AC*","","AT*","CG*","CT*","GT*"}};	
	MGCustomRateBiasTermsN = MGCustomRateBiasTermsS;	
}
else
{
	done = 0;
	while (!done)
	{
		fprintf (stdout,"\nHow many rate classes should there be ([3-8] is a good range, must be >=2):");
		fscanf  (stdin,"Number", resp);
		if (resp>=2)
		{	
			done = 1;
		}
	}			

	global betaNSP = 1;
	global betaNSQ = 1;
	betaNSP:>0.05;betaNSP:<85;
	betaNSQ:>0.05;betaNSQ:<85;
	category _pd = (resp-1, EQUAL, MEAN, 
					_x_^(betaNSP-1)*(1-_x_)^(betaNSQ-1)/Beta(betaNSP,betaNSQ), /* density */
					IBeta(_x_,betaNSP,betaNSQ), /*CDF*/
					0, 				   /*left bound*/
					1, 			   /*right bound*/
					IBeta(_x_,betaNSP+1,betaNSQ)*betaNSP/(betaNSP+betaNSQ)
	);

	global beta = .5;
	beta:>0.01;beta:<100;
	category NC = (resp, _pd, MEAN, 
					GammaDist(_x_,beta,beta), 
					CGammaDist(_x_,beta,beta), 
					0 , 
			  	    1e25,
			  	    CGammaDist(_x_,beta+1,beta)
			  	 );
		  	 
	if (rateType == 1)
	{
		MGCustomRateBiasTermsS = {{"AC*","","AT*","CG*","CT*","GT*"}};	
		MGCustomRateBiasTermsN = {{"AC*NC*","NC*","AT*NC*","CG*NC*","CT*NC*","GT*NC*"}};
	}
	else
	{
		global alpha = 1;
		alpha:>0.01;
		alpha:<100;

		global betaP = 1;
		global betaQ = 1;
		betaP:>0.05;betaP:<85;
		betaQ:>0.05;betaQ:<85;
		
		category _pc = (resp-1, EQUAL, MEAN, 
						_x_^(betaP-1)*(1-_x_)^(betaQ-1)/Beta(betaP,betaQ), /* density */
						IBeta(_x_,betaP,betaQ), /*CDF*/
						0, 				   /*left bound*/
						1, 			   /*right bound*/
						IBeta(_x_,betaP+1,betaQ)*betaP/(betaP+betaQ)
		);

		global alpha = .5;
		alpha:>0.01;alpha:<100;
		category SC = (resp, _pc, MEAN, 
						GammaDist(_x_,alpha,alpha), 
						CGammaDist(_x_,alpha,alpha), 
						0 , 
				  	    1e25,
				  	    CGammaDist(_x_,alpha+1,alpha)
				  	 );
		  	 
		MGCustomRateBiasTermsS = {{"AC*SC*","SC*","AT*SC*","CG*SC*","CT*SC*","GT*SC*"}};
		MGCustomRateBiasTermsN = {{"AC*NC*","NC*","AT*NC*","CG*NC*","CT*NC*","GT*NC*"}};	
	}
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

global AC = 1;
global AT = 1;
global CG = 1;
global CT = 1;
global GT = 1;
		
				
ChoiceList	(modelSpec,"Nucleotide Substitution Model",1,SKIP_NONE,
					   "Default", "MG94xHKY85. Only corrects for transition/transversion bias.",
					   "Custom",  "MG94xarbitrary nucleotide bias corrections, from none (xF81) to full (xREV).");
					   
if (modelSpec < 0)
{
	return;
}

if (modelSpec == 1)
{
	done = 0;
	while (!done)
	{
		fprintf (stdout,"\nPlease enter a 6 character model designation (e.g:010010 defines HKY85):");
		fscanf  (stdin,"String", modelDesc);
		if (Abs(modelDesc)==6)
		{	
			done = 1;
		}
	}			
}
else
{
	modelDesc = "010010";
}

paramCount	  = 0;
_nucBiasTerms = {4,4};
_nucBiasTerms[0][0] = "";

_nucBiasTerms2 = {4,4};
_nucBiasTerms2[0][0] = "";


if (modelDesc[0]==modelDesc[1])
{
	MGCustomRateBiasTermsS[0] = MGCustomRateBiasTermsS[1];
	MGCustomRateBiasTermsN[0] = MGCustomRateBiasTermsN[1];
}

_nucBiasTerms[1][0] = MGCustomRateBiasTermsS[0];
_nucBiasTerms[0][1] = MGCustomRateBiasTermsS[0];
_nucBiasTerms[2][0] = MGCustomRateBiasTermsS[1];
_nucBiasTerms[0][2] = MGCustomRateBiasTermsS[1];

_nucBiasTerms2[1][0] = MGCustomRateBiasTermsN[0];
_nucBiasTerms2[0][1] = MGCustomRateBiasTermsN[0];
_nucBiasTerms2[2][0] = MGCustomRateBiasTermsN[1];
_nucBiasTerms2[0][2] = MGCustomRateBiasTermsN[1];


h = 0;
v = 3;

for (customLoopCounter2=2; customLoopCounter2<6; customLoopCounter2=customLoopCounter2+1)
{
	for (customLoopCounter=0; customLoopCounter<customLoopCounter2; customLoopCounter=customLoopCounter+1)
	{
		if (modelDesc[customLoopCounter]==modelDesc[customLoopCounter2])
		{
			_nucBiasTerms[h][v] = MGCustomRateBiasTermsS[customLoopCounter];
			_nucBiasTerms[v][h] = MGCustomRateBiasTermsS[customLoopCounter];
			_nucBiasTerms2[h][v] = MGCustomRateBiasTermsN[customLoopCounter];
			_nucBiasTerms2[v][h] = MGCustomRateBiasTermsN[customLoopCounter];
			break;
		}
	}
	if (customLoopCounter == customLoopCounter2)
	{
		_nucBiasTerms[h][v]  = MGCustomRateBiasTermsS[customLoopCounter2];
		_nucBiasTerms[v][h]  = MGCustomRateBiasTermsS[customLoopCounter2];
		_nucBiasTerms2[h][v] = MGCustomRateBiasTermsN[customLoopCounter2];
		_nucBiasTerms2[v][h] = MGCustomRateBiasTermsN[customLoopCounter2];
	}
	
	v = v+1;
	if (v==4)
	{
		h=h+1;
		v=h+1;
	}
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

incFileName = HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"MGwAA.ibf";
ExecuteCommands  ("#include \""+incFileName+"\";");

HarvestFrequencies (observedFreq,filteredData,3,1,1);
MG94custom = 0;
MULTIPLY_BY_FREQS = PopulateModelMatrix ("MG94custom", observedFreq);
vectorOfFrequencies = BuildCodonFrequencies (observedFreq);
Model MG94customModel = (MG94custom,vectorOfFrequencies,0);
USE_POSITION_SPECIFIC_FREQS = 1;

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

incFileName = HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"queryTree.bf";
ExecuteCommands  ("#include \""+incFileName+"\";");

Tree		   givenTree = treeString;

internalNodes = BranchCount(givenTree);
leafNodes     = TipCount(givenTree);
choiceMatrix  = {internalNodes+leafNodes,2};

for (bc=0; bc<internalNodes; bc=bc+1)
{
	choiceMatrix[bc][0] = BranchName(givenTree,bc);
	choiceMatrix[bc][1] = "Internal Branch Rooting " + givenTree[bc];
}

for (bc=0; bc<leafNodes; bc=bc+1)
{
	choiceMatrix[bc+internalNodes][0] = TipName(givenTree,bc);
	choiceMatrix[bc+internalNodes][1] = "Terminal branch endin in " + choiceMatrix[bc+internalNodes][0];
}

mxTreeSpec = {5,1};
mxTreeSpec [0] = "givenTree";
mxTreeSpec [1] = "8240";
mxTreeSpec [2] = "10,40,-10,-175,1";
mxTreeSpec [3] = "";
mxTreeSpec [4] = "";
OpenWindow (TREEWINDOW, mxTreeSpec);

ChoiceList  (bOption,"Choose the branch to test",0,NO_SKIP,choiceMatrix);

if (bOption[0] < 0)
{
	return -1;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

nsRateClasses = Rows(aaRateClassIDs);
nsChoiceMatrix = {Columns(nsRateClasses),2};
defineGlobals  = "";
defineGlobals * 1024;

for (bc=0; bc<Columns(nsRateClasses); bc=bc+1)
{
	nsChoiceMatrix[bc][0] = nsRateClasses[bc];
	nsChoiceMatrix[bc][1] = "Non-synonymous rate parameter nsClass" + nsRateClasses[bc];
	defineGlobals * ("global SharedNS"+nsRateClasses[bc]+"=1;ReplicateConstraint(\"this1.?.nsClass"+nsRateClasses[bc]+":=SharedNS"+nsRateClasses[bc]+"\",givenTree);\n");
}

defineGlobals * 0;

ExecuteCommands (defineGlobals);
							

ChoiceList	(paramSpec,"Parameters to test:",0,SKIP_NONE,nsChoiceMatrix);

if (paramSpec[0]<0)
{
	return 0;
}

LikelihoodFunction lf = (filteredData, givenTree);

fprintf (stdout, "\n1. Working on the null hypothesis.\n");
Optimize (res_null, lf);
fprintf (stdout, lf);

fprintf (stdout, "\n\n2. Working on the alternative hypothesis.\n");

for (h=0; h<Columns(bOption)*Rows(bOption); h=h+1)
{
	v = bOption[h];
	for (bc=0; bc<Columns(paramSpec)*Rows(paramSpec); bc=bc+1)
	{
		k = paramSpec[bc];
		ExecuteCommands ("givenTree."+choiceMatrix[v][0]+".nsClass"+nsRateClasses[k]+"="+"givenTree."+choiceMatrix[v][0]+".nsClass"+nsRateClasses[k]+";");
	}
}

som = OPTIMIZATION_METHOD;
sulr = USE_LAST_RESULTS;
OPTIMIZATION_METHOD = 0;
USE_LAST_RESULTS = 1;
Optimize (res_alt, lf);
fprintf (stdout, lf);
OPTIMIZATION_METHOD = som;
USE_LAST_RESULTS = sulr;

fprintf (stdout, "\nBranch parameter estimates\n\n");

for (h=0; h<Columns(bOption)*Rows(bOption); h=h+1)
{
	v = bOption[h];
	for (bc=0; bc<Columns(paramSpec)*Rows(paramSpec); bc=bc+1)
	{
		k = paramSpec[bc];
		paramName = "givenTree."+choiceMatrix[v][0]+".nsClass"+nsRateClasses[k];
		COVARIANCE_PRECISION = 0.95;
		COVARIANCE_PARAMETER = paramName;
		CovarianceMatrix (cmx, lf);
		fprintf (stdout, paramName, " = ", cmx[1], ", 95% CI: ", cmx[0], "-", cmx[2], "\n");
	}
}

lnLikDiff = -2(res_null[1][0]-res_alt[1][0]);
degFDiff = res_alt[1][1]-res_null[1][1];

fprintf (stdout, "\n\n-2(Ln Likelihood Difference)=",lnLikDiff,"\n","Difference in number of parameters:",Format(degFDiff,0,0));
fprintf (stdout, "\np-value:",1-CChi2(lnLikDiff,degFDiff));





