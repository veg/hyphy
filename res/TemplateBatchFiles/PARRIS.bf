ModelMatrixDimension = 0;

baselineOutput   = "";

/*---------------------------------------------------------------------------------------------------------------------------------------------------*/

function computeExpSubWeights (dum)
{
	rateTypeWeights     = {6,2};
	rateTypeWeightsGY94 = {6,2};
		
	hshift = 0;

	for (h=0; h<64; h=h+1)
	{
		if (_Genetic_Code[h]==10) 
		{
			continue; 
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
					
					mxIndex = 0;
					if (transition<transition2)
					{
						t1 = transition;
						t2 = transition2;
					}
					else
					{
						t1 = transition2;
						t2 = transition;
					}
					
					if (t1 == 0)
					{
						mxIndex = t2-1;
					}
					else
					{
						if (t1 == 1)
						{
							mxIndex = 1+t2;
						}
						else
						{
							mxIndex = 5;
						}
					}
					
					t1 = (_Genetic_Code[0][h]!=_Genetic_Code[0][v]);
					rateTypeWeights[mxIndex][t1] = rateTypeWeights[mxIndex][t1]+
												   vectorOfFrequencies[h-hshift]*observedFreq[transition][nucPosInCodon]+
												   vectorOfFrequencies[v-vshift]*observedFreq[transition2][nucPosInCodon];
					rateTypeWeightsGY94[mxIndex][t1] = rateTypeWeightsGY94[mxIndex][t1]+
												   2*vectorOfFrequencies[h-hshift]*vectorOfFrequencies[v-vshift];
				}
		   }
	    }		
	}
	return 0;		
}

/*---------------------------------------------------------------------------------------------------------------------------------------------------*/

/* echoCatVar: calculate category mean and variance; output distribution info */ 
function echoCatVar (distrInfo)
{
	DD = Columns(distrInfo);
	EE = 0.0;
	sampleVar = 0.0;
	cumulative_weight = 0.0;
	median = -1;
	for (k=0; k<DD; k=k+1)
	{
		T = distrInfo[0][k]*distrInfo[1][k];
		EE = EE+T;
		sampleVar = T*distrInfo[0][k]+sampleVar;
		cumulative_weight = cumulative_weight + distrInfo[1][k];
		if (cumulative_weight >= 0.5 && median < 0)
		{
		    median = distrInfo[0][k];
		}
	}
		
	sampleVar = sampleVar-EE*EE;
	
	fprintf  (distribOutput,"\n------------------------------------------------\n\nSample mean = ",EE, " (sample variance = ",sampleVar,", median = ",median,")\n");

	for (k=0; k<DD; k=k+1)
	{
		fprintf (distribOutput,"\nRate[",Format(k,0,0),"]=",Format(distrInfo[0][k],12,8), " (weight=", 
						  Format(distrInfo[1][k],9,7),")");
	}
	fprintf  (distribOutput,"\n------------------------------------------------\n\n");
	return EE;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------------*/

function echoRatio (distrInfo, distrInfo2)
{
	D1 = Columns(distrInfo);
	D2 = Columns(distrInfo2);
		
	ratioInfo = {3,D1*D2};
	
	for (k=0; k<D1; k=k+1)
	{
		EE = k*D2;
		for (k2 = 0; k2<D2; k2=k2+1)
		{
			ratioInfo [0][EE+k2] = distrInfo[0][k]/distrInfo2[0][k2];
			ratioInfo [1][EE+k2] = distrInfo[1][k]*distrInfo2[1][k2];
			ratioInfo [2][EE+k2] = k2*D1+k;
		}
	}
	
	done = 0;
	EE = D1*D2;
	while (!done)
	{
		done = 1;
		for (k=1; k<EE; k=k+1)
		{
			if (ratioInfo [0][k] < ratioInfo[0][k-1])
			{
				DD = ratioInfo [0][k];
				ratioInfo [0][k] = ratioInfo [0][k-1];
				ratioInfo [0][k-1] = DD;
				DD = ratioInfo [1][k];
				ratioInfo [1][k] = ratioInfo [1][k-1];
				ratioInfo [1][k-1] = DD;
				DD = ratioInfo [2][k];
				ratioInfo [2][k] = ratioInfo [2][k-1];
				ratioInfo [2][k-1] = DD;
				done = 0;
			}
		}
	}
	return echoCatVar (ratioInfo);
}

/*---------------------------------------------------------------------------------------------------------------------------------------------------*/

/* F3x4 empirical frequency estimation; obsF is a 4x3 matrix containing position-specific nucleotide frequencies: */
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

/*---------------------------------------------------------------------------------------------------------------------------------------------------*/
fprintf (stdout, "\n+--------------------------------+\n",
                   "|     PARRIS SELECTION TEST      |\n",
                   "|          written by            |\n",
                   "|       Konrad Scheffler         |\n",
                   "|             and                |\n",
                   "|        Cathal Seoighe          |\n",
                   "|                                |\n",
                   "|    If you use this analysis    |\n",
                   "| in a publication, please cite  |\n",
                   "|  Bioinformatics, 22:2493Ð2499  |\n",
                   "+--------------------------------+\n");


#include "TemplateModels/chooseGeneticCode.def";
#include "_MFReader_.ibf";

/* _MFReader_ stores position-specific nucleotide frequencies in positionFrequencies (4x3)
   and overall nucleotide frequencies in overallFrequencies (4x1). */
observedFreq = positionFrequencies;

/* Input whether to optimise branch lengths in codon model or use pre-optimised nucleotide model branch lengths: */
ChoiceList (branchLengths,"Branch Lengths",1,SKIP_NONE,
			"Codon Model","Jointly optimize rate parameters and branch lengths (slow and thorough)",
			"Nucleotide Model","Estimate branch lengths once, using an appropriate nucleotide model (quick and dirty)."
		    );
		    
if (branchLengths<0)
{
	return 0;
}

/* Input type of underlying nucleotide model (optionally with "multi") and set up constraints: */
/* (In earlier versions these options were input as modelChoice using a single ChoiceList.) */
modelConstraintString = "";

ChoiceList (MGGYChoice, "Options for handling equilibrium frequencies",1,SKIP_NONE,
	                "MG", "Muse-Gaut 94: rates are proportional to position-specific target nucleotide frequency.",
	                "GY", "Goldman-Yang 94: rates are proportional to target codon frequency (F3x4 estimate)."
           );

if (MGGYChoice<0)
{
	return 0;
}

ChoiceList (nucModelChoice, "Nucleotide Rate Matrix Options",1,SKIP_NONE,
			"Jukes-Cantor",	 "Single rate type (all substitutions equally likely).",
			"HKY85","Two rate types (transition/transverion ratio parameter kappa).",
			"REV","Fully reversible model (6 rate types - recommended rate matrix)",
			"Custom","Arbitrary nucleotide reversible model, except F81."
           );

if (nucModelChoice<0)
{
	return 0;
}

ChoiceList (AAModelChoice, "Options for multiple classes of non-synonymous substitutions",1,SKIP_NONE,
                 	    "Single","Only one class of non-synonymous substitutions (standard model).",
	                    "Multi","Multiple classes of non-synonymous substitutions.",
	                    "NMulti","Multi with numerical bias corrections for various amino-acid substitutions."
           );

if (AAModelChoice<0)
{
	return 0;
}

modelInFile = "2RatesAnalyses/MG94GY94xREV_PARRIS_syn3.mdl";
if (MGGYChoice == 0)
{
    ModelTitle = "MG94";
}
else
{
    ModelTitle = "GY94";
}

if (nucModelChoice == 0)
{
	modelConstraintString = "AC:=1;AT:=1;CG:=1;CT:=1;GT:=1"; 
}
else
{
	if (nucModelChoice == 1)
	{
		modelConstraintString = "CT:=1;AT:=AC;CG:=AC;GT:=AC"; 
		ModelTitle = ModelTitle+"xHKY85";
	}
	else
	{
		if (nucModelChoice == 3)
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
			ModelTitle = ModelTitle+modelDesc[0];
						
			rateBiasTerms = {{"AC","1","AT","CG","CT","GT"}};
			paramCount	  = 0;

			modelConstraintString = "";

			for (customLoopCounter2=1; customLoopCounter2<6; customLoopCounter2=customLoopCounter2+1)
			{
				for (customLoopCounter=0; customLoopCounter<customLoopCounter2; customLoopCounter=customLoopCounter+1)
				{
					if (modelDesc[customLoopCounter2]==modelDesc[customLoopCounter])
					{
						ModelTitle  = ModelTitle+modelDesc[customLoopCounter2];	
						if (rateBiasTerms[customLoopCounter2] == "1")
						{
							modelConstraintString = modelConstraintString + rateBiasTerms[customLoopCounter]+":="+rateBiasTerms[customLoopCounter2]+";";
						}
						else
						{
							modelConstraintString = modelConstraintString + rateBiasTerms[customLoopCounter2]+":="+rateBiasTerms[customLoopCounter]+";";			
						}
						break;
					}
				}
				if (customLoopCounter==customLoopCounter2)
				{
					ModelTitle = ModelTitle+modelDesc[customLoopCounter2];	
				}
			}	
		}
		else
		{
		    ModelTitle = ModelTitle+"xREV";
		}
	}
			
	if (AAModelChoice > 0)
	{
	    fprintf(stdout, "Warning: running untested multi code!");
	    if (AAModelChoice == 2)
	    {
		ModelTitle = ModelTitle+"xNMulti";
		_AA_RM_NUMERIC = 2;
	    }
	    else
	    {
		ModelTitle = ModelTitle+"xMulti";
	    }
            #include "TemplateModels/MGwAA.ibf";
	    rateKeys = Rows (aaRateClassIDs);
	    userAARateMultipliers = {21,21};
	    if (AAModelChoice == 2)
	    {
		_AA_RM_NUMERIC = 0;
		for (h=0; h<21;h=h+1)
		{
		    for (v=0; v<21; v=v+1)
		    {
			userAARateMultipliers[h][v] = "R*" + aaRateMultipliers[h][v] + "*";
			userAARateMultipliers[v][h] = "R*" + aaRateMultipliers[v][h] + "*";
		    }
		}								
	    }
	    else
	    {
		for (aaIndex = 0; aaIndex < Abs(aaRateClassIDs); aaIndex = aaIndex+1)
		{
		    ExecuteCommands ("global DN_"+rateKeys[aaIndex]+" = 1;");
		}
		for (h=0; h<21;h=h+1)
		{
		    for (v=0; v<21; v=v+1)
		    {
			userAARateMultipliers[h][v] = "DN_" + aaRateMultipliers[h][v] + "*";
			userAARateMultipliers[v][h] = "DN_" + aaRateMultipliers[v][h] + "*";
		    }
		}								
	    }
	}
}

ChoiceList (rateVarModelChoice, "Rate Variation Models",1,SKIP_NONE,
			     "Constant","Constant Rate Model: no rate variation across sites.", /* index 0 */
			     "Proportional","Proportional Variable Rates Model: dS and dN vary along the sequence, but dN = R*dS for every site. Recommended model.",
			     "Nonsynonymous","Non-synonymous Variable Rates Model: dS = 1 for every site, while dN is drawn from a given distribution.",
			     "Dual","Dual Variable Rates Model: dS and dN are drawn from a bivariate distribution (independent or correlated components).",
			     "Lineage Dual","Lineage Dual Variable Rates Model:  dS and dN are drawn from a bivariate distribution (independent or correlated components), plus each lineage has an adjustment factor for the E[dN]/E[dS]."
			     );
			    
if (rateVarModelChoice<0)
{
    return 0;
}
	
/* For models with dual rate variation, input whether NonSynRate should be independent or a multiplicative funtion of synRate, and whether to model synonymous rate variation at the codon level (syn1) or the nucleotide level (syn3): */
multiplicativeNonSynRate = 0;
nucSynVar = 0;
if (rateVarModelChoice >= 3)
{
    ChoiceList (multiplicativeNonSynRate, "Independent or multiplicative nonsynonymous rate",1,SKIP_NONE,
		"Independent","Nonsynonymous rate is independent of synonymous rate.",
		"Multiplicative","Nonsynonymous rate is omega multiplied by synonymous rate."
	       );
    if (multiplicativeNonSynRate<0)
    {
	return 0;
    }

    ChoiceList (nucSynVar, "Codon or nucleotide level synonymous rate variation",1,SKIP_NONE,
		"Codon (syn1)","1 synonymous rate per codon.",
		"Nucleotide (syn3)","3 synonymous rates per codon.",
		/*"syn2", "2 independent synonymous rates per codon (2nd pos is avg of 1st and 3rd)"*/
	    );
    if (nucSynVar<0)
    {
	return 0;
    }    
}

/* Input list of models to be run (out of various rate variation options): */
/* (chosenModelList will contain nonzero entries for selected models) */
num_models = 8;
chosenModelList = {num_models,1};

ChoiceList (modelChoice,"Distribution Options",1,SKIP_NONE,
			"Run All","Run all available models.",
			"Run Custom","Choose some of the available models.",
	                "PARRIS","Run PARRIS model comparison.");
			
if (modelChoice<0)
{
	return 0;
}

if (modelChoice==2)
{
    for (mi = 0; mi<5; mi=mi+1)
    {
		chosenModelList[mi][0] = 0;
    }
    chosenModelList[5][0] = 1;
    chosenModelList[6][0] = 1;
    resp = 3; /* nr of synonymous rate classes (not used for PARRIS models except in output and postprocessing) */
    resp2 = 3; /* nr of nonsynonymous rate classes (in alternative model) (not used for PARRIS models except in output and postprocessing) */
	USE_LAST_RESULTS = 1;
}
else
{

/* Input nrs of rate classes: (presumably unused when they don't apply) */
    resp  = 0;
    resp2 = 0;

    while (resp<2)
    {
	fprintf (stdout,"Number of synonymous (and single variable rate models) rate classes (>=2):");
	fscanf  (stdin, "Number", resp);
    }

    while (resp2<2)
    {
	fprintf (stdout,"Number of non-synonymous rate classes (>=2):");
	fscanf  (stdin, "Number", resp2);
    }
			
    if (modelChoice==0)
    {
	for (mi = 0; mi<num_models; mi=mi+1)
	{
		chosenModelList[mi][0] = 1;
	}
    }
    else
    {

/* Input distribution options: */
/* M3 is similar to Independent Discrete, differences being that it fixes the nr of rate categories and that it
   allow use of syn3 files. Independent Discrete is the older (KP) implementation and does not support syn3. */
 
	ChoiceList (distribChoices, "Distribution Option",0,SKIP_NONE,
			"Syn:Gamma, Non-syn:Gamma",	 "Both syn and non-syn rates are drawn from the gamma distributions for all models.",
			"Syn:Gamma, Non-syn:Inv+Gamma","Rates are drawn from the gamma distributions for Proportional and Nonsynonymous. For Dual and Local Dual, syn rates are drawn from the gamma distribution, and non-syn rates - from Inv+Gamma.",
			"Independent Discrete", "Independent General Discrete Distributions (Recommended setting)",
			"Correlated Discrete", "Correlated General Discrete Distributions",
			"Non-positive Discrete", "General Discrete Distribution for dS, and dN, but constrained so that dN<=dS. Useful to perform a LRT for presence of selection in an alignment",
	                "M1a", "General Discrete Distribution for dS (3 cat), M1a omega distribution",
	                "M2a", "General Discrete Distribution for dS (3 cat), M2a omega distribution",
		        "M3", "General Discrete Distribution for dS (3 cat), M3 omega distribution (3 cat)");	

			
	if (distribChoices[0] < 0)
	{
	    return 0;
	}
	
	for (mi = 0; mi < Rows(distribChoices)*Columns(distribChoices); mi = mi + 1)
	{
	    modelChoice = distribChoices[mi];
	    chosenModelList[modelChoice] = 1;
	}	
    }
}
modelNamesShort = {{"GamGam","GamIGam","Discrete","CorrDiscrete","NPDiscrete","M1a","M2a", "M3"}};


/* Input whether to use default or randomised initial values for rate distribution parameters: */
ChoiceList (randomizeInitValues, "Initial Value Options",1,SKIP_NONE,
			"Default",	 "Use default inital values for rate distribution parameters.",
			"Randomized",	 "Select initial values for rate distribution parameters at random.");


if (randomizeInitValues < 0)
{
	return 0;
}

/* This variable is dereferenced in the unused function SetCodonNorm in the MG and GY model definition files, where it is used to scale codonFactor; keep it set to 1 and it will have no effect even if SetCodonNorm is used: */
fudgeFactor = 1.0;


SetDialogPrompt ("Save summary result file to:");
fprintf (PROMPT_FOR_FILE,CLEAR_FILE);
baselineOutput = LAST_FILE_PATH;

distribOutput = baselineOutput + ".distributions";
fprintf (distribOutput,CLEAR_FILE);


/* Don't know why this is negated: */
if (branchLengths)
{
	branchLengths = -branchLengths;
}


/* Output stuff: */
separator = "+--------------------+----------------+---------------+-----------------+------------------+-----------+-----+-----------+\n";

fprintf (stdout, "\n\nRUNNING ",ModelTitle," MODEL COMPARISONS on ",fileCount, " files\n\n",
	 "\n\n########### ",resp,"x",resp2," CLASSES ###########\n\n", separator,
	    "|       Model        | Log Likelihood | Synonymous CV |  NS Exp and CV  |  N/S Exp and CV  |  P-Value  | Prm |    AIC    |\n", separator);

fprintf (baselineOutput,  "\n\nRUNNING ",ModelTitle," MODEL COMPARISONS on ",dataFilePath, "\n\n",
	 "\n\n########### ",resp,"x",resp2," CLASSES ###########\n\n", separator,
	 "|       Model        | Log Likelihood | Synonymous CV |  NS Exp and CV  |  N/S Exp and CV  |  P-Value  | Prm |    AIC    |\n", separator);


modelNames = {{"| Gamma, Gamma       |",
	       "| Gamma, Inv+Gamma   |",
	       "| Independent Discr  |",
	       "| Correlated Discr   |",
	       "| Non-positive Discr |",
	       "| discr(3), M1a      |",
	       "| discr(3), M2a      |",
               "| discr(3), discr(3) |"}};
			   

if (MPI_NODE_COUNT>1)
/* Check to see if we are running with more than one MPI node.
   If not - bail to single CPU execution */			   
{
	MPINodeState = {MPI_NODE_COUNT-1,2};
	
	/* One of the nodes (0) will be the dispatcher, and
	   for every other node, we need to store two values:
	   whether it is busy or not (0 or 1), and 
	   an integer (from 0 to 5) to indicate which model
	   (M0-M5) the node is working on */
	   
	OPTIMIZE_SUMMATION_ORDER = 0;
	
	/* the master node will be creating likelihood functions, 
	but not evaluating them - thus there is no need to perform
	the extra step of optimizing data column ordering for faster
	likelihood evaluations */
}

doNucFit = 1;

for (midx = 0; midx<num_models; midx=midx+1)
{
	if (chosenModelList[midx])
	{

/* Depending on distribution option, read category definitions from file: */
/* for gamma models, read from file into categDef1 (for syn distribution) and categDef2 
   (for nonsyn distribution) variables;
   for discrete models, execute generator file which will create these variables (randomizeInitValues is used).
   The category definitions are string variables to be executed in the mdl files. */
	    if (midx<2) /* gamma rate models */
	    {
		fscanf ("2RatesAnalyses/gamma1.def","Raw",categDef1);

		if (midx == 0)
		{
		    fscanf ("2RatesAnalyses/gamma2.def","Raw",categDef2);
		}
		else
		{
		    fscanf ("2RatesAnalyses/gamma2+Inv.def","Raw",categDef2);
		}
	    }
	    else
	    {
		if (midx < 5) /* general discrete rate models */
		{
		    if (midx < 4)
		    {
			correlationOn = (distribChoice>2);
			fscanf ("2RatesAnalyses/discreteGenerator.bf","Raw",tempstr);
		    }
		    else
		    {
			fscanf ("2RatesAnalyses/discreteGeneratorNoPS.bf","Raw",tempstr);	
		    }
		    ExecuteCommands (tempstr);
		}
		else /* PARRIS rate models */
		{
		    if (nucSynVar)
		    {
			fscanf ("2RatesAnalyses/PARRIS_syn3.def","Raw",categDef1);
		    }
		    else
		    {
			fscanf ("2RatesAnalyses/PARRIS_synvar.def","Raw",categDef1);
		    }

		    if (midx == 5) /* null (M1a omega distrib) */
		    {
			fscanf ("2RatesAnalyses/PARRIS_M1.def","Raw",categDef2);
		    }
		    if (midx == 6) /* alternative (M2a omega distrib) */
		    {
			fscanf ("2RatesAnalyses/PARRIS_M2.def","Raw",categDef2);
		    }		    
		    if (midx == 7) /* PARRIS discrete (M3 omega distrib) */
		    {
			fscanf ("2RatesAnalyses/PARRIS_M3.def","Raw",categDef2);
		    }
		}
	    }

/* execute model file: */
	    ExecuteAFile (modelInFile);
	    
	    if (Abs(modelConstraintString))
	    {
			ExecuteCommands (modelConstraintString);
	    }

	    
	    theRateMatrix = 0;
	    
	    MULTIPLY_BY_FREQS = PopulateModelMatrix ("theRateMatrix", observedFreq, rateVarModelChoice);
	    vectorOfFrequencies = BuildCodonFrequencies (observedFreq);
	    
	    if (doNucFit) /* Fit nucleotide model to optimise branch lengths if this hasn't been done already: */
	    {
		ExecuteCommands (nucModelString+"\nModel nucModel = (nucModelMatrix,overallFrequencies);");
		
		populateTrees   ("nucTree", fileCount);
		ExecuteCommands (constructLF ("nuc_lf", "nucData", "nucTree", fileCount));
		Optimize (nuc_res, nuc_lf);
		
		computeExpSubWeights (0);		
		global codonFactor = 0.33;
		
		/* Count nr of branch lengths optimised: */
		totNumBranchLengths = 0;
		for (fileID = 1; fileID <= fileCount; fileID = fileID + 1)
		{
		    ExecuteCommands ("totNumBranchLengths = totNumBranchLengths + TipCount (nucTree_"+fileID+") + BranchCount (nucTree_"+fileID+");");
		}

		doNucFit = 0;
	    }
	    
	    Model theModel = (theRateMatrix,vectorOfFrequencies,MULTIPLY_BY_FREQS);
	    populateTrees     ("givenTree", fileCount);
	    
	    
	    for (fileID = 1; fileID <= fileCount; fileID = fileID + 1)
	    {
			if (branchLengths)
			{
				ExecuteCommands ("ClearConstraints (givenTree_" + fileID + ");");
				ExecuteCommands ("ReplicateConstraint (\"this1.?.synRate:=this2.?.t__/codonFactor\",givenTree_" + fileID + ",nucTree_" + fileID + ");");
			}
			else
			{
				bnames = Eval ("BranchName (givenTree_"+fileID+",-1)");
				for (lc = 0; lc < Columns (bnames) - 1; lc = lc+1)
				{
					ExecuteCommands ("givenTree_" + fileID + "." + bnames[lc] + ".synRate=nucTree_" + fileID + "." + bnames[lc] + ".t/codonFactor;");
				}
			}
	    }
	    
	    ExecuteCommands (constructLF ("lf", "filteredData", "givenTree", fileCount));
	    
	    if (midx >= 4) /* nonpositive discrete distribution and PARRIS models */
	    {
			R := 1;
	    }	
		
	    if (MPI_NODE_COUNT>1)
	    {
		/* look for an idle node */
		for (mpiNode = 0; mpiNode < MPI_NODE_COUNT-1; mpiNode = mpiNode+1)
		{
		    if (MPINodeState[mpiNode][0]==0)
		    {
			break;	
		    }
		}
		
		if (mpiNode==MPI_NODE_COUNT-1)
		    /* all nodes busy */
		{
		    /* wait for some node to complete and send out current job */
		    mpiNode = ReceiveJobs (1);
		}
		else
		{
		    /* send the job to an idle node; update node state */
		    MPISend (mpiNode+1,lf);
		    MPINodeState[mpiNode][0] = 1;
		    MPINodeState[mpiNode][1] = midx;
		}
	    }
	    else
	    {
		/* Non-MPI execution */
			Optimize (res,lf);
			modelIndex = midx;
			ReceiveJobs (0);
	    }
	}
}

if (MPI_NODE_COUNT>1)
/* wait for all the jobs to finish, process their results */
{
	while (1)
	{
		for (nodeCounter = 0; nodeCounter < MPI_NODE_COUNT-1; nodeCounter = nodeCounter+1)
		{
			if (MPINodeState[nodeCounter][0]==1)
			{
				fromNode = ReceiveJobs (0);
				break;	
			}
		}
		if (nodeCounter == MPI_NODE_COUNT-1)
		{
			break;
		}
	}	
	OPTIMIZE_SUMMATION_ORDER = 1;
}

/* ____________________________________________________________________________________________________________________*/

function ReceiveJobs (sendOrNot)

/* This function receieves and processes 
   model results. The parameter is a boolean,
   set to 1 if there are jobs waiting to be sent to
   an MPI node */
{
	if (MPI_NODE_COUNT>1)
	{

	    /* Fetch completed job (waiting if necessary); send new job if one is waiting; update node status; 
	       reconstruct likelihood function and result matrix from completed job as if it had run locally. */

		MPIReceive (-1, fromNode, result_String);
		modelIndex = MPINodeState[fromNode-1][1];
		
		if (sendOrNot)
		{
			/* send the likelihood function to the node which just finished */
			MPISend (fromNode,lf);
			MPINodeState[fromNode-1][1] = midx;
		}
		else
		{
			/* mark the node as idle */
			MPINodeState[fromNode-1][0] = 0;
			MPINodeState[fromNode-1][1] = 0;
		}
		
		/* reset category variables */
		/* (we are about to build the likelihood function of the completed job, 
		   which may be incompatible with constraints still in place from the previous lf) */

		if (rateVarModelChoice)
		{
			ClearConstraints (c);
		}
		if (rateVarModelChoice>2)
		{
			ClearConstraints (d);
			if (nucSynVar)
			{
			    ClearConstraints (c1,c2,c3);
			}
		}
		
		/* build the likelihood function lf with MLE parameter values
		   just returned from an MPI node. MLE matrix is returned 
		   in lf_MLES. */
		   		
		ExecuteCommands (result_String);
				
		/* store the result matrix */
		res = lf_MLES;
	}
	

	/* res[1][1] contains the number of independent parameters optimised. If branchlengths were fixed but
	   previously optimised on a nucleotide model, we need to add these, minus 1 for the codonFactor 
	   (tree length scaling) parameter: */
	if (branchLengths)
	{
		res[1][1] = res[1][1] + totNumBranchLengths - 1;
	}

	/* Initialise output files, output whole-sequence results: */
	LIKELIHOOD_FUNCTION_OUTPUT = 6;

	fprintf (stdout, modelNames[modelIndex], " ", Format (res[1][0],14,5));
	fprintf (baselineOutput, modelNames[modelIndex], " ", Format (res[1][0],14,5));
	
	fitOutFile = baselineOutput + "." + modelNamesShort[modelIndex] + ".fit";
       	fprintf (fitOutFile,CLEAR_FILE,lf);

	if (rateVarModelChoice) /* Create output file for marginals if a category variable exists */
	{
		marginalOutFile = baselineOutput + "." + modelNamesShort[modelIndex] + ".marginals";
		fprintf (marginalOutFile,CLEAR_FILE);
	}
	
	fprintf (distribOutput, modelNames[modelIndex]);
	if (rateVarModelChoice==0) /* No site-to-site rate variation: just output mean omega (variance of category variable is zero when only a single category exists */
	{
		EE = 0;
		sampleVar = 0;
		fprintf (stdout, " |      N/A      | ", Format (R,7,5) , ",", Format (sampleVar,7,5), " | ", Format (R,7,5) , ",", Format (sampleVar,8,5), " | ");
		fprintf (baselineOutput, " |      N/A      | ", Format (R,7,5) , ",", Format (sampleVar,7,5), " | ", Format (R,7,5) , ",", Format (sampleVar,8,5), " | ");
	}
	else
	{
	    /* For models with category variables, calculate and output category mean and variance information: */
		if (rateVarModelChoice==2) /* Nonsynonymous rate variation only */
		{
		        GetInformation(NSdistrInfo,d);
			DD = Columns (NSdistrInfo);
			for (EE=0; EE<DD; EE=EE+1)
			{
				NSdistrInfo[0][EE] = R*NSdistrInfo[0][EE];
			}
			
			EE  = echoCatVar (NSdistrInfo);
			fprintf (marginalOutFile, NSdistrInfo);
			
			fprintf (stdout, " |      N/A      | ", Format (EE,7,5) , ",", Format (Sqrt(sampleVar)/EE,7,5), " | ", Format (EE,7,5) , ",", Format (Sqrt(sampleVar)/EE,8,5), " | ");
			fprintf (baselineOutput, " |      N/A      | ", Format (EE,7,5) , ",", Format (Sqrt(sampleVar)/EE,7,5), " | ", Format (EE,7,5) , ",", Format (Sqrt(sampleVar)/EE,8,5), " | ");
		}  
		else /* Model has synonymous rate variation */
		{
		    if (nucSynVar) /*syn3 */
		    {
			GetInformation(synDistrInfo,c1); /* for syn3 all 3 category variables are identical, so still only need to call this for one of them */
		    }
		    else /* syn1 */
		    {
			GetInformation(synDistrInfo,c);
		    }
		    fprintf (marginalOutFile, synDistrInfo);
		    EE  = echoCatVar (synDistrInfo);
		    synmedian = median;
		    fprintf (stdout, " | ",Format (Sqrt(sampleVar),13,8)," | ");
		    fprintf (baselineOutput, " | ",Format (Sqrt(sampleVar),13,8)," | ");
			
		    if (rateVarModelChoice>=3) /* Synonymous and nonsynonymous rate variation */
		    {
			GetInformation(NSdistrInfo,d);
			DD = Columns (NSdistrInfo);
			for (EE=0; EE<DD; EE=EE+1)
			{
			    NSdistrInfo[0][EE] = R*NSdistrInfo[0][EE];
			}
			fprintf (marginalOutFile, NSdistrInfo);
			EEN  = echoCatVar (NSdistrInfo);
			varN = sampleVar;
			EER  = echoRatio  (NSdistrInfo,synDistrInfo);
			varR = sampleVar;
			
			fprintf (stdout,  Format (EEN,7,5) , ",", Format (Sqrt(varN)/EEN,7,5), " | ",Format (EER,7,5) , ",", Format (Sqrt(varR)/EER,8,5), " | ");	
			fprintf (baselineOutput, Format (EEN,7,5) , ",", Format (Sqrt(varN)/EEN,7,5), " | ",Format (EER,7,5) , ",", Format (Sqrt(varR)/EER,8,5), " | ");	
					
		    }
		    else /* Proportional model: Nonsynonymous rate variation is just a scaled version of synonymous rate variation. */
		    {
			fprintf (stdout, Format (EE*R,7,5) , ",", Format (Sqrt(sampleVar)/EE,7,5), " | ", Format (R,7,5) , ",", Format (0,8,5), " | ");
			fprintf (baselineOutput, Format (EE*R,7,5) , ",", Format (Sqrt(sampleVar)/EE,7,5), " | ", Format (R,7,5) , ",", Format (0,8,5), " | ");
		    }
		}
	}
		
	/* Likelihood ratio tests: */
	if (modelIndex == 5)
	{
	    nullLL = res[1][0];
	    nullPC = res[1][1];
	}
	/* we don't know that alternative returned before null, so we won't print the p-value in the MPI mode: */
	if ((modelIndex == 6)&&(MPI_NODE_COUNT<=1))
	{
	    fprintf (stdout, Format (1-CChi2(2*(res[1][0]-nullLL),res[1][1]-nullPC),9,7), " |");
	    fprintf (baselineOutput, Format (1-CChi2(2*(res[1][0]-nullLL),res[1][1]-nullPC),9,7), " |");
	}
	else
	{
	    fprintf (stdout,"   N/A    |");
	    fprintf (baselineOutput, "   N/A    |");
	}

	/* Output parameter count and AIC scores: */
	fprintf (stdout, Format (res[1][1],5,0), "|", Format (2*(res[1][1]-res[1][0]),11,2),"|");
	fprintf (baselineOutput, Format (res[1][1],5,0), "|", Format (2*(res[1][1]-res[1][0]),11,2),"|");
	fprintf (stdout, "\n",separator);
	fprintf (baselineOutput, "\n",separator);

	/* Site-specific posterior analysis: */
	if (rateVarModelChoice  > 1 && nucSynVar < 2) /* Models with independent nonsynonymous rate variation */
	{
	    sigLevel = 0.95;
	    D = Columns(NSdistrInfo);
	    for (k=0; k<D; k=k+1)
	    {
		if (NSdistrInfo[0][k]>1) break;
	    }
	    ConstructCategoryMatrix(site_likelihoods,lf,COMPLETE);
	    CC = Columns (site_likelihoods);
	    DD = Rows (site_likelihoods);
	    /* set up variable order information:
	       (ordernum[n] is the number of the nth variable in site_likelihoods, 
	       where omega (d) is variable 0 and c (synRate, for syn1) or 
	       c1, c2, c3 (for syn3) are variables 1 or 1-3): */
	    GetInformation (order_ID, lf);
	    if (rateVarModelChoice < 3)
	    {
		ordernum = {{0}};
		order0 = 0;
	    } else {
		if (!nucSynVar) /* syn1 */
		{
		    if (order_ID[0] == "d")
		    {
			ordernum = {{0,1}}; /* first variable is nonsyn (0), second is syn (1) */
		    } else {
			ordernum = {{1,0}}; /* first variable is syn (1), second is nonsyn (0) */
		    }
		    order0 = ordernum[0];
		    order1 = ordernum[1];
		}
		else /* syn3 */
		{
		    if (multiplicativeNonSynRate)
		    {
			ordernum = {1,4};
			for (posnum = 0; posnum < 4; posnum = posnum+1)
			{
			    if (order_ID[posnum] == "d")
			    {
				ordernum[posnum] = 0;
			    }
			    if (order_ID[posnum] == "c1")
			    {
				ordernum[posnum] = 1;
			    }
			    if (order_ID[posnum] == "c2")
			    {
				ordernum[posnum] = 2;
			    }
			    if (order_ID[posnum] == "c3")
			    {
				ordernum[posnum] = 3;
			    }
			}
			order0 = ordernum[0];
			order1 = ordernum[1];
			order2 = ordernum[2];
			order3 = ordernum[3];
		    }
		    else /* Independent syn3: c2 is not used. */
		    {
			ordernum = {1,3};
			for (posnum = 0; posnum < 3; posnum = posnum+1)
			{
			    if (order_ID[posnum] == "d")
			    {
				ordernum[posnum] = 0;
			    }
			    if (order_ID[posnum] == "c1")
			    {
				ordernum[posnum] = 1;
			    }
			    if (order_ID[posnum] == "c3")
			    {
				ordernum[posnum] = 2;
			    }
			}
			order0 = ordernum[0];
			order1 = ordernum[1];
			order2 = ordernum[2];
		    }
		}
	    }

	    if (k<D)
		/* have rates > 1 */
	    {
		fprintf  (marginalOutFile,"\n\n------------------------------------------------\n\n Sites with dN/dS > 1 (Posterior cutoff = ",sigLevel,")\n\n");
	    }
	    else
	    {
		fprintf  (marginalOutFile,"\n\n------------------------------------------------\n\n No rate classes with dN/dS>1.");
	    }

	    /* allocate memory for data structures: */
	    site_posteriors = {DD,CC};
	    negprobs = {1,CC}; /* vector for storing site-specific negative selection probs */
	    negsynprobs = {1,CC}; /* vector for storing site-specific negative synonymous selection probs */
	    negsynprobs_1 = {1,CC};
	    negsynprobs_2 = {1,CC};
	    negsynprobs_3 = {1,CC};
	    omega_marginals = {D,CC};
	    synrate_marginals = {resp,CC};
	    synrate_1_marginals = {resp,CC};
	    synrate_2_marginals = {resp,CC};
	    synrate_3_marginals = {resp,CC};

	    for (v=0; v<CC; v=v+1) /* loop over sites */
	    {
		totProb = 0;
		positiveProb = 0;
		negativeProb = 0;
		negSynProb = 0;
		negSynProb_1 = 0;
		negSynProb_2 = 0;
		negSynProb_3 = 0;

		/* initialise index and dimensions variables: */
		indexes = {{0, 0, 0, 0}};
		if (rateVarModelChoice < 3)
		{
		    dimensions = {{D}};
		} 
		else 
		{
		    if (!nucSynVar) /* syn1 */
		    {
			dimensions = {{D, resp}};
		    }
		    else /* syn3 */
		    {
			dimensions = {{D, resp, resp, resp}};
		    }
		}   

		/* loop through all category combinations */
		for (h=0; h<DD; h=h+1) 
		{
		    idx0 = indexes[0]; /* index for nonsynonymous rate */
		    idx1 = indexes[1]; /* index for synonymous rate 1, etc */
		    idx2 = indexes[2];
		    idx3 = indexes[3];

		    if (rateVarModelChoice < 3)
		    {
			site_posteriors[h][v] = NSdistrInfo[1][h]*site_likelihoods[h][v];
		    }
		    else
		    {
			if (!nucSynVar) /* syn1 */
			{
			    site_posteriors[h][v] = NSdistrInfo[1][idx0]*synDistrInfo[1][idx1]*site_likelihoods[h][v];
			}
			else /* syn3 */
			{
			    if (multiplicativeNonSynRate)
			    {
				site_posteriors[h][v] = NSdistrInfo[1][idx0]*synDistrInfo[1][idx1]*synDistrInfo[1][idx2]*synDistrInfo[1][idx3]*site_likelihoods[h][v];
			    }
			    else
			    {
				site_posteriors[h][v] = NSdistrInfo[1][idx0]*synDistrInfo[1][idx1]*synDistrInfo[1][idx2]*site_likelihoods[h][v];
			    }				
			}
		    }
		    totProb = totProb + site_posteriors[h][v]; /* accumulate conditional posteriors, summing to probability of data at this site */

		    /* Accumulate all marginals of interest: */
		    if (NSdistrInfo[0][idx0]>1) /* positive AA selection category */
		    {
			positiveProb = positiveProb + site_posteriors[h][v];
		    }
		    if (NSdistrInfo[0][idx0]<1) /* purifying AA selection category */
		    {
			negativeProb = negativeProb + site_posteriors[h][v];
		    }
		    if (rateVarModelChoice >= 3)
		    {
			if (!nucSynVar) /* syn1 */
			{
			    if (synDistrInfo[0][idx1]<synmedian) /* purifying nuc-level (synonymous) selection category (whole codon) */
			    {
				negSynProb = negSynProb + site_posteriors[h][v];
			    }
			    omega_marginals[idx0][v] = omega_marginals[idx0][v] + site_posteriors[h][v];
			    synrate_marginals[idx1][v] = synrate_marginals[idx1][v] + site_posteriors[h][v];
			}
			else /* syn3 */
			{
			    if (synDistrInfo[0][idx1]<synmedian) /* purifying nuc-level (synonymous) selection category (position 1) */
			    {
				negSynProb_1 = negSynProb_1 + site_posteriors[h][v];
			    }
			    if (synDistrInfo[0][idx2]<synmedian) /* purifying nuc-level (synonymous) selection category (position 2) */
			    {
				negSynProb_2 = negSynProb_2 + site_posteriors[h][v];
			    }
			    if (multiplicativeNonSynRate)
			    {
				if (synDistrInfo[0][idx3]<synmedian) /* purifying nuc-level (synonymous) selection category (position 3) */
				{
				    negSynProb_3 = negSynProb_3 + site_posteriors[h][v];
				}
			    }
			    omega_marginals[idx0][v] = omega_marginals[idx0][v] + site_posteriors[h][v];
			    synrate_1_marginals[idx1][v] = synrate_1_marginals[idx1][v] + site_posteriors[h][v];
			    synrate_2_marginals[idx2][v] = synrate_2_marginals[idx2][v] + site_posteriors[h][v];
			    if (multiplicativeNonSynRate)
			    {
				synrate_3_marginals[idx3][v] = synrate_3_marginals[idx3][v] + site_posteriors[h][v];
			    }
			}
		    }

		    /* keep track of individual category variable indices: */
		    if (rateVarModelChoice < 3)
		    {
			indexes[0] = indexes[0]+1;
		    } else {
			if (!nucSynVar) /* syn1 */
			{
			    indexes[order1] = indexes[order1]+1; /* increment index for second variable in category matrix */
			    if (indexes[order1] == dimensions[order1])
			    {
				indexes[order1] = 0;
				indexes[order0] = indexes[order0]+1;
			    }
			}
			else /* syn3 */
			{
			    if (multiplicativeNonSynRate)
			    {
				indexes[order3] = indexes[order3]+1;
				if (indexes[order3] == dimensions[order3])
				{
				    indexes[order3] = 0;
				    indexes[order2] = indexes[order2]+1;
				    if (indexes[order2] == dimensions[order2])
				    {
					indexes[order2] = 0;
					indexes[order1] = indexes[order1]+1;
					if (indexes[order1] == dimensions[order1])
					{
					    indexes[order1] = 0;
					    indexes[order0] = indexes[order0]+1;
					}
				    }
				}
			    }
			    else /* independent syn3: only 3 category variables */
			    {
				indexes[order2] = indexes[order2]+1;
				if (indexes[order2] == dimensions[order2])
				{
				    indexes[order2] = 0;
				    indexes[order1] = indexes[order1]+1;
				    if (indexes[order1] == dimensions[order1])
				    {
					indexes[order1] = 0;
					indexes[order0] = indexes[order0]+1;
				    }
				}
			    }
			}
		    }
		} /* for h */

		positiveProb = positiveProb/totProb;
		negprobs[0][v] = negativeProb/totProb;
		if (rateVarModelChoice >= 3)
		{
		    if (!nucSynVar) /* syn1 */
		    {
			negsynprobs[0][v] = negSynProb/totProb;
			for (h = 0; h < D; h=h+1)
			{
			    omega_marginals[h][v] = omega_marginals[h][v]/totProb;
			}
			for (h = 0; h < resp; h=h+1)
			{
			    synrate_marginals[h][v] = synrate_marginals[h][v]/totProb;
			}
		    }
		    else /* syn3 */
		    {
			negsynprobs_1[0][v] = negSynProb_1/totProb;
			negsynprobs_2[0][v] = negSynProb_2/totProb;
			if (multiplicativeNonSynRate)
			{
			    negsynprobs_3[0][v] = negSynProb_3/totProb;
			}
			for (h = 0; h < D; h=h+1)
			{
			    omega_marginals[h][v] = omega_marginals[h][v]/totProb;
			}
			for (h = 0; h < resp; h=h+1)
			{
			    synrate_1_marginals[h][v] = synrate_1_marginals[h][v]/totProb;
			    synrate_2_marginals[h][v] = synrate_2_marginals[h][v]/totProb;
			    if (multiplicativeNonSynRate)
			    {
				synrate_3_marginals[h][v] = synrate_3_marginals[h][v]/totProb;
			    }
			}
		    }
		}

		for (h=0; h<DD; h=h+1)
		{
		    site_posteriors[h][v] = site_posteriors[h][v]/totProb;
		}
		if (positiveProb>=sigLevel)
		{
		    fprintf (marginalOutFile,Format (v+1,0,0)," (",positiveProb,")\n");
		}
	    } /* for v */

	    fprintf  (marginalOutFile,"\n\n------------------------------------------------\n\n Sites with dN/dS < 1 (Posterior cutoff = ",sigLevel,")\n\n");
	    for (v=0; v<CC; v=v+1)
	    {
		if (negprobs[0][v]>=sigLevel)
		{
		    fprintf (marginalOutFile,Format (v+1,0,0)," (",negprobs[0][v],")\n");
		}
	    }

	    if (rateVarModelChoice >= 3)
	    {
		if (!nucSynVar) /* syn1 */
		{
		    fprintf  (marginalOutFile,"\n\n------------------------------------------------\n\n Codon sites under purifying synonymous selection (Posterior cutoff = ",sigLevel,")\n\n");
		    for (v=0; v<CC; v=v+1)
		    {
			if (negsynprobs[0][v]>=sigLevel)
			{
			    fprintf (marginalOutFile,Format (v+1,0,0)," (",negsynprobs[0][v],")\n");
			}
		    }
		}
		else /* syn3 */
		{
		    fprintf  (marginalOutFile,"\n\n------------------------------------------------\n\n Nucleotide sites under purifying synonymous selection (Posterior cutoff = ",sigLevel,")\n\n");
		    fprintf  (marginalOutFile,"Position 1:\n\n");
		    for (v=0; v<CC; v=v+1)
		    {
			if (negsynprobs_1[0][v]>=sigLevel)
			{
			    fprintf (marginalOutFile,Format (v+1,0,0)," (",negsynprobs_1[0][v],")\n");
			}
		    }
		    if (multiplicativeNonSynRate)
		    {
			fprintf  (marginalOutFile,"\nPosition 2:\n\n");
		    }
		    else
		    {
			fprintf  (marginalOutFile,"\nPosition 3:\n\n");
		    }
		    for (v=0; v<CC; v=v+1)
		    {
			if (negsynprobs_2[0][v]>=sigLevel)
			{
			    fprintf (marginalOutFile,Format (v+1,0,0)," (",negsynprobs_2[0][v],")\n");
			}
		    }
		    if (multiplicativeNonSynRate)
		    {
			fprintf  (marginalOutFile,"\nPosition 3:\n\n");
			for (v=0; v<CC; v=v+1)
			{
			    if (negsynprobs_3[0][v]>=sigLevel)
			    {
				fprintf (marginalOutFile,Format (v+1,0,0)," (",negsynprobs_3[0][v],")\n");
			    }
			}
		    }
		}
	    }


/* tabulate posteriors for omega and synRate rate classes */

	    fprintf  (marginalOutFile,"\n\n------------------------------------------------\n\nTabulated Posteriors for Each Site\nSites in Columns, Rate Classes in Rows\nImport the following part into a data processing program\nfor further analysis\n\n");
	    if (rateVarModelChoice < 3)
	    {
		outString = "Omega/Site";
		for (v=0; v<D; v=v+1)
		{
		    outString = outString + "\tOmega=" + NSdistrInfo[0][v];
		}	
	    
		for (v=0; v<CC; v=v+1)
		{
		    outString = outString + "\n" + (v+1);
		    for (h=0; h<D; h=h+1)
		    {
			outString = outString + "\t" + site_posteriors[h][v];
		    }
		}
	    }	    

	    if (rateVarModelChoice >= 3)
	    {
		outString = "Rates/Site";
		for (v=0; v<D; v=v+1)
		{
		    outString = outString + "\tOmega=" + NSdistrInfo[0][v];
		}	
		if (!nucSynVar) /* syn1 */
		{
		    for (v=0; v<resp; v=v+1)
		    {
			outString = outString + "\tS=" + synDistrInfo[0][v];
		    }	
		    
		    for (v=0; v<CC; v=v+1)
		    {
			outString = outString + "\n" + (v+1);
			for (h=0; h<D; h=h+1)
			{
			    outString = outString + "\t" + omega_marginals[h][v];
			}
			for (h=0; h<resp; h=h+1)
			{
			    outString = outString + "\t" + synrate_marginals[h][v];
			}
		    }
		}
		else /* syn3 */
		{
		    for (v=0; v<resp; v=v+1)
		    {
			outString = outString + "\tS_pos1=" + synDistrInfo[0][v];
		    }	
		    if (multiplicativeNonSynRate)
		    {
			for (v=0; v<resp; v=v+1)
			{
			    outString = outString + "\tS_pos2=" + synDistrInfo[0][v];
			}	
		    }
		    for (v=0; v<resp; v=v+1)
		    {
			outString = outString + "\tS_pos3=" + synDistrInfo[0][v];
		    }	
	    
		    for (v=0; v<CC; v=v+1)
		    {
			outString = outString + "\n" + (v+1);
			for (h=0; h<D; h=h+1)
			{
			    outString = outString + "\t" + omega_marginals[h][v];
			}
			for (h=0; h<resp; h=h+1)
			{
			    outString = outString + "\t" + synrate_1_marginals[h][v];
			}
			for (h=0; h<resp; h=h+1)
			{
			    outString = outString + "\t" + synrate_2_marginals[h][v];
			}
			if (multiplicativeNonSynRate)
			{
			    for (h=0; h<resp; h=h+1)
			    {
				outString = outString + "\t" + synrate_3_marginals[h][v];
			    }
			}
		    }
		}
	    }	    
	    fprintf  (marginalOutFile,"\n", outString,"\n");

	    site_posteriors = 0;
	} /* if rateVarModelChoice */
	
	return fromNode-1;
}
