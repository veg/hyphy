/* 1. include a file to define the genetic code
   Note the use of base directory and path forming variables to make this analysis 
   independent of directory placement
 */

ExecuteAFile  (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"chooseGeneticCode.def");
ExecuteAFile  (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Utility"+DIRECTORY_SEPARATOR+"GrabBag.bf");

/* 2. load a codon partition  */

SetDialogPrompt 			("Please locate a coding alignment:");
DataSet 	  ds		   = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);
coding_path = LAST_FILE_PATH;

fprintf (stdout, "\nLoaded a ", filteredData.species, " sequence alignment with ", filteredData.sites, " codons from\n",coding_path,"\n");

/* 3. include a file to prompt for a tree */

incFileName = HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"queryTree.bf";
ExecuteCommands  ("#include \""+incFileName+"\";");

/* 4. Compute nucleotide counts by position for the F3x4 estimator */

COUNT_GAPS_IN_FREQUENCIES = 0;
HarvestFrequencies (baseFreqs,filteredData,3,1,1);
COUNT_GAPS_IN_FREQUENCIES = 1;

fprintf (stdout, "\nBase composition:\n\tA: ", Format (baseFreqs[0][0],10,5),",",Format (baseFreqs[0][1],10,5),",",Format (baseFreqs[0][2],10,5),
								    "\n\tC: ", Format (baseFreqs[1][0],10,5),",",Format (baseFreqs[1][1],10,5),",",Format (baseFreqs[1][2],10,5), 
									"\n\tG: ", Format (baseFreqs[2][0],10,5),",",Format (baseFreqs[2][1],10,5),",",Format (baseFreqs[2][2],10,5), 
									"\n\tT: ", Format (baseFreqs[3][0],10,5),",",Format (baseFreqs[3][1],10,5),",",Format (baseFreqs[3][2],10,5), "\n");
										  


/* 6. define the GY94 rate matrix; for now each branch will have it's own
   dS and dN, we will constrain them later */

global kappa_inv = 1;

ModelMatrixDimension = 64;
for (h = 0; h<64; h=h+1) 
{
	if (_Genetic_Code[h]==10) /* stop codon */
	{
		ModelMatrixDimension = ModelMatrixDimension-1;
	}
}

GY_Matrix = {ModelMatrixDimension,ModelMatrixDimension};
hshift = 0;
for (h=0; h<64; h=h+1)
{
	if (_Genetic_Code[h]==10) 
	{
		hshift = hshift+1;
	}
	else
	{
		vshift = hshift;
		for (v = h+1; v<64; v=v+1)
		{
			diff = v-h;
			if (_Genetic_Code[v]==10) 
			{
				vshift = vshift+1;
			}
			else
			{
			  	if ((h$4==v$4)||((diff%4==0)&&(h$16==v$16))||(diff%16==0)) /* one step */
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
			  			}
			  			else
			  			{
			  				transition = v%16$4;
			  				transition2= h%16$4;
			  			}
			  		}
			  		if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) /* synonymous */
			  		{
			  			if (Abs(transition-transition2)%2) /* transversion */
			  			{
			  				GY_Matrix[h-hshift][v-vshift] := kappa_inv*synRate;
			  				GY_Matrix[v-vshift][h-hshift] := kappa_inv*synRate;
			  			}
			  			else
			  			{
			  				GY_Matrix[h-hshift][v-vshift] := synRate;
			  				GY_Matrix[v-vshift][h-hshift] := synRate;			  			
			  			}
				  	}
			  		else
			  		{
			  			if (Abs(transition-transition2)%2) /* transversion */
			  			{
			  				GY_Matrix[h-hshift][v-vshift] := kappa_inv*nonSynRate;
			  				GY_Matrix[v-vshift][h-hshift] := kappa_inv*nonSynRate;
			  			}
			  			else
			  			{
			  				GY_Matrix[h-hshift][v-vshift] := nonSynRate;
			  				GY_Matrix[v-vshift][h-hshift] := nonSynRate;			  			
			  			}
		  			}
			  	}
			 }
		 }
	}	
}

/*8. build codon frequencies (use the F3x4 estimator) */

PIStop = 1.0;
codonFreqs = {ModelMatrixDimension,1};
hshift = 0;

for (h=0; h<64; h=h+1)
{
	first  = h$16;
	second = h%16$4;
	third  = h%4;
	if (_Genetic_Code[h]==10) 
	{
		hshift = hshift+1;
		PIStop = PIStop-baseFreqs[first][0]*baseFreqs[second][1]*baseFreqs[third][2];
		continue; 
	}
	codonFreqs[h-hshift]=baseFreqs[first][0]*baseFreqs[second][1]*baseFreqs[third][2];
}

codonFreqs = codonFreqs*(1.0/PIStop);

/*9. define the codon model */

Model 		CodonModel = (GY_Matrix,codonFreqs,1);

/*10. Approximate kappa and branch lengths using an HKY85 fit */

HKY85_Matrix = {{*,t*kappa_inv,t,t*kappa_inv}
				{t*kappa_inv,*,kappa_inv*t,t}
				{t,t*kappa_inv,*,kappa_inv*t}
				{t*kappa_inv,t,kappa_inv*t,*}};
			
HarvestFrequencies (nucFreqs,ds,1,1,1);
Model HKY85_Model = (HKY85_Matrix,nucFreqs,1);

Tree		  nucTree = treeString;
DataSetFilter nucData = CreateFilter (ds,1);

fprintf 					 (stdout, "Obtaining nucleotide branch lengths and kappa to be used as starting values...\n");
LikelihoodFunction	nuc_lf = (nucData,nucTree);
Optimize					 (nuc_mle,nuc_lf);
fprintf 					 (stdout, "\n", Format (nucTree,1,1), "\nkappa=", Format (1/kappa_inv,8,3), "\n");
nucBL					   = BranchLength (nucTree,-1);
nucBN					   = BranchName	  (nucTree,-1);
nucBC					   = Columns (nucBN) - 1;

sampleCount				   = filteredData.sites*filteredData.species;
baseParamCount			   = nuc_mle[1][1] + 9;

fixGlobalParameters		  ("nuc_lf");

/* work out the contributions of synRate and non-synRate to branch lengths */

synMultiplier = computeScaling (1,0);
nsMultiplier  = computeScaling (0,1);

fprintf						 (stdout, "Fitting (c-AIC) step-up local models to determine the appropriate number (and proportions) of rate classes\n");

bestBIC					   = 1e100;
bestBIC_classes			   = 0;
currentRate				   = 1;

currentLFDef			   = "LikelihoodFunction codonLF = (";

UseModel				   (CodonModel);
USE_LAST_RESULTS		   = 1;
logLByRateClassC		   = {};
p_thresh				   = 0.01;

while(1)
{
	fprintf 			(stdout, "Model with ", currentRate, " classes\n");
	/* define the new tree */
	
	ExecuteCommands 	("Tree codon_tree_" + currentRate + " = " + treeString);
	/* copy nucleotide branch lengths to this tree, using a random [0.05-1.95] dN/dS assignment */
	
	
	if (currentRate>1)
	{
		currentLFDef	  = currentLFDef+",";
		for (k = 0; k < nucBC; k = k+1)
		{
			randomProp = Random (0.05, 1.95);
			ExecuteCommands ("codon_tree_" + currentRate+ "." + nucBN[k] + ".synRate       =  3*" + codonBL[k] + "/(synMultiplier + nsMultiplier)");
			ExecuteCommands ("codon_tree_" + currentRate+ "." + nucBN[k] + ".nonSynRate    := ((codon_tree_1." + nucBN[k] + ".synRate - codon_tree_" + currentRate + "." + nucBN[k] + ".synRate)* synMultiplier__ + codon_tree_1." + nucBN[k] + ".nonSynRate * nsMultiplier__)/nsMultiplier__;");
		}	
	}
	else
	{
		for (k = 0; k < nucBC; k = k+1)
		{
			randomProp = Random (0.05, 1.95);
			ExecuteCommands ("codon_tree_1." + nucBN[k] + ".synRate    = 3*" + nucBL[k] + "/(synMultiplier + randomProp*nsMultiplier)");
			ExecuteCommands ("codon_tree_1." + nucBN[k] + ".nonSynRate = codon_tree_1." + nucBN[k] + ".synRate * randomProp");
		}	
	}
	currentLFDef 		  = currentLFDef + "filteredData,codon_tree_" + currentRate;
	if (currentRate > 1)
	{
		generate_gdd_freqs (currentRate, "mixingFreq", "lfMix", "PS", 1);
		ExecuteCommands (currentLFDef + ",\"" + lfMix + "\");");
	}
	else
	{
		ExecuteCommands (currentLFDef + ");");
		mixingFreq = {{"1"}};
	}
	Optimize (resCurrent, codonLF);
	codonBL = BranchLength (codon_tree_1, -1);
	myDF	= baseParamCount+nucBC*currentRate;
	myBIC 	= 2*(myDF*sampleCount/(sampleCount-myDF-1)-resCurrent[1][0]);
	fprintf (stdout, codonLF);
	
	logLByRateClassC [currentRate] = resCurrent[1][0];
	
	if (currentRate > 1)
	{
		LRT = 2(logLByRateClassC [currentRate]-logLByRateClassC [currentRate-1]);
		pp  = 1-CChi2(LRT,nucBC+1);
		fprintf (stdout, "\nAdvance rate class count p-value = ", pp, "\n");
		if (pp > p_thresh)
		{
			break;	
		}
	}

	bestBIC_classes = currentRate;
	bestFreqs		= mixingFreq;
	Export 			  (bestLF, codonLF);
	currentRate 		  = currentRate+1;	
}

ExecuteCommands (bestLF);

LFCompute (codonLF,LF_START_COMPUTE);
LFCompute (codonLF,alt_res);
LFCompute (codonLF,LF_DONE_COMPUTE);

hasPS    = {};
numFreqs = {bestBIC_classes,1};

SetDialogPrompt ("Save to:");

fprintf  		(PROMPT_FOR_FILE, CLEAR_FILE, KEEP_OPEN);
summaryOut = 	LAST_FILE_PATH;

fprintf (summaryOut, "Branch,LogL,p1,p2");

for (branchID   = 0; branchID < bestBIC_classes; branchID = branchID + 1)
{
	ExecuteCommands ("numFreqs[branchID] = " + bestFreqs[branchID] + ";");
	fprintf (summaryOut, ",synRate_",branchID+1,",nonSynRate_",branchID+1);
}

fprintf (summaryOut, "\nSummary,",alt_res,",,");

for (branchID   = 0; branchID < bestBIC_classes; branchID = branchID + 1)
{
	ExecuteCommands ("numFreqs[branchID] = " + bestFreqs[branchID] + ";");
	fprintf (summaryOut, ",", numFreqs[branchID], ",");
}




nullLogLByBranch = {};


for (branchID   = 0; branchID < nucBC; branchID = branchID + 1)
{
	thisBranchName = nucBN [branchID];

	fprintf (summaryOut, "\n", thisBranchName);

	fprintf (stdout, "Branch ", thisBranchName, "\n");
	
	rateString = "";
	rateString * 128;
	
	for (classCount = 1; classCount <= bestBIC_classes; classCount = classCount + 1)
	{
		ExecuteCommands ("beta = codon_tree_" + classCount + "." + thisBranchName + ".nonSynRate;alpha = codon_tree_" + classCount + "." + thisBranchName + ".synRate;");
		if (beta>alpha)
		{
			hasPS[thisBranchName] = hasPS[thisBranchName]+1;
		}
		rateString * ("," + alpha + "," + beta);
		fprintf (stdout, "\tClass ", classCount, " (p=", numFreqs[classCount-1], ")\n\t\tdN/dS = ", beta/alpha, "\n");
	}
	rateString * 0;
	if (hasPS[thisBranchName])
	{
		fprintf (stdout, "Testing for significance...\n");
		global branchRatio_1;
		branchRatio_1 :< 1; 
		ExecuteCommands ("branchRatio_1 = codon_tree_1." + thisBranchName + ".nonSynRate / codon_tree_1." + thisBranchName + ".synRate");
		ExecuteCommands ("codon_tree_1." + thisBranchName + ".nonSynRate := branchRatio_1 * codon_tree_1." + thisBranchName + ".synRate");

		for (classCount = 2; classCount <= bestBIC_classes; classCount = classCount + 1)
		{
			ExecuteCommands ("global branchRatio_" + classCount + ":< 1;");
			ExecuteCommands ("branchRatio_" + classCount + " = codon_tree_" + classCount+ "." + thisBranchName + ".nonSynRate / codon_tree_" + classCount+ "." + thisBranchName + ".synRate");
			ExecuteCommands ("codon_tree_" + classCount+ "." + thisBranchName + ".nonSynRate := branchRatio_" + classCount +" * codon_tree_" + classCount+ "." + thisBranchName + ".synRate");
			ExecuteCommands ("codon_tree_" + classCount+ "." + thisBranchName + ".synRate    := (codon_tree_1." + thisBranchName + ".synRate*synMultiplier__ + codon_tree_1." + thisBranchName + ".nonSynRate*nsMultiplier__)/(synMultiplier__+nsMultiplier__*branchRatio_" + classCount + ");");
		}
		Optimize (res_null, codonLF);
		nullLogLByBranch [thisBranchName] = res_null[1][0];
		LRT = 2*(alt_res-res_null[1][0]);
		p1  = 1-CChi2(LRT,hasPS[thisBranchName]);
		p2  = 1-CChi2(LRT,bestBIC_classes);
		
		fprintf (summaryOut, ",", res_null[1][0],",", p1 , ",", p2);
		fprintf (stdout, "\n\tLR test statistic = ", LRT, "\n\tp-value (less conservative) = ", p1, "\n\tp-value (more conservative) = ", p2, "\n");
		ExecuteCommands (bestLF);
	}
	else
	{
		fprintf (summaryOut, ",,");
	}
	fprintf (summaryOut, rateString);
}

fprintf (summaryOut, CLOSE_FILE);






/*---------------------------------*/

function computeScaling (sR,nsR)
{
	sc = 0;
	synRate    = sR;
	nonSynRate = nsR;
	GYMN	   = GY_Matrix;
	for (h=0; h<ModelMatrixDimension; h=h+1)
	{
		lsum = 0;
		for (v=0;v<ModelMatrixDimension; v=v+1)
		{
			lsum = lsum + GYMN[h][v] * (v!=h) * codonFreqs[v];
		}
		sc = sc + lsum * codonFreqs[h];
	}
	return sc;
}

	
