ModelMatrixDimension = 0;

baselineOutput   = "";

/*---------------------------------------------------------------------------------------------------------------------------------------------------*/

function computeModelBranchLengthExpression (freqVector)
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
												   freqVector[h-hshift]*freqVector[transition][nucPosInCodon]+
												   freqVector[v-vshift]*freqVector[transition2][nucPosInCodon];
					rateTypeWeightsGY94[mxIndex][t1] = rateTypeWeightsGY94[mxIndex][t1]+
												   2*freqVector[h-hshift]*freqVector[v-vshift];
				}
		   }
	    }		
	}
	return 0;		
}

/*---------------------------------------------------------------------------------------------------------------------------------------------------*/

function echoCatVar (distrInfo)
{
	DD = Columns(distrInfo);
	EE = 0.0;
	sampleVar = 0.0;
	for (k=0; k<DD; k=k+1)
	{
		EE = distrInfo[0][k]*distrInfo[1][k]+EE;
		sampleVar = sampleVar+distrInfo[0][k]*distrInfo[0][k]*distrInfo[1][k];
	}
		
	sampleVar = sampleVar-EE*EE;
	
	fprintf  (distribOutput,"\n\n------------------------------------------------\n\nSample mean = ",EE, " (sample variance = ",sampleVar,")\n");
	for (k=0; k<DD; k=k+1)
	{
		fprintf (distribOutput,"\nRate[",Format(k,0,0),"]=",Format(distrInfo[0][k],12,8), " (weight=", 
						  Format(distrInfo[1][k],9,7),")");
	}
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


SetDialogPrompt ("Please specify a codon data file:");

DataSet ds = ReadDataFile (PROMPT_FOR_FILE);
dataFilePath = LAST_FILE_PATH;
fprintf (stdout,"\n______________READ THE FOLLOWING DATA______________\n",ds);

#include "TemplateModels/chooseGeneticCode.def";

DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);

HarvestFrequencies (observedFreq,	   filteredData,3,1,1);

#include "queryTree.bf";

ChoiceList (branchLengths,"Branch Lengths",1,SKIP_NONE,
			"Codon Model","Jointly optimize rate parameters and branch lengths (slow and thorough)",
			"Nucleotide Model","Estimate branch lengths once, using an appropriate nucleotide model (quick and dirty)."
		    );

if (branchLengths<0)
{
	return;
}

modelConstraintString = "";

ChoiceList (modelChoice, "Rate Matrix Options",1,SKIP_NONE,
			"MG94",	 "Standard Muse-Gaut 94 model.",
			"MG94xHKY85","MG94 with the transition/transverion ratio parameter kappa.",
			"MG94xREV","MG94 with 5 additional parameters for each type of nucleotide substitution ratio. (Recommended Rate Matrix)",
			"MG94xCustom","MG94 crossed with an arbitrary nucelotide reversible model, except F81.",
			/*"MG94xCustomxAA","MG94 crossed with an arbitrary nucelotide reversible model and an arbitrary amino-acid bias profile.",*/
			"GY94","Goldman-Yang, 94 model (similar to MG94xHKY85)"
);

if (modelChoice<0)
{
	return;
}

marginalOutFile = "2RatesAnalyses/MG94xREV.mdl";
if (modelChoice == 0)
{
	modelConstraintString = "AC:=1;AT:=1;CG:=1;CT:=1;GT:=1"; 
	ModelTitle = "MG94";
}
else
{
	if (modelChoice == 1)
	{
		modelConstraintString = "CT:=1;AT:=AC;CG:=AC;GT:=AC"; 
		ModelTitle = "MG94xHKY85";
	}
	else
	{
		if (modelChoice == 3)
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
			ModelTitle = "MG94x"+modelDesc[0];
						
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
			if (modelChoice == 4)
			{
				marginalOutFile = "2RatesAnalyses/GY94.mdl";
				ModelTitle = "GY94";
			}
			else
			{
				ModelTitle = "MG94xREV";
			}
		}			
	}
}

chosenModelList = {5,1};

ChoiceList (modelChoice,"Rate Variation Options",1,SKIP_NONE,
			"Run All","Run all available models.",
			"Run Custom","Choose some of the available models.");
			
if (modelChoice<0)
{
	return;
}

if (modelChoice==0)
{
	for (mi = 0; mi<5; mi=mi+1)
	{
		chosenModelList[mi][0] = 1;
	}
}
else
{
	ChoiceList (modelTypes,"Rate Variation Models",0,SKIP_NONE,
				"Constant","Constant Rate Model: no rate variation across sites.", /* index 0 */
				"Proportional","Proportional Variable Rates Model: dS and dN vary along the sequence, but dN = R*dS for every site",
				"Nonsynonymous","Non-synonymous Variable Rates Model: dS = 1 for every site, while dN is drawn from a given distribution.",
				"Dual","Dual Variable Rates Model: dS and dN are drawn from a bivariate distribution (independent or correlated components). Recommened model.",
				"Lineage Dual","Lineage Dual Variable Rates Model:  dS and dN are drawn from a bivariate distribution (independent or correlated components), plus each lineage has an adjustment factor for the E[dN]/E[dS]."
			    );
			    
	if (modelTypes[0]<0)
	{
		return;
	}
	
	for (mi = 0; mi < Rows(modelTypes)*Columns(modelTypes); mi = mi + 1)
	{
		modelChoice = modelTypes[mi];
		chosenModelList[modelChoice] = 1;
	}	
}

modelNamesShort = {{"Constant","Proportional","Nonsynonymous","Dual","LineageDual"}};


ChoiceList (modelChoice, "Distribution Option",1,SKIP_NONE,
			"Syn:Gamma, Non-syn:Gamma",	 "Both syn and non-syn rates are drawn from the gamma distributions for all models.",
			"Syn:Gamma, Non-syn:Inv+Gamma","Syn and non-syn rates are drawn from the gamma distributions for Proportional and Nonsynonymous. For Dual and Local Dual, syn rates are drawn from the gamma distribution, and non-syn rates - from Inv+Gamma.",
			"Independent Discrete", "Independent General Discrete Distributions (Recommended setting)",
			"Correlated Discrete", "Correlated General Discrete Distributions");
			
			
if (modelChoice < 0)
{
	return;
}

ChoiceList (randomizeInitValues, "Initial Value Options",1,SKIP_NONE,
			"Default",	 "Use default inital values for rate distribution parameters.",
			"Randomized",	 "Select initial values for rate distribution parameters at random.");


if (randomizeInitValues < 0)
{
	return;
}

resp  = 0;
resp2 = 0;

while (resp<2)
{
	fprintf (stdout,"Number of synonymous (and single variable rate modles) rate classes (>=2):");
	fscanf  (stdin, "Number", resp);
}

while (resp2<2)
{
	fprintf (stdout,"Number of non-synonymous rate classes (>=2):");
	fscanf  (stdin, "Number", resp2);
}
			
fudgeFactor = 1.0;

if (modelChoice<2)
{
	fscanf ("2RatesAnalyses/gamma1.def","Raw",categDef1);

	if (modelChoice == 0)
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
	correlationOn = (modelChoice>2);
	fscanf ("2RatesAnalyses/discreteGenerator.bf","Raw",mi);
	ExecuteCommands (mi);
}

ExecuteCommands ("#include \""+marginalOutFile+"\";");

if (Abs(modelConstraintString))
{
	ExecuteCommands (modelConstraintString);
}


SetDialogPrompt ("Save summary result file to:");
fprintf (PROMPT_FOR_FILE,CLEAR_FILE);
baselineOutput = LAST_FILE_PATH;

distribOutput = baselineOutput + ".distributions";
fprintf (distribOutput,CLEAR_FILE);



if (branchLengths)
{
	branchLengths = -branchLengths;
}

separator =      "+---------------------+---------------------+---------------+-------------------+-------------------+---------------+-----+-----------+\n";

fprintf (stdout, "\n\nRUNNING ",ModelTitle," MODEL COMPARISONS on ",dataFilePath, "\n\n",
				 "\n\n########### ",resp,"x",resp2," CLASSES ###########\n\n",
				 separator,
				 "|       Model         |   Log Likelihood    | Synonymous CV |  NS Exp and CV    |  N/S Exp and CV   |    P-Value    | Prm |    AIC    |\n",
				 separator);

fprintf (baselineOutput,  "\n\nRUNNING ",ModelTitle," MODEL COMPARISONS on ",dataFilePath, "\n\n",
						   "\n\n########### ",resp,"x",resp2," CLASSES ###########\n\n",
				 separator,
				 "|       Model         |   Log Likelihood    | Synonymous CV |  NS Exp and CV    |  N/S Exp and CV   |    P-Value    | Prm |    AIC    |\n",
				 separator);

lastLikValue = 0;
R = 1;

modelNames = {{"| Constant Rates      |",
			   "| Prop. Var. Rates    |",
			   "| Var. N.Syn. Rates   |",
			   "| Dual Variable Rates |",
			   "| Dual V.R. + Lineage |"}};
			   


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

for (mi = 0; mi<5; mi=mi+1)
{
	if (chosenModelList[mi])
	{
		theRateMatrix = 0;

		MULTIPLY_BY_FREQS = PopulateModelMatrix ("theRateMatrix", observedFreq, mi);
		vectorOfFrequencies = BuildCodonFrequencies (observedFreq);
		
		if (doNucFit)
		{
			HarvestFrequencies 					   (observedFreqSingle,filteredData,1,1,1);
			DataSetFilter nucFilter = CreateFilter (filteredData,1);
			ExecuteCommands (nucModelString+"\nModel nucModel = (nucModelMatrix,observedFreqSingle);");

			Tree  nucTree = treeString;
			LikelihoodFunction nuc_lf = (nucFilter,nucTree);
			Optimize (nuc_res, nuc_lf);
			
			computeExpSubWeights (0);		
			global codonFactor = 0.33;
						
			doNucFit = 0;
		}
			
		Model MG94model = (theRateMatrix,vectorOfFrequencies,MULTIPLY_BY_FREQS);
		
		Tree  givenTree = treeString;
		if (branchLengths)
		{
			ClearConstraints (givenTree);
			ReplicateConstraint ("this1.?.synRate:=this2.?.t__/codonFactor",givenTree,nucTree);
		}
		else
		{
			lbc = TipCount (givenTree);
			initString = "";
			for (lc = 0; lc < lbc; lc = lc+1)
			{
				bn = TipName (givenTree,lc);
				initString = initString + "givenTree." + bn + ".synRate=nucTree." + bn + ".t/codonFactor;";
			}	
			lbc = BranchCount (givenTree);
			for (lc = 0; lc < lbc; lc = lc+1)
			{
				bn = BranchName (givenTree,lc);
				initString = initString + "givenTree." + bn + ".synRate=nucTree." + bn + ".t/codonFactor;";
			}	
			ExecuteCommands (initString);
			initString = "";
			
		}

		LikelihoodFunction lf = (filteredData,givenTree);
		
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
				MPINodeState[mpiNode][1] = mi;
			}
		}
		else
		{
			/* Non-MPI execution */
			Optimize (res,lf);
			modelIndex = mi;
			dummy = ReceiveJobs (0);
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
		MPIReceive (-1, fromNode, result_String);
		modelIndex = MPINodeState[fromNode-1][1];
		
		if (sendOrNot)
		{
			/* send the likelihood function to the node which just finished */
			MPISend (fromNode,lf);
			MPINodeState[fromNode-1][1] = mi;
		}
		else
		{
			/* mark the node as idle */
			MPINodeState[fromNode-1][0] = 0;
			MPINodeState[fromNode-1][1] = 0;
		}
		
		/* reset category variables */
		if (modelIndex)
		{
			ClearConstraints (c);
		}
		if (modelIndex>2)
		{
			ClearConstraints (d);
		}
		
		/* build the likelihood function lf with MLE parameter values
		   just returned from an MPI node. MLE matrix is returned 
		   in lf_MLES. */
		   		
		ExecuteCommands (result_String);
				
		/* store the result matrix */
		res = lf_MLES;
	}
	
	
	if (branchLengths)
	{
		res[1][1] = res[1][1] + TipCount (givenTree) + BranchCount (givenTree) - 1;
	}

	LIKELIHOOD_FUNCTION_OUTPUT = 6;

	fprintf (stdout, modelNames[modelIndex], " ", Format (res[1][0],19,5));
	fprintf (baselineOutput, modelNames[modelIndex], " ", Format (res[1][0],19,5));
	
	marginalOutFile = baselineOutput + "." + modelNamesShort[modelIndex] + ".fit";
	
	fprintf (marginalOutFile,CLEAR_FILE,lf);

	if (modelIndex)
	{
		marginalOutFile = baselineOutput + "." + modelNamesShort[modelIndex] + ".marginals";
		fprintf (marginalOutFile,CLEAR_FILE);
	}
	
	fprintf (distribOutput, "\n", modelNames[modelIndex]);
	if (modelIndex==0)
	{
		EE = 0;
		sampleVar = 0;
		fprintf (stdout, " |      N/A      | ", Format (R,8,5) , ",", Format (sampleVar,8,5), " | ", Format (R,8,5) , ",", Format (sampleVar,8,5), " | ");
		fprintf (baselineOutput, " |      N/A      | ", Format (R,8,5) , ",", Format (sampleVar,8,5), " | ", Format (R,8,5) , ",", Format (sampleVar,8,5), " | ");
	}
	else
	{
		GetInformation(dI,c);
		if (modelIndex==2)
		{
			DD = Columns (dI);
			for (EE=0; EE<DD; EE=EE+1)
			{
				dI[0][EE] = R*dI[0][EE];
			}
			
			NVRMLL = res[1][0];
			NVRMPC = res[1][1];
			
			EE  = echoCatVar (dI);
			fprintf (marginalOutFile, dI);
			
			fprintf (stdout, " |      N/A      | ", Format (EE,8,5) , ",", Format (Sqrt(sampleVar)/EE,8,5), " | ", Format (EE,8,5) , ",", Format (Sqrt(sampleVar)/EE,8,5), " | ");
			fprintf (baselineOutput, " |      N/A      | ", Format (EE,8,5) , ",", Format (Sqrt(sampleVar)/EE,8,5), " | ", Format (EE,8,5) , ",", Format (Sqrt(sampleVar)/EE,8,5), " | ");
		}  
		else
		{
			fprintf (marginalOutFile, dI);
			EE  = echoCatVar (dI);
			fprintf (stdout, " | ",Format (Sqrt(sampleVar),13,8)," | ");
			fprintf (baselineOutput, " | ",Format (Sqrt(sampleVar),13,8)," | ");
			
			if (modelIndex>=3)
			{
				GetInformation(dI2,d);
				if (modelIndex!=5)
				{
					DD = Columns (dI2);
					for (EE=0; EE<DD; EE=EE+1)
					{
						dI2[0][EE] = R*dI2[0][EE];
					}
				}
				fprintf (marginalOutFile, dI2);
				EEN  = echoCatVar (dI2);
				varN = sampleVar;
				EER  = echoRatio  (dI2,dI);
				varR = sampleVar;
				
				fprintf (stdout,  Format (EEN,8,5) , ",", Format (Sqrt(varN)/EEN,8,5), " | ",Format (EER,8,5) , ",", Format (Sqrt(varR)/EER,8,5), " | ");	
				fprintf (baselineOutput, Format (EEN,8,5) , ",", Format (Sqrt(varN)/EEN,8,5), " | ",Format (EER,8,5) , ",", Format (Sqrt(varR)/EER,8,5), " | ");	
					
				if (modelIndex==3)
				{
					DVRMLL = res[1][0];
					DVRMPC = res[1][1];
					/* we don't know that DVRM returned before M4, so we won't print the p-value in the MPI mode */
					if ((chosenModelList[2])&&(MPI_NODE_COUNT<=1))
					{
						pVal = 1-CChi2(2*(DVRMLL-NVRMLL),DVRMPC-NVRMPC);
						fprintf (stdout, Format (pVal,13,8), " |");
						fprintf (baselineOutput, Format (pVal,13,8), " |");
					}
					else
					{
						fprintf (stdout,"     N/A      |");
						fprintf (baselineOutput, "     N/A      |");
					}
				}
			}
			else
			{
				fprintf (stdout, Format (EE*R,8,5) , ",", Format (Sqrt(sampleVar)/EE,8,5), " | ", Format (R,8,5) , ",", Format (0,8,5), " | ");
				fprintf (baselineOutput, Format (EE*R,8,5) , ",", Format (Sqrt(sampleVar)/EE,8,5), " | ", Format (R,8,5) , ",", Format (0,8,5), " | ");
			}
		}
	}
		
	if (modelIndex==4) /* local rates */
	{
		/* we don't know that DVRM returned before M4, so we won't print the p-value in the MPI mode */
		if ((chosenModelList[3])&&(MPI_NODE_COUNT<=1))
		{
			fprintf (stdout, Format (1-CChi2(2*(res[1][0]-DVRMLL),res[1][1]-DVRMPC),13,8), " |");
			fprintf (baselineOutput, Format (1-CChi2(2*(res[1][0]-DVRMLL),res[1][1]-DVRMPC),13,8), " |");
		}
		else
		{
			fprintf (stdout,"     N/A      |");
			fprintf (baselineOutput, "     N/A      |");
		}
	}
	else
	{
		if (modelIndex < 3)
		{
			fprintf (stdout,"     N/A      |");
			fprintf (baselineOutput,"     N/A      |");
		}
	}
	
	fprintf (stdout, Format (res[1][1],5,0), "|", Format (2*(res[1][1]-res[1][0]),11,2),"|");
	fprintf (baselineOutput, Format (res[1][1],5,0), "|", Format (2*(res[1][1]-res[1][0]),11,2),"|");
	fprintf (stdout, "\n",separator);
	fprintf (baselineOutput, "\n",separator);
	if (modelIndex)
	{
		ConstructCategoryMatrix(marginals,lf,COMPLETE);
		
		if (modelIndex>=3)
		{
			GetInformation (categVarIDs,lf);
			if (categVarIDs[0]!="c")
			{
				marginalsCorrected = marginals;
				fprintf (MESSAGE_LOG,"Adjusting marginal matrix rows.\n"); 
				for (h=0; h<Columns(marginals); h=h+1)
				{
					transition = 0;
					for (diff=0; diff<resp; diff = diff+1)
					{
						for (v=diff; v<Rows(marginals); v=v+resp)
						{
							marginalsCorrected[transition][h] = marginals[v][h];
							transition = transition+1;
						}
					}
				}
				marginals = marginalsCorrected;
			}
		}

		fprintf (marginalOutFile,marginals);
	}
	lastLikValue = res[1][0];	
	
	return fromNode-1;
}
