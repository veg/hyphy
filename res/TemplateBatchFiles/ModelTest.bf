/* first define global variables for the rate matrix */

	global AC = 1;
	global AT = 1;
	global CG = 1;
	global CT = 1;
	global GT = 1;

/* this is the matrix which does not include rate variation */

	Q_Matrix = {{*,AC*t,t,AT*t}
				{AC*t,*,CG*t,CT*t}
				{t,CG*t,*,GT*t}
				{AT*t,CT*t,GT*t,*}};
			
/* this is the matrix which does not include rate variation */

	Q_Matrix_RV = {{*,AC*t*c,t*c,AT*t*c}
				  {AC*t*c,*,CG*t*c,CT*t*c}
				  {t*c,CG*t*c,*,GT*t*c}
				  {AT*t*c,CT*t*c,GT*t*c,*}};
			
/* a utility matrix used to prepare variable constraints */

	RateBiasTerms = {{"AC","1","AT","CG","CT","GT"}};	

/* model fit cache: an associative array index by the full model name;
   each entry is an associative array with three keys:
   
   	"LogL" 		- the log-likelihood score of the model
   	"DF"   		- the number of model parameters
   	"FitInfo"	- serialized likelihood function for the model
*/

	modelCache = {};

/* 
	model names and corresponding constraint strings;
	the first column are names with equal base frequencies,
	the second column are names with unequal base frequencies
*/

	 modelNameStrings = {{"JC69","F81"},
						 {"K80","HKY85"},
						 {"TrNef","TrN"},
						 {"K81","K81uf"},
						 {"TIMef","TIM"},
						 {"TVMef","TVM"},
						 {"SYM","GTR"}};

	 modelStringsMap 	 = {};
	 
	 modelStringsMap ["000000"] = 0;
	 modelStringsMap ["010010"] = 1;
	 modelStringsMap ["010020"] = 2;
	 modelStringsMap ["012210"] = 3;
	 modelStringsMap ["012230"] = 4;
	 modelStringsMap ["012310"] = 5;
	 modelStringsMap ["012345"] = 6;
	 
	 modelStrings	= Rows (modelStringsMap); /* keys of modelStringsMap */
                
/* include this for directory functions */

	ExecuteAFile  (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Utility"+DIRECTORY_SEPARATOR+"GrabBag.bf");

/**************************************************************************	
Generate a function name from the model string, gamma type and base frequencies
**************************************************************************/

function generateModelDescriptor (modelSpecString, gammaType, freqsType)
{
	modelName = modelNameStrings[modelStringsMap[modelSpecString]][freqsType];
	if (gammaType == 1)
	{
		return modelName + "+G";
	}
	if (gammaType == 2)
	{
		return modelName + "+G+I";
	}
	if (gammaType == 3)
	{
		return modelName + "+I";
	}
	return modelName;
}

 
/**************************************************************************/
/* run a model with the rate matrix given in modelSpecString,
   rate variation options given by gammaType (0: none, 1: gamma, 2: gamma+Inv, 3: Inv)
   freqsVector contains the 4x1 vector of base frequencies to use
   
   returns a model description string
*/
	
/**************************************************************************/
function runAModel (modelSpecString, gammaType, freqsType)
{
	/* first check to see if the model has already been fitted */
	
	modelDescriptionString = generateModelDescriptor (modelSpecString, gammaType, freqsType);
	
	if (Abs (modelCache[modelDescriptionString]))
	/* already done */
	{
		return modelDescriptionString;
	}

	/* 
		clear parameter constraints 1st and then set them according 
	    to the 6-char modelSpecString (as in 010010 for HKY85)
	*/
	
	ClearConstraints (AC,AT,CG,CT,GT);
	
	/* generate a constraint string based on modelSpecString */
	
	CustomModelConstraintString = "";
	CustomModelConstraintString * 128;

	for (customLoopCounter2=1; customLoopCounter2<6; customLoopCounter2=customLoopCounter2+1)
	{
		for (customLoopCounter=0; customLoopCounter<customLoopCounter2; customLoopCounter=customLoopCounter+1)
		{
			if (modelSpecString[customLoopCounter2]==modelSpecString[customLoopCounter])
			{
				if (RateBiasTerms[customLoopCounter2] == "1")
				{
					CustomModelConstraintString * (RateBiasTerms[customLoopCounter]+":="+RateBiasTerms[customLoopCounter2]+";");
				}
				else
				{
					CustomModelConstraintString * (RateBiasTerms[customLoopCounter2]+":="+RateBiasTerms[customLoopCounter]+";");			
				}
				break;
			}
		}
	}	
	CustomModelConstraintString * 0;

	/* define rate variation */
	if (gammaType==1)
	{
		global alpha = .5;
		alpha:>0.01;alpha:<100;
		category c = (classCount, EQUAL, MEAN, 
						GammaDist(_x_,alpha,alpha), 
						CGammaDist(_x_,alpha,alpha), 
						0 , 
				  	    1e25,
				  	    CGammaDist(_x_,alpha+1,alpha)
				  	 );
	}
	else
	{
		if (gammaType==2)
		{
			global alpha   = .5;
			alpha:>0.01;alpha:<100;

			global P;
			P :< 1-(1e-10);
			P = 1/(classCount+1);
			
			global beta = 1;
			beta:>0.01;
			beta:<200;

			beta := (1-P)*alpha;
			
			catFreqMatrix = {classCount+1,1};
			catFreqMatrix [0] := P;
			
			for (h=1; h<=classCount;h=h+1)
			{
				catFreqMatrix [h] := (1-P)/classCount__;
			}
			
			category c = (classCount+1, catFreqMatrix, MEAN, 
					(1-P)*GammaDist(_x_,alpha,beta)*(_x_>0), 
					(1-P)*CGammaDist(_x_,alpha,beta)*(_x_>0)+P, 
					0, 
			  	    1e25,
			  	    (1-P)*CGammaDist(_x_,alpha+1,beta)*(alpha/beta)*(_x_>0)
			  	 );
		}	
		else
		{
			if (gammaType==3)
			{
				global P;
				P :< 1-(1e-10);
				P = 1/2;
				
				catFreqMatrix = {2,1};
				catFreqMatrix [0] := P;
				catFreqMatrix [1] := (1-P);
				
				catRateMatrix = {{0,1/(1-P)}};
				
				category c = (2, catFreqMatrix, MEAN, 
								, 
								catRateMatrix, 
								0, 
						  	    1e25,
						  	 );
			}	
		}
	}

	if (freqsType == 0)
	{
		freqsVector = equal_EFV;
	}
	else
	{
		freqsVector = PI_EFV;
	}
	
	if (gammaType)
	{
		Model currentModel = (Q_Matrix_RV, freqsVector);
	}
	else
	{
		Model currentModel = (Q_Matrix, freqsVector);	
	}
	
	if (Abs(CustomModelConstraintString)) /* have constraints */
	{
		ExecuteCommands (CustomModelConstraintString);
	}
	
	Tree MTTree 		  = treeString;
	LikelihoodFunction lf = (filteredData, MTTree);
	Optimize 				(res,lf);
	
	/* cache model results */
	
	modelCache[modelDescriptionString] = {};
	(modelCache[modelDescriptionString])["LogL"]    = res[1][0];
	(modelCache[modelDescriptionString])["DF"]      = res[1][1] + 3*freqsType; 
							/* add extra parameters for the frequencies if needed */
	
 	/* export a serialized version of the likelihood function */
 	   
	Export (lfInfo, lf);
	(modelCache[modelDescriptionString])["FitInfo"] = lfInfo;
	
	return   modelDescriptionString;
}

/*************************************************************************

	Compare two models, using the chi_1^2 or 0.5(chi_1^2+chi_0^2) test 
	statistic, returning the p-value
	
	model2 is assumed to be the alternative, model1 is assumed to be the null
	
	If reportFlag is set, print the result of model comparison

*************************************************************************/

function compareTwoModels (model1, gamma1, freqs1, model2, gamma2, freqs2, useChiBar, reportFlag)
{	
	modelString1 = runAModel (model1, gamma1, freqs1);
	modelString2 = runAModel (model2, gamma2, freqs2);

	LRT	   = 2*((modelCache[modelString2])["LogL"]-(modelCache[modelString1])["LogL"]);
	DEGF   = (modelCache[modelString2])["DF"]-(modelCache[modelString1])["DF"];
		
	if (useChiBar)
	{
		pV =  .5-.5*CChi2(LRT,DEGF);	
	}
	else
	{
		pV = 1-CChi2(LRT,DEGF);
	}
	
	if (reportFlag)
	{
		reportTest (model1, gamma1, freqs1, model2, gamma2, freqs2, pV);
	}
	return pV;
}

/*************************************************************************
	Report the results of model testing
*************************************************************************/

function reportTest (model1, gamma1, freqs1, model2, gamma2, freqs2, pValue)
{
	modelName1 = generateModelDescriptor (model1, gamma1, freqs1);
	modelName2 = generateModelDescriptor (model2, gamma2, freqs2);

	fprintf (stdout, "\n\tNull:", modelName1, "\n\t\tLog-likelihood = ", Format((modelCache[modelString1])["LogL"],12,4),
														", ", (modelCache[modelString1])["DF"], " parameters."); 
														
	fprintf (stdout, "\n\tAlt :", modelName2, "\n\t\tLog-likelihood = ", Format((modelCache[modelString2])["LogL"],12,4), 
														", ", (modelCache[modelString2])["DF"], " parameters."); 
														
	fprintf (stdout, "\n\n\tLR statistic       : ", 2*((modelCache[modelString2])["LogL"]-(modelCache[modelString1])["LogL"]),
					 "\n\tDegrees of freedom : ",(modelCache[modelString2])["DF"]-(modelCache[modelString1])["DF"],
					 "\n\tp-Value            : ", pValue);
	
	if (pValue<rejectAt)
	{
		fprintf (stdout, "\n\n\tNull hypothesis rejected\n\t", modelName2, " chosen.\n");
		return modelName2;
	}
	else
	{
		fprintf (stdout, "\n\n\tNull hypothesis accepted\n\t", modelName1, " chosen.\n");
		return modelName1;
	}	
	return 1;
}     

/*************************************************************************
print model information
*************************************************************************/

function printModelDescription (modelString, gammaRate, modelFreq)
{
	sep =            "+---+-------------+-------------+-------------+-------------+\n";
	fprintf (stdout, "Rate matrix\n", 
					 sep,
					 "|   |      A      |      C      |      G      |      T      |\n",
					 sep,
					 "| A |      *      | ",Format(AC,11,5)," | ", Format(1.0,11,5), " | ", Format(AT,11,5), " |\n",
					 sep,
					 "| C | ",Format(AC,11,5)," |      *      | ", Format(CG,11,5), " | ", Format(CT,11,5), " |\n",
					 sep,
					 "| G | ",Format(1,11,5), " | ",  Format(CG,11,5), " |      *      | ", Format(GT,11,5), " |\n",
					 sep,
					 "| T | ",Format(AT,11,5), " | ",  Format(CT,11,5), " | ",  Format(GT,11,5)," |      *      |\n",
					 sep);
						 
	fprintf (stdout, "\nBase frequencies",
					 "\nA = ", Format (freqsVector[0],10,4),
					 "\nC = ", Format (freqsVector[1],10,4),
					 "\nG = ", Format (freqsVector[2],10,4),
					 "\nT = ", Format (freqsVector[3],10,4), "\n");
					 
	if (gammaRate)
	{
		fprintf (stdout, "\nRate variation");
		GetInformation 		(catInfo, c);
		for (idx = 0; idx < Columns (catInfo); idx = idx + 1)
		{
			fprintf (stdout, "\nRate ",Format (idx+1,5,0), " = ", Format (catInfo[0][idx], 10,5), " (weight = ", Format (catInfo[1][idx], 10,5), ")");
		}
		fprintf (stdout, "\n");
	}
	
	return 1;
}

/*************************************************************************
echo instructions of how to choose the model using the custom model 
*************************************************************************/


function echoModelSpec (modelString, modelRate, modelFreq)
{
	fprintf (stdout, "\n\nTo define this model (if it isn't predefined) for analyses,\nUse 'Custom' model with the following options:\n",
					 "Model String:", modelString, "\n");
					 
	fprintf (stdout, "Model Options: ");
	if (modelRate==0)
	{
		fprintf (stdout, "Global");
	}
	else
	{
		if (modelRate==1)
		{
			fprintf (stdout, "Rate Variation, then choose the Gamma distribution");
		}
		else
		{
			if (modelRate==2)
			{
				fprintf (stdout, "Rate Variation, then choose the Gamma+Inv distribution");
			}
			else
			{
				if (modelRate==3)
				{
					fprintf (stdout, "Rate Variation, then choose the Inv Class distribution");
				}
			}
		}
	}
	fprintf (stdout, "\nEquilibrium Frequencies Option: ");
	if (modelFreq)
	{
		fprintf (stdout, "Observed\n\n");	
	}
	else
	{
		fprintf (stdout, "Equal\n\n");	
	}
	return 1;
}   

                    

/**************************************************************************/

function doHierarchicalTest (echoFlag)
{
	if (echoFlag)
	{
		fprintf (stdout, "\n\n****** RUNNING HIERARCHICAL MODEL TESTING ******\n\n");
	}

	chosenMatrix = 0;
	chosenRates  = 0;
	chosenFreqs  = 0;
	chosen2		 = 0;
	chosen3      = 0;
	chosen4		 = 0;
	
	/* do the frequencies test */
	
	if (echoFlag)
	{
		fprintf (stdout,"\n1). Checking for equilibrium frequencies equality.\n");
	}

	pVal 		=  compareTwoModels (modelStrings[0],0,0,modelStrings[0],0,1,0, echoFlag);
	chosenFreqs = (pVal<rejectAt);
	
	/* do the tr=tv test */
	
	if (echoFlag)
	{
		fprintf (stdout,"\n2). Checking for equality of transition and transversion rates.\n");
	}
	
	pVal 		=  compareTwoModels (modelStrings[0],0,chosenFreqs,modelStrings[1],0,chosenFreqs,0, echoFlag);
	
	chosen2 = (pVal<rejectAt);
		
	if (chosen2)
	/* finetune the matrix */
	{
		if (echoFlag)
		{
			fprintf (stdout,"\n3). Checking for equality of A<->G and C<->T (transition) rates.\n");
		}
		pVal =  compareTwoModels (modelStrings[1],0,chosenFreqs,modelStrings[2],0,chosenFreqs,0, echoFlag);
		chosen3 = (pVal<rejectAt);
				
		if (echoFlag)
		{
			fprintf (stdout,"\n4). Checking whether there are 1 or 2 transversion rates.\n");
		}
		
		pVal =  compareTwoModels (modelStrings[1+chosen3],0,chosenFreqs,modelStrings[3+chosen3],0,chosenFreqs,0, echoFlag);
		
		chosen4      = (pVal<rejectAt);
		chosenMatrix = 1+chosen3+2*chosen4;
		if (chosen4)
		/* fine tune some more */
		{
			if (echoFlag)
			{
				fprintf (stdout,"\n5). Checking whether there are 2 or 4 transversion rates.\n");
			}
			pVal =  compareTwoModels (modelStrings[3+chosen3],0,chosenFreqs,modelStrings[5+chosen3],0,chosenFreqs,0, echoFlag);
			chosenMatrix = 3+chosen3+2*(pVal<rejectAt);
		}
		else
		{
			if (echoFlag)
			{
				fprintf (stdout,"\n...Skipping step 5...\n");
			}
		}
	}
	else
	{
		if (echoFlag)
		{
			fprintf (stdout,"\n...Skipping steps 3 through 6...\n");
		}
		chosenMatrix = 0;
	}
	
	if (echoFlag)
	{
		fprintf (stdout,"\n6). Checking for evidence of rate variation.\n");
	}
	pVal =  compareTwoModels (modelStrings[chosenMatrix],0,chosenFreqs,modelStrings[chosenMatrix],1,chosenFreqs,1, echoFlag);
	
	if (echoFlag)
	{
			fprintf (stdout,"\n7). Checking for evidence of an invariant rate class.\n");
	}
	if (pVal<rejectAt)
	/* accept rate variation */
	{
		pVal =  compareTwoModels (modelStrings[chosenMatrix],1,chosenFreqs,modelStrings[chosenMatrix],2,chosenFreqs,1, echoFlag);
		if (pVal<rejectAt)
		{
			chosenRates = 2;
		}
		else
		{
			chosenRates = 1;
		}
	}
	else
	{
		pVal =  compareTwoModels (modelStrings[chosenMatrix],0,chosenFreqs,modelStrings[chosenMatrix],3,chosenFreqs,1, echoFlag);
		if (pVal<rejectAt)
		{
			chosenRates = 3;
		}
		else
		{
			chosenRates = 0;
		}	
	}
	return 0;
}

/**************************************************************************/

function doAICTest (echoFlag)
{
	if (echoFlag)
	{
		fprintf (stdout, "\n\n****** RUNNING AIC BASED MODEL SELECTION ******\n\n");
	}
	aicScore = 1e100;
	
	if (echoFlag)
	{	
		fprintf (stdout,"\n|  Model     | # prm |    lnL    |    AIC     |\n",
				          "|------------|-------|-----------|------------|");   
	}		
	for (k=0; k<Columns(modelStrings); k=k+1)
	{
		for (i=0; i<2; i=i+1) /* freqOption */
		{
			for (j=0; j<4; j=j+1) /* gammaOption */
			{
				thisModelName = runAModel (modelStrings[k],j,i);
				
				modelLL	 = (modelCache[thisModelName])["LogL"];
				modelDF  = (modelCache[thisModelName])["DF"];
				
				modelAIC = 2*(modelDF-modelLL);
								
				if (echoFlag)
				{
					while (Abs(thisModelName)<11)
					{
						thisModelName = thisModelName+" ";						
					}
					fprintf (stdout, "\n| ", thisModelName, "| ", Format (modelDF,5,0),
										   " | ", Format (modelLL,9,3), " | ", Format (modelAIC,10,3), " |");
				}
				if (modelAIC < aicScore)
				/* best model so far */
				{
					aicFreqs   = i;
					aicRates   = j;
					aicMatrix  = k;
					aicScore = modelAIC;
					if (echoFlag)
					{
						fprintf (stdout, " *");
					}
				}
			}
		}
	}
	if (echoFlag)
	{
		fprintf (stdout,"\n|------------|-------|-----------|------------|");   
	}
	
	return 0;
}

/**************************************************************************/
/* MAIN BODY 															  */
/**************************************************************************/

OPTIMIZE_SUMMATION_ORDER = 1;

fprintf (stdout, "\n+--------------------------------+\n",
				 "| RUNNING MODEL TESTING ANALYSIS |\n",
				 "| Based on the program ModelTest |\n",
				 "|              by                |\n",
				 "|        David  Posada           |\n",
				 "|             and                |\n",
				 "|        Keith Crandall          |\n",
				 "|                                |\n",
				 "|    If you use this analysis,   |\n",
				 "| be sure to cite the original   |\n",
				 "| reference, which can be found  |\n",
				 "| in Bioinformatics (1998)       |\n",
				 "| Vol 14, ppg. 817-818           |\n",
				 "+--------------------------------+\n");

fprintf			  (stdout, "\nTesting space includes ", Abs(modelStringsMap)*8, " models\n");

SetDialogPrompt   ("Please specify a nucleotide data file:");
DataSet ds    	= ReadDataFile (PROMPT_FOR_FILE);

dataFilePath	= LAST_FILE_PATH;
DataSetFilter 	  filteredData = CreateFilter (ds,1);
fprintf (stdout, "\n\nData read from \"",dataFilePath,"\"\n",ds,"\n");

filePathComponents = splitFilePath (dataFilePath);

HarvestFrequencies (PI_EFV,filteredData,1,1,1);
equal_EFV = {{.25}{.25}{.25}{.25}};

#include "queryTree.bf";

classCount = 0;
while (classCount<=0)
{
	fprintf (stdout, "\nNumber of rate classes in rate variation models (e.g. 4):");
	fscanf  (stdin,"Number", classCount);
	classCount = classCount $ 1; /* truncate the number */
}

KEEP_OPTIMAL_ORDER = 1;
MESSAGE_LOGGING    = 0;
USE_LAST_RESULTS   = 1;

ChoiceList (runType,"Model Selection Method",1,SKIP_NONE,
			"Hierarchical Test","Perform a series of nested model comparisons to select the model.",
			"AIC Test","Obtain MLEs for each model, and select the best one using Akaike Information Criterion.",
			"Both","Run the AIC test, followed by the Hierarchical test");
			
if (runType<0)
{
	return;
}

modelNum		   = 0;
rejectCount 	   = 0;

if (runType != 1) /* need to run hierarchical testing */
{
	rejectAt = 0;

	while ((rejectAt<=0)||(rejectAt>=1))
	{
		fprintf (stdout, "\nModel rejection level (e.g. 0.05):");
		fscanf  (stdin,"Number", rejectAt);
	}
}

if (runType)
{
	doAICTest (1);
	aicMdl = generateModelDescriptor (modelStrings[aicMatrix],aicRates,aicFreqs);
	fprintf (stdout, "\n\nAIC-based model: ",aicMdl,", AIC = ", aicScore);
}

if (runType!=1)
{
	hMdl = doHierarchicalTest (1);
}



OPTIMIZE_SUMMATION_ORDER = 0;

if (runType!=1)
{
	hMdl = generateModelDescriptor (modelStrings[chosenMatrix],chosenRates,chosenFreqs);
	ExecuteCommands ((modelCache[hMdl])["FitInfo"]);
	DEFAULT_FILE_SAVE_NAME = filePathComponents["FILENAME"] + "_hierarchicalFit.bf";
	SetDialogPrompt ("Save hierarchical fit to:");
	fprintf (PROMPT_FOR_FILE,CLEAR_FILE,(modelCache[hMdl])["FitInfo"]);
	fprintf (stdout, "\n\n ******* Hierarchical Testing based model (",hMdl,") ******** \n\n");
	printModelDescription (modelStrings[chosenMatrix],chosenRates,chosenFreqs);
	echoModelSpec	 	  (modelStrings[chosenMatrix], chosenRates, chosenFreqs);
}


if (runType)
{
	if (runType == 2)
	{
		if (aicMdl != hMdl)
		{
			ExecuteCommands ((modelCache[aicMdl])["FitInfo"]);
			DEFAULT_FILE_SAVE_NAME = filePathComponents["FILENAME"] + "_aicFit.bf";
			SetDialogPrompt ("Save AIC-based fit to:");
			fprintf (PROMPT_FOR_FILE,CLEAR_FILE,(modelCache[aicMdl])["FitInfo"]);
		}
		else
		{
			fprintf (stdout, "\nAIC based testing concurs with hierarchical testing\n");
			OPTIMIZE_SUMMATION_ORDER = 1;
			return 0;
		}
	}
	else
	{
		ExecuteCommands ((modelCache[aicMdl])["FitInfo"]);
		DEFAULT_FILE_SAVE_NAME = filePathComponents["FILENAME"] + "_aicFit.bf";
		SetDialogPrompt ("Save AIC-based fit to:");
		fprintf (PROMPT_FOR_FILE,CLEAR_FILE,(modelCache[aicMdl])["FitInfo"]);
	}
	fprintf (stdout, "\n\n ******* AIC based model (",aicMdl,") ******** \n\n");
	printModelDescription (modelStrings[aicMatrix],aicRates,aicFreqs);
	echoModelSpec	 	  (modelStrings[aicMatrix],aicRates,aicFreqs);
}

OPTIMIZE_SUMMATION_ORDER = 1;
