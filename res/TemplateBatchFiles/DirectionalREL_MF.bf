BayesFactorCutoff = 20;

AAString    = "ACDEFGHIKLMNPQRSTVWY";
AACharToIdx = {};
for (k=0; k<20; k=k+1)
{
	AACharToIdx [AAString[k]] = k;
}

SKIP_MODEL_PARAMETER_LIST = 0;

#include "AddABias.ibf";

test_p_values = {20,2};

/*--------------------------------------------------------------------------------------------*/

function GetEqFreqs (ModelMatrixName&, baseFreqs)
{
	t = 1;
	numRateMx = ModelMatrixName;
	for (ri = 0; ri < 20; ri = ri+1)
	{
		for (ci = 0; ci < 20; ci = ci+1)
		{
			if (ri != ci)
			{
				numRateMx[ri][ci] = numRateMx[ri][ci] * baseFreqs[ci];
				numRateMx[ri][ri] = numRateMx[ri][ri] - numRateMx[ri][ci];
			}
		}
	}
	for (ri = 0; ri < 20; ri = ri+1)
	{
		numRateMx [ri][19] = 1;
	}
	numRateMxI = Inverse (numRateMx);
	return numRateMxI [19][-1];
}

/*--------------------------------------------------------------------------------------------*/


ChoiceList						(reloadFlag, "Reload/New", 1, SKIP_NONE, "New analysis","Start a new analysis",
																	      "Reload","Reload a baseline protein fit.");

ACCEPT_ROOTED_TREES 			= 1;

if (reloadFlag == 0)
{

	SetDialogPrompt 			("File to examine:");
	_doAAcids					= 1;
	ExecuteAFile    			("_MFReader_.ibf");
	ExecuteAFile				("TreeTools.ibf");
	basePath 					 = LAST_FILE_PATH;

	UseModel	    	(USE_NO_MODEL);
	roots	= {};
	for (k=0; k<fileCount; k=k+1)
	{
		Tree		  T = treeStrings[k+1];
		roots[k]	    = selectATreeBranch ("T", "Choose a root for partition :" + (k+1));
		
		treeStrings [k+1] = RerootTree (T, roots[k]);
	}
	
	SKIP_FREQ_HARVESTING = 1;
	vectorOfFrequencies = overallFrequencies;
	promptModel (0);
	SKIP_FREQ_HARVESTING = 0;
	populateTrees ("givenTree", fileCount);	
	
	ChoiceList						(fixFB, "Fix Branch", 1, SKIP_NONE, "Unknown root","The character at the root of the tree is drawn from the stationary distribution",
																		"Fix rooting sequence as root","The rooting sequence in the file is used to populate the root sequences.");
	if (fixFB>0)
	{
		for (k=0; k<fileCount; k=k+1)
		{
			ExecuteCommands ("givenTree_"+(k+1)+"."+roots[k]+".t:=0");
		}
	}
	else
	{
		if (fixFB < 0)
		{
			return 0;
		}
	}
	
	ExecuteCommands (constructLF ("lf", "nucData", "givenTree", fileCount));
	
	fprintf							(stdout, "[PHASE 0.1] Standard model fit\n"); 
	
	VERBOSITY_LEVEL				= 1;
	AUTO_PARALLELIZE_OPTIMIZE	= 1;
	Optimize 						(res0,lf);
	AUTO_PARALLELIZE_OPTIMIZE	= 0;
	VERBOSITY_LEVEL				= -1;
	
	
	LIKELIHOOD_FUNCTION_OUTPUT = 7;
	outPath = basePath + ".base";
	fprintf (outPath, CLEAR_FILE, lf);
	
	baselineLogL				= res0[1][0];
	
}
else
{
	SetDialogPrompt ("Locate an existing fit:");
	ExecuteAFile 				(PROMPT_FOR_FILE);
	basePath 					 = LAST_FILE_PATH;
	GetString					 (lfInfo, lf, -1);
	fileCount 					= Columns(lfInfo["Trees"]);
	ExecuteCommands 			 ("GetString	(modelInfo,"+(lfInfo["Models"])[0]+",-2);");
	modelNameString = modelInfo["Matrix"];
	treeStrings = {}; 
	totalCodonCount = 0;
	for (k=0; k<fileCount; k=k+1)
	{
		ExecuteCommands ("treeStrings [k+1] = Format(givenTree_"+(k+1)+",1,1);");
		ExecuteCommands ("totalCodonCount = totalCodonCount + nucData_" + (k+1) + ".sites");
	}
	treeString = Format (givenTree,0,0);
	LFCompute (lf,LF_START_COMPUTE);
	LFCompute (lf,baselineLogL);
	LFCompute (lf,LF_DONE_COMPUTE);
	_runAsFunctionLibrary = 1;
	ExecuteAFile    			("_MFReader_.ibf");
	_runAsFunctionLibrary = 0;

}

siteToPartititon = {totalCodonCount,2};
shifter			 = 0;
for (k=0; k<fileCount; k=k+1)
{
	ExecuteCommands ("mySize = nucData_" + (k+1) + ".sites");
	for (k2 = 0; k2 < mySize; k2 = k2+1)
	{
		siteToPartititon[k2+shifter][0] = k;
		siteToPartititon[k2+shifter][1] = k2;
	}
	shifter = shifter + mySize;
}

	
referenceL 	= {fileCount, 1};	
	
for (k=0; k<fileCount; k=k+1)
{
	ExecuteCommands			 		("baselineBL = BranchLength (givenTree_"+(k+1)+",-1);");
	referenceL[k]					= (baselineBL * (Transpose(baselineBL)["1"]))[0];
}

summaryPath					   = basePath+".summary";
substitutionsPath			   = basePath+"_subs.csv";
siteReportMap				   = basePath+"_bysite.csv";
fprintf 						(summaryPath, CLEAR_FILE, KEEP_OPEN);
fprintf							(stdout,      "[PHASE 0.2] Standard model fit. Log-L = ",baselineLogL,".\n");
fprintf							(summaryPath, "[PHASE 0.2] Standard model fit. Log-L = ",baselineLogL,".\n");

for (k=0; k<fileCount; k=k+1)
{
	fprintf (stdout, "\tTree ",k+1," length = ",referenceL[k], " subs/site \n"); 
	fprintf (summaryPath, "\tTree ",k+1," length = ",referenceL[k], " subs/site \n"); 
}

ExecuteAFile 					(HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "Utility" + DIRECTORY_SEPARATOR + "GrabBag.bf");
ExecuteAFile 					(HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "Utility" + DIRECTORY_SEPARATOR + "AncestralMapper.bf");


fixGlobalParameters ("lf");

byResidueSummary = {};
bySiteSummary	 = {};



/*------------------------------------------------------------------------------*/

if (MPI_NODE_COUNT > 1)
{
	MPINodeStatus = {MPI_NODE_COUNT-1,1}["-1"];
}



for (residue = 0; residue < 20; residue = residue + 1)
{
	AddABiasREL 					(modelNameString,"biasedMatrix",residue);
	Model							biasedModel = (biasedMatrix, vectorOfFrequencies, 1);
	
	populateTrees 					("biasedTree", fileCount);	
	copyTreeBL	 					("biasedTree", "givenTree", "treeScaler", fileCount);	
	ExecuteCommands				     (constructLF ("lfb", "nucData", "biasedTree", fileCount));
	
	
	if (MPI_NODE_COUNT > 1)
	{
		SendAJob (residue);
	}
	else
	{
		Optimize 						(lfb_MLES,lfb);
		DoResults 						(residue);
	}
}

/*------------------------------------------------------------------------------*/

if (MPI_NODE_COUNT > 1)
{
	jobsLeft = ({1,MPI_NODE_COUNT-1}["1"] * MPINodeStatus["_MATRIX_ELEMENT_VALUE_>=0"])[0];

	for (nodeID = 0; nodeID < jobsLeft; nodeID = nodeID + 1)
	{
		MPIReceive (-1, fromNode, theJob);
		oldRes = MPINodeStatus[fromNode-1];
		ExecuteCommands (theJob);
		DoResults (oldRes);	
	}
}

/*------------------------------------------------------------------------------*/

fprintf							(substitutionsPath, CLEAR_FILE, KEEP_OPEN, "Site,From,To,Count");
fprintf							(siteReportMap,     CLEAR_FILE, KEEP_OPEN, "Site");
for (k=0; k<20; k=k+1)
{
	fprintf (siteReportMap, ",", AAString[k]);
}
fprintf (siteReportMap, "\nLRT p-value"); 

test_p_values       = test_p_values % 0;
rejectedHypotheses   = {};

for (k=0; k<20; k=k+1)
{
	pv      = (byResidueSummary[AAString[k]])["p"];
	fprintf (siteReportMap, ",", pv);
}

fprintf (stdout, 	  "\nResidues (and p-values) for which there is evidence of directional selection\n");
fprintf (summaryPath, "\nResidues (and p-values) for which there is evidence of directional selection");

for (k=0; k<20; k=k+1)
{
	if (test_p_values[k][0] < (0.05/(20-k)))
	{
		rejectedHypotheses  [test_p_values[k][1]]           = 1;
		rejectedHypotheses  [AAString[test_p_values[k][1]]] = 1;
		fprintf (stdout, 		"\n\t", AAString[test_p_values[k][1]], " : ",test_p_values[k][0] );
		fprintf (summaryPath, 	"\n\t", AAString[test_p_values[k][1]], " : ",test_p_values[k][0] );
	}
	else
	{
		break;
	}
}

fprintf (stdout, 	  "\n");
fprintf (summaryPath, "\n");

ancCacheID = {fileCount,1};
for (k=0; k<fileCount; k=k+1)
{
	ancCacheID 	[k]					= _buildAncestralCache ("lf", k);
}
outputcount						= 0;

for (k=0; k<totalCodonCount; k=k+1)
{
	thisSite = _substitutionsBySite (ancCacheID[siteToPartititon[k][0]],siteToPartititon[k][1]);
	
	for (char1 = 0; char1 < 20; char1 = char1+1)
	{
		for (char2 = 0; char2 < 20; char2 = char2+1)
		{
			if (char1 != char2 && (thisSite["COUNTS"])[char1][char2])
			{	
				ccount = (thisSite["COUNTS"])[char1][char2];
				fprintf (substitutionsPath, "\n", k+1, ",", AAString[char1], ",", AAString[char2], "," , ccount);
			}
		}
	}
	

	if (Abs(bySiteSummary[k]))
	{
		fprintf (siteReportMap, "\n", k+1);
		
		didSomething = 0;
		for (k2=0; k2<20; k2=k2+1)
		{
			if (Abs((byResidueSummary[AAString[k2]])["BFs"]) == 0 || rejectedHypotheses[k2] == 0)
			{
				fprintf (siteReportMap, ",N/A");
			}
			else
			{
				pv = ((byResidueSummary[AAString[k2]])["BFs"])[k];
				fprintf (siteReportMap, ",", pv);			
				if (pv > BayesFactorCutoff)
				{
					didSomething = 1;
				}
			}
		}
		
		if (!didSomething)
		{
			continue;
		}
		
		if (outputcount == 0)
		{
			outputcount = 1;
			fprintf (stdout, 		"\nThe list of sites which show evidence of directional selection (Bayes Factor > 20)\n",
							 		"together with the target residues and inferred substitution counts\n");
			fprintf (summaryPath, 	"\nThe list of sites which show evidence of directional selection (Bayes Factor > 20)\n",
							 		"together with the target residues and inferred substitution counts\n");
		}	
		fprintf (stdout,      "\nSite ", Format (k+1,8,0), "\n Preferred residues: ");
		fprintf (summaryPath, "\nSite ", Format (k+1,8,0), "\n Preferred residues: ");
		
		
		for (k2 = 0; k2 < Abs (bySiteSummary[k]); k2=k2+1)
		{
			thisChar = (bySiteSummary[k])[k2];
			if (rejectedHypotheses[thisChar])
			{
				fprintf (stdout,      thisChar);
				fprintf (summaryPath, thisChar);
			}
		}

		fprintf (stdout,      	   "\n Substitution counts:");
		fprintf (summaryPath,      "\n Substitution counts:");

		for (char1 = 0; char1 < 20; char1 = char1+1)
		{
			for (char2 = char1+1; char2 < 20; char2 = char2+1)
			{
				ccount  = (thisSite["COUNTS"])[char1][char2];
				ccount2 = (thisSite["COUNTS"])[char2][char1];
				if (ccount+ccount2)
				{	
					fprintf (stdout, 	  "\n\t", AAString[char1], "->", AAString[char2], ":", Format (ccount, 5, 0), "/",
											 AAString[char2], "->", AAString[char1], ":", Format (ccount2, 5, 0));
					fprintf (summaryPath, "\n\t", AAString[char1], "->", AAString[char2], ":", Format (ccount, 5, 0), "/",
											 AAString[char2], "->", AAString[char1], ":", Format (ccount2, 5, 0));
				}
			}
		}

	}
}	

for (k=0; k<fileCount; k=k+1)
{
	_destroyAncestralCache 			(ancCacheID[k]);
}
fprintf (substitutionsPath, CLOSE_FILE);
fprintf (summaryPath, 		CLOSE_FILE);
fprintf (siteReportMap, 	CLOSE_FILE);
fprintf (stdout, "\n");

/*--------------------------------------------------------------------------------------------*/

function computeDelta (ModelMatrixName&, efv, t_0, which_cat)
{
	t   	= t_0;
	c   	= 1;
	catVar  = which_cat;
	rmx 	= ModelMatrixName;
	for (r=0; r<20; r=r+1)
	{	
		diag = 0;
		for (c=0; c<20; c=c+1)
		{
			rmx[r][c] = rmx[r][c] * efv[c];
			diag = diag - rmx[r][c];
		}
		rmx[r][r] = diag;
	}
	return Transpose(efv)*(Exp (rmx) - {20,20}["_MATRIX_ELEMENT_ROW_==_MATRIX_ELEMENT_COLUMN_"]);
}

/*------------------------------------------------------------------------------*/

function SendAJob (residueIn)
{
	for (nodeID = 0; nodeID < MPI_NODE_COUNT -1; nodeID = nodeID + 1)
	{
		if (MPINodeStatus[nodeID] < 0)
		{
			MPINodeStatus[nodeID] = residueIn;
			MPISend (nodeID+1,lfb);
			break;
		}
	}
	if (nodeID == MPI_NODE_COUNT - 1)
	{
		MPIReceive (-1, fromNode, theJob);
		oldRes = MPINodeStatus[fromNode-1];
		MPISend (fromNode,lfb);
		MPINodeStatus[fromNode-1] = residueIn;
		ExecuteCommands (theJob);
		DoResults (oldRes);
	}
	return 0;
}

/*------------------------------------------------------------------------------*/

function DoResults (residueIn)
{
	residueC 					= 	AAString[residueIn];
	fprintf							(stdout, "[PHASE ",residueIn+1,".1] Model biased for ",residueC,"\n"); 
	fprintf							(summaryPath, "[PHASE ",residueIn+1,".1] Model biased for ",residueC,"\n"); 

	pv							=   1-CChi2(2(lfb_MLES[1][0]-baselineLogL),2);
	fprintf							(stdout, "[PHASE ",residueIn+1,".2] Finished with the model biased for ",residueC,". Log-L = ",Format(lfb_MLES[1][0],8,3),"\n"); 
	fprintf							(summaryPath, "[PHASE ",residueIn+1,".2] Finished with the model biased for ",residueC,". Log-L = ",Format(lfb_MLES[1][0],8,3),"\n"); 
	
	fr1 						= 	P_bias;
	
	
	
	fprintf							(stdout, "\n\tBias term                = ", Format(rateBiasTo,8,3),
											 "\n\tproportion               = ", Format(fr1,8,3),
											 "\n\tp-value         = ", Format(pv,8,3));

	fprintf							(summaryPath, "\n\tBias term                = ", Format(rateBiasTo,8,3),
											 	  "\n\tproportion               = ", Format(fr1,8,3),
											      "\n\tp-value         = ", Format(pv,8,3));

	for		(rc = 0; rc < fileCount; rc = rc + 1)
	{
		rateAccel1					=   (computeDelta("biasedMatrix",vectorOfFrequencies,referenceL,1))[residueIn];
		fprintf (stdout, 			"\n\tExp freq increase (tree ", Format(rc+1,3,0), ")  = ", Format(rateAccel1*100,8,3), "%");
		fprintf (summaryPath, 		"\n\tExp freq increase (tree ", Format(rc+1,3,0), ")  = ", Format(rateAccel1*100,8,3), "%");
	}
	
	fprintf							(stdout,			"\n");
	fprintf							(summaryPath, 		"\n");
		
		
											 											
	LIKELIHOOD_FUNCTION_OUTPUT = 7;
	outPath = basePath + "." + residueC;
	fprintf (outPath, CLEAR_FILE, lfb);

	byResidueSummary [residueC] = {};
	(byResidueSummary [residueC])["p"] = pv;		

	test_p_values [residueIn][0] = pv;
	test_p_values [residueIn][1] = residueIn;

	/*if (pv < 0.0025)*/
	{
		(byResidueSummary [residueC])["sites"] = {};		
		(byResidueSummary [residueC])["BFs"]   = {};		
		
		ConstructCategoryMatrix (mmx,lfb,COMPLETE);
		GetInformation			(catOrder, lfb);		
		dim = Columns (mmx);
		_MARGINAL_MATRIX_	= {2, dim};
		
		GetInformation 				(cInfo, c);
		GetInformation 				(_CATEGORY_VARIABLE_CDF_, catVar);
		
		ccc	= Columns (cInfo);
		
		_CATEGORY_VARIABLE_CDF_ = _CATEGORY_VARIABLE_CDF_[1][-1];
		if (catOrder [0] == "c")
		{
			for (k=0; k<dim; k=k+1)
			{
				for (k2 = 0; k2 < ccc; k2=k2+1)
				{
					_MARGINAL_MATRIX_ [0][k] = _MARGINAL_MATRIX_ [0][k] + mmx[2*k2][k]  *cInfo[1][k2];
					_MARGINAL_MATRIX_ [1][k] = _MARGINAL_MATRIX_ [1][k] + mmx[2*k2+1][k]*cInfo[1][k2];
				}
			}
		}
		else
		{
			for (k=0; k<dim; k=k+1)
			{
				for (k2 = 0; k2 < ccc; k2=k2+1)
				{
					_MARGINAL_MATRIX_ [0][k] = _MARGINAL_MATRIX_ [0][k] + mmx[k2][k]*cInfo[1][k2];
					_MARGINAL_MATRIX_ [1][k] = _MARGINAL_MATRIX_ [1][k] + mmx[ccc+k2][k]*cInfo[1][k2];
				}
			}
		}
		ExecuteAFile 					(HYPHY_LIB_DIRECTORY + "ChartAddIns" + DIRECTORY_SEPARATOR + "DistributionAddIns" + DIRECTORY_SEPARATOR + "Includes" + DIRECTORY_SEPARATOR + "posteriors.ibf");
		
		prior = (_CATEGORY_VARIABLE_CDF_[1])/(1-_CATEGORY_VARIABLE_CDF_[1]);
				
		for (k=0; k<dim; k=k+1)
		{
			bayesF = _MARGINAL_MATRIX_[1][k]/_MARGINAL_MATRIX_[0][k]/prior;
			((byResidueSummary [residueC])["BFs"])[k] = bayesF;
			if (bayesF > BayesFactorCutoff)
			{
				((byResidueSummary [residueC])["sites"])[Abs((byResidueSummary [residueC])["sites"])] = k+1;
				if (Abs(bySiteSummary[k]) == 0)
				{
					bySiteSummary[k] = {};
				}
				(bySiteSummary[k])[Abs(bySiteSummary[k])] = residueC;
			}
		}
		
	}	
	return 0;
}


