
/* define an associative array with key = amino acid and value = integer
	index for AA in alphabetical order */
AAString    = "ACDEFGHIKLMNPQRSTVWY";

AACharToIdx = {};		/* this is never used :-P */
for (k=0; k<20; k=k+1)
{
	AACharToIdx [AAString[k]] = k;
}

SKIP_MODEL_PARAMETER_LIST = 0;

#include "AddABias.ibf";
#include "Utility/GrabBag.bf";

test_p_values = {20,2};

/*--------------------------------------------------------------------------------------------*/
/*	
	returns stationary distribution vector for amino acids
	UNUSED?!?
*/
function GetEqFreqs (ModelMatrixName&, baseFreqs)
{
	t = 1;
	
	/* duplicate instantaneous rate parameterization from template  model */
	numRateMx = ModelMatrixName;
	
	for (ri = 0; ri < 20; ri = ri+1)
	{
		for (ci = 0; ci < 20; ci = ci+1)
		{
			if (ri != ci)
			{
				/* multiply entries by base frequencies */
				numRateMx[ri][ci] = numRateMx[ri][ci] * baseFreqs[ci];
				
				/* diagonal entry is sum of off-diagonal row entries (Markov process) */
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

SetDialogPrompt ("Protein file to examine:");

ChoiceList						(reloadFlag, "Reload/New", 1, SKIP_NONE, "New analysis","Start a new analysis",
																	      "Reload","Reload a baseline protein fit.",
																	      "Reprocess the results","Regenerate sites from model fits.");


ACCEPT_ROOTED_TREES 			= 1;

if (reloadFlag == 0)
{
	/* new analysis, fit baseline model */
	
	/* this should be a protein alignment */
	DataSet			ds 			 = ReadDataFile (PROMPT_FOR_FILE);
	basePath 					 = LAST_FILE_PATH;
	DataSetFilter   filteredData = CreateFilter (ds,1);
	
	GetDataInfo		(checkCharacters, filteredData, "CHARACTERS");
	if (Columns (checkCharacters) != 20)
	{
		fprintf (stdout, "ERROR: please ensure that the input alignment contains protein sequences");
		return 0;
	}
	
	/* from AddABias.ibf; select from one of several AA rate matrices estimated from
		curated alignments (e.g., HIV within); by default, set rate variation to
		adaptive 4-bin beta-gamma -- argument 0 skips dialog to set [pickATarget] */
	
	promptModel (0);
	
	
	ExecuteAFile 					(HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "queryTree.bf");
	ChoiceList						(fixFB, "Fix Branch", 1, SKIP_NONE, "Unknown root","The character at the root of the tree is drawn from the stationary distribution",
																		"Fix 1st sequence as root","The 1st sequence in the file (assumed to be a direct descendant of the root) is used to populate the root sequences.");
	
	
	/* check if the tree is rooted */
	
	treeAVL  = givenTree^0;
	rootNode = treeAVL[(treeAVL[0])["Root"]];
	if (Abs(rootNode["Children"]) != 2)
	{
		fprintf (stdout, "ERROR: please ensure that the tree is rooted");
		return 0;
	}
	root_left  = "givenTree." + (treeAVL[(rootNode["Children"])[0]])["Name"] + ".t";
	root_right = "givenTree." + (treeAVL[(rootNode["Children"])[1]])["Name"] + ".t";
	
	if (fixFB>0)
	{
		/* branch to first sequence is collapsed to length zero;
			enforcing identifiability with root sequence */
		ExecuteCommands ("givenTree."+TipName(givenTree,0)+".t:=0");
	}
	else
	{
		if (fixFB < 0)
		{
			return 0;
		}
	}
	
	/* only the sum of branch lengths leading to two immediate descendants
		of root can be estimated; this measure prevents one of the branches
		from collapsing to zero length */
	ExecuteCommands					(root_left + ":=" + root_right); 
	
	LikelihoodFunction lf 		= 	(filteredData, givenTree);
	fprintf							(stdout, "[PHASE 0.1] Standard model fit\n"); 
	
	
	VERBOSITY_LEVEL				= 1;
	AUTO_PARALLELIZE_OPTIMIZE	= 1;
	Optimize 						(res0,lf);
	AUTO_PARALLELIZE_OPTIMIZE	= 0;
	VERBOSITY_LEVEL				= -1;
	
	
	/* export baseline model LF */
	LIKELIHOOD_FUNCTION_OUTPUT = 7;
	outPath = basePath + ".base";
	fprintf (outPath, CLEAR_FILE, lf);
	
	baselineLogL				= res0[1][0];
	
}
else
{
	/* import baseline LF from file */
	modelNameString = "_customAAModelMatrix";
	SetDialogPrompt ("Locate an existing fit:");
	ExecuteAFile (PROMPT_FOR_FILE);
	GetString (lfInfo,lf,-1);
	if ((lfInfo["Models"])[0] == "mtREVModel")
	{
		modelNameString = "mtREVMatrix";	
	}
	bpSplit						 = splitFilePath (LAST_FILE_PATH);
	basePath					 = bpSplit["DIRECTORY"] + bpSplit["FILENAME"];
	outPath						 = basePath + ".base";
	treeString 					 = Format (givenTree,0,0);	/* topology only, no branch lengths */

	treeAVL  = givenTree^0;
	rootNode = treeAVL[(treeAVL[0])["Root"]];
	if (Abs(rootNode["Children"]) != 2)
	{
		fprintf (stdout, "ERROR: please ensure that the tree is rooted");
		return 0;
	}


	ChoiceList						(fixFB, "Fix Branch", 1, SKIP_NONE, "Unknown root","The character at the root of the tree is drawn from the stationary distribution",
																		"Fix 1st sequence as root","The 1st sequence in the file (assumed to be a direct descendant of the root) is used to populate the root sequences.");
	if (fixFB>0)
	{
		ExecuteCommands ("givenTree."+TipName(givenTree,0)+".t:=0");
	}
	else
	{
		if (fixFB < 0)
		{
			return 0;
		}
	}
	
	LFCompute (lf,LF_START_COMPUTE);
	LFCompute (lf,baselineLogL);
	LFCompute (lf,LF_DONE_COMPUTE);
}


/* when the heck does biasedTree get called?  It's tucked inside runAFit() in AddABias.ibf... */
root_left  = "biasedTree." + (treeAVL[(rootNode["Children"])[0]])["Name"] + ".t";
root_right = "biasedTree." + (treeAVL[(rootNode["Children"])[1]])["Name"] + ".t";


/* vector of branch lengths for entire tree */
baselineBL						= BranchLength (givenTree,-1);

referenceL						= (baselineBL * (Transpose(baselineBL)["1"]))[0];

summaryPath					   = basePath+".summary";
substitutionsPath			   = basePath+"_subs.csv";
siteReportMap				   = basePath+"_bysite.csv";
fprintf 						(summaryPath, CLEAR_FILE, KEEP_OPEN);
fprintf							(stdout,      "[PHASE 0.2] Standard model fit. Log-L = ",baselineLogL,". Tree length = ",referenceL, " subs/site \n"); 
fprintf							(summaryPath, "[PHASE 0.2] Standard model fit. Log-L = ",baselineLogL,". Tree length = ",referenceL, " subs/site \n"); 


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
	if (reloadFlag == 2)	/* regenerate sites from results files */
	{
		lfb_MLES = {2,1};
		ExecuteAFile (basePath + "." + AAString[residue]);
		LFCompute (lfb,LF_START_COMPUTE);
		LFCompute (lfb,res);
		LFCompute (lfb,LF_DONE_COMPUTE);	
		lfb_MLES [1][0] = res;
		DoResults 						(residue);
	}
	else
	{
		AddABiasREL 					(modelNameString,"biasedMatrix",residue);	/* returns biasedMatrix object */
		global P_bias2 					:= 1;
		global relBias					:= 1;
		
		/* vectorOfFrequencies comes from Custom_AA_empirical.mdl, in turn imported from a file such as "HIVWithin" */
		/* rate matrix is multiplied by this vector (third argument) */
		Model							biasedModel = (biasedMatrix, vectorOfFrequencies, 1);
		Tree							biasedTree = treeString;
		global							treeScaler = 1;
		
		/* constrains tree to be congruent to that estimated under baseline model */
		ReplicateConstraint 			("this1.?.?:=treeScaler*this2.?.?__",biasedTree,givenTree);
		ExecuteCommands					(root_left + "=" + root_left);
		ExecuteCommands					(root_right + "=" + root_right);
		LikelihoodFunction lfb 		= 	(filteredData, biasedTree);
		
		
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


ancCacheID 						= _buildAncestralCache ("lf", 0);
outputcount						= 0;

for (k=0; k<filteredData.sites; k=k+1)
{
	thisSite = _substitutionsBySite (ancCacheID,k);
	
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
		pv			 = 0;
		for (k2=0; k2<20; k2=k2+1)
		{
			if (Abs((byResidueSummary[AAString[k2]])["BFs"]) == 0 || rejectedHypotheses[k2] == 0)
			{
				fprintf (siteReportMap, ",N/A");
			}
			else
			{
				thisSitePV = ((byResidueSummary[AAString[k2]])["BFs"])[k];
				pv = Max(pv,thisSitePV);
				fprintf (siteReportMap, ",", thisSitePV);			
				if (pv > 100)
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
		fprintf (stdout,      "\nSite ", Format (k+1,8,0), " (max BF = ", pv, ")\n Preferred residues: ");
		fprintf (summaryPath, "\nSite ", Format (k+1,8,0), " (max BF = ", pv, ")\n Preferred residues: ");
		
		
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

_destroyAncestralCache 			(ancCacheID);
fprintf (substitutionsPath, CLOSE_FILE);
fprintf (summaryPath, 		CLOSE_FILE);
fprintf (siteReportMap, 	CLOSE_FILE);
fprintf (stdout, "\n");

/*--------------------------------------------------------------------------------------------*/
/* 
	Compute the difference vector between the stationary distribution at the root (efv) and the expected
	distribution of residues after time t0.
*/
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
		if (MPINodeStatus[nodeID] < 0)	/* semaphore indicates available node */
		{
			MPINodeStatus[nodeID] = residueIn;
			MPISend (nodeID+1,lfb);
			break;
		}
	}
	if (nodeID == MPI_NODE_COUNT - 1)	/* looped through semaphore with no idle nodes */
	{
		MPIReceive (-1, fromNode, theJob);	/* wait until one of the node completes its job */
		oldRes = MPINodeStatus[fromNode-1];
		MPISend (fromNode,lfb);
		MPINodeStatus[fromNode-1] = residueIn;
		ExecuteCommands (theJob);	/* load serialized LF into master node memory */
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

	pv							=   1-CChi2(2(lfb_MLES[1][0]-baselineLogL),3);	/* approximate p-value */
	fprintf							(stdout, "[PHASE ",residueIn+1,".2] Finished with the model biased for ",residueC,". Log-L = ",Format(lfb_MLES[1][0],8,3),"\n"); 
	fprintf							(summaryPath, "[PHASE ",residueIn+1,".2] Finished with the model biased for ",residueC,". Log-L = ",Format(lfb_MLES[1][0],8,3),"\n"); 
	
	fr1 						= 	P_bias;
	
	rateAccel1					=   (computeDelta("biasedMatrix",vectorOfFrequencies,referenceL,1))[residueIn];
	
	fprintf							(stdout, "\n\tBias term           = ", Format(rateBiasTo,8,3),
											 "\n\tproportion          = ", Format(fr1,8,3),
											 "\n\tExp freq increase   = ", Format(rateAccel1*100,8,3), "%",
											 "\n\tp-value    = ", Format(pv,8,3),"\n");
											 
	fprintf							(summaryPath, "\n\tBias term           = ", Format(rateBiasTo,8,3),
											 	  "\n\tproportion          = ", Format(fr1,8,3),
											      "\n\tExp freq increase   = ", Format(rateAccel1*100,8,3), "%",
											      "\n\tp-value    = ", Format(pv,8,3),"\n");
	if (reloadFlag != 2)
	{
		LIKELIHOOD_FUNCTION_OUTPUT = 7;
		outPath = basePath + "." + residueC;
		fprintf (outPath, CLEAR_FILE, lfb);
	}

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
			if (bayesF > 100)
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


