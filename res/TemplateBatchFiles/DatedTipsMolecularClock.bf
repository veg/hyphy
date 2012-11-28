RequireVersion ("2.0020110101");
VERBOSITY_LEVEL = -1;


LoadFunctionLibrary ("TreeTools");
ExecuteAFile 		("_tipDater.ibf");


/*-----------------------------------------------------------*/

function generateDatedTipConstraints (treeNameID, parameterToConstrain, tipDateAVL, rateAVL, initialGuesses)
{
	DT_String = "";
	DT_String * 8192;
	

	treePostOrderAVL = Eval(treeNameID+"^0");
	treePreOrderAVL  = Eval(treeNameID+"^1");
	internalNodeCount = Eval("BranchCount(`treeNameID`)");
	nodeCount        =  Abs(treePostOrderAVL);
	
	rateAVL          = {};
	doneAssignment   = {};

    // define separate global variables for each rate class (currently not used)
	for (nodeIndex = 0; nodeIndex < nodeCount-1; nodeIndex += 1)
	{
		aRate = rateAVL[nodeIndex];
		if (doneAssignment[aRate] == 0)
		{
			doneAssignment[aRate] = 1;
			DT_String * ("global "+treeNameID+"_scaler_"+aRate+" = 1.0;\n"+treeNameID+"_scaler_"+aRate+"_ :> 0.0;\n");
		}
	}

	timeStops				  = {};
    
    

	for (nodeIndex = 1; nodeIndex < nodeCount; nodeIndex += 1)
	{
		nodeInfo 	= treePostOrderAVL[nodeIndex];
		nodeNameS	= nodeInfo["Name"];
		timeStops[nodeNameS] = 1e100;
	
		if (Abs(nodeInfo["Children"]))
		{
			DT_String * ("\n\n`treeNameID`.`nodeNameS`.T = 1;\n`treeNameID`.`nodeNameS`.T:>(0);\n");
			if (Abs(nodeInfo["Parent"]) == 0)
			{
				
				if (minV > 0)
				{
					minV = minV/2;
				}
				else
				{
					minV = minV*2;
				}
                
				DT_String * ("`treeNameID`.`nodeNameS`.T = " + minV + ";");
				timeStops[nodeNameS] = minV;
			}
			else
			{
				DT_String * ("`treeNameID`.`nodeNameS`.BL = 0.0001;\n`treeNameID`.`nodeNameS`.BL :> 0;\n");
			}
		}
		else
		{
			treePostOrderAVL[nodeIndex] = nodeInfo;
		}
	}
	
	descendantsList 		  = {};
	
    PARAMETER_GROUPING = {internalNodeCount,1};
    PARAMETER_GROUPING [0] = "";
    internalsHit = 0;
	for (nodeIndex = 1; nodeIndex < nodeCount; nodeIndex += 1)
	{
		nodeInfo 	= treePostOrderAVL[nodeIndex];
		nodeNameS	= nodeInfo["Name"];
		nodeParent = nodeInfo["Parent"];
		if (nodeParent)
		{
			parentNodeID = nodeParent;
			nodeParent = treePostOrderAVL[nodeParent];
		}
		if (Abs(nodeParent))
		{
			treePostOrderAVL[parentNodeID] = nodeParent;
			pName    = nodeParent["Name"];
			insIndex = Abs (descendantsList[pName]);
			if (insIndex == 0)
			{
				descendantsList[pName] = {};
			}
			
			rateClass = rateAVL[nodeIndex-1];
						
			if (Abs(nodeInfo["Children"]))
			{
			    PARAMETER_GROUPING [internalsHit] =  "`treeNameID`.`nodeNameS`.BL";
			    internalsHit += 1;
				DT_String * ("\n`treeNameID`.`nodeNameS`.`parameterToConstrain`:=`treeNameID`_scaler_"+rateClass+"*("+
							 "`treeNameID`.`nodeNameS`.BL);\n`treeNameID`.`nodeNameS`.T:=`treeNameID`."+nodeParent["Name"]+".T+`treeNameID`."+nodeNameS+".BL;\n");
							 
				(descendantsList[pName])[insIndex] = treeNameID+"."+nodeNameS+".T";
				timeStops[pName] = Min (timeStops[pName],timeStops[nodeNameS]);
			}
			else
			{
				DT_String * ("\n" + treeNameID+"."+nodeNameS+"."+parameterToConstrain+":="+treeNameID+"_scaler_"+rateClass+"*("+
							 tipDateAVL[nodeNameS]+"-"+treeNameID+"."+nodeParent["Name"]+".T);\n");			
				(descendantsList[pName])[insIndex] = tipDateAVL[nodeNameS];
				timeStops[pName] = Min (timeStops[pName],tipDateAVL[nodeNameS]);
			}
		}
	}
	
	fprintf (stdout, PARAMETER_GROUPING, "\n");
	PARAMETER_GROUPING = {"0": PARAMETER_GROUPING};
	
	for (nodeIndex = 1; nodeIndex < nodeCount; nodeIndex = nodeIndex+1)
	{
		nodeInfo 	= treePreOrderAVL[nodeIndex];
		if (Abs(nodeInfo["Children"]) && Abs (nodeInfo["Parent"]))
		{
			nodeNameS	= nodeInfo["Name"];
			pName		= timeStops[(treePreOrderAVL[(nodeInfo["Parent"])])["Name"]];
    
    
            rateClass = initialGuesses[nodeNameS] - pName;
            
            fprintf (stdout, nodeNameS, " initial guess = ", initialGuesses[nodeNameS], ":", rateClass, "\n");

            if (rateClass < 0 || initialGuesses[nodeNameS] >= timeStops[nodeNameS])
            {
                rateClass = Random(0.2,0.8)*(timeStops[nodeNameS] - pName);
            }
			DT_String * (treeNameID+"."+nodeNameS+".BL = " + rateClass + ";\n");
			timeStops[nodeNameS] = pName + rateClass;			
		}
	}


	rateClass 			= Rows (descendantsList); 
	minTimeByNode		= {};
	
	for (nodeIndex = 0; nodeIndex < Columns (rateClass); nodeIndex = nodeIndex + 1)
	{
		pName 				= rateClass[nodeIndex];
		nodePT				= descendantsList[pName];
		doneAssignment		= 1e100;
		
		for (nodeInfo = 0; nodeInfo< Abs (nodePT); nodeInfo = nodeInfo + 1)
		{
			nodeNameS = nodePT[nodeInfo];
			if (Type (nodeNameS) == "Number")
			{
				doneAssignment = Min (doneAssignment, nodeNameS);
			}
			else
			{
				break;
			}
		}
		minTimeByNode[pName]= doneAssignment;
	}

	for (nodeIndex = 1; nodeIndex < nodeCount; nodeIndex = nodeIndex+1)
	{
		nodeInfo 	= treePostOrderAVL[nodeIndex];
		if (Abs(nodeInfo["Children"]) )
		{
			nodeNameS	= nodeInfo["Name"];
			if (Abs (nodeInfo["Parent"]))
			{
				pName		= (treePostOrderAVL[nodeInfo["Parent"]])["Name"];
				minTimeByNode[pName] = Min(minTimeByNode[pName], minTimeByNode[nodeNameS]);
			}
			
			DT_String * (treeNameID+"."+nodeNameS+".T:<"+minTimeByNode[nodeNameS]+";\n");
		}
	}
	
	if (initialGuesses["Root"] > 0.0)
	{
	    DT_String * ("\n`treeNameID`.Node0.T = 0.5*" + initialGuesses["Root"] + ";");
 	}
 	else
 	{
	    DT_String * ("\n`treeNameID`.Node0.T = 2.0*" + initialGuesses["Root"] + ";");

 	}
 	DT_String * 0;
 	fprintf (stdout, DT_String, "\n");
 	 	
	return DT_String;
}

/*-----------------------------------------------------------*/

function generateBLVector (treeNameID)
{
	ExecuteCommands ("treePostOrderAVL = "+treeNameID+"^0;");
	_blVector = {};
		
	for (nodeIndex = 1; nodeIndex < nodeCount; nodeIndex = nodeIndex+1)
	{
		nodeInfo 	= treePostOrderAVL[nodeIndex];
		nodeNameS	= nodeInfo["Name"];
		nodeParent = nodeInfo["Parent"];
		if (nodeParent)
		{
			parentNodeID = nodeParent;
			nodeParent = treePostOrderAVL[nodeParent];
		}
		if (Abs(nodeParent))
		{
			pName    = nodeParent["Name"];
			if (Abs(nodeInfo["Children"]))
			{
				ExecuteCommands ("_thisBL = `treeNameID`.`nodeNameS`.BL;tipDateAVL[\"`nodeNameS`\"]=`treeNameID`.`nodeNameS`.T;");
			}
			else
			{
				ExecuteCommands ("_thisBL = "+ tipDateAVL[nodeNameS] +" -" + treeNameID+"."+pName+".T;");			
			}
			_blVector [nodeNameS] = _thisBL;
		}
	}
	
	
	return _blVector;
}

/*-----------------------------------------------------------*/


_DO_TREE_REBALANCE_ = 0;


fprintf(stdout,"\n ---- RUNNING MOLECULAR CLOCK ANALYSIS WITH DATED TIPS---- \n",
			   "The underlying methodology developed by Andrew Rambaut is discussed \nin Bioinformatics 16(4): 395-399\n" );

ChoiceList (dataType,"Data type",1,SKIP_NONE,"Nucleotide/Protein","Nucleotide or amino-acid (protein).",
				     "Codon","Codon (several available genetic codes).");

if (dataType<0) 
{
	return;
}

if (dataType == 1)
{
	LoadFunctionLibrary ("chooseGeneticCode");
}

SetDialogPrompt ("Choose the data file:");
DataSet ds = ReadDataFile (PROMPT_FOR_FILE);
fprintf (stdout,"The following data were read:\n",ds,"\n");

byPosition = 0;

if (dataType)
{
	DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);
}
else
{
	DataSetFilter filteredData = CreateFilter (ds,1);
	GetDataInfo (_charHandles, filteredData, "CHARACTERS");
	if (Columns (_charHandles) == 4)
	{
		ChoiceList (byPosition,"By codon position",1,SKIP_NONE,"Single Partition","Do not split the alignment into codon positions.",
				    		   "Three Partitions","Allow each codon position to evolve at its own rates");

		if (byPosition<0) 
		{
			return;
		}
	}
}


LoadFunctionLibrary ("queryTree");

UseModel 		 		 (USE_NO_MODEL);
Tree			 		 aTree = treeString;


ChoiceList			   (df, "Date Format", 1, SKIP_NONE,
							"TipDate","Taxon names include sampling dates, e.g. ((s1_85,s2_90),t_83)",
							"Branch Lengths", "Branch lengths are the sampling dates, e.g. ((s1:85,s2:90),t:83)");

if (df<0)
{
	return 0;
}			
else
{	
    tipDateAVL = Eval ("getTipDatesFromNames" + (df+1) + "(\"aTree\")");
}

if (Abs(tipDateAVL) == 0)
{
	fprintf (stdout, "\nERROR: Input tree must contain date information in the chosen format for all sequences.\n",
					  "Branch lengths for internal nodes are ignored.\n");
	return 0;
}

fprintf (stdout, "What units are the dates measured in (e.g. months; this is only used for reporting the results):");
fscanf  (stdin,   "String", dateUnit);

fprintf (stdout, "Read the following dates: \n");
seqNames = Rows (tipDateAVL);
for (sid = 0; sid < Columns (seqNames); sid += 1)
{
    fprintf (stdout, seqNames[sid], ":\t", retAVL[seqNames[sid]]/maxV,"\t",dateUnit,"\n");
}



SelectTemplateModel		(filteredData);
Tree					 givenTree = treeString;

if (byPosition)
{
	DataSetFilter filteredData  = CreateFilter (ds,1,"<100>");
	DataSetFilter filteredData2 = CreateFilter (ds,1,"<010>");
	DataSetFilter filteredData3 = CreateFilter (ds,1,"<001>");

	Tree givenTree2 = treeString;
	Tree givenTree3 = treeString;
	global			T_Scale_2 = 1;
	global			T_Scale_3 = 1;
	ReplicateConstraint ("this1.?.?:=T_Scale_2*this2.?.?", givenTree2, givenTree);
	ReplicateConstraint ("this1.?.?:=T_Scale_3*this2.?.?", givenTree3, givenTree);
}


parameter2Constrain = 0;

ACCEPT_ROOTED_TREES = 1;
Tree clockTree = treeString;

if (byPosition)
{
	Tree clockTree2 = treeString;
	Tree clockTree3 = treeString;
	global			TC_Scale_2 = 1;
	global			TC_Scale_3 = 1;
	ReplicateConstraint ("this1.?.?:=TC_Scale_2*this2.?.?", clockTree2, clockTree);
	ReplicateConstraint ("this1.?.?:=TC_Scale_3*this2.?.?", clockTree3, clockTree);

}
ACCEPT_ROOTED_TREES = 0;

if (Rows("LAST_MODEL_PARAMETER_LIST")>1)
{
	ChoiceList (parameter2Constrain, "Parameter(s) to constrain:",1,SKIP_NONE,LAST_MODEL_PARAMETER_LIST);

	if (parameter2Constrain<0)
	{
		return;
	}
	if (parameter2Constrain==0)
	{
		parameter2ConstrainString = "";
		for (parameter2Constrain=Rows("LAST_MODEL_PARAMETER_LIST")-1; parameter2Constrain; parameter2Constrain = parameter2Constrain-1)
		{
			GetString (funnyString,LAST_MODEL_PARAMETER_LIST,parameter2Constrain);
			parameter2ConstrainString = parameter2ConstrainString + funnyString + ",";
		}
		GetString (funnyString,LAST_MODEL_PARAMETER_LIST,0);
		parameter2ConstrainString = parameter2ConstrainString + funnyString;
	}
	else
	{
		GetString (parameter2ConstrainString,LAST_MODEL_PARAMETER_LIST,parameter2Constrain-1);
	}
}
else
{
	GetString (parameter2ConstrainString,LAST_MODEL_PARAMETER_LIST,0);	
}



timer = Time(0);
if (byPosition)
{
	LikelihoodFunction lf = (filteredData,givenTree,filteredData2,givenTree2,filteredData3,givenTree3);
}
else
{
	LikelihoodFunction lf = (filteredData,givenTree);
}
Optimize (res,lf);



separator = "*-----------------------------------------------------------*";

fprintf (stdout, "\nCPU time taken: ", Time(0)-timer, " seconds.\n");
fprintf (stdout, "\n", separator, "\nRESULTS WITHOUT THE CLOCK CONSTRAINT:\n",lf);

fullModelLik = res[1][0];
fullVars 	 = res[1][1];

returnAVL = {"Free Log(L)": fullModelLik,
			 "Free degrees of freedom": fullVars};


/* now specify the constraint */

_initialGuesses = PathDistanceToRoot (givenTree^0, "");

// SLKP -- this needs to be updated to deal with partitions 

linearXY = {filteredData.species, 2};

//--------------------------------------------------------------------------

function popXY (key,value)
{
    linearXY [_k][1] = value;
    linearXY [_k][0] = _initialGuesses[key];
    _k += 1;
    return 0;
}

//--------------------------------------------------------------------------

_k = 0;
tipDateAVL["popXY"][""];

LoadFunctionLibrary ("ProbabilityDistributions");

lfit           = linearFit (linearXY);
ExecuteCommands ("predictor:="+lfit["Slope"]+"*d+"+lfit["Intercept"]);

for (k = 0; k < Abs (_initialGuesses); k += 1)
{
    d = _initialGuesses [_initialGuesses["INDEXORDER"][k]];
    _initialGuesses [_initialGuesses["INDEXORDER"][k]] = predictor;
}

d = 0;
_initialGuesses["Root"] = predictor;

ExecuteCommands (generateDatedTipConstraints  ("clockTree",parameter2ConstrainString,tipDateAVL,0,_initialGuesses));

if (byPosition)
{
	LikelihoodFunction lfConstrained = (filteredData,clockTree,filteredData2,clockTree2,filteredData3,clockTree3);
}
else
{
	LikelihoodFunction lfConstrained = (filteredData, clockTree);
}


sud  	= USE_DISTANCES;
sulr 	= USE_LAST_RESULTS;
slfo	= LIKELIHOOD_FUNCTION_OUTPUT;

USE_DISTANCES 	 				= 0;
USE_LAST_RESULTS 				= 1;
LIKELIHOOD_FUNCTION_OUTPUT		= 0;

timer = Time(0);

VERBOSITY_LEVEL = 10;
OPTIMIZATION_METHOD = 0;

Optimize (res1,lfConstrained);


fprintf (stdout, "\n", separator,"\n\nRESULTS WITH DATED TIPS CLOCK:\nLog-likelihood: ",lfConstrained);
lnLikDiff = 2(fullModelLik-res1[1][0]);
degFDiff = fullVars - res1[1][1];
fprintf (stdout, "\nCPU time taken: ", Time(0)-timer, " seconds.\n");

fprintf (stdout, "\n", separator,"\n\n-2(Ln Likelihood Ratio)=",lnLikDiff,"\n","Constrained parameters:",Format(degFDiff,0,0));
P = 1-CChi2(lnLikDiff,degFDiff);
fprintf (stdout, "\nAsymptotic p-value:",P);


returnAVL ["Clock Log(L)"] = res1[1][0];
returnAVL ["Clock degrees of freedom"] = res1[1][1];
returnAVL ["P"] = P;

ExecuteCommands ("rateEstimate=clockTree_scaler_0*maxV;");
COVARIANCE_PARAMETER = "clockTree_scaler_0";
COVARIANCE_PRECISION = 0.95;
CovarianceMatrix (cmx,lfConstrained);


/* convert to branch lengths */

bls 	= BranchLength (clockTree, -1);
tbc		= Columns 	   (bls);
bns 	= BranchName   (clockTree, -1);
btms	= {tbc,1};
ttl		= (bls*Transpose (bls["1"]))[0];
trl		= 0;
blv		= generateBLVector ("clockTree");

bScaleFactor = 0;
scaledBT	 = {tbc-1,1};
chartLabels  = {1,2};
DT_String * 128;
DT_String * "Node";
for (k=0; k<tbc-1; k=k+1)
{
	nodeNameS 			  = bns[k];
	blv[nodeNameS]  	  = blv[nodeNameS]/maxV;
	trl 	  			  = trl + blv[nodeNameS];
	scaledBT[k]			  = tipDateAVL[nodeNameS]/maxV;
	DT_String * (";" + nodeNameS);
}

DT_String * 0;
chartLabels [0] = "Date (" + dateUnit + ")";
chartLabels [1] = DT_String;
bScaleFactor 	= ttl/trl;

OpenWindow (CHARTWINDOW,{{"Inferred Dates"}
		{"chartLabels"}
		{"scaledBT"}
		{"None"}
		{"Index"}
		{"None"}
		{""}
		{""}
		{""}
		{"0"}
		{""}
		{"0;0"}
		{"10;1.309;0.785398"}
		{"Times:12:0;Times:10:0;Times:12:2"}
		{"0;0;13816530;16777215;0;0;6579300;11842740;13158600;14474460;0;3947580;16777215;15670812;6845928;16771158;2984993;9199669;7018159;1460610;16748822;11184810;14173291"}
		{"16,0,0"}
		},
		"735;540;70;70");
		
datedTree = PostOrderAVL2StringDistances (treePostOrderAVL,blv);

returnAVL ["Dated Tree"] = datedTree;

UseModel (USE_NO_MODEL);
ACCEPT_ROOTED_TREES = 1;
Tree	dT = datedTree;
ACCEPT_ROOTED_TREES = 0;

COVARIANCE_PARAMETER = "clockTree.Node0.T";
COVARIANCE_PRECISION = 0.95;
CovarianceMatrix (cmxT0,lfConstrained);
cmxT0 = cmxT0 * (1/maxV);

fprintf (stdout, "\nTree with branch lengths scaled in ",dateUnit,":\n",datedTree, "\n");
fprintf (stdout, "\nRoot of the tree placed at ", cmxT0[1], "(95% CI: ", cmxT0[0],",", cmxT0[2],") ", dateUnit, "\n");

returnAVL ["Root"] = cmxT0[1];
returnAVL ["Root lower bound"] = cmxT0[0];
returnAVL ["Root upper bound"] = cmxT0[2];

if (byPosition)
{
	fprintf (stdout, "\nSubstitution rate for the first codon position estimated at ", bScaleFactor, "(95% CI: ", cmx[0]/rateEstimate*bScaleFactor,",", cmx[2]/rateEstimate*bScaleFactor,") subs per ",dateUnit," per site\n");
	fprintf (stdout, "Substitution rate for the second codon position estimated at ", TC_Scale_2*bScaleFactor, "(95% CI: ", TC_Scale_2*cmx[0]/rateEstimate*bScaleFactor,",", TC_Scale_2*cmx[2]/rateEstimate*bScaleFactor,") subs per ",dateUnit," per site\n");
	fprintf (stdout, "Substitution rate for the third codon position estimated at ", TC_Scale_3*bScaleFactor, "(95% CI: ", TC_Scale_3*cmx[0]/rateEstimate*bScaleFactor,",", TC_Scale_3*cmx[2]/rateEstimate*bScaleFactor,") subs per ",dateUnit," per site\n");
}
else
{
	fprintf (stdout, "\nSubstitution rate estimated at ", bScaleFactor, "(95% CI: ", maxV*cmx[0]/rateEstimate*bScaleFactor,",", maxV*cmx[2]/rateEstimate*bScaleFactor,") subs per ",dateUnit," per site\n");
	returnAVL ["Rate"] = bScaleFactor;
	returnAVL ["Rate lower bound"] = maxV*cmx[0]/rateEstimate*bScaleFactor;
	returnAVL ["Rate upper bound"] = maxV*cmx[2]/rateEstimate*bScaleFactor;

}
mxTreeSpec = {5,1};
mxTreeSpec [0] = "dT";
mxTreeSpec [3] = "";
mxTreeSpec [4] = "";
mxTreeSpec [1] = "8211";
mxTreeSpec [2] = "100,100,-300,-300,1";
OpenWindow (TREEWINDOW, mxTreeSpec );

VERBOSITY_LEVEL 				= 0;
USE_DISTANCES 	 				= sud;
USE_LAST_RESULTS 				= sulr;
LIKELIHOOD_FUNCTION_OUTPUT		= slfo;

return returnAVL;
