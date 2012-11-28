RequireVersion  ("2.0020110101");

timer = Time(1);

LoadFunctionLibrary ("GrabBag");
LoadFunctionLibrary ("WriteDelimitedFiles");

ExecuteAFile 		("qndhelper1_mf.ibf");

SHORT_MPI_RETURN 	= 1;
SAVE_FIT_TO_FILE = "";

rOptions = 2;
ExecuteAFile 		("qndhelper2_mf.ibf");

fprintf (stdout, "\n");
_in_dNdSPValue = prompt_for_a_value ("Significance level for Likelihood Ratio Tests",0.05,0,1,0);
fprintf (stdout, "MEME will use test p-value of ", _in_dNdSPValue, "\n");

saveNucs = {{AC__,AT__,CG__,CT__,GT__}};

OPTIMIZATION_METHOD = 0;
USE_LAST_RESULTS   = 1;

ClearConstraints (AC,AT,CG,CT,GT);  

for (k = 1; k <= fileCount; k+=1)
{
	ExecuteCommands ("ClearConstraints (codonTree_" + k + ");");
	ExecuteCommands ("ReplicateConstraint(\"this1.?.nonSynRate:=dNdS*this2.?.synRate\",codonTree_"+k+",codonTree_"+k+")");
}

fprintf (stdout, "\n[MEME PHASE 0] Retuning branch lengths and nucleotide rates under the codon model...\n");

finishedPatterns = 0;


T0 = Time(1);
Optimize (codonLF, lf);
OPTIMIZATION_TIME_HARD_LIMIT = (Time(1)-T0)*4;

AUTO_PARALLELIZE_OPTIMIZE = 0;
fprintf (stdout, "Improved Log(L) BY ", codonLF[1][0]-resC[1][0], " points\n");
OPTIMIZATION_METHOD = 4;

LoadFunctionLibrary ("BranchSiteTemplate");

global      mixingP     =         0.5; mixingP :< 1-1e-9; mixingP :> 1e-9;

PopulateModelMatrix              ("MGMatrix1",  positionFrequencies, "alpha", "", "beta1");
PopulateModelMatrix              ("MGMatrix2",  positionFrequencies, "alpha", "", "beta2");
AC := saveNucs__[0];
AT := saveNucs__[1];
CG := saveNucs__[2];
CT := saveNucs__[3];
GT := saveNucs__[4];

Model 		MG1		=		  ("Exp(MGMatrix1)*mixingP+Exp(MGMatrix2)*(1-mixingP)",codonFrequencies,EXPLICIT_FORM_MATRIX_EXPONENTIAL);
Model       MGFEL   =         (MGMatrix2,codonFrequencies,0);

global      sFactor   =  1;
global      nsFactor1 =  1;
nsFactor1               :< 1;
global      nsFactor2 =  1;
global      omega2    =  1;
			omega2    :< 1;

omega2    :> 0;
nsFactor2 :> 0;
nsFactor1 :> 0;
sFactor   :> 0;


doneSites    = {totalUniqueSites,8};
fullSites    = {totalCodonCount,9};						
labels       = {{"Beta-","Weight-","Beta+","Weight+","alpha","LRT","p-value","Full Log(L)"}};

fprintf (stdout, "[MEME PHASE 1]. Fitting a codon model site-by-site\n",
"| Codon |     Beta- |   Weight- |    Beta+  |   Weight+ |     alpha |      LRT  |   p-value |Full Log(L)|   Est. remaining time (secs)\n");


if (MPI_NODE_COUNT > 1) { 								
    MPINodeState 		= {MPI_NODE_COUNT-1,5};
    haz_MPI             = 1;
} else {
    MPINodeState 		= {1,5};
    haz_MPI             = 0;
}

treeLengths 		= {fileCount,1};
MEME_RUN_TIMER    	= Time(1);

vOffset  = 0;
vuOffset = 0;
alreadyDone   = {totalUniqueSites,1};
timesPerSite  = {totalUniqueSites,4} ["(_MATRIX_ELEMENT_ROW_+1)*(_MATRIX_ELEMENT_COLUMN_==3)"];
bySiteCache   = {totalUniqueSites, 3};

LoadFunctionLibrary ("CodonTools.def");

GetString (funcInfo, obtainBranchWiseEBEstimatesMPI, -1);
funcText = "function " + funcInfo["ID"] + "(" + Join (",", funcInfo["Arguments"]) + ") { " + funcInfo ["Body"] + "}";


_memeExtraNull = "SAVE_FIT_TO_FILE = \"`SAVE_FIT_TO_FILE`\";
if (Abs (SAVE_FIT_TO_FILE)) {
    LIKELIHOOD_FUNCTION_OUTPUT = 7;
    saveLFTo = SAVE_FIT_TO_FILE + \".\" + fileID + \".\" + siteID + jobSuffix;
    fprintf (saveLFTo, CLEAR_FILE, siteLikelihood);
    LIKELIHOOD_FUNCTION_OUTPUT = 2;
}
";

 _memeExtra = _memeExtraNull + funcText + 
"_OBSERVED_S_ = " + _OBSERVED_S_ + ";\n" +
"_OBSERVED_NS_ = " + _OBSERVED_NS_ + ";\n" +
"MPI_NEXUS_FILE_RETURN = {};
 MPI_NEXUS_FILE_RETURN [\"MLES\"]     = siteLikelihood_MLES;
 _vtr = {{\"nsFactor1\", \"nsFactor2\", \"sFactor\", \"mixingP\"}};
 MPI_NEXUS_FILE_RETURN [\"VALUES\"]   = {};
 for (_k = 0; _k < Columns (_vtr); _k += 1) {
    (MPI_NEXUS_FILE_RETURN [\"VALUES\"])[_vtr[_k]] = Eval (_vtr[_k]);
 }
 MPI_NEXUS_FILE_RETURN [\"BRANCHES\"] = obtainBranchWiseEBEstimatesMPI (sFactor, nsFactor1, nsFactor2, mixingP);
";

bySiteBranchReports = {};
saveTreesForReports = {};

for (fileID = 1; fileID <= fileCount; fileID = fileID+1)
{
	ClearConstraints (siteTree);
	ClearConstraints (felTree);
	UseModel		  (MG1);
	Tree		   	  siteTree = treeStrings[fileID];
	saveTreesForReports [fileID] = Eval("Format (codonTree_" + fileID +",1,1)");
	UseModel		  (MGFEL);
	Tree		   	  felTree = treeStrings[fileID];
	
	ExecuteCommands ("GetDataInfo  (dupInfo, filteredData_"+fileID+");");			
	ExecuteCommands ("thisFilterSize  = filteredData_"+fileID+".sites;");			
	ExecuteCommands ("thisFilterSizeU = filteredData_"+fileID+".unique_sites;");
	
	ExecuteCommands ("ReplicateConstraint (\"this1.?.alpha:=sFactor*this2.?.synRate__\",siteTree,codonTree_"+fileID+");");
 	ExecuteCommands ("ReplicateConstraint (\"this1.?.beta1:=nsFactor1*sFactor*this2.?.synRate__\",siteTree,codonTree_"+fileID+");");
 	ExecuteCommands ("ReplicateConstraint (\"this1.?.beta2:=nsFactor2*this2.?.synRate__\",siteTree,codonTree_"+fileID+");");
                                
 	ExecuteCommands ("ReplicateConstraint (\"this1.?.alpha:=sFactor*this2.?.synRate__\",felTree,codonTree_"+fileID+");");
 	ExecuteCommands ("ReplicateConstraint (\"this1.?.beta2:=nsFactor2*this2.?.synRate__\",felTree,codonTree_"+fileID+");");

	treeLengths [fileID-1] = + Eval("BranchLength(codonTree_"+fileID+",-1)");
	
	lfSpawnDone 	  = 0;
	//debugVerboseFlag  = 1;
	toDoList      	  = {};
	
	// populate the initial queue of things to do
	
	for (siteCount = 0; siteCount < thisFilterSize; siteCount += 1)
	{
		siteMap = dupInfo[siteCount] + vuOffset;
		if (alreadyDone[siteMap] == 0)
		{
			alreadyDone[siteMap] = 1;				
			filterString = "" + (siteCount*3) + "-" + (siteCount*3+2);
			ExecuteCommands ("DataSetFilter siteFilter = CreateFilter (ds_"+fileID+",3,filterString,\"\",GeneticCodeExclusions);");
			ExecuteCommands ("DataSetFilter felFilter = CreateFilter (ds_"+fileID+",3,filterString,\"\",GeneticCodeExclusions);");

			HarvestFrequencies (f1, siteFilter, 3, 3, 0);
			m1 = +(f1["_MATRIX_ELEMENT_VALUE_>0"]); // how many unique characters?
			
			if (lfSpawnDone == 0)
			{
				LikelihoodFunction siteLikelihood = (siteFilter, siteTree);	
				LikelihoodFunction felLikelihood  = (felFilter, felTree);	
				lfSpawnDone = 1;
			}


			if (m1>1)
			{
				toDoList ["FEL_" + siteCount] = {{siteCount__,siteMap__,siteCount__+vOffset__,2}}; // 2 in the 2nd column means FEL
			}
			else
			{
				doneSites[siteMap][4] = 1;											
			}
			
		}
		
	}
	
	if (debugVerboseFlag)
	{
		fprintf (stdout, toDoList);
	}   

	while (SendJobMEME ()) 
	{
	
	}
					
	vOffset 	= vOffset  + thisFilterSize;
	vuOffset 	= vuOffset + thisFilterSizeU;
}

		
vOffset  = 0;
vuOffset = 0;

alreadyDone = {totalUniqueSites,1};

posSelected = 0;
qValues = {totalCodonCount,2};

codon_labels = {};

for (siteCount = 0; siteCount < totalCodonCount; siteCount += 1){
    qValues [siteCount][0] = siteCount;
    qValues [siteCount][1] = fullSites[siteCount][6];
    codon_labels + (siteCount+1);
}

qValues = qValues % 1;


for (siteCount = 0; siteCount < totalCodonCount; siteCount += 1){
    qValues [siteCount][1] = Min(1,qValues [siteCount][1] * (totalCodonCount) / (1+siteCount));
}

qValues = qValues % 0;

for (siteCount = 0; siteCount < totalCodonCount; siteCount += 1){
    fullSites[siteCount][8] = qValues[siteCount][1];			 
}

fprintf (stdout, "\nMEME analysis done\n");

fprintf (stdout, "\n\n",
"Codon |     Beta- |   Weight- |    Beta+  |   Weight+ |     alpha |      LRT  |   p-value |Full Log(L)|   q-value\n");


for (fileID = 1; fileID <= fileCount; fileID = fileID+1)
{
	ExecuteCommands ("GetDataInfo  (dupInfo, filteredData_"+fileID+");");			
	ExecuteCommands ("thisFilterSize  = filteredData_"+fileID+".sites;");			
	ExecuteCommands ("thisFilterSizeU = filteredData_"+fileID+".unique_sites;");			
	for (siteCount = 0; siteCount < thisFilterSize; siteCount = siteCount+1)
	{
		siteMap 	= dupInfo[siteCount];
		ReportSiteMEME (siteCount+vOffset, siteMap+vuOffset, 0);	
		if (fullSites[siteCount+vOffset][6] < =_in_dNdSPValue && fullSites[siteCount+vOffset][2]>fullSites[siteCount+vOffset][4]) {
			posSelected += 1;
		}
	}
	vOffset 	= vOffset  + thisFilterSize;
	vuOffset 	= vuOffset + thisFilterSizeU;
}


_csv_result_labels = {{"Codon","beta1","weight1","beta2","weight2","alpha","LRT","pvalue","LogL","qvalue"}};

SetDialogPrompt ("Save .csv file to");
WriteSeparatedTable ("", _csv_result_labels, fullSites, codon_labels, ",");

tree_results = LAST_FILE_PATH + ".branches";

fprintf (tree_results, CLEAR_FILE, "Tree,Site,Branch,PosteriorProbability,EmpiricalBayesFactor,SynSubs,NonsynSubs");

vOffset  = 0;
vuOffset = 0;

for (fileID = 1; fileID <= fileCount; fileID = fileID+1)
{
	ExecuteCommands ("GetDataInfo  (dupInfo, filteredData_"+fileID+");");			
	ExecuteCommands ("thisFilterSize  = filteredData_"+fileID+".sites;");			
	ExecuteCommands ("thisFilterSizeU = filteredData_"+fileID+".unique_sites;");			
	for (siteCount = 0; siteCount < thisFilterSize; siteCount += 1)
	{
		siteMap =  dupInfo[siteCount]+vuOffset;
        if (Abs (bySiteBranchReports [siteMap])) {
                bNames = Rows (bySiteBranchReports [siteMap]);
                for (branchCount = 0; branchCount < Columns(bNames); branchCount += 1) {
                    fprintf (tree_results, "\n", fileID, ",", vOffset + siteCount+1, ",", bNames[branchCount], ",", Join (",",(bySiteBranchReports [siteMap])[bNames[branchCount]]));
                }
            }		
	}
	vOffset 	+= thisFilterSize;
	vuOffset 	+= thisFilterSizeU;
}

fprintf (tree_results, CLOSE_FILE);

/*------------------------------------------------------------------------*/

function ReportSiteMEME (siteI, siteM, doEstimatedTime)
{
    // labels       = {{"beta1","beta2","weight1","weight2","alpha","LRT","p-value","Full Log(L)"}};
	fullSites[siteI][0] = doneSites[siteM][0];   // beta1
	fullSites[siteI][1] = doneSites[siteM][2];   // q1
	fullSites[siteI][2] = doneSites[siteM][1];   // beta2
	fullSites[siteI][3] = 1-doneSites[siteM][2]; // q2
	fullSites[siteI][4] = doneSites[siteM][7];   // alpha
	fullSites[siteI][5] = doneSites[siteM][3];   // LRT
	fullSites[siteI][6] = doneSites[siteM][4];   // p-value
	fullSites[siteI][7] = doneSites[siteM][5];   // Log L
	
	// index 8 is a q-value


	/*fprintf (stdout, "| Codon: ", 		Format(siteI+1,4,0),
					 "| Beta1: ", 		Format(fullSites[siteI][0],10,2),
					 "| P(Beta1): ", 	Format(fullSites[siteI][2],5,2),
					 "| Beta2: ", 		Format(fullSites[siteI][1],10,2),
					 "| P(Beta2): ",	Format(fullSites[siteI][3],5,2),
					 "| alpha: ",		Format(fullSites[siteI][4],10,2),
					 "| LRT: ",			Format(fullSites[siteI][5],6,2),
					 "| p: ",			Format(fullSites[siteI][6],5,2),		
					 "| Log(L): ",		Format(fullSites[siteI][7],5,2));*/		

	fprintf (stdout, "|", Format (siteI+1, 6,0));
	for (mxI = 0; mxI < Columns (labels); mxI += 1) {
		fprintf (stdout, " | ",Format(fullSites[siteI][mxI],8,2)," ");
	}

    if (doEstimatedTime) {
	    fprintf (stdout, " | ",Format ((Time(1) - MEME_RUN_TIMER+1)*((totalUniqueSites-finishedPatterns+1)/(1+finishedPatterns)),10,1));
	} else {
		fprintf (stdout, " | ",Format(fullSites[siteI][8],8,2)," ");	
	}
	if (fullSites[siteI][6]<_in_dNdSPValue && fullSites[siteI][2]>fullSites[siteI][4]){
		fprintf (stdout, " *P");		
	}
	fprintf (stdout, "\n");	
	
	return 0;
}

/*------------------------------------------------------------------------*/

function fakeMPIReturnFromLF (lfID, res) {
    ExecuteCommands ("GetString(lfInfo,`lfID`,-1)");
    _fakeMPIAVL = {"VALUES" : {}};
    _fakeMPIAVL["MLES"] = res;
    avlKeys = {{"Local Independent","Global Independent"}};
    for (_kID = 0; _kID < Columns (avlKeys); _kID += 1) {
        key = avlKeys[_kID];    
        for (_vID = 0; _vID < Columns (lfInfo[key]); _vID += 1) {
            _vName = (lfInfo[key])[_vID];
            (_fakeMPIAVL["VALUES"])[_vName] = Eval (_vName);
        }
    }
    return _fakeMPIAVL;
}

/*------------------------------------------------------------------------*/

function ReceiveJobsMEME ()
{
    if (haz_MPI) {
	    MPIReceive (-1, fromNode, result_String);
	} else {
	    fromNode = 1;
	}
	
	siteIndex 		= MPINodeState[fromNode-1][1];
	siteNAF	  		= MPINodeState[fromNode-1][2];
	siteIndexMap	= MPINodeState[fromNode-1][3];
	siteFilterMap   = MPINodeState[fromNode-1][4];
	
 	timesPerSite [siteIndexMap][siteNAF] = Time(1)-timesPerSite [siteIndexMap][siteNAF];

    MPINodeState[fromNode-1][0] = 0;
    MPINodeState[fromNode-1][1] = -1;		
	
	
	
    if (siteNAF == 1)
	{	
	    res                                 = Eval (result_String);
	    siteLikelihood_MLES                 = res["MLES"];
	    siteLikelihood_MLE_VALUES           = res["VALUES"];
	    bySiteBranchReports [siteIndexMap]  = res["BRANCHES"];
	}
	else
	{
	    if (haz_MPI) {
	        ExecuteCommands (result_String);
	    } else {
            res                                 = Eval (result_String);
	        if (siteNAF < 2) {
                siteLikelihood_MLES                 = res["MLES"];
                siteLikelihood_MLE_VALUES           = res["VALUES"];
	        } else {
                felLikelihood_MLES                 = res["MLES"];
                felLikelihood_MLE_VALUES           = res["VALUES"];
	        
	        }
	    }
	}
		
    if (siteNAF < 2)
    {
        nsf1V   = siteLikelihood_MLE_VALUES ["nsFactor1"];
        nsf2V   = siteLikelihood_MLE_VALUES ["nsFactor2"];
        omega2F = siteLikelihood_MLE_VALUES ["omega2"];

        mixingF = siteLikelihood_MLE_VALUES ["mixingP"];
        sFValue = siteLikelihood_MLE_VALUES ["sFactor"];
    }
    else
    {
        nsf2V   = felLikelihood_MLE_VALUES ["nsFactor2"];
        sFValue = felLikelihood_MLE_VALUES ["sFactor"];
    }
	
	if (siteNAF == 1) // alternative
	{
        doneSites[siteIndexMap][0] = nsf1V*sFValue;
        doneSites[siteIndexMap][1] = nsf2V;
        doneSites[siteIndexMap][2] = mixingF;
        doneSites[siteIndexMap][6] = 1-mixingF;                                        
        doneSites[siteIndexMap][7] = sFValue;
		
		doneSites[siteIndexMap][3] = doneSites[siteIndexMap][3]+2*siteLikelihood_MLES[1][0];
		doneSites[siteIndexMap][5] = siteLikelihood_MLES[1][0];
        
        bySiteCache[siteIndexMap][0] = sFValue;
        bySiteCache[siteIndexMap][1] = nsf1V;
        bySiteCache[siteIndexMap][2] = mixingF;
        
        if (debugVerboseFlag)
        {
            fprintf (stdout, "[DEBUG: Received MEME alternative fit of site ", siteIndex, " from node ", fromNode, "]");
            fprintf (stdout, "\n\talpha  = ", doneSites[siteIndexMap][7],
                             "\n\tbeta1  = ", doneSites[siteIndexMap][0],
                             "\n\tbeta2  = ", doneSites[siteIndexMap][1],
                             "\n\tmixing = ", doneSites[siteIndexMap][2],
                             "\n");
        }


       if (nsf2V > sFValue)
        {
            toDoList["MEME_NULL_" + siteIndex] = {{siteFilterMap__,siteIndexMap__,siteIndex__,0}};
            if (debugVerboseFlag)
            {
                fprintf (stdout, "[DEBUG: Added null model fit for site ", siteIndex, " to the queue]\n");
            }       
        }
        else
        {
			finishedPatterns 		  += 1;
            doneSites[siteIndexMap][3] = 0;
            doneSites[siteIndexMap][4] = -1;
        }
	}
	else 
    {
        if (siteNAF == 2) // FEL
        {
            bySiteCache[siteIndexMap][0] = sFValue;
            bySiteCache[siteIndexMap][1] = nsf2V;
            toDoList ["MEME_ALT_" + siteIndex] = {{siteFilterMap__,siteIndexMap__,siteIndex__,1}};
            if (debugVerboseFlag)
            {
                fprintf (stdout, "[DEBUG: Received FEL fit of site ", siteIndex, " from node ", fromNode, "]");
                fprintf (stdout, "\n\talpha  = ", sFValue,
                                 "\n\tbeta  = ", nsf2V,
                                 "\n");
            }
            
        }
        else // null
        {
            doneSites[siteIndexMap][3] = doneSites[siteIndexMap][3]-2*siteLikelihood_MLES[1][0];	
			finishedPatterns 		  += 1;
            if (debugVerboseFlag)
            {
                fprintf (stdout, "[DEBUG: Received MEME NULL fit of site ", siteIndex, " from node ", fromNode, "]");
                fprintf (stdout, "\n\talpha  = ",  sFValue,
                                 "\n\tbeta1  = ",  nsf1V,
                                 "\n\tomega2  = ", omega2F,
                                 "\n\tmixing = ",  mixingF,
                                 "\n");
            }
        }
    }

    if (siteNAF < 2)
    {
        if (doneSites[siteIndexMap][4] == 0)
        {
            doneSites[siteIndexMap][4] = -1;
        }
        else
        {
            if (doneSites[siteIndexMap][4] == (-1))
            {
                doneSites[siteIndexMap][4] = 0.67-0.67*(0.45*CChi2(doneSites[siteIndexMap][3],1)+0.55*CChi2(doneSites[siteIndexMap][3],2));						
                ReportSiteMEME (siteIndex, siteIndexMap,1);
            }
        }
    }
	
	return fromNode-1;
}
/*------------------------------------------------------------------------*/

function SendJobMEME ()
{
    for (mpiNode = 0; mpiNode < Max(1,MPI_NODE_COUNT-1); mpiNode += 1) {
        if (MPINodeState[mpiNode][0]==0)
        {
            break;	
        }
    }
    
    if (mpiNode == Max(1,MPI_NODE_COUNT-1)) {
        mpiNode = ReceiveJobsMEME ();
    }
    
    lastKey = Abs(toDoList);
    
    if (lastKey > 0)
    {
        lastKey  = toDoList ["INDEXORDER"][0];
        theJob   = toDoList [lastKey];
        toDoList - lastKey;

        filterString = "" + (theJob[0]*3) + "-" + (theJob[0]*3+2);
		ExecuteCommands ("DataSetFilter siteFilter = CreateFilter (ds_"+fileID+",3,filterString,\"\",GeneticCodeExclusions);");
		ExecuteCommands ("DataSetFilter felFilter = CreateFilter (ds_"+fileID+",3,filterString,\"\",GeneticCodeExclusions);");
        
        if (theJob[3] == 2)
        {
            sFactor   = 1;
            ClearConstraints (nsFactor2);
            nsFactor2 :> 0;
            nsFactor2 = dNdS;
            OPTIMIZATION_METHOD = 0;
            
            if (haz_MPI) {
                MPISend (mpiNode+1, felLikelihood);
            } else {
                Optimize (res, felLikelihood);
                result_String = ""+fakeMPIReturnFromLF ("felLikelihood", res);
            }
            if (debugVerboseFlag) {
                fprintf (stdout, "[DEBUG: Sending FEL fit of site ", theJob[0], " to node ", mpiNode+1, "]\n");
            }
        }
        else
        {
            if (theJob[3] == 1)
            {
                sFactor   = bySiteCache[theJob[1]][0];
                ClearConstraints (nsFactor2);
                nsFactor2 :> 0;
                nsFactor2 = bySiteCache[theJob[1]][1];
                if (nsFactor2 > sFactor)
                {
                    nsFactor1    = 1;
                    mixingP      = 0.25;
                }
                else
                {
                    nsFactor1    = nsFactor2/sFactor;
                    if (nsFactor2 == 0)
                    {
                        nsFactor2    = sFactor*0.5;
                        mixingP      = 0.05;
                    }
                    else
                    {
                        nsFactor2    = sFactor*1.5;
                        mixingP      = 0.75;
                    }
                }
                OPTIMIZATION_METHOD = 0;
                if (haz_MPI) {
                    LF_NEXUS_EXPORT_EXTRA = "fileID = " + fileID + ";siteID = " + theJob[0] + "; jobSuffix = \".alt\";" + _memeExtra;
                    MPISend (mpiNode+1, siteLikelihood);
                    LF_NEXUS_EXPORT_EXTRA = "";
                } else {
                    Optimize (res, siteLikelihood);
                    result_String = fakeMPIReturnFromLF ("siteLikelihood", res);
                    result_String ["BRANCHES"] = obtainBranchWiseEBEstimates (sFactor, nsFactor1, nsFactor2, mixingP,filterString);
                    result_String = "" + result_String;
                }
                if (debugVerboseFlag)
                {
                    fprintf (stdout, "[DEBUG: Sending MEME fit of site ", theJob[0], " to node ", mpiNode+1, "]");
                    fprintf (stdout, "\n\talpha  = ", sFactor,
                                     "\n\tbeta1  = ", nsFactor1*sFactor,
                                     "\n\tbeta2  = ", nsFactor2, 
                                     "\n\tmixing = ", mixingP,
                                     "\n");
                }
            }
            else
            {
                sFactor      = bySiteCache[theJob[1]][0];
                nsFactor1    = bySiteCache[theJob[1]][1];
                mixingP      = bySiteCache[theJob[1]][2];
                omega2       = 1;
                nsFactor2    := omega2 * sFactor;
                if (sFactor == 0) {
                    sFactor = 0.001;
                }
                OPTIMIZATION_METHOD = 0;
                if (haz_MPI) {
                    LF_NEXUS_EXPORT_EXTRA = "fileID = " + fileID + ";siteID = " + theJob[0] + "; jobSuffix = \".null\";" + _memeExtraNull;
                    MPISend (mpiNode+1, siteLikelihood);
                    LF_NEXUS_EXPORT_EXTRA = "";
                } else {
                    Optimize (res, siteLikelihood);
                    result_String = "" + fakeMPIReturnFromLF ("siteLikelihood", res);
                }
                if (debugVerboseFlag)
                {
                    fprintf (stdout, "[DEBUG: Sending MEME NULL fit of site ", theJob[0], " to node ", mpiNode+1, "]");
                    fprintf (stdout, "\n\talpha  = ", sFactor,
                                     "\n\tbeta1  = ", nsFactor1*sFactor,
                                     "\n\tbeta2  = ", nsFactor2, 
                                     "\n\tmixing = ", mixingP,
                                     "\n");
                }
            }
        }

		timesPerSite [theJob[1]][1] = Time(1);

 		MPINodeState[mpiNode][0] = 1;
		MPINodeState[mpiNode][1] = theJob[2];		
		MPINodeState[mpiNode][2] = theJob[3];		
		MPINodeState[mpiNode][3] = theJob[1];		
        MPINodeState[mpiNode][4] = theJob[0];
        
        /*if (debugVerboseFlag)
        {
            fprintf (fileOut, CLEAR_FILE, MPI_LAST_SENT_MSG);
        }*/
    }
    else
    {
        if (+(MPINodeState [-1][0]))
        {
            ReceiveJobsMEME();
        }
    }
    
    return +(MPINodeState [-1][0]) || Abs(toDoList);
}

// ------------------------------------------------------------------------------------
// ------ MEME helper function --------------------------------------------------------

function obtainBranchWiseEBEstimatesMPI (_sFactor,_nsFactor1,_nsFactor2,_mixingP) {
        if (_nsFactor2 <= _sFactor || _mixingP == 1 || _mixingP == 0)
        {
            return {};
        }
        
        sFactor   = _sFactor;
        nsFactor1 = _nsFactor1;
        nsFactor2 = _nsFactor2;
        mixingP   = _mixingP;
        
        treeString = Format (siteTree,1,1);
        
        LoadFunctionLibrary ("AncestralMapper");
        ancID = _buildAncestralCache ("siteLikelihood",0);
        subMap = _tabulateSubstitutionsAtSiteByBranch (ancID,0);
        _destroyAncestralCache (ancID);

        Model 		MGLocalMix		=		  ("Exp(MGMatrix1)*lmp+Exp(MGMatrix2)*(1-lmp)",codonFrequencies,EXPLICIT_FORM_MATRIX_EXPONENTIAL);
        Tree        perBranchTree 			= treeString;
        ClearConstraints    (perBranchTree);
        ReplicateConstraint ("this1.?.alpha:=this2.?.alpha__",perBranchTree,siteTree);
        ReplicateConstraint ("this1.?.beta1:=this2.?.beta1__",perBranchTree,siteTree);
        ReplicateConstraint ("this1.?.beta2:=this2.?.beta2__",perBranchTree,siteTree);
        ReplicateConstraint ("this1.?.lmp:=_mixingP", perBranchTree);
        
        _bn = BranchName (perBranchTree, -1);
        
		LikelihoodFunction siteLikelihoodLoc = (siteFilter, perBranchTree);
        LFCompute (siteLikelihoodLoc,LF_START_COMPUTE);
        LFCompute (siteLikelihoodLoc,baseline);
        
        _totalBranchCount  = Columns (_bn) - 1;
        posteriorEstimates = {};
        
        _priorOdds = (1-_mixingP)/_mixingP;
        
        for (k = 0; k < _totalBranchCount; k+=1)
        {
             _pname = "perBranchTree." + _bn[k] + ".lmp";
            ExecuteCommands ("`_pname`=1");
            LFCompute (siteLikelihoodLoc,LOGL0);
            
            MaxL     = -Max (LOGL0,baseline);
            
            baseline += MaxL;
            LOGL0 = Exp(MaxL+LOGL0);
            LOGL1 = (Exp(baseline) - _mixingP * LOGL0) / (1-_mixingP);
            
            ExecuteCommands ("`_pname`=_mixingP");
            _posteriorProb = {{LOGL0 * _mixingP, LOGL1 * (1-_mixingP)}};
            _posteriorProb = _posteriorProb * (1/(+_posteriorProb));
            if ( _priorOdds != 0) {
                eBF = _posteriorProb[1] / (1 - _posteriorProb[1]) / _priorOdds;
            } else {
                eBF = 1;
            }
            posteriorEstimates [_bn[k]] = {1,4};
            (posteriorEstimates [_bn[k]])[0] = _posteriorProb[1];
            (posteriorEstimates [_bn[k]])[1] = eBF;
            (posteriorEstimates [_bn[k]])[2] = (subMap[_bn[k]])[0];
            (posteriorEstimates [_bn[k]])[3] = (subMap[_bn[k]])[1];
            baseline += -MaxL;
        } 
        
        LFCompute (siteLikelihoodLoc,LF_DONE_COMPUTE);
        
         
        return posteriorEstimates;
}

// ------ MEME helper function --------

function obtainBranchWiseEBEstimates (_sFactor,_nsFactor1,_nsFactor2, _mixingP,filterString) {

        ClearConstraints    (perBranchTree);
        ReplicateConstraint ("this1.?.alpha:=_sFactor*this2.?.synRate__",perBranchTree,codonTree);
        ReplicateConstraint ("this1.?.beta1:=_nsFactor1*sFactor*this2.?.synRate__",perBranchTree,codonTree);
        ReplicateConstraint ("this1.?.beta2:=_nsFactor2*this2.?.synRate__",perBranchTree,codonTree);
        ReplicateConstraint ("this1.?.lmp:=_mixingP", perBranchTree);
        
        
        LoadFunctionLibrary ("AncestralMapper");
        ancID = _buildAncestralCache ("siteLikelihood",0);
        subMap = _tabulateSubstitutionsAtSiteByBranch (ancID,0);
        _destroyAncestralCache (ancID);

       _bn = BranchName (perBranchTree, -1);
        
        DataSetFilter locSiteFilter = CreateFilter (ds,3,filterString,"",GeneticCodeExclusions);
        
		LikelihoodFunction siteLikelihoodLoc = (locSiteFilter, perBranchTree);
        LFCompute (siteLikelihoodLoc,LF_START_COMPUTE);
        LFCompute (siteLikelihoodLoc,baseline);
        
        _totalBranchCount  = Columns (_bn) - 1;
        posteriorEstimates = {};
        
        if (_mixingP != 1 && _mixingP != 0) {
            _priorOdds = (1-_mixingP)/_mixingP;
        } else {
            _priorOdds = 0;
        }
        
        for (k = 0; k < _totalBranchCount; k+=1)
        {
             _pname = "perBranchTree." + _bn[k] + ".lmp";
            ExecuteCommands ("`_pname`=1");
            LFCompute (siteLikelihoodLoc,LOGL0);
            
            MaxL     = -Max (LOGL0,baseline);
            
            baseline += MaxL;
            LOGL0 = Exp(MaxL+LOGL0);
            LOGL1 = (Exp(baseline) - _mixingP * LOGL0) / (1-_mixingP);
            
            ExecuteCommands ("`_pname`=_mixingP");
            _posteriorProb = {{LOGL0 * _mixingP, LOGL1 * (1-_mixingP)}};
            _posteriorProb = _posteriorProb * (1/(+_posteriorProb));
            if ( _priorOdds != 0) {
                eBF = _posteriorProb[1] / (1 - _posteriorProb[1]) / _priorOdds;
            } else {
                eBF = 1;
            }
            posteriorEstimates [_bn[k]] = {1,4};
            (posteriorEstimates [_bn[k]])[0] = _posteriorProb[1];
            (posteriorEstimates [_bn[k]])[1] = eBF;
            (posteriorEstimates [_bn[k]])[2] = (subMap[_bn[k]])[0];
            (posteriorEstimates [_bn[k]])[3] = (subMap[_bn[k]])[1];
            baseline += -MaxL;
        } 
        LFCompute (siteLikelihoodLoc,LF_DONE_COMPUTE);
        
        return posteriorEstimates;
}
