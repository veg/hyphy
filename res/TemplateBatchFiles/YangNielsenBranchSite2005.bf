/* 1. include a file to define the genetic code
   Note the use of base directory and path forming variables to make this analysis 
   independent of directory placement
 */

incFileName = HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"chooseGeneticCode.def";
ExecuteCommands  ("#include \""+incFileName+"\";");

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
										  

/* 5. prompt for the type of model to run */

ChoiceList (modelKind, "Model", 1, SKIP_NONE,
					   "Alternative", "Fit Model A with 4 rate classes. Class 1: Negative selection in FG and BG. Class 2: Neutral evolution in FG and BG. Class 3: Negative selection in BG, Positive in FG. Class 4: Neutral evolution in BG, Positive in FG",
					   "Null for Test 1", "Fit a model with 2 rate classes. Class 1: Negative selection in FG and BG. Class 2: Neutral evolution in FG and BG.",
					   "Null for Test 2", "Fit Model A with 3 rate classes. Class 1: Negative selection in FG and BG. Class 2: Neutral evolution in FG and BG. Class 3: Negative selection in BG, Neutral in FG.",
);					   

if (modelKind < 0)
{
	return 0;
}


/* 6. define the 'site_kind' variable as a discrete category variable; 
	  it decides which class a site belongs to, but does not
	  determine omega ratios directly (see below for this) */

global P_0     = 0.5;
P_0:<1;
P_0:>0;

if (modelKind == 1)
{
	rateClasses = 2;

	categFreqMatrix = {{P_0,1-P_0}} ;
	categRateMatrix = {{1,2}};
}
else
{
	P_0        = 1/4;
	global 		P_1_aux = 1/4;
	global 		P_1     := Min(P_1_aux,1-P_0);
	P_1:<1;
	P_1:>0;
	
	if (modelKind == 0)
	{
		rateClasses = 4;
		categFreqMatrix = {{P_0,P_1,(1-P_0-P_1)/(P_0+P_1)*P_0,(1-P_0-P_1)/(P_0+P_1)*P_1}} ;
		categRateMatrix = {{1,2,3,4}};
	}
	else
	{
		rateClasses = 3;
		categFreqMatrix = {{P_0,P_1/(P_0+P_1),(1-P_0-P_1)/(P_0+P_1)*P_0}} ;
		categRateMatrix = {{1,2,3}};
	}
}

category site_kind = (rateClasses, categFreqMatrix , MEAN, ,categRateMatrix, 1, 4);

/* 7. define the GY94 rate matrix; for now each branch will have it's own
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

Model GY_Model = (GY_Matrix,codonFreqs,1);

/*10. Define the tree and pick the foreground branch, displaying a tree window to facilitate selection; 
the latter step is executed for 2 of 3 model choices */

Tree 	  givenTree = treeString;

USE_LAST_RESULTS    = 0;
OPTIMIZATION_METHOD = 4;

/* Approximate kappa and branch lengths using an HKY85 fit */

HKY85_Matrix = {{*,t*kappa_inv,t,t*kappa_inv}
				{t*kappa_inv,*,kappa_inv*t,t}
				{t,t*kappa_inv,*,kappa_inv*t}
				{t*kappa_inv,t,kappa_inv*t,*}};
			
HarvestFrequencies (nucFreqs,ds,1,1,1);
Model HKY85_Model = (HKY85_Matrix,nucFreqs,1);

Tree		  nucTree = treeString;
DataSetFilter nucData = CreateFilter (ds,1);

fprintf (stdout, "Obtaining nucleotide branch lengths and kappa to be used as starting values...\n");
LikelihoodFunction	nuc_lf = (nucData,nucTree);
Optimize(nuc_mle,nuc_lf);
fprintf (stdout, "\n", Format (nucTree,1,1), "\nkappa=", Format (1/kappa_inv,8,3), "\n");

USE_LAST_RESULTS = 1;

if (modelKind != 1)
{

	mxTreeSpec = {5,1};
	mxTreeSpec [0] = "givenTree";
	mxTreeSpec [1] = "8240";
	mxTreeSpec [2] = "10,40,-10,-175,1";
	mxTreeSpec [3] = "";
	mxTreeSpec [4] = "";
	OpenWindow 		(TREEWINDOW, mxTreeSpec,"(SCREEN_WIDTH-50)/2;(SCREEN_HEIGHT-50)/2;30+(SCREEN_WIDTH-30)/2;45");				

	leafNodes	  = TipCount   (givenTree);			
	internalNodes = BranchCount(givenTree);

	choiceMatrix = {internalNodes+leafNodes,2};
	for (bc=0; bc<internalNodes; bc=bc+1)
	{
		choiceMatrix[bc][0] = BranchName(givenTree,bc);
		choiceMatrix[bc][1] = "Internal Branch Rooting " + givenTree[bc];
	}
	for (bc=0; bc<leafNodes; bc=bc+1)
	{
		choiceMatrix[bc+internalNodes][0] = TipName(givenTree,bc);
		choiceMatrix[bc+internalNodes][1] = "Leaf node " + choiceMatrix[bc+internalNodes][0];
	}

	ChoiceList  (stOption,"Choose the foreground branch",0,NO_SKIP,choiceMatrix);

	if (stOption[0] < 0)
	{
		return -1;
	}

	fprintf (stdout, "\n\n", Columns (stOption)," foreground branch(es) set to: ", "\n");
	for (bc = 0; bc < Columns (stOption); bc = bc + 1)
	{
		fprintf (stdout, choiceMatrix[stOption[bc]][0], "\n");
	}
	OpenWindow (CLOSEWINDOW, "Tree givenTree");
}
	
/* 15. Constrain dS and dN in the tree to based upon different models */

global omega_0 = 0.25;
omega_0 :< 1;
ClearConstraints (givenTree);

if (modelKind == 1)
{
	global omega := ((site_kind==1)*omega_0+(site_kind==2));
	/*  will evaluate to omega_0 for sites in class site_kind=1 and to 1 for sites in class site_kind=2 */
	ReplicateConstraint ("this1.?.nonSynRate:=omega*this2.?.synRate",givenTree,givenTree);
}
else
{
	if (modelKind == 2)
	{
		global omega_FG := ((site_kind==1)*omega_0+(site_kind>1)); 						/* foreground model */
		global omega_BG := (((site_kind==1)+(site_kind==3))*omega_0+(site_kind==2));    /* background model */
	}
	else
	{
		global omega_2 = 2.0;
		omega_2:>1;
		global omega_FG := ((site_kind==1)*omega_0+(site_kind==2)+(site_kind>2)*omega_2); 			/* foreground model */
		global omega_BG := (((site_kind==1)+(site_kind==3))*omega_0+(site_kind==2)+(site_kind==4)); /* background model */
	}
	
	/* constrain the foreground branch first */
	for (bc = 0; bc < Columns (stOption); bc = bc + 1)
	{
		ExecuteCommands ("givenTree."+choiceMatrix[stOption[bc]][0]+".nonSynRate:=omega_FG*givenTree."+choiceMatrix[stOption[bc]][0]+".synRate;");
	}

	/* constrain other branches next */
	ReplicateConstraint ("this1.?.nonSynRate:=omega_BG*this2.?.synRate",givenTree,givenTree);

}

/* 16. define and optimize the likelihood function */

bNames = BranchName   (givenTree,-1);
nucBL  = BranchLength (nucTree,-1);

for (bc=0; bc<Columns(bNames)-1; bc=bc+1)
{
	ExecuteCommands ("givenTree."+bNames[bc]+".synRate=nucTree."+bNames[bc]+".t;");
}

codBL  = BranchLength (givenTree,-1);

for (bc=0; bc<Columns(bNames)-1; bc=bc+1)
{
	if (nucBL[bc]>0)
	{
		ExecuteCommands ("givenTree."+bNames[bc]+".synRate=nucTree."+bNames[bc]+".t*"+nucBL[bc]/codBL[bc]+";");
	}
}

OPTIMIZATION_PRECISION = 0.001;
LikelihoodFunction lf = (filteredData, givenTree);


while (1)
{
	Optimize 		   (mles,lf);
	
	fprintf (stdout, lf);

	GetString 		  (lfParameters, lf, -1);
	glV 		 	= lfParameters["Local Independent"];
	stashedValues 	= {};
	for (glVI = 0; glVI < Columns (glV); glVI = glVI + 1)
	{
		ExecuteCommands ("stashedValues[\""+glV[glVI]+"\"] = " + glV[glVI] + ";\n");
	}
	glV 		 	= lfParameters["Global Independent"];
	for (glVI = 0; glVI < Columns (glV); glVI = glVI + 1)
	{
		ExecuteCommands ("stashedValues[\""+glV[glVI]+"\"] = " + glV[glVI] + ";\n");
	}
	
	mlBL		 	= BranchLength (givenTree,-1);
	
	samples = 500;
	fprintf (stdout, "\nChecking for convergence by Latin Hypercube Sampling (this may take a bit of time...)\n");
		
	steps 		  = 50;
	if (modelKind == 0)
	{
		vn			  = {{"P_0","P_1_aux","omega_0", "omega_2"}};
		ranges		  = {{0.0001,1}{0.0001,1}{0.0001,1}{1,10}};
	}
	else
	{
		if (modelKind==2)
		{
			vn			  = {{"P_0","P_1_aux","omega_0"}};
			ranges		  = {{0.0001,1}{0.0001,1}{0.0001,1}};
		}
		else
		{
			if (modelKind==1)
			{
				vn			  = {{"P_0","omega_0"}};
				ranges		  = {{0.0001,1}{0.0001,1}};
			}		
		}
	}
		
	LFCompute (lf,LF_START_COMPUTE);
	
	for (sample = 0; sample < samples; sample = sample + 1)
	{
		rv = Random({1,steps}["_MATRIX_ELEMENT_COLUMN_"],0);
		for (vid = 0; vid < Columns (vn); vid = vid + 1)
		{
			ctx = vn[vid] + "=" + (ranges[vid][0] + (ranges[vid][1]-ranges[vid][0])/steps*rv[vid]);
			ExecuteCommands (ctx);
		}
		currentBL = BranchLength (givenTree,-1);
		for (bc=0; bc<Columns(bNames)-1; bc=bc+1)
		{
			if (currentBL[bc]>0)
			{
				ExecuteCommands ("givenTree."+bNames[bc]+".synRate=givenTree."+bNames[bc]+".synRate*"+mlBL[bc]/currentBL[bc]+";");
			}
		}		

		LFCompute (lf,sample_value);	
		if (sample_value>mles[1][0])
		{
			fprintf (stdout, "\nFound a better likelihood score. Restarting the optimization routine.\n");
			break;
		}	
	}
	LFCompute (lf,LF_DONE_COMPUTE);
	
	if (sample < samples)
	{
		continue;
	}
	storedV = Rows (stashedValues);
	for (k=0; k<Columns (storedV); k=k+1)
	{
		ExecuteCommands (storedV[k] + "=" + stashedValues[storedV[k]]);
	}
	fprintf (stdout, "\nThe estimation procedure appears to have converged.\n");
	break;
}


/* 17. Report inferred rate distribition to screen */

if (modelKind == 0)
{
	fprintf (stdout, "\nInferred rate distribution:",
					 "\n\tClass 0.  omega_0 = ", Format (omega_0, 5,3), " weight = ", Format (P_0,5,3),
					 "\n\tClass 1.  omega  := ", Format (1, 5,3), " weight = ", Format (P_1,5,3),
					 "\n\tClass 2a. Background omega_0 = ", Format (omega_0, 5,3), " foreground omega_2 = ", Format (omega_2, 5,3), " weight = ", Format (P_0(1-P_0-P_1)/(P_0+P_1),5,3), 
					 "\n\tClass 2b. Background omega  := ", Format (1, 5,3), " foreground omega_2 = ", Format (omega_2, 5,3), " weight = ", Format (P_1(1-P_0-P_1)/(P_0+P_1),5,3), "\n");
}
if (modelKind == 1)
{
	fprintf (stdout, "\nInferred rate distribution:",
					 "\n\tClass 0.  omega_0 = ", Format (omega_0, 5,3), " weight = ", Format (P_0,5,3),
					 "\n\tClass 1.  omega  := ", Format (1, 5,3), " weight = ", Format (1-P_0,5,3), "\n");
}
if (modelKind == 2)
{
	fprintf (stdout, "\nInferred rate distribution:",
					 "\n\tClass 0.  omega_0 = ", Format (omega_0, 5,3), " weight = ", Format (P_0,5,3),
					 "\n\tClass 1.  omega  := ", Format (1, 5,3), " weight = ", Format (P_1,5,3),
					 "\n\tClass 2a. Background omega_0 = ", Format (omega_0, 5,3), " foreground omega_2 := ", Format (1, 5,3), " weight = ", Format (P_0(1-P_0-P_1)/(P_0+P_1),5,3), 
					 "\n\tClass 2b. Background omega  := ", Format (1, 5,3), " foreground omega_2 := ", Format (1, 5,3), " weight = ", Format (P_1(1-P_0-P_1)/(P_0+P_1),5,3), "\n");
}

/* 18. Prepare and open a window of conditional probabilities at every site (this requires a GUI but will still run - 
- just not open any windows in a console build */

ConstructCategoryMatrix (posteriorMatrix, lf, COMPLETE);
posteriorMatrix = Transpose (posteriorMatrix);

GetInformation			(siteProfile, site_kind); 
						/* this call returns the distribution function {{value1,..., valueN}{prob1, ..., probN}} for site_kind */

headers = {1, rateClasses};

for (k=1; k<=rateClasses;k=k+1)
{
	headers [k-1] = "Class " + k;
}

/* convert distribution info to the form expected by OpenWindow */
disributionInfo = "site_kind";

for (k = 0; k < rateClasses; k=k+1)
{
	disributionInfo = disributionInfo + ":" + siteProfile[1][k];
}
for (k = 0; k < rateClasses; k=k+1)
{
	disributionInfo = disributionInfo + ":" + siteProfile[0][k];
}

OpenWindow (DISTRIBUTIONWINDOW,{{"Conditional probabilities by site"}
		{"headers"}
		{"posteriorMatrix"}
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
		{disributionInfo}
		},
		"600;600;50;50");
		
_MARGINAL_MATRIX_ = Transpose(posteriorMatrix);
_CATEGORY_VARIABLE_CDF_ = siteProfile[1][-1];

ExecuteAFile(HYPHY_LIB_DIRECTORY+"ChartAddIns"+DIRECTORY_SEPARATOR+"DistributionAddIns"+DIRECTORY_SEPARATOR+"Includes"+DIRECTORY_SEPARATOR+"posteriors.ibf");
ExecuteAFile(HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Utility"+DIRECTORY_SEPARATOR+"WriteDelimitedFiles.bf");
siteCount = Columns(_MARGINAL_MATRIX_);
siteCounter = {};
for (k=0; k<siteCount; k=k+1)
{
	siteCounter[k] = k + 1;
}
classCount = Columns(_CATEGORY_VARIABLE_CDF_);
columnHeaders = {1,classCount+1};
columnHeaders[0] = "Site";
for (k=1; k<=classCount; k=k+1)
{
	columnHeaders[k] = "Class "+k;
}

SetDialogPrompt ("Write site-by-site conditional probabilities to this file:");
WriteSeparatedTable("",columnHeaders,Transpose(_MARGINAL_MATRIX_), siteCounter,",");
