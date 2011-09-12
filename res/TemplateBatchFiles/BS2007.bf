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

Model GY_Model = (GY_Matrix,codonFreqs,1);

/*10. Define the tree and pick the foreground branch, displaying a tree window to facilitate selection; 
the latter step is executed for 2 of 3 model choices */

Tree 	  givenTree1 	= treeString;
Tree 	  givenTree2	= treeString;
Tree 	  givenTree3	= treeString;
Tree	  givenTree4	= treeString;

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

fprintf 					 (stdout, "Obtaining nucleotide branch lengths and kappa to be used as starting values...\n");
LikelihoodFunction	nuc_lf = (nucData,nucTree);
Optimize					 (nuc_mle,nuc_lf);
fprintf 					 (stdout, "\n", Format (nucTree,1,1), "\nkappa=", Format (1/kappa_inv,8,3), "\n");

USE_LAST_RESULTS = 1;

mxTreeSpec = {5,1};
mxTreeSpec [0] = "nucTree";
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
OpenWindow (CLOSEWINDOW, "Tree nucTree");
	
/* 15. Constrain dS and dN in the tree to based upon different models */

global omega_1 = 0.25;
global omega_2 = 0.5;
global omega_3 = 1.5;

omega_1 		:< 1;
omega_2 		:< 1;
omega_3 		:> 1;

ClearConstraints (givenTree1);
ClearConstraints (givenTree2);
ClearConstraints (givenTree3);
ClearConstraints (givenTree4);

/* 16. define and optimize the likelihood function */

bNames = BranchName   (givenTree,-1);
nucBL  = BranchLength (nucTree,-1);

for (bc=0; bc<Columns(bNames)-1; bc=bc+1)
{
	ExecuteCommands ("givenTree1."+bNames[bc]+".synRate		=nucTree."+bNames[bc]+".t;");
	ExecuteCommands ("givenTree1."+bNames[bc]+".nonSynRate	=nucTree."+bNames[bc]+".t;");
	ExecuteCommands ("givenTree2."+bNames[bc]+".synRate		=nucTree."+bNames[bc]+".t;");
	ExecuteCommands ("givenTree2."+bNames[bc]+".nonSynRate	=nucTree."+bNames[bc]+".t;");
	ExecuteCommands ("givenTree3."+bNames[bc]+".synRate		=nucTree."+bNames[bc]+".t;");
	ExecuteCommands ("givenTree3."+bNames[bc]+".nonSynRate	=nucTree."+bNames[bc]+".t;");
	ExecuteCommands ("givenTree4."+bNames[bc]+".synRate		=nucTree."+bNames[bc]+".t;");
	ExecuteCommands ("givenTree4."+bNames[bc]+".nonSynRate	=nucTree."+bNames[bc]+".t;");
}

codBL  = BranchLength (givenTree1,-1);

for (bc=0; bc<Columns(bNames)-1; bc=bc+1)
{
	if (nucBL[bc]>0)
	{
		scalingFactor = nucBL[bc]/codBL[bc];
		ExecuteCommands ("givenTree1."+bNames[bc]+".synRate=nucTree."+bNames[bc]+".t*"+scalingFactor+";");
		ExecuteCommands ("givenTree1."+bNames[bc]+".nonSynRate=nucTree."+bNames[bc]+".t*"+scalingFactor+";");
		ExecuteCommands ("givenTree2."+bNames[bc]+".synRate=nucTree."+bNames[bc]+".t*"+scalingFactor+";");
		ExecuteCommands ("givenTree2."+bNames[bc]+".nonSynRate=nucTree."+bNames[bc]+".t*"+scalingFactor+";");
		ExecuteCommands ("givenTree3."+bNames[bc]+".synRate=nucTree."+bNames[bc]+".t*"+scalingFactor+";");
		ExecuteCommands ("givenTree3."+bNames[bc]+".nonSynRate=nucTree."+bNames[bc]+".t*"+scalingFactor+";");
		ExecuteCommands ("givenTree4."+bNames[bc]+".synRate		=nucTree."+bNames[bc]+".t;");
		ExecuteCommands ("givenTree4."+bNames[bc]+".nonSynRate	=nucTree."+bNames[bc]+".t;");
	}
}

for (bc = 0; bc < Columns (stOption); bc = bc + 1)
{
	bName = choiceMatrix[stOption[bc]][0];
	ExecuteCommands ("givenTree1."+bName+".nonSynRate:=omega_1*givenTree1."+bName+".synRate");
	ExecuteCommands ("givenTree2."+bName+".nonSynRate:=omega_2*givenTree2."+bName+".synRate");
	ExecuteCommands ("givenTree3."+bName+".nonSynRate:=givenTree3."+bName+".synRate");
	ExecuteCommands ("givenTree4."+bName+".nonSynRate:=omega_3*givenTree4."+bName+".synRate");
}


OPTIMIZATION_PRECISION = 0.001;

global P_1 = 1/4;
global P_2 = 1/3;
global P_3 = 1/2;
P_1 :< 1;
P_2 :< 1;
P_3 :< 1;


fprintf (stdout, "\nFitting the model with selection...\n");

LikelihoodFunction lf  = (filteredData, givenTree1, filteredData, givenTree2,filteredData, givenTree3,filteredData, givenTree4,
						  "Log(P_1*SITE_LIKELIHOOD[0]+(1-P_1)*P_2*SITE_LIKELIHOOD[1]+(1-P_1)(1-P_2)*P_3*SITE_LIKELIHOOD[2]+(1-P_1)(1-P_2)(1-P_3)*SITE_LIKELIHOOD[3])");
						  

Optimize 		    (mles,lf);
fprintf 			(stdout, lf);

	
