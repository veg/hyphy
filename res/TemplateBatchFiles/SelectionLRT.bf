/*---------------------------------------------------------------------------------------------------------------------------------------------*/


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

ModelMatrixDimension = 0;

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function PopulateModelMatrix (ModelMatrixName&, EFV)
{
	if (!ModelMatrixDimension)
	{
		ModelMatrixDimension = 64;
		for (h = 0 ;h<64; h=h+1)
		{
			if (_Genetic_Code[h]==10)
			{
				ModelMatrixDimension = ModelMatrixDimension-1;
			}
		}
	}
	
	ModelMatrixName = {ModelMatrixDimension,ModelMatrixDimension}; 

	hshift = 0;
	
	modelDefString = "";
	modelDefString*16384;
	
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
				hs = Format(h-hshift,0,0);
				vs = Format(v-vshift,0,0);
				ts = Format(transition,0,0);
				ts2= Format(transition2,0,0);
				ps = Format(nucPosInCodon,0,0);
				aa1 = _Genetic_Code[0][h];
				aa2 = _Genetic_Code[0][v];
				if (aa1==aa2) 
				{
					modelDefString*("ModelMatrixName["+hs+"]["+vs+"] := "+_nucBiasTerms[transition][transition2]+"synRate*EFV__["+ts+"]["+ps+"];\n"+
													 "ModelMatrixName["+vs+"]["+hs+"] := "+_nucBiasTerms[transition][transition2]+"synRate*EFV__["+ts2+"]["+ps+"];\n");
				}
				else
				{
					modelDefString*("ModelMatrixName["+hs+"]["+vs+"] := "+_nucBiasTerms[transition][transition2]+"nonSynRate*EFV__["+ts+"]["+ps+"];\n"+
													 "ModelMatrixName["+vs+"]["+hs+"] := "+_nucBiasTerms[transition][transition2]+"nonSynRate*EFV__["+ts2+"]["+ps+"];\n");	
				}
			}
	    }
    }		
	modelDefString*0;
	ExecuteCommands (modelDefString);
	return 0;
}




incFileName = HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"chooseGeneticCode.def";
ExecuteCommands  ("#include \""+incFileName+"\";");

SetDialogPrompt 	  		  ("Locate a codon file:");
DataSet 		ds  		 = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter	filteredData = CreateFilter(ds,3,"","",GeneticCodeExclusions);



/*---------------------------------------------------------------------------------------------------------------------------------------------*/

global AC = 1;
global AT = 1;
global CG = 1;
global CT = 1;
global GT = 1;

catCounterAL = {totalBranches,1};

MGCustomRateBiasTerms = {{"AC*","","AT*","CG*","CT*","GT*"}};	

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

		
paramCount	  = 0;
_nucBiasTerms = {4,4};
_nucBiasTerms[0][0] = "";


if (modelDesc[0]==modelDesc[1])
{
	MGCustomRateBiasTerms[0] = MGCustomRateBiasTerms[1];
}

_nucBiasTerms[1][0] = MGCustomRateBiasTerms[0];
_nucBiasTerms[0][1] = MGCustomRateBiasTerms[0];
_nucBiasTerms[2][0] = MGCustomRateBiasTerms[1];
_nucBiasTerms[0][2] = MGCustomRateBiasTerms[1];

h = 0;
v = 3;

for (customLoopCounter2=2; customLoopCounter2<6; customLoopCounter2=customLoopCounter2+1)
{
	for (customLoopCounter=0; customLoopCounter<customLoopCounter2; customLoopCounter=customLoopCounter+1)
	{
		if (modelDesc[customLoopCounter]==modelDesc[customLoopCounter2])
		{
			_nucBiasTerms[h][v] = MGCustomRateBiasTerms[customLoopCounter];
			_nucBiasTerms[v][h] = MGCustomRateBiasTerms[customLoopCounter];
			break;
		}
	}
	if (customLoopCounter == customLoopCounter2)
	{
		_nucBiasTerms[h][v] = MGCustomRateBiasTerms[customLoopCounter2];
		_nucBiasTerms[v][h] = MGCustomRateBiasTerms[customLoopCounter2];
	}
	
	v = v+1;
	if (v==4)
	{
		h=h+1;
		v=h+1;
	}
}

HarvestFrequencies (observedFreq,filteredData,3,1,1);

MG94plain = 0;
MULTIPLY_BY_FREQS 				= PopulateModelMatrix ("MG94", observedFreq);
vectorOfFrequencies 			= BuildCodonFrequencies (observedFreq);
Model MG94_Model		 		= (MG94,vectorOfFrequencies,0);

global 		  dNdS0 = 1;

_DO_TREE_REBALANCE_ = 0;
incFileName = HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"queryTree.bf";
ExecuteCommands  ("#include \""+incFileName+"\";");

internalNodes = BranchCount(givenTree);

choiceMatrix = {internalNodes,2};

for (bc=0; bc<internalNodes; bc=bc+1)
{
	choiceMatrix[bc][0] = BranchName(givenTree,bc);
	choiceMatrix[bc][1] = "Internal Branch Rooting " + givenTree[bc];
}


ChoiceList  (bOption,"Choose root of subtree",1,NO_SKIP,choiceMatrix);

if (bOption < 0)
{
	return -1;
}

breakBranch = "givenTree."+BranchName (givenTree,bOption);

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

global dNdS_A 	   = 1;
global dNdS_B 	   = 1;
global dNdS_Branch = 1;

ClearConstraints(givenTree);
ExecuteCommands ("ReplicateConstraint(\"this1.?.nonSynRate:=dNdS_A*this2.?.synRate\","+breakBranch+","+breakBranch+");");
ReplicateConstraint("this1.?.nonSynRate:=dNdS_B*this2.?.synRate",givenTree,givenTree);
ExecuteCommands (breakBranch+".nonSynRate:=dNdS_Branch*"+breakBranch+".synRate;");

dNdS_B:=dNdS_A;
dNdS_Branch:=dNdS_A;

LikelihoodFunction lfNull = (filteredData, givenTree);
fprintf (stdout, "\nPHASE 1: One dN/dS model\n\n");

Optimize (res000, lfNull);

fprintf (stdout, lfNull, "\n");

COVARIANCE_PARAMETER = "dNdS_A";
COVARIANCE_PRECISION = 0.95;
CovarianceMatrix (cmx, lfNull);

fprintf (stdout, "\nError bounds on global dN/dS: ",Format (cmx[0],10,4),"<",Format (cmx[1],10,4),"<",Format (cmx[2],10,4),"\n");


/*---------------------------------------------------------------------------------------------------------------------------------------------*/

dNdS_Branch = dNdS_Branch;

fprintf (stdout, "\nPHASE 2: Separating Branch Versus Two Clades\n\n");
Optimize (res100, lfNull); 
fprintf (stdout, lfNull, "\n");

COVARIANCE_PARAMETER = "dNdS_A";
COVARIANCE_PRECISION = 0.95;
CovarianceMatrix (cmx, lfNull);

fprintf (stdout, "\nError bounds on clades dN/dS: ",Format (cmx[0],10,4),"<",Format (cmx[1],10,4),"<",Format (cmx[2],10,4),"\n");

COVARIANCE_PARAMETER = "dNdS_Branch";
COVARIANCE_PRECISION = 0.95;
CovarianceMatrix (cmx, lfNull);

fprintf (stdout, "Error bounds on branch dN/dS: ",Format (cmx[0],10,4),"<",Format (cmx[1],10,4),"<",Format (cmx[2],10,4),"\n");

p1 = 1-CChi2(2*(res100[1][0]-res000[1][0]),1);
fprintf (stdout, "\nLRT p-value vs the single rate model = ", p1, "\n\n");

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

fprintf (stdout, "\nPHASE 3: Clade A (subtree) and branch vs Clade B\n\n");

dNdS_Branch:=dNdS_A;
dNdS_B = dNdS_B;

Optimize (res110, lfNull); 
fprintf (stdout, lfNull, "\n");


COVARIANCE_PARAMETER = "dNdS_A";
COVARIANCE_PRECISION = 0.95;
CovarianceMatrix (cmx, lfNull);
fprintf (stdout, "\nError bounds on clade A + branch dN/dS: ",Format (cmx[0],10,4),"<",Format (cmx[1],10,4),"<",Format (cmx[2],10,4),"\n");

COVARIANCE_PARAMETER = "dNdS_B";
COVARIANCE_PRECISION = 0.95;
CovarianceMatrix (cmx, lfNull);
fprintf (stdout, "Error bounds on clade B dN/dS: ",Format (cmx[0],10,4),"<",Format (cmx[1],10,4),"<",Format (cmx[2],10,4),"\n");

p1 = 1-CChi2(2*(res110[1][0]-res000[1][0]),1);
fprintf (stdout, "\nLRT p-value vs the single rate model = ", p1, "\n\n");

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

dNdS_Branch:=dNdS_B;

fprintf (stdout, "\nPHASE 4: Clade A (subtree) vs branch and Clade B\n\n");

Optimize (res101, lfNull); 
fprintf (stdout, lfNull, "\n");

COVARIANCE_PARAMETER = "dNdS_A";
COVARIANCE_PRECISION = 0.95;
CovarianceMatrix (cmx, lfNull);
fprintf (stdout, "\nError bounds on clade A dN/dS: ",Format (cmx[0],10,4),"<",Format (cmx[1],10,4),"<",Format (cmx[2],10,4),"\n");

COVARIANCE_PARAMETER = "dNdS_B";
COVARIANCE_PRECISION = 0.95;
CovarianceMatrix (cmx, lfNull);
fprintf (stdout, "Error bounds on clade B + branch dN/dS: ",Format (cmx[0],10,4),"<",Format (cmx[1],10,4),"<",Format (cmx[2],10,4),"\n");

p1 = 1-CChi2(2*(res101[1][0]-res000[1][0]),1);
fprintf (stdout, "\nLRT p-value vs the single rate model = ", p1, "\n\n");

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

dNdS_Branch = dNdS_Branch;

fprintf (stdout, "\nPHASE 5: Clade A, Clade B and branch\n\n");

Optimize (res111, lfNull); 
fprintf (stdout, lfNull, "\n");

COVARIANCE_PARAMETER = "dNdS_A";
COVARIANCE_PRECISION = 0.95;
CovarianceMatrix (cmx, lfNull);

fprintf (stdout, "\nError bounds on clade A (subtree) dN/dS: ",Format (cmx[0],10,4),"<",Format (cmx[1],10,4),"<",Format (cmx[2],10,4),"\n");

COVARIANCE_PARAMETER = "dNdS_B";
COVARIANCE_PRECISION = 0.95;
CovarianceMatrix (cmx, lfNull);

fprintf (stdout, "Error bounds on clade B dN/dS: ",Format (cmx[0],10,4),"<",Format (cmx[1],10,4),"<",Format (cmx[2],10,4),"\n");

COVARIANCE_PARAMETER = "dNdS_Branch";
COVARIANCE_PRECISION = 0.95;
CovarianceMatrix (cmx, lfNull);

fprintf (stdout, "Error bounds on branch dN/dS: ",Format (cmx[0],10,4),"<",Format (cmx[1],10,4),"<",Format (cmx[2],10,4),"\n");

p1 = 1-CChi2(2*(res111[1][0]-res000[1][0]),2);
fprintf (stdout, "\nLRT p-value vs the single rate model = ", p1, "\n\n");

/*---------------------------------------------------------------------------------------------------------------------------------------------*/


fprintf (stdout, "\n+-------+",
				 "\n|SUMMARY|",
				 "\n+-------+\n");

fprintf (stdout, "\nPHASE 1 AIC: ", Format(2*(res000[1][1]-res000[1][0]),20,5));
fprintf (stdout, "\nPHASE 2 AIC: ", Format(2*(res100[1][1]-res100[1][0]),20,5));
fprintf (stdout, "\nPHASE 3 AIC: ", Format(2*(res110[1][1]-res110[1][0]),20,5));
fprintf (stdout, "\nPHASE 4 AIC: ", Format(2*(res101[1][1]-res101[1][0]),20,5));
fprintf (stdout, "\nPHASE 5 AIC: ", Format(2*(res111[1][1]-res111[1][0]),20,5),"\n");


