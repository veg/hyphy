global		R;
global     AC;
global 	   AT;
global     CG;
global	   CT;
global     GT;

R  = 1;
AT = 1;
CG = 1;
CT = 1;
GT = 1;
AC = 1;

ModelMatrixDimension = 0;

/*--------------------------------------------------------------------------------------*/

function testLRT (vec1, vec2)
{
	size1 = Columns(vec1);
	
	sumVec1 = {size1,1};
	jvec	= {2,size1};
	resMx1	= {khSamples,1};
	resMx2	= {khSamples,1};
	
	for (k=0; k<size1; k=k+1)
	{
		sumVec1 [k]	   = 1;
		jvec	[0][k] = Log(vec1[k]);
		jvec	[1][k] = Log(vec2[k]);
	}
	
	
	for (k=0; k<khSamples; k=k+1)
	{
		resampled = Random(jvec,1);
		resampled = resampled*sumVec1;
		resMx1[k] = resampled[0];
		resMx2[k] = resampled[1];
	}
	
	return (resMx1-resMx2)*2;
}


/*--------------------------------------------------------------------------------------------------------------------------------------*/

function BuildCodonFrequencies (obsF)
{
	if (!ModelMatrixDimension)
	{
		ModelMatrixDimension = 64;
		for (h = 0 ; h<64; h=h+1)
		{
			if (_Genetic_Code[h]==10)
			{
				ModelMatrixDimension = ModelMatrixDimension-1;
			}
		}
	}
	
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

/*--------------------------------------------------------------------------------------------------------------------------------------*/

function BuildCodonFrequencies2 (obsF)
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
			PIStop = PIStop-obsF[first]*obsF[second]*obsF[third];
			continue; 
		}
		result[h-hshift][0]=obsF[first]*obsF[second]*obsF[third];
	}
	return result*(1.0/PIStop);
}

/*--------------------------------------------------------------------------------------------------------------------------------------*/

nucModelString = "nucModelMatrix = {{*,AC*t,t,AT*t}{AC*t,*,CG*t,CT*t}{t,CG*t,*,GT*t}{AT*t,CT*t,GT*t,*}};";

/*--------------------------------------------------------------------------------------------------------------------------------------*/

function PopulateModelMatrix (ModelMatrixName&, EFV, modelType)
{
	
	ModelMatrixName = {ModelMatrixDimension,ModelMatrixDimension}; 
	hshift = 0;

	if (modelType == 0)
	{
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
					if (transition<transition2)
					{
						trSM = transition;
						trLG = transition2;
					}
					else
					{
						trSM = transition2;
						trLG = transition;
					}
					
					if (trSM==0)
					{
						if (trLG==1)
						{
							if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
							{
								ModelMatrixName[h-hshift][v-vshift] := AC*synRate*EFV__[transition__][nucPosInCodon__];
								ModelMatrixName[v-vshift][h-hshift] := AC*synRate*EFV__[transition2__][nucPosInCodon__];
							}
							else
							{
								ModelMatrixName[h-hshift][v-vshift] := AC*R*synRate*EFV__[transition__][nucPosInCodon__];
								ModelMatrixName[v-vshift][h-hshift] := AC*R*synRate*EFV__[transition2__][nucPosInCodon__];
							}
						}
						else
						{
							if (trLG==2)
							{
								if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
								{
									ModelMatrixName[h-hshift][v-vshift] := synRate*EFV__[transition__][nucPosInCodon__];
									ModelMatrixName[v-vshift][h-hshift] := synRate*EFV__[transition2__][nucPosInCodon__];
								}
								else
								{
									ModelMatrixName[h-hshift][v-vshift] := R*synRate*EFV__[transition__][nucPosInCodon__];
									ModelMatrixName[v-vshift][h-hshift] := R*synRate*EFV__[transition2__][nucPosInCodon__];
								}							
							}
							else
							{
								if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
								{
									ModelMatrixName[h-hshift][v-vshift] := AT*synRate*EFV__[transition__][nucPosInCodon__];
									ModelMatrixName[v-vshift][h-hshift] := AT*synRate*EFV__[transition2__][nucPosInCodon__];
								}
								else
								{
									ModelMatrixName[h-hshift][v-vshift] := AT*R*synRate*EFV__[transition__][nucPosInCodon__];
									ModelMatrixName[v-vshift][h-hshift] := AT*R*synRate*EFV__[transition2__][nucPosInCodon__];
								}							
							}
						}
					}
					else
					{
						if (trSM==1)
						{
							if (trLG==2)
							{
								if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
								{
									ModelMatrixName[h-hshift][v-vshift] := CG*synRate*EFV__[transition__][nucPosInCodon__];
									ModelMatrixName[v-vshift][h-hshift] := CG*synRate*EFV__[transition2__][nucPosInCodon__];
								}
								else
								{
									ModelMatrixName[h-hshift][v-vshift] := CG*R*synRate*EFV__[transition__][nucPosInCodon__];
									ModelMatrixName[v-vshift][h-hshift] := CG*R*synRate*EFV__[transition2__][nucPosInCodon__];
								}
							}
							else
							{
								if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
								{
									ModelMatrixName[h-hshift][v-vshift] := CT*synRate*EFV__[transition__][nucPosInCodon__];
									ModelMatrixName[v-vshift][h-hshift] := CT*synRate*EFV__[transition2__][nucPosInCodon__];
								}
								else
								{
									ModelMatrixName[h-hshift][v-vshift] := CT*R*synRate*EFV__[transition__][nucPosInCodon__];
									ModelMatrixName[v-vshift][h-hshift] := CT*R*synRate*EFV__[transition2__][nucPosInCodon__];
								}							
							}
						}
						else
						{
							if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
							{
								ModelMatrixName[h-hshift][v-vshift] := GT*synRate*EFV__[transition__][nucPosInCodon__];
								ModelMatrixName[v-vshift][h-hshift] := GT*synRate*EFV__[transition2__][nucPosInCodon__];
							}
							else
							{
								ModelMatrixName[h-hshift][v-vshift] := GT*R*synRate*EFV__[transition__][nucPosInCodon__];
								ModelMatrixName[v-vshift][h-hshift] := GT*R*synRate*EFV__[transition2__][nucPosInCodon__];
							}							
						}
					}
				}
		     }
	    }		
	}
	else
	{
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
					if (transition<transition2)
					{
						trSM = transition;
						trLG = transition2;
					}
					else
					{
						trSM = transition2;
						trLG = transition;
					}
					
					if (trSM==0)
					{
						if (trLG==1)
						{
							if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
							{
								ModelMatrixName[h-hshift][v-vshift] := AC*synRate;
								ModelMatrixName[v-vshift][h-hshift] := AC*synRate;
							}
							else
							{
								ModelMatrixName[h-hshift][v-vshift] := AC*R*synRate;
								ModelMatrixName[v-vshift][h-hshift] := AC*R*synRate;
							}
						}
						else
						{
							if (trLG==2)
							{
								if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
								{
									ModelMatrixName[h-hshift][v-vshift] := synRate;
									ModelMatrixName[v-vshift][h-hshift] := synRate;
								}
								else
								{
									ModelMatrixName[h-hshift][v-vshift] := R*synRate;
									ModelMatrixName[v-vshift][h-hshift] := R*synRate;
								}							
							}
							else
							{
								if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
								{
									ModelMatrixName[h-hshift][v-vshift] := AT*synRate;
									ModelMatrixName[v-vshift][h-hshift] := AT*synRate;
								}
								else
								{
									ModelMatrixName[h-hshift][v-vshift] := AT*R*synRate;
									ModelMatrixName[v-vshift][h-hshift] := AT*R*synRate;
								}							
							}
						}
					}
					else
					{
						if (trSM==1)
						{
							if (trLG==2)
							{
								if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
								{
									ModelMatrixName[h-hshift][v-vshift] := CG*synRate;
									ModelMatrixName[v-vshift][h-hshift] := CG*synRate;
								}
								else
								{
									ModelMatrixName[h-hshift][v-vshift] := CG*R*synRate;
									ModelMatrixName[v-vshift][h-hshift] := CG*R*synRate;
								}
							}
							else
							{
								if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
								{
									ModelMatrixName[h-hshift][v-vshift] := CT*synRate;
									ModelMatrixName[v-vshift][h-hshift] := CT*synRate;
								}
								else
								{
									ModelMatrixName[h-hshift][v-vshift] := CT*R*synRate;
									ModelMatrixName[v-vshift][h-hshift] := CT*R*synRate;
								}							
							}
						}
						else
						{
							if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
							{
								ModelMatrixName[h-hshift][v-vshift] := GT*synRate;
								ModelMatrixName[v-vshift][h-hshift] := GT*synRate;
							}
							else
							{
								ModelMatrixName[h-hshift][v-vshift] := GT*R*synRate;
								ModelMatrixName[v-vshift][h-hshift] := GT*R*synRate;
							}							
						}
					}
				}
		     }
	    }		
		return 1;
	}
	return 0;
}


#include "TemplateModels/chooseGeneticCode.def";

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
ModelTitle = modelDesc[0];
			
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

if (Abs(modelConstraintString))
{
	ExecuteCommands (modelConstraintString);
}

khSamples = -1;
while (khSamples < 0)
{
	fprintf (stdout,"\nHow many samples should be generated for the KH test (0 to skip test):");
	fscanf  (stdin,"Number", khSamples);
}			

ChoiceList (branchLengths,"Branch Lengths",1,SKIP_NONE,
			"Codon Model","Jointly optimize rate parameters and branch lengths (slow and thorough)",
			"Nucleotide Model","Estimate branch lengths once, using an appropriate nucleotide model (quick and dirty)."
		    );

if (branchLengths<0)
{
	return;
}

SetDialogPrompt ("Please specify a codon data file:");

DataSet					 ds = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter filteredData  = CreateFilter (ds,3,"","",GeneticCodeExclusions);
HarvestFrequencies  (observedFreq ,filteredData,3,1,1);
HarvestFrequencies  (observedFreqE,filteredData,1,1,1);

vectorOfFrequencies 	= BuildCodonFrequencies  (observedFreq);
vectorOfFrequencies1x4  = BuildCodonFrequencies2 (observedFreqE);
HarvestFrequencies  (empiricalCodon,filteredData,3,3,0);

cfm = {ModelMatrixDimension, 3};

for (k=0; k<ModelMatrixDimension; k=k+1)
{
	cfm [k][0] = vectorOfFrequencies[k];
	cfm [k][1] = vectorOfFrequencies1x4[k];
	cfm [k][2] = empiricalCodon[k];
}

cL 	   = "Codon";
labels = {{"F3x4","F1x4", "Empirical",""}};

nucs = {{"A","C","G","T"}};

k=0;
for (p1=0; p1<4; p1=p1+1)
{
	for (p2=0; p2<4; p2=p2+1)
	{
		for (p3=0; p3<4; p3=p3+1)
		{
			if (_Genetic_Code[k] != 10)
			{
				cL = cL + ";" + nucs[p1] + nucs[p2] + nucs[p3];
			}
			k=k+1;
		}
	}
}

labels[3] = cL;

OpenWindow (CHARTWINDOW,
	{{"Codon Frequencies"}
	{"labels"}
	{"cfm"}
	{"Bar Chart"}
	{"Index"}
	{labels[0]+";"+labels[1]+";"+labels[2]}
	{"Codon Position"}
	{""}
	{"Proportion"}
	{"3"}
	{""}
	{""}
	{"10;1.309;0.785398"}
	{"Times:12:0;Times:10:0;Times:14:2"}
	{"0;0;16777215;5000268;0;0;6750054;11842740;13158600;14474460;0;3947580;13421772;8388608;32768;4194432;1644825;7018159;1460610;16748822;11184810;14173291;14173291"}
	{"16"}
	},
	"(SCREEN_WIDTH-30);(SCREEN_HEIGHT-80)/2;20;SCREEN_HEIGHT/2");
	

fm = Transpose (observedFreq);
labels = {{"A","C", "G","T","Position;1;2;3"}};
OpenWindow (CHARTWINDOW,
	{{"Nucleotide Frequencies"}
	{"labels"}
	{"fm"}
	{"Bar Chart"}
	{"Index"}
	{"A;C;G;T"}
	{"Codon Position"}
	{""}
	{"Proportion"}
	{"3"}
	{""}
	{""}
	{"10;1.309;0.785398"}
	{"Times:12:0;Times:10:0;Times:14:2"}
	{"0;0;16777215;5000268;0;0;6750054;11842740;13158600;14474460;0;3947580;13421772;8388608;32768;4194432;1644825;7018159;1460610;16748822;11184810;14173291;14173291"}
	{"16"}
	},
	"(SCREEN_WIDTH-30)/2;(SCREEN_HEIGHT-80)/2;20;40");

fprintf (stdout,"\n0). Read the data\n",ds);
MG94Matrix = 0;
PopulateModelMatrix ("MG94Matrix", observedFreq, 0);
Model MG94Model = (MG94Matrix, vectorOfFrequencies, 0);

_DO_TREE_REBALANCE_ = 1;

#include "queryTree.bf";

if (branchLengths)
{
	DataSetFilter nucFilter = CreateFilter (filteredData,1);
	ExecuteCommands (nucModelString+"\nModel nucModel = (nucModelMatrix,observedFreqE);");

	Tree  nucTree = treeString;
	LikelihoodFunction nuc_lf = (nucFilter,nucTree);
	Optimize (nuc_res, nuc_lf);
	global codonFactor = 0.33;
	
}

GY94Matrix = 0;
PopulateModelMatrix ("GY94Matrix", observedFreq, 0);
Model GY94Model = (GY94Matrix, vectorOfFrequencies, 1);

Tree givenTreeGY = treeString;

LikelihoodFunction lf_MG = (filteredData,givenTree);

if (branchLengths)
{
	ClearConstraints (givenTree);
	ReplicateConstraint ("this1.?.synRate:=this2.?.t__/codonFactor",givenTree,nucTree);
}

Optimize (resMG,lf_MG);
ConstructCategoryMatrix (MGconditionals,lf_MG,COMPLETE);


fprintf (stdout,"\n1). MG94x", ModelTitle, " fit\n", lf_MG, "\n");

COVARIANCE_PARAMETER = "R";
COVARIANCE_PRECISION = 0.95;
CovarianceMatrix (MGR, lf_MG);

LikelihoodFunction lf_GY = (filteredData,givenTreeGY);

if (branchLengths)
{
	ClearConstraints (givenTree);
	ReplicateConstraint ("this1.?.synRate:=this2.?.t__/codonFactor",givenTreeGY,nucTree);
}

Optimize (resGY,lf_GY);
ConstructCategoryMatrix (GYconditionals,lf_GY,COMPLETE);
CovarianceMatrix (GYR, lf_GY);

fprintf (stdout,"\n2). GY94x", ModelTitle, " fit\n", lf_GY, "\n");

fprintf (stdout, "\n3). Profile likelihood intervals on dN/dS\n",	
				  "\tMG94: ", Format(MGR[0],10,5), " <= dN/dS (", Format(MGR[1],10,5), ") <= ", Format(MGR[2],10,5),"\n",
				  "\tGY94: ", Format(GYR[0],10,5), " <= dN/dS (", Format(GYR[1],10,5), ") <= ", Format(GYR[2],10,5),"\n");

siteWiseLikelihoods	= {filteredData.sites,3};

MGbetter = 0;

for (k=0; k < filteredData.sites; k=k+1)
{
	 siteWiseLikelihoods [k][0] = Log (MGconditionals[k]);
	 siteWiseLikelihoods [k][1] = Log (GYconditionals[k]);
	 siteWiseLikelihoods [k][2] = siteWiseLikelihoods [k][0]-siteWiseLikelihoods [k][1];
	 if (siteWiseLikelihoods [k][2]>=0)
	 {
	 	MGbetter = MGbetter + 1;
	 }
}

labels = {{"MG94 Log L","GY94 Log L", "MG94-GY94 Log L difference"}};
OpenWindow (CHARTWINDOW,
	{{"Site by site comparison"}
	{"labels"}
	{"siteWiseLikelihoods"}
	{"Bar Chart"}
	{"Index"}
	{labels[2]}
	{"Codon"}
	{"LRT"}
	{"Log L difference"}
	{""}
	{""}
	{""}
	{"10;1.309;0.785398"}
	{"Times:12:0;Times:10:0;Times:14:2"}
	{"0;0;16777215;5000268;0;0;6750054;11842740;13158600;14474460;0;3947580;13421772;5000268;11776947;10066329;9199669;7018159;1460610;16748822;11184810;14173291;14173291"}
	{"16"}
	},
	"(SCREEN_WIDTH-30)/2;(SCREEN_HEIGHT-80)/2;20+(SCREEN_WIDTH)/2;40");
}

fprintf (stdout,"\n4). Site-by-site fits\n\tMG94 is better on ", MGbetter, "/", filteredData.sites, " codon columns\n");

if (khSamples)
{
	fprintf (stdout,"\n5). Kishino-Hasegawa resampling test\n");

	LRT1 = 2*(resMG[1][0]-resGY[1][0]);

	if (LRT1>0)
	{
		fprintf (stdout, "\n\nSUMMARY:\n\n\tMG94 is better than GY94. LRT = ",LRT1);
		test1 = testLRT (MGconditionals,GYconditionals)%0;
	}
	else
	{
		fprintf (stdout, "\n\nSUMMARY:\n\n\tGY94 is better than MG94. LRT = ",LRT1);
		test1 = testLRT (GYconditionals,MGconditionals)%0;
	}

	for (k=0; k<Columns(GYconditionals); k=k+1)
	{
		if (test1[k]>0)
		{
			break;
		}
	}	
	fprintf (stdout, "\n\nEstimated p-value for  (based on ", khSamples," replicates): ", Format(k/khSamples,10,4),"\n");

}
