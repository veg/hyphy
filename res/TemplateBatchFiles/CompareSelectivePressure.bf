SHORT_MPI_RETURN=1;

/*------------------------------------------------------------------------*/

function BuildCodonFrequencies (obsF)
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

/*------------------------------------------------------------------------*/

function MakeModelConstraints (dummy)
{
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
	return 0;
}	


/*------------------------------------------------------------------------*/

function ReportSite2 (siteI, siteM)
{
	fullSites[siteI][0] = doneSites[siteM][0];
	fullSites[siteI][1] = doneSites[siteM][1];
	fullSites[siteI][2] = doneSites[siteM][2];
	fullSites[siteI][3] = doneSites[siteM][3];
	fullSites[siteI][4] = doneSites[siteM][1]/doneSites[siteM][0];
	fullSites[siteI][5] = doneSites[siteM][3]/doneSites[siteM][2];
	fullSites[siteI][6] = doneSites[siteM][4];
	fullSites[siteI][7] = doneSites[siteM][5];
	fullSites[siteI][8] = doneSites[siteM][6];
	fullSites[siteI][9] = doneSites[siteM][7];
	fullSites[siteI][10]= doneSites[siteM][8];


	fprintf (stdout, "| Codon: ", 		Format(siteI+1,4,0),
					 "| Sample 1 (dS,dN,dN/dS): ", 		
					 					Format(fullSites[siteI][0],6,2),
					 					Format(fullSites[siteI][1],6,2),
					 					Format(fullSites[siteI][4],8,2),
					 "| Sample 2 (dS,dN,dN/dS): ", 		
					 					Format(fullSites[siteI][2],6,2),
					 					Format(fullSites[siteI][3],6,2),
					 					Format(fullSites[siteI][5],8,2),
					 					
					 "| LRT: ",			Format(fullSites[siteI][9],7,2),
					 "| p: ",			Format(fullSites[siteI][10],5,2), "\n");		

	return 0;
}


/*------------------------------------------------------------------------*/

function ReceiveJobs2 (sendOrNot, nullAlt)
{
	MPIReceive (-1, fromNode, result_String);
	
	siteIndex = MPINodeState[fromNode-1][1];
	siteNA	  = MPINodeState[fromNode-1][2];
	
	if (sendOrNot)
	{
		MPISend (fromNode,siteLikelihood);
		MPINodeState[fromNode-1][1] = siteCount;			
		MPINodeState[fromNode-1][2] = nullAlt;			
	}
	else
	{
		MPINodeState[fromNode-1][0] = 0;
		MPINodeState[fromNode-1][1] = -1;		
	}
	
	siteMap = siteIndex;
	
	/*fprintf ("../dump",CLEAR_FILE,result_String);*/
		
	ExecuteCommands (result_String);
	
	sRateV = siteLikelihood_MLE_VALUES["sRate"];
	dNdSV = siteLikelihood_MLE_VALUES["dNdS"];
	sRate_2V = siteLikelihood_MLE_VALUES["sRate_2"];
	dNdS_2V = siteLikelihood_MLE_VALUES["dNdS_2"];

	if (siteNA)
	{
		doneSites[siteMap][0] = sRateV;
		doneSites[siteMap][1] = dNdSV;
		doneSites[siteMap][2] = sRate_2V;
		doneSites[siteMap][3] = dNdS_2V;
		doneSites[siteMap][7] = doneSites[siteMap][7]+2*siteLikelihood_MLES[1][0];
	}
	else
	{
		doneSites[siteMap][4] = sRateV;
		doneSites[siteMap][5] = sRate_2V;
		doneSites[siteMap][6] = dNdSV;
		doneSites[siteMap][7] = doneSites[siteMap][7]-2*siteLikelihood_MLES[1][0];	
	}

	if (doneSites[siteMap][8] == 0)
	{
		doneSites[siteMap][8] = -1;
	}
	else
	{
		if (doneSites[siteMap][8] == (-1))
		{
			doneSites[siteMap][8] = 1-CChi2(doneSites[siteMap][7],1);						
			dummy = ReportSite2 (siteIndex, siteMap);
		}
	}
	
	return fromNode-1;
}

/*------------------------------------------------------------------------*/

#include   "TemplateModels/chooseGeneticCode.def";


global AC 	= 1;
global AT 	= 1;
global CG 	= 1;
global CT 	= 1;
global GT 	= 1;
global dNdS = 1;		

global AC_2 	= 1;
global AT_2 	= 1;
global CG_2 	= 1;
global CT_2 	= 1;
global GT_2 	= 1;
global dNdS_2   = 1;		


NucleotideMatrix	 = {{*,AC*t,t,AT*t}{AC*t,*,CG*t,CT*t}{t,CG*t,*,GT*t}{AT*t,CT*t,GT*t,*}};
NucleotideMatrix_2	 = {{*,AC_2*t,t,AT_2*t}{AC_2*t,*,CG_2*t,CT_2*t}{t,CG_2*t,*,GT_2*t}{AT_2*t,CT_2*t,GT_2*t,*}};

SetDialogPrompt ("Please specify the first codon data file:");

DataSet 	  ds 		= ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);
DataSetFilter nucData 	= CreateFilter (filteredData,1);
HarvestFrequencies (overallFrequencies,  nucData,1,1,0);
HarvestFrequencies (positionFrequencies, nucData,3,1,1);
codonFrequencies = BuildCodonFrequencies(positionFrequencies);

#include	 "queryTree.bf";
treeString_1 = ""+givenTree;

SetDialogPrompt ("Please specify the second codon data file:");

DataSet 	  ds_2		= ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter filteredData_2 = CreateFilter (ds_2,3,"","",GeneticCodeExclusions);
DataSetFilter nucData_2 	= CreateFilter (filteredData_2,1);
HarvestFrequencies (overallFrequencies_2,  nucData_2,1,1,0);
HarvestFrequencies (positionFrequencies_2, nucData_2,3,1,1);
codonFrequencies_2 = BuildCodonFrequencies(positionFrequencies_2);

if (filteredData.sites != filteredData_2.sites)
{
	fprintf (stdout, "\n\nBoth data sets must have the same number of sites.\n\n");
	return 0;
}

#include	 "queryTree.bf";
treeString_2 = ""+givenTree;

ChoiceList (modelChoice, "Model Options",1,SKIP_NONE,
			"Default","Use HKY85 and MG94xHKY85.",
			"Custom", "Use any reversible nucleotide model crossed with MG94.");
			
if (modelChoice < 0)
{
	return;
}

modelDesc = "";

if (modelChoice)
{
	done = 0;
	while (!done)
	{
		fprintf (stdout,"\nPlease enter a 6 character model designation for the first data set(e.g:010010 defines HKY85):");
		fscanf  (stdin,"String", modelDesc);
		if (Abs(modelDesc)==6)
		{	
			done = 1;
		}
	}			

	done = 0;
	while (!done)
	{
		fprintf (stdout,"\nPlease enter a 6 character model designation for the second data set(e.g:010010 defines HKY85):");
		fscanf  (stdin,"String", modelDesc2);
		if (Abs(modelDesc2)==6)
		{	
			done = 1;
		}
	}			
}
else
{
	modelDesc  = "010010";
	modelDesc2 = "010010";
}

ModelTitle = ""+modelDesc[0];
rateBiasTerms = {{"AC","1","AT","CG","CT","GT"}};
paramCount	  = 0;

dummy = MakeModelConstraints (0);

if (Abs(modelConstraintString))
{
	ExecuteCommands (modelConstraintString);
}

modelDesc = modelDesc2;
rateBiasTerms = {{"AC_2","1","AT_2","CG_2","CT_2","GT_2"}};
paramCount	  = 0;

dummy = MakeModelConstraints (0);

if (Abs(modelConstraintString))
{
	ExecuteCommands (modelConstraintString);
}

Model NucModel	    = (NucleotideMatrix, overallFrequencies, 1);
Tree  givenTree 	= treeString_1;

LikelihoodFunction 	nucLF = (nucData,givenTree);
Optimize (res, nucLF);

fprintf (stdout, "\nNucleotide fit for the first data set\n\n",nucLF,"\n\n");

Model NucModel_2 	= (NucleotideMatrix_2, overallFrequencies_2, 1);
Tree  givenTree_2 	= treeString_2;

LikelihoodFunction 	nucLF_2 = (nucData_2,givenTree_2);
Optimize (res_2, nucLF_2);

fprintf (stdout, "\nNucleotide fit for the second data set\n\n",nucLF_2,"\n\n");
				
CodonMatrix = {ModelMatrixDimension,ModelMatrixDimension};

hshift = 0;

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
						CodonMatrix[h-hshift][v-vshift] := AC__*synRate*positionFrequencies__[transition__][nucPosInCodon__];
						CodonMatrix[v-vshift][h-hshift] := AC__*synRate*positionFrequencies__[transition2__][nucPosInCodon__];
					}
					else
					{
						CodonMatrix[h-hshift][v-vshift] := AC__*nonSynRate*positionFrequencies__[transition__][nucPosInCodon__];
						CodonMatrix[v-vshift][h-hshift] := AC__*nonSynRate*positionFrequencies__[transition2__][nucPosInCodon__];
					}
				}
				else
				{
					if (trLG==2)
					{
						if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
						{
							CodonMatrix[h-hshift][v-vshift] := synRate*positionFrequencies__[transition__][nucPosInCodon__];
							CodonMatrix[v-vshift][h-hshift] := synRate*positionFrequencies__[transition2__][nucPosInCodon__];
						}
						else
						{
							CodonMatrix[h-hshift][v-vshift] := nonSynRate*positionFrequencies__[transition__][nucPosInCodon__];
							CodonMatrix[v-vshift][h-hshift] := nonSynRate*positionFrequencies__[transition2__][nucPosInCodon__];
						}							
					}
					else
					{
						if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
						{
							CodonMatrix[h-hshift][v-vshift] := AT__*synRate*positionFrequencies__[transition__][nucPosInCodon__];
							CodonMatrix[v-vshift][h-hshift] := AT__*synRate*positionFrequencies__[transition2__][nucPosInCodon__];
						}
						else
						{
							CodonMatrix[h-hshift][v-vshift] := AT__*nonSynRate*positionFrequencies__[transition__][nucPosInCodon__];
							CodonMatrix[v-vshift][h-hshift] := AT__*nonSynRate*positionFrequencies__[transition2__][nucPosInCodon__];
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
							CodonMatrix[h-hshift][v-vshift] := CG__*synRate*positionFrequencies__[transition__][nucPosInCodon__];
							CodonMatrix[v-vshift][h-hshift] := CG__*synRate*positionFrequencies__[transition2__][nucPosInCodon__];
						}
						else
						{
							CodonMatrix[h-hshift][v-vshift] := CG__*nonSynRate*positionFrequencies__[transition__][nucPosInCodon__];
							CodonMatrix[v-vshift][h-hshift] := CG__*nonSynRate*positionFrequencies__[transition2__][nucPosInCodon__];
						}
					}
					else
					{
						if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
						{
							CodonMatrix[h-hshift][v-vshift] := CT__*synRate*positionFrequencies__[transition__][nucPosInCodon__];
							CodonMatrix[v-vshift][h-hshift] := CT__*synRate*positionFrequencies__[transition2__][nucPosInCodon__];
						}
						else
						{
							CodonMatrix[h-hshift][v-vshift] := CT__*nonSynRate*positionFrequencies__[transition__][nucPosInCodon__];
							CodonMatrix[v-vshift][h-hshift] := CT__*nonSynRate*positionFrequencies__[transition2__][nucPosInCodon__];
						}							
					}
				}
				else
				{
					if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
					{
						CodonMatrix[h-hshift][v-vshift] := GT__*synRate*positionFrequencies__[transition__][nucPosInCodon__];
						CodonMatrix[v-vshift][h-hshift] := GT__*synRate*positionFrequencies__[transition2__][nucPosInCodon__];
					}
					else
					{
						CodonMatrix[h-hshift][v-vshift] := GT__*nonSynRate*positionFrequencies__[transition__][nucPosInCodon__];
						CodonMatrix[v-vshift][h-hshift] := GT__*nonSynRate*positionFrequencies__[transition2__][nucPosInCodon__];
					}							
				}
			}
		}
   }
}		

CodonMatrix_2 = {ModelMatrixDimension,ModelMatrixDimension};

hshift = 0;

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
						CodonMatrix_2[h-hshift][v-vshift] := AC_2__*synRate*positionFrequencies_2__[transition__][nucPosInCodon__];
						CodonMatrix_2[v-vshift][h-hshift] := AC_2__*synRate*positionFrequencies_2__[transition2__][nucPosInCodon__];
					}
					else
					{
						CodonMatrix_2[h-hshift][v-vshift] := AC_2__*nonSynRate*positionFrequencies_2__[transition__][nucPosInCodon__];
						CodonMatrix_2[v-vshift][h-hshift] := AC_2__*nonSynRate*positionFrequencies_2__[transition2__][nucPosInCodon__];
					}
				}
				else
				{
					if (trLG==2)
					{
						if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
						{
							CodonMatrix_2[h-hshift][v-vshift] := synRate*positionFrequencies_2__[transition__][nucPosInCodon__];
							CodonMatrix_2[v-vshift][h-hshift] := synRate*positionFrequencies_2__[transition2__][nucPosInCodon__];
						}
						else
						{
							CodonMatrix_2[h-hshift][v-vshift] := nonSynRate*positionFrequencies_2__[transition__][nucPosInCodon__];
							CodonMatrix_2[v-vshift][h-hshift] := nonSynRate*positionFrequencies_2__[transition2__][nucPosInCodon__];
						}							
					}
					else
					{
						if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
						{
							CodonMatrix_2[h-hshift][v-vshift] := AT_2__*synRate*positionFrequencies_2__[transition__][nucPosInCodon__];
							CodonMatrix_2[v-vshift][h-hshift] := AT_2__*synRate*positionFrequencies_2__[transition2__][nucPosInCodon__];
						}
						else
						{
							CodonMatrix_2[h-hshift][v-vshift] := AT_2__*nonSynRate*positionFrequencies_2__[transition__][nucPosInCodon__];
							CodonMatrix_2[v-vshift][h-hshift] := AT_2__*nonSynRate*positionFrequencies_2__[transition2__][nucPosInCodon__];
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
							CodonMatrix_2[h-hshift][v-vshift] := CG_2__*synRate*positionFrequencies_2__[transition__][nucPosInCodon__];
							CodonMatrix_2[v-vshift][h-hshift] := CG_2__*synRate*positionFrequencies_2__[transition2__][nucPosInCodon__];
						}
						else
						{
							CodonMatrix_2[h-hshift][v-vshift] := CG_2__*nonSynRate*positionFrequencies_2__[transition__][nucPosInCodon__];
							CodonMatrix_2[v-vshift][h-hshift] := CG_2__*nonSynRate*positionFrequencies_2__[transition2__][nucPosInCodon__];
						}
					}
					else
					{
						if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
						{
							CodonMatrix_2[h-hshift][v-vshift] := CT_2__*synRate*positionFrequencies_2__[transition__][nucPosInCodon__];
							CodonMatrix_2[v-vshift][h-hshift] := CT_2__*synRate*positionFrequencies_2__[transition2__][nucPosInCodon__];
						}
						else
						{
							CodonMatrix_2[h-hshift][v-vshift] := CT_2__*nonSynRate*positionFrequencies_2__[transition__][nucPosInCodon__];
							CodonMatrix_2[v-vshift][h-hshift] := CT_2__*nonSynRate*positionFrequencies_2__[transition2__][nucPosInCodon__];
						}							
					}
				}
				else
				{
					if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
					{
						CodonMatrix_2[h-hshift][v-vshift] := GT_2__*synRate*positionFrequencies_2__[transition__][nucPosInCodon__];
						CodonMatrix_2[v-vshift][h-hshift] := GT_2__*synRate*positionFrequencies_2__[transition2__][nucPosInCodon__];
					}
					else
					{
						CodonMatrix_2[h-hshift][v-vshift] := GT_2__*nonSynRate*positionFrequencies_2__[transition__][nucPosInCodon__];
						CodonMatrix_2[v-vshift][h-hshift] := GT_2__*nonSynRate*positionFrequencies_2__[transition2__][nucPosInCodon__];
					}							
				}
			}
		}
   }
}

Model MGModel   = (CodonMatrix,codonFrequencies,0);
Tree  codonTree = treeString_1;

Model MGModel_2   = (CodonMatrix_2,codonFrequencies_2,0);
Tree  codonTree_2 = treeString_2;

global	sRate   		= 1;
global  sRate_2 		= 1;

DataSet jointDS = Combine (ds, ds_2);


doneSites    = {nucData.sites/3,9};
fullSites    = {nucData.sites/3,11};

smipv = MAXIMUM_ITERATIONS_PER_VARIABLE;
MAXIMUM_ITERATIONS_PER_VARIABLE = 5000;

labels = {{"dS 1","dN 1","dS 2","dN 2", "dN1/dS1", "dN2/dS2", "Joint dS 1","Joint dS 2","Joint dN/dS", "LRT","p-value"}};

ClearConstraints (codonTree);
ReplicateConstraint("this1.?.synRate:=sRate*this2.?.t__",codonTree,givenTree);
ReplicateConstraint("this1.?.nonSynRate:=dNdS*this2.?.t__",codonTree,givenTree);
ClearConstraints (codonTree_2);
ReplicateConstraint("this1.?.synRate:=sRate_2*this2.?.t__",codonTree_2,givenTree_2);
ReplicateConstraint("this1.?.nonSynRate:=dNdS_2*this2.?.t__",codonTree_2,givenTree_2);

GetDataInfo    (dupInfo, filteredData);
GetDataInfo    (dupInfo_2, filteredData_2);

donePairs = {};

if (MPI_NODE_COUNT<=1)
{
	for (siteCount = 0; siteCount < nucData.sites/3; siteCount = siteCount+1)
	{
		siteMap   = dupInfo[siteCount];
		siteMap_2 = dupInfo_2[siteCount];
		
		mapKey = ""+siteMap+","+siteMap_2;
		
		mapIndex = donePairs[mapKey];
		
		if (mapIndex == 0)
		{
			siteMap = siteCount;
			
			filterString = "";
			filterString = filterString + (siteCount*3) + "-" + (siteCount*3+2);
			filterString = "";
			filterString = filterString + (siteCount*3) + "-" + (siteCount*3+2);
			DataSetFilter filteredData   = CreateFilter (jointDS,3,filterString,speciesIndex<nucData.species,GeneticCodeExclusions);
			DataSetFilter filteredData_2 = CreateFilter (jointDS,3,filterString,speciesIndex>=nucData.species,GeneticCodeExclusions);

			/* check to see if the site is constant in the first data set */
			
			HarvestFrequencies (f1, filteredData, 3, 3, 0);
			m1 = 0;
			for (mpiNode=0; mpiNode < 64; mpiNode=mpiNode+1)
			{
				if (f1[mpiNode]>0)
				{
					m1=m1+1;
				}
			}
			
			/* check to see if the site is constant in the first data set */
			HarvestFrequencies (f1, filteredData_2, 3, 3, 0);
			m2 = 0;
			for (mpiNode=0; mpiNode < 64; mpiNode=mpiNode+1)
			{
				if (f1[mpiNode]>0)
				{
					m2=m2+1;
				}
			}
			
			if (m1>1 || m2>1)
			/* at least one site is not constant */
			{
				sRate   = 1;
				sRate_2 = 1;
				dNdS    = 1;
				dNdS_2  = 1;
				if (m1 > 1)
				{
					if (m2 > 1)
					{
						ClearConstraints (codonTree);
						ReplicateConstraint("this1.?.synRate:=sRate*this2.?.t__",codonTree,givenTree);
						ReplicateConstraint("this1.?.nonSynRate:=dNdS*this2.?.t__",codonTree,givenTree);
						ClearConstraints (codonTree_2);
						ReplicateConstraint("this1.?.synRate:=sRate_2*this2.?.t__",codonTree_2,givenTree_2);
						ReplicateConstraint("this1.?.nonSynRate:=dNdS_2*this2.?.t__",codonTree_2,givenTree_2);
						LikelihoodFunction siteLikelihood = (filteredData, codonTree,filteredData_2,codonTree_2);				
						Optimize (site_res, siteLikelihood);
					}
					else
					{
						ClearConstraints (codonTree);
						ReplicateConstraint("this1.?.synRate:=sRate*this2.?.t__",codonTree,givenTree);
						ReplicateConstraint("this1.?.nonSynRate:=dNdS*this2.?.t__",codonTree,givenTree);
						LikelihoodFunction siteLikelihood = (filteredData, codonTree);				
						Optimize (site_res, siteLikelihood);
						sRate_2 = 0;
						dNdS_2  = 0;
					}
				}
				else
				{
					ClearConstraints (codonTree_2);
					ReplicateConstraint("this1.?.synRate:=sRate_2*this2.?.t__",codonTree_2,givenTree_2);
					ReplicateConstraint("this1.?.nonSynRate:=dNdS_2*this2.?.t__",codonTree_2,givenTree_2);
					LikelihoodFunction siteLikelihood = (filteredData_2,codonTree_2);				
					Optimize (site_res, siteLikelihood);				
					sRate   = 0;
					dNdS    = 0;
				}
				doneSites[siteMap][0] = sRate;
				doneSites[siteMap][1] = dNdS;
				doneSites[siteMap][2] = sRate_2;
				doneSites[siteMap][3] = dNdS_2;

				sRate   = 1;
				sRate_2 = 1;
				dNdS    = 1;
				ClearConstraints (codonTree);
				ReplicateConstraint("this1.?.synRate:=sRate*this2.?.t__",codonTree,givenTree);
				ReplicateConstraint("this1.?.nonSynRate:=dNdS*sRate*this2.?.t__",codonTree,givenTree);
				ClearConstraints (codonTree_2);
				ReplicateConstraint("this1.?.synRate:=sRate_2*this2.?.t__",codonTree_2,givenTree_2);
				ReplicateConstraint("this1.?.nonSynRate:=dNdS*sRate_2*this2.?.t__",codonTree_2,givenTree_2);

				Optimize (site_resN, siteLikelihood);
				doneSites[siteMap][4] = sRate;
				doneSites[siteMap][5] = sRate_2;
				doneSites[siteMap][6] = dNdS;
				doneSites[siteMap][7] = 2*(site_res[1][0]-site_resN[1][0]);
				doneSites[siteMap][8] = 1-CChi2(doneSites[siteMap][7],1);		
			}
			else
			{
				doneSites[siteMap][7] = 0;
				doneSites[siteMap][8] = 1;					
			}					
			donePairs [mapKey] = siteCount+1;
			ReportSite2 (siteCount, siteCount);	
		}
		else
		{
			dummy = ReportSite2 (siteCount, mapIndex-1);				 
		}
	}	
}
else
{
	MPINodeState = {MPI_NODE_COUNT-1,3};
	for (siteCount = 0; siteCount < nucData.sites/3; siteCount = siteCount+1)
	{
		siteMap   = dupInfo[siteCount];
		siteMap_2 = dupInfo_2[siteCount];
		mapKey = ""+siteMap+","+siteMap_2;		
		mapIndex = donePairs[mapKey];
		
		if (mapIndex == 0)
		{
			
			filterString = "";
			filterString = filterString + (siteCount*3) + "-" + (siteCount*3+2);
			DataSetFilter filteredData   = CreateFilter (jointDS,3,filterString,speciesIndex<nucData.species,GeneticCodeExclusions);
			DataSetFilter filteredData_2 = CreateFilter (jointDS,3,filterString,speciesIndex>=nucData.species,GeneticCodeExclusions);
			
			HarvestFrequencies (f1, filteredData, 3, 3, 0);
			m1 = 0;
			for (mpiNode=0; mpiNode < 64; mpiNode=mpiNode+1)
			{
				if (f1[mpiNode]>0)
				{
					m1=m1+1;
				}
			}
			
			/* check to see if the site is constant in the first data set */
			HarvestFrequencies (f1, filteredData_2, 3, 3, 0);
			m2 = 0;
			for (mpiNode=0; mpiNode < 64; mpiNode=mpiNode+1)
			{
				if (f1[mpiNode]>0)
				{
					m2=m2+1;
				}
			}
			
			if (m1>1 || m2>1)
			{
				siteMap = siteCount;
				sRate   = 1;
				sRate_2 = 1;
				dNdS    = 1;
				dNdS_2  = 1;
				
				ClearConstraints (codonTree);
				ReplicateConstraint("this1.?.synRate:=sRate*this2.?.t__",codonTree,givenTree);
				ReplicateConstraint("this1.?.nonSynRate:=dNdS*this2.?.t__",codonTree,givenTree);
				ClearConstraints (codonTree_2);
				ReplicateConstraint("this1.?.synRate:=sRate_2*this2.?.t__",codonTree_2,givenTree_2);
				ReplicateConstraint("this1.?.nonSynRate:=dNdS_2*this2.?.t__",codonTree_2,givenTree_2);
				LikelihoodFunction siteLikelihood = (filteredData, codonTree,filteredData_2,codonTree_2);				

				for (mpiNode = 0; mpiNode < MPI_NODE_COUNT-1; mpiNode = mpiNode+1)
				{
					if (MPINodeState[mpiNode][0]==0)
					{
						break;	
					}
				}
				
				if (mpiNode==MPI_NODE_COUNT-1)
				{
					mpiNode = ReceiveJobs2 (1,1);
				}
				else
				{
					MPISend (mpiNode+1,siteLikelihood);
					MPINodeState[mpiNode][0] = 1;
					MPINodeState[mpiNode][1] = siteCount;
					MPINodeState[mpiNode][2] = 1;
				}
				
				sRate   = 1;
				sRate_2 = 1;
				dNdS    = 1;
				
				ClearConstraints (codonTree);
				ReplicateConstraint("this1.?.synRate:=sRate*this2.?.t__",codonTree,givenTree);
				ReplicateConstraint("this1.?.nonSynRate:=dNdS*sRate*this2.?.t__",codonTree,givenTree);
				ClearConstraints (codonTree_2);
				ReplicateConstraint("this1.?.synRate:=sRate_2*this2.?.t__",codonTree_2,givenTree_2);
				ReplicateConstraint("this1.?.nonSynRate:=dNdS*sRate_2*this2.?.t__",codonTree_2,givenTree_2);

				DataSetFilter filteredData   = CreateFilter (jointDS,3,filterString,speciesIndex<nucData.species,GeneticCodeExclusions);
				DataSetFilter filteredData_2 = CreateFilter (jointDS,3,filterString,speciesIndex>=nucData.species,GeneticCodeExclusions);
				LikelihoodFunction siteLikelihood = (filteredData, codonTree,filteredData_2,codonTree_2);				

				for (mpiNode = 0; mpiNode < MPI_NODE_COUNT-1; mpiNode = mpiNode+1)
				{
					if (MPINodeState[mpiNode][0]==0)
					{
						break;	
					}
				}
				if (mpiNode==MPI_NODE_COUNT-1)
				{
					mpiNode = ReceiveJobs2 (1,0);
				}
				else
				{
					MPISend (mpiNode+1,siteLikelihood);
					MPINodeState[mpiNode][0] = 1;
					MPINodeState[mpiNode][1] = siteCount;
					MPINodeState[mpiNode][2] = 0;
				}
				
			}
			else
			{
				doneSites[siteCount][7] = 0;
				doneSites[siteCount][8] = 1;					
				ReportSite2 (siteCount, siteCount);	
			}					
			donePairs [mapKey] = siteCount+1;
		}
		else
		{
			doneSites[siteCount][0] = -mapIndex;
		}
	}					
	while (1)
	{
		for (nodeCounter = 0; nodeCounter < MPI_NODE_COUNT-1; nodeCounter = nodeCounter+1)
		{
			if (MPINodeState[nodeCounter][0]==1)
			{
				fromNode = ReceiveJobs2 (0,0);
				break;	
			}
		}
		if (nodeCounter == MPI_NODE_COUNT-1)
		{
			break;
		}
	}			
			

	fprintf (stdout, "\n\n\n");

	for (siteCount = 0; siteCount < nucData.sites/3; siteCount = siteCount+1)
	{
		nodeCounter = doneSites[siteCount][0];
		if (nodeCounter < (-0.5))
		{
			ReportSite2 (siteCount, -nodeCounter-1);	
		} 
		else
		{
			ReportSite2 (siteCount, siteCount);	
		}			 
	}
}

OpenWindow (CHARTWINDOW,{{"Pairwise Comparison Results"}
					   {"labels"},
					   {"fullSites"},
					   {"Bar Chart"},
					   {"Index"},
					   {labels[10]},
					   {"Site Index"},
					   {""},
					   {labels[10]},
					   {"0"}},
					   "SCREEN_WIDTH-60;SCREEN_HEIGHT-50;30;50");
					   
SetDialogPrompt ("Save site-by-site LRT results to:");

fprintf (PROMPT_FOR_FILE,CLEAR_FILE,"dS 1,dN 1,dS 2, dN 2, dN1/dS1, dN2/dS2, Joint dS1, Joint dS2, Joint dN/dS, LRT, p-value");

dummy = Columns (fullSites);

for (nodeCounter=0; nodeCounter < Rows (fullSites); nodeCounter = nodeCounter+1)
{
	fprintf (LAST_FILE_PATH,"\n",fullSites[nodeCounter][0]);
	for (siteCount = 1; siteCount < dummy; siteCount = siteCount + 1)
	{
		fprintf (LAST_FILE_PATH,",",fullSites[nodeCounter][siteCount]);
	}
}

MAXIMUM_ITERATIONS_PER_VARIABLE = smipv;

donePairs = 0;
