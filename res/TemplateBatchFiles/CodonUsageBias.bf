#include   "TemplateModels/chooseGeneticCode.def";

/*___________________________________________________________________________________________________________*/

function PadString (padLength,padChar)
{
	for (padCounter=0;padCounter<padLength;padCounter=padCounter+1)
	{
		fprintf (stdout,padChar);
	}
	return padLength;
}

/*___________________________________________________________________________________________________________*/

function	PrintASCIITable (dataMatrix, titleMatrix)
{
	columnWidths = {1,Columns(titleMatrix)};
	for (counter1=0; counter1<Columns(titleMatrix); counter1 = counter1+1)
	{
		counter2 = Abs (titleMatrix[0][counter1])+2;
		if (counter2<12)
		{
			counter2 = 12;
		}
		columnWidths[0][counter1] = counter2;
	}
	fprintf (stdout, "\n");
	for (counter2 = 0; counter2 < Columns (titleMatrix); counter2 = counter2+1)
	{
		fprintf (stdout,"+-");
		dummy = PadString (columnWidths[0][counter2],"-");
		fprintf (stdout,"-");
	}
	fprintf (stdout,"+\n| ");
	
	for (counter1=0; counter1<Columns(titleMatrix); counter1 = counter1+1)
	{
		fprintf (stdout, titleMatrix[counter1]);
		dummy = PadString (columnWidths[0][counter1]-Abs(titleMatrix[counter1])," ");
		fprintf (stdout, " | ");
	}
	
	fprintf (stdout, "\n");
	
	for (counter1=-1; counter1<Rows(dataMatrix); counter1 = counter1 + 1)
	{
		if (counter1>=0)
		{
			fprintf (stdout,"| ");
			fprintf (stdout,Format(counter1+1,columnWidths[0][0],0));
			for (counter2 = 1; counter2 < Columns (titleMatrix); counter2 = counter2+1)
			{
				fprintf (stdout," | ");
				fprintf (stdout,Format(dataMatrix[counter1][counter2-1],columnWidths[0][counter2],-1));
			}
			fprintf (stdout," ");
			fprintf (stdout, "|\n");
		}
		for (counter2 = 0; counter2 < Columns (titleMatrix); counter2 = counter2+1)
		{
			fprintf (stdout,"+-");
			dummy = PadString (columnWidths[0][counter2],"-");
			fprintf (stdout,"-");
		}
		fprintf (stdout, "+\n");
	}
	
	return 1;
}


/*---------------------------------------------------------------------------*/

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

/*---------------------------------------------------------------------------*/

global	P1_1 = 1/4;
global  P1_2 = 1/3;
global  P1_3 = 1/2;

P1_1:<1;
P1_2:<1;
P1_3:<1;

global	F1_A := P1_1;
global	F1_C := (1-P1_1)*P1_2;
global  F1_G := (1-P1_1)*(1-P1_2)*P1_3;
global  F1_T := (1-P1_1)*(1-P1_2)*(1-P1_3);


global	P2_1 = 1/4;
global  P2_2 = 1/3;
global  P2_3 = 1/2;

P2_1:<1;
P2_2:<1;
P2_3:<1;

global	F2_A := P2_1;
global	F2_C := (1-P2_1)*P2_2;
global  F2_G := (1-P2_1)*(1-P2_2)*P2_3;
global  F2_T := (1-P2_1)*(1-P2_2)*(1-P2_3);


global	P3_1 = 1/4;
global  P3_2 = 1/3;
global  P3_3 = 1/2;

P3_1:<1;
P3_2:<1;
P3_3:<1;

global	F3_A := P3_1;
global	F3_C := (1-P3_1)*P3_2;
global  F3_G := (1-P3_1)*(1-P3_2)*P3_3;
global  F3_T := (1-P3_1)*(1-P3_2)*(1-P3_3);

nucCharacters = "ACGT";


/*---------------------------------------------------------------------------*/

function BuildCodonFrequenciesEst (freqNames)
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

	result = {ModelMatrixDimension,1};
	
	/* first pass to find the normalizing factor */
	
	constraintString = "1";
	
	
	for (h=0; h<64; h=h+1)
	{
		if (_Genetic_Code[h]==10) 
		{
			first = h$16;
			second = h%16$4;
			third = h%4;
			{
				constraintString = constraintString + "-F1_" + nucCharacters[first] + "*F2_" + nucCharacters[second] + "*F3_" + nucCharacters[third];
			}
		}
	}
	
	hshift = 0;
	
	ExecuteCommands ("global F_NORM:="+constraintString+";");
	
	dummy=constraintString*1024; 
	dummy=constraintString*(freqNames+"={"+ModelMatrixDimension+",1};");

	for (h=0; h<64; h=h+1)
	{
		if (_Genetic_Code[h]==10) 
		{
			hshift = hshift+1;
		}
		else
		{
			first = h$16;
			second = h%16$4;
			third = h%4;
			dummy = constraintString*(freqNames + "[" +(h-hshift)+ "]:=(F1_" + nucCharacters[first] + "*F2_" + nucCharacters[second] + "*F3_" + nucCharacters[third] + ")/F_NORM;\n");
		}
	}
	dummy=constraintString*0; 
	ExecuteCommands (constraintString);
	constraintString = 0;
	return 0;
}

/*---------------------------------------------------------------------------*/

global AC 	= 1;
global AT 	= 1;
global CG 	= 1;
global CT 	= 1;
global GT 	= 1;
global dNdS = 1;		

NucleotideMatrix	 = {{*,AC*t,t,AT*t}{AC*t,*,CG*t,CT*t}{t,CG*t,*,GT*t}{AT*t,CT*t,GT*t,*}};

SetDialogPrompt ("Please specify a codon data file:");

DataSet 	  ds 		= ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);


ChoiceList (modelChoice, "Model Options",1,SKIP_NONE,
			"Default","Use MG94xHKY85.",
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
		fprintf (stdout,"\nPlease enter a 6 character model designation (e.g:010010 defines HKY85):");
		fscanf  (stdin,"String", modelDesc);
		if (Abs(modelDesc)==6)
		{	
			done = 1;
		}
	}			
}
else
{
	modelDesc = "010010";
}

ModelTitle = ""+modelDesc[0];
			
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


HarvestFrequencies (positionFrequencies, ds,3,1,1);

codonFrequencies = BuildCodonFrequencies(positionFrequencies);

#include	 "queryTree.bf";

CodonMatrix    = {ModelMatrixDimension,ModelMatrixDimension};
CodonMatrixEst = {ModelMatrixDimension,ModelMatrixDimension};

constraintString="";
dummy=constraintString*65536; 

global alpha_s  = 1;     
alpha_s         :>0.01; 
alpha_s         :<100;
 
category S = (4, EQUAL, MEAN, GammaDist(_x_,alpha_s,alpha_s),  	/* density */
							  CGammaDist(_x_,alpha_s,alpha_s), 	/* CDF */
							  0,1e25, 						   	/* support */
							  CGammaDist(_x_,alpha_s+1,alpha_s) /* E[alpha_s * x] - conditional mean */
			 );


global alpha_ns  = 1;     
alpha_ns         :>0.01; 
alpha_ns         :<100;
 
category NS = (4, EQUAL, MEAN, GammaDist(_x_,alpha_ns,alpha_ns),  	/* density */
							  CGammaDist(_x_,alpha_ns,alpha_ns), 	/* CDF */
							  0,1e25, 						   	/* support */
							  CGammaDist(_x_,alpha_ns+1,alpha_ns) /* E[alpha_ns * x] - conditional mean */
			 );


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
						CodonMatrix[h-hshift][v-vshift] := AC*S*synRate*positionFrequencies__[transition__][nucPosInCodon__];
						CodonMatrix[v-vshift][h-hshift] := AC*S*synRate*positionFrequencies__[transition2__][nucPosInCodon__];
						dummy = constraintString*("CodonMatrixEst["+(h-hshift)+"]["+(v-vshift)+"]:=AC*S*synRate*F"+(nucPosInCodon+1)+"_"+nucCharacters[transition]+";");
						dummy = constraintString*("CodonMatrixEst["+(v-vshift)+"]["+(h-hshift)+"]:=AC*S*synRate*F"+(nucPosInCodon+1)+"_"+nucCharacters[transition2]+";");						
					}
					else
					{
						CodonMatrix[h-hshift][v-vshift] := AC*dNdS*NS*synRate*positionFrequencies__[transition__][nucPosInCodon__];
						CodonMatrix[v-vshift][h-hshift] := AC*dNdS*NS*synRate*positionFrequencies__[transition2__][nucPosInCodon__];
						dummy = constraintString*("CodonMatrixEst["+(h-hshift)+"]["+(v-vshift)+"]:=AC*dNdS*NS*synRate*F"+(nucPosInCodon+1)+"_"+nucCharacters[transition]+";");
						dummy = constraintString*("CodonMatrixEst["+(v-vshift)+"]["+(h-hshift)+"]:=AC*dNdS*NS*synRate*F"+(nucPosInCodon+1)+"_"+nucCharacters[transition2]+";");						
					}
				}
				else
				{
					if (trLG==2)
					{
						if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
						{
							CodonMatrix[h-hshift][v-vshift] := S*synRate*positionFrequencies__[transition__][nucPosInCodon__];
							CodonMatrix[v-vshift][h-hshift] := S*synRate*positionFrequencies__[transition2__][nucPosInCodon__];
							dummy = constraintString*("CodonMatrixEst["+(h-hshift)+"]["+(v-vshift)+"]:=S*synRate*F"+(nucPosInCodon+1)+"_"+nucCharacters[transition]+";");
							dummy = constraintString*("CodonMatrixEst["+(v-vshift)+"]["+(h-hshift)+"]:=S*synRate*F"+(nucPosInCodon+1)+"_"+nucCharacters[transition2]+";");						
						}
						else
						{
							CodonMatrix[h-hshift][v-vshift] := dNdS*NS*synRate*positionFrequencies__[transition__][nucPosInCodon__];
							CodonMatrix[v-vshift][h-hshift] := dNdS*NS*synRate*positionFrequencies__[transition2__][nucPosInCodon__];
							dummy = constraintString*("CodonMatrixEst["+(h-hshift)+"]["+(v-vshift)+"]:=dNdS*NS*synRate*F"+(nucPosInCodon+1)+"_"+nucCharacters[transition]+";");
							dummy = constraintString*("CodonMatrixEst["+(v-vshift)+"]["+(h-hshift)+"]:=dNdS*NS*synRate*F"+(nucPosInCodon+1)+"_"+nucCharacters[transition2]+";");						
						}							
					}
					else
					{
						if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
						{
							CodonMatrix[h-hshift][v-vshift] := S*AT*synRate*positionFrequencies__[transition__][nucPosInCodon__];
							CodonMatrix[v-vshift][h-hshift] := S*AT*synRate*positionFrequencies__[transition2__][nucPosInCodon__];
							dummy = constraintString*("CodonMatrixEst["+(h-hshift)+"]["+(v-vshift)+"]:=S*AT*synRate*F"+(nucPosInCodon+1)+"_"+nucCharacters[transition]+";");
							dummy = constraintString*("CodonMatrixEst["+(v-vshift)+"]["+(h-hshift)+"]:=S*AT*synRate*F"+(nucPosInCodon+1)+"_"+nucCharacters[transition2]+";");						
						}
						else
						{
							CodonMatrix[h-hshift][v-vshift] := AT*dNdS*NS*synRate*positionFrequencies__[transition__][nucPosInCodon__];
							CodonMatrix[v-vshift][h-hshift] := AT*dNdS*NS*synRate*positionFrequencies__[transition2__][nucPosInCodon__];
							dummy = constraintString*("CodonMatrixEst["+(h-hshift)+"]["+(v-vshift)+"]:=AT*dNdS*NS*synRate*F"+(nucPosInCodon+1)+"_"+nucCharacters[transition]+";");
							dummy = constraintString*("CodonMatrixEst["+(v-vshift)+"]["+(h-hshift)+"]:=AT*dNdS*NS*synRate*F"+(nucPosInCodon+1)+"_"+nucCharacters[transition2]+";");						
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
							CodonMatrix[h-hshift][v-vshift] := S*CG*synRate*positionFrequencies__[transition__][nucPosInCodon__];
							CodonMatrix[v-vshift][h-hshift] := S*CG*synRate*positionFrequencies__[transition2__][nucPosInCodon__];
							dummy = constraintString*("CodonMatrixEst["+(h-hshift)+"]["+(v-vshift)+"]:=S*CG*synRate*F"+(nucPosInCodon+1)+"_"+nucCharacters[transition]+";");
							dummy = constraintString*("CodonMatrixEst["+(v-vshift)+"]["+(h-hshift)+"]:=S*CG*synRate*F"+(nucPosInCodon+1)+"_"+nucCharacters[transition2]+";");						
						}
						else
						{
							CodonMatrix[h-hshift][v-vshift] := CG*dNdS*NS*synRate*positionFrequencies__[transition__][nucPosInCodon__];
							CodonMatrix[v-vshift][h-hshift] := CG*dNdS*NS*synRate*positionFrequencies__[transition2__][nucPosInCodon__];
							dummy = constraintString*("CodonMatrixEst["+(h-hshift)+"]["+(v-vshift)+"]:=CG*dNdS*NS*synRate*F"+(nucPosInCodon+1)+"_"+nucCharacters[transition]+";");
							dummy = constraintString*("CodonMatrixEst["+(v-vshift)+"]["+(h-hshift)+"]:=CG*dNdS*NS*synRate*F"+(nucPosInCodon+1)+"_"+nucCharacters[transition2]+";");						
						}
					}
					else
					{
						if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
						{
							CodonMatrix[h-hshift][v-vshift] := S*CT*synRate*positionFrequencies__[transition__][nucPosInCodon__];
							CodonMatrix[v-vshift][h-hshift] := S*CT*synRate*positionFrequencies__[transition2__][nucPosInCodon__];
							dummy = constraintString*("CodonMatrixEst["+(h-hshift)+"]["+(v-vshift)+"]:=S*CT*synRate*F"+(nucPosInCodon+1)+"_"+nucCharacters[transition]+";");
							dummy = constraintString*("CodonMatrixEst["+(v-vshift)+"]["+(h-hshift)+"]:=S*CT*synRate*F"+(nucPosInCodon+1)+"_"+nucCharacters[transition2]+";");						
						}
						else
						{
							CodonMatrix[h-hshift][v-vshift] := CT*dNdS*NS*synRate*positionFrequencies__[transition__][nucPosInCodon__];
							CodonMatrix[v-vshift][h-hshift] := CT*dNdS*NS*synRate*positionFrequencies__[transition2__][nucPosInCodon__];
							dummy = constraintString*("CodonMatrixEst["+(h-hshift)+"]["+(v-vshift)+"]:=CT*dNdS*NS*synRate*F"+(nucPosInCodon+1)+"_"+nucCharacters[transition]+";");
							dummy = constraintString*("CodonMatrixEst["+(v-vshift)+"]["+(h-hshift)+"]:=CT*dNdS*NS*synRate*F"+(nucPosInCodon+1)+"_"+nucCharacters[transition2]+";");						
						}							
					}
				}
				else
				{
					if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
					{
						CodonMatrix[h-hshift][v-vshift] := S*GT*synRate*positionFrequencies__[transition__][nucPosInCodon__];
						CodonMatrix[v-vshift][h-hshift] := S*GT*synRate*positionFrequencies__[transition2__][nucPosInCodon__];
						dummy = constraintString*("CodonMatrixEst["+(h-hshift)+"]["+(v-vshift)+"]:=S*GT*synRate*F"+(nucPosInCodon+1)+"_"+nucCharacters[transition]+";");
						dummy = constraintString*("CodonMatrixEst["+(v-vshift)+"]["+(h-hshift)+"]:=S*GT*synRate*F"+(nucPosInCodon+1)+"_"+nucCharacters[transition2]+";");						
					}
					else
					{
						CodonMatrix[h-hshift][v-vshift] := GT*dNdS*NS*synRate*positionFrequencies__[transition__][nucPosInCodon__];
						CodonMatrix[v-vshift][h-hshift] := GT*dNdS*NS*synRate*positionFrequencies__[transition2__][nucPosInCodon__];
						dummy = constraintString*("CodonMatrixEst["+(h-hshift)+"]["+(v-vshift)+"]:=GT*dNdS*NS*synRate*F"+(nucPosInCodon+1)+"_"+nucCharacters[transition]+";");
						dummy = constraintString*("CodonMatrixEst["+(v-vshift)+"]["+(h-hshift)+"]:=GT*dNdS*NS*synRate*F"+(nucPosInCodon+1)+"_"+nucCharacters[transition2]+";");						
					}							
				}
			}
		}
   }
}		

dummy = constraintString*0;
ExecuteCommands (constraintString);
constraintString = 0;

Model MGModelObs   = (CodonMatrix,codonFrequencies,0);
Tree	obsTree	   = treeString;

dummy 			   = BuildCodonFrequenciesEst("codonFrequenciesEst");
Model MGModelEst   = (CodonMatrixEst,codonFrequenciesEst,0);
Tree	estTree	   = treeString;

ReplicateConstraint ("this1.?.synRate:=this2.?.synRate",estTree,obsTree);

global			  MP = .5;
MP:<1;

LikelihoodFunction lf = (filteredData,obsTree,filteredData,estTree,"Log(MP*SITE_LIKELIHOOD[0]+(1-MP)*SITE_LIKELIHOOD[1])");
Optimize (res,lf);

fprintf (stdout, "\n", lf, "\n");

labelMatrix = {{"Position","A","C","G","T"}};


fprintf (stdout,"\n\nObserved nucleotide frequency estimates\n\n");
dummy = PrintASCIITable (Transpose(positionFrequencies),labelMatrix);

estFreqs = {{F1_A,F1_C,F1_G,F1_T}{F2_A,F2_C,F2_G,F2_T}{F3_A,F3_C,F3_G,F3_T}};

fprintf (stdout,"\nEstimated nucleotide frequency estimates (Prop=",1-MP,")\n\n");
dummy = PrintASCIITable (estFreqs,labelMatrix);

assignments = {filteredData.sites,2};

ConstructCategoryMatrix (marginals,lf,COMPLETE);

idx2 = filteredData.sites;

for (idx=0; idx<filteredData.sites; idx=idx+1)
{
	PrObs = MP*marginals[idx];
	PrEst = (1-MP)*marginals[idx2];
	
	assignments [idx][0] = PrObs/(PrObs+PrEst);
	assignments [idx][1] = PrEst/(PrObs+PrEst);
	
	idx2 = idx2 + 1;
}

labelMatrix = {{"Observed Frequencies","Estimated Frequencies"}};
specString  = labelMatrix[0]+";"+labelMatrix[1];

OpenWindow    (CHARTWINDOW,{{"Model Posteriors"}
						   {"labelMatrix"},
						   {"assignments"},
						   {"Stacked Bars"},
						   {"Index"},
						   {specString},
						   {"Codon"},
						   {""},
						   {"Model Posterior"},
						   {"3"}},
						   "(SCREEN_WIDTH-30);(SCREEN_HEIGHT-50);10;30");




