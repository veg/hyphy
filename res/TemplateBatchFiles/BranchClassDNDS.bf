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
					if (rvOptions == 2)
					{
						modelDefString*("ModelMatrixName["+hs+"]["+vs+"] := "+_nucBiasTerms[transition][transition2]+"c*nonSynRate*EFV__["+ts+"]["+ps+"];\n"+
														 "ModelMatrixName["+vs+"]["+hs+"] := "+_nucBiasTerms[transition][transition2]+"c*nonSynRate*EFV__["+ts2+"]["+ps+"];\n");	
					}
					else
					{
						modelDefString*("ModelMatrixName["+hs+"]["+vs+"] := "+_nucBiasTerms[transition][transition2]+"nonSynRate*EFV__["+ts+"]["+ps+"];\n"+
														 "ModelMatrixName["+vs+"]["+hs+"] := "+_nucBiasTerms[transition][transition2]+"nonSynRate*EFV__["+ts2+"]["+ps+"];\n");						
					}
				}
			}
	    }
    }		
	modelDefString*0;
	ExecuteCommands (modelDefString);
	return 0;
}


/*---------------------------------------------------------------------------------------------------------------------------------------------*/



ExecuteAFile  (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"chooseGeneticCode.def");

SetDialogPrompt 	  		  ("Locate a codon file:");
DataSet 		ds  		 = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter	filteredData = CreateFilter(ds,3,"","",GeneticCodeExclusions);

numberOfSamples = filteredData.sites;

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

global AC = 1;
global AT = 1;
global CG = 1;
global CT = 1;
global GT = 1;

catCounterAL = {totalBranches,1};

ChoiceList  (rvOptions,"Site-to-site rate variation",1,NO_SKIP,
					   "None","No site-to-site rate variation",
					   "Mutation Rate","Vary mutation rate site-to-site",
					   "dN/dS","Vary dN/dS ratio site-to-site",
					   "Both","Vary dS and dN across sites");
					 

MGCustomRateBiasTerms = {{"AC*","","AT*","CG*","CT*","GT*"}};	

if (rvOptions < 0)
{
	return;
}

if (rvOptions > 0)
{
	if (rvOptions < 3)
	{
		ExecuteAFile  (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"defineGamma.mdl");
	}
	else
	{
		
		ExecuteAFile  (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Utility"+DIRECTORY_SEPARATOR+"GrabBag.bf");	
		ExecuteAFile  (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TreeTools.ibf");	
		partCount = prompt_for_a_value ("How many site classes", 3, 2, 8, 1);
		
		generate_gdd_freqs (partCount, "freqs", "lfMix", "P", 0);
		
		
	}
	
	if (rvOptions == 1)
	{
		MGCustomRateBiasTerms = {{"AC*c*","c*","AT*c*","CG*c*","CT*c*","GT*c*"}};	
	}
}


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

SetDialogPrompt ("Locate tree and branch splits file:");
fscanf (PROMPT_FOR_FILE,"Raw",treeDef);
ExecuteCommands (treeDef);

bclasses = Abs(_branchClasses);

if (bclasses == 0)
{
	fprintf (stdout, LAST_FILE_PATH, " does not appear to contain a valid branch partition information\n");
	return 0;
}

classAllocations = {};
fprintf (stdout, "\nLoaded ", bclasses, " branch partitions\n");

doAllocation = 1;

if (Abs (_rateClasses) == bclasses)
{
	fprintf (stdout, "Would you like to use allocation information loaded from the file (y/n)?");
	fscanf (stdin,"String",doAllocation);
	if (doAllocation[0] == "n")
	{
		doAllocation = 0;
	}
	else
	{
		doAllocation = 1;
	}
}

if (bclasses > 1)
{
	cclasses = 0;
	while (cclasses < 1 || cclasses > bclasses)
	{
		fprintf (stdout, "How many dN/dS branch classes would you like to use [1-", bclasses,"]?");
		fscanf (stdin,"Number",cclasses);
	}
}
else
{
	cclasses = 1;
}

keys = Rows (_branchClasses);
classAllocations[0] = 1;

fprintf (stdout, "\nUsing ", cclasses, " dN/dS branch classes\n");
fprintf (stdout, "\n\nAllocated branch set ", keys[0], " to dN/dS class 1\n");

if (cclasses > 1)
{
	allocatedBC   = 1;
	allocatedRC   = 1;

	while (cclasses-allocatedRC < bclasses-allocatedBC)
	{
		allocMeTo = -1;
		maxRate = Min (allocatedRC+1,cclasses);
		while (allocMeTo < 1 || allocMeTo > maxRate)
		{
			fprintf (stdout, "Allocate branch set ", keys[allocatedBC], " to rate class [1-", maxRate,"]:");
			fscanf (stdin, "Number", allocMeTo);
		}
		classAllocations[allocatedBC] = allocMeTo;
		allocatedBC = allocatedBC + 1;
		if (allocMeTo > allocatedRC)
		{
			allocatedRC = allocatedRC + 1; 
		}
	}
	allocatedRC = allocatedRC+1;
	for (k=allocatedBC; k<bclasses; k=k+1)
	{
		fprintf (stdout, "\nAllocated branch set ", keys[k], " to rate class ",allocatedRC);
		classAllocations[k] = allocatedRC;
		allocatedRC = allocatedRC + 1;
	}
	fprintf (stdout, "\n");
}
else
{
	for (k=1; k<bclasses; k=k+1)
	{
		fprintf (stdout, "\nAllocated branch set ", keys[k], " to rate class 1");
		classAllocations [k] = 1;
	}
}

fprintf (stdout, "\n\n");

if (rvOptions == 3)
{
	treeString = Format (theTree,1,0);
	
	likeFunc = {partCount+1,1};
	
	scalerExpr = {};
	
	for (k2 = 0; k2 < bclasses; k2 += 1)
	{
		scalerExpr[k2] = {partCount,1};
	}

	for (k = 0; k < partCount; k+=1)
	{
		ExecuteCommands ("Tree T_" + k + " = " + treeString + ";\n");
		likeFunc[k] = "filteredData,T_"+k;
		for (k2 = 0; k2 < bclasses; k2 += 1)
		{
			(scalerExpr[k2])[k] = "S_S" + k + "_B" + k2;
			ExecuteCommands ("global " + (scalerExpr[k2])[k] + " = 1;");
			(scalerExpr[k2])[k] = (scalerExpr[k2])[k] + "*" + freqs[k];
		}
	}
	likeFunc[k] = "\"" + lfMix + "\"";
	
	for (k = 0; k < bclasses; k+=1)
	{
		ExecuteCommands ("global syn_scaler_" + k +" := " + Join("+",scalerExpr[k]) + "\n");
	}
		
		
	for (k=0; k<bclasses;k=k+1)
	{
		for (k3 = 0; k3 < partCount; k3+=1)
		{
			class = classAllocations [k]-1;
			rateParmID =  "omega_B"+ class + "_S" + k3;
			fprintf (stdout, "Branch set ", keys[k], " has been mapped to ",rateParmID,"\n");
				
			ExecuteCommands ("global `rateParmID`=1;");
			stringMatrix = _branchClasses[keys[k]];
			for (k2=0; k2<Columns (stringMatrix); k2=k2+1)
			{
				stringCommand = "T_"+k3+"."+stringMatrix[k2]+".nonSynRate:="+rateParmID+"*T_"+k3+"."+stringMatrix[k2]+".synRate;" + 
								"T_"+k3+"."+stringMatrix[k2]+".synRate:=S_S"+k3+"_B"+class+"/syn_scaler_"+class+"*theTree."+stringMatrix[k2]+".synRate;";

				ExecuteCommands (stringCommand);
			}
		}
	}
	
	ExecuteCommands ("LikelihoodFunction lf = (" + Join (",",likeFunc) + ");");
}
else
{
	for (k=0; k<bclasses;k=k+1)
	{
		rateParmID =  "dNdS_"+classAllocations [k];
		fprintf (stdout, "Branch set ", keys[k], " has been mapped to ",rateParmID,"\n");
			
		ExecuteCommands ("global `rateParmID`=1;");
		stringMatrix = _branchClasses[keys[k]];
		for (k2=0; k2<Columns (stringMatrix); k2=k2+1)
		{
			stringCommand = "theTree."+stringMatrix[k2]+".nonSynRate:="+rateParmID+"*theTree."+stringMatrix[k2]+".synRate;";
			/*fprintf (stdout, stringCommand, "\n");*/
			ExecuteCommands (stringCommand);
		}
	}
	LikelihoodFunction lf = (filteredData,theTree);
}

Optimize (res,lf);

fprintf (stdout, "\nFitting summary:\n\n",
				 "   Log-L = ", res[1][0], "\n",
				 "   d.f   = ", res[1][1], "\n",
				 "   sites = ", filteredData.sites, "\n",
				 "   AIC   = ", 2*(res[1][1]-res[1][0]), "\n",
				 "   c-AIC = ", 2*(res[1][1]*(filteredData.sites/(filteredData.sites-res[1][1]-1))-res[1][0]), "\n\n", lf);
				 
if (rvOptions == 3)
{
	weights = {partCount,1}["Eval(freqs[_MATRIX_ELEMENT_ROW_])"];
	fprintf (stdout, "\nSubstitution rate summary\n");

	for (k3 = 0; k3 < bclasses; k3+=1)
	{
		fprintf (stdout, "\n\tBranch class ", k3+1);
		rateMatrix = {3,partCount}["Eval(\"S_S\"+_MATRIX_ELEMENT_COLUMN_+\"_B\"+k3+\"/syn_scaler_\"+k3)"];
		rateMatrix = rateMatrix["Eval(\"omega_B\"+k3+\"_S\"+_MATRIX_ELEMENT_COLUMN_)"]["_MATRIX_ELEMENT_ROW_==1"];
		rateMatrix = rateMatrix["weights[_MATRIX_ELEMENT_COLUMN_]"]["_MATRIX_ELEMENT_ROW_==2"];
		
		fprintf (stdout, "\n\tRate class | Synonymous Rate |   Omega   |  Weight "); 
		
		for (k2 = 0; k2 < partCount; k2 += 1)
		{
			fprintf (stdout, "\n\t ", Format (k2+1,9,0), " | ", Format (rateMatrix[0][k2], 15, 3), " | ", Format (rateMatrix[1][k2], 9, 3), " | ", Format (rateMatrix[2][k2], 6, 2));
		}
		
		fprintf (stdout, "\n");
	}
	
	BL = BranchLength(T_0,-1)*weights[0];
	for (k2 = 1; k2 < partCount; k2 += 1)
	{
		BL += Eval ("BranchLength(T_"+k2+",-1)*weights[k2]");
	}
	
	TAVL = T_0^0;
	
	for (k2 = 1; k2 < Abs(TAVL); k2+=1)
	{
		(TAVL[k2])["Length"] = BL[k2-1];
	}
	fprintf (stdout, "\nScaled mixture-model tree:\n", PostOrderAVL2StringDL(TAVL,1), "\n");
	
}
				 

