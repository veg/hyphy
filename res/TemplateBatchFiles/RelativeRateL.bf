#include "relrateBootstrapL.bf";

fprintf(stdout,"\n ---- RUNNING dN/dS RELATIVE RATE ANALYSIS ---- \n");

NICETY_LEVEL = 3;
SetDialogPrompt ("Please choose a codon data file:");
#include "TemplateModels/chooseGeneticCode.def";

fileCount = 0;
speciesCount = 0;

if ((OPTIMIZATION_PRECISION==0)||(OPTIMIZATION_PRECISION>0.001))
{
	OPTIMIZATION_PRECISION = 0.001;
}

while (speciesCount < 3)
{
	if (fileCount>=1)
	{
		SetDialogPrompt ("Please choose an additional data file:");
		DataSet ds2 = ReadDataFile (PROMPT_FOR_FILE);
		speciesCount = speciesCount + ds2.species;
		DataSet ds = Combine (ds,ds2);
	}
	else
	{
		DataSet ds = ReadDataFile (PROMPT_FOR_FILE);
		speciesCount = speciesCount + ds.species;
	}
	
	fileCount = fileCount+1;
}

firstSpec = 0;
secondSpec = 1;
thirdSpec = 2;
fprintf (stdout,"The following data was read:\n",ds,"\n");

if (speciesCount>3) 
{
	speciesNames = {1,3};
	ChoiceList (twoSpec, "Choose 2 taxa to share rates:",2,SKIP_NONE,ds);
	if (twoSpec[0]<0)
	{
		return ;
	}
	firstSpec = twoSpec[0];
	secondSpec = twoSpec[1];
	ChoiceList (thirdSpec, "Choose the outgroup:",1,twoSpec,ds);
	
}
else
	{		
		ChoiceList (thirdSpec,"Choose the outgroup:",1,SKIP_NONE,ds);
		
		if (thirdSpec == 1)
		{
			firstSpec = 0;
			secondSpec = 2;
		}
		else
		{
			if (thirdSpec == 2)
			{
				firstSpec = 0;
				secondSpec = 1;
			}
			else
			{
				firstSpec = 1;
				secondSpec = 2;
				thirdSpec = 0;
			}

		}
			
	}

		
DataSetFilter filteredData = CreateFilter (ds,3,"",(speciesIndex==firstSpec)||(speciesIndex==secondSpec)||(speciesIndex==thirdSpec),GeneticCodeExclusions);

ChoiceList (modelType,"Select Model Frequencies",1,SKIP_NONE,
			"1x4","MG94 with 3 frequency parameters",
			"3x4","MG94 with 9 frequency parameters");

if (modelType<0)
{
	return;
}

/* Explicit model definition */
ModelMatrixDimension = 0;

function BuildCodonFrequencies (obsF)
{
	PIStop = 1.0;
	result = {ModelMatrixDimension,1};
	hshift = 0;

	if (modelType)
	/* 9 freqs */
	{
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
	else
	/* 3 freqs */
	{

		for (h=0; h<64; h=h+1)
		{
			first = h$16;
			second = h%16$4;
			third = h%4;
			if (_Genetic_Code[h]==10) 
			{
				hshift = hshift+1;
				PIStop = PIStop-obsF[first][0]*obsF[second][0]*obsF[third][0];
				continue; 
			}
			result[h-hshift][0]=obsF[first][0]*obsF[second][0]*obsF[third][0];
		}
		return result*(1.0/PIStop);
	}
}

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
			  			}
			  			else
			  			{
			  				transition = v%16$4;
			  				transition2= h%16$4;
			  			}
			  		}
			  		if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
			  		{
			  			ModelMatrixName[h-hshift][v-vshift] := synRate*EFV__[transition__][0];
			  			ModelMatrixName[v-vshift][h-hshift] := synRate*EFV__[transition2__][0];
				  	}
			  		else
			  		{
				  		ModelMatrixName[h-hshift][v-vshift] := nonSynRate*EFV__[transition__][0];
			  			ModelMatrixName[v-vshift][h-hshift] := nonSynRate*EFV__[transition2__][0];
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
			  		if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
			  		{
			  			ModelMatrixName[h-hshift][v-vshift] := synRate*EFV__[transition__][nucPosInCodon__];
			  			ModelMatrixName[v-vshift][h-hshift] := synRate*EFV__[transition2__][nucPosInCodon__];
				  	}
			  		else
			  		{
				  		ModelMatrixName[h-hshift][v-vshift] := nonSynRate*EFV__[transition__][nucPosInCodon__];
			  			ModelMatrixName[v-vshift][h-hshift] := nonSynRate*EFV__[transition2__][nucPosInCodon__];
		  			}
			  	}
			  }
		}
	}
	return 0;
}

if (modelType)
{
	HarvestFrequencies (observedFreq,filteredData,3,1,1);
}
else
{
	HarvestFrequencies (observedFreq,filteredData,1,1,0);
}

MG94 = 0;

MULTIPLY_BY_FREQS = PopulateModelMatrix ("MG94", observedFreq);

vectorOfFrequencies = BuildCodonFrequencies (observedFreq);

Model MG94model = (MG94,vectorOfFrequencies,0);

/* define the tree with fixed species names */

if ((firstSpec<thirdSpec)&&(secondSpec<thirdSpec))
{
	treeString = "(Ingroup1,Ingroup2,OutGroup)";
}
else
{
	if ((firstSpec>thirdSpec)&&(secondSpec>thirdSpec))
	{
		treeString = "(OutGroup,Ingroup1,Ingroup2)";
	}
	else
	{
		treeString = "(Ingroup1,OutGroup,Ingroup2)";
	}
}

Tree threeTaxaTree = treeString;

LikelihoodFunction lf = (filteredData,threeTaxaTree);

GetString (specName,ds,thirdSpec);
fprintf (stdout,"\nOutgroup = ",specName);
GetString (specName,ds,firstSpec);
fprintf (stdout,"\nIngroup1 = ",specName);
GetString (specName,ds,secondSpec);
fprintf (stdout,"\nIngroup2 = ",specName);


Optimize (res,lf);

fprintf (stdout, "\n_________________RESULTS_________________\n\nFULL MODEL RESULTS\n******************\n\n",lf);

Tree constrained3TaxaTree = treeString;

/* now specify the constraint */

LikelihoodFunction lfConstrained = (filteredData,constrained3TaxaTree);

global		R;

constrained3TaxaTree.Ingroup1.nonSynRate := R*constrained3TaxaTree.Ingroup1.synRate;
constrained3TaxaTree.Ingroup2.nonSynRate := R*constrained3TaxaTree.Ingroup2.synRate;

Optimize (res1,lfConstrained);

fprintf (stdout, "\n\nCONSTRAINED MODEL RESULTS\n*************************\n\n",lfConstrained);

lnLikDiff = -2(res1[1][0]-res[1][0]);

degFDiff = res[1][1]-res1[1][1];

fprintf (stdout, "\n\nLik. ratio statistic = ",lnLikDiff,"\n",
					 "Degrees of freedom   = ",Format(degFDiff,0,0));

fprintf (stdout, "\nChi2 Derived P-Value = ",1-CChi2(lnLikDiff,degFDiff));
