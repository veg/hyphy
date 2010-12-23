RequireVersion ("2.0020101201");

/*----------Branch-Site Directional Selection Analysis---------*/
/*Usage: Takes an alignment and a tree (with foreground branches tagged with "{FG}"), a nucleotide model
string ("010010" for HKY85, "012345" for REV), a range of sites to test, and a site reference shift,
which allows for alignments that don't start where they should*/

/*--Analysis Setup--*/


fprintf (stdout, "\n[RUNNING MEDS (Models for Episodic Directional Selection). For help please refer to http://www.hyphy.org/wiki/MEDS ]\n");
nucModelString 	= "012345";

SetDialogPrompt ("Load a coding alignment");
/*---------Loading alignment and tree files-------------------------------------*/

DataSet 		myData = ReadDataFile (PROMPT_FOR_FILE);
fprintf (stdout, "Loaded ", myData.species, " sequences with ", myData.sites, " sites from ",LAST_FILE_PATH,"\n");

SetDialogPrompt ("Load an annotated tree file"); 
fscanf 			(PROMPT_FOR_FILE, "Raw", treeString);

SetDialogPrompt ("Specify the output (.csv) file"); 
fprintf 		(PROMPT_FOR_FILE, CLEAR_FILE);

outputFile 		= LAST_FILE_PATH;

siteShift 		= -1;  /*Used to standardize codon positions. Remember HyPhy indexes from 0, so siteShift = -1 will report first codon as 1*/


/*--Code Overview:--*/
/*
1) Estimate branch lengths by fitting a custom nuc model. These are used throughout.
2) Fit foreground and background codon models, with all parameters except the nonsyn rate tied
	between fg and bg. This is to get estimates of the nucleotide transition rates from a codon model.
3) Site by site FEL directional selection analysis: Nuc transition rates are fixed from 2).
	For each site, fit a positive selection null model that allows freely variable syn rate,
	and nonsynBG rates but with nonsynFG constrained to equal syn. Fit a second positive selection
	model, with nonsynFG unconstrained. An LRT is used to determine whether the unconstrained
	positive selection model significantly outperformed the null. So far this is only positive
	selection. The unconstrained positive selection model will be used as the null for tests of
	directional selection. For each of 20 AAs, fit a seperate model that allows a rate multiplier
	favouring substitutions towards that AA. This produces 20 nested models, so a set Bonferroni
	corrected LRTs can be used to identify evidence of direction selection. This batch file outputs
	p-values BEFORE Bonferroni correction.
*/

/*--Code--*/
ACCEPT_ROOTED_TREES=1;
_DO_TREE_REBALANCE_ = 0;
COUNT_GAPS_IN_FREQUENCIES = 0;

LoadFunctionLibrary ("chooseGeneticCode");
countSenseCodons = + ({64,1}["IsStop(_MATRIX_ELEMENT_ROW_,_Genetic_Code)==0"]);


/*---------Estimate Branch Lengths Using Nucleotide Model-----------------------*/
DataSetFilter myFilter = CreateFilter (myData,1);
HarvestFrequencies (obsNucFreqs, myFilter, 1, 1, 1);

/*---------Begin Setting up custom nuc model---------------*/
/*The parameters that the nucModelString indexes into*/
nucBiasMult = {{"AC*","","AT*","CG*","CT*","GT*"}};

/*Initialises and populates a matrix of string values nuc multipliers*/
customRateString = {{"*","","",""}
					{"","*","",""}
					{"","","*",""}
					{"","","","*"}
				 };
				 
for (i=0; i<3; i+=1)
{
	shift = -(i == 0);
	for(j=i+1;j<4;j += 1)
	{
		customRateString[i][j] = Eval ("nucBiasMult["+nucModelString[j+i+shift]+"]");
		customRateString[j][i] = customRateString[i][j];
	}
}

global AC = 1; global AT = 1; global CG = 1; global CT = 1; global GT = 1; 

/*To set up a nucleotide model, populate a string version of the rate matrix*/
modelDefString = "";
modelDefString * 16384;

modelDefString* "{" ;
for (i=0; i<4; i += 1)
{
	modelDefString*("{");
	for(j=0;j<4;j=j+1)
	{
		if(j>0)
		{
			modelDefString*(",");
		}
		if(i==j)
		{
			modelDefString*("*");
		}
		else
		{
			modelDefString*(customRateString[i][j]+"t");
		}
	}
	modelDefString*"}";
}
modelDefString*"}";
modelDefString*0;
ExecuteCommands("nucModel = " + modelDefString);
Model RT = (nucModel, obsNucFreqs);

/*This is in case the treestring came with model assignments!*/
Model FG = (nucModel, obsNucFreqs);

/*-----------End of defining nuc model------------*/


/*FIDDLING WITH THE ROOTING. This (commented out) code reroots on an internal branch and is only certain to work on trees with external branches marked FG. This is removed and now only rooted trees are accepted.*/
/*
Tree tempTree = treeString;
treeString = RerootTree(treeString,BranchName(tempTree,0));
*/
Tree givenTree = treeString;

fprintf (stdout, "\n\n[PHASE 1. Estimating Branch Lengths using a Nucleotide Model]\n");
LikelihoodFunction theLikFun = (myFilter, givenTree, obsFreqs);
Optimize (paramValues, theLikFun);
fprintf (stdout, theLikFun);

/*---------end nuc model fit----------*/

/*---------Labels all branches containing branchID with the FG model. This code segment does nothing when the branches are explicitly tagged, and branchMatch doesn't match anything-------------*/
/*--First populate a list of branches containing branchID. branchID could be regExp text--*/
branchMatch = "DONTNAMEANYTAXATHISUNLESSYOUREALLYMEANIT"; /*This is a different (old) way to tag foreground branches - All terminal branches with any substring = branchMatch will get tagged. This is untested and explicit tree tagging is preffered*/
branchID = branchMatch;
BOInames = {};
tips = TipCount(givenTree);
loc = 0;
for (i=0; i<tips; i=i+1)
{
	tempTipName = TipName(givenTree,i);
	tempResult = tempTipName$branchID;
	if (tempResult[1][1]>=0)
		{
			BOInames[loc] = TipName(givenTree,i);
			loc=loc+1;
		}
}
/*--use regular expressions to append the FG model text--*/
numBOI = loc;
newTreeString = treeString;

fprintf (stdout, "The following branches were labeled as foreground: \n");

for (i=0; i< numBOI; i=i+1)
{
	newTreeString = newTreeString^{{BOInames[i]}{BOInames[i]+"{FG}"}};
	/*Just to check what I label*/
	//fprintf (stdout, BOInames[i], "\n");
}
treeString = newTreeString;
/*-----End branch labeling code-----*/

/*---------------------Harvesting Frequencies-----------------------*/
/*---A function that converts a nucFreqMatrix to a vector of freqs--*/
function BuildCodonFrequencies (nucFreqMatrix)
{
	
	PIStop = 1.0; 		/* denominator */
	result = {countSenseCodons,1};    /* resulting codon frequencies */
	hshift = 0;         /* how many stop codons have been counted so far */

	for (h=0; h<64; h=h+1) /* loop over all possible codons */
	{
		first  = h$16;    /* Decompose a codon into 3 nucleotides. 
							 The index of the first nucleotide (A=0,C=1,G=2,T=3) is found here,
							 by doing integer division by 16  */
		second = h%16$4;  /* The index of the second nucleotide. 
							 First take the remainder of division by 16, i.e. positions 2 and 3
							 and then extract position 2 by integer division by 4*/
		third  = h%4;     /* The index of the third nucleotide.
							 Remainder of integer division by 4*/
							 
						  /* in the end: h = 16*first + 4*second + third */
							 
		if (_Genetic_Code[h]==10) /* stop codon */ 
		{
			hshift = hshift+1; 
			PIStop = PIStop-nucFreqMatrix[first][0]*nucFreqMatrix[second][1]*nucFreqMatrix[third][2]; /* adjust the denominator */
		}
		else
		{
			result[h-hshift] = nucFreqMatrix[first][0]*nucFreqMatrix[second][1]*nucFreqMatrix[third][2]; 
											/* store the frequency for codon h. Notice the substraction of hshift to compensate
											  for the absense of stop codons. The first codon affected by it is 
											  TAC (h=49), which gets stored in result[48], because TAA (a stop codon) was skipped. */
		}
	}
	return result*(1.0/PIStop);
}
DataSetFilter 		codonFilter  = CreateFilter	(myData,3,"","","TAA,TAG,TGA"); /* define the codon filter, excluding stop codons */
HarvestFrequencies  (nuc3by4,myData,3,1,1); 			     /* collect position specific nucleotide frequencies */
estimatedCodonFreqs = BuildCodonFrequencies(nuc3by4);
/*----------------------Done with freqs-----------------------------------------*/


/*----------------------Defines a function for creating a custom codon model--------------*/
/*-------Usage: "PopulateModelMatrix ("ModelVarName",Freq3x4,targetAA,customRateString,"nonSynRateTag");"---*/
/*specify targetAA = 21 if you don't want directional selection. NonSynRateTag allows one to control the
name of the "nonsyn" rate variable by appending nonSynRateTag to the end. This allows having seperate foreground
and background nonsyn rates. "customRateString" is NOT a PAUP specifier, but rather a matrix of string valued
multiplers, derived above in the setup of the nuc model*/
/*Part copypasta from MG94customCF3x4.mdl*/

function PopulateModelMatrix (ModelMatrixName&, EFV, targetAA,customRateString,nonSynRateTag)
{
	ModelMatrixDimension = countSenseCodons;
	_localNucBiasMult = customRateString;
	ModelMatrixName = {ModelMatrixDimension,ModelMatrixDimension}; 
	modelDefString = "";
	modelDefString*16384;
	hshift = 0;
	for (h=0; h<64; h=h+1)
	{
		if (_Genetic_Code[h]==10) 
		{
			hshift = hshift+1;
			continue; 
		}
		vshift = 0;
		for (v = 0; v<64; v=v+1)
		{
			if(h==v)
			{
				continue;
			}
			diff = v-h;
			if (_Genetic_Code[v]==10) 
			{
				vshift = vshift+1;
				continue; 
			}
			nucPosInCodon = 2;
			if ((h$4==v$4)||((diff%4==0)&&(h$16==v$16))||(diff%16==0)) /* differ by one subsitution only */
			{
				if (h$4==v$4) /* third position */
				{
					transition = v%4;
					transition2= h%4;
				}
				else
				{
					if(diff%16==0) /* first position */
					{
						transition = v$16;
						transition2= h$16;
						nucPosInCodon = 0;
					}
					else   /* second position */
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
				aa1 = _Genetic_Code[h];
				aa2 = _Genetic_Code[v];
				
				synOrNon = "nonsyn"+ nonSynRateTag +"*";
				if (aa1==aa2) {synOrNon = "syn*";}
				targetAAmult = "";
				if (aa2==targetAA) {targetAAmult = "((1/(1-omegaT))-1)*";} /*This is a different way of parameterizing the omegaT multiplier that improves optimization speed. Needs to be transformed back at the end*/
				modelDefString*("ModelMatrixName["+hs+"]["+vs+"] := " + targetAAmult + synOrNon + "t*"+_localNucBiasMult[transition][transition2]+"EFV__["+ts+"]["+ps+"];\n");   /*EFV__["+ts+"]["+ps+"] multiplies the 3x4 equilibrium frequency of the column codon*/
			}
	    }
    }		
	modelDefString*0;
	ExecuteCommands (modelDefString);
	return 0;
}
/*--------------------End of custom codon model setup function-------------*/

/* --------------------Set up the background model--------------------*/
global syn = 1;
global nonsyn=1;
/*--Populate the transition matrix for the background codon model--*/
PopulateModelMatrix ("MG94xCustomRateMatrix",nuc3by4,21,customRateString,"");   
Model MG94xCustom = (MG94xCustomRateMatrix,estimatedCodonFreqs,0);

/* -------------------Set up the foreground model-------------------------*/
global nonsynFG=1;
/*Populate the transition matrix for the foreground codon model*/
PopulateModelMatrix ("MG94xCustomRateMatrixFG",nuc3by4,21,customRateString,"FG");  
Model FG = (MG94xCustomRateMatrixFG,estimatedCodonFreqs,0);

/*------------------Assign models to tree------------------------*/
UseModel(MG94xCustom); /*This sets assigns background model to all unlabeled branches*/

/*Just to check my {FG} tags!*/
fprintf (stdout, treeString);

Tree myTreeFG = treeString; /*treeString has "{FG}" appended to all foreground branches*/

/*Forces all branch lengths to be those estimated by the nuc model*/
ReplicateConstraint ("this1.?.t:=this2.?.t__",myTreeFG,givenTree);

/*----------------Optimise likelihood function--------------------*/

fprintf (stdout, "\n\n[PHASE 2. Estimating Substitution Parameters using the Global Codon Model]\n");

LikelihoodFunction lf = (codonFilter, myTreeFG);
Optimize (myRes, lf);
fprintf (stdout, "\n", lf, "\n");

/*--------Forever constrain nuc rates-------*/
AC := AC__; AT := AT__; CG := CG__; CT := CT__; GT := GT__; 

/*--------------Setting up the output file--------------------------*/
fprintf (outputFile,CLEAR_FILE,"Site,NullLL,NullBgNonSyn,NullSyn,DivLL,Div_p,DivFgNonSyn,DivBgNonSyn,DivSyn,0LL,0p,0w,0FgNonSyn,0BgNonSyn,0Syn,1LL,1p,1w,1FgNonSyn,1BgNonSyn,1Syn,2LL,2p,2w,2FgNonSyn,2BgNonSyn,2Syn,3LL,3p,3w,3FgNonSyn,3BgNonSyn,3Syn,4LL,4p,4w,4FgNonSyn,4BgNonSyn,4Syn,5LL,5p,5w,5FgNonSyn,5BgNonSyn,5Syn,6LL,6p,6w,6FgNonSyn,6BgNonSyn,6Syn,7LL,7p,7w,7FgNonSyn,7BgNonSyn,7Syn,8LL,8p,8w,8FgNonSyn,8BgNonSyn,8Syn,9LL,9p,9w,9FgNonSyn,9BgNonSyn,9Syn,10LL,10p,10w,10FgNonSyn,10BgNonSyn,10Syn,11LL,11p,11w,11FgNonSyn,11BgNonSyn,11Syn,12LL,12p,12w,12FgNonSyn,12BgNonSyn,12Syn,13LL,13p,13w,13FgNonSyn,13BgNonSyn,13Syn,14LL,14p,14w,14FgNonSyn,14BgNonSyn,14Syn,15LL,15p,15w,15FgNonSyn,15BgNonSyn,15Syn,16LL,16p,16w,16FgNonSyn,16BgNonSyn,16Syn,17LL,17p,17w,17FgNonSyn,17BgNonSyn,17Syn,18LL,18p,18w,18FgNonSyn,18BgNonSyn,18Syn,19LL,19p,19w,19FgNonSyn,19BgNonSyn,19Syn,20LL,20p,20w,20FgNonSyn,20BgNonSyn,20Syn");

/*----------------------------For loop over sites-----------------------------------*/

fprintf (stdout, "\n\n[PHASE 3. Testing for Directional Selection on Foreground Branches]\n");

for(siteIn=1;siteIn<=codonFilter.sites;siteIn += 1)
{
	fprintf (stdout, "Working on site ", siteIn, "\n");
	
	/*-------------Allow user to select sites and options-------------*/
	site = siteIn +siteShift;
	siteString = "" + (site*3) + "-" + (site*3+2);
	DataSetFilter siteFilter = CreateFilter (myData,3,siteString,"","TAA,TAG,TGA");
	
	/*-----Count Site Specific Codon and AA Freqs. This is just to exclude certain sites------*/
	HarvestFrequencies (siteCodFreqs, siteFilter, 3, 3, 0);
	AAfreqs = {21,1};
	for (h=0; h<64; h=h+1) /* loop over all possible codons */
	{
		AAfreqs[ _Genetic_Code[h]] += siteCodFreqs[h];
	}
	
	/*Count how many AA have frequencies greater than 0*/
	numGrtZero = 0;
	for (h=0; h<21; h=h+1) /* loop over all possible codons */
	{
		if(AAfreqs[h]>0)
		{
			numGrtZero = numGrtZero+1;
		}
	}
	/*-------------Only test a site if there is more than one observed AA-------------*/
	if(numGrtZero>1) 
	{
		AAlower = 0; AAupper = 20;
	
		AAlikes = {};
		AAlikesOmegaTs = {};
		AAlikesSyn = {};
		AAlikesBgNonSyn = {};
		AAlikesFgNonSyn = {};
		
		/*-------------------------loop over amino acid targets--------------------------------*/
		for (AAcount=AAlower; AAcount<AAupper+1; AAcount=AAcount+1)
		{
			/*------Only test for directional selection towards observed AAs-------*/
			if(AAfreqs[AAcount]>0) /*Stop codons are implicitly excluded*/
			{
				fprintf (stdout, "\tTesting target residue ", _hyphyAAOrdering[AAcount], "\n");
				targetAA = AAcount;
			
				/* --------------------Construct Directional Model---------------------------------*/
				global nonsynDIR=1;
				global omegaT = 0.5; /*reparameterized. 0.5 = 1*/
				PopulateModelMatrix ("MG94xCustomRateMatrixDIR",nuc3by4,targetAA,customRateString,"DIR");
				Model FG = (MG94xCustomRateMatrixDIR,estimatedCodonFreqs,0);

				UseModel(MG94xCustom); /*This assigns background model to all unlabeled branches*/
				Tree myTreeFG = treeString;
				
				/*Forces all branch lengths to be those estimated by the nuc model*/
				ReplicateConstraint ("this1.?.t:=this2.?.t__",myTreeFG,givenTree);
				omegaT :<1; /*For the reparameterization*/
				syn = 1; /*this is shared between foreground and background*/
				nonsyn=1; /*background nonsyn rate*/
				
				LikelihoodFunction lf = (siteFilter, myTreeFG);
				Optimize (mySiteRes, lf); 
				/*fprintf (stdout, "\n", lf, "\n");*/
				
				unConstrainedRatio = mySiteRes[1][0];
				AAlikes[AAcount] = unConstrainedRatio;
				AAlikesOmegaTs[AAcount] = omegaT;
				AAlikesSyn[AAcount] = syn;
				AAlikesBgNonSyn[AAcount] = nonsyn;
				AAlikesFgNonSyn[AAcount] = nonsynDIR;
			}
			else
			{
				fprintf (stdout, "\tSkipping unobserved target residue ", _hyphyAAOrdering[AAcount], "\n");
				AAlikes[AAcount] = -99999999;
				AAlikesOmegaTs[AAcount] = -99;
				AAlikesSyn[AAcount] = -99;
				AAlikesBgNonSyn[AAcount] = -99;
				AAlikesFgNonSyn[AAcount] = -99;
			}
		}
		
		/*---------------------Set up model allowing non-neutral selection on FG---------------*/
		syn = 1;
		nonsyn=1;
		nonsynDIR=1; /*Allow positive selection on FG*/
		omegaT :=0.5; /*Constrain omegaT to 1 for null model - reparameterized*/
		LikelihoodFunction lf = (siteFilter, myTreeFG);
		Optimize (mySiteRes, lf);
		constrainedRatioPos = mySiteRes[1][0];
		DivFgNonSyn = nonsynDIR;
		DivBgNonSyn = nonsyn;
		DivSyn = syn;
		/*--End non-neutral model--*/
		
		/*---------------------Set up null model forcing neutral selection on FG----------------*/
		syn = 1;
		nonsyn=1;
		nonsynDIR:=syn; /*Force neutral selection on FG*/
		omegaT :=0.5; /*Constrain omegaT to 1 for null model - reparameterized*/
		LikelihoodFunction lf = (siteFilter, myTreeFG);
		Optimize (mySiteRes, lf); 
		/*fprintf (stdout, "\n", lf, "\n");*/
		constrainedRatioNoPos = mySiteRes[1][0];
		/*--End null model--*/
		
		/*-----------------------Display and Write Results------------------------*/
		/*---"Site,NullLL,NullBgNonSyn,NullSyn,DivLL,Div_p,DivFgNonSyn,DivBgNonSyn,DivSyn,0LL,0p,0w,0FgNonSyn,0BgNonSyn,0Syn"---*/
		outputString = ""+siteIn;
		outputString = outputString + "," + constrainedRatioNoPos + "," + nonsyn + "," + syn;
		fprintf (stdout, "\nTests for site ", siteIn,"\n");
		
		/*--First a test for general positive selection at the site--*/
		lrtPosVsNoPos = 2*(constrainedRatioPos-constrainedRatioNoPos);
		pValPosVsNoPos = 1-CChi2 (lrtPosVsNoPos, 1);
		fprintf (stdout, "\tLikelihood Ratio Test for diversifying selection");
		fprintf (stdout, ": ", lrtPosVsNoPos);
		fprintf (stdout, "  p-value: ", pValPosVsNoPos, "\n");
		
		outputString = outputString + "," + constrainedRatioNoPos + "," + pValPosVsNoPos + "," + DivFgNonSyn + "," + DivBgNonSyn + "," + DivSyn;
		
		/*--Now test for directional selection vs positive selection--*/
		fprintf (stdout, "\tTesting for directional selection againts null allowing positive selection in foreground", "\n");
		for (AAcount=AAlower; AAcount<AAupper+1; AAcount=AAcount+1)
		{
			lrtScore = 2*(AAlikes[AAcount]-constrainedRatioPos);
			pValues = 1-CChi2 (2*(AAlikes[AAcount]-constrainedRatioPos), 1);
			if (AAfreqs[AAcount] > 0)
			{
				fprintf (stdout, "\tLikelihood Ratio Test for: ", _hyphyAAOrdering[AAcount]);
				fprintf (stdout, "  : ", lrtScore);
				fprintf (stdout, "  p-value: ", pValues, "\n");
			}
			outputString = outputString + "," + AAlikes[AAcount] + "," + pValues + "," + AAlikesOmegaTs[AAcount] + "," + AAlikesFgNonSyn[AAcount] + "," + AAlikesBgNonSyn[AAcount] + "," + AAlikesSyn[AAcount];
		}
		/*Testing for directional vs neutral selection can be done in post-processesing*/
		fprintf (outputFile,"\n",outputString);
		/*Clear constraints to test another site*/
		ClearConstraints(omegaT)
		ClearConstraints(nonsynDIR);
	}
	else
	{
		fprintf (stdout,"Skipped site ",siteIn," because it is invariable\n");
	}
	
} /*End main loop*/
