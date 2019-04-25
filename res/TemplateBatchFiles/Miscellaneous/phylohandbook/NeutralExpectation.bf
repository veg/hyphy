/* 

THIS HyPhy BATCH LANGUAGE FILE IS PROVIDED AS AN EXAMPLE/TOOL FOR 
THE 2nd EDITION OF 'The Phylogenetic Handbook'

Written by 
Sergei L Kosakovsky Pond (spond@ucsd.edu)
Art FY Poon				 (apoon@ucsd.edu)
Simon DW Frost			 (sdfrost@ucsd.edu)

*/

RequireVersion ("0.9920060830");

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

/*---------------------------------------------------------------------------------------------------------------------------------------------*/


ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"chooseGeneticCode.def");

ModelMatrixDimension = 64;
for (k=0; k<64; k=k+1)
{
	if (_Genetic_Code[k] == 10)
	{
		ModelMatrixDimension = ModelMatrixDimension -1;
	}
}

ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"2RatesAnalyses"+DIRECTORY_SEPARATOR+"MG94xREV.mdl");

ChoiceList (freqs, "Base composition", 1, SKIP_NONE, "Equal", 		"All nucleotides are equally likely",
													 "HIV-1 pol", 	"Base frequences from the HIV-1 polymerase gene",
													 "Custom", 		"Read in an alignment and use its base compostition");
													 
if (freqs < 0)
{
	return 0;
}
													 
if (freqs == 0)
{
	equalFreqs = {{0.25,0.25,0.25}
				  {0.25,0.25,0.25}			  
				  {0.25,0.25,0.25}			  
				  {0.25,0.25,0.25}};
}
else
{
	if (freqs == 1)
	{
		equalFreqs = {{    0.333103130755,    0.354898324512,    0.451571333823}
					  {    0.182965975093,     0.19593620008,    0.137179004237}
					  {    0.322726306842,    0.182924120079,    0.179560425762}
					  {     0.16120458731,     0.26624135533,    0.231689236178}};
	}
	else
	{
		SetDialogPrompt     ("Choose a nucleotide alignment");
		DataSet ds        = ReadDataFile (PROMPT_FOR_FILE);
		HarvestFrequencies  (equalFreqs,ds,3,1,1);
	}
}
			  			  
fprintf (stdout, "Transition/Transversion (kappa) ratio?");
fscanf 	(stdin, "Number", kappa);
if (kappa < 0)
{
	kappa = 0.001;
}
			  			  
ecf		   		  = BuildCodonFrequencies (equalFreqs);

PopulateModelMatrix ("MGREVMX", equalFreqs, 0);
Model MGREV 	  = (MGREVMX, ecf, 0);

Tree T 				= (1,2);
T.2.synRate 		= 1;
multFactor		    = BranchLength (T,0);

T.2.synRate			= 1./multFactor;
R 					= 1.0;
AC					= 1/kappa;
AT					= 1/kappa;
CG					= 1/kappa;
CT					= 1;
GT					= 1/kappa;

fprintf				(stdout, "\nTransition/Transversion ratio = ", kappa, "\n",
							 "Base frequencies (%)\n",
							 "Codon Position       A      C      G      T\n",
							 "           1st ", Format (100*equalFreqs[0][0],7,2),Format (100*equalFreqs[1][0],7,2),Format (100*equalFreqs[2][0],7,2),Format (100*equalFreqs[3][0],7,2),"\n",
							 "           2nd ", Format (100*equalFreqs[0][1],7,2),Format (100*equalFreqs[1][1],7,2),Format (100*equalFreqs[2][1],7,2),Format (100*equalFreqs[3][1],7,2),"\n",
						 	 "           3rd ", Format (100*equalFreqs[0][2],7,2),Format (100*equalFreqs[1][2],7,2),Format (100*equalFreqs[2][2],7,2),Format (100*equalFreqs[3][2],7,2),"\n");
							 

ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"dSdNTreeTools.ibf");

cl  = ReturnVectorsOfCodonLengths 	(ComputeScalingStencils		(0), "T");
fprintf (stdout, "\nUnder a Muse-Gaut 94 (9 frequency parameters) Markov model, the expected proportions of substitutions are\n", 
				 "\nSynonymous        : ", Format((cl["Syn"])[0]/(cl["Total"])[0]*100,10,2),"%",
				 "\nNon-synonymous    : ", Format((cl["NonSyn"])[0]/(cl["Total"])[0]*100,10,2), "%\n");