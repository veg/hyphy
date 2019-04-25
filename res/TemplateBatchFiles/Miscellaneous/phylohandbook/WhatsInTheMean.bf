* 

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


fprintf (stdout, "\n\nThis example analysis fits the global MG94xHKY85 model with to a codon dataset,\n",
				 "then generates a number of simulated datasets, assuming all the same\n",
				 "model parameters as fitted to the data, except the distribution of dN/dS\n",
				 "which is drawn from the different distributions all sharing the same mean\n",
				 "and finally fits the global MG94 model to each replicate\n",
				 "to demonstrate that a number of selective pressure profiles can result\n",
				 "in the same average dN/dS.\n");

ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Utility"+DIRECTORY_SEPARATOR+"DescriptiveStatistics.bf");
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

ChoiceList (freqs, "Alignment", 1, SKIP_NONE, "Flu H5N1 HA", 		"Use the example Influenza H5N1 hemagglutinin alignment (5 sequences)",
											  "Drospophila adh", 	"Use the example Drosophila ADH alignment (6 sequences).",
											  "Custom", 			"Load your own alignment and tree.");
													 
if (freqs < 0)
{
	return 0;
}
													 
if (freqs == 2)
{
	SetDialogPrompt     ("Choose a nucleotide alignment");
	DataSet ds        = ReadDataFile (PROMPT_FOR_FILE);
}
else
{
	if (freqs == 1)
	{
		GetURL			(dataString, "http://www.hyphy.org/phylohandbook/data/Drosophila_adh.nex");	
	}
	else
	{
		GetURL			(dataString, "http://www.hyphy.org/phylohandbook/data/H5N1_HA_5.nex");
	}	
	DataSet ds		 = ReadFromString (dataString);
}

HarvestFrequencies  (baseFreqs,ds,3,1,1);
DataSetFilter	  filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);
ecf		   		  = BuildCodonFrequencies (baseFreqs);
PopulateModelMatrix ("MGREVMX", baseFreqs, 0);
Model MGREV 	  = (MGREVMX, ecf, 0);

AC					= 1;
AT					:= AC;
CG					:= AC;
CT					:= 1;
GT					:= AC;

ExecuteAFile 		(HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"queryTree.bf");

LikelihoodFunction  baseLikelihoodFunction = (filteredData, givenTree);

fprintf (stdout, "Step 1. Fitting global MG94xHKY85 to observed data. \n");
Optimize (res, baseLikelihoodFunction);
fprintf	 (stdout, "Done. Estimated mean dN/dS = ", Format (R, 5,2), " and ts/tv ratio = ", 1/AC, "\n");

savedAC = AC;
savedR  = R;

ClearConstraint (c);
global 	c := 1;
categDef1 = "AC = savedAC;R  = savedR;";
simulateAndFitADistribution (2," a point mass distribution (global model) with dN/dS = " + R);

category c = (2,{{0.5}{0.5}} ,MEAN, ,{{0.5}{1.5}}, 0,2);
simulateAndFitADistribution (3," a bimodal distribution with 50% at " + R*0.5 + " and 50% at " + R*1.5);

category c = (3,{{0.6}{0.3}{0.1}} ,MEAN, ,{{0.1}{1}{6.4}}, 0,10);
simulateAndFitADistribution (4," a trimodal distribution with 60% at " + R*0.1 + ", 30% at " + R + " and 10% at " + R*6.4);

global alpha = 0.1;

category catVar =  (4,		/* number of rates */
					EQUAL,  /* probs. of rates */
					MEAN,	/* sampling method */
					GammaDist(_x_,alpha,alpha), /* density */
					CGammaDist(_x_,alpha,alpha), /*CDF*/
					0, 				   /*left bound*/
					1e25, 			   /*right bound*/
					CGammaDist(_x_,alpha+1,alpha)
						/* "antiderivative" of x f(x) */
				   );

simulateAndFitADistribution (5," a gamma distribution with variance "+ R^2*10);

/*---------------------------------------------------------------------------------------------------------------------------------------------*/
function simulateAndFitADistribution (stepID, distributionString)
{
	fprintf (stdout, "\nStep ", stepID, "a. Simulating under ", distributionString, "\n");
	PopulateModelMatrix ("MGREVMX_RV", baseFreqs, 2);
	Model MGREV_RV 	  = (MGREVMX_RV, ecf, 0);
	Tree  simT		  = treeString;
	ReplicateConstraint ("this1.?.?:=this2.?.?", simT, givenTree);
	LikelihoodFunction  simLikelihoodFunction = (filteredData, simT);
	DataSet simmedData = SimulateDataSet (simLikelihoodFunction);
	HarvestFrequencies  (baseFreqs_sim,simmedData,3,1,1);
	DataSetFilter	  sim_filter = CreateFilter (simmedData,3,"","",GeneticCodeExclusions);
	ecf_sim		   	  = BuildCodonFrequencies (baseFreqs_sim);
	PopulateModelMatrix ("MGREVMX_Sim", baseFreqs_sim, 0);
	Model MGREV_Sim 	  = (MGREVMX_Sim, ecf_sim, 0);
	Tree simT			  = treeString;
	LikelihoodFunction  simLikelihoodFunction = (sim_filter, simT);
	fprintf (stdout, "Step ", stepID, "b. Fitting a global model...\n");
	Optimize (res, simLikelihoodFunction);
	fprintf	 (stdout, "Done. Estimated mean dN/dS = ", Format (R, 5,2), " and ts/tv ratio = ", 1/AC, "\n");
	ClearConstraint (c);
	return 0;
}
