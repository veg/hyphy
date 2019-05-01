/* 

THIS HyPhy BATCH LANGUAGE FILE IS PROVIDED AS AN EXAMPLE/TOOL FOR 
THE 2nd EDITION OF 'The Phylogenetic Handbook'

Written by 
Sergei L Kosakovsky Pond (spond@ucsd.edu)
Art FY Poon				 (apoon@ucsd.edu)
Simon DW Frost			 (sdfrost@ucsd.edu)

*/


RequireVersion ("0.9920060830");
VERBOSITY_LEVEL = -1;

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

fprintf (stdout, "\n\nThis example analysis fits the Muse-Gaut 94 model with the \nuser selected nucleotide correction to a small codon dataset",
				 "\ntests the hypothesis that global omega=beta/alpha!=1 using LRT",
				 "\nand performs the parametric bootstrap to illustrate the distribution",
				 "\nof the likelihood ratio under the null model\n");

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

ChoiceList (freqs, "Alignment", 1, SKIP_NONE, "Flu H5N1 HA", 		"Use the example Influenza H5N1 heamagglutinin alignment (5 sequences)",
											  "Drospophila adh", 	"Use the example Drosophila ADH alignment (6 sequences).",
											  "Custom", 			"Load your own coding alignment.");
													 
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
		DataSet ds = ReadDataFile (PATH_TO_CURRENT_BF + '/datasets/Drosophilia_adh.nex');
	}
	else
	{
		DataSet ds = ReadDataFile (PATH_TO_CURRENT_BF + '/datasets/H5N1_HA_5.nex');
	}
}

			  			  			  			  
DataSetFilter	  	filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);

SKIP_MODEL_PARAMETER_LIST = 1;
done 					  = 0;
while (!done)
{
	fprintf (stdout,"\nPlease enter a 6 character model designation (e.g:010010 defines HKY85):");
	fscanf  (stdin,"String", modelDesc);
	if (Abs(modelDesc)==6)
	{	
		done = 1;
	}
}	
modelType 				  = 1;
ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"MG94custom.mdl");
SKIP_MODEL_PARAMETER_LIST = 0;

ExecuteAFile 		(HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"queryTree.bf");
LikelihoodFunction  theLnLik = (filteredData, givenTree);


R				:= 1;
fprintf 					   (stdout, "Fitting the null model (omega=1) to the data...\n");
Optimize 					   (res, theLnLik);		
null_L			= res[1][0];
null_DF			= res[1][1] + 9;
fprintf (stdout, "H_0 has log L of ", null_L, " and ", null_DF, " model parameters.\n");
/* save all parameter values for later simulations */

stashedValues = {};
for (k=0; k<res[1][1]; k=k+1)
{
	GetString (param, theLnLik, k);
	ExecuteCommands ("stashedValues[\""+param+"\"]="+param+";");
}

R				= 1;
fprintf 		(stdout, "Fitting the alternative model (omega is estimated) to the data...\n");
Optimize 		(res, theLnLik);		
alt_L			= res[1][0];
alt_DF			= res[1][1] + 9;
fprintf (stdout, "H_A has log L of ", alt_L, " and ", alt_DF, " model parameters. omega MLE = ", R, "\n");
LR				= 2 (alt_L-null_L);
DF				= alt_DF - null_DF;
fprintf (stdout, "The likelihood ratio test statistic (LR) = ", Format (LR,10,5), "\n");

fprintf (stdout, "Critical values to reject H_0 are");
FindRoot (root, CChi2(x,1)-0.95 ,x, 0, 1000);
fprintf (stdout,"\n\tp=0.05 :", Format (root, 10, 5));
FindRoot (root, CChi2(x,1)-0.99 ,x, 0, 1000);
fprintf (stdout,"\n\tp=0.01 :", Format (root, 10, 5));
FindRoot (root, CChi2(x,1)-0.999 ,x, 0, 1000);
fprintf (stdout,"\n\tp=0.001:", Format (root, 10, 5));
fprintf (stdout, "\nThe null hypothesis can be rejected at the alpha-level (p-value) of ", Format(1-CChi2(LR,DF),20,10), "\n");

iterates = 100;
fprintf (stdout, "\nRunning ", iterates," parametric simulations to estimate the distribution of LR ...\n");
fprintf (stdout, "Iterate | omega |   LR   | Prop (LR>",Format(LR, 5,2), ")  | Est Remaining Time, secs\n");
simulatedLR = {iterates,1};

timer = Time (1);
sim_p = 0;

paramNames = Rows (stashedValues);
maxHist	   = 8;
bins	   = 20;
step	   = maxHist/bins;
LR_counts  = {bins,2}["_MATRIX_ELEMENT_ROW_*(_MATRIX_ELEMENT_COLUMN_==0)*step"];

for (it = 0; it < iterates; it = it+1)
{
	for (k=0; k<Abs(stashedValues); k=k+1)
	{
		ExecuteCommands (paramNames[k]+"=stashedValues[\""+paramNames[k]+"\"];");
	}
	R = 1;
	DeleteObject (simLnLik, :shallow);
	DataSet 			simData 		= 	SimulateDataSet (theLnLik);
	DataSetFilter	  	simFilter 	= 	CreateFilter (simData,3,"","",GeneticCodeExclusions);
	HarvestFrequencies (simFreq,simFilter,3,1,1);
   	PopulateModelMatrix ("MG94sim", simFreq);
	simCodonF = BuildCodonFrequencies (simFreq);
	
	Tree simTree = treeString;
	
	R:=1;
	LikelihoodFunction  simLnLik = (simFilter, simTree);
	Optimize 		(simRes, simLnLik);		
	sim_null_L		= simRes[1][0];
	sim_null_DF		= simRes[1][1] + 9;
	R=1;
	Optimize 		(simRes, simLnLik);		
	sim_alt_L		= simRes[1][0];
	sim_alt_DF		= simRes[1][1] + 9;
	sim_LR			= 2 (sim_alt_L-sim_null_L);
	if (sim_LR>LR)
	{
		sim_p=sim_p+1;
	}
	fprintf (stdout, Format (it+1, 7,0), " |", Format (R,6,3)," |",Format (sim_LR,7,3)," |", Format (sim_p/(it+1),17,2), " | ",
					 Format((Time(1)-timer)*(iterates-it-1)/(it+1),8,2), " seconds\n");
	simulatedLR [it] = sim_LR;
	sim_LR = Min((sim_LR/step)$1,bins-1);
	LR_counts	[sim_LR][1] = LR_counts	[sim_LR][1] + 1/iterates;
}

VERBOSITY_LEVEL = 0;

LR_counts[0][0] = 1e-25;

dsc    = GatherDescriptiveStats(simulatedLR);

fprintf (stdout, "\n\nLR report:",
				 "\n\tMean                       = ", dsc ["Mean"], 
				 "\n\tMin                        = ", dsc ["Min"],
				 "\n\tMax                        = ", dsc ["Max"],"\n");