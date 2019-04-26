* 

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


fprintf (stdout, "\n\nThis example analysis fits a series of global MG94 models \n",
				 "with different nucleotide biases \n",
				 "and estimates the effect of these corrections on alpha and beta\n");
				 
ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Utility"+DIRECTORY_SEPARATOR+"DescriptiveStatistics.bf");
ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"chooseGeneticCode.def");


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

DataSetFilter	  	filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);
dNdSEstimates 		= {203,3};

modelVector			= {6,1};
modelCounter		= 0;
modelStrings		= {};

produceAModel 		(0,0,0);


dsc    = GatherDescriptiveStats(dNdSEstimates[-1][0]);
sorted = dNdSEstimates%1;

VERBOSITY_LEVEL = 0;
				 
bestAIC = sorted[0][1];
sorted2 = (sorted[-1][1])["_MATRIX_ELEMENT_VALUE_ - bestAIC__"];
sorted2 = sorted2["Exp(-_MATRIX_ELEMENT_VALUE_)"];
sorted2 = sorted2 * (1/(Transpose(sorted2)*sorted2["1"])[0]);
maic	= (Transpose(sorted[-1][0])*sorted2)[0];

fprintf (stdout, "\n\ndN/dS report:",
				 "\n\tMean                       = ", dsc ["Mean"], 
				 "\n\tMin                        = ", dsc ["Min"],
				 "\n\tMax                        = ", dsc ["Max"],
				 "\n\tBest model (",    
				 				 modelStrings[sorted[0][2]]  ,           
				 				      ")   	     = ",sorted[0][0],
				 "\n\tWorst model (",    
				 				 modelStrings[sorted[202][2]]           , 
				 				      ")   	     = ",sorted[202][0],
				 "\n\tREV                        = ", dNdSEstimates[202][0],
				 "\n\tModel Aver.                = ", maic, "\n");

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function produceAModel (currentPosition, currentValue, currentMax)
{
	for (currentValue = 0; currentValue <= currentMax; currentValue = currentValue + 1)
	{
		modelVector     [currentPosition] = currentValue;
		if (currentPosition < 5)
		{
			ExecuteCommands ("produceAModel (currentPosition+1,0,(currentValue==currentMax)*(currentMax+1)+(currentValue<currentMax)*currentMax);");
		}
		else
		{
			SKIP_MODEL_PARAMETER_LIST = 1;
			modelDesc 				  = ""+modelVector[0]+modelVector[1]+modelVector[2]+modelVector[3]+modelVector[4]+modelVector[5];
			modelType 				  = 1;
			ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"MG94custom.mdl");
			SKIP_MODEL_PARAMETER_LIST = 0;
			if (modelCounter == 0)
			{
				ExecuteAFile 		(HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"queryTree.bf");
				fprintf (stdout, "Model#  Model   beta/alpha    LogL    AC    AG    AT    CG    CT    GT   c-AIC");
			}
			
			LikelihoodFunction  theLnLik = (filteredData, givenTree);
			Optimize (res, theLnLik);
			
			logL = res[1][0];
			df = res[1][1]+9;
			cAIC = 2(-res[1][0] + df*(filteredData.sites/(filteredData.sites-df-1)));
			
			fprintf  (stdout, "\n", Format (modelCounter+1,6,0), 
							 "  ", modelDesc, 
							 Format (R,11,3), " ",
							 Format (logL,7,2), " ", 
							 Format (AC, 5,2), " ",
							 Format (1, 5,2), " ",
							 Format (AT, 5,2), " ",
							 Format (CG, 5,2), " ",
							 Format (CT, 5,2), " ",
							 Format (GT, 5,2), " ", Format (cAIC,7,2));
							  
			dNdSEstimates[modelCounter][0] = R;
			dNdSEstimates[modelCounter][1] = cAIC;
			dNdSEstimates[modelCounter][2] = modelCounter;
			modelStrings [modelCounter]    = modelDesc;
			modelCounter 				   = modelCounter + 1;
		}
	}
	return 0;
}