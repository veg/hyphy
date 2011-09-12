incFileName = HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"chooseGeneticCode.def";
ExecuteCommands  ("#include \""+incFileName+"\";");

incFileName = HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"2RatesAnalyses"+DIRECTORY_SEPARATOR+"MG94xREVxBivariate.mdl";
ExecuteCommands  ("#include \""+incFileName+"\";");

s_rateMin    = 0.01;
s_rateMax    = 1.00;

ns_rateMin    = 0.01;
ns_rateMax    = 5.00;

SetDialogPrompt ("Choose a model fit:");
fscanf (PROMPT_FOR_FILE, "Raw", inFile);
outFile = LAST_FILE_PATH;

ExecuteCommands (inFile);
GetInformation(vars,"^P_[0-9]+");
rateCount = Columns (vars)+2;

fprintf (stdout, "\nLoaded a fit with ", rateCount-1, " rates\n");

ChoiceList (searchKind, "Search directions:", 1, SKIP_NONE,
						"All", "Search all directions",
						"Negative", "Only dS > dN directions",
						"Neutral", "Only dS = dN directions",
						"Positive", "Only dS < dN directions");
						
if (searchKind < 0)
{
	return 0;
}

timerCap = 0;
while (timerCap < 1)
{
	fprintf (stdout, "\nSeconds (x3) to spend looking for directional improvement in fit (>=1, or -1 to search a fixed grid)?");
	fscanf	(stdin, "Number",timerCap);
	if (timerCap == (-1))
	{
		break;
	}
}

if (timerCap < 0)
{
	fprintf (stdout, "Sampling 355 directions...\n");
}
else
{
	fprintf (stdout, "Allocated ", timerCap*3, " seconds for approximate search...\n");
}

LFCompute (lf,LF_START_COMPUTE);
LFCompute (lf,res);
LFCompute (lf,LF_DONE_COMPUTE);

fprintf (stdout, "\nBaseline likelihood:", res);

HarvestFrequencies (observedFreq,	   filteredData,3,1,1);
ExecuteCommands ("global S_"+(rateCount-1)+"=0.5; global NS_"+(rateCount-1)+"=0.5;global P_"+(rateCount-1)+"=0.5;P_"+(rateCount-1)+":<1;");
ExecuteCommands("PopulateModelMatrix(\"rate_matrix_"+(rateCount-1)+"\",observedFreq,\"S_\"+(rateCount-1)+\"/c_scale\",\"NS_\"+(rateCount-1)+\"/c_scale\");");
ExecuteCommands("Model MG94model_"+(rateCount-1)+"= (rate_matrix_"+(rateCount-1)+",vectorOfFrequencies,0);");
ExecuteCommands ("Tree tree_"+(rateCount-1)+"="+tree_0+";");
ExecuteCommands ("ReplicateConstraint (\"this1.?.synRate:=this2.?.synRate\",tree_"+(rateCount-1)+",tree_0);global P_"+(rateCount-1)+"=1/"+filteredData.sites+";");

lfDef1 	  = "";
lfDef1 	  * 1024;
lfDef1 	  * "\"Log(";

ps 			 = {rateCount,1};

for (mi = 0; mi < rateCount -2; mi=mi+1)
{
	ExecuteCommands ("ps["+mi+"]="+"P_"+(mi+1)+";");
}

scaler		 = 1-1/filteredData.sites;
freqValueMatrix    = {rateCount,1};
freqStrMatrix      = {rateCount,1};
freqValueMatrix    = freqValueMatrix["scaler"];

for (mi=0; mi<rateCount-2; mi=mi+1)
{
	freqStrMatrix[mi]   = "";
	for (mi2 = 0; mi2 < mi; mi2=mi2+1)
	{
		freqValueMatrix[mi] = freqValueMatrix[mi] * (1-ps[mi2]);
		freqStrMatrix[mi]   = freqStrMatrix[mi]+"(1-P_"+(mi2+1)+")*";
	}
	freqValueMatrix[mi] = freqValueMatrix[mi] * ps[mi];
	freqStrMatrix[mi]   = freqStrMatrix[mi]+"P_"+(mi+1)+"";
}

freqStrMatrix[rateCount-2]   = "";
freqStrMatrix[rateCount-1]   = "";


for (mi2 = 0; mi2 < mi; mi2=mi2+1)
{
	freqValueMatrix[mi] = freqValueMatrix[mi] * (1-ps[mi2]);
	freqStrMatrix[mi]   = freqStrMatrix[mi]+"(1-P_"+(mi2+1)+")*";
}


freqStrMatrix[rateCount-1] = freqStrMatrix[rateCount-2]+"(1-P_"+(rateCount-1)+")";
freqStrMatrix[rateCount-2] = freqStrMatrix[rateCount-2]+"P_"+(rateCount-1);
freqValueMatrix[rateCount-1] = 1-scaler;

/* now solve for P_.. */

P_1 = freqValueMatrix[0];
runningTerm = 1-P_1;
for (mi2 = 2; mi2 < rateCount-1; mi2=mi2+1)
{
	ExecuteCommands ("P_"+mi2+"=freqValueMatrix[mi2-1]/runningTerm;runningTerm=runningTerm*(1-P_"+mi2+");");
}
ExecuteCommands ("P_"+(rateCount-1)+"=1-freqValueMatrix[rateCount-1]/runningTerm;");

lfParts = "LikelihoodFunction lf = (filteredData,tree_0";
c_scalerS = "";
for (mi=0; mi<rateCount; mi=mi+1)
{
	if (mi)
	{
		lfDef1 * "+";
		c_scalerS = c_scalerS+"+";
	}
	c_scalerS = c_scalerS + freqStrMatrix[mi] + "*S_"+ mi;
	lfDef1*(""+freqStrMatrix[mi]+"*SITE_LIKELIHOOD["+mi+"]");
	if (mi)
	{
		lfParts = lfParts + ",filteredData,tree_"+mi;
	}
}

lfDef1*0;

ExecuteCommands ("c_scale:="+c_scalerS+";");
ExecuteCommands (lfParts + "," + lfDef1 + ")\");");

timer = Time(0);

LFCompute (lf,LF_START_COMPUTE);
bestDiff = 0;
bestS    = 0;
bestNS   = 0;
bestAlpha = 0;
bestBeta  = 0;

sampleCounter = 0;

if (searchKind == 0 || searchKind == 1)
{
	fprintf (stdout, "\nSampling negatively selected directions\n");

	if (timerCap > 0)
	{
		while (Time(0)-timer <= timerCap)
		{
			v1 = Random(s_rateMin,s_rateMax);
			v2 = Random(ns_rateMin,v1);
			
			checkASample (0);
		}
	}
	else
	{
		step = 1/16;
		v1 = 0;
		for (v1c = 0; v1c < 15; v1c = v1c + 1)
		{
			v1 = v1+step;
			v2 = step/2; 
			for (v2c = 0; v2c < v1c; v2c = v2c+1)
			{
				checkASample (0);	
				v2 = v2 + step;
			}
		} 
	}
}

if (searchKind == 0 || searchKind == 2)
{
	timer = Time(0);
	fprintf (stdout, "\nSampling neutral directions\n");
	if (timerCap > 0)
	{
		while (Time(0)-timer <= timerCap)
		{
			v1 = Random(s_rateMin,s_rateMax);
			v2 = v1;
			
			checkASample (0);
		}
	}
	else
	{
		step = 0.02;
		v1 = 0;
		for (v1c = 0; v1c < 50; v1c = v1c + 1)
		{
			v1 = v1 + step;
			v2 = v1;
			
			checkASample (0);
		}	
	}	
}

if (searchKind == 0 || searchKind == 3)
{
	timer = Time(0);
	fprintf (stdout, "\nSampling positively selected directions\n");
	if (timerCap > 0)
	{
		while (Time(0)-timer <= timerCap)
		{
			v1 = Random(s_rateMin,s_rateMax);
			v2 = Random(v1,ns_rateMax);
			checkASample (0);
		}
	}
	else
	{
		step  = 0.1;
		step2 = 0.25;
		
		v1 = 0.05;
		for (v1c = 0; v1c < 10; v1c = v1c + 1)
		{
			v2 = v1+step; 
			for (v2c = 0; v2c < 20; v2c = v2c+1)
			{
				checkASample (0);	
				v2 = v2 + step2;
			}
		} 
		v1 = v1+step;
	}
}

LFCompute (lf,LF_DONE_COMPUTE);

fprintf (stdout, "\nExamined ", sampleCounter, " directions...\n");

if (bestDiff > 0)
{
	fprintf (stdout, "\nFound a likelihood improvement in the direction (", bestS, ",", bestNS, ")\n");
	slfo = LIKELIHOOD_FUNCTION_OUTPUT;
	LIKELIHOOD_FUNCTION_OUTPUT = 7;
	outFile = outFile+"."+rateCount;
	ExecuteCommands ("S_"+(rateCount-1)+"=bestAlpha;NS_"+(rateCount-1)+"=bestBeta;");
	if (searchKind == 1 || searchKind == 3)
	{
		ExecuteCommands ("global R_"+(rateCount-1)+
						 "=NS_"+(rateCount-1)+
						 "/S_"+(rateCount-1)+";NS_"+(rateCount-1)+":=R_"+(rateCount-1)+"*S_"+(rateCount-1)+";");
		if (searchKind == 1)
		{
			ExecuteCommands ("R_"+(rateCount-1)+":<1.0;\n")
		}
		else
		{
			ExecuteCommands ("R_"+(rateCount-1)+":>1.0;\n")	
		}
	}
	else
	{
		if (searchKind == 2)
		{
			ExecuteCommands ("NS_"+(rateCount-1)+":=S_"+(rateCount-1)+";");
		}
	}
	fprintf(outFile,CLEAR_FILE,lf);
	LIKELIHOOD_FUNCTION_OUTPUT = slfo;
}
else
{
	fprintf (stdout, "\nFailed to find an improvement...\n");
}


/*------------------------------------------------------------------------------------------------------------*/

function checkASample (dummy)
{
	ExecuteCommands ("S_"+(rateCount-1)+"=v1;NS_"+(rateCount-1)+"=v2;");
	LFCompute (lf,res_n);
	if (res_n-res > bestDiff)
	{
		fprintf (stdout, "synRate = ", Format (v1,10,5), " nonSynRate = ", Format (v2,10,5), " Log Likelihood Diff: ", Format (res_n-res,10,5), "\n");
		bestDiff 	= res_n-res;
		bestS 		= v1/c_scale;
		bestNS 		= v2/c_scale;
		bestAlpha 	= v1;
		bestBeta  	= v2;
	}
	sampleCounter = sampleCounter + 1;
	return 0;
}
