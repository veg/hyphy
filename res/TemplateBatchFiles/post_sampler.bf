likelihoodFnChoice = 0;

if (Rows("LikelihoodFunction")>1)
{
	ChoiceList  (likelihoodFnChoice,"Choose a Likelihood Function",1,NO_SKIP,LikelihoodFunction);
}		

if (likelihoodFnChoice<0)
{
	return;
} 
	
GetString (LF_NAME,LikelihoodFunction,likelihoodFnChoice);

ChoiceList  (sMethod,"Sampling Method",1,NO_SKIP,
			 "Normal + SiR","Sample parameters using multivariate normal with the estimated variance-covariance matrix",
			 "Latin Hypercube + SiR","Sample parameters using latin hypercube sampling scheme with intervals derived using profile likelihood (possibly inflated)");
			 
if (sMethod<0)
{
	return;
}

ExecuteCommands("GetString (LF_VARS,"+LF_NAME+",-1);");

catVars  = LF_VARS["Categories"];
globVars = LF_VARS["Global Independent"];
locVars  = LF_VARS["Local Independent"];

lvc = Columns (locVars);
gvc = Columns (globVars);

sMarginals = 0;

SetDialogPrompt ("Save results to:");
fprintf (PROMPT_FOR_FILE,CLEAR_FILE);
baseResPath = LAST_FILE_PATH;

if (Columns(catVars)>0)
{
	ChoiceList  (sMarginals,"Save marginals",1,NO_SKIP,
				 "Skip","Do not save site-by-site conditional probabilities",
				 "Save","Save site-by-site conditional probabilities to a file for later processing",
				 "Save all","Save even the samples which were rejected");
				 
	if (sMarginals < 0)
	{
		return;
	}	
	else
	{
		if (sMarginals)
		{
			marginalOutFile = baseResPath+".marginals";
			fprintf (marginalOutFile,CLEAR_FILE);
			ExecuteCommands ("GetInformation (catVarList,"+LF_NAME+");");
			if (sMarginals > 1)
			{
				marginalOutFileAll = baseResPath+".marginals.all";
				marginalOutFileLF  = baseResPath+".marginals.all.logL";
				fprintf (marginalOutFileAll,CLEAR_FILE);
				fprintf (marginalOutFileLF ,CLEAR_FILE);
			}
		}
	}
}

if (lvc && gvc)
{
ChoiceList  (pMethod,"Parameters to Sample",1,NO_SKIP,
			 "All","Sample all independent parameters",
			 "Global","Sample global model parameters only",
			 "Local","Sample local model parameter only",
			 "Choose","Choose which parameters to sample",
			 "RegExp","Sample model parameters which match a regular expression");
			 
	if (pMethod < 0) {
		return;
	}
}
else
{
	if (Columns(globVars))
	{
		ChoiceList  (pMethod,"Parameters to Sample",1,NO_SKIP,
					 "All","Sample all independent parameters",
					 "Global","Sample global model parameters only",
					 "Choose","Choose which parameters to sample",
			         "RegExp","Sample model parameters which match a regular expression");
		if (pMethod < 0) {
			return;
		}
		else {
			if (pMethod >= 2) {
				pMethod += 1;
			}
		}
	}
	else
	{
		ChoiceList  (pMethod,"Parameters to Sample",1,NO_SKIP,
					 "All","Sample all independent parameters",
					 "Local","Sample locl model parameters only",
					 "Choose","Choose which parameters to sample",
			         "RegExp","Sample model parameters which match a regular expression");
		if (pMethod < 0) {
			return;
		}
		else {
			if (pMethod) {
				pMethod += 1;
			}
		}	
	}
}

COVARIANCE_PARAMETER = {};

choiceList  = {lvc + gvc, 2};
chosenParms = {lvc + gvc, 1}; 

for (k=0; k < gvc; k=k+1)
{
	choiceList[k][0] = globVars[k];
	choiceList[k][1] = "Global variable " + globVars[k];
}

for (k=0; k < lvc; k=k+1)
{
	choiceList[gvc+k][0] = locVars[k];
	choiceList[gvc+k][1] = "Local variable " + locVars[k];
}

if (pMethod == 3)
{
	ChoiceList  (vMethod,"Parameters to Sample",0,NO_SKIP,choiceList);
	if (vMethod[0] < 0)
	{
		return;
	}
	gvc = Columns(vMethod);
	for (k=0; k<gvc; k=k+1)
	{
		lvc = vMethod[k];
		chosenParms[lvc] = 1;
	}
}
else
{
	if (pMethod == 0) {
	    chosenParms = chosenParms["1"];
	}
	else {
		if (pMethod == 1) {
			for (k=0; k<gvc; k=k+1) {
				chosenParms[k] = 1;
			}
		}
		else {
		    if (pMethod == 4) {
		        fprintf (stdout, "Please enter a regular expression to filter the parameter names with:");
		        fscanf  (stdin,"String",filtering_regexp);   
		        chosenParms = chosenParms["(choiceList[_MATRIX_ELEMENT_ROW_][0]$filtering_regexp)[0]>=0"];
		        for (k = 0; k < Rows (chosenParms); k+=1) {
		            if (chosenParms[k]) {
		                fprintf (stdout, "Selected ", choiceList[k][0], " for sampling\n");
		            }
		        }
		    } else {
		        chosenParms = chosenParms["1"];
			}
		}
	}
}

for (k=Columns(locVars)+Columns(globVars)-1; k>=0; k=k-1)
{
	if (chosenParms[k])
	{
		lvc = choiceList[k][0];
		COVARIANCE_PARAMETER [lvc] = 1;
	}
}

SAMPLE_N = 0;

while (SAMPLE_N<1)
{
	fprintf (stdout, "How many samples should be drawn? (>=1)");
	fscanf  (stdin,"Number",SAMPLE_N);
}	

SAMPLE_N = SAMPLE_N$1;
SAMPLE_M = 0;

while (SAMPLE_M<1 || SAMPLE_M > SAMPLE_N)
{
	fprintf (stdout, "How many re-samples should be drawn? [1,", SAMPLE_N, "]");
	fscanf  (stdin,"Number",SAMPLE_M);
}	
SAMPLE_M = SAMPLE_M$1;

if (sMethod)
{
	ExecuteAFile("Samplers/lhc.bf");
}
else
{
	ExecuteAFile("Samplers/sir.bf");
}
