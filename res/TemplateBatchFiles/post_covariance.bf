likelihoodFnChoice = 0;

if (Rows("LikelihoodFunction")>1)
{
	ChoiceList  (likelihoodFnChoice,"Choose a Likelihood Function",1,NO_SKIP,LikelihoodFunction);
}		

if (likelihoodFnChoice<0)
{
	return;
} 
			 
ChoiceList  (response,"Covariance Matrix",1,NO_SKIP,
			 "Crude covariances","1st order derivative approximation. This could be slow, since it entails ~p(p-1)/2 calculations of the likelihood function, where p = # independent parameters.",
			 "Finer covariances","2nd order derivative approximation. This could be very slow, since it entails ~p(p-1)*2 calculations of the likelihood function, where p = # of independent parameters.",
			 "Likelihood Profile","95% CI based on the quadratic approximation to the likelihood surface. Should be faster than full covariance matrices for large numbers of parameters and may be more robust to statistical errors.");

			 
if (response<0)
{
	return;
}

if (response == 2)
{
	COVARIANCE_PRECISION = 0.95;
}
else
{	
	COVARIANCE_PRECISION = response+1;
}

GetString (likelihoodFunctionName,LikelihoodFunction,likelihoodFnChoice);

if (RESTORE_GLOBALS)
{
	RestoreGlobalValues (likelihoodFnChoice);
}

CovarianceMatrix (covMatrix, *likelihoodFunctionName);

SetDialogPrompt ("Save Covariance Matrix to:");

fprintf (PROMPT_FOR_FILE,CLEAR_FILE,"Order of parameters:\n");

covDim = Rows(covMatrix);

for (covCounter = 0; covCounter<covDim; covCounter=covCounter+1)
{
	GetString (argName,*likelihoodFunctionName,covCounter);
	fprintf (LAST_FILE_PATH,"\t",argName);
}

fprintf (LAST_FILE_PATH,covMatrix);
