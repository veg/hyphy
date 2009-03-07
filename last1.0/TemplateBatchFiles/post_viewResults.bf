likelihoodFnChoice = 0;
if (Rows("LikelihoodFunction")>1)
{
	ChoiceList  (likelihoodFnChoice,"Choose a Likelihood Function",1,NO_SKIP,LikelihoodFunction);
}		

if (likelihoodFnChoice<0)
{
	return;
} 
			 
ChoiceList  (response,"Likelihood Function Display",1,NO_SKIP,
			 "Value Only","Display Only the ln-likelihood value.",
			 "List","Display the ln-likelihood value followed by a list of parameter estimates.",
			 "Standard Tree","Display the ln-likelihood value followed by a tree string with branch lengths representing the expected number of substitutions per site.",
			 "Parameters only","Display the list of independent parameters and constraints.",
			 "Parameters+Values","Display the list of parameter names, constraints and MLEs.",
			 "Parameters+Values+Trees","Display the list of parameter names, constraints and MLEs followed by tree string(s) with branch lengths representing the expected number of substitutions per site.");
			 
			 
if (response<0)
{
	return;
}

LIKELIHOOD_FUNCTION_OUTPUT = response;

if (RESTORE_GLOBALS)
{
	dumb = RestoreGlobalValues (likelihoodFnChoice);
}

GetString (likelihoodFunctionName,LikelihoodFunction,likelihoodFnChoice);
fprintf   (stdout, "\n", likelihoodFunctionName__ ,"\n");
