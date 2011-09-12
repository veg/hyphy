likelihoodFnChoice = 0;
if (Rows("LikelihoodFunction")>1)
{
	ChoiceList  (likelihoodFnChoice,"Choose a Likelihood Function",1,NO_SKIP,LikelihoodFunction);
}		

if (likelihoodFnChoice<0)
{
	return;
} 
			 
LIKELIHOOD_FUNCTION_OUTPUT = response;

if (RESTORE_GLOBALS)
{
	dumb = RestoreGlobalValues (likelihoodFnChoice);
}

GetString(likelihoodFunctionName,LikelihoodFunction,likelihoodFnChoice);

DataSet 		ancestralSequences = ReconstructAncestors(likelihoodFunctionName__);

DataSetFilter 	ancestralSequencesDF = CreateFilter (ancestralSequences,1);

SetDialogPrompt ("Save ancestral sequences to:");

fprintf (PROMPT_FOR_FILE,CLEAR_FILE,ancestralSequencesDF);
