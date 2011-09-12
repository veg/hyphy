likelihoodFnChoice = 0;

if (Rows("LikelihoodFunction")>1)
{
	ChoiceList  (likelihoodFnChoice,"Choose a Likelihood Function",1,NO_SKIP,LikelihoodFunction);
}		

if (likelihoodFnChoice<0)
{
	return;
} 

GetString (likelihoodFunctionName,LikelihoodFunction,likelihoodFnChoice);

if (RESTORE_GLOBALS)
{
	dumb = RestoreGlobalValues (likelihoodFnChoice);
}
