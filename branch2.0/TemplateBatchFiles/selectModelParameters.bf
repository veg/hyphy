parameter2Constrain = 0;

if (Rows("LAST_MODEL_PARAMETER_LIST")>1)
{
	ChoiceList (parameter2Constrain, "Parameter(s) to constrain:",1,SKIP_NONE,LAST_MODEL_PARAMETER_LIST);

	if (parameter2Constrain<0)
	{
		return;
	}
}

if (parameter2Constrain==0)
{
	funnyString = "?";
}
else
{
	GetString (funnyString,LAST_MODEL_PARAMETER_LIST,parameter2Constrain-1);
}


constraintString = "this1.?."+funnyString+relationString+"this2.?."+funnyString;
