TemplateNames = {{
	"Transitions",
	"Transversions"
}};

SubClassNames = {{
	"Purines",
	"Pyrimidines"
}};

function CreateTemplates (classID)
{
	if (classID==0)
	{	
		MatrixTemplate = {{2,7,8,13}};
	}
	else
	{
		if (classID==1)
		{
			MatrixTemplate = {{1,3,4,6,9,11,12,14}};
		}
	}
	return 0;
}

function CreateSubClasses (classID)
{
	if (classID==0)
	{
		MatrixTemplate = {{0,2}};
	}
	else
	{
		if (classID==1)
		{
			MatrixTemplate = {{1,3}};
		}
	}
	return 0;
}
