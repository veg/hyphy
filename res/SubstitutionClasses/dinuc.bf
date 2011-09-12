TemplateNames = {{
	"One-step",
	"Two-step",
	"1-step Transitions",
	"1-step Transversions"
}};


function CreateTemplates (classID)
{
	MatrixTemplate = {240,1};
	matchCount = 0;

	if (classID < 2)
	{
		for (h=0; h<15; h=h+1)
		{
			for (v=h+1; v<16; v=v+1)
			{
				posCount = (h%4!=v%4)+(h$4!=v$4)-1;
								
				if (posCount == classID)
				{
					MatrixTemplate[matchCount] = h*16 + v;
					matchCount = matchCount+1;
					MatrixTemplate[matchCount] = v*16 + h;
					matchCount = matchCount+1;
				}
			}
		}
	}
	else
	{
		if (classID < 4)
		{
			for (h=0; h<15; h=h+1)
			{
				for (v=h+1; v<16; v=v+1)
				{
					diff = v-h;
					if ((h$4==v$4)||(h%4==v%4))
					{
						if (diff>=4)
						{
							diff = diff$4;
						}
					}
					else
					{
						continue;
					}
					
					if ((diff%2)==classID-2)
					{
						MatrixTemplate[matchCount] = h*16 + v;
						matchCount = matchCount+1;
						MatrixTemplate[matchCount] = v*16 + h;
						matchCount = matchCount+1;
					}
				}
			}
		}
	}
	return 0;
}

SubClassNames = {{
	"Purines-1",
	"Purines-2",
	"Pyrimidines-1",
	"Pyrimidines-2"
}};

function CreateSubClasses (classID)
{
	if (classID==0)
	{
		MatrixTemplate = {{0,1,2,3,8,9,10,11}};
	}
	else
	{
		if (classID==2)
		{
			MatrixTemplate = {{4,5,6,7,12,13,14,15}};
		}
		else
		{
			if (classID==1)
			{
				MatrixTemplate = {{0,2,4,6,8,10,12,14}};
			}
			else
			{
				if (classID==3)
				{
					MatrixTemplate = {{1,3,5,7,9,11,13,15}};
				}
			}
		}
	}
	return 0;
}

