/* 
0  - A (Alanine)
1  - C (Cysteine)
2  - D (Aspartic Acid)
3  - E (Glutamin Acid)
4  - F (Phenylalanine)
5  - G (Glycine)
6  - H (Histidine)
7  - I (Isoleucine)
8  - K (Lysine)
9  - L (Leucine)
10 - M (Methionine)
11 - N (Asparagine)
12 - P (Prolone)
13 - Q (Glutamine)
14 - R (Arginine)
15 - S (Serine)
16 - T (Threonine)
17 - V (Valine)
18 - W (Tryptophan)
19 - Y (Tyrosine)
*/


SubClassNames = {{
	"Hydrophobic",
	"Hydrophilic",
	"Charged[+]",
	"Charged[-]",
	"Neutral",
	"Polar",
	"Aliphatic",
	"Aromatic",
	"All"
}};

function CreateSubClasses (classID)
{
	if (classID==0) /*hydrophobic*/
	{
		MatrixTemplate = {{0,1,4,7,9,10,12,17,18,19}};
	}
	else
	{
		if (classID==1) /*hydrophilic*/
		{
			MatrixTemplate = {{2,3,6,8,11,13,14,15,16}};
		}
		else
		{
			if (classID==2) /* charged (+) */
			{
				MatrixTemplate = {{6,8,14}};
			}
			else
			{
				if (classID==3) /* charged (-) */
				{
					MatrixTemplate = {{2,3}};
				}
				else
				{
					if (classID==4) /* neutral */
					{
						MatrixTemplate = {{0,1,4,5,7,9,10,11,12,13,15,16,17,18}};
					}
					else
					{
						if (classID==5) /* polar */
						{
							MatrixTemplate = {{1,2,3,6,8,11,13,14,15,16,19}};
						}
						else
						{
							if (classID==6) /* alpihatic */
							{
								MatrixTemplate = {{0,5,7,9,17}};
							}
							else
							{
								if (classID==7) /* aromatic */
								{
									MatrixTemplate = {{4,6,18,19}};
								}
								else
								{
									if (classID==8) /* neutral */
									{
										MatrixTemplate = {{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19}};
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return 0;
}

/*SubClassNames = {{
	"Hydrophobic",
	"Charged",
	"Polar",
	"All"
}};

function CreateSubClasses (classID)
{
	if (classID==0)
	{
		MatrixTemplate = {{0,4,7,9,10,12,17}};
	}
	else
	{
		if (classID==1)
		{
			MatrixTemplate = {{2,3,8,14}};
		}
		else
		{
			if (classID==2)
			{
				MatrixTemplate = {{1,6,11,13,15,16,18,19}};
			}
			else
			{
				if (classID==3)
				{
					MatrixTemplate = {{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19}};
				}
			}
		}
	}
	return 0;
}*/
