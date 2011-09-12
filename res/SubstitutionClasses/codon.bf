ModelGeneticCode = {
		{14, 13, 14,  13,
		  7,  7,  7,   7, 
		 19,  5, 19,   5, 
		  2,  2, 	3,   2, 
		 12, 11, 12,  11,
		  6,  6,  6,   6, 
		 19, 19, 19,  19,
		  1,  1,  1,   1, 
		 16, 15, 16,  15,
		  8,  8,  8,   8, 
		 20, 20, 20,  20,
		  4,  4,  4,   4, 
		 10,  9, 10,   9, 
		  5,  5,  5,   5, 
		 10, 17, 18,  17,
		  1,  0,  1,   0 }
};

ModelGeneticCodeSize = 64;

for (h=0; h<64; h=h+1)
{
	if (ModelGeneticCode[h]==10)
	{
		ModelGeneticCodeSize=ModelGeneticCodeSize-1;
	}
}

TemplateNames = {{
	"One-step",
	"Two-step",
	"Three-step",
	"Synonymous",
	"Non-synonymous",
	"1-step Transitions",
	"1-step Transversions",
	"A-C",
	"A-G",
	"A-T",
	"C-G",
	"C-T",
	"G-T"
}};

SubClassNames = {{
	"Alanin",
	"Arginine",
	"Asparagine",
	"Aspartic Acid",
	"Cysteine",
	"Glutamin Acid",
	"Glutamine",
	"Glycine",
	"Histidine",
	"Isoleucine",
	"Leucine",
	"Lysine",
	"Methionine",
	"Phenylalanine",
	"Proline",
	"Serine",
	"Threonine",
	"Tryptophan",
	"Tyrosine",
	"Valine",
	"Hydrophobic",
	"Charged",
	"Polar"
}};

function CreateTemplates (classID)
{
	vShift = 0;
	hShift = 0;
	MatrixTemplate = {ModelGeneticCodeSize*ModelGeneticCodeSize,1};
	matchCount = 0;

	if (classID < 3)
	{
		for (h=0; h<63; h=h+1)
		{
			if (ModelGeneticCode[h]==10)
			{
				hShift = hShift+1;
				continue;
			}
			
			vShift = hShift;
			
			for (v=h+1; v<64; v=v+1)
			{
				if (ModelGeneticCode[v]==10)
				{
					vShift = vShift+1;
					continue;
				}
				
				posCount = (h%4!=v%4)+(h$16!=v$16)+(h%16$4!=v%16$4)-1;
								
				if (posCount == classID)
				{
					MatrixTemplate[matchCount] = (h-hShift)*ModelGeneticCodeSize + v-vShift;
					matchCount = matchCount+1;
					MatrixTemplate[matchCount] = (v-vShift)*ModelGeneticCodeSize + h-hShift;
					matchCount = matchCount+1;
				}
			}
		}
	}
	else
	{
		if (classID<5)
		{
			classID = 4-classID;
			for (h=0; h<63; h=h+1)
			{
				if (ModelGeneticCode[h]==10)
				{
					hShift = hShift+1;
					continue;
				}
				vShift = hShift;
				for (v=h+1; v<64; v=v+1)
				{
					if (ModelGeneticCode[v]==10)
					{
						vShift = vShift+1;
						continue;
					}
					if ((ModelGeneticCode[v]==ModelGeneticCode[h])==classID)
					{
						MatrixTemplate[matchCount] = (h-hShift)*ModelGeneticCodeSize + v-vShift;
						matchCount = matchCount+1;
						MatrixTemplate[matchCount] = (v-vShift)*ModelGeneticCodeSize + h-hShift;
						matchCount = matchCount+1;
					}
				}
			}
		}
		else
		{
			if (classID < 7)
			{
				for (h=0; h<63; h=h+1)
				{
					if (ModelGeneticCode[h]==10)
					{
						hShift = hShift+1;
						continue;
					}
					
					vShift = hShift;
					for (v=h+1; v<64; v=v+1)
					{
						if (ModelGeneticCode[v]==10)
						{
							vShift = vShift+1;
							continue;
						}
						diff = v-h;
						if ((h$4==v$4)||((diff%4==0)&&(h$16==v$16))||(diff%16==0))
						{
							if (diff>=4)
							{
								if (diff<16)
								{
									diff = diff$4;
								}
								else
								{
									diff = diff$16;
								}
							}
						}
						else
						{
							continue;
						}
						
						if ((diff%2)==classID-5)
						{
							MatrixTemplate[matchCount] = (h-hShift)*ModelGeneticCodeSize + v-vShift;
							matchCount = matchCount+1;
							MatrixTemplate[matchCount] = (v-vShift)*ModelGeneticCodeSize + h-hShift;
							matchCount = matchCount+1;
						}
					}
				}
			}
			else
			{
				if (classID < 13)
				{
					
					m11 = 0;
					if (classID>=10)
					{
						if (classID<12)
						{
							m11 = 1;
						}
						else
						{
							m11 = 2;
						}
					}
					
					m22 = 3;
					if (classID==7)
					{
						m22 = 1;
					}
					else
					{
						if ((classID==8)||(classID==10))
						{
							m22 = 2;
						}
					}
					
					
					for (h=0; h<63; h=h+1)
					{
						if (ModelGeneticCode[h]==10)
						{
							hShift = hShift+1;
							continue;
						}
						
						vShift = hShift;
						for (v=h+1; v<64; v=v+1)
						{
							if (ModelGeneticCode[v]==10)
							{
								vShift = vShift+1;
								continue;
							}
							
							diff = v-h;
							p11  = h%4;
							p21  = v%4;
							p12  = (h%16)$4;
							p22  = (v%16)$4;
							p13  = h$16;
							p23  = v$16;
							
							if (((p11==m11)&&(p21==m22))||((p12==m11)&&(p22==m22))||((p13==m11)&&(p23==m22)))
							{
								MatrixTemplate[matchCount] = (h-hShift)*ModelGeneticCodeSize + v-vShift;
								matchCount = matchCount+1;
								MatrixTemplate[matchCount] = (v-vShift)*ModelGeneticCodeSize + h-hShift;
								matchCount = matchCount+1;
							}
						}
					}
				}
			}
		}
		
	}
	return 0;
}

NumberToIndex = {{
8,
19,
13,
15,
17,
16,
12,
20,
11,
2,
1,
14,
3,
0,
6,
5,
7,
18,
9,
4
}};

function CreateSubClasses (classID)
{
	hShift = 0;
	MatrixTemplate = {ModelGeneticCodeSize,1};
	matchCount = 0;
	if (classID<20)
	{
		v = NumberToIndex [classID];
		for (h=0; h<64; h=h+1)
		{
			if (ModelGeneticCode[h-hShift]==10)
			{
				hShift = hShift+1;
				continue;
			}
			if (v==ModelGeneticCode[h-hShift])
			{
				MatrixTemplate[matchCount] = h-hShift;
				matchCount = matchCount+1;
			}
		}
	}
	else
	{
		if (classID==20)
		{
			for (h=0; h<64; h=h+1)
			{
				vShift = ModelGeneticCode[h-hShift];
				if (vShift==10)
				{
					hShift = hShift+1;
					continue;
				}
				if ((vShift==8)||(vShift==4)||(vShift==0)||(vShift==6)||(vShift==2)||(vShift==1))
				{
					MatrixTemplate[matchCount] = h-hShift;
					matchCount = matchCount+1;
				}
			}
		
		}
		else
		{
			if (classID==21)
			{
				for (h=0; h<64; h=h+1)
				{
					vShift = ModelGeneticCode[h-hShift];
					if (vShift==10)
					{
						hShift = hShift+1;
						continue;
					}
					if ( (vShift==13)||(vShift==16)||(vShift==14)||(vShift==19))
					{
						MatrixTemplate[matchCount] = h-hShift;
						matchCount = matchCount+1;
					}
				}
			
			}
			else
			{
				if (classID==22)
				{
					for (h=0; h<64; h=h+1)
					{
						vShift = ModelGeneticCode[h-hShift];
						if (vShift==10)
						{
							hShift = hShift+1;
							continue;
						}
						if ( (vShift==5)||(vShift==7)||(vShift==9)||(vShift==11)||(vShift==17)||(vShift==13)||(vShift==12)||(vShift==18))
						{
							MatrixTemplate[matchCount] = h-hShift;
							matchCount = matchCount+1;
						}
					}
				
				}
			}
		}
	}
	return 0;
}
