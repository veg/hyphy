/* 

THIS HyPhy BATCH LANGUAGE FILE IS PROVIDED AS AN EXAMPLE/TOOL FOR 
THE 2nd EDITION OF 'The Phylogenetic Handbook'

Written by 
Sergei L Kosakovsky Pond (spond@ucsd.edu)
Art FY Poon				 (apoon@ucsd.edu)
Simon DW Frost			 (sdfrost@ucsd.edu)

*/

RequireVersion ("0.9920060901");
ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"chooseGeneticCode.def");

function addACodon (codon, reference, index)
{
	if (codon!=reference)
	{
		diff = Abs (codon-reference);
		if (diff > 16)
		{
			diff = diff $ 16;
		}
		else
		{
			if (diff > 4)
			{
				diff = (diff % 16) $ 4;
			}	
		}
		
		codon = _Genetic_Code[codon];
		
		if (codon != 10)
		{
			if (myCode != codon)
			{
				idx = 2;
			}
			else
			{
				idx = 0;
			}
			if (diff % 2 == 1)
			{
				idx = idx + 1;
			}
			(countMatrix[index])[idx] = (countMatrix[index])[idx] + 1;
		}
		else
		{
			(countMatrix[index])[4] = (countMatrix[index])[4] + 1;
		}
	}
	return 0;
}

nonStopCount = 0;

for (codonIndex=0; codonIndex < 64; codonIndex = codonIndex + 1)
{
	if (_Genetic_Code[codonIndex] != 10)
	{
		nonStopCount = nonStopCount + 1;
	}
}

fprintf (stdout, "This genetic code has ", nonStopCount, " sense (non-termination) codons\n");

countMatrix 	= {};
countMatrix [0] = {5,1};
countMatrix [1] = {5,1};
countMatrix [2] = {5,1};

for (codonIndex=0; codonIndex < 64; codonIndex = codonIndex + 1)
{
	myCode = _Genetic_Code[codonIndex];
	if (myCode != 10)
	{
		/* 1st position */
		remdr = codonIndex%16;
		addACodon (remdr,    codonIndex, 0);
		addACodon (remdr+16, codonIndex, 0);
		addACodon (remdr+32, codonIndex, 0);
		addACodon (remdr+48, codonIndex, 0);
		/* 2nd position */
		remdr = codonIndex - ((codonIndex$4)%4)*4;
		addACodon (remdr,   codonIndex, 1);
		addACodon (remdr+4, codonIndex, 1);
		addACodon (remdr+8, codonIndex, 1);
		addACodon (remdr+12, codonIndex,1);
		/* 3rd position */
		remdr = codonIndex - codonIndex%4;
		addACodon (remdr,   codonIndex, 2);
		addACodon (remdr+1, codonIndex, 2);
		addACodon (remdr+2, codonIndex, 2);
		addACodon (remdr+3, codonIndex, 2);
	}
}

s0 = (countMatrix[2])[0]+(countMatrix[1])[0]+(countMatrix[0])[0];
s1 = (countMatrix[2])[1]+(countMatrix[1])[1]+(countMatrix[0])[1];
s2 = (countMatrix[2])[2]+(countMatrix[1])[2]+(countMatrix[0])[2];
s3 = (countMatrix[2])[3]+(countMatrix[1])[3]+(countMatrix[0])[3];

fprintf (stdout, "\nSubstitution types\n",
				 "                              Synonymous                     Non-synonymous            To a stop codon\n", 
				 "                   Transitions Transversions Total | Transitions Transversions Total |      Total",
				 "\n1st position:    ", Format((countMatrix[0])[0],13,0),
				 					 Format((countMatrix[0])[1],14,0),
				 					 Format((countMatrix[0])[0]+(countMatrix[0])[1],6,0),
				 					 Format((countMatrix[0])[2],14,0),
				 					 Format((countMatrix[0])[3],14,0),
				 					 Format((countMatrix[0])[2]+(countMatrix[0])[3],6,0),
				 					 Format((countMatrix[0])[4],13,0),
				 
				 "\n2nd position:    ", Format((countMatrix[1])[0],13,0),
				 					 Format((countMatrix[1])[1],14,0),
				 					 Format((countMatrix[1])[0]+(countMatrix[1])[1],6,0),
				 					 Format((countMatrix[1])[2],14,0),
				 					 Format((countMatrix[1])[3],14,0),
				 					 Format((countMatrix[1])[2]+(countMatrix[1])[3],6,0),
				 					 Format((countMatrix[1])[4],13,0),

				 "\n3rd position:    ", Format((countMatrix[2])[0],13,0),
				 					 Format((countMatrix[2])[1],14,0),
				 					 Format((countMatrix[2])[0]+(countMatrix[2])[1],6,0),
				 					 Format((countMatrix[2])[2],14,0),
				 					 Format((countMatrix[2])[3],14,0),
				 					 Format((countMatrix[2])[2]+(countMatrix[2])[3],6,0),
				 					 Format((countMatrix[2])[4],13,0),"\n",
				 "                   ------------------------------------------------------------------------------\n",
				 "Total            ", Format(s0,13,0),
				 					  Format(s1,14,0),
				 					  Format(s1+s0,6,0),
				 					  Format(s2,14,0),
				 					  Format(s3,14,0),
				 					  Format(s2+s3,6,0),
				 					  Format((countMatrix[2])[4]+(countMatrix[1])[4]+(countMatrix[0])[4],13,0),"\n");