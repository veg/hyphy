
/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function BuildCodonFrequencies (obsF, freqMx&)
{
	result = {ModelMatrixDimension,1};
	hshift = 0;
	normString = "1";

	for (h=0; h<64; h=h+1)
	{
		first = h$16;
		second = h%16$4;
		third = h%4;
		if (_Genetic_Code[h]==10) 
		{
			hshift = hshift+1;
			normString = normString + "-" + obsF[first][0] + "*" + obsF[second][1]+ "*" +obsF[third][2];
			continue; 
		}
		result[h-hshift]=obsF[first][0] + "*" + obsF[second][1]+ "*" +obsF[third][2];
	}
	
	ExecuteCommands ("global _mlFreqNormalizer := " + normString);
	
	freqMx = {ModelMatrixDimension,1};
	
	hshift = 0;
	for (h=0; h<64; h=h+1)
	{
		first = h$16;
		second = h%16$4;
		third = h%4;
		if (_Genetic_Code[h]==10) 
		{
			hshift = hshift+1;
		}
		else
		{
			ExecuteCommands ("freqMx[h-hshift]:="+result[h-hshift]+"/_mlFreqNormalizer");
		}
	}
	return 0;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function DefineFrequencyParameters (observedFreq)
{
	global AF1:<1;
	global AF2:<1;
	global AF3:<1;
	global C1:<1;
	global C2:<1;
	global C3:<1;
	global G1:<1;
	global G2:<1;
	global G3:<1;
	
	AF1 = observedFreq[0][0];
	AF2 = observedFreq[0][1];
	AF3 = observedFreq[0][2];
	C1 = observedFreq[1][0]/(1-AF1);
	C2 = observedFreq[1][1]/(1-AF2);
	C3 = observedFreq[1][2]/(1-AF3);
	G1 = observedFreq[2][0]/(1-AF1)/(1-C1);
	G2 = observedFreq[2][1]/(1-AF2)/(1-C2);
	G3 = observedFreq[2][2]/(1-AF3)/(1-C3);
	
	global CF1 := (1-AF1)*C1;
	global CF2 := (1-AF2)*C2;
	global CF3 := (1-AF3)*C3;
	
	global GF1 := (1-AF1)*(1-C1)*G1;
	global GF2 := (1-AF2)*(1-C2)*G2;
	global GF3 := (1-AF3)*(1-C3)*G3;

	global TF1 := (1-AF1)*(1-C1)*(1-G1);
	global TF2 := (1-AF2)*(1-C2)*(1-G2);
	global TF3 := (1-AF3)*(1-C3)*(1-G3);

	return  {{"AF1","AF2","AF3"}
		     {"CF1","CF2","CF3"}
		     {"GF1","GF2","GF3"}
			 {"TF1","TF2","TF3"}};

}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

if (!USE_LAST_RESULTS)
{
	global AC = 1;
	global AT = 1;
	global CG = 1;
	global CT = 1;
	global GT = 1;
}
else
{
	global AC;
	global AT;
	global CG;
	global CT;
	global GT;
}

MGCustomRateBiasTerms = {{"AC*","","AT*","CG*","CT*","GT*"}};	

_nucBiasTerms = {4,4};
_nucBiasTerms[0][0] = "";

hv = 0;

for (h=0; h<4; h=h+1)
{
	for (v=h+1; v<4; v=v+1)
	{
		_nucBiasTerms[h][v] = MGCustomRateBiasTerms[hv];
		_nucBiasTerms[v][h] = MGCustomRateBiasTerms[hv];
		hv = hv + 1;	
	}
}

h=0;
v=0;

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function PopulateModelMatrix (ModelMatrixName&, EFV)
{
	if (!ModelMatrixDimension)
	{
		ModelMatrixDimension = 64 - (_Genetic_Code["_MATRIX_ELEMENT_VALUE_ == 10"]*Transpose(_Genetic_Code["1"]))[0];
	}
	
	_localNucBiasMult = _nucBiasTerms;
		
	ModelMatrixName = {ModelMatrixDimension,ModelMatrixDimension}; 
	
	if (modelType >= 0)
	{
		synCatRateMult 	  = "synRate*";
		nonsynCatRateMult = "nonSynRate*";
	
		if (modelType)
		{
			ExecuteCommands ("global R"+_MG94GlobalSuffix+"=1;");
			nonsynCatRateMult = synCatRateMult + "R"+_MG94GlobalSuffix+"*";
			if (modelType > 1)
			{
				synCatRateMult 	      = synCatRateMult    + "c"+_MG94GlobalSuffix+"*";
				nonsynCatRateMult 	  = nonsynCatRateMult + "c"+_MG94GlobalSuffix+"*";
			}
		}
	}
	
	
		
	modelDefString = "";
	modelDefString*16384;
	
	hshift = 0;
	
	for (h=0; h<64; h=h+1)
	{
		if (_Genetic_Code[h]==10) 
		{
			hshift = hshift+1;
			continue; 
		}
		vshift = hshift;
		for (v = h+1; v<64; v=v+1)
		{
			diff = v-h;
			if (_Genetic_Code[v]==10) 
			{
				vshift = vshift+1;
				continue; 
			}
			nucPosInCodon = 2;
			if ((h$4==v$4)||((diff%4==0)&&(h$16==v$16))||(diff%16==0))
			{
				if (h$4==v$4)
				{
					transition = v%4;
					transition2= h%4;
				}
				else
				{
					if(diff%16==0)
					{
						transition = v$16;
						transition2= h$16;
						nucPosInCodon = 0;
					}
					else
					{
						transition = v%16$4;
						transition2= h%16$4;
						nucPosInCodon = 1;
					}
				}
				hs = Format(h-hshift,0,0);
				vs = Format(v-vshift,0,0);
				ts = Format(transition,0,0);
				ts2= Format(transition2,0,0);
				ps = Format(nucPosInCodon,0,0);
				aa1 = _Genetic_Code[0][h];
				aa2 = _Genetic_Code[0][v];
				
				if (aa1==aa2) 
				{
					modelDefString*("ModelMatrixName["+hs+"]["+vs+"] := "+synCatRateMult+_localNucBiasMult[transition][transition2]+EFV[transition][nucPosInCodon]+";\n"+
													 "ModelMatrixName["+vs+"]["+hs+"] := "+synCatRateMult+_localNucBiasMult[transition][transition2]+EFV[transition2][nucPosInCodon]+";\n");
				}
				else
				{
					modelDefString*("ModelMatrixName["+hs+"]["+vs+"] := "+nonsynCatRateMult+_localNucBiasMult[transition][transition2]+EFV[transition][nucPosInCodon]+";\n"+
													 "ModelMatrixName["+vs+"]["+hs+"] := "+nonsynCatRateMult+_localNucBiasMult[transition][transition2]+EFV[transition2][nucPosInCodon]+";\n");						
				}
			}
	    }
    }		
	modelDefString*0;
	ExecuteCommands (modelDefString);

	if (Abs(MGCustomModelConstraintString))
	{
		ExecuteCommands (MGCustomModelConstraintString);
	}
	return 0;
}
