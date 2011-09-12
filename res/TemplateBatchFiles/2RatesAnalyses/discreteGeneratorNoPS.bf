categDef1 = "";
categDef2 = "";

categDef1 * 1024;
categDef2 * 1024;

defTail = "";

for (mi=1; mi<resp; mi=mi+1)
{
	if (randomizeInitValues)
	{
		categDef1*("global PS_"+mi+" = Random(0.05,0.95);\nPS_"+mi+":<1;\n");	
	}
	else
	{
		categDef1*("global PS_"+mi+" = 1/"+(resp+1-mi)+";\nPS_"+mi+":<1;\n");
	}
}

if (randomizeInitValues)
{
	categDef1*( "\n\nglobal RS_1 = Random(0.05,0.95);\nRS_1:<1;RS_1:>0.000000001;\n");
}
else
{
	categDef1*( "\n\nglobal RS_1 = .3;\nRS_1:<1;RS_1:>0.000000001;\n");
}

for (mi=3; mi<=resp; mi=mi+1)
{
	if (randomizeInitValues)
	{
		categDef1*("global RS_" + mi + " = Random (1.05,5);\nRS_" + mi + ":>1;RS_" + mi + ":<100000;\n");	
	}
	else
	{
		categDef1*("global RS_" + mi + " = 1.5;\nRS_" + mi + ":>1;RS_" + mi + ":<100000;\n");
	}
} 

rateStrMx    = {resp,1};
rateStrMx[0] = "RS_1";
rateStrMx[1] = "1";

for (mi=3; mi<=resp; mi=mi+1)
{
	rateStrMx[mi-1] = rateStrMx[mi-2]+"*RS_"+mi;
} 	

freqStrMx    = {resp,1};
freqStrMx[0] = "PS_1";

for (mi=1; mi<resp-1; mi=mi+1)
{
	freqStrMx[mi] = "";
	for (mi2=1;mi2<=mi;mi2=mi2+1)
	{
		freqStrMx[mi] = freqStrMx[mi]+"(1-PS_"+mi2+")";		
	}
	freqStrMx[mi] = freqStrMx[mi]+"PS_"+(mi+1);	
}	

freqStrMx[mi] = "";
for (mi2=1;mi2<mi;mi2=mi2+1)
{
	freqStrMx[mi] = freqStrMx[mi]+"(1-PS_"+mi2+")";		
}
freqStrMx[mi] = freqStrMx[mi]+"(1-PS_"+mi+")";	


categDef1*( "\n\nglobal c_scale:=" + rateStrMx[0] + "*" + freqStrMx[0]);

for (mi=1; mi<resp; mi=mi+1)
{
	categDef1*( "+" + rateStrMx[mi] + "*" + freqStrMx[mi]);
}

categDef1*( ";\ncategFreqMatrix={{"+freqStrMx[0] );

for (mi=1; mi<resp; mi=mi+1)
{
	categDef1*( "," + freqStrMx[mi]);
}

categDef1*( "}};\ncategRateMatrix={{" + rateStrMx[0] + "/c_scale");

for (mi=1; mi<resp; mi=mi+1)
{
	categDef1*( "," + rateStrMx[mi] + "/c_scale");
}

categDef1*( "}};\n\ncategory c      = (" + resp + ", categFreqMatrix , MEAN, ,categRateMatrix, 0, 1e25"+defTail+");\n\n");

/* begin non-syn */

if (randomizeInitValues)
{
	categDef2*( "\n\nglobal RN_1 = Random(0.05,0.95);\nRN_1:<1;\n");
}
else
{
	categDef2*( "\n\nglobal RN_1 = .0.9;\nRN_1:<1;\n");
}


for (mi=2; mi<=resp2; mi=mi+1)
{
	if (randomizeInitValues)
	{
		categDef2*("global RN_" + mi + " = Random(0.05,0.95);\nRN_" + mi + ":<1;\n");
	}
	else
	{
		categDef2*("global RN_" + mi + " = 0.75;\nRN_" + mi + ":<1;\n");
	}
} 

multFactor = rateStrMx[0];

rateStrMx    	   = {resp2,1};
rateStrMx[0] = "";
rateStrMx[resp2-1] = multFactor+"/c_scale*RN_" + resp2;

for (mi=resp2-1; mi>=1; mi=mi-1)
{
	rateStrMx[mi-1] = rateStrMx[mi]+"*RN_"+mi;
} 	


for (mi=1; mi<=resp2; mi=mi+1)
{
	if (randomizeInitValues)
	{
		categDef2*("global PN_" + mi + " = Random(0.05,0.95);\nPN_" + mi + ":<1;\n");	
	}
	else
	{
		categDef2*("global PN_" + mi + " = 1/" + (resp2+1-mi) + ";\nPN_" + mi + ":<1;\n");	
	}
}

freqStrMx    = {resp2,1};
freqStrMx[0] = "PN_1";

for (mi=1; mi<resp2-1; mi=mi+1)
{
	freqStrMx[mi] = "";
	for (mi2=1;mi2<=mi;mi2=mi2+1)
	{
		freqStrMx[mi] = freqStrMx[mi]+"(1-PN_"+mi2+")";		
	}
	freqStrMx[mi] = freqStrMx[mi]+"PN_"+(mi+1);	
}	

freqStrMx[mi] = "";
for (mi2=1;mi2<mi;mi2=mi2+1)
{
	freqStrMx[mi] = freqStrMx[mi]+"(1-PN_"+mi2+")";		
}
freqStrMx[mi] = freqStrMx[mi]+"(1-PN_"+mi+")";	


categDef2*( "\ncategRateMatrixN={{" + rateStrMx[0]);

for (mi=1; mi<resp2; mi=mi+1)
{
	categDef2*( "," + rateStrMx[mi]);
}

categDef2*( "}};\n");

categDef2*( "\n\ncategFreqMatrixN={{" + freqStrMx[0] );

for (mi=1; mi<resp2; mi=mi+1)
{
	categDef2*( "," + freqStrMx[mi]);
}

categDef2*( "}};\n\ncategory d     = (" + resp2 + ", categFreqMatrixN , MEAN, ,categRateMatrixN, 0, 1e25"+defTail+");\n\n");


categDef1 * 0;
categDef2 * 0;
