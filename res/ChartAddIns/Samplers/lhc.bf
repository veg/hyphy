/*---------------------------------------------------------------------------------------------------------------------*/

varCount = Abs(COVARIANCE_PARAMETER);
usedVars = Rows(COVARIANCE_PARAMETER);

inflationFactor = -1;
while (inflationFactor <= 0)
{
	fprintf (stdout, "Inflation factor for profile likelihood bounds (>0.0):");
	fscanf (stdin,"Number",inflationFactor);
}

fprintf (stdout, "\nObtaining profile likeihood bounds...\n");

svpc = COVARIANCE_PRECISION;
COVARIANCE_PRECISION = 0.95;

ExecuteCommands ("CovarianceMatrix (covMx, "+LF_NAME+");");

COVARIANCE_PRECISION = svpc;
stashedValues = {varCount,4};

assignmentString = "";
assignmentString * 256;

for (k=0; k<varCount; k=k+1)
{
	aKey = usedVars[k];
	stashedValues[k][0] = covMx[k][1];
	stashedValues[k][1] = covMx[k][1]-(covMx[k][1]-covMx[k][0])*inflationFactor;
	if (stashedValues[k][1] < 0)
	{
		stashedValues[k][1] = 0;
	}
	stashedValues[k][2] = covMx[k][1]+(covMx[k][2]-covMx[k][1])*inflationFactor;
	stashedValues[k][3] = (covMx[k][2]-covMx[k][0])*inflationFactor/SAMPLE_N;
	assignmentString * (aKey+ "=generatedSamples[itCount]["+k+"];\n");
}

assignmentString * 0;

varCount	= Rows (covMx);
generatedSamples = {SAMPLE_N,varCount};
indexVector		 = {1,SAMPLE_N};

for (k=0; k<SAMPLE_N; k=k+1)
{
	indexVector[k] = k;
}

fprintf (stdout, "\nDoing Latin hypercube sampling...\n");

for (k=0; k<varCount; k=k+1)
{
	permVector = Random (indexVector,1);
	lb = stashedValues[k][1];
	st = stashedValues[k][3];
	for (m = 0; m<SAMPLE_N; m=m+1)
	{
		generatedSamples[m][k] = lb+st*permVector[m];
	}
}


#include "srs.ibf";

