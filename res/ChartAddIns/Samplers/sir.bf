/*---------------------------------------------------------------------------------------------------------------------*/

pi_const = 3.141592653589793;

/*---------------------------------------------------------------------------------------------------------------------*/

function x2ymap (x,l,u, m)
{	
	if (m==1)
		return Log(x-l);
	else
	{
		x = (x-l)/(u-l);
		return Log(x/(1-x));
	}
}

/*---------------------------------------------------------------------------------------------------------------------*/

function  y2xmap (y,l,u, m)
{
	if (m==1)
		return Exp(y)+l;
	else
	{
		x = Exp(y)/(1+Exp(y));
		return x*(u-l)+l;
	}
}

/*---------------------------------------------------------------------------------------------------------------------*/

function CholeskyDecomposition (aMatrix)
{
	matrixDim = Rows(aMatrix);
	_CholeskyL_ 	 =     {matrixDim,matrixDim};
	for (i=0; i<matrixDim; i=i+1)
	{
		for (j=i; j<matrixDim; j=j+1)
		{
			k = i-1;
			for (sum=aMatrix[i][j]; k>=0; k=k-1)
			{
				sum = sum - aMatrix[i][k]*aMatrix[j][k];
			}
			if (i==j)
			{
				if (sum<=0)
				{
					fprintf (stdout, "\nMatrix passed to CholeskyDecomposition was not positive definite\n");
					return 0;
				}
				_CholeskyL_[i][i] = Sqrt (sum);
			}
			else
			{
				_CholeskyL_[j][i] = sum/_CholeskyL_[i][i];
			}
		}
	}
	return _CholeskyL_;
}

/*---------------------------------------------------------------------------------------------------------------------*/

varCount = Abs(COVARIANCE_PARAMETER);
stashedValues = {varCount,4};
MeanVector    = {varCount,1};

usedVars = Rows(COVARIANCE_PARAMETER);


smif = USE_INTERVAL_MAPPING;
svpc = COVARIANCE_PRECISION;
USE_INTERVAL_MAPPING = 1;
COVARIANCE_PRECISION = 2;

ExecuteCommands ("CovarianceMatrix (covMx, "+LF_NAME+");");

COVARIANCE_PRECISION = svpc;
USE_INTERVAL_MAPPING = smif;

assignmentString = "";
assignmentString * 256;

for (k=0; k<varCount; k=k+1)
{
	aKey = usedVars[k];
	GetInformation(varInfo,aKey);
	stashedValues[k][0] = varInfo[0];
	stashedValues[k][1] = varInfo[1];
	stashedValues[k][2] = varInfo[2];
	stashedValues[k][3] = INTERVAL_MAPPING_METHOD[aKey];
	MeanVector[k] = x2ymap(varInfo[0],varInfo[1],varInfo[2],stashedValues[k][3]);
	assignmentString * (aKey+ "=generatedSamples[itCount]["+k+"];\n");
}

assignmentString * 0;

varCount	= Rows (covMx);
generatedSamples = {SAMPLE_N,varCount};

LT = CholeskyDecomposition(covMx);


if (Abs(LT))
{
	fprintf (stdout, "Generating approximate normal samples...\n");
	timer = Time (0);
	for (itCount = 0; itCount < SAMPLE_N; itCount = itCount + 1)
	{
		sampleVector = {varCount,1};
		
		for (coord = 0; coord < varCount; coord = coord+1)
		{
			_cdfCut = Random(0,1);
			FindRoot (stdNormal,ZCDF(x)-_cdfCut,x,-10,10);
			sampleVector[coord] = stdNormal;
		}
		
		sampleVector = LT*sampleVector+MeanVector;

		for (coord = 0; coord < varCount; coord = coord+1)
		{
			generatedSamples[itCount][coord] = y2xmap(sampleVector[coord],stashedValues[coord][1],stashedValues[coord][2],stashedValues[coord][3]);
		}
		
		if ((1+itCount) % 500 == 0)
		{
			fprintf (stdout, itCount+1, "/", SAMPLE_N, " samples done. Estimated remaining time: ",Format (((SAMPLE_N-itCount-1)/(itCount+1))*(Time(0)-timer),5,2)," seconds \n");
		}
	}
	
    #include "srs.ibf";
}
