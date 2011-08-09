/*---------------------------------------------------------------------------------------------------------------------*/
LoadFunctionLibrary ("GrabBag");

varCount = Abs (COVARIANCE_PARAMETER);
usedVars = Rows(COVARIANCE_PARAMETER);

inflationFactor = prompt_for_a_value ("Inflation factor for profile likelihood bounds", 2.0, 0.0, 1000.0, 0);

fprintf (stdout, "\nObtaining profile likeihood bounds...\n");

svpc = COVARIANCE_PRECISION;
COVARIANCE_PRECISION = 0.95;

useMPIFlag = MPI_NODE_COUNT > 2 && MPI_NODE_ID == 0 ;

if (useMPIFlag && varCount > 1)
{
	availableComputeNodes = Min(MPI_NODE_COUNT - 1,varCount);
	
	currentIndex = 0;
	
	mapper		 = {};
	dispatched   = {};
	
	for (k = 0; k < varCount; k+=1)
	{
		mapper [usedVars[k]] = k;
	}
	
	for (nodeIndex = 0; nodeIndex < availableComputeNodes; nodeIndex += 1)
	{
		localCovarianceParameters = {};
		for (currentIndex = nodeIndex; currentIndex < varCount; currentIndex += availableComputeNodes)
		{
			localCovarianceParameters[usedVars[currentIndex]] = 1;
		}
		
		LF_NEXUS_EXPORT_EXTRA = "COVARIANCE_PARAMETER = " + localCovarianceParameters + 
							  ";\nCOVARIANCE_PRECISION = 0.95;\nCovarianceMatrix (MPI_NEXUS_FILE_RETURN, `LF_NAME`);fprintf (stdout, MPI_NEXUS_FILE_RETURN);";
		
		ExecuteCommands ("Export (lfString, `LF_NAME`);");
		dispatched [nodeIndex+1] = localCovarianceParameters;
		MPISend (nodeIndex+1, lfString);
	}
	
	covMx = {varCount, 3};
	
	for (nodeIndex = 0; nodeIndex < availableComputeNodes; nodeIndex += 1)
	{
		MPIReceive  (-1, fromNode, res);
		Eval ("nodeMx="+res);
		varNames = Rows(dispatched [fromNode]);
		for (k = 0; k < Columns (varNames); k+=1)
		{
			k2			 = mapper[varNames[k]];
			covMx[k2][0] = nodeMx[k][0];
			covMx[k2][1] = nodeMx[k][1];
			covMx[k2][2] = nodeMx[k][2];
		}
	}
}
else
{
	ExecuteCommands ("CovarianceMatrix (covMx, `LF_NAME`);");

}


profileString = "";
profileString * 128;
profileString * "Parameter, Lower Bound, MLE, Upper Bound";
for (k=0; k<varCount; k=k+1)
{
    profileString * ("\n"+usedVars[k]+","+covMx[k][0]+","+covMx[k][1]+","+covMx[k][2]);
}
profileString * 0;

k = baseResPath+".profile";
fprintf (k,CLEAR_FILE,profileString);
profileString = 0;


COVARIANCE_PRECISION = svpc;
stashedValues = {varCount,4};

assignmentString = "";
assignmentString * 256;

variableValues     = {};

fscanf ("lhc_supp.ibf", REWIND, "Raw", sampleInstructions);

availableComputeNodes = MPI_NODE_COUNT - 1;

if (useMPIFlag)
{
	variablesToExport = {{"SAMPLE_N",
	                      "SAMPLE_M",
						  "SAMPLE_LEFT",
						  "SAMPLE_RIGHT",
						  "leftSpan",
						  "leftSTD",
						  "rightSpan",
						  "rightSTD",
						  "currentVarValue"}};
	MPI_NODE_STATUS = {availableComputeNodes, 1}["-1"];
}

for (k=0; k<varCount; k=k+1)
{
	aKey 			    = usedVars[k];
	ExecuteCommands 	("GetInformation(varRange,"+aKey+");");
	stashedValues[k][0] = covMx[k][1];
	/* compute the variance of the approximate normal to the left */
	leftSTD 		= (covMx[k][1]-covMx[k][0])/1.96 * inflationFactor;
	rightSTD		= (covMx[k][2]-covMx[k][1])/1.96 * inflationFactor;
	if (leftSTD)
	{
		SAMPLE_LEFT = Min(1,(covMx[k][1]-varRange[1])/(1.96*leftSTD));
		leftSpan	= Min(covMx[k][1]-varRange[1],1.96*leftSTD);
	}
	else
	{
		SAMPLE_LEFT = 0;
	}
	if (rightSTD)
	{
		SAMPLE_RIGHT = Min(1,(varRange[2]-covMx[k][1])/(1.96*rightSTD));
		rightSpan	 = Min(varRange[2]-covMx[k][1],1.96*rightSTD);
	}
	else
	{
		SAMPLE_RIGHT = 0;
	}
    
	currentVarValue     = covMx[k][1];
	stashedValues[k][1] = Max(covMx[k][1]-(covMx[k][1]-covMx[k][0])*inflationFactor,varRange[1]);
	stashedValues[k][2] = Min(covMx[k][1]+(covMx[k][2]-covMx[k][1])*inflationFactor,varRange[2]);
	stashedValues[k][3] = (stashedValues[k][2]-stashedValues[k][1])/SAMPLE_N;
	assignmentString * (aKey+ "=generatedSamples[itCount]["+k+"];\n");

	if (useMPIFlag && varCount > 1)
	{
		stringToExport = exportVarList (variablesToExport) + 
						 sampleInstructions + 
						 "return thisVarValues;";
						 
						 
		for (nodeIndex = 0; nodeIndex < availableComputeNodes; nodeIndex += 1)
		{
			if (MPI_NODE_STATUS [nodeIndex] < 0)
			{
				MPISend (nodeIndex+1, stringToExport);
				MPI_NODE_STATUS [nodeIndex] = k;
				fprintf (stdout, "Sending variable ", aKey, " to node ", nodeIndex + 1, "\n");
				break;
			}
		}	
		if (nodeIndex == availableComputeNodes)
		{
			ReceiveLHC (1);
		}
	}
	else
	{
		fprintf			(stdout, "Preparing sampling grid for ", aKey, " over [", leftSpan, "-", rightSpan, "]\n");		 
		ExecuteCommands (sampleInstructions);
	}
	

	variableValues[k]   = thisVarValues;
}

if (availableComputeNodes)
{
	leftOver = +MPI_NODE_STATUS["_MATRIX_ELEMENT_VALUE_>=0"];

	for (nodeIndex = 0; nodeIndex < leftOver; nodeIndex += 1)
	{
		ReceiveLHC (0);
	}
}

assignmentString * 0;

//fprintf		(stdout, varCount, "\n", Rows (variableValues), "\n");

varCount	= Rows (covMx);
generatedSamples = {SAMPLE_N,varCount};
indexVector		 = {1,SAMPLE_N}["_MATRIX_ELEMENT_COLUMN_"];

fprintf (stdout, "\nDrawing the LHC sample points...\n");

generatedSamples = generatedSamples["(variableValues[_MATRIX_ELEMENT_COLUMN_])[_MATRIX_ELEMENT_ROW_]"];
generatedSamples = Random(generatedSamples,"LHS");

//fprintf ("/tmp/cubedump.bf", CLEAR_FILE, "generatedSamples = ", generatedSamples, ";\nnames = ", usedVars, ";\n", assignmentString, "\n");

#include "srs.ibf";

//------------------------------------------------------------------------------

function ReceiveLHC (sendOnwards)
{
	MPIReceive									(-1, fromNode, res);
	
	Eval										("thisVarValue="+res);
	variableValues								[MPI_NODE_STATUS[fromNode-1]] = thisVarValue;
	
	fprintf (stdout, "Received variable ", usedVars[MPI_NODE_STATUS[fromNode-1]], " from node ", fromNode, "\n");
	
	if (sendOnwards)
	{
		fprintf (stdout, "Sending variable ", aKey, " to node ", fromNode, "\n");
		MPISend (fromNode, stringToExport);
		MPI_NODE_STATUS [fromNode-1] = k;	
	}
	else
	{
		MPI_NODE_STATUS								[fromNode-1] = -1;
	}
	
	return 0;
}

