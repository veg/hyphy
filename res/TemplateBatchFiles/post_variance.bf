likelihoodFnChoice = 0;

if (Rows("LikelihoodFunction")>1) {
	ChoiceList  (likelihoodFnChoice,"Choose a Likelihood Function",1,NO_SKIP,LikelihoodFunction);
}		

if (likelihoodFnChoice<0) {
	return;
} 
			 
ChoiceList  (response,"Covariance Matrix",1,NO_SKIP,
			 "Crude covariances","1st order derivative approximation. This could be slow, since it entails ~p(p-1)/2 calculations of the likelihood function, where p = # independent parameters.",
			 "Finer covariances","2nd order derivative approximation. This could be very slow, since it entails ~p(p-1)*2 calculations of the likelihood function, where p = # of independent parameters.",
			 "Likelihood Profile","95% CI based on the quadratic approximation to the likelihood surface. Should be faster than full covariance matrices for large numbers of parameters and may be more robust to statistical errors.");

			 
if (response<0) {
	return 1;
}

if (response == 2) {
	COVARIANCE_PRECISION = 0.95;
} else {	
	COVARIANCE_PRECISION = response+1;
}

GetString (likelihoodFunctionName,LikelihoodFunction,likelihoodFnChoice);

if (RESTORE_GLOBALS) {
	RestoreGlobalValues (likelihoodFnChoice);
}

CovarianceMatrix (covMatrix, *likelihoodFunctionName);

fprintf   (stdout, "\n\n\t\tVARIANCE ESTIMATES\n\n");

nameWidth = 10;

for (covCounter = 0; covCounter<Rows(covMatrix); covCounter += 1) {
	GetString (argName,*likelihoodFunctionName,covCounter);
	nameWidth = Max (nameWidth, Abs (argName));
}



fprintf (stdout,"+");
for (counter = 0; counter<nameWidth+1; counter=counter+1)
{
	fprintf (stdout,"-");
}
fprintf (stdout,"+");
for (counter = counter+1; counter<nameWidth+11; counter=counter+1)
{
	fprintf (stdout,"-");
}
fprintf (stdout,"+");
counter = counter+1;
for (; counter<nameWidth+24; counter=counter+1)
{
	fprintf (stdout,"-");
}
fprintf (stdout,"+");
counter = counter+1;
for (; counter<nameWidth+45; counter=counter+1)
{
	fprintf (stdout,"-");
}
fprintf (stdout,"+\n| Parameter ");
for (covCounter = 0; covCounter<nameWidth-10; covCounter=covCounter+1)
{
	fprintf (stdout," ");
}
fprintf (stdout,"|  Value  |  Std. Dev. | 95% Conf. Interval |");
for (covCounter = 0; covCounter<Rows(covMatrix); covCounter=covCounter+1)
{
	fprintf (stdout,"\n+");
	for (counter = 0; counter<nameWidth+1; counter=counter+1)
	{
		fprintf (stdout,"-");
	}
	fprintf (stdout,"+");
	for (counter = counter+1; counter<nameWidth+11; counter=counter+1)
	{
		fprintf (stdout,"-");
	}
	fprintf (stdout,"+");
	counter = counter+1;
	for (; counter<nameWidth+24; counter=counter+1)
	{
		fprintf (stdout,"-");
	}
	fprintf (stdout,"+");
	counter = counter+1;
	for (; counter<nameWidth+45; counter=counter+1)
	{
		fprintf (stdout,"-");
	}
	fprintf (stdout,"+\n| ");
	GetString (argName,*likelihoodFunctionName,covCounter);
	fprintf (stdout, argName);
	argLength = Abs(argName);
	for (counter = 0; counter<nameWidth-argLength; counter=counter+1)
	{
		fprintf (stdout," ");
	}
	GetInformation (parValue,argName);
	if (response<2)
	{
		parVar = covMatrix[covCounter][covCounter];
		parValue=parValue[0];
		parVar = Sqrt(Abs(parVar));
		fprintf (stdout,"|",Format(parValue,9,4),"|",Format(parVar,12,8),"|", "(", Format(Max(parValue-1.96*parVar,0),8,4),",",Format(parValue+1.96*parVar,8,4),") |");
	}
	else
	{
		fprintf (stdout,"|",Format(covMatrix[covCounter][1],9,4),"|    N/A     |", "(", Format(Max(covMatrix[covCounter][0],0),8,4),",",Format(covMatrix[covCounter][2],8,4),") |");	
	}
}
fprintf (stdout,"\n+");
for (counter = 0; counter<nameWidth+1; counter=counter+1)
{
	fprintf (stdout,"-");
}
fprintf (stdout,"+");
for (counter = counter+1; counter<nameWidth+11; counter=counter+1)
{
	fprintf (stdout,"-");
}
fprintf (stdout,"+");
counter = counter+1;
for (; counter<nameWidth+24; counter=counter+1)
{
	fprintf (stdout,"-");
}
fprintf (stdout,"+");
counter = counter+1;
for (; counter<nameWidth+45; counter=counter+1)
{
	fprintf (stdout,"-");
}
fprintf (stdout,"+\n| ");
