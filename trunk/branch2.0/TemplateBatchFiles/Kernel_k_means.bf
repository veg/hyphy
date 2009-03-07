ExecuteAFile 	("Kernel_support.ibf");
read_kernel_matrix (0);

fprintf (stdout, "How many clusters?");
fscanf (stdin,"Number",rateClassesCount);
	
DEFAULT_FILE_SAVE_NAME = "K-means.csv";

SetDialogPrompt ("Write CSV clusters to:");
fprintf 		 (PROMPT_FOR_FILE,CLEAR_FILE,KEEP_OPEN,"Gene,Cluster\n");


autoStepFlag	 = 0;

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

produceOffspring		= 16;
populationSize  		= 2*produceOffspring;
incestDistance  		= 0;
generationCount		  	= 5000;
maxSampleTries			= populationSize*10;
mutationThreshhold		= 0.0001;
mutationProb			= 0.10;
mutationProbDecrease	= 1.0;
annealingPhase			= 100;
totalSampleCounter		= 0;
localMutationRate		= 0.05;
localMutationInterval	= 20;

stoppingCriterion		= 20;
sampleCount				= 0;
familyControlSize		= produceOffspring$6;

predef   				= {};

MasterList				= {};
verboseFlag				= 0;
stateVectorDimension    = points;
branchUpdateThreshold   = generationCount;

k_trace					= mx_diag (kernel_matrix);

currentPopulation	  = {};
currentPopulation [0] = {stateVectorDimension,1};
crapAIC			   	  = dualMeans(kernel_matrix,1,currentPopulation[0]);
sortedScores	      = {populationSize,2};
sortedScores[0][0]    = crapAIC;
sortedScores[0][1]    = 0;

ExecuteAFile ("GA_CHC_kernel.ibf");

dualMeans 	 (kernel_matrix, rateClassesCount, currentPopulation[populationSize-1]);
allocation = distances[-1][0];

byCluster  = {};

for (k=0; k<points; k=k+1)
{
	fprintf (LAST_FILE_PATH,names[k+1],",",allocation[k],"\n");
	byCluster [allocation[k]] = byCluster [allocation[k]]+1;
}

columnHeaders = {{"Cluster ID"}};

OpenWindow (CHARTWINDOW,{{"Kernel k-means cluster allocation"}
						{"columnHeaders"}
						{"allocation"}
						{"Scatterplot"}
						{"Index"}
						{"Cluster ID"}
						{"Data set ID"}
						{""}
						{"Cluster ID"}
						{"0"}
						{""}
						{"-1;-1"}
						{"10;1.309;0.785398"}
						{"Times:12:0;Times:10:0;Times:12:2"}
						{"0;0;16777215;16777215;0;0;6579300;11842740;13158600;14474460;0;3947580;16777215;0;6845928;16771158;2984993;9199669;7018159;1460610;16748822;11184810;14173291"}
						{"16,0,0"}
						},
						"572;516;37;51");
		
fprintf (stdout, "Cluster size report:\n");

for (k=0; k<rateClassesCount; k=k+1)
{
	fprintf (stdout, "\tCluster ", Format(k+1,3,0), " has ", Format (byCluster[k],5,0), " points\n"); 
}

fprintf (LAST_FILE_PATH,CLOSE_FILE);

/*------------------------------------------------------------------------------*/
	
function dualMeans (K,N,initVector)
{
	K_dim = Rows(K);
	A	  = {K_dim,N};
	/* random initial assignment */
	distances = {K_dim,2};
	for (_k = 0; _k < K_dim; _k=_k+1)
	{
		distances[_k][0] = initVector[_k];
		A[_k][distances[_k][0]] = 1;
	}
	
	goOn = 1;

	while (goOn)
	{
	  goOn = 0;
	  /*E = A * diag(1./sum(A));*/
	  E = A * diag ((sum(A))["1./_MATRIX_ELEMENT_VALUE_"]);
	  Z = {K_dim,1}["1"] * Transpose(mx_diag (Transpose(E)*K*E)) - K*E*2.0;
	  
	  for (k = 0; k < K_dim; k=k+1)
	  {
	  		min     = 1e100;
	  		min_idx = 0;
	  		
	  		for (k2 = 0; k2 < N; k2=k2+1)
	  		{
	  			if (Z[k][k2] < min)
	  			{
	  				min = Z[k][k2];
	  				min_idx = k2;
	  			}
	  		}
	  		
	  		
	  		distances[k][1] = min;
	  		
	  		if (min_idx != distances[k][0])
	  		{
	  			A[k][min_idx] 		  = 1;
	  			A[k][distances[k][0]] = 0;
	  			distances[k][0] = min_idx;
	  			goOn = 1;
	  		}
	  }
	}
	
	return  -Abs(k_trace+distances[-1][1]);
}


/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function MatrixToString (rateMatrix)
{
	outString = "";
	outString * 256;
	for (h=0; h<stateVectorDimension; h=h+1)
	{
		outString * (","+rateMatrix[h]);				
	}
	outString * 0;
	return outString;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function SpawnRandomString (clsCnt)
{
	rModel = {stateVectorDimension,1};
	for (h=0; h<stateVectorDimension; h=h+1)
	{
		rModel[h] = Random(0,clsCnt)$1;
	}
	
	return MakeStringCanonical(rModel,clsCnt);
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function MakeStringCanonical (randomModel, classCount)
{
	compressedString = randomModel%0;
	return randomModel;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function UpdateBL (dummy)
{
	return 0;
}
/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function RunASample (modelDF, jobIndex)
{
	sampleString = MatrixToString (cString);
	myAIC 		 = MasterList[sampleString];
	
	if (myAIC > 0.1)
	{
		if (resultProcessingContext==0)
		{
			sortedScores[jobIndex][0] = myAIC;
		}
		else
		{
			intermediateProbs[jobIndex][0] = myAIC;	
		}			
		return 0;
	}
	
	myAIC   = dualMeans (kernel_matrix,rateClassesCount,cString);
	ReceiveJobs (1, jobIndex);
	return 0;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function IsChildViable (v)
{
	return v;
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function ReceiveJobs (sendOrNot, ji)
{
	if (resultProcessingContext==0)
	{
		sortedScores[ji][0] = myAIC;
	}
	else
	{
		intermediateProbs[ji][0] = myAIC;	
	}
	if (ji>=0)
	{
		if (resultProcessingContext==0)
		{
			MasterList [MatrixToString (currentPopulation[ji])] = myAIC;
		}
		else
		{
			MasterList [MatrixToString (children[ji-populationSize])] = myAIC;
		}
	}
	return fromNode-1;
}

