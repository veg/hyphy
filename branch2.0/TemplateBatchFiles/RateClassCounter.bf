REPLACE_TREE_STRUCTURE = 1;
/*****************************************************************************/
function echoCatVar (distrInfo)
{
	D = Columns(distrInfo);
	E = 0.0;
	for (k=0; k<D; k=k+1)
	{
		E = distrInfo[0][k]*distrInfo[1][k]+E;
	}
	sampleVar = 0.0;
	for (k=0; k<D; k=k+1)
	{
		sampleVar = sampleVar+(distrInfo[0][k]-E)*(distrInfo[0][k]-E);
	}
	sampleVar = sampleVar/(D-1);
	fprintf  (stdout,"\n\nSample mean = ",E, " (sample variance = ",sampleVar,")\n");
	for (k=0; k<D; k=k+1)
	{
		fprintf (stdout,"\nRate[",Format(k,0,0),"]=",Format(distrInfo[0][k],12,8), " (weight=", 
						  Format(distrInfo[1][k],9,7),")");
	}
	fprintf (stdout,"\n------------------------------------------------\n");
	return 0.0;
}

/*****************************************************************************/

function fitASite (t)
{
	return (treePath < rateMatrix*t) - t;
}


/*****************************************************************************/
function doNCatLFit (useApprox, resp)
{
	gdDefString = "";
	gdDefString * 1024;
	for (mi=1; mi<=resp; mi=mi+1)
	{
		gdDefString*("global PS_"+mi+" = 1/"+((resp+1)-mi)+";\nPS_"+mi+":<1;\n");
	}
	
	gdDefString*("\n\nglobal RS_1 = Random(0,0.5);\nRS_1:<1;RS_1:>0.000000001;\n");

	for (mi=3; mi<=resp; mi=mi+1)
	{
		gdDefString*("global RS_"+mi+" = Random(1,2);"+"\nRS_"+mi+":>1;RS_"+mi+":<100000;\n");
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


	gdDefString*("\n\nglobal c_scale:="+rateStrMx[0]+"*"+freqStrMx[0]);

	for (mi=1; mi<resp; mi=mi+1)
	{
		gdDefString*("+"+rateStrMx[mi]+"*"+freqStrMx[mi]);
	}

	gdDefString*(";\ncategFreqMatrix={{"+freqStrMx[0]);

	for (mi=1; mi<resp; mi=mi+1)
	{
		gdDefString*(","+freqStrMx[mi]);
	}

	gdDefString*("}};\ncategRateMatrix={{"+rateStrMx[0]+"/c_scale");

	for (mi=1; mi<resp; mi=mi+1)
	{
		gdDefString*(","+rateStrMx[mi]+"/c_scale");
	}

	gdDefString*("}};\n\ncategory c  = ("+resp+", categFreqMatrix , MEAN, ,categRateMatrix, 0, 1e25);\n\n");
	gdDefString*0;
	
	ExecuteCommands (gdDefString);

	if (Abs(useApprox)>0)
	{
		mi = Columns (useApprox);
		
		RS_1 = useApprox[0][0]/useApprox[0][1];
		
		if (mi>2)
		{
			RS_3 = useApprox[0][2]/useApprox[0][1];
		}
		
		for (mpiNode = 4; mpiNode < mi; mpiNode = mpiNode+1)
		{
			ExecuteCommands ("RS_"+(mpiNode+1)+"="+useApprox[0][mpiNode]/useApprox[0][mpiNode-1]+";");
		}

		mi2 = 1-useApprox[1][0];

		for (mpiNode = 1; mpiNode < mi; mpiNode = mpiNode+1)
		{
			useApprox[1][mpiNode] = useApprox[1][mpiNode]/mi2;
			mi2 = mi2*(1-useApprox[1][mpiNode]);
			if (mi2 == 0)
			{
				break;
			}
		}

		for (mpiNode = 0; mpiNode < mi; mpiNode = mpiNode+1)
		{
			ExecuteCommands ("PS_"+(mpiNode+1)+"="+useApprox[1][mpiNode]+";");
		}

	}

	Tree		   testTree = treeString;
	if (branchLengths)
	{
		for (mpiNode = 0; mpiNode < totalBranchCount; mpiNode = mpiNode+1)
		{
			eCommand = "testTree."+branchNames[mpiNode]+".t:="+Format(stashedLengths [mpiNode],20,12)+"/totalFactor";
			ExecuteCommands (eCommand);
		}	
	}

	fprintf (stdout, "\n\nFitting ", resp, " rate classes\n");
	/*fprintf (stdout, gdDefString, "\n")
	fprintf (stdout, LIST_ALL_VARIABLES, "\n");
	fscanf (stdin,"String", d);*/
	LikelihoodFunction lfR = (filteredData,testTree);
	Optimize (resR,lfR);
	
	slfo = LIKELIHOOD_FUNCTION_OUTPUT;
	LIKELIHOOD_FUNCTION_OUTPUT = 7;
	
	fName = pathPrefix + "."+ resp;
	fprintf (fName, CLEAR_FILE, lfR);

	LIKELIHOOD_FUNCTION_OUTPUT = slfo;
	
	if (branchLengths)
	{
		AICScore = 2*(resR[1][1]+res[1][1]-res[1][2]-resR[1][0]);
	}
	else
	{
		AICScore = 2*(resR[1][1]-resR[1][0]);	
	}
	
	fprintf (stdout, "\nLikelihood = ", resR[1][0], "\nAIC Score:", AICScore);
	GetInformation (cI, c);
	echoCatVar (cI);
	ClearConstraints (c);
	return 0;

}

/*****************************************************************************/

function  DoNCategoryFitK (weightMatrix)
{
	breakPointMatrix = (SELECTED_CHART_DATA == weightMatrix);
	return breakPointMatrix[3][0];
}

/*****************************************************************************/

function  DoNCategoryFit (weightMatrix)
{
	weightClasses	 = Rows(weightMatrix);
	breakPointMatrix = {3,weightClasses};
	
	runningSum 		= 0;
	targetSum  		= weightMatrix[0];
	currentIndex 	= 0;
	currentSlider	= 0;
	runningSize		= 1;
	valueSum		= 0;
	currentSpan		= SELECTED_CHART_DATA[1][currentIndex];
	
	measCount 		= Columns (SELECTED_CHART_DATA);
	
	
	while (currentIndex < measCount - 1)
	{
		runningSum = runningSum + SELECTED_CHART_DATA[1][currentIndex]/dataPoints;
		
		if ((runningSum >= targetSum)||(measCount-currentIndex <= weightClasses - currentSlider))
		{
			breakPointMatrix[0][currentSlider] = currentIndex;
			breakPointMatrix[1][currentSlider] = runningSize;
			breakPointMatrix[2][currentSlider] = (SELECTED_CHART_DATA[0][currentIndex]*SELECTED_CHART_DATA[1][currentIndex]+valueSum)/(currentSpan+SELECTED_CHART_DATA[1][currentIndex]);	
			runningSize   = 1;	
			valueSum	  = 0;
			currentSlider = currentSlider + 1;
			targetSum 	  = targetSum + weightMatrix[currentSlider];
			currentSpan	  = 0;
		}	
		else
		{
			valueSum	= valueSum	  + SELECTED_CHART_DATA[0][currentIndex]*SELECTED_CHART_DATA[1][currentIndex];
			runningSize = runningSize + 1;
			currentSpan = currentSpan + SELECTED_CHART_DATA[1][currentIndex];
		}
		currentIndex = currentIndex + 1;
	}

	currentSpan	  = currentSpan + SELECTED_CHART_DATA[1][currentIndex];
	valueSum = valueSum	+ SELECTED_CHART_DATA[0][currentIndex]*SELECTED_CHART_DATA[1][currentIndex];
	breakPointMatrix[0][currentSlider] = currentIndex;
	breakPointMatrix[1][currentSlider] = runningSize;	
	breakPointMatrix[2][currentSlider] = valueSum/currentSpan;	
	
		
	currentIndex 	= 0;
	currentSlider	= 0;
	runningOffset	= 0;
	
	logLikelihood	= 0;
	valueSum		= Rows(weightMatrix);
	
	while (currentSlider < valueSum)
	{
		classSize 	= breakPointMatrix[1][currentSlider];
		classWeight = weightMatrix[currentSlider];
			
		if (classWeight > 0.0)
		{
			if (classSize == 1)
			{
				logLikelihood = logLikelihood + SELECTED_CHART_DATA[1][runningOffset] * Log (weightMatrix[currentSlider]);
			}
			else
			{
				classMean 		= breakPointMatrix[2][currentSlider];
				classNorm 		= 0;
				currentIndex 	= runningOffset+classSize;
					
				REWEIGHTED_MATRIX = {measCount,1};
				for (reslider = runningOffset; reslider < currentIndex; reslider = reslider+1)
				{
					targetSum = SELECTED_CHART_DATA[0][reslider];
					targetSum = Exp(-(targetSum-classMean)^2/(2*Abs(classMean)));		
					if (targetSum==0)
					{
						targetSum = 1e-50;
					}
					REWEIGHTED_MATRIX[reslider] = targetSum;
					classNorm = classNorm + targetSum;
				}
				
				REWEIGHTED_MATRIX = Log (REWEIGHTED_MATRIX * (classWeight/classNorm));								
				for (reslider = runningOffset; reslider < currentIndex; reslider = reslider+1)
				{
					logLikelihood = logLikelihood + REWEIGHTED_MATRIX[reslider]*SELECTED_CHART_DATA[1][reslider];
				}
			}
		}
		else
		{
			if (classSize>0)
			{
				return -(1e100);
			}
		}
		runningOffset = runningOffset + classSize;
		currentSlider = currentSlider + 1;
	}
	
	return logLikelihood;
}

/*****************************************************************************/

SetDialogPrompt ("Please specify a nucleotide file:");

DataSet ds = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter filteredData = CreateFilter (ds,1);

fprintf (stdout,"\n______________READ THE FOLLOWING DATA______________\n",ds);

global AC 	= 1;
global AT 	= 1;
global CG 	= 1;
global CT 	= 1;
global GT 	= 1;
global c;

NucleotideMatrix	 = {{*,AC*t,t,AT*t}{AC*t,*,CG*t,CT*t}{t,CG*t,*,GT*t}{AT*t,CT*t,GT*t,*}};
NucleotideMatrixC	 = {{*,AC*t*c,t*c,AT*t*c}{AC*t*c,*,CG*t*c,CT*t*c}{t*c,CG*t*c,*,GT*t*c}{AT*t*c,CT*t*c,GT*t*c,*}};

modelDesc = "";

done = 0;
while (!done)
{
	fprintf (stdout,"\nPlease enter a 6 character model designation (e.g:010010 defines HKY85):");
	fscanf  (stdin,"String", modelDesc);
	if (Abs(modelDesc)==6)
	{	
		done = 1;
	}
}			

rateBiasTerms = {{"AC","1","AT","CG","CT","GT"}};
paramCount	  = 0;

modelConstraintString = "";

for (customLoopCounter2=1; customLoopCounter2<6; customLoopCounter2=customLoopCounter2+1)
{
	for (customLoopCounter=0; customLoopCounter<customLoopCounter2; customLoopCounter=customLoopCounter+1)
	{
		if (modelDesc[customLoopCounter2]==modelDesc[customLoopCounter])
		{
			if (rateBiasTerms[customLoopCounter2] == "1")
			{
				modelConstraintString = modelConstraintString + rateBiasTerms[customLoopCounter]+":="+rateBiasTerms[customLoopCounter2]+";";
			}
			else
			{
				modelConstraintString = modelConstraintString + rateBiasTerms[customLoopCounter2]+":="+rateBiasTerms[customLoopCounter]+";";			
			}
			break;
		}
	}
}	

if (Abs(modelConstraintString))
{
	ExecuteCommands (modelConstraintString);
}
	
ChoiceList (branchLengths,"Estimate Branch Lengths",1,SKIP_NONE,
			"Every Time","Branch lengths are reestimated for every model.",
			"Once","Branch lenghts obtained from the single rate model are reused for subsequent models."
	       );

if (branchLengths<0)
{
	return;
}

SetDialogPrompt ("Save Analysis Results to");
fprintf (PROMPT_FOR_FILE, CLEAR_FILE);
pathPrefix = LAST_FILE_PATH;

HarvestFrequencies (overallFrequencies,  filteredData,1,1,0);
Model NucModel = (NucleotideMatrix, overallFrequencies, 1);

pia = overallFrequencies[0];
pic = overallFrequencies[1];
pig = overallFrequencies[2];
pit = overallFrequencies[3];

	
_DO_TREE_REBALANCE_ = 1;
#include "queryTree.bf";
LikelihoodFunction lf = (filteredData,givenTree);
Optimize (res,lf);
fprintf (stdout, "\nSingle Rate fit\n",lf,"\n");
fprintf (pathPrefix, "\nSingle Rate fit\n",lf,"\n");

DataSet dsA					= ReconstructAncestors (lf);
DataSetFilter filteredDataA = CreateFilter (dsA,1);

seqToBranchMap = {stateCharCount,1};

DataSet		   dsJoint = Combine(dsA,ds);
DataSetFilter filteredDataJ = CreateFilter (dsJoint,1);

branchNames = BranchName (givenTree,-1);
h = Columns (branchNames);

seqToBranchMap 	= {h, 2};

for (k=0; k<filteredData.species; k=k+1)
{
	GetString (seqName, filteredData, k);
	seqToBranchMap[k][0] = -1;
	for (v=0; v<h; v=v+1)
	{
		if (branchNames[v] % seqName)
		{
			seqToBranchMap[k][0] = v;
			seqToBranchMap[v][1] = k;
			break;
		}
	}
}



seqToBranchMap[filteredData.species][0] = h-1;
seqToBranchMap[h-1][1] = filteredData.species;

for (k=1; k<filteredDataA.species; k=k+1)
{
	GetString (seqName, filteredDataA, k);
	seqToBranchMap[filteredData.species+k][0] = -1;
	for (v=0; v<h; v=v+1)
	{
		if (branchNames[v] % seqName)
		{
			seqToBranchMap[k+filteredData.species][0] = v;
			seqToBranchMap[v][1] = k+filteredData.species;
			break;
		}
	}
}

treePath = {3,Columns(branchNames)-1};
for (v=Columns(branchNames)-2; v>=0; v=v-1)
{
	ExecuteCommands ("treePath[2]["+v+"]=givenTree."+branchNames[v]+".t;");
}

/* get codon matrix */

codonInfo  = {filteredData.species, filteredData.unique_sites};
codonInfo2 = {filteredDataA.species, filteredDataA.unique_sites};

GetDataInfo    (dupInfo, filteredData);
GetDataInfo	   (dupInfoA, filteredDataA);
GetDataInfo	   (dupInfoJ, filteredDataJ);

doneSites = {1,filteredDataJ.unique_sites};

for (h=filteredDataJ.unique_sites-1; h>=0; h=h-1)
{
	doneSites  [h] = -1;
}

stateCharCount = 4;

matrixTrick  = {1,stateCharCount};
matrixTrick2 = {1,stateCharCount};

for (h=Columns(matrixTrick)-1; h>=0; h=h-1)
{
	matrixTrick  [h] = h;
	matrixTrick2 [h] = 1;
}

for (v=0; v<filteredData.unique_sites;v=v+1)
{
	for (h=0; h<filteredData.species;h=h+1)
	{
		GetDataInfo (siteInfo, filteredData, h, v);
		_SITE_ES_COUNT = matrixTrick2 * siteInfo; 
		if (_SITE_ES_COUNT[0] == 1)
		{
			siteInfo = matrixTrick * siteInfo;
			codonInfo[h][v] = siteInfo[0];
		}
		else
		{
			codonInfo[h][v] = -1;
		}
	}
}

for (v=0; v<filteredDataA.unique_sites;v=v+1)
{
	for (h=0; h<filteredDataA.species;h=h+1)
	{
		GetDataInfo (siteInfo, filteredDataA, h, v);
		siteInfo = matrixTrick * siteInfo;
		codonInfo2[h][v] = siteInfo[0];
	}
}

_SITE_RESULTS = {4,filteredData.sites};
flatTreeRep	  = Abs (givenTree);

resultMatrix = {filteredData.sites,2};

t = 1;
rateMatrix = NucleotideMatrix*{4,4{0,0,pia}{1,1,pic}{2,2,pig}{3,3,pit}};

rateMatrix[0][0] = -rateMatrix[0][1]-rateMatrix[0][2]-rateMatrix[0][3];
rateMatrix[1][1] = -rateMatrix[1][0]-rateMatrix[1][2]-rateMatrix[1][3];
rateMatrix[2][2] = -rateMatrix[2][0]-rateMatrix[2][1]-rateMatrix[2][3];
rateMatrix[3][3] = -rateMatrix[3][0]-rateMatrix[3][1]-rateMatrix[3][2];

siteRate:>0;

svl = VERBOSITY_LEVEL;
VERBOSITY_LEVEL = -1;

rateSum = 0;

/* estimate the second derivative on the overall rate */

LFCompute (lf,LF_START_COMPUTE);
LFCompute (lf,c1);
h = 0.00001;
for (v=Columns(branchNames)-2; v>=0; v=v-1)
{
	ExecuteCommands ("givenTree."+branchNames[v]+".t = givenTree."+branchNames[v]+".t*(1+h);");
}
LFCompute (lf,c2);
for (v=Columns(branchNames)-2; v>=0; v=v-1)
{
	ExecuteCommands ("givenTree."+branchNames[v]+".t = givenTree."+branchNames[v]+".t*(1-h)/(1+h);");
}
LFCompute (lf,c3);
for (v=Columns(branchNames)-2; v>=0; v=v-1)
{
	ExecuteCommands ("givenTree."+branchNames[v]+".t = givenTree."+branchNames[v]+".t/(1-h);");
}
LFCompute (lf,LF_DONE_COMPUTE);

ce = -(c2+c3-2*c1)/(4*h*h)/filteredData.sites;

fprintf (stdout, "\nApproximate Information per site: ", ce , "\n");
fprintf (pathPrefix, "\nApproximate Information per site: ", ce , "\n");


for (v=0; v<filteredData.sites;v=v+1)
{
	c3 = dupInfoJ[v];
	if (doneSites[c3]>=0)
	{
		resultMatrix[v][0] = doneSites[c3];
		resultMatrix[v][1] = doneSites[c3];
		rateSum = rateSum+doneSites[c3];
		
	}
	else
	{
		k = filteredData.species+1;
		c1 = dupInfoA[v];
		for (h=1; h<filteredDataA.species; h=h+1)
		{
			p1 = seqToBranchMap[k][0];
			p2 = flatTreeRep[p1];
			p2 = seqToBranchMap[p2][1]-filteredData.species;
			
			treePath[0][p1] = codonInfo2[h] [c1];
			treePath[1][p1] = codonInfo2[p2][c1];;
			
			k=k+1;
		}
		
		c2 = dupInfo[v];
		for (h=0; h<filteredData.species; h=h+1)
		{
			p1 = seqToBranchMap[h][0];
			p2 = flatTreeRep[p1];
			p2 = seqToBranchMap[p2][1]-filteredData.species;
			
			cd2 = codonInfo2[p2][c1];
			cd1 =  codonInfo[h] [c2];
			
			treePath[0][p1] = cd2;
			treePath[1][p1] = cd1;
			
			if (cd1<0)
			{
				GetDataInfo    (ambInfo, filteredData, h, c2);	
				weightFactor = 0;
				tempMx = -1;
				ambInfo  	 = ambInfo$overallFrequencies;
				for (k=0; k<stateCharCount; k=k+1)
				{
					if (ambInfo[k]>weightFactor)
					{
						weightFactor = ambInfo[k];
						tempMx = k;
					}
				}
				if (tempMx>=0)
				{
					treePath[1][p1] = tempMx;
				}
				else
				{
					treePath[1][p1] = cd2;			
				}
			}
		}
		
		Optimize (tFit,fitASite(siteRate));
		
		h = tFit[0][0];
		
		resultMatrix[v][0] = h;
		resultMatrix[v][1] = h;
		
		rateSum = rateSum+h;
		doneSites[c3] = h;
	}
}

resultMatrix = resultMatrix*(filteredData.sites/rateSum);
resultMatrix = resultMatrix%1;

dataPoints = Rows(resultMatrix);
temp_data_vector = {2, dataPoints};
currentIndex = 0;

temp_data_vector[0][0] = resultMatrix[0][1];
temp_data_vector[1][0] = 1;

for (nextIndex = 1; nextIndex < dataPoints; nextIndex = nextIndex + 1)
{
	if (resultMatrix[nextIndex][1]!=resultMatrix[nextIndex-1][1])
	{
		currentIndex = currentIndex+1;
		temp_data_vector[0][currentIndex] = resultMatrix[nextIndex][1];
	}
	temp_data_vector[1][currentIndex] = temp_data_vector[1][currentIndex] + 1;
}

SELECTED_CHART_DATA = {2, currentIndex+1};

for (nextIndex = 0; nextIndex <= currentIndex; nextIndex = nextIndex + 1)
{
	SELECTED_CHART_DATA[0][nextIndex] = temp_data_vector[0][nextIndex];
	SELECTED_CHART_DATA[1][nextIndex] = temp_data_vector[1][nextIndex];
}


temp_data_vector = 0;
normal_sigma 	 = .25;
lastMax 		 = -100000000000;

fprintf (stdout, "\nDoing a very approximate fit.\n");
fprintf (pathPrefix, "\nDoing a very approximate fit.\n");

PROFILE_MEAN_VAR_MULT = 0.25/ce;

for (resp = 2; resp <= Columns(SELECTED_CHART_DATA); resp = resp+1)
{
	freqStrMx    = {resp,1};
	freqStrMx[0] = "APS_1";

	for (mi=1; mi<resp-1; mi=mi+1)
	{
		freqStrMx[mi] = "";
		for (mi2=1;mi2<=mi;mi2=mi2+1)
		{
			freqStrMx[mi] = freqStrMx[mi]+"(1-APS_"+mi2+")";		
		}
		freqStrMx[mi] = freqStrMx[mi]+"APS_"+(mi+1);	
	}	

	freqStrMx[mi] = "";
	for (mi2=1;mi2<mi;mi2=mi2+1)
	{
		freqStrMx[mi] = freqStrMx[mi]+"(1-APS_"+mi2+")";		
	}
	freqStrMx[mi] = freqStrMx[mi]+"(1-APS_"+mi+")";	
	mx = {resp,1};
	
	execString = "";
	for (resp2 = 0; resp2 < resp; resp2 = resp2 + 1)
	{
		if (resp2)
		{
			/*execString = execString + "APS_"+resp2+":>0;APS_"+resp2+":<1;APS_"+resp2+"=1/"+(resp-resp2+1)+";";*/
			execString = execString + "global APS_"+resp2+":>0;APS_"+resp2+":<1;";
		}
		execString = execString + "mx[" + resp2 + "]:=" + freqStrMx[resp2] + ";";
	}
	ExecuteCommands (execString);
	
	
	bF = {{0}{-1e100}};
	
	for (k=0; k<10; k=k+1)
	{
		execString = "";
		rM = {resp-1,1};
		for (resp2 = 1; resp2 < resp; resp2 = resp2 + 1)
		{	
			rM [resp2-1] = Random(0,1);
			execString = execString +"APS_"+resp2+"=rM["+(resp2-1)+"];";
		}
		ExecuteCommands (execString);		
		Optimize (bestFit, DoNCategoryFitK(mx));
		if (bestFit[1][0]>bF[1][0])
		{
			bF = bestFit;
			smx = mx;
			sbp = breakPointMatrix;
		}
	}
	
	bestFit = bF;
	mx = smx;
	breakPointMatrix = sbp;
	
	fprintf (stdout, Format (resp,3,0), " rate classes: Log(L) = ", bestFit[1][0], " \n");
	fprintf (pathPrefix, Format (resp,3,0), " rate classes: Log(L) = ", bestFit[1][0], " \n");
	mi = 0;
	for (resp2 = 0; resp2 < resp; resp2 = resp2+1)
	{
		mi = mi + breakPointMatrix[2][resp2]*mx[resp2];
	}

	cI_R = {2,resp};

	for (resp2 = 0; resp2 < resp; resp2 = resp2+1)
	{
		cI_R[0][resp2] = breakPointMatrix[2][resp2]/mi;
		cI_R[1][resp2] = mx[resp2];
		fprintf (stdout, "Class ", Format (resp2+1,3,0), " : ",  Format (cI_R[0][resp2], 8, 5), " weight = ", Format (cI_R[1][resp2], 8, 5), "\n");
		fprintf (pathPrefix, "Class ", Format (resp2+1,3,0), " : ",  Format (cI_R[0][resp2], 8, 5), " weight = ", Format (cI_R[1][resp2], 8, 5), "\n");
	}	
	
	if (bestFit[1][0] - lastMax < 2)
	{
		break;
	}
	saveBPM = breakPointMatrix;
	saveMx	= mx;
	lastMax = bestFit[1][0];
}

resp = resp-1;

fprintf (stdout, "\n\nFitted ", resp, " categories to the data\n\n");
fprintf (pathPrefix, "\n\nFitted ", resp, " categories to the data\n\n");

mi = 0;
for (resp2 = 0; resp2 < resp; resp2 = resp2+1)
{
	mi = mi + saveBPM[2][resp2]*saveMx[resp2];
}

cI_R = {2,resp};

for (resp2 = 0; resp2 < resp; resp2 = resp2+1)
{
	cI_R[0][resp2] = saveBPM[2][resp2]/mi;
	cI_R[1][resp2] = saveMx[resp2];
	fprintf (stdout, "Class ", Format (resp2+1,3,0), " : ",  Format (cI_R[0][resp2], 8, 5), " weight = ", Format (cI_R[1][resp2], 8, 5), "\n");
	fprintf (pathPrefix, "Class ", Format (resp2+1,3,0), " : ",  Format (cI_R[0][resp2], 8, 5), " weight = ", Format (cI_R[1][resp2], 8, 5), "\n");
}

VERBOSITY_LEVEL = svl;

if (branchLengths)
{
	totalBranchCount     = TipCount(givenTree) + BranchCount (givenTree);
	stashedLengths		 = {totalBranchCount,1};
	
	branchNames = BranchName (givenTree,-1);
	
	
	global totalFactor := AC*(2*pia__*pic__)+2*pia__*pig__+(2*pia__*pit__)*AT+
							 (2*pic__*pig__)*CG+(2*pic__*pit__)*CT+(2*pig__*pit__)*GT;
							 
	for (v2 = 0; v2 < totalBranchCount; v2 = v2+1)
	{
		stashedLengths [v2] = BranchLength(givenTree,v2);
	}
}

Model NucModelC   = (NucleotideMatrixC, overallFrequencies, 1);

doNCatLFit (cI_R, resp);

/*labelMatrix =     {1,2};
labelMatrix[0] = "Observed Changes";
labelMatrix[1] = "Scaled Changes";


OpenWindow (CHARTWINDOW,{{"Site Subs"}
						   {"labelMatrix"},
						   {"resultMatrix"},
						   {"Bar Chart"},
						   {"Index"},
						   {labelMatrix[1]},
						   {"Site Index"},
						   {""},
						   {"Scaled Subs"},
						   {"0"}},
						   "SCREEN_WIDTH-60;SCREEN_HEIGHT-50;30;50");*/


go_up   = 1;

lastAICScore = AICScore;

while (resp>2)
{
	resp = resp-1;
	doNCatLFit (0, resp);
	fprintf (stdout, "\nLikelihood = ", resR[1][0], "\nAIC Score:", AICScore);
	fprintf (pathPrefix, "\nLikelihood = ", resR[1][0], "\nAIC Score:", AICScore);
	
	if (lastAICScore < AICScore)
	{
		resp = resp+1;
		break;
	}
	else
	{
		lastAICScore = AICScore;
		go_up = 0;
	}			
}


while (go_up && resp>=2)
{
	resp = resp+1;
	doNCatLFit (0, resp);
		
	fprintf (stdout, "\nLikelihood = ", resR[1][0], "\nAIC Score:", AICScore);
	fprintf (pathPrefix, "\nLikelihood = ", resR[1][0], "\nAIC Score:", AICScore);
	
	if (lastAICScore < AICScore)
	{
		break;
	}
	else
	{
		lastAICScore = AICScore;
	}				
}

if (resp>2)
{
	resp = resp-1;
}

global betaP = 1;
global betaQ = 1;
betaP:>0.05;betaP:<85;
betaQ:>0.05;betaQ:<85;
category pc = (resp-1, EQUAL, MEAN, 
				_x_^(betaP-1)*(1-_x_)^(betaQ-1)/Beta(betaP,betaQ), 
				IBeta(_x_,betaP,betaQ), 
				0, 				   
				1, 			  
				IBeta(_x_,betaP+1,betaQ)*betaP/(betaP+betaQ)
);

global alpha = .5;
alpha:>0.01;alpha:<100;
category c = (resp, pc, MEAN, 
				GammaDist(_x_,alpha,alpha), 
				CGammaDist(_x_,alpha,alpha), 
				0 , 
		  	    1e25,
		  	    CGammaDist(_x_,alpha+1,alpha)
		  	 );

Tree		   testTree = treeString;
if (branchLengths)
{
	for (mpiNode = 0; mpiNode < totalBranchCount; mpiNode = mpiNode+1)
	{
		eCommand = "testTree."+branchNames[mpiNode]+".t:="+Format(stashedLengths [mpiNode],20,12)+"/totalFactor";
		ExecuteCommands (eCommand);
	}	
}


fprintf (stdout, "\n\nFitting gamma mod beta with ", resp, " rate classes\n");
fprintf (pathPrefix, "\n\nFitting gamma mod beta with ", resp, " rate classes\n");
LikelihoodFunction lfR = (filteredData,testTree);


Optimize (resR,lfR);

slfo = LIKELIHOOD_FUNCTION_OUTPUT;
LIKELIHOOD_FUNCTION_OUTPUT = 7;

fName = pathPrefix + ".gammaModBeta."+ resp;
fprintf (fName, CLEAR_FILE, lfR);

LIKELIHOOD_FUNCTION_OUTPUT = slfo;

if (branchLengths)
{
	AICScore = 2*(resR[1][1]+res[1][1]-res[1][2]-resR[1][0]);
}
else
{
	AICScore = 2*(resR[1][1]-resR[1][0]);	
}

fprintf (stdout, "\nLikelihood = ", resR[1][0], "\nAIC Score:", AICScore);
fprintf (pathPrefix, "\nLikelihood = ", resR[1][0], "\nAIC Score:", AICScore);
GetInformation (cI, c);
echoCatVar (cI);
