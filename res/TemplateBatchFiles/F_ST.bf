RequireVersion ("0.9920060901");

/*------------------------------------------------------------------------------*/
function doSNN (vpart1, vpart2, vsize1, vsize2)
{
	for (k=0; k<vsize1;k=k+1)
	{
		idx1 = vpart1[k];
		minD = 1e100;
		X_j	 = 0;
		W_j  = 0;
		
		for (k2=0; k2<vsize1; k2=k2+1)
		{
			if (k!=k2)
			{
				idx2 = vpart1[k2];
				cd   = distanceMatrix[idx1][idx2];
				if (cd < minD)
				{
					W_j  = 1;
					X_j  = 1;
					minD = cd;
				}
				else
				{
					if (cd == minD)
					{
						W_j = W_j+1;
						X_j = X_j+1;
					}
				}
			}
		}

		for (k2=0; k2<vsize2; k2=k2+1)
		{
			if (k!=k2)
			{
				idx2 = vpart2[k2];
				cd   = distanceMatrix[idx1][idx2];
				if (cd < minD)
				{
					W_j  = 1;
					X_j  = 0;
					minD = cd;
				}
				else
				{
					if (cd == minD)
					{
						W_j = W_j+1;
					}
				}
			}
		}
		s_nn = s_nn + X_j/W_j;
	}
	return 0;
}
/*------------------------------------------------------------------------------*/

function computeCompartmentValues (part1, part2)
{
	resMatrix = {4,1};
	template1 = {ds.species, ds.species};
	template2 = {ds.species, ds.species};
	for (k=0; k<clASize;k=k+1)
	{
		idx1 = part1[k];
		for (k2=k+1; k2<clASize;k2=k2+1)
		{
			idx2 = part1[k2];
			template1 [idx1][idx2] = 1;
			template1 [idx2][idx1] = 1;
		}
	}

	for (k=0; k<clBSize;k=k+1)
	{
		idx1 = part2[k];
		for (k2=k+1; k2<clBSize;k2=k2+1)
		{
			idx2 = part2[k2];
			template2 [idx1][idx2] = 1;
			template2 [idx2][idx1] = 1;
		}
	}
	
	template3 = totalUnitSqr-template1-template2;
	
	count1 = totalUnitRow*(template1*totalUnitCol);
	count2 = totalUnitRow*(template2*totalUnitCol);
	count3 = totalUnitRow*(template3*totalUnitCol);
	
	sum1   = totalUnitRow*(template1$distanceMatrix*totalUnitCol);
	sum2   = totalUnitRow*(template2$distanceMatrix*totalUnitCol);
	sum3   = totalUnitRow*(template3$distanceMatrix*totalUnitCol);
	
	dd1 = sum1[0]/count1[0]*f_1^2+sum2[0]/count2[0]*f_2^2;
	
	resMatrix[2] = sum3[0]/(count3[0]-bothCladesSize);	/* pi_B */
	resMatrix[1] = dd1/(f_1^2+f_2^2);			        /* pi_S */
	resMatrix[0] = dd1 + 2*f_1*f_2*resMatrix[2];        /* pi_T */
	
	/* now compute S_nn */
	
	s_nn = 0;
	
	doSNN (part1, part2, clASize, clBSize);
	doSNN (part2, part1, clBSize, clASize);
	
	resMatrix[3] = s_nn / (clASize+clBSize);
	
	return    resMatrix;
}

/*------------------------------------------------------------------------------*/

function reportBSTRP (_ds, obs)
{
	fprintf (stdout, "Observed value  : ", Format(obs,8,3), "\n",
					 "Bootst. mean    : ", Format(_ds["Mean"],8,3), "\n",
					 "Bootst. median  : ", Format(_ds["Median"],8,3), "\n",
					 "Bootst. st. dev.: ", Format(_ds["Std.Dev"],8,3), "\n",
					 "Bootst. 95% CI  : ", Format(_ds["2.5%"],8,3), " - ", Format(_ds["97.5%"],8,3), "\n");
	return 0;
}

/*------------------------------------------------------------------------------*/

ChoiceList (distanceChoice, "Distance Computation",1,SKIP_NONE,
			"Distance formulae","Use one of the predefined distance measures based on data comparisons. Fast.",
			"Full likelihood","Estimate distances using pairwise MLE. More choices but slow.",
			"Load Matrix","Re-load a previously computed matrix in HyPhy format");
			
if (distanceChoice < 0) {
	return 0;
}


if (distanceChoice == 2) {
    ChoiceList (useSeqData,"Sequence names",1,SKIP_NONE,"Included","Sequence names are in the matrix file preceding the distance matrix.",
                           "Read from an alignment","Gather sequence names from a separate alignment file");

    if (useSeqData < 0) {
        return 0;
    }
} else {
    useSeqData = 1;
}


if (useSeqData) 
{
    ChoiceList (dataType,"Data type",1,SKIP_NONE,"Nucleotide/Protein","Nucleotide or amino-acid (protein).",
                         "Codon","Codon (several available genetic codes).");
    
    if (dataType<0) {
        return;
    }

    if (dataType) {
        SetDialogPrompt ("Please choose a codon data file:");
        ExecuteAFile ("TemplateModels/chooseGeneticCode.def");
    }
    else {
        SetDialogPrompt ("Please choose a nucleotide or amino-acid data file:");
    }
    
    DataSet ds = ReadDataFile (PROMPT_FOR_FILE);

    fprintf (stdout, "\nRead the following data:", ds,"\n\n");
} else {
	SetDialogPrompt ("Load the names/distance matrix");
	fscanf          (PROMPT_FOR_FILE, "Matrix,NMatrix", namesMatrix,distanceMatrix);
	dim_names = Rows(namesMatrix)*Columns(namesMatrix);
	assert (dim_names == Columns (distanceMatrix), "Dimension mismatch between the names vector and the distance matrix");
	fakeDS = ""; fakeDS * 128;
	for (_i = 0; _i < dim_names; _i+=1) {
	    fakeDS * (">" + namesMatrix [_i] + "\nA\n");
	}
	fakeDS * 0;
	DataSet ds = ReadFromString (fakeDS);
}

promptFor2ndRegExp = 1;
ExecuteAFile ("partitionSequences.ibf");
promptFor2ndRegExp = 0;

bothCladesSize = clASize + clBSize; 
f_1 = clASize/bothCladesSize;
f_2 = clBSize/bothCladesSize;

fprintf (stdout, "\nProportion of sequence in population 1: ", f_1,
				 "\nProportion of sequence in population 2: ", f_2, "\n");

resultAVL = {"Proportion 1": f_1,
			 "Proportion 2": f_2};

p1vector = {};
p2vector = {};
overallSample2 = {1,bothCladesSize};
overallSample  = {};
sequencesIn	   = {};

for (specIndex = 0; specIndex < ds.species; specIndex = specIndex + 1)
{
	GetString (specName, ds, specIndex);
	
	if (cladeA[specName]>0)
	{
		p1vector			[Abs(p1vector)] 		= specIndex;
		overallSample2		[Abs(overallSample)] 	= 1;
		overallSample		[Abs(overallSample)] 	= specIndex;
		sequencesIn			[specIndex] 			= 1;
	}
	else
	{
		if (cladeB[specName]>0)
		{
			p2vector		[Abs(p2vector)] 		= specIndex;
			overallSample2	[Abs(overallSample)] 	= 2;
			overallSample	[Abs(overallSample)] 	= specIndex;
			sequencesIn		[specIndex] 			= 1;
		}
	}
}

if (dataType)
{
	DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);
}
else
{
	DataSetFilter filteredData = CreateFilter (ds,1,"","");
}

if (useSeqData) {
    distanceMatrix = {ds.species, ds.species};
}

if (distanceChoice == 1)
{
	SelectTemplateModel(filteredData);
	
	fprintf (stdout,"\nHYPHY is computing pairwise maximum likelihood distance estimates. A total of ", Format(ds.species*(ds.species-1)/2,0,0),
				    " estimations will be performed.\n");

	MESSAGE_LOGGING = 0;

	for (i = 0; i<ds.species-1; i=i+1)
	{
		for (j = 0; j<=i; j = j+1)
		{
			if (dataType)
			{
				DataSetFilter twoSpecFilter = CreateFilter (filteredData,3,"",(speciesIndex==i+1)||(speciesIndex==j),GeneticCodeExclusions);
			}
			else
			{
				DataSetFilter twoSpecFilter = CreateFilter (filteredData,1,"",(speciesIndex==i+1)||(speciesIndex==j));
			}
			if (FREQUENCY_SENSITIVE)
			{
				if (USE_POSITION_SPECIFIC_FREQS)
				{
					HarvestFrequencies (vectorOfFrequencies,filteredData,3,1,1);
				}
				else
				{
					HarvestFrequencies (vectorOfFrequencies,twoSpecFilter,1,1,0);
				}
			}
			if (FREQUENCY_SENSITIVE)
			{
				MULTIPLY_BY_FREQS = PopulateModelMatrix ("modelMatrix",vectorOfFrequencies);
				if (dataType)
				{
					codonFrequencies = BuildCodonFrequencies (vectorOfFrequencies);
					Model pairModel = (modelMatrix, codonFrequencies, MULTIPLY_BY_FREQS);
				}
				else
				{
					Model pairModel = (modelMatrix, vectorOfFrequencies, MULTIPLY_BY_FREQS);
				}
			}		
			else
			{
				if (i+j==0)
				{
					MULTIPLY_BY_FREQS = PopulateModelMatrix ("modelMatrix",equalFreqs);
					Model pairModel = (modelMatrix, equalFreqs, MULTIPLY_BY_FREQS);
				}
			}
			
			Tree Inferred_Tree = (1,2);
			LikelihoodFunction lf = (twoSpecFilter,Inferred_Tree);
			Optimize (res,lf);
			k = BranchLength (Inferred_Tree,0);
			distanceMatrix[j][i+1] = k;
			distanceMatrix[i+1][j] = k;
		}
	}
}
else
{
	if (distanceChoice == 2)
	{
        if (useSeqData) {
            SetDialogPrompt ("Load the distance matrix");
            fscanf (PROMPT_FOR_FILE,"NMatrix",distanceMatrix);
        }
		if ((Rows(distanceMatrix) != ds.species)||(Columns(distanceMatrix) != ds.species))
		{
			fprintf (stdout, "\nThe dimensions of the distance matrix are incompatible with the data set.\n");
			return  0;
		}		
	}
	else
	{
		DISTANCE_PROMPTS = 1;
		#include "chooseDistanceFormula.def";
		dummy = InitializeDistances (0);
					    
		tdc = 0;
		tdp = 0;
		ldp = 0;
		togo = ds.species*(ds.species-1)/2;

		fprintf (stdout,"\nHYPHY is computing pairwise distance estimates. A total of ", Format(togo,0,0),
					    " estimations will be performed.\n");
		
		for (i = 0; i<ds.species-1; i=i+1)
		{
			for (j = 0; j<=i; j = j+1)
			{
				k = ComputeDistanceFormula (i+1,j);
				distanceMatrix[j][i+1] = k;
				distanceMatrix[i+1][j] = k;
			}
			tdc = tdc+i+1;
			tdp = (tdc/togo * 100)$1;
			if (tdp>ldp)
			{
				ldp = tdp;
				fprintf (stdout, ldp, "% done\n");
			}
		}

		DISTANCE_PROMPTS = 0;
	}
}

/* make some unit matrices */

totalUnitRow = {1,ds.species};
totalUnitCol = {ds.species,1};
totalUnitSqr = {ds.species,ds.species};

for (k=0; k<ds.species; k=k+1)
{
	totalUnitRow[k] = 1;
	totalUnitCol[k] = 1;
	if (sequencesIn[k])
	{
		for (k2=0; k2<ds.species; k2=k2+1)
		{
			if (sequencesIn[k2])
			{
				totalUnitSqr [k][k2] = 1;
			}
		}
	}
}

resMx = computeCompartmentValues (p1vector,p2vector);
pi_D = resMx[2]-resMx[1];

F_ST_OBS = {{pi_D/(resMx[1]+pi_D),pi_D/(resMx[1]+resMx[2]),1-resMx[1]/resMx[0], resMx[3]}};

fprintf (stdout, "\n\nPopulation characterisitcs:",
				 "\nMetapopulation diversity (pi_T)       = ", resMx[0],
				 "\nMean subpopulation diversity (pi_S)   = ", resMx[1],
				 "\nMean interpopulation diversity (pi_B) = ", resMx[2],
				 "\n\nF_ST\n",
				 "\nHudson, Slatkin and Maddison (Genetics     132:583-589): ",F_ST_OBS[0],
				 "\nSlatkin                      (Evolution    47:264-279) : ",F_ST_OBS[1],
				 "\nHudson, Boos and Kaplan	     (Mol Bio Evol 9: 138-151) : ",F_ST_OBS[2],
				 "\nHudson (S_nn)                (Genetics     155:2011-14): ",F_ST_OBS[3], "\n");
				 
F_ST_OBS = F_ST_OBS;


resultAVL ["pi_T"] = resMx[0];
resultAVL ["pi_S"] = resMx[1];
resultAVL ["pi_B"] = resMx[2];
resultAVL ["Hudson, Slatkin and Maddison"] = F_ST_OBS[0];
resultAVL ["Slatkin"] = F_ST_OBS[1];
resultAVL ["Hudson, Boos and Kaplan"] = F_ST_OBS[2];
resultAVL ["Hudson (S_nn)"] = F_ST_OBS[3];

				 
ChoiceList (resample,"Bootstrap Estimators",1,SKIP_NONE,
					 "Skip","Do not perform a permutation test.",
				     "Sure","Resample with replacement within populations to estimate sampling properties of the estimators.");
				 

if (resample < 0)
{
	return 0;
}

if (resample)
{
	sampleCount = 0;
	while (sampleCount < 1)
	{
		fprintf (stdout, "\nHow many permutations would you like to perform?");
		fscanf (stdin,"Number",sampleCount);
	}

	F_ST_1 = {sampleCount,1}; 
	F_ST_2 = {sampleCount,1}; 
	F_ST_3 = {sampleCount,1}; 
	F_ST_4 = {sampleCount,1}; 
	step  = sampleCount$20; 

	fprintf (stdout, "Running the simulations...\n"); 

	saveDM = distanceMatrix; 

	p1Size = Abs (p1vector); 
	p2Size = Abs (p2vector); 
	basePartition1 = {1,p1Size}; 
	basePartition2 = {1,p2Size}; 

	for (k=0; k<p1Size;k=k+1) 
	{ 
		basePartition1[k] = p1vector[k]; 
	} 

	for (k=0; k<p2Size;k=k+1) 
	{ 
		basePartition2[k] = p2vector[k]; 
	} 
	 
	maxSampleCount =  10*sampleCount;
	 
	for (sampleID = 0; sampleID < sampleCount; sampleID = sampleID + 1) 
	{ 
		 /* repeat the following code for each replicate */ 
		 resampleP1 = Random (basePartition1,1); /* resample vector of indices with replacement */ 
		 resampleP2 = Random (basePartition2,1);
		 distanceMatrix = {Rows(saveDM),Rows(saveDM)}; 
		 
		 for (k=0; k<p1Size;k=k+1) /* repopulate distances for population 1 */ 
		 { 
			ki = resampleP1[k];
			ii = p1vector[k]; 
			for (k2 = k+1; k2 < p1Size; k2=k2+1) 
			{ 
				k2i = resampleP1[k2]; 
				ii2 = p1vector[k2];
				distanceMatrix [ii][ii2] = saveDM[ki][k2i]; 
				distanceMatrix [ii2][ii] = saveDM[k2i][ki]; 
			} 
		 } 

		 for (k=0; k<p2Size;k=k+1) /* repopulate distances for population 2 */ 
		 { 
		  	ki = resampleP2[k]; 
			ii = p2vector  [k]; 
		  	for (k2 = k+1; k2 < p2Size; k2=k2+1) 
			{ 
				k2i = resampleP2[k2]; 
				ii2 = p2vector  [k2]; 
				distanceMatrix  [ii][ii2] = saveDM[ki][k2i]; 
				distanceMatrix  [ii2][ii] = saveDM[k2i][ki]; 
			} 
		 } 

		for (k=0; k<p1Size;k=k+1) /* repopulate interpopulation distances */ 
		{ 
			ki = resampleP1[k]; 
			ii = p1vector[k]; 
			for (k2 = 0; k2 < p2Size; k2=k2+1) 
			{ 
				k2i = resampleP2[k2]; 
				ii2 = p2vector  [k2]; 
				distanceMatrix [ii][ii2] = saveDM[ki][k2i]; 
				distanceMatrix [ii2][ii] = saveDM[k2i][ki]; 
			} 
		} 

		if (Abs(distanceMatrix) == 0)
		{
			if (maxSampleCount)
			{
				maxSampleCount = maxSampleCount - 1;
				sampleID = sampleID-1;
				continue;
			}
			else
			{
				fprintf (stdout, "[ERROR: TOO MANY IDENTICAL SEQUENCES; CAN'T RESAMPLE WITHOUT OBTAINING ZERO DISTANCE MATRICES IN THE ALLOCATED NUMBER OF TRIES]\n");
			}
		}
	
		resMx = computeCompartmentValues (basePartition1,basePartition2); 
		pi_D = resMx[2]-resMx[1]; 
	
		F_ST_1[sampleID] = pi_D/(resMx[1]+pi_D); 
		F_ST_2[sampleID] = pi_D/(resMx[1]+resMx[2]); 
		F_ST_3[sampleID] = 1-resMx[1]/resMx[0]; 
		F_ST_4[sampleID] = resMx[3]; 

		if ((sampleID+1)%step == 0) 
		{ 
			fprintf (stdout, Format ((sampleID+1)*100/sampleCount, 6, 2), "% done\n"); 
		} 
		
		maxSampleCount = maxSampleCount - 1;
	} 

	ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "Utility" + DIRECTORY_SEPARATOR + 
										 "DescriptiveStatistics.bf");
										 
	fprintf (stdout, "\n\nBootstrapped estimator statistics.\n",
					 "\nHudson, Slatkin and Madison (Genetics     132:583-589)\n");
	_stats = GatherDescriptiveStats (F_ST_1);				 
	reportBSTRP (_stats, F_ST_OBS[0]);
	resultAVL ["Bootstrap (Hudson, Slatkin and Maddison)"] = _stats;

	fprintf (stdout, "\nSlatkin                     (Evolution    47:264-279)\n");
	_stats = GatherDescriptiveStats (F_ST_2);				 
	reportBSTRP (_stats, F_ST_OBS[1]);
	resultAVL ["Bootstrap (Slatkin)"] = _stats;

	fprintf (stdout, "\nHudson, Boos and Kaplan	    (Mol Bio Evol 9: 138-151)\n");
	_stats = GatherDescriptiveStats (F_ST_3);				 
	reportBSTRP (_stats, F_ST_OBS[2]);
	resultAVL ["Bootstrap (Hudson, Boos and Kaplan)"] = _stats;

	fprintf (stdout, "\nHudson (S_nn)               (Genetics     155:2011-14)\n");
	_stats = GatherDescriptiveStats (F_ST_4);				 
	reportBSTRP (_stats, F_ST_OBS[3]);
	resultAVL ["Bootstrap (Hudson (S_nn))"] = _stats;
	
	distanceMatrix = saveDM;
}

ChoiceList (resample,"Permutation Test",1,SKIP_NONE,
					 "Skip","Do not perform a permutation test.",
				     "But of course","Randomly allocate sequences into subpopulations and tabulate the distribution of various F_ST statistics.");
				 
if (resample > 0)
{

	sampleCount = 0;
	while (sampleCount < 1)
	{
		fprintf (stdout, "\nHow many permutations would you like to perform?");
		fscanf (stdin,"Number",sampleCount);
	}
	
	F_ST_1 = {sampleCount,2};
	F_ST_2 = {sampleCount,2};
	F_ST_3 = {sampleCount,2};
	F_ST_4 = {sampleCount,2};
	step   = sampleCount$20;
	
	fprintf (stdout, "Running the simulations...\n");
	
	for (sampleID = 0; sampleID < sampleCount; sampleID = sampleID + 1)
	{
		aSample = Random(overallSample2,0);
		p1_1 	= {};
		p2_1 	= {};
		for (k=0; k<bothCladesSize;k=k+1)
		{
			k2 = overallSample[k];
			if (aSample[k] == 2)
			{
				p2_1[Abs(p2_1)] = k2;
			}
			else
			{
				p1_1[Abs(p1_1)] = k2;
			}
		}
		resMx = computeCompartmentValues (p1_1,p2_1);
		pi_D = resMx[2]-resMx[1];	
		F_ST_1[sampleID][1] = pi_D/(resMx[1]+pi_D);
		F_ST_2[sampleID][1] = pi_D/(resMx[1]+resMx[2]);
		F_ST_3[sampleID][1] = 1-resMx[1]/resMx[0];
		F_ST_4[sampleID][1] = resMx[3];
		
		if ((sampleID+1)%step == 0)
		{
			fprintf (stdout, Format ((sampleID+1)*100/sampleCount, 6, 2), "% done\n");
		}
	}
	
	F_ST_1 = F_ST_1%1;
	F_ST_2 = F_ST_2%1;
	F_ST_3 = F_ST_3%1;
	F_ST_4 = F_ST_4%1;
	
	pv = {{sampleCount,sampleCount,sampleCount,sampleCount}};
	
	for (sampleID = 0; sampleID < sampleCount; sampleID = sampleID + 1)
	{
		step = sampleID/sampleCount;
		F_ST_1[sampleID][0] = step;
		F_ST_2[sampleID][0] = step;
		F_ST_3[sampleID][0] = step;
		F_ST_4[sampleID][0] = step;
		if (pv[0]==sampleCount)
		{
			if (F_ST_1[sampleID][1] > F_ST_OBS[0])
			{
				pv[0] = sampleID;
			}
		}
		if (pv[1]==sampleCount)
		{
			if (F_ST_2[sampleID][1] > F_ST_OBS[1])
			{
				pv[1] = sampleID;
			}
		}
		if (pv[2]==sampleCount)
		{
			if (F_ST_3[sampleID][1] > F_ST_OBS[2])
			{
				pv[2] = sampleID;
			}
		}
		if (pv[3]==sampleCount)
		{
			if (F_ST_4[sampleID][1] > F_ST_OBS[3])
			{
				pv[3] = sampleID;
			}
		}
	}
	
	
	fprintf (stdout, "\n\nProb {Random F_ST > Observed F_ST}\n",
					 "\nHudson, Slatkin and Maddison : ",(sampleCount-pv[0])/sampleCount,
					 "\nSlatkin                      : ",(sampleCount-pv[1])/sampleCount,
					 "\nHudson, Boos and Kaplan	     : ",(sampleCount-pv[2])/sampleCount,
					 "\nHudson, S_nn                 : ",(sampleCount-pv[3])/sampleCount, "\n");
	
	resultAVL ["P (Hudson, Slatkin and Maddison)"] = (sampleCount-pv[0])/sampleCount;
	resultAVL ["P (nSlatkin)"] = (sampleCount-pv[1])/sampleCount;
	resultAVL ["P (Hudson, Boos and Kaplan)"] = (sampleCount-pv[2])/sampleCount;
	resultAVL ["P (Hudson (S_nn))"] = (sampleCount-pv[3])/sampleCount;
	
	labels = {{"Cumulative Weight","F_ST"}};
	
	OpenWindow (CHARTWINDOW,{{"Hudson, Slatkin and Madison F_ST"}
			{"labels"}
			{"F_ST_1"}
			{"Step Plot"}
			{"F_ST"}
			{"Cumulative Weight"}
			{"F_ST"}
			{""}
			{"Cumulative Probability"}
			{"3"}
			{""}
			{"-1;-1"}
			{"10;1.309;0.785398"}
			{"Times:14:0;Times:12:0;Times:14:1"}
			{"16777215;16777215;16512;11776947;0;16777215;16711680;11842740;13158600;14474460;0;3947580;79;16744448;16777215;2984993;9199669;7018159;1460610;16748822;11184810;14173291;14173291"}
			{"16"}
			},
			"SCREEN_WIDTH/2-50;SCREEN_HEIGHT/2-100;50;50");
			
	OpenWindow (CHARTWINDOW,{{"Slatkin F_ST"}
			{"labels"}
			{"F_ST_2"}
			{"Step Plot"}
			{"F_ST"}
			{"Cumulative Weight"}
			{"F_ST"}
			{""}
			{"Cumulative Probability"}
			{"3"}
			{""}
			{"-1;-1"}
			{"10;1.309;0.785398"}
			{"Times:14:0;Times:12:0;Times:14:1"}
			{"16777215;16777215;16512;11776947;0;16777215;16711680;11842740;13158600;14474460;0;3947580;79;16744448;16777215;2984993;9199669;7018159;1460610;16748822;11184810;14173291;14173291"}
			{"16"}
			},
			"SCREEN_WIDTH/2-50;SCREEN_HEIGHT/2-100;50+SCREEN_WIDTH/2;50");
	
	OpenWindow (CHARTWINDOW,{{"Hudson, Boos and Kaplan F_ST"}
			{"labels"}
			{"F_ST_3"}
			{"Step Plot"}
			{"F_ST"}
			{"Cumulative Weight"}
			{"F_ST"}
			{""}
			{"Cumulative Probability"}
			{"3"}
			{""}
			{"-1;-1"}
			{"10;1.309;0.785398"}
			{"Times:14:0;Times:12:0;Times:14:1"}
			{"16777215;16777215;16512;11776947;0;16777215;16711680;11842740;13158600;14474460;0;3947580;79;16744448;16777215;2984993;9199669;7018159;1460610;16748822;11184810;14173291;14173291"}
			{"16"}
			},
			"SCREEN_WIDTH/2-50;SCREEN_HEIGHT/2-100;50+SCREEN_WIDTH/2;100+SCREEN_HEIGHT/2");
	
	labels = {{"Cumulative Weight","S_NN"}};
	
	OpenWindow (CHARTWINDOW,{{"Hudson S_NN"}
			{"labels"}
			{"F_ST_4"}
			{"Step Plot"}
			{"S_NN"}
			{"Cumulative Weight"}
			{"S_NN"}
			{""}
			{"Cumulative Probability"}
			{"3"}
			{""}
			{"-1;-1"}
			{"10;1.309;0.785398"}
			{"Times:14:0;Times:12:0;Times:14:1"}
			{"16777215;16777215;16512;11776947;0;16777215;16711680;11842740;13158600;14474460;0;3947580;79;16744448;16777215;2984993;9199669;7018159;1460610;16748822;11184810;14173291;14173291"}
			{"16"}
			},
			"SCREEN_WIDTH/2-50;SCREEN_HEIGHT/2-100;50;100+SCREEN_HEIGHT/2");
}

return resultAVL;