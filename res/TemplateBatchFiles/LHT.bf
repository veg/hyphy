
/* A Likelihood Ratio Test to detect conflicting phylogentic signal
  Huelsenbeck and Bull, 1996 */
/* implemented by Olivier Fedrigo Feb. 23 2006 */


fprintf(stdout,"\n\n+----------------------------------------------------------------------------------------------+\n");
fprintf(stdout,"| A Likelihood Ratio Test to detect conflicting phylogenetic signal Huelsenbeck and Bull, 1996 |\n");
fprintf(stdout,"|            implemented by Olivier Fedrigo Sep. 23 2006             |\n");
fprintf(stdout,"+----------------------------------------------------------------------------------------------+\n\n");
nreplicates=10;
fprintf(stdout,"\nHow many replicates?");
fscanf(stdin,"Number",nreplicates);
fprintf(stdout,"\n",nreplicates," bootstrap replicates will be performed\n");

SetDialogPrompt ("Please locate a partitioned NEXUS file:");
DataSet ds = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter filteredData = CreateFilter (ds,1);
fprintf (stdout, "\nLoaded a ", filteredData.species, " sequence alignment with ", filteredData.sites," sites\n");

SelectTemplateModel(filteredData);

/*check if there are any partitions.*/

partitionsFound = partitionsFound;

fprintf (stdout,"\n\n",partitionsFound," character sets found\n");
if (partitionsFound==0)
{
	fprintf (stdout, "\n\nPlease add character sets to your NEXUS file...");
	return 0;
}
else
{
	for (part=0; part<partitionsFound; part=part+1)
	{
	    fprintf (stdout,DATA_FILE_PARTITION_MATRIX[0][part],"\n");
	    fprintf (stdout,DATA_FILE_PARTITION_MATRIX[1][part],"\n\n");
	}
}

/*check if there are any trees*/
fprintf (stdout,"\n\n",Rows(NEXUS_FILE_TREE_MATRIX)," trees found\n");
if (partitionsFound>Rows(NEXUS_FILE_TREE_MATRIX))
{
    fprintf (stdout, "\n\nPlease add trees to your nexus file...");
    return 0;
}

ntree=Rows(NEXUS_FILE_TREE_MATRIX);
listOfTree	=	{ntree,1};
likList		=	{partitionsFound,1};
listOfLk	=	{partitionsFound,ntree};
bestLk		=	{partitionsFound,1};
bestTree	=	{partitionsFound,1};
bestLk2		=	{partitionsFound,1};
globalLk	=	{ntree,1};

/*GET LIKELIHOOD FOR TREES FOR ALL PARTITIONS AND CHOOSE THE BEST TREE*/    
for (part=0; part<partitionsFound; part=part+1)
{
  DataSetFilter dsf=CreateFilter (ds,1,DATA_FILE_PARTITION_MATRIX[1][part]);
  if (FREQUENCY_SENSITIVE) /*get the base frequency independently for each partition*/
  {
      HarvestFrequencies (baseFreqs, dsf, 1,1,1);
      MBF = PopulateModelMatrix ("filterModelMatrix", baseFreqs);
      Model filterModel = (filterModelMatrix, baseFreqs, MBF);
  }
  
  fprintf(stdout,"Analysis for partition ",DATA_FILE_PARTITION_MATRIX[0][part]," given trees:\n");
  for (i=0; i<ntree; i=i+1)
  {
		Tree uniqueTree=NEXUS_FILE_TREE_MATRIX[i][1];
		LikelihoodFunction lf = (dsf,uniqueTree);
		Optimize (mles,lf);
		if (part==0) 
			{fprintf(stdout,lf);}
		listOfLk[part][i]=mles[1][0];
		if (i==0)
		{
			bestLk[part]=mles[1][0];
			bestTree[part]=i;
		}
		else
		{
			if (bestLk[part]<mles[1][0])
			{
				bestLk[part]=mles[1][0];
				bestTree[part]=i;
			}
		}
	}
}

/* GET LK1 and LK2 AND CHOOSE THE BEST TREE*/
lk2=0;  /*H1: unconstrained different tree for different partition*/
for (part=0; part<partitionsFound; part=part+1) {lk2=lk2+bestLk[part];}

lk1=0;  /*H0: constrained same tree for different partition*/
for (i=0; i<ntree; i=i+1)
{
	globalLk[i]=0;
	for (part=0; part<partitionsFound; part=part+1) 
	{
		globalLk[i]=globalLk[i]+listOfLk[part][i];
	}
	if (i==0)
	{
		lk1=globalLk[i];
		theBestTree=i;
	}
	else
	{
		if (lk1<globalLk[i])
		{
			lk1=globalLk[i];
			theBestTree=i;
		}
	}
}

fprintf (stdout,"\n| Tree n.  |");
for (i=0;i<partitionsFound;i=i+1) 
{
	fprintf(stdout," ",DATA_FILE_PARTITION_MATRIX[0][i],"\t|");
}
fprintf (stdout," Total\n");
for (j=0; j<ntree; j=j+1)
{
	fprintf (stdout,"| ",NEXUS_FILE_TREE_MATRIX[j][0]," |");
	for (i=0;i<partitionsFound;i=i+1)
	{
		fprintf(stdout," ",Format(listOfLk[i][j],10,10));
		if (listOfLk[i][j]==bestLk[i]) {fprintf(stdout,"*\t|");} else {fprintf(stdout,"\t|");}
	}
	fprintf(stdout," ",Format(globalLk[j],10,10));
	if (globalLk[j]==lk1) {fprintf(stdout,"*\n");} else {fprintf(stdout,"\n");}
}

fprintf(stdout,"\n");
fprintf(stdout,"L0= ",Format(lk1,10,10),"\n");
fprintf(stdout,"L1= ",Format(lk2,10,10),"\n");
lnlikDelta= 2*(lk2-lk1);
fprintf(stdout,"LRT= ",Format(lnlikDelta,10,10),"\n\n");

/*PARAMETRIC BOOTSTRAP*/    
count=0;

for (replicate=0;replicate<nreplicates;replicate=replicate+1)
{
	lk1=0;
	lk2=0;
	for (part=0; part<partitionsFound; part=part+1)
	{
		DataSetFilter dsf=CreateFilter (ds,1,DATA_FILE_PARTITION_MATRIX[1][part]);
		if (FREQUENCY_SENSITIVE)
		{
			HarvestFrequencies (baseFreqs, dsf, 1,1,1);
			MBF = PopulateModelMatrix ("filterModelMatrix", baseFreqs);
			Model filterModel = (filterModelMatrix, baseFreqs, MBF);
		}
		Tree uniqueTree=NEXUS_FILE_TREE_MATRIX[theBestTree][1];
		LikelihoodFunction lf = (dsf,uniqueTree);
		Optimize (mles,lf);

		/*simulate constrained different partition*/
		DataSet simData = SimulateDataSet(lf);
		DataSetFilter simDataf=CreateFilter (simData,1);     
		if (FREQUENCY_SENSITIVE)
		{
			HarvestFrequencies (baseFreqs, simDataf, 1,1,1);
			MBF = PopulateModelMatrix ("filterModelMatrix", baseFreqs);
			Model filterModel = (filterModelMatrix, baseFreqs, MBF);
		}

		/*H0: constrained same tree for different partition*/ 
		Tree globalTree=NEXUS_FILE_TREE_MATRIX[theBestTree][1];
		LikelihoodFunction lfsim = (simDataf,globalTree);
		Optimize (mles,lfsim);
		lk1=lk1+mles[1][0];  

		/*H1: unconstrained difffeent tree for different partition*/
		for (i=0; i<ntree; i=i+1)
		{
			Tree uniqueTree=NEXUS_FILE_TREE_MATRIX[i][1];
			LikelihoodFunction lfsim2 = (simDataf,uniqueTree);
			Optimize (mles,lfsim2);
			if (i==0) {bestLk=mles[1][0];} else {bestLk=Max(bestLk,mles[1][0]);}
		}
		lk2=lk2+bestLk;
	}
	simLRT = 2*(lk2-lk1);
	if (simLRT>0)
	{
		fprintf(stdout,"replicate ",replicate,"- ",simLRT,"\n");
	}
	else
	{
		fprintf(stdout,"replicate ",replicate,"- 0.00000 *",simLRT,"\n");
	}
	if (lnlikDelta>simLRT) 
	{
		count = count+1;
	}
}
/*If your observed lnlikDelta is greater than 95% of the simulated LRT's, then reject the null hypothesis (that the partitions share the same tree) at the 5% level. negative LRTs are forced to be 0*/
fprintf(stdout,"\n\n*** P-value (Parametric BS)= ",Format(1.0-(count/nreplicates),10,10),"\n");
