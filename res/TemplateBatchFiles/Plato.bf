/************************************************************************************/
/* A Likelihood Method for the Detection of Selection and Recombination using Nucleotide Sequences 
	-Grassly and Holmes, 1997
  HBL implementation by Olivier Fedrigo (ofedrigo@duke.edu) 
  October 2006
*/

RequireVersion ("0.9920060901");

VERBOSITY_LEVEL		  = -1; /* HELPS REDUCE GUI REFRESH OVERHEAD */
ACCEPT_BRANCH_LENGTHS = 1;
ACCEPT_ROOTED_TREES   = 1;
nreps=100;

fprintf			(stdout, "Partial Likelihood Anomalies detected Through Optimization - PLATO\nVersion 2.0\n(c) Copyright, 1998 Nick Grassly and Andrew Rambaut\nDepartment of Zoology, University of Oxford\nSouth Parks Road, Oxford OX1 3PS, U.K.\n");

SetDialogPrompt ("Please locate an alignment file:");
DataSet 			nucleotideSequences = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter       filteredData 	    = CreateFilter (nucleotideSequences,1);
HarvestFrequencies (observedFreqs, filteredData, 1, 1, 1);

fprintf (stdout, "\nLoaded a ", filteredData.species, " sequence alignment with ", filteredData.sites,"\n");
fprintf (stdout, "\nBase composition:\n\tA: ", Format (observedFreqs[0],10,5),
									"\n\tC: ", Format (observedFreqs[1],10,5),
									"\n\tG: ", Format (observedFreqs[2],10,5),
									"\n\tT: ", Format (observedFreqs[3],10,5), "\n\n");
									
ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"queryTree.bf");
SelectTemplateModel(filteredData);

/*CALCULATE ML PARAMETERS AND LIKELIHOOD PER SITES*/
Tree givenTree				= treeString;
LikelihoodFunction theLnLik = (filteredData, givenTree);
Optimize 					  (paramValues, theLnLik);

/* SLKP: NEED TO SAVE GLOBAL VARIABLES; OTHERWISE THEY WILL BE OVERWRITTEN
   DURING OPTIMIZATIONS FROM SIMULATED DATA */
   
GetString 			  (likelihoodInfo, theLnLik, -1);   
globalVariableList 	= likelihoodInfo["Global Independent"];
globalVariableCount	= Columns (globalVariableList);

if (globalVariableCount)
{
	stashed_GV = {globalVariableCount,1};
	for (vc = 0; vc < globalVariableCount; vc = vc + 1)
	{
		ExecuteCommands ("stashed_GV[vc]="+globalVariableList[vc]+";");
	}
}


ConstructCategoryMatrix		  (L, theLnLik, COMPLETE);

/* SLKP: If there are category variables involved in the model, one needs to collapse site likelihoods 
   conditional on the value of the category into an averaged value at a site */
 
if (Rows(L)>1)
{
	catVars 	= likelihoodInfo["Categories"];
	catVarCount = Columns 		(catVars);
	
	if (catVarCount > 1)
	{
		fprintf (stdout, "ERROR: only one category variable can be handled by this code\n");
		return 0;
	}

	ExecuteCommands ("GetInformation (catWeight,"+catVars[vc]+");");
	L = catWeight[1][-1]*L;
}

/* END SLKP */

fprintf						  (stdout,"\nML parameters estimated. Log(L) = ", paramValues[1][0], "\n");

/*CALCULATE THE RATIO: (SUM OF LIKELIHOOD INSIDE THE WINDOW PER SITE)/(SUM OF LIKELIHOOD OUTSITE THE WINDOW PER SITE)*/
/*FOR EACH WINDOW SIZE (FROM 5bp TO HALF OF THE DATASET*/
fprintf(stdout,"Actual surface \n");

smin			=	5;
smax			=	filteredData.sites$2;
liktable		=	getSurface(smin,smax,filteredData.sites,paramValues[1][0],L);

/*DO THE SAME THING FOR SIMULATIONS. 100 SIMULATION, FOR EACH ITERATE TAKE THE MAXIMAL VALUE FOR A GIVEN WINDOW SIZE*/
fprintf			(stdout,"Simulating \n");
winLikList		=	{nreps,smax};
for (simCounter=0;simCounter<nreps;simCounter=simCounter+1)
{
	Tree simTree=treeString;
	ClearConstraints(simTree);
	/* SLKP: NEED TO RESTORE GLOBAL VARIABLES BEFORE SIMULATION */
	if (globalVariableCount)
	{
		for (vc = 0; vc < globalVariableCount; vc = vc + 1)
		{
			ExecuteCommands (globalVariableList[vc]+"=stashed_GV[vc];");
		}
	}	 
	DataSet simData				=	SimulateDataSet(theLnLik);
	DataSetFilter simFilter		=	CreateFilter(simData,1);
	HarvestFrequencies				(simFreqs,simFilter,1,1,1);
	LikelihoodFunction simLik 	= (simFilter,simTree,simFreqs);
	Optimize (simParamValues,simLik);
	ConstructCategoryMatrix (Lsim,simLik, COMPLETE);
	if (catVarCount == 1)
	{
		ExecuteCommands ("GetInformation (catWeight,"+catVars[vc]+");");
		Lsim = catWeight[1][-1]*Lsim;
	}
	fprintf(stdout,"Replicate #",simCounter,"\nTree -ln Likelihood = ",simParamValues[1][0],"\n");
	liktableTemp=getSurface(smin,smax,simFilter.sites,simParamValues[1][0],Lsim);
	for (sp=0;sp<simFilter.sites-smin;sp=sp+1)
	{
		for (s=0;s<smax-smin+1;s=s+1) {winLikList[simCounter][s]=Max(liktableTemp[s][sp],winLikList[simCounter][s]);}
	}
}

fprintf(stdout,"\nSimulations done\n");

/*CALCULATE MEAN AND VARIANCE OF SIMUALTIONS FOR EACH WINDOW SIZE*/
fprintf(stdout,"Calculate mean and variance \n");
mean		=	{smax,1};
variance	=	{smax,1};
sumsq=0.0; sums=0.0;

for (i=0;i<(smax-smin+1);i=i+1)
{
	sumsq=0.0; sums=0.0;
	for(j=0;j<nreps;j=j+1)
	{
		sums=sums+winLikList[j][i];
		sumsq=sumsq+((winLikList[j][i])*(winLikList[j][i]));
	}
	variance[i]=((sumsq-((sums*sums)/nreps))/(nreps-1));
	mean[i]=sums/nreps;
	windowsize=i+smin;
}

/*CALCULATE Z-VALUES*/       
fprintf(stdout,"Results \n");
alpha=0.05/(smax-smin+1);
Z={smax,filteredData.sites};
for(i=0;i<(smax-smin+1);i=i+1)
  {
     for(j=0;j<(filteredData.sites+1-i-smin);j=j+1)
          {
        Z[i][j]=(liktable[i][j]-mean[i])/(Sqrt(variance[i]));
        }
      }
fprintf(stdout,"\nZ values calculated\n");

/*DETERMINE CUTOFF*/
fprintf(stdout, "\nBonferroni-corrected significance level for alpha=0.05: ", alpha,"\n");
Z_cutOff=-normZval(alpha);
fprintf(stdout, "Z-values greater than ",Z_cutOff," are significant\n\n");

/* GET THE BEST SCORE*/
temp=0.0;
for (i=smin;i<=smax;i=i+1)
{
	for (j=0;j<(filteredData.sites+1-i);j=j+1)
	{
		if (Z[i-smin][j]>temp && Z[i-smin][j]>Z_cutOff)
		{
			temp=Z[i-smin][j];
			size=i;
			sp=j;
		}
	}

}

map={filteredData.sites,1};
if (temp>Z_cutOff) /*IF THE BEST Z-SCORE IS SMALLER THAN THE CUTOFF THEN NO NEED TO DO IT*/
{
	while (temp>Z_cutOff)
	{
		for (i=sp;i<(sp+size);i=i+1) {map[i]=1;}
		fprintf(stdout,Format(sp+1,5,0)," - ",Format(sp+size,5,0)," : ",Format(Z[size-smin][site],6,2),"\n");
		temp=0.0;
		for (i=smin;i<=smax;i=i+1)
		{
			for (j=0;j<(filteredData.sites-i+1);j=j+1)
			{
			if (Z[i-smin][j]>temp && Z[i-smin][j]>Z_cutOff)
			{
				check=0;
				for (n=j;n<(j+i);n=n+1) {check=check+map[n];}
				if (check==0) /*RECORD A GOOD SCORE ONLY IF NONE OF ITS SITES HAS BEEN INCLUDED IN AN PREVIOUS BETTER SCORE*/
				{
					temp=Z[i-smin][j];
					size=i;
					sp=j;
				}
			}
			}
		}
	}
}
else
{
    fprintf(stdout, "No anomalous regions found\n");
}

fprintf(stdout,"\nFinished\n");

        
/************************************************************************************/

function normZval(P)
  {
     /*modified from Fortran algorithm AS241
    APPL. STATIST. (1988) VOL. 37, NO. 3, 477-484.
 Produces the normal deviate Z corresponding to a given lower
       tail area of P; Z is accurate to about 1 part in 10e7.*/    
  SPLIT1 = 0.425;
  SPLIT2 = 5.0;
  CONST1 = 0.180625;
  CONST2 = 1.6;
  /*Coefficients for P close to 0.5*/
  A0 = 3.3871327179;
     A1 = 50.434271938;
  A2 = 159.29113202;
  A3 = 59.109374720;
  B1 = 17.895169469;
  B2 = 78.757757664;
  B3 = 67.187563600;
  /*Coefficients for P not close to 0, 0.5 or 1.*/
  C0 = 1.4234372777;
   C1 = 2.7568153900;
  C2 = 1.3067284816;
  C3 = 0.17023821103;
  D1 = 0.73700164250;
  D2 = 0.12021132975;
  /*Coefficients for P near 0 or 1.*/
      E0 = 6.6579051150;
      E1 = 3.0812263860;
  E2 = 0.42868294337;
  E3 = 0.017337203997;
  F1 = 0.24197894225;
  F2 = 0.012258202635;
      Q=P-0.5;
        if(Abs(Q)<=SPLIT1)
          {
              Rx=CONST1-(Q*Q);
                x=Q*((((((A3*Rx)+A2)*Rx)+A1)*Rx)+A0)/((((((B3*Rx)+B2)*Rx)+B1)*Rx)+1.0);
                return x;
         }
  else {if(Q<0.0) {Rx=P;} else {Rx=1.0-P;}}
  if(Rx<=0.0) {return 0.0;}    
         Rx=Sqrt(-Log(Rx));
        if(Rx<=SPLIT2)
         {
        Rx=Rx-CONST2;
        x=(((C3*Rx+C2)*Rx+C1)*Rx+C0)/((D2*Rx+D1)*Rx+1.0);
        }
       else
         {
              Rx=Rx-SPLIT2;
    x=(((E3*Rx+E2)*Rx+E1)*Rx+E0)/((F2*Rx+F1)*Rx+1.0);
      }
       if(Q<0.0) {x = -x;}
   return x;
  }

/************************************************************************************/

function getSurface(windowMin,windowMax,nsites,MLS,Likelihood)
{
    /*CALCULATE THE RATIO: (SUM OF Likelihood INSIDE THE WINDOW PER SITE)/(SUM OF Likelihood OUTSITE THE WINDOW PER SITE)*/
  	/*FOR EACH WINDOW SIZE (FROM windowMin TO windowMax FOR A DATASET WITH nsites AND MLS*/
  	
  	timer 		= 	Time(1);
    temp		=	{windowMax,nsites};
    loopUB  	= 	nsites-windowMin+1;
    windowSpan 	=   windowMax-windowMin+1;
    loggedL		=   Log (Likelihood);
    
    for 	(sp=0;sp<loopUB;sp=sp+1)
	{
		/*
		window	=	0;
		loopUB2 = windowMin+sp;
		for (s=sp;s<loopUB2;s=s+1) 
		{	
			window=window+Log(Likelihood[s]);
		}*/
		
		/* SLKP: matrix hackery to do the same as above; 
		   speeds things up quite a bit
		   see Examples/BatchLanguage/MatrixIndexing.bf */
		   
		logExtract = loggedL[{{0,sp}}][{{0,windowMin+sp-1}}];
		window  = (logExtract*Transpose(logExtract["1"]))[0];
		
		temp[0][sp]= 	window/windowMin / ((MLS-window) / (nsites-windowMin));
		
		localWS = Min (windowSpan, nsites-sp-windowMin-1);
		i		   =	1;
		while (i<=localWS)
		{
			cww			=   windowMin+i;
			window		=	window+loggedL[sp+cww];
			temp[i][sp]	=	(window/cww)/((MLS-window)/(nsites-cww));
			i			=	i+1;
		}
	}
	
	fprintf(stdout,"\nML surface done in ", Time(1)-timer, " seconds\n");      
	return temp;
}
/************************************************************************************/
