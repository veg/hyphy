/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2006  
Primary Development:
  Sergei L Kosakovsky Pond (sergeilkp@mac.com)
Significant contributions from:
  Spencer V Muse (muse@stat.ncsu.edu)
  Simon DW Frost (sdfrost@ucsd.edu)
  Art FY Poon    (apoon@biomail.ucsd.edu)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "likefunc.h"
#include "calcnode.h"
#include "site.h"
#include "batchlan.h"
#include "category.h"


#ifdef __WINDOZE__
	#include     "windows.h"
#endif

#ifdef __MAC__
	#include	 "timer.h"
#endif

#ifdef		__MACPROFILE__
	#include "profiler.h"
#endif

#ifdef	_SLKP_LFENGINE_REWRITE_
	#include "scfg.h"
#endif

#if !defined __UNIX__ || defined __HEADLESS__
	#ifndef __HEADLESS__
		#include 	 "HYTreePanel.h"
		extern		 _HYTreePanel*  feedbackTreePanel;
	#endif
	
	_Parameter 	 nicetyLevel 	= 	1.0;
	_String nicetyMacLevel   ("NICETY_LEVEL");	
	
	long	 	 siteEvalCount	=	0, 
				 divideBy 		= 	10000000;
				 
	 void	DecideOnDivideBy (_LikelihoodFunction*);

#endif



#ifdef __MP__
	#ifndef __MACHACKMP__
		#include <pthread.h>
	#else
		#include "mypthread.h"
	#endif
	struct	 ThreadReleafTask {
		_TheTree 	*    t;
		
		 long		*siteList,
					*leafList;
		_SimpleList	*freqs;
					
		_DataSetFilter*
					df;
					
		long		startAt,
					endAt,
					threadIndex;
					
		_Parameter  result,
				  * resultsStorage;
	};
	
	struct	 WancReleafTask {
		_TheTree 	*tree;
		
		long		startAt,
					endAt,
					*doneSites,
					*lastDone,
					totalUniqueSites,
					threadIndex;
					
		_DataSetFilter*
					dsf;
					
		_List	   *dupList;
		_Formula   *fla;
										
	};

	pthread_mutex_t wancMutex = PTHREAD_MUTEX_INITIALIZER;

	pthread_t*		theThreads = nil;
	ThreadReleafTask* theTasks = nil;
	void*			ThreadReleafFunctionAA (void*);
#endif

#ifdef	__HYPHYMPI__

	_List			parallelOptimizerTasks;
	
	_Matrix			varTransferMatrix,
					resTransferMatrix,
					catFreqMatrix;
				
	long			MPICategoryCount,
					transferrableVars;
					
	bool			lfAutoMPI = false;

	extern			int		  _hy_mpi_node_rank;

#endif

extern	_Parameter  machineEps;

#define SQR(A) (A)*(A)
#define GOLDEN_RATIO 1.618034

#define	  PERTURBATION_OF_ZERO 	  0.0

long	  likeFuncEvalCallCount = 0,
		  systemCPUCount 		= 1,
		  lockedLFID	 		= -1;

#ifndef  __HYALTIVEC__
#define	 STD_GRAD_STEP 1.0e-7
#else
#define	 STD_GRAD_STEP 5.0e-6
#endif

#ifdef 	  __HYPHYDMALLOC__
	#include "dmalloc.h"
#endif


// some string constants

_String
	
	globalStartingPoint 			("GLOBAL_STARTING_POINT"),
	randomStartingPerturbations 	("RANDOM_STARTING_PERTURBATIONS"),
	optimizationPrecision 			("OPTIMIZATION_PRECISION"),
	startingPrecision 				("STARTING_PRECISION"),
	optimizationMethod 				("OPTIMIZATION_METHOD"),
	useLastResults 					("USE_LAST_RESULTS"),
	allowBoundary 					("ALLOW_BOUNDARY"),
	bracketingPersistence			("BRACKETING_PERSISTENCE"),
	intermediatePrecision			("INTERMEDIATE_PRECISION"),
	keepOptimalOrder 				("KEEP_OPTIMAL_ORDER"),
	skipOmissions 					("SKIP_OMISSIONS"),
	optimizeSummationOrder 			("OPTIMIZE_SUMMATION_ORDER"),
	optimizePartitionSize 			("OPTIMIZE_SUMMATION_ORDER_PARTITION"),
	maximumIterationsPerVariable 	("MAXIMUM_ITERATIONS_PER_VARIABLE"),
	optimizationPrecisionMethod 	("OPTIMIZATION_PRECISION_METHOD"),
	relativePrecision 				("RELATIVE_PRECISION"),
	likefuncOutput					("LIKELIHOOD_FUNCTION_OUTPUT"),
	dataFilePrintFormat				("DATA_FILE_PRINT_FORMAT"), 
	dataFileDefaultWidth			("DATA_FILE_DEFAULT_WIDTH"),
	dataFileGapWidth				("DATA_FILE_GAP_WIDTH"),
	categorySimulationMethod		("CATEGORY_SIMULATION_METHOD"),
	useInitialDistanceGuess			("USE_DISTANCES"),
	randomSeed						("RANDOM_SEED"),
	assignedSeed					("ASSIGNED_SEED"),
	covariancePrecision				("COVARIANCE_PRECISION"),
	cacheSubtrees 					("CACHE_SUBTREES"),
	likeFuncCountVar 				("LF_CALL_COUNT"),
	doShuffleOrder					("SHUFFLE_ORDER_OF_PARAMETERS"),
	forceDistanceEstimates			("FORCE_DISTANCE_ESTIMATES"),
	useDuplicateMatrixCaching		("USE_DUPLICATE_MATRIX_CACHING"),
	siteWiseMatrix					("SITE_LIKELIHOOD"),
	blockWiseMatrix					("BLOCK_LIKELIHOOD"),
	useFullMST						("USE_MST_HEURISTIC"),
	stateCountMatrix				("STATE_COUNT_MATRIX"),
	wStateCountMatrix				("WSTATE_COUNT_MATRIX"),
	tryNumericSequenceMatch			("TRY_NUMERIC_SEQUENCE_MATCH"),
	allowSequenceMismatch			("ALLOW_SEQUENCE_MISMATCH"),
	shortMPIReturn					("SHORT_MPI_RETURN"),
	mpiPrefixCommand				("MPI_PREFIX_COMMAND"),
	skipConjugateGradient			("SKIP_CONJUGATE_GRADIENT"),
	useIntervalMapping				("USE_INTERVAL_MAPPING"),
	intervalMappingMethod			("INTERVAL_MAPPING_METHOD"),
	useAdaptiveVariableStep			("USE_ADAPTIVE_VARIABLE_STEP"),
	storeRootSupportFlag			("STORE_ROOT_SUPPORT"),
	supportMatrixVariable			("SUPPORT_MATRIX_LIST"),
	optimizationStatusFile			("SAVE_OPT_STATUS_TO"),
	autoParalellizeLF				("AUTO_PARALLELIZE_OPTIMIZE");
	
	

extern _String useNexusFileData,
			   VerbosityLevelString;


void		StashMSTStuff 			 (_Parameter**, long **, char **, long, long);
void		countingTraverse 		 (node<long>*, long&, long, long&, bool);
void		countingTraverseArbRoot  (node<long>*, node<long>*, long&, long, long&);
long		findAvailableSlot 		 (_SimpleList&, long&);
void		setComputingArrays 		 (node<long>*, node<long>*, _SimpleList&, _SimpleList&, _SimpleList &, _SimpleList&, _SimpleList&, long&);

#ifdef 		__MP__
	void*	StateCounterMP 			 (void*);
#endif	
_SimpleList	Fibonacci;

_Parameter	go2Bound = 0, 
			precision,
			optimizationPrecMethod,
			maxItersPerVar,
			relPrec,
			dFPrintFormat,
			dFDefaultWidth,
			assignedSeedVal = -1.0,
			categorySimMethod;
			
bool		forceRecomputation = false, 
			isInOptimize 	   = false,
			usedCachedResults  = false;

long		bracketFCount 	= 0, 
			bracketCount 	= 0, 
			oneDFCount 		= 0, 
			oneDCount 		= 0,
			categID 		= 0, 
			offsetCounter 	= 1;

extern		long	
			matrixExpCount, 
			taylorTermsCount, 
			squaringsCount;
			
extern		_List	
			dataSetList,
			likeFuncList;

bool		CheckOneDStep 	(_Parameter&,_Parameter,_Parameter);
bool 		CheckEqual 		(_Parameter, _Parameter);

_Variable*  siteWiseVar     = nil,
		 *  blockWiseVar	= nil;
		 
		 
_String  *  progressFileString = nil;
		 
node<long>* DepthWiseStepTraverserLevel  (long&, node<long>* root);
_Parameter  myLog (_Parameter);
_Parameter	SumUpHiddenMarkov (_Matrix&, _Matrix&, _Matrix&, _SimpleList &, long);

		 
#ifdef 		__UNIX__

void		UpdateOptimizationStatus (_Parameter, long, char, bool, _String * fileName = nil);

//__________________________________________________________________________________

void		UpdateOptimizationStatus (_Parameter max, long pdone, char init, bool optimization, _String* fileName)
{
	static long    lCount;
	static time_t  startedAt;
	static clock_t userTimeStart;
	static long	   lastDone;
	FILE		   *outFile = fileName?doFileOpen (fileName->sData,"w"):nil;
	
	if (init==0)
	{
		lCount = likeFuncEvalCallCount;
		time    (&startedAt);
		setvbuf (stdout,nil, _IONBF,1);
		lastDone = 0;
		userTimeStart = clock();
	}
	else
	if (init==1)
	{
		time_t timeNow;
		time (&timeNow);
		if (timeNow == startedAt) timeNow++;
		if (pdone<0)
			pdone = lastDone;
		if (optimization)
		{
			if (outFile)
				fprintf (outFile,"Likelihood function optimization status\nCurrent Maximum: %-14.8g (%d %% done)\nLikelihood Function evaluations/second: %-8.4g", (double)max, pdone, 
						(likeFuncEvalCallCount-lCount)/(double) (timeNow-startedAt));			
			else
			{
				printf ("\033\015 Current Max: %-14.8g (%d %% done) LF Evals/Sec: %-8.4g", (double)max, pdone, 
						(likeFuncEvalCallCount-lCount)/(double) (timeNow-startedAt));
	
				if (timeNow-startedAt)
					printf ("CPU Load: %-8.4g", (clock()-userTimeStart)/((_Parameter)CLOCKS_PER_SEC*(timeNow-startedAt)));
			}
		}
		else
			printf ("\033\015 Sites done: %g (%d %% done)", (double)max, pdone);		

		lastDone = pdone;
	}
	else
	{
		if (outFile)
			fprintf (outFile,"DONE");
		else
		{
			printf ("\033\015 ");
			setvbuf (stdout,nil,_IOLBF,1024);
		}	
	}
	if (outFile)
		fclose (outFile);
}

#else

//__________________________________________________________________________________


void		 DecideOnDivideBy (_LikelihoodFunction* lf)
{
	long		 alterIndex = 0;
	
	if		(lf->HasComputingTemplate())
	{
		for (long k=0; k<lf->GetIndependentVars().lLength;k++)
			if (!LocateVar (lf->GetIndependentVars().lData[k])->IsGlobal())
			{
				alterIndex = k;
				break;
			}
	}
	
	
#ifndef _SLKP_LFENGINE_REWRITE_	
	long				   stash1 = siteEvalCount;
#endif
#ifdef	_OPENMP
	lf->SetThreadCount (1);
#endif	
	TimerDifferenceFunction (false);
	lf->SetIthIndependent (alterIndex,lf->GetIthIndependent(alterIndex));
	lf->Compute			  ();

	
#ifdef _SLKP_LFENGINE_REWRITE_	
	_Parameter			  tdiff = TimerDifferenceFunction(true);
#ifdef	_OPENMP
	if (systemCPUCount > 1)
	{
		_Parameter			minDiff = tdiff;
		long				bestTC  = 1;
		
		for (long k = 2; k <= systemCPUCount; k++)
		{
			lf->SetThreadCount				(k);
			TimerDifferenceFunction			(false);
			lf->SetIthIndependent			(alterIndex,lf->GetIthIndependent(alterIndex));
			lf->Compute						();
			tdiff = TimerDifferenceFunction (true);
			if (tdiff < minDiff)
			{
				minDiff = tdiff;
				bestTC	= k;
			}
			else
				break;
		}
		lf->SetThreadCount (bestTC);
		divideBy			  = 0.5 / minDiff;		
		ReportWarning		(_String("Auto-benchmarked an optimal number (") & bestTC & ") of threads.");
	}
	else
#endif
		divideBy			  = 0.5 / tdiff;
	ReportWarning		(_String("Set GUI update interval to every ") & divideBy & "-th LF evaluation.");
		
#else	
	divideBy			  = (siteEvalCount-stash1) * 0.5 / TimerDifferenceFunction(true);
#endif
	
}

#endif

//__________________________________________________________________________________

_Parameter myLog (_Parameter arg)
{
	return (arg>0.0)?log(arg):-1000000.;
}

//__________________________________________________________________________________

_Parameter			SumUpHiddenMarkov (_Matrix& hmc, _Matrix& hmm, _Matrix& hmf, _SimpleList & duplicateMap, long bl)
{
	long      		   ni = hmm.GetHDim(),
					   k,
					   m,
					   mi;
					   					   
	_Matrix	  		   temp  (ni,1,false,true),
					   temp2 (ni,1,false,true);
					  
	_Parameter		   *p1 = hmc.fastIndex(),
					   scrap,
					   min;
			
	long		       correctionFactor = 0; // correction factor
	
	p1 = hmc.fastIndex();
	mi = duplicateMap.lData[duplicateMap.lLength-1];
	
	for (m=0; m<ni; m++,mi += bl)
		temp2.theData[m] = p1[mi];
	
	for (long i=duplicateMap.lLength-2; i>=0; i--)
	{
		min = 1.e300;
		for (k=0; k<ni;k++)
		{
			mi    = duplicateMap.lData[i];
			scrap = 0.0;
			
			for (m=0; m<ni; m++,mi += bl)
				scrap += hmm.theData[k*ni+m] * p1[mi] * temp2.theData[m];
				
			temp.theData[k] = scrap;
			
			if (scrap<min)
				min = scrap;
		}
		
		if (min <= 0.0)
			return -1000000.;
			
		k = log (min)-1;
		correctionFactor -= k;
		scrap = exp ((_Parameter)(-k));
		for (k=0; k<ni; k++)
			temp.theData[k] *= scrap;
		
		_Parameter* swap = temp.theData;
		temp.theData     = temp2.theData;
		temp2.theData    = swap;
	}
	scrap = 0.0;
	
	for (k=0; k<ni; k++)
		scrap += temp2.theData[k] * hmf.theData[k];
	
	return myLog(scrap) - correctionFactor;
}

//_______________________________________________________________________________________

_LikelihoodFunction::_LikelihoodFunction (void)
{
	Init();
}

//_______________________________________________________________________________________

void _LikelihoodFunction::Init (void)
{
	siteResults 		= nil;
	bySiteResults 		= nil;
	hasBeenOptimized 	= false;
	hasBeenSetUp 		= 0;
	templateKind        = false;
	computingTemplate   = nil;
	mstCache			= nil;
	nonConstantDep      = nil;
	
#ifdef	_SLKP_LFENGINE_REWRITE_
	conditionalInternalNodeLikelihoodCaches = nil;
	conditionalTerminalNodeStateFlag		= nil;
	siteScalingFactors						= nil;
	branchCaches							= nil;
#ifdef	_OPENMP
	SetThreadCount		(systemCPUCount);
#endif
#endif	
	
}

//_______________________________________________________________________________________

_LikelihoodFunction::_LikelihoodFunction (_String& s) // from triplets
//format: datasetfilter name, tree name, frequency matrix name; etc...
{
	Init();
	Construct(s,nil);
}

//_______________________________________________________________________________________

bool	_LikelihoodFunction::MapTreeTipsToData (long f, bool leafScan) // from triplets
{
	_TheTree*  		t = (_TheTree*)LocateVar (theTrees.lData[f]);
	_CalcNode* 		travNode = t->StepWiseTraversal(true);
	_DataSetFilter* df = (_DataSetFilter*)dataSetFilterList.lData[theDataFilters.lData[f]];
	long 			dfDim = df->GetDimension(true);
	
	_List 	   		tips;
	
	while (travNode)
	{
		if (t->IsCurrentNodeATip())
		{
			_String tipName (travNode->GetName()->Cut(travNode->GetName()->FindBackwards('.',0,-1)+1,-1));
			tips&& &tipName;
		}
		if (!t->IsCurrentNodeTheRoot())
		{
			if (!travNode->GetModelMatrix())
			{
				_String warnMsg ("Model is not associated with the node:");
				warnMsg = warnMsg&*travNode->GetName();
				acknError (warnMsg);
				return false;
			}				
			else
				if (travNode->GetModelMatrix()->GetHDim()!=dfDim)
				{
					_String warnMsg ("The dimension of the transition matrix at node ");
					warnMsg = warnMsg&travNode->GetName()->Cut(travNode->GetName()->FindBackwards('.',0,-1)+1,-1)
						&" is not equal to the state count in the data filter associated with the tree.";
					acknError (warnMsg);
					return false;
				}	
		}			
		travNode = t->StepWiseTraversal(false);		
	}
	// now that "tips" contains all the names of tree tips we can 
	// scan thru the names in the datafilter and check whether there is a 1-1 match
	if ((t->IsDegenerate()?2:tips.lLength)!=df->NumberSpecies())
	{
		_String warnMsg ("The number of tree tips in ");
		warnMsg = warnMsg&*t->GetName()& " (" & _String((long)tips.lLength) & ") is not equal to the number of species in the data filter associated with the tree " &
				 '(' & _String((long)df->NumberSpecies()) & ")." ;
		acknError (warnMsg);
		return false;
	}
	
	if (!t->IsDegenerate())
	{
		long		j,
					k;
		
		_SimpleList tipMatches;
		// produce a sorted list of sequence names
		
		j = df->FindSpeciesName (tips, tipMatches);		
		
		if (j!=tips.lLength) // try numeric match
		{
			_Parameter 		doNum = 0.0;
			checkParameter (tryNumericSequenceMatch, doNum, 0.0);
			if (doNum>0.5)
			{
				long sj = j;
				for (j=0;j<tips.lLength;j++)
				{
					_String *thisName = (_String*)tips(j);
					k = atoi (thisName->sData);
					_String tryAgain (k);
					if ((tryAgain.Equal(thisName))&&(k<=tips.lLength))
						tipMatches<<k;
					else
						break;
				}			
				if (j==tips.lLength)
				{
					if (tipMatches.Find(0)==-1) // map to indexing from 0
						for (j=0;j<tips.lLength;j++)
							tipMatches.lData[j]--;
							
					_SimpleList *dfMap = (_SimpleList*)df->GetMap();
					if (dfMap)
					{
						for (sj = 0; sj < tips.lLength; sj++)
							tipMatches.lData[sj] = dfMap->lData[tipMatches.lData[sj]];
					}	
				}
				else
					j=sj;
			}
		}
		
		if (j==tips.lLength) // all matched
		{
			_String warnMsg ("The tips of the tree:"),*s = (_String*)tipMatches.toStr();
			warnMsg = warnMsg & *t->GetName() &" were matched with the species names from the data as follows "&*s;
			DeleteObject(s);
			ReportWarning (warnMsg);
			/*for (j=0; j<tips.lLength;j++)
			{
				tipMatches.lData[j] = df->theNodeMap.lData[tipMatches.lData[j]];
			}*/
			df->SetMap(tipMatches);
		}
		else
		{
			_String warnMsg ("The leaf of the tree:");
			warnMsg = warnMsg & *t->GetName() &" labeled " &*(_String*)tips(j)&" had no match in the data. Please make sure that all leaf names correspond to a sequence name in the data file.";
			_Parameter asmm = 0.0;
			checkParameter (allowSequenceMismatch, asmm, 0.0);
			if (asmm<.5)
			{
				WarnError (warnMsg);
				return false;
			}
			ReportWarning (warnMsg);
		}
	}
	if (leafScan)
	{
		((_SimpleList*)leafSkips(f))->Clear();
		df->MatchStartNEnd(*(_SimpleList*)optimalOrders(f),*(_SimpleList*)leafSkips(f));
		t->BuildINodeDependancies();
	}
	return true;
}

//_______________________________________________________________________________________

bool	_LikelihoodFunction::UpdateFilterSize (long f) // from triplets
{
	_TheTree*		t = (_TheTree*)LocateVar (theTrees.lData[f]);
	_CalcNode* 		travNode = t->StepWiseTraversal(true);
	_DataSetFilter* df = (_DataSetFilter*)dataSetFilterList.lData[f];
	
	_List 	   tips;
	while (travNode)
	{
		if (t->IsCurrentNodeATip())
		{
			_String tipName (travNode->GetName()->Cut(travNode->GetName()->FindBackwards('.',0,-1)+1,-1));
			tips&& &tipName;
		}
		travNode = t->StepWiseTraversal(false);		
	}
	
	if (!t->IsDegenerate())
	{
		long	j;
		
		_SimpleList tipMatches;
		_List*		specNames = &df->GetData()->GetNames();
		
		for (j=0;j<tips.lLength;j++)
		{
			long k = specNames->Find((_String*)tips(j));
			if	 (k==-1) break;
			tipMatches<<k;
		}
		if (j==tips.lLength) // all matched
		{
			_SimpleList  sortedList, 
						 vertPart, 
						 theExclusions;
			
			long						unitSize = df->GetUnitLength();
			sortedList.Duplicate		(&tipMatches);
			sortedList.Sort				();
			vertPart.Duplicate			(&df->theOriginalOrder);
			theExclusions.Duplicate		(&df->theExclusions);
			df->SetFilter				(df->GetData(),unitSize,sortedList,vertPart,false);
			df->SetMap				    (tipMatches);
			df->FilterDeletions		    (&theExclusions);
			df->theExclusions.Duplicate (&theExclusions);
			df->SetupConversion		    ();
			sortedList.Clear		    ();
			
			_SimpleList*				theOO = (_SimpleList*)optimalOrders (f),
					   *				theLS = (_SimpleList*)leafSkips (f);
			
			theOO->Clear();
			theLS->Clear();
			theOO->Populate (df->theMap.lLength/unitSize,0,1);
			df->MatchStartNEnd (*theOO,*theLS);
		}
		else
		{
			return false;
		}
	}
	return true;
}
//_______________________________________________________________________________________

void	 _LikelihoodFunction::Rebuild (void)
{
	blockDependancies.Clear();
	computationalResults.Clear();
	hasBeenSetUp     = 0;
	hasBeenOptimized = false;
	Cleanup();
	RescanAllVariables();
	Setup();
} 
//_______________________________________________________________________________________

void	 _LikelihoodFunction::Clear (void)
{
	DeleteCaches  ();
	theTrees.Clear();
	theDataFilters.Clear();
	theProbabilities.Clear();
	indexInd.Clear();
	indexDep.Clear();
	indexCat.Clear();
	blockDependancies.Clear();
	computationalResults.Clear();
	
	optimalOrders.Clear();
	leafSkips.Clear();
	hasBeenSetUp			= 0;
	hasBeenOptimized		= false;
	if (computingTemplate)
	{
		delete computingTemplate;
		computingTemplate = nil;
		templateKind	  = false;
	}
	if (mstCache)
	{
		delete (mstCache);
		mstCache = nil;
	}
	
#ifdef _SLKP_LFENGINE_REWRITE_	
	treeTraversalMasks.Clear();
#ifdef	_OPENMP
	SetThreadCount		(systemCPUCount);
#endif
#endif
}


//_______________________________________________________________________________________

bool	 _LikelihoodFunction::Construct(_String& s, _VariableContainer* theP) // from triplets
//format: datasetfilter name, tree name, frequency matrix name; etc...
{	
	
	Clear ();
	
	_List Triplets(&s,';');
	
	for (int i=0; i<Triplets.lLength; i++)
	{
		_List	TripletItems (Triplets(i),',');
		if (TripletItems.lLength==3) // good triplet - continue
		{
			_String* s;
			long 	 f;
			
			// add datasetfilter
			s = (_String*)TripletItems(0);
			f = FindDataSetFilterName(AppendContainerName(*s,theP));
			if (f==-1)
			{
				WarnError (_String("\nCould not locate a data filter named: ")&*s);
				return false;			
			}
			else
				theDataFilters<<f;

			// add the tree
			
			s = (_String*)TripletItems(1);
			f = LocateVarByName(AppendContainerName(*s,theP));
			_Variable* theVar = FetchVar(f);
			if (!theVar || theVar->ObjectClass()!=TREE)
			{
				WarnError (_String("\nCould not locate a tree variable named: ")&*s);
				return false;
			}
			else
				theTrees<<variableNames.GetXtra(f);
			
			// add the matrix of probabilities
			s = (_String*)TripletItems(2);
			f = LocateVarByName(AppendContainerName(*s,theP));
			theVar = FetchVar(f);
			if ( !theVar|| theVar->ObjectClass()!=MATRIX) 
			{	
				WarnError (_String("\nCould not locate a frequencies matrix named: ")&*s);
				return false;
			}
			else
				theProbabilities<<variableNames.GetXtra(f);
				
			// at this stage also check to see whether tree tips match to species names in the dataset filter and
			// if they do - then remap
			_SimpleList remap;
			_DataSetFilter* df = ((_DataSetFilter*)dataSetFilterList(theDataFilters(i)));
			long dfDim = df->GetDimension(true),
				 efvDim = ((_Matrix*)theVar->Compute())->GetHDim();
			if ( efvDim!=dfDim)
			{	
				_String warnMsg ("The dimension of the equilibrium frequencies vector ");
				warnMsg = warnMsg & *s & " (" & efvDim & ") doesn't match the number of states in the dataset filter (" & dfDim & ')';
				&*(_String*)TripletItems(0);
				acknError (warnMsg);
				//theDataFilters.Delete(theDataFilters.lLength-1);
				//theTrees.Delete(theTrees.lLength-1);
				return false;
			}
			// first - produce the list of tip node names
			if (!MapTreeTipsToData (theProbabilities.lLength-1))	
				return false;
		}
		else
		{
			if ((i==Triplets.lLength-1)&&(TripletItems.lLength==1)&&theDataFilters.lLength) // formula spec
			{
				_String flaString (ProcessLiteralArgument((_String*)TripletItems(0),theP));
				
				if (flaString.sLength)
				{
					siteWiseVar  = CheckReceptacle (&siteWiseMatrix,empty),
					blockWiseVar = CheckReceptacle (&blockWiseMatrix,empty);
							 
					_String    flaString2 (flaString);
					_Formula   theFla (flaString,theP);
					
					
					if (theFla.IsEmpty()|| terminateExecution)
					{
						flaString = flaString2 & " is not a valid formula specification in call to LikelihoodFunction";
						WarnError (flaString);
						Clear ();
						return false;
					}
					
					bool  hasSiteMx = false,
						  hasBlkMx  = false;
						  
					_List*fList		= &theFla.GetList ();
					
					// check on proper variable dependence
					
					for (long f = 0; f<fList->lLength; f++)
					{
						_Operation * anOp = (_Operation*)(*fList)(f);
						if (anOp->IsAVariable())
						{
							long kk = anOp->GetAVariable();
							if (kk == siteWiseVar->GetAVariable())
								hasSiteMx = true;
							else
								if (kk == blockWiseVar->GetAVariable())
									hasBlkMx = true;
						}
					}
					
					templateKind = 0;
					
					if ( hasBlkMx&&hasSiteMx || (!(hasBlkMx||hasSiteMx)) )
					{
						if (!(hasBlkMx||hasSiteMx))
						{
							if (theFla.NumberOperations()==1) // hidden Markov variable
							{
								_Operation * firstOp = (_Operation*)(theFla.GetList())(0);
								if (firstOp->IsAVariable()&& !firstOp->GetANumber())
								{
									_Variable * daVar = LocateVar(firstOp->GetAVariable());
									if (daVar->IsCategory())
									{
										if (((_CategoryVariable*)daVar)->IsHiddenMarkov())
										{
											templateKind = -daVar->GetAVariable()-1;
											hasSiteMx = true;
										}
									}
								}
							}
							else
								if (theFla.NumberOperations()>=2) // user function
								{
									_Operation * lastOp = (_Operation*)(theFla.GetList())(theFla.NumberOperations()-1);
									long		 nOps = lastOp->GetNoTerms();
									if (nOps<0) // user defined function
									{
										templateKind = -nOps+1;
										hasSiteMx = true;
									}
								}
						}
						if (templateKind==0)
						{
							flaString = flaString2 & " must depend either on " & siteWiseMatrix & " or " & blockWiseMatrix & " (but not on both). Alternatively, it could be a reference to a HM category variable or a user defined BL function." ;
							WarnError (flaString);
							Clear ();
							return false;
						}
					}
					
					// verify appropriate data filter dimension for the "per-site" case
					
					_DataSetFilter*  firstFilter = (_DataSetFilter*)dataSetFilterList(theDataFilters(0));
					long			 firstFilterSize = firstFilter->GetFullLengthSpecies()/firstFilter->GetUnitLength();
					
					if (hasSiteMx)
					{
						if (templateKind==0)
							templateKind   = 1;
						
						for (long f=1; f<theDataFilters.lLength; f++)
						{
							_DataSetFilter* currentFilter     = (_DataSetFilter*)dataSetFilterList(theDataFilters(f));
							long			currentFilterSize = currentFilter->GetFullLengthSpecies()/currentFilter->GetUnitLength();
							
							if (currentFilterSize > firstFilterSize)
							{
								firstFilterSize = currentFilterSize;
								/*flaString = _String("For use of site-wise computational template all data filters in the likelihood function must have the same number of sites.")
											& *(_String*)dataSetFilterNamesList (theDataFilters(0)) & " and " & *(_String*)dataSetFilterNamesList (theDataFilters(f)) & " did not.";
								WarnError (flaString);
								Clear ();
								return false;*/
								
							}
						}
					}
					else
						templateKind = 0;
						
					// now test evaluate the formula
					
					if (templateKind==1)
					{
						_Matrix			testMx (theTrees.lLength,1,false,true);
						for (long testCount = 0; testCount<25; testCount++)
						{
							for (long di = 0; di < theTrees.lLength; di++)
								testMx.theData[di] = genrand_int32 ()/(_Parameter)RAND_MAX_32;
								
							siteWiseVar->SetValue (&testMx);
							blockWiseVar->SetValue (&testMx);
							
							_PMathObj    testResult = theFla.Compute();
							flaString = empty;
							
							if ((!testResult)|| terminateExecution)
								flaString = _String ("Failed to evaluate computation template formula (")& flaString2 & ")in LikelihoodFunction constructor.";
							else
							{
								if (testResult->ObjectClass() == NUMBER)
								{
									if (testResult->Value()>0.0)
										flaString = _String ("Computation template formula (")& flaString2 & ")in LikelihoodFunction constructor evaluated to a positive value (as a log-likelihood - must be non-positive).";
								}
								else
									flaString = _String ("Computation template formula (")& flaString2 & ")in LikelihoodFunction constructor evaluated to a non-scalar value.";
							}
							
							if (flaString.sLength)
							{
								WarnError (flaString);
								Clear ();
								return false;
							}
						}
					}
					
					computingTemplate = (_Formula*)theFla.makeDynamic();
					
					if (templateKind && templateKind<2)
					{
						bySiteResults = new _Matrix (theTrees.lLength,firstFilterSize,false,true);
						checkPointer (bySiteResults);
					}
					else
						bySiteResults = nil;
						
				}
			}
		}
	}
	
	if (theTrees.lLength>1)
	{
		_SimpleList	     checkDupTreeIDs (theTrees);
		checkDupTreeIDs.Sort ();
		for (long i=1; i<checkDupTreeIDs.lLength; i++)
			if (checkDupTreeIDs.lData[i] == checkDupTreeIDs.lData[i-1])
			{
				_String errMsg = _String("The same tree - ")&*LocateVar(checkDupTreeIDs.lData[i])->GetName() & " - can't be used in multiple partitions in the same likelihood function. You should create an independent tree for each partition, and constrain the parameters instead.";
				WarnError (errMsg);
				return false;
			}
	}
			
	
	ScanAllVariables();
	
	Setup();
	
	return true;
}

//_______________________________________________________________________________________

_LikelihoodFunction::_LikelihoodFunction (_LikelihoodFunction& lf) // stack copy
{	
	Clear();
	
	hasBeenOptimized    = lf.hasBeenOptimized;
	templateKind        = lf.templateKind;
	
	if (lf.computingTemplate)
		computingTemplate   = (_Formula*)lf.computingTemplate->makeDynamic();
	else
		computingTemplate   = nil;
		
	mstCache		= nil;
	nonConstantDep  = nil;
	
	Duplicate (&lf);
}
	
//_______________________________________________________________________________________

BaseRef _LikelihoodFunction::makeDynamic (void)  // dynamic copy of this object
{
	_LikelihoodFunction * res = new _LikelihoodFunction;
	checkPointer(res);
	memcpy ((char*)res, (char*)this, sizeof (_LikelihoodFunction));
	if (!res)
	{
		isError(0);
		return nil;
	}
	res->Duplicate(this);
	return res;
}	

//_______________________________________________________________________________________

void	_LikelihoodFunction::Duplicate (BaseRef obj) // duplicate an object into this one
{
	_LikelihoodFunction* lf = (_LikelihoodFunction*)obj;
	theTrees.Duplicate(&lf->theTrees);
	theProbabilities.Duplicate(&lf->theProbabilities);
	theDataFilters.Duplicate(&lf->theDataFilters);
	indexInd.Duplicate(&lf->indexInd);
	indexDep.Duplicate(&lf->indexDep);
	indexCat.Duplicate(&lf->indexCat);
	blockDependancies.Duplicate(&lf->blockDependancies);
	computationalResults.Duplicate(&lf->computationalResults);
	siteResults = nil;
	
	optimalOrders.Duplicate(&lf->optimalOrders);
	leafSkips.Duplicate (&lf->leafSkips);
	templateKind        = lf->templateKind;
	
	if (lf->computingTemplate)
		computingTemplate   = (_Formula*)lf->computingTemplate->makeDynamic();
	else
		computingTemplate   = nil;
		
	if (lf->mstCache)
	{
		mstCache = new MSTCache;
		checkPointer (mstCache);
		
		mstCache->computingOrder.Duplicate(&lf->mstCache->computingOrder);
		mstCache->storageOrder.Duplicate  (&lf->mstCache->storageOrder);
		mstCache->referenceOrder.Duplicate(&lf->mstCache->referenceOrder);
		mstCache->parentOrder.Duplicate (&lf->mstCache->parentOrder);
		mstCache->resultCache.Duplicate (&lf->mstCache->resultCache);
		mstCache->statesCache.Duplicate (&lf->mstCache->statesCache);
		mstCache->cacheSize.Duplicate(&lf->mstCache->cacheSize);
	}
		
	if (lf->bySiteResults)
		bySiteResults = (_Matrix*)lf->bySiteResults->makeDynamic();
	else
		bySiteResults = nil;
	
	if (lf->nonConstantDep)
		nonConstantDep = (_SimpleList*)lf->nonConstantDep->makeDynamic();
	else
		nonConstantDep = nil;
}	
	

//_______________________________________________________________________________________
_SimpleList&	_LikelihoodFunction::GetIndependentVars (void) 
{
	return indexInd;
}

//_______________________________________________________________________________________
_SimpleList&	_LikelihoodFunction::GetDependentVars (void) 
{
	return indexDep;
}

//_______________________________________________________________________________________
_SimpleList&	_LikelihoodFunction::GetCategoryVars (void) 
{
	return indexCat;
}

//_______________________________________________________________________________________
void	_LikelihoodFunction::GetGlobalVars (_SimpleList& rec) 
{
	_Variable* 		thisV;
	long			k;
	
	for (k=0; k<indexInd.lLength; k++)
	{
		thisV = LocateVar (indexInd.lData[k]);
		if (thisV->IsGlobal())
			rec << indexInd.lData[k];
	}
	for (k=0; k<indexDep.lLength; k++)
	{
		thisV = LocateVar (indexDep.lData[k]);
		if (thisV->IsGlobal())
			rec << indexDep.lData[k];
	}
}

//_______________________________________________________________________________________
_Parameter	_LikelihoodFunction::GetIthIndependent (long index) 
{
	return ((_Constant*) LocateVar (indexInd.lData[index])->Compute())->Value();
}

//_______________________________________________________________________________________
_Parameter	_LikelihoodFunction::GetIthDependent (long index) 
{
	return ((_Constant*) LocateVar (indexDep.lData[index])->Compute())->Value();
}

//_______________________________________________________________________________________
void	_LikelihoodFunction::SetIthIndependent (long index, _Parameter p) 
{
	_Variable * v =(_Variable*) LocateVar (indexInd.lData[index]);
	_Constant c(p);
	v->SetValue (&c);
}

//_______________________________________________________________________________________
bool	_LikelihoodFunction::CheckAndSetIthIndependent (long index, _Parameter p) 
{
	_Variable * v =(_Variable*) LocateVar (indexInd.lData[index]);

	bool set;
	if (p!=0.0)
		set = (fabs(v->Value()-p)/p)>machineEps;
	else
		set = fabs(v->Value()-p)>machineEps;

	if (set)
	{
		_Constant c(p);
		v->SetValue (&c);
	}
	
	return set;
}

//_______________________________________________________________________________________
void	_LikelihoodFunction::SetIthDependent (long index, _Parameter p) 
{
	_Variable * v =(_Variable*) LocateVar (indexDep.lData[index]);
	_Constant c(p);
	v->SetValue (&c);
}
	
//_______________________________________________________________________________________
_Matrix*	_LikelihoodFunction::RemapMatrix(_Matrix* source)
{
	_DataSetFilter* df;
	long hDim = source->GetHDim(), vDim = 0,i,offset=0,bigOffset=0;
	for (i=0; i<theTrees.lLength;i++)
	{
		df = ((_DataSetFilter*)dataSetFilterList(theDataFilters.lData[i]));
		vDim+=df->GetFullLengthSpecies()/df->GetUnitLength();
	}
	_Matrix* res = new _Matrix (hDim,vDim,false,true);
	checkPointer(res);
	for (i=0; i<theTrees.lLength;i++)
	{
		df = ((_DataSetFilter*)dataSetFilterList(theDataFilters(i)));
		long j,k,dl = df->GetFullLengthSpecies()/df->GetUnitLength();
		if (HasHiddenMarkov(blockDependancies.lData[i])>=0)
		{
			for (j=bigOffset+dl-1;j>=bigOffset;j--)
			{
				for (long l=0; l<hDim; l++)
					res->Store(l,j,(*source)(l,j+offset));
			}
			offset+=dl;
			bigOffset+=dl;
		}
		else
		{
			for (j=bigOffset+dl-1;j>=bigOffset;j--)
			{
				k = df->duplicateMap.lData[j-bigOffset];
				for (long l=0; l<hDim; l++)
				{
					_Parameter value = (*source)(l,k+offset);
					if (value)
						res->Store(l,j,value);
				}
			}
			offset+=BlockLength(i);
			bigOffset+=dl;
		}
	}
	res->AmISparse();
	return res;

}	

//_______________________________________________________________________________________

#ifdef __HYPHYMPI__
void		_LikelihoodFunction::MPI_LF_Compute (long senderID, bool partMode)
#else
void		_LikelihoodFunction::MPI_LF_Compute (long, bool)
#endif
{
	#ifdef __HYPHYMPI__
		long i, 
			 hDim = 0;
			 
		for (i=0;i<theTrees.lLength;i++)
			hDim += BlockLength (i);
		
		if (!partMode)
		{
			siteResults = new _Matrix   (1,hDim,false,true);	
			checkPointer(siteResults);			
		}
		
		_Matrix		  variableStash (indexInd.lLength,1,false,true);
		_Parameter	  siteLL = 0.0;
		
		MPI_Status	  status;
		//ReportWarning (_String ("Waiting on:") & (long) indexInd.lLength & " MPI_DOUBLES");
		ReportMPIError(MPI_Recv(variableStash.theData, indexInd.lLength, MPI_DOUBLE, senderID, HYPHY_MPI_VARS_TAG, MPI_COMM_WORLD,&status),false);
	
		long  loopie = 0;
		while (variableStash.theData[0]>=-1e100)
		{
			//ReportWarning (_String("In at step  ") & loopie);
			bool	doSomething = false;
			for (i = 0; i<indexInd.lLength; i++)
			{
				_Variable *anInd = LocateVar(indexInd.lData[i]);
				//ReportWarning (*anInd->GetName() & " = " & variableStash.theData[i]);
				if (!CheckEqual(anInd->Value(), variableStash.theData[i]))
				{
					doSomething = true;
					SetIthIndependent (i,variableStash.theData[i]);
				}
			}
			if (doSomething)
			{
				if (partMode)
					siteLL = Compute();
				else
				{
					categID = -1;
					long    shifter = 0;
					
					for (i=0; i<theTrees.lLength; i++)
					{
						ComputeBlock (i,siteResults->theData+shifter);
						shifter += BlockLength (i);
					}
						
					for (i=0;i<indexDep.lLength; i++)
						LocateVar (indexDep.lData[i])->Compute();
						
					for (i=0;i<indexInd.lLength; i++)
						LocateVar (indexInd.lData[i])->MarkDone();					
				}
		 	}	
					
			if (partMode)
				ReportMPIError(MPI_Send(&siteLL, 1, MPI_DOUBLE, senderID, HYPHY_MPI_DATA_TAG, MPI_COMM_WORLD),true);				
			else
				ReportMPIError(MPI_Send(siteResults->theData, hDim, MPI_DOUBLE, senderID, HYPHY_MPI_DATA_TAG, MPI_COMM_WORLD),true);
			//ReportWarning (_String("Sent at step  ") & loopie++ & " " & siteLL);
			//ReportWarning (_String("Sent Back... step  ") & loopie++);
			ReportMPIError(MPI_Recv(variableStash.theData, indexInd.lLength, MPI_DOUBLE, senderID, HYPHY_MPI_VARS_TAG, MPI_COMM_WORLD,&status),false);
		}
		//ReportWarning (_String("Exiting slave loop after step  ") & loopie);
		if (!partMode)
		{
			DeleteObject (siteResults);
			siteResults = nil;			
		}
	#endif
}
									  
//_______________________________________________________________________________________
_Matrix*	_LikelihoodFunction::ConstructCategoryMatrix (bool complete, bool remap) 
// 
{
	long		 				i,
				 				hDim,
				 				vDim,
				 				currentOffset = 0;
				 				
	_DataSetFilter			    *df;
	
	PrepareToCompute();
	if (!siteResults)
		AllocateSiteResults ();
	
	if (complete) 
	// just return the maximum conditional likelihood category
	{
		vDim = 0;
		hDim = 1;
		// compute the horizontal dimension of the matrix 
		for (i=0;i<theTrees.lLength;i++)
			if (HasHiddenMarkov(blockDependancies.lData[i])>=0)
			{
				df = (_DataSetFilter*)dataSetFilterList(theDataFilters.lData[i]);
				vDim+=df->GetFullLengthSpecies()/df->GetUnitLength();
			}
			else
				vDim+=BlockLength(i);
			
		_Matrix *result = new _Matrix (hDim,vDim,false,true);
		checkPointer (result);
		// now proceed to compute each of the blocks
		for (i=0;i<theTrees.lLength;i++)
		{
			categID = 0;
			offsetCounter = 1;
			if (!blockDependancies.lData[i])
			{
				long u = currentOffset+BlockLength(i);
				for (long j=currentOffset; j<u; j++)
					(*result)[j]=0.0;
				currentOffset+=BlockLength(i);
			}
			else
			{
				long	hmmID = HasHiddenMarkov (blockDependancies.lData[i]);
				if (hmmID>=0)
				{
						// do hidden markov
						df = (_DataSetFilter*)dataSetFilterList(theDataFilters.lData[i]);
						_CategoryVariable* thisC = (_CategoryVariable*)LocateVar(indexCat.lData[hmmID]);
						thisC->Refresh();
						
						long      		   ni = thisC->GetNumberOfIntervals(),
										   bl = BlockLength (i),
										   k,
										   m,
										   mi,
										   mi2,
										   ul = df->GetUnitLength(),
										   tl = df->theOriginalOrder.lLength,
										   ml = tl/ul,
										   j,
										   hmmShifter = 1;
										   
						_Matrix	  		   hmc (ni,bl,false,true),
										  *hmm = thisC->ComputeHiddenMarkov(),
										   temp(ni,ml,false,true),
										  *hmf = thisC->ComputeHiddenMarkovFreqs();
										  
						_Parameter		  *p1 = hmc.fastIndex(),
										  *p2,
										  scrap,
										  max,
										  min;
								
						for (m=HighestBit(blockDependancies.lData[i]);m>=0;m--)
							if (CheckNthBit (blockDependancies.lData[i],m))
							{
								_CategoryVariable * aC = (_CategoryVariable*)LocateVar(indexCat.lData[i]);
								if (! (aC->IsHiddenMarkov() || aC->IsConstantOnPartition()))
									hmmShifter *= aC->GetNumberOfIntervals();
							}
						
						offsetCounter = 1;
						for (m=0; m<ni; m++,categID++)
						{
							categID = m*hmmShifter;
							thisC->SetIntervalValue (m);
							RecurseCategory (i,0,blockDependancies.lData[i],HighestBit(blockDependancies.lData[i]),1);
							p2 = siteResults->fastIndex();
							for   (k = 0; k<bl; k++,p2++,p1++)
								*p1 = *p2;
						}
						
						categID = 0; // correction factor
						p1 = hmc.fastIndex();
						mi = df->duplicateMap.lData[tl/ul-1];
						for (m=ml-1; m<ml*ni; m+=ml,mi += bl)
							temp.theData[m] = p1[mi];
						
						for (j=ml-2; j>=0; j--)
						{
							for (k=0; k<ni;k++)
							{
								mi  = df->duplicateMap.lData[j];
								mi2 = j+1;
								max = 0.0;
								min = 1.e300;
								for (m=0; m<ni; m++,mi += bl,mi2 += ml)
								{
									scrap = hmm->theData[k*ni+m] * p1[mi] * temp.theData[mi2];
									if (scrap > max)
										max = scrap;
									if (scrap < min)
										min = scrap;
									
								}
								temp.theData[j+k*ml] = max;
							}
							k = log (min)-1;
							//categID -= k;
							scrap = exp ((_Parameter)(-k));
							mi2 = j;
							for (; mi2<ni*ml; mi2 += ml)
								temp.theData[mi2] *= scrap;
							
						}
						max = 0.;
						mi2 = 0;
						for (k=0,mi=0; k<ni; k++,mi+=ml)
						{
							scrap = temp.theData[mi] * hmf->theData[k];
							if (scrap > max)
							{
								mi2 = mi;
								max = scrap;
							}
						}
						mi2 /= ml;
						result->theData[currentOffset] = mi2;
						for (k=1; k<ml; k++)
						{
							mi = k;
							max = 0.;
							long mi3  = df->duplicateMap.lData[k],
								 mi4  = 0;
							for (m=0; m<ni; m++, mi += ml, mi3 += bl)
							{
								scrap = hmm->theData[mi2*ni+m] * p1[mi3] * temp.theData[mi];
								if (scrap > max)
								{
									max = scrap;
									mi4 = m;
								}
							}
							mi2 = mi4;
							result->theData[currentOffset+k] = mi2;
						}
						currentOffset+=ml;
					}
					else
					{
						long ff = HasHiddenMarkov(blockDependancies.lData[i],false);
						if (ff<0)
							FindMaxCategory(i,  HighestBit( blockDependancies.lData[i]), blockDependancies.lData[i], LowestBit( blockDependancies.lData[i]), currentOffset, *result);						
						else
						{
							long	mxDim = 1,
								    bl = BlockLength (i),
								    j;

							for (j=0; j<=ff; j++)
								if (CheckNthBit (blockDependancies.lData[i],j))
								{
									_CategoryVariable* thisC = (_CategoryVariable*)LocateVar(indexCat.lData[j]);
									mxDim *= thisC->GetNumberOfIntervals();
								}
							
							_Matrix ccache (mxDim,1,false,true);
							
							categID = 0;
							RecurseConstantOnPartition    (i,0,blockDependancies.lData[j],ff,1,ccache);
							categID = 0;
							
							/*
								FindMaxCategory(i,  HighestBit( blockDependancies.lData[i]), blockDependancies.lData[i], ff+1, currentOffset, result);						
							*/
							
							_Parameter maxValue = -1.e300;
							long	   mxIndex  = -1;
							
							for (j=0; j<mxDim; j++)
								if (ccache.theData[j] > maxValue)
								{
									maxValue = ccache.theData[j];
									mxIndex = j;
								}
								
							mxDim = 1;
							for (j=HighestBit (blockDependancies.lData[i]); j>ff; j--)
								if (CheckNthBit (blockDependancies.lData[i],j))
								{
									_CategoryVariable* thisC = (_CategoryVariable*)LocateVar(indexCat.lData[j]);
									mxDim *= thisC->GetNumberOfIntervals();
								}
							
							mxIndex *= mxDim;
							
							for   (long kk = 0; kk<bl; kk++)
								result->theData[currentOffset+kk] += mxIndex;							
							
						}					
						currentOffset+=BlockLength(i);
					}
			}
		}
		DoneComputing();
		if (remap)
		{
			_Matrix * retObj = RemapMatrix(result);
			DeleteObject (result);
			return retObj;
		}
		return result;
	}
	else
	{
		vDim = 0;
		hDim = 1;
		for (i=0;i<theTrees.lLength;i++)
		{
			vDim+=BlockLength(i);
			long dim = 1, checkByte = blockDependancies.lData[i],m;
			for (m=sizeof(long)*4-1;m>=0;m--)
				if (CheckNthBit(checkByte,m))
					dim *= ((_CategoryVariable*)LocateVar(indexCat.lData[m]))-> GetNumberOfIntervals() ;
			if (dim>hDim) hDim = dim;
		}
		//decide what the vertical dimension should be
	
		_Matrix result (hDim,vDim,false,true);
		// now proceed to compute each of the blocks
		for (i=0;i<theTrees.lLength;i++)
		{
			categID = 0;
			offsetCounter = 1;
			if (!blockDependancies.lData[i])
			{
				categID = -1;
				forceRecomputation = true;
				ComputeBlock(i,result.fastIndex()+currentOffset);
				forceRecomputation = false;
			}
			else
				WriteAllCategories(i,  HighestBit( blockDependancies.lData[i]), blockDependancies.lData[i], LowestBit( blockDependancies.lData[i]), currentOffset, result);
			currentOffset+=BlockLength(i);
		}
		DoneComputing();
		//return (_Matrix*)result.makeDynamic();
		if (remap)
			return RemapMatrix(&result);
		return (_Matrix*)result.makeDynamic();
	}
	return nil;
}

//_______________________________________________________________________________________
void	_LikelihoodFunction::AllocateSiteResults 		(void)
{
	long i, 
		 hDim = 0;
		 
	for (i=0;i<theTrees.lLength;i++)
	{
		//if ((BlockLength(i)>hDim)&& blockDependancies (i))	/* mod 05/01/2003 */
		long bl = BlockLength (i);
		if (bl > hDim)
			hDim = bl;
	}
	siteResults = new _Matrix (hDim,2,false,true);
	checkPointer(siteResults);
}

//_______________________________________________________________________________________

#ifndef __HYPHYMPI__
void	  _LikelihoodFunction::SendOffToMPI		  (long)
{
#else
void	  _LikelihoodFunction::SendOffToMPI		  (long index)
{
	long 			vi  = 0;
	_SimpleList*	vim = (_SimpleList*)parallelOptimizerTasks (index);
	
	for (long k=0; k<indexCat.lLength; k++,vi++)
		varTransferMatrix.theData[vim->lData[vi]] = LocateVar(indexCat.lData[k])->Value();
	for (long k=0; k<indexInd.lLength; k++,vi++)
	{
		long idx = vim->lData[vi];
		if (idx>=0)
			varTransferMatrix.theData[idx] = LocateVar(indexInd.lData[k])->Value();
	}
	ReportMPIError(MPI_Send(varTransferMatrix.theData, transferrableVars, MPI_DOUBLE, index+1 , HYPHY_MPI_VARS_TAG, MPI_COMM_WORLD),true);
#endif //__HYPHYMPI__
}

//_______________________________________________________________________________________
		
#ifdef __HYPHYMPI__
void	  _LikelihoodFunction::RecurseCategoryMPI (long index, long lowestIndex, _Parameter weight)
#else
void	  _LikelihoodFunction::RecurseCategoryMPI (long, long, _Parameter)
#endif
{
	 #ifdef __HYPHYMPI__
	_CategoryVariable* thisC = (_CategoryVariable*)LocateVar(indexCat.lData[index]);
	if (index>lowestIndex)
	{
		thisC->Refresh();
		long nI = thisC->GetNumberOfIntervals ();
		offsetCounter *= nI;
		for (long k = 0; k<nI; k++)
		{
			thisC->SetIntervalValue(k);
			RecurseCategoryMPI (index-1,lowestIndex,weight*thisC->GetIntervalWeight(k));
			categID+=offsetCounter/nI;
		}
		offsetCounter/=nI;
		if (offsetCounter>1)
			categID-=nI*offsetCounter;
	}
	else
	{
		//ReportWarning (_String("RecurseCategoryMPI: In ") & categID);
		long 		nI = thisC->GetNumberOfIntervals ();
							
		thisC->Refresh();

		_Matrix*	cws= thisC->GetWeights();
		for (long k = 0; k<nI; k++)
		{
			thisC->SetIntervalValue(k,!k);
			catFreqMatrix.theData[categID] = cws->theData[k]*weight ;
			SendOffToMPI (categID);
			categID+=offsetCounter;
		}
		if (offsetCounter>1)
			categID-=nI*offsetCounter;
		//ReportWarning (_String("RecurseCategoryMPI: Out ") & categID);
	}
	#endif
}

//_______________________________________________________________________________________
		
_Parameter	_LikelihoodFunction::ComputeMasterMPI (void)
{
	_Parameter res = 0.0;
	 
	 #ifdef __HYPHYMPI__
		
		if (mpiPartitionOptimizer)
		{
			//printf ("In ComputeMasterMPI\n");
			long	totalSent = 0;
			
			for (long blockID = 0; blockID < parallelOptimizerTasks.lLength; blockID ++)
			{
				bool sendToSlave = (computationalResults.lLength < parallelOptimizerTasks.lLength);
				_SimpleList * slaveParams = (_SimpleList*)parallelOptimizerTasks(blockID);
				
				for (long varID = 0; varID < slaveParams->lLength; varID++)
				{
					_Variable * aVar = LocateVar (slaveParams->lData[varID]);
					varTransferMatrix.theData[varID] = aVar->Value();
					sendToSlave = sendToSlave || aVar->HasChanged();
				}
				//printf ("Block %d with %d parameters (send flag %d)\n", blockID, slaveParams->lLength, sendToSlave);
				if (sendToSlave)
				{
					ReportMPIError(MPI_Send(varTransferMatrix.theData, slaveParams->lLength, MPI_DOUBLE, blockID+1 , HYPHY_MPI_VARS_TAG, MPI_COMM_WORLD),true);
					totalSent++;
				}
				else
					res += ((_Constant*)computationalResults(blockID))->Value();
			}
			
			if (computationalResults.lLength < parallelOptimizerTasks.lLength)
			{
				computationalResults.Clear();
				for (long k=0; k<parallelOptimizerTasks.lLength;k=k+1)
				{
					_Constant * c = new _Constant (0.0);
					checkPointer (c);
					computationalResults << c;
					DeleteObject (c); 
				}
			}
					
		 	while (totalSent)
		 	{
		 		long 		index = 0;
		 		MPI_Status	status;
		 		_Parameter  blockRes;
				ReportMPIError(MPI_Recv (&blockRes, 1, MPI_DOUBLE, MPI_ANY_SOURCE , HYPHY_MPI_DATA_TAG, MPI_COMM_WORLD,&status),true);
				
				index = status.MPI_SOURCE;
				//printf ("Got %g from block %d \n", blockRes, index);

				res += blockRes;
				((_Constant*)computationalResults(index-1))->SetValue (blockRes);
				totalSent--;
		 	}
		}
		else
		{
			categID = 0;
		 	RecurseCategoryMPI (indexCat.lLength-1,0,1.0);
		 	
		 	for (long k=resTransferMatrix.GetVDim(); k<2*resTransferMatrix.GetVDim(); k++)
		 		resTransferMatrix.theData[k] = 0.0;

			if (computingTemplate&& templateKind>1)
			{
				_Matrix ** receptacles = new _Matrix* [theTrees.lLength];
				checkPointer (receptacles);
				
				for (long i=0; i<theTrees.lLength; i++)
				{
					_DataSetFilter	*df 			= ((_DataSetFilter*)dataSetFilterList(theDataFilters.lData[i]));
					long 			blockWidth 		= df->GetFullLengthSpecies()/df->GetUnitLength();
						
					_Matrix 		*blockResult 	= nil;
					_String 		mxName 			= siteWiseMatrix & (long)(i+1); 
					_Variable* 		matrixStorage 	= nil;
						
					categID = LocateVarByName (mxName);
						
					if (categID<0)
					{
						_Variable dummy (mxName);
						matrixStorage = LocateVar (dummy.GetAVariable());
					}
					else
					{
						matrixStorage = FetchVar (categID);
						if (matrixStorage->ObjectClass() == MATRIX)
							blockResult = (_Matrix*)matrixStorage->GetValue();
					}
					
					bool			madeMatrix 		= (!blockResult);	
					long			nCat = ((_TheTree*)LocateVar(theTrees.lData[i]))->categoryCount;

					if (!madeMatrix)
					{
						if ((blockResult->GetHDim() != blockWidth) || (blockResult->GetVDim() != nCat))
							madeMatrix = true;
						else
							blockResult->CheckIfSparseEnough(true);
					}

					if (madeMatrix)
					{
						blockResult 		= new _Matrix (nCat,blockWidth,false,true);
						checkPointer	(blockResult);
						matrixStorage->SetValue (blockResult);
						DeleteObject  	(blockResult);
						receptacles [i] = (_Matrix*)matrixStorage->GetValue();
					}
					else
						receptacles [i] = blockResult;

				}	
					
			 	for (long k=0; k<MPICategoryCount; k++)
			 	{

			 		long index = 0;
			 		MPI_Status	status;
					ReportMPIError(MPI_Recv (resTransferMatrix.theData, resTransferMatrix.GetVDim(), MPI_DOUBLE, MPI_ANY_SOURCE , HYPHY_MPI_DATA_TAG, MPI_COMM_WORLD,&status),true);
					
					index = status.MPI_SOURCE-1;

					_Parameter *p1 	  = resTransferMatrix.theData;
					
					for (long i=0; i<theTrees.lLength; i++)
					{
						_DataSetFilter* df = ((_DataSetFilter*)dataSetFilterList(theDataFilters(i)));
						
						long *dupInfo	  = df->duplicateMap.lData;
						
						long upTo 		  = receptacles [i]->GetVDim();
						_Parameter *p2	  = receptacles [i]->theData+upTo*index;
						
						for (long jj = 0; jj < upTo; jj++)
							*(p2++) = p1[dupInfo[jj]];
							
						p1 += BlockLength(i);
							
					}
			 	}
			 	
				delete receptacles;
				if (terminateExecution)
					return -1e100;
				res = computingTemplate->Compute()->Value();
			}
			else
			{
			 	for (long k=0; k<MPICategoryCount; k++)
			 	{
			 		long index = 0;
			 		MPI_Status	status;
					ReportMPIError(MPI_Recv (resTransferMatrix.theData, resTransferMatrix.GetVDim(), MPI_DOUBLE, MPI_ANY_SOURCE , HYPHY_MPI_DATA_TAG, MPI_COMM_WORLD,&status),true);
					
					index = status.MPI_SOURCE;

					_Parameter weight = catFreqMatrix.theData[index-1],
							  *p1 	  = resTransferMatrix.theData,
							  *p2	  = p1+resTransferMatrix.GetVDim();
					
			 		for (long kk=resTransferMatrix.GetVDim()-1; kk>=0; kk--)
			 			*(p2++) += weight**(p1++);
			 	}
			 			
			 	_Parameter	  *sR;
				sR = resTransferMatrix.theData+(2*resTransferMatrix.GetVDim()-1);
				
				for (long j=0;j<theTrees.lLength;j++)
				{
					_DataSetFilter* df = ((_DataSetFilter*)dataSetFilterList(theDataFilters(j)));
					for (long i=BlockLength(j)-1;i>=0;i--,sR--)
						res+=myLog(*sR)*df->theFrequencies.lData[i];

				}
			}			
		}

		blockComputed = -1;
		likeFuncEvalCallCount++;
	 #endif
	 
	 return res;
}

//_______________________________________________________________________________________
		
bool	_LikelihoodFunction::PreCompute 		(void)
{
	
	//useGlobalUpdateFlag = true; 
	// mod 20060125 to only update large globals once
	long i = 0;
	
	_SimpleList * arrayToCheck = nonConstantDep?nonConstantDep:&indexDep;
	
	for (;i < arrayToCheck->lLength;i++)
	{
		_Variable* cornholio = LocateVar(arrayToCheck->lData[i]);
		_Parameter result 	 = ((_Constant*) cornholio->Compute())->Value();
		
		if (result<cornholio->GetLowerBound() || result>cornholio->GetUpperBound())
			break;
	}

	//useGlobalUpdateFlag = false; 
	// mod 20060125 to only update large globals once

	for (long j=0;j<arrayToCheck->lLength;j++)
	{
		_Variable* cornholio = LocateVar(arrayToCheck->lData[j]);
		if (cornholio->varFlags&HY_DEP_V_COMPUTED)
			cornholio->varFlags -= HY_DEP_V_COMPUTED;
	}

	return (i==arrayToCheck->lLength);
}


//_______________________________________________________________________________________
		
void	_LikelihoodFunction::PostCompute 		(void)
{
	_SimpleList * arrayToCheck = nonConstantDep?nonConstantDep:&indexDep;

	for (long i=0;i<arrayToCheck->lLength; i++)
		//LocateVar (indexDep.lData[i])->PostMarkChanged();
		LocateVar (arrayToCheck->lData[i])->Compute();
	// mod 20060125 comment out the compute loop; seems redundant	
	{
	for (long i=0;i<indexInd.lLength; i++)
		LocateVar (indexInd.lData[i])->MarkDone();
	}
}



//_______________________________________________________________________________________
		
void	_LikelihoodFunction::ComputeBlockForTemplate 		(long i, bool f)
{
	long 	blockWidth = bySiteResults->GetVDim();
	_Parameter	 * resStore = bySiteResults->theData+i*blockWidth;
#ifdef _SLKP_LFENGINE_REWRITE_
	categID = -1;
#endif
	ComputeBlock (i,resStore);
#ifdef _SLKP_LFENGINE_REWRITE_
	// multiply scaling factors into the expression
	if (! usedCachedResults)
	{
		long	*   siteCorrectors	= ((_SimpleList**)siteCorrections.lData)[i]->lData,
					upto			= ((_SimpleList**)siteCorrections.lData)[i]->lLength;
		for (long s = 0; s < upto; s++)
		{
		//	printf ("%d\t%d\t%g\t%d\n", i, s, resStore[s], siteCorrectors[s]);
			resStore[s] *= acquireScalerMultiplier(siteCorrectors[s]);
		}
	}
#endif
	if (f || !usedCachedResults)
	// remap compressed sites to full length
		ComputeBlockForTemplate2 (i, resStore, resStore, blockWidth);
}


//_______________________________________________________________________________________
		
void	_LikelihoodFunction::ComputeBlockForTemplate2 		(long i, _Parameter* resTo, _Parameter* resFrom, long blockWidth)
{
	_DataSetFilter		*df = ((_DataSetFilter*)dataSetFilterList(theDataFilters.lData[i]));
	long*				dupMap = df->duplicateMap.lData,
						dupL   = df->duplicateMap.lLength;
	
	if (resTo == resFrom)
	{
		_Matrix 			temp (1,blockWidth,false,true);
		for (long k1 = 0; k1 < dupL; k1++)
			temp.theData[k1] = resFrom[dupMap[k1]];
			
		for (long k2 = 0; k2< dupL; k2++)
			resTo[k2] = temp.theData[k2];	
			
		for (long k3 = dupL; k3 < blockWidth; k3++)
			resTo[k3] = 1.;
	}
	else
	{
		for (long k1 = 0; k1 < dupL; k1++)
			resTo[k1] = resFrom[dupMap[k1]];	
		
		for (long k3 = dupL; k3 < blockWidth; k3++)
			resTo[k3] = 1.;
	}
}

//_______________________________________________________________________________________
		
_Parameter	_LikelihoodFunction::Compute 		(void)
{
		
	_Parameter result = 0.;
	long	   i;
	
	if (!PreCompute())
		return -A_LARGE_NUMBER;
	
	if (!isInOptimize && hasBeenOptimized)
		for (i=0;i<indexInd.lLength; i++)
			if (LocateVar (indexInd.lData[i])->HasChanged())
			{
				hasBeenOptimized = false;
				break;
			}
			
	#ifdef __HYPHYMPI__
		if (parallelOptimizerTasks.lLength)
		{
			return ComputeMasterMPI ();
		}
		else
		{
	#endif
	
	if (computingTemplate&& templateKind>1) // spawn a matrix of partial likelihoods (one for each partition - with appropriate dimensions)
	{
		if (!siteResults)
			AllocateSiteResults ();
			
		for (i=0; i<theTrees.lLength; i++)
		{
			if (HasBlockChanged(i))
			{
				_DataSetFilter	*df 			= ((_DataSetFilter*)dataSetFilterList(theDataFilters.lData[i]));
				long 			blockWidth 		= df->GetFullLengthSpecies()/df->GetUnitLength();
				
				_Matrix 		*blockResult 	= nil;
				_String 		mxName 			= siteWiseMatrix & (long)(i+1); 
				_Variable* 		matrixStorage 	= CheckReceptacle (&mxName, empty, false);
				
				if (matrixStorage->ObjectClass() == MATRIX)
					blockResult = (_Matrix*)matrixStorage->GetValue();
				
				bool			madeMatrix 		= !blockResult;
				
				
				if (blockDependancies.lData[i])
				{
					categID 		= 0;
					offsetCounter	= 1;
					long			nCat = ((_TheTree*)LocateVar(theTrees.lData[i]))->categoryCount;
																			
					if (!madeMatrix)
					{
						if ((blockResult->GetHDim() != blockWidth) || (blockResult->GetVDim() != nCat))
							madeMatrix = true;
						else
							blockResult->CheckIfSparseEnough(true);
					}

					if (madeMatrix)
					{
						blockResult 		= new _Matrix (nCat,blockWidth,false,true);
						checkPointer	(blockResult);
					}
					
					WriteAllCategoriesTemplate(i,  HighestBit( blockDependancies.lData[i]), blockDependancies.lData[i], LowestBit( blockDependancies.lData[i]), blockWidth,blockResult->theData,df->duplicateMap.lData);
				}
				else
				{
					categID = 0; /* mod 07/23/2004 - was = -1 */
					ComputeBlock (i,siteResults->theData);
					if (!usedCachedResults)
					// remap compressed sites to full length
					{
						if (!madeMatrix)
						{
							if ((blockResult->GetHDim() != blockWidth) || (blockResult->GetVDim() != 1))
								madeMatrix = true;
							else
								blockResult->CheckIfSparseEnough(true);
						}
						
						if (madeMatrix)
							checkPointer (blockResult 		= new _Matrix (1,blockWidth,false,true));
						
						long*				dupMap = df->duplicateMap.lData;
						
						for (long k1 = 0; k1 < blockWidth; k1++)
							blockResult->theData[k1] = siteResults->theData[dupMap[k1]];
					}
				}
				
				if (madeMatrix)
				{
					matrixStorage->SetValue (blockResult);
					DeleteObject (blockResult);
				}
			}
		}
		
		blockComputed = -1;
		likeFuncEvalCallCount++;
		if (terminateExecution)
			return -1e100;
		result = computingTemplate->Compute()->Value(); // mod 20060125 added the PostCompute call here
		PostCompute();
		return result;
	}
	else
	{
		if (indexCat.lLength==0)
		{
			if (computingTemplate && templateKind)
			{
				categID = 0; /* mod 07/23/2004 - was = -1 */
				for (i=0; i<theTrees.lLength; i++)
					ComputeBlockForTemplate (i);
			}
			else
				for (i=0; i<theTrees.lLength; i++)
					result+=ComputeBlock (i);				
		}
		else
		{
			long hDim = 0;				 
				 
			if (!siteResults)
				AllocateSiteResults ();
			
			hDim = siteResults->GetHDim();
			
			_Parameter *sR = siteResults->fastIndex();

			for (long j=0;j<theTrees.lLength;j++)
			{
				for (i=0;i<2*hDim;i++)
					sR[i] = 0.;
				
				categID = 0;
				
				long	blockLength = BlockLength(j);
				
				if (blockDependancies.lData[j])
				{
					if ( computationalResults.lLength<=j || HasBlockChanged(j))
					{
						offsetCounter = HasHiddenMarkov(blockDependancies.lData[j]);
						_DataSetFilter* df = ((_DataSetFilter*)dataSetFilterList(theDataFilters(j)));
						if (offsetCounter<0)
						{
							offsetCounter = 1;
							_Parameter blockResult = 0.0;
							
							long ff = HasHiddenMarkov(blockDependancies.lData[j],false);
							if (ff<0)
							{
#ifdef _SLKP_LFENGINE_REWRITE_
								_SimpleList	scalerTabs (blockLength, 0, 0);
								long		cumulativeCorrection = 0;
#endif
								RecurseCategory (j,0,blockDependancies.lData[j],HighestBit(blockDependancies.lData[j]),1
#ifdef _SLKP_LFENGINE_REWRITE_
								,&scalerTabs
#endif
								);
								
								if (computingTemplate && templateKind)
									ComputeBlockForTemplate2 (j, bySiteResults->theData+j*bySiteResults->GetVDim(), sR, bySiteResults->GetVDim());
								else
									for (long s = 0; s < blockLength; s++)
									{
										blockResult+=myLog(sR[s])*df->theFrequencies.lData[s];
#ifdef _SLKP_LFENGINE_REWRITE_
										cumulativeCorrection += scalerTabs.lData[s];
#endif
									}
#ifdef _SLKP_LFENGINE_REWRITE_
								blockResult -= cumulativeCorrection*_logLFScaler;
								//printf ("%d,%g\n", likeFuncEvalCallCount, blockResult);
#endif
								
							}
							else
							// do constant on partition
							{
								long	mxDim = 1;
								for (i=0; i<=ff; i++)
									if (CheckNthBit (blockDependancies.lData[j],i))
									{
										_CategoryVariable* thisC = (_CategoryVariable*)LocateVar(indexCat.lData[i]);
										mxDim *= thisC->GetNumberOfIntervals();
									}
								
								_Matrix ccache (mxDim,1,false,true);
								
								RecurseConstantOnPartition (j,0,blockDependancies.lData[j],ff,1,ccache);
								
								_Parameter maxValue = -1.e300;
																	
								for (i=0; i<mxDim; i++)
									if (ccache.theData[i]>maxValue)
										maxValue = ccache.theData[i];

								for (i=0; i<mxDim; i++)
									blockResult += exp (ccache.theData[i]-maxValue);
									
								blockResult = myLog (blockResult)+maxValue;
							}
							
							if (computationalResults.lLength>j)
								((_Constant*)computationalResults(j))->SetValue (blockResult);
							else
								computationalResults.AppendNewInstance(new _Constant(blockResult));

							result  += blockResult;
							categID  = 0;
						}
						else
						{
							// do hidden markov
							_CategoryVariable* thisC = (_CategoryVariable*)LocateVar(indexCat.lData[offsetCounter]);
							thisC->Refresh();
							
							long      		   ni 				= thisC->GetNumberOfIntervals(),
											   bl 				= BlockLength (j),
											   hmmShifter		= 1;
											   
							_Matrix	  		   hmc (ni,bl,false,true),
											  *hmm = thisC->ComputeHiddenMarkov(),
											  *hmf = thisC->ComputeHiddenMarkovFreqs();
											  
							_Parameter		   *p1 = hmc.fastIndex(),
											   *p2;
									
							for (i=HighestBit(blockDependancies.lData[j]);i>=0;i--)
								if (CheckNthBit (blockDependancies.lData[j],i))
								{
									_CategoryVariable * aC = (_CategoryVariable*)LocateVar(indexCat.lData[i]);
									if (! (aC->IsHiddenMarkov() || aC->IsConstantOnPartition()))
										hmmShifter *= aC->GetNumberOfIntervals();
								}
									
							offsetCounter = 1;
							categID 	  = 0;
							// SLKP mod: 20070301; fixed category indexing bug when HMM AND 
								// another category is used
							for (i=0; i<ni; i++)
							{
								categID		= hmmShifter * i;
								thisC->SetIntervalValue (i,!i);
								RecurseCategory (j,0,blockDependancies.lData[j],HighestBit(blockDependancies.lData[j]),1);
								p2 = siteResults->fastIndex();
								for   (long k = 0; k<bl; k++,p2++,p1++)
									*p1 = *p2;
							}
							
							result += SumUpHiddenMarkov (hmc, *hmm, *hmf, df->duplicateMap, bl);
						}
					}
					else
						result+=((_Constant*)computationalResults(j))->Value();
				}
				else
				{ 
					if (computingTemplate && templateKind == 1)
					// 20060808: need to fill out site-by-site probabilities 
						//categID = 0;
						ComputeBlockForTemplate (j);
					else
						result+=ComputeBlock(j);
				}
					
			}
		}
	}
	#ifdef __HYPHYMPI__
	}
	#endif

	PostCompute ();
	blockComputed = -1;
	likeFuncEvalCallCount++;
	
	if (terminateExecution)
		return -1e100;

	
	if (computingTemplate)
	{
		siteArrayPopulated = true;
		_Matrix				  blockValues (theTrees.lLength,1,false,true);
		siteWiseVar->SetValue (&blockValues);
		
		_Matrix* varMxValue = (_Matrix*)siteWiseVar->GetValue();
		
		if (templateKind>0)
		{
			long 	blockWidth = bySiteResults->GetVDim();
			result = 0.0;
			for (i=0; i<blockWidth; i++)
			{
				for (long k=0; k<theTrees.lLength; k++)
					varMxValue->theData[k] = bySiteResults->theData[k*blockWidth+i];
					
				result += computingTemplate->Compute()->Value();
			}
			return result;
		}
		else
		{
			if (templateKind<0)
			{
				_CategoryVariable * thisC = (_CategoryVariable*)LocateVar (-templateKind-1);							
								   
				_Matrix	  		  *hmm = thisC->ComputeHiddenMarkov(),
								  *hmf = thisC->ComputeHiddenMarkovFreqs();
								  
				_SimpleList		  rm;
				
				for (i=0; i<bySiteResults->GetVDim(); i++)
					rm << i;
								  				
				result 			  += SumUpHiddenMarkov (*bySiteResults, *hmm, *hmf, rm, bySiteResults->GetVDim());
				
			}
			else
			{
				for (i=0; i<theTrees.lLength; i++)
					blockValues.theData[i] = ((_Constant*)computationalResults(i))->Value();
					
				blockWiseVar->SetValue (&blockValues);
				return 	   computingTemplate->Compute()->Value();
			}
		}
		
	}

//	if (likeFuncEvalCallCount>=10)
//	{
//		ProfilerDump("\pProfile");
//		ProfilerTerm();
//		exit(0);
//	}
	return result;
}
//_______________________________________________________________________________________

long		_LikelihoodFunction::BlockLength(long index) 
{
	return ((_DataSetFilter*)dataSetFilterList(theDataFilters(index)))->NumberDistinctSites();
}
//_______________________________________________________________________________________

bool		_LikelihoodFunction::HasBlockChanged(long index) 
{
	return ((_TheTree*)LocateVar(theTrees(index)))->HasChanged2(); 
}

//_______________________________________________________________________________________
		
void	  _LikelihoodFunction::RecurseConstantOnPartition (long blockIndex, long index, long dependance, long highestIndex, _Parameter weight, _Matrix& cache)
{
	_CategoryVariable* thisC = (_CategoryVariable*)LocateVar(indexCat.lData[index]);
	if (index<highestIndex)
	{
		if ((!CheckNthBit(dependance,index))||thisC->IsHiddenMarkov())
			RecurseCategory (blockIndex, index+1, dependance,highestIndex,weight);
		else
		{
			thisC->Refresh();
			long nI = thisC->GetNumberOfIntervals ();
			offsetCounter *= nI;
			for (long k = 0; k<nI; k++)
			{
				thisC->SetIntervalValue(k);
				RecurseConstantOnPartition(blockIndex,index+1,dependance, highestIndex,weight*thisC->GetIntervalWeight(k),cache);
				categID+=offsetCounter/nI;
			}
			offsetCounter/=nI;
			if (offsetCounter>1)
				categID-=nI*offsetCounter;
		}
	}
	else
	{
		long 		nI 			  = thisC->GetNumberOfIntervals (), 
					currentOffset = BlockLength(blockIndex),
					hBit 		  = HighestBit(blockDependancies.lData[blockIndex]);
					
		thisC->Refresh();
		
		_Parameter* sR = siteResults->fastIndex();
		_Matrix*	cws= thisC->GetWeights();
	
		_DataSetFilter* df = ((_DataSetFilter*)dataSetFilterList(theDataFilters(blockIndex)));

		for (long k = 0; k<nI; k++)
		{
			thisC->SetIntervalValue(k,!k);
						
			for   (long kk = 0; kk<currentOffset; kk++)
				sR[kk] = 0.0;
				
			//sR -= currentOffset;
			
			if (hBit>index)
			{
				offsetCounter *= nI;
				RecurseCategory (blockIndex,index+1,blockDependancies.lData[blockIndex],hBit,1);
				offsetCounter /= nI;	
			}
			else
				ComputeBlock (blockIndex,sR);	
				
			_Parameter prod = 0.0;
			{
			for   (long kk = 0; kk<currentOffset; kk++,sR++)
				prod +=  myLog (*sR)*df->theFrequencies.lData[kk];
			}
			prod += myLog (cws->theData[k]*weight);
			cache.theData[categID] = prod;
			
			categID+=offsetCounter;
			sR -= currentOffset;
		}
			
		if (offsetCounter>1)
			categID-=nI*offsetCounter;	
	}
}

//_______________________________________________________________________________________
		
void	  _LikelihoodFunction::RecurseCategory(long blockIndex, long index, long dependance, long highestIndex, _Parameter weight
#ifdef _SLKP_LFENGINE_REWRITE_
											,_SimpleList* siteMultipliers  
#endif				
	)
{
	_CategoryVariable* thisC = (_CategoryVariable*)LocateVar(indexCat.lData[index]);
	if (index<highestIndex)
	{
		if ((!CheckNthBit(dependance,index))||thisC->IsHiddenMarkov())
			RecurseCategory (blockIndex, index+1, dependance,highestIndex,weight
#ifdef _SLKP_LFENGINE_REWRITE_
			,siteMultipliers  
#endif				
			);
		else
		{
			thisC->Refresh();
			long nI = thisC->GetNumberOfIntervals ();
			offsetCounter *= nI;
			for (long k = 0; k<nI; k++)
			{
				thisC->SetIntervalValue(k);
				RecurseCategory(blockIndex,index+1,dependance, highestIndex,weight*thisC->GetIntervalWeight(k)
#ifdef _SLKP_LFENGINE_REWRITE_
								,siteMultipliers  
#endif				
				);
				categID+=offsetCounter/nI;
			}
			offsetCounter/=nI;
			if (offsetCounter>1)
				categID-=nI*offsetCounter;
		}
	}
	else
	{
		if (thisC->IsHiddenMarkov())
		{
			if (offsetCounter == 1) // this is the only categ for the block
				ComputeBlock (blockIndex,siteResults->fastIndex());						
		}
		else
		{
			long 		hDim 			= siteResults->GetHDim(),
						nI 				= thisC->GetNumberOfIntervals (), 
						currentOffset 	= BlockLength(blockIndex);
						
			thisC->Refresh();
			
			_Parameter* sR 	= siteResults->fastIndex();
			_Matrix*	cws	= thisC->GetWeights();
			
#ifdef _SLKP_LFENGINE_REWRITE_
			long	*   siteCorrectors	= ((_SimpleList**)siteCorrections.lData)[blockIndex]->lLength?
										  (((_SimpleList**)siteCorrections.lData)[blockIndex]->lData) + categID * currentOffset
										  :nil;
#endif			
			
			for (long k = 0; k<nI; k++)
			{
				thisC->SetIntervalValue		(k,!k);
				ComputeBlock				(blockIndex,sR+hDim);						
				_Parameter					localWeight = cws->theData[k]*weight ;
				for (long r1 = 0, r2 = hDim; r1 < currentOffset; r1++,r2++)
				{
					#ifdef _SLKP_LFENGINE_REWRITE_
					if (siteCorrectors)
					{
						long scv = *siteCorrectors;
						
						if (scv < siteMultipliers->lData[r1]) // this has a _smaller_ scaling factor
						{
							siteMultipliers->lData[r1] = scv;
							sR[r1] = localWeight * sR[r2] + sR[r1] * acquireScalerMultiplier (siteMultipliers->lData[r1] - scv);
						}
						else
						{
							if (scv > siteMultipliers->lData[r1]) // this is a _larger_ scaling factor
								sR[r1] += localWeight * sR[r2] * acquireScalerMultiplier (scv - siteMultipliers->lData[r1]);							
							else // same scaling factors
								sR[r1] += localWeight * sR[r2];
						}
						
						siteCorrectors++;
					}
					else
					#endif
						sR[r1] += localWeight * sR[r2];
				}
				
				categID+=offsetCounter;
				#ifdef _SLKP_LFENGINE_REWRITE_
				if (offsetCounter > 1)
					siteCorrectors += (offsetCounter-1) * currentOffset;
				#endif
			}
			if (offsetCounter>1)
				categID-=nI*offsetCounter;
		}
	}
}

//_______________________________________________________________________________________
		
void	  _LikelihoodFunction::FindMaxCategory(long blockIndex, long index, long dependance, long lowestIndex, long matrixOffset, _Matrix&storage)
{
	_CategoryVariable* thisC = (_CategoryVariable*)LocateVar(indexCat.lData[index]);
	
	if (index>lowestIndex)
	{
		if (!CheckNthBit(dependance,index)) 
			FindMaxCategory (blockIndex, index-1, dependance,lowestIndex,matrixOffset, storage);
		
		thisC->Refresh();
		long nI = thisC->GetNumberOfIntervals ();
		offsetCounter *= nI;
		
		for (long k = 0; k<nI; k++)
		{
			thisC->SetIntervalValue(k);
			FindMaxCategory(blockIndex,index-1,dependance, lowestIndex,matrixOffset,storage);
			categID+=offsetCounter/nI;
		}
		
		offsetCounter/=nI;
		if (offsetCounter>1)
			categID-=nI*offsetCounter;
	}
	else
	{
		thisC->Refresh();
		
		long 	hDim 			= siteResults->GetHDim(),
				nI 				= thisC->GetNumberOfIntervals (),
			 	currentOffset 	= BlockLength(blockIndex);

		_Parameter* sR = siteResults->fastIndex();

		if (categID==0)
			for (long i=0;i<currentOffset;i++)
				sR[i] = 0.;
		
		for (long k = 0; k<nI; k++)
		{
			thisC->SetIntervalValue(k);
			ComputeBlock (blockIndex, sR+hDim);		
			for (long i=0, l=hDim; i<currentOffset; i++,l++)
				if (sR[l] > sR[i] && !CheckEqual (sR[l], sR[i]))
				{
					sR[i] = sR[l];
					storage[matrixOffset+i] = categID;
				}

			
			categID+=offsetCounter;
		}
		if (offsetCounter>1)
			categID-=nI*offsetCounter;
	}
}


//_______________________________________________________________________________________
		
void	  _LikelihoodFunction::WriteAllCategories(long blockIndex, long index, long dependance, long lowestIndex, long matrixOffset, _Matrix&storage)
{
	_CategoryVariable* thisC = (_CategoryVariable*)LocateVar(indexCat(index));
	if (index>lowestIndex)
	{
		if ((!CheckNthBit(dependance,index))||(thisC->IsHiddenMarkov()))
			WriteAllCategories (blockIndex, index-1, dependance,lowestIndex,matrixOffset,storage);
		else
		{
			thisC->Refresh();
			long nI = thisC->GetNumberOfIntervals ();
			offsetCounter *= nI;
			for (long k = 0; k<nI; k++)
			{
				thisC->SetIntervalValue(k);
				WriteAllCategories(blockIndex,index-1,dependance, lowestIndex,matrixOffset,storage);
				categID+=offsetCounter/nI;
			}
			offsetCounter/=nI;
			if (offsetCounter>1)
				categID-=nI*offsetCounter;
		}
	}
	else
	{
		thisC->Refresh();
		long vDim = storage.GetVDim(), nI = thisC->GetNumberOfIntervals (),
			 currentOffset = BlockLength(blockIndex);
		_Parameter* sR = siteResults->fastIndex();
		for (long k = 0; k<nI; k++)
		{
			thisC->SetIntervalValue(k);
			ComputeBlock (blockIndex,sR);		
			long j=currentOffset-1;
			for (; j>=0;j--)
			{
				storage.theData[categID*vDim+matrixOffset+j]=sR[j];
			}
			categID+=offsetCounter;
		}
		if (offsetCounter>1)
			categID-=nI*offsetCounter;
	}
}

//_______________________________________________________________________________________
		
void	  _LikelihoodFunction::WriteAllCategoriesTemplate (long blockIndex, long index, long dependance, long lowestIndex, long vDim, _Parameter* storage, long* dupMap)
{
	_CategoryVariable* thisC = (_CategoryVariable*)LocateVar(indexCat(index));
	if (index>lowestIndex)
	{
		if ((!CheckNthBit(dependance,index))||(thisC->IsHiddenMarkov()))
			WriteAllCategoriesTemplate (blockIndex, index-1, dependance,lowestIndex,vDim,storage,dupMap);
		else
		{
			thisC->Refresh();
			long nI = thisC->GetNumberOfIntervals ();
			offsetCounter *= nI;
			for (long k = 0; k<nI; k++)
			{
				thisC->SetIntervalValue(k);
				WriteAllCategoriesTemplate(blockIndex,index-1,dependance, lowestIndex,vDim,storage,dupMap);
				categID+=offsetCounter/nI;
			}
			offsetCounter/=nI;
			if (offsetCounter>1)
				categID-=nI*offsetCounter;
		}
	}
	else
	{
		thisC->Refresh();
		long nI            = thisC->GetNumberOfIntervals ();
			 
		_Parameter* sR = siteResults->fastIndex();
				  
		for (long k = 0; k<nI; k++)
		{
			thisC->SetIntervalValue(k);
			ComputeBlock (blockIndex,sR);	
			
			_Parameter * rR = storage+categID*vDim;

			for (long k1 = 0; k1 < vDim; k1++)
				rR[k1] = sR[dupMap[k1]];

			categID+=offsetCounter;
		}
		if (offsetCounter>1)
			categID-=nI*offsetCounter;
	}
}

//_______________________________________________________________________________________
void		_LikelihoodFunction::RandomizeList (_SimpleList& orderList, long elements)
{
	long divisor = RAND_MAX_32/elements-1,i,n;
	if (divisor<=0) divisor = 1;
	orderList.Clear();
	for (i=0; i<elements; i++)
		orderList<<-1;
	for (i=0; i<elements; i++)
	{
		do{
			n = genrand_int32()/divisor;
			if (n>elements) n = elements;
		}
		while (orderList(n)>=0);
		orderList[n]=i;
	}
}

//_______________________________________________________________________________________
void		_LikelihoodFunction::CheckFibonacci (_Parameter shrinkFactor)
{
	long n=Fibonacci.lLength-1;
	if (n<0)
	{
		Fibonacci<<1;
		Fibonacci<<1;
		n+=2;
	}
	while (Fibonacci(n)<shrinkFactor)
	{
		Fibonacci<<Fibonacci(n)+Fibonacci(n-1);
		n++;
	}
}
//_______________________________________________________________________________________
long	 _LikelihoodFunction::HasPrecisionBeenAchieved (_Parameter funcVal, bool cleanup)
{
	static _Parameter lastValue = -A_LARGE_NUMBER, callcount = likeFuncEvalCallCount;
	static _Parameter* oldVarValues = nil;
	
	if (cleanup)
	{
		lastValue = 0;
		callcount = likeFuncEvalCallCount;
		if (oldVarValues)
			delete oldVarValues;
		oldVarValues = nil;
		return 0;
	}
	if (funcVal >= A_LARGE_NUMBER) // reset 
	{
		lastValue = 0;
		callcount = likeFuncEvalCallCount;
		if (oldVarValues)
			delete oldVarValues;
		oldVarValues = new _Parameter[indexInd.lLength];
		checkPointer(oldVarValues);
		for (long i=indexInd.lLength-1; i>=0; i--)
			oldVarValues[i]=0.0;
		return 0;
	}
	if (likeFuncEvalCallCount-callcount>maxItersPerVar)
	// too many computations - return "precision has been reached and spit a warning
	{
		_String warningMsg ("Optimization routines returning before requested precision goal met. The maximum iteration number specified by MAXIMUM_ITERATIONS_PER_VARIABLE has been reached");
		ReportWarning (warningMsg);
		// return precision reached
		warningMsg = _String("Last absolute error in ln-lik function was:")&_String ((_Parameter)fabs(funcVal-lastValue));
		ReportWarning (warningMsg);
		if (optimizationPrecMethod>0.5) // return precision in parameters
		{
			long maxIndex;
			_Parameter average = 0.0, max = 0.0, vVal;
			for (long i=0; i<=indexInd.lLength-1; i--)
			{
				vVal = fabs(GetIthIndependent(i)-oldVarValues[i]);
				if (vVal>max)
				{
					max = vVal;
					maxIndex = i;
				}
				average+=vVal;
			}
			warningMsg = _String("Average last step = ")&_String(average/indexInd.lLength)&
						", with maximum occurring at "&(*LocateVar(indexInd(maxIndex))->GetName())&_String(" =")&_String (max);
			ReportWarning (warningMsg);
		}
		return 1;
	} 
	else // check whether the precision has been reached
	{
		bool done = false;
		if (optimizationPrecMethod<0.5) 
		{
			if (relPrec>0.5)
			{
				if (fabs((funcVal-lastValue)/funcVal)<precision)
				{
					done = true;
				}
				else
					lastValue = funcVal;
			}
			else
			{
				if (fabs(funcVal-lastValue)<precision)
				{
					done = true;
				}
				else
					lastValue = funcVal;
			}
			if (done)
			{
				long maxIndex = 0;
				_Parameter average = 0.0, max = 0.0, vVal;
				for (long i=0; i<=indexInd.lLength-1; i--)
				{
					vVal = fabs(GetIthIndependent(i)-oldVarValues[i]);
					if (vVal>max)
					{
						max = vVal;
						maxIndex = i;
					}
					average+=vVal;
				}
				_String warningMsg = _String("Average last step = ")&_String((_Parameter)(average/indexInd.lLength))
				&", with maximum occurring at "&(*LocateVar(indexInd(maxIndex))->GetName())&_String("=")&_String (max);
				ReportWarning (warningMsg);
				return 1;
			}

			for (long i=0; i<=indexInd.lLength-1; i++)
			{
				_Variable* currentVar = LocateVar(indexInd[i]);
				oldVarValues[i]=currentVar->Value();
			}
			
			return 0;
		}
		else
		{
			done = true;
			if (relPrec>0.5)
			{
				for (long i=0; i<=indexInd.lLength-1; i++)
				{
					_Variable* currentVar = LocateVar(indexInd[i]);
					if (done)
						done = fabs ((currentVar->Value()-oldVarValues[i])/currentVar->Value())<precision;
					oldVarValues[i]=currentVar->Value();
				}
			}
			else
			{
				for (long i=0; i<=indexInd.lLength-1; i++)
				{
					_Variable* currentVar = LocateVar(indexInd[i]);
					if (done)
						done = fabs (currentVar->Value()-oldVarValues[i])<precision;
					oldVarValues[i]=currentVar->Value();
				}
			}
			if (done)
			{
					_String warningMsg = _String("Last absolute error in ln-lik was:")&_String ((_Parameter)fabs(lastValue-funcVal));
					ReportWarning(warningMsg);
					return 1;
			}
			lastValue = funcVal;
			return 0;
		}
	}
	return 0;
}
//_______________________________________________________________________________________

void 	_LikelihoodFunction::CheckDependentBounds (void)
// this function makes sure that a constrained optimization starts within the domain 
// of allowed parameter values
{
	if 	(!indexDep.lLength) // nothing to do here
		return; 
	
	long index,
		 i,
		 j,
		 badConstraint;
		 
	_Matrix currentValues (indexDep.lLength,1,false,true),
			lowerBounds   (indexDep.lLength,1,false,true),
			upperBounds   (indexDep.lLength,1,false,true);
			
	bool	ohWell = false;
			
	_Variable* 	cornholio;
	_SimpleList badIndices; //indices of dependent variables which are out of bounds
		
	nonConstantDep = new _SimpleList;

	for (index = 0; index<indexDep.lLength && !ohWell;index++) 
	// check whether any of the dependent variables are out of bounds
	{
		cornholio 						= 	LocateVar(indexDep.lData[index]);
		
		currentValues.theData[index]	=	cornholio->Compute()->Value();
		lowerBounds.theData[index]		=	cornholio->GetLowerBound();
		upperBounds.theData[index]		=	cornholio->GetUpperBound();
		
		bool badApple = currentValues.theData[index]<lowerBounds.theData[index] || currentValues.theData[index]>upperBounds.theData[index];
		if (badApple)	
			badIndices<<index;
		
		if (cornholio->IsConstant())
		{
			badConstraint = indexDep.lData[index];
			ohWell = badApple;
			j	   = index; // for error reporting at the bottom
		}
		else
			(*nonConstantDep) << indexDep.lData[index];
	}

	if (badIndices.lLength && !ohWell) // one of the variables has left its prescribed bounds
	// build a table of dependancies
	{
		_Matrix 	dependancies (indexDep.lLength,indexInd.lLength,true,true);
		
		// element (i,j) represents the dependance of i-th dep var on the j-th ind var
		// 0 - no dep, 
		// 1 -> monotone increase, 
		// -1 -> monotone decrease
		
		for (index = 0; index<indexInd.lLength;index++) 
		{
			_Parameter 			temp = GetIthIndependent(index);
			SetIthIndependent  (index,temp*1.000000000001);
			
			for (j=indexDep.lLength-1;j>-1;j--)
			{
				_Parameter temp1 = GetIthDependent(j);
				if (temp1>currentValues[j])
					dependancies.Store(j,index,1.0);
				else
					if (temp1<currentValues[j])
						dependancies.Store(j,index,-1.0);					
			}
			SetIthIndependent (index,temp);
		}

		// now we can go through the dependant variables which are out of bounds one at a time
		// and attempt to move them back in.
		for (index = badIndices.lLength-1; index>-1;index--) 
		{
			_Parameter temp = GetIthDependent(badIndices(index));
			
			currentValues[badIndices(index)]=temp;
			
			if (temp<lowerBounds[badIndices(index)]) // too low
			{
				// look for an independent which will increase this variable w/o decreasing the others
				// this is equivalent to searching for the matrix of dependancies for a column w/o negative
				// entries,such that row "index" has "1" in it
				for (j = indexInd.lLength-1; j>-1; j--)
				{
					if (dependancies(badIndices(index),j)==1.0) 
					{
						for (i=indexDep.lLength-1; i>-1;i--)
						{
							if (dependancies(i,j)<0) break;
						}
						if (i==-1) break; //match found
					}
				}
				if (j==-1) // no suitable variable found - try random permutations
					break;
				// now nudge the found variable (indexed by "j") upward, until var "index" is within limits again
				// try this by a trivial bisection
				cornholio = LocateVar (indexInd(j));
				_Parameter left = cornholio->Value(), right = cornholio->GetUpperBound(), 
						   middle, lb = cornholio->GetLowerBound();
				temp=right-left>0.001?0.00001:(right-left)/1000.0;
				while (right-left>temp)
				{
					middle = (left+right)/2.0;
					SetIthIndependent(j,middle);
					if (GetIthDependent  (badIndices(index))>lb) // in range
					{
						right = middle;
					}
					else
					{
						left = middle;
					}
				}
				// take "right" as the new value for the ind
				SetIthIndependent(j,right);
			}
			else
			if (temp>upperBounds[badIndices(index)]) // too high
			{
				// look for an independent which will increase this variable w/o decreasing the others
				// this is equivalent to searching for the matrix of dependancies for a column w/o negative
				// entries,such that row "index" has "1" in it
				for (j = indexInd.lLength-1; j>-1; j--)
				{
					if (dependancies(badIndices(index),j)==-1.0) 
					{
						for (i=indexDep.lLength-1; i>-1;i--)
						{
							if (dependancies(i,j)>0) break;
						}
						if (i==-1) break; //match found
					}
				}
				if (j==-1) // no suitable variable found - try random permutations
					break;
				// now nudge the found variable (indexed by "j") upward, until var "index" is within limits again
				// try this by a trivial bisection
				cornholio = LocateVar (indexInd(j));
				_Parameter left = cornholio->Value(), right = cornholio->GetUpperBound(), 
						   middle, ub = cornholio->GetUpperBound();
				temp=right-left>0.001?0.00001:(right-left)/1000.0;
				while (right-left>temp)
				{
					middle = (left+right)/2.0;
					SetIthIndependent(j,middle);
					if (GetIthDependent  (badIndices(index))<ub) // in range
					{
						right = middle;
					}
					else
					{
						left = middle;
					}
				}
				// take "right" as the new value for the ind
				SetIthIndependent(j,right);
			}

		}
		
		// now we check to see whether all is well
		
		if (index==-1)
		{
			// verify that all dependent vars are now within allowed domains
			for (index = indexDep.lLength-1; index>-1;index--) 
			// check whether any of the dependent variables are out of bounds
			{
				cornholio = LocateVar(indexDep(index));
				currentValues[index]=cornholio->Compute()->Value();
				if ((currentValues[index]<lowerBounds[index])||(currentValues[index]>upperBounds[index]))
				{	
					break;
				}		
			}	
			if (index == -1) // we did it!
				return;
		}
		
		// if we get here, then some of the dependent variables couldn't be moved withing bound
		// we will try 10,000 random assignments of independent variables
		// hoping that one of them will do the trick
		// reuse the first two rows of "dependancies" to store
		// the lower and upper bounds for the dependent vars
		
		for (index = 0; index<indexInd.lLength; index++)
		{
			cornholio 			= LocateVar(indexInd.lData[index]);
			dependancies.Store (0,index,cornholio->GetLowerBound());
			dependancies.Store (1,index,(cornholio->GetUpperBound()>10?10:cornholio->GetUpperBound())-dependancies(0,index));
		}
		
		for (i=0; i<10000; i++)
		{
			for (index = 0; index < indexInd.lLength; index++)
			{
				SetIthIndependent	(index,dependancies(0,index)+genrand_real2()*dependancies(1,index));
				for (j = 0; j < nonConstantDep->lLength; j++) 
				// check whether any of the dependent variables are out of bounds
				{
					currentValues.theData[j]	=	LocateVar(nonConstantDep->lData[j])->Compute()->Value();
					if (currentValues.theData[j]<lowerBounds.theData[j] || currentValues.theData[j]>upperBounds.theData[j])
					{	
						badConstraint = nonConstantDep->lData[j];
						break;
					}		
				}	
				if (j == nonConstantDep->lLength) 
					break;
			}
			if(index < indexInd.lLength) 
				break;
		}
		ohWell = i==10000;
	}
	if (ohWell)
	{
		cornholio 			   = LocateVar(badConstraint);
		_String 		* cStr = (_String*)cornholio->GetFormulaString(),
						badX   = (*cornholio->GetName()) & ":=" & *cStr & " must be in [" & lowerBounds[j] & "," & upperBounds[j] &"]. Current value = " & currentValues[j] & ".";
		
		DeleteObject 			 (cStr);

		_String errMsg = _String("Constrained optimization failed, since a starting point within the domain specified for the variables couldn't be found. Set it by hand, or check your constraints for compatibility. Failed constraint:") & badX;
		acknError(errMsg);
	}
}
//_______________________________________________________________________________________
inline _Parameter sqr (_Parameter x)
{
	return x*x;
}


//_______________________________________________________________________________________
void		_LikelihoodFunction::GetInitialValues (void)
// compute pairwise distances for block index
{
	for (long index=0;index<theTrees.lLength; index++)
	{
		_DataSetFilter* df = (_DataSetFilter*)dataSetFilterList(theDataFilters(index));
		_TheTree *t 	   = (_TheTree*)LocateVar(theTrees.lData[index]);
		long mDim 		   = df->NumberSpecies();
		
		if (df->GetData()->GetTT()->IsStandardNucleotide() && df->IsNormalFilter() && mDim<150 && LocateVar(theProbabilities.lData[index])->IsIndependent())
		{
			// if not - use distance estimates
				
			if (t->IsDegenerate()) 
				continue;
				
			long		branchCount=-1,
						i,
						j,
						k,
						m,
						eqCount=0;
					
			_SimpleList history, 
						history2;
						
			_CalcNode*  dumpkopf = t->StepWiseTraversal(true);
			
			bool two = false;
			
			while (dumpkopf)
			{
				branchCount++;
				history.Clear();
				{
					_AVLList hA (&history);
					dumpkopf->ScanForVariables(hA,hA);
					hA.ReorderList ();
				}
				if (history.lLength>=2)
					two = true;
				dumpkopf = t->StepWiseTraversal();
			}
			
			// first of all, construct the matrix of branch sums,
			// which will be topology dependent
			// the branches are associated with the node they terminate with
			_Matrix  theSystem ((mDim-1)*mDim/2,branchCount,true,true), dichotomy (branchCount-mDim,mDim,true,true);
			node<long>* theRoot = &t->GetRoot();
			k = NodePathTraverser (history, theRoot);
			while (k!=-1)
			{
				for (m=history.lLength-1;m>=0;m--)
				{
					dichotomy[history.lData[m]*mDim+k]=1.0;
				}
				k = NodePathTraverser (history, (node<long>*)nil);
			}
			
			for (i=0;i<mDim-1;i++)
			{
				for (j=i+1;j<mDim;j++)
				{
					theSystem.Store(eqCount,i,1.0);
					theSystem.Store(eqCount,j,1.0);
					for (k=0; k<branchCount-mDim; k++)
						if (dichotomy(k,i)!=dichotomy(k,j))
						{
							theSystem.Store(eqCount,mDim+k,1.0);
						}
					eqCount++;
				}
			}
			
			dichotomy.Clear();
			_Matrix transversions (theSystem.GetHDim(),1,false,true),
					transitions   (theSystem.GetHDim(),1,false,true),
					jc (theSystem.GetHDim(),1,false,true),
					transpose(theSystem);
			
			transpose.Transpose();
			_Matrix AAst (transpose);
			AAst*=theSystem;
			_Matrix* LUSystem = (_Matrix*)AAst.LUDecompose();
			if (!LUSystem)
				return;
			// bingo - the system is now formulated
			// now we build the right-hand sides.
			AAst.Clear();
			_Matrix *freq;
			if (df->GetUnitLength()==1)
			 	freq = (_Matrix*)LocateVar(theProbabilities(index))->GetValue();
			else
				freq = df->HarvestFrequencies (1,1,false);
				
 			_Variable* cornholio;
			_Matrix diffs (1,7,false,true);

			long h = df->NumberSpecies();
			// compute the universal frequency weights
			_Parameter tmp,t2,
					   cR = (*freq)[0]+(*freq)[2],
					   cY = (*freq)[1]+(*freq)[3],
					   c1=2*(*freq)[0]*(*freq)[2]/cR,
					   c2=2*(*freq)[1]*(*freq)[3]/cY,
					   c3=2*((*freq)[0]*(*freq)[2]*cY/cR
					   +(*freq)[1]*(*freq)[3]*cR/cY), comps, fP = 1.0-sqr((*freq)[0])-sqr((*freq)[1])-sqr((*freq)[2])-sqr((*freq)[3]);
					   
			if (df->GetUnitLength()>1)
				DeleteObject(freq);
					   
			eqCount = 0;
			for (i=0;i<h-1;i++)
				for (j=i+1; j<h; j++)
				{
					// Tamura - Nei distance formula
					df->ComputePairwiseDifferences (diffs,i,j) ;
					_Parameter Q = diffs.theData[1]+diffs.theData[3]+
									diffs.theData[4]+diffs.theData[6],P1=diffs.theData[2],
									P2=diffs.theData[5];
					comps = Q+P1+P2+diffs.theData[0];
					if (comps==0.0) continue;
					if (two)
					{
						Q/=comps;
						P1/=comps;
						P2/=comps;
						tmp = 1-Q/(2*cR*cY);
						if (tmp>0.0)
							transversions[eqCount] = -2*cR*cY*log(tmp);
						else
							transversions[eqCount] = 200*cR*cY;
						
						tmp = 1.0-1.0/c1*P1-.5/cR*Q;
						if (tmp>0.0)
							t2 = -c1*log(tmp);
						else
							t2 = c1*100;
							
						tmp = 1.0-1.0/c2*P2-.5/cY*Q;
						if (tmp>0.0)
							t2 -= c2*log(tmp);
						else
							t2 += c2*100;
							
						tmp = 1-.5/(cR*cY)*Q;
						if (tmp>0.0)
							t2 += c3*log(tmp);
						else
							t2 += c3*100;
						
						transitions[eqCount]= t2;
					}
					_Parameter P = (Q+P1+P2)/comps;
					P = 1.0-P/fP;
					if (P>0)
						jc[eqCount]=-fP*log(P);
					else
						jc[eqCount]=20.0;
					eqCount++;
				}
			//DeleteObject(diffs);
				
			_Matrix *trstEst=nil,
			        *trvrEst=nil,
			        *jcEst=nil;
			
			theSystem=transpose;
			theSystem*=jc;
			jc.Clear();
			jcEst = (_Matrix*)LUSystem->LUSolve(&theSystem);
			if (two)
			{
				theSystem=transpose;
				transpose*=transversions;
				theSystem*=transitions;
				transitions.Clear();
				transversions.Clear();
				trstEst = (_Matrix*)LUSystem->LUSolve(&theSystem);
				trvrEst = (_Matrix*)LUSystem->LUSolve(&transpose);
			}
			DeleteObject (LUSystem);
			// now that we have the estimates set the values of branch parameters 
			// the numbering is as follows
			// the tips : 0 to # tips in the same order as in the tree string
			// the branches: # tips + 1 to # tips + #branches
			// within each branch:
				// if there is one ind variable present: set to sum of both
				// if two parameters present - set first one to transition and the 2nd one to transversion
				// for the third parameter and afterwards use the average of the two
			// produce a list of depthwise traversed nodes
			dumpkopf = t->StepWiseTraversal (true);
			dumpkopf = t->StepWiseTraversal(false);
			eqCount = 0; // used to count tips
			h = mDim; // used to count internal nodes
			
			while (dumpkopf)
			{
				history.Clear();
				history2.Clear();
				{
					_AVLList h1 (&history),
							 h2 (&history2);
							 
					dumpkopf->ScanForVariables(h1,h2);
					
					h1.ReorderList();
					h2.ReorderList();
				}
				if (t->IsCurrentNodeATip())
				{
					c3 = (*jcEst)[eqCount];
					if (two)
					{
						c1 = (*trstEst)[eqCount];
						c2 = (*trvrEst)[eqCount++];
					}
					else 	
						eqCount++;
				}
				else
				{
					c3 = (*jcEst)[h];
					if (two)
					{
						c1 = (*trstEst)[h];
						c2 = (*trvrEst)[h++];
					}
					else 	
						h++;
				}
				if (history.lLength==1)
				{
					cornholio = LocateVar(history.lData[0]);
					if (!cornholio->HasChanged())
						cornholio->CheckAndSet (c3*.66667);
				}
				else
				if (history.lLength>=2)
				{
					cornholio = LocateVar(history.lData[0]);
					if (!cornholio->HasChanged())
						cornholio->CheckAndSet (c1);
					cornholio = LocateVar(history.lData[1]);
					if (!cornholio->HasChanged())
						cornholio->CheckAndSet (c2);
					c1 = (c1+c2)/2;
					for (i=2; i<history.lLength;i++)
					{
						cornholio = LocateVar(history.lData[i]);
						if (!cornholio->HasChanged())
							cornholio->CheckAndSet (c1);
					}
				}
				dumpkopf = t->StepWiseTraversal(false);
			}
				
			DeleteObject(trstEst);
			DeleteObject(trvrEst);
			DeleteObject(jcEst);
		}
		else
		{
			_Parameter  	initValue = 0.1;
			checkParameter  (globalStartingPoint, initValue, 0.1);
			
			_SimpleList		indeps;
			_AVLList 		iavl (&indeps);
			
			t->ScanForVariables (iavl, iavl);

			for (long vc = 0; vc < indeps.lLength; vc++)
			{
				//char buf[512];
				_Variable * localVar = LocateVar(indeps.lData[vc]);
				if (localVar->IsIndependent() && !localVar->HasChanged() && !localVar->IsGlobal())
				{
						localVar->CheckAndSet (initValue);
				//		sprintf (buf,"[PRESET]%s = %g\n", localVar->GetName()->sData, localVar->Compute()->Value());
				}
				//else
				//	sprintf (buf,"[PRESET]%s = %g\n", localVar->GetName()->sData, localVar->Compute()->Value());
				//BufferToConsole(buf);
			}


		}
	}
}
//_______________________________________________________________________________________
	
void		_LikelihoodFunction::SetReferenceNodes (void)
{
	_Parameter			cv;
	
	checkParameter 		(useDuplicateMatrixCaching, cv, 0.);
	
	if (cv>0.5)
	{
		_List   			mappedNodes;
		_SimpleList			mappedTo,
							canMap; 
							
		long				i;
							
		for (i=0; i<theTrees.lLength; i++)
		{
			_TheTree * cT = ((_TheTree*)(LocateVar(theTrees(i))));
			
			_CalcNode * aNode = cT->DepthWiseTraversal (true);
			
			while (aNode)
			{
				long rV = aNode->CheckForReferenceNode ();
				if (rV >= 0)
				{
					mappedNodes << aNode;
					mappedTo	<< rV;
				}
				else
					canMap << aNode->GetAVariable();
					
				aNode = cT->DepthWiseTraversal (false);
			}
		}
		
		if (mappedNodes.lLength)
		{
			canMap.Sort();
			for (i=0; i<mappedNodes.lLength; i++)
			{
				if (canMap.BinaryFind (mappedTo.lData[i])>=0)
				{
					_CalcNode * travNode = (_CalcNode*)mappedNodes(i);
					travNode->SetRefNode (mappedTo.lData[i]);
					((_CalcNode*)LocateVar(mappedTo.lData[i]))->AddRefNode();
					_String msg = _String ("Matrix for node ") & travNode->GetName()->getStr() & " mapped to " & 
								   LocateVar(mappedTo.lData[i])->GetName()->getStr();
					ReportWarning (msg);
				}
			}
		}
	}
}

//_______________________________________________________________________________________

void	_LikelihoodFunction::InitMPIOptimizer (void)
{
	#ifdef __HYPHYMPI__
		  if (mpiParallelOptimizer)
		  {
			  parallelOptimizerTasks.Clear();
			  int 	size;	   			   			 
			  MPI_Comm_size(MPI_COMM_WORLD, &size);
			  MPICategoryCount = 1;
			  
			  long		cacheSize  = 0;
			  transferrableVars    = 0;
			  
			  _SimpleList cVarList;
			  {
			  	  _AVLList cAVL (&cVarList);
				  for (long i=0; i<theTrees.lLength; i++)
				  {
					  _TheTree * cT = ((_TheTree*)(LocateVar(theTrees(i))));
					   cT->ScanForCVariables (cAVL);
					   cacheSize += BlockLength (i);
				  }
				  cAVL.ReorderList ();
			  }
			  
			  
			  for (long j=0; j<cVarList.lLength; j++)
				   MPICategoryCount *= ((_CategoryVariable*)LocateVar(cVarList (j)))->GetNumberOfIntervals();
			  
			  if ((MPICategoryCount>1)&&(size>MPICategoryCount))
			  {
				    if (computingTemplate)
				    {
				    	if (templateKind>1)
				    	{
							  for (long i=1; i<theTrees.lLength; i++)
								   if (((_TheTree*)(LocateVar(theTrees(i))))->categoryCount !=
								   	   ((_TheTree*)(LocateVar(theTrees(i-1))))->categoryCount)
								   {
					  	   				WarnError ("MPI Parallel Optimizer does not presently work with functional computational templates for partitions with unequal number of rate classes");
					  	   				return;
					  	   		   }
				    	}
				    	else
				  	    {
					  	   	WarnError ("MPI Parallel Optimizer does not presently work with inline or HM computational templates");
					  	   	return;
					    }
				    }
				    
				    for (long i=0; i<indexCat.lLength; i++)
				    {
				  	   _CategoryVariable* thisCV = (_CategoryVariable*) LocateVar (indexCat.lData[i]);
				  	   if (thisCV->IsHiddenMarkov()||thisCV->IsConstantOnPartition())
				  	   {
					  	   	WarnError ("MPI Parallel Optimizer does not presently work with Hidden Markov or Constant on Partition category variables");
					  	   	return;
					   }
				    }

			  		ReportWarning    (_String ("InitMPIOptimizer with:") & MPICategoryCount & " categories on " & (long)size & " MPI nodes");
			  		_String	  		 sLF (8192L, true),	
			  						 errMsg;
			  						 
					SerializeLF		 (sLF,4);
					sLF.Finalize 	 ();
					
									
					long			 senderID = 0;
					
					for (long i = 1; i<=MPICategoryCount; i++)
						MPISendString (sLF,i);
						
					for (long i = 1; i<=MPICategoryCount; i++)
					{
						_String 	*mapString  = MPIRecvString (i,senderID);
						
						//ReportWarning	 (_String("Got:\n")& *mapString & "\n from node " & (long)senderID);
						
						_List 		*varNames	= mapString->Tokenize (";");
						_SimpleList varMap;
						
						long		mappedVars = indexCat.lLength;
														
						for (long vi = 0; vi < indexCat.lLength; vi++)
						{
							_String * searchTerm = LocateVar(indexCat.lData[vi])->GetName();
							
							long f = varNames->Find (searchTerm);
							if (f<0)
							{
								errMsg = _String ("InitMPIOptimizer: Failed to map cat. var ") & *searchTerm & " for MPI node " & i &". Had variable string:" & *mapString;
								FlagError (errMsg);
							}
							varMap << f;
						}

						for (long vi = 0; vi < indexInd.lLength; vi++)
						{
							_String * searchTerm = LocateVar(indexInd.lData[vi])->GetName();
							
							long f = varNames->Find (searchTerm);
							if ((f<0)&&(!LocateVar(indexInd.lData[vi])->IsGlobal()))
							{
								errMsg = _String ("InitMPIOptimizer: Failed to map independent. var ") & *searchTerm & " for MPI node " & i &". Had variable string:" & *mapString;
								FlagError (errMsg);
							}
							varMap << f;
							if (f>=0)
								mappedVars ++;
						}
						
						parallelOptimizerTasks && & varMap;
						
						if (transferrableVars > 0)
						{
							if (transferrableVars != mappedVars)
							{
								errMsg = _String ("InitMPIOptimizer: Inconsistent variable counts between spawned instances.");
								FlagError (errMsg);
							}
						}
						else
							transferrableVars = mappedVars;
						
						DeleteObject (varNames);
						DeleteObject (mapString);
						
						//ReportWarning (((_String*)parallelOptimizerTasks.toStr())->getStr());
					}
					
					CreateMatrix (&varTransferMatrix, 1, transferrableVars, false, true, false);
					CreateMatrix (&resTransferMatrix, 2, cacheSize, false, true, false);
					CreateMatrix (&catFreqMatrix,1,MPICategoryCount,false,true, false);
					
					ReportWarning(_String("InitMPIOptimizer:Finished with the setup. ") & transferrableVars & " transferrable parameters, " & cacheSize & " sites to be cached.");
			  }
		}
		else
		{
		  int 	size,
		  		slaveNodes;
		  			   			   			 
		  MPI_Comm_size(MPI_COMM_WORLD, &size);
		  slaveNodes = size-1;

		  _Parameter 	  aplf = 0.0;
		  checkParameter (autoParalellizeLF,aplf,0.0);
		  
		  if (theDataFilters.lLength != 1 || ((_DataSetFilter*)dataSetFilterList(theDataFilters.lData[0]))->NumberDistinctSites() < size)
		  	  aplf = 0.0;
		  
		  if (mpiPartitionOptimizer || aplf > 0.5)
		  {
		  	  
			  parallelOptimizerTasks.Clear();
			  
		  	  if (!mpiPartitionOptimizer)
		  	  {
		  	  	  _String css ("_CONTEXT_SWITCH_MPIPARTITIONS_");
		  	  	  // send a context switch signal
		  	  	  for (long ni = 1; ni < size; ni++)
		  	  	  	  MPISendString (css, ni);
		  	  	  // receive confirmation of successful switch
		  	  	  for (long ni = 1; ni < size; ni++)
		  	  	  {
		  	  	  	  long fromNode = ni;
		  	  	  	  
		  	  	  	  _String t (MPIRecvString (ni,fromNode));
		  	  	  	  if (!t.Equal (&css))
		  	  	  	  {
					  	 WarnError (_String("Failed to confirm MPI mode switch at node ") & ni);
					  	 return;
		  	  	  	  }
		  	  	  	  else
		  	  	  	  	 ReportWarning (_String("Successful mode switch to MPI partition optimizer confirmed from node ") & ni);
		  	  	  }
		  	  	  	
		  	  	  	  
		  	  	  lfAutoMPI = true;
		  	  	  mpiPartitionOptimizer = true;
		  	  }

			  if (theDataFilters.lLength>1 && slaveNodes<=theDataFilters.lLength || lfAutoMPI)
			  {
				    if (computingTemplate)
				    {
					  	 WarnError ("MPI Partition Parallel Optimizer does not presently work with computational templates");
					  	 return;
				    }
				    
					long			 senderID = 0,
									 fromPart = 0,
									 toPart   = 0;
					
					_String	  		 errMsg;
					
			  						 
					transferrableVars = 0;		 


			 		_SimpleList	 	  copiedIvars 	(indexInd);
			  		copiedIvars.Sort();				   									 

					if (lfAutoMPI)
					{
						_SimpleList    * optimalOrder = (_SimpleList*)(optimalOrders (0));
						_DataSetFilter * aDSF 	      = (_DataSetFilter*)dataSetFilterList(theDataFilters.lData[0]);
						
						long 		  sitesPerNode	=  optimalOrder->lLength / slaveNodes,
									  overFlow 		=  optimalOrder->lLength % slaveNodes,
									  unitLength	=  aDSF->GetUnitLength();
						
						_SimpleList * dfMap  	= &aDSF->theMap,
									* dupMap 	= &aDSF->duplicateMap,
									* orOrder	= &aDSF->theOriginalOrder;

						if (overFlow)
							overFlow = slaveNodes/overFlow;
							
				  		ReportWarning    (_String ("InitMPIOptimizer (autoLF) with:") & (long)optimalOrder->lLength & " site patterns on " & (long)slaveNodes 
				  						   									 & " MPI computational nodes. " & sitesPerNode & " site patterns per node (+1 every " 
				  						   									 & overFlow & "-th node)");
	
						for (long i = 1; i<size; i++)
						{
							toPart = sitesPerNode;
							if (overFlow && i%overFlow == 0)
								toPart++;
							
							if (fromPart+toPart > optimalOrder->lLength || i == size-1)
								toPart = optimalOrder->lLength-fromPart;
							
							_SimpleList		map		(optimalOrder->lLength, 0, 0),
											subset;
							
							for (long i2 = 0; i2 < toPart; i2++)
								map.lData[optimalOrder->lData[fromPart+i2]] = 1;
								
							//ReportWarning    (_String((_String*)map.toStr()));

							for (long i2 = 0; i2 < dupMap->lLength; i2++)
								if (map.lData[dupMap->lData[i2]])
									for (long cp = 0; cp < unitLength; cp++ )
										subset << orOrder->lData[i2*unitLength+cp];
							
							ReportWarning (_String((_String*)subset.toStr()));
							fromPart += toPart;
							
							_String	  		 sLF (8192L, true);				  						 
							SerializeLF		 (sLF,0,nil,&subset);
							sLF.Finalize 	 ();
							
							//ReportWarning (sLF);
							
							MPISendString (sLF,i);
							_SimpleList 	  dummy;
							parallelOptimizerTasks && &dummy; 
						}
					}
					else
					{
						long	 perNode		=  theDataFilters.lLength / slaveNodes, 
								 overFlow   	=  theDataFilters.lLength % slaveNodes;
								 
								 
						if (overFlow)
							overFlow = slaveNodes/overFlow;
							
				  		ReportWarning    (_String ("InitMPIOptimizer with:") & theDataFilters.lLength & " partitions on " & (long)slaveNodes 
				  						   									 & " MPI computational nodes. " & perNode & " partitions per node (+1 every " 
				  						   									 & overFlow & "-th node)");
				  						   									 
				  						   									 
						for (long i = 1; i<size; i++)
						{
							toPart = perNode;
							if (overFlow && i%overFlow == 0)
								toPart++;
							
							if (fromPart+toPart > theDataFilters.lLength || i == size-1)
								toPart = theDataFilters.lLength-fromPart;
							
							_SimpleList		subset (toPart,fromPart,1);
							
							_String pStr ((_String*)subset.toStr());
							ReportWarning (pStr);
							fromPart += toPart;
							
							_String	  		 sLF (8192L, true);				  						 
							SerializeLF		 (sLF,0,&subset);
							sLF.Finalize 	 ();
							
							MPISendString (sLF,i);
							_SimpleList 	  dummy;
							parallelOptimizerTasks && &dummy; 
						}
					}
					
					for (long i = 1; i<size; i++)
					{						
						_String 		*mapString  	= 	MPIRecvString (-1,senderID);
						_List 			*varNames		= 	mapString->Tokenize (";");
						_SimpleList*	indexedSlaveV   = 	(_SimpleList*)parallelOptimizerTasks(senderID-1);
						
						for (long i2 = 0; i2 < varNames->lLength; i2++)
						{
							_Variable   *aVar = FetchVar(LocateVarByName (*(_String*)(*varNames)(i2)));
							long vi2	= (aVar?copiedIvars.BinaryFind (aVar->GetAVariable()):-1);
								 
							if (aVar == nil || vi2 < 0)
							{
								errMsg = _String ("InitMPIOptimizer: Failed to map independent. var ") 
												& *(_String*)(*varNames)(i2) 
												& " for MPI node " & senderID &". Had variable string:" 
												& *mapString;
												
								FlagError (errMsg);
							}
							
							
							(*indexedSlaveV) << aVar->GetAVariable();
						}
														
						if (indexedSlaveV->lLength > transferrableVars)
							transferrableVars = indexedSlaveV->lLength;
						
						DeleteObject (varNames);
						DeleteObject (mapString);
						
					}
					
					CreateMatrix (&varTransferMatrix, 1, transferrableVars, false, true, false);
					ReportWarning(_String("InitMPIOptimizer:Finished with the setup. Maximum of ") & transferrableVars & " transferrable parameters.");
			  }
		  }
		}
	#endif //__HYPHYMPI__
}

//_______________________________________________________________________________________

void	_LikelihoodFunction::CleanupMPIOptimizer (void)
{
	#ifdef __HYPHYMPI__
		if (parallelOptimizerTasks.lLength)
		{
			 varTransferMatrix.theData[0] = -1.e101;
			 if (mpiPartitionOptimizer)
				 for (long i=0; i<parallelOptimizerTasks.lLength; i++)
					 ReportMPIError(MPI_Send(varTransferMatrix.theData, ((_SimpleList*)parallelOptimizerTasks(i))->lLength , MPI_DOUBLE, i+1, HYPHY_MPI_VARS_TAG, MPI_COMM_WORLD),true);			 
			 else
			 {
				 for (long i=0; i<MPICategoryCount; i++)
					 ReportMPIError(MPI_Send(varTransferMatrix.theData, varTransferMatrix.GetVDim(), MPI_DOUBLE, i+1, HYPHY_MPI_VARS_TAG, MPI_COMM_WORLD),true);
					 
				 resTransferMatrix.Clear();
				 catFreqMatrix.Clear();			 	
			 }
			 varTransferMatrix.Clear();
			 parallelOptimizerTasks.Clear();
			 if (lfAutoMPI)
			 {	
				int size;
			    MPI_Comm_size(MPI_COMM_WORLD, &size);
		  	    for (long ni = 1; ni < size; ni++)
		  	   		MPISendString (empty, ni);
		  	   	  
			    lfAutoMPI = false;
			    mpiPartitionOptimizer = false;
			 }
		}
	#endif
}
	
	//_______________________________________________________________________________________
#ifdef _SLKP_LFENGINE_REWRITE_
void			_LikelihoodFunction::SetupLFCaches				(void)
{
	// need to decide which data represenation to use, 
	// large trees short alignments 
	// an acceptable cache size etc
	categID = 0;
	checkPointer(conditionalInternalNodeLikelihoodCaches = new _Parameter*   [theTrees.lLength]);
	checkPointer(branchCaches							 = new _Parameter*   [theTrees.lLength]);
	checkPointer(siteScalingFactors						 = new _Parameter*   [theTrees.lLength]);
	checkPointer(conditionalTerminalNodeStateFlag		 = new long*		 [theTrees.lLength]);
	overallScalingFactors.Populate						  (theTrees.lLength, 0,0);

	for (long i=0; i<theTrees.lLength; i++)
	{
		_TheTree * cT = ((_TheTree*)(LocateVar(theTrees(i))));
		_DataSetFilter *theFilter = ((_DataSetFilter*)dataSetFilterList(theDataFilters(i)));
		
		if (!theFilter->IsNormalFilter())
		{
			conditionalInternalNodeLikelihoodCaches[i] = nil;
			conditionalTerminalNodeStateFlag	   [i] = nil;
			siteScalingFactors					   [i] = nil;
			branchCaches						   [i] = nil;
			siteCorrections.AppendNewInstance	   (new _SimpleList);
			continue;
		}
		long patternCount	= theFilter->NumberDistinctSites(),
		stateSpaceDim	= theFilter->GetDimension (),
		leafCount		= cT->GetLeafCount(),
		iNodeCount		= cT->GetINodeCount(),
		atomSize		= theFilter->GetUnitLength();
		
		checkPointer (conditionalInternalNodeLikelihoodCaches[i] = new _Parameter [patternCount*stateSpaceDim*iNodeCount*cT->categoryCount]);
		checkPointer (conditionalTerminalNodeStateFlag[i]		 = new long		  [patternCount*leafCount]);
		checkPointer (siteScalingFactors[i]						 = new _Parameter [patternCount*iNodeCount*cT->categoryCount]);
		checkPointer (branchCaches[i]							 = new _Parameter [2*patternCount*stateSpaceDim*cT->categoryCount]);
		
		cachedBranches.AppendNewInstance (new _SimpleList (cT->categoryCount,-1,0));
		if (cT->categoryCount == 1)
			siteCorrections.AppendNewInstance (new _SimpleList (patternCount,0,0));
		else
			siteCorrections.AppendNewInstance (new _SimpleList (cT->categoryCount*patternCount,0,0));
		
		for (long k = 0; k < patternCount*iNodeCount*cT->categoryCount; (siteScalingFactors[i])[k] = 1., k++) ; 
		
		// now process filter characters by site / column
		
		_List		 foundCharactersAux; 
		_AVLListX	 foundCharacters (&foundCharactersAux);
		_String		 aState ((unsigned long)atomSize);
		
		char			** columnBlock		= new char*[atomSize]; checkPointer (columnBlock);
		_Parameter		* translationCache	= new _Parameter [stateSpaceDim]; checkPointer (translationCache);
		_GrowingVector  * ambigs			= new _GrowingVector();
		
		for (long siteID = 0; siteID < patternCount; siteID ++)
		{
			siteScalingFactors[i][siteID] = 1.;
			for (long k = 0; k < atomSize; k++)
				columnBlock[k] = theFilter->GetColumn(siteID*atomSize+k); 
			
			for (long leafID = 0; leafID < leafCount; leafID ++)
			{
				long mappedLeaf  = theFilter->theNodeMap.lData[leafID],
				translation;
				
				for (long k = 0; k < atomSize; k++)
					aState.sData[k] = columnBlock[k][mappedLeaf];
				
				translation = foundCharacters.Find (&aState);
				if (translation < 0)
				{
					translation = theFilter->Translate2Frequencies (aState, translationCache);
					if (translation < 0)
					{
						for (long j = 0; j < stateSpaceDim; j++)
							ambigs->Store(translationCache[j]);
						translation = -ambigs->GetUsed()/stateSpaceDim;
					}
					foundCharacters.Insert (new _String(aState), translation);
				}
				else
					translation = foundCharacters.GetXtra (translation);
				conditionalTerminalNodeStateFlag [i][leafID*patternCount + siteID] = translation;
			}
		}
		conditionalTerminalNodeLikelihoodCaches.AppendNewInstance (ambigs);
		delete columnBlock; delete translationCache;
	}
}
	
#endif

//extern long marginalLFEvals, marginalLFEvalsAmb;

//_______________________________________________________________________________________
	
_Matrix*		_LikelihoodFunction::Optimize ()
{	
	char		   buffer [1024];
	
 	if (lockedLFID != -1)
 	{
 		_String nrerr ("Optimize() could not be executed, because another optimization is already in progress.");
 		WarnError (nrerr);
		return new _Matrix (1,1,false,true);
 	}
 	
 	if (indexInd.lLength == 0)
 	{
		_Matrix result (2,indexDep.lLength<3?3:indexDep.lLength, false, true);	
		PrepareToCompute();	
		result.Store (1,0,Compute());
		result.Store (1,1,indexInd.lLength);
		result.Store (1,2,0);
		for (long i=0; i<indexDep.lLength; i++)
		{
			_PMathObj pm = (_PMathObj)(LocateVar (indexDep(i)))->Compute();
			result.Store(0,i,pm->Value());
		}
		DoneComputing();	
 		return (_Matrix*)result.makeDynamic();
 	}
 	
	#ifdef		__MACPROFILE__
	 	ProfilerInit(collectDetailed,bestTimeBase,1000,500);
	#endif
	
	_Parameter  intermediateP,
				wobble 				= 0.,
				maxSoFar 			= -A_LARGE_NUMBER, 
				bestVal, 
				lastMax, 
			    currentPrecision 	= 0.1,
			    keepStartingPoint,
			    bP,
			    optMethodP, 
			    percentDone 		= 0.0, 
			    bigLastMax;
				
	long 		i,
		 		j, 
		 		fnDim 		= MaximumDimension(), 
		 		categCount  = 1,
		 		evalsIn 	= likeFuncEvalCallCount;


	TimerDifferenceFunction (false);
							  	
	isInOptimize = true;
	lockedLFID   = likeFuncList._SimpleList::Find ((long)this);

	RankVariables();
	VerbosityLevel();
	
	#if defined __UNIX__ && ! defined __HEADLESS__
		#ifdef __HYPHYMPI__
			if (_hy_mpi_node_rank == 0)
			{
		#endif
		_Variable 	  * progressFile = CheckReceptacle (&optimizationStatusFile, empty, false);
		progressFileString = nil;
		
		if (progressFile->ObjectClass () == STRING)
			progressFileString = ((_FString*)progressFile->Compute())->theString;
		
		if (verbosityLevel==1)
			UpdateOptimizationStatus (0,0,0,true,progressFileString);
		#ifdef __HYPHYMPI__
			}
		#endif
	#endif
	

	for (i=0; i<theTrees.lLength; i++)
	{
		_TheTree * cT = ((_TheTree*)(LocateVar(theTrees(i))));
		
		if (cT->CountTreeCategories() >categCount)
			categCount = cT->categoryCount;
	}
#ifdef _SLKP_LFENGINE_REWRITE_	
	SetupLFCaches();
#endif

#ifdef __HYPHYMPI__
	InitMPIOptimizer ();
	if (!mpiPartitionOptimizer || !parallelOptimizerTasks.lLength)
	{
		
#endif
	
	checkParameter (cacheSubtrees,precision,1);
	SetReferenceNodes ();
	checkParameter (useFullMST,intermediateP,0.0);
	
	for (i=0; i<theTrees.lLength; i++)
	{
		_TheTree * cT = ((_TheTree*)(LocateVar(theTrees(i))));
		cT->SetUpMatrices(cT->categoryCount);
		#if USE_SCALING_TO_FIX_UNDERFLOW
			cT->AllocateUnderflowScalers (((_DataSetFilter*)dataSetFilterList (theDataFilters(i)))->NumberDistinctSites());
		#endif
		if (mstCache&&(intermediateP>0.5))
		{
			long  cacheSize = mstCache->cacheSize[i];
			if (cacheSize)
			{
				j = cT->GetLeafCount()+cT->GetINodeCount();
				
				_Parameter**		mstResultCacheIndex = new _Parameter* [cacheSize+1];
				checkPointer		(mstResultCacheIndex);
				
				for (long kk=0; kk<cacheSize; kk++)
				{
					_Parameter*		cacheVector = new _Parameter [j*cT->GetCodeBase()];
					checkPointer (cacheVector);
					mstResultCacheIndex[kk] = cacheVector;
				}
				
				mstResultCacheIndex[cacheSize] = nil;
				mstCache->resultCache << (long)mstResultCacheIndex;
				
				
				long		**		mstStateCacheIndex = new long* [cacheSize+1];
				checkPointer		(mstStateCacheIndex);
				
				for (long kk2=0; kk2<cacheSize; kk2++)
				{
					long*		cacheVector = new long [j];
					checkPointer (cacheVector);
					mstStateCacheIndex[kk2] = cacheVector;
				}
				
				mstStateCacheIndex[cacheSize] = nil;
				mstCache->statesCache << (long)mstStateCacheIndex;
				
				
				
				char		**		mstStateNCacheIndex = new char* [cacheSize+1];
				checkPointer		(mstStateNCacheIndex);
				
				j = cT->GetINodeCount();
				
				for (long kk3=0; kk3<cacheSize; kk3++)
				{
					char*		cacheVector = new char [j];
					checkPointer (cacheVector);
					mstStateNCacheIndex[kk3] = cacheVector;
				}
				
				mstStateCacheIndex[cacheSize] = nil;
				mstCache->statesNCache << (long)mstStateNCacheIndex;
				
				mstCache->stashedLeafOrders && ((_SimpleList*)leafSkips(i));
				((_SimpleList*)leafSkips(i))->Clear();
				((_DataSetFilter*)dataSetFilterList.lData[theDataFilters.lData[i]])->MatchStartNEnd(*(_SimpleList*)mstCache->computingOrder(i),*(_SimpleList*)leafSkips(i),(_SimpleList*)mstCache->parentOrder(i));
				ReportWarning (_String("Using Full MST heurisic on block ") & i & " of likelihood function " & ((_String*)likeFuncNamesList (lockedLFID))->getStr());
			}
			else
			{
				mstCache->resultCache  << 0;
				mstCache->statesCache  << 0;
				mstCache->statesNCache << 0;
				_SimpleList tl;
				mstCache->stashedLeafOrders	&& & tl;
			}
		}		
		//if (dupTrees.lData[i] == 0)
		//{
			_DataSetFilter* dsf = (_DataSetFilter*)dataSetFilterList (theDataFilters.lData[i]);
			_Constant*  tC = (_Constant*)cT->TipCount();
			if ((precision>.1)&&((dsf->GetUnitLength()>1)||(tC->Value()>7)))
			{
				cT->BuildTopLevelCache();
				cT->AllocateResultsCache(dsf->NumberDistinctSites());
			}
			DeleteObject (tC);
		//}
	}

	hasBeenSetUp = 1;
	
#ifdef __HYPHYMPI__
	}
#endif
	
	_Matrix variableValues(indexInd.lLength,1,false,true);
	computationalResults.Clear();
	
	#if !defined __UNIX__ || defined __HEADLESS__
		SetStatusBarValue (0,maxSoFar,0);
	#endif


	#ifdef __HYPHYMPI__
		if (_hy_mpi_node_rank == 0)
		{
	#endif
	
	
	#if !defined __UNIX__ || defined __HEADLESS__
		checkParameter (nicetyMacLevel,nicetyLevel,1);
		nicetyLevel  = nicetyLevel>4?4:nicetyLevel;
		if ((fnDim>=16)&&(nicetyLevel<2))
			nicetyLevel = 2;
		if ((fnDim>25)&&(nicetyLevel<3))
			nicetyLevel = 3;
				
		divideBy     = nicetyLevel>0.01?10000000./(long)(exp(log(10.0)*nicetyLevel)):1000000;
		#if defined __MAC__ || defined __WINDOZE__
			DecideOnDivideBy (this);
		#endif
	#endif
	
	#ifdef __UNIX__
		Compute();
	#endif
	
	#ifdef __HYPHYMPI__
		}
	#endif

	bool			skipCG = false;
	checkParameter (skipConjugateGradient,wobble,0.0);
	skipCG = wobble>0.5;
	checkParameter (randomStartingPerturbations,wobble,0.0);
	checkParameter (useLastResults,keepStartingPoint,0.0);
	checkParameter (allowBoundary,go2Bound,1.0);
	checkParameter (useInitialDistanceGuess,precision,1);

	if (floor(keepStartingPoint) == 0.0 && precision>0.1)
		GetInitialValues();
	
	checkParameter  (globalStartingPoint,precision,0.1);
	_Constant 	  c (precision);

	if ((long)keepStartingPoint==0)
	{
		if ((long)wobble==0)
		{
			for (i=0; i<indexInd.lLength; i++)
			{
				_Variable *iv = LocateVar (indexInd.lData[i]);
				if (iv->IsGlobal())
				{
					if (iv->Value() < 1.e-10)
						iv->SetValue(&c);
				}
				else
				{
					if (!iv->HasChanged())
						iv->SetValue(&c);
				}
			}
		}
		else
		{		
			for (i=0; i<indexInd.lLength; i++)
			{
				_Variable *iv = LocateVar (indexInd(i));
				if (!iv->HasChanged())
				{
					_Parameter newV = precision*(1.0+genrand_int32()/(_Parameter)RAND_MAX_32);
					c.SetValue(newV);
					iv->SetValue(&c);
				}
			}
		}	
	}
	
	#if !defined __UNIX__ || defined __HEADLESS__
		SetStatusBarValue (5,maxSoFar,(likeFuncEvalCallCount-evalsIn)/TimerDifferenceFunction(true));
	#endif
	
	CheckDependentBounds();
	
	checkParameter (startingPrecision,currentPrecision,0.1);
	checkParameter (optimizationMethod,optMethodP,4.0);
	checkParameter (optimizationPrecisionMethod,optimizationPrecMethod,0.0);
	checkParameter (optimizationPrecision,precision,0.001);
	checkParameter (maximumIterationsPerVariable,maxItersPerVar,500+1000/indexInd.lLength);
	
	maxItersPerVar *= indexInd.lLength;
	checkParameter (relativePrecision,relPrec,0.0);
	
	#if !defined __UNIX__ || defined __HEADLESS__
		#ifdef __HYPHYMPI__
			if (_hy_mpi_node_rank == 0)
		#endif
		if (terminateExecution) 
		{
			CleanUpOptimize();
			return new _Matrix (1,1,false,true);
		}
	#endif

	int optMethod = optMethodP;
	
	for (j=0; j<indexInd.lLength; j++)
		variableValues[j]=GetIthIndependent(j);
	
	/*for (long pc = 0; pc < indexInd.lLength; pc++)
	{
		_String pv = _String (pc) & ' ' & GetIthIndependent (pc) &'\n';
		StringToConsole (pv);
	}*/

	if (optMethod == 3) // Powell's Method
	{
		_Matrix 		bestSoFar (indexInd.lLength,1,false,true), 
						savedStartingPoint(indexInd.lLength,1,false,true),
						extrapolatePoint(indexInd.lLength,1,false,true),
						temp(indexInd.lLength,1,false,true);
		_List			startingDirections;
		
		_Parameter		saveItMax;
		
		long	indOfMaxIncrease;
		currentPrecision = 0.01;
		for (i=0; i<indexInd.lLength; i++)	
		{
			_Matrix *ithDir = new _Matrix (indexInd.lLength,1,false,true);
			ithDir->theData[i]=1.0;
			startingDirections<<ithDir;
			DeleteObject(ithDir);
			bestSoFar.theData[i]=GetIthIndependent(i);
		}
		indOfMaxIncrease = 1;
		
		HasPrecisionBeenAchieved(2.*A_LARGE_NUMBER);
		
					
		maxSoFar = Compute();

		while (!HasPrecisionBeenAchieved (maxSoFar))
		{
			bP = 0;
			saveItMax = maxSoFar;
			savedStartingPoint.Duplicate(&bestSoFar);
			keepStartingPoint = Compute();
			for (i=0; i<indexInd.lLength; i++)
			{
				lastMax = maxSoFar;
				GradientLocateTheBump (currentPrecision, maxSoFar, bestSoFar, *((_Matrix*)startingDirections(i)));	
				//((_Matrix*)startingDirections(i))->Duplicate(&bestSoFar);
				//*((_Matrix*)startingDirections(i))-=temp;
				#if !defined __UNIX__ || defined __HEADLESS__
					#ifdef __HYPHYMPI__
						if (_hy_mpi_node_rank == 0)
					#endif
					if (terminateExecution) 
					{
						CleanUpOptimize();
						return new _Matrix (1,1,false,true);
					}
				#endif
				if (fabs(lastMax-maxSoFar)>bP)
				{
					indOfMaxIncrease = i;
					bP = fabs(lastMax-maxSoFar);
				}
				if (verbosityLevel>=5)
				{
			 	 	sprintf (buffer,"\nPowell's direction %d current Max = %g\n", i, maxSoFar);
			 	 	BufferToConsole (buffer);
			 	}

			}
			
			/*_Constant * normConst;
			extrapolatePoint.Duplicate(&bestSoFar);
			extrapolatePoint*=2.0;
			extrapolatePoint-=savedStartingPoint;
			normConst = (_Constant*)extrapolatePoint.Abs();
			extrapolatePoint*=1./normConst->Value();
			DeleteObject (normConst);
			
			for (i=0; i<indexInd.lLength; i++)
				SetIthIndependent (i,extrapolatePoint(0,i));
				
			wobble = Compute();
			
			for (i=0; i<indexInd.lLength; i++)
				SetIthIndependent (i,bestSoFar(0,i));

			if (wobble > saveItMax)
			{
				//if (2.*(saveItMax-2.* maxSoFar+wobble) * SQR (
				temp = bestSoFar;
				temp -= savedStartingPoint;
				normConst = (_Constant*)temp.Abs();
				temp*=1./normConst->Value();
				DeleteObject (normConst);
				GradientLocateTheBump (currentPrecision, maxSoFar, bestSoFar, temp);	
				
				if (indOfMaxIncrease != indexInd.lLength-1)
					((_Matrix*)startingDirections(indOfMaxIncrease))->Duplicate(startingDirections(indexInd.lLength-1));
					
				((_Matrix*)startingDirections(indOfMaxIncrease))->Duplicate(&temp);
			}*/
			
			if (verbosityLevel>=5)
			{
			 	 sprintf (buffer,"\nAt Powell's Precision %g  current Max = %g", currentPrecision, maxSoFar);
			 	 BufferToConsole (buffer);
			}

			if (currentPrecision>1e-8)
				currentPrecision/=10;
		}	
	}
	
	if ((optMethod == 4)&&(indexInd.lLength == 1))
		optMethod = 0;

	if (optMethod)
		checkParameter (bracketingPersistence,bP,2.5);
	else
		checkParameter (bracketingPersistence,bP,3);
	if ((optMethod == 4)||(optMethod == 6)||(optMethod == 7)) // gradient descent
	{
		_Matrix bestSoFar (indexInd.lLength,1,false,true);
		for (i=0; i<indexInd.lLength; i++)	
			bestSoFar.theData[i]=GetIthIndependent(i);
		if (fnDim<21)
			checkParameter (intermediatePrecision,intermediateP,.1);
		else
			checkParameter (intermediatePrecision,intermediateP,.1);
			
		if (verbosityLevel>20)
		{
		  sprintf (buffer,"\nGradient Precision %g  Opt Precision = %g", intermediateP, precision);
		  BufferToConsole (buffer);
		}

		if (optMethod!=7)
			ConjugateGradientDescent (0.05, bestSoFar);
		else
			ConjugateGradientDescent(machineEps, bestSoFar);	
		#if !defined __UNIX__ || defined __HEADLESS__
			#ifdef __HYPHYMPI__
				if (_hy_mpi_node_rank == 0)
				{
			#endif
				if (terminateExecution) 
				{
					CleanUpOptimize();
					return new _Matrix (1,1,false,true);
				}
				#ifndef __HEADLESS__
					if (feedbackTreePanel)
						if (windowObjectRefs.Find ((long)feedbackTreePanel) >= 0)
						{
							feedbackTreePanel->BuildTree (true);
							feedbackTreePanel->RenderTree();
						}
						else
							feedbackTreePanel = nil;
				#endif
			#ifdef __HYPHYMPI__
				}
		#endif
		#endif
		maxSoFar = Compute();	
		if (optMethod != 7)
			optMethod = optMethod!=6?0:5;
		currentPrecision = optMethod==7?sqrt(precision):intermediateP;
		percentDone = 10.0;
		HasPrecisionBeenAchieved(2.*A_LARGE_NUMBER);
	}

	if (optMethod==5) // randomized bracketing method
	{
		_SimpleList passOrder;
		_Matrix optimizationStats (5, indexInd.lLength, false,true);
		
		// initialize the statsMatrix
		
		for (j=0; j<indexInd.lLength; j++)
		{
			optimizationStats.Store(0,j,0); // number of passes
			optimizationStats.Store(1,j,0); // average X change
			optimizationStats.Store(2,j,0); // average Y change
			optimizationStats.Store(3,j,currentPrecision); // bracketing precision
			optimizationStats.Store(4,j,1); // probability of optimizing over this leaf
		}
			
		long	totalPasses, passesDone = 0, lastPassed = -1, counterR = 0;
		totalPasses = indexInd.lLength*bP;
		while (passesDone<totalPasses)
		{
			counterR++;
			RandomizeList (passOrder,indexInd.lLength);
			_Parameter passThresh = (_Parameter)(genrand_int32())/RAND_MAX_32;
			if (passesDone>=indexInd.lLength)
				passThresh/=exp((_Parameter)indexInd.lLength*log((_Parameter)(passesDone+1)/(_Parameter)indexInd.lLength));
			for (j=0; j<indexInd.lLength; j++)
			{
				if (optimizationStats(4,passOrder[j])>=passThresh)
				{
					long jj = passOrder[j];
					if (passOrder[j]==lastPassed) continue;
					_Parameter oldXValue = GetIthIndependent(passOrder[j]),
							   oldYValue = passesDone?lastMax:Compute();
					if (!passesDone) maxSoFar = oldYValue;
					bestVal = oldXValue;
					lastMax = maxSoFar;
					LocateTheBump (jj,optimizationStats(3,jj), maxSoFar, bestVal);
					#ifdef __HYPHYMPI__
						if (_hy_mpi_node_rank == 0)
					#endif
					if (terminateExecution) 
					{
						CleanUpOptimize();
						return new _Matrix (1,1,false,true);
					}
					if (verbosityLevel>1)
					{
						_String *s = (_String*)LocateVar(indexInd(passOrder[j]))->GetName();
						BufferToConsole ("\nAt ");
						StringToConsole (*s);
						sprintf (buffer," with bracketing precision %g current Max = %g. Prior value = %g, current value = %g, %ld", optimizationStats(3,jj), maxSoFar, oldXValue, GetIthIndependent(passOrder[j]), counterR);
						BufferToConsole (buffer);
					}
					passesDone++;
					_Parameter currentThresh = optimizationStats(4,jj);
					if (optimizationStats(1,jj)&&(optimizationStats(1,jj)<fabs(oldXValue-GetIthIndependent(jj))))
					{
						currentThresh*=2;
					}
					if (optimizationStats(2,passOrder[j])&&(optimizationStats(2,passOrder[j])<fabs(oldYValue-maxSoFar)))
					{
						currentThresh*=2;
					}
					if (optimizationStats(0,jj)>totalPasses/passesDone) currentThresh/=indexInd.lLength;
					currentThresh/=indexInd.lLength;
					if (currentThresh>.5) currentThresh = .5;
					optimizationStats.Store(0,jj,optimizationStats(0,jj)+1);
					optimizationStats.Store(1,jj,(optimizationStats(1,jj)*(optimizationStats(0,jj)-1)+fabs(oldXValue-GetIthIndependent(jj)))/optimizationStats(0,jj));
					optimizationStats.Store(2,jj,(optimizationStats(1,jj)*(optimizationStats(0,jj)-1)+fabs(oldYValue-maxSoFar))/optimizationStats(0,jj));
					optimizationStats.Store(4,jj,currentThresh);
					optimizationStats.Store(3,jj,currentPrecision/exp(log((_Parameter)10.0)*optimizationStats(0,jj)));
					lastPassed = jj;

					break;
				}
			
			}	
		}
		
		//control pass
		RandomizeList (passOrder,indexInd.lLength);
		for (j=0; j<indexInd.lLength; j++)
		{
			bestVal = GetIthIndependent(j);
			lastMax = maxSoFar;
			LocateTheBump (j,precision, maxSoFar, bestVal);
			if (verbosityLevel>1)
			{
				sprintf (buffer,"\nControl Pass At Var#%ld with precision %g current Max = %g. Prior value = %g, current value = %g", j, currentPrecision, maxSoFar, bestVal, GetIthIndependent(j));
				BufferToConsole (buffer);
			}
		}
		

	}
	
	if (optMethod == 0)
	{
		bool 	  forward = false;
		
		_Parameter averageChange = 0, 
				   divFactor	 = -1.0 , 
				   oldAverage 	 = -1.0, 
				   stdFactor 	 = 10, 
				   loopCounter 	 = 0.0,
				   sY 			 = 0.0, 
				   sXY 		     = 0.0,
				   s 			 = 0.0, 
				   sX 			 = 0.0, 
				   sXX 			 = 0.0, 
				   nPercentDone, 
				   sigma, 
				   lastMaxValue,
				   averageChange2 = 0.0, 
				   //currentPrecision2 = precision,
				   doShuffle,
				   useAdaptiveStep = 0.0;
				   
		long	   stayPut 		= 0, 
				   inCount 		= 0, 
				   termFactor, 
				   lfCount 		= likeFuncEvalCallCount;
				   
		_SimpleList noChange,
					glVars, 
					shuffledOrder;
					
		_List		*stepHistory = nil;
		
		checkParameter (useAdaptiveVariableStep, useAdaptiveStep, 0.0);
		if (useAdaptiveStep>0.5)
		{		
			stepHistory = new _List;
			for (j=0; j<indexInd.lLength; j++)
			{
				_List 		*varHistory = new _List;
				_Constant 	*aVal = new _Constant;
				aVal->theValue = fabs (variableValues[j])*0.1;
				if (aVal->theValue < 0.0001) 
					aVal->theValue = 0.0001;
					
				(*varHistory) << aVal;
				DeleteObject (aVal);
				(*stepHistory) << varHistory;
				DeleteObject (varHistory);
			}
		}
		
		noChange << -1;		   

		#if !defined __UNIX__ || defined __HEADLESS__
			SetStatusBarValue (percentDone,maxSoFar,(likeFuncEvalCallCount-evalsIn)*TimerDifferenceFunction(true));
		#endif
		
		if (indexInd.lLength<8)
			stdFactor = 8.0;
		else
			if (indexInd.lLength<16)
				stdFactor = 4.0;
			else
				if (indexInd.lLength<32)
					stdFactor = 2.5;
				else	
					stdFactor = .5*(1.+sqrt(5.));
					
		for (j=0; j<indexInd.lLength; j++)
			averageChange+=fabs (variableValues[j]-GetIthIndependent(j));
		
		if (averageChange == 0.0)
			averageChange = currentPrecision;

		currentPrecision = precision>.1?precision:.1;	
			
		lastMaxValue = Compute();

		termFactor = stdFactor+1;
		if (termFactor>indexInd.lLength/2)
			termFactor = indexInd.lLength/2;
			
		if (termFactor<3) termFactor = 3;
		
		checkParameter (doShuffleOrder, doShuffle, 0.0);
		
		for (j=0; j<indexInd.lLength; j++)
			if (LocateVar (indexInd.lData[j])->IsGlobal())
				glVars << j;
 
		while (inCount<termFactor)
		{				
			if (likeFuncEvalCallCount-lfCount>maxItersPerVar)
			{
				_String warningMsg ("Optimization routines returning before requested precision goal met. The maximum iteration number specified by MAXIMUM_ITERATIONS_PER_VARIABLE has been reached");
				ReportWarning (warningMsg);
				DeleteObject (stepHistory);
				break;
			}
			if (averageChange <1e-20)
			{
				averageChange = fabs(oldAverage)<1?fabs(oldAverage*oldAverage):.5;
				if (averageChange<1e-20)
					averageChange = 1e-8;
			}

			if (oldAverage<0.0 || stdFactor*averageChange>oldAverage)
				divFactor = stdFactor;
			else
				divFactor = oldAverage/averageChange;
			
			oldAverage = averageChange;

			averageChange  = 1e-25;
			averageChange2 = 1e-25;

			bigLastMax = maxSoFar;
			
			if (verbosityLevel>1)
			{
				sprintf (buffer,"\n\nOptimization Pass %ld (%ld)\n", (long)loopCounter, inCount);
				BufferToConsole (buffer);
			}
			#if defined __UNIX__ && ! defined __HEADLESS__
			else
				if (verbosityLevel==1)
					UpdateOptimizationStatus (maxSoFar,-1,1,true,progressFileString);
			#endif
			
			_SimpleList nc2;
			
			long		ncp = 0,
						jjj;
			
			averageChange2 = 0.0;
			
			if (doShuffle > 0.1)
			{
				for (j=0; j<indexInd.lLength; j++)
					nc2 << j;
					
				shuffledOrder.Subtract (nc2, glVars);
				shuffledOrder.Permute  (1);
				
				for (j=0; j<glVars.lLength; j++)
					shuffledOrder << glVars.lData[j];
					
				nc2.Clear();
				shuffledOrder.Flip();
			}
			
			for (jjj=forward?0:indexInd.lLength-1; forward?(jjj<indexInd.lLength):jjj>=0; forward?jjj++:jjj--)
			{
				if (doShuffle > 0.1)
					j = shuffledOrder.lData[jjj];
				else
					j = jjj;
				
				bool amIGlobal = LocateVar (indexInd.lData[j])->IsGlobal();
				
				#ifdef __HYPHYMPI__
					if (_hy_mpi_node_rank == 0)
				#endif
				if (terminateExecution) 
				{
					CleanUpOptimize();
					return new _Matrix (1,1,false,true);
				}

				bestVal = GetIthIndependent(j);
				lastMax = maxSoFar;
				if (j==noChange.lData[ncp])
				{
					if (ncp<noChange.lLength-1)
						ncp++;
					averageChange2+=bestVal;
					continue;
				}
				
				_List *vH = nil;
				_Parameter brackStep;
				
				if (useAdaptiveStep>0.5)
				{
					vH  = (_List*)(*stepHistory)(j);
					if (vH->lLength>1)
					{
						brackStep = ((_Constant*)(*vH)(vH->lLength-1))->theValue/((_Constant*)(*vH)(vH->lLength-2))->theValue;
						if (brackStep > 0.5)
							brackStep = 0.5;
						brackStep = ((_Constant*)(*vH)(vH->lLength-1))->theValue*brackStep;
						
					}
					else
						brackStep = ((_Constant*)(*vH)(vH->lLength-1))->theValue*0.5;
						
					if (brackStep < 1e-6)
						brackStep = 1e-6;
				}
				else
				{
					if (amIGlobal)
						brackStep = pow(currentPrecision/**(bestVal>1.?pow(e,long(log(bestVal))):1.)*/,.5+loopCounter/(indexInd.lLength*(loopCounter/indexInd.lLength+1)));
					else
						brackStep = currentPrecision;
				}
				
				LocateTheBump (j,brackStep, maxSoFar, bestVal);
				
				_Parameter  cj = GetIthIndependent(j),
							ch = fabs(bestVal-cj);
	
				averageChange  += ch;
				averageChange2 += cj;
					
				if (vH)
				{
					_Constant tc;
					
					if (ch > 1.e-7)
						tc.theValue = ch;
					else
						tc.theValue = brackStep*0.5;
						
					(*vH) && & tc;
				}

				if (ch<currentPrecision/indexInd.lLength)
					nc2 << j;
					
				variableValues[j] = ch;
				
				if (verbosityLevel>1)
				{
					sprintf (buffer,"\n%ld\t%14.10g\t%10.6g (%10.6g, %10.6g) ", j, maxSoFar, cj,  ch, brackStep);
					BufferToConsole (buffer);
					StringToConsole (*LocateVar(indexInd.lData[j])->GetName());
				}
				#if defined __UNIX__ && ! defined __HEADLESS__
					else
						if (verbosityLevel==1)
							UpdateOptimizationStatus (maxSoFar,-1,1,true,progressFileString);
				#endif
				
			}
			
			if (!forward || doShuffle > 0.1)
				nc2.Sort();
				
			averageChange2/=indexInd.lLength;									// mean of parameter values
			averageChange/=(_Parameter)(indexInd.lLength-nc2.lLength+1);		// mean of change in parameter values during the last step
			
			if (glVars.lLength == indexInd.lLength)
			{
				noChange.Clear();
				noChange.Duplicate (&nc2);
			}
			else
				noChange.Subtract	   (nc2,glVars);
			
			if (noChange.lLength==0)
				noChange << -1;
			

			#if !defined __UNIX__ && ! defined __HEADLESS__
				#ifdef __HYPHYMPI__
					if (_hy_mpi_node_rank == 0)
				#endif
				if (feedbackTreePanel)
					if (windowObjectRefs.Find ((long)feedbackTreePanel) >= 0)
					{
						feedbackTreePanel->BuildTree (true);
						feedbackTreePanel->RenderTree();
					}
					else
						feedbackTreePanel = nil;
			#endif
				
			forward = true;
			loopCounter+=1;
			nPercentDone = lastMax-bigLastMax;
			if (nPercentDone == 0)
				nPercentDone = machineEps;

			sigma = 2.0/nPercentDone;
			s+= sigma;
			sX+= sigma*(loopCounter);
			sXX += sigma*(loopCounter)*(loopCounter);
			sY+=nPercentDone*sigma;
			sXY+=(loopCounter)*nPercentDone*sigma;
			
			if (loopCounter>1.0)
			{
				nPercentDone = (sX*sY-s*sXY);
				if (nPercentDone != 0.0)
					nPercentDone = (sXX*sY-sX*sXY)/nPercentDone;
				if (nPercentDone != 0.0)
					nPercentDone = loopCounter/nPercentDone;
				nPercentDone -= .05;
				if ((nPercentDone>0.0)&&(nPercentDone<.9999))
				{
					nPercentDone = pow (nPercentDone,4)*100.0;
					if (nPercentDone>percentDone)
					{
						stayPut = 0;
						percentDone = nPercentDone;
					}
					else
					{
						if (percentDone<50.0)
							percentDone+=1.0;
						stayPut++;
					}
				}
				else
					if (percentDone<50.0)
					{
						stayPut++;
						percentDone+=1.0;
					}
				#if !defined __UNIX__ || defined __HEADLESS__
				SetStatusBarValue (percentDone,maxSoFar,(likeFuncEvalCallCount-evalsIn)/TimerDifferenceFunction(true));
				#else
					if (verbosityLevel==1)
						UpdateOptimizationStatus (maxSoFar,percentDone,1,true,progressFileString);
				#endif
			}
			
			if (verbosityLevel>5)
			{
				sprintf (buffer,"\nAverage Variable Change: %g %g %g %g %ld", averageChange, nPercentDone,divFactor,oldAverage/averageChange,stayPut);
				BufferToConsole (buffer);
			}
				
			currentPrecision = averageChange/divFactor;
			
			if(currentPrecision<1e-6)
				currentPrecision = 1e-6;
			else
			 	if (currentPrecision>averageChange2)
			 		currentPrecision = averageChange2;
					
			if (maxSoFar-lastMaxValue<precision/termFactor)
				inCount++;
			else
				inCount = 0;
				
			lastMaxValue = maxSoFar;
			
			if (!skipCG&&loopCounter&&(indexInd.lLength>1)&&((((long)loopCounter)%indexInd.lLength)==0))
			{
				_Matrix bestMSoFar (indexInd.lLength,1,false,true);
				for (i=0; i<indexInd.lLength; i++)	
					bestMSoFar[i]=GetIthIndependent(i);
				ConjugateGradientDescent (currentPrecision, bestMSoFar);
			}			
		}

		if (optMethod == 7)
		{
			_Matrix bestMSoFar (indexInd.lLength,1,false,true);
			for (i=0; i<indexInd.lLength; i++)	
				bestMSoFar[i]=GetIthIndependent(i);
			ConjugateGradientDescent (currentPrecision*.01, bestMSoFar);
		}
		
		DeleteObject (stepHistory);

	}
	if (optMethod==1)
		SimplexMethod (precision);
	if (optMethod==2)
		Anneal (precision);
	
	_Matrix result (2,indexInd.lLength+indexDep.lLength<3?3:indexInd.lLength+indexDep.lLength, false, true);
	
	//forceRecomputation = true;
	result.Store (1,0,Compute());
	//forceRecomputation = false;
	result.Store (1,1,indexInd.lLength);
	i = 0;
	while (LocateVar(indexInd.lData[i])->IsGlobal())
	{
		i++;
		if (i==indexInd.lLength)
			break;
	}
	result.Store (1,2,i);
	_PMathObj pm;
	for (i=0; i<indexInd.lLength; i++)
	{
		pm = (LocateVar (indexInd(i)))->Compute();
		result.Store(0,i,pm->Value());
	}
	for (i=0; i<indexDep.lLength; i++)
	{
		pm = (LocateVar (indexDep(i)))->Compute();
		result.Store(0,i+indexInd.lLength,pm->Value());
	}
	
	CleanUpOptimize();
	#if !defined __UNIX__ || defined __HEADLESS__
		SetStatusBarValue (-1,maxSoFar,(likeFuncEvalCallCount-evalsIn)/TimerDifferenceFunction(true));
	#endif
	
	#ifdef		__MACPROFILE__
		ProfilerDump("\pProfile");
		ProfilerTerm();
	#endif
	
	return (_Matrix*)result.makeDynamic();
}

//_______________________________________________________________________________________

void _LikelihoodFunction::CleanUpOptimize (void)
{
	
	categID = 0;
	//printf ("Done OPT LF eval %d MEXP %d\n", likeFuncEvalCallCount, matrixExpCount);
#ifdef __HYPHYMPI__
	if (!mpiPartitionOptimizer || !parallelOptimizerTasks.lLength)
	{
#endif
	for (long i=0; i<theTrees.lLength; i++)
	{
		_TheTree * cT = ((_TheTree*)(LocateVar(theTrees(i))));
		cT->CleanUpMatrices();
		cT->KillTopLevelCache();
	}

#ifdef	_SLKP_LFENGINE_REWRITE_
		DeleteCaches (false);
#endif		
		
#ifdef __HYPHYMPI__
	}
	CleanupMPIOptimizer(); 
	if (!mpiPartitionOptimizer || !parallelOptimizerTasks.lLength)
	{
#endif

	if (mstCache)
	{
		_Parameter 		umst = 0.0;
		checkParameter (useFullMST,umst,0.0);
		if (umst>.5)
			for (long kk=0; kk<mstCache->cacheSize.lLength; kk++)
			{
				long cS = mstCache->cacheSize.lData[kk];
				if ((cS>0)&&(mstCache->resultCache[kk]))
				{
					_Parameter ** c1 = (_Parameter**)mstCache->resultCache[kk];
					for (long k2 = 0; k2<cS; k2++)
						delete c1[k2];
					delete c1;

					long ** c2 = (long**)mstCache->statesCache[kk];
					for (long k3 = 0; k3<cS; k3++)
						delete c2[k3];
					delete c2;
					
					char ** c3 = (char**)mstCache->statesNCache[kk];
					for (long k4 = 0; k4<cS; k4++)
						delete c3[k4];
					delete c3;

					((_SimpleList*)leafSkips(kk))->Clear();
					((_SimpleList*)leafSkips(kk))->Duplicate (mstCache->stashedLeafOrders(kk));
					//printf ("\n%s\n", ((_String*)((_SimpleList*)leafSkips(kk))->toStr())->getStr());
				}
			}
		mstCache->resultCache.Clear();
		mstCache->statesCache.Clear();
		mstCache->statesNCache.Clear();
		mstCache->stashedLeafOrders.Clear();
	}
#ifdef __HYPHYMPI__
	}
#endif
	
	
#if defined __UNIX__ && ! defined __HEADLESS__
	if (verbosityLevel==1)
		UpdateOptimizationStatus (0,0,2,true,progressFileString);
#endif

	setParameter (likeFuncCountVar,likeFuncEvalCallCount);
	HasPrecisionBeenAchieved (0,true);
	isInOptimize = false;
	//DoneComputing();
	hasBeenOptimized = true;
	hasBeenSetUp = 0;
	//printf ("\n%d\n",matrixExpCount);
	lockedLFID	 = -1;
	/*printf ("\nOptimization spool:\n");
	for (i=0; i<indexInd.lLength; i++)
	{	
		_Variable* vspool = LocateVar (indexInd(i));
		printf ("%d %s %g\n", i, vspool->GetName()->getStr(), vspool->Compute()->Value());
	}*/
	DeleteObject (nonConstantDep);
	nonConstantDep = nil;
	#ifndef __UNIX__
		divideBy = 0x0fffffff;
	#endif
}

//_______________________________________________________________________________________
	
//_LikelihoodFunction*		_LikelihoodFunction::SearchTopologies ()
// implement a star decomposition type topology search returning the best tree and 
// optimization results for the best tree.

inline	bool CheckOneDStep (_Parameter& val, _Parameter lB, _Parameter uB)
{
	if (val<lB)
	{
		val = lB;
	}
	else
		if (val>uB)
			val = uB;
	return true;
}

//_______________________________________________________________________________________
bool CheckEqual (_Parameter a, _Parameter b)
{
	if (a!=0.0)
	{
		a = (a>b)?(a-b)/a:(b-a)/a;
		return ((a>0.0)?(a<=machineEps):(a>=-machineEps));
	}
	return (b<=machineEps)&&(b>=-machineEps);
}

//_______________________________________________________________________________________
long 	_LikelihoodFunction::BracketOneVar (long index, _Parameter& left, _Parameter& middle, _Parameter& right,  _Parameter& leftValue, _Parameter& middleValue, _Parameter& rightValue, _Parameter& initialStep)
{
	_Variable* curVar = LocateVar (indexInd (index));
	bool	   movingLeft = false, 
			   first = true;
	_Parameter lowerBound = curVar->GetLowerBound(), 
			   upperBound = curVar->GetUpperBound(), 
			   practicalUB = upperBound>DEFAULTPARAMETERUBOUND?DEFAULTPARAMETERUBOUND:upperBound,
			   magR = 1.61803,
			   //r,q,u,d,
			   leftStep = initialStep*.5, 
			   rightStep = initialStep*.5, 
			   lastLStep = -1.0, 
			   lastRStep = -1.0,
			   saveL = -1.,
			   saveM = -1.,
			   saveR = -1.,
			   saveLV = 0.0,
			   saveMV = 0.0,
			   saveRV = 0.0;
	
	long	   funcCounts = likeFuncEvalCallCount;
	
	//if (index==0)
	//{
		//printf ("Index 0\n");
	//}
	
	middle  =  GetIthIndependent (index);
	
	if ((lowerBound>middle)||(upperBound<middle))
		middle = (lowerBound+practicalUB)/2;
	
	if (middle==lowerBound)
		middle = lowerBound+initialStep/5;
	
	if (middle==upperBound)
		middle = upperBound-initialStep/5;

	while (1)
	{
		
		while ((middle-leftStep)<lowerBound)
		{
			leftStep/=8;
			if (leftStep<initialStep/10) 
			{
				if (!first)
				{
					if (go2Bound>.1)
					{
						middle=lowerBound==0.0?PERTURBATION_OF_ZERO:lowerBound;
						SetIthIndependent (index,middle);
						middleValue=Compute();
					}
					return -2;
				}
				else
					middle=lowerBound+initialStep/10;
			}
		}
			
		while ((rightStep+middle)>upperBound)
		{
			rightStep/=8;
			if (rightStep<initialStep/10) 
			{
				if (!first)
				{
					if (go2Bound>.1)
					{
						middle=upperBound;
						SetIthIndependent (index,upperBound);
						middleValue=Compute();
					}
					return -2;
				}
				else
					middle=upperBound-initialStep/10;
			}

		}
		/*if ((lastLStep==leftStep)&&(lastRStep == rightStep))
		{
			if (movingLeft)
			{
				right = left+leftStep+rightStep;
				rightValue=middleValue;
				middleValue = leftValue;
				left = middle-leftStep;
				SetIthIndependent (index, left);
				leftValue = Compute();
			}
			else
			{
				left = right-leftStep-rightStep;
				leftValue=middleValue;
				middleValue = rightValue;
				right = middle+rightStep;
				SetIthIndependent (index, right);
				rightValue = Compute();
			}
		}
		else
		{
			if (lastLStep==leftStep)
			{	
				if(movingLeft)
				{
					middleValue = leftValue;
					left = middle-leftStep;
					SetIthIndependent (index, left);
					leftValue = Compute();	
		 
				}
				else
				{
					left += leftStep;
					leftValue = middleValue;
					middleValue = rightValue;
				}
				
	 			right = middle+rightStep;
				SetIthIndependent (index, right);
				rightValue = Compute();
			}
			else
				if (lastRStep == rightStep)
				{	
					if(movingLeft)
					{
						rightValue = middleValue;
						right -= rightStep;
						middleValue = leftValue;
					}
					else
					{
						middleValue = rightValue;
						right = middle+rightStep;
						SetIthIndependent (index, right);
						rightValue = Compute();	
					}
					
					left = middle-leftStep;
					SetIthIndependent (index, left);
					leftValue = Compute();	
				}
				else
				{
					CheckAndSetIthIndependent (index, middle);
					middleValue = Compute();
					
					left = middle-leftStep;
					SetIthIndependent (index, left);
					leftValue = Compute();	
					
		 			right = middle+rightStep;
					SetIthIndependent (index, right);
					rightValue = Compute();
					
				}
		}*/
		
		left = middle-leftStep;
		if (CheckEqual(left,saveL))
			leftValue = saveLV;
		else
			if (CheckEqual(left,saveR))
				leftValue = saveRV;
			else
			if (CheckEqual(left,saveM))
				leftValue = saveMV;
			else
			{
				SetIthIndependent (index,left);
				leftValue = Compute();
			}
			
		if (CheckEqual(middle,saveL))
			middleValue = saveLV;
		else
			if (CheckEqual(middle,saveR))
				middleValue = saveRV;
			else
			if (CheckEqual(middle,saveM))
				middleValue = saveMV;
			else
			{
				CheckAndSetIthIndependent (index,middle);
				middleValue = Compute();
			}
		
		right = middle+rightStep;
		if (CheckEqual(right,saveL))
			rightValue = saveLV;
		else
			if (CheckEqual(right,saveR))
				rightValue = saveRV;
			else
			if (CheckEqual(right,saveM))
				rightValue = saveMV;
			else
			{
				SetIthIndependent (index,right);
				rightValue = Compute();
			}
		
		saveL = left;
		saveM = middle;
		saveR = right;
		saveLV = leftValue;
		saveMV = middleValue;
		saveRV = rightValue;

		lastLStep = leftStep;
		lastRStep = rightStep;

		if ((rightValue<=middleValue)&&(leftValue<=middleValue))
			break;

		// case 1: /
		if ((rightValue>=middleValue)&&(middleValue>=leftValue))
		{
			if (movingLeft) 
				rightStep/=magR;
			else 
			{
				leftStep = rightStep;
				/*r = (middle-left)*(middleValue-rightValue);
				q = (middle-right)*(middleValue-leftValue);
				if (fabs(q-r)<1e-20)
					d = q<r?-1e-20:1e-20;
				else
					d = q-r;
				u = middle - ((middle-right)*q-(middle-left)*r)/(2.*d);
				if (u>right)
				{
					if ((u-right)/rightStep<8)
						rightStep=u-right;
					else
						rightStep*=magR*magR;
				}
				else
					if ((u>left)&&(u<right))
					{
						leftStep = u-left;
						rightStep = right-u;
						right = u;
					}
					else*/
				rightStep*=magR;
			}
			middle = right;
			movingLeft = false;
		}
		else // case 2
		
			if ((rightValue<=middleValue)&&(middleValue<=leftValue))
			{
				if ((!movingLeft)&&(!first)) 
					leftStep/=magR;
				else
				{
					rightStep = leftStep;
					/*q = (middle-left)*(middleValue-rightValue);
					r = (middle-right)*(middleValue-leftValue);
					if (fabs(q-r)<1e-20)
						d = q<r?-1e-20:1e-20;
					else
						d = q-r;
					u = middle - ((middle-left)*q-(middle-right)*r)/(2.*d);
					if (u<left)
					{
						if ((left-u)/leftStep<8)
							leftStep=left-u;
						else
							leftStep*=magR*magR;
					}
					else
						if ((u>left)&&(u<right))
						{
							leftStep = u-left;
							rightStep = right-u;
							left = u;
						}
						else*/
					leftStep*=magR;
				}
				middle = left;						
				movingLeft = true;
			}
			else 
			{
				if (movingLeft)
					middle = left;
				else
					middle = right;
			}
						
		if (middle>=practicalUB)
		{
			middle = practicalUB;
			SetIthIndependent (index, middle);
			middleValue = Compute();
			break;
		}
		first = false;
		
	}	
	
	if (!CheckEqual(GetIthIndependent (index),middle))
	{
		SetIthIndependent (index, middle);
		middleValue = Compute();
	}
	
	
	
	bracketFCount+=likeFuncEvalCallCount-funcCounts;
	bracketCount++;
	return 0;	
}
//_______________________________________________________________________________________

void	_LikelihoodFunction::CheckStep (_Parameter& tryStep, _Matrix vect, _Matrix* selection)
{
	for (long index = 0; index<indexInd.lLength; index++)
	{
		_Variable*  curVar = LocateVar (indexInd.lData[index]);
		
		_Parameter  Bound ,currentValue,
					locValue = vect.theData[index];
		
		if (fabs(locValue)<1e-14)
		{
			Bound = curVar->GetUpperBound();
			locValue = 0.0;
		}
		else
			if (locValue<0)
				Bound = curVar->GetLowerBound();
			else
				Bound = curVar->GetUpperBound();
				
		
		if (selection) 
			currentValue = selection->theData[index];
		else
			currentValue = GetIthIndependent (index);

		if (Bound>1000.) Bound = 1000.;
		if (locValue>=0)
		{
			while (currentValue+locValue*tryStep>Bound-1e-8)
			{
				tryStep/=5;
				if (tryStep<1e-8) 
				{
					tryStep = 0.;
					return;
				}
			}
		}
		else
		{
			while (currentValue+locValue*tryStep<Bound+1e-8)
			{
				tryStep/=5;
				if (tryStep<1e-8) 
				{
					tryStep = 0.;
					return;
				}
			}
		}
	}
}

//_______________________________________________________________________________________
long 	_LikelihoodFunction::GradientBracketOneVar (_Matrix& gradient, _Matrix& left, _Matrix& middle, _Matrix& right,  _Parameter& leftValue, _Parameter& middleValue, _Parameter& rightValue, _Parameter& initialStep, bool retry)
{
	_Parameter  leftStep 	= initialStep/1.75, 
			  	rightStep 	= initialStep/1.75, 
				lastLStep 	= -1.0,
				lastRStep 	= -1.0;

	bool	    movingLeft = false, 
				first = true;
		
	_Matrix     unit (middle);
	
	long	    index;
	
	for			(index=0; index<indexInd.lLength; index++)
	{
		_Variable*  thisVar = LocateVar(indexInd.lData[index]);
		if ((thisVar->Value()-thisVar->GetLowerBound()<1e-8)||(thisVar->GetUpperBound()-thisVar->Value()<1e-8))
			gradient.theData[index] = 0.;
		else
			first = false;
	}
	
	if (first) return -1;
	
	first = true;
	// now try to obtain a bracket
	while (1)
	{
		
		_Matrix	  negGrad (gradient);
		negGrad *= -1.;
		CheckStep (leftStep, negGrad, nil);
		CheckStep (rightStep, gradient, nil);
		if ((!leftStep)||(!rightStep)) 
		{
			if (!retry)
				if (!leftStep)
				{
					//left = right;
					//middle = right;
					left = gradient;
					left *= rightStep;
					middle += left;
					return GradientBracketOneVar (gradient, left,middle,right,leftValue,middleValue,rightValue, initialStep, true);
				}
				else
				{
					if (!rightStep)
					{
						//right = left;
						left = gradient;
						left *= leftStep;
						middle -= left;
						return GradientBracketOneVar (gradient, left,middle,right,leftValue,middleValue,rightValue, initialStep, true);
					}
					else 
						return -1;
			}
				
			return -1;
		}
		else
		 	if (leftStep+rightStep<initialStep*.5)
		 		return -2;
		
		
		if ((lastLStep==leftStep)&&(lastRStep == rightStep))
		{
			if (movingLeft)
			{
				right = gradient;
				right*=(leftStep+rightStep);
				right += left;
				rightValue=middleValue;
				middleValue = leftValue;
				left = gradient;
				left*=(-leftStep);
				left += middle;
				for (index=0;index<indexInd.lLength; index++)
					SetIthIndependent (index, left.theData[index]);
				leftValue = Compute();
			}
			else
			{
			    left = gradient;
			    left*=-(leftValue+rightStep);
				left += right;
				leftValue=middleValue;
				middleValue = rightValue;
				right = gradient;
				right*=rightStep;
				right += middle;
				for (index=0;index<indexInd.lLength; index++)
					SetIthIndependent (index, right.theData[index]);
				rightValue = Compute();
			}
		}
		else
		{
			if (lastLStep==leftStep)
			{	
				if(movingLeft)
				{
					middleValue = leftValue;
					left = gradient;
					left*=(-leftStep);
					left += middle;
					for (index=0;index<indexInd.lLength; index++)
						SetIthIndependent (index, left.theData[index]);
					leftValue = Compute();	
		 
				}
				else
				{
					_Matrix dummy;
					dummy = gradient;
					dummy*=leftStep;
					left += dummy;
					leftValue = middleValue;
					middleValue = rightValue;
				}
				
	 			right = gradient;
	 			right*=rightStep;
	 			right += middle;
				for (index=0;index<indexInd.lLength; index++)
					SetIthIndependent (index, right.theData[index]);
				rightValue = Compute();
			}
			else
				if (lastRStep == rightStep)
				{	
					if(movingLeft)
					{
						rightValue = middleValue;
						_Matrix dummy;
						dummy =  gradient;
						dummy*=(-rightStep);
						right += dummy;
						middleValue = leftValue;
					}
					else
					{
						middleValue = rightValue;
						right = gradient;
						right*=rightStep;
						right += middle;
						for (index=0;index<indexInd.lLength; index++)
							SetIthIndependent (index, right.theData[index]);
						rightValue = Compute();	
					}
					left = gradient;
					left*=(-leftStep);
					left += middle;
					for (index=0;index<indexInd.lLength; index++)
						SetIthIndependent (index, left.theData[index]);
					leftValue = Compute();	
				}
				else
				{
					for (index=0;index<indexInd.lLength; index++)
						SetIthIndependent (index, middle.theData[index]);
					middleValue = Compute();
					
					left = gradient;
					left*=(-leftStep);
					left += middle;
					for (index=0;index<indexInd.lLength; index++)
						SetIthIndependent (index, (left.theData)[index]);
					leftValue = Compute();	
					
		 			right = gradient;
		 			right*=rightStep;
		 			right += middle;
					for (index=0;index<indexInd.lLength; index++)
						SetIthIndependent (index, right.theData[index]);
					rightValue = Compute();
					
				}
		}


		lastLStep = leftStep;
		lastRStep = rightStep;
		
		if (verbosityLevel>1)
		{
			char    buffer[128];
			sprintf (buffer,"%12.8g(%10.5g)\t%12.8g\t%12.8g(%10.5g)\n",leftValue,leftStep,middleValue,rightValue,rightStep);
			BufferToConsole (buffer);
		}
		#if defined  __UNIX__ && !defined __HEADLESS__
			else
				if (verbosityLevel==1)
					UpdateOptimizationStatus (middleValue,-1,1,true,progressFileString);
		#endif

	
		if ((rightValue<=middleValue)&&(leftValue<=middleValue))
		{
			break;
		}
		
		// case 1: /
		
		if ((rightValue>=middleValue)&&(middleValue>=leftValue))
		{
			if (movingLeft) rightStep/=4;
			middle = right;
			if (!movingLeft) 
			{
				leftStep = rightStep;
				lastLStep = leftStep;
				rightStep*=4;
			}
			movingLeft = false;
		}
		else // case 2
		
			if ((rightValue<=middleValue)&&(middleValue<=leftValue))
			{
				if ((!movingLeft)&&(!first)) 
				{
					leftStep/=4;
				}
				middle = left;
				if (movingLeft) 
				{
					rightStep = leftStep;
					lastRStep = rightStep;
					leftStep*=4;
				}
				movingLeft = true;
			}
			else // non-parabolic behavior!
				{
					// try to skip the point in the direction that we were moving
					middle = movingLeft?left:right;
				}
						
		first = false;
		
	}	
	
	bool taint = false;
	for (index=0; index<indexInd.lLength; index++)
		if (!CheckEqual(GetIthIndependent (index),middle.theData[index]))
		{
			taint = true;
			SetIthIndependent (index, middle.theData[index]);
		}
	if (taint) 			
		middleValue = Compute();
	
	if (verbosityLevel>1)
	{
		char 	buffer[64];
		sprintf (buffer,"%12.8g\n",middleValue);
		BufferToConsole (buffer);
	}

	return 0;	
}
//_______________________________________________________________________________________

_PMathObj	_LikelihoodFunction::CovarianceMatrix (_SimpleList* parameterList)
{
	if (indexInd.lLength==0)
		return nil;
		
	_Parameter	h = 1.e-5,//STD_GRAD_STEP, 
				functionValue, 
				t1,
				t2,
				cm; // small h and L(x_opt)
	
	checkParameter (covariancePrecision,cm,1.0);
	
	long		i = parameterList?parameterList->lLength:indexInd.lLength,
				j;

	PrepareToCompute();

	functionValue = Compute();
	_Variable      *thisVar;

	bool		useIndirectIndexing = false;

	if (parameterList)
		useIndirectIndexing = true;
	else
		parameterList = &indexInd;

	if (cm<1.)
	// use likelihood profile with the appropriate signifcance level
	{
		_Matrix 	sigLevels  (i,3,false,true);
		
		// find the appropriate significance level for chi2
		
		{
			_String xxc ("_xx_");
			thisVar = CheckReceptacle (&xxc, empty);
			thisVar->SetBounds (-1.e30,1.e30);
		}
		
		if (cm>0.0)
		{
			
			_String		fString = _String ("CChi2(_xx_,1)-") & cm;
			_Formula 	CChi2Fla (fString,nil,true);
			t1 = CChi2Fla.Brent (thisVar,0,0);
			if (fabs(CChi2Fla.Compute()->Value())>1.e-6)
			{
				_String errMsg ("Failed to compute chi-square significance level in call to CovarianceMatrix");
				WarnError (errMsg);
				DoneComputing();
				return	  (_PMathObj)sigLevels.makeDynamic();
			}
		}
		else
		{
			if (cm<0.0)
				t1 = -2.*cm;
			else
			{
				_String errMsg ("Must have a non-zero COVARIANCE_PRECISION in call to CovarianceMatrix.");
				WarnError (errMsg);
				DoneComputing();
				return	  (_PMathObj)sigLevels.makeDynamic();
			}
		}

		#ifndef __UNIX__
			_Parameter totalCount    = 2*parameterList->lLength;
					   
			long	   finishedCount = 0;
					   
			TimerDifferenceFunction	 (false);		   			
			DecideOnDivideBy (this);
		#endif

		
		thisVar->SetNumericValue (t1);

		
		/* ____________________________________ NEW CODE BY AFYP _____________________________________ */
		long		mastodon	= likeFuncList._SimpleList::Find((long)this);
		_String *	myName;
		
		if (mastodon < 0)
		{
			mastodon	= scfgList._SimpleList::Find((long)this);
			myName		= (_String *) scfgNamesList (mastodon);
			
		}
		else
		{
			// generate HBL
			myName		= (_String*)likeFuncNamesList (mastodon);
		}
		
		_String		fString		= _String("function _profileFit(_xxv_,_variableIndex){SetParameter(")&*myName&",_variableIndex,_xxv_);LFCompute("
		//	&*myName&(",_xxres);  fprintf (stdout,\"\\n\",_xxv_,\" \",_xxres);  return _xxres;}");
			&*myName&(",_xxres);return _xxres;}");
		
		
		/* ___________________________________ ! NEW CODE BY AFYP ____________________________________ */
		
						  
		t1 = functionValue-.5*t1;

		_ExecutionList exL (fString);
		exL.Execute ();
				
		for (i=0;i<parameterList->lLength;i++)
		{
			j = useIndirectIndexing?parameterList->lData[i]:i;	
			t2 = GetIthIndependent (j);
			_Variable* thisVar2 = LocateVar (indexInd.lData[j]);
			thisVar->SetBounds (thisVar2->GetLowerBound(),thisVar2->GetUpperBound());
			sigLevels.Store (i,1,t2);
			
			char buffer[255];
			sprintf (buffer,"%.14g",t1);
			
			fString = _String("_profileFit(_xx_,") & j & ")-(" & buffer& ')';
			_Formula 	FitFla (fString,nil,true);
			if (CheckEqual(t2,thisVar2->GetLowerBound()))
				sigLevels.Store (i,0,t2);
			else
			{	
				h = FitFla.Brent (thisVar,t2+1,t2,t2*0.0001+0.000001);
				sigLevels.Store (i,0,MAX(h,thisVar2->GetLowerBound()));
			}
				
			#ifndef __UNIX__
				finishedCount += 1;
				if (TimerDifferenceFunction(true)>1.)
				{
					SetStatusBarValue (finishedCount/totalCount*100.,1,0); 
					TimerDifferenceFunction (false);
				}
			#endif


			if (CheckEqual(t2,thisVar2->GetUpperBound()))
				sigLevels.Store (i,2,t2);
			else
			{
				h = FitFla.Brent (thisVar,t2,t2,t2*0.0001+0.000001);
				sigLevels.Store (i,2,MIN(thisVar2->GetUpperBound(),h));
			}

			#ifndef __UNIX__
				finishedCount += 1;
				if (TimerDifferenceFunction(true)>1.)
				{
					SetStatusBarValue (finishedCount/totalCount*100.,1,0); 
					TimerDifferenceFunction (false);
				}
			#endif

			SetIthIndependent (j,t2);
		}
		
		DoneComputing();
		return  (_PMathObj)sigLevels.makeDynamic();
	}

	_Matrix		hessian    (i,i,false,true),
				funcValues (i,5,false,true),
				*iMap = nil;
	
	_AssociativeList * mapMethod = nil;
	
	_Parameter  	uim = 0.0;
	checkParameter  (useIntervalMapping, uim, 0.0);
	
	if (uim > 0.5)
	{
		iMap = new _Matrix (i,3,false,true);
		mapMethod = new _AssociativeList;
		
		checkPointer (iMap);
		checkPointer (mapMethod);
		
	}
		// y,x',x''
	// first check for boundary values and move the parameter values a bit if needed 
	for (i=0;i<parameterList->lLength;i++)
	{
		long 	 dIndex = useIndirectIndexing?parameterList->lData[i]:i;
		thisVar = LocateVar (indexInd.lData[dIndex]);
		t1 = thisVar->Value();
		funcValues.Store (i,3,t1);
		
		_Parameter locH = 1./131072.; //*(t1>10.?exp(log(10.)*((long)log(t1)/log(10.))):1.);
		
		if (locH<1e-7)
			locH = 1e-7;
			
		/*_Parameter locH = 1./1024.,
				   tryH = locH * 0.25,
				   dv1,
				   dv2,
				   lastErr = 1e100;
				   
		SetIthIndependent (dIndex, t1+locH);
		dv1 = Compute();
		SetIthIndependent (dIndex, t1-locH);
		dv1 += Compute()-2.*functionValue;
		dv1 /= 4.*locH*locH;

		while (tryH >= 1e-12)				   
		{
			SetIthIndependent (dIndex, t1+tryH);
			dv2 = Compute ();
			SetIthIndependent (dIndex, t1-tryH);
			dv2 += Compute ()-2.*functionValue;
			dv2 /= 4.*tryH*tryH;
			_Parameter err = fabs (dv2-dv1) / MAX (fabs(dv2), fabs(dv1));
			if (err > lastErr)
				break;
			else
				lastErr = err;
			dv1 = dv2;
			locH = tryH;
			tryH *= 0.25;
		}*/
		
		
		funcValues.Store (i,4,(locH+t1)-t1);
		
		if (t1+locH > thisVar->GetUpperBound())
			SetIthIndependent (dIndex,thisVar->GetUpperBound()-2.0*locH);
		else
			if (t1-locH < thisVar->GetLowerBound())
				SetIthIndependent (dIndex,thisVar->GetLowerBound()+2.0*locH);
				
		if (uim > 0.5)
		{
			_Parameter 		lb  = thisVar->GetLowerBound(),
					   		ub  = thisVar->GetUpperBound(),
					   		y ,
					   		dy2,
					   		dy;
					   		
			_FString		varKeyName (*thisVar->GetName());
			_Constant		varMapMethod;
					   		
			if (CheckEqual (lb, DEFAULTPARAMETERLBOUND) && CheckEqual (ub, DEFAULTPARAMETERUBOUND))
			// use exp<->log map
			{
				y = log (t1-lb);
				dy = exp(y);
				dy2 = dy;
				varMapMethod.SetValue(1.);
			}
			// use logistic map
			else
			{
				y   = (t1-lb)/(ub-lb);
				t1	= y/(1-y);
				y   = log (t1);
				dy  = (ub-lb)*(t1/(t1+1)/(t1+1));
				dy2 = dy*((1-t1*t1)/(t1+1)/(t1+1));
				varMapMethod.SetValue(2.);
			}					   		   
			iMap->Store(i,0,y);
			iMap->Store(i,1,dy);
			iMap->Store(i,2,dy2);
			mapMethod->MStore (&varKeyName, &varMapMethod, true);
		}
	}
	
	
	
	#ifndef __UNIX__
		_Parameter totalCount    = 2*parameterList->lLength;
				   				   
		long	   finishedCount = 0;
				   
		if (cm>1.1)
			totalCount += 2*parameterList->lLength*(parameterList->lLength-1);
		else
			totalCount += (parameterList->lLength*(parameterList->lLength-1)/2);
		
		DecideOnDivideBy (this);
	#endif
	
	// fill in funcValues with L(...,x_i\pm h,...) and 1st derivatives and get 2nd derivatives
	for (i=0;i<parameterList->lLength;i++)
	{
		long			  pIdx = useIndirectIndexing?parameterList->lData[i]:i;
		
		_Parameter		  pVal = GetIthIndependent (pIdx),
						  d1,
						  locH = funcValues (i,4);
		
		SetIthIndependent (pIdx,pVal-locH); // - step
		t1 = Compute();
		funcValues.Store (i,0,t1);
		SetIthIndependent (pIdx,pVal+locH); // + step
		t2 = Compute();
		funcValues.Store (i,1,t2);          // reset value
		SetIthIndependent (pIdx,pVal);
		d1 = (t2-t1)/(2.0*locH);
		// central 1st derivative
		funcValues.Store (i,2,d1);	
		
		t1  = ((t1-functionValue)+(t2-functionValue))/(locH*locH); 
		// Standard central second derivative
		
		if (uim < 0.5)
			hessian.Store (i,i,-t1);			
		else
			hessian.Store (i,i,-(t1*(*iMap)(i,1)*(*iMap)(i,1)+(*iMap)(i,2)*d1));	
			
		#ifndef __UNIX__
			finishedCount += 2;
			if (TimerDifferenceFunction(true)>1.)
			{
				SetStatusBarValue (finishedCount/totalCount*100.,1,0); 
				TimerDifferenceFunction (false);
			}
		#endif
	}
	

	if (cm>1.1)
	{
		// fill in off-diagonal elements using the f-la
		// f_xy = 1/4h^2 (f(x+h,y+h)-f(x+h,y-h)+f(x-h,y-h)-f(x-h,y+h))

		for (i=0;i<parameterList->lLength-1;i++)
		{
			long	    iidx = useIndirectIndexing?parameterList->lData[i]:i;
			
			_Parameter	ival  = GetIthIndependent(iidx),
						locHi = 1/8192.;//funcValues (i,4);
					   
			for (j=i+1; j<parameterList->lLength; j++)
			{
				long	    jidx = useIndirectIndexing?parameterList->lData[j]:j;
				
				_Parameter	jval  = GetIthIndependent(jidx),
							locHj = locHi, //funcValues (j,4),
							a, // f (x+h,y+h)
					   		b, // f (x+h,y-h)
					   		c, // f (x-h,y-h)
					   		d; // f (x-h,y+h)

				SetIthIndependent (iidx,ival+locHi);
				SetIthIndependent (jidx,jval+locHj);
				a = Compute();
				SetIthIndependent (jidx,jval-locHj);
				b = Compute();
				SetIthIndependent (iidx,ival-locHi);
				c = Compute();
				SetIthIndependent (jidx,jval+locHj);
				d = Compute();
				
				t2 = (a-b-d+c)/(4*locHi*locHj);

				if (uim > 0.5)
					t2 *= (*iMap)(i,1)*(*iMap)(j,1);

				hessian.Store (i,j,-t2);
				hessian.Store (j,i,-t2);				
				SetIthIndependent (iidx,ival);
				SetIthIndependent (jidx,jval);
				#ifndef __UNIX__
					finishedCount += 4;
					if (TimerDifferenceFunction(true)>1.)
					{
						SetStatusBarValue (finishedCount/totalCount*100.,1,0); 
						TimerDifferenceFunction (false);
					}
				#endif
			}
		}
		
	}
	else
	{
		// fill in off-diagonal elements using the f-la
		// f_xy = 1/h^2 (f(x+h,y+h)-f(x)-f_x h -f_y h -.5h^2(f_xx+f_yy))
		
		if (CheckEqual(cm,1.))
		{
			for (i=0;i<parameterList->lLength-1;i++)
			{
				_Parameter t3 = GetIthIndependent(useIndirectIndexing?parameterList->lData[i]:i),
						   t5 = hessian(i,i),
						   t6 = funcValues(i,2);
						   
				SetIthIndependent (useIndirectIndexing?parameterList->lData[i]:i,t3+h);
				for (j=i+1; j<parameterList->lLength; j++)
				{
					_Parameter t4 = GetIthIndependent(useIndirectIndexing?parameterList->lData[j]:j);
					SetIthIndependent (useIndirectIndexing?parameterList->lData[j]:j,t4+h);
					t1 = Compute();
					t2 = (t1-functionValue-(t6+funcValues(j,2)-.5*(t5+hessian(j,j))*h)*h)/(h*h);
					hessian.Store (i,j,-t2);
					hessian.Store (j,i,-t2);
					SetIthIndependent (useIndirectIndexing?parameterList->lData[j]:j,t4);
					#ifndef __UNIX__
						finishedCount ++;
						if (TimerDifferenceFunction(true)>1.)
						{
							SetStatusBarValue (finishedCount/totalCount*100.,1,0); 
							TimerDifferenceFunction (false);
						}
					#endif
				}
				SetIthIndependent (useIndirectIndexing?parameterList->lData[i]:i,t3);
			}		
		}
	}
	// undo changes to var values if needed
	
	DoneComputing();
	
	if (iMap)
	{
		DeleteObject (iMap);
		setParameter (intervalMappingMethod, mapMethod);
		DeleteObject (mapMethod);
	}
	
	for (i=0;i<parameterList->lLength;i++)
	{
		t1 = funcValues(i,3);
		t2 = GetIthIndependent (useIndirectIndexing?parameterList->lData[i]:i);
		if (!CheckEqual (t1,t2))
			SetIthIndependent (useIndirectIndexing?parameterList->lData[i]:i,t1);
	}
	return hessian.Inverse();
	//return (_Matrix*)hessian.makeDynamic();
}

//_______________________________________________________________________________________

void	_LikelihoodFunction::ComputeGradient (_Matrix& gradient, _Matrix&unit,  _Parameter& gradientStep, _Matrix& values,_SimpleList& freeze, long order, bool normalize)
{
	long index;
	_Parameter funcValue;
	
	CheckStep (gradientStep,unit,&values);
	if (order>1) 
	{
		_Matrix nG (unit);
		nG*=-1;
		CheckStep (gradientStep,nG,&values);
	}
	if (gradientStep==0)
	{
		return;
	}
	#ifndef __HYALTIVEC__
	if (order==1)
	{
		funcValue = Compute();
		for (index=0;index<indexInd.lLength;index++)
		{
			if (freeze.Find(index)!=-1)
			{
				gradient[index]=0;
				continue;
			}
			SetIthIndependent(index,GetIthIndependent(index)+gradientStep);
			gradient[index]=(Compute()-funcValue)/gradientStep;
			SetIthIndependent(index,GetIthIndependent(index)-gradientStep);
		}
	}
	else
	#endif
	{
		for (index=0;index<indexInd.lLength;index++)
		{
			if (freeze.Find(index)!=-1)
			{
				gradient[index]=0;
				continue;
			
			}
			SetIthIndependent(index,GetIthIndependent(index)-gradientStep);
			_Parameter temp = Compute();
			SetIthIndependent(index,GetIthIndependent(index)+2*gradientStep);
			gradient[index]=(Compute()-temp)/gradientStep/2;
			SetIthIndependent(index,GetIthIndependent(index)-gradientStep);
		}
	
	}
	
	
	// normalize the gradient 
	if (normalize)
	{
		funcValue = 0;
		for (index=0;index<indexInd.lLength;index++)
		{
			funcValue+=gradient.theData[index]*gradient.theData[index];
		}
		funcValue = sqrt(funcValue);
		for (index=0;index<indexInd.lLength;index++)
		{
			gradient[index]/=funcValue;
		}	
	}
}
//_______________________________________________________________________________________

bool	_LikelihoodFunction::SniffAround (_Matrix& values, _Parameter& bestSoFar, _Parameter& step)
{
	for (long index = 0; index<indexInd.lLength; index++)
	{
		
		_Variable*  curVar = LocateVar (indexInd (index));
		_Parameter lowerBound = curVar->GetLowerBound(),tryStep = step, funcValue, 
				   upperBound = curVar->GetUpperBound(), practicalUB = upperBound>1000?1000:upperBound, val = GetIthIndependent(index);
		// try moving backwards
		while (val-tryStep<lowerBound+1e-8)
		{
			tryStep/=2;
			if (tryStep<1e-8)
				break;
		}
		if (tryStep>=1e-8)
		{
			SetIthIndependent(index,val-tryStep);
			funcValue = Compute();
			if (funcValue>bestSoFar)
			{
				bestSoFar = funcValue;
				values[index]=val-tryStep;
				return true;
			}
		
		}
		tryStep = step	;
		while (val+tryStep>practicalUB-1e-8)
		{
			tryStep/=2;
			if (tryStep<1e-8)
				break;
		}
		if (tryStep>=1e-8)
		{
			SetIthIndependent(index,val+tryStep);
			funcValue = Compute();
			if (funcValue>bestSoFar)
			{
				bestSoFar = funcValue;
				values[index]=val-tryStep;
				return true;
			}
		
		}
		SetIthIndependent(index,val);
	
	
	}
	return false;
}

//_______________________________________________________________________________________

void	_LikelihoodFunction::ConjugateGradientDescent (_Parameter precision, _Matrix& bestVal)
{

	_Parameter  gradientStep	 = STD_GRAD_STEP, 
				temp, 
				maxSoFar 		 = Compute(), 
				currentPrecision = .01;
				
	_SimpleList	freeze;
	_Matrix 	unit (bestVal), 
				gradient (bestVal);
	
	long		vl = verbosityLevel, 
				index;
	
	char	buffer[256];
	
	for (index=0; index<unit.GetHDim(); index++)
		unit[index]=1;
		
	if (vl>1)
	{
		sprintf (buffer,"\nConjugate Gradient Pass %d, max so far %g\n",0,maxSoFar);
		BufferToConsole (buffer);
	}
	
	_Matrix 	G (bestVal), 
				H (bestVal), 
				S(bestVal);
				
	_Parameter  gradL, 
				hL, 
				stepL;
				
	ComputeGradient (gradient, unit, gradientStep, bestVal, freeze, 1, false);
	
	G.Duplicate(&gradient);
	H.Duplicate(&gradient);
	_Constant* absV;
	
	for (index = 0; index<200; index++, currentPrecision/=4)
	{
		temp = maxSoFar;
		absV = (_Constant*)gradient.Abs();
		gradL = absV->Value();
		DeleteObject(absV);
	
		if (CheckEqual(gradL,0.0))
			break;
		
		gradient *= 1.0/gradL;
			
		if (currentPrecision < 0.00001)
			currentPrecision = 0.00001;
		
		GradientLocateTheBump(currentPrecision, maxSoFar, bestVal, gradient);
		if (vl>1)
		{
			sprintf (buffer,"Conjugate Gradient Pass %d, max so far %g\n",index+1,maxSoFar);
			BufferToConsole (buffer);		
		}
		if (fabs((maxSoFar-temp)/maxSoFar)<=precision)
				break;
		
		ComputeGradient (gradient, unit, gradientStep, bestVal, freeze, 1, false);
		absV = (_Constant*)G.Abs();
		gradL = absV->Value();
		DeleteObject(absV);
		if (gradL==0.0) 
			break;
			
		//absV = (_Constant*)gradient.Abs();
		//hL = absV->Value();
		//DeleteObject(absV);
		//stepL = hL/gradL;
		
		hL = 0.;
		for (long j=0; j<indexInd.lLength; j++)
		{
			stepL = gradient.theData[j];
			hL+=(stepL-G.theData[j])*stepL;
		}
		stepL = hL/(gradL*gradL);
		//stepL = SQR(stepL);
		G.Duplicate(&gradient);
		H*=stepL;
		H+=G;
		gradient.Duplicate(&H);
		if (terminateExecution)
			return;
		
		#if !defined __UNIX__ && !defined __HEADLESS__
			if (feedbackTreePanel)
				if (windowObjectRefs.Find ((long)feedbackTreePanel) >= 0)
				{
					feedbackTreePanel->BuildTree (true);
					feedbackTreePanel->RenderTree();
				}
				else
					feedbackTreePanel = nil;
		#endif		
	}
	for (index=0;index<indexInd.lLength;index++)
		SetIthIndependent(index,bestVal(0,index));
	if (vl>1)
		BufferToConsole("\n");
	
}
	
//_______________________________________________________________________________________

void	_LikelihoodFunction::GradientDescent (_Parameter& gPrecision, _Matrix& bestVal)
{

	_Parameter  currentPrecision = 0.1, gradientStep=STD_GRAD_STEP, temp, tryStep, maxSoFar = Compute(), bestTry;
	_SimpleList	leastChange, freeze, countLC;
	_Matrix unit (bestVal), gradient (bestVal);
	long	vl = verbosityLevel, index;
	
	for (index=0; index<unit.GetHDim(); index++)
		unit[index]=1;
		
	while ((currentPrecision>=gPrecision)&&(freeze.lLength<indexInd.lLength))
	{
			
		gradientStep = STD_GRAD_STEP;
		if (terminateExecution)
			return;
		
		char	buffer[128];
		ComputeGradient (gradient, unit, gradientStep, bestVal, freeze, 1);
		if (gradientStep==0) break;
		bool done = false;
		bestTry = tryStep = currentPrecision;
		while(!done)
		{
			CheckStep (tryStep,gradient, &bestVal);
			if (tryStep == 0.0)
			{
				long wereFrozen = freeze.lLength;
				for (index=0;index<indexInd.lLength;index++)
				{
					if (freeze.Find(index)>=0) continue;
					_Variable* thisV = LocateVar(indexInd.lData[index]);
					if ((GetIthIndependent(index)-thisV->GetLowerBound()<1.0e-20)
						&& (gradient(0,index)<0.0))
					{
						freeze<<index;
						break;
					}
					if ((-GetIthIndependent(index)+thisV->GetUpperBound()<1.0e-20)
						&& (gradient(0,index)>0.0))
					{
						freeze<<index;
						break;
					}
				}
				tryStep = currentPrecision;
				if (freeze.lLength==wereFrozen)
					return;
				break;
			}
			for (index=0;index<indexInd.lLength;index++)
			{
				SetIthIndependent(index,bestVal(0,index)+tryStep*gradient(index,0));
			}
			temp = Compute();
			if (temp>maxSoFar)
			{
				if (vl>=5)
				{
					sprintf (buffer,"\nMoving down along the gradient with step %g value %g", tryStep, temp);
					BufferToConsole (buffer);
				}
				_Matrix delta;
				delta = gradient;
				delta*=tryStep;
				if ((temp-maxSoFar)/fabs(maxSoFar)<gPrecision)
				{
					maxSoFar = temp;
					bestVal += delta;
					return;
				}
				
				maxSoFar = temp;
				bestVal += delta;
				//see which variable changed the least
				temp = A_LARGE_NUMBER;
				long  suspect,f;
				for (long i = 0; i<indexInd.lLength; i++)
				{
					if (fabs(delta(i,0))<temp)
					{
						if (freeze.Find(i)!=-1) continue;
						temp = fabs(delta(i,0));
						suspect = i;
					}
				
				}
				
				f = leastChange.Find(suspect);
				if (f==-1)
				{
					leastChange<<suspect;
					countLC<<1;
				}
				
				if (tryStep<currentPrecision/10)
					currentPrecision/=10;
				else
					tryStep*=2;
				done = true;
			}
			else
			{
				tryStep/=10;
				if (tryStep<gPrecision) 
				{
					for (index=0;index<indexInd.lLength;index++)
					{
						SetIthIndependent(index,bestVal(0,index));
					}
					if (leastChange.lLength)
					{
						SortLists(&leastChange, &countLC);
						if (vl>=5)
						{
							sprintf (buffer,"\nFreezing Variable %d",leastChange(leastChange.lLength-1));
							BufferToConsole (buffer);
						}
						freeze<<leastChange(leastChange.lLength-1);
						leastChange.Delete (leastChange.lLength-1);
						countLC.Delete (countLC.lLength-1);
						currentPrecision = sqrt(gPrecision);
						done = true;
						break;
					}
					currentPrecision=0;
					break;
				}
				if (vl>=5)
				{
					sprintf (buffer,"\nShrinking step to %g (%g %g)", tryStep, tryStep, gPrecision);
					BufferToConsole (buffer);
				}
			}
		}
 	}
}
			
		

//_______________________________________________________________________________________

void	_LikelihoodFunction::GradientLocateTheBump (_Parameter gPrecision, _Parameter& maxSoFar, _Matrix& bestVal, _Matrix& gradient)
{
	_Parameter  leftValue = maxSoFar, middleValue = maxSoFar, rightValue = maxSoFar, bp = gPrecision*2.0;
	_Matrix	    left (indexInd.lLength,1,false,true), 
				right (indexInd.lLength,1,false,true), middle(indexInd.lLength,1,false,true),
				newMiddle(indexInd.lLength,1,false,true);
	_Parameter  gRatio = (3.0-sqrt(5.0))*.5;
	long		index;
							   
	middle = bestVal;
	
	int  outcome = GradientBracketOneVar(gradient, left,middle,right,leftValue, middleValue, rightValue,bp);
	
	bool reset = false;

	if (outcome!=-1) // successfull bracket
	{
		// set up left, right, middle
		_Parameter lV, rV;
		for (outcome = 0; outcome<indexInd.lLength; outcome++)
		{
			if (fabs(gradient.theData[outcome])>machineEps)
			{
				lV = (left.theData[outcome]-middle.theData[outcome])/gradient.theData[outcome];
				rV = (right.theData[outcome]-middle.theData[outcome])/gradient.theData[outcome];
				break;
			}
		}
		if (outcome == indexInd.lLength)
			reset = true;
		else
		{
			_Parameter U,V,W,X=0.0,E=0,FX,FW,FV,XM,R,Q,P,ETEMP,D,FU;
			W = .0;
			V = .0;
			FX = middleValue;
			FV = FX;
			FW = FX;
			outcome = 0;
			bool pFitGood;
			while (outcome < 20)
			{
				pFitGood = false;
				XM = .5*(lV+rV);
				if (fabs(X-XM)<gPrecision*2) break;
				if (fabs(E)>machineEps)
				{
					R = (X-W)*(FX-FV);
					Q = (X-V)*(FX-FW);
					P = (X-V)*Q-(X-W)*R;
					Q = 2.0 * (Q-R);
					if (Q>0)
						P = -P;
					Q = fabs(Q);
					ETEMP = E;
					E = D;
					if (!((fabs(P)>=fabs(.5*Q*ETEMP))||(P<=Q*(lV-X))||(P>=Q*(rV-X))))
					{
						D = P/Q;
						U = X+D;
						pFitGood = true;
					}
				}
				if (!pFitGood)
				{
						if (X<XM)
							E = lV-X;
						else
							E = rV - X;
						D = gRatio*E;
				}
				U = X+D;
				for (index = 0; index < indexInd.lLength; index++)
					SetIthIndependent (index,middle.theData[index]+U*gradient.theData[index]);
				FU = Compute();
				
				if (FU>=FX)
				{
					if (U>=X)
						lV = X;
					else
						rV = X;
					V = W;
					FV = FW;
					W = X;
					FW = FX;
					X = U;
					FX = FU;
				}
				else
				{
					if (U<=X)
						lV = U;
					else
						rV = U;
					if ((FU>=FW)||(W==X))
					{
						V = W;
						FV = FW;
						W = U;
						FW = FU;
					}
					else
					{
						if ((FU>=FV)||(V==X)||(V==W))
						{
							V = U;
							FV = FU;
						}
					}
					
				}
				outcome++;
			
			}
			
			middleValue = FX;
			for (outcome = 0; outcome < indexInd.lLength; outcome++)
			{
				middle.theData[outcome] += X*gradient.theData[outcome];
				if (!CheckEqual(GetIthIndependent(outcome),middle.theData[outcome]))
					SetIthIndependent (outcome,middle.theData[outcome]);
			}
			bestVal = middle;
			maxSoFar = middleValue;
		}
		//middle = X;
	}
			
	if (outcome == -1)
	{			
		reset = true;
		if (verbosityLevel>1)
		{
			BufferToConsole ("Line optimization unsuccessful\n");
		}
		if (leftValue>middleValue)
		{
			middleValue = leftValue;
			middle = left;
		}
		if (rightValue>middleValue)
		{
			middleValue = rightValue;
			middle = right;
		}
			
		if (middleValue>maxSoFar)
		{
			for (long i=0;i<indexInd.lLength;i++)
			{
				bestVal[i]=middle(i,0);
				if (!CheckEqual(GetIthIndependent(i),middle(i,0)))
					SetIthIndependent (i,middle(i,0));
			}
			maxSoFar = middleValue;
			reset = false;
		}
	}
	if (reset)
		for (long i=0;i<indexInd.lLength;i++)
			SetIthIndependent(i,bestVal.theData[i]);
}

//_______________________________________________________________________________________

void	_LikelihoodFunction::LocateTheBump (long index,_Parameter gPrecision, _Parameter& maxSoFar, _Parameter& bestVal)
{
	_Parameter left, 
			   right, 
			   middle,
			   leftValue, 
			   middleValue, 
			   rightValue,  
			   bp = 2.*gPrecision,
			   gRatio = (3.0-sqrt(5.0))*.5;
	
#ifdef _SLKP_LFENGINE_REWRITE_
	DetermineLocalUpdatePolicy			 ();
	//printf ("%s\n", LocateVar (indexInd.lData[index])->GetName()->sData);
#endif	
	
	int outcome = BracketOneVar(index,left,middle,right,leftValue, middleValue, rightValue,bp);
	long		fCount = likeFuncEvalCallCount;
	if (outcome!=-1) // successfull bracket
	{
		if (outcome == -2)
		{
			if (middleValue>maxSoFar)
			{
				maxSoFar = middleValue;
				bestVal = middle;
				if (GetIthIndependent(index)!=middle)
					SetIthIndependent (index,middle);
			}
			else
				SetIthIndependent (index,bestVal);
#ifdef _SLKP_LFENGINE_REWRITE_
			FlushLocalUpdatePolicy			  ();
#endif		
			return;	
		}
			
		{
			_Parameter U,V,W,X=middle,E=0,FX,FW,FV,XM,R,Q,P,ETEMP,D,FU;
			W = middle;
			V = middle;
			FX = middleValue;
			FV = FX;
			FW = FX;
			outcome = 0;
			bool pFitGood;
			while (outcome < 20)
			{
				pFitGood = false;
				XM = .5*(left+right);
				if (fabs(X-XM)<gPrecision*2.) break;
				if (fabs(E)>1.0e-10)
				{
					R = (X-W)*(FX-FV);
					Q = (X-V)*(FX-FW);
					P = (X-V)*Q-(X-W)*R;
					Q = 2.0 * (Q-R);
					if (Q>0)
						P = -P;
					Q = fabs(Q);
					ETEMP = E;
					E = D;
					if (!((fabs(P)>=fabs(.5*Q*ETEMP))||(P<=Q*(left-X))||(P>=Q*(right-X))))
					{
						D = P/Q;
						U = X+D;
						pFitGood = true;
					}
				}
				if (!pFitGood)
				{
					if (X<XM)
						E = left-X;
					else
						E = right - X;
					D = gRatio*E;
				}
				U = X+D;
				SetIthIndependent (index,U);
				FU = Compute();
				
				if (FU>=FX)
				{
					if (U>=X)
						left = X;
					else
						right = X;
					V = W;
					FV = FW;
					W = X;
					FW = FX;
					X = U;
					FX = FU;
				}
				else
				{
					if (U<=X)
						left = U;
					else
						right = U;
					if ((FU>=FW)||(W==X))
					{
						V = W;
						FV = FW;
						W = U;
						FW = FU;
					}
					else
					{
						if ((FU>=FV)||(V==X)||(V==W))
						{
							V = U;
							FV = FU;
						}
					}
					
				}
				outcome++;
			
			}
		
			middleValue = FX;
			middle = X;
		}
		
		//printf ("/nMax at (%g,%g)", X, FX);
		if (middleValue<maxSoFar)
			SetIthIndependent(index,bestVal);
		else
		{	
			if (!CheckEqual(GetIthIndependent(index),middle))
				SetIthIndependent (index,middle);
			maxSoFar = middleValue;
		}	 	
	}
	
	oneDFCount += likeFuncEvalCallCount-fCount;
	oneDCount ++;
#ifdef _SLKP_LFENGINE_REWRITE_
	FlushLocalUpdatePolicy			  ();
#endif		
}

//_______________________________________________________________________________________
bool		_LikelihoodFunction::checkPermissibility (_Matrix&m, long row)
{
	for (long j=0; j<indexInd.lLength; j++)
	{
		_Parameter junk = m (row,j);
		
		_Variable *v = LocateVar (indexInd(j));

		if (junk<v->GetLowerBound())
		{
			return FALSE;
		}
		else
			if (junk>v->GetUpperBound())
			{
				return FALSE;
			}			
	}		
	
	return true;
}

//_______________________________________________________________________________________
_Parameter		_LikelihoodFunction::computeAtAPoint (_Matrix&m, long row)
{

	if (!checkPermissibility(m,row)) return -1e300;

	for (long k=0; k < indexInd.lLength; k++)
		SetIthIndependent (k, m(row,k));
	
	return Compute();
}

//_______________________________________________________________________________________
_Parameter		_LikelihoodFunction::replaceAPoint (_Matrix&m, long row, _Matrix&p, _Parameter& nV, _Matrix& fv)
{
	for (long k=0; k < indexInd.lLength; k++)
		m.Store (row,k,p(0,k));
	fv.Store (0,row,nV);
	return nV;
}


//_______________________________________________________________________________________
_Parameter		_LikelihoodFunction::SimplexMethod (_Parameter& gPrecision)
{
	_Matrix		points (indexInd.lLength+1, indexInd.lLength, false, true), //the matrix of points
				functionalEvaluations (1, indexInd.lLength+1,false, true),
				scratch1 (1, indexInd.lLength,false, true),
				scratch2 (1, indexInd.lLength,false, true),
				bestSoFar (1, indexInd.lLength+1,false, true),
				temp;
			
	_Parameter  lastBError 		= 0.,
				lastError 		= 0.0,
				simplexAlpha 	= 1.0,
				simplexBeta 	= 0.5, 
				simplexGamma 	= 2.0, 
				testValue, 
				minError 		= 1.e300,
				bumpingFactor 	= 1.;
				
	long		iterationsCount = 0, 
				nBumps 			= 0;
	
	
	// first populate the matrix of initial points
	long j,k;
	
	for (k=0;k<indexInd.lLength;k++)
	{
		_Variable *v = LocateVar (indexInd (k));
		
		_Parameter lowP 	= v->GetLowerBound(), 
				   highP 	= v->GetUpperBound(), 
				   span 	= highP-lowP;
				   
		if (span>1)
		{
			lowP += 0.1;
			highP = lowP+1;
		}
		else
		{
			lowP += span/10;
			highP += lowP+span/2;
		}
		
		for (j=0; j<k; j++)
			points.Store (j,k,lowP);

		for (j=k+1; j<=indexInd.lLength; j++)
			points.Store (j,k,lowP);

		points.Store (k,k,highP);
	}
	
	
	for (j=0; j<=indexInd.lLength; j++)
		functionalEvaluations.Store (0,j,computeAtAPoint(points,j));		
	
	do
	{
		// compute the reflection
		long 			indexMin = 0, 
						indexMax = 0, 
						index2Min;

		_Parameter 		max  = -1e300, 
						min  = 1e300, 
						min2 = 1e300,
						error = 0;
		
		for (j=0; j<=indexInd.lLength; j++)
		{
			if (functionalEvaluations(0,j)>max)
			{
				indexMax = j;
				max = functionalEvaluations(0,j);
			}
			if (functionalEvaluations(0,j)<min)
			{
				if (min<min2)
				{
					min2 = min;
					index2Min = indexMin;
				}
				indexMin = j;
				min = functionalEvaluations(0,j);
			}
			else
			{
				if (functionalEvaluations(0,j)<min2)
				{
					index2Min = j;
					min2 = functionalEvaluations(0,j);
				}
			}
		}
		
		
		
		
		error = functionalEvaluations(0,indexMax)-functionalEvaluations(0,indexMin);
		if (verbosityLevel>1)
		{
			char	buffer[64];
			sprintf (buffer,"\n Error = %15.15g", error);
			BufferToConsole (buffer);
		}
		if (error<minError)
		{
			minError = error;
			for (k=0; k < indexMax; k++)
				bestSoFar.Store(0,k,points(indexMax,k));
		}
			
		// find the centroid of all the points excluding the worst!	
		for (j=0; j<indexInd.lLength; j++)
		{
			_Parameter junk = 0;
			for (k=0; k < indexMin; k++)
				junk+= points(k,j);
			
			for (k=indexMin+1; k <= indexInd.lLength; k++)
				junk+= points(k,j);
			_Variable *v = LocateVar (indexInd(j));
			
			junk/=indexInd.lLength;
			if (junk<v->GetLowerBound())
			{
				scratch1.Store (0,j,v->GetLowerBound());
			}
			else
				if (junk>v->GetUpperBound())
				{
					scratch1.Store (0,j,v->GetUpperBound());
				}
				else
					scratch1.Store (0,j,junk);
				
		}
	
			
		
		if (error<gPrecision) break;
		if (fabs(lastError-error)<gPrecision*error)
		{
			iterationsCount++;
			if (iterationsCount%10==0)
			{
				if ((error>lastBError)||(iterationsCount%(30)==0))
				{
					//perturb the coodinates of simplex vertices to escape a vicios loop
					nBumps++;
					for (k=0;k<indexInd.lLength;k++)
					{

						_Variable *v = LocateVar (indexInd (k));
						_Parameter lowP = v->GetLowerBound() , highP = v->GetUpperBound(), trial;
						
						if (fabs(lastBError-error)>gPrecision)
							bumpingFactor*=10;
						else
							bumpingFactor = 1;
						//trial = bestSoFar(0,k)+bumpingFactor*(genrand_int32()-RAND_MAX_32)/(_Parameter)RAND_MAX_32*error;
//						else
//							trial = scratch1(0,k)+100*(genrand_int32()-RAND_MAX_32)/(_Parameter)RAND_MAX_32*error;
							trial = lowP+(highP-lowP)*(genrand_int32()-RAND_MAX_32)/(_Parameter)RAND_MAX_32*error;
							
						if ((trial<(lowP+gPrecision))||(trial>(highP-gPrecision))) 
							{
								trial = lowP+ genrand_int32()*(highP-lowP)/RAND_MAX_32;
							}
						scratch1.Store (0,k,trial);
						lastBError = error;
						
					}
					if (verbosityLevel>1)
					{
						char	buffer[64];
						sprintf (buffer,"\nBUMPING... with factor %g", bumpingFactor);
						BufferToConsole (buffer);
					}
				}
				
			}
			if (nBumps>25)
			{
				char str[255];
				sprintf (str,"Simplex Method Failed to Converge to Desired Precision. Precision attained: %15.15g", minError);
				ReportWarning (_String(str));
				break;
			}
		}
		else
		{
			iterationsCount = 0;
		}
		lastError = error;

		for (k=0; k < indexInd.lLength; k++)
			scratch2.Store (0,k,points(indexMin,k));
		
		// reflect the worst point across the centroid
		_Matrix dummy;
		temp = scratch2;
		dummy = scratch1;
		dummy -= temp;
		dummy *=(1.0+simplexAlpha); 
		temp += dummy;
			
		testValue = computeAtAPoint (temp);
		
		
		if ((testValue>functionalEvaluations(0,indexMin))&&(testValue<functionalEvaluations(0,indexMax)))
		{
				replaceAPoint (points, indexMin, temp, testValue, functionalEvaluations);
		}
		else
			if (testValue>functionalEvaluations(0,indexMax))
			{
					replaceAPoint (points, indexMin, temp, testValue, functionalEvaluations);
					scratch2 =temp;
					scratch2*=simplexGamma;
					scratch1*=(1-simplexGamma);
					scratch2+=scratch1;
					_Parameter testValue3 = computeAtAPoint (scratch2);
					
					if (testValue3>functionalEvaluations(0,index2Min))
						replaceAPoint (points, index2Min, scratch2,testValue3, functionalEvaluations);
			}
			else
			{
				for (k=0; k < indexInd.lLength; k++)
					scratch2.Store (0,k,points(indexMin,k));
				temp=scratch2;
				temp*=simplexBeta;
				scratch1*=(1-simplexBeta);
				temp+=scratch1;

				_Parameter testValue2 = computeAtAPoint (temp);
				if (testValue2>=testValue)
				{
					replaceAPoint (points, indexMin, temp, testValue, functionalEvaluations);
				}
				else
				{
					for (j=0; j<indexInd.lLength; j++)
						for (k=0; k<=indexInd.lLength; k++)
						{
							points.Store (j,k,(points(j,k)+points(indexMax,k))/2);
						}
				}
			}
	}
	while (1);
	
	for (k=0; k < indexInd.lLength; k++)
		SetIthIndependent (k, bestSoFar(0,k));
	
	return true;
						
}
		
													
//_______________________________________________________________________________________

void	_LikelihoodFunction::Anneal (_Parameter&)
// simple simulated annealing approach
{
	/*_Parameter jumpMagnitude , decreaseRate, startingRate, currentValue = -1e300, trialValue, bestValue = -1e300;
	_Parameter iterationsPerDecrease, decisionOnConvergence, terminateAfter;
	_Matrix	jumpVector (1, indexInd.lLength,false, true),
			currentPoint (1, indexInd.lLength,false, true),
			trialPoint (1, indexInd.lLength,false, true),
			bestPoint (1, indexInd.lLength,false, true);
	long i,count = 0, smCount = 0;
			
	checkParameter (_String("SA_JUMP_MAGNITUDE"),jumpMagnitude,0.1);
	checkParameter (_String("SA_DECREASE_RATE"),decreaseRate,1.25);
	checkParameter (_String("SA_STARTING_RATE"),startingRate,1);
	checkParameter (_String("SA_ITERATIONS_PER_DECREASE"),iterationsPerDecrease,500);
	checkParameter (_String("SA_CONVERGENCE_THERESHOLD"),decisionOnConvergence,100);
	checkParameter (_String("SA_TERMINATION_THERESHOLD"),terminateAfter,1000000);
	
	
	for (i=0; i<indexInd.lLength; i++)
	{
		trialPoint.Store(0,i,GetIthIndependent(0));
	}
		
	
	currentValue = computeAtAPoint (currentPoint);
	while (1)
	{
		count++;
		if (count>terminateAfter)
		{
			computeAtAPoint (bestPoint);
			ReportWarning (_String("Simulated annealing terminated before convergence was attained"));
			return;
		}
		if (count%(long)iterationsPerDecrease==0) 
		{
			startingRate *= decreaseRate;
			iterationsPerDecrease/=decreaseRate;
		}
		trialValue = computeAtAPoint (trialPoint);
		// compare the trial value and the current value and decide what to do;
		if (fabs(trialValue-currentValue)<gPrecision)
		{
			smCount++;
		}
		else
			smCount=0;
			
		if (smCount > decisionOnConvergence)
		{
			computeAtAPoint (bestPoint);
			return;
		}
		
		if (trialValue > currentValue)
		{
			if (VerbosityLevel()>1)
			{
				printf ("\nJump Up to %g accepted", trialValue);
			}
			currentValue = trialValue;
			currentPoint = trialPoint;
		}
		else
		{
			_Parameter decisionProb = exp((trialValue-currentValue)*startingRate);
			if (((_Parameter)genrand_int32())/RAND_MAX_32<decisionProb)
			{
				currentValue = trialValue;
				currentPoint = trialPoint;
				if (VerbosityLevel()>1)
				{
					printf ("\nJump down to %g accepted", currentValue);
				}
			}
			else
				if (VerbosityLevel()>900)
				{
					printf ("\nJump to %g rejected", trialValue);
				}
			
		}
		
		// see if we have the max
		
		if (currentValue>bestValue)
		{
			bestValue = currentValue;
			bestPoint = currentPoint;
		}
		
		// now generate a random jump
		_Parameter distance = jumpMagnitude*genrand_int32()/(_Parameter)RAND_MAX_32, dist1 = 0;
		
		// generate an arbitrary point in the allowed domain
		
		
		for (i=0; i<indexInd.lLength; i++)
		{
			_Variable*v = LocateVar (indexInd(i));
			jumpVector.Store (0,i,v->GetLowerBound()+(v->GetUpperBound()-v->GetLowerBound())*genrand_int32()/(_Parameter)RAND_MAX_32-currentPoint(0,i));
			dist1+=jumpVector(0,i)*jumpVector(0,i);
		}
		
		jumpVector *= (distance/sqrt(dist1)); // normalize the jump;
		
		trialPoint = jumpVector;
		trialPoint += currentPoint;
		
		
	}		*/
	_String errMsg ("Simulated Annealing is yet to be implemented. Sorry about that.");
	WarnError (errMsg);
}

//_______________________________________________________________________________________
void	_LikelihoodFunction::RescanAllVariables (void)
{
	indexCat.Clear();
	indexDep.Clear();
	indexInd.Clear();
	computationalResults.Clear();
	ScanAllVariables();
}

//_______________________________________________________________________________________
long	_LikelihoodFunction::DependOnTree (_String& treeName)
{
	long f = LocateVarByName (treeName);
	if (f>=0)
		return theTrees.Find(variableNames.GetXtra(f));		
	
	return -1;
}

//_______________________________________________________________________________________
long	_LikelihoodFunction::DependOnDS (long ID)
{
	for (long k = 0; k<theDataFilters.lLength; k++)
		if (dataSetList._SimpleList::Find 
			((long)((_DataSetFilter*)dataSetFilterList (theDataFilters.lData[k]))->GetData()) == ID)
			return k;
	return -1;
}

//_______________________________________________________________________________________
long	_LikelihoodFunction::DependOnModel (_String& modelTitle)
{
	long f = modelNames.Find (&modelTitle);
	if ((f>=0)&&(modelTitle.sLength))
	{
		for (long k=0; k<theTrees.lLength; k++)
		{
			_TheTree * t = (_TheTree*)LocateVar (theTrees.lData[k]);
			_CalcNode* cN = t->DepthWiseTraversal (true);
			while (cN)
			{
				if (cN->GetModelIndex()==f)
					return k;
				cN = t->DepthWiseTraversal (false);
			}
		}
	}
	return -1;
}

//_______________________________________________________________________________________
_CategoryVariable*	_LikelihoodFunction::FindCategoryVar (long index)
{
	_CategoryVariable* res = nil;
	if ((index>=0)&&(index<blockDependancies.lLength))
		res = (_CategoryVariable*)LocateVar(indexCat(HighestBit(blockDependancies.lData[index])));
	return res;
}

//_______________________________________________________________________________________
void	_LikelihoodFunction::ScanAllVariables (void)
{
	_SimpleList	  allVariables,
				  covCat,
				  cpCat;
				  
	long		  i;
	
	{
		_AVLList avl (&allVariables);
		
		for (i=0; i<theProbabilities.lLength; i++)
			(LocateVar(theProbabilities(i)))->ScanForVariables (avl,true);
			
		if (computingTemplate)
			computingTemplate->ScanFForVariables (avl,true,false,true);

		avl.ReorderList();
	}

	if (templateKind<0)
		allVariables.Delete (allVariables.Find(-templateKind-1));
		
	{
		_AVLList iia (&indexInd),
				 iid (&indexDep);
	
		for (i=0; i<allVariables.lLength; i++)
		{
			_Variable* theV = ((_Variable*)LocateVar(allVariables(i)));
			if (theV->IsCategory())
			{
				if (((_CategoryVariable*)theV)->IsUncorrelated())
				{
					if (((_CategoryVariable*)theV)->IsConstantOnPartition())
						indexCat << allVariables(i);
					else
						cpCat << allVariables(i);
				}
				else
					covCat<<allVariables(i);
				continue;
			}
			if (theV->IsIndependent())
				iia.Insert ((BaseRef)allVariables(i));
			else
				iid.Insert ((BaseRef)allVariables(i));
		}
						
		for (i=0; i<theTrees.lLength; i++)
			((_TheTree*)(LocateVar(theTrees(i))))->ScanForGVariables (iia, iid);
		
		for (i=0; i<theTrees.lLength; i++)
		{
			_TheTree * cT = ((_TheTree*)(LocateVar(theTrees(i))));
			cT->ScanForVariables  (iia, iid);
			cT->ScanForDVariables (iid, iia);
			cT->SetUp();
		}
		
		iia.ReorderList ();
		iid.ReorderList ();
	}

	/* mod 10/28/2005; Need to make two passes on trees; one to check for all category variables defined on the tree
	   and second to lay out dependency bits (after all the category variables have been concatenated) */

	for (i=0; i<theTrees.lLength; i++)
	{
		_SimpleList localCategVars;
		{
			_AVLList ca (&localCategVars);
			((_TheTree*)(LocateVar(theTrees(i))))->ScanForCVariables (ca);
			ca.ReorderList();
		}
		
		for (long i2 = 0; i2 < localCategVars.lLength; i2++)
		{
			_CategoryVariable * cvRef = (_CategoryVariable*)LocateVar (localCategVars.lData[i2]);
			if (cvRef->IsUncorrelated())
			{
				if (cvRef->IsConstantOnPartition())
					indexCat >> cvRef->GetAVariable();
				else
					cpCat >> cvRef->GetAVariable();
			}
			else
				covCat >> cvRef->GetAVariable();		
		}
	}

	indexCat << cpCat;
	indexCat << covCat;

	if (indexCat.lLength>sizeof(long)*8)
	{
		_String errMsg ("The number of category variables exceeded ");
		errMsg = errMsg&_String((long)sizeof(long)*8);
		WarnError (errMsg);
		return;
	}

	for (i=0; i<theTrees.lLength; i++)
	{
		_SimpleList categVars;
		{
			_AVLList ca (&categVars);
			((_TheTree*)(LocateVar(theTrees(i))))->ScanForCVariables (ca);
			ca.ReorderList();
		}
		
		if (categVars.lLength)
		{
			long dep = 0;
			for (long k=categVars.lLength-1; k>=0; k--)
				SetNthBit(dep,indexCat.Find (categVars(k)));
			blockDependancies<<dep;
		}
		else
			blockDependancies<<0;
	}
	
	if (indexCat.lLength) // scan global variables upon which category variables depend 
						  // which may have been missed 
						  // indexInd will be sorted at this point
	{
		_SimpleList sortedCats (indexCat);
		sortedCats.Sort();
		if (sortedCats.CountCommonElements(indexInd,true))
		{
			_SimpleList newList;
			newList.Subtract (indexInd,sortedCats);
			indexInd.Clear();
			indexInd.Duplicate(&newList);
		}
		
		{
			_SimpleList   auxL;
			_AVLList      auxa (&auxL);
			
			for (i=0; i<indexCat.lLength; i++)
			{
				_CategoryVariable *theCV = (_CategoryVariable*)LocateVar(indexCat(i));
				theCV->ScanForGVariables (auxa);
			}
			
			auxa.ReorderList();
			if (auxL.lLength)
			{
				_SimpleList nl;
				nl.Union (indexInd,auxL);
				if (nl.lLength > indexInd.lLength)
				{
					indexInd.Clear();
					indexInd.Duplicate (&nl);
				}
			}
		}
	}
	
	_Parameter l = DEFAULTLOWERBOUND*(1.0-machineEps),
			   u = DEFAULTUPPERBOUND*(1.0-machineEps);
			   
	for (i=0;i<indexInd.lLength;i++)
	{
		_Variable *_cv = LocateVar(indexInd(i));
		if (_cv->GetLowerBound()<=l)
			_cv->SetBounds(DEFAULTPARAMETERLBOUND,_cv->GetUpperBound());
		if (_cv->GetUpperBound()>=u)
			_cv->SetBounds(_cv->GetLowerBound(),DEFAULTPARAMETERUBOUND);
	}
	for (i=0;i<indexDep.lLength;i++)
	{
		_Variable *_cv = LocateVar(indexDep(i));
		if (_cv->GetLowerBound()<=l)
			_cv->SetBounds(DEFAULTPARAMETERLBOUND,_cv->GetUpperBound());
		if (_cv->GetUpperBound()>=u)
			_cv->SetBounds(_cv->GetLowerBound(),DEFAULTPARAMETERUBOUND);
	}
}

//_______________________________________________________________________________________
void	_LikelihoodFunction::ScanAllVariablesOnPartition (_SimpleList& pidx, _SimpleList& iind, _SimpleList& idep, _SimpleList &icat)
{
	_SimpleList	  allVariables,
				  covCat,
				  cpCat;
				  
	{
		_AVLList avl (&allVariables);
		for (long i = 0; i < pidx.lLength; i++)
			LocateVar(theProbabilities(pidx(i)))->ScanForVariables (avl,true);
			
		if (computingTemplate)
			computingTemplate->ScanFForVariables (avl,true,false,true);

		avl.ReorderList();
	}

	if (templateKind<0)
		allVariables.Delete (allVariables.Find(-templateKind-1));
		
	{
		_AVLList iia (&iind),
				 iid (&idep);
	
		for (long i=0; i<allVariables.lLength; i++)
		{
			_Variable* theV = ((_Variable*)LocateVar(allVariables(i)));
			if (theV->IsCategory())
			{
				if (((_CategoryVariable*)theV)->IsUncorrelated())
				{
					if (((_CategoryVariable*)theV)->IsConstantOnPartition())
						icat << allVariables(i);
					else
						cpCat << allVariables(i);
				}
				else
					covCat<<allVariables(i);
				continue;
			}
			if (theV->IsIndependent())
				iia.Insert ((BaseRef)allVariables(i));
			else
				iid.Insert ((BaseRef)allVariables(i));
		}
		
		indexCat << cpCat;
		indexCat << covCat;
				
		for (long i2=0; i2<pidx.lLength; i2++)
		{
			_TheTree * cT = (_TheTree*)(LocateVar(theTrees.lData[pidx.lData[i2]]));
			cT->ScanForGVariables (iia, iid);
		}
		{
		for (long i=0; i<pidx.lLength; i++)
		{
			_TheTree * cT = (_TheTree*)(LocateVar(theTrees.lData[pidx.lData[i]]));
			cT->ScanForVariables (iia,iid);
			cT->ScanForDVariables (iid, iia);
		}
		}
		iia.ReorderList ();
		iid.ReorderList ();
	}


	for (long i=0; i<pidx.lLength; i++)
	{
		_SimpleList categVars;
		_AVLList ca (&categVars);
		((_TheTree*)(LocateVar(theTrees.lData[pidx.lData[i]])))->ScanForCVariables (ca);
		ca.ReorderList();
		
		if (categVars.lLength)
			for (long k=categVars.lLength-1; k>=0; k--)
			{
				long f = icat.Find (categVars(k));
				if (f==-1)
					icat<<categVars(k);
			}
	}
	
	if (icat.lLength)
	{
		for (long i=0; i<iind.lLength; i++)
		{
			long t = icat.Find (iind.lData[i]);
			if (t>=0)
			{
				iind.Delete(i);
				i--;
			}
		}
		{
			_SimpleList   auxL;
			_AVLList      auxa (&auxL);
			
			for (long i=0; i<icat.lLength; i++)
			{
				_CategoryVariable *theCV = (_CategoryVariable*)LocateVar(icat(i));
				theCV->ScanForGVariables (auxa);
			}
			
			auxa.ReorderList();
			if (auxL.lLength)
			{
				_SimpleList nl;
				nl.Union (iind,auxL);
				if (nl.lLength > iind.lLength)
				{
					iind.Clear();
					iind.Duplicate (&nl);
				}
			}
		}
	}
}

//_______________________________________________________________________________________
void	_LikelihoodFunction::UpdateIndependent (long index, bool purgeResults)
{
	long f = indexInd.Find (index);
	if (f!=-1)
	{
		indexInd.Delete(f);
		indexDep<<index;
		_SimpleList   newVars;
		{
			_AVLList  al (&newVars);
			LocateVar(index)->ScanForVariables(al,true);
			al.ReorderList();
		}
		for (f=0; f<newVars.lLength; f++)
		{
			_Variable* cv = LocateVar(newVars.lData[f]);
			if (cv->IsIndependent()&&(indexInd.Find(newVars.lData[f])==-1))
				indexInd << newVars.lData[f];
		}
		if (purgeResults)
			computationalResults.Clear();
	}

}

//_______________________________________________________________________________________
void	_LikelihoodFunction::UpdateDependent (long index)
{
	long f = indexDep.Find (index);
	if (f!=-1)
	{
		indexDep.Delete(f);
		indexInd<<index;
	}

}


//_______________________________________________________________________________________

void	_LikelihoodFunction::Cleanup (void)
{
	DeleteCaches();
}
	

//_______________________________________________________________________________________

void	_LikelihoodFunction::DeleteCaches (bool all)
{
	if (all)
	{
		DeleteObject (siteResults);   siteResults = nil;
		DeleteObject (bySiteResults); bySiteResults = nil;
	}
	
#ifdef	_SLKP_LFENGINE_REWRITE_
	conditionalTerminalNodeLikelihoodCaches.Clear();
	cachedBranches.Clear();
	siteCorrections.Clear();
	canUseReversibleSpeedups.Clear();
//	computedLocalUpdatePolicy.Clear();
//	treeTraversalMasks.Clear();
//	matricesToExponentiate.Clear();
//	overallScalingFactors.Clear();
//  localUpdatePolicy.Clear();
	
	if (conditionalInternalNodeLikelihoodCaches)
	{
		for (long k = 0; k < theTrees.lLength; k++)
			if (conditionalInternalNodeLikelihoodCaches[k]) delete (conditionalInternalNodeLikelihoodCaches[k]);
		delete (conditionalInternalNodeLikelihoodCaches);
		conditionalInternalNodeLikelihoodCaches = nil;
	}
	if (branchCaches)
	{
		for (long k = 0; k < theTrees.lLength; k++)
			if (branchCaches[k]) delete (branchCaches[k]);
		delete (branchCaches);
		branchCaches = nil;
	}
	if (conditionalTerminalNodeStateFlag)
	{
		for (long k = 0; k < theTrees.lLength; k++)
			if (conditionalTerminalNodeStateFlag[k]) delete (conditionalTerminalNodeStateFlag[k]);
		delete (conditionalTerminalNodeStateFlag);
		conditionalTerminalNodeStateFlag = nil;
	}
	if (siteScalingFactors)
	{
		for (long k = 0; k < theTrees.lLength; k++)
			if (siteScalingFactors[k]) delete (siteScalingFactors[k]);
		delete (siteScalingFactors);
		siteScalingFactors = nil;
	}
	#endif
}
	

//_______________________________________________________________________________________

void	_LikelihoodFunction::Setup (void)
{
	_Parameter kp		= 0.0;
	blockComputed		= -1;
	//RankVariables();
	checkParameter		(useFullMST,kp,1.0);
	
	if (kp>.5 && !mstCache)
		checkPointer (mstCache = new MSTCache);

	if (theTrees.lLength==optimalOrders.lLength) 
		//check to see if we need to recompute the
		// optimal summation order
	{
		checkParameter (keepOptimalOrder,kp,0.0);
		if (kp)
		{			
			for (int i=0; i<theTrees.lLength; i++)
			{
				_SimpleList*	s = (_SimpleList*)optimalOrders(i),
						   *	l = (_SimpleList*)leafSkips(i);
						   
				_DataSetFilter* df			= ((_DataSetFilter*)dataSetFilterList(theDataFilters(i)));
				_Matrix		  *glFreqs		= (_Matrix*)LocateVar(theProbabilities.lData[i])->GetValue();
				_TheTree	  *t			= ((_TheTree*)LocateVar(theTrees.lData[i]));
				
				t->InitializeTreeFrequencies (glFreqs, true);
				if (s->lLength!=df->NumberDistinctSites())
				{
					s->Clear();
					l->Clear();
					OptimalOrder (i,*s);
					df->MatchStartNEnd (*s,*l);
				}
			}
			return;
		}
	}
			
	optimalOrders.Clear();
	leafSkips.Clear();
	
#ifdef	_SLKP_LFENGINE_REWRITE_
	treeTraversalMasks.Clear();
	canUseReversibleSpeedups.Clear();
	_SimpleList alreadyDoneModelsL;
	_AVLListX   alreadyDoneModels (&alreadyDoneModelsL);
#endif
	
	for (long i=0; i<theTrees.lLength; i++)
	{
		_Matrix			*glFreqs = (_Matrix*)LocateVar(theProbabilities.lData[i])->GetValue();
		_DataSetFilter* df		= ((_DataSetFilter*)dataSetFilterList(theDataFilters(i)));
		_TheTree		*t		 = ((_TheTree*)LocateVar(theTrees.lData[i]));
		t->InitializeTreeFrequencies (glFreqs, true);
		_SimpleList		   *s = new _SimpleList,
						   *l = new _SimpleList;
		
#ifdef	_SLKP_LFENGINE_REWRITE_
		treeTraversalMasks.AppendNewInstance(new _SimpleList (t->GetINodeCount() * df->NumberDistinctSites() / _HY_BITMASK_WIDTH_ + 1,0,0));			
#endif
		
		OptimalOrder	  (i,*s);
		df->MatchStartNEnd(*s,*l);
		optimalOrders.AppendNewInstance(s);
		leafSkips.AppendNewInstance(l);
#ifdef	_SLKP_LFENGINE_REWRITE_
		_SimpleList treeModels;
		t->CompileListOfModels (treeModels);
		bool isReversiblePartition = true;
		for (long m = 0; m < treeModels.lLength && isReversiblePartition; m++)
		{
			long alreadyDone = alreadyDoneModels.Find ((BaseRef)treeModels.lData[m]);
			if (alreadyDone>=0)
				alreadyDone = alreadyDoneModels.GetXtra (alreadyDone);
			else
			{
				alreadyDone = IsModelReversible (treeModels.lData[m]);
				alreadyDoneModels.Insert ((BaseRef)treeModels.lData[m], alreadyDone);
			}
			isReversiblePartition = isReversiblePartition && alreadyDone;
		}
		ReportWarning (_String ("Partition ") & i & " reversible model flag computed as " & (long)isReversiblePartition);
		canUseReversibleSpeedups << isReversiblePartition;
		
#endif
	}
	
}

//_______________________________________________________________________________________

#ifdef 		__MP__
	void*	ThreadReleafFunctionAA (void* arg)
	{
		ThreadReleafTask * theTask = (ThreadReleafTask *)arg;
		
		theTask->siteList += theTask->startAt;
		theTask->leafList += theTask->startAt;
		
		theTask->result =  myLog(theTask->t->ThreadReleafTreeCharCache (theTask->df, *(theTask->siteList), -1, 0, 
							   theTask->t->GetLeafCount()-1,theTask->startAt,theTask->threadIndex));
							
		theTask->result*=(_Parameter)theTask->freqs->lData[*(theTask->siteList++)];
		for (long i = theTask->startAt+1; i<theTask->endAt; i++,theTask->siteList++,theTask->leafList++)
			theTask->result += myLog(theTask->t->ThreadReleafTreeCharCache (theTask->df, *theTask->siteList, 
				   *(theTask->siteList-1), (*(theTask->leafList))&0x0000ffff, (*(theTask->leafList))>>16,i,theTask->threadIndex))
				   *(_Parameter)theTask->freqs->lData[*(theTask->siteList)];
		return nil;
	}
	
	void*	ThreadReleafFunctionCodon (void* arg)
	{
		ThreadReleafTask * theTask = (ThreadReleafTask *)arg;
		
		theTask->siteList += theTask->startAt;
		theTask->leafList += theTask->startAt;
		
		theTask->result =  myLog(theTask->t->ThreadReleafTreeCache (theTask->df, *(theTask->siteList), -1, 0, 
							   theTask->t->GetLeafCount()-1,theTask->startAt,theTask->threadIndex));
							
		theTask->result*=(_Parameter)theTask->freqs->lData[*(theTask->siteList++)];
		for (long i = theTask->startAt+1; i<theTask->endAt; i++,theTask->siteList++,theTask->leafList++)
			theTask->result += myLog(theTask->t->ThreadReleafTreeCache (theTask->df, *theTask->siteList, 
				   *(theTask->siteList-1), (*(theTask->leafList))&0x0000ffff, (*(theTask->leafList))>>16,i,theTask->threadIndex))
				   *(_Parameter)theTask->freqs->lData[*(theTask->siteList)];
		return nil;
	}

	void*	ThreadReleafFunctionNuc (void* arg)
	{
		ThreadReleafTask * theTask = (ThreadReleafTask *)arg;
		
		theTask->siteList += theTask->startAt;
		theTask->leafList += theTask->startAt;
		
		theTask->result =  myLog(theTask->t->ThreadReleafTreeChar4 (theTask->df, *(theTask->siteList), -1, 0, 
							   theTask->t->GetLeafCount()-1,theTask->startAt,theTask->threadIndex));
							
		theTask->result*=(_Parameter)theTask->freqs->lData[*(theTask->siteList++)];
		for (long i = theTask->startAt+1; i<theTask->endAt; i++,theTask->siteList++,theTask->leafList++)
			theTask->result += myLog(theTask->t->ThreadReleafTreeChar4 (theTask->df, *theTask->siteList, 
				   *(theTask->siteList-1), (*(theTask->leafList))&0x0000ffff, (*(theTask->leafList))>>16,i,theTask->threadIndex))
				   *(_Parameter)theTask->freqs->lData[*(theTask->siteList)];
		return nil;
	}

	void*	ThreadReleafFunctionMNuc (void* arg)
	{
		ThreadReleafTask * theTask = (ThreadReleafTask *)arg;
		
		theTask->siteList += theTask->startAt;
		theTask->leafList += theTask->startAt;
		
		theTask->resultsStorage[*(theTask->siteList)] =  theTask->t->ThreadReleafTreeChar4 (theTask->df, *(theTask->siteList), -1, 0, 
							         theTask->t->GetLeafCount()-1,theTask->startAt*theTask->t->categoryCount+categID,theTask->threadIndex);
							
		theTask->siteList++;
		for (long i = theTask->startAt+1; i<theTask->endAt; i++,theTask->siteList++,theTask->leafList++)
			theTask->resultsStorage[*(theTask->siteList)] = theTask->t->ThreadReleafTreeChar4 (theTask->df, *theTask->siteList, 
				   *(theTask->siteList-1), (*(theTask->leafList))&0x0000ffff, (*(theTask->leafList))>>16,
				   i*theTask->t->categoryCount+categID,theTask->threadIndex);
		return nil;
	}

	void*	ThreadReleafFunctionMAA (void* arg)
	{
		ThreadReleafTask * theTask = (ThreadReleafTask *)arg;
		
		theTask->siteList += theTask->startAt;
		theTask->leafList += theTask->startAt;
		
		theTask->resultsStorage[*(theTask->siteList)] =  theTask->t->ThreadReleafTreeCharCache (theTask->df, *(theTask->siteList), -1, 0, 
							         theTask->t->GetLeafCount()-1,theTask->startAt*theTask->t->categoryCount+categID,theTask->threadIndex);
							
		theTask->siteList++;
		for (long i = theTask->startAt+1; i<theTask->endAt; i++,theTask->siteList++,theTask->leafList++)
			theTask->resultsStorage[*(theTask->siteList)] = theTask->t->ThreadReleafTreeCharCache (theTask->df, *theTask->siteList, 
				   *(theTask->siteList-1), (*(theTask->leafList))&0x0000ffff, (*(theTask->leafList))>>16,
				   i*theTask->t->categoryCount+categID,theTask->threadIndex);
		return nil;
	}

	void*	ThreadReleafFunctionMCodon (void* arg)
	{
		ThreadReleafTask * theTask = (ThreadReleafTask *)arg;
		
		theTask->siteList += theTask->startAt;
		theTask->leafList += theTask->startAt;
		
		theTask->resultsStorage[*(theTask->siteList)] =  theTask->t->ThreadReleafTreeCache (theTask->df, *(theTask->siteList), -1, 0, 
							         theTask->t->GetLeafCount()-1,theTask->startAt*theTask->t->categoryCount+categID,theTask->threadIndex);
							
		theTask->siteList++;
		for (long i = theTask->startAt+1; i<theTask->endAt; i++,theTask->siteList++,theTask->leafList++)
			theTask->resultsStorage[*(theTask->siteList)] = theTask->t->ThreadReleafTreeCache (theTask->df, *theTask->siteList, 
				   *(theTask->siteList-1), (*(theTask->leafList))&0x0000ffff, (*(theTask->leafList))>>16,
				   i*theTask->t->categoryCount+categID,theTask->threadIndex);
		return nil;
	}

#endif


//_______________________________________________________________________________________

void	StashMSTStuff (_Parameter** mstCI, long **mstSI, char ** mstNI, long index, long index2)
{
	_Parameter* temp = mstCI[index];
	mstCI[index] = mstCI[index2];
	mstCI[index2] = temp;
	
	long	  * t2   = mstSI[index];
	mstSI[index] = mstSI[index2];
	mstSI[index2] = t2;

	char	  * t3   = mstNI[index];
	mstNI[index] = mstNI[index2];
	mstNI[index2] = t3;
}

//_______________________________________________________________________________________

void	_LikelihoodFunction::ComputeBlockInt1 (long index, _Parameter& res, _TheTree* t, _DataSetFilter* df, char type)
{
	long i;

	_SimpleList	*sl = (_SimpleList*)optimalOrders.lData[index],
				*ls = (_SimpleList*)leafSkips.lData[index];

	_SimpleList* co = (_SimpleList*)mstCache->computingOrder.lData[index],
			   * so = (_SimpleList*)mstCache->storageOrder.lData[index],
			   * ro = (_SimpleList*)mstCache->referenceOrder.lData[index],
			   * po = (_SimpleList*)mstCache->parentOrder.lData[index];
			   
	_Parameter*  stashMLC = t->marginalLikelihoodCache,
			  ** mstMLC	  = (_Parameter**)mstCache->resultCache.lData[index];
			  
	long	  ** mstNC	  = (long **)mstCache->statesCache.lData[index];
	char	  ** mstNMC	  = (char **)mstCache->statesNCache.lData[index];
			  
			   
	long		copySize2  = t->GetLeafCount()+t->GetINodeCount(),
				copySize1  = copySize2 * t->GetCodeBase() * sizeof (_Parameter),
				copySize3  = t->GetINodeCount(),
				soI	   	   = so->lData[0],
				cS		   = mstCache->cacheSize.lData[index],
			 	* leafList = ls->lData;				   
				
	
	copySize2 *= sizeof (long);

	
	if (type==0)
		res =  myLog(t->ReleafTreeAndCheckChar4 (df, co->lData[0], true));
	else
		if (type == 1)
			res =  myLog(t->ReleafTreeAndCheckChar (df, co->lData[0], true));		
		else
			res =  myLog(t->ReleafTreeAndCheck (df, co->lData[0], true));

	res *= (_Parameter)df->theFrequencies.lData[co->lData[0]];
	
	mstMLC [cS] = t->marginalLikelihoodCache;
	mstNC  [cS] = t->nodeStates;
	mstNMC [cS] = t->nodeMarkers;

	if (soI>=0)
		StashMSTStuff (mstMLC, mstNC, mstNMC, soI, cS);																
											
	for (i = 1; i<sl->lLength; i++, leafList++)
	{
		#if !defined __UNIX__ || defined __HEADLESS__
			siteEvalCount ++;
			if (divideBy && (siteEvalCount%divideBy == 0) || siteEvalCount%1000 == 0)
				yieldCPUTime();
		#endif
		
		soI = ro->lData[i];
		if (soI>=0)
		{
			if (soI!=so->lData[i])
			{
				t->marginalLikelihoodCache = mstMLC [cS];
				t->nodeStates			   = mstNC  [cS];
				t->nodeMarkers			   = mstNMC [cS];
				
				_Parameter *p0 = t->marginalLikelihoodCache,
						   *p1 = mstMLC[soI],
						   *pS = t->marginalLikelihoodCache+copySize1/sizeof (_Parameter);
						   
				for (;p0!=pS;)
					*(p0++)=*(p1++);

				long	   *l0 = t->nodeStates,
						   *l1 = mstNC[soI],
						   *lS = t->nodeStates+copySize2/sizeof (long);
						   
				for (;l0!=lS;)
					*(l0++)=*(l1++);
						   
				char	   *c0 = t->nodeMarkers,
						   *c1 = mstNMC[soI],
						   *cS = t->nodeMarkers+copySize3;
						   
				for (;c0!=cS;)
					*(c0++)=*(c1++);

				//memcpy (t->marginalLikelihoodCache,mstMLC[soI],copySize1);
				//memcpy (t->nodeStates, mstNC[soI], copySize2);
				//memcpy (t->nodeMarkers,mstNMC[soI], copySize3);
			}
			else
			{
				t->marginalLikelihoodCache = mstMLC [soI];
				t->nodeStates			   = mstNC  [soI];
				t->nodeMarkers			   = mstNMC [soI];										
			}
		}
		
		if (type==0)
			res += myLog(t->ThreadReleafTreeChar4 (df, co->lData[i], 
				   po->lData[i], (*leafList)&0x0000ffff, (*leafList)>>16,i))
				   *(_Parameter)df->theFrequencies.lData[co->lData[i]];
		else
			if (type == 1)
				res += myLog(t->ThreadReleafTreeCharCache (df, co->lData[i], 
					   po->lData[i], (*leafList)&0x0000ffff, (*leafList)>>16,i))
					   *(_Parameter)df->theFrequencies.lData[co->lData[i]];
			else
				res += myLog(t->ThreadReleafTreeCache (df, co->lData[i], 
					   po->lData[i], (*leafList)&0x0000ffff, (*leafList)>>16,i))
					   *(_Parameter)df->theFrequencies.lData[co->lData[i]];
			   
			   
		if ((soI!=so->lData[i])&&(so->lData[i]>=0))
			StashMSTStuff (mstMLC, mstNC, mstNMC, so->lData[i], cS);							
	}		
	
	for (i = 0; i<cS; i++)
		if (mstMLC[i] == stashMLC)
			break;
			
	if (i<cS)
	{
		StashMSTStuff (mstMLC, mstNC, mstNMC, i, cS);							
		t->marginalLikelihoodCache = mstMLC [cS];
		t->nodeStates			   = mstNC  [cS];
		t->nodeMarkers			   = mstNMC [cS];
	}
}

//_______________________________________________________________________________________

_Parameter	_LikelihoodFunction::ComputeBlock (long index, _Parameter* siteRes)
// compute likelihood over block index i
{
	blockComputed				= index;
	// set up global matrix frequencies
	_Matrix 			*glFreqs = (_Matrix*)LocateVar(theProbabilities.lData[index])->GetValue();
	SetGlobalFrequencyMatrix (glFreqs);
	
	_SimpleList 		*sl = (_SimpleList*)optimalOrders.lData[index],
						*ls = (_SimpleList*)leafSkips.lData[index];
						
	_TheTree 			*t = ((_TheTree*)LocateVar(theTrees.lData[index]));
	t->InitializeTreeFrequencies (GlobalFrequenciesMatrix);
	
	_DataSetFilter		*df = ((_DataSetFilter*)dataSetFilterList(theDataFilters.lData[index]));
	
	if (!forceRecomputation)
	{
		forceRecomputation = glFreqs->HasChanged();
//		if ((!forceRecomputation)&&computingTemplate)
//			forceRecomputation = computingTemplate->HasChanged();
	}
		
	usedCachedResults = false;
	
	if (computingTemplate&&templateKind)
	{
		if (!(forceRecomputation||t->HasChanged()||(!siteArrayPopulated)))
		{		
			usedCachedResults = true;
			return			  -1e300;
		}	
	}
	else
	{
		if ((!forceRecomputation)&&computationalResults.lLength&&(computationalResults.lLength==optimalOrders.lLength)&&(siteRes)&&(!t->HasChanged()))
		{		
			usedCachedResults = true;
			return ((_Constant*)computationalResults.lData[index])->Value();
		}
	}
	
#ifdef _SLKP_LFENGINE_REWRITE_
	if (conditionalInternalNodeLikelihoodCaches && conditionalInternalNodeLikelihoodCaches[index])
	{
		_SimpleList changedBranches, *branches;
		_List		changedModels,   *matrices;
		long		doCachedComp	 = 0;	  // whether or not to use a cached branch calculation when only one
											  // local tree parameter is being adjusted at a time
		
		long		catID			 = siteRes?categID:-1,
					ciid			 = MAX(0,categID),
					*cbid			 = &(((_SimpleList*)cachedBranches(index))->lData[ciid]);
		
		if (computedLocalUpdatePolicy.lLength)
		{
			branches = (_SimpleList*)(*((_List*)localUpdatePolicy(index)))(ciid);
			matrices = (_List*)      (*((_List*)matricesToExponentiate(index)))(ciid) ;
			
			long nodeID = ((_SimpleList*)computedLocalUpdatePolicy(index))->lData[ciid];
			if (nodeID == 0 || nodeID == 1)
			{
				long snID = -1;
				
				if (nodeID == 1)
				{
					if (matrices->lLength == 2)
					{
						branches->Clear();
						matrices->Clear();
						
						snID = t->DetermineNodesForUpdate		   (*branches, matrices,catID,*cbid);			
					}
				}
				else
					snID = t->DetermineNodesForUpdate		   (*branches, matrices,catID,*cbid);
				
				*cbid = -1;
				if (snID >= 0 && canUseReversibleSpeedups.lData[index])
				{
					((_SimpleList*)computedLocalUpdatePolicy(index))->lData[ciid] = snID+3;
					doCachedComp = -snID-1;
					//doCachedComp = 0;
				}
				else
					((_SimpleList*)computedLocalUpdatePolicy(index))->lData[ciid] = nodeID + 1;
			}
			else
				doCachedComp = nodeID;
		}
		else
		{
			t->DetermineNodesForUpdate		   (changedBranches,&changedModels,catID,*cbid);
			*cbid						= -1;
			branches					= &changedBranches;
			matrices					= &changedModels;
		}
		
		if (matrices->lLength)
			t->ExponentiateMatrices(*matrices, GetThreadCount(),catID);
	
		//printf ("%d %d\n",likeFuncEvalCallCount, matrixExpCount);
#if !defined __UNIX__ || defined __HEADLESS__
		if (divideBy && (likeFuncEvalCallCount % divideBy == 0))
			yieldCPUTime();
	#endif


		long blockID    = df->NumberDistinctSites()*t->GetINodeCount();
		
		_SimpleList* tcc = (_SimpleList*)treeTraversalMasks(index);

		_Parameter*inc  = (categID<1)?conditionalInternalNodeLikelihoodCaches[index]:
									  conditionalInternalNodeLikelihoodCaches[index] + categID*df->GetDimension()*blockID,
				  *ssf  = (categID<1)?siteScalingFactors[index]: siteScalingFactors[index] + categID*blockID,
				  *bc   = (categID<1)?branchCaches[index]: (branchCaches[index] + categID*df->NumberDistinctSites()*df->GetDimension()*2);
		
		long *scc = nil;
		
		if (siteRes)
			scc = ((_SimpleList*)siteCorrections(index))->lData + ((categID<1)?0:df->NumberDistinctSites()*categID);
		
		_Parameter sum  = 0.;
		if (doCachedComp >= 3)
		{
			sum = t->ComputeLLWithBranchCache (*sl,
												doCachedComp-3,
												bc,
												df,
												0,
												df->NumberDistinctSites (),
												catID,
												siteRes)
												-	_logLFScaler * overallScalingFactors[index];
			return sum;
		}

		long np = 1, sitesPerP = df->NumberDistinctSites() / np + 1;
#ifdef _OPENMP
		np			 = MIN(GetThreadCount(),omp_get_max_threads());
		sitesPerP    = df->NumberDistinctSites() / np + 1;		
#endif		
			
			
		#pragma omp  parallel for default(shared) schedule(static,1) private(blockID) num_threads (np) reduction(+:sum) if (np>1)
		for (blockID = 0; blockID < np; blockID ++)
		{
			sum += t->ComputeTreeBlockByBranch (*sl, 
												*branches, 
												tcc,
												df, 
												inc,
												conditionalTerminalNodeStateFlag[index],
												ssf,
												(_GrowingVector*)conditionalTerminalNodeLikelihoodCaches(index),
												overallScalingFactors[index],
												blockID * sitesPerP,
												(1+blockID) * sitesPerP,
												catID,
												siteRes,
												scc);
		}
		
		sum -= _logLFScaler * overallScalingFactors[index];


		if (doCachedComp < 0)
		{
			//printf ("Cache check in %d %d\n", doCachedComp, overallScalingFactors[index]);
			doCachedComp = -doCachedComp-1;
			//printf ("Set up %d\n", doCachedComp);
			*cbid = doCachedComp;
#pragma omp  parallel for default(shared) schedule(static,1) private(blockID) num_threads (np) if (np>1)
			for (blockID = 0; blockID < np; blockID ++)
			{
				t->ComputeBranchCache (*sl,doCachedComp, bc, inc, df, 
									   conditionalTerminalNodeStateFlag[index], 
									   ssf, 
									   scc,
									   (_GrowingVector*)conditionalTerminalNodeLikelihoodCaches(index),
									   overallScalingFactors[index],
									   blockID * sitesPerP,
									   (1+blockID) * sitesPerP,
									   catID,tcc);
			}
						
			/*_Parameter check2 =  t->ComputeLLWithBranchCache (*sl,
																   doCachedComp,
																   bc,
																   df,
																   0,
																   df->NumberDistinctSites (),
																   catID,
																   siteRes)
			-	_logLFScaler * overallScalingFactors[index];
			
			if (fabs(check2-sum)>0.1)
				printf ("Cache check failed LF-eval %d node %d scaled difference %g %d\n", likeFuncEvalCallCount, doCachedComp, fabs(sum-check2)/_logLFScaler, overallScalingFactors[index]);
			*/
			 //else
			//	printf ("Cache check OK\n");
			//printf ("Cache compute %d\n", catID);
			
		}
		return sum;
	}
#endif	
	
	if (1)
	{
		long * siteList = sl->lData, 
			 * leafList = ls->lData,				   
			   filterUnit = df->GetUnitLength(),
			   i;
			   
		#ifdef __MP__
			#if !defined __UNIX__ || defined __HEADLESS__
				bool	yieldTime = false;
			#endif
			long tstep = sl->lLength/systemCPUCount;
			if ((t->HasCache())&&(systemCPUCount>1)&&(tstep>5))
			{
				theThreads = new pthread_t [systemCPUCount-1];
				theTasks = new ThreadReleafTask [systemCPUCount-1];
				for (long tc = 0; tc<systemCPUCount-1; tc++)
				{
					theTasks[tc].t = t;
					theTasks[tc].siteList = siteList;
					theTasks[tc].leafList = leafList;
					theTasks[tc].startAt = tstep*(tc+1);
					theTasks[tc].freqs = &df->theFrequencies;
					theTasks[tc].endAt = tstep*(tc+2);
					theTasks[tc].df = df;
					if (tc == systemCPUCount-2)
						theTasks[tc].endAt = sl->lLength;
					theTasks[tc].threadIndex = tc+1;
				}
			}
			else
				theTasks = nil;
		#endif
		
		if (!siteRes)
		{
			_Parameter res = 0.0;
			if ( filterUnit==1 && t->GetCodeBase()==4)
			{
				if (t->HasCache())
				{
					#ifdef __MP__
					if (theTasks)
					{
						res = myLog(t->ReleafTreeAndCheckChar4 (df, *siteList, true));
																								 		
						for (long tc = 0; tc<systemCPUCount-1; tc++)
						{
							#ifdef __MACHACKMP__
							if ( pthread_create( theThreads+tc, NULL, ThreadReleafFunctionNucHook,(void*)(theTasks+tc)))
							#else
							if ( pthread_create( theThreads+tc, NULL, ThreadReleafFunctionNuc,(void*)(theTasks+tc)))							
							#endif
							{
								FlagError("Failed to initialize a POSIX thread in ComputeBlock()");
								exit(1);
							}		
						}
												
						res*=(_Parameter)df->theFrequencies.lData[*(siteList++)];
						
						for (i = 1; i<tstep; i++,siteList++,leafList++)
						{
							#ifndef __UNIX__
								siteEvalCount ++;
								if (divideBy && (siteEvalCount%divideBy == 0) || siteEvalCount%1000 == 0)
									yieldTime = true;
							#endif
							res += myLog(t->ThreadReleafTreeChar4 (df, *siteList, 
								   *(siteList-1), (*leafList)&0x0000ffff, (*leafList)>>16,i))
								   *(_Parameter)df->theFrequencies.lData[*siteList];
						}
					}
					else
					{
					#endif					
						if (mstCache&&(mstCache->resultCache.lLength>index)&&(mstCache->cacheSize.lData[index]>0))
						{
							ComputeBlockInt1 (index, res, t, df, 0);
						}
						else
						{
							res=myLog(t->ReleafTreeAndCheckChar4 (df, *siteList, true));
							
							res*=(_Parameter)df->theFrequencies.lData[*(siteList++)];
							for (i = 1; i<sl->lLength; i++,siteList++,leafList++)
							{
								#if !defined __UNIX__ || defined __HEADLESS__
									siteEvalCount ++;
									if (divideBy && (siteEvalCount%divideBy == 0) || siteEvalCount%1000 == 0)
									#ifdef __MP__
										yieldTime = true;
									#else
										yieldCPUTime();
									#endif
								#endif
								//res += myLog(t->ReleafTreeChar4 (df, *siteList, 
								res += myLog(t->ThreadReleafTreeChar4 (df, *siteList, 
									   *(siteList-1), (*leafList)&0x0000ffff, (*leafList)>>16,i))
									   *(_Parameter)df->theFrequencies.lData[*siteList];
						}
							   
					}
					#ifdef __MP__
					}
					#endif
				}
				else // no cache
					if (t->IsDegenerate())
					{
						res=myLog(t->ReleafTreeAndCheckChar4 (df, *siteList, false));
						res*=(_Parameter)df->theFrequencies.lData[*(siteList++)];
						for (i = 1; i<sl->lLength; i++,siteList++)
						{
							#if !defined __UNIX__ || defined __HEADLESS__
								siteEvalCount ++;
								if (divideBy && (siteEvalCount%divideBy == 0) || siteEvalCount%1000 == 0)
									yieldCPUTime();
							#endif
							res += myLog(t->ReleafTreeChar4Degenerate (df, *siteList))
								   *(_Parameter)df->theFrequencies.lData[*siteList];
						}											
					}
					else
					{
						res=myLog(t->ReleafTreeAndCheckChar4 (df, *siteList, false));
						res*=(_Parameter)df->theFrequencies.lData[*(siteList++)];
						if (df->IsNormalFilter())
							for (i = 1; i<sl->lLength; i++,siteList++,leafList++)
							{
								#if !defined __UNIX__ || defined __HEADLESS__
									siteEvalCount ++;
									if (divideBy && (siteEvalCount%divideBy == 0) || siteEvalCount%1000 == 0)
										yieldCPUTime();
								#endif
								res += myLog(t->ReleafTreeChar4 (df, *siteList, 
									   *(siteList-1), (*leafList)&0x0000ffff, (*leafList)>>16))
									   *(_Parameter)df->theFrequencies.lData[*siteList];
							}	
						else
						{
                            _Parameter accumulator = 1.0;
                            for (i = 1; i<sl->lLength; i++,siteList++)
                            {
								//res += myLog(t->ReleafTreeCharNumFilter4Tree3 ((_DataSetFilterNumeric*)df, *siteList)) \
								//  *(_Parameter)df->theFrequencies.lData[*siteList];                                    \
                                                                                                                        
                                long thisSiteFreq = df->theFrequencies.lData[*siteList];
                                _Parameter thisSiteScore = t->ReleafTreeCharNumFilter4Tree3 ((_DataSetFilterNumeric*)df, *siteList);
	                            if (thisSiteScore>0.0)
									for (; thisSiteFreq; thisSiteFreq--)
									{
										_Parameter testAcc = accumulator * thisSiteScore;
										if (testAcc>1.e-300)
											accumulator = testAcc;
										else
										{
											res += log (accumulator);
											accumulator = thisSiteScore;
										}
									}
	                             else
                            		res -= thisSiteFreq*1000000.;
                             }
                           	 res += log (accumulator);
                          }							
					}
			}
			else
			{
				if (filterUnit>1)
				{
					if (t->HasCache())
					{
						#ifdef __MP__
						if (theTasks)
						{
							res = myLog(t->ReleafTreeAndCheck (df, *siteList, true));
							
							for (long tc = 0; tc<systemCPUCount-1; tc++)
							{
								#ifdef __MACHACKMP__
								if ( pthread_create( theThreads+tc, NULL, ThreadReleafFunctionCodonHook,(void*)(theTasks+tc)))
								#else
								if ( pthread_create( theThreads+tc, NULL, ThreadReleafFunctionCodon,(void*)(theTasks+tc)))
								#endif
								{
									FlagError("Failed to initialize a POSIX thread in ComputeBlock()");
									exit(1);
								}		
							}
													
							res*=(_Parameter)df->theFrequencies.lData[*(siteList++)];
							
							for (i = 1; i<tstep; i++,siteList++,leafList++)
							{
								#ifndef __UNIX__
									siteEvalCount ++;
									if (divideBy && (siteEvalCount%divideBy == 0) || siteEvalCount%1000 == 0)
										yieldTime = true;
								#endif
								res += myLog(t->ThreadReleafTreeCache (df, *siteList, 
									   *(siteList-1), (*leafList)&0x0000ffff, (*leafList)>>16,i))
									   *(_Parameter)df->theFrequencies.lData[*siteList];
							}
						}
						else
						{
						#endif
							if (mstCache&&(mstCache->resultCache.lLength>index)&&(mstCache->cacheSize.lData[index]>0))
							{
								ComputeBlockInt1 (index, res, t, df, 2);
							}
							else
							{
								res=myLog(t->ReleafTreeAndCheck (df, *siteList, true));
								res*=(_Parameter)df->theFrequencies.lData[*(siteList++)];
								
								for (i = 1; i<sl->lLength; i++,siteList++,leafList++)
								{
									#if !defined __UNIX__ || defined __HEADLESS__
										siteEvalCount ++;
										if (!(siteEvalCount%divideBy)) 
											#ifdef __MP__
												yieldTime = true;
											#else
												yieldCPUTime();
											#endif
										
									#endif
									//res += myLog(t->ReleafTreeCache (df, *siteList, 
									res += myLog(t->ThreadReleafTreeCache (df, *siteList, 
										   *(siteList-1), (*leafList)&0x0000ffff, (*leafList)>>16,i))
										   *(_Parameter)df->theFrequencies.lData[*siteList];
									/*res += myLog(t->ReleafTree (df, *siteList, 
										   *(siteList-1), (*leafList)&0x0000ffff, (*leafList)>>16))
										   *(_Parameter)df->theFrequencies.lData[*siteList];*/
										   
										   
								}
							}
						#ifdef __MP__
						}
						#endif
					} // no cache
					else
						if (t->IsDegenerate())
						{
							res=myLog(t->ReleafTreeAndCheck (df, *siteList, false));
							res*=(_Parameter)df->theFrequencies.lData[*(siteList++)];
							for (i = 1; i<sl->lLength; i++,siteList++)
							{
								#if !defined __UNIX__ || defined __HEADLESS__
									siteEvalCount ++;
									if (divideBy && (siteEvalCount%divideBy == 0) || siteEvalCount%1000 == 0)
										yieldCPUTime();
								#endif
								res += myLog(t->ReleafTreeDegenerate (df, *siteList))
									   *(_Parameter)df->theFrequencies.lData[*siteList];
							}											
						}
						else
						{
							res=myLog(t->ReleafTreeAndCheck (df, *siteList, false));
							res*=(_Parameter)df->theFrequencies.lData[*(siteList++)];
							for (i = 1; i<sl->lLength; i++,siteList++,leafList++)
							{
								#if !defined __UNIX__ || defined __HEADLESS__
									siteEvalCount ++;
									if (divideBy && (siteEvalCount%divideBy == 0) || siteEvalCount%1000 == 0)
										yieldCPUTime();
								#endif
								res += myLog(t->ReleafTree (df, *siteList, 
									   *(siteList-1), (*leafList)&0x0000ffff, (*leafList)>>16))
									   *(_Parameter)df->theFrequencies.lData[*siteList];
							}	
						}
				}
				else
				{
					if (t->HasCache())
					{
						#ifdef __MP__
						long tstep = sl->lLength/systemCPUCount;
						if (theTasks)
						{
							
							res=myLog(t->ReleafTreeAndCheckChar (df, *siteList, true));
							for (long tc = 0; tc<systemCPUCount-1; tc++)
								#ifdef __MACHACKMP__
								if ( pthread_create( theThreads+tc, NULL, ThreadReleafFunctionAAHook, (void*)(theTasks+tc)) ) 
								#else
								if ( pthread_create( theThreads+tc, NULL, ThreadReleafFunctionAA, (void*)(theTasks+tc)) ) 
								#endif
								{
									FlagError("Failed to initialize a POSIX thread in ComputeBlock()");
									exit(1);
								}		
													
							res*=(_Parameter)df->theFrequencies.lData[*(siteList++)];
							for (i = 1; i<tstep; i++,siteList++,leafList++)
							{
								#ifndef __UNIX__
									siteEvalCount ++;
									if (divideBy && (siteEvalCount%divideBy == 0) || siteEvalCount%1000 == 0)
										yieldTime = true;
								#endif
								res += myLog(t->ThreadReleafTreeCharCache (df, *siteList, 
									   *(siteList-1), (*leafList)&0x0000ffff, (*leafList)>>16,i))
									   *(_Parameter)df->theFrequencies.lData[*siteList];
							}
						}
						else
						{
						#endif
							if (mstCache&&(mstCache->resultCache.lLength>index)&&(mstCache->cacheSize.lData[index]>0))
							{
								ComputeBlockInt1 (index, res, t, df, 1);
							}
							else
							{
								res=myLog(t->ReleafTreeAndCheckChar (df, *siteList, true));
								res*=(_Parameter)df->theFrequencies.lData[*(siteList++)];
			   					for (i = 1; i<sl->lLength; i++,siteList++,leafList++)
								{
									#if !defined __UNIX__ || defined __HEADLESS__
										siteEvalCount ++;
										if (divideBy && (siteEvalCount%divideBy == 0) || siteEvalCount%1000 == 0)
											#ifdef __MP__
												yieldTime = true;
											#else
												yieldCPUTime();
											#endif

									#endif
									res += myLog(t->ThreadReleafTreeCharCache (df, *siteList, 
										   *(siteList-1), (*leafList)&0x0000ffff, (*leafList)>>16,i))
										   *(_Parameter)df->theFrequencies.lData[*siteList];
								}
							}
							
						#ifdef __MP__
						}
						#endif
					}
					else
						if (t->IsDegenerate())
						{
							res=myLog(t->ReleafTreeAndCheck (df, *siteList, false));
							res*=(_Parameter)df->theFrequencies.lData[*(siteList++)];
							for (i = 1; i<sl->lLength; i++,siteList++)
							{
								#if !defined __UNIX__ || defined __HEADLESS__
									siteEvalCount ++;
									if (divideBy && (siteEvalCount%divideBy == 0) || siteEvalCount%1000 == 0)
										yieldCPUTime();
								#endif
								res += myLog(t->ReleafTreeCharDegenerate (df, *siteList))
									   *(_Parameter)df->theFrequencies.lData[*siteList];
							}											
						}
						else
						{
							res=myLog(t->ReleafTreeAndCheckChar (df, *siteList, false));
							res*=(_Parameter)df->theFrequencies.lData[*(siteList++)];
							for (i = 1; i<sl->lLength; i++,siteList++,leafList++)
							{
								#if !defined __UNIX__ || defined __HEADLESS__
									siteEvalCount ++;
									if (divideBy && (siteEvalCount%divideBy == 0) || siteEvalCount%1000 == 0)
										yieldCPUTime();
								#endif
								res += myLog(t->ReleafTreeChar (df, *siteList, 
									   *(siteList-1), (*leafList)&0x0000ffff, (*leafList)>>16))
									   *(_Parameter)df->theFrequencies.lData[*siteList];
							}					
						}
				}
			}
			#ifdef __MP__
				if (theTasks)
				{
					for (long wc = 0; wc<systemCPUCount-1; wc++)
					{
						if ( pthread_join ( theThreads[wc], NULL ) ) 
						{
							FlagError("Failed to join POSIX threads in ComputeBlock()");
							exit(1);
						}
						res += theTasks[wc].result;
					}
					delete theTasks;
					delete theThreads;				
				}
			#endif 
			
			#if USE_SCALING_TO_FIX_UNDERFLOW
				// add scaling corrections
				_Parameter subTimes = t->GetLeafCount();
				//_String    outString = _String("Likelihood evaluation ") & likeFuncEvalCallCount & ": before " & res;
				//StringToConsole (outString);
				for (i=0; i<df->theFrequencies.lLength; i++)
				{
					if (t->scalingForUnderflow->theData[i] < 0.)
					{
						_String errMsg = _String ("Failed to scale site ") & i & " at likelihood evaluation " & likeFuncEvalCallCount;
						t->scalingForUnderflow->theData[i] = -t->scalingForUnderflow->theData[i];
						ReportWarning (errMsg);
						//StringToConsole(errMsg);
						//NLToConsole();
						res -= 1.e7;
					}
					else
						if (t->scalingForUnderflow->theData[i] != 1.0)
						{
							//char buffer [2048];
							//sprintf (buffer,"Scaled (LF Eval %d) at site %d scaling factor %g\n",likeFuncEvalCallCount, i, t->scalingForUnderflow->theData[i]);
							//BufferToConsole (buffer);
							res -= log (t->scalingForUnderflow->theData[i])*subTimes*df->theFrequencies.lData[i];
						}

				}
				//outString = _String (". After :") & res;
				//StringToConsole (outString); 
				//NLToConsole();
			#endif
			
			if (computationalResults.lLength>index)
				((_Constant*)computationalResults.lData[index])->SetValue (res);
			else
			{
				_Constant c_res (res);
				computationalResults&&(&c_res);
			}
				#if (!defined __UNIX__ || defined __HEADLESS__) && defined __MP__
				if (yieldTime)
					yieldCPUTime();
			#endif

			return res;
		}
		else
		// categories on
		{
			#ifdef __MP__
				if (theTasks)
					for (long tc = 0; tc<systemCPUCount-1; tc++)
						theTasks[tc].resultsStorage = siteRes;
			#endif
			if (filterUnit==1 && t->GetCodeBase()==4) // nucleotide data
			{
				if (t->HasCache())
				{
					#ifdef __MP__
					if (theTasks)
					{
						siteRes[*siteList] = t->ReleafTreeAndCheckChar4 (df, *siteList, true, categID);				
						for (long tc = 0; tc<systemCPUCount-1; tc++)
							#ifdef __MACHACKMP__
							if ( pthread_create( theThreads+tc, NULL, ThreadReleafFunctionMNucHook,(void*)(theTasks+tc)))
							#else
							if ( pthread_create( theThreads+tc, NULL, ThreadReleafFunctionMNuc,(void*)(theTasks+tc)))
							#endif
							{
								FlagError("Failed to initialize a POSIX thread in ComputeBlock()");
								exit(1);
							}		
												
						siteList++;
						for (i = 1; i<tstep; i++,siteList++,leafList++)
						{
							#if !defined __UNIX__ || defined __HEADLESS__
								siteEvalCount ++;
								if (divideBy && (siteEvalCount%divideBy == 0) || siteEvalCount%1000 == 0)
									yieldTime = true;
							#endif
							siteRes[*siteList]=t->ThreadReleafTreeChar4 (df, *siteList, *(siteList-1), (*leafList)&0x0000ffff, (*leafList)>>16, 
																		i*t->categoryCount+categID);
						}
						for (long wc = 0; wc<systemCPUCount-1; wc++)
							if ( pthread_join ( theThreads[wc], NULL ) ) 
							{
								FlagError("Failed to join POSIX threads in ComputeBlock()");
								exit(1);
							}
						delete theTasks;
						delete theThreads;	
					}
					else
					{
					#endif
						siteRes[*siteList] = t->ReleafTreeAndCheckChar4 (df, *siteList, true, categID);				
						siteList++;
						for (int i = 1; i<sl->lLength; i++,	siteList++,	leafList++)
						{
							#if !defined __UNIX__ || defined __HEADLESS__
							siteEvalCount ++;
							if (divideBy && (siteEvalCount%divideBy == 0) || siteEvalCount%1000 == 0)
								#ifdef __MP__
								yieldTime = true;
								#else
								yieldCPUTime();
								#endif
							#endif
							siteRes[*siteList]=t->ThreadReleafTreeChar4 (df, *siteList, *(siteList-1), (*leafList)&0x0000ffff, (*leafList)>>16, i*t->categoryCount+categID);
						}
					
						
					#ifdef __MP__
					}
					#endif
				}
				else
				{
					if (t->IsDegenerate())
					{
						siteRes[*siteList] = t->ReleafTreeAndCheckChar4 (df, *siteList, false, categID);				
						siteList++;
						for (int i = 1; i<sl->lLength; i++,siteList++,leafList++)
						{
							#if !defined __UNIX__ || defined __HEADLESS__
								siteEvalCount ++;
								if (divideBy && (siteEvalCount%divideBy == 0) || siteEvalCount%1000 == 0)
									yieldCPUTime();
							#endif
							siteRes[*siteList]=t->ReleafTreeChar4Degenerate (df, *siteList);
						}
					}
					else
					{
						siteRes[*siteList] = t->ReleafTreeAndCheckChar4 (df, *siteList, false, categID);				
						siteList++;
						if (df->IsNormalFilter())
						{
							for (int i = 1; i<sl->lLength; i++,siteList++,leafList++)
							{
								#if !defined __UNIX__ || defined __HEADLESS__
									siteEvalCount ++;
									if (divideBy && (siteEvalCount%divideBy == 0) || siteEvalCount%1000 == 0)
										yieldCPUTime();
								#endif
								siteRes[*siteList]=t->ReleafTreeChar4 (df, *siteList, *(siteList-1), (*leafList)&0x0000ffff, (*leafList)>>16);
							}
						}
						else
						{
							/*for (i = 0; i < 3; i++)
							{
								_Matrix *mp = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[t->theRoot->nodes.data[2]->in_object])->GetCompExp();
								_String t ((long)mp);
								StringToConsole (t);
								NLToConsole ();
							}*/
							for (i = 1; i<sl->lLength; i++,siteList++)
								siteRes[*siteList] = t->ReleafTreeCharNumFilter4Tree3 ((_DataSetFilterNumeric*)df, *siteList,categID);
						}
				
					}
				}
			}
			else
			{
				if (df->GetUnitLength()>1)
				{
					if (t->HasCache())
					{
						#ifdef __MP__
						if (theTasks)
						{
							siteRes[*siteList] = t->ReleafTreeAndCheck (df, *siteList, true, categID);				
							for (long tc = 0; tc<systemCPUCount-1; tc++)
								#ifdef __MACHACKMP__
								if ( pthread_create( theThreads+tc, NULL, ThreadReleafFunctionMCodonHook,(void*)(theTasks+tc)))
								#else
								if ( pthread_create( theThreads+tc, NULL, ThreadReleafFunctionMCodon,(void*)(theTasks+tc)))
								#endif
								{
									FlagError("Failed to initialize a POSIX thread in ComputeBlock()");
									exit(1);
								}		
													
							siteList++;
							for (i = 1; i<tstep; i++,siteList++,leafList++)
							{
								#if !defined __UNIX__ || defined __HEADLESS__
									siteEvalCount ++;
									if (divideBy && (siteEvalCount%divideBy == 0) || siteEvalCount%1000 == 0)
										yieldTime = true;
								#endif

								siteRes[*siteList]=t->ThreadReleafTreeCache (df, *siteList, *(siteList-1), (*leafList)&0x0000ffff, (*leafList)>>16, 
																			i*t->categoryCount+categID);
							}
							for (long wc = 0; wc<systemCPUCount-1; wc++)
								if ( pthread_join ( theThreads[wc], NULL ) ) 
								{
									FlagError("Failed to join POSIX threads in ComputeBlock()");
									exit(1);
								}
							delete theTasks;
							delete theThreads;				
						}
						else
						{
						#endif
							siteRes[*siteList] = t->ReleafTreeAndCheck (df, *siteList, true, categID);								
							siteList++;
							for (int i = 1; i<sl->lLength; i++,siteList++,leafList++)
							{
								#if !defined __UNIX__ || defined __HEADLESS__
									siteEvalCount ++;
									if (divideBy && (siteEvalCount%divideBy == 0) || siteEvalCount%1000 == 0)
										#ifdef __MP__
											yieldTime = true;
										#else
											yieldCPUTime();
										#endif
									
								#endif
								siteRes[*siteList]=t->ThreadReleafTreeCache (df, *siteList, *(siteList-1), (*leafList)&0x0000ffff, (*leafList)>>16, i*t->categoryCount+categID);
							}
						#ifdef __MP__
						}
						#endif
					}
					else
					{
						if (t->IsDegenerate())
						{
							siteRes[*siteList] = t->ReleafTreeAndCheck (df, *siteList, false, categID);				
							siteList++;
							for (int i = 1; i<sl->lLength; i++,siteList++,leafList++)
							{
								#if !defined __UNIX__ || defined __HEADLESS__
									siteEvalCount ++;
									if (divideBy && (siteEvalCount%divideBy == 0) || siteEvalCount%1000 == 0)
										yieldCPUTime();
								#endif
								siteRes[*siteList]=t->ReleafTreeDegenerate (df, *siteList);
							}
						}
						else
						{
							siteRes[*siteList] = t->ReleafTreeAndCheck (df, *siteList, false, categID);				
							siteList++;
							for (int i = 1; i<sl->lLength; i++,siteList++,leafList++)
							{
								#if !defined __UNIX__ || defined __HEADLESS__
									siteEvalCount ++;
									if (divideBy && (siteEvalCount%divideBy == 0) || siteEvalCount%1000 == 0)
										yieldCPUTime();
								#endif
								siteRes[*siteList]=t->ReleafTree (df, *siteList, *(siteList-1), (*leafList)&0x0000ffff, (*leafList)>>16);
							}
						}
					}
				}
				else
				{
					if (t->HasCache())
					{
					#ifdef __MP__
						if (theTasks)
						{
							siteRes[*siteList] = t->ReleafTreeAndCheckChar (df, *siteList, true, categID);				
							for (long tc = 0; tc<systemCPUCount-1; tc++)
								#ifdef __MACHACKMP__
								if ( pthread_create( theThreads+tc, NULL, ThreadReleafFunctionMAAHook,(void*)(theTasks+tc)))
								#else
								if ( pthread_create( theThreads+tc, NULL, ThreadReleafFunctionMAA,(void*)(theTasks+tc)))
								#endif
								{
									FlagError("Failed to initialize a POSIX thread in ComputeBlock()");
									exit(1);
								}		
													
							siteList++;
							for (i = 1; i<tstep; i++,siteList++,leafList++)
							{
								#if !defined __UNIX__ || defined __HEADLESS__
									siteEvalCount ++;
									if (divideBy && (siteEvalCount%divideBy == 0) || siteEvalCount%1000 == 0)
										yieldTime = true;
								#endif
								siteRes[*siteList]=t->ThreadReleafTreeCharCache (df, *siteList, *(siteList-1), (*leafList)&0x0000ffff, (*leafList)>>16, 
																			i*t->categoryCount+categID);
							}
							for (long wc = 0; wc<systemCPUCount-1; wc++)
								if ( pthread_join ( theThreads[wc], NULL ) ) 
								{
									FlagError("Failed to join POSIX threads in ComputeBlock()");
									exit(1);
								}
							delete theTasks;
							delete theThreads;				
						}
						else
						{
						#endif
							siteRes[*siteList] = t->ReleafTreeAndCheckChar (df, *siteList, true, categID);				
							siteList++;
							for (int i = 1; i<sl->lLength; i++,siteList++,leafList++)
							{
								#if !defined __UNIX__ || defined __HEADLESS__
									siteEvalCount ++;
									if (divideBy && (siteEvalCount%divideBy == 0) || siteEvalCount%1000 == 0)
									#ifdef __MP__
											yieldTime = true;
										#else
											yieldCPUTime();
										#endif
										
								#endif
								siteRes[*siteList]=t->ThreadReleafTreeCharCache (df, *siteList, *(siteList-1), (*leafList)&0x0000ffff, (*leafList)>>16, i*t->categoryCount+categID);
							}
						#ifdef __MP__
						}
						#endif
					}
					else
					{
						if (t->IsDegenerate())
						{
							siteRes[*siteList] = t->ReleafTreeAndCheckChar (df, *siteList, false, categID);				
							siteList++;
							for (int i = 1; i<sl->lLength; i++,siteList++,leafList++)
							{
								#if !defined __UNIX__ || defined __HEADLESS__
									siteEvalCount ++;
									if (divideBy && (siteEvalCount%divideBy == 0) || siteEvalCount%1000 == 0)
										yieldCPUTime();
								#endif
								siteRes[*siteList]=t->ReleafTreeCharDegenerate (df, *siteList);
							}
						}
						else
						{
							siteRes[*siteList] = t->ReleafTreeAndCheckChar (df, *siteList, false, categID);				
							siteList++;
							for (int i = 1; i<sl->lLength; i++,siteList++,leafList++)
							{
								#if !defined __UNIX__ || defined __HEADLESS__
									siteEvalCount ++;
									if (divideBy && (siteEvalCount%divideBy == 0) || siteEvalCount%1000 == 0)
										yieldCPUTime();
								#endif
								siteRes[*siteList]=t->ReleafTreeChar (df, *siteList, *(siteList-1), (*leafList)&0x0000ffff, (*leafList)>>16);
							}
						}
					}			
				}
			}
			#if !defined __UNIX__ || defined __HEADLESS__
				#ifdef __MP__
					if (yieldTime)
						yieldCPUTime();
				#endif
			#endif

		}
	}
	return 0.0;

}

//_______________________________________________________________________________________
long		_LikelihoodFunction::CostOfPath	 (_DataSetFilter* df, _TheTree* t, _SimpleList& sl, _SimpleList* tcc)
{
	long res = 0;
	for (long i=0; i<(long)sl.lLength-1; i++)
		res+=t->ComputeReleafingCost (df,sl.lData[i],sl.lData[i+1], tcc, i+1);
	return res;
}
//_______________________________________________________________________________________

/*void		_LikelihoodFunction::DumpingOrder	 (long index, _SimpleList& sl)
{
// 
	_DataSetFilter* df = (_DataSetFilter*)(dataSetFilterList(theDataFilters (index)));
	_TheTree		*t = (_TheTree*)LocateVar(theTrees(index));
	t->DumpingOrder(df,sl);
}*/
			


//_______________________________________________________________________________________

//#define __SORTING_DEBUG
#ifdef 	 __SORTING_DEBUG
	#include "HYChartWindow.h"
#endif

void	countingTraverse (node<long>* startingNode, long& totalLength, long currentSize, long& maxSize, bool add2Size)
{
	if (startingNode->parent)
		totalLength+=startingNode->in_object;
		
	if (add2Size)
		currentSize++;
		
	if (currentSize>maxSize)
		maxSize = currentSize;
		
	for (long k=1; k<startingNode->nodes.length; k++)
		countingTraverse (startingNode->go_down(k), totalLength, currentSize, maxSize, true);
		
	if (startingNode->nodes.length)
		countingTraverse (startingNode->go_down(startingNode->nodes.length), totalLength, currentSize, maxSize, false);
}

//_______________________________________________________________________________________

void	countingTraverseArbRoot (node<long>* startingNode, node<long>* childNode, long& totalLength, long currentSize, long& maxSize)
{
	if (childNode)
		for (long k=1; k<=startingNode->nodes.length; k++)
		{
			node<long>* cNode = startingNode->go_down(k);
			if (cNode!=childNode)
				countingTraverse (cNode, totalLength, currentSize, maxSize, true);
		}
	else
		for (long k=1; k<=startingNode->nodes.length; k++)
			countingTraverse (startingNode->go_down(k), totalLength, currentSize, maxSize, true);
			
	if (startingNode->parent)
	{
		totalLength+=startingNode->in_object;
		countingTraverseArbRoot (startingNode->parent, startingNode, totalLength, currentSize, maxSize);
	}
}
//_______________________________________________________________________________________

long	findAvailableSlot (_SimpleList& cs, long& ff)
{
	for (long k = ff; k<cs.lLength; k++)
		if (cs.lData[k] == -1)
		{
			ff = k+1;
			return k;
		}
	{
	for (long k = 0; k<ff; k++)
		if (cs.lData[k] == -1)
		{
			ff = k+1;
			return k;
		}
	}
	cs << -1;
	ff = 0;
	return cs.lLength-1;
}

//_______________________________________________________________________________________

void	setComputingArrays (node<long>* startingNode, node<long>* childNode, _SimpleList& co, _SimpleList& so, _SimpleList &cs, _SimpleList& ro, _SimpleList& po, long& ff)
{
	bool isFirstNode = (co.lLength == 0);
	
	co << startingNode->in_object;
	
	if (isFirstNode || startingNode->nodes.length)
	{
		long id = findAvailableSlot(cs,ff);
		cs.lData[id] = startingNode->in_object;
		startingNode->in_object = id;
		so << id;
		if (isFirstNode)
		{
			ro << -1;
			po << -1;
		}
	}
	else
		so << -1;
		
	if (childNode) // above root node
	{
		ro << cs.lData[childNode->in_object];
		po << (long)childNode;
		
		if (startingNode->parent)
		{
			for (long k=1; k<=startingNode->nodes.length; k++)
			{
				node<long>* cNode = startingNode->go_down(k);
				if (cNode!=childNode)
					setComputingArrays (startingNode->go_down(k), nil, co,so,cs,ro,po,ff);
			}		
			cs.lData[childNode->in_object] = -1;
			setComputingArrays (startingNode->parent, startingNode, co,so,cs,ro,po,ff);
		}
		else // actual tree root
		{
			bool passedUpPath = false;
			
			for (long k=1; k<=startingNode->nodes.length; k++)
			{
				node<long>* cNode = startingNode->go_down(k);
				if (cNode!=childNode)
				{
					if ((k==startingNode->nodes.length)||((k==startingNode->nodes.length-1)&&(!passedUpPath)))
						cs.lData[startingNode->in_object] = -1;
					setComputingArrays (startingNode->go_down(k), nil, co,so,cs,ro,po,ff);
				}
				else	
					passedUpPath = true;
			}		
		}
	}
	else // "below" root node
	{
		if ((!isFirstNode)&&startingNode->parent)
		{
			ro << startingNode->parent->in_object;
			po << (long)startingNode->parent;
		}
			
		for (long k=1; k<startingNode->nodes.length; k++)
			setComputingArrays (startingNode->go_down(k), nil, co,so,cs,ro,po,ff);
			
		if (!isFirstNode)
			if (startingNode->nodes.length)
				cs.lData[startingNode->in_object] = -1;

		if (startingNode->nodes.length)
			setComputingArrays (startingNode->go_down(startingNode->nodes.length), nil, co,so,cs,ro,po,ff);

		if (isFirstNode&&startingNode->parent)
		{
			cs.lData[startingNode->in_object] = -1;
			setComputingArrays (startingNode->parent,startingNode, co,so,cs,ro,po,ff);
		}
		
	}
}

//_______________________________________________________________________________________

void		_LikelihoodFunction::OptimalOrder	 (long index, _SimpleList& sl)
{
	
	_DataSetFilter* df = (_DataSetFilter*)(dataSetFilterList(theDataFilters (index)));
		
	long 			partition 			= -1, 
					totalSites 			= 0, 
					completedSites 		= 0,
	 				startpt, 
	 				endpt, 
	 				j, 
	 				max 				= -1, 
	 				k, 
	 				intI,
	 				totalLength;
	 				
					
	 				
	_Parameter		skipo = 1.0;
	_TheTree		*t = (_TheTree*)LocateVar(theTrees(index));
	checkParameter	(optimizeSummationOrder,skipo,1.0);
	
	if (!skipo || df->NumberDistinctSites()==1 || t->IsDegenerate() || !df->IsNormalFilter ()) // do not optimize 
	{
		for (k = 0; k < df->NumberDistinctSites(); k++)
		    sl<<k;
		return;
	}
	SetStatusLine ("Optimizing data ordering");
	checkParameter (optimizePartitionSize,skipo,0.0);
	totalSites = df->NumberDistinctSites();
	if (skipo) //  partition the sequence into smaller subseqs. for optimization
	{
		partition = (long)skipo;
		if ((partition<=0)||(partition>totalSites))
			partition = totalSites;
	}
	else // revert to default partitioning
		partition = totalSites>1500?1500:totalSites;
		
	long	vLevel 		 = VerbosityLevel(),
			globalLength = 0;

	if (vLevel>5)
	{
		char buffer [128];
		sprintf(buffer,"\nOptimizing Column Order for block %ld", index);
		BufferToConsole (buffer);
	}
	
	_SimpleList   partitionSites, distances, edges;
	
	
	while (completedSites<totalSites)
	{
		 if (totalSites-completedSites<partition) 
		 	partition = totalSites-completedSites;	
		 
		 intI = 0; // internal index for partition
		 
		 // populate sites allowed
		 
		 if (df->GetUnitLength()==1)
		 {
		 	for (k=completedSites+1; k<completedSites+partition; k++)
		 	{
		 		partitionSites<<k;
		 		distances<<t->ComputeReleafingCostChar (df, completedSites, k);
		 		edges<<completedSites;
		 	}
		 }
		 else
		 {
		 	for (k=completedSites+1; k<completedSites+partition; k++)
		 	{
		 		partitionSites<<k;
		 		distances<<t->ComputeReleafingCost (df, completedSites, k);
		 		edges<<completedSites;
		 	}
		 }
		 	
	  	 node<long>*   spanningTreeRoot;
		 _SimpleList   spanningTreePointers,
		 			   spanningTreeSites;
		 
		 if (mstCache)
		 {
			 spanningTreeRoot = new node<long>;
			 spanningTreeSites << 0;
			 for (k=completedSites+1; k<completedSites+partition; k++)
			 {
			  	spanningTreePointers << (long)spanningTreeRoot;
			 	spanningTreeSites << 0; 
			 }
			 spanningTreeSites.lData[0] = (long)spanningTreeRoot;
		 }
		
		 sl << completedSites;
		 	
		 while (intI<partition-1)
		 {
		 	// search for the shortest branch to add to the tree.
		 	max = 0x0fffffff;
		 	long * pl  = distances.lData;
		 		 	    
		 	for (k=0;k<distances.lLength;k++,pl++)
		 		if (*pl<=max)
		 		{
		 			max = *pl;
		 			endpt = partitionSites.lData[k];
		 			startpt = k;
		 		}
		 	
		 	// delete k from allowed sites
		 	partitionSites.Delete(startpt);
		 	distances.Delete(startpt);
		 	// insert the endpt into the tree before the pt that the edge came from
		 	
			node<long>*   spanningTreeNode;
			if (mstCache)
			{
			  	spanningTreeNode = new node<long>;
				spanningTreeNode->in_object = max;
			 	((node<long>*)spanningTreePointers.lData[startpt])->add_node (*spanningTreeNode);
			 	spanningTreeSites.lData[endpt-completedSites] = (long)spanningTreeNode;
			 	spanningTreePointers.Delete(startpt);
			 }
		 	
		 	j = sl.Find(edges.lData[startpt],completedSites);
		 	sl.InsertElement ((BaseRef)endpt,j+1,false,false);
		 	edges.Delete(startpt);
		 	
		 	// make one more pass and update the distances if needed
		 	if (df->GetUnitLength()==1)
		 	{
			 	for (k=0;k<distances.lLength;k++)
			 	{
			 		j = t->ComputeReleafingCostChar (df,endpt, partitionSites.lData[k]);
			 		if (j<distances.lData[k])
			 		{
			 			distances.lData[k]=j;
			 			edges.lData[k] = endpt;
			 			
						if (mstCache)
							spanningTreePointers.lData[k] = (long)spanningTreeNode;
			 		}
			 	}
			 }
			 else
		 	 {
			 	for (k=0;k<distances.lLength;k++)
			 	{
			 		j = t->ComputeReleafingCost (df,endpt, partitionSites.lData[k]);
			 		if (j<distances.lData[k])
			 		{
			 			distances.lData[k]=j;
			 			edges.lData[k] = endpt;

						if (mstCache)
			 				spanningTreePointers.lData[k] = (long)spanningTreeNode;
			 		}
			 	}
			 }
			 
			 intI++;
			 if (intI%50==0)
			 {
				#if !defined __UNIX__ || defined __HEADLESS__
					SetStatusBarValue ((intI+completedSites)*100/totalSites,1.0, 0.0);
				#endif
				//yieldCPUTime();
			 }
		 }
		 
		 
		 
		 partitionSites.Clear();
		 distances.Clear();
		 edges.Clear();
		 
		 if (mstCache)
		 {
		 	// print out the info on the spanning tree
		 	// keep track of the "deepest path"
		 	long	    maxLevel    = 0;
		 	
		 	totalLength = 0;					 				
		 	countingTraverse (spanningTreeRoot,totalLength,0,maxLevel,true);
		 	
		 	globalLength += totalLength;
		 	
			_String		* res = new _String((unsigned long)128,true);
		
			long 		level,
				 		nc			=   0;
			 		

			node<long>* curNode= DepthWiseStepTraverser(spanningTreeRoot),
					  *	bestRoot = spanningTreeRoot;
			
			while (curNode!=spanningTreeRoot)
			{
				long maxLevel2 = 0; 
				totalLength = 0;
				countingTraverseArbRoot (curNode, nil, totalLength, 1, maxLevel2);
				if (maxLevel2<maxLevel)
				{
					maxLevel = maxLevel2;
					bestRoot = curNode;	
				}
				curNode= DepthWiseStepTraverser((node <long>*) nil);
			}
			
			_SimpleList	  computingOrder,
						  storageOrder,
						  cacheSlots,
						  referenceOrder,
						  parentOrder;
						  
			for (level=0; level<maxLevel; level++)
				cacheSlots << -1;
				
			for (level=0; level<spanningTreeSites.lLength; level++)
				((node<long>*)spanningTreeSites.lData[level])->in_object = level+completedSites;
				
			setComputingArrays (spanningTreeRoot, nil, computingOrder, storageOrder, cacheSlots, referenceOrder, parentOrder, nc);
			
			for (level=0; level<spanningTreeSites.lLength; level++)
				((node<long>*)spanningTreeSites.lData[level])->in_object = level+completedSites;
				
			for (level=1; level<parentOrder.lLength; level++)
				parentOrder.lData[level] = ((node<long>*)parentOrder.lData[level])->in_object;

			if (completedSites)
			{
				level = mstCache->computingOrder.lLength-1;
				(*(_SimpleList*)mstCache->computingOrder(level)) << computingOrder;
				(*(_SimpleList*)mstCache->storageOrder(level)) << storageOrder;
				(*(_SimpleList*)mstCache->referenceOrder(level)) << referenceOrder;
				(*(_SimpleList*)mstCache->parentOrder(level)) << parentOrder;
				if (cacheSlots.lLength > mstCache->cacheSize.lData[level])
					mstCache->cacheSize.lData[level] = cacheSlots.lLength;
			}
			else
			{
				mstCache->computingOrder 	&& & computingOrder;
				mstCache->storageOrder 		&& & storageOrder;
				mstCache->referenceOrder 	&& & referenceOrder;
				mstCache->parentOrder 		&& & parentOrder;
				mstCache->cacheSize			<< cacheSlots.lLength;			
			}
						
			#ifdef __SORTING_DEBUG
				char		buffer[512];
				_Matrix 	dataMx (computingOrder.lLength,4,false,true);
				
				nc = 0;
				for (level = 0; level < computingOrder.lLength; level++)
				{
					dataMx.theData[nc++] = computingOrder.lData[level];
					dataMx.theData[nc++] = storageOrder.lData[level];
					dataMx.theData[nc++] = referenceOrder.lData[level];
					dataMx.theData[nc++] = parentOrder.lData[level];
				}
				

				_List		titles;
				titles  && &_String("Computing Order");
				titles  && &_String("Storage Order");
				titles  && &_String("Storage Reference Order");
				titles  && &_String("Parent Order");

				_HYChartWindow* ncw = new _HYChartWindow ("MST Debug", titles, dataMx, nil);
				ncw->Activate();

				for (level=0; level<spanningTreeSites.lLength; level++)
					((node<long>*)spanningTreeSites.lData[level])->in_object = level;
					
				nc = 0;
				for (level=1; level<computingOrder.lLength; level++)
					nc += t->ComputeReleafingCost (df,computingOrder.lData[level],parentOrder.lData[level]);
					
				sprintf (buffer,"\nDouble check cost %d\n", nc);
				BufferToConsole (buffer);

				nc = 0;
				curNode= DepthWiseStepTraverserLevel(myLevel,spanningTreeRoot);						
				node<long>*	 nextNode;
				level 		= myLevel;
				

				nextNode	=  DepthWiseStepTraverserLevel(myLevel,nil) ;

				while (nextNode)
				{
					if (level>lastLevel)
					{
						if (lastLevel)
							(*res)<<',';
						for (j=0;j<level-lastLevel;j++)
							(*res)<<'(';
					}
					else
						if (level<lastLevel)
						{
							for (j=0;j<lastLevel-level;j++)
								(*res)<<')';
						}
						else
							(*res)<<',';
					
					(*res)<<_String(nc++);
					
					lastLevel = level;
					level 	  = myLevel;
					curNode   = nextNode;
					nextNode	=  DepthWiseStepTraverserLevel(myLevel,nil) ;
				
				}
				for (j=0;j<lastLevel-level;j++)
					(*res)<<')';

				(*res).Finalize();
				
				_String		 stree ("MST_STRING");
				_FString fr  (res);
				setParameter (stree, &fr);
				
				//DeleteObject  (res);
			#else
				res->Finalize();
				DeleteObject (res);
			#endif
			 spanningTreeRoot->delete_tree();
			 delete (spanningTreeRoot); // dmalloc fix 06162005
		 	
		 }
				 
		 completedSites+=partition;
		 if (vLevel>5)
		 {
		 	char   buffer[64];
			sprintf(buffer,"\n%ld %% done", (long)(completedSites*100/totalSites));
			BufferToConsole (buffer);		
		 }
	}
		
	_SimpleList straight (sl.lLength, 0, 1),
				* tcc = nil;
	
#ifdef _SLKP_LFENGINE_REWRITE_
	 if (treeTraversalMasks.lLength > index)
		 tcc =  (_SimpleList*) treeTraversalMasks(index);
	
#endif
	
	_Parameter strl = CostOfPath (df,t,straight), 
			   optl = CostOfPath (df,t,sl,tcc);
	
	if (vLevel>500)
	{
		_String* opPath = (_String*)sl.toStr();
		BufferToConsole("\nSite ordering:");
		StringToConsole(*opPath);
		DeleteObject(opPath);
	}
	
	char buffer [512];
	sprintf (buffer,"\nPseudo-optimal path's cost %ld vs %ld for 1..k=> a %g x improvement",(long)optl,(long)strl,strl/(double)optl );
	ReportWarning (buffer);
	if (vLevel>5)
		BufferToConsole (buffer);		

	if (mstCache)
	{
		long     memOverhead = mstCache->cacheSize.lData[mstCache->cacheSize.lLength-1];
		if (memOverhead)
		{
			memOverhead *= (t->GetINodeCount()*(sizeof(_Parameter)*t->GetCodeBase()+sizeof(long)+sizeof (char))+t->GetLeafCount()*(sizeof(_Parameter)+sizeof(long)))/1024;
			sprintf (buffer,"\nIf using full MST heurisitcs: %ld vs %ld for 1..k=> a %g x (relative %g x) improvement with %d KB memory overhead",globalLength,(long)strl,strl/(double)globalLength,optl/(double)globalLength,memOverhead);
			ReportWarning (buffer);
			if (vLevel>5)
				BufferToConsole (buffer);		
		}
	}
		
	#ifdef __SORTING_DEBUG
		printf ("\nTheortical cost lower bound %ld",t->GetLowerBoundOnCost (df));		
		printf ("\nSave all cost lower bound %ld",t->GetLowerBoundOnCostWithOrder (df,&sl));
	#endif
}
//_______________________________________________________________________________________

void	_LikelihoodFunction::ComputePruningEfficiency (long& full, long& saved)
{
	full = 0;
	saved = 0;
	for (long i=0; i<theTrees.lLength; i++)
	{
		_TheTree 	*cT = ((_TheTree*)(LocateVar(theTrees(i))));
		_SimpleList *l  = (_SimpleList*)leafSkips (i);
		_PMathObj	lc  = cT->TipCount();
		
		long		 leafCount = lc->Value(),
					 iCount;
					 
		DeleteObject (lc);
		lc = cT->BranchCount ();
		iCount = lc->Value();
		DeleteObject (lc);
		
		saved += leafCount+iCount;
		full  += (leafCount+iCount) * (l->lLength+1);
		
		for (long k=0; k<l->lLength;k++)
		{
			unsigned long j = l->lData[k],
				 		  p1 = j&0xffff,
				 		  p2 = ((j>>16)&0xffff);
				 		  
			saved += leafCount - p1 - (leafCount - 1 - p2);
			saved += iCount - cT->GetLeftINodes().lData[p1];
		}
	}
}	
		
//_______________________________________________________________________________________

long	_LikelihoodFunction::CountObjects (char flag)
{
	switch (flag)
	{
		case 1:
			{
				long res = 0;
				for (long k=0; k<indexInd.lLength; k++)
				{
					_Variable *v = LocateVar (indexInd.lData[k]);
					if (v->IsGlobal()) res++;
				}
				return res;
			}
		case 2:
			return indexInd.lLength - CountObjects (1);
		case 3:
			return indexDep.lLength;
	}
	return theTrees.lLength;
}

//_______________________________________________________________________________________

void	_LikelihoodFunction::SerializeLF (_String& rec, char opt, _SimpleList * partitionList, _SimpleList* exportPart)
{
	// first - check how many data sets the LF depends on;
	// if only one - then spool it out into a NEXUS string;
	// if more than one - output a NEXUS file where
	// data sets after the first one are embedded as strings.
	
	if (partitionList)
	// validate the list
	{
		_String errMsg; 
		if (partitionList->lLength == 0 || partitionList->lLength > theDataFilters.lLength)
			errMsg = "The partition list for sub-export passed to SerializeLF was either empty or too long";
		else
		{
			// check for duplicates and index overrun
			partitionList->Sort();
			if (partitionList->lData[partitionList->lLength-1] >= theDataFilters.lLength)
				errMsg = "The partition list contained invalid partition indices (too high) in call to SerializeLF";
			else
				for (long k2 = 1; k2 < partitionList->lLength; k2=k2+1)
					if (partitionList->lData[k2] == partitionList->lData[k2-1])
					{
						errMsg = "The partition list contained duplicate partition indices in call to SerializeLF";
						break;
					}
				
		}
		
		if (errMsg.sLength)
		{
			WarnError (errMsg);
			return;
		}
	}
	
	_String *lfName = (_String*)likeFuncNamesList (likeFuncList._SimpleList::Find((long)this));

	if (opt == 5)
	{
		_String resVarName = *lfName & "_MLES";
		_Variable* resVar = FetchVar(LocateVarByName (resVarName));
		if (resVar)
		{
			rec << & resVarName;
			rec << '=';
			_String* resVarVal = (_String*)resVar->Compute()->toStr();
			rec << resVarVal; 
			rec << ";\n";
			DeleteObject (resVarVal);
		}
		
		resVarName = *lfName & "_MLE_VALUES";
		rec << & resVarName;
		rec << "={};\n";
		
		for (long k=0; k<indexInd.lLength; k++)
		{
			_Variable *tiv = LocateVar(indexInd.lData[k]);
			rec << &resVarName;
			rec << "[\"";
			rec << *tiv->GetName();
			rec << "\"]=";
			rec << _String(tiv->Compute()->Value());
			rec << ";\n";
		}
		{
		for (long k=0; k<indexDep.lLength; k++)
		{
			_Variable *tiv = LocateVar(indexDep.lData[k]);
			rec << &resVarName;
			rec << "[\"";
			rec << tiv->GetName();
			rec << "\"]=";
			rec << _String(tiv->Compute()->Value());
			rec << ";\n";
		}
        }
		return;
	}

				
	_SimpleList * redirector = nil;
	if (partitionList)
	{
		redirector = new _SimpleList;
		checkPointer (redirector);
		for (long pidx = 0; pidx < partitionList->lLength; pidx = pidx+1)
			(*redirector) << theDataFilters.lData[partitionList->lData[pidx]];
	}
	else
		redirector = &theDataFilters;
				
				
				
	_SimpleList	  taggedDataSets,
				  dataSetsByFilter;
				  
	_AVLListX	  taggedDS (&taggedDataSets);
	
	for (long idx = 0; idx < redirector->lLength; idx++)
	{
		long tIdx = dataSetList._SimpleList::Find ((long)(((_DataSetFilter*)dataSetFilterList (redirector->lData[idx]))->GetData()));
		tIdx = taggedDS.Insert ((BaseRef)tIdx, taggedDS.countitems());
		
		if (tIdx < 0)
			tIdx = -tIdx-1;
			
		dataSetsByFilter << tIdx;
	}
	
	if (taggedDS.countitems()>1 && exportPart)
	{
		_String errMsg = _String("Can't represent likelihood function ") & *(_String*)likeFuncNamesList(likeFuncList._SimpleList::Find ((long)this)) & 
						" as a single file with a prespecified partition, because it depends on multiple datasets. This is an internal error.";
		WarnError (errMsg);
		return;
	}
	
	_SimpleList indexedDataSets (taggedDS.countitems(),0,0);
	for (long idx2 = 0; idx2 < taggedDataSets.lLength; idx2++)
		indexedDataSets.lData[taggedDS.xtraD.lData[idx2]] = taggedDataSets.lData[idx2];
	
	_List		 involvedSites,
				 avlSupport,
				 dV;
	
	_SimpleList	 unique_sites;
				 
	for (long idx3 = 0; idx3 < indexedDataSets.lLength; idx3++)
	{
		unique_sites << 0;
		_SimpleList * esl  = new _SimpleList;
		checkPointer (esl);
		_AVLListX	* invS = new _AVLListX (esl);
		avlSupport << esl;
		involvedSites << invS;
		DeleteObject (esl);
		checkPointer (esl  = new _SimpleList);
		dV << esl;
		DeleteObject (esl);
		DeleteObject (invS);
	}
	
	// create a dummy data filter and spool the data to a file
	stashParameter (skipOmissions,0.0,true);

	if (partitionList || !exportPart)
	{
		for (long idx = 0; idx < redirector->lLength; idx++)
		{
			_SimpleList* originalOrderFF = &((_DataSetFilter*)dataSetFilterList (redirector->lData[idx]))->theOriginalOrder;
			_AVLListX  * involvedSitesL  = (_AVLListX*)(involvedSites (dataSetsByFilter.lData[idx]));
			long	   * unique_sitesL   = unique_sites.lData + dataSetsByFilter.lData[idx];
			
			for (long idx2 = 0; idx2 < originalOrderFF->lLength; idx2++)
				if (involvedSitesL->Insert ((BaseRef)(originalOrderFF->lData[idx2]),*unique_sitesL) >= 0)
					(*unique_sitesL)++;
		}
		
		for (long idx4 = 0; idx4 < indexedDataSets.lLength; idx4++)
		{
			_AVLListX  * involvedSitesL  = (_AVLListX*)involvedSites (idx4);

			_SimpleList tcache,
						* dVL =  (_SimpleList*)dV (idx4);
			
			long		iv,
						k = involvedSitesL->Traverser (tcache, iv, involvedSitesL->GetRoot());
						
			for (; k>=0; k = involvedSitesL->Traverser (tcache, iv))
				(*dVL) << (long)involvedSitesL->Retrieve (k);
		}
	}
	else
		if (exportPart)
			 ((_SimpleList*) dV (dataSetsByFilter(0)))->Duplicate (exportPart);
	
	
	for (long idx5 = 0; idx5 < indexedDataSets.lLength; idx5++)
	{
		_DataSetFilter 	* entireSet = new _DataSetFilter ();
		checkPointer (entireSet);
		_SimpleList dH; 
		stashParameter (skipOmissions,0.0,true);
		entireSet->SetFilter ((_DataSet*)dataSetList(indexedDataSets.lData[idx5]),1,dH,*(_SimpleList*)dV (idx5));
		stashParameter (skipOmissions,0.0,false);
		if (idx5)
			stashParameter (dataFilePrintFormat,0.0,true);
		else
			stashParameter (dataFilePrintFormat,4.0,true);
		
		_String * cs =(_String*) entireSet->toStr ();
		
		if (idx5)
		{
			if (idx5 == 1)
				rec << "\n\nBEGIN HYPHY;\n\n";
			rec  << "_tdsstring_ = \"";
			rec.EscapeAndAppend (*cs);
			rec << "\";\nDataSet ";
			rec << (_String*)dataSetNamesList (indexedDataSets.lData[idx5]);
			rec << " = ReadFromString (_tdsstring_);\n_tdsstring_=0;\n";
		}
		else
		{
			rec << *cs;
			DeleteObject (cs);
		}
		DeleteObject (entireSet);
		stashParameter (dataFilePrintFormat,0.0,false);
	}
	
	if (indexedDataSets.lLength == 1)
		rec << "\n\nBEGIN HYPHY;\n\n";
	
	if (opt==1)
	{
		_Parameter     p1,p2;
		checkParameter (useLastResults,p1,0.0);
		checkParameter (useInitialDistanceGuess,p2,1.0);
		if (CheckEqual (p1,0.0) && p2>.1)
		{
			for (long i=0;i<indexInd.lLength;i++)
				LocateVar(indexInd.lData[i])->MarkDone();
			
			GetInitialValues();	
		}
		long	f = LocateVarByName (mpiPrefixCommand);
		
		if (f>=0)
		{
			_PMathObj iv = FetchVar (f)->Compute();
			if (iv->ObjectClass () == STRING)
				rec << *((_FString*)iv)->theString;
		}	
	}
		
	// write out all globals
	
	_String		glVars (1024L,true),
				locVars(1024L,true);
								
	char		str[4096];
		
	_SimpleList * indepVarList = nil,
				* depVarList   = nil,
				* catVarList   = nil;
	
	if (partitionList)
	{
		indepVarList 	= new _SimpleList;
		depVarList 		= new _SimpleList;
		catVarList 		= new _SimpleList;
		
		checkPointer ((Ptr)(indepVarList && depVarList && catVarList));
		ScanAllVariablesOnPartition (*partitionList, *indepVarList, *depVarList, *catVarList);
	}
	else
	{
		indepVarList = &indexInd;
		depVarList	 = &indexDep;
		catVarList   = &indexCat;
	}
	
	ExportIndVariables (glVars,locVars,indepVarList);
	ExportDepVariables (glVars,locVars,depVarList);
	
	glVars.Finalize();
	locVars.Finalize();
	
	rec << glVars;
	
	// write out all categs
	
	if (opt == 4)
		for (long idx = 0; idx < catVarList->lLength; idx++)
		{
			_CategoryVariable* theC = (_CategoryVariable*)LocateVar(catVarList->lData[idx]);
			sprintf (str, "\nglobal %s;", theC->GetName()->getStr());
			rec<<str;			
		}
	else
	{
		ExportCatVariables (rec, catVarList);
	}
	
	//write out all the models

	
	_SimpleList * redirectorT = nil,
				  dH,
				  dV2,
				  dU;
				  
	if (partitionList)
	{
		redirectorT = new _SimpleList;
		checkPointer (redirectorT);
		for (long pidx = 0; pidx < partitionList->lLength; pidx = pidx+1)
			(*redirectorT) << theTrees.lData[partitionList->lData[pidx]];
	}
	else
		redirectorT = &theTrees;
	{
	for (long idx = 0; idx < redirectorT->lLength; idx++)
	{
		_SimpleList dT;
		((_TheTree*)LocateVar(redirectorT->lData[idx]))->CompileListOfModels(dT);
		
		if (dT.lLength == 1)
			dV2<<dT.lData[0];
		else
			dV2<<-1;
			
		dH.Union (dU,dT);
		dU.Duplicate (&dH);
	}
	}
	
	{
		_SimpleList modelAux;
		_AVLList	modelDupChecker (&modelAux);
		for (long idx = 0; idx<dH.lLength; idx++)
			SerializeModel (rec, dH.lData[idx],&modelDupChecker);
		
	}
		
		
	// write out all the trees, including model definitions if needed
	
	_Parameter     stashIM = 0.0;
	checkParameter (includeModelSpecs,stashIM,0.0);
	{
	for (long idx = 0; idx < redirectorT->lLength; idx++)
	{
		_String* cStr;
		
		if (dV2.lData[idx]>=0)
		{
			rec << "\nUseModel (";
			rec << *((_String*)modelNames (dV2.lData[idx]));
			rec << ");\n";
			setParameter (includeModelSpecs,0.0);
		}
		else
			setParameter (includeModelSpecs,1.0);
			
		rec << "Tree ";
		rec << LocateVar(redirectorT->lData[idx])->GetName();
		rec << '=';
		cStr = (_String*)((_TheTree*)LocateVar(redirectorT->lData[idx]))->toStr();
		rec << *cStr;
		rec << '\n';
		DeleteObject (cStr);
	}
    }
	setParameter (includeModelSpecs,stashIM);
	
	rec 	<< locVars;
	
	// spool out the specs for datafilter
	
	rec  	<< "\nDataSet ";
	rec 	<< ((_String*)dataSetNamesList(indexedDataSets.lData[0]))->getStr();	
	rec		<< " = ReadDataFile(";
	rec     << useNexusFileData;
	rec     << ");\n";
	
	dU.Clear();
	_AVLList writtenDF (&dU);
	{
	for (long idx = 0; idx < redirector->lLength; idx++)
	{
		if (writtenDF.Insert ((BaseRef)redirector->lData[idx])>=0)
		{
			_DataSetFilter* theDF = (_DataSetFilter*)dataSetFilterList (redirector->lData[idx]);
			_String		  * horPart;
			
			 if (partitionList || !exportPart)
			 {
				 _SimpleList remappedOO;
				 _AVLListX	 * involvedSitesL = (_AVLListX*)involvedSites (dataSetsByFilter.lData[idx]);
				 for (long idx2 = 0; idx2 < theDF->theOriginalOrder.lLength; idx2++)
					remappedOO << involvedSitesL->GetXtra(involvedSitesL->Find ((BaseRef)theDF->theOriginalOrder.lData[idx2]));
				 horPart = (_String*)remappedOO.ListToPartitionString();
			 }
			 else
				 if (exportPart)
					 horPart = new _String(_String("0-") & (long)exportPart->lLength-1);
				// else
					// horPart = (_String*)theDF->theOriginalOrder.ListToPartitionString();
			 
			 rec    << "DataSetFilter ";
			 rec    << ((_String*)dataSetFilterNamesList(redirector->lData[idx]))->getStr();
			 rec    << " = CreateFilter(";
			 rec 	<< ((_String*)dataSetNamesList(indexedDataSets.lData[dataSetsByFilter.lData[idx]]))->getStr();	
			 rec    << ',';
			 rec    << _String((long)theDF->GetUnitLength());
			 rec    << ',';
			 rec    << '"';
			 rec    << *horPart;			
			 rec    << '"';
			 DeleteObject (horPart);
			 horPart = (_String*)theDF->theNodeMap.ListToPartitionString();
			 rec    << ',';
			 rec    << '"';
			 rec    << *horPart;			
			 rec    << '"';
			 DeleteObject (horPart);

			 horPart = theDF->GetExclusions();
			 if (horPart->sLength)
			 {
				rec    << ',';
				rec    << '"';
				rec    << *horPart;			
				rec    << '"';
				DeleteObject (horPart);
			 }
			 
			 rec << ");\n";	
		}
	}
	}
	
	rec << "LikelihoodFunction ";
	rec << *lfName;
	rec << " = (";
	
	long dsID = 0;
	{
	for (long idx = 0; idx < redirector->lLength; idx++)
	{
		if (dsID)
			rec << ',';
		else
			dsID = 1;
			
		rec << (_String*)dataSetFilterNamesList (redirector->lData[idx]);
		rec << ',';
		rec << *LocateVar(redirectorT->lData[idx])->GetName();
	}
	}
	if (computingTemplate&&templateKind == 1)
	{
		rec << ",\"";
		rec << (_String*)computingTemplate->toStr();
		rec << '"';
	}

	if (opt==1)
	{
		_Parameter pv;
		rec << ");\n";
		rec << optimizationPrecision;
		checkParameter (optimizationPrecision, pv ,0.001);
		rec << '=';
		rec << _String(pv);
		rec << ";\n";
		/*rec << optimizationMethod;
		checkParameter (optimizationMethod, pv ,7);
		rec << '=';
		rec << _String(pv);
		rec << ";\n";
		rec << useInitialDistanceGuess;
		checkParameter (useInitialDistanceGuess, pv ,1);
		rec << '=';
		rec << _String(pv);*/
		rec << ";\nOptimize(";
		rec << lfName;
		rec << "_MLES,";
		rec << lfName;
		rec << ");\n";
		rec << mpiMLELFValue;
		rec << '=';
		rec << lfName;
		rec << "_MLES[1][0];\n";
		rec << lf2SendBack;
		rec << "=\"";
		rec << lfName;
		rec << "\";\n";
		checkParameter (shortMPIReturn, pv ,0);
		rec << shortMPIReturn;
		rec << '=';
		rec << _String(pv);
		rec << ";\n";
	}
	else
		if (opt==2)
		{
			rec << ");\n";
			rec << mpiMLELFValue;
			rec << '=';
			rec << _String(FetchVar(LocateVarByName(mpiMLELFValue))->Compute()->Value());
			rec << ";\n";
			_String resVarName = *lfName & "_MLES";
			_Variable* resVar = FetchVar(LocateVarByName (resVarName));
			if (resVar)
			{
				rec << resVarName;
				rec << '=';
				_String* resVarVal = (_String*)resVar->Compute()->toStr();
				rec << *resVarVal; 
				rec << ";\n";
				DeleteObject (resVarVal);
			}
		}
		else
			rec << ");";
		
	rec << "\n\nEND;";
	
	if (partitionList)
	{
		DeleteObject (redirector);
		DeleteObject (indepVarList);
		DeleteObject (depVarList);
		DeleteObject (catVarList);
		DeleteObject (redirectorT);
	}
}

//_______________________________________________________________________________________

BaseRef	_LikelihoodFunction::toStr (void)
{
	_Parameter longOrShort,
			   value = 0.0;
			   
	checkParameter(likefuncOutput,longOrShort,2.0);
	
	if (longOrShort < 4.0)
	{
		PrepareToCompute();
		value = Compute();		
	}
	
	_String res(1,true);
	char str [2048];
	
	if (longOrShort>5.5) // serialized self-contained LF
	{
		_String	  *sLF = new _String (8192L, true);
		checkPointer 	(sLF);
		//_SimpleList 	dummySL (8,0,1);
		SerializeLF		(*sLF);
		sLF->Finalize	();
		return 		 	sLF;
	}
	
	if (longOrShort>=4.0) // just spool out the names of parameters
	{
		_Variable * thisVar;
		
		for (long i=0;i<indexInd.lLength;i++)
		{
			thisVar = LocateVar(indexInd.lData[i]);
			if (thisVar->IsGlobal())
				sprintf (str, "\nglobal %s=%.16g;", thisVar->GetName()->getStr(),(double)GetIthIndependent(i));
			else
				sprintf (str, "\n%s=%.16g;", thisVar->GetName()->getStr(),(double)GetIthIndependent(i));
				
			res<<str;

			if (!CheckEqual(thisVar->GetLowerBound(),DEFAULTPARAMETERLBOUND))
			{
				sprintf (str, "\n%s:>%.16g;", thisVar->GetName()->getStr(),(double)thisVar->GetLowerBound());
				res<<str;
			}
			if (!CheckEqual(thisVar->GetUpperBound(),DEFAULTPARAMETERUBOUND))
			{
				sprintf (str, "\n%s:<%.16g;", thisVar->GetName()->getStr(),(double)thisVar->GetUpperBound());
				res<<str;
			}
		}
		if (indexDep.lLength>0)
		{
			for (long i=0;i<indexDep.lLength;i++)
			{
				thisVar = LocateVar(indexDep.lData[i]);
				if (thisVar->IsGlobal())
					sprintf (str, "\nglobal %s",thisVar->GetName()->getStr());
				else
					sprintf (str, "\n%s",thisVar->GetName()->getStr());
				res<<str;
				res<<":=";
				_String* s = thisVar->GetFormulaString();
				res<<s;
				DeleteObject(s);
				res<<';';
				if (!CheckEqual(thisVar->GetLowerBound(),DEFAULTPARAMETERLBOUND))
				{
					sprintf (str, "\n%s:>%.16g;", thisVar->GetName()->getStr(),(double)thisVar->GetLowerBound());
					res<<str;
				}
				if (!CheckEqual(thisVar->GetUpperBound(),DEFAULTPARAMETERUBOUND))
				{
					sprintf (str, "\n%s:<%.16g;", thisVar->GetName()->getStr(),(double)thisVar->GetUpperBound());
					res<<str;
				}
			}
		}
		if (longOrShort>4.0)
			longOrShort = 2.5;
		else
		{
			DoneComputing();
			res.Finalize();
			return res.makeDynamic();
		}

	}
	if (longOrShort==3.0) // just spool out the names of parameters
	{
		sprintf (str, "Log Likelihood = %g;\n\n",(double)value);
		res<<str;
		sprintf (str,"Independent Parameters List\n");
		res<<str;
		for (long i=0;i<indexInd.lLength;i++)
		{
			sprintf (str, "\n  Parameter %4d : %s", i+1, LocateVar(indexInd.lData[i])->GetName()->getStr());
			res<<str;
		}
		if (indexDep.lLength>0)
		{
			sprintf (str,"\n\nConstrained Parameters List\n");
			res<<str;
			for (long i=0;i<indexDep.lLength;i++)
			{
				sprintf (str, "\n  Parameter %4d : %s", i+1, LocateVar(indexDep.lData[i])->GetName()->getStr());
				res<<str;
				res<<'=';
				_String* s = LocateVar(indexDep.lData[i])->GetFormulaString();
				res<<s;
				DeleteObject(s);
			}
		}
	}
	else
	if (longOrShort>1.1)
	{
		if (longOrShort!=2.5)
		{
			sprintf (str, "Log Likelihood = %15.15g;",(double)value);
			res<<str;
			bool globals = false;
			for (long i = 0; i<indexInd.lLength+indexDep.lLength; i++)
			{
				bool doDep = (i>=indexInd.lLength);
				_Variable* v = LocateVar (doDep?indexDep(i-indexInd.lLength):indexInd(i));
				if (v->IsGlobal())
				{
					if (!globals)
					{
						globals = true;
						res<< "\nShared Parameters:\n";
					}
					res<<v->GetName();
					res<<'=';
					if (doDep)
					{
						_String* fs = v->GetFormulaString();
						res<<fs;
						res<<'=';
						DeleteObject(fs);
					}
					_String* s = (_String*)v->toStr();
					res<<s;
					res<<'\n';
					DeleteObject(s);
				}
			}
		}
		// traverse the trees now
		for (long treeCounter = 0; treeCounter<theTrees.lLength; treeCounter++)
		{
			_TheTree* currentTree = (_TheTree*)LocateVar(theTrees(treeCounter));
			long level = 0,myLevel=0,lastLevel = 0,l1,l2,j;
			res<<"\nTree ";
			res<<currentTree->GetName();
			l1 = currentTree->GetName()->Length();
			res<<'=';
			_CalcNode* currentNode=currentTree->DepthWiseTraversalLevel(myLevel,true),*nextNode;
			level = myLevel;
			nextNode=currentTree->DepthWiseTraversalLevel(myLevel) ;
			_SimpleList iV, dV2;
			// decide if we can use expected substituion measure
			bool	useExpectedSubstitutions = longOrShort<3;
			/*if (blockDependancies.lData[i])
				useExpectedSubstitutions = false;*/
					
			while (nextNode)
			{
				iV.Clear();
				dV2.Clear();
				if (level>lastLevel)
				{
					if (lastLevel)
						res<<',';
					for (j=0;j<level-lastLevel;j++)
						res<<'(';
				}
				else
					if (level<lastLevel)
					{
						for (j=0;j<lastLevel-level;j++)
							res<<')';
					}
					else
						res<<',';
				_String num = currentNode->GetName()->Cut(l1+1,-1);
				l2 = l1+num.Length();
				res<<&num;
				if (!useExpectedSubstitutions)
				{
					res<<'{';
					{
						_AVLList 	iVA (&iV),
									dVA (&dV2);
									
						currentNode->ScanForVariables(iVA,dVA);
						currentNode->ScanForDVariables(dVA,iVA);
						
						iVA.ReorderList ();
						dVA.ReorderList ();
					}
					_String* s;
					_Variable* currentVar;
					
					for (j=0;j<iV.lLength;j++)
					{
						currentVar = LocateVar (iV(j));
						num = currentVar->GetName()->Cut(l2+2,-1);
						res<<&num;
						res<<'=';
						s = (_String*)currentVar->toStr();
						res<<s;
						if (j<iV.lLength-1)
							res<<',';
						DeleteObject(s);
					}
					
					for (j=0;j<dV2.lLength;j++)
					{
						if ((iV.lLength)||(j>0))
							res<<',';
						currentVar = LocateVar (dV2(j));
						num = currentVar->GetName()->Cut(l2+2,-1);
						res<<&num;
						res<<'=';
						s = currentVar->GetFormulaString();
						res<<s;
						res<<'=';
						DeleteObject(s);					
						s = (_String*)currentVar->toStr();
						res<<s;
						DeleteObject(s);
					}
					res<<'}';
				}
				else
				{
					res<<':';
					_Matrix*	theMx = currentNode->ComputeModelMatrix();
					
					_String expSubStr;
					if (!theMx)
						expSubStr = "NAN";
					else
					{
						_Parameter expSubs = currentNode->BranchLength();
						expSubStr = (_Parameter)fabs(expSubs);
					}
					res<<&expSubStr;
					
				}
				lastLevel = level;
				level = myLevel;
				currentNode = nextNode;
				nextNode=currentTree->DepthWiseTraversalLevel(myLevel) ;
			
			}
			for (j=0;j<lastLevel-level;j++)
			res<<')';
			res<<';';
		}
		res << '\n';
	}
	else
		if (longOrShort>.1)
		{
			sprintf (str, "Likelihood Function's Current Value = %15.15g\n",(double)value);
			_String num (str);
			res<<&num;
				
			for (long i = 0; i<indexInd.lLength+indexDep.lLength; i++)
			{
				bool doDep = (i>=indexInd.lLength);
				_Variable* v = LocateVar (doDep?indexDep(i-indexInd.lLength):indexInd(i));
				res<<v->GetName();
				if (doDep)
				{
					_String* fs = v->GetFormulaString();
					res<<fs;
					res<<'=';
					DeleteObject (fs);
				}
				res<<'=';
				_String* s = (_String*)v->toStr();
				res<<s;
				res<<'\n';
				DeleteObject(s);
			}
		}
		else
		{
			sprintf (str, "%15.15g",(double)value);		
			_String num (str);
			res<<&num;
		}
	
	if (longOrShort < 4.0)
		DoneComputing();
		
	res.Finalize();
	return res.makeDynamic();
	
}
	
//_______________________________________________________________________________________

void	_LikelihoodFunction::ReconstructAncestors (_DataSet &target, bool sample)
{
	long  i,
		  j,
		  k,
		  m,
		  n,
		  sc,
		  offset = 0,
		  offset2= 0;
	
	_DataSetFilter *dsf = (_DataSetFilter*)dataSetFilterList (theDataFilters(0));	
	target.SetTranslationTable (dsf->GetData());	

	_SimpleList 	thisSite,
					localExclusions;
	
	_TheTree    	*firstTree = (_TheTree*)LocateVar(theTrees(0));
	
	_Matrix			*rateAssignments = nil;
	
	_List			categoryVars;
	
	if  (indexCat.lLength>0)
	{
		rateAssignments = ConstructCategoryMatrix(false,false);
		checkPointer (rateAssignments);
		// use Bayes' formula now to obtain posterior probabilities of rates at sites
		// first, compute the probabilities for each class
		_Matrix		  probs (rateAssignments->GetHDim(),1,false,true);
		for (i = 0; i<theTrees.lLength; i++)
		{
			unsigned    long p2 = 0x80000000; // set high bit
			_SimpleList partitionClasses,
						partitionVars;
						
			_List       partitionProbs;
						
			m = 1;
						
			for (j=31; j>=0; j--)
			{
				if (blockDependancies.lData[i] & p2)
				{
					_CategoryVariable * cV = ((_CategoryVariable*)LocateVar(indexCat.lData[j]));
					
					partitionVars    << indexCat.lData[j];
					partitionClasses << (n=cV->GetNumberOfIntervals());
					partitionProbs   && cV->GetWeights ();
					m*=n;
				}
				p2 /= 2;
			}
			
			categoryVars && &partitionVars;
	
			// now compute the probs 
			
			for (j=0; j<m; j++)
			{
				k = j;
				
				_Parameter theProb = 1.0;
				
				for (n=0; n<partitionClasses.lLength; n++)
				{
					theProb *= (*(_Matrix*)partitionProbs (n)) [k%partitionClasses.lData[n]];
					k/=partitionClasses.lData[n];
				}
				
				probs[j] = theProb;
			}
			
			dsf = (_DataSetFilter*)dataSetFilterList (theDataFilters(i));
			
			
			for (j=offset; j<offset+dsf->theFrequencies.lLength; j++)
			{
				_Parameter maxProb  = 0.0;
						   
				long 	   maxIndex = 0;
				
				for (k=0; k<m; k++)
				{
					_Parameter tProb = (*rateAssignments)(k,j)*probs[k];
					if (tProb>maxProb)
					{
						maxProb  = tProb;
						maxIndex = k;
					}
				}
				
				// now that we found that maximum index, unroll it and 
				// put the data in the matrix

				rateAssignments->Store (0,j,maxIndex);
			}
			offset += dsf->theFrequencies.lLength;
		}
		offset = 0;
	}

	PrepareToCompute();

	computationalResults.Clear();
	Compute(); // set up all the matrices
	
	_Parameter			doRootSupportFlag = 0.0;
				
	checkParameter 		(storeRootSupportFlag, doRootSupportFlag, 0.0);
	_AssociativeList *  supportAVL = nil;
	
	if (doRootSupportFlag > 0.5)
		supportAVL = new _AssociativeList;

	for (i = 0; i<theTrees.lLength; i++)
	{
		_TheTree   *tree  = (_TheTree*)LocateVar(theTrees(i));
		long	   treeNameL = tree->GetName()->sLength;
		dsf = (_DataSetFilter*)dataSetFilterList (theDataFilters(i));
		
		bool  doSetRates = (rateAssignments && ((_SimpleList*)categoryVars(i))->lLength);
		if (i==0)
		{
			_CalcNode* travNode = tree->StepWiseTraversal (true);
			while (travNode)
			{
				if (!tree->IsCurrentNodeATip ())
				{
					_String nodeName (*travNode->GetName(),treeNameL+1,-1);
					target.GetNames() && & nodeName;
				}
				
				travNode = tree->StepWiseTraversal (false);
			}
			sc = target.GetNames().lLength;
		}
		else
			if (!tree->Equal(firstTree)) // incompatible likelihood function
			{
				ReportWarning ((_String("Ancestor reconstruction had to ignore part ")&_String(i+1)&" of the likelihood function since it has a different tree topology than the first part."));
				continue;
			}
			
		_List      	thisSet;
		
		_Parameter	verbLevel = VerbosityLevel();
		
		long		sitesDone = 0,
					filterStates = dsf->GetDimension (true); 

		_Matrix*    rootSupportMatrix = supportAVL&&(!sample) ? new _Matrix (dsf->theFrequencies.lLength,filterStates,false,true) : nil;
		
		for (j = 0; j<dsf->theFrequencies.lLength; j++)
		{
			if (doSetRates)
				tree->SetCompMatrices ((*rateAssignments)(0,offset2+j));
				
			if (sample)
			{
				_List * recStrings = new _List;
				checkPointer (recStrings);
				tree->ReleafTree (dsf, j, -1, 0, tree->flatCLeaves.lLength-1);
				for (long jj = 0; jj<dsf->theFrequencies.lData[j]; jj++)
				{
					_List * recoveredStrings = tree->SampleAncestors (dsf, tree->theRoot);
					(* recStrings) << recoveredStrings;
					DeleteObject (recoveredStrings);
					if (verbLevel > 9.999)
					{	
						char buffer[64];
						sprintf (buffer,"Site %d/%d\n", sitesDone++, dsf->theMap.lLength/dsf->GetUnitLength());
						BufferToConsole (buffer);		
					}

				}
				thisSet << recStrings;
				DeleteObject (recStrings);
			}
			else
			{
				_List *  recoveredStrings = tree->RecoverAncestralSequences (dsf,j,j-1, rootSupportMatrix?rootSupportMatrix->theData+j*filterStates: nil);
				for (m = 0; m<recoveredStrings->lLength;m++)
					thisSet << (*recoveredStrings)(m);
					
				DeleteObject (recoveredStrings);
			}
		}
		
		m = dsf->GetUnitLength();
		
		if (sample)
		{
			_SimpleList patMatched (dsf->theFrequencies.lLength); 

			for (j = 0; j<dsf->theOriginalOrder.lLength/m; j++)
			{
				k = dsf->duplicateMap.lData[j];
				
				_List*   thisList = (_List*)(*((_List*)thisSet(k)))(patMatched.lData[k]++);
				
				for (n = 0; n<m; n++)
					target.AddSite (((_String*)(*thisList)(n))->sData[0]);
					
				for (long p = 1; p<sc; p++)
					for (n = 0; n<m; n++)
						target.Write2Site (offset+j*m+n,((_String*)thisList->lData[n])->sData[p]);	
						
				for (n = 0; n<m; n++)
					target.Compact (offset+j*m+n);
						
				target.ResetIHelper();
			}
		}
		else
		{
			
			if (rootSupportMatrix)
			{
				long 	  k 				  = dsf->theOriginalOrder.lLength/m;
				_Matrix * mappedSupportMatrix = new _Matrix (k, filterStates,false,true);
				
				for (j = 0; j<dsf->theOriginalOrder.lLength/m; j++)
				{
					_Parameter * sourceVec = rootSupportMatrix->theData+dsf->duplicateMap.lData[j]*filterStates;
					
					for (long j2 = filterStates * j; j2 < filterStates * (j+1); j2++, sourceVec++)	
						mappedSupportMatrix->theData[j2] = *sourceVec;
						
				}	
				
				DeleteObject (rootSupportMatrix);
				_FString	 aKey;
				*aKey.theString = i;
				supportAVL->MStore (&aKey, mappedSupportMatrix, false);
			}
			
			for (j = 0; j<dsf->theOriginalOrder.lLength/m; j++)
			{
				k = dsf->duplicateMap.lData[j] * m;
				for (n = 0; n<m; n++)
					target.AddSite (((_String*)thisSet(k+n))->sData[0]);
				
			}
			
			for (long p = 1; p<sc; p++)
				for (j = 0; j<dsf->theOriginalOrder.lLength/m; j++)
				{
					k = dsf->duplicateMap.lData[j] * m;
					for (n = 0; n<m; n++)
						target.Write2Site (offset+j*m+n,((_String*)(thisSet.lData[k+n]))->sData[p]);	
				}
		}

		offset  += dsf->theOriginalOrder.lLength;
		offset2 += dsf->theFrequencies.lLength;
	}
	
	if (supportAVL)
	{
		_Variable * vid = CheckReceptacle (&supportMatrixVariable, "Ancestral State Reconstruction", false);
		if (vid)
			vid->SetValue (supportAVL, false);
		else
			DeleteObject (supportAVL);
	}
	
	for (i=0; i<theTrees.lLength; i++)
	{
		_TheTree * cT = ((_TheTree*)(LocateVar(theTrees(i))));
		cT->CleanUpMatrices();
	}
	target.Finalize();
	target.SetNoSpecies(target.GetNames().lLength);
	
	if (rateAssignments)
		DeleteObject (rateAssignments);
		
	DoneComputing ();
	
}

//_______________________________________________________________________________________

#ifdef 		__MP__
	void*	StateCounterMP (void* arg)
	{
		WancReleafTask * theTask = (WancReleafTask *)arg;
		
		_Parameter		vls = VerbosityLevel();
		
		for (long sites = theTask->startAt; sites < theTask->endAt; sites++)
		{
			_Parameter siteLikelihood = theTask->tree->ThreadReleafTreeCache (theTask->dsf, sites, ((sites>theTask->startAt)?sites-1:-1), 0, theTask->tree->flatCLeaves.lLength-1,sites,theTask->threadIndex);
			pthread_mutex_lock(&wancMutex);
			
			_Matrix	 res1(theTask->tree->GetCodeBase() , theTask->tree->GetCodeBase(), false, true),
					 res2(theTask->tree->GetCodeBase(),  theTask->tree->GetCodeBase(), false, true);

			pthread_mutex_unlock(&wancMutex);	
			
			if (vls>9.99)
			{
				char     buffer[64];
				sprintf (buffer,"WeightedCharacterDifferences at site %d\n", sites);
				BufferToConsole (buffer);		
			}
				
			theTask->tree->WeightedCharacterDifferences (siteLikelihood, &res1, &res2, theTask->threadIndex);
			
			_SimpleList * dSites = (_SimpleList*)(*theTask->dupList)(sites);
			StateCounterResultHandler (*theTask->fla,dSites, *theTask->doneSites, *theTask->lastDone, theTask->totalUniqueSites, res1, res2);
		}
		return nil;
	}
	
#endif

//_______________________________________________________________________________________

void	StateCounterResultHandler (_Formula& fla, _SimpleList* dSites, long & doneSites, long & lastDone, long totalUniqueSites, _Matrix& res1, _Matrix& res2)
{
	#ifdef __MP__
		if (systemCPUCount>1)
			pthread_mutex_lock(&wancMutex);
	#endif

	setParameter (stateCountMatrix,  &res1);
	setParameter (wStateCountMatrix, &res2);
	
	for (long dupSites = 0; dupSites < dSites->lLength; dupSites++)
	{
		_Operation tempO (new _Constant(dSites->lData[dupSites]));
		fla.GetList().InsertElement(&tempO,1,true);
		fla.Compute();
		fla.GetList().Delete(1);
	}
	doneSites ++;
	if (((doneSites-lastDone)*100.)/totalUniqueSites>1.)
	{
		lastDone = doneSites;
		#if !defined __UNIX__ || defined __HEADLESS__
			SetStatusBarValue ((doneSites*100.)/totalUniqueSites,1.,0.);
		#endif
		#ifdef __MAC__
			handleGUI(true);
		#endif
		if (terminateExecution)
			return;
	}
	#if defined __UNIX__ && ! defined __HEADLESS__
		if (VerbosityLevel()==1)
			UpdateOptimizationStatus (doneSites, (doneSites*100.)/totalUniqueSites, 1, false);
	#endif


	#ifdef __MP__
		if (systemCPUCount>1)
			pthread_mutex_unlock(&wancMutex);
	#endif
}

//_______________________________________________________________________________________

void	_LikelihoodFunction::StateCounter (long functionCallback)
{
	PrepareToCompute();
	computationalResults.Clear();
	
	_Operation functionCallbackOp;
	
	functionCallbackOp.SetTerms (-functionCallback-1);
	functionCallbackOp.TheCode() = functionCallback;
	
	_Formula		   fla;
	
	fla.GetList() && & functionCallbackOp;
	
	long  totalUniqueSites = 0,
		  doneSites = 0,
		  lastDone	= 0;
    {
	for (long i = 0; i<theTrees.lLength; i++)
	{
		_DataSetFilter * dsf  = (_DataSetFilter*)dataSetFilterList (theDataFilters(i));
		totalUniqueSites += dsf->NumberDistinctSites();
	}
	}
	for (long i = 0; i<theTrees.lLength; i++)
	{
		long 			 offset = -1;
		
		_TheTree       * tree = (_TheTree*)LocateVar(theTrees(i));
		_DataSetFilter * dsf  = (_DataSetFilter*)dataSetFilterList (theDataFilters(i));
		
		_Matrix 			*glFreqs = (_Matrix*)LocateVar(theProbabilities.lData[i])->GetValue();
		SetGlobalFrequencyMatrix (glFreqs);
	
		tree->InitializeTreeFrequencies (GlobalFrequenciesMatrix);
	
		_List		   duplicateMatches;
		
		while (fla.GetList().countitems()>1)
			fla.GetList().Delete(0);
			
		{
			_Operation tempO (new _Constant(i+1));
			fla.GetList().InsertElement(&tempO,0,true);
		}
		
		for (long mapper = 0; mapper < dsf->duplicateMap.lLength; mapper++)
		{
			long matched = dsf->duplicateMap.lData[mapper];
			
			if (duplicateMatches.lLength<=matched)
			{
				_SimpleList tlist;
				duplicateMatches && & tlist;
			}
				
			(*((_SimpleList*)duplicateMatches(matched))) << mapper;
		} 
		
		_Parameter	 	  scaler = 0.;
		
		_CalcNode*		  rambler = tree->DepthWiseTraversal (true);
		
		while (!tree->IsCurrentNodeTheRoot())
		{
			_Parameter bLength = rambler->BranchLength();
			_Constant tConst (bLength);
			rambler->SetValue (&tConst);
			scaler += bLength;
			rambler = tree->DepthWiseTraversal (false);
		}
		
		rambler = tree->DepthWiseTraversal (true);
		
		while (!tree->IsCurrentNodeTheRoot())
		{
			_Parameter bLength = rambler->Value();
			_Constant tConst (bLength/scaler);
			rambler->SetValue (&tConst);
			rambler = tree->DepthWiseTraversal (false);
		}
		
		_SimpleList * dSites = (_SimpleList*)duplicateMatches(0);
		SetStatusLine 	  ("Weighted ancestor counting...Computing transition matrices.");
		#ifdef __MAC__
			handleGUI(true);
		#endif
		#if defined __UNIX__ && ! defined __HEADLESS__
			if (VerbosityLevel()==1)
				UpdateOptimizationStatus (0,0,0,true);
		#endif
		#ifdef __MP__
			long tstep = (dsf->NumberDistinctSites()-1)/systemCPUCount;
			if ((systemCPUCount>1)&&tstep)
			{
				tree->BuildTopLevelCache();
				tree->AllocateResultsCache(dsf->NumberDistinctSites());
				offset = 0;
				for (long idx = 0; idx < tree->flatCLeaves.lLength; idx++)
					((_CalcNode*)tree->flatCLeaves(idx))->theProbs[0] = idx;
				for (long idx = 0; idx < tree->flatTree.lLength; idx++)
					((_CalcNode*)tree->flatTree(idx))->theProbs[0] = idx+tree->flatCLeaves.lLength;
			}
			_Parameter siteLikelihood = tree->ReleafTreeAndCheck (dsf,0,tree->HasCache());
		#else
			_Parameter siteLikelihood = tree->ReleafTreeAndCheck (dsf,0,false);		
		#endif
		#if !defined __UNIX__ || defined __HEADLESS__
			SetStatusLine 	  ("Weighted ancestor counting...Doing the counting.");
			SetStatusBarValue (0,1.,0.);
		#endif
		#ifdef __MAC__
			handleGUI(true);
		#endif
		if (terminateExecution)
			return;

		if (1)
		{
			_Matrix	 res1(tree->GetCodeBase() , tree->GetCodeBase(), false, true),
					 res2(tree->GetCodeBase(),  tree->GetCodeBase(), false, true);
			tree->WeightedCharacterDifferences (siteLikelihood, &res1, &res2, offset);					
			StateCounterResultHandler (fla, dSites,doneSites,lastDone,totalUniqueSites,res1,res2);
		}
		
		// mp stuff begins here
		
		#ifdef __MP__
			
			pthread_t 		*countingThreads = nil;
			WancReleafTask	*wancTasks		 = nil;
			
			if (tree->HasCache())
			{
				countingThreads = new pthread_t 	 [systemCPUCount-1];
				wancTasks	    = new WancReleafTask [systemCPUCount-1];
				
				for (long tc = 0; tc<systemCPUCount-1; tc++)
				{
					wancTasks[tc].tree 	 =  tree;
					wancTasks[tc].dupList = 	&duplicateMatches;

					wancTasks[tc].startAt = 1+tstep*(tc+1);
					wancTasks[tc].endAt   = 1+tstep*(tc+2);
					wancTasks[tc].dsf = dsf;

					wancTasks[tc].totalUniqueSites = totalUniqueSites;
					wancTasks[tc].lastDone 		  = &lastDone;
					wancTasks[tc].doneSites		  = &doneSites;

					wancTasks[tc].fla			  = &fla;
					wancTasks[tc].threadIndex	  = tc+1;
					
					if (tc == systemCPUCount-2)
						wancTasks[tc].endAt = dsf->NumberDistinctSites();

					#ifdef __MACHACKMP__
					if ( pthread_create( countingThreads+tc, NULL, StateCounterMPHook,(void*)(wancTasks+tc)))
					#else
					if ( pthread_create( countingThreads+tc, NULL, StateCounterMP,(void*)(wancTasks+tc)))							
					#endif
					{
						FlagError("Failed to initialize a POSIX thread in StateCounter()");
						exit(1);
					}		
				}
				
				_Parameter		vls = 0.0;
				checkParameter (VerbosityLevelString, vls, 0.0);

				for (long sites = 1; sites < 1+tstep; sites++)
				{
					pthread_mutex_lock(&wancMutex);
					_Matrix	 res1(tree->GetCodeBase() , tree->GetCodeBase(), false, true),
							 res2(tree->GetCodeBase(),  tree->GetCodeBase(), false, true);
							 
					pthread_mutex_unlock(&wancMutex);

					dSites = (_SimpleList*)duplicateMatches(sites);					
					siteLikelihood = tree->ThreadReleafTreeCache (dsf, sites, sites-1, 0, tree->flatCLeaves.lLength-1,sites);									
					if (vls>9.99)
					{
						char   buffer[64];
						sprintf (buffer,"WeightedCharacterDifferences at site %d\n", sites);
						BufferToConsole (buffer);		
					}
					tree->WeightedCharacterDifferences (siteLikelihood, &res1, &res2,0);								
					StateCounterResultHandler (fla, dSites,doneSites,lastDone,totalUniqueSites,res1,res2);
				}
				
				for (long wc = 0; wc<systemCPUCount-1; wc++)
					if ( pthread_join ( countingThreads[wc], NULL ) ) 
					{
						FlagError("Failed to join POSIX threads in StateCounter()");
						exit(1);
					}

				tree->KillTopLevelCache();
				
				delete countingThreads;
				delete wancTasks;				
			}
			else
			{
		#endif


		for (long sites = 1; sites < dsf->theFrequencies.lLength; sites++)
		{
			dSites = (_SimpleList*)duplicateMatches(sites);
			
			// populate tips
			
			siteLikelihood = tree->ReleafTree (dsf, sites, sites-1, 0, tree->flatCLeaves.lLength-1);
							
			_Matrix	 res1(tree->GetCodeBase() , tree->GetCodeBase(), false, true),
					 res2(tree->GetCodeBase(),  tree->GetCodeBase(), false, true);
					 
			tree->WeightedCharacterDifferences (siteLikelihood, &res1, &res2);			
			
			StateCounterResultHandler (fla, dSites,doneSites,lastDone,totalUniqueSites,res1,res2);
		}
	}
	
	#ifdef __MP__
	}
	#endif

	#if !defined __UNIX__ || defined __HEADLESS__
		SetStatusBarValue (-1,1.,0.);
		SetStatusLine ("Idle");
	#endif
	#ifdef __MAC__
		handleGUI(true);
	#endif

	#if defined __UNIX__ && ! defined __HEADLESS__
		if (VerbosityLevel()==1)
			UpdateOptimizationStatus (0, 0, 2, false);
	#endif

	DoneComputing ();	
}		
//_______________________________________________________________________________________

void	_LikelihoodFunction::Simulate (_DataSet &target, _List& theExclusions, _Matrix* catValues, _Matrix* catNames, _Matrix* spawnValues, _String* storeIntermediates)
{
	// will step thru multiple trees of the project and simulate  a dataset from the likelihood function
	
	long i,
		 j,
		 k;
		 
	char catSim = 0;
		// 0 - not at all
		// 1 - discrete, 
		// 2 - continuous
				
	/*for (i=0; i<theTrees.lLength; i++)
	{
	}*/

	_SimpleList 	catCont, 
					catDiscrete, 
					catHMM, 
					HMMState;
					
	bool  			columnWise = false;
	
	if (indexCat.lLength)
	{
		checkParameter (categorySimulationMethod,categorySimMethod,2.0);
		// check to see if continuous method can actually be implemented
		catSim = (char)categorySimMethod;
		if (categorySimMethod >1.5)
		{
			for (i=0; i<indexCat.lLength; i++)
			{
				_CategoryVariable* thisC = (_CategoryVariable*)LocateVar (indexCat.lData[i]);
				if (thisC->IsHiddenMarkov())
				{
					catHMM << i;
					HMMState << 0;
				}
				else
					if (thisC->GetCumulative().IsEmpty())
					{
						_String warnMsg ("SimulateDataSet will treat category variable ");
						warnMsg = warnMsg & *thisC->GetName() & " as discrete since no cumulative dist'n is available.";
						ReportWarning(warnMsg);
						catDiscrete<<i;
						HMMState<<i;
					}
					else
						catCont<<i;
			}
		}	
		else		
			for (i=0; i<indexCat.lLength; i++)
			{
				((_CategoryVariable*)LocateVar (indexCat.lData[i]))->Refresh(true); // mod Jan 3rd 2006
				catDiscrete<<i;
			}
				
		if ((catNames)&&(indexCat.lLength))
		{
			catNames->Clear();
			CreateMatrix (catNames,indexCat.lLength,1,false,true,false);
			catNames->Convert2Formulas();
			_String   thisFla;
			for (i=0; i<catHMM.lLength; i++)
			{
				thisFla = _String ("\"") & *LocateVar(indexCat.lData[catHMM.lData[i]])->GetName() & _String ("\"");
				_Formula f (thisFla);
				catNames->StoreFormula (i,0,f);
			}
			for (i=0; i<catDiscrete.lLength; i++)
			{
				thisFla = _String ("\"") & *LocateVar(indexCat.lData[catDiscrete.lData[i]])->GetName() & _String ("\"");
				_Formula f (thisFla);
				catNames->StoreFormula (i+catHMM.lLength,0,f);
			}
			for (i=0; i<catCont.lLength; i++)
			{
				thisFla = _String ("\"") & *LocateVar(indexCat.lData[catCont.lData[i]])->GetName() & _String ("\"");
				_Formula f (thisFla);
				catNames->StoreFormula (i+catHMM.lLength+catDiscrete.lLength,0,f);
			}
		}
	}
	/*else
	{
		Compute(); 
	}*/
	
	_DataSetFilter *dsf = (_DataSetFilter*)dataSetFilterList (theDataFilters(0));
	//_String duh ("Simulated Species #");
	
	i = dsf->NumberSpecies(); 
	
	target.GetNames ().RequestSpace (i);
	
	{
		_List 		* dsfNames = &dsf->GetData()->GetNames();
		_SimpleList * mapList  = (_SimpleList*)dsf->GetMap();
		
		for (long i2=0; i2< i; i2++)
			target.GetNames() <<  (*dsfNames)(mapList->lData[i2]);
	}
		
	long	countIntermediates = target.GetNames().lLength;

	if (storeIntermediates && storeIntermediates->sLength == 0)
	{
		_TheTree *datree = (_TheTree*)LocateVar(theTrees(0));
		datree->AddNodeNamesToDS (&target,false,true,false);
		countIntermediates = target.GetNames().lLength - countIntermediates;
	}
	else
		countIntermediates = 0;
		
	target.SetTranslationTable (dsf->GetData());
	
	_List* 		localExc;
	_SimpleList localExclusions;
	long   DSOffset 	  = 0,
		   siteOffset	  = 0,
		   totalPositions = 0;
	
	for (i = 0; i<theTrees.lLength; i++)
	{
		dsf = (_DataSetFilter*)dataSetFilterList (theDataFilters(i));
		totalPositions += dsf->theOriginalOrder.lLength/dsf->GetUnitLength();
	}

	if ((catValues)&&(indexCat.lLength))
	{
		catValues->Clear();
		CreateMatrix (catValues,indexCat.lLength,totalPositions,false,true,false);
	}
	
	_Constant lastTime (0.0),
			  currentTime (0.0);
	_PMathObj tValue = lastTime.Time();
	lastTime.SetValue (tValue->Value());
	DeleteObject (tValue);
	
	for (i = 0; i<theTrees.lLength; i++) // loop thru the trees one at a time
	{
		// will propagate a complete column thru the tree
		// first generate a root sequence
		dsf = (_DataSetFilter*)dataSetFilterList (theDataFilters(i));
		if (dsf->NumberSpecies()+countIntermediates!=target.GetNames().lLength) // incompatible likelihood function
		{
			ReportWarning ((_String("DataSet Simulation had to ignore part ")&_String(i)&" of the likelihood function since it has a different length than the first part."));
			continue;
		}
		_Parameter *freqs = ((_Matrix*)LocateVar(theProbabilities(i))->Compute())->fastIndex();
		_TheTree *tree = (_TheTree*)LocateVar(theTrees(i));
		// how large will the propagating vector be?
		long	vecSize = 0,leafCount = 0, good = 0;
		
		localExclusions.Clear();
		
		vecSize = dsf->theOriginalOrder.lLength/dsf->GetUnitLength();
			
		if (theExclusions.lLength>i)
		{
			localExc = (_List*)theExclusions(i);
			if (localExc->lLength)
			{
				_Parameter* translatedState = new _Parameter [dsf->GetDimension(false)];
				checkPointer(translatedState);
				for (j=0; j<localExc->lLength; j++)
				{
					dsf->Translate2Frequencies(*(_String*)(*localExc)(j),translatedState,true);
					k = 0;
					long l = -1;
					
					for (long h=dsf->GetDimension(false)-1; h>=0; h--)
						if (translatedState[h]!=0)
						{
							k++;
							l = h;
						}
					
					if (k != 1)
					{
						_String warnMsg ("Excluded character ");
						warnMsg = warnMsg & *(_String*)(*localExc)(j) &" does not represent a unique state and will therefore be ignored in data simulation";
						ReportWarning (warnMsg);
					}				
					
					localExclusions<<l;
				}
				delete (translatedState);
				columnWise = true;
			}
		}
		
		if (catSim != 0) 
			columnWise = true;
		
		if (columnWise)
		{
			_TheTree * cT = ((_TheTree*)(LocateVar(theTrees(i))));
			cT->SetUpMatrices(1);
			while (good<vecSize)
			{
				_Parameter randVal = genrand_real2(), 
						   sumSoFar = 0.0;
						   
				long rootState = 0;
				
				if (catSim)
				{				
					while (randVal<1.0e-30)
						randVal = genrand_real2();
						
					for (j=0; j<catHMM.lLength; j++) // deal with HMM
					{
						_CategoryVariable* thisC = (_CategoryVariable*)LocateVar(indexCat.lData[catHMM.lData[j]]);
						_Parameter* iWeights;
						
						if (good==0)
							iWeights = thisC->GetWeights()->fastIndex();
						else
						{
							_Matrix * hmm = thisC->ComputeHiddenMarkov();
							iWeights = hmm->theData+hmm->GetVDim()*HMMState.lData[j];
						}
						
						while (sumSoFar<randVal)
							sumSoFar+=iWeights[rootState++];

						if (rootState)
							rootState--;
						thisC->SetIntervalValue(rootState);
						if (good>0)
							HMMState.lData[j] = rootState;
													
						randVal = genrand_real2();
						rootState = 0;
						sumSoFar = 0.0;
						
						if (randVal==0.0)
						{	
							j--;
							continue;
						}
						if (catValues)
							catValues->Store (j,siteOffset+good,thisC->Compute()->Value());
					}

					for (j=0; j<catDiscrete.lLength; j++) // use discrete values here
					{
						_CategoryVariable* thisC = (_CategoryVariable*)LocateVar(indexCat.lData[catDiscrete.lData[j]]);
						_Parameter* iWeights = thisC->GetWeights()->fastIndex();
						
						while (sumSoFar<randVal)
							sumSoFar+=iWeights[rootState++];
						
						if (rootState)
							rootState--;
						thisC->SetIntervalValue(rootState);
						randVal = genrand_real2();
						rootState = 0;
						sumSoFar = 0.0;
						if (randVal==0.0)
						{	
							j--;
							continue;
						}
						if (catValues)
							catValues->Store (j+catHMM.lLength,siteOffset+good,thisC->Compute()->Value());
					}
					
					for (j=0; j<catCont.lLength; j++) // use continuous values here
					{
						_CategoryVariable* thisC = (_CategoryVariable*)LocateVar(indexCat.lData[catCont.lData[j]]);
						randVal = thisC->GetCumulative().Newton(thisC->GetDensity(),randVal,thisC->GetMinX(),thisC->GetMaxX(),_x_);
						_Constant randCat (randVal);
						thisC->SetValue(&randCat);
						if (catValues)
							catValues->Store (j+catHMM.lLength+catDiscrete.lLength,siteOffset+good,randVal);
						randVal = 0.0;
						while (randVal<1.0e-30)
							randVal = genrand_real2();
						
					}
				}
				
				if (randVal==0.0)
					continue;
				
				if (spawnValues)
					rootState = spawnValues->theData[siteOffset+good];
				else
				{
					k=0;
					while (randVal>sumSoFar) 
					{
						sumSoFar+=freqs[k];
						k++;
					}
					if (k==0)
						rootState = 0;
					else
						rootState=k-1;
				}
					
				_SimpleList intermediates,
							thisSite;
				
				if (SingleBuildLeafProbs	(tree->GetRoot(), rootState, thisSite, localExclusions,tree, true, dsf, storeIntermediates?&intermediates:nil))
				{
					good++;
					//add this site to the simulated dataset
					rootState = dsf->GetUnitLength();
					
					_String topState(dsf->ConvertCodeToLetters (dsf->CorrectCode(thisSite.lData[0]),rootState));
					for (k=0; k<topState.sLength; k++)
					{
						target.AddSite(topState.sData[k]);
						leafCount++;
					}
					for (j=1; j<dsf->NumberSpecies();j++)
					{
						topState = dsf->ConvertCodeToLetters (dsf->CorrectCode(thisSite.lData[j]),rootState);
						for (k=0; k<topState.sLength; k++)
							target.Write2Site(DSOffset+leafCount-rootState+k,topState.sData[k]);
					}
					
					target.ResetIHelper();
					for (k=0; k<topState.sLength; k++)
						target.Compact(DSOffset+leafCount-rootState+k);
					
					if (storeIntermediates)
					{
						for (j=0; j<intermediates.lLength;j++)
						{
							topState = dsf->ConvertCodeToLetters (dsf->CorrectCode(intermediates.lData[j]),rootState);
							for (k=0; k<topState.sLength; k++)
								target.Write2Site(DSOffset+leafCount-rootState+k,topState.sData[k]);
							target.ResetIHelper();
							for (k=0; k<topState.sLength; k++)
								target.Compact(DSOffset+leafCount-rootState+k);
						}
					}
					//for (k=0; k<rootState; k++)
						//target.CheckMapping (DSOffset+leafCount-rootState+k);
					
				}
				currentTime.SetValue (0.0);
				tValue = currentTime.Time();
				currentTime.SetValue (tValue->Value());
				DeleteObject (tValue);					
				if (currentTime.Value()-lastTime.Value () > .25)
				{
					#if !defined __UNIX__ || defined __HEADLESS__
						SetStatusBarValue (100.*(DSOffset+good)/totalPositions, 1, 0);
					#endif
					lastTime.SetValue (currentTime.Value());
					#ifndef __UNIX__
						yieldCPUTime();
					#endif
				}
			}
			DSOffset   += dsf->theOriginalOrder.lLength;
			siteOffset += dsf->theOriginalOrder.lLength/dsf->GetUnitLength();
			cT->CleanUpMatrices();
			continue;
		}
		
		long*   baseVector = (long*)MemAllocate (sizeof(long)*(vecSize));
		
		// generate a random "spawning vector"
		
		if (spawnValues) // use supplied starting values
			for (j=0;  j<vecSize; j++)
				baseVector[j] = spawnValues->theData[j+DSOffset];
				
		else
			for (j=0;  j<vecSize; j++)
			{
				_Parameter randVal  = genrand_real2(), 
						   sumSoFar = 0;
						   
				if (randVal==0.0)
				{
					j--;
					continue;
				}
				k=0;
				while (randVal>sumSoFar) 
				{
					sumSoFar+=freqs[k];
					k++;
				}
				if (k==0)
					baseVector[j]=0;
				else
					baseVector[j]=k-1;
			}
					
		// now proceed down the tree branches to get the values of the probabilities at the leaves
		// this is done recursively
		
		_DataSet * intermediates = nil;
		
		if (storeIntermediates)
		{
			if (storeIntermediates->sLength)
			{
				FILE * iff = doFileOpen (storeIntermediates->sData,"w");
				if (!iff)
				{
					_String errMsg = _String ("Failed to open ") & *storeIntermediates & " for writing.";
					WarnError (errMsg);
					target.Finalize();
					return;
				}
				else
				{
					intermediates = new _DataSet (iff);
					_TheTree *datree = (_TheTree*)LocateVar(theTrees(0));
					datree->AddNodeNamesToDS (intermediates,false,true,false);
				}
			}	
			else
				intermediates = new _DataSet (vecSize);

			checkPointer (intermediates);
			intermediates->SetTranslationTable (dsf->GetData());
		}
		
		BuildLeafProbs (tree->GetRoot(), baseVector, vecSize, target, tree, leafCount, true, dsf->GetUnitLength(), dsf, DSOffset,intermediates);
		
		if (intermediates)
		{
			intermediates->Finalize();
			if (storeIntermediates->sLength == 0)
			{
				long	  siteCount = intermediates->GetTheMap().lLength;
				
				for (long p = 0; p<countIntermediates; p++)
					for (long c=DSOffset; c<siteCount+DSOffset; c++)
					{
						_Site*  aSite = intermediates->GetSite(c-DSOffset);
						target.Write2Site (c,aSite->sData[p]);
					}
			}	
			DeleteObject (intermediates);
		}
		
		free (baseVector);
		DSOffset += dsf->theOriginalOrder.lLength;
		
	}
	
	/*for (i=0; i<theTrees.lLength; i++)
	{
		_TheTree * cT = ((_TheTree*)(LocateVar(theTrees(i))));
		cT->CleanUpMatrices();
	}*/
	
	target.Finalize();
	target.SetNoSpecies(target.GetNames().lLength);
}

//_______________________________________________________________________________________

void	_LikelihoodFunction::BuildLeafProbs (node<long>& curNode, long* baseVector, long& vecSize, _DataSet& target, _TheTree* curTree, long& leafCount, bool isRoot, long baseLength, _DataSetFilter* dsf, long DSOffset, _DataSet* intNodes)
{
	long* curVector = nil,
		  i,
		  k,
		  m;
		  
	_CalcNode* ccurNode = (_CalcNode*)LocateVar (curNode.get_data());
	
	if (!isRoot)
	{	
	
		curVector = (long*)MemAllocate (vecSize*sizeof(long));
		// first "mutate" the parent vector 	

		if (ccurNode->NeedToExponentiate(-1))
			ccurNode->RecomputeMatrix(0,1);
			
		_Parameter * baseI = ccurNode->GetCompExp()->fastIndex();
		
		m = ccurNode->GetCompExp()->GetVDim();
	
		for (i = 0; i<vecSize; i++)
		{
			_Parameter randVal = genrand_int32()/(_Parameter)RAND_MAX_32, 
					   sumSoFar = 0,
					  *fastI = baseI + baseVector[i]*m;
			k=0;
			while ((randVal>sumSoFar)&&(k<m))
			{
				sumSoFar+=fastI[k];
				k++;
			}
			if (k==0)
				curVector[i]=0;
			else
				curVector[i]=k-1;	
		}
	}
	else
	{
		// handle the degenerate tree case
		if (curNode.nodes.length == 1)
		{
			for (k = 0; k<vecSize; k++)
			{
				_String letterValue = dsf->ConvertCodeToLetters (dsf->CorrectCode(baseVector[k]), baseLength);
				for (m = 0; m<letterValue.sLength; m++)
					target.AddSite (letterValue.sData[m]);
			}
			leafCount++;
			BuildLeafProbs (*curNode.go_down(1), baseVector, vecSize, target, curTree, leafCount, false, baseLength, dsf,DSOffset,intNodes);
			return;
		}
	}
	
	// now scan the "children" and pass on the parameters as needed
	
	if (curNode.nodes.length) 
	{
		if (intNodes)
		{
			bool writeOrAdd = intNodes->lLength;
			for (k = 0; k<vecSize; k++)
			{
				_String letterValue = dsf->ConvertCodeToLetters (dsf->CorrectCode(curVector?curVector[k]:baseVector[k]), baseLength);
				if (writeOrAdd)
					for (m = 0; m<letterValue.sLength; m++)
						intNodes->Write2Site (letterValue.sLength*k+m, letterValue.sData[m]);
				else
					for (m = 0; m<letterValue.sLength; m++)
						intNodes->AddSite (letterValue.sData[m]);
			}
		}
		for (k = 1; k<=curNode.get_num_nodes(); k++)
			BuildLeafProbs (*curNode.go_down(k), curVector?curVector:baseVector, vecSize, target, curTree, leafCount, false, baseLength, dsf,DSOffset,intNodes);
	}
	else
	// reached a leaf
	{
		// attach a row to the new data set
		long siteCount = DSOffset;
		if (!leafCount) 
		{
			for (k = 0; k<vecSize; k++)
			{
				_String letterValue = dsf->ConvertCodeToLetters (dsf->CorrectCode(curVector[k]), baseLength);
				for (m = 0; m<letterValue.sLength; m++)
					target.AddSite (letterValue.sData[m]);
			}
			leafCount++;
		}
		else
		{
			for (k = 0; k<vecSize; k++)
			{
				_String letterValue = dsf->ConvertCodeToLetters (dsf->CorrectCode(curVector[k]), baseLength);
				for (m = 0; m<letterValue.sLength; m++)
					target.Write2Site (siteCount++, letterValue.sData[m]);
			}
		}
	}
	
	if (!isRoot)
		ccurNode->FreeUpMemory (0);
	
	if (curVector)
		free (curVector);
}

//_______________________________________________________________________________________

bool	_LikelihoodFunction::SingleBuildLeafProbs (node<long>& curNode, long parentState, _SimpleList& target, _SimpleList& theExc, _TheTree* curTree, bool isRoot, _DataSetFilter* dsf, _SimpleList * iNodes)
{
	long myState = 0;
	if (!isRoot)
	{	
	
		// first "mutate" the parent vector 	
		_CalcNode* ccurNode = (_CalcNode*)LocateVar (curNode.get_data());
					
		if (ccurNode->NeedToExponentiate(-1))
			ccurNode->RecomputeMatrix(0,1);
	
		_Parameter* fastI = ccurNode->GetCompExp()->fastIndex()+parentState*ccurNode->GetCompExp()->GetVDim(),
					randVal = genrand_int32()/(_Parameter)RAND_MAX_32, 
					sumSoFar = 0.0;
		
		long   k=0, 
			   n = ccurNode->GetCompExp()->GetVDim();
			   
		while  ((randVal>sumSoFar)&&(k<n))
		{
			sumSoFar+=fastI[k];
			k++;
		}
		
		if (k==0)
			myState = 0;
		else
			myState=k-1;
			
		if (curNode.nodes.length) 
		{
			if (iNodes)
			{
				if (theExc.Find(myState)!=-1) 
					return false;
				(*iNodes)<<myState;
				//return true;
			}
		}
		else// reached a leaf
		{
			// attach a row to the new data set
			if (theExc.Find(myState)!=-1) 
				return false;
			target<<myState;
			return true;
			
		}		
	}	
	else
	{
		if (curNode.nodes.length == 1)
			target << parentState;
		else
			if (iNodes)
				(*iNodes)<<parentState;
		
	}
			
	// now scan the "children" and pass on the parameters as needed

	for (long k = 1; k<=curNode.get_num_nodes(); k++)
	{
		if(!SingleBuildLeafProbs (*curNode.go_down(k), isRoot?parentState:myState, target, theExc,curTree, false, dsf, iNodes))
		   return false;
	}	
	return true;
}

//_______________________________________________________________________________________

_AssociativeList*	_LikelihoodFunction::SimulateCodonNeutral (_Matrix* synCost, _Matrix* nsCost, long countPerState)
{
	_AssociativeList * resList = new _AssociativeList;
	
	if (indexCat.lLength || theTrees.lLength != 1)
	{
		_String errMsg ("SimulateCodonNeutral works only with likelihood functions which do not include rate variation and contain exactly one partition.");
		WarnError (errMsg);
	}
	else
	{
		PrepareToCompute ();
		Compute();
		_TheTree * tree = (_TheTree*)LocateVar(theTrees(0));
		
		long	   stateCount = nsCost->GetVDim();
		_FString   aKey;
		
		long	   maxSubCount = 3*(tree->GetLeafCount()+tree->GetINodeCount()),
				   mxDim 	   = 1+3*((1+maxSubCount)*maxSubCount);
				   
		SetStatusLine ("Simulating the null distribution");
		
		long	  totalSimCount = stateCount * countPerState / 100,
				  doneCount		= 0;
				  
		for (long k=0; k<stateCount; k++)
		{
			/*_Matrix			   simmedSites (countPerState,2,false,true),
							       *sortedStates = nil;*/
							       
			_Matrix			   	  simmedStates (mxDim,1,false,true);
			
			for (long it = 0; it < countPerState; it++)
			{
				_Parameter s  = 0.0,
						   ns = 0.0;
						   
				
				 doneCount++;
				 #if !defined __UNIX__ || defined __HEADLESS__
				 if (doneCount%totalSimCount==0)
					SetStatusBarValue (doneCount/totalSimCount,1.0, 0.0);
				 #endif

				CodonNeutralSimulate	(tree->GetRoot(), k, true, synCost, nsCost, s, ns);

				#ifndef __BORLAND_HACK__
				long nsi = round (s+ns);
				#else
				long nsi = (long)(s+ns+0.5);
            	#endif
				if (ns <= maxSubCount)
				{
					if (nsi)
					{
						#ifndef __BORLAND_HACK__
						long si  = round(s*6.);
						#else
						long si  = (long)(s*6.+0.5);
						#endif

						simmedStates.theData[1+3*((nsi-1)*nsi)+si] += 1.0;
					}
					else
						simmedStates.theData[0] += 1.0;
				}
				
				/*simmedSites.theData[2*it]   = floor(s+ns);
				simmedSites.theData[2*it+1] = s;*/
			}
			
			/* resort the matrix */
			
			_AssociativeList 	  *stateAVL = new _AssociativeList;
			
			for (long it2 = 0; it2 < maxSubCount; it2++)
			{
				long	mxDim2 = it2?2+6*it2:2,
						offset = (it2>0)+3*((it2-1)*it2);
						
				_Matrix * conditionalDistribution = new _Matrix (mxDim2,2,false,true);
				
				_Parameter total = 0.0,
						   mb 	 = 1./6.;
				
				for (long idx2 = 0; idx2 < mxDim2-1; idx2++)
				{
					conditionalDistribution->theData[2*idx2+2] = idx2*mb;
					total += (conditionalDistribution->theData[2*idx2+3] = simmedStates.theData[offset+idx2]);
				}
				
				if (total>0.0)
				{
					conditionalDistribution->theData[0] = total;
					total = 1./total;
					conditionalDistribution->theData[3] *= total;
					for (long idx2 = 5; idx2 < 2*mxDim2; idx2+=2)
						conditionalDistribution->theData[idx2] = conditionalDistribution->theData[idx2]*total+conditionalDistribution->theData[idx2-2];

					*aKey.theString = it2;
					stateAVL->MStore (&aKey, conditionalDistribution, false);
				}
				else
					DeleteObject (conditionalDistribution);
			}
			
			*aKey.theString = k;
			resList->MStore (&aKey, stateAVL, false);
			
			/* now resort the values */
			
			/*sortedStates = (_Matrix*)simmedSites.SortMatrixOnColumn (&_Constant (0.0));
			
			_Parameter		lastV = sortedStates->theData[0];
			long			lastI = 0,
							curI  = 0,
							mD	  = 2*countPerState;
							
			while (1)
			{
				while ( curI < mD && sortedStates->theData[curI] == lastV)
					curI += 2;
				
				long 		mxD = (curI - lastI)/2;
				
				_Matrix	    conditionalSample (mxD,1,false,true),
							* scs = nil;
				
				for (long i2 = 0; i2 < mxD; i2++)
					conditionalSample.theData[i2] = sortedStates->theData[lastI+2*i2+1];
				
				
				_Parameter 		mxDI = 1./mxD;
				scs = (_Matrix*)conditionalSample.SortMatrixOnColumn (&_Constant (0.0));
				
				long	lastI2 = 0,
						curI2  = 0;
						
				_Parameter lastV2 = scs->theData[0];
				
				_AssociativeList condASL;
				
				while (1)
				{
					while (curI2 < mxD && scs->theData[curI2] == lastV2)
						curI2++;
					
					long mxD2 = curI2-lastI2;
					*aKey.theString = _String (lastV2);
					
					condASL.MStore (&aKey, new _Constant (mxD2*mxDI), false);	
					lastI2 = curI2;
					
					if (curI2 < mxD)
						lastV2 = scs->theData[curI2];
					else
						break; 
				}
				
				curI2 = condASL.avl.countitems();
				_Matrix * condProbs = new _Matrix (curI2+1, 2, false, true);
				
				condProbs->theData[0] = mxD;
				
				_SimpleList tcache;
				
				long		iv,
							k2 = condASL.avl.Traverser (tcache, iv, condASL.avl.GetRoot());
							
				lastI2 = 2;			
				for (; k2>=0; k2 = condASL.avl.Traverser (tcache, iv))
				{
					condProbs->theData[lastI2] = ((_String*)condASL.avl.Retrieve (k2))->toNum();
					condProbs->theData[lastI2+1] = ((_Constant*)condASL.avl.GetXtra (k2))->Value()+condProbs->theData[lastI2-1];
					lastI2+=2;
				}
						
				DeleteObject (scs);
				
				*aKey.theString = lastV;
				stateK->MStore (&aKey, condProbs, false);
				
				lastI = curI;
				if (curI < mD)
					lastV = sortedStates->theData[curI];
				else
					break; 
					
			}
			
			DeleteObject (sortedStates);
			*aKey.theString = k;
			
			resList->MStore (&aKey, stateK,false);*/
			
		}
		DoneComputing();
	}
	SetStatusLine ("Idle");
	return resList;
}

//_______________________________________________________________________________________

void	_LikelihoodFunction::CodonNeutralSimulate (node<long>& curNode, long parentState, bool isRoot, _Matrix* costMatrixS, _Matrix* costMatrixNS, _Parameter& synCount, _Parameter& nsCount)
{
	long myState = 0;
	
	if (!isRoot)
	{	
		_CalcNode* ccurNode = (_CalcNode*)LocateVar (curNode.get_data());
		
		_Matrix  * compExpMatrix = ccurNode->GetCompExp();
					
		long   k=0, 
			   n = compExpMatrix->GetVDim();

		_Parameter* fastI = compExpMatrix->theData+parentState*n,
					randVal = genrand_int32()/(_Parameter)RAND_MAX_32, 
					sumSoFar = 0.0;
		
			   
		while  ( randVal>sumSoFar && k<n)
		{
			sumSoFar+=fastI[k];
			k++;
		}
		
		if (k==0)
			myState = 0;
		else
			myState=k-1;
			
		n = parentState * n + myState;
		synCount += costMatrixS->theData[n];
		nsCount  += costMatrixNS->theData[n];
			
	}	
			
	// now scan the "children" and pass on the parameters as needed

	for (long k = curNode.get_num_nodes(); k ; k--)
		CodonNeutralSimulate (*curNode.go_down(k), isRoot?parentState:myState, false, costMatrixS, costMatrixNS, synCount, nsCount);
}

//_______________________________________________________________________________________


void	_LikelihoodFunction::SetNthBit (long& reference, char n)
{
	long bitshifter = 1;
	bitshifter = bitshifter<<n;
	reference = reference|bitshifter;
}			

//_______________________________________________________________________________________
	
bool	_LikelihoodFunction::CheckNthBit (long& reference, char n)
{
	unsigned long bitshifter = 1;
	bitshifter = bitshifter<<n;
	return reference&bitshifter;
}	

//_______________________________________________________________________________________
	
char	_LikelihoodFunction::HighestBit (long reference)
{
	unsigned long bitshifter = 1,count = sizeof(long)*8-1;
	bitshifter = bitshifter<<(sizeof(long)*8-1);
	while(!(reference&bitshifter))
	{
		bitshifter = bitshifter>>1;
		count--;
	}
	return count;
}

//_______________________________________________________________________________________
	
long	_LikelihoodFunction::HasHiddenMarkov (long reference, bool hmm)
{
	unsigned long bitshifter = 1, count = sizeof(long)*8-1;
	long	 hMarkov = -1;
	bitshifter = bitshifter<<(sizeof(long)*8-1);
	while(bitshifter)
	{
		if (bitshifter&reference)
		{
			_CategoryVariable* thisC = (_CategoryVariable*)LocateVar(indexCat.lData[count]);
			if (hmm)
			{
				if (thisC->IsHiddenMarkov ())
					hMarkov = count;
			}
			else
				if (thisC->IsConstantOnPartition ())
					return count;
		}
		bitshifter = bitshifter>>1;
		count--;
	}
	return hMarkov;
}

//_______________________________________________________________________________________
	
char	_LikelihoodFunction::LowestBit (long reference)
{
	unsigned long bitshifter = 1,count = 0;
	while(!(reference&bitshifter))
	{
		bitshifter = bitshifter<<1;
		count++;
	}
	return count;
}

//_______________________________________________________________________________________
	
void	_LikelihoodFunction::BuildIncrements (long reference, _SimpleList& incList)
{
	unsigned long currentIncr = 1;
	for (long i=0; i<indexCat.lLength; i++)
	{
		if (CheckNthBit(reference,i))
		{
			incList<<currentIncr;
			currentIncr*=((_CategoryVariable*)LocateVar(indexCat(i)))->GetNumberOfIntervals();
		}
		else 	
			incList<<0;
	}
}
	
//_______________________________________________________________________________________
	
long	_LikelihoodFunction::MaximumDimension (void)
{
	long maxDim = 0;
	for (long i=0; i<theTrees.lLength; i++)
	{
		_Matrix *cM = (_Matrix*)LocateVar(theProbabilities.lData[i])->GetValue();
		long myDim = cM->GetHDim()>cM->GetVDim()?cM->GetHDim():cM->GetVDim();
		if (myDim>maxDim)
			maxDim = myDim;
	}
	return maxDim;
}		

//_______________________________________________________________________________________
	
long	_LikelihoodFunction::SequenceCount (long partID)
{
	if ((partID>=0)&&(partID<theTrees.lLength))
	{
		_TheTree * cT = ((_TheTree*)(LocateVar(theTrees(partID))));
		_PMathObj  seqs = cT->TipCount();
		long res = seqs->Value();
		DeleteObject (seqs);
		return res;
	}
	return -1;
}	

//_______________________________________________________________________________________
	
long	_LikelihoodFunction::SiteCount (void)
{
	long res = 0;
	
	for (long i=0; i<theDataFilters.lLength; i++)
	{
		_DataSetFilter * df = (_DataSetFilter*) dataSetFilterList (theDataFilters.lData[i]);
		res += df->theOriginalOrder.lLength;
	}

	return res;
}		

//_______________________________________________________________________________________

	
	
//_______________________________________________________________________________________
	
void	_LikelihoodFunction::PrepareToCompute (bool disableClear)
{
	if (hasBeenSetUp == 0)
	{
		long categCount = 1;
								
		for (long i=0; i<theTrees.lLength; i++)
		{
			_TheTree			* cT = ((_TheTree*)(LocateVar(theTrees(i))));
						
			long locCC = cT->CountTreeCategories();
			if   (locCC >categCount) categCount = locCC;
			cT->SetUpMatrices (locCC);
			
			#if USE_SCALING_TO_FIX_UNDERFLOW
				cT->AllocateUnderflowScalers (((_DataSetFilter*)dataSetFilterList (theDataFilters(i)))->NumberDistinctSites());
			#endif
			
		}
		
		for (long i2=0; i2<theProbabilities.lLength; i2++)
			((_Matrix*)LocateVar(theProbabilities.lData[i2])->GetValue())->MakeMeSimple();
#ifdef _SLKP_LFENGINE_REWRITE_
		SetupLFCaches	  ();
#endif
		SetReferenceNodes ();
		
		if (disableClear)
			hasBeenSetUp = 0x6FFFFFFF;
		else
			hasBeenSetUp ++;
		siteArrayPopulated = false;
	}	
	else
		hasBeenSetUp++;
}		

//_______________________________________________________________________________________
	
void	_LikelihoodFunction::DoneComputing (bool force)
{
	if (hasBeenSetUp == 1 || (hasBeenSetUp > 0 && force))
	{
		for (long i=0; i<theTrees.lLength; i++)
		{
			_TheTree * cT = ((_TheTree*)(LocateVar(theTrees(i))));
			cT->CleanUpMatrices();
			#if USE_SCALING_TO_FIX_UNDERFLOW
				cT->DeallocateUnderflowScalers ();
			#endif
		}
		if (mstCache)
		{
			mstCache->resultCache.Clear();
			mstCache->statesCache.Clear();
		}
		{
		for (long i=0; i<theProbabilities.lLength; i++)
			((_Matrix*)LocateVar(theProbabilities.lData[i])->GetValue())->MakeMeGeneral();
        }
		DeleteCaches		(false);
		hasBeenSetUp 	   = 0;
		siteArrayPopulated = false;
	}
	else
		hasBeenSetUp --;
}

//_______________________________________________________________________________________
void	_LikelihoodFunction::RankVariables(void)
{
	_SimpleList	varRank, holder;
	long		k,j,f;
	
	for (k=0; k<indexInd.lLength; k++)
	{
		if (LocateVar(indexInd.lData[k])->IsGlobal())
		{
			varRank<<10000;
		}
		else
		{
			varRank<<10050;
		}
	}

	for (k=0; k<indexDep.lLength;k++)
	{
		holder.Clear();
		{
			_AVLList   al (&holder);
			LocateVar (indexDep.lData[k])->ScanForVariables(al,true);
			al.ReorderList ();
		}
		for (j=0;j<holder.lLength; j++)
		{
			f = indexInd.Find(holder.lData[j]);
			if (f>=0)
			{
				varRank.lData[f]--;
			}
		}
	}
	SortLists (&varRank,&indexInd);
}

//_______________________________________________________________________________________

_CustomFunction::_CustomFunction (_String* arg)
{
	_String	body (*arg);
	Parse (&myBody, body,nil,nil,true);
	_SimpleList myVars;
	{
		_AVLList al (&myVars);
		myBody.ScanFForVariables(al,true,false,false);
		al.ReorderList();
	}
	for (long k=0; k<myVars.lLength; k++)
		if (LocateVar(myVars.lData[k])->IsIndependent())
			GetIndependentVars() << myVars.lData[k];
}

//_______________________________________________________________________________________

_Parameter _CustomFunction::Compute (void)
{
	likeFuncEvalCallCount++;
	_SimpleList * iv = &GetIndependentVars ();
	for (long i=0;i<iv->lLength;i++)
	{
		_Variable* cornholio = LocateVar(iv->lData[i]);
		_Parameter result = GetIthIndependent(i);
		
		if (result<cornholio->GetLowerBound() || result>cornholio->GetUpperBound())
			return -A_LARGE_NUMBER;
	}

	_PMathObj res = myBody.Compute();
	if (res)
		return res->Value();
	return 0.0;	
}