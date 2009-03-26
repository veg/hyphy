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

#ifndef __LIKELIHOODF__
#define	__LIKELIHOODF__

//#pragma once
#include "category.h"
#include "calcnode.h"

#ifdef __HYALTIVEC__
#define	  A_LARGE_NUMBER		  1.e35
#else
#define	  A_LARGE_NUMBER		  1.e100
#endif

#define	  DEFAULTPARAMETERLBOUND  0.0
#define   DEFAULTPARAMETERUBOUND  10000.0

/* various run modes for PopulateConditionalProbabilities */

#define	  _hyphyLFConditionProbsRawMatrixMode		0
#define	  _hyphyLFConditionProbsScaledMatrixMode	1
#define	  _hyphyLFConditionProbsWeightedSum			2
#define	  _hyphyLFConditionProbsMaxProbClass		3
#define	  _hyphyLFConditionProbsClassWeights		4

/* computational template kinds for the likelihood function */

#define	  _hyphyLFComputationalTemplateNone			0
#define	  _hyphyLFComputationalTemplateBySite		1
#define	  _hyphyLFComputationalTemplateByPartition  2

//_______________________________________________________________________________________

struct	MSTCache
{
	_List				computingOrder,
						storageOrder,
						referenceOrder,
						parentOrder,
						stashedLeafOrders;
						
	_SimpleList			statesNCache,
						resultCache,
						statesCache,
						cacheSize;
	
};

//_______________________________________________________________________________________

class	_LikelihoodFunction: public BaseObj {

	public:
	
		// constructors
		
		_LikelihoodFunction (void); // default - doesn't do much
		_LikelihoodFunction (_String&, _VariableContainer*); // from triplets
		_LikelihoodFunction (_LikelihoodFunction&); // stack copy
		void		Init	  (void);
		
		bool		Construct (_List&,_VariableContainer*);
		
virtual ~_LikelihoodFunction (void) 
					{Cleanup();} // destructor

virtual	BaseRef 	toStr (void);

virtual	BaseRef		makeDynamic (void); 	 // dynamic copy of this object

virtual	void		Duplicate (BaseRef);		 // duplicate an object into this one

		_SimpleList&GetIndependentVars (void); // return a list of all indepenent variables 
		_SimpleList&GetDependentVars   (void); // get all dependent vars of this object	
		_SimpleList&GetCategoryVars	   (void); // get all category variables
		void		GetGlobalVars	   (_SimpleList&);

		long		CountObjects 	   (char);
												   // return object count
												   // 0 - partitions
												   // 1 - global variables
												   // 2 - local independents
												   // 3 - dependents
		
		
		_Parameter	GetIthIndependent (long); 	// get the value of i-th independent variable 
		_Parameter	GetIthDependent   (long); 	// get the value of i-th dependent variable
		
		void		SetIthIndependent (long, _Parameter); 	// set the value of i-th independent variable 
		bool		CheckAndSetIthIndependent (long, _Parameter); 	// sget the value of i-th independent variable 
		void		SetIthDependent (long, _Parameter); 	// set the value of i-th dependent variable 
		
		void		UpdateIndependent (long,bool);
		void		UpdateDependent (long);
		_AssociativeList*	
					SimulateCodonNeutral (_Matrix*, _Matrix*, long);
		
		bool		PreCompute		(void);
		void		PostCompute		(void);
virtual		
		_Parameter	Compute 		(void);
		
		void		PrepareToCompute (bool = false);
		void		DoneComputing 	 (bool = false);
virtual		
		_Matrix*	Optimize ();
		_Matrix* 	ConstructCategoryMatrix		(const _SimpleList&, char, bool = true, _String* = nil);
		
		_Parameter	SimplexMethod				(_Parameter& precision);
		void		Anneal						(_Parameter& precision);
			
		void		Simulate					(_DataSet &,_List&, _Matrix* = nil, _Matrix* = nil, _Matrix* = nil, _String* = nil);
	
		void		ReconstructAncestors		(_DataSet &, _SimpleList&, _String&, bool = false, bool = false);
					// 20090224: added an argument to allow 
					// the marginal state reconstruction
	
		long		MaximumDimension			(void);
		
virtual	_PMathObj	CovarianceMatrix			(_SimpleList* = nil);
		
		// compute  covariance matrix  based on the Hessian
		// optional list of parameters to estimate the conditional covariance for


		void		RescanAllVariables		(void);
													
		long		DependOnTree	   		(_String&);
		long		DependOnModel	   		(_String&);
		long		DependOnDS	   	   		(long);
		bool		DependOnDF	       		(long ID) {return theDataFilters.Find(ID)>=0;}
		bool		MapTreeTipsToData  		(long, bool leafScan = false);
		bool		UpdateFilterSize   		(long);
		void		VoidOldResults 	   		(void) {computationalResults.ZeroUsed();}
		_CategoryVariable*		
					FindCategoryVar	  		(long);
					// return the category variable for a given partition 
		void		RankVariables 			(void);
		_SimpleList&GetTheTrees   			(void) { return theTrees;}
		_SimpleList&GetTheFilters 			(void) { return theDataFilters;}
		_SimpleList&GetBaseFreqs 			(void) { return theProbabilities;}
		void		Setup 		  			(void); 
		bool&		HasBeenOptimized (void) { return hasBeenOptimized; }
		void		ComputePruningEfficiency(long&, long&);
		
		long		SequenceCount			(long);
		long		SiteCount				(void);
		void		Rebuild					(void);
		void		SerializeLF				(_String&, char=0, _SimpleList* = nil, _SimpleList* = nil);
		_Formula*	HasComputingTemplate	(void) 
											{return computingTemplate;}
		void		StateCounter 		  	(long);
		void		MPI_LF_Compute 			(long, bool = false);
												
#if defined	_SLKP_LFENGINE_REWRITE_ 
	#if defined _OPENMP
		void		SetThreadCount			  (long tc) { lfThreadCount = tc;}
		long		GetThreadCount			  (void) { return lfThreadCount;}
	#else
		long		GetThreadCount			  (void) { return 1;}
	#endif
#endif 

		bool			ProcessPartitionList		(_SimpleList&, _Matrix*, _String);	
						// given a matrix argument (argument 2; can be nil to include all)
						// populate a sorted list (argument 1)
						// of partitions indexed in the matrix (e.g. {{1,3}} would include the 2-nd and 4-th partitions
						// the third argument is the ID of the calling function for error reporting
						// returns true if at least one of the requested partitions is valid
						// otherwise returns false
	
		long			TotalRateClassesForAPartition (long);
						// given a partition index (assumping that category caches have been set up
						// returns how many total rate classes there are for this partition
		
	
	
protected:
	
		_Matrix* 		PairwiseDistances 		(long index);
		void			CheckDependentBounds	(void);
		void			AllocateSiteResults 	(void);
		void			ZeroSiteResults			(void);
		// 20090211: A utility function to reset site results.
		void			GetInitialValues 		(void);
		bool			checkPermissibility 	(_Matrix&m, long row);
		_Parameter		computeAtAPoint 		(_Matrix&m, long row = 0);
		_Parameter		replaceAPoint 			(_Matrix&m, long row, _Matrix&p, _Parameter& nV, _Matrix& fv);
		
		void			ScanAllVariablesOnPartition 
												(_SimpleList&, _SimpleList&, _SimpleList&, _SimpleList&);  
							// internal function to scan all the on a set of partitions variables in

virtual	void			ScanAllVariables 		(void);  
							// internal function to scan all the variables in
		
		void			OptimalOrder	 		(long, _SimpleList&);	
							// determine the optimal order of compuation for a block
		
		_Parameter		ComputeBlock 			(long, _Parameter* siteResults = nil, long currentRateClass = -1, long = -1, _SimpleList* = nil);
		// 20090224: SLKP
		// added the option to pass an interior branch (referenced by the 3rd argument in the same order as flatTree)
		// and a set of values for each site pattern (indexed left to right) in the 4th argument
	
		void			SetReferenceNodes		(void); 
							// compute likelihood over block index i
		
							// internal Setup, called from within constructors
		
//		void			DumpingOrder	 		(long index, _SimpleList& sl);
							// used for FreeUpMemory
		
		void			Cleanup 				(void); 
							// internal cleanup, called from within destructors 
							
		void			Clear					(void);
		

		long 			BracketOneVar 		  (long , _Parameter& , _Parameter& , _Parameter& ,  
												_Parameter& , _Parameter& , _Parameter& , _Parameter&);
		long 			GradientBracketOneVar (_Matrix&, _Matrix& , _Matrix& , _Matrix&,  _Parameter& ,
											 	_Parameter&, _Parameter&, _Parameter&, bool retry = false);
		void			LocateTheBump 		  (long,_Parameter , _Parameter& , _Parameter&);
		void			GradientLocateTheBump (_Parameter, _Parameter&, _Matrix&, _Matrix&);
		void			GradientDescent 	  (_Parameter& , _Matrix& );
		void			ConjugateGradientDescent 
											  (_Parameter , _Matrix& );

		long			CostOfPath	 		  (_DataSetFilter*, _TheTree* , _SimpleList&, _SimpleList* = nil);

		void			BuildLeafProbs 		  (node<long>& , long*, long&, _DataSet&, _TheTree*, long&, bool, long, _DataSetFilter*, long, _DataSet* = nil);
		bool			SingleBuildLeafProbs  (node<long>&, long, _SimpleList&, _SimpleList&, _TheTree*, bool,_DataSetFilter*, _SimpleList* = nil);
		void			CodonNeutralSimulate  (node<long>&, long, bool,_Matrix*,_Matrix*, _Parameter&, _Parameter&);		
					    
		bool			HasBlockChanged		  (long);
		long			BlockLength			  (long);
		void			PartitionCatVars	  (_SimpleList&, long);
		// 20090210: extract variable indices for category variables in i-th partition 
		// and append them to _SimpleList
	
		
static	void			RandomizeList				(_SimpleList&, long);		
static	void			CheckFibonacci				(_Parameter);
	
		long			PartitionLengths		    (char = 0);
						/* 
							SLKP: 20090317 
							
							argument = 0 
								return the length (sites) of the longest part in this
								likelihood function
						    argument = 1
								return the cumulative length (sites) of all parts in this
								likelihood function
						 
						*/

		
	private: 	
	
		void	  		SendOffToMPI		  		(long);
		void			RecurseCategoryMPI 			(long,  long, _Parameter);
		_Parameter		ComputeMasterMPI			(void);
		void			InitMPIOptimizer			(void);
		void			CleanupMPIOptimizer			(void);
		void			ComputeBlockInt1 			(long,_Parameter&,_TheTree*,_DataSetFilter*, char);
		void		 	CheckStep 					(_Parameter&, _Matrix, _Matrix* selection = nil);
		void		 	ComputeGradient 			(_Matrix&, _Matrix&,  _Parameter&, _Matrix&, _SimpleList&, 
										 			long, bool normalize = true);
		bool		 	SniffAround 				(_Matrix& , _Parameter& , _Parameter&);
		long		 	HasPrecisionBeenAchieved 	(_Parameter funcValue = 2.*A_LARGE_NUMBER, bool = false);
		void		 	RecurseCategory				(long,long,long,long,_Parameter
#ifdef _SLKP_LFENGINE_REWRITE_
													,_SimpleList* = nil, char = 0, _Parameter* = nil,
													long = -1, _SimpleList* = nil
#endif
													 );
		void	  		RecurseConstantOnPartition  (long, long, long, long, _Parameter, _Matrix&);	

	
		
		void		 	SetNthBit 					(long&,char);
		bool		 	CheckNthBit 				(long&,char);
		void		 	BuildIncrements 			(long, _SimpleList&);
		char		 	HighestBit 					(long);
		char			LowestBit  					(long);
		long			HasHiddenMarkov				(long, bool hmm = true);
		_Matrix*	 	RemapMatrix					(_Matrix*, const _SimpleList&);
						/* given a matrix where each column corresponds to a site pattern
						 and a list of partitions that are covered
						 this function returns a new matrix where each pattern is mapped
						 to the corresponding sites in the alignment
						*/
		void			CleanUpOptimize				(void);
		void			ComputeBlockForTemplate		(long, bool = false);
		void			ComputeBlockForTemplate2	(long, _Parameter*, _Parameter*, long);
		void			DeleteCaches				(bool = true);
		void			PopulateConditionalProbabilities	
													(long index, char runMode, _Parameter* buffer, _SimpleList& scalers, long = -1, _SimpleList* = nil);
		void			ComputeSiteLikelihoodsForABlock
													(long, _Parameter*, _SimpleList&, long = -1, _SimpleList* = nil);
	
						// this function computes a list of site probabilities for the i-th block (1st parameter)
						// stores them in pattern (left to right) order (2nd argument)
						// and writes scaling factors for each site into the third argument
						// arguments 4 (branch index) and 5 (pattern assignments)
						// allows the calculation of the probability vector while setting a specific interior branch
						// to a given sequence
						
		_Parameter		SumUpSiteLikelihoods		(long, const _Parameter*, const _SimpleList&);
					    /* 
							SLKP 20090318 
							
							given a partition index (argument 1),
							a vector of _pattern_ likelihoods (argument 2) and	
							a list of pattern scaling factors (argument 3)
						 
							compute the log likelihood of the partition
						 
						 */
	
		void			UpdateBlockResult			(long, _Parameter);
						/*
							SLKP 20090318 
						 
							update cached log-likelihood value for a given 
							partition
						*/
	
	
		_List*			RecoverAncestralSequencesMarginal
													(long, _Matrix&,_List&);
		void			DetermineLocalUpdatePolicy  (void);
		void			FlushLocalUpdatePolicy		(void);
		void			RestoreScalingFactors		(long, long, long, long*, long *);
		void			SetupLFCaches				(void);
		void			SetupCategoryCaches			(void);

		_SimpleList	 	theTrees, 
						theDataFilters, 
						theProbabilities, 
						indexInd, 
						indexDep, 
						indexCat,
						*nonConstantDep,
					 	blockDependancies;

		_GrowingVector  computationalResults;
	
		_List			optimalOrders,
						leafSkips,
						categoryTraversalTemplate;
						
/*SLKP: 20090225 
						This list contains as many entries (themselves of type _List) as there are partitions
						The entry will be empty for a partition without category variables
						For a partition with N category variables the entry will contain a
							1). _List of references to category variables themselves in the order that 
								they appear in the blockDependancies
							2). _SimpleList of category counts for each variable C_1, C_2, ... C_N + a final entry with
								C_1 * C_2 * ... * C_N -- the total number of rate categories for this partition
							3). _SimpleList of incremental loop offsets for each variable, e.g. if the
								_SimpleList in 2). is 2,3,4, then this _SimpleList is
								12,4,1
*/
						
		long		 	templateKind,
	
						/* 
							_hyphyLFComputationalTemplateNone: 
								no computational template: simply return the summed up log-likelihood
								for all partitions in the likelihood function
						 
							_hyphyLFComputationalTemplateBySite:
								all partitons have the same number of sites
								the template is in terms of SITE_LIKELIHOOD
								meant to implement model/tree mixture constructs
						 
							_hyphyLFComputationalTemplateByPartition
								the template is in terms of BLOCK_LIKELIHOOD
								meant to implement random effects on the level of partitions 
								models
							
							<0											: (-1-variable index) if a hidden markov process on site likelihoods
							>_hyphyLFComputationalTemplateByPartition	: a _hyphyLFComputationalTemplateByPartition shifted index of a user
																		  function that will assemble the likelihood
																		  from site conditionals, written into matrices
																		  SITE_LIKELIHOODS_0, SITE_LIKELIHOODS_1, .. 
						*/
						hasBeenSetUp;
		
		_Matrix			*siteResults,
						*bySiteResults;
		
		bool			hasBeenOptimized,
						siteArrayPopulated;
																		
		_Formula*		computingTemplate;
		MSTCache*		mstCache;
	
#ifdef	_SLKP_LFENGINE_REWRITE_
		/* 
		   these variables store conditional likelihoods for every paritition
		   internal node (post-order traversal)
		   and leaves;
		 
		   -conditionalInternalNodeLikelihoodCaches: 3D vector
		 
				1st coordinate - the node index (in post-order traversal)
				2nd coordinate - the site (unique pattern) index (left-to-right)
				3rd coordiante - i-th marginal (for the i-th character; 0-filterCharDimension)
		 
				stores the probability for the subtree below this node given that i-th character
				is present at the node
		 
				for partitions with category variables this vector is copied enough times to 
				store every combination of category values
		
			-conditionalTerminalNodeStateFlag: 2D vector
				1st coordinate - the leaf index (in post-order traversal)
				2nd coordinate - the site (unique pattern) index (left-to-right) 
				
				stores a non-negative integer for a leaf with a simple (non-ambiguous) character
				stores a negative value to index (-value-1) into the vector of conditionalTerminalNodeLikelihoodCaches
				and read off filterCharDimension characters from there
		*/
	
		_Parameter**		conditionalInternalNodeLikelihoodCaches,
				  **		siteScalingFactors,
				  **		branchCaches;
	
		_List				conditionalTerminalNodeLikelihoodCaches;	
		long	  **		conditionalTerminalNodeStateFlag;
	
		_SimpleList			overallScalingFactors,
							overallScalingFactorsBackup,
							// this (and siteCorrectionsBackup)
							// is needed to store site/partition scaling factors when treeTraversalMasks
							// and branch caching are used, otherwise scaling done by ComputeBranchCaches
							// will fubar the scaling for ComputeTreeBlockByBranch depending on site
							// patterns
							canUseReversibleSpeedups,
							// a partition will be tagged with '1' if its tree has only
							// time-reversible models
							siteScalerBuffer;
							// used for LF with category variables to 
							// store site-by-site scaling factors
		 
		_List				localUpdatePolicy, 
							matricesToExponentiate,
							treeTraversalMasks,
							computedLocalUpdatePolicy,
							siteCorrections,
							siteCorrectionsBackup,
							cachedBranches,
							// for models with categories, a list of site by site scaling operation counts
							partScalingCache
							// used to store site by site scalers in computations that are performed
							// on a site-by-site basis; includes scratch cache for remapping
							;
	
#ifdef	_OPENMP
		long				lfThreadCount;
#endif	
	
#endif	
};

//_______________________________________________________________________________________

class	_CustomFunction: public _LikelihoodFunction {

public:

						_CustomFunction 		(_String*);
	
virtual		_Parameter	Compute 				(void);
	
	

	_Formula myBody;
};


//_______________________________________________________________________________________

extern	bool forceRecomputation,
			 isInOptimize;
			 
extern	long lockedLFID;

extern  _String	// declare shared global keys/settings
		globalStartingPoint 			,
		randomStartingPerturbations 	,
		optimizationPrecision 			,
		startingPrecision 				,
		optimizationMethod 				,
		useLastResults 					,
		allowBoundary 					,
		bracketingPersistence			,
		intermediatePrecision			,
		keepOptimalOrder 				,
		skipOmissions 					,
		optimizeSummationOrder 			,
		optimizePartitionSize 			,
		maximumIterationsPerVariable 	,
		optimizationPrecisionMethod 	,
		relativePrecision 				,
		likefuncOutput					,
		dataFilePrintFormat				, 
		dataFileDefaultWidth			,
		dataFileGapWidth				,
		categorySimulationMethod		,
		useInitialDistanceGuess			,
		randomSeed						,
		assignedSeed					,
		covariancePrecision				,
		cacheSubtrees 					,
		likeFuncCountVar 				,
		doShuffleOrder					,
		forceDistanceEstimates			,
		useDuplicateMatrixCaching		,
		siteWiseMatrix					,
		blockWiseMatrix					,
		useFullMST						,
		stateCountMatrix				,
		wStateCountMatrix				,
		tryNumericSequenceMatch			,
		allowSequenceMismatch			,
		shortMPIReturn					,
		mpiPrefixCommand				,
		skipConjugateGradient			,
		useIntervalMapping				,
		intervalMappingMethod			,
		useAdaptiveVariableStep			,
		storeRootSupportFlag			,
		supportMatrixVariable			,
		optimizationStatusFile			,
		autoParalellizeLF				;



bool 	CheckEqual 					(_Parameter, _Parameter);
void	StateCounterResultHandler 	(_Formula&, _SimpleList*,long&,long&,long,_Matrix&,_Matrix&);

_LikelihoodFunction*
		FindLikeFuncByName		 	 (_String&);
extern _Parameter			_lfScalerUpwards,
							_lfScalingFactorThreshold,
							_logLFScaler;	

extern	_GrowingVector		_scalerMultipliers, 
							_scalerDividers;

_Parameter					acquireScalerMultiplier (long);
_Parameter					myLog (_Parameter);

#endif
		
		
 
		
		
		
		
		

	