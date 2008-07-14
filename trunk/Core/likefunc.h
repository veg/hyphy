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
		_LikelihoodFunction (_String&); // from triplets
		_LikelihoodFunction (_LikelihoodFunction&); // stack copy
		
		bool		Construct (_String&,_VariableContainer*);
		
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
		_Matrix* 	ConstructCategoryMatrix (bool complete = false, bool = true);
		_Parameter	SimplexMethod (_Parameter& precision);
		void		Anneal (_Parameter& precision);
			
		void		Simulate (_DataSet &,_List&, _Matrix* = nil, _Matrix* = nil, _Matrix* = nil, _String* = nil);
		void		ReconstructAncestors (_DataSet &, bool = false);
		
		long		MaximumDimension (void);
		
virtual	_PMathObj	CovarianceMatrix (_SimpleList* = nil);
		
		// compute  covariance matrix  based on the Hessian
		// optional list of parameters to estimate the conditional covariance for


		void		RescanAllVariables (void);
													
		long		DependOnTree	   		(_String&);
		long		DependOnModel	   		(_String&);
		long		DependOnDS	   	   		(long);
		bool		DependOnDF	       		(long ID) {return theDataFilters.Find(ID)>=0;}
		bool		MapTreeTipsToData  		(long, bool leafScan = false);
		bool		UpdateFilterSize   		(long);
		void		VoidOldResults 	   		(void) {computationalResults.Clear();}
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
												
	protected:
	
		_Matrix* 		PairwiseDistances 		(long index);
		void			CheckDependentBounds	(void);
		void			AllocateSiteResults 	(void);
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
		
		_Parameter		ComputeBlock 			(long, _Parameter* siteResults = nil);	 
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

		long			CostOfPath	 		  (_DataSetFilter*, _TheTree* , _SimpleList&);

		void			BuildLeafProbs 		  (node<long>& , long*, long&, _DataSet&, _TheTree*, long&, bool, long, _DataSetFilter*, long, _DataSet* = nil);
		bool			SingleBuildLeafProbs  (node<long>&, long, _SimpleList&, _SimpleList&, _TheTree*, bool,_DataSetFilter*, _SimpleList* = nil);
		void			CodonNeutralSimulate  (node<long>&, long, bool,_Matrix*,_Matrix*, _Parameter&, _Parameter&);		
					    
		bool			HasBlockChanged		  (long);
		long			BlockLength			  (long);
		
static	void			RandomizeList (_SimpleList&, long);		
static	void			CheckFibonacci (_Parameter);

		
		
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
		void		 	RecurseCategory				(long,long,long,long,_Parameter);
		void	  		RecurseConstantOnPartition  (long, long, long, long, _Parameter, _Matrix&);
		void		 	FindMaxCategory				(long,long,long,long,long,_Matrix&);
		void		 	WriteAllCategories			(long,long,long,long,long,_Matrix&);
		void			WriteAllCategoriesTemplate  (long, long, long, long, long, _Parameter*, long*);
		void		 	SetNthBit 					(long&,char);
		bool		 	CheckNthBit 				(long&,char);
		void		 	BuildIncrements 			(long, _SimpleList&);
		char		 	HighestBit 					(long);
		char			LowestBit  					(long);
		long			HasHiddenMarkov				(long, bool hmm = true);
		_Matrix*	 	RemapMatrix					(_Matrix*);
		void			CleanUpOptimize				(void);
		void			ComputeBlockForTemplate		(long, bool = false);
		void			ComputeBlockForTemplate2	(long, _Parameter*, _Parameter*, long);
		void			DeleteCaches				(bool = true);
		_SimpleList	 	theTrees, 
						theDataFilters, 
						theProbabilities, 
						indexInd, 
						indexDep, 
						indexCat,
						*nonConstantDep,
					 	blockDependancies;
					 	
		_List			optimalOrders,
						leafSkips,
						computationalResults;
						
		long		 	blockComputed,
						templateKind,
						// 1 if template is per site, otherwise = 0. - variable index if it is hidden markov
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
	
		_Parameter**		conditionalInternalNodeLikelihoodCaches;
		_GrowingVector*		conditionalTerminalNodeLikelihoodCaches;
		long	  **		conditionalTerminalNodeStateFlag;
	
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

#endif
		
		
 
		
		
		
		
		

	