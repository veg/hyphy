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

#ifndef		__CALCNODE__
#define		__CALCNODE__

#include "parser.h"
#include "classes.h"
#include "site.h"


#define UNROOTED 						0
#define ROOTED_LEFT 					1
#define ROOTED_RIGHT 					2


#ifdef	USE_SCALING_TO_FIX_UNDERFLOW
	#undef USE_SCALING_TO_FIX_UNDERFLOW
#endif

#define USE_SCALING_TO_FIX_UNDERFLOW	0
#define	UNDERFLOW_SCALING_MAX			16




//_______________________________________________________________________________________________
					   
class		_CalcNode: public _VariableContainer {

	public:
	
												// constructors
	
	_CalcNode 			(void); 				// default constructor, doesn't do much
	_CalcNode 			(_String, _String, int  = 4, _VariableContainer* = nil, _AVLListXL * = nil);
												// construct a node from a string of the form
												// codeBase specifies the number of distinct states (4 for nucleotides, 61 for codons etc)
												// matrix name, <optional comma separated variable declarations, inititalizations>
												// also should be passed the pointer to a container tree
												
	virtual				~_CalcNode 		(void);

	virtual	long		ObjectClass 	(void) 			
							{ return TREE_NODE; }
	
	virtual void		Duplicate 		(BaseRef); 
	
	virtual	long 		FreeUpMemory 	(long);
	
	void				InitializeCN 	( _String&, int, _VariableContainer*, _AVLListXL * = nil);
	
	virtual BaseRef		makeDynamic 	(void); 	
						// creates a dynamic copy of this object
	
	virtual BaseRef 	toStr 			(void);		
						// converts this object to string
	
	_Parameter&			operator[]		(unsigned long); 	
						// access the i-th element of the
						// probabilities (i = 0..codeBase-1)
												
	void				SetCodeBase      (int);		 	
						// change the codeBase value for this node
						// this will resize the vector used to handle frequencies
		
	void				RecomputeMatrix  (long = 0, long = 1,_Matrix* = nil);		
						// reexponentiate the transition matrix and 
						// store it in compExp.
												
	virtual	bool	    HasChanged 		 (void);
	virtual	bool	    NeedToExponentiate(long = -1);
												
	bool				IsFlagged 		 (void)  
								{return theProbs[0]==-3.1415296;}
	void				SetFlag  		 (void)  
								{theProbs[0]=-3.1415296; }
	
	void				SetSummedFlag    (void)  
								{ if (theProbs[0]>=0) theProbs[0] -= 2.0; }
	bool				IsSummedFlagged (void)  
								{return theProbs[0]<0.0;}
	void				RemoveSummedFlag (void)  
								{ if (theProbs[0]<0) theProbs[0]+= 2.0;}
	
	_Parameter		    GetProbs 		(long k) 
							{return theProbs[k];}
	_Parameter* 	    GetProbs 		(void)   
							{return theProbs;}
	
	void			    SetCompExp 		(_Matrix* mx) 
							{ mx->nInstances++; compExp = mx; }
	void			    SetCompMatrix 	(long);
	_Matrix*    		GetCompExp 		(long catID = -1);
							
	_Formula*		    RecurseMC 		(long , node<long>* , bool first = false, char rooted = UNROOTED);
	
	long			    GetCodeBase 	(void) 
							{ return cBase;  }
	
	_Parameter  	    BranchLength 	(void); 
	virtual	long		SetDependance   (long);
	
	node<long>* 	    LocateMeInTree  (void);
				// return the tree structure node corresponing to this one...
	long			    ConvertToSimpleMatrix (void);
	void			    ConvertFromSimpleMatrix (void);
	_Matrix*		    ComputeModelMatrix(bool expMe=false);
	long			    GetTheModelID 	(void)
							{ return theModel;}
	bool				MatchSubtree	(_CalcNode*);
	virtual	void		RemoveModel		(void);
	virtual long		CheckForReferenceNode
										(void);
										
			void		SetRefNode		(long rn) 
										{referenceNode = rn;slaveNodes = 0;}	
			void		AddRefNode		(void)									
										{referenceNode --;}
							
	friend	class 	    _TheTree;
	
	public:
		_Parameter*	 	theProbs; 	 	// list of transitional probabilities
		long		 	lastState;
		
	protected:
	
		_SimpleList  	categoryVariables, 
						categoryIndexVars;
		

		_Matrix	  *  	compExp;	  	// matrix exponential computed previously
		_Matrix	  ** 	matrixCache; 	// only meaningful for category computations
								 
	    long     		cBase,			// dimension of theProbs
	    				nodeIndex,
	    				referenceNode,
	    				slaveNodes;
};	

//_______________________________________________________________________________________________

#define		HY_BRANCH_SELECT		0x01
#define		HY_BRANCH_DESELECT		0xFFFFFFFE

struct		nodeCoord {

	_Parameter	h,
				v,
				auxD, 
				bL,
				label1,
				label2;
				
	long		varRef,
				auxL,
				textWidth,
				color,
				labelColor,
				flags;
	
	_String		branchName,
				branchTag;
	
}; // used for tree imaging

//_______________________________________________________________________________________________
 
class _TheTree; // forward declaration for xlc

//_______________________________________________________________________________________________

class _TreeTopology: public _CalcNode {
	
	protected:
	
		virtual	void			PreTreeConstructor 					(bool);
		virtual bool			MainTreeConstructor					(_String&,bool = true);
		virtual void			PostTreeConstructor					(bool);
 				node<long>* 	prepTree4Comparison	 				(_List&, _SimpleList&, node<long>* = nil);
 				void			destroyCompTree						(node<long>*);
				_List*			SplitTreeIntoClustersInt 			(node<long>*, _List*, _AVLListX&, long, long);
 				char			internalTreeCompare 				(node<long>*, node<long>*, _SimpleList*, char, long, node<long>*, bool = false);
 				char			internalNodeCompare 				(node<long>*, node<long>*, _SimpleList&, _SimpleList*, bool, long, node<long>*, bool = false);
		virtual	_PMathObj 		FlatRepresentation  				(void);
				void			FindCOTHelper 						(node<long>*, long, _Matrix&, _Matrix&, _Matrix&, _List&, _AVLListX&, _Parameter);
				void 			FindCOTHelper2 						(node<long>*, _Matrix&, _Matrix&, _AVLListX&, node<long>*, _Parameter);

	public:

 				node<long>		*theRoot, 
 								*currentNode;

				_List	    	flatTree, 
								flatCLeaves;

				char			rooted;

		virtual void			toFileStr 							(FILE*);
		virtual	BaseRef 		toStr								(void);
				void 			RerootTreeInternalTraverser 		(long, bool,_String&, long  = -1, bool = false);
	
					   			_TreeTopology 						(void);	
					   			_TreeTopology 						(_String, _String&, bool = true);
					   			_TreeTopology 						(_String*);
					   			_TreeTopology						(_TheTree*);
					   			
		virtual		  			~_TreeTopology 						(void);

		virtual	 bool			Equal 								(_PMathObj);
		virtual	 _PMathObj		Compute 							(void);
		virtual  BaseRef		makeDynamic							(void);
				 node<long>*	CopyTreeStructure   				(node<long>*, bool);	
		virtual	 bool	 		FinalizeNode 						(node<long>*, long, _String&, _String&, _String&);


				 bool			IsCurrentNodeATip					(void);
				 bool			IsCurrentNodeTheRoot				(void);
				 bool			IsDegenerate 						(void);

		virtual _PMathObj 		Execute 							(long, _PMathObj = nil , _PMathObj = nil);
		virtual	_PMathObj 		TipCount  	 						(void);
		virtual	_PMathObj 		BranchCount 						(void);
		virtual	_PMathObj 		AVLRepresentation 					(_PMathObj);
		virtual	long			ObjectClass 						(void) 			
																		{ return TOPOLOGY; } 
		virtual _AssociativeList*	
								FindCOT								(_PMathObj);

				void			DepthWiseT							(bool = false); 
				void	 		DepthWiseTRight 					(bool = false); 
				void			DepthWiseTLevel 					(long& level, bool = false);
				void		 	StepWiseT 							(bool = false); 
				void		 	StepWiseTLevel  					(long&, bool = false); 
				void		 	LeafWiseT 							(bool = false); 
		
		virtual	void			GetNodeName							(node<long> *, _String&, bool = false);
		virtual	void			GetBranchLength						(node<long> *, _String&);
		virtual	void			GetBranchLength						(node<long> *, _Parameter&);
		virtual	void			GetBranchValue						(node<long> *, _String&);
		virtual	void			GetBranchVarValue					(node<long> *, _String&, long);
		virtual void			PasteBranchLength 					(node<long> *, _String&, long, _Parameter factor = 1.);

				node<long>& 	GetRoot 							(void) 				
																		{ return *theRoot;}
				void			SetRoot 							(node<long>* r) 	
																		{ theRoot = r;}
				node<long>& 	GetCurrentNode 						(void) 
																		{return *currentNode;}
			    void			SubTreeString 						(_String&, bool = false, long = -1, _AVLListXL* = nil);
			    
 				_String			CompareTrees						(_TreeTopology*);
				_String			MatchTreePattern					(_TreeTopology*);
		virtual	_PMathObj 		TipName	 							(_PMathObj);
		virtual	_PMathObj 		BranchName	 						(_PMathObj, bool = false);
		virtual	_PMathObj 		BranchLength	 					(_PMathObj);
		virtual	_PMathObj 		RerootTree	 						(_PMathObj);
				_List*			SplitTreeIntoClusters 				(unsigned long, unsigned long);
				void			SetLeafName						    (long, _String*);
};

#if USE_SCALING_TO_FIX_UNDERFLOW
	extern _Parameter scalingLogConstant;
#endif

//_______________________________________________________________________________________________

class _TheTree: public _TreeTopology {

// theModel matrix of _TheTree contains the column matrix of probabilities, which is computed
// based on the DataSetFilter passed on to the tree at the initialization stage
	
 public:

	_TheTree ();												// default constructor - doesn't do much
	_TheTree (_String name, _String& parms, bool = true);		// builds a tree from a string
	_TheTree (_String name, _TreeTopology*);					// builds a tree from a tree topology
	
	
	virtual			   		~_TheTree 					(void);
	virtual	bool	   		HasChanged 					(void);
	virtual	void	   		MarkDone					(void);
			bool	   		HasChanged2 				(void);
	
			_CalcNode* 	 	DepthWiseTraversal 			(bool = false); 
																//performs a post-order traversal
			_CalcNode* 	 	DepthWiseTraversalRight 	(bool = false); 
																//performs a post-order tree traversal going right first
			_CalcNode* 	 	DepthWiseTraversalLevel 	(long&, bool = false); 
																//performs a post-order tree traversal
																//storing current node depth
			_CalcNode* 	 	StepWiseTraversal 			(bool = false); 
																//performs a pre-order  tree traversal
			_CalcNode* 	 	StepWiseTraversalLevel  	(long&, bool = false); 
																//performs a pre-order wise tree traversal
																//storing current node depth
																
			_CalcNode* 	 	LeafWiseTraversal 			(bool = false); 
																//iterate through the leaves (left-to-right)
	
	virtual	 bool	 		FinalizeNode 				(node<long>*, long, _String&, _String&, _String&);
	virtual  BaseRef		makeDynamic					(void);

	virtual  BaseRef		makeDynamicCopy				(_String*);
			 node<long>*	DuplicateTreeStructure 		(node<long>*, _String*, bool);
	virtual	 BaseRef 		toStr						(void);
	virtual	 long			ObjectClass 				(void) 			
															{ return TREE; } 

	virtual  _PMathObj 		Execute 					(long, _PMathObj = nil , _PMathObj = nil);
	virtual	 _PMathObj 		TEXTreeString  				(_PMathObj);
	virtual	 _PMathObj 		PlainTreeString				(_PMathObj,_PMathObj);
	
	virtual	 void			GetNodeName					(node<long> *, _String&, bool = false);
	virtual	 void			GetBranchLength				(node<long> *, _String&);
	virtual	 void			GetBranchLength				(node<long> *, _Parameter&);
	virtual	 void			GetBranchValue				(node<long> *, _String&);
	virtual  _String*		GetBranchSpec				(node<long> *);
	virtual	 void			GetBranchVarValue			(node<long> *, _String&, long);
	
			 void	 		InitializeTreeFrequencies 	(_Matrix *, bool = false);
	
			_Parameter		ReleafTreeAndCheck 			(_DataSetFilter*, long, bool, long categID = -1);				 
			_Parameter		ReleafTreeAndCheckChar  	(_DataSetFilter*, long, bool, long categID = -1);	
			_Parameter		ReleafTreeAndCheckChar4 	(_DataSetFilter*, long, bool, long categID = -1);	

			_Parameter		ReleafTree 					(_DataSetFilter*,long,long,long,long);	
			_Parameter		ReleafTreeDegenerate 		(_DataSetFilter*,long);	

			_Parameter		ReleafTreeCache 			(_DataSetFilter*,long,long,long,long,long);	
#if USE_SCALING_TO_FIX_UNDERFLOW
			_Parameter		ThreadReleafTreeCache		(_DataSetFilter*,long,long,long,long,long,long offset = 0,long fixAttempt = 0, _Parameter = 690.);
			_Parameter		doScaling					(_DataSetFilter*,long,long,long,long,_Parameter, bool, bool);
#else
			_Parameter		ThreadReleafTreeCache		(_DataSetFilter*,long,long,long,long,long,long offset = 0);
#endif

			void			ThreadMatrixUpdate			(long, bool);
			void			SerialMatrixUpdate			(long, bool);
			void			MatrixCacheUpdate			(void);
	
			_Parameter  	ReleafTreeChar 				(_DataSetFilter*,long,long,long,long);	
			_Parameter		ReleafTreeCharCache 		(_DataSetFilter*,long,long,long,long,long);	
			_Parameter		ThreadReleafTreeCharCache 	(_DataSetFilter*,long,long,long,long,long,long offset = 0);	
			_Parameter		ReleafTreeCharDegenerate 	(_DataSetFilter*,long);	
			_Parameter		ReleafTreeChar4 			(_DataSetFilter*,long,long,long,long,long);
			_Parameter		ReleafTreeChar4Degenerate   (_DataSetFilter*,long);	

#if USE_SCALING_TO_FIX_UNDERFLOW
			_Parameter		ThreadReleafTreeChar4		(_DataSetFilter*,long,long,long,long,long,long offset = 0,long fixAttempt = 0, _Parameter = 690.);	
			_Parameter		ReleafTreeChar4 			(_DataSetFilter*,long,long,long,long,long fixAttempt = 0, _Parameter = 690.);
			_Parameter		doChar4Scaling				(_DataSetFilter*,long,long,long,long, _Parameter, bool, bool);
			_Parameter		doChar4Scaling_nc			(_DataSetFilter*,long,long,_Parameter, bool, bool);
#else
			_Parameter		ThreadReleafTreeChar4		(_DataSetFilter*,long,long,long,long,long,long offset = 0);
			_Parameter		ReleafTreeChar4 			(_DataSetFilter*,long,long,long,long);	
#endif
			_Parameter	Probij 							(long, long, _CalcNode*);
			_Parameter	ReleafTreeCharNumFilter4Tree3	(_DataSetFilterNumeric*, long, long = 0);
			_Parameter	PruneTree	  					(long categID = -1);
			_Parameter	PruneTreeChar 					(long categID = -1);
			_Parameter	PruneTreeCharCache 				(long categID = -1);
			_Parameter	PruneTreeChar4					(long categID = -1);
			_Parameter	PruneTreeChar4Cache 			(long categID = -1);
			
			_List*		RecoverAncestralSequences 		(_DataSetFilter*, long, long, _Parameter* = nil);
			void		RecoverNodeSupportStates 		(_DataSetFilter*, long, long, _Matrix&);
			void		RecoverNodeSupportStates2 		(node<long>*,_Parameter*,_Parameter*,long);
			_List*		SampleAncestors 				(_DataSetFilter*, node<long>*);
			
			void		PurgeTree						(void); 		
			
			long	 	ComputeReleafingCost    		(_DataSetFilter*, long, long);
			long	 	ComputeReleafingCostChar 		(_DataSetFilter*, long, long);
			void	 	DumpingOrder 					(_DataSetFilter*, _SimpleList&);
			void	 	SetTreeCodeBase 				(long);
			long		IsLinkedToALF					(long&);
			
			bool		HasCache 						(void)	
															{return topLevelNodes.lLength>0;}
							
			long		GetLeafCount 					(void) 
															{return flatLeaves.lLength;}
							
			long		GetINodeCount 				(void) 
															{return flatNodes.lLength	;}

			void 		ScanForVariables 				(_AVLList& l, _AVLList& l2);
			void 		ScanForDVariables 				(_AVLList& l, _AVLList& l2);
			void 		ScanForGVariables 				(_AVLList&, _AVLList&);
			void 		ScanForCVariables 				(_AVLList&);
			void		MolecularClock 					(_String&, _List&);
			
			void		SetUp 							(void);
			void		SetUpMatrices 					(long);
			void		CleanUpMatrices 				(void);
			void		BuildTopLevelCache 				(void);
			void		KillTopLevelCache 				(void);
			void		SetCompMatrices					(long);				

	virtual	void 		ClearConstraints 				(void);
	
			bool		FindScalingVariables 			(_SimpleList&);
			bool		HaveStringBranchLengths 		(void);
		 	void 		AssignLabelsToBranches 			(node<nodeCoord>*, _String*, bool);
			
			node<nodeCoord>* 
						AlignedTipsMapping 				(bool first = false, bool respectRoot = true);
						
			void		AlignNodes 						(node<nodeCoord>*);
						
			node<nodeCoord>* 
						ScaledBranchMapping				(node<nodeCoord>* , _String*, long, long&);

			void 
						ScaledBranchReMapping			(node<nodeCoord>*, _Parameter);
			char&		RootedFlag						(void) 
															{ return rooted; }

			nodeCoord	TreeTEXRecurse 		  			(node<nodeCoord>*,_String&,_Parameter,_Parameter,long,long);
			void		TreePSRecurse 		  			(node<nodeCoord>*,_String&,_Parameter,_Parameter,long,long,long,long,_AssociativeList* = nil);
			
			bool		AllBranchesHaveModels 			(long);
			void		ScanSubtreeVars	 				(_List&, char, _CalcNode*); 
			void		BuildINodeDependancies			(void);
			void		AllocateResultsCache			(long);
			long		CountTreeCategories				(void);
			void		CompileListOfModels				(_SimpleList&);
			
			void 		MarkMatches						(_DataSetFilter*,long,long);
			long 		GetLowerBoundOnCost				(_DataSetFilter*);
			long 		GetLowerBoundOnCostWithOrder	(_DataSetFilter*,_SimpleList*);
			_SimpleList&GetLeftINodes					(void) 
															{return leftiNodes;}
			bool		MatchLeavesToDF					(_SimpleList&, _DataSetFilter*, bool);
			virtual void
						RemoveModel						(void);
			_String*	TreeUserParams					(void);

		 	
		 	_String		CompareSubTrees					(_TheTree*, node<long>*);
		 	_String		FindMaxCommonSubTree			(_TheTree*, long&, _List*);
			void		WeightedCharacterDifferences	(_Parameter, _Matrix*, _Matrix*, long = -1);
			void		AddNodeNamesToDS				(_DataSet*, bool, bool, bool);
			_Parameter	PSStringWidth					(_String&);

		 #if USE_SCALING_TO_FIX_UNDERFLOW
			void		AllocateUnderflowScalers		(long); 				
			void		DeallocateUnderflowScalers		(void); 				
		 #endif
		 
#ifdef	_SLKP_LFENGINE_REWRITE_
		_Parameter		ComputeTreeBlockByBranch		(_SimpleList&, _SimpleList&, _DataSetFilter*, _Parameter*, long*, _GrowingVector*, long = -1);
#endif			

	// --------------------------
	

	long	  * nodeStates;
	char      * nodeMarkers;

	_Parameter* rootIChildrenCache,
			  * marginalLikelihoodCache;
			  
	_AVLListXL* aCache;

	long		categoryCount;

 	#if USE_SCALING_TO_FIX_UNDERFLOW
 		_Matrix * scalingForUnderflow;
 	#endif
 protected:
 
 				
	bool		IntPopulateLeaves	(_DataSetFilter*, long, long);						
 	 	
	_Parameter	ConditionalBranchLikelihood   
									(node<long>* , node<long>* , _Parameter* , _Parameter* , long, long);
									
	_Parameter	ConditionalNodeLikelihood   
									(node<long>* , node<long>* , _Parameter* , _Parameter* , long ,long);
									
	_List*		MapCBaseToCharacters(_DataSetFilter*, bool = true);
					
	
	virtual		void				PreTreeConstructor 					(bool);
	virtual 	void				PostTreeConstructor					(bool);
			  
				
 	// all of the following members exist to speed-up the pruning algorithm
 	// the are created by calling the function set up
 	_SimpleList flatLeaves, 
 				flatNodes, 
 				leftiNodes, 
 	#ifdef	_SLKP_LFENGINE_REWRITE_
				flatParents,
	#endif
				topLevelNodes,
 				topLevelLeftL,
 				topLevelRightL;
 				
};

extern char 	isDefiningATree;
extern _String	expectedNumberOfSubs,
				stringSuppliedLengths,
				includeModelSpecs;
				
#endif