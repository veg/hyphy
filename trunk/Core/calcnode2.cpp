/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2008 
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

#include "calcnode.h"
#include "likefunc.h"
#include "scfg.h"
#include <math.h>

#ifdef 	  __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

extern	long likeFuncEvalCallCount;

#ifdef	_SLKP_LFENGINE_REWRITE_

_Parameter			_lfScalerUpwards		  = pow(2.,200.),
					_lfScalingFactorThreshold = 1./_lfScalerUpwards,
					_logLFScaler			  = 200. *log(2.);

_GrowingVector		_scalerMultipliers, 
					_scalerDividers;

/*----------------------------------------------------------------------------------------------------------*/

_Parameter	acquireScalerMultiplier (long s)
{
	if (s>0)
	{
		if (s >= _scalerMultipliers.used)
			for (long k = _scalerMultipliers.used; k <= s; k++)
				_scalerMultipliers.Store (exp (-_logLFScaler * k));
		return _scalerMultipliers.theData[s];
	}
	s = -s;
	if (s >= _scalerDividers.used)
		for (long k = _scalerDividers.used; k <= s; k++)
			_scalerDividers.Store (exp (_logLFScaler * k));
	return _scalerDividers.theData[s];
}

/*----------------------------------------------------------------------------------------------------------*/
void		_TheTree::ExponentiateMatrices	(_List& expNodes, long tc, long catID)
{
	_List	  matrixQueue,
			  nodesToDo;
		
	for (long nodeID = 0; nodeID < expNodes.lLength; nodeID++)
	{
		long didIncrease = matrixQueue.lLength;
		_CalcNode* thisNode = (_CalcNode*) expNodes(nodeID);
		thisNode->RecomputeMatrix (catID, categoryCount, nil, &matrixQueue);
		if (matrixQueue.lLength - didIncrease)
			nodesToDo << thisNode;
	}
	
	long matrixID;
	
#ifdef _OPENMP
	long nt = cBase<20?1:(MIN(tc, matrixQueue.lLength / 3 + 1));
#endif
	
#pragma omp parallel for default(shared) private (matrixID) schedule(static) if (nt>1)  num_threads (nt)
	for  (matrixID = 0; matrixID < matrixQueue.lLength; matrixID++)
		((_CalcNode*) nodesToDo(matrixID))->SetCompExp (((_Matrix*)matrixQueue(matrixID))->Exponentiate(), catID);	
}

/*----------------------------------------------------------------------------------------------------------*/

long		_TheTree::DetermineNodesForUpdate	(_SimpleList& updateNodes, _List* expNodes, long catID, long addOne)
{
	nodesToUpdate.Populate (flatLeaves.lLength + flatTree.lLength - 1, 0, 0); 
	_CalcNode		*currentTreeNode;
	long			lastNodeID = -1;

	// look for nodes with model changes and mark the path up to the root as needing an update
	
	
	if (addOne >= 0)
		nodesToUpdate.lData[addOne] = 1;
			
	if (forceRecalculationOnTheseBranches.lLength)
	{
		for (long markedNode = 0; markedNode < forceRecalculationOnTheseBranches.lLength; markedNode++)
			nodesToUpdate.lData[forceRecalculationOnTheseBranches.lData[markedNode]] = 1;
		
		forceRecalculationOnTheseBranches.Clear();
	}
	
	for (long nodeID = 0; nodeID < nodesToUpdate.lLength; nodeID++)
	{
		bool	isLeaf	   = nodeID < flatLeaves.lLength;

		currentTreeNode = isLeaf? (((_CalcNode**) flatCLeaves.lData)[nodeID]):
								  (((_CalcNode**) flatTree.lData)  [nodeID - flatLeaves.lLength]);
		
		if (currentTreeNode->NeedToExponentiate (catID))
		{
			if (expNodes)
			{
				(*expNodes) << currentTreeNode;
				//printf ("EXP>%s\n", currentTreeNode->GetName()->sData);
				lastNodeID = nodeID;
			}
			else
				currentTreeNode->RecomputeMatrix (catID, categoryCount, nil);
				
			nodesToUpdate.lData[nodeID] = 1;
		}
				
		if (nodesToUpdate.lData[nodeID])
			nodesToUpdate.lData[flatParents.lData[nodeID]+flatLeaves.lLength] = 1;
	}
	
	
	// one more pass to pick up all descendants of changed internal nodes
	
	for (long nodeID = 0; nodeID < nodesToUpdate.lLength; nodeID++)
		if (nodesToUpdate.lData[flatLeaves.lLength+flatParents.lData[nodeID]] && nodesToUpdate.lData[nodeID] == 0)
			nodesToUpdate.lData[nodeID] = 1;
	
	// write out all changed nodes
	
	for (long nodeID = 0; nodeID < nodesToUpdate.lLength; nodeID++)
		if (nodesToUpdate.lData[nodeID])
			updateNodes << nodeID;
	
	if (expNodes && expNodes->lLength == 1)
		return lastNodeID;
	
	return -1;
}

/*----------------------------------------------------------------------------------------------------------*/

void		_TheTree::FillInConditionals		(_DataSetFilter*		theFilter, _Parameter*	iNodeCache,  _SimpleList*	tcc)
// this utility function will simply fill in all the conditional probability vectors for internal nodes,
// including those that were skipped due to column sorting optimization
// this is useful to avoid code duplication for other functions (e.g. ancestral sampling) that 
// make use of conditional probability vectors, but would not benefit from subtree caching
{
	if (!tcc)
		return;
	
	long			alphabetDimension	  =			theFilter->GetDimension(),
					siteCount			  =			theFilter->NumberDistinctSites();

	for  (long nodeID = 0; nodeID < flatTree.lLength; nodeID++)
	{
		_Parameter * conditionals		= iNodeCache +(nodeID  * siteCount) * alphabetDimension;
		long		currentTCCIndex		= siteCount * nodeID,
					currentTCCBit		= currentTCCIndex % _HY_BITMASK_WIDTH_;
		
		currentTCCIndex /= _HY_BITMASK_WIDTH_;
		for (long siteID = 0; siteID < siteCount; siteID++, conditionals += alphabetDimension)
		{
			if (siteID  && (tcc->lData[currentTCCIndex] & bitMaskArray.masks[currentTCCBit]) > 0)
			{
				for (long k = 0; k < alphabetDimension; k++)
					conditionals[k] = conditionals[k-alphabetDimension];
			}
			if (++currentTCCBit == _HY_BITMASK_WIDTH_)
			{
				currentTCCBit   = 0;
				currentTCCIndex ++;
			}
			
		}
	}
}

/*----------------------------------------------------------------------------------------------------------*/

_Parameter		_TheTree::ComputeTreeBlockByBranch	(					_SimpleList&		siteOrdering, 
																		_SimpleList&		updateNodes, 
																		_SimpleList*		tcc,
																		_DataSetFilter*		theFilter,
																		_Parameter*			iNodeCache,
																		long	  *			lNodeFlags,
																		_Parameter*			scalingAdjustments,
																		_GrowingVector*		lNodeResolutions,
																		long&				overallScaler,
																		long				siteFrom,
																		long				siteTo,
																		long				catID,
																		_Parameter*			storageVec,
																		long*				siteCorrectionCounts,
																		long				setBranch,
																		long*				setBranchTo
																 )
// the updateNodes flags the nodes (leaves followed by inodes in the same order as flatLeaves and flatNodes)
// that must be recomputed
{
	// process the leaves first 
	
	_SimpleList		taggedInternals					(flatNodes.lLength, 0, 0);
	long			alphabetDimension	  =			theFilter->GetDimension(),
					siteCount			  =			theFilter->NumberDistinctSites(),
					alphabetDimensionmod4  =		alphabetDimension-alphabetDimension%4;
	
	_CalcNode		*currentTreeNode;
	long			localScalerChange	  =			0;
	
	if (siteTo	> siteCount)	siteTo = siteCount;
	
	for  (long nodeID = 0; nodeID < updateNodes.lLength; nodeID++)
	{
		long	nodeCode   = updateNodes.lData [nodeID],
				parentCode = flatParents.lData [nodeCode];
		
		bool	isLeaf	   = nodeCode < flatLeaves.lLength;
		
		if (!isLeaf)
			nodeCode -=  flatLeaves.lLength;
		
		_Parameter * parentConditionals = iNodeCache +			  (siteFrom + parentCode  * siteCount) * alphabetDimension;
		if (taggedInternals.lData[parentCode] == 0)
		// mark the parent for update and clear its conditionals if needed
		{
			taggedInternals.lData[parentCode]	  = 1;
			_Parameter		_hprestrict_ *localScalingFactor	  = scalingAdjustments + parentCode*siteCount;
			
			bool	matchSet   = (parentCode == setBranch);

			if (alphabetDimension == 4)
			{
				long k3		= 0;
				if (matchSet)
					for (long k = siteFrom; k < siteTo; k++, k3+=4)
					{
						parentConditionals [k3]   = 0.;
						parentConditionals [k3+1] = 0.;
						parentConditionals [k3+2] = 0.;
						parentConditionals [k3+3] = 0.;
						parentConditionals [k3+setBranchTo[siteOrdering.lData[k]]] = localScalingFactor[k];
					}
				else					
					for (long k = siteFrom; k < siteTo; k++, k3+=4)
					{
						_Parameter scaler = localScalingFactor[k];
						parentConditionals [k3]   = scaler;
						parentConditionals [k3+1] = scaler;
						parentConditionals [k3+2] = scaler;
						parentConditionals [k3+3] = scaler;
					}
			}
			else
			{
				long k3		= 0;
				if (matchSet)
				{
					for (long k = siteFrom; k < siteTo; k++)
					{
						for (long k2 = 0; k2 < alphabetDimension; k2++)
							parentConditionals [k3+k2] = 0.;
						
						parentConditionals[k3 + setBranchTo[siteOrdering.lData[k]]] = localScalingFactor[k];
						k3				   +=	alphabetDimension;
					}					
				}
				else
				{
					for (long k = siteFrom; k < siteTo; k++)
					{
						_Parameter scaler = localScalingFactor[k];
						for (long k2 = 0; k2 < alphabetDimension; k2++, k3++)
							parentConditionals [k3] = scaler;
					}
				}
			}
		}
		
		currentTreeNode = isLeaf? ((_CalcNode*) flatCLeaves (nodeCode)):
								  ((_CalcNode*) flatTree    (nodeCode));
						
		_Parameter  *		_hprestrict_ transitionMatrix = currentTreeNode->GetCompExp(catID)->theData;
		_Parameter  *	    childVector, 
					*		lastUpdatedSite;
		
		if (!isLeaf)
			childVector = iNodeCache + (siteFrom + nodeCode * siteCount) * alphabetDimension;
		
		long currentTCCIndex		,
			 currentTCCBit			,
			 parentTCCIIndex		,
			 parentTCCIBit			;
		
		if (tcc)
		{
			parentTCCIIndex = siteCount * parentCode + siteFrom;
			parentTCCIBit	= parentTCCIIndex % _HY_BITMASK_WIDTH_;
			parentTCCIIndex = parentTCCIIndex / _HY_BITMASK_WIDTH_;
			if (! isLeaf)
			{
				currentTCCIndex = siteCount * nodeCode + siteFrom;
				currentTCCBit	= currentTCCIndex % _HY_BITMASK_WIDTH_;
				currentTCCIndex /= _HY_BITMASK_WIDTH_;				
			}
		}
		
		//long successiveSkips = 0;
		
		for (long siteID = siteFrom; siteID < siteTo; siteID++, parentConditionals += alphabetDimension)
		{		
			if (tcc)
			{
				if (parentTCCIBit == _HY_BITMASK_WIDTH_)
				{
					parentTCCIBit   = 0;
					parentTCCIIndex ++;
				}
				
				if (siteID > siteFrom && (tcc->lData[parentTCCIIndex] & bitMaskArray.masks[parentTCCIBit]) > 0)
				{
					if (!isLeaf)
					{
						childVector     += alphabetDimension;
						if (++currentTCCBit == _HY_BITMASK_WIDTH_)
						{
							currentTCCBit   = 0;
							currentTCCIndex ++;
						}
					}
					parentTCCIBit++;
					continue;
				}
				parentTCCIBit++;
			}
			
			_Parameter  *tMatrix = transitionMatrix,
						sum		 = 0.0;

			char		didScale = 0;
			
			if (isLeaf)
			{
				long siteState = lNodeFlags[nodeCode*siteCount + siteOrdering.lData[siteID]] ;
				if (siteState >= 0)
				// a single character state; sweep down the appropriate column 
				{
					tMatrix  +=  siteState;
					if (alphabetDimension == 4)
					{
						parentConditionals[0] *= tMatrix[0];
						parentConditionals[1] *= tMatrix[4];
						parentConditionals[2] *= tMatrix[8];
						parentConditionals[3] *= tMatrix[12];
					}
					else
					{
						long k = 0;
						for (; k < alphabetDimensionmod4; k+=4, tMatrix += alphabetDimension+alphabetDimension+alphabetDimension+alphabetDimension)  
						{
								parentConditionals[k]   *= tMatrix[0];
								parentConditionals[k+1] *= tMatrix[alphabetDimension];
								parentConditionals[k+2] *= tMatrix[alphabetDimension+alphabetDimension];
								parentConditionals[k+3] *= tMatrix[alphabetDimension+alphabetDimension+alphabetDimension];
						}
						for (; k < alphabetDimension; k++, tMatrix += alphabetDimension)  
							parentConditionals[k] *= *tMatrix;
					}
					continue;
				}
				else
					childVector = lNodeResolutions->theData + (-siteState-1) * alphabetDimension;
			}
			else
			{
				if (tcc)
				{
					if ((tcc->lData[currentTCCIndex] & bitMaskArray.masks[currentTCCBit]) > 0 && siteID > siteFrom)
					{
						for (long k = 0; k < alphabetDimension; k++) 
							childVector[k] = lastUpdatedSite[k];
					}			
					if (++currentTCCBit == _HY_BITMASK_WIDTH_)
					{
						currentTCCBit   = 0;
						currentTCCIndex ++;
					}
					lastUpdatedSite = childVector;
				}
			}
			
			if (alphabetDimension == 4) // special case for nuc data 
			{
				_Parameter t1 = childVector[0] - childVector[3],
						   t2 = childVector[1] - childVector[3],
						   t3 = childVector[2] - childVector[3],
						   t4 = childVector[3];
				
				
				parentConditionals [0] *= tMatrix[0]  * t1 + tMatrix[1] * t2 + tMatrix[2] * t3 + t4;
				parentConditionals [1] *= tMatrix[4]  * t1 + tMatrix[5] * t2 + tMatrix[6] * t3 + t4;
				parentConditionals [2] *= tMatrix[8]  * t1 + tMatrix[9] * t2 + tMatrix[10] * t3 + t4;
				parentConditionals [3] *= tMatrix[12] * t1 + tMatrix[13] * t2 + tMatrix[14] * t3 + t4;
				
				
				// handle scaling if necessary 
				// the check for sum > 0.0 is necessary for 'degenerate' log-L functions (-infinity)
				
				sum		= parentConditionals [0] + parentConditionals [1] + parentConditionals [2] + parentConditionals [3];
				if (sum < _lfScalingFactorThreshold && sum > 0.0)
				{
					_Parameter tryScale									= scalingAdjustments [parentCode*siteCount + siteID] * _lfScalerUpwards;
					if (tryScale < HUGE_VAL)
					{
						parentConditionals [0]							   *= _lfScalerUpwards;
						parentConditionals [1]							   *= _lfScalerUpwards;
						parentConditionals [2]							   *= _lfScalerUpwards;
						parentConditionals [3]							   *= _lfScalerUpwards;
						localScalerChange								   += theFilter->theFrequencies [siteOrdering.lData[siteID]];
						scalingAdjustments [parentCode*siteCount + siteID]  = tryScale;
						didScale											= 1;
					}
				}
				else
				{
					if (sum > _lfScalerUpwards)
					{
						scalingAdjustments [parentCode*siteCount + siteID] *= _lfScalingFactorThreshold;
						parentConditionals [0]							   *= _lfScalingFactorThreshold;
						parentConditionals [1]							   *= _lfScalingFactorThreshold;
						parentConditionals [2]							   *= _lfScalingFactorThreshold;
						parentConditionals [3]							   *= _lfScalingFactorThreshold;
						localScalerChange								   -= theFilter->theFrequencies [siteOrdering.lData[siteID]];
						didScale											= -1;
					}
				}
				
				childVector += 4;
			}
			else
			{
				_Parameter sum = 0.0;
				
				if (alphabetDimension > alphabetDimensionmod4)
					for (long p = 0; p < alphabetDimension; p++)
					{
						_Parameter		accumulator = 0.0;
						
						for (long c = 0; c < alphabetDimensionmod4; c+=4) // 4 - unroll the loop
							accumulator +=  tMatrix[c]   * childVector[c] + 
											tMatrix[c+1] * childVector[c+1] +
											tMatrix[c+2] * childVector[c+2] +
											tMatrix[c+3] * childVector[c+3];
						
						for (long c = alphabetDimensionmod4; c < alphabetDimension; c++) 
							accumulator +=  tMatrix[c] * childVector[c];
						
						tMatrix				  += alphabetDimension;
						sum += (parentConditionals[p] *= accumulator);
					}
				else
					for (long p = 0; p < alphabetDimension; p++)
					{
						_Parameter		accumulator = 0.0;
						
						for (long c = 0; c < alphabetDimensionmod4; c+=4) // 4 - unroll the loop
							accumulator +=  tMatrix[c]   * childVector[c] + 
							tMatrix[c+1] * childVector[c+1] +
							tMatrix[c+2] * childVector[c+2] +
							tMatrix[c+3] * childVector[c+3];
						
						tMatrix				  += alphabetDimension;
						sum += (parentConditionals[p] *= accumulator);
					}
				
				if (sum < _lfScalingFactorThreshold && sum > 0.0)
				{
					_Parameter tryScale									= scalingAdjustments [parentCode*siteCount + siteID] * _lfScalerUpwards;
					if (tryScale < HUGE_VAL)
					{
						scalingAdjustments [parentCode*siteCount + siteID] = tryScale;
						for (long c = 0; c < alphabetDimension; c++) 
							parentConditionals [c] *= _lfScalerUpwards;
						
						localScalerChange									   += theFilter->theFrequencies [siteOrdering.lData[siteID]];
						didScale											    = 1;
					}
				}
				else
				{
					if (sum > _lfScalerUpwards)
					{
						scalingAdjustments [parentCode*siteCount + siteID] *= _lfScalingFactorThreshold;
						for (long c = 0; c < alphabetDimension; c++) 
							parentConditionals [c] *= _lfScalingFactorThreshold;
						localScalerChange								   -= theFilter->theFrequencies [siteOrdering.lData[siteID]];
						didScale											= -1;
					}
				}
				childVector	   += alphabetDimension;
			}

			if (didScale)
			{
				if (siteCorrectionCounts)
					siteCorrectionCounts [siteOrdering.lData[siteID]] += didScale; 
				
				//printf ("NS: site %d node %d \n", siteOrdering.lData[siteID], parentCode);

				if (tcc)
				{
					// look ahead to see if we need to correct for downstream cached nodes  
					long cparentTCCIIndex	=	parentTCCIIndex,
						  cparentTCCIBit	=	parentTCCIBit;
					
					_Parameter				scM;
					if (didScale < 0)
					{	scM = _lfScalingFactorThreshold; }
					else
					{	scM = _lfScalerUpwards;  }
						
					
					for (long sid = siteID + 1; sid < siteTo; sid++,cparentTCCIBit++)
					{
						if (cparentTCCIBit == _HY_BITMASK_WIDTH_)
						{
							cparentTCCIBit   = 0;
							cparentTCCIIndex ++;
						}
						
						if ((tcc->lData[cparentTCCIIndex] & bitMaskArray.masks[cparentTCCIBit]) > 0)
						{
							if (siteCorrectionCounts)
								siteCorrectionCounts [siteOrdering.lData[sid]] += didScale; 
							scalingAdjustments   [parentCode*siteCount + sid] *= scM;
							localScalerChange								+= didScale * theFilter->theFrequencies [siteOrdering.lData[sid]];
						}
						else
							break;
					}
				}
			}
		}
	}
		
	// assemble the entire likelihood
	
	_Parameter _hprestrict_	* rootConditionals = iNodeCache + alphabetDimension * (siteFrom + (flatTree.lLength-1)  * siteCount);
	_Parameter				  result = 0.0;
	

	for (long siteID = siteFrom; siteID < siteTo; siteID++)
	{
		_Parameter accumulator = 0.;
		if (setBranch == flatTree.lLength-1)
		{
			long				rootState = setBranchTo[siteOrdering.lData[siteID]];
			accumulator      += rootConditionals[rootState] * theProbs[rootState];
			rootConditionals += alphabetDimension;
		}
		else
			for (long p = 0; p < alphabetDimension; p++,rootConditionals++)
				accumulator += *rootConditionals * theProbs[p];
							   		
		if (storageVec)
			storageVec [siteOrdering.lData[siteID]] = accumulator;
		else
		{
			if (accumulator <= 0.0)
			{
				result = -A_LARGE_NUMBER;
				break;
			}
			result += log(accumulator) * theFilter->theFrequencies [siteOrdering.lData[siteID]];
		}
	}
	
	if (!storageVec && localScalerChange) 
	{
#pragma omp atomic
		overallScaler += localScalerChange; 
	}
	
	return result;
}

#ifdef _SLKP_DEBUG_

/*----------------------------------------------------------------------------------------------------------*/
void echoNodeList (_SimpleList& theNodes, _SimpleList& leaves, _SimpleList& iNodes)
{
	for (long n = 0; n < theNodes.lLength; n++)
	{
		node<long>* nd = (node<long>*)(theNodes(n)<leaves.lLength?leaves(theNodes(n)):iNodes(theNodes(n)-leaves.lLength));
		printf ("%d %d %s\n", n, theNodes(n), LocateVar(nd->in_object)->GetName()->sData);
	}
}

#endif

/*----------------------------------------------------------------------------------------------------------*/

void			_TheTree::ComputeBranchCache	(	 
													  _SimpleList&			siteOrdering,
													  long					brID,
													  _Parameter*			cache,
													  _Parameter*			iNodeCache,
													  _DataSetFilter*		theFilter,
													  long		   *		lNodeFlags,
													  _Parameter*			scalingAdjustments,
													  long		*			siteCorrectionCounts,
													  _GrowingVector*		lNodeResolutions,
													  long&					overallScaler,
													  long					siteFrom,
													  long					siteTo,
													  long					catID,
													  _SimpleList*			tcc,
													  _Parameter*			siteRes
												)
{	
	
	_SimpleList	taggedNodes (flatLeaves.lLength + flatNodes.lLength, 0, 0),
				nodesToProcess,
				rootPath;
	
	long		myParent			   = brID		-flatLeaves.lLength,
				alphabetDimension	   =			theFilter->GetDimension(),
				alphabetDimensionmod4  =			alphabetDimension - alphabetDimension % 4,
				siteCount			   =			theFilter->NumberDistinctSites();
	
	if (siteTo	> siteCount)	siteTo = siteCount;
	
	do
	{
		taggedNodes.lData[myParent+flatLeaves.lLength] = 1;
		myParent = flatParents.lData[myParent+flatLeaves.lLength];
	}
	while (myParent >= 0);
	

	for (long k = 0; k <  flatLeaves.lLength+flatNodes.lLength; k++)
	{
		myParent = flatParents.lData[k];
		if (taggedNodes.lData[myParent+flatLeaves.lLength] == 1 && taggedNodes.lData[k] == 0)
		{
			if (myParent != brID - flatLeaves.lLength)
				nodesToProcess << k;
		}
		if (taggedNodes.lData[k])
			rootPath << k;
	}
	
	//printf ("ComputeBranchCache at branch %d; siteOdering %s\n", 
	//		brID, _String((_String*)siteOrdering.toStr()).sData);
	
	//echoNodeList (rootPath,flatLeaves,flatNodes );
	//echoNodeList (nodesToProcess,flatLeaves,flatNodes);

	_Parameter * state = cache + alphabetDimension * siteFrom,
			   * childVector;
	
	long		localScalerChange = 0;

	// first populate the downward looking vector of conditionals

	if (brID < flatLeaves.lLength) // a leaf
	{
		for (long siteID = siteFrom; siteID < siteTo; siteID ++, state += alphabetDimension)
		{
			long siteState = lNodeFlags[brID*siteCount + siteOrdering.lData[siteID]] ;
			if (siteState >= 0)
				// a single character state; sweep down the appropriate column 
			{
				for (long s = 0; s < alphabetDimension; s++)
					state[s] = 0.;
				state[siteState] = 1.;
			}
			else
			{
				childVector = lNodeResolutions->theData + (-siteState-1) * alphabetDimension;
				for (long s = 0; s < alphabetDimension; s++)
					state[s] = childVector[s];
			}
		}
	}
	else // an internal branch
	{
		long		nodeCode = brID - flatLeaves.lLength;
		_Parameter *lastUpdated = iNodeCache + (nodeCode * siteCount + siteFrom) * alphabetDimension;
		
		long currentTCCIndex		,
			 currentTCCBit			;	
		
		if (tcc)
		{
				currentTCCIndex = siteCount * nodeCode + siteFrom;
				currentTCCBit	= currentTCCIndex % _HY_BITMASK_WIDTH_;
				currentTCCIndex /= _HY_BITMASK_WIDTH_;				
		}
		
		for (long siteID = siteFrom; siteID < siteTo; siteID ++, state += alphabetDimension)
		{
			if (tcc)
			{
				if ((tcc->lData[currentTCCIndex] & bitMaskArray.masks[currentTCCBit]) == 0)
					lastUpdated = iNodeCache + (nodeCode * siteCount + siteID) * alphabetDimension;
			}
			
			for (long s = 0; s < alphabetDimension; s++)
				state[s] = lastUpdated[s];
			
			if (tcc)
			{
				if (++currentTCCBit == _HY_BITMASK_WIDTH_)
				{
					currentTCCBit   = 0;
					currentTCCIndex ++;
				}
			}
			else
				lastUpdated += alphabetDimension;
		}
	}
	
	state = cache + alphabetDimension * siteCount; 
	
	taggedNodes.Populate (flatTree.lLength, 0, 0);
	rootPath.Flip ();
	
	for  (long nodeID = 0; nodeID < nodesToProcess.lLength + rootPath.lLength - 2; nodeID++)
	{
		bool	notPassedRoot = nodeID<nodesToProcess.lLength;

		long	nodeCode   = notPassedRoot?nodesToProcess.lData [nodeID]:rootPath.lData[nodeID-nodesToProcess.lLength],
				parentCode = notPassedRoot?flatParents.lData [nodeCode]:(rootPath.lData[nodeID-nodesToProcess.lLength+1] - flatLeaves.lLength);
		
		bool	isLeaf	   = nodeCode < flatLeaves.lLength;
		
		if (!isLeaf)
			nodeCode -=  flatLeaves.lLength;
		
		_Parameter * parentConditionals = iNodeCache +			  (siteFrom + parentCode  * siteCount) * alphabetDimension;
		if (taggedNodes.lData[parentCode] == 0)
			// mark the parent for update and clear its conditionals if needed
		{
			taggedNodes.lData[parentCode]	  = 1;
			_Parameter		_hprestrict_ *localScalingFactor	  = scalingAdjustments + parentCode*siteCount;
			if (alphabetDimension == 4)
			{
				long k3		= 0;
				for (long k = siteFrom; k < siteTo; k++, k3+=4)
				{
					_Parameter scaler = localScalingFactor[k];
					parentConditionals [k3]   = scaler;
					parentConditionals [k3+1] = scaler;
					parentConditionals [k3+2] = scaler;
					parentConditionals [k3+3] = scaler;
				}
			}
			else
			{
				long k3		= 0;
				for (long k = siteFrom; k < siteTo; k++)
				{
					_Parameter scaler = localScalingFactor[k];
					for (long k2 = 0; k2 < alphabetDimension; k2++, k3++)
						parentConditionals [k3] = scaler;
				}
			}
		}
		
		_CalcNode    * currentTreeNode = isLeaf? ((_CalcNode*) flatCLeaves (nodeCode)):
											     ((_CalcNode*) flatTree    (notPassedRoot?nodeCode:parentCode));
		
		_Parameter  *		_hprestrict_ transitionMatrix = currentTreeNode->GetCompExp(catID)->theData;
		
		_Parameter  *	    childVector,
					*		lastUpdatedSite;
		
		if (!isLeaf)
			lastUpdatedSite = childVector = iNodeCache + (siteFrom + nodeCode * siteCount) * alphabetDimension;
		
		long currentTCCIndex		,
			 currentTCCBit			,
			 parentTCCIIndex		,
			 parentTCCIBit			;
		
		if (tcc)
		{
			parentTCCIIndex = siteCount * parentCode + siteFrom;
			parentTCCIBit	= parentTCCIIndex % _HY_BITMASK_WIDTH_;
			parentTCCIIndex = parentTCCIIndex / _HY_BITMASK_WIDTH_;
			if (! isLeaf)
			{
				currentTCCIndex = siteCount * nodeCode + siteFrom;
				currentTCCBit	= currentTCCIndex % _HY_BITMASK_WIDTH_;
				currentTCCIndex /= _HY_BITMASK_WIDTH_;				
			}
		}
		
		for (long siteID = siteFrom; siteID < siteTo; siteID++, parentConditionals += alphabetDimension)
		{		 
			_Parameter  *tMatrix = transitionMatrix;			
			
			char canScale = !notPassedRoot;
			
			if (isLeaf)
			{
				long siteState = lNodeFlags[nodeCode*siteCount + siteOrdering.lData[siteID]] ;
				if (siteState >= 0)
					// a single character state; sweep down the appropriate column 
				{
					if (alphabetDimension == 4) // special case for nuc data
					{
						parentConditionals[0] *= tMatrix[siteState];
						parentConditionals[1] *= tMatrix[siteState+4];
						parentConditionals[2] *= tMatrix[siteState+8];
						parentConditionals[3] *= tMatrix[siteState+12];						
					}
					else
					{
						tMatrix  +=  siteState;
						for (long k = 0; k < alphabetDimension; k++, tMatrix += alphabetDimension)  
							parentConditionals[k] *= *tMatrix;
						
					}
					continue;
				}
				else
					childVector = lNodeResolutions->theData + (-siteState-1) * alphabetDimension;
				canScale = 0;
			}
			else
			{
				if (tcc&&notPassedRoot)
				{			
					if ((tcc->lData[currentTCCIndex] & bitMaskArray.masks[currentTCCBit]) > 0 && siteID > siteFrom)
						// the value of this conditional vector needs to be copied from a previously stored site
						// subtree duplication
							for (long k = 0; k < alphabetDimension; k++) 
								childVector[k] = lastUpdatedSite[k];
						else
							lastUpdatedSite = childVector;
					
					if (++currentTCCBit == _HY_BITMASK_WIDTH_)
					{
						currentTCCBit   = 0;
						currentTCCIndex ++;
					}
					if (++parentTCCIBit == _HY_BITMASK_WIDTH_)
					{
						parentTCCIBit   = 0;
						parentTCCIIndex ++;
					}
				}	
			}

			_Parameter sum		= .0;
			char	   didScale = 0;
			
			if (alphabetDimension == 4) // special case for nuc data 
			{
				_Parameter t1 = childVector[0] - childVector[3],
						   t2 = childVector[1] - childVector[3],
						   t3 = childVector[2] - childVector[3],
						   t4 = childVector[3];
				
				
				parentConditionals [0] *= tMatrix[0]  * t1 + tMatrix[1] * t2 + tMatrix[2] * t3 + t4;
				parentConditionals [1] *= tMatrix[4]  * t1 + tMatrix[5] * t2 + tMatrix[6] * t3 + t4;
				parentConditionals [2] *= tMatrix[8]  * t1 + tMatrix[9] * t2 + tMatrix[10] * t3 + t4;
				parentConditionals [3] *= tMatrix[12] * t1 + tMatrix[13] * t2 + tMatrix[14] * t3 + t4;
				
				if (canScale)
				{
					sum		= parentConditionals [0] + parentConditionals [1] + parentConditionals [2] + parentConditionals [3];
					if (sum < _lfScalingFactorThreshold && sum > 0.0)
					{
						_Parameter tryScale									= scalingAdjustments [nodeCode*siteCount + siteID] * _lfScalerUpwards;
						if (tryScale < HUGE_VAL)
						{
								parentConditionals[0]							  *= _lfScalerUpwards;
								parentConditionals[1]							  *= _lfScalerUpwards;
								parentConditionals[2]							  *= _lfScalerUpwards;
								parentConditionals[3]							  *= _lfScalerUpwards;
							
								localScalerChange								   += theFilter->theFrequencies [siteOrdering.lData[siteID]];
								didScale											= 1;
						}
					}
					else
					{
						if (sum > _lfScalerUpwards)
						{
							parentConditionals [0]							   *= _lfScalingFactorThreshold;
							parentConditionals [1]							   *= _lfScalingFactorThreshold;
							parentConditionals [2]							   *= _lfScalingFactorThreshold;
							parentConditionals [3]							   *= _lfScalingFactorThreshold;
							
							localScalerChange								   -= theFilter->theFrequencies [siteOrdering.lData[siteID]];
							didScale											= -1;
						}
					}
				}
				childVector += 4;
			}
			else
			{
				for (long p = 0; p < alphabetDimension; p++)
				{
					_Parameter		accumulator = 0.0;
					
					for (long c = 0; c < alphabetDimensionmod4; c+=4) // 4 - unroll the loop
						accumulator +=  tMatrix[c]   * childVector[c] + 
						tMatrix[c+1] * childVector[c+1] +
						tMatrix[c+2] * childVector[c+2] +
						tMatrix[c+3] * childVector[c+3];
					
					for (long c = alphabetDimensionmod4; c < alphabetDimension; c++) 
						accumulator +=  tMatrix[c] * childVector[c];
					
					tMatrix				  += alphabetDimension;
					sum += (parentConditionals[p] *= accumulator);
				}
				if (canScale)
				{
					if (sum < _lfScalingFactorThreshold && sum > 0.0)
					{
						_Parameter tryScale									= scalingAdjustments [nodeCode*siteCount + siteID] * _lfScalerUpwards;
						if (tryScale < HUGE_VAL)
						{
							for (long c = 0; c < alphabetDimension; c++) 
								parentConditionals [c] *= _lfScalerUpwards;
							
							localScalerChange									   += theFilter->theFrequencies [siteOrdering.lData[siteID]];
							didScale											    = 1;
						}
					}
					else
					{
						if (sum > _lfScalerUpwards)
						{
							for (long c = 0; c < alphabetDimension; c++) 
								parentConditionals [c] *= _lfScalingFactorThreshold;
							
							localScalerChange								   -= theFilter->theFrequencies [siteOrdering.lData[siteID]];
							didScale											= -1;
						}
					}			
				}
				childVector	   += alphabetDimension;
			}
			
			if (didScale&&siteCorrectionCounts)
			{
				siteCorrectionCounts [siteOrdering.lData[siteID]] += didScale;		
				siteRes[siteOrdering.lData[siteID]] *= didScale<0?_lfScalingFactorThreshold:_lfScalerUpwards;
			}
		}
	}
	
	
	_Parameter _hprestrict_	*rootConditionals   = iNodeCache +	(rootPath.lData[rootPath.lLength-2] - flatLeaves.lLength)  * siteCount * alphabetDimension;
	
	for (long ii = siteFrom * alphabetDimension; ii < alphabetDimension*siteTo; ii++)
		state[ii] = rootConditionals[ii];

	if (!siteCorrectionCounts && localScalerChange) 
	{
#pragma omp atomic
		overallScaler += localScalerChange; 
	}
}

/*----------------------------------------------------------------------------------------------------------*/

_Parameter			_TheTree::ComputeLLWithBranchCache (	 
												     _SimpleList&			siteOrdering,
												     long					brID,
												     _Parameter*			cache,
													 _DataSetFilter*		theFilter,
													 long					siteFrom,
													 long					siteTo,
													 long					catID,
													 _Parameter*			storageVec
													 )
{
	long		alphabetDimension	   = theFilter->GetDimension(),
				//alphabetDimensionmod4  = alphabetDimension - alphabetDimension % 4,
				siteCount			   =  theFilter->NumberDistinctSites();
					
	if (siteTo	> siteCount)	siteTo = siteCount;

	_Parameter _hprestrict_	*branchConditionals	= cache				 + siteFrom * alphabetDimension;
	_Parameter _hprestrict_	*rootConditionals   = branchConditionals + siteCount * alphabetDimension;
	_Parameter	result = 0.0;
							
	
	//printf ("ComputeLLWithBranchCache @ %d catID = %d branchID = %d\n", likeFuncEvalCallCount, catID, brID);
	
	_CalcNode *givenTreeNode = brID < flatLeaves.lLength ? (((_CalcNode**) flatCLeaves.lData)[brID]):
														   (((_CalcNode**) flatTree.lData)   [brID - flatLeaves.lLength]);
	
	_Parameter  *	_hprestrict_ transitionMatrix = givenTreeNode->GetCompExp(catID)->theData;
	
	for (long siteID = siteFrom; siteID < siteTo; siteID++)
	{
		_Parameter accumulator = 0.,
					*rmx	       = transitionMatrix;
				

		if (alphabetDimension == 4)
		{
			  accumulator =    rootConditionals[0] * theProbs[0] * 
								 (branchConditionals[0] *  rmx[0] + branchConditionals[1] *  rmx[1] + branchConditionals[2] *  rmx[2] + branchConditionals[3] *  rmx[3]) + 
								 rootConditionals[1] * theProbs[1] * 
								 (branchConditionals[0] *  rmx[4] + branchConditionals[1] *  rmx[5] + branchConditionals[2] *  rmx[6] + branchConditionals[3] *  rmx[7]) + 	
								 rootConditionals[2] * theProbs[2] * 
								 (branchConditionals[0] *  rmx[8] + branchConditionals[1] *  rmx[9] + branchConditionals[2] *  rmx[10] + branchConditionals[3] *  rmx[11]) + 	
								 rootConditionals[3] * theProbs[3] * 
								 (branchConditionals[0] *  rmx[12] + branchConditionals[1] *  rmx[13] + branchConditionals[2] *  rmx[14] + branchConditionals[3] *  rmx[15]);
			 rootConditionals += 4;
		}
		else
		{
			for (long p = 0; p < alphabetDimension; p++,rootConditionals++)
			{
				_Parameter	   r2 = 0.;
				for (long c = 0; c < alphabetDimension; c++, rmx++)
					r2 += branchConditionals[c] *  *rmx;
				accumulator += *rootConditionals * theProbs[p] * r2;
			}
			
		}
		
		branchConditionals += alphabetDimension;
		
		if (storageVec)
			storageVec [siteOrdering.lData[siteID]] = accumulator;
		else
		{
			if (accumulator <= 0.0)
			{
				result = -A_LARGE_NUMBER;
				break;
			}
			result += log(accumulator) * theFilter->theFrequencies [siteOrdering.lData[siteID]];
		}
	}
	return result;
}

/*----------------------------------------------------------------------------------------------------------*/

_Parameter		_TheTree::ComputeTwoSequenceLikelihood	
													(			
													 _SimpleList   &	siteOrdering,
													 _DataSetFilter*	theFilter,
													 long	  *			lNodeFlags,
													 _GrowingVector*	lNodeResolutions,
													 long				siteFrom,
													 long				siteTo,
													 long				catID,
													 _Parameter*		storageVec
													 )
// the updateNodes flags the nodes (leaves followed by inodes in the same order as flatLeaves and flatNodes)
// that must be recomputed
{
	// process the leaves first 
	
	long			alphabetDimension	   =			theFilter->GetDimension(),
					siteCount			   =			theFilter->NumberDistinctSites(),
					alphabetDimensionmod4  =			alphabetDimension-alphabetDimension%4;
	
	_CalcNode		*theNode			   =			((_CalcNode*) flatCLeaves (0));
	_Parameter  *	_hprestrict_ transitionMatrix 
										   =			theNode->GetCompExp(catID)->theData,
					result				   =			0.;

	if (siteTo	> siteCount)	siteTo = siteCount;
	
	for (long siteID = siteFrom; siteID < siteTo; siteID++)
	{		
		_Parameter  *tMatrix = transitionMatrix,
					sum	 = 0.;
		
		long siteState1 = lNodeFlags[siteOrdering.lData[siteID]],
			 siteState2 = lNodeFlags[siteCount + siteOrdering.lData[siteID]];
		
		if (siteState1 >= 0)
			// a single character state; sweep down the appropriate column 
		{
			if (siteState2 >= 0) // both completely resolved;
			{
				sum = tMatrix[siteState1*alphabetDimension + siteState2];
			}
			else // first resolved, second is not
			{
				_Parameter* childVector = lNodeResolutions->theData + (-siteState2-1) * alphabetDimension;
				tMatrix	  +=  siteState1*alphabetDimension;
				if (alphabetDimension == 4) // special case for nuc data
				{
					sum = tMatrix[0] * childVector[0] +
						  tMatrix[1] * childVector[1] +
						  tMatrix[2] * childVector[2] +
						  tMatrix[3] * childVector[3];
				}
				else
				{
					for (long c = 0; c < alphabetDimensionmod4; c+=4) // 4 - unroll the loop
						sum +=  tMatrix[c]   * childVector[c] + 
								tMatrix[c+1] * childVector[c+1] +
								tMatrix[c+2] * childVector[c+2] +
								tMatrix[c+3] * childVector[c+3];
					
					for (long c = alphabetDimensionmod4; c < alphabetDimension; c++) 
						sum +=  tMatrix[c] * childVector[c];
				}
			}
			sum *= theProbs[siteState1];
		}
		else
		{
			if (siteState2 >=0 ) // second resolved, but not the first
			{
				_Parameter* childVector = lNodeResolutions->theData + (-siteState1-1) * alphabetDimension;
				tMatrix	               +=  siteState2;
				if (alphabetDimension == 4) // special case for nuc data
				{
					sum = tMatrix[0] * childVector[0]  * theProbs[0]+
						  tMatrix[4] * childVector[1]  * theProbs[1]+
						  tMatrix[8] * childVector[2]  * theProbs[2]+
						  tMatrix[12] * childVector[3] * theProbs[3];
					
				}
				else
				{
					for (long c = 0; c < alphabetDimensionmod4; c+=4, tMatrix += 4*alphabetDimension) // 4 - unroll the loop
						sum +=  tMatrix[0]   * childVector[c] * theProbs[c]+ 
							    tMatrix[alphabetDimension] * childVector[c+1] * theProbs[c+1]+
								tMatrix[alphabetDimension+alphabetDimension] * childVector[c+2] * theProbs[c+2]+
							    tMatrix[alphabetDimension+alphabetDimension+alphabetDimension] * childVector[c+3] * theProbs[c+3];
					
					for (long c = alphabetDimensionmod4; c < alphabetDimension; c++, tMatrix += alphabetDimension) 
						sum +=  tMatrix[0] * childVector[c] * theProbs[c];
				}
			}
			else
			// both unresolved
			{
				_Parameter *childVector1 = lNodeResolutions->theData + (-siteState1-1) * alphabetDimension,
						   *childVector2 = lNodeResolutions->theData + (-siteState2-1) * alphabetDimension;

				if (alphabetDimension == 4) // special case for nuc data
				{
					sum = (tMatrix[0] * childVector2[0] + tMatrix[1] * childVector2[1] + tMatrix[2] * childVector2[2] + tMatrix[3] * childVector2[3])     * childVector1[0] * theProbs[0]+ 
						  (tMatrix[4] * childVector2[0] + tMatrix[5] * childVector2[1] + tMatrix[6] * childVector2[2] + tMatrix[7] * childVector2[3])     * childVector1[1] * theProbs[1]+
					      (tMatrix[8] * childVector2[0] + tMatrix[9] * childVector2[1] + tMatrix[10] * childVector2[2] + tMatrix[11] * childVector2[3])   * childVector1[2] * theProbs[2] +
					      (tMatrix[12] * childVector2[0] + tMatrix[13] * childVector2[1] + tMatrix[14] * childVector2[2] + tMatrix[15] * childVector2[3]) * childVector1[3] * theProbs[3];
					
				}
				else
				{
					for (long r = 0; r < alphabetDimension; r++) // 4 - unroll the loop
					{
						if (childVector1[r] > 0.0)
						{
							_Parameter sum2 = 0.0;
							for (long c = 0; c < alphabetDimensionmod4; c+=4, tMatrix += 4) // 4 - unroll the loop
								sum2   +=  tMatrix[0] * childVector2[c] + 
										   tMatrix[1] * childVector2[c+1]+
										   tMatrix[2] * childVector2[c+2]+
										   tMatrix[3] * childVector2[c+3];
					
							for (long c = alphabetDimensionmod4; c < alphabetDimension; c++, tMatrix ++) 
								sum2 +=  tMatrix[0] * childVector2[c];
							
							sum += sum2 * childVector1[r] * theProbs[r];
						}
						else
							tMatrix += alphabetDimension;
					}
				}
			}
		}
		if (storageVec)
			storageVec [siteOrdering.lData[siteID]] = sum;
		else
		{
			if (sum <= 0.0)
				return -A_LARGE_NUMBER;
			else
			{
				//printf ("%d: %g\n", siteID, sum);
				result += log(sum) * theFilter->theFrequencies [siteOrdering.lData[siteID]];
			}
		}
	}
		
	return result;
}

//_______________________________________________________________________________________________

void	 _TheTree::SampleAncestorsBySequence (_DataSetFilter* dsf, _SimpleList& siteOrdering, node<long>* currentNode, _AVLListX* nodeToIndex, _Parameter* iNodeCache, 
																   _List& result, _SimpleList* parentStates, _List& expandedSiteMap, _Parameter* catAssignments, long catCount)	

// must be called initially with the root node 


// dsf:							the filter to sample from
// siteOrdering:				the map from cache ordering to actual pattern ordering
// currentNode:					the node index to sample for
// nodeToIndex:					an AVL that maps the address of an internal node pointed to by node<long> to its order in the tree postorder traversal
// iNodeCache:					internal node likelihood caches
// results:						the list that will store sampled strings
// parentStates:				sampled states for the parent of the current node
// expandedSiteMap:				a list of simple lists giving site indices for each unique column pattern in the alignment
// catAssignments:				a vector assigning a (partition specific) rate category to each site (nil if no rate variation)
// catCount:					the number of rate classes

// this needs to be updated to deal with traversal caches!
{	
	long					  childrenCount		= currentNode->get_num_nodes();
	
	if (childrenCount)
	{
		long			siteCount						= dsf->NumberDistinctSites	(),
						alphabetDimension				= dsf->GetDimension			(),
						nodeIndex						= nodeToIndex->GetXtra (nodeToIndex->Find ((BaseRef)currentNode)),
						unitLength						= dsf->GetUnitLength(),
						catBlockShifter					= catAssignments?(dsf->NumberDistinctSites()*GetINodeCount()):0;

		
		_CalcNode *		currentTreeNode	= ((_CalcNode*)	flatTree (nodeIndex));
		_SimpleList		sampledStates	  (dsf->GetSiteCount (), 0, 0);
		
		_Parameter  *		_hprestrict_ transitionMatrix = (catAssignments|| !parentStates)?nil:currentTreeNode->GetCompExp()->theData;
		_Parameter  *		_hprestrict_ conditionals	  = catAssignments?nil:(iNodeCache + nodeIndex  * siteCount * alphabetDimension);
		_Parameter	*		_hprestrict_ cache			  = new _Parameter [alphabetDimension];
		
		for (long			pattern = 0; pattern < siteCount; pattern++)
		{
			_SimpleList*	patternMap = (_SimpleList*) expandedSiteMap (siteOrdering.lData[pattern]);
			if (catAssignments)
			{
				long localCatID = catAssignments[siteOrdering.lData[pattern]];
				if (parentStates)
					transitionMatrix = currentTreeNode->GetCompExp(localCatID)->theData;
					
				conditionals	 = iNodeCache + localCatID*alphabetDimension*catBlockShifter + (pattern + nodeIndex  * siteCount) * alphabetDimension;
			}
			
			for (long site = 0; site < patternMap->lLength; site++)
			{
				long		siteID =   patternMap->lData[site];
				
				_Parameter	randVal  = genrand_real2(),
							totalSum = 0.;
				
				_Parameter  *		_hprestrict_  matrixRow; 
				
				if  (parentStates == nil)
					matrixRow = theProbs;
				else
					matrixRow = transitionMatrix + parentStates->lData[siteID] * alphabetDimension;
				
				for (long i = 0; i<alphabetDimension; i++)
					totalSum += (cache[i] = matrixRow[i]*conditionals[i]);
				
				randVal *= totalSum;
				totalSum	= 0.0;
				long	    sampledChar = -1;
				while		(totalSum < randVal)
				{
					sampledChar ++;
					totalSum += cache[sampledChar];
				}
				
				sampledStates.lData[siteID] = sampledChar;
			}
			
			if (catAssignments == nil)
				conditionals += alphabetDimension;
		}
		
		delete (cache);
		
		_SimpleList  conversion;
		_AVLListXL	 conversionAVL (&conversion);

		_String * sampledSequence = new _String (siteCount*unitLength, true);
		_String  letterValue (unitLength,false);
		for (long charIndexer = 0; charIndexer < sampledStates.lLength; charIndexer++)
		{
			dsf->ConvertCodeToLettersBuffered (dsf->CorrectCode(sampledStates.lData[charIndexer]), unitLength, letterValue.sData, &conversionAVL);
			(*sampledSequence) << letterValue;
		}
		sampledSequence->Finalize();
		result.AppendNewInstance(sampledSequence);
		//printf ("%d: %s\n", nodeIndex, sampledSequence->sData);
		
		for (long child = 1; child <= childrenCount; child ++)
			SampleAncestorsBySequence (dsf, siteOrdering, currentNode->go_down(child), nodeToIndex, iNodeCache, result, &sampledStates, expandedSiteMap, catAssignments, catCount);
	}
}


//_______________________________________________________________________________________________

_List*	 _TheTree::RecoverAncestralSequences (_DataSetFilter* dsf, 
											  _SimpleList& siteOrdering,  
											  _List& expandedSiteMap, 
											  _Parameter* iNodeCache, 
											  _Parameter* catAssignments, 
											  long catCount,
											  long* lNodeFlags,
											  _GrowingVector* lNodeResolutions
											  )	


// dsf:							the filter to sample from
// siteOrdering:				the map from cache ordering to actual pattern ordering
// nodeToIndex:					an AVL that maps the address of an internal node pointed to by node<long> to its order in the tree postorder traversal
// expandedSiteMap:				a list of simple lists giving site indices for each unique column pattern in the alignment
// iNodeCache:					internal node likelihood caches
// catAssignments:				a vector assigning a (partition specific) rate category to each site
// catCount:					the number of rate classes

{	
	
	long			patternCount					= dsf->NumberDistinctSites	(),
					alphabetDimension				= dsf->GetDimension			(),
					unitLength						= dsf->GetUnitLength		(),
					iNodeCount						= GetINodeCount				(),
					siteCount						= dsf->GetSiteCount			(),
					allNodeCount					= 0,
					*stateCache						= new long [patternCount*(iNodeCount-1)*alphabetDimension],
					*leafBuffer						= new long [alphabetDimension];	
	
					// a Patterns x Int-Nodes x CharStates integer table
					// with the best character assignment for node i given that its parent state is j for a given site
	
	_Parameter			*buffer							= new _Parameter [alphabetDimension];
					// iNodeCache will be OVERWRITTEN with conditional pair (i,j) conditional likelihoods
	
	checkPointer	(stateCache);
	checkPointer	(leafBuffer);
	
	_SimpleList		taggedInternals (iNodeCount, 0, 0),
					postToIn;
	
	MapPostOrderToInOderTraversal (postToIn);
	// all nodes except the root
	
	allNodeCount = iNodeCount + GetLeafCount () - 1;
	
	for  (long nodeID = 0; nodeID < allNodeCount; nodeID++)
	{
		long	parentCode = flatParents.lData [nodeID],
				nodeCode   = nodeID;
		
		bool	isLeaf	   = nodeID < flatLeaves.lLength;
		

		if (!isLeaf)
		{
			nodeCode -=  flatLeaves.lLength;
			AddBranchToForcedRecomputeList (nodeCode);
		}
		
		_Parameter * parentConditionals = iNodeCache + parentCode * alphabetDimension * patternCount;
		
		if (taggedInternals.lData[parentCode] == 0)
			// mark the parent for update and clear its conditionals if needed
		{
			taggedInternals.lData[parentCode]	  = 1;
			long k3		= 0;
			for (long k = 0; k < patternCount; k++)
			{
				for (long k2 = 0; k2 < alphabetDimension; k2++, k3++)
					parentConditionals [k3] = 1.;
			}
		}
		
		_CalcNode *			 currentTreeNode = isLeaf? ((_CalcNode*) flatCLeaves (nodeCode)):((_CalcNode*) flatTree    (nodeCode));
		_Parameter  *		_hprestrict_ transitionMatrix = catAssignments?nil:currentTreeNode->GetCompExp()->theData;
		// this will need to be toggled on a per site basis
		_Parameter  *	    childVector;
		
		if (!isLeaf)
			childVector = iNodeCache + (nodeCode * patternCount) * alphabetDimension;
				
		for (long siteID = 0; siteID < patternCount; siteID++, parentConditionals += alphabetDimension)
		{					
			if (catAssignments)
				transitionMatrix = currentTreeNode->GetCompExp(catAssignments[siteOrdering.lData[siteID]])->theData;

			_Parameter  _hprestrict_ *tMatrix = transitionMatrix;
			if (isLeaf)
			{
				long siteState = lNodeFlags[nodeCode*patternCount + siteOrdering.lData[siteID]] ;
				if (siteState >= 0)
					// a fully resolved leaf
				{
					tMatrix  +=  siteState;
					for (long k = 0; k < alphabetDimension; k++, tMatrix += alphabetDimension)  
						parentConditionals[k] *= *tMatrix;
					continue;
				}
				else
					// an ambiguous leaf
					childVector = lNodeResolutions->theData + (-siteState-1) * alphabetDimension;
					
			}
			
			// now repopulate this vector as necessary -- if we are here this means 
			// that the subtree below has been completely processed, 
			// the i-th cell of childVector contains the likelihood of the _optimal_
			// assignment in the subtree below given that the character at the current
			// node is i.
			
			// hence, given parent state 'p', we optimize 
			// max_i pr (p->i) childVector [i] and store it in the p cell of vector childVector
			
			_Parameter overallMax					  = 0.0;
			
			long	   *stateBuffer					  = isLeaf?leafBuffer:stateCache;
			
			// check for degeneracy
			
			long howManyOnes = 0;
			for (long k = 0; k < alphabetDimension; k++)
				howManyOnes += childVector[k]==1.;
			
			if (howManyOnes == alphabetDimension)
			{
				for (long k = 0; k < alphabetDimension; k++)
					stateBuffer[k] = -1;
			}
			else
			{
				for (long p = 0; p < alphabetDimension; p++)
				{
					_Parameter max_lik = 0.;
					long	   max_idx = 0;
					
					for (long c = 0; c < alphabetDimension; c++)
					{
						_Parameter thisV = tMatrix[c] * childVector[c];
						if (thisV > max_lik)
						{
							max_lik = thisV;
							max_idx = c;
						}
					}
					
					stateBuffer [p] = max_idx;
					buffer [p]      = max_lik;
					
					if (max_lik > overallMax)
						overallMax = max_lik;

					tMatrix += alphabetDimension;
				}
				
				if (overallMax > 0.0 && overallMax < _lfScalingFactorThreshold)
				{
					for (long k = 0; k < alphabetDimension; k++)
						buffer[k] *= _lfScalerUpwards;
				}
				
				// buffer[p] now contains the maximum likelihood of the tree 
				// from this point forward given that parent state is p
				// and stateBuffer[p] stores the maximimizing assignment 
				// for this node
				
				for (long k = 0; k < alphabetDimension; k++)
				{
					long stateValue = stateBuffer[k];
					if (stateValue >= 0)
						parentConditionals[k] *= buffer[k]; 
				}
			}

			if (!isLeaf)
				stateCache += alphabetDimension;
			
			childVector += alphabetDimension;
		}
	}
	
	_List	   *result = new _List;
	for (long k = 0; k < iNodeCount; k++)
		result->AppendNewInstance (new _String(siteCount*unitLength,false));
	
	_Parameter _hprestrict_	* rootConditionals = iNodeCache + alphabetDimension * ((iNodeCount-1)  * patternCount);
	_SimpleList  parentStates (iNodeCount,0,0),
				 conversion;
	
	stateCache -= patternCount*(iNodeCount-1)*alphabetDimension;
	
	_AVLListXL	  conversionAVL (&conversion);
	_String		  codeBuffer	(unitLength, false);
	
	for (long siteID = 0; siteID < patternCount; siteID++, rootConditionals += alphabetDimension)
	{
		_Parameter max_lik = 0.;
		long	   max_idx = 0;
		
		long howManyOnes = 0;
		for (long k = 0; k < alphabetDimension; k++)
			howManyOnes += rootConditionals[k]==1.;
		
		_SimpleList*	patternMap = (_SimpleList*) expandedSiteMap (siteOrdering.lData[siteID]);

		if (howManyOnes != alphabetDimension)
		{
			for (long c = 0; c < alphabetDimension; c++)
			{
				_Parameter thisV = theProbs[c] * rootConditionals[c];
				if (thisV > max_lik)
				{
					max_lik = thisV;
					max_idx = c;
				}
			}
		
			parentStates.lData[iNodeCount-1] = max_idx;
			for  (long nodeID = iNodeCount-2; nodeID >=0 ; nodeID--)
			{
				long parentState = parentStates.lData[flatParents.lData [nodeID+flatLeaves.lLength]];
				if (parentState == -1)
					parentStates.lData[nodeID] = -1;
				else
					parentStates.lData[nodeID] = stateCache[(patternCount*nodeID+siteID)*alphabetDimension + parentState];
			}
		}
		else
			parentStates.Populate(iNodeCount,-1,0);
		
		for  (long nodeID = 0; nodeID < iNodeCount ; nodeID++)
		{
			dsf->ConvertCodeToLettersBuffered (dsf->CorrectCode(parentStates.lData[nodeID]), unitLength, codeBuffer.sData, &conversionAVL);
			_String  *sequence   = (_String*) (*result)(postToIn.lData[nodeID]);
			
			for (long site = 0; site < patternMap->lLength; site++)
			{
				char* storeHere = sequence->sData + patternMap->lData[site]*unitLength;
				for (long charS = 0; charS < unitLength; charS ++)
					storeHere[charS] = codeBuffer.sData[charS];
			}
		}
	}
		
	delete stateCache; delete leafBuffer;
	delete buffer;
	
	return result;
}
//_______________________________________________________________________________________________

void	 _TheTree::SetupCategoryMapsForNodes (_List& containerVariables, _SimpleList& classCounter, _SimpleList& multipliers)
{
	_CalcNode* iterator = DepthWiseTraversal (true);
	while (iterator)
	{
		iterator->SetupCategoryMap (containerVariables,classCounter,multipliers);
		iterator = DepthWiseTraversal(false);
	}
}
//_______________________________________________________________________________________________

void	 _CalcNode::SetupCategoryMap (_List& containerVariables, _SimpleList& classCounter, _SimpleList& multipliers)
{
	
	long	totalCategories = classCounter.Element(-1),
			globalCatCount  = containerVariables.lLength-1,
			localCategories = 1,
			catCount		= categoryVariables.lLength-1,
			entriesPerCat   = 2+catCount; 

	//for (long k = 0; k<containerVariables.lLength;k++)
	//	printf ("%d %s\n", k, ((_Variable*)containerVariables(k))->GetName()->sData);

	if (catCount<0)
		remapMyCategories.Clear();
	else
	{
		
		remapMyCategories.Populate (totalCategories*entriesPerCat,0,0);
		
		_SimpleList		remappedIDs,
						rateMultiplers (categoryVariables.lLength,1,0),
						categoryValues (globalCatCount+1,0,0);
		
		for (long myCatID = 0; myCatID <= catCount; myCatID++)
		{
			long coordinate = containerVariables.FindPointer(LocateVar(categoryVariables.lData[myCatID]));
			if (coordinate < 0)
				WarnError ("Internal error in SetupCategoryMap. Please report to spond@ucsd.edu");
			localCategories *= classCounter.lData[coordinate];
			//printf ("%d %s\n", myCatID, LocateVar(categoryVariables.lData[myCatID])->GetName()->sData);
			remappedIDs << coordinate;
		}
		
		for (long myCatID = catCount-1; myCatID >= 0; myCatID--)
			rateMultiplers.lData[myCatID] = rateMultiplers.lData[myCatID]*classCounter.lData[remappedIDs.lData[myCatID+1]];
		
		for	(long currentRateCombo  = 0; currentRateCombo < totalCategories; currentRateCombo++)
		{
			long copyRateCombo = currentRateCombo;
			for (long variableID = 0; variableID <= globalCatCount; variableID++)
			{
				categoryValues.lData[variableID] = copyRateCombo / multipliers.lData[variableID];
				copyRateCombo = copyRateCombo%multipliers.lData[variableID];
				//printf ("%d %d %d %d\n", currentRateCombo, variableID, multipliers.lData[variableID], categoryValues.lData[variableID]);
			}
			
			long localCatID = 0;
			
			for  (long localVariableID = 0; localVariableID<=catCount; localVariableID++)
				localCatID += rateMultiplers.lData[localVariableID] * categoryValues.lData[remappedIDs.lData[localVariableID]];
			
			long offset = currentRateCombo * entriesPerCat;
			remapMyCategories.lData[offset] = localCatID;

			offset++;
			for  (long localVariableID = 0; localVariableID<=catCount; localVariableID++)
				remapMyCategories[offset++] = categoryValues.lData[remappedIDs.lData[localVariableID]];

		}
	}
	
	//printf ("Node remap at %s yielded %s\n", GetName()->sData, _String((_String*)remapMyCategories.toStr()).sData);
	
}

//_______________________________________________________________________________________________

_Parameter	 _TheTree::Process3TaxonNumericFilter (_DataSetFilterNumeric* dsf, long catID)	
{
	
	_Parameter *l0 =  dsf->probabilityVectors.theData + 
					  dsf->categoryShifter * catID + dsf->theNodeMap.lData[0]*dsf->shifter,
				*l1 = dsf->probabilityVectors.theData + 
					  dsf->categoryShifter * catID + dsf->theNodeMap.lData[1]*dsf->shifter,
				*l2 = dsf->probabilityVectors.theData + 
					  dsf->categoryShifter * catID + dsf->theNodeMap.lData[2]*dsf->shifter,
				* matrix0 = ((_CalcNode*)(LocateVar(theRoot->nodes.data[0]->in_object)))->GetCompExp(catID)->theData,
				* matrix1 = ((_CalcNode*)(LocateVar(theRoot->nodes.data[1]->in_object)))->GetCompExp(catID)->theData,
				* matrix2 = ((_CalcNode*)(LocateVar(theRoot->nodes.data[2]->in_object)))->GetCompExp(catID)->theData,
				  overallResult = 0.;
	
	long		patternCount =  dsf->NumberDistinctSites();
	
	_Parameter  currentAccumulator = 1.;
	
	for (long patternIndex = 0; patternIndex < patternCount; patternIndex ++, l0+=4, l1+=4, l2+=4)
	{
		_Parameter rp0 = l0[0] * matrix0[0]+ l0[1]  * matrix0[1]  + l0[2] * matrix0[2]  + l0[3] * matrix0[3];
		_Parameter rp1 = l0[0] * matrix0[4]+ l0[1]  * matrix0[5]  + l0[2] * matrix0[6]  + l0[3] * matrix0[7];
		_Parameter rp2 = l0[0] * matrix0[8]+ l0[1]  * matrix0[9]  + l0[2] * matrix0[10] + l0[3] * matrix0[11];
		_Parameter rp3 = l0[0] * matrix0[12]+ l0[1] * matrix0[13] + l0[2] * matrix0[14] + l0[3] * matrix0[15];
		
		rp0 *= l1[0] * matrix1[0] + l1[1] * matrix1[1]  + l1[2] * matrix1[2]  + l1[3] * matrix1[3];
		rp1 *= l1[0] * matrix1[4] + l1[1] * matrix1[5]  + l1[2] * matrix1[6]  + l1[3] * matrix1[7];
		rp2 *= l1[0] * matrix1[8] + l1[1] * matrix1[9]  + l1[2] * matrix1[10] + l1[3] * matrix1[11];
		rp3 *= l1[0] * matrix1[12]+ l1[1] * matrix1[13] + l1[2] * matrix1[14] + l1[3] * matrix1[15];
		
		rp0 *= l2[0] * matrix2[0] + l2[1] * matrix2[1]  + l2[2] * matrix2[2]  + l2[3] * matrix2[3];
		rp1 *= l2[0] * matrix2[4] + l2[1] * matrix2[5]  + l2[2] * matrix2[6]  + l2[3] * matrix2[7];
		rp2 *= l2[0] * matrix2[8] + l2[1] * matrix2[9]  + l2[2] * matrix2[10] + l2[3] * matrix2[11];
		rp3 *= l2[0] * matrix2[12]+ l2[1] * matrix2[13] + l2[2] * matrix2[14] + l2[3] * matrix2[15];
	
		_Parameter 	result = theProbs[0]*rp0+
							 theProbs[1]*rp1+
							 theProbs[2]*rp2+
							 theProbs[3]*rp3;
	
	
		if (result<=0.0) 
		{
			return -A_LARGE_NUMBER;
		}
		
		long patternFreq = dsf->theFrequencies[patternIndex];
		for  (long freqIterator = 0; freqIterator < patternFreq; freqIterator++)
		{
			_Parameter tryMultiplication = currentAccumulator*result;
			if (tryMultiplication > 1.e-300)
				currentAccumulator = tryMultiplication;
			else
			{
				overallResult += myLog (currentAccumulator);
				currentAccumulator = result;
			}
		}
	}
	return overallResult + myLog (currentAccumulator);
}
//__________________________________________________________________________________

//__________________________________________________________________________________

void		 _treePSWRepresentation (node<long>* aNode, Ptr extraData)
{
	
}

//__________________________________________________________________________________

_Matrix*	 _TreeTopology::SplitsIdentity (_PMathObj p)
// compare tree topologies
{
	_Matrix * result = (_Matrix*) checkPointer(new _Matrix (2,1,false,true));
	_Constant * bc = (_Constant*) BranchCount ();
	result->theData[0] = bc->Value();
	result->theData[1] = -1;
	
	if (p && (p->ObjectClass() == TOPOLOGY || p->ObjectClass() == TREE))
	{
		DepthWiseT (true, _treePSWRepresentation, nil);
		while (currentNode)
			DepthWiseT (true, _treePSWRepresentation, nil);
	}
	
	DeleteObject (bc);
	return result;
}


#endif