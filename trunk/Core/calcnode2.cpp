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

#ifdef	_SLKP_LFENGINE_REWRITE_

_Parameter			_lfScalerUpwards		  = pow(2.,100.),
					_lfScalingFactorThreshold = 1./_lfScalerUpwards,
					_logLFScaler			  = 100.*log(2.);

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
	long nt = MIN(tc, matrixQueue.lLength / 19 + 1);
#endif
	
#pragma omp parallel for default(shared) private (matrixID) schedule(static) if (cBase >= 20)  num_threads (nt)
	for  (matrixID = 0; matrixID < matrixQueue.lLength; matrixID++)
		((_CalcNode*) nodesToDo(matrixID))->SetCompExp (((_Matrix*)matrixQueue(matrixID))->Exponentiate(), catID);	
}

/*----------------------------------------------------------------------------------------------------------*/

void		_TheTree::DetermineNodesForUpdate	(_SimpleList& updateNodes, _List* expNodes, long catID)
{
	_SimpleList		nodesToUpdate (flatLeaves.lLength + flatTree.lLength - 1, 0, 0); 
	_CalcNode		*currentTreeNode;

	// look for nodes with model changes and mark the path up to the root as needing an update
	
	for (long nodeID = 0; nodeID < nodesToUpdate.lLength; nodeID++)
	{
		bool	isLeaf	   = nodeID < flatLeaves.lLength;

		currentTreeNode = isLeaf? (((_CalcNode**) flatCLeaves.lData)[nodeID]):
								  (((_CalcNode**) flatTree.lData)  [nodeID - flatLeaves.lLength]);
		
		if (currentTreeNode->NeedToExponentiate (catID))
		{
			if (expNodes)
				(*expNodes) << currentTreeNode;
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
																		_Parameter&			overallScaler,
																		long				siteFrom,
																		long				siteTo,
																		long				catID,
																		_Parameter*			storageVec
																 )
// the updateNodes flags the nodes (leaves followed by inodes in the same order as flatLeaves and flatNodes)
// that must be recomputed
{
	// process the leaves first 
	
	_SimpleList		taggedInternals					(flatNodes.lLength, 0, 0);
	long			alphabetDimension	  =			theFilter->GetDimension(),
					siteCount			  =			theFilter->NumberDistinctSites(),
					alphabetDimensionmod4  =		alphabetDimension-alphabetDimension%8;
	
	_CalcNode		*currentTreeNode;
	_Parameter		localScalerChange	  =			0.;
	
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
			_Parameter		
#ifdef __GNUC__ 
__restrict 
#endif
 *localScalingFactor	  = scalingAdjustments + parentCode*siteCount;
			
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
		
		currentTreeNode = isLeaf? ((_CalcNode*) flatCLeaves (nodeCode)):
								  ((_CalcNode*) flatTree    (nodeCode));
						
		_Parameter  *		
#ifdef __GNUC__ 
__restrict 
#endif
 transitionMatrix = currentTreeNode->GetCompExp(catID)->theData;
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

			bool		didScale = false;
			
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
				parentConditionals [0] *= tMatrix[0] * childVector[0]		+tMatrix[1] * childVector[1]		+tMatrix[2] * childVector[2]		+tMatrix[3] * childVector[3];
				parentConditionals [1] *= tMatrix[4] * childVector[0]		+tMatrix[5] * childVector[1]		+tMatrix[6] * childVector[2]		+tMatrix[7] * childVector[3];
				parentConditionals [2] *= tMatrix[8] * childVector[0]		+tMatrix[9] * childVector[1]		+tMatrix[10] * childVector[2]		+tMatrix[11] * childVector[3];
				parentConditionals [3] *= tMatrix[12] * childVector[0]		+tMatrix[13] * childVector[1]		+tMatrix[14] * childVector[2]		+tMatrix[15] * childVector[3];
				
				// handle scaling if necessary 
				// the check for sum > 0.0 is necessary for 'degenerate' log-L functions (-infinity)
				
				sum		= parentConditionals [0] + parentConditionals [1] + parentConditionals [2] + parentConditionals [3];
				if (sum < _lfScalingFactorThreshold && sum > 0.0)
				{
					scalingAdjustments [parentCode*siteCount + siteID] *= _lfScalerUpwards;
					parentConditionals [0]							   *= _lfScalerUpwards;
					parentConditionals [1]							   *= _lfScalerUpwards;
					parentConditionals [2]							   *= _lfScalerUpwards;
					parentConditionals [3]							   *= _lfScalerUpwards;
					localScalerChange								   += _logLFScaler * theFilter->theFrequencies [siteOrdering.lData[siteID]];
					didScale											= true;
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
						localScalerChange								   -= _logLFScaler * theFilter->theFrequencies [siteOrdering.lData[siteID]];
						didScale											= true;
					}
				}
				
				childVector += 4;
			}
			else
			{
				_Parameter sum = 0.0;
				
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
				if (sum < _lfScalingFactorThreshold && sum > 0.0)
				{
					scalingAdjustments [parentCode*siteCount + siteID] *= _lfScalerUpwards;
					for (long c = 0; c < alphabetDimension; c++) 
						parentConditionals [c] *= _lfScalerUpwards;
					localScalerChange									   += _logLFScaler * theFilter->theFrequencies [siteOrdering.lData[siteID]];
					didScale											    = true;
			}
				else
				{
					if (sum > _lfScalerUpwards)
					{
						scalingAdjustments [parentCode*siteCount + siteID] *= _lfScalingFactorThreshold;
						for (long c = 0; c < alphabetDimension; c++) 
							parentConditionals [c] *= _lfScalingFactorThreshold;
						localScalerChange								   -= _logLFScaler * theFilter->theFrequencies [siteOrdering.lData[siteID]];
						didScale											= true;
					}
				}
				childVector	   += alphabetDimension;
			}

			if (tcc && didScale)
			{
				// look ahead to see if we need to correct for downstream cached nodes  
				long cparentTCCIIndex	=	parentTCCIIndex,
				cparentTCCIBit		=	parentTCCIBit;
				
				_Parameter				scM, lsA;
				if (sum > _lfScalerUpwards)
				{	scM = _lfScalingFactorThreshold; lsA = -_logLFScaler; }
				else
				{	scM = _lfScalerUpwards; lsA = _logLFScaler; }
					
				
				for (long sid = siteID + 1; sid < siteTo; sid++,cparentTCCIBit++)
				{
					if (cparentTCCIBit == _HY_BITMASK_WIDTH_)
					{
						cparentTCCIBit   = 0;
						cparentTCCIIndex ++;
					}
					
					if ((tcc->lData[cparentTCCIIndex] & bitMaskArray.masks[cparentTCCIBit]) > 0)
					{
						scalingAdjustments [parentCode*siteCount + sid] *= scM;
						localScalerChange								+= lsA * theFilter->theFrequencies [siteOrdering.lData[sid]];
					}
					else
						break;
				}
			}
		}
	}
		
	// assemble the entire likelihood
	
	_Parameter 
#ifdef __GNUC__ 
__restrict 
#endif
 * rootConditionals = iNodeCache + alphabetDimension * (siteFrom + (flatTree.lLength-1)  * siteCount),
			   result = 0.0;
	

	for (long siteID = siteFrom; siteID < siteTo; siteID++)
	{
		_Parameter accumulator = 0.;
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
	
	if (localScalerChange > 0.1 || localScalerChange < -0.1) 
	{
#pragma omp atomic
		overallScaler += localScalerChange; 
	}
	
	return result;
}

#endif