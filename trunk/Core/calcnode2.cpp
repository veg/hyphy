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
#include "math.h"

#ifdef	_SLKP_LFENGINE_REWRITE_

_Parameter		_TheTree::ComputeTreeBlockByBranch	(
																		_SimpleList& siteOrdering, 
																		_SimpleList& updateNodes, 
																		_DataSetFilter* theFilter,
																		_Parameter* iNodeCache,
																		long	  * lNodeFlags,
																		_GrowingVector*
																				    lNodeResolutions,
																		long		catID
																 )
// the updatePolicy flags the nodes (leaves followed by inodes in the same order as flatLeaves and flatNodes)
// that must be recomputed
{
	// process the leaves first 
	_SimpleList		taggedInternals			(flatNodes.lLength, 0, 0);
	long			alphabetDimension =		theFilter->GetDimension(),
					siteCount		  =		theFilter->NumberDistinctSites();
					
	_CalcNode		*currentTreeNode;
	
	for  (long nodeID = 0; nodeID < updateNodes.lLength; nodeID++)
	{
		long	nodeCode   = updateNodes.lData [nodeID],
				parentCode = flatParents.lData [nodeCode];
				
		bool	isLeaf	   = nodeCode < flatLeaves.lLength;
		
		if (!isLeaf)
			nodeCode -=  flatLeaves.lLength;
				
		_Parameter		* parentConditionals = iNodeCache + parentCode * alphabetDimension * siteCount;
		if (taggedInternals.lData[parentCode] == 0)
		// mark the parent for update and clear its conditionals if needed
		{
			taggedInternals.lData[parentCode] = 1;
			for (long k = 0; k < alphabetDimension*siteCount; parentConditionals[k] = 1., k++) ;
		}
		
		currentTreeNode = isLeaf? ((_CalcNode*) flatCLeaves (nodeCode)):
								  ((_CalcNode*) flatTree    (nodeCode));
						
		currentTreeNode->RecomputeMatrix (catID, categoryCount, nil);
		_Parameter *		transitionMatrix = currentTreeNode->GetCompExp(catID)->theData;
		_Parameter *		childVector;
		
		for (long siteID = 0; siteID < siteCount; siteID++)
		{
			_Parameter		*tMatrix = transitionMatrix;
			if (isLeaf)
			{
				long siteState = lNodeFlags[nodeCode*siteCount + siteID] ;
				if (siteState >= 0)
				/* a single character state; sweep down the appropriate column */
				{
					tMatrix  +=  siteState;
					for (long k = 0; k < alphabetDimension; parentConditionals[k] *= *tMatrix, k++, tMatrix += alphabetDimension)  ;
					parentConditionals += alphabetDimension;
					continue;
				}
				else
					childVector = lNodeResolutions->theData + (-siteState-1) * alphabetDimension;
			}
			else
				childVector = iNodeCache + (nodeCode * siteCount + siteID ) * alphabetDimension;
						
			for (long p = 0; p < alphabetDimension; p++)
			{
				_Parameter		accumulator = 0.0;
				for (long c = 0; c < alphabetDimension; c++, tMatrix++)
					accumulator +=  *tMatrix * childVector[c];
				
				parentConditionals[p] *= accumulator;
			}
			parentConditionals += alphabetDimension;
		}
	}
	
	// assemble the entire likelihood
	
	_Parameter * rootConditionals = iNodeCache + (flatTree.lLength-1) * alphabetDimension * siteCount,
			   result = 0.0;
	
	for (long siteID = 0; siteID < siteCount; siteID++)
	{
		_Parameter accumulator = 0.;
		for (long p = 0; p < alphabetDimension; p++,rootConditionals++)
			accumulator += *rootConditionals * theProbs[p];
		
		if (accumulator <= 0.0)
			return -A_LARGE_NUMBER;
		result += log(accumulator) * theFilter->theFrequencies [siteID];
	}
	
	return result;
}

#endif