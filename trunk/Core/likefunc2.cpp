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

#include "likefunc.h"

#ifdef	_SLKP_LFENGINE_REWRITE_

_Parameter		_LikelihoodFunction::ComputeTreeBlockByBranch	(long index, _TheTree* tree, _SimpleList& siteOrdering, _SimpleList& updateLeaves, _SimpleList& updateInternals)
// the updatePolicy flags the nodes (leaves followed by inodes in the same order as flatLeaves and flatNodes)
// that must be recomputed
{
	// process the leaves first 
	long			nodeCount	   = tree->GetLeafCount();
	_SimpleList		taggedInternals (tree->GetINodeCount(), 0, 0);
	
	for  (long leafID = 0; leafID < updateLeaves.lLength; leafID++)
	{
		long parentCode = tree->flatParents.lData[updateLeaves.lData[leafID]];
	}
}

#endif