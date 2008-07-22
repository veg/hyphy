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

/*--------------------------------------------------------------------------------------------------*/

void	_LikelihoodFunction::DetermineLocalUpdatePolicy (void)
{
	for (long k = 0; k < theTrees.lLength; k ++)
	{
		long catCount = ((_TheTree*)LocateVar(theTrees(k)))->categoryCount;
		_List * lup = new _List,
			  * mte = new _List;
		
		computedLocalUpdatePolicy.AppendNewInstance (new _SimpleList (catCount,0,0));
		for (long l = 0; l < catCount; l++)
		{
			lup->AppendNewInstance (new _SimpleList);
			mte->AppendNewInstance (new _List);
		}
			
		localUpdatePolicy.AppendNewInstance		 (lup);
		matricesToExponentiate.AppendNewInstance (mte);
	}
}

/*--------------------------------------------------------------------------------------------------*/

void	_LikelihoodFunction::FlushLocalUpdatePolicy (void)
{
	computedLocalUpdatePolicy.Clear();
	localUpdatePolicy.Clear();
	matricesToExponentiate.Clear();
}

/*--------------------------------------------------------------------------------------------------*/

#endif