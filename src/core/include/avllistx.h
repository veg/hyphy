/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2011
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

#ifndef _AVLLISTX_
#define _AVLLISTX_
//#pragma once
#include "list.h"
#include "avllist.h"

#define  MEMORYSTEP 8

//_____________________________________________________________________________

class _AVLListX: public _AVLList
{

    public:
        /* SLKP: 20090817
            add key: index values from the list of strings
         */

        //Data Members
        _SimpleList xtraD;

        //Methods
        _AVLListX(_SimpleList*);

        virtual ~_AVLListX(void){}
        virtual BaseRef toStr(void);

        virtual void Clear(bool = false);
        virtual void DeleteXtra(long);
        virtual void PopulateFromList(_List&);

        virtual long InsertData(BaseRef, long, bool);

        void SetXtra(long,long);
        long GetXtra(long);

};

#endif
