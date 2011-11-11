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

#ifndef _AVLLIST_
#define _AVLLIST_
//#pragma once
#include "simplelist.h"

#define  MEMORYSTEP 8

//_____________________________________________________________________________
class _AVLList: public BaseObj
{

    public:

        //Data Members
        _SimpleList *dataList,
                    leftChild,
                    rightChild,
                    balanceFactor,
                    emptySlots;

        long root;

        //Methods
        _AVLList(_SimpleList*);

        virtual ~_AVLList(void){}
        virtual void Clear(bool = false);
        virtual bool HasData(long);

        virtual void ReorderList(_SimpleList* = nil);
        virtual long InsertData(BaseRef, long, bool);
        virtual BaseRef toStr(void);
        virtual long Traverser(_SimpleList&, long &, long = -1);
        virtual long GetRoot(void){return root;}
        virtual void DeleteXtra(long){};
        virtual void DeleteAll(bool cL){
            Clear(cL);        
            DeleteObject(dataList);
        }

        unsigned long countitems(void);

        // 20100623: a shortcut function to look for integers only
        // avoids calling ::Compare

        long Find(BaseRef);
        long Find(BaseRef,_SimpleList&);
        long FindLong(long);
        char FindBest(BaseRef, long&);

        long Next(long,_SimpleList&);
        long Prev(long,_SimpleList&);
        long First(void);
        long Last(void);

        long GetByIndex(const long);

        // the 1st bool flag is to say whether to dup the object being inserted
        // the 2nd bool flag (if the first flag is false) if set to true,
        // will cause failed inserts (key already exists) to delete the key
        long Insert(BaseRef,long = 0,bool=true,bool=false);

        BaseRef Retrieve(long);

        void Delete(BaseRef,bool=false);
        void ConsistencyCheck(void);

};

#endif
