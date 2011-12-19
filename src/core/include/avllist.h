/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-now
Core Developers:
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon    (apoon@cfenet.ubc.ca)
  Steven Weaver (sweaver@ucsd.edu)
  
Module Developers:
	Lance Hepler (nlhepler@gmail.com)
	Martin Smith (martin.audacis@gmail.com)

Significant contributions from:
  Spencer V Muse (muse@stat.ncsu.edu)
  Simon DW Frost (sdf22@cam.ac.uk)

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

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
