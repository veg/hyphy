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
        virtual long UpdateValue (BaseRef, long, long);

        void        SetXtra(long,long);
        long        GetXtra(long);

};

#endif
