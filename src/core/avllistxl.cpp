/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2009
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon              (apoon@cfenet.ubc.ca)

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

_AVLList    structure inspired by the excellent documentation of
GNU libavl 2.0.1 by Ben Pfaff (http://www.msu.edu/~pfaffben/avl/index.html)

*/

#include "hy_strings.h"
#include "errorfns.h"
#include "parser.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

//______________________________________________________________

_AVLListXL::_AVLListXL (_SimpleList* d):_AVLList(d)
{
}

//______________________________________________________________

BaseRef _AVLListXL::GetXtra (long d)
{
    return xtraD(d);
}

//______________________________________________________________

void    _AVLListXL::SetXtra (long i, BaseRef d, bool dup)
{
    xtraD.Replace (i,d, dup);
}


//______________________________________________________________

BaseRef _AVLListXL::toStr (void)
{
    _String * str = new _String (128L, true);
    checkPointer (str);

    if (countitems() == 0) {
        (*str) << "Empty Associative List";
    } else {
        _SimpleList  hist;
        long         ls, cn;

        cn = Traverser (hist,ls,root);

        while (cn>=0) {
            _String * keyVal = (_String*)Retrieve (cn);
            (*str) << keyVal;
            (*str) << " : ";
            (*str) << (_String*)GetXtra (cn);
            (*str) << '\n';
            cn = Traverser (hist,ls);
        }
    }

    str->Finalize();
    return str;
}


//______________________________________________________________

void _AVLListXL::Clear (bool cL)
{
    xtraD.Clear();
    _AVLList::Clear(cL);
}



//______________________________________________________________

long  _AVLListXL::InsertData (BaseRef b, long xl, bool cp)
{
    long w = (long)emptySlots.lLength - 1,
         n;

    BaseRef x = (BaseRef)xl;

    if (w>=0) {
        n = emptySlots.lData[w];
        emptySlots.Delete (w);
        leftChild.lData[n] = -1;
        rightChild.lData[n] = -1;
        balanceFactor.lData[n] = 0;
        ((BaseRef*)xtraD.lData)[n] = x;
        if (cp) {
            x->nInstances++;
        }
        ((BaseRef*)dataList->lData)[n] = b;
    } else {
        n = dataList->lLength;
        dataList->InsertElement (b,-1,false,false);
        leftChild  << -1;
        rightChild << -1;
        balanceFactor << 0;
        xtraD << x;
        if (!cp) {
            x->nInstances--;
        }
    }
    return n;
}



//______________________________________________________________

void _AVLListXL::DeleteXtra (long i)
{
    DeleteObject (((BaseRef*)xtraD.lData)[i]);
    (((BaseRef*)xtraD.lData)[i]) = nil;
}

