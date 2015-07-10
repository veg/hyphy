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

BaseRef _AVLListXL::GetDataByKey(BaseRef key)
{
    long f = Find (key);
    if (f < 0) {
        return nil;
    }
    return GetXtra(f);
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

long  _AVLListXL::UpdateValue(BaseRef b, BaseRef d, bool do_copy, bool copy_key) {
    long exists = Find (b);
    if (exists >= 0) {
        SetXtra (exists, d, do_copy);       
    } else {
        Insert (copy_key?b->makeDynamic():b,(long)d, do_copy);
    }
    return exists;
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

