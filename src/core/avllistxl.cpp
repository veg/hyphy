/*
 
 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
 Art FY Poon    (apoon42@uwo.ca)
 Steven Weaver (sweaver@temple.edu)
 
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
#include "hy_string_buffer.h"
#include "parser.h"

#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>

//______________________________________________________________

_AVLListXL::_AVLListXL (_SimpleList* d):_AVLList(d){
}

//______________________________________________________________

_AVLListXL::_AVLListXL (_SimpleList* d, _List&& payload):_AVLList(d) {
    for (int i = 0; i < payload.countitems(); i+=2) {
        payload.GetItem(i)->AddAReference();
        payload.GetItem(i+1)->AddAReference();
        PushPair((_String*)payload.GetItem(i), payload.GetItem(i+1));
    }
}


//______________________________________________________________

BaseRef _AVLListXL::GetXtra (long d) const {
    return xtraD.GetItem(d);
}


//______________________________________________________________

BaseRef _AVLListXL::GetDataByKey(BaseRefConst key) const {
    long f = Find (key);
    if (f < 0L) {
        return nil;
    }
    return GetXtra(f);
}


//______________________________________________________________

void    _AVLListXL::SetXtra (long i, BaseRef d, bool dup) {
    xtraD.Replace (i,d, dup);
}


//______________________________________________________________

BaseRef _AVLListXL::toStr (unsigned long)
{
    _StringBuffer * str = new _StringBuffer (128L);

    if (countitems() == 0) {
        (*str) << "Empty Associative List";
    } else {
        _SimpleList  hist;
        long         ls, cn;

        cn = Traverser (hist,ls,root);

        while (cn>=0) {
            (*str) << (_String*)Retrieve (cn)
                   << " : ";
            str->AppendNewInstance((_String*)GetXtra (cn)->toStr());
            (*str) << '\n';
          
            cn = Traverser (hist,ls);
        }
    }

    return str;
}
//______________________________________________________________

_AVLListXL&  _AVLListXL::PushPairCopyKey  (_String const key, BaseRef data) {
    UpdateValue (new _String (key), data, false, false);
    return *this;
}
//______________________________________________________________

_AVLListXL&  _AVLListXL::PushPair          (_String* key, BaseRef data) {
    if (UpdateValue (key, data, false, false) >= 0) {
        DeleteObject (key);
    }
    return *this;
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
        n = emptySlots.list_data[w];
        emptySlots.Delete (w);
        leftChild.list_data[n] = -1;
        rightChild.list_data[n] = -1;
        balanceFactor.list_data[n] = 0;
        ((BaseRef*)xtraD.list_data)[n] = x;
        if (cp) {
            x->AddAReference();
        }
        ((BaseRef*)dataList->list_data)[n] = b;
    } else {
        n = dataList->lLength;
        dataList->InsertElement (b,-1,false,false);
        leftChild  << -1;
        rightChild << -1;
        balanceFactor << 0;
        xtraD << x;
        if (!cp) {
            x->RemoveAReference();
        }
    }
    return n;
}



//______________________________________________________________

void _AVLListXL::DeleteXtra (long i) {
    DeleteObject (((BaseRef*)xtraD.list_data)[i]);
    (((BaseRef*)xtraD.list_data)[i]) = nil;
}

