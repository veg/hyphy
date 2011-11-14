/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2009
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon              (apoon@cfenet.ubc.ca)

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

_AVLList    structure inspired by the excellent documentation of
GNU libavl 2.0.1 by Ben Pfaff (http://www.msu.edu/~pfaffben/avl/index.html)

*/

#include "avllistx.h"
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

_AVLListX::_AVLListX (_SimpleList* d):_AVLList(d)
{
}

//______________________________________________________________

long    _AVLListX::GetXtra (long d)
{
    return xtraD.lData[d];
}

//______________________________________________________________

void    _AVLListX::SetXtra (long i, long d)
{
    xtraD.lData[i] = d;
}

//______________________________________________________________

void _AVLListX::Clear (bool cL)
{
    xtraD.Clear();
    _AVLList::Clear(cL);
}

//______________________________________________________________

BaseRef _AVLListX::toStr (void)
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
            (*str) << _String(GetXtra (cn));
            (*str) << '\n';
            cn = Traverser (hist,ls);
        }
    }

    str->Finalize();
    return str;
}

//______________________________________________________________

long  _AVLListX::InsertData (BaseRef b, long d, bool)
{
    long w = (long)emptySlots.lLength - 1,
         n;

    if (w>=0) {
        n = emptySlots.lData[w];
        emptySlots.Delete (w);
        leftChild.lData[n] = -1;
        rightChild.lData[n] = -1;
        balanceFactor.lData[n] = 0;
        xtraD.lData[n] = d;
        ((BaseRef*)dataList->lData)[n] = b;
    } else {
        n = dataList->lLength;
        dataList->InsertElement (b,-1,false,false);
        leftChild  << -1;
        rightChild << -1;
        balanceFactor << 0;
        xtraD << d;
    }
    return n;
}

//______________________________________________________________

void _AVLListX::DeleteXtra (long i)
{
    xtraD.lData[i] = -1;
}

//______________________________________________________________

void _AVLListX::PopulateFromList (_List& src)
{
    Clear(false);
    for (long k = 0; k < src.lLength; k++) {
        Insert (src(k)->makeDynamic(),k,false);
    }
}

