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

#include "list.h"
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


/*
==============================================================
Constructors
==============================================================
*/

// Does nothing
_List::_List ()
{
}

// Length constructor
_List::_List (unsigned long l):_SimpleList(l)
{
}

// Stack Copy contructor
_List::_List (const _List& l, long from, long to)
{
    if (from == 0 && to == -1) { // copy the whole thing
        BaseRef   br = (BaseRef)&l;
        Duplicate (br);
    } else {
        Initialize           ();
        NormalizeCoordinates (from, to, l.lLength);

        for (long i = from; i <= to; i++) {
            (*this) << ((BaseRef*)l.lData)[i];
        }
    }
}

// Construct a list of substrings from the original string separated by char
_List::_List (BaseRef ss, char sep)
{
    _String* s = (_String*)ss;
    if (s->Length()!=0) {
        long cp=0,cpp;
        while ((cpp = s->Find(sep,cp,-1))!=-1) {
            if (cpp>cp) {
                AppendNewInstance (new _String(*s,cp,cpp-1));
            } else {
                AppendNewInstance (new _String);
            }
            cp=cpp+1;
        }

        AppendNewInstance (new _String(*s,cp,-1));
    }
}

// Data constructor (1 member list)
_List::_List (BaseRef br)
{
    lLength = 1;
    laLength = MEMORYSTEP;
    lData = (long*)MemAllocate (laLength * sizeof(Ptr));
    ((BaseRef*)lData)[0]= br->makeDynamic();
}

//Destructor
_List::~_List(void)
{
    if (nInstances<=1) {
        for (unsigned long i = 0; i<lLength; i++) {
            BaseRef t = ((BaseRef*)lData)[i];
            if (t) {
                if (t->nInstances<=1) {
                    DeleteObject(t);
                } else {
                    t->nInstances--;
                }
            }
        }
    }


}

/*
==============================================================
Operator Overloads
==============================================================
*/

// Element location functions (0,llength - 1)
BaseRef& _List::operator [] (long i)
{
    BaseRef t = BaseRef(_SimpleList::operator[] (i));
    if (t)
        if (t->nInstances>1) {
            t->nInstances--;
            ((BaseRef*)(lData))[i]=t->makeDynamic();
        }

    return ((BaseRef*)(lData))[i];
}

// Element location functions (0,llength - 1)
BaseRef _List::operator () (unsigned long i)
{
    return ((BaseRef*)lData)[i];
}

// Assignment operator
_List _List::operator = (_List& l)
{
    this->~_List();
    lLength = l.lLength;
    laLength = l.laLength;
    lData = l.lData;
    l.nInstances++;
    for (unsigned long i = 0; i<lLength; i++) {
        ((BaseRef*)(lData))[i]->nInstances++;
    }
    return *this;
}

// Assignment compare 
bool _List::operator == (_List& l2)
{
    return this->Equal(l2);
}

// Append operator
_List _List::operator & (_List& l)
{
    _List res (l.lLength + lLength);
    if (!res.laLength) {
        return res;
    }

    if (lData&&lLength) {
        memcpy(res.lData,lData,lLength*sizeof(void*));
    }
    if (l.lData&&l.lLength) {
        memcpy((char*)res.lData+lLength*sizeof(void*),l.lData, l.lLength*sizeof(void*));
    }
    res.lLength = l.lLength + lLength;
    unsigned long i;
    for (i = 0; (i<lLength); i++) {
        ((BaseRef*)lData)[i]->nInstances++;

    }
    for (i=0; i<l.lLength; i++) {
        ((BaseRef*)l.lData)[i]->nInstances++;
    }
    return res;
}

// Append operator
_List _List::operator & (BaseRef br)
{
    _List res (lLength+1);
    if (!res.laLength) {
        return res;
    }

    if (lData) {
        memcpy(res.lData,lData,lLength*sizeof(void*));
    }
    for (unsigned long i = 0; (i<lLength); i++) {
        ((BaseRef*)lData)[i]->nInstances++;
    }
    res.lLength=lLength+1;
    ((BaseRef*)res.lData)[lLength]=br->makeDynamic();
    return res;
}

void _List::operator && (BaseRef br)
{
    InsertElement (br);
}

void _List::operator && (const char * buffer)
{
    _String*        newString = new _String (buffer);
    checkPointer    (newString);
    InsertElement   (newString,-1,false);
    DeleteObject    (newString);
}

void _List::operator << (BaseRef br)
{
//  InsertElement (br, -1, false);
    lLength++;
    if (lLength>laLength) {
        unsigned long incBy = (MEMORYSTEP*5 > lLength)? MEMORYSTEP: lLength/5;

        laLength+=incBy;

        if (lData) {
            checkPointer (lData = (long*)MemReallocate((char*)lData, laLength*sizeof(void*)));
        } else {
            checkPointer (lData = (long*)MemAllocate(laLength*sizeof(void*)));
        }
    }

    ((BaseRef*)lData)[lLength-1]=br;
    br->nInstances++;
}

void _List::operator << (_List& source)
{
    for (long k=0; k<source.lLength; k++) {
        (*this) << ((BaseRef*)source.lData)[k];
    }
}


/*
==============================================================
Methods
==============================================================
*/

void _List::AppendNewInstance (BaseRef br)
{
    if (br) {
        (*this)<<br;
        br->nInstances--;
    } else {
        checkPointer (br);
    }
}

long  _List::BinaryFind (BaseRef s)
{
    _String * st = (_String*)s;
    long top=lLength-1, bottom=0, middle;

    if (top==-1) {
        return -1;
    }
    while (top>bottom) {
        middle = (top+bottom)/2;
        _String* stp = (_String*)(((BaseRef*)lData)[middle]->toStr());

        int      cres = st->Compare (stp);
        DeleteObject (stp);

        if (cres < 0) {
            top = middle==top?top-1:middle;
        } else if (cres == 0) {
            return middle;
        } else {
            bottom = middle==bottom?bottom+1:middle;
        }

    }
    middle = top;
    _String* stp=(_String*)(((BaseRef*)lData)[middle]->toStr());
    if (st->Equal(stp)) {
        DeleteObject(stp);
        return middle;
    }
    DeleteObject(stp);
    return -middle-2;
}

long  _List::BinaryInsert (BaseRef s)
{

    if (!lLength) {
        InsertElement (s,0,true);
        return 0;
    }

    long pos = -BinaryFind (s)-2;
    if (pos<0) {
        return -pos+2;
    }
    _String *s1 = (_String*)s->toStr(), *s2 =(_String*) ((*this)(pos))->toStr();
    if (*s2<*s1) {
        pos++;
    }
    DeleteObject(s1);
    DeleteObject(s2);
    InsertElement (s,pos,true);
    return pos>=lLength?lLength-1:pos;

}

void    _List::bumpNInst (void)
{
    for (unsigned long i = 0; i<lLength; i++) {
        ((BaseRef*)lData)[i]->nInstances++;
    }
}

void  _List::Clear (bool completeClear)
{
    if (nInstances<=1) {
        for (unsigned long i = 0; i<lLength; i++) {
            DeleteObject (((BaseRef*)lData)[i]);
        }
        _SimpleList::Clear(completeClear);

    } else {
        nInstances--;
    }
}

long  _List::Compare (long i, long j)
{
    _String             *si = (_String*)lData[i],
                         *sj = (_String*)lData[j];

    return  si->Compare(sj);
}

long  _List::Compare (BaseRef i, long j)
{
    _String             *sj = (_String*)lData[j],
                         *si = (_String*)i;

    return  si->Compare(sj);
}

unsigned long _List::Count()
{
    return lLength;
}

// Delete item at index (>=0)
void  _List::Delete (long index)
{
    if ((index>=0)&&(index<lLength)) {
        BaseRef theObj = ((BaseRef*)lData)[index];
        DeleteObject (theObj);
        lLength--;
        if (lLength-index)
            for (unsigned long i = index; i < lLength; i++) {
                lData[i] = lData[i+1];
            }
        //memcpy ((Ptr)lData+sizeof(BaseRef)*(index),(Ptr)lData+sizeof(BaseRef)*(index+1),sizeof(BaseRef)*(lLength-index));
    }
    if (laLength-lLength>MEMORYSTEP) {
        laLength -= ((laLength-lLength)/MEMORYSTEP)*MEMORYSTEP;
        lData = (long*)MemReallocate ((char*)lData, laLength*sizeof(Ptr));
    }

}

// Delete item at index (>=0)
void  _List::DeleteList (const _SimpleList& toDelete)
{
    if (toDelete.lLength) {
        long k = 0;
        for (long i = 0; i<lLength; i++) {
            if (k<toDelete.lLength && i==toDelete.lData[k]) {
                DeleteObject (((BaseRef*)lData)[i]);
                //if (k<toDelete.lLength)
                k++;
            } else {
                ((BaseRef*)lData)[i-k] = ((BaseRef*)lData)[i];
            }
        }
        lLength -= toDelete.lLength;
        if (laLength-lLength>MEMORYSTEP) {
            laLength -= ((laLength-lLength)/MEMORYSTEP)*MEMORYSTEP;
            lData = (long*)MemReallocate ((char*)lData, laLength*sizeof(Ptr));
        }
    }
}

void    _List::Duplicate (const BaseRef theRef)
{
    _SimpleList::Duplicate (theRef);
    if (lData)
        for (unsigned long i = 0; i<lLength; i++) {
            if (((BaseRef*)lData)[i]) {
                (((BaseRef*)lData)[i])->nInstances++;
            }
        }
}

bool _List::Equal(_List& l2)
{
    if (lLength!=l2.lLength) {
        return false;
    }

    for (long i=0; i<lLength; i++)
        if (!((_String*)lData[i])->Equal ((_String*)l2.lData[i])) {
            return false;
        }

    return true;
}

long  _List::Find (BaseRef s, long startat)
{
    _String * st = (_String*)s;
    for (unsigned long i = startat; i<lLength; i++) {
        _String * sp = (_String*)(((BaseRef*)lData)[i]->toStr());

        if (st->Equal(sp)) {
            DeleteObject(sp);
            return i;
        }
        DeleteObject(sp);
    }
    return -1;
}

long  _List::FindString (BaseRef s, long startat, bool caseSensitive, long upTo)
{
    char * s1, *s2;
    long t = ((_String*)s)->sLength;

    if (upTo < 0 || upTo >= lLength) {
        upTo = lLength-1;
    }

    for (long i = startat; i<= upTo; i++) {
        s1 = ((_String*)s)->sData;
        if (((_String*)(((BaseRef*)lData)[i]))->sLength==t) {
            s2 = ((_String*)(((BaseRef*)lData)[i]))->sData;
            long j = 0;
            if (caseSensitive)
                for (j=0; (*s1==*s2)&&(j<t); j++,s1++,s2++) ;
            else
                for (j=0; (toupper(*s1)==toupper(*s2))&&(j<t); j++,s1++,s2++) ;

            if (j==t) {
                return i;
            }
        }
    }
    return -1;
}

long  _List::FreeUpMemory (long requestedBytes)
{
    long freed = 0;
    for (unsigned long i = 0; i<lLength; i++) {
        BaseRef t = ((BaseRef*)lData)[i];
        freed+=t->FreeUpMemory(requestedBytes-freed);
        if (freed>=requestedBytes) {
            return freed;
        }
    }
    return freed;
}

// Append & store operator
void _List::InsertElement (BaseRef br, long insertAt, bool store)
{
    _SimpleList::InsertElement (br, insertAt, store);
}

void    _List::Intersect (_List& l1, _List& l2, _SimpleList* idx, _SimpleList* idx2)
{
    // compute the union of two sorted lists
    // each repeat appears exactly once
    if (lLength) {
        Clear();
    }

    long  c1 = 0,
          c2 = 0;

    while (c1<l1.lLength && c2<l2.lLength) {
        while (c1<l1.lLength && ((_String*)l1(c1))->Compare((_String*)l2(c2))<0) {
            c1++;
        }
        if (c1==l1.lLength) {
            break;
        }

        while (c1<l1.lLength && c2<l2.lLength &&((_String*)l1(c1))->Equal((_String*)l2(c2))) {
            if (idx) {
                (*idx) << c1;
            }
            if (idx2) {
                (*idx2) << c2;
            }

            (*this) << l1(c1++);
            c2++;
        }
        if (c1==l1.lLength || c2==l2.lLength) {
            break;
        }
        while (c2<l2.lLength && ((_String*)l2(c2))->Compare((_String*)l1(c1))<0) {
            c2++;
        }
    }
}

BaseRef  _List::Join (BaseRef spacer)
{
    _String *joined = new _String (256L,true);

    for (long k = 0; k < lLength; k++) {
        if (k) {
            (*joined) << *(_String*)spacer;
        }
        joined->AppendNewInstance((_String*) ((BaseRef*)lData)[k]->toStr());
    }

    joined->Finalize();
    return joined;
}

void _List::Place (BaseRef br)
{
//  InsertElement (br, -1, false);
    lLength++;
    if (lLength>laLength) {
        laLength+=MEMORYSTEP;
        if (lData) {
            checkPointer (lData = (long*)MemReallocate((char*)lData, laLength*sizeof(void*)));
        } else {
            checkPointer (lData = (long*)MemAllocate(laLength*sizeof(void*)));
        }
    }
    ((BaseRef*)lData)[lLength-1]=br;
}

//TODO: makeDynamic should be MakeDynamic to follow convention.
BaseRef _List::makeDynamic(void)
{
    _List * Res = new _List;
    checkPointer(Res);
    //lData = nil;
    memcpy ((char*)Res, (char*)this, sizeof (_List));
    Res->nInstances = 1;
    Res->lData = nil;
    Res->Duplicate (this);
    return Res;
}

void  _List::Replace (long index, BaseRef newObj, bool dup)
{
    if ((index>=0)&&(index<lLength)) {
        BaseRef theObj = ((BaseRef*)lData)[index];
        DeleteObject (theObj);
        ((BaseRef*)lData)[index] = dup?newObj->makeDynamic():newObj;
    }
}

// Char* conversion
//TODO: toFileStr should be ToFileStr to follow convention.
void _List::toFileStr(FILE* dest)
{
    fprintf (dest,"{");

    for (unsigned long i = 0; (i<lLength); i++) {
        ((BaseRef*)lData)[i]->toFileStr(dest);
        if (i<lLength-1) {
            fprintf (dest,",");
        }
    }
    fprintf (dest,"}");
}

// Char* conversion
//TODO: toFileStr should be ToStr to follow convention.
BaseRef _List::toStr(void)
{
    _String * s = new _String((unsigned long)20*(lLength+1),true);

    checkPointer (s);

    (*s)<<'{';

    for (unsigned long i = 0; (i<lLength); i++) {
        _String* t = (_String*)(((BaseRef*)lData)[i]->toStr());
        if (t) {
            (*s)<<t;
            DeleteObject (t);
        }
        if (i<lLength-1) {
            (*s)<<',';
        }
    }
    (*s)<<'}';
    s->Finalize();
    return s;
}

