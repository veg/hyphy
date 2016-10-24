/*
 
 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
 Art FY Poon    (apoon@cfenet.ubc.ca)
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
#include <stdarg.h>

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
_List::_List (BaseRef br) {
    lLength = 1;
    laLength = MEMORYSTEP;
    lData = (long*)MemAllocate (laLength * sizeof(Ptr));
    BaseRef   object_copy = br->makeDynamic();
    ((BaseRef*)lData)[0]= object_copy;
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
BaseRef _List::operator () (const unsigned long i)
{
    return ((BaseRef*)lData)[i];
}

// Element location functions (0,llength - 1)
BaseRef _List::GetItem (const unsigned long i) const {
    return ((BaseRef*)lData)[i];
}

  // Element location functions (0,llength - 1)
BaseRef _List::GetItemRangeCheck(const unsigned long i) const {
  return i < lLength ? ((BaseRef*)lData)[i] : nil;
}


// Assignment operator
const _List _List::operator = (const _List& l) {

    Clear(true);
  
    lLength = l.lLength;
    RequestSpace(laLength);
    for (unsigned long i = 0UL; i<lLength; i++) {
        ((BaseRef*)(lData))[i] = l.GetItem (i);
        ((BaseRef*)(lData))[i] -> AddAReference();
    }

    return *this;

}

// Assignment compare 
bool _List::operator == (_List const& l2) const
{
    return this->Equal(l2);
}

// Append operator
const _List _List::operator & (_List const & l) const
{
    _List res (l.lLength + lLength);
  
    if (!res.laLength) {
        return res;
    }

    BaseRef * res_data = (BaseRef*) res.lData;
  
    if (lData) {
      BaseRef * base_data = (BaseRef*) lData;
      
      for (res.lLength = 0UL; res.lLength < lLength; res.lLength++) {
        res_data [res.lLength] = base_data[res.lLength];
        res_data [res.lLength] -> AddAReference();
      }
    }
  
    if (l.lData) {
      BaseRef * l_data = (BaseRef*) l.lData;
      for (unsigned long li = 0UL; li < l.lLength; li++, res.lLength++) {
        res_data [res.lLength] = l_data [li];
        l_data [li] -> AddAReference();
      }
    }

    return res;
}

// Append operator
const _List _List::operator & (BaseRef br) const
{
    _List res (lLength+1UL);
    if (res.laLength == 0UL) {
        return res;
    }
  
    BaseRef * res_data = (BaseRef*) res.lData;
    if (lData) {
      BaseRef * base_data = (BaseRef*) lData;
      
      for (res.lLength = 0UL; res.lLength < lLength; res.lLength++) {
        res_data [res.lLength] = base_data[res.lLength];
        res_data [res.lLength] -> AddAReference();
      }
      
    }
  
    res_data [res.lLength++]=br->makeDynamic();
    return res;
}

_List& _List::operator && (BaseRef br) {
    InsertElement (br);
    return *this;
}

_List& _List::operator && (const char * buffer) {
    InsertElement   (new _String (buffer),-1,false,false);
    return *this;
}

_List& _List::operator < (BaseRef br) {
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
  return *this;
}

_List& _List::operator << (BaseRef br) {
    br->AddAReference();
    return (*this) < br;
}

_List& _List::operator << (_List const& source) {
  for (unsigned long k=0UL; k<source.lLength; k++) {
      (*this) << ((BaseRef*)source.lData)[k];
  }
  return *this;
}

_List& _List::operator < (_List const& source) {
  for (unsigned long k=0UL; k<source.lLength; k++) {
    (*this) << ((BaseRef*)source.lData)[k];
  }
  return *this;
}


_List& _List::operator < (const char * text) {
  return (*this) < new _String (text);
}

/*
==============================================================
Methods
==============================================================
*/

void _List::AppendNewInstance (BaseRef br) {
    if (br) {
        (*this)<<br;
        br->RemoveAReference();
    } else {
        checkPointer (br);
    }
}




long  _List::BinaryFindObject (BaseObj const * s, long startAt) const {
    _String const * st = (_String const*)s;
    long top    = lLength-1,
         bottom = startAt,
         middle;

    if (top < 0L) {
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

    long pos = -BinaryFindObject (s)-2;
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
        for (unsigned long i = 0UL; i<lLength; i++) {
            DeleteObject (((BaseRef*)lData)[i]);
        }
        _SimpleList::Clear(completeClear);

    } else {
        nInstances--;
    }
}

long  _List::Compare (long i, long j) const {
    _String             *si = (_String*)lData[i],
                         *sj = (_String*)lData[j];

    return  si->Compare(sj);
}

long  _List::Compare (BaseObj const * i, long j) const {
    _String const       *sj = (_String const*)lData[j],
                        *si = (_String const*)i;

    return  si->Compare(sj);
}

unsigned long _List::Count()
{
    return lLength;
}

// Delete item at index (>=0)
void  _List::Delete (long index, bool delete_object)
{
    if (index>=0 && index<lLength) {
        if (delete_object) {
            DeleteObject (((BaseRef*)lData)[index]);
        }
        lLength--;
        if (lLength-index)
            for (unsigned long i = index; i < lLength; i++) {
                lData[i] = lData[i+1];
            }
        //memcpy ((Ptr)lData+sizeof(BaseRef)*(index),(Ptr)lData+sizeof(BaseRef)*(index+1),sizeof(BaseRef)*(lLength-index));
    }
    if (laLength-lLength>MEMORYSTEP) {
        laLength -= ((laLength-lLength)/MEMORYSTEP)*MEMORYSTEP;
        if (laLength > 0)
          lData = (long*)MemReallocate ((char*)lData, laLength*sizeof(Ptr));
        else {
          free (lData); lData = nil;
        }
    }

}

// Delete item at index (>=0)
void  _List::DeleteList (const _SimpleList& toDelete)
{
    if (toDelete.lLength) {
        long k = 0;
        for (unsigned long i = 0; i<lLength; i++) {
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
            if (laLength > 0)
              lData = (long*)MemReallocate ((char*)lData, laLength*sizeof(Ptr));
            else {
              free (lData); lData = nil;
            }
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

bool _List::Equal(_List const & l2) const
{
    if (lLength!=l2.lLength) {
        return false;
    }

    for (unsigned long i=0; i<lLength; i++)
        if (!((_String*)lData[i])->Equal ((_String*)l2.lData[i])) {
            return false;
        }

    return true;
}

long  _List::FindObject (BaseRefConst s, long startat) const {
    _String const * st = (_String const*)s;
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
void _List::InsertElement (BaseRef br, long insertAt, bool store, bool pointer) {
    _SimpleList::InsertElement (br, insertAt, store, pointer);
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

BaseRef  _List::Join (BaseRef spacer, long startAt, long endAt)
{
    _String *joined = new _String (256L,true);
    if (endAt < 0) { 
        endAt = lLength; 
    } else if (endAt > lLength) {
        endAt = lLength;
    }

    for (unsigned long k = MAX(0,startAt); k < endAt; k++) {
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
    if (index>=0 && index<lLength) {
        DeleteObject (((BaseRef*)lData)[index]);
        ((BaseRef*)lData)[index] = dup?newObj->makeDynamic():newObj;
    }
}

// Char* conversion
//TODO: toFileStr should be ToFileStr to follow convention.
void _List::toFileStr(FILE* dest, unsigned long)
{
    fprintf (dest,"{");
  
    if (lLength) {
      GetItem(0)->toFileStr(dest);
      for (unsigned long i = 1UL; i<lLength; i++) {
          fprintf (dest,", ");
          GetItem(i)->toFileStr(dest);
      }
    }
    fprintf (dest,"}");
}

// Char* conversion
//TODO: toFileStr should be ToStr to follow convention.
BaseRef _List::toStr(unsigned long) {
    _String * s = new _String((unsigned long)20*(lLength+1),true);

    (*s)<<'{';
  
    if (lLength) {
      s->AppendNewInstance ((_String*)GetItem(0)->toStr());
      for (unsigned long i = 1UL; i<lLength; i++) {
          (*s)<< ", ";
          s->AppendNewInstance ((_String*)GetItem(i)->toStr());
      }
    }
  
    (*s)<<'}';
    s->Finalize();
    return s;
}

void    _List::Map (_List& target, _SimpleList& mapping) {
    mapping.Clear();
    if (lLength == 0) {
        return ;
    }
    _List     aux;
    _AVLListX theMapping (&aux);
    for (unsigned long t = 0; t < target.lLength; t++) {
        _String * s = (_String *)target.GetItem(t)->toStr();
        theMapping.Insert (s, t);
        //printf ("Target %ld:%s\n", t, s->sData);
    }
    
    mapping.Clear();
    for (unsigned long s = 0; s < lLength; s++) {
        _String * s_object = (_String*)GetItem (s)->toStr();
        //printf ("Source %ld : %s\n", s, s_object->sData);
        long      idx = theMapping.Find (s_object);
        if (idx >= 0) {
            mapping << theMapping.GetXtra (idx);
        } else {
            mapping << -1;
        }
        DeleteObject (s_object);
    }
}

const _String _List::GenerateUniqueNameForList (_String const& base, bool sorted) const {
  _String try_name (base);
  long    suffix = 1L;
  
  if (sorted)
    while (BinaryFindObject (&try_name)>=0) {
      suffix++;
      try_name = base & "_" & suffix;
    }
  else
    while (FindObject(&try_name)>=0) {
      suffix++;
      try_name = base & "_" & suffix;
    }
  
  return try_name;

}



