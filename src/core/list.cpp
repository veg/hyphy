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

#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <stdarg.h>

#include "list.h"
#include "hy_strings.h"
#include "hy_string_buffer.h"
#include "parser.h"

using namespace hy_global;


/*
==============================================================
Constructors
==============================================================
*/

// Does nothing
_List::_List () {
}

// Length constructor
_List::_List (unsigned long l):_SimpleList(l) {
}

// Stack Copy contructor
_List::_List (const _List& l, long from, long to) {
    if (from == 0 && to == -1) { // copy the whole thing
        BaseRef   br = (BaseRef)&l;
        _List::Duplicate (br);
    } else {
        Initialize           ();
        NormalizeCoordinates (from, to, l.lLength);

        for (long i = from; i <= to; i++) {
            (*this) << ((BaseRef*)l.list_data)[i];
        }
    }
}

// Construct a list of substrings from the original string separated by char
_List::_List (BaseRef ss, char sep) {
    _String* s = (_String*)ss;
    if (s->empty() == false) {
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
_List::_List (BaseRef br) : _SimpleList (0L) {
    ((BaseRef*)list_data)[0]= br->makeDynamic();
}


//Destructor
_List::~_List(void) {
    if (CanFreeMe()) {
        BaseRef * references =(BaseRef*)list_data;
    
        for (unsigned long i = 0UL; i<lLength; i++) {
            BaseRef an_object = references[i];
            if (an_object) {
                if (an_object->CanFreeMe()) {
                    DeleteObject(an_object);
                } else {
                    an_object->RemoveAReference();
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
BaseRef& _List::operator [] (long i) {
    // TODO 20170426: SLKP, maybe deprecate
    BaseRef ith_object = GetItem(i);
    if (ith_object)
        if (!ith_object->CanFreeMe()) {
            ith_object->RemoveAReference();
            ((BaseRef*)(list_data))[i]=ith_object->makeDynamic();
        }

    return ((BaseRef*)(list_data))[i];
}

// Element location functions (0,llength - 1)
BaseRef _List::operator () (const unsigned long i) {
    return ((BaseRef*)list_data)[i];
}


  // Element location functions (0,llength - 1)
BaseRef _List::GetItemRangeCheck(const unsigned long i) const {
  return i < lLength ? ((BaseRef*)list_data)[i] : nil;
}


// Assignment operator
const _List & _List::operator = (const _List& l) {

    Clear(true);
    lLength = l.lLength;
    RequestSpace(laLength);
    for (unsigned long i = 0UL; i<lLength; i++) {
        (((BaseRef*)(list_data))[i] = l.GetItem (i))-> AddAReference();
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
  
    BaseRef * res_data = (BaseRef*) res.list_data;
  
    if (list_data) {
      BaseRef * base_data = (BaseRef*) list_data;
      
      for (res.lLength = 0UL; res.lLength < lLength; res.lLength++) {
        res_data [res.lLength] = base_data[res.lLength];
        res_data [res.lLength] -> AddAReference();
      }
    }
  
    if (l.list_data) {
      BaseRef * l_data = (BaseRef*) l.list_data;
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
  
    BaseRef * res_data = (BaseRef*) res.list_data;
    if (list_data) {
      BaseRef * base_data = (BaseRef*) list_data;
      
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
    unsigned long incBy = ((MEMORYSTEP >> 2) > lLength)? MEMORYSTEP: (lLength << 2);
    laLength+=incBy;
    _EnsureCorrectStorageType();
    /*
    if (list_data) {
      list_data = (long*)MemReallocate((char*)list_data, laLength*sizeof(void*));
    } else {
      list_data = (long*)MemAllocate(laLength*sizeof(void*));
    }*/
  }
  
  ((BaseRef*)list_data)[lLength-1]=br;
  return *this;
}

_List& _List::operator << (BaseRef br) {
    br->AddAReference();
    return (*this) < br;
}

_List& _List::operator << (_List const& source) {
  for (unsigned long k=0UL; k<source.lLength; k++) {
      (*this) << ((BaseRef*)source.list_data)[k];
  }
  return *this;
}

_List& _List::operator < (_List const& source) {
  for (unsigned long k=0UL; k<source.lLength; k++) {
    (*this) < ((BaseRef*)source.list_data)[k];
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

long _List::AppendNewInstance (BaseRef br) {
    if (br) {
        (*this)<br;
    } else {
        HandleApplicationError(_String ("Passed a null reference to ") & __PRETTY_FUNCTION__);
    }
    return lLength;
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
        _String* stp = (_String*)(((BaseRef*)list_data)[middle]->toStr());
        hyComparisonType      cres = st->Compare (*stp);
        DeleteObject (stp);

        if (cres == kCompareLess) {
            top = middle==top?top-1:middle;
        } else if (cres == kCompareEqual) {
            return middle;
        } else {
            bottom = middle==bottom?bottom+1:middle;
        }

    }
    middle = top;
    _String* stp=(_String*)(((BaseRef*)list_data)[middle]->toStr());
    if (st->Equal(*stp)) {
        DeleteObject(stp);
        return middle;
    }
    DeleteObject(stp);
    return -middle-2;
}

long  _List::BinaryInsert (BaseRef s) {

    if (!lLength) {
        InsertElement (s,0,true);
        return 0;
    }

    long pos = -BinaryFindObject (s)-2;
    if (pos<0) {
        return -pos+2;
    }
    
    if (pos < lLength) {
        _String *s1 = (_String*)s->toStr(),
                *s2 = (_String*) ((*this)(pos))->toStr();
        if (*s2<*s1) {
            pos++;
        }
        DeleteObject(s1);
        DeleteObject(s2);
    }
    InsertElement (s,pos,true);
    return pos>=lLength?lLength-1:pos;

}

void    _List::bumpNInst (void) {
    for (unsigned long i = 0UL; i<lLength; i++) {
        ((BaseRef*)list_data)[i]->AddAReference();
    }
}

void  _List::Clear (bool completeClear) {
    if (CanFreeMe()) {
        for (unsigned long i = 0UL; i<lLength; i++) {
            DeleteObject (((BaseRef*)list_data)[i]);
        }
        _SimpleList::Clear(completeClear);

    } else {
        RemoveAReference();
    }
}

hyComparisonType  _List::Compare (long i, long j) const {
    _String             *si = (_String*)list_data[i],
                         *sj = (_String*)list_data[j];

    return  si->Compare(*sj);
}

hyComparisonType  _List::Compare (BaseObj const * i, long j) const {
    _String const       *sj = (_String const*)list_data[j],
                        *si = (_String const*)i;

    return  si->Compare(*sj);
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
            DeleteObject (((BaseRef*)list_data)[index]);
        }
        lLength--;
        if (lLength-index)
            for (unsigned long i = index; i < lLength; i++) {
                list_data[i] = list_data[i+1];
            }
        //memcpy ((hyPointer)list_data+sizeof(BaseRef)*(index),(hyPointer)list_data+sizeof(BaseRef)*(index+1),sizeof(BaseRef)*(lLength-index));
    }
    
    _UpdateStorageType();

}

// Delete item at index (>=0)
void  _List::DeleteTail (long from, bool delete_object) {
    if (from>=0 && from<lLength) {
        if (delete_object)
            for (long k = lLength-1; k >= from; k--) {
                     DeleteObject (((BaseRef*)list_data)[k]);
            }
        lLength = from;
        //memcpy ((hyPointer)list_data+sizeof(BaseRef)*(index),(hyPointer)list_data+sizeof(BaseRef)*(index+1),sizeof(BaseRef)*(lLength-index));
    }
    
    _UpdateStorageType();

}

// Delete item at index (>=0)
void  _List::DeleteList (const _SimpleList& toDelete)
{
    if (toDelete.lLength) {
        long k = 0;
        for (unsigned long i = 0; i<lLength; i++) {
            if (k<toDelete.lLength && i==toDelete.list_data[k]) {
                DeleteObject (((BaseRef*)list_data)[i]);
                //if (k<toDelete.lLength)
                k++;
            } else {
                ((BaseRef*)list_data)[i-k] = ((BaseRef*)list_data)[i];
            }
        }
        lLength -= k;
        _UpdateStorageType();

    }
}

void    _List::Duplicate (BaseRefConst theRef) {
    _List::Clear();
    _SimpleList::Duplicate (theRef);
    for (unsigned long i = 0UL; i<lLength; i++) {
        if (((BaseRef*)list_data)[i]) {
            (((BaseRef*)list_data)[i])->AddAReference();
        }
    }
}

bool _List::Equal(_List const & l2) const
{
    if (lLength!=l2.lLength) {
        return false;
    }

    for (unsigned long i=0; i<lLength; i++) {
        if (!((_String*)list_data[i])->Equal (*(_String*)l2.list_data[i])) {
            return false;
        }
    }

    return true;
}

long  _List::FindObject (BaseRefConst s, long startat) const {
    _String const * st = (_String const*)s;
    for (unsigned long i = startat; i<lLength; i++) {
      
        _String * sp = (_String*)(((BaseRef*)list_data)[i]->toStr());

        if (*st == *sp) {
            DeleteObject(sp);
            return i;
        }
        DeleteObject(sp);
    }
    return kNotFound;
}


// Append & store operator
void _List::InsertElement (BaseRef br, long insertAt, bool store, bool pointer) {
    _SimpleList::InsertElement (br, insertAt, store, pointer);
}

void    _List::Intersect (_List& l1, _List& l2, _SimpleList* idx, _SimpleList* idx2) {
    // compute the union of two sorted lists
    // each repeat appears exactly once
    if (nonempty()) {
        Clear();
    }

    long  c1 = 0L,
          c2 = 0L;

    while (c1<l1.lLength && c2<l2.lLength) {
        while (c1<l1.lLength && ((_String*)l1(c1))->Compare(*(_String*)l2(c2)) == kCompareLess) {
            c1++;
        }
      
        if (c1==l1.lLength) {
            break;
        }

        while (c1<l1.lLength && c2<l2.lLength &&((_String*)l1(c1))->Equal(*(_String*)l2(c2))) {
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
        while (c2<l2.lLength && ((_String*)l2(c2))->Compare(*(_String*)l1(c1)) == kCompareLess) {
            c2++;
        }
    }
}

_String*  _List::Join (_String const & spacer, long startAt, long endAt) const {
    _StringBuffer *joined = new _StringBuffer (256L);
    if (endAt < 0) { 
        endAt = lLength; 
    } else if (endAt > lLength) {
        endAt = lLength;
    }

    for (unsigned long k = MAX(0L,startAt); k < endAt; k++) {
        if (k > startAt) {
            (*joined) << spacer;
        }
        joined->AppendNewInstance((_String*) ((BaseRef*)list_data)[k]->toStr());
    }

    return joined;
}

void _List::Place (BaseRef br) {
//  InsertElement (br, -1, false);
    lLength++;
    if (lLength>laLength) {
        laLength+=MEMORYSTEP;
        _EnsureCorrectStorageType();
    }
    ((BaseRef*)list_data)[lLength-1]=br;
}

//TODO: makeDynamic should be MakeDynamic to follow convention.
BaseRef _List::makeDynamic(void) const {
    _List * Res = new _List;
    Res->Duplicate (this);
    return Res;
}

void  _List::Replace (long index, BaseRef newObj, bool dup) {
    if (index>=0 && index<lLength) {
        BaseRef payload = dup?newObj->makeDynamic():newObj;
        // important to do this BEFORE calling DeleteObject in case newObj == existing object
        DeleteObject (((BaseRef*)list_data)[index]);
        ((BaseRef*)list_data)[index] = payload;
    }
}

// Char* conversion
//TODO: toFileStr should be ToFileStr to follow convention.
void _List::toFileStr(FILE* dest, unsigned long) {
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
    _StringBuffer * s = new _StringBuffer((unsigned long)20*(lLength+1));

    (*s)<<'{';
  
    if (lLength) {
      s->AppendNewInstance ((_String*)GetItem(0)->toStr());
      for (unsigned long i = 1UL; i<lLength; i++) {
          (*s)<< ", ";
          s->AppendNewInstance ((_String*)GetItem(i)->toStr());
      }
    }
  
    (*s)<<'}';
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



