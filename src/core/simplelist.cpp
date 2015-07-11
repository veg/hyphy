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
#include "list.h"
#include "simplelist.h"
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
_SimpleList::_SimpleList ()
{
    Initialize(false);
}

//Data constructor (1 member list)
_SimpleList::_SimpleList (long br)
{
    lLength = 1;
    laLength = MEMORYSTEP;
    lData = (long*)MemAllocate (laLength * sizeof(Ptr));
    ((long*)lData)[0]= br;
}

//Length constructor and populator
_SimpleList::_SimpleList (long l, long start, long step)
{
    Initialize (false);
    Populate   (l,start,step);
}

//Length constructor
_SimpleList::_SimpleList (unsigned long l)
{
    lLength = 0;
    laLength = (l/MEMORYSTEP + 1)*MEMORYSTEP;
    lData = (long*)MemAllocate (laLength * sizeof(Ptr));
    memset (lData,0,laLength * sizeof(Ptr));
}

//Stack copy contructor
_SimpleList::_SimpleList (_SimpleList& l, long from, long to)
{
    if (from == 0 && to == -1) { // copy the whole thing
        Duplicate (&l);
    } else {
        Initialize           ();
        NormalizeCoordinates (from, to, l.lLength);
        RequestSpace(to-from);
        long upto = to-from ; 
        for (long k = 0; k < upto; k++) {
            lData[k] = l.lData[from+k];
        }
        /*
        for (long i = from; i < to; i++) {
            (*this) << l.lData[i];
        }
        */
    }
}

// Data constructor (variable number of long constants)
_SimpleList::_SimpleList (const long value1, const unsigned long number, ...)
{
    Initialize (true);
    va_list vl;
    
    (*this) << value1;
    
    va_start(vl,number);
    for (unsigned long arg_id =0;arg_id<number;arg_id++) {
        const long this_arg =va_arg(vl,long);
        (*this) << this_arg;
    }
    va_end(vl);
}

//Destructor
_SimpleList::~_SimpleList(void)
{
    if (nInstances<=1) {
        if (lData) {
            free (lData);
        }
    } else {
        nInstances--;
    }


}

/*
==============================================================
Operator Overloads
==============================================================
*/

//Element location functions (0,llength - 1)
long& _SimpleList::operator [] (const long i)
{
    if (lLength == 0) {
        return lData[0];
    }

    const unsigned long in = (const unsigned long)i;
    if (in>lLength-1) {
        return lData[lLength-1];
    }

    return lData[in];
}
//Element location functions (0,llength - 1)
long _SimpleList::operator () (const unsigned long i)
{
    //if (lLength == 0) return 0;
    //Is there a reason why this is commented out?
    //if (i>=lLength) i = lLength-1;
    if(i < lLength) {
        return lData[i];
    }
    warnError("List index out of range");
    return -1;
}

//Assignment operator
_SimpleList _SimpleList::operator = (_SimpleList l)
{
    Clear();
    lLength  = l.lLength;
    laLength = l.laLength;
    if (laLength) {
        checkPointer (lData = (long*)MemAllocate (laLength*sizeof (Ptr)));
        if (lLength) {
            memcpy (lData,l.lData,lLength*sizeof (Ptr));
        }
    }

    return *this;
}

//Append operator
_SimpleList _SimpleList::operator & (_SimpleList l)
{
    _SimpleList res (l.lLength + lLength);
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

    return res;
}

void _SimpleList::operator << (long br)
{
    InsertElement ((BaseRef)br, -1, false, false);
}

bool _SimpleList::operator >> (long br)
{
    if (Find(br) == -1) {
        InsertElement ((BaseRef)br, -1, false, false);
        return true;
    }
    return false;
}

void _SimpleList::operator << (_SimpleList& source)
{
    for (unsigned long k=0; k<source.lLength; k++) {
        (*this) << source.lData[k];
    }
}

/*
==============================================================
Methods
==============================================================
*/


//Element location functions (0,llength - 1), negative values return 
// elements from the end of the list

long _SimpleList::GetElement (const long index)
{
    if (index >= 0) {
        if ((const unsigned long) index < lLength) {
            return lData [index];
        }
    } 
    if ((const unsigned long) (-index) <= lLength) {
        return lData[lLength + index];
    }
    warnError(_String("List index '") & (long)((const unsigned long) (-index)) & "' out of range in _SimpleList::GetElement on list of length " & long (lLength));    
    return 0;
}


long  _SimpleList::BinaryFind (long s, long startAt)
{
    long top    =   lLength-1,
         bottom  =   startAt,
         middle;

    if (top == -1) {
        return -2;
    }

    while (top>bottom) {
        middle = (top+bottom)/2;
        if (s<((long*)lData)[middle]) {
            top = middle==top?top-1:middle;
        } else if (s>((long*)lData)[middle]) {
            bottom = middle==bottom?bottom+1:middle;
        } else {
            return middle;
        }
    }

    middle     = top;
    long comp  = ((long*)lData)[middle]-s;
    if (!comp) {
        return middle;
    }

    return comp<0?-middle-3:-middle-2;
}

long  _SimpleList::BinaryInsert (long n)
{
    if (!lLength) {
        (*this) << n;
        return 0;
    }

    long pos = -BinaryFind (n)-2;

    if (pos<0) {
        return -pos+2;
    }

    if (lData[pos]<n) {
        pos++;
    }

    InsertElement ((BaseRef)n,pos,false,false);

    return pos>=lLength?lLength-1:pos;
}

void _SimpleList::ClearFormulasInList(void)
{
    for (unsigned long k = 0; k < lLength; k++)
        if (lData[k]) {
            delete (_Formula*)lData[k];
        }
}

long _SimpleList::Sum (void) {
    long sum = 0;
    for (unsigned long k = 0; k < lLength; k++) {
       sum += lData[k];
    }
    return sum;
}

long  _SimpleList::Compare (long i, long j)
{
    long    v1 = ((long*)lData)[i],
            v2 = ((long*)lData)[j];


    if (v1<v2) {
        return -1;
    } else if (v1==v2) {
        return 0;
    } else {
        return 1;
    }


    //return ((long*)lData)[i]-((long*)lData)[j];
}

long  _SimpleList::Compare (BaseRef i, long j)
{
    long    v1 = (long)i,
            v2 = ((long*)lData)[j];


    if (v1<v2) {
        return -1;
    } else if (v1==v2) {
        return 0;
    } else {
        return 1;
    }

    //return (long)i-((long*)lData)[j];
}

// Compute the number of shared of two sorted lists
long    _SimpleList::CountCommonElements (_SimpleList& l1, bool yesNo)
{
    long  c1    = 0,
          c2    = 0,
          res   = 0;


    while (c1<l1.lLength && c2<lLength) {
        while (l1.lData[c1]<lData[c2]) {
            c1++;
            if (c1==l1.lLength) {
                break;
            }
        }
        if (c1==l1.lLength) {
            break;
        }

        while (l1.lData[c1]==lData[c2]) {
            c2++;
            if (yesNo) {
                return 1;
            } else {
                res++;
            }
            if (c1==l1.lLength || c2==lLength) {
                break;
            }
        }
        if (c1==l1.lLength || c2==lLength) {
            break;
        }
        while (lData[c2]<l1.lData[c1]) {
            c2++;
            if (c2==lLength) {
                break;
            }
        }
    }

    return res;
}

//List length
unsigned long _SimpleList::countitems(void)
{
    return lLength;
}


_SimpleList*  _SimpleList::CountingSort (long upperBound, _SimpleList* ordering)
{
    if (ordering) {
        ordering->Clear();
    }

    if (lLength) {
        if (upperBound < 0) {
            upperBound = Max()+1;
        }

        _SimpleList buffer      (upperBound, 0, 0),
                    * result    =  new _SimpleList (lLength);
        for (long pass1 = 0; pass1 < lLength; pass1 ++) {
            buffer.lData[lData[pass1]] ++;
        }
        for (long pass2 = 1; pass2 < upperBound; pass2 ++) {
            buffer.lData[pass2] += buffer.lData[pass2-1];
        }
        if (ordering) {
            ordering->Populate (lLength, 0, 0);
            for (long pass3 = lLength-1; pass3 >=0; pass3--) {
                result->lData[--buffer.lData[lData[pass3]]] = lData[pass3];
                ordering->lData[buffer.lData[lData[pass3]]] = pass3;
            }
        } else
            for (long pass3 = lLength-1; pass3 >=0; pass3--) {
                result->lData[--buffer.lData[lData[pass3]]] = lData[pass3];
            }
        result->lLength = lLength;

        return result;
    }
    return new _SimpleList;
}

void  _SimpleList::Clear (bool completeClear)
{
    if (nInstances<=1) {
        lLength = 0;
        if (completeClear) {
            laLength = 0;
            if (lData) {
                free (lData);
            }
            lData = nil;
        }
    } else {
        nInstances--;
    }
}

void _SimpleList::DebugVarList(void)
{
    printf ("\nVariable list dump:\n");
    for  (unsigned long e = 0; e < lLength; e++) {
        if (lData[e] >= 0) {
            _Variable * theV = LocateVar (lData[e]);
            if (theV) {
                printf ("[%s]\n", theV->GetName()->getStr());
                continue;
            }
        }
        printf ("[Empty]\n");
    }
}

//Delete item at index (>=0)
void  _SimpleList::Delete (long index, bool compact)
{
    if (index>=0 && index<lLength) {
        lLength--;
        if (lLength-index) {
            memmove ((Ptr)lData+sizeof(BaseRef)*(index),(Ptr)lData+sizeof(BaseRef)*(index+1),sizeof(BaseRef)*(lLength-index));
        }
    }
    if (compact && laLength-lLength>MEMORYSTEP) {
        laLength -= ((laLength-lLength)/MEMORYSTEP)*MEMORYSTEP;
        if (laLength) {
            lData = (long*)MemReallocate ((char*)lData, laLength*sizeof(Ptr));
        } else {
            free (lData);
            lData = nil;
        }
    }

}

//Delete duplicates from a sorted list
void  _SimpleList::DeleteDuplicates (void)
{
    if (lLength>1) {
        _SimpleList noDups;

        long    lastValue = lData[0]+1;
        for (unsigned long k=0; k<lLength; k++) {
            long thisValue = lData[k];
            if (thisValue!=lastValue) {
                noDups << thisValue;
                lastValue = thisValue;
            }
        }

        if (noDups.lLength!=lLength) {
            Duplicate (&noDups);
        }
    }

}

//Delete items from a sorted list
void  _SimpleList::DeleteList (const _SimpleList& toDelete)
{
    if (toDelete.lLength) {
        unsigned long k = 0;
        for (unsigned long i = 0; i<lLength; i++) {
            if (k<toDelete.lLength && i==toDelete.lData[k])
                //if (k<toDelete.lLength)
            {
                k++;
            } else {
                lData[i-k] = lData[i];
            }
        }
        lLength -= toDelete.lLength;
    }
    if (laLength-lLength>MEMORYSTEP) {
        laLength -= ((laLength-lLength)/MEMORYSTEP)*MEMORYSTEP;
        if (laLength) {
            lData = (long*)MemReallocate ((char*)lData, laLength*sizeof(Ptr));
        } else {
            free (lData);
            lData = nil;
        }
    }
}

//Shift the range from start to end
void  _SimpleList::Displace (long start, long end, long delta)
{
    if (start<0) {
        start = 0;
    } else if (start>=lLength) {
        start = lLength-1;
    }

    if (end<0) {
        end = lLength-1;
    } else if (end>=lLength) {
        end = lLength-1;
    }

    if ((end-start>=0) && delta && (end-start<lLength-1))
    {
        if (delta>0 && lLength-end <= delta) { // shift up
            delta = lLength-end-1;
        } else if (start-delta<0){
            delta = start;
        }

        if (delta) {
            long i,j,delta2 = end-start+1;
            _SimpleList swapList ((unsigned long)(end-start+1));

            for (i=start; i<=end; i++) {
                swapList << lData[i];
            }

            if (delta>0) {

                for (i=end+1; i<=end+delta; i++) {
                    lData[i-delta2] = lData[i];
                }

            } else {
                for (i=start-1; i>=start+delta; i--) {
                    lData[i+delta2] = lData[i];
                }
            }

            for (i=start+delta, j=0; i<=end+delta; i++,j++) {
                lData[i] = swapList.lData[j];
            }
        }

    }
}

void _SimpleList::Duplicate(BaseRef theRef)
{
    _SimpleList* l  = (_SimpleList*)theRef;
    lLength         = l->lLength;
    laLength        = l->laLength;
    lData           = l->lData;
    if (lData) {
        checkPointer (lData = (long*)MemAllocate (laLength*sizeof (Ptr)));
        memcpy ((char*)lData, (char*)l->lData, lLength*sizeof (Ptr));
    }
}

//Element location functions (0,llength - 1)
//Negative indices return offsets from the end of the list
long _SimpleList::Element(long index)
{
    if (index >= 0 && index < lLength) {
        return lData[index];
    }

    else if(-index <= lLength) {
        return lData[(long)lLength+index];
    }

    return 0;
}

bool _SimpleList::Equal(_SimpleList& l2)
{
    if (lLength!=l2.lLength) {
        return false;
    }

    for (long i=0; i<lLength; i++)
        if (lData[i] != l2.lData[i]) {
            return false;
        }

    return true;
}

long  _SimpleList::Find (long s, long startAt)
{
    for (unsigned long i = startAt; i<lLength; i++) {
        if ( ((long*)(lData))[i] == s ) {
            return i;
        }
    }
    return -1;
}

long  _SimpleList::FindStepping (long s, long step, long startAt)
{
    for (unsigned long i = startAt; i<lLength; i+=step)
        if (lData[i] == s) {
            return i;
        }

    return -1;
}

void  _SimpleList::FilterRange (long lb, long ub)
{
    if (ub <= lb) {
        Clear();
    } else {
        _SimpleList toDelete;
        for (long k = 0; k < lLength; k++)
            if (lData[k] <= lb || lData[k] >= ub) {
                toDelete << k;
            }
        DeleteList (toDelete);
    }
}

void _SimpleList::Flip ()
{
    for (long k=0, l=lLength-1; k<l; k++,l--) {
        void * pt = ((void**)lData)[k];
        ((void**)lData)[k] = ((void**)lData)[l];
        ((void**)lData)[l] = pt;
    }
}

void _SimpleList::Initialize(bool doMemAlloc)
{
    BaseObj::Initialize();
    lLength = 0;
    if (doMemAlloc) {
        laLength = MEMORYSTEP;
        lData = (long*)MemAllocate (laLength * sizeof(Ptr));
    } else {
        laLength = 0;
        lData    = nil;
    }
}

//Append & store operator
void _SimpleList::InsertElement (BaseRef br, long insertAt, bool store, bool pointer)
{
    lLength++;
    if (lLength>laLength) {
        unsigned long incBy = (MEMORYSTEP*5 > lLength)? MEMORYSTEP: lLength/5;

        laLength+=incBy;

        //memAlloc += sizeof(Ptr)*incBy;

        if (lData) {
            lData = (long*)MemReallocate((char*)lData, laLength*sizeof(void*));
        } else {
            lData = (long*)MemAllocate(laLength*sizeof(void*));
        }

        if (!lData) {
            checkPointer (lData);
        }
    }
    if (insertAt==-1) {
        if (store) {
            ((BaseRef*)lData)[lLength-1]=br->makeDynamic();
        } else {
            ((BaseRef*)lData)[lLength-1]=br;
            if (pointer) {
                br->nInstances++;
            }
        }
    } else {
        //insertAt = insertAt>=lLength?lLength:insertAt;
        insertAt = insertAt>=lLength?lLength-1:insertAt;
        long     moveThisMany = (laLength-insertAt-1);
        if (moveThisMany < 32)
            for (long k=insertAt+moveThisMany; k> insertAt ; k--) {
                lData[k] = lData[k-1];
            }
        else {
            memmove (((char**)lData)+(insertAt+1), ((char**)lData)+insertAt, moveThisMany*sizeof(void*));
        }

        if (store) {
            ((BaseRef*)lData)[insertAt]=br->makeDynamic();
        } else {
            ((BaseRef*)lData)[insertAt]=br;
            if (pointer) {
                br->nInstances++;
            }
        }
    }


}

//Convert a list into a partition style string
BaseRef _SimpleList::ListToPartitionString ()
{
    _String *result = new _String ((unsigned long)64,true),
    conv;

    for (long k=0; k<lLength; k++) {
        long m;
        for (m=k+1; m<lLength; m++)
            if (lData[m]-lData[m-1]!=1) {
                break;
            }
        if (m>k+2) {
            conv = lData[k];
            (*result) << & conv;
            (*result) << '-';
            conv = lData[m-1];
            (*result) << & conv;
            if (m<lLength) {
                (*result) << ',';
            }
            k = m-1;
        } else {
            conv = lData[k];
            (*result) << &conv;
            if (k<lLength-1) {
                (*result) << ',';
            }
        }
    }
    (*result).Finalize();
    return result;
}

BaseRef _SimpleList::makeDynamic(void)
{
    _SimpleList * Res = new _SimpleList;
    checkPointer(Res);
    memcpy ((char*)Res, (char*)this, sizeof (_SimpleList));
    Res->nInstances = 1;
    Res->lData = nil;
    Res->Duplicate (this);
    return Res;
}

long _SimpleList::Max(void)
{
    long res = LONG_MIN;
    for  (long e = 0; e < lLength; e++)
        if (lData[e] > res) {
            res = lData[e];
        }
    return res;
}

long _SimpleList::Min(void)
{
    long res = LONG_MAX;
    for  (long e = 0; e < lLength; e++)
        if (lData[e] < res) {
            res = lData[e];
        }
    return res;
}

//Merge 2 lists (sorted)

void _SimpleList::Merge(_SimpleList& l1, _SimpleList& l2, _SimpleList* mergeResults1, _SimpleList* mergeResults2) {
  
  this->Clear();
  if (mergeResults1) {
    mergeResults1->Clear();
  }
  if (mergeResults2) {
    mergeResults2->Clear();
  }

  enum    machine_states {INIT, ADVANCE1, ADVANCE2, FLUSH1, FLUSH2} advancing = INIT;

  unsigned long   pos1 = 0,
                  pos2 = 0,
                  nt1 = l1.countitems(),
                  nt2 = l2.countitems();


  bool  keep_going = true;

  while (keep_going) {
    switch (advancing) {
        
      case ADVANCE1: { // advancing in the 1st list
        pos1++;
        if (pos1==nt1) {
          advancing = FLUSH2;
          continue;
        }
        
        long cmp = l1.lData[pos1] - l2.lData [pos2];
        //if (l1lData[pos1] <= l2->lData[pos2]) {
        if ( cmp <= 0L) {
          if (mergeResults1) {
            //printf("MERGE line number %d in file %s\n", __LINE__, __FILE__);
            (*mergeResults1)<<this->lLength;
          }
          
          (*this) << (l1.lData[pos1]);
          
          //        if (mergeResults2 && l1->lData[pos1] < l2->lData[pos2]) {
          if (cmp < 0L) {
            if (mergeResults2 && mergeResults2->countitems () == 0UL && pos1 >= pos2) {
              for (unsigned long i=0UL; i<pos2; i++) {
                (*mergeResults2) << i;
              }
            }
            continue;
          }
        }
        
        if (cmp > 0L) {
          //if (l1->lData[pos1] > l2->lData[pos2]) {
          advancing = ADVANCE2;
          if (mergeResults1 && mergeResults1->countitems () == 0UL ) {
            for (unsigned long i=0UL; i<pos1; i++) {
              //printf("MERGE line number %d in file %s\n", __LINE__, __FILE__);
              (*mergeResults1) << i;
            }
          }
          
          if (mergeResults2) {
            (*mergeResults2)<<this->lLength;
          }
          
          //printf ("this->append (l2->AtIndex(pos2))\n");
          (*this)  << l2.lData[pos2];
          continue;
        }
        break;
      }
        
      case ADVANCE2: { // advancing in the 2nd list
        pos2++;
        if (pos2==nt2) {
          advancing = FLUSH1;
          continue;
        }
        
        long cmp = l2.lData[pos2] - l1.lData[pos1];
        
        //if (l2->lData[pos2] <= l1->lData[pos1]) {
        if (cmp <= 0L) {
          if (mergeResults2) {
            (*mergeResults2)<<this->lLength;
          }
          
          (*this) << l2.lData[pos2];
          
          //if (l2->lData[pos2] < l1->lData[pos1] && mergeResults1 && !doMerge1 && pos2>=pos1) {
          if (cmp < 0L ) {
            if(mergeResults1 && mergeResults1->countitems() == 0 && pos2>=pos1) {
              for (unsigned long i=0UL; i<pos1; i++) {
                //printf("MERGE line number %d in file %s\n", __LINE__, __FILE__);
                (*mergeResults1)<<i;
              }
            }
            continue;
          }
        }
        
        //if (l2->lData[pos2] > l1->lData[pos1]) {
        if (cmp > 0L) {
          advancing = ADVANCE1;
          if (mergeResults2 && mergeResults2->countitems() == 0) {
            for (unsigned long i=0UL; i<pos2; i++) {
              (*mergeResults2)<<i;
            }
          }
          if (mergeResults1) {
            //printf("MERGE line number %d in file %s\n", __LINE__, __FILE__);
            (*mergeResults1)<<this->lLength;
          }
          (*this) << l1.lData[pos1];
          continue;
        }
        break;
      }
      case FLUSH2: { // flush out the 2nd list
        if (mergeResults1 && pos2<nt2 ) {
          for (unsigned long i=pos1; i<nt1; i++) {
            //printf("MERGE line number %d in file %s\n", __LINE__, __FILE__);
            (*mergeResults1)<<i;
          }
        }
        if (mergeResults2)
          while (pos2<nt2) {
            (*mergeResults2)<<this->lLength;
            (*this) << (l2.lData[pos2++]);
          }
        else
          while (pos2<nt2) {
            (*this) << (l2.lData[pos2++]);
          }
        keep_going = false;
        break;
      }
      case FLUSH1: { // flush out the 1st list
        if (mergeResults2 && pos1<nt1) {
          for (unsigned long i=pos2; i<nt2; i++) {
            (*mergeResults2)<<i;
          }
        }
        if (mergeResults1)
          while (pos1<nt1) {
            //printf("MERGE line number %d in file %s\n", __LINE__, __FILE__);
            (*mergeResults1)<<this->lLength;
            (*this) << (l1.lData[pos1++]);
          }
        else
          while (pos1<nt1) {
            (*this) << (l1.lData[pos1++]);
          }
        keep_going = false;
        break;
      }
      case INIT: { // just starting
        if (!nt1) { // first list is empty!
          advancing = FLUSH2;
          continue;
        }
        if (!nt2) { // second list is empty!
          advancing = FLUSH1;
          continue;
        }
        
        //if (l1->lData[pos1] <= l2->lData[pos2]) { // begin with the first list
        
        long cmp = l1.lData[pos1] - l2.lData[pos2];
        
        if (cmp <= 0L) { // begin with the first list
          if (mergeResults1) {
            (*mergeResults1)<<this->lLength;
          }
          (*this) << (l1.lData[0]);
          advancing = ADVANCE1;
          if (cmp != 0L) {
            continue;
          }
        } else {
          if (mergeResults2) {
            (*mergeResults2)<<this->lLength;
          }
          (*this) << (l2.lData[0]);
          advancing = ADVANCE2;
          continue;
        }
        break;
      }
    } // end SWITCH
    
    if (keep_going) {
      if (advancing == ADVANCE1) { // moving up in the second term
        pos1++;
        if (pos1==nt1) {
          pos2++;
          if (mergeResults2) {
            (*mergeResults2)<<this->lLength-1UL;
          }
          advancing = FLUSH2;
          continue;
        } else {
          advancing = ADVANCE2;
          if (mergeResults2) {
            (*mergeResults2)<<this->lLength-1UL;
          }
        }
      } else {
        pos2++;
        if (pos2==nt2) {
          pos1++;
          if (mergeResults1) {
            (*mergeResults1)<<this->lLength-1UL;
          }
          advancing = FLUSH1;
          continue;
        } else {
          advancing = ADVANCE1;
          if (mergeResults1) {
            (*mergeResults1)<<this->lLength-1UL;
          }
        }
      }
    }
  }
}

// Together with the next function
// Implements algorithm NEXKSB from p.27 of http://www.math.upenn.edu/~wilf/website/CombinatorialAlgorithms.pdf
bool  _SimpleList::NChooseKInit (_SimpleList& state, _SimpleList& store, unsigned long stride, bool algorithm)
{
    if (stride <= lLength && lLength) {
        state.Clear();
        state.RequestSpace (stride+3);
        state << stride;
        store.Clear();
        store.RequestSpace (stride);
        return true;
    }
    return false;
}

// Implements algorithm NEXKSB from p.27 of http://www.math.upenn.edu/~wilf/website/CombinatorialAlgorithms.pdf
bool  _SimpleList::NChooseK (_SimpleList& state, _SimpleList& store)
{
    if (state.lLength == 1) { // first pass
        state << 0;              // m
        state << state.lData[0]; // h
        state.lLength = state.lData[0]+3;
        store.lLength = state.lData[0];
        if (store.lLength == 0) {
            return false;
        }
    } else {
        if (state.lData[1] < lLength - state.lData[2]) {
            state.lData[2] = 0;
        }
        state.lData[2] ++;
        state.lData[1] = state.lData[3 + state.lData[0] - state.lData[2]] + 1;
    }
    for (long j=1; j <= state.lData[2]; j++) {
        long  anIndex   = j+state.lData[0]-state.lData[2],
              anIndex2    = state.lData[1]+j;
        state.lData[anIndex+2]   = anIndex2-1;
        store.lData[anIndex-1]   = lData[anIndex2-1];
    }
    return state.lData[3] < lLength-state.lData[0];
}

//Coordinate normalizer
void _SimpleList::NormalizeCoordinates (long& from, long& to, const unsigned long refLength)
{
    if (to < 0) {
        to = refLength+to;
    } else {
        to = to < refLength-1 ? to : refLength - 1;
    }
    if (from < 0) {
        from = refLength+from;
    }
}

void _SimpleList::Offset (long shift)
{
    if (lData) {
        for (long k=0; k<lLength; k++) {
            lData[k] += shift;
        }
    }
}

_SimpleList* _SimpleList::Subset (unsigned long size, bool replacement)
{
    _SimpleList* result = new _SimpleList;
    if (size > 0) {
        size = MIN(size, lLength);
        if (replacement) {
            for (long k = 0; k < size; k++) {
                (*result) << lData[genrand_int32()%lLength];
            }
        } else {
            (*result) << (*this);
            for (long k = 0; k < size; k++) {
                long idx = lData[genrand_int32()%(lLength-k)];
                long t = result->lData[k];
                result->lData[k] = result->lData[idx];
                result->lData[idx] = t;
            }
            result->lLength = size;
            result->TrimMemory();
        }
    }
    return result;
}

// Create a permutation of the list's elements
void  _SimpleList::Permute (long blockLength)
{
    unsigned long blockCount = lLength/blockLength;

    if (blockLength>1) {
        /*_SimpleList     result ((unsigned long)(blockCount*blockLength));
        while (blockCount)
        {
            unsigned long sample = (unsigned long)(genrand_real2()*blockCount);
            sample *= blockLength;
            for (long j = 0; j<blockLength; j++)
            {
                result<<lData[sample];
                Delete(sample);
            }
            blockCount --;
        }
        Duplicate(&result); */

        for (unsigned long k=0; k<blockCount-1; k=k+1) {
            unsigned long k2 = genrand_real2()*(blockCount-k);
            if (k2) {
                k2 += k;
                k2 *= blockLength;

                for (long j = 0; j<blockLength; j++) {
                    long t = lData[k2+j];
                    lData[k2+j] = lData[k*blockLength+j];
                    lData[k*blockLength+j] = t;
                }
            }
        }

    } else {
        for (unsigned long k=0; k<blockCount-1; k=k+1) {
            unsigned long k2 = genrand_real2()*(blockCount-k);
            if (k2) {
                k2+=k;
                long t = lData[k2];
                lData[k2] = lData[k];
                lData[k] = t;
            }
        }
        /*{
            while (blockCount)
            {
                unsigned long sample = genrand_real2()*blockCount;
                result<<lData[sample];
                Delete(sample);
                blockCount --;
            }
        }*/
    }
}

// Create a permutation of the list's elements with possible repetitions
void  _SimpleList::PermuteWithReplacement (long blockLength)
{
    unsigned long blockCount = lLength/blockLength;
    _SimpleList   result ((unsigned long)(blockCount*blockLength));
    if (blockLength>1)
        for (long i = 0; i<blockCount; i++) {
            unsigned long sample = (unsigned long)(genrand_real2()*blockCount);
            sample *= blockLength;
            for (long j = 0; j<blockLength; j++,sample++) {
                result<<lData[sample];
            }
        } else {
            for (long i = 0; i<blockCount; i++) {
                unsigned long sample = genrand_real2()*blockCount;
                result<<lData[sample];
            }
        }

    Clear();
    Duplicate(&result);

}

long _SimpleList::Pop (void)
{
    if (lLength > 0) {
        lLength --;
        return lData[lLength];
    }

    return 0;
}

//Length constructor and populator
void _SimpleList::Populate (long l, long start, long step)
{
    RequestSpace (l);
    for (long k = 0; k < l; k++, start+=step) {
        lData[k] = start;
    }

    lLength = l;
}

void _SimpleList::RecursiveIndexSort (long from, long to, _SimpleList* index)
{
    long middle = (from+to)/2, middleV = lData[middle],
         bottommove = 1, topmove = 1, temp,i, imiddleV = (*index)(middle);
    long *idata = (*index).quickArrayAccess();
    if (middle)
        while ((middle-bottommove>=from)&&(Compare(middle-bottommove,middle)>=0)) {
            bottommove++;
        }

    if (from<to)
        while ((middle+topmove<=to)&&(Compare(middle+topmove,middle)<=0)) {
            topmove++;
        }

    // now shuffle
    for (i=from; i<middle-bottommove; i++) {
        if (Compare(i,middle)>=0) {
            temp = lData[middle-bottommove];
            lData[middle-bottommove] = lData[i];
            lData[i]=temp;
            temp = idata[middle-bottommove];
            idata[middle-bottommove] = idata[i];
            idata[i]=temp;
            bottommove++;
            while ((middle-bottommove>=from)&&(Compare(middle-bottommove,middle)>=0)) {
                bottommove++;
            }
        }
    }

    for (i=middle+topmove+1; i<=to; i++) {
        if (Compare(i,middle)<=0) {
            temp = lData[middle+topmove];
            lData[middle+topmove] = lData[i];
            lData[i]=temp;
            temp = idata[middle+topmove];
            idata[middle+topmove] = idata[i];
            idata[i]=temp;
            topmove++;
            while ((middle+topmove<=to)&&(Compare(middle+topmove,middle)<=0)) {
                topmove++;
            }
        }
    }

    if (topmove==bottommove) {
        for (i=1; i<bottommove; i++) {
            temp = lData[middle+i];
            lData[middle+i] = lData[middle-i];
            lData[middle-i]=temp;
            temp = idata[middle+i];
            idata[middle+i] = idata[middle-i];
            idata[middle-i]=temp;
        }
    } else if (topmove>bottommove) {
        long shift = topmove-bottommove;
        for (i=1; i<bottommove; i++) {
            temp = lData[middle+i+shift];
            lData[middle+i+shift] = lData[middle-i];
            lData[middle-i]=temp;
            temp = idata[middle+i+shift];
            idata[middle+i+shift] = idata[middle-i];
            idata[middle-i]=temp;
        }
        for (i=0; i<shift; i++) {
            lData[middle+i]=lData[middle+i+1];
            idata[middle+i]=idata[middle+i+1];
        }
        middle+=shift;
        lData[middle]=middleV;
        idata[middle]=imiddleV;
    } else {
        long shift = bottommove-topmove;
        for (i=1; i<topmove; i++) {
            temp = lData[middle-i-shift];
            lData[middle-i-shift] = lData[middle+i];
            lData[middle+i]=temp;
            temp = idata[middle-i-shift];
            idata[middle-i-shift] = idata[middle+i];
            idata[middle+i]=temp;
        }
        for (i=0; i<shift; i++) {
            lData[middle-i]=lData[middle-i-1];
            idata[middle-i]=idata[middle-i-1];
        }
        middle-=shift;
        lData[middle]=middleV;
        idata[middle]=imiddleV;
    }
    if (to>middle+1) {
        RecursiveIndexSort (middle+1,to, index);
    }
    if (from<middle-1) {
        RecursiveIndexSort (from,middle-1, index);
    }
}

//Append & store operator
void _SimpleList::RequestSpace (long slots)
{
    if (slots>laLength) {
        laLength=(slots/MEMORYSTEP+1)*MEMORYSTEP;
        if (lData) {
            checkPointer (lData = (long*)MemReallocate((char*)lData, laLength*sizeof(void*)));
        } else {
            checkPointer (lData = (long*)MemAllocate(laLength*sizeof(void*)));
        }
    }
}

// Compute the union of two sorted lists
// Each repeat appears exactly once
void    _SimpleList::Subtract (_SimpleList& l1, _SimpleList& l2)
{
    if (lLength) {
        Clear();
    }

    long  c1 = 0,
          c2 = 0;

    while (c1<l1.lLength && c2<l2.lLength) {
        while (c1<l1.lLength && l1.lData[c1]<l2.lData[c2]) {
            (*this) << l1.lData[c1++];
        }
        if (c1==l1.lLength) {
            break;
        }
        while ( c1<l1.lLength && c2<l2.lLength  && l1.lData[c1]==l2.lData[c2] ) {
            c1++;
            c2++;
        }
        if (c1==l1.lLength || c2==l2.lLength) {
            break;
        }
        while (c2<l2.lLength && l2.lData[c2]<l1.lData[c1]) {
            c2++;
        }
    }

    while (c1<l1.lLength) {
        (*this) << l1.lData[c1++];
    }
}

void _SimpleList::Swap (long i, long j)
{
    if ( i>=lLength || j>=lLength ) {
        return;
    }

    void * pt = ((void**)lData)[j];
    ((void**)lData)[j] = ((void**)lData)[i];
    ((void**)lData)[i] = pt;
}

//Char* conversion
BaseRef _SimpleList::toStr(void)
{
    if (lLength) {
        unsigned long ssi = _String::storageIncrement,
                      ma  = lLength*(1+log10((double)lLength));

        if ( ma > ssi) {
            _String::storageIncrement = ma;
        }

        _String * s = new _String (10L, true);

        checkPointer (s);

        (*s) << "{";

        for (unsigned long i = 0; i<lLength; i++) {
            char c[32];
            snprintf (c, sizeof(c),"%ld",((long*)lData)[i]),
                    (*s) << c;
            if (i<lLength-1) {
                (*s) << ',';
            }
        }

        (*s) << '}';

        s->Finalize();
        _String::storageIncrement = ssi;
        return s;
    } else {
        return new _String ("{}");
    }
}


//Delete item at index (>=0)
void  _SimpleList::TrimMemory (void)
{
    if (laLength>lLength) {
        laLength = lLength;
        if (laLength) {
            if (lData) {
                lData = (long*)MemReallocate ((char*)lData, laLength*sizeof(Ptr));
            } else {
                lData = (long*)MemAllocate (laLength*sizeof(Ptr));
            }
            if (!lData) {
                checkPointer (lData);
            }
        } else {
            if (lData) {
                free (lData);
                lData = nil;
            }
        }
    }
}

/*
==============================================================
Sort Methods
==============================================================
*/

void  _SimpleList::Sort (bool ascending)
{
    if (lLength<10) { // use bubble sort
        BubbleSort();
    } else {
        QuickSort(0,lLength-1);
    }

    if (!ascending) {
        long swap,i,j;
        for (i=0, j=lLength-1; i<j; i++,j--) {
            swap = ((long*)lData)[i];
            ((long*)lData)[i]=((long*)lData)[j];
            ((long*)lData)[j]=swap;
        }
    }
}

void  _SimpleList::BubbleSort (void)
{
    bool done = false;
    long swap,i,j;
    while (!done) {
        done = true;
        for (i=lLength-1,j=i-1; i>0; i--,j--) {
            if (Compare(i,j)<0) {
                done = false;
                swap = ((long*)lData)[i];
                ((long*)lData)[i]=((long*)lData)[j];
                ((long*)lData)[j]=swap;
            }
        }
    }
}

void  _SimpleList::QuickSort (long from, long to)
{
    long middle = (from+to)/2,
         middleV = ((long*)lData)[middle],
         top = to,
         bottommove = 1,
         topmove = 1,
         temp,
         i;

    if (middle)
        //while ((middle-bottommove>=from)&&(((long*)lData)[middle-bottommove]>middleV))
        while ((middle-bottommove>=from)&&(Compare(middle-bottommove,middle)>0)) {
            bottommove++;
        }

    if (from<to)
        //while ((middle+topmove<=to)&&(((long*)lData)[middle+topmove]<middleV))
        while ((middle+topmove<=to)&&(Compare(middle+topmove,middle)<0)) {
            topmove++;
        }
    // now shuffle
    for (i=from; i<middle-bottommove; i++) {
        if (Compare(i,middle)>0) {
            temp = ((long*)lData)[middle-bottommove];
            ((long*)lData)[middle-bottommove] = ((long*)lData)[i];
            ((long*)lData)[i]=temp;
            bottommove++;

            //while ((middle-bottommove>=from)&&(((long*)lData)[middle-bottommove]>middleV))
            while ((middle-bottommove>=from)&&(Compare(middle-bottommove,middle)>0)) {
                bottommove++;
            }
        }
    }

    for (i=middle+topmove+1; i<=top; i++) {
        if (Compare(i,middle)<0) {
            temp = ((long*)lData)[middle+topmove];
            ((long*)lData)[middle+topmove] = ((long*)lData)[i];
            ((long*)lData)[i]=temp;
            topmove++;

            //while ((middle+topmove<=to)&&(((long*)lData)[middle+topmove]<middleV))
            while ((middle+topmove<=to)&&(Compare(middle+topmove,middle)<0)) {
                topmove++;
            }
        }
    }

    if (topmove==bottommove) {
        for (i=1; i<bottommove; i++) {
            temp = ((long*)lData)[middle+i];
            ((long*)lData)[middle+i] = ((long*)lData)[middle-i];
            ((long*)lData)[middle-i]=temp;
        }
    } else if (topmove>bottommove) {
        long shift = topmove-bottommove;
        for (i=1; i<bottommove; i++) {
            temp = ((long*)lData)[middle+i+shift];
            ((long*)lData)[middle+i+shift] = ((long*)lData)[middle-i];
            ((long*)lData)[middle-i]=temp;
        }
        for (i=0; i<shift; i++) {
            ((long*)lData)[middle+i]=((long*)lData)[middle+i+1];
        }
        middle+=shift;
        ((long*)lData)[middle]=middleV;
    } else {
        long shift = bottommove-topmove;
        for (i=1; i<topmove; i++) {
            temp = ((long*)lData)[middle-i-shift];
            ((long*)lData)[middle-i-shift] = ((long*)lData)[middle+i];
            ((long*)lData)[middle+i]=temp;
        }
        for (i=0; i<shift; i++) {
            ((long*)lData)[middle-i]=((long*)lData)[middle-i-1];
        }
        middle-=shift;
        ((long*)lData)[middle]=middleV;
    }
    if (to>middle+1) {
        QuickSort (middle+1,top);
    }
    if (from<middle-1) {
        QuickSort (from,middle-1);
    }
}

//TODO: This is a global. Should it be here?
void SortLists (_SimpleList* ref, _SimpleList* index)
{
    if ((*ref).lLength!=index->lLength) {
        return;
    }
    if ((*ref).lLength<=10) {
        bool done = false;

        while (!done) {
            done = true;
            for (long i=1; i<(*ref).lLength; i++) {
                if (ref->Compare(i-1,i)>0) {
                    long swap;
                    swap = ((long*)ref->lData)[i];
                    ((long*)ref->lData)[i]=((long*)ref->lData)[i-1];
                    ((long*)ref->lData)[i-1]=swap;
                    swap = ((long*)index->lData)[i];
                    ((long*)index->lData)[i]=((long*)index->lData)[i-1];
                    ((long*)index->lData)[i-1]=swap;
                    done = false;
                }
            }
        }
    } else {
        (*ref).RecursiveIndexSort (0, (*ref).lLength-1,index);
    }
}

/*
==============================================================
Set Methods
==============================================================
*/

// Compute the intersection of two sorted lists
void    _SimpleList::Intersect (_SimpleList& l1, _SimpleList& l2)
{
    if (lLength) {
        Clear();
    }

    long  c1 = 0,
          c2 = 0;

    while (c1<l1.lLength && c2<l2.lLength) {
        while (l1.lData[c1]<l2.lData[c2]) {
            c1++;
            if (c1==l1.lLength) {
                break;
            }
        }
        if (c1==l1.lLength) {
            break;
        }

        while (l1.lData[c1]==l2.lData[c2]) {
            (*this) << l1.lData[c1++];
            c2++;
            if (c1==l1.lLength || c2==l2.lLength) {
                break;
            }
        }
        if (c1==l1.lLength || c2==l2.lLength) {
            break;
        }
        while (l2.lData[c2]<l1.lData[c1]) {
            c2++;
            if (c2==l2.lLength) {
                break;
            }
        }
    }
}

// Compute the union of two sorted lists
// Each repeat appears exactly once
void    _SimpleList::Union (_SimpleList& l1, _SimpleList& l2)
{
    if (lLength) {
        Clear();
    }

    long  c1 = 0,
          c2 = 0;

    while (c1<l1.lLength && c2<l2.lLength) {
        while (l1.lData[c1]<l2.lData[c2]) {
            (*this) << l1.lData[c1++];
            if (c1==l1.lLength) {
                break;
            }
        }

        if (c1==l1.lLength) {
            break;
        }

        while (l1.lData[c1]==l2.lData[c2]) {
            (*this) << l1.lData[c1++];
            c2++;
            if (c1==l1.lLength || c2==l2.lLength) {
                break;
            }
        }

        if (c1==l1.lLength || c2==l2.lLength) {
            break;
        }

        while (l2.lData[c2]<l1.lData[c1]) {
            (*this) << l2.lData[c2++];
            if (c2==l2.lLength) {
                break;
            }
        }
    }

    while (c1<l1.lLength) {
        (*this) << l1.lData[c1++];
    }
    while (c2<l2.lLength) {
        (*this) << l2.lData[c2++];
    }

}

// Compute the union of two sorted lists
// Each repeat appears exactly once
void    _SimpleList::XOR (_SimpleList& l1, _SimpleList& l2)
{
    if (lLength) {
        Clear();
    }

    long  c1 = 0,
          c2 = 0;

    while ((c1<l1.lLength)&&(c2<l2.lLength)) {
        while (c1<l1.lLength && l1.lData[c1]<l2.lData[c2]) {
            (*this) << l1.lData[c1++];
        }
        if (c1==l1.lLength) {
            break;
        }
        while (c1<l1.lLength && c2<l2.lLength && l1.lData[c1]==l2.lData[c2] ) {
            c1++;
            c2++;
        }
        if (c1==l1.lLength||c2==l2.lLength) {
            break;
        }

        while (c2<l2.lLength && l2.lData[c2]<l1.lData[c1]) {
            (*this) << l2.lData[c2++];
        }
    }

    while (c1<l1.lLength) {
        (*this) << l1.lData[c1++];
    }
    while (c2<l2.lLength) {
        (*this) << l2.lData[c2++];
    }
}

