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

/*
 ==============================================================
 Constructors
 ==============================================================
 */ // Does nothing
template<typename PAYLOAD>
_hyListOrderable<PAYLOAD>::_hyListOrderable() {
}

  //Data constructor (1 member list)
template<typename PAYLOAD>
_hyListOrderable<PAYLOAD>::_hyListOrderable(const PAYLOAD item) : _hyList<PAYLOAD> (item) {
}


  //Stack copy contructor
template<typename PAYLOAD>
_hyListOrderable<PAYLOAD>::_hyListOrderable(const _hyListOrderable <PAYLOAD> &l, const long from, const long to) {
  this->Initialize ();
  this->Clone (&l, from, to);
}

  // Data constructor (variable number of long constants)
template<typename PAYLOAD>
_hyListOrderable<PAYLOAD>::_hyListOrderable(const unsigned long number, const PAYLOAD items[]) : 
_hyList<PAYLOAD> (number, items) {
}

/*
==============================================================
Search/Insert Functions
==============================================================
*/

template<typename PAYLOAD>
long _hyListOrderable<PAYLOAD>::BinaryFind(const PAYLOAD item, const long startAt) const {

  
  long top = this->lLength - 1L,
       bottom = startAt,
       middle;

  if (top == -1L) {
    return -2L;
  }

  while (top > bottom) {
    middle = (top + bottom) >> 1;
    if (CompareToValue (middle, item) == HY_COMPARE_GREATER) {
      top = middle == top ? top - 1L : middle;
    } else if (CompareToValue (middle, item) == HY_COMPARE_LESS) {
      bottom = middle == bottom ? bottom + 1 : middle;
    } else {
      return middle;
    }
  }

  middle = top;
  if (this->ItemEqualToValue (middle,item)) {
    return middle;
  } else {
    if (CompareToValue(middle,item) == HY_COMPARE_LESS) {
      return -middle - 3;
    }
  }

  return -middle - 2;
}

template<typename PAYLOAD>
long _hyListOrderable<PAYLOAD>::BinaryInsert(const PAYLOAD item) {
  if (this->lLength == 0UL) {
    this->append (item);
    return 0L;
  }

  long pos = -BinaryFind(item) - 2;

  if (pos < 0L) {
    return -pos + 2;
  }

  if (CompareToValue (pos,item) == HY_COMPARE_LESS) {
    pos++;
  }

  this->InsertElement(item, pos);

  return pos >= this->lLength ? this->lLength - 1 : pos;
}


/*
==============================================================
Summary Functions 
==============================================================
*/

template<typename PAYLOAD>
PAYLOAD _hyListOrderable<PAYLOAD>::Max(void) const{
  PAYLOAD res = this->lData[0];
  for (unsigned long e = 1UL; e < this->lLength; e++)
    if (CompareToValue (e, res) == HY_COMPARE_GREATER) {
      res = this->lData[e];
    }
  return res;
}

template<typename PAYLOAD>
PAYLOAD _hyListOrderable<PAYLOAD>::Min(void) const{
  PAYLOAD res = this->lData[0];
  for (unsigned long e = 1UL; e < this->lLength; e++)
    if (CompareToValue (e, res) == HY_COMPARE_LESS) {
      res = this->lData[e];
    }
  return res;
}




/*
==============================================================
Sorting Functions 
==============================================================
*/


template<typename PAYLOAD>
bool _hyListOrderable<PAYLOAD>::IsSorted (void) const {
  
  for (unsigned long i = 0UL; 1UL + i < this->countitems(); i++) {
    if (this->Compare (i, i+1) == HY_COMPARE_GREATER) {
      return false;
    }
  }
  
  return true;
}

template<typename PAYLOAD>
long _hyListOrderable<PAYLOAD>::Compare(const long i, const long j) const {
  /*
  if (i < 0 || j < 0 || i > (long)this->lLength || j > (long)this->lLength) {
    printf ("Comparing out of bounds %ld %ld %ld\n", i, j, this->lLength);
  }
  */
  if (this->lData[i] < this->lData[j]) {
    return HY_COMPARE_LESS;
  } else if (this->lData[i] == this->lData[j]) {
    return HY_COMPARE_EQUAL;
  }
  return HY_COMPARE_GREATER;
}

template<typename PAYLOAD>
long _hyListOrderable<PAYLOAD>::CompareToValue(const long item, const PAYLOAD& value) const {
  
  if (this->lData[item] < value) {
    return HY_COMPARE_LESS;
  } else if (this->lData[item] > value) {
    return HY_COMPARE_GREATER;
  }
  return HY_COMPARE_EQUAL;
}


template<typename PAYLOAD>
void _hyListOrderable<PAYLOAD>::BubbleSort(_hyListOrderable <long> * index_list) {
  bool done = false;
  while (!done) {
    done = true;
    for (unsigned long i = 0UL, j = i + 1UL; i + 1UL < this->countitems(); i++, j++) {
      if (this->Compare(i, j) == HY_COMPARE_GREATER) {
        done = false;
        this->Swap (i,j);
        if (index_list) {
          index_list->Swap (i,j);
        }
      }
    }
  }
}

template<typename PAYLOAD>
void _hyListOrderable<PAYLOAD>::QuickSort(const unsigned long from, const unsigned long to, _hyListOrderable <long> * index_list) {

  if (from >= to || from >= this->countitems()) {
    return;
  } 
  if (to >= this->countitems() && from < this->countitems()) {
     this->QuickSort (from, this->countitems() - 1, index_list); 
     return;
  }

  unsigned long middle = (from + to) >> 1, 
                top = to,
                bottommove = 1UL, 
                topmove = 1UL;
                
  long          imiddleV = index_list?index_list->AtIndex(middle):0L;
         
  PAYLOAD middleV = this->AtIndex(middle);
  
  if (middle) {
    while (middle >= from + bottommove &&
           this->Compare(middle - bottommove, middle) == HY_COMPARE_GREATER) {
      bottommove++;
    }
  }

  if (from < to) {
    while (middle + topmove <= to &&
           this->Compare(middle + topmove, middle) == HY_COMPARE_LESS) {
      topmove++;
    }
  }
  
  // now shuffle
  for (unsigned long i = from; i + bottommove < middle; i++) {
    if (this->Compare(i, middle) == HY_COMPARE_GREATER) {
      this->Swap (i, middle-bottommove);
      if (index_list) {
        index_list->Swap (i, middle-bottommove);
      }
      bottommove++;

      while (middle >= from + bottommove &&
             Compare(middle - bottommove, middle) == HY_COMPARE_GREATER) {
        bottommove++;
      }
    }
  }

  for (unsigned long i = middle + topmove + 1; i <= top; i++) {
    if (Compare(i, middle) == HY_COMPARE_LESS) {
      this->Swap (i, middle+topmove);
      if (index_list) {
        index_list->Swap (i,middle+topmove);
      }
      topmove++;

      while (middle + topmove <= to &&
             Compare(middle + topmove, middle) == HY_COMPARE_LESS) {
        topmove++;
      }
    }
  }

  if (topmove == bottommove) {
    for (unsigned long i = 1UL; i < bottommove; i++) {
      this->Swap (middle + i, middle -i);
      if (index_list) {
        index_list->Swap (middle + i, middle -i);
      }
    }
  } else if (topmove > bottommove) {
    unsigned long shift = topmove - bottommove;
    for (unsigned long i = 1UL; i < bottommove; i++) {
      this->Swap (middle-i, middle+i+shift);
      if (index_list) {
        index_list->Swap (middle-i, middle+i+shift);
      }
    }
    for (unsigned long i = 0UL; i < shift; i++) {
      this->lData[middle+i] = this->lData[middle+i+1];
      if (index_list) {
        (*index_list)[middle+i] = index_list->AtIndex(middle+i+1);
      }
    }
    middle += shift;
    this->lData[middle] = middleV;
    if (index_list) {
      (*index_list)[middle] = imiddleV;
    }
  } else {
    unsigned long shift = bottommove - topmove;
    for (unsigned long i = 1UL; i < topmove; i++) {
      this->Swap (middle+i, middle-i-shift);
      if (index_list) {
        index_list->Swap (middle+i, middle-i-shift);
      }
    }
    for (unsigned long i = 0UL; i < shift; i++) {
      this->lData[middle-i] = this->lData[middle-i-1];
      if (index_list) {
        (*index_list)[middle-i] = index_list->AtIndex(middle-i-1);
      }
    }
    middle -= shift;
    this->lData[middle] = middleV;
    if (index_list) {
      (*index_list)[middle] = imiddleV;
    }
  }
  if (to > middle + 1) {
    QuickSort(middle + 1, top, index_list);
  }
  if (from + 1 < middle) {
    QuickSort(from, middle - 1,index_list);
  }
}

template<typename PAYLOAD>
void _hyListOrderable<PAYLOAD>::Sort(bool ascending, _hyListOrderable <long> * index_list) {
  if (this->countitems() < 10UL) { 
    this->BubbleSort(index_list);
  } else {
    this->QuickSort(0, this->lLength - 1, index_list);
  }

  if (!ascending) {
    this->Flip();
    if (index_list) {
      index_list->Flip();
    }
  }
}

/*
==============================================================
Set Operations
==============================================================
*/

//Delete duplicates from a sorted list
template<typename PAYLOAD>
void _hyListOrderable<PAYLOAD>::DeleteDuplicates(void) {
  if (this->lLength > 1) {
    _hyListOrderable <PAYLOAD> noDups;

    PAYLOAD lastValue = this->lData[0];
    noDups << lastValue;
    
    for (unsigned long k = 1; k < this->lLength; k++) {
      PAYLOAD thisValue = this->lData[k];
      if (!this->EqualToValue (k,lastValue)) {
        noDups << thisValue;
        lastValue = thisValue;
      }
    }

    if (noDups.lLength != this->lLength) {
      Duplicate(&noDups);
    }
  }
}

template<typename PAYLOAD>
long _hyListOrderable<PAYLOAD>::CountCommonElements(const _hyListOrderable &l1, bool at_least_one) const {
  long c1 = 0L, c2 = 0L, res = 0L;

  while (c1 < l1.lLength && c2 < this->lLength) {
    while (l1.CompareToValue (c1,this->lData[c2]) == HY_COMPARE_LESS) {
      c1++;
      if (c1 == l1.lLength) {
        break;
      }
    }
    if (c1 == l1.lLength) {
      break;
    }
    
    while (l1.EqualToValue (c1, this->lData[c2])) {
      c2++;
      if (at_least_one) {
        return 1;
      } else {
        res++;
      }
      if (c1 == l1.lLength || c2 == this->lLength) {
        break;
      }
    }
    if (c1 == l1.lLength || c2 == this->lLength) {
      break;
    }
    while (CompareToValue (c2, l1.lData[c1]) == HY_COMPARE_LESS) {
      c2++;
      if (c2 == this->lLength) {
        break;
      }
    }
  }
  return res;
}

template<typename PAYLOAD>
void _hyListOrderable<PAYLOAD>::FilterRange(const PAYLOAD lb, const PAYLOAD ub) {
  if (ub <= lb) {
    this->Clear();
  } else {
    _hyListOrderable <long> toDelete;
    for (unsigned long k = 0UL; k < this->lLength; k++)
      if (CompareToValue (k,lb) != HY_COMPARE_GREATER || CompareToValue (k,ub) != HY_COMPARE_LESS) {
        toDelete << k;
      }
    this->DeleteList(&toDelete);
  }
}

template<typename PAYLOAD>
void _hyListOrderable<PAYLOAD>::Union(const _hyListOrderable <PAYLOAD> & l1,const _hyListOrderable <PAYLOAD> & l2) {
  this->Clear();

  long c1 = 0, c2 = 0;

  while (c1 < l1.lLength && c2 < l2.lLength) {
    while (l1.lData[c1] < l2.lData[c2]) {
      (*this) << l1.lData[c1++];
      if (c1 == l1.lLength) {
        break;
      }
    }

    if (c1 == l1.lLength) {
      break;
    }

    while (l1.lData[c1] == l2.lData[c2]) {
      (*this) << l1.lData[c1++];
      c2++;
      if (c1 == l1.lLength || c2 == l2.lLength) {
        break;
      }
    }

    if (c1 == l1.lLength || c2 == l2.lLength) {
      break;
    }

    while (l2.lData[c2] < l1.lData[c1]) {
      (*this) << l2.lData[c2++];
      if (c2 == l2.lLength) {
        break;
      }
    }
  }

  while (c1 < l1.lLength) {
    (*this) << l1.lData[c1++];
  }
  while (c2 < l2.lLength) {
    (*this) << l2.lData[c2++];
  }

}

template<typename PAYLOAD>
void _hyListOrderable<PAYLOAD>::XOR(const _hyListOrderable <PAYLOAD> & l1,const _hyListOrderable <PAYLOAD> & l2) {
  this->Clear();

  long c1 = 0, c2 = 0;

  while (c1 < l1.lLength && c2 < l2.lLength) {
    while (c1 < l1.lLength && l1.lData[c1] < l2.lData[c2]) {
      (*this) << l1.lData[c1++];
    }
    if (c1 == l1.lLength) {
      break;
    }
    while (c1 < l1.lLength && c2 < l2.lLength && l1.lData[c1] == l2.lData[c2]) {
      c1++;
      c2++;
    }
    if (c1 == l1.lLength || c2 == l2.lLength) {
      break;
    }

    while (c2 < l2.lLength && l2.lData[c2] < l1.lData[c1]) {
      (*this) << l2.lData[c2++];
    }
  }

  while (c1 < l1.lLength) {
    (*this) << l1.lData[c1++];
  }
  while (c2 < l2.lLength) {
    (*this) << l2.lData[c2++];
  }
}


template<typename PAYLOAD>
void _hyListOrderable<PAYLOAD>::Intersect(const _hyListOrderable <PAYLOAD> & l1,const _hyListOrderable <PAYLOAD> & l2) {
  this->Clear();

  long c1 = 0, c2 = 0;

  while (c1 < l1.lLength && c2 < l2.lLength) {
    while (l1.lData[c1] < l2.lData[c2]) {
      c1++;
      if (c1 == l1.lLength) {
        break;
      }
    }
    if (c1 == l1.lLength) {
      break;
    }

    while (l1.lData[c1] == l2.lData[c2]) {
      (*this) << l1.lData[c1++];
      c2++;
      if (c1 == l1.lLength || c2 == l2.lLength) {
        break;
      }
    }
    if (c1 == l1.lLength || c2 == l2.lLength) {
      break;
    }
    while (l2.lData[c2] < l1.lData[c1]) {
      c2++;
      if (c2 == l2.lLength) {
        break;
      }
    }
  }
}

template<typename PAYLOAD>
void _hyListOrderable<PAYLOAD>::Subtract(const _hyListOrderable <PAYLOAD> & l1,const _hyListOrderable <PAYLOAD> & l2) {
  this->Clear();

  long c1 = 0, c2 = 0;

  while (c1 < l1.lLength && c2 < l2.lLength) {
    while (c1 < l1.lLength && l1.lData[c1] < l2.lData[c2]) {
      (*this) << l1.lData[c1++];
    }
    if (c1 == l1.lLength) {
      break;
    }
    while (c1 < l1.lLength && c2 < l2.lLength && l1.lData[c1] == l2.lData[c2]) {
      c1++;
      c2++;
    }
    if (c1 == l1.lLength || c2 == l2.lLength) {
      break;
    }
    while (c2 < l2.lLength && l2.lData[c2] < l1.lData[c1]) {
      c2++;
    }
  }

  while (c1 < l1.lLength) {
    (*this) << l1.lData[c1++];
  }
}


template<typename PAYLOAD>
void _hyListOrderable<PAYLOAD>::Merge(const _hyListOrderable <PAYLOAD> &l1, const _hyListOrderable <PAYLOAD> &l2, 
  _hyListOrderable <long> *mergeResults1, _hyListOrderable <long> *mergeResults2)
{
  this->Clear();
  if (mergeResults1) {
    mergeResults1->Clear();
  }
  if (mergeResults2) {
    mergeResults2->Clear();
  }
  
  enum    machine_states {INIT, ADVANCE1, ADVANCE2, FLUSH1, FLUSH2} advancing = INIT;

  unsigned long    pos1 = 0,
  pos2 = 0,
  nt1 = l1.lLength,
  nt2 = l2.lLength;
  
  
  bool   doMerge1 = false,
  doMerge2 = false;
  
  while (1) { // stuff left to do
    if (advancing == ADVANCE1) { // advancing in the 1st list
      pos1++;
      if (pos1==nt1) {
        advancing = FLUSH2;
        continue;
      }

      if (l1->lData[pos1] <= l2->lData[pos2]) {
        if (mergeResults1 && doMerge1) {
          (*mergeResults1)<<this->lLength;
        }
        
        this->append (l1->lData[pos1]);
          
        if (mergeResults2 && l1->lData[pos1] < l2->lData[pos2]) {
          if (!doMerge2) {
            if (pos1>=pos2) {
              doMerge2 = true;
              for (unsigned long i=0UL; i<pos2; i++) {
                mergeResults2->append (i);
              }
            }
          }
          continue;
        }
      }
      
      if (l1->lData[pos1] > l2->lData[pos2]) {
        advancing = ADVANCE1;
        if (mergeResults1 && !doMerge1) {
          for (unsigned long i=0UL; i<pos1; i++) {
            mergeResults1->append (i);
          }
          doMerge1 = true;
        }
        
        if (mergeResults2 && doMerge2) {
          (*mergeResults2)<<this->lLength;
        }
        
        this->append (l2->lData[pos2]);
        continue;
      }
      
    } else if (advancing == ADVANCE2) { // advancing in the 2nd list
      pos2++;
      if (pos2==nt2) {
        advancing = FLUSH1;
        continue;
      }

      if (l2->lData[pos2] <= l1->lData[pos1]) {
        if (mergeResults2 && doMerge2) {
          (*mergeResults2)<<this->lLength;
        }
        this->append (l2->lData[pos2]);

        if (l2->lData[pos2] < l1->lData[pos1]) {
          if (mergeResults1 && !doMerge1) {
            if (pos2>=pos1) {
              doMerge1 = true;
              for (unsigned long i=0UL; i<pos1; i++) {
                (*mergeResults1)<<i;
              }
            }
          }
          continue;
        }
      }
      
      if (l2->lData[pos2] > l1->lData[pos1]) {
        advancing = ADVANCE1;
        if (mergeResults2 && !doMerge2) {
          for (unsigned long i=0UL; i<pos2; i++) {
            (*mergeResults2)<<i;
          }
          doMerge2 = true;
        }
        if (mergeResults1 && doMerge1) {
          (*mergeResults1)<<this->lLength;
        }
        this->append (l1->lData[pos1]);
        continue;
      }
    } else if (advancing == FLUSH2) { // flush out the 2nd list
      if (mergeResults1 && !doMerge1&& pos2<nt2 ) {
        for (unsigned long i=0UL; i<nt1; i++) {
          (*mergeResults1)<<i;
        }
      }
      if (mergeResults2 && doMerge2)
        while (pos2<nt2) {
          (*mergeResults2)<<this->lLength;
          this->append (l2->lData[pos2++]);
        }
      else
        while (pos2<nt2) {
          this->append (l2->lData[pos2++]);
        }
      break;
    } else if (advancing == FLUSH1) { // flush out the 1st list
      if (mergeResults2 && !doMerge2 && pos1<nt1) {
        for (unsigned long i=0UL; i<nt2; i++) {
          (*mergeResults2)<<i;
        }
      }
      if (mergeResults1 && doMerge1)
        while (pos1<nt1) {
          (*mergeResults1)<<this->lLength;
          this->append (l1->lData[pos1++]);
        }
      else
        while (pos1<nt1) {
          this->append (l1->lData[pos1++]);
        }
      break;
    } else if (advancing == INIT) { // just starting
      if (!nt1) { // first list is empty!
        advancing = FLUSH2;
        continue;
      }
      if (!nt2) { // second list is empty!
        advancing = FLUSH1;
        continue;
      }
      
      if (l1->lData[pos1] <= l2->lData[pos2]) { // begin with the first list
        this->append (l1->lData[pos1]);
        advancing = ADVANCE1;
        if (mergeResults2 && l1->lData[pos1] != l2->lData[pos2]) {
          doMerge2 = true;
          continue;
        }
      } else {
        this->append (l1->lData[pos1]);
        advancing = ADVANCE2;
        doMerge1 = true;
        continue;
      }
      
    }
    
    if (advancing == ADVANCE1) { // moving up in the second term
      pos1++;
      if (pos1==nt1) {
        pos2++;
        if (mergeResults2 && doMerge2) {
          (*mergeResults2)<<this->lLength-1UL;
        }
        advancing = FLUSH2;
        continue;
      } else {
        advancing = ADVANCE2;
        if (mergeResults2 && doMerge2) {
          (*mergeResults2)<<this->lLength-1UL;
        }
      }
    } else {
      pos2++;
      if (pos2==nt2) {
        pos1++;
        if (mergeResults1 && doMerge1) {
          (*mergeResults1)<<this->lLength-1UL;
        }
        advancing = FLUSH1;
        continue;
      } else {
        advancing = ADVANCE2;
        if (mergeResults1 && doMerge1) {
          (*mergeResults1)<<this->lLength-1UL;
        }
      }
    }
  }
}


