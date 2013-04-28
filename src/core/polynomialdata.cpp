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

#include "string.h"
#include "defines.h"
#include "math.h"
#include "polynoml.h"
#include "polynomialdata.h"

extern _Parameter dropPrecision, drop2Precision, topPolyCap, dropTerms,
    enforcePolyCap, *varCheckArray, dropThreshold, maximumPolyTermsPerVariable,
    maxPolynomialExpIterates, polynomialExpPrecision;

extern long varCheckAllocated, polyTermCap;

extern bool checkReset;

//______________________________________________________________________________
_PolynomialData::_PolynomialData(void) {
  theCoeff = nil;
  thePowers = nil;
  numberVars = 0;
  actTerms = 0;
  allocTerms = 0;
}

//______________________________________________________________________________
_PolynomialData::_PolynomialData(long vars) {
  numberVars = vars >= 0 ? vars : 0;
  theCoeff =
      (_Parameter *)MemAllocate(sizeof(_Parameter) * POLY_DATA_INCREMENT);
  if (numberVars) {
    thePowers = (long *)MemAllocate(sizeof(long) * POLY_DATA_INCREMENT * vars);
  } else {
    thePowers = nil;
  }
  allocTerms = POLY_DATA_INCREMENT;
  actTerms = 0;
}

//______________________________________________________________________________
_PolynomialData::_PolynomialData(_PolynomialData &source) {
  Duplicate(&source);
}

//______________________________________________________________________________
_PolynomialData::_PolynomialData(long vars, long terms, _Parameter *theCoeffs) {
  numberVars = vars >= 0 ? vars : 0;
  allocTerms = (terms / POLY_DATA_INCREMENT + 1) * POLY_DATA_INCREMENT;
  actTerms = terms;
  theCoeff = (_Parameter *)MemAllocate(sizeof(_Parameter) * allocTerms);
  memcpy(theCoeff, theCoeffs, sizeof(_Parameter) * terms);
  thePowers = nil;
}

//______________________________________________________________________________
_PolynomialData::~_PolynomialData(void) {
  if (nInstances <= 1) {
    if (theCoeff) {
      free(theCoeff);
    }
    if (thePowers) {
      free(thePowers);
    }
    allocTerms = 0;
    actTerms = 0;
  }
}

//______________________________________________________________________________
BaseRef _PolynomialData::makeDynamic(void) {
  _PolynomialData *ret = new _PolynomialData;
  checkPointer(ret);
  nInstances++;
  ret->allocTerms = allocTerms;
  ret->actTerms = actTerms;
  ret->numberVars = numberVars;
  ret->theCoeff = theCoeff;
  ret->thePowers = thePowers;
  return ret;
}

//______________________________________________________________________________
void _PolynomialData::Duplicate(BaseRef source) {
  _PolynomialData *s = (_PolynomialData *)source;
  nInstances = 1;
  allocTerms = s->allocTerms;
  actTerms = s->actTerms;
  numberVars = s->numberVars;
  if (actTerms) {
    theCoeff = (_Parameter *)MemAllocate(sizeof(_Parameter) * allocTerms);
    memcpy(theCoeff, s->theCoeff, sizeof(_Parameter) * actTerms);
    if (numberVars) {
      thePowers = (long *)MemAllocate(sizeof(long) * allocTerms * numberVars);
      memcpy(thePowers, s->thePowers, sizeof(long) * actTerms * numberVars);
    } else {
      thePowers = nil;
    }

  } else {
    thePowers = nil;
    theCoeff = nil;
  }
}

//______________________________________________________________________________
bool _PolynomialData::checkMe(void) {
  if (actTerms > 1) {
    long *t1 = GetTerm(0), *t2;
    for (long i = 1; i < actTerms; i++) {
      t2 = GetTerm(i);
      if (CompareTerms(t1, t2) >= 0) {
        WarnError("\n Internal polynomial error!");
        return false;
      }
      t1 = t2;
    }
  }
  return true;
}

//______________________________________________________________________________
long *_PolynomialData::GetTerm(long index) {
  return thePowers + index * numberVars;
}

//______________________________________________________________________________
void _PolynomialData::AddTerm(long *theTerm, _Parameter theC) {
  if (!(allocTerms - actTerms)) { // no space left - reallocate
    allocTerms += POLY_DATA_INCREMENT;
    if (theCoeff) {
      theCoeff = (_Parameter *)MemReallocate((char *)theCoeff,
                                             allocTerms * sizeof(_Parameter));
    } else {
      theCoeff = (_Parameter *)MemAllocate(allocTerms * sizeof(_Parameter));
    }
    if (numberVars) {
      if (thePowers) {
        thePowers = (long *)MemReallocate(
            (char *)thePowers, allocTerms * sizeof(long) * numberVars);
      } else {
        thePowers = (long *)MemAllocate(allocTerms * sizeof(long) * numberVars);
      }
    }
  }
  theCoeff[actTerms] = theC;
  if (numberVars > 2) {
    long *stTerm = thePowers + actTerms * numberVars;
    for (long i = 0; i < numberVars; i++, theTerm++, stTerm++) {
      *stTerm = *theTerm;
    }
    //      memcpy (thePowers+actTerms*numberVars, theTerm,
    // numberVars*sizeof(long));
  } else {
    if (numberVars == 2) {
      *(thePowers + actTerms *numberVars) = *theTerm;
      *(thePowers + actTerms *numberVars + 1) = *(theTerm + 1);
    } else if (numberVars == 1) {
      *(thePowers + actTerms *numberVars) = *theTerm;
    }
  }
  actTerms++;
}

//______________________________________________________________________________
void _PolynomialData::WriteTerm(long *theTerm, long index) {
  if (numberVars > 2) {
    long *stTerm = thePowers + index * numberVars;
    for (long i = 0; i < numberVars; i++, theTerm++, stTerm++) {
      *stTerm = *theTerm;
    }
    // memcpy (thePowers+index*numberVars, theTerm,
    // numberVars*sizeof(long));
  } else {
    if (numberVars == 2) {
      *(thePowers + index *numberVars) = *theTerm;
      *(thePowers + index *numberVars + 1) = *(theTerm + 1);
    } else if (numberVars == 1) {
      *(thePowers + index *numberVars) = *theTerm;
    }
  }
  actTerms++;
}

//______________________________________________________________________________
void _PolynomialData::AddTerm(long *theTerm, _Parameter theC, long *reindexer,
                              long actLength) {

  if (!(allocTerms - actTerms)) { 
    // no space left - reallocate
    allocTerms += POLY_DATA_INCREMENT;
    if (theCoeff) {
      theCoeff = (_Parameter *)MemReallocate((char *)theCoeff,
                                             allocTerms * sizeof(_Parameter));
    } else {
      theCoeff = (_Parameter *)MemAllocate(allocTerms * sizeof(_Parameter));
    }
    if (numberVars) {
      if (thePowers) {
        thePowers = (long *)MemReallocate(
            (char *)thePowers, allocTerms * sizeof(long) * numberVars);
      } else {
        thePowers = (long *)MemAllocate(allocTerms * sizeof(long) * numberVars);
      }
    }
  }

  theCoeff[actTerms] = theC;

  if (numberVars > 2) {
    long *newTerm = thePowers + actTerms * numberVars, i;
    for (i = 0; i < numberVars; i++, newTerm++) {
      *newTerm = 0;
    }

    newTerm -= numberVars;
    //!     memset (newTerm, 0, numberVars*sizeof(long));
    for (i = 0; i < actLength; i++, reindexer++, theTerm++) {
      newTerm[*reindexer] = *theTerm;
    }

  } else {
    *(thePowers + actTerms *numberVars) = 0;
    *(thePowers + actTerms *numberVars + 1) = 0;
    *(thePowers + actTerms *numberVars + *reindexer) = *theTerm;
  }
  actTerms++;

}

//______________________________________________________________________________
void _PolynomialData::AddTerm(_Parameter theC) {
  if (numberVars == 0) {
    AddTerm(nil, 0);
  } else {
    if (!(allocTerms - actTerms)) { 
      // no space left - reallocate
      allocTerms += POLY_DATA_INCREMENT;
      if (theCoeff) {
        theCoeff = (_Parameter *)MemReallocate((char *)theCoeff,
                                               allocTerms * sizeof(_Parameter));
      } else {
        theCoeff = (_Parameter *)MemAllocate(allocTerms * sizeof(_Parameter));
      }
      if (numberVars) {
        if (thePowers) {
          thePowers = (long *)MemReallocate(
              (char *)thePowers, allocTerms * sizeof(long) * numberVars);
        } else {
          thePowers =
              (long *)MemAllocate(allocTerms * sizeof(long) * numberVars);
        }
      }
    }
    memmove(thePowers + numberVars, thePowers, numberVars * actTerms * sizeof(long));
    memmove(theCoeff + 1, theCoeff, sizeof(_Parameter) * actTerms);
    *theCoeff = theC;
    for (long i = 0; i < numberVars; *thePowers = 0, i++, thePowers++) {
    }
    thePowers -= numberVars;
    // memset (thePowers,0,numberVars*sizeof(long));
    actTerms++;
  }

}

//______________________________________________________________________________
void _PolynomialData::DeleteTerm(long index) {
  actTerms--;
  if (index != actTerms) { // shift stuff
    memmove((void *)(theCoeff + index), (void *)(theCoeff + index + 1),
            (actTerms - index) * sizeof(_Parameter));
    if (numberVars)
      memmove((void *)(thePowers + index * numberVars),
              (void *)(thePowers + (1 + index) * numberVars),
              (actTerms - index) * numberVars * sizeof(long));
  }
  if (allocTerms - actTerms > POLY_DATA_INCREMENT) {
    allocTerms -= POLY_DATA_INCREMENT;
    theCoeff = (_Parameter *)MemReallocate((char *)theCoeff,
                                           allocTerms * sizeof(_Parameter));
    if (numberVars) {
      thePowers = (long *)MemReallocate((char *)thePowers,
                                        allocTerms * sizeof(long) * numberVars);
    }
  }
}

//______________________________________________________________________________
long _PolynomialData::SumOfPowers(long index) {
  if (numberVars) {
    long *theTerm = GetTerm(index), res = 0;
    for (long i = 0; i < numberVars; i++, theTerm++) {
      res += *theTerm;
    }
    return res;
  }
  return 0;
}

//______________________________________________________________________________
long _PolynomialData::WeightedSumOfPowers(long index, _Parameter *w) {
  if (numberVars) {
    long *theTerm = GetTerm(index), res = 0;
    for (long i = 0; i < numberVars; i++, theTerm++, w++) {
      res += ((_Parameter) * theTerm) * *w;
    }
    return res;
  }
  return 0;
}

//______________________________________________________________________________
bool _PolynomialData::IsFirstANumber(void) {
  long *fst = GetTerm(0);
  for (long i = 0; i < numberVars; i++)
    if (fst[i]) {
      return false;
    }
  return true;
}

//______________________________________________________________________________
void _PolynomialData::MultiplyTerms(long *target, long *s1, long *s2) {
  for (long i = 0; i < numberVars; i++, target++, s1++, s2++) {
    *target = *s1 + *s2;
  }
}

//______________________________________________________________________________
void _PolynomialData::RaiseTerm(long *target, long power) {
  for (long i = 0; i < numberVars; i++, target++) {
    *target *= power;
  }
}

//______________________________________________________________________________
_Parameter _PolynomialData::BinaryRaise(_Parameter base, long pwr) {
  _Parameter result = 1;
  char bits[sizeof(long) * 8];
  unsigned char nLength = 0;
  while (pwr) {
    bits[nLength] = pwr % 2;
    pwr /= 2;
    nLength++;
  }
  while (nLength) {
    nLength--;
    if (bits[nLength]) {
      result *= base;
    }
    if (nLength) {
      result *= result;
    }
  }
  return result;
}

//______________________________________________________________________________
long _PolynomialData::FindTerm(long *theTerm, long *reIndex, long start) {
  long top = actTerms - 1, bottom = start, middle;
  char comp;

  if (top == -1) {
    return -2;
  }

  while (top > bottom) {
    middle = (top + bottom) / 2;
    // compare the two
    comp = CompareTerms(GetTerm(reIndex[middle]), theTerm);
    if (comp == 1) {
      top = middle == top ? top - 1 : middle;
    } else if (comp == -1) {
      bottom = middle == bottom ? bottom + 1 : middle;
    } else {
      return middle;
    }

  }
  middle = top;
  comp = CompareTerms(GetTerm(reIndex[middle]), theTerm);
  if (!comp) {
    return middle;
  }
  return comp < 0 ? -middle - 3 : -middle - 2;
}

//______________________________________________________________________________
void _PolynomialData::ResortTerms(long *reIndex) {
  long i, j, *source, *target, deleted = 0;

  _Parameter *newCoeff =
      (_Parameter *)MemAllocate(allocTerms * sizeof(_Parameter));
  long *newPowers = (long *)MemAllocate(allocTerms * numberVars * sizeof(long));

  // pass 1 to check for deletions
  for (i = 0; i < actTerms; i++, reIndex++) {
    if (checkTerm(theCoeff[*reIndex], *reIndex)) {
      newCoeff[i] = theCoeff[*reIndex];
    } else {
      newCoeff[i] = 0.0;
    }
  }

  reIndex -= actTerms;
  // pass 2 copy
  for (i = 0; i < actTerms; i++, newCoeff++, reIndex++) {
    if (*newCoeff != 0.0) {
      if (deleted) {
        *(newCoeff - deleted) = *newCoeff;
      }
      target = newPowers + numberVars * (i - deleted);
      source = thePowers + numberVars * (*reIndex);
      for (j = 0; j < numberVars; j++, source++, target++) {
        *target = *source;
      }
    } else {
      deleted++;
    }
  }
  free(theCoeff);
  free(thePowers);
  theCoeff = newCoeff - actTerms;
  thePowers = newPowers;
  actTerms -= deleted;
  if (allocTerms - actTerms > POLY_DATA_INCREMENT) {
    long theCut =
        ((allocTerms - actTerms) / POLY_DATA_INCREMENT) * POLY_DATA_INCREMENT;
    allocTerms -= theCut;
    theCoeff = (_Parameter *)MemReallocate((char *)theCoeff,
                                           allocTerms * sizeof(_Parameter));
    if (numberVars) {
      thePowers = (long *)MemReallocate((char *)thePowers,
                                        allocTerms * sizeof(long) * numberVars);
    }
  }

}

//______________________________________________________________________________
void _PolynomialData::ChopTerms(void) {
  long maxAllowedTerms = maximumPolyTermsPerVariable * numberVars;
  if (actTerms <= maxAllowedTerms) {
    return; // nothing to do
  } else {
    _SimpleList terms, index;
    long i, *qa, k, j;
    _Parameter logTop = log(topPolyCap), *qc;
    for (i = 0; i < actTerms; i++, theCoeff++) {
      index << i;
      terms << (long)(log(fabs(*theCoeff)) + logTop * SumOfPowers(i));
    }
    SortLists(&terms, &index);
    terms.Clear();
    theCoeff -= actTerms;
    qa = index.quickArrayAccess() + maxAllowedTerms;
    for (i = maxAllowedTerms; i < actTerms; i++, qa++) {
      theCoeff[*qa] = 0.0;
    }

    allocTerms =
        (maxAllowedTerms / POLY_DATA_INCREMENT + 1) * POLY_DATA_INCREMENT;
    _Parameter *newCoeff =
                   (_Parameter *)MemAllocate(allocTerms * sizeof(_Parameter)),
               *nP = newCoeff;
    long *newPowers =
             (long *)MemAllocate(allocTerms * numberVars * sizeof(long)),
         *target, *source;

    target = newPowers;
    source = thePowers;
    k = 0;
    qc = GetCoeff();
    for (i = 0; i < actTerms; i++, qc++, source += numberVars) {
      if (*qc != 0.0) {
        *nP = *qc;
        nP++;
        for (j = 0; j < numberVars; j++, source++, target++) {
          *target = *source;
        }
      } else {
        k++;
      }
    }
    free(theCoeff);
    free(thePowers);
    theCoeff = newCoeff;
    thePowers = newPowers;
    actTerms -= k;
  }

}

//______________________________________________________________________________
void _PolynomialData::RearrangeTerm(long *target, long *source, long *markup,
                                    long items) {
  for (long i = 0; i < items; i++, source++, markup++) {
    target[*markup] = *source;
  }
}

//______________________________________________________________________________
char _PolynomialData::CompareTerms(long *s1, long *s2) {
  for (long i = 0; i < numberVars; i++) {
    long comp = s1[i] - s2[i];
    if (comp > 0) {
      return 1;
    }
    if (comp < 0) {
      return -1;
    }
  }
  return 0;
}

//______________________________________________________________________________
char _PolynomialData::CompareTerms(long *s1, long *s2, long *secondReindex,
                                   long actLength) {
  long comp, i;
  for (i = 0; i < actLength; i++, s1++, s2++, secondReindex++) {
    comp = *secondReindex - i;
    if (comp > 0) {
      if (*s1) {
        return 1;
      } else {
        s2--;
        secondReindex--;
        actLength++;
        continue;
      }
    }
    comp = *s1 - *s2;
    if (comp > 0) {
      return 1;
    }
    if (comp < 0) {
      return -1;
    }
  }
  for (; i < numberVars; i++, s1++) {
    if (*s1) {
      return 1;
    }
  }
  return 0;
}

//______________________________________________________________________________
char _PolynomialData::CompareTerms(long *s1, long *s2, long *firstReindex,
                                   long *secondReindex, long actLength1,
                                   long actLength2) {
  bool secondLonger =
      actLength1<actLength2, firstLonger = actLength1> actLength2;

  long minLength = actLength1 < actLength2 ? actLength1 : actLength2, i;

  for (i = 0; i < minLength; i++, s1++, s2++, firstReindex++, secondReindex++) {
    long comp = *firstReindex - *secondReindex;
    if (comp < 0) {
      if (*s1) {
        return 1;
      }
      secondReindex--;
      s2--;

      if (firstLonger && minLength != actLength1) {
        minLength++;
      }

      continue;
    }
    if (comp > 0) {
      if (*s2) {
        return -1;
      }
      firstReindex--;
      s1--;

      if (minLength != actLength2 && secondLonger) {
        minLength++;
      }
      continue;
    }

    comp = *s1 - *s2;
    if (comp > 0) {
      return 1;
    }
    if (comp < 0) {
      return -1;
    }
  }
  if (actLength1 > minLength)
    for (; i < actLength1; i++, s1++) {
      if (*s1) {
        return 1;
      }
    }
  else
    for (; i < actLength2; i++, s2++) {
      if (*s2) {
        return -1;
      }
    }
  return 0;
}

//______________________________________________________________________________
bool _PolynomialData::checkTerm(_Parameter myCoeff, long myIndex) {
  if (myCoeff == 0.0) { // delete in any case
    return false;
  }
  if (checkReset) {
    checkReset = false;
    dropThreshold = dropPrecision + log(fabs(myCoeff));
    if (dropThreshold < drop2Precision) {
      dropThreshold = drop2Precision;
      if (enforcePolyCap) {
        dropThreshold += SumOfPowers(myIndex) * log(topPolyCap);
      } else {
        dropThreshold += WeightedSumOfPowers(myIndex, varCheckArray);
      }
      return false;
    }
    if (enforcePolyCap) {
      dropThreshold += SumOfPowers(myIndex) * log(topPolyCap);
    } else {
      dropThreshold += WeightedSumOfPowers(myIndex, varCheckArray);
    }
    return true;

  }
  if (dropTerms) {
    if (enforcePolyCap) {
      if (log(fabs(myCoeff)) +
              log(topPolyCap) * (_Parameter)(SumOfPowers(myIndex)) <
          dropThreshold) {
        return false;
      }
    } else if (enforcePolyCap) {
      if (log(fabs(myCoeff)) + WeightedSumOfPowers(myIndex, varCheckArray) <
          dropThreshold) {
        return false;
      }
    }
  }
  return true;
}
