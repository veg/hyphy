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

#include <math.h>
#include <string.h>

#include "function_templates.h"
#include "global_things.h"
#include "polynoml.h"

using namespace hy_global;

hyFloat dropPrecision = log(1e-9), drop2Precision = log(1e-18),
        topPolyCap = 5.0, *varCheckArray = nil, dropThreshold = 0.0,
        polynomialExpPrecision = 1e-10;

bool dropTerms = true, enforcePolyCap = true;

long varCheckAllocated = 0, polyTermCap = 0x1fffffff;

bool checkReset = true;

const unsigned long maximumPolyTermsPerVariable = 500UL,
                    maxPolynomialExpIterates = 20UL;

//__________________________________________________________________________________

HBLObjectRef _returnPolyOrUseCache(HBLObjectRef cache) {
  if (cache && cache->ObjectClass() == HY_UNDEFINED) {
    return cache;
  }
  return new _MathObject;
}
//__________________________________________________________________________________
/**
 * @brief Construct a new _PolynomialData object.
 *
 */
_PolynomialData::_PolynomialData(void)
    : BaseObj(), theCoeff(nullptr), thePowers(nullptr), numberVars(0),
      actTerms(0), allocTerms(0) {}

//__________________________________________________________________________________
/**
 * @brief Construct a new _PolynomialData object.
 * @param vars The number of variables.
 */
_PolynomialData::_PolynomialData(const long vars) {
  numberVars = vars >= 0 ? vars : 0;
  theCoeff = (hyFloat *)MemAllocate(sizeof(hyFloat) * POLY_DATA_INCREMENT);
  if (numberVars) {
    thePowers = (long *)MemAllocate(sizeof(long) * POLY_DATA_INCREMENT * vars);
  } else {
    thePowers = nil;
  }
  allocTerms = POLY_DATA_INCREMENT;
  actTerms = 0;
}

//__________________________________________________________________________________

/**
 * @brief Construct a new _PolynomialData object.
 * @param source The _PolynomialData object to copy.
 */
_PolynomialData::_PolynomialData(_PolynomialData const &source) {
  _PolynomialData::Duplicate(&source);
}

//__________________________________________________________________________________
/**
 * @brief Assign a _PolynomialData object to this object.
 * @param source The _PolynomialData object to assign.
 * @return This object.
 *
 */
_PolynomialData const &
_PolynomialData::operator=(_PolynomialData const &source) {
  if (this != &source) {
    Duplicate(&source);
  }
  return *this;
}

//__________________________________________________________________________________
/**
 * @brief Construct a new _PolynomialData object.
 * @param vars The number of variables.
 * @param terms The number of terms.
 * @param theCoeffs The coefficients.
 */
_PolynomialData::_PolynomialData(const long vars, const long terms,
                                 hyFloat const *theCoeffs) {
  numberVars = vars >= 0 ? vars : 0;
  allocTerms = (terms / POLY_DATA_INCREMENT + 1) * POLY_DATA_INCREMENT;
  actTerms = terms;
  theCoeff = (hyFloat *)MemAllocate(sizeof(hyFloat) * allocTerms);
  memcpy(theCoeff, theCoeffs, sizeof(hyFloat) * terms);
  thePowers = nil;
}

//__________________________________________________________________________________

/**
 * @brief Destroy the _PolynomialData object.
 */
_PolynomialData::~_PolynomialData(void) {
  if (CanFreeMe()) {
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

//__________________________________________________________________________________
/**
 * @brief Make a dynamic copy of this object.
 * @return A pointer to the new object.
 */
BaseRef _PolynomialData::makeDynamic(void) const {
  _PolynomialData *ret = new _PolynomialData;
  ret->allocTerms = allocTerms;
  ret->actTerms = actTerms;
  ret->numberVars = numberVars;
  ret->theCoeff = theCoeff;
  ret->thePowers = thePowers;
  return ret;
}

//__________________________________________________________________________________
/**
 * @brief Duplicate this object.
 * @param source The object to duplicate.
 */
void _PolynomialData::Duplicate(BaseRefConst source) {
  _PolynomialData const *s = (_PolynomialData const *)source;
  allocTerms = s->allocTerms;
  actTerms = s->actTerms;
  numberVars = s->numberVars;
  if (actTerms) {
    theCoeff = (hyFloat *)MemAllocate(sizeof(hyFloat) * allocTerms);
    memcpy(theCoeff, s->theCoeff, sizeof(hyFloat) * actTerms);
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

//__________________________________________________________________________________
/**
 * @brief Check the polynomial data.
 * @return True if the polynomial data is valid, false otherwise.
 */
bool _PolynomialData::checkMe(void) {
  if (actTerms > 1) {
    long *t1 = GetTerm(0), *t2;
    for (long i = 1; i < actTerms; i++) {
      t2 = GetTerm(i);
      if (CompareTerms(t1, t2) >= 0) {
        HandleApplicationError("Internal polynomial error!", false);
        return false;
      }
      t1 = t2;
    }
  }
  return true;
}

//__________________________________________________________________________________
long *_PolynomialData::GetTerm(long index) const {
  return thePowers + index * numberVars;
}

//__________________________________________________________________________________
void _PolynomialData::EnsureSpaceIsAvailable(void) {
  if (allocTerms == actTerms) { // no space left - reallocate
    allocTerms += POLY_DATA_INCREMENT;
    theCoeff = (hyFloat *)MemReallocate(theCoeff, allocTerms * sizeof(hyFloat));
    if (numberVars) {
      thePowers = (long *)MemReallocate(thePowers,
                                        allocTerms * sizeof(long) * numberVars);
    }
  }
}

//__________________________________________________________________________________
void _PolynomialData::ShrinkExtraSpace(void) {
  if (allocTerms - actTerms > POLY_DATA_INCREMENT) {
    allocTerms -= POLY_DATA_INCREMENT;
    theCoeff = (hyFloat *)MemReallocate((char *)theCoeff,
                                        allocTerms * sizeof(hyFloat));
    if (numberVars) {
      thePowers = (long *)MemReallocate((char *)thePowers,
                                        allocTerms * sizeof(long) * numberVars);
    }
  }
}

//__________________________________________________________________________________
/**
 * @brief Add a term to the polynomial.
 * @param theTerm The term to add.
 * @param theC The coefficient of the term.
 */
void _PolynomialData::AddTerm(long *theTerm, hyFloat theC) {
  EnsureSpaceIsAvailable();
  theCoeff[actTerms] = theC;
  if (theTerm) {
    memcpy(thePowers + actTerms * numberVars, theTerm,
           sizeof(long) * numberVars);
  } else {
    memset(thePowers + actTerms * numberVars, 0, sizeof(long) * numberVars);
  }
  actTerms++;
}

//__________________________________________________________________________________
/**
 * @brief Write a term to the polynomial.
 * @param theTerm The term to write.
 * @param index The index of the term.
 */
void _PolynomialData::WriteTerm(long *theTerm, long index) {
  if (numberVars > 0) {
    memcpy(thePowers + index * numberVars, theTerm, numberVars * sizeof(long));
  }
  actTerms++;
}

//__________________________________________________________________________________
/**
 * @brief Add a term to the polynomial.
 * @param theTerm The term to add.
 * @param theC The coefficient of the term.
 * @param reindexer A pointer to the reindexer.
 * @param actLength The active length.
 */
void _PolynomialData::AddTerm(long *theTerm, hyFloat theC, long *reindexer,
                              long actLength) {
  EnsureSpaceIsAvailable();
  theCoeff[actTerms] = theC;
  if (numberVars > 0) {
    long *newTerm = thePowers + actTerms * numberVars;
    memset(newTerm, 0, numberVars * sizeof(long));
    for (long i = 0; i < actLength; i++)
      newTerm[reindexer[i]] = theTerm[i];
  }
  actTerms++;
}

//__________________________________________________________________________________
/**
 * @brief Add a constant  term to the polynomial.
 * @param theC The coefficient of the term.
 */
void _PolynomialData::AddTerm(hyFloat theC) {
  EnsureSpaceIsAvailable();
  memmove(thePowers + numberVars, thePowers,
          numberVars * actTerms * sizeof(long));
  memmove(theCoeff + 1, theCoeff, sizeof(hyFloat) * actTerms);
  *theCoeff = theC;
  memset(thePowers, 0, numberVars * sizeof(long));
  actTerms++;
}

//__________________________________________________________________________________
/**
 * @brief Delete a term from the polynomial.
 * @param index The index of the term to delete.
 */
void _PolynomialData::DeleteTerm(long index) {
  actTerms--;
  if (index != actTerms) { // shift stuff
    memmove((void *)(theCoeff + index), (void *)(theCoeff + index + 1),
            (actTerms - index) * sizeof(hyFloat));
    if (numberVars)
      memmove((void *)(thePowers + index * numberVars),
              (void *)(thePowers + (1 + index) * numberVars),
              (actTerms - index) * numberVars * sizeof(long));
  }
  ShrinkExtraSpace();
}

//__________________________________________________________________________________
/**
 * @brief Get the sum of the powers of a term.
 * @param index The index of the term.
 * @return The sum of the powers.
 */
long _PolynomialData::SumOfPowers(long index) const {
  if (numberVars) {
    long *theTerm = GetTerm(index), res = 0;
    for (long i = 0; i < numberVars; i++) {
      res += theTerm[i];
    }
    return res;
  }
  return 0;
}

//__________________________________________________________________________________
/**
 * @brief Get the weighted sum of the powers of a term.
 * @param index The index of the term.
 * @param w The weights.
 * @return The weighted sum of the powers.
 */
hyFloat _PolynomialData::WeightedSumOfPowers(long index, hyFloat *w) const {
  if (numberVars) {
    long *theTerm = GetTerm(index);
    hyFloat res = 0.;
    for (long i = 0; i < numberVars; i++) {
      res += w[i] * theTerm[i];
    }
    return res;
  }
  return 0.;
}

//__________________________________________________________________________________
/**
 * @brief Check if the first term is a number.
 * @return True if the first term is a number, false otherwise.
 */
bool _PolynomialData::IsFirstANumber(void) {
  return !ArrayAny(GetTerm(0), numberVars,
                   [](long p, unsigned long) { return p > 0; });
}

//__________________________________________________________________________________
/**
 * @brief Multiply two terms.
 * @param target The target term.
 * @param s1 The first source term.
 * @param s2 The second source term.
 */
void _PolynomialData::MultiplyTerms(long *target, long *s1, long *s2) {
  for (long i = 0; i < numberVars; i++) {
    target[i] = s1[i] + s2[i];
  }
}

//__________________________________________________________________________________

/**
 * @brief Raise a term to a power.
 * @param target The term to raise.
 * @param power The power.
 */
void _PolynomialData::RaiseTerm(long *target, long power) {
  for (long i = 0L; i < numberVars; i++) {
    target[i] *= power;
  }
}

//__________________________________________________________________________________
/**
 * @brief Raise a number to a power.
 * @param base The base.
 * @param pwr The power.
 * @return The result.
 */
hyFloat _PolynomialData::BinaryRaise(hyFloat base, long pwr) {
  if (pwr >= 0) {
    return ComputePower(base, pwr);
  }
  return 1. / ComputePower(base, -pwr);
}

/**
 * @brief Find a term using binary search.
 * @param theTerm The term to find.
 * @param reIndex The reindex array providing sorted access to terms.
 * @param start The start index for the search.
 * @return The index of the term if found. If not found, returns a negative
 * value encoding the insertion point based on a custom convention.
 */
long _PolynomialData::FindTerm(long *theTerm, long *reIndex, long start) const {
  long top = actTerms - 1;
  long bottom = start;

  if (top < bottom) { // No terms to search in the given range.
    return -bottom - 2;
  }

  char comp;

  while (bottom < top) {
    // Using `bottom + (top - bottom) / 2` prevents potential integer overflow.
    long middle = bottom + (top - bottom) / 2;
    comp = CompareTerms(GetTerm(reIndex[middle]), theTerm);

    if (comp < 0) {
      // The term at `middle` is smaller, so we search in the upper half.
      bottom = middle + 1;
    } else {
      if (comp == 0) {
        return middle;
      }
      // The term at `middle` is greater or equal, so we search in the lower
      // half, including `middle`.
      top = middle;
    }
  }

  // After the loop, `bottom` (which equals `top`) is the index of the first
  // element that is not less than `theTerm`. We check if it's an exact match.

  comp = CompareTerms(GetTerm(reIndex[bottom]), theTerm);
  if (comp == 0) {
    return bottom; // Exact match found.
  }

  // The term was not found.
  // return the inserion point
  return comp < 0 ? -bottom - 3 : -bottom - 2;
}

//__________________________________________________________________________________
void _PolynomialData::ResortTerms(long *reIndex) {
  long *source, *target, deleted = 0;

  hyFloat *newCoeff = (hyFloat *)MemAllocate(allocTerms * sizeof(hyFloat));
  long *newPowers = (long *)MemAllocate(allocTerms * numberVars * sizeof(long));

  // pass 1 to check for deletions
  for (long i = 0; i < actTerms; i++) {
    if (KeepTerm(theCoeff[reIndex[i]], reIndex[i])) {
      newCoeff[i] = theCoeff[reIndex[i]];
    } else {
      newCoeff[i] = 0.0;
    }
  }

  // pass 2 copy
  for (long i = 0; i < actTerms; i++) {
    if (newCoeff[i] != 0.0) {
      if (deleted) {
        newCoeff[i - deleted] = newCoeff[i];
      }
      target = newPowers + numberVars * (i - deleted);
      source = thePowers + numberVars * (reIndex[i]);
      for (long j = 0; j < numberVars; j++) {
        target[j] = source[j];
      }
    } else {
      deleted++;
    }
  }
  free(theCoeff);
  free(thePowers);
  theCoeff = newCoeff;
  thePowers = newPowers;
  actTerms -= deleted;
  ShrinkExtraSpace();
}
//__________________________________________________________________________________
/**
 * @brief Chop the terms.
 */
void _PolynomialData::ChopTerms(void) {
  long maxAllowedTerms = maximumPolyTermsPerVariable * numberVars;
  if (actTerms <= maxAllowedTerms) {
    return; // nothing to do
  } else {
    _SimpleList terms, index;
    long i, *qa, k, j;
    hyFloat logTop = log(topPolyCap), *qc;
    for (i = 0; i < actTerms; i++) {
      index << i;
      terms << (long)(log(fabs(theCoeff[i])) + logTop * SumOfPowers(i));
    }
    SortLists(&terms, &index);
    terms.Clear();
    qa = index.quickArrayAccess() + maxAllowedTerms;
    for (i = maxAllowedTerms; i < actTerms; i++) {
      theCoeff[qa[i]] = 0.0;
    }

    allocTerms =
        (maxAllowedTerms / POLY_DATA_INCREMENT + 1) * POLY_DATA_INCREMENT;
    hyFloat *newCoeff = (hyFloat *)MemAllocate(allocTerms * sizeof(hyFloat)),
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

//__________________________________________________________________________________
/**
 * @brief Rearrange a term.
 * @param target The target term.
 * @param source The source term.
 * @param markup The markup.
 * @param items The number of items.
 */
void _PolynomialData::RearrangeTerm(long *target, long *source, long *markup,
                                    long items) {
  for (long i = 0; i < items; i++) {
    target[markup[i]] = source[i];
  }
}
//__________________________________________________________________________________
/**
 * @brief Compare two terms.
 * @param s1 The first term.
 * @param s2 The second term.
 * @return 0 if the terms are equal, 1 if the first term is greater, -1 if the
 * second term is greater.
 */
char _PolynomialData::CompareTerms(long *s1, long *s2) const {
  for (long i = 0; i < numberVars; i++) {
    long comp = s1[i] - s2[i];
    if (comp > 0L) {
      return 1;
    }
    if (comp < 0L) {
      return -1;
    }
  }
  return 0;
}

//__________________________________________________________________________________
/**
 * @brief Compare two terms.
 * @param s1 The first term.
 * @param s2 The second term.
 * @param secondReindex The reindex for the second term.
 * @param actLength The active length.
 * @return 0 if the terms are equal, 1 if the first term is greater, -1 if the
 * second term is greater.
 */
char _PolynomialData::CompareTerms(long *s1, long *s2, long *secondReindex,
                                   long actLength) const {
  long second_index = 0;
  for (long k = 0L; k < numberVars; k++) {
    long v1 = s1[k], v2 = 0;
    if (second_index < actLength && k == secondReindex[second_index]) {
      v2 = s2[second_index];
      second_index++;
    }

    if (v1 != v2) {
      return v1 > v2 ? 1 : -1;
    }
  }

  return 0;
}

//__________________________________________________________________________________
/**
 * @brief Compare two terms.
 * @param s1 The first term.
 * @param s2 The second term.
 * @param firstReindex The reindex for the first term.
 * @param secondReindex The reindex for the second term.
 * @param actLength1 The active length of the first term.
 * @param actLength2 The active length of the second term.
 * @return 0 if the terms are equal, 1 if the first term is greater, -1 if the
 * second term is greater.
 */
char _PolynomialData::CompareTerms(long *s1, long *s2, long *firstReindex,
                                   long *secondReindex, long actLength1,
                                   long actLength2) const {
  long first_index = 0L, second_index = 0L, k;

  for (k = 0; k < numberVars; k++) {

    long v1 = 0, v2 = 0;

    if (first_index < actLength1) {
      if (k == firstReindex[first_index]) {
        v1 = s1[first_index];
        first_index++;
      }
    }

    if (second_index < actLength2) {
      if (k == secondReindex[second_index]) {
        v2 = s2[second_index];
        second_index++;
      }
    }

    if (v1 != v2) {
      return v1 < v2 ? -1 : 1;
    }
  }

  return 0;
}

//__________________________________________________________________________________
/**
 * @brief Check a term.
 * @param myCoeff The coefficient of the term.
 * @param myIndex The index of the term.
 * @return True if the term is valid, false otherwise.
 */

bool _PolynomialData::KeepTerm(hyFloat myCoeff, long myIndex) const {
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
              log(topPolyCap) * (hyFloat)(SumOfPowers(myIndex)) <
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

//__________________________________________________________________________________

_Polynomial::_Polynomial(void) { theTerms = new _PolynomialData; }

//__________________________________________________________________________________

_Polynomial::_Polynomial(_SimpleList &vars) {
  variableIndex.Duplicate(&vars);
  theTerms = new _PolynomialData(vars.countitems());
}

//__________________________________________________________________________________

_Polynomial::_Polynomial(hyFloat value) {
  // a constant polynomial
  theTerms = new _PolynomialData;
  theTerms->AddTerm(nil, value);
}

//__________________________________________________________________________________

_Polynomial::_Polynomial(_Variable &v)
// a monic monomial
{
  variableIndex << v.get_index();
  theTerms = new _PolynomialData(1);
  long vIndex = 1;
  theTerms->AddTerm(&vIndex, 1.0);
}

//__________________________________________________________________________________

_Polynomial::_Polynomial(_Polynomial const &p) { *this = p; }

//__________________________________________________________________________________

_Polynomial const &_Polynomial::operator=(_Polynomial const &p) {
  if (this != &p) {
    variableIndex.Duplicate(&p.variableIndex);
    theTerms = new _PolynomialData;
    if (p.theTerms) {
      theTerms->Duplicate(p.theTerms);
    } else {
      theTerms->numberVars = variableIndex.countitems();
    }
    compList1.Duplicate(&p.compList1);
    compList2.Duplicate(&p.compList2);
  }
  return *this;
}
//__________________________________________________________________________________

_Polynomial::~_Polynomial() {
  if (theTerms) {
    DeleteObject(theTerms);
  }
}

//__________________________________________________________________________________

BaseObj *_Polynomial::makeDynamic(void) const {
  _Polynomial *res = new _Polynomial;

  res->variableIndex.Duplicate(&variableIndex);
  res->compList1.Duplicate(&compList1);
  res->compList2.Duplicate(&compList2);

  if (theTerms) {
    res->theTerms->Duplicate(theTerms);
  } else {
    DeleteObject(res->theTerms);
    res->theTerms = nil;
  }
  return res;
}

//__________________________________________________________________________________

void _Polynomial::Duplicate(BaseRefConst tp) {
  _Polynomial const *p = (_Polynomial const *)tp;
  variableIndex.Clear();
  variableIndex.Duplicate(&p->variableIndex);
  compList1.Duplicate(&p->compList1);
  compList2.Duplicate(&p->compList2);
  DeleteAndZeroObject(theTerms);
  if (p->theTerms) {
    theTerms = new _PolynomialData(*(p->theTerms));
  }
}
//__________________________________________________________________________________

HBLObjectRef _Polynomial::ExecuteSingleOp(long opCode, _List *arguments,
                                          _hyExecutionContext *context,
                                          HBLObjectRef cache) {
  switch (opCode) {     // first check operations without arguments
  case HY_OP_CODE_TYPE: // Type
    return Type(cache);
  case HY_OP_CODE_EVAL: // Eval
    return _returnConstantOrUseCache(ComputePolynomial(), cache);
  case HY_OP_CODE_POLYNOMIAL: // Polynomial
    return new _Polynomial(*this);
  }

  _MathObject *arg0 = _extract_argument(arguments, 0UL, false);

  switch (
      opCode) { // next check operations without arguments or with one argument
  case HY_OP_CODE_SUB: // -
    if (arg0) {
      return Sub(arg0, cache);
    } else {
      return Minus(cache);
    }
  case HY_OP_CODE_ADD: // +
    if (arg0) {
      return Add(arg0, cache);
    } else {
      return Sum(cache);
    }
  }

  if (arg0) {
    switch (opCode) {
    case HY_OP_CODE_MUL: //*
      return Mult(arg0, cache);
    case HY_OP_CODE_POWER: // ^
      return Raise(arg0, cache);
    case HY_OP_CODE_EQ: // ==
      return _returnConstantOrUseCache(Equal(arg0), cache);
    case HY_OP_CODE_NEQ: // !=
      return _returnConstantOrUseCache(!Equal(arg0), cache);
    }
  }

  switch (opCode) {
  case HY_OP_CODE_MUL:
  case HY_OP_CODE_POWER:
  case HY_OP_CODE_EQ:
  case HY_OP_CODE_NEQ:
    WarnWrongNumberOfArguments(this, opCode, context, arguments);
    break;
  default:
    WarnNotDefined(this, opCode, context);
  }

  return new _MathObject;
}

//__________________________________________________________________________________
_MathObject *_Polynomial::IsANumber(bool returnLeading) {
  long nV = variableIndex.countitems();
  if (!nV) {
    if (theTerms->NumberOfTerms() > 0) {
      return new _Constant(theTerms->theCoeff[0]);
    } else {
      return new _Constant(0.0);
    }
  }
  if (theTerms->NumberOfTerms() == 0) {
    return new _Constant(0.0);
  }
  if (theTerms->NumberOfTerms() == 1) {
    if (theTerms->IsFirstANumber() || returnLeading) {
      return new _Constant(theTerms->theCoeff[0]);
    }
  }
  return nil;
}

//__________________________________________________________________________________
_MathObject *_Polynomial::Add(_MathObject *m, HBLObjectRef cache) {
  return Plus(m, false, cache);
}
//__________________________________________________________________________________

void ResetPolynomialCheck(_Polynomial *p) {
  if (dropTerms) {
    if (!enforcePolyCap) {
      if (varCheckAllocated != (long)p->variableIndex.countitems()) {
        if (varCheckArray) {
          free(varCheckArray);
        }
        varCheckAllocated = p->variableIndex.countitems();
        varCheckArray =
            (hyFloat *)MemAllocate(varCheckAllocated * sizeof(hyFloat));
        hyFloat lb, ub;
        for (long j = varCheckAllocated - 1; j >= 0; j--) {
          _Variable *theV = LocateVar(p->variableIndex(j));
          lb = fabs(theV->GetLowerBound());
          ub = fabs(theV->GetUpperBound());
          varCheckArray[j] = log(lb > ub ? lb : ub);
        }
      }
    }
  }
  checkReset = true;
}

//__________________________________________________________________________________
void _Polynomial::CheckTerm(void) {
  long myIndex = theTerms->actTerms - 1;
  hyFloat myCoeff = theTerms->GetCoeff(myIndex);
  if (!theTerms->KeepTerm(myCoeff, myIndex)) {
    theTerms->DeleteTerm(myIndex);
  }
}

//__________________________________________________________________________________

bool _Polynomial::Equal(_MathObject *m) {
  bool result = false;
  if (m->ObjectClass() == POLYNOMIAL || m->ObjectClass() == NUMBER) {
    _Polynomial *diff = (_Polynomial *)Sub(m, nil);
    if (diff) {
      _Constant *v = (_Constant *)diff->IsANumber(true);
      if (v) {
        result = CheckEqual(v->Value(), 0.0);
        DeleteObject(v);
      }
      DeleteObject(diff);
    }
  }
  return result;
}
//__________________________________________________________________________________
void _Polynomial::DropZeroCoefficient(void) {
  if (theTerms->GetCoeff(0) == 0.0) {
    theTerms->DeleteTerm(0);
  }
}

//__________________________________________________________________________________
_Polynomial *_Polynomial::NegateCoefficients(void) {
  hyFloat *invC = theTerms->theCoeff;
  for (long inv = 0; inv < theTerms->actTerms; inv++) {
    invC[inv] = -invC[inv];
  }
  return this;
}

//__________________________________________________________________________________
_MathObject *_Polynomial::Plus(_MathObject *m, bool subtract, HBLObjectRef) {
  long objectT = m->ObjectClass();

  if (objectT == NUMBER) { // a number
    Convert2OperationForm();
    _Polynomial *result = new _Polynomial(*this);
    hyFloat mV = subtract ? -m->Value() : m->Value();

    if (variableIndex.empty()) { // constant poly
      if (theTerms->NumberOfTerms()) {
        result->theTerms->GetCoeff(0) += mV;
        result->DropZeroCoefficient();
      } else {
        if (mV != 0.0) {
          result->theTerms->AddTerm(nil, mV);
        }
      }
    } else { // variable containing poly
      if (theTerms->NumberOfTerms()) {
        if (theTerms->IsFirstANumber()) {
          result->theTerms->GetCoeff(0) += mV;
          result->DropZeroCoefficient();
        } else {
          result->theTerms->AddTerm(mV);
        }
      } else {
        result->theTerms->AddTerm(mV);
      }
    }
    return result;
  }

  if (objectT == POLYNOMIAL) { // another polynomial

    Convert2OperationForm();
    _Polynomial *p2 = (_Polynomial *)m;
    p2->Convert2OperationForm();
    if (variableIndex.empty()) {       // this polynomial has no variables
      if (theTerms->NumberOfTerms()) { // this polynomial is just a number
        _Constant coef1(theTerms->GetCoeff(0));
        if (!subtract) {
          return p2->Plus(&coef1, false, nil);
        } else {
          return ((_Polynomial *)p2->Plus(&coef1, true, nil))
              ->NegateCoefficients();
        }
      } else {
        if (!subtract) { // this polynomial is 0
          return new _Polynomial(*p2);
        } else {
          return (new _Polynomial(*p2))->NegateCoefficients();
        }
      }
    }

    p2->Convert2OperationForm();

    if (p2->variableIndex.empty()) {       // the second term is a constant
      if (p2->theTerms->NumberOfTerms()) { // the second term is NOT empty
        _Constant coef2(p2->theTerms->GetCoeff(0));
        return Plus(&coef2, subtract, nil);
      } else {
        return new _Polynomial(*this);
      }
    }

    long nt2 = p2->theTerms->NumberOfTerms(), nt1 = theTerms->NumberOfTerms(),
         pos1 = 0, pos2 = 0, *term1, *term2;

    hyFloat *coeff1 = theTerms->GetCoeff(), *coeff2 = p2->theTerms->GetCoeff();

    _Polynomial *res;

    _SimpleList mergedVariables, merge1, merge2;
    mergedVariables.Merge(variableIndex, p2->variableIndex, &merge1, &merge2);

    res = new _Polynomial(mergedVariables); // create a blank new result holder
    ResetPolynomialCheck(res);

    while (pos1 < nt1 && pos2 < nt2) {
      term1 = theTerms->GetTerm(pos1);
      term2 = p2->theTerms->GetTerm(pos2);
      char comp = res->theTerms->CompareTerms(
          term1, term2, merge1.quickArrayAccess(), merge2.quickArrayAccess(),
          merge1.countitems(), merge2.countitems());
      if (comp < 0) {
        res->theTerms->AddTerm(term1, coeff1[pos1], merge1.quickArrayAccess(),
                               merge1.countitems());
        pos1++;
      } else if (comp > 0) {
        res->theTerms->AddTerm(term2, subtract ? -coeff2[pos2] : coeff2[pos2],
                               merge2.quickArrayAccess(), merge2.countitems());
        pos2++;
      } else {
        res->theTerms->AddTerm(
            term1, coeff1[pos1] + (subtract ? -coeff2[pos2] : coeff2[pos2]),
            merge1.quickArrayAccess(), merge1.countitems());
        pos1++;
        pos2++;
      }
      res->CheckTerm();
    }

    while (pos1 < nt1) {
      term1 = theTerms->GetTerm(pos1);
      res->theTerms->AddTerm(term1, coeff1[pos1++], merge1.quickArrayAccess(),
                             merge1.countitems());
      res->CheckTerm();
    }

    while (pos2 < nt2) {
      term2 = p2->theTerms->GetTerm(pos2);
      if (subtract) {
        res->theTerms->AddTerm(term2, -coeff2[pos2++],
                               merge2.quickArrayAccess(), merge2.countitems());
      } else {
        res->theTerms->AddTerm(term2, coeff2[pos2++], merge2.quickArrayAccess(),
                               merge2.countitems());
      }
      res->CheckTerm();
    }

    if (!res->theTerms->checkMe()) {
      return nil;
    }
    //      res->theTerms->ChopTerms();
    if (res->theTerms->GetNoTerms() == 0) {
      DeleteObject(res);
      return new _Polynomial(0.0);
    }
    return res;
  }

  HandleApplicationError("An incompatible operand was supplied to polynomial "
                         "addition/subtraction");
  return nil;
}

//__________________________________________________________________________________
_MathObject *_Polynomial::Sub(_MathObject *m, HBLObjectRef cache) {
  return Plus(m, true, cache);
}

//__________________________________________________________________________________
_MathObject *_Polynomial::Minus(HBLObjectRef) {
  return (new _Polynomial(*this))->NegateCoefficients();
}

//__________________________________________________________________________________
_MathObject *_Polynomial::Mult(_MathObject *m, HBLObjectRef) {
  long objectT = m->ObjectClass();

  if (objectT == NUMBER) { // a number or a monomial
    Convert2OperationForm();
    hyFloat nb = ((_Constant *)m)->Value();
    if (nb == 0.0) {
      return new _Polynomial(0.);
    }

    _Polynomial *result = new _Polynomial(*this);
    for (long i = theTerms->NumberOfTerms() - 1; i >= 0; i--) {
      result->theTerms->GetCoeff(i) *= nb;
    }
    return result;
  }

  if (objectT == POLYNOMIAL) { // another polynomial
    Convert2OperationForm();
    _Polynomial *p2 = (_Polynomial *)m;
    p2->Convert2OperationForm();
    if (!variableIndex.countitems()) {
      if (theTerms->NumberOfTerms()) {
        _Constant coef1(theTerms->GetCoeff(0));
        return p2->Mult(&coef1, nil);
      } else {
        return new _Polynomial();
      }
    }
    if (!p2->variableIndex.countitems()) {
      if (p2->theTerms->NumberOfTerms()) {
        _Constant coef2(p2->theTerms->GetCoeff(0));
        return Mult(&coef2, nil);
      } else {
        return new _Polynomial();
      }
    }
    long nt2 = p2->theTerms->NumberOfTerms(), nt1 = theTerms->NumberOfTerms(),
         ref = (variableIndex.countitems() + p2->variableIndex.countitems()) *
               maximumPolyTermsPerVariable,
         *idx1;
    _Polynomial *res;
    _SimpleList reIndex, terms1, terms2;

    if (nt1 * nt2 > ref)
    // too many terms expected - chop b4 multiplying
    {
      _SimpleList index;
      long i;
      hyFloat logTop = log(topPolyCap), *theCoeff = theTerms->GetCoeff();
      for (i = 0; i < nt1; i++, theCoeff++) {
        index << i;
        terms1 << (long)(log(fabs(*theCoeff)) +
                         logTop * theTerms->SumOfPowers(i));
      }
      theCoeff = p2->theTerms->GetCoeff();
      for (i = 0; i < nt2; i++, theCoeff++) {
        index << i + nt1;
        terms1 << (long)(log(fabs(*theCoeff)) +
                         logTop * p2->theTerms->SumOfPowers(i));
      }
      SortLists(&terms1, &index);
      terms1.Clear();
      terms2.Clear();
      idx1 = index.quickArrayAccess() + nt1 + nt2 - 1;

      for (i = 0; i < ref; i++, idx1--) {
        if (*idx1 >= nt1) {
          terms2 << *idx1;
        } else {
          terms1 << *idx1;
        }
      }
      if (terms1.countitems() == 0) {
        terms1 << 0;
      }
      if (terms2.countitems() == 0) {
        terms2 << 0;
      }
    }

    if (!terms1.countitems()) {
      if (variableIndex.Equal(
              p2->variableIndex)) { // same variable arrays - proceed to the
                                    // operation directly
        res = new _Polynomial(variableIndex);
        ResetPolynomialCheck(res);
        long *newTerm = new long[res->variableIndex.countitems()], *term1, f;
        for (long i = 0; i < nt1; i++) {
          ref = 0;
          hyFloat c1 = theTerms->GetCoeff(i);
          term1 = theTerms->GetTerm(i);
          for (long j = 0; j < nt2; j++) {
            res->theTerms->MultiplyTerms(newTerm, term1,
                                         p2->theTerms->GetTerm(j));
            f = res->theTerms->FindTerm(newTerm, reIndex.quickArrayAccess(),
                                        ref);
            if (f >= 0) { // existing term - add or sub
              res->theTerms->GetCoeff(reIndex(f)) +=
                  c1 * p2->theTerms->GetCoeff(j);
              ref = f + 1;
            } else { // new term - insert
              reIndex.InsertElement((BaseRef)res->theTerms->actTerms, -f - 2,
                                    false, false);
              res->theTerms->AddTerm(newTerm, c1 * p2->theTerms->GetCoeff(j));
              ref = -f - 1;
            }
          }
        }

        delete[] newTerm;
      } else {
        _SimpleList merge1, merge2, joint;
        joint.Merge(variableIndex, p2->variableIndex, &merge1, &merge2);
        res = new _Polynomial(joint);
        ResetPolynomialCheck(res);
        long f, nv = res->variableIndex.countitems(), *scratch1 = new long[nv],
                *scratch3, k, *scratch2;
        bool do1 = merge1.countitems() != 0, do2 = merge2.countitems() != 0;
        if (do1) {
          scratch2 = new long[nv];
        }
        if (do2) {
          scratch3 = new long[nv];
        }
        for (long i = 0; i < nt1; i++) {
          hyFloat c1 = theTerms->GetCoeff(i);
          ref = 0;
          if (do1) {
            for (k = 0; k < nv; k++, scratch2++) {
              *scratch2 = 0;
            }
            scratch2 -= nv;
            theTerms->RearrangeTerm(scratch2, theTerms->GetTerm(i),
                                    merge1.quickArrayAccess(),
                                    merge1.countitems());
          } else {
            scratch2 = theTerms->GetTerm(i);
          }
          for (long j = 0; j < nt2; j++) {
            if (do2) {
              for (k = 0; k < nv; k++, scratch3++) {
                *scratch3 = 0;
              }
              scratch3 -= nv;
              theTerms->RearrangeTerm(scratch3, p2->theTerms->GetTerm(j),
                                      merge2.quickArrayAccess(),
                                      merge2.countitems());
            } else {
              scratch3 = p2->theTerms->GetTerm(j);
            }
            res->theTerms->MultiplyTerms(scratch1, scratch2, scratch3);
            f = res->theTerms->FindTerm(scratch1, reIndex.quickArrayAccess(),
                                        ref);
            if (f >= 0) { // existing term - add or sub
              res->theTerms->GetCoeff(reIndex(f)) +=
                  c1 * p2->theTerms->GetCoeff(j);
              ref = f + 1;
            } else { // new term - insert
              reIndex.InsertElement((BaseRef)res->theTerms->actTerms, -f - 2,
                                    false, false);
              res->theTerms->AddTerm(scratch1, c1 * p2->theTerms->GetCoeff(j));
              ref = -f - 1;
            }
          }
        }
        delete[] scratch1;
        if (do1) {
          delete[] scratch2;
        }
        if (do2) {
          delete[] scratch3;
        }
      }

    } else {
      terms1.Sort();
      terms2.Sort();
      long *onOff = (long *)MemAllocate(nt2 * sizeof(long));
      onOff += (nt2 - 1);
      for (long i = nt2 - 1; i >= 0; i--, onOff--)
        if (terms2.BinaryFind(i) >= 0) {
          *onOff = 1;
        } else {
          *onOff = 0;
        }
      onOff++;
      if (variableIndex.Equal(
              p2->variableIndex)) { // same variable arrays - proceed to the
                                    // operation directly
        res = new _Polynomial(variableIndex);
        ResetPolynomialCheck(res);
        long *newTerm = new long[res->variableIndex.countitems()], *term1, f;
        for (long i = 0; i < nt1; i++) {
          if (terms1.BinaryFind(i) < 0) {
            continue;
          }
          ref = 0;
          hyFloat c1 = theTerms->GetCoeff(i);
          term1 = theTerms->GetTerm(i);
          for (long j = 0; j < nt2; j++) {
            if (onOff[j] == 0) {
              continue;
            }
            res->theTerms->MultiplyTerms(newTerm, term1,
                                         p2->theTerms->GetTerm(j));
            f = res->theTerms->FindTerm(newTerm, reIndex.quickArrayAccess(),
                                        ref);
            if (f >= 0) { // existing term - add or sub
              res->theTerms->GetCoeff(reIndex(f)) +=
                  c1 * p2->theTerms->GetCoeff(j);
              ref = f + 1;
            } else { // new term - insert
              reIndex.InsertElement((BaseRef)res->theTerms->actTerms, -f - 2,
                                    false, false);
              res->theTerms->AddTerm(newTerm, c1 * p2->theTerms->GetCoeff(j));
              ref = -f - 1;
            }
          }
        }

        delete[] newTerm;
      } else {
        _SimpleList merge1, merge2, joint;
        joint.Merge(variableIndex, p2->variableIndex, &merge1, &merge2);
        res = new _Polynomial(joint);
        ResetPolynomialCheck(res);
        long f, nv = res->variableIndex.countitems(), *scratch1 = new long[nv],
                *scratch3, k, *scratch2;
        bool do1 = merge1.countitems() != 0, do2 = merge2.countitems() != 0;
        if (do1) {
          scratch2 = new long[nv];
        }
        if (do2) {
          scratch3 = new long[nv];
        }
        for (long i = 0; i < nt1; i++) {
          hyFloat c1 = theTerms->GetCoeff(i);
          ref = 0;
          if (do1) {
            for (k = 0; k < nv; k++, scratch2++) {
              *scratch2 = 0;
            }
            scratch2 -= nv;
            theTerms->RearrangeTerm(scratch2, theTerms->GetTerm(i),
                                    merge1.quickArrayAccess(),
                                    merge1.countitems());
          } else {
            scratch2 = theTerms->GetTerm(i);
          }
          for (long j = 0; j < nt2; j++) {
            if (do2) {
              for (k = 0; k < nv; k++, scratch3++) {
                *scratch3 = 0;
              }
              scratch3 -= nv;
              theTerms->RearrangeTerm(scratch3, p2->theTerms->GetTerm(j),
                                      merge2.quickArrayAccess(),
                                      merge2.countitems());
            } else {
              scratch3 = p2->theTerms->GetTerm(j);
            }
            res->theTerms->MultiplyTerms(scratch1, scratch2, scratch3);
            f = res->theTerms->FindTerm(scratch1, reIndex.quickArrayAccess(),
                                        ref);
            if (f >= 0) { // existing term - add or sub
              res->theTerms->GetCoeff(reIndex(f)) +=
                  c1 * p2->theTerms->GetCoeff(j);
              ref = f + 1;
            } else { // new term - insert
              reIndex.InsertElement((BaseRef)res->theTerms->actTerms, -f - 2,
                                    false, false);
              res->theTerms->AddTerm(scratch1, c1 * p2->theTerms->GetCoeff(j));
              ref = -f - 1;
            }
          }
        }
        delete[] scratch1;
        if (do1) {
          delete[] scratch2;
        }
        if (do2) {
          delete[] scratch3;
        }
      }
      free(onOff);
    }

    res->theTerms->ResortTerms(reIndex.quickArrayAccess());
    res->theTerms->checkMe();
    //      res->theTerms->ChopTerms();
    return res;
  }

  HandleApplicationError(
      "An incompatible operand was supplied to polynomial multiplication");
  return nil;
}

//__________________________________________________________________________________
_MathObject *_Polynomial::Raise(_MathObject *m, HBLObjectRef) {
  long objectT = m->ObjectClass();
  bool del = false;
  if (objectT == POLYNOMIAL) {
    m = ((_Polynomial *)m)->IsANumber();
    if (!m) {
      return nil;
    }
    del = true;
    objectT = m->ObjectClass();
  }
  if (objectT == NUMBER) { // a number
    Convert2OperationForm();
    _Polynomial *result;
    if (theTerms->NumberOfTerms() == 1) { // just a monomial
      long power = (long)m->Value();
      result = new _Polynomial(*this);
      result->theTerms->RaiseTerm(result->theTerms->GetTerm(0), power);
      result->theTerms->GetCoeff(0) =
          _PolynomialData::BinaryRaise(result->theTerms->GetCoeff(0), power);

    } else { // binary raise
      result = new _Polynomial(1.0);

      _Polynomial *oldR;
      unsigned char bits[sizeof(long) * 8];
      long pwr = (long)m->Value(), nLength = 0L;
      while (pwr) {
        bits[nLength] = pwr % 2;
        pwr /= 2;
        nLength++;
      }
      while (nLength) {
        oldR = result;
        nLength--;
        if (bits[nLength]) {
          result = (_Polynomial *)result->Mult(this, nil);
          DeleteObject(oldR);
        }
        oldR = result;
        if (nLength) {
          result = (_Polynomial *)result->Mult(result, nil);
          DeleteObject(oldR);
        }
      }
    }
    if (del) {
      DeleteObject(m);
    }
    return result;
  }
  HandleApplicationError(
      "An incompatible operand was supplied to polynomial raise to power");
  return nil;
}
//__________________________________________________________________________________
void _Polynomial::Convert2ComputationForm(_SimpleList *c1, _SimpleList *c2,
                                          _SimpleList *termsToInclude) {
  // compList1 and two one will contain pairs of numbers to be interpreted as
  // follows for pair (m,n), where m is in compList1 and n is compList2
  // m>=0, n>=0, m!=last var
  // m<0, n<0 - mult by -n-th power of var -m-1, do NOT add yet
  // m<0, n>0 - mult by -n-th power of var -m-1, DO add
  // m>=0, n<0 - reset all vars after m and set the power of m-the var to -n, do
  // not compute yet m>=0, n>0 - reset all vars after m and set the power of
  // m-the var to -n, do add last var, n>0 - add n successive powers of last var

  if ((theTerms->NumberOfTerms() != 0) &&
      (compList1.countitems() == 0)) { // stuff to do
    _SimpleList *op_var_index, *op_power_code, ti;

    const long variable_count = variableIndex.countitems();
    long i = 0, n = variable_count - 1, count_of_items_to_process,
         accumulator_counter = 0;

    if (c1 && c2 && termsToInclude) {
      op_var_index = c1;
      op_power_code = c2;
      ti.Duplicate(termsToInclude);
    } else {
      op_var_index = &compList1;
      op_power_code = &compList2;
      ti.Populate(theTerms->actTerms, 0, 1);
      theTerms->actTerms = 0;
    }

    count_of_items_to_process = ti.countitems();
    // vi = ti.quickArrayAccess();
    op_var_index->Clear();
    op_power_code->Clear();

    auto _TERM_add_coefficient_to_result = [&]() -> void {
      (*op_var_index) << n;
      (*op_power_code) << 0;
    };

    auto _TERM_last_var_jump = [&](long power) -> void {
      (*op_var_index) << n;
      (*op_power_code) << -power;
    };

    auto _TERM_raise_var_and_add = [&](long var_index, long power) -> void {
      (*op_var_index) << -var_index - 1;
      (*op_power_code) << power;
    };

    auto _TERM_raise_var_deferred = [&](long var_index, long power) -> void {
      if (power) {
        (*op_var_index) << -var_index - 1;
        (*op_power_code) << -power;
      }
    };

    auto _TERM_raise_var_deferred_and_reset = [&](long var_index,
                                                  long power) -> void {
      (*op_var_index) << var_index;
      (*op_power_code) << -power;
    };

    auto _TERM_raise_variable_deferred = [&](long var_index,
                                             long power) -> void {
      if (power) {
        (*op_var_index) << -var_index - 1;
        (*op_power_code) << -power;
      }
    };

    auto _TERM_last_var_succession = [&](long count) -> void {
      (*op_var_index) << n;
      (*op_power_code) << count;
    };

    auto _TERM_remove_deferral = [&](void) -> void {
      (*op_power_code)[op_power_code->countitems() - 1] =
          -op_power_code->get(op_power_code->countitems() - 1);
    };

    if (!theTerms->IsFirstANumber()) {
      long *first_term = theTerms->GetTerm(ti.get(0));
      // set up the first term
      // check later why the ugly add then delete
      _TERM_raise_variable_deferred(n, first_term[n]);

      for (long j = variable_count - 2; j >= 0; j--) {
        long power = first_term[j];
        if (power) {
          _TERM_raise_variable_deferred(j, power);
        }
      }

      _TERM_remove_deferral();
      //(*op_power_code)[op_power_code->countitems() - 1] = -op_power_code->get
      //(op_power_code->countitems() - 1);
      // remove deferral; add this term

      if (op_power_code->countitems() > 1) {
        if (first_term[n] == 0) {
          op_var_index->Delete(0);
          op_power_code->Delete(0);
        }
      }
    } else {
      _TERM_add_coefficient_to_result();
    }
    i++;

    long *powerDiff = (long *)alloca(sizeof(long) * (n + 1)); // (why n+2)?

    for (; i < count_of_items_to_process; i++) { // even more stuff to do

      long *current_term_powers = theTerms->GetTerm(ti.get(i)),
           *previous_term_powers = theTerms->GetTerm(ti.get(i - 1));

      long k = 0, ch = -1;
      bool reset = false;

      for (long j = 0; j < n; j++) {
        powerDiff[j] = current_term_powers[j] - previous_term_powers[j];
        if (powerDiff[j]) {
          if (ch < 0) {
            ch = j;
          }
          k--;
          reset = reset || (powerDiff[n] < 0);
        }
      }

      powerDiff[n] = current_term_powers[n] - previous_term_powers[n];
      reset = reset || (powerDiff[n] < 0);

      if (k == 0) {
        k = powerDiff[n];
      } else if (powerDiff[n]) {
        k--;
      }
      // analyze the difference
      if (k == 1) {
        accumulator_counter++; // the only difference is a +1 power in the last
                               // variable
        continue;
      } else {
        if (accumulator_counter > 0) {
          // counting the number of x + x^2 + x^3 type tems
          _TERM_last_var_succession(accumulator_counter);
          accumulator_counter = 0;
        }
        if (k > 0) { // only changed the last variable, but the power step > 1
          _TERM_last_var_jump(k);
          continue;
        }
        if (k < 0) {
          // power decreased from the previous term
          if (k == -1) { // change in exactly one variable
            _TERM_raise_var_and_add(ch, powerDiff[ch]);
            continue;
          }
          // change in more than one variable

          if (reset) {
            _TERM_raise_var_deferred_and_reset(ch, powerDiff[ch]);
          } else {
            _TERM_raise_var_deferred(ch, powerDiff[ch]);
          }

          for (long c = ch + 1; c <= n; c++) {
            if (powerDiff[c] > 0) {
              _TERM_raise_var_deferred(c, reset ? previous_term_powers[c] +
                                                      powerDiff[c]
                                                : powerDiff[c]);
              //(*op_var_index) << -c - 1;
              //(*op_power_code) << (reset ? -(previous_term_powers[c] +
              // powerDiff[c]) : -powerDiff[c]);
            } else {
              _TERM_raise_var_deferred(c,
                                       previous_term_powers[c] + powerDiff[c]);
            }
          }
          _TERM_remove_deferral();
        }
      }
    }
    if (accumulator_counter > 0) {
      _TERM_last_var_succession(accumulator_counter);
      //(*op_var_index) << n;
      //(*op_power_code) << accumulator_counter;
    }

    if (!(c1 && c2)) {
      free(theTerms->thePowers);
      theTerms->thePowers = nil;
    }
  }
}

//__________________________________________________________________________________
void _Polynomial::Convert2OperationForm(void) {
  // see if anything needs to be done at all
  if (compList1.nonempty() && !theTerms->thePowers) {
    long n = variableIndex.countitems() + 1, m = compList1.countitems();
    long lp = n - 2, *scratch = nil, index = 0;
    if (n > 1) {
      theTerms->thePowers =
          (long *)MemAllocate(theTerms->allocTerms * sizeof(long) * (n - 1));
      scratch = (long *)alloca(sizeof(long) * (n - 1));
      memset(scratch, 0, sizeof(long) * (n - 1));
      memset(theTerms->thePowers, 0,
             theTerms->allocTerms * sizeof(long) * (n - 1));
    }

    for (long i = 0; i < m; i++) { // loop over all commands
      long i1 = compList1.get(i);
      long i2 = compList2.get(i);
      if (i1 == lp) { // operations with the last var
        if (i2 > 0) { // do several iterations
          for (long k = 0; k < i2; k++, index++) {
            scratch[lp]++;
            theTerms->WriteTerm(scratch, index);
          }
        } else {
          if (!i2) {
            theTerms->WriteTerm(scratch, index);
          } else {
            scratch[lp] -= i2;
            theTerms->WriteTerm(scratch, index);
          }
          index++;
        }
      } else {
        if (i1 < 0) {
          bool compute = i2 < 0;
          i1 = -i1 - 1;
          if (i2 < 0) {
            i2 = -i2;
          }
          if (i2 == 1) {
            scratch[i1]++;
          } else {
            scratch[i1] += i2;
          }
          if (compute) {
            continue;
          }
          theTerms->WriteTerm(scratch, index);
          index++;
        } else {
          bool compute = i2 < 0;
          for (long k = i1 + 1; k <= lp; k++) {
            scratch[k] = 0;
          }
          if (i2 < 0) {
            i2 = -i2;
          }
          if (i2 == 1) {
            scratch[i1]++;
          } else {
            scratch[i1] += i2;
          }

          if (compute) {
            continue;
          }
          theTerms->WriteTerm(scratch, index);
          index++;
        }
      }
    }
    compList1.Clear();
    compList2.Clear();
  }
}
//__________________________________________________________________________________
_MathObject *_Polynomial::Compute(void) {
  return this; // new _Constant(ComputePolynomial());
}

//__________________________________________________________________________________
hyFloat _Polynomial::ComputePolynomial(void)
// assumed that the poly is already in the comp form
{
  Convert2ComputationForm();
  long n = variableIndex.countitems() + 1;
  hyFloat *varValues = (hyFloat *)alloca(n * sizeof(hyFloat));
  for (long i = 0; i < n - 1; i++) {
    varValues[i] = LocateVar(variableIndex(i))->Compute()->Value();
  }
  hyFloat result =
      ComputeP(varValues, theTerms->GetCoeff(), n, compList1.countitems(),
               compList1.quickArrayAccess(), compList2.quickArrayAccess());
  return result;
}
//__________________________________________________________________________________

hyFloat _Polynomial::ComputeP(hyFloat *varValues, hyFloat *compCoeff,
                              long variable_count, long array_length, long *c1,
                              long *c2) {
  hyFloat *holder = (hyFloat *)alloca(sizeof(hyFloat) * variable_count),
          term = 1., result = 0., lv;

  long lp = variable_count -
            2, // array index of the last variable, variable_count is one extra
      current_coefficient = 0;

  for (long var_index = 0; var_index < variable_count - 1; var_index++) {
    holder[var_index] = 1.;
  }

  lv = variable_count > 1 ? varValues[variable_count - 2] : 1.;

  /*
        Opcode meanings
        LVI = index of the last variable (varcount - 1)
        VV (i) : value of variable indexed i
        TERM   : current term
        RESULT : overall result
        MADD   : multiply by the current coefficient, then add to the RESULT,
     advacing the current coefficient HOLDER : an array of N powers (one for
     each variable + 1)

        opcode1 == LVI
            opcode2 == 0 : add the current coefficient to the result
            opcode2 < 0  : multiply TERM by VV (LVI) ^ (-opcode2 - 1), then
     current coefficient, and MADD this to the RESULT opcode2 > 0  : repeat
     opcode2 times multiply TERM by VV (LVI) the MADD to the result opcode1 < 0
            CVI (current variable index) -opcode1 - 1
            POWER = Abs (opcode2)

                multiply HOLDER (CVI) by VV (CVI) ^ POWER
                multiply TERM by VV (CVI) ^ POWER

            opcode2 > 0:
                MADD TERM to RESULT

        opcode1 > 0
            POWER = Abs (opcode2)
               Reset values opcode + 1 to the end in HOLDER to 1
               multiply HOLDER (opcode1) by VV (opcode1) ^ POWER

            opcode2 < 0
                MADD TERM TO RESULT

   */

  for (long array_index = 0; array_index < array_length;
       array_index++) { // loop over all commands
    long opcode1 = c1[array_index];
    long opcode2 = c2[array_index];
    if (opcode1 == lp) { // operations with the last var
      if (opcode2 > 0) { // do several iterations
        for (long k = 0; k < opcode2; k++, current_coefficient++) {
          term *= lv;
          result += term * compCoeff[current_coefficient];
        }
      } else {
        if (!opcode2) {
          result += compCoeff[current_coefficient];
        } else {
          term *= _PolynomialData::BinaryRaise(lv, -opcode2);
          result += term * compCoeff[current_coefficient];
        }
        current_coefficient++;
      }
    } else {
      if (opcode1 < 0) {
        bool compute = opcode2 < 0;
        opcode1 = -opcode1 - 1;
        if (compute) {
          opcode2 = -opcode2;
        }
        if (opcode2 == 1) {
          holder[opcode1] *= varValues[opcode1];
          term *= varValues[opcode1];
        } else {
          hyFloat p2 =
              _PolynomialData::BinaryRaise(varValues[opcode1], opcode2);
          holder[opcode1] *= p2;
          term *= p2;
        }
        if (compute) {
          continue;
        }
        result += term * compCoeff[current_coefficient++];
      } else { // opcode1 < 0
        bool compute = opcode2 < 0;
        for (long k = opcode1 + 1; k <= lp; k++) {
          holder[k] = 1.;
        }
        if (opcode2 < 0) {
          opcode2 = -opcode2;
        }
        if (opcode2 == 1) {
          holder[opcode1] *= varValues[opcode1];
        } else {
          holder[opcode1] *=
              _PolynomialData::BinaryRaise(varValues[opcode1], opcode2);
        }
        term = 1.;
        for (long k = 0; k <= opcode1; k++) {
          term *= holder[k];
        }

        if (compute) {
          continue;
        }
        result += term * compCoeff[current_coefficient++];
      }
    }
  }
  return result;
}
//__________________________________________________________________________________

void _Polynomial::RankTerms(_SimpleList *receptacle) {
  receptacle->Clear();
  hyFloat logTop = log(topPolyCap), *theCoeff = theTerms->GetCoeff();

  for (long i = 0; i < theTerms->actTerms; i++, theCoeff++) {
    (*receptacle) << (long)(log(fabs(*theCoeff)) +
                            logTop * theTerms->SumOfPowers(i));
  }
}

//__________________________________________________________________________________

BaseObj *_Polynomial::toStr(unsigned long padding) {
  _StringBuffer *result = new _StringBuffer(32UL);

  /*
        Opcode meanings
        LVI = index of the last variable (varcount - 1)
        VV (i) : value of variable indexed i
        TERM   : current term
        RESULT : overall result
        MADD   : multiply by the current coefficient, then add to the RESULT,
     advacing the current coefficient HOLDER : an array of N powers (one for
     each variable + 1)

        opcode1 == LVI
            opcode2 == 0 : add the current coefficient to the result
            opcode2 < 0  : multiply TERM by VV (LVI) ^ (-opcode2 - 1), then
     current coefficient, and MADD this to the RESULT opcode2 > 0  : repeat
     opcode2 times multiply TERM by VV (LVI) the MADD to the result opcode1 < 0
            CVI (current variable index) -opcode1 - 1
            POWER = Abs (opcode2)

                multiply HOLDER (CVI) by VV (CVI) ^ POWER
                multiply TERM by VV (CVI) ^ POWER

            opcode2 > 0:
                MADD TERM to RESULT

        opcode1 > 0
            POWER = Abs (opcode2)
               Reset values opcode + 1 to the end in HOLDER to 1
               multiply HOLDER (opcode1) by VV (opcode1) ^ POWER

            opcode2 < 0
                MADD TERM TO RESULT

   */

  if (compList1.nonempty() && !theTerms->thePowers) {
    long N = variableIndex.countitems() - 1;

    for (long i = 0; i <= N; i++) {
      (*result) << "Variable  " << _String(i) << " = "
                << LocateVar(variableIndex.get(i))->GetName()->Enquote('"')
                << "\n";
    }

    for (unsigned long i = 0; i < compList1.countitems(); i++) {
      long i1 = compList1.get(i), i2 = compList2.get(i);

      (*result) << "\nStep " << _String(i + 1) << ": ";
      if (N >= 0 && i1 == N) {
        (*result) << " last var "
                  << LocateVar(variableIndex.get(N))->GetName()->Enquote('"');
        if (i2 < 0) {
          (*result) << " raise to power " << _String(-i2 - 1);
        } else {
          (*result) << " " << _String(i2) << " consecutive powers";
        }
      } else {
        if (i1 < 0) {
          (*result)
              << "Operation on "
              << LocateVar(variableIndex.get(-i1 - 1))->GetName()->Enquote('"')
              << " power " << _String(fabs(i2));
          if (i2 < 0) {
            (*result) << " add to results";
          }
        } else {
          (*result)
              << "Reset from "
              << LocateVar(variableIndex.get(i1 + 1))->GetName()->Enquote('"')
              << " power " << _String(fabs(i2));
          if (i2 < 0) {
            (*result) << " add to results";
          }
        }
      }
    }
    (*result) << "\n";
  }
  Convert2OperationForm();
  if (theTerms->NumberOfTerms()) {
    unsigned long i;
    _List _varNames;
    for (i = 0; i < variableIndex.countitems(); i++) {
      _varNames << LocateVar(variableIndex(i))->GetName();
    }

    bool firstN = theTerms->IsFirstANumber();
    for (i = 0; i < theTerms->NumberOfTerms(); i++) {
      char number[100];

      hyFloat coeff = theTerms->GetCoeff(i);

      bool print_coeff = false;

      if ((i == 0 && firstN) || !CheckEqual(coeff, 1)) {
        print_coeff = true;
      }

      snprintf(number, sizeof(number), PRINTF_FORMAT_STRING, coeff);
      if (i > 0 && coeff >= 0.) {
        *result << '+';
      }

      if (print_coeff) {
        *result << number;
      }

      if (i > 0 || !firstN) {
        if (print_coeff) {
          *result << '*';
        }
        long *cT = theTerms->GetTerm(i);
        bool printedFirst = false;
        for (unsigned long k = 0; k < variableIndex.countitems(); k++, cT++) {
          if (*cT != 0) {
            if (printedFirst) {
              *result << '*';
            } else {
              printedFirst = true;
            }

            *result << (_String *)_varNames(k);
            if (*cT > 1) {
              *result << '^' << _String(*cT);
            } else {
              if (*cT < 0) {
                *result << '^' << _String(*cT).Enquote('(', ')');
              }
            }
          }
        }
      }
    }
  } else {
    *result << _String((_String *)compList1.toStr(padding)) << '\n'
            << _String((_String *)compList2.toStr(padding)) << '\n';
  }
  result->TrimSpace();
  return result;
}

//__________________________________________________________________________________

void _Polynomial::toFileStr(hyFile *f, unsigned long padding) {
  Convert2OperationForm();
  if (theTerms->NumberOfTerms() && theTerms->thePowers) {
    f->puts("p(");
    _List _varNames;
    unsigned long i;
    for (i = 0; i < variableIndex.countitems(); i++) {
      _varNames << LocateVar(variableIndex(i))->GetName();
      f->puts(((_String *)_varNames(i))->get_str());
      if (i + 1 < variableIndex.countitems()) {
        f->puts(",");
      }
    }
    f->puts(")=");
    for (i = 0; i < theTerms->NumberOfTerms(); i++) {
      char number[100];
      snprintf(number, sizeof(number), PRINTF_FORMAT_STRING,
               theTerms->GetCoeff(i));
      if ((i > 0) && (number[0] != '-')) {
        f->puts("+");
      }
      f->puts(number);
      if ((i > 0) || !theTerms->IsFirstANumber()) {
        f->puts("*");
        long *cT = theTerms->GetTerm(i);
        for (unsigned long k = 0; k < variableIndex.countitems(); k++, cT++) {
          if (*cT > 0) {
            f->puts(((_String *)_varNames(k))->get_str());
            if (*cT > 1) {
              f->puts("^");
              snprintf(number, sizeof(number), "%ld", *cT);
              f->puts(number);
            }
          }
        }
      }
    }
  } else {
    compList1.toFileStr(f, padding);
    compList2.toFileStr(f, padding);
  }
}

//__________________________________________________________________________________
void _Polynomial::ScanForVariables(_AVLList &l, bool globals, _AVLListX *tagger,
                                   long weight) const {
  for (unsigned long i = 0; i < variableIndex.lLength; i++) {
    long vi = variableIndex(i);

    _Variable *v = LocateVar(vi);
    if (v->IsGlobal()) {
      if (globals) {
        l.Insert((BaseRef)vi);
        if (tagger) {
          tagger->UpdateValue((BaseRef)vi, weight, 0);
        }
      }
    } else {
      l.Insert((BaseRef)vi);
      if (tagger) {
        tagger->UpdateValue((BaseRef)vi, weight, 0);
      }
    }
  }
}

//__________________________________________________________________________________
bool _Polynomial::IsObjectEmpty(void) {
  if (compList1.countitems()) {
    return false;
  }
  if (theTerms->NumberOfTerms()) {
    if (theTerms->NumberOfTerms() == 1) {
      if (theTerms->IsFirstANumber()) {
        return theTerms->theCoeff[0] == 0.0;
      }
    }
    return false;
  }
  return true;
}
//__________________________________________________________________________________
bool _Polynomial::HasChanged(bool gl, _AVLListX *cache) {
  for (long k = variableIndex.countitems() - 1; k >= 0; k--) {
    if (LocateVar(variableIndex(k))->HasChanged(gl, cache)) {
      return true;
    }
  }
  return false;
}

//__________________________________________________________________________________
bool _Polynomial::IsMaxElement(hyFloat bench) {
  hyFloat *tc = theTerms->GetCoeff();
  for (long k = 0; k < theTerms->actTerms; k++, tc++) {
    if (fabs(*tc) >= bench) {
      return true;
    }
  }
  return false;
}

//__________________________________________________________________________________
void SetPolyTermCap(long t) { polyTermCap = t; }
