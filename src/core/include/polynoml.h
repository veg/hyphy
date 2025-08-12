/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-now
Core Developers:
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon    (apoon42@uwo.ca)
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

#ifndef __POLY__
#define __POLY__

#include "parser.h"

#define POLY_DATA_INCREMENT 8
#define GLOBAL_VARIABLE 1
#define CATEGORY_VARIABLE 2
#define RANDOM_VARIABLE 3

//__________________________________________________________________________________

/**
 * @class _PolynomialData
 * @brief A class to store a single polynomial term
 * @description This object stores terms term of the form C * x_1^p_1 * x_2^p_2
 * * ... * x_n^p_n C is the numeric coeffieint (stored in theCoeff) The number
 * of variables is stored in numberVars Variables are tracked by their index in
 * the array of powers (thePowers) Some powers can be 0 (variable not present)
 *      Offsets within this array (thePowers) are multiples of numberVars
 *      actStored and allocTerms allow storage increments > 1 when arrays need
 * to be resized The terms are stored lexicographically in terms of powers, so
 * that (p_1, ..., p_n) vectors are in increasing order
 */
class _PolynomialData : public BaseObj {

public:
  /**
   * @brief Construct a new _PolynomialData object.
   */
  _PolynomialData(void);
  /**
   * @brief Construct a new _PolynomialData object.
   * @param vars The number of variables.
   */
  _PolynomialData(const long);
  /**
   * @brief Construct a new _PolynomialData object.
   * @param p The _PolynomialData object to copy.
   */
  _PolynomialData(_PolynomialData const &);
  /**
   * @brief Construct a new _PolynomialData object.
   * @param vars The number of variables.
   * @param terms The number of terms.
   * @param coeffs The coefficients.
   */
  _PolynomialData(const long, const long, hyFloat const *);
  /**
   * @brief Assign a _PolynomialData object to this object.
   * @param p The _PolynomialData object to assign.
   * @return This object.
   */
  _PolynomialData const &operator=(_PolynomialData const &);

  /**
   * @brief Destroy the _PolynomialData object.
   */
  virtual ~_PolynomialData();

  /**
   * @brief Make a dynamic copy of this object.
   * @return A pointer to the new object.
   */
  virtual BaseObj *makeDynamic(void) const;
  /**
   * @brief Duplicate this object.
   * @param p The object to duplicate.
   */
  virtual void Duplicate(BaseRefConst);

  /**
   * @brief Get the coefficients.
   * @return A pointer to the coefficients.
   */
  inline hyFloat *GetCoeff(void) { return theCoeff; }
  /**
   * @brief Get a coefficient.
   * @param index The index of the coefficient.
   * @return The coefficient.
   */
  inline hyFloat &GetCoeff(long index) { return theCoeff[index]; }

  /**
   * @brief Get a term.
   * @param index The index of the term.
   * @return A pointer to the term.
   */
  long *GetTerm(long) const;
  /**
   * @brief Get the number of terms.
   * @return The number of terms.
   */
  long GetNoTerms(void) { return actTerms; }
  /**
   * @brief Add a term to the polynomial.
   * @param term The term to add.
   * @param coeff The coefficient of the term.
   */
  void AddTerm(long *, hyFloat);
  /**
   * @brief Add a term to the polynomial.
   * @param term The term to add.
   * @param coeff The coefficient of the term.
   * @param reindexer A pointer to the reindexer.
   * @param actLength The active length.
   */
  void AddTerm(long *, hyFloat, long *, long);
  /**
   * @brief Add a term to the polynomial.
   * @param coeff The coefficient of the term.
   */
  void AddTerm(hyFloat);
  /**
   * @brief Write a term to the polynomial.
   * @param term The term to write.
   * @param index The index of the term.
   */
  void WriteTerm(long *, long);
  /**
   * @brief Delete a term from the polynomial.
   * @param index The index of the term to delete.
   */
  void DeleteTerm(long);
  /**
   * @brief Check if the first term is a number.
   * @return True if the first term is a number, false otherwise.
   */
  bool IsFirstANumber(void);
  /**
   * @brief Get the number of terms.
   * @return The number of terms.
   */
  inline unsigned long NumberOfTerms(void) { return actTerms; }
  /**
   * @brief Get the sum of the powers of a term.
   * @param index The index of the term.
   * @return The sum of the powers.
   */
  long SumOfPowers(long) const;
  /**
   * @brief Get the weighted sum of the powers of a term.
   * @param index The index of the term.
   * @param w The weights.
   * @return The weighted sum of the powers.
   */
  hyFloat WeightedSumOfPowers(long, hyFloat *) const;

  /**
   * @brief Check the polynomial data.
   * @return True if the polynomial data is valid, false otherwise.
   */
  bool checkMe(void);

  friend class _Polynomial;

  /**
   * @brief Multiply two terms.
   * @param target The target term.
   * @param s1 The first source term.
   * @param s2 The second source term.
   */
  void MultiplyTerms(long *, long *, long *);
  /**
   * @brief Raise a term to a power.
   * @param target The term to raise.
   * @param power The power.
   */
  void RaiseTerm(long *, long);
  /**
   * @brief Raise a number to a power.
   * @param base The base.
   * @param power The power.
   * @return The result.
   */
  static hyFloat BinaryRaise(hyFloat, long);
  /**
   * @brief Rearrange a term.
   * @param target The target term.
   * @param source The source term.
   * @param markup The markup.
   * @param items The number of items.
   */
  static void RearrangeTerm(long *, long *, long *, long);
  /**
   * @brief Compare two terms.
   * @param s1 The first term.
   * @param s2 The second term.
   * @return 0 if the terms are equal, 1 if the first term is greater, -1 if the
   * second term is greater.
   */
  char CompareTerms(long *, long *) const;
  /**
   * @brief Compare two terms.
   * @param s1 The first term.
   * @param s2 The second term.
   * @param secondReindex The reindex for the second term.
   * @param actLength The active length.
   * @return 0 if the terms are equal, 1 if the first term is greater, -1 if the
   * second term is greater.
   */
  char CompareTerms(long *, long *, long *, long) const;
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
  char CompareTerms(long *, long *, long *, long *, long, long) const;
  /**
   * @brief Find a term.
   * @param theTerm The term to find.
   * @param reIndex The reindex.
   * @param start The start index.
   * @return The index of the term, or a negative number if not found.
   */
  long FindTerm(long *, long *, long start = 0) const;
  /**
   * @brief Resort the terms.
   * @param reIndex The reindex.
   */
  void ResortTerms(long *);
  /**
   * @brief Chop the terms.
   */
  void ChopTerms(void);
  /**
   * @brief Check a term.
   * @param myCoeff The coefficient of the term.
   * @param myIndex The index of the term.
   * @return True if the term is valid, false otherwise.
   */
  bool KeepTerm(hyFloat, long) const;

protected:
  hyFloat *theCoeff;
  long *thePowers;
  long numberVars, actTerms, allocTerms;

private:
  /**
   * @brief Check if there's room to add a new term, and if not increase storage
   */
  void EnsureSpaceIsAvailable(void);
  /**
   * @brief Check if there's too much allocation buffer and shrink if necessary
   */
  void ShrinkExtraSpace(void);
};

//__________________________________________________________________________________

/**
 * @class _Polynomial
 * @brief A class to represent a polynomial.
 */
class _Polynomial : public _MathObject {

public:
  /**
   * @brief Construct a new _Polynomial object.
   */
  _Polynomial(void);
  /**
   * @brief Construct a new _Polynomial object.
   * @param vars The variables of the polynomial.
   */
  _Polynomial(_SimpleList &);
  /**
   * @brief Construct a new _Polynomial object.
   * @param p The _Polynomial object to copy.
   */
  _Polynomial(_Polynomial const &);
  /**
   * @brief Assign a _Polynomial object to this object.
   * @param p The _Polynomial object to assign.
   * @return This object.
   */
  _Polynomial const &operator=(_Polynomial const &);
  /**
   * @brief Construct a new _Polynomial object.
   * @param value The value of the polynomial.
   */
  _Polynomial(hyFloat);
  /**
   * @brief Construct a new _Polynomial object.
   * @param v The variable of the polynomial.
   */
  _Polynomial(_Variable &);
  /**
   * @brief Destroy the _Polynomial object.
   */
  virtual ~_Polynomial();
  /**
   * @brief Execute a single operation.
   * @param opCode The operation code.
   * @param arguments The arguments.
   * @param context The execution context.
   * @param cache The cache.
   * @return The result of the operation.
   */
  virtual _MathObject *ExecuteSingleOp(
      long opCode, _List *arguments = nil,
      _hyExecutionContext *context = _hyDefaultExecutionContext,
      HBLObjectRef cache = nil); // execute this operation with the list of Args

  /**
   * @brief Make a dynamic copy of this object.
   * @return A pointer to the new object.
   */
  virtual BaseObj *makeDynamic(void) const;
  /**
   * @brief Duplicate this object.
   * @param p The object to duplicate.
   */
  virtual void Duplicate(BaseRefConst);

  /**
   * @brief Add a math object to this polynomial.
   * @param m The math object to add.
   * @param cache The cache.
   * @return The result of the addition.
   */
  virtual _MathObject *Add(_MathObject *, HBLObjectRef cache);
  /**
   * @brief Add a math object to this polynomial.
   * @param m The math object to add.
   * @param subtract True to subtract, false to add.
   * @param cache The cache.
   * @return The result of the addition.
   */
  virtual _MathObject *Plus(_MathObject *, bool subtract, HBLObjectRef cache);
  /**
   * @brief Subtract a math object from this polynomial.
   * @param m The math object to subtract.
   * @param cache The cache.
   * @return The result of the subtraction.
   */
  virtual _MathObject *Sub(_MathObject *, HBLObjectRef cache);
  /**
   * @brief Raise this polynomial to a power.
   * @param m The power.
   * @param cache The cache.
   * @return The result of the exponentiation.
   */
  virtual _MathObject *Raise(_MathObject *, HBLObjectRef cache);
  /**
   * @brief Negate this polynomial.
   * @param cache The cache.
   * @return The negated polynomial.
   */
  virtual _MathObject *Minus(HBLObjectRef cache);
  /**
   * @brief Multiply this polynomial by a math object.
   * @param m The math object to multiply by.
   * @param cache The cache.
   * @return The result of the multiplication.
   */
  virtual _MathObject *Mult(_MathObject *, HBLObjectRef cache);
  /**
   * @brief Compute the value of the polynomial.
   * @return The value of the polynomial.
   */
  virtual _MathObject *Compute(void);
  /**
   * @brief Check if this polynomial is equal to another math object.
   * @param m The math object to compare to.
   * @return True if the objects are equal, false otherwise.
   */
  virtual bool Equal(_MathObject *);
  /**
   * @brief Compute the value of the polynomial.
   * @return The value of the polynomial.
   */
  hyFloat ComputePolynomial(void);

  /**
   * @brief Compute the value of the polynomial.
   * @param c1 The first list of coefficients.
   * @param c2 The second list of coefficients.
   * @param nc1 The number of coefficients in the first list.
   * @param nc2 The number of coefficients in the second list.
   * @param p1 The first list of powers.
   * @param p2 The second list of powers.
   * @return The value of the polynomial.
   */
  hyFloat ComputeP(hyFloat *, hyFloat *, long, long, long *, long *);
  /**
   * @brief Check if this polynomial is a number.
   * @param returnLeading True to return the leading coefficient, false
   * otherwise.
   * @return A pointer to a _Constant object if the polynomial is a number, nil
   * otherwise.
   */
  _MathObject *IsANumber(bool = false);
  /**
   * @brief Check if this polynomial is empty.
   * @return True if the polynomial is empty, false otherwise.
   */
  virtual bool IsObjectEmpty(void);

  /**
   * @brief Get the object class.
   * @return The object class.
   */
  virtual unsigned long ObjectClass(void) const { return POLYNOMIAL; }
  /**
   * @brief Get the value of the polynomial.
   * @return The value of the polynomial.
   */
  virtual hyFloat Value(void) { return ComputePolynomial(); }

  /**
   * @brief Convert the polynomial to a string.
   * @param format The format of the string.
   * @return A pointer to the string.
   */
  virtual BaseObj *toStr(unsigned long = 0UL);
  /**
   * @brief Check a term of the polynomial.
   */
  void CheckTerm(void);

  /**
   * @brief Convert the polynomial to a file string.
   * @param file The file to write to.
   * @param format The format of the string.
   */
  virtual void toFileStr(hyFile *, unsigned long = 0UL);

  /**
   * @brief Get the number of variables.
   * @return The number of variables.
   */
  long GetNoVariables(void) const { return variableIndex.countitems(); }

  /**
   * @brief Get the Ith variable.
   * @param i The index of the variable.
   * @return A pointer to the variable.
   */
  _Variable *GetIthVariable(unsigned long i) const {
    return LocateVar(variableIndex.get(i));
  }

  /**
   * @brief Get the terms of the polynomial.
   * @return A pointer to the terms.
   */
  _PolynomialData *GetTheTerms(void) const { return theTerms; }
  /**
   * @brief Set the terms of the polynomial.
   * @param td A pointer to the terms.
   */
  void SetTheTerms(_PolynomialData *td) { theTerms = td; }
  /**
   * @brief Set the computation lists.
   * @param c1 The first computation list.
   * @param c2 The second computation list.
   */
  void SetCLists(_SimpleList &c1, _SimpleList &c2) {
    compList1.Duplicate(&c1);
    compList2.Duplicate(&c2);
  }
  /**
   * @brief Scan for variables.
   * @param l The list of variables.
   * @param globals True to scan for global variables, false otherwise.
   * @param tagger The tagger.
   * @param weight The weight.
   */
  virtual void ScanForVariables(_AVLList &l, bool globals = false,
                                _AVLListX *tagger = nil, long weight = 0) const;
  /**
   * @brief Check if the polynomial has changed.
   * @param ignore A boolean value.
   * @param p A pointer to a list of variables.
   * @return True if the polynomial has changed, false otherwise.
   */
  virtual bool HasChanged(bool = false, _AVLListX * = nil);
  /**
   * @brief Reset the polynomial check.
   * @param p The polynomial to reset.
   */
  friend void ResetPolynomialCheck(_Polynomial *);
  /**
   * @brief Get the computational size.
   * @return The computational size.
   */
  long ComputationalSize(void) { return compList1.countitems(); }
  /**
   * @brief Check if the polynomial is the maximum element.
   * @param v The value to check.
   * @return True if the polynomial is the maximum element, false otherwise.
   */
  bool IsMaxElement(hyFloat);
  /**
   * @brief Convert the polynomial to computation form.
   * @param c1 The first computation list.
   * @param c2 The second computation list.
   * @param termsToInclude The terms to include.
   */
  void Convert2ComputationForm(_SimpleList *c1 = nil, _SimpleList *c2 = nil,
                               _SimpleList *termsToInclude = nil);
  /**
   * @brief Rank the terms of the polynomial.
   * @param p The list of terms to rank.
   */
  void RankTerms(_SimpleList *);

protected:
  /**
   * @brief Check if the numeric term is 0, and drop it if so
   */
  void DropZeroCoefficient(void);

  /**
   * @brief Set all coefficients to -coefficient; return self
   */
  _Polynomial *NegateCoefficients(void);

  /**
   * @brief Drop small terms from the polynomial.
   */
  void DropSmallTerms(void);
  /**
   * @brief Convert the polynomial to operation form.
   */
  void Convert2OperationForm(void);

  _SimpleList variableIndex, compList1, compList2;

  _PolynomialData *theTerms;
};

extern hyFloat dropPrecision, topPolyCap, polynomialExpPrecision;

extern bool dropTerms, enforcePolyCap;

extern const unsigned long maximumPolyTermsPerVariable,
    maxPolynomialExpIterates;
/**
 * @brief Set the polynomial term cap.
 * @param cap The new cap.
 */
void SetPolyTermCap(long);

#endif
