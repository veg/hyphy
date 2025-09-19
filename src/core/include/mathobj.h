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

#ifndef __MATHOBJ__
#define __MATHOBJ__

#include "_hyExecutionContext.h"
#include "avllistx.h"
#include "baseobj.h"
#include "defines.h"
#include "hy_strings.h"

class _MathObject : public BaseObj { // abstract math operations class

protected:
  _MathObject *_extract_argument(_List *arguments, unsigned long index,
                                 bool fill_in) const;

public:
  /**
   * @brief Add two math objects
   *
   * @param o The object to add
   * @param cache The cache to use
   * @return _MathObject* The result of the addition
   */
  virtual _MathObject *Add(_MathObject *, _MathObject *cache = nil);
  /**
   * @brief Subtract two math objects
   *
   * @param o The object to subtract
   * @param cache The cache to use
   * @return _MathObject* The result of the subtraction
   */
  virtual _MathObject *Sub(_MathObject *, _MathObject *cache = nil);
  /**
   * @brief Negate a math object
   *
   * @param cache The cache to use
   * @return _MathObject* The negated object
   */
  virtual _MathObject *Minus(_MathObject *cache = nil);
  /**
   * @brief Sum a math object
   *
   * @param cache The cache to use
   * @return _MathObject* The sum
   */
  virtual _MathObject *Sum(_MathObject *cache = nil);
  /**
   * @brief Multiply two math objects
   *
   * @param o The object to multiply by
   * @param cache The cache to use
   * @return _MathObject* The result of the multiplication
   */
  virtual _MathObject *Mult(_MathObject *, _MathObject *cache = nil);
  /**
   * @brief Divide two math objects
   *
   * @param o The object to divide by
   * @param cache The cache to use
   * @return _MathObject* The result of the division
   */
  virtual _MathObject *Div(_MathObject *, _MathObject *cache = nil);
  /**
   * @brief Left divide two math objects
   *
   * @param o The object to divide by
   * @param cache The cache to use
   * @return _MathObject* The result of the division
   */
  virtual _MathObject *lDiv(_MathObject *, _MathObject *cache = nil);
  /**
   * @brief Long divide two math objects
   *
   * @param o The object to divide by
   * @param cache The cache to use
   * @return _MathObject* The result of the division
   */
  virtual _MathObject *longDiv(_MathObject *, _MathObject *cache = nil);
  /**
   * @brief Raise a math object to a power
   *
   * @param o The power to raise to
   * @param cache The cache to use
   * @return _MathObject* The result of the exponentiation
   */
  virtual _MathObject *Raise(_MathObject *, _MathObject *cache = nil);
  /**
   * @brief Check if two math objects are equal
   *
   * @param o The object to compare with
   * @return true if the objects are equal, false otherwise
   */
  virtual bool Equal(_MathObject *);
  /**
   * @brief Get the absolute value of a math object
   *
   * @param cache The cache to use
   * @return _MathObject* The absolute value of the object
   */
  virtual _MathObject *Abs(_MathObject *cache = nil);
  /**
   * @brief Get the sine of a math object
   *
   * @param cache The cache to use
   * @return _MathObject* The sine of the object
   */
  virtual _MathObject *Sin(_MathObject *cache = nil);
  /**
   * @brief Get the cosine of a math object
   *
   * @param cache The cache to use
   * @return _MathObject* The cosine of the object
   */
  virtual _MathObject *Cos(_MathObject *cache = nil);
  /**
   * @brief Get the tangent of a math object
   *
   * @param cache The cache to use
   * @return _MathObject* The tangent of the object
   */
  virtual _MathObject *Tan(_MathObject *cache = nil);
  /**
   * @brief Get the exponential of a math object
   *
   * @param cache The cache to use
   * @return _MathObject* The exponential of the object
   */
  virtual _MathObject *Exp(_MathObject *cache = nil);
  /**
   * @brief Get the natural logarithm of a math object
   *
   * @param cache The cache to use
   * @return _MathObject* The natural logarithm of the object
   */
  virtual _MathObject *Log(_MathObject *cache = nil);
  /**
   * @brief Get the square root of a math object
   *
   * @param cache The cache to use
   * @return _MathObject* The square root of the object
   */
  virtual _MathObject *Sqrt(_MathObject *cache = nil);
  /**
   * @brief Get the gamma function of a math object
   *
   * @param cache The cache to use
   * @return _MathObject* The gamma function of the object
   */
  virtual _MathObject *Gamma(_MathObject *cache = nil);
  /**
   * @brief Get the error function of a math object
   *
   * @param cache The cache to use
   * @return _MathObject* The error function of the object
   */
  virtual _MathObject *Erf(_MathObject *cache = nil);
  /**
   * @brief Get the natural logarithm of the gamma function of a math object
   *
   * @param cache The cache to use
   * @return _MathObject* The natural logarithm of the gamma function of the object
   */
  virtual _MathObject *LnGamma(_MathObject *cache = nil);
  /**
   * @brief Get the beta function of two math objects
   *
   * @param o The second object
   * @param cache The cache to use
   * @return _MathObject* The beta function of the two objects
   */
  virtual _MathObject *Beta(_MathObject *, _MathObject *cache = nil);
  /**
   * @brief Get the incomplete gamma function of two math objects
   *
   * @param o The second object
   * @param cache The cache to use
   * @return _MathObject* The incomplete gamma function of the two objects
   */
  virtual _MathObject *IGamma(_MathObject *, _MathObject *cache = nil);
  /**
   * @brief Get the cumulative chi-squared distribution of two math objects
   *
   * @param o The second object
   * @param cache The cache to use
   * @return _MathObject* The cumulative chi-squared distribution of the two objects
   */
  virtual _MathObject *CChi2(_MathObject *, _MathObject *cache = nil);
  /**
   * @brief Get the incomplete beta function of three math objects
   *
   * @param o1 The second object
   * @param o2 The third object
   * @param cache The cache to use
   * @return _MathObject* The incomplete beta function of the three objects
   */
  virtual _MathObject *IBeta(_MathObject *, _MathObject *,
                             _MathObject *cache = nil);
  /**
   * @brief Perform a simplex optimization on a math object
   *
   * @param cache The cache to use
   * @return _MathObject* The result of the simplex optimization
   */
  virtual _MathObject *Simplex(_MathObject *cache = nil);
  /**
   * @brief Convert a math object to a polynomial
   *
   * @param cache The cache to use
   * @return _MathObject* The polynomial representation of the object
   */
  virtual _MathObject *ConvertToPolynomial(_MathObject *cache = nil);

  /**
   * @brief Simplify a math object
   *
   * @param cache The cache to use
   * @return _MathObject* The simplified object
   */
  virtual _MathObject *Simplify(_MathObject *cache = nil);

  /**
   * @brief Get the minimum of two math objects
   *
   * @param o The second object
   * @param cache The cache to use
   * @return _MathObject* The minimum of the two objects
   */
  virtual _MathObject *Min(_MathObject *, _MathObject *cache = nil);
  /**
   * @brief Get the maximum of two math objects
   *
   * @param o The second object
   * @param cache The cache to use
   * @return _MathObject* The maximum of the two objects
   */
  virtual _MathObject *Max(_MathObject *, _MathObject *cache = nil);
  /**
   * @brief Get the inverse chi-squared distribution of two math objects
   *
   * @param o The second object
   * @param cache The cache to use
   * @return _MathObject* The inverse chi-squared distribution of the two objects
   */
  virtual _MathObject *InvChi2(_MathObject *, _MathObject *cache = nil);
  /**
   * @brief Get the standard normal cumulative distribution function of a math object
   *
   * @param cache The cache to use
   * @return _MathObject* The standard normal cumulative distribution function of the object
   */
  virtual _MathObject *ZCDF(_MathObject *cache = nil);
  /**
   * @brief Get the current time
   *
   * @param cache The cache to use
   * @return _MathObject* The current time
   */
  virtual _MathObject *Time(_MathObject *cache = nil);
  /**
   * @brief Get the arctangent of a math object
   *
   * @param cache The cache to use
   * @return _MathObject* The arctangent of the object
   */
  virtual _MathObject *Arctan(_MathObject *cache = nil);
  /**
   * @brief Check if a math object is less than another
   *
   * @param o The object to compare with
   * @param cache The cache to use
   * @return _MathObject* 1 if less, 0 otherwise
   */
  virtual _MathObject *Less(_MathObject *, _MathObject *cache = nil);
  /**
   * @brief Generate a random number
   *
   * @param o The type of random number to generate
   * @param cache The cache to use
   * @return _MathObject* The random number
   */
  virtual _MathObject *Random(_MathObject *, _MathObject *cache = nil);
  /**
   * @brief Check if a math object is greater than another
   *
   * @param o The object to compare with
   * @param cache The cache to use
   * @return _MathObject* 1 if greater, 0 otherwise
   */
  virtual _MathObject *Greater(_MathObject *, _MathObject *cache = nil);
  /**
   * @brief Check if a math object is less than or equal to another
   *
   * @param o The object to compare with
   * @param cache The cache to use
   * @return _MathObject* 1 if less or equal, 0 otherwise
   */
  virtual _MathObject *LessEq(_MathObject *, _MathObject *cache = nil);
  /**
   * @brief Check if a math object is greater than or equal to another
   *
   * @param o The object to compare with
   * @param cache The cache to use
   * @return _MathObject* 1 if greater or equal, 0 otherwise
   */
  virtual _MathObject *GreaterEq(_MathObject *, _MathObject *cache = nil);
  /**
   * @brief Check if two math objects are equal
   *
   * @param o The object to compare with
   * @param cache The cache to use
   * @return _MathObject* 1 if equal, 0 otherwise
   */
  virtual _MathObject *AreEqual(_MathObject *, _MathObject *cache = nil);
  /**
   * @brief Check if two math objects are not equal
   *
   * @param o The object to compare with
   * @param cache The cache to use
   * @return _MathObject* 1 if not equal, 0 otherwise
   */
  virtual _MathObject *NotEqual(_MathObject *, _MathObject *cache = nil);
  /**
   * @brief Logical AND of two math objects
   *
   * @param o The second object
   * @param cache The cache to use
   * @return _MathObject* The result of the logical AND
   */
  virtual _MathObject *LAnd(_MathObject *, _MathObject *cache = nil);
  /**
   * @brief Logical OR of two math objects
   *
   * @param o The second object
   * @param cache The cache to use
   * @return _MathObject* The result of the logical OR
   */
  virtual _MathObject *LOr(_MathObject *, _MathObject *cache = nil);
  /**
   * @brief Get the gamma distribution of three math objects
   *
   * @param o1 The second object
   * @param o2 The third object
   * @param cache The cache to use
   * @return _MathObject* The gamma distribution of the three objects
   */
  virtual _MathObject *GammaDist(_MathObject *, _MathObject *,
                                 _MathObject *cache = nil);
  /**
   * @brief Get the cumulative gamma distribution of three math objects
   *
   * @param o1 The second object
   * @param o2 The third object
   * @param cache The cache to use
   * @return _MathObject* The cumulative gamma distribution of the three objects
   */
  virtual _MathObject *CGammaDist(_MathObject *, _MathObject *,
                                  _MathObject *cache = nil);
  /**
   * @brief Logical NOT of a math object
   *
   * @param cache The cache to use
   * @return _MathObject* The result of the logical NOT
   */
  virtual _MathObject *LNot(_MathObject *cache = nil);
  /**
   * @brief Get the number of tips in a tree
   *
   * @param cache The cache to use
   * @return _MathObject* The number of tips
   */
  virtual _MathObject *TipCount(_MathObject *cache = nil);
  /**
   * @brief Get the number of branches in a tree
   *
   * @param cache The cache to use
   * @return _MathObject* The number of branches
   */
  virtual _MathObject *BranchCount(_MathObject *cache = nil);
  /**
   * @brief Get the name of a tip in a tree
   *
   * @param o The index of the tip
   * @param cache The cache to use
   * @return _MathObject* The name of the tip
   */
  virtual _MathObject *TipName(_MathObject *, _MathObject *cache = nil);
  /**
   * @brief Get the name of a branch in a tree
   *
   * @param o The index of the branch
   * @param cache The cache to use
   * @return _MathObject* The name of the branch
   */
  virtual _MathObject *BranchName(_MathObject *, _MathObject *cache = nil);
  /**
   * @brief Get the length of a branch in a tree
   *
   * @param o The index of the branch
   * @param cache The cache to use
   * @return _MathObject* The length of the branch
   */
  virtual _MathObject *BranchLength(_MathObject *, _MathObject *cache = nil);
  /**
   * @brief Reroot a tree
   *
   * @param o The new root
   * @param cache The cache to use
   * @return _MathObject* The rerooted tree
   */
  virtual _MathObject *RerootTree(_MathObject *, _MathObject *cache = nil);
  /**
   * @brief Get the Newick string representation of a tree
   *
   * @param o The tree
   * @param cache The cache to use
   * @return _MathObject* The Newick string
   */
  virtual _MathObject *TEXTreeString(_MathObject *,
                                     _MathObject *cache = nil) const;
  /**
   * @brief Get the type of a math object
   *
   * @param cache The cache to use
   * @return _MathObject* The type of the object
   */
  virtual _MathObject *Type(_MathObject *cache = nil);
  /**
   * @brief Get the plain string representation of a tree
   *
   * @param o1 The first object
   * @param o2 The second object
   * @param cache The cache to use
   * @return _MathObject* The plain string representation of the tree
   */
  virtual _MathObject *PlainTreeString(_MathObject *, _MathObject *,
                                       _MathObject *cache = nil);
  /**
   * @brief Format a number as a string
   *
   * @param o1 The number to format
   * @param o2 The format string
   * @param cache The cache to use
   * @return _MathObject* The formatted string
   */
  virtual _MathObject *FormatNumberString(_MathObject *, _MathObject *,
                                          _MathObject *cache = nil);
  /**
   * @brief Get the value of a math object
   *
   * @return hyFloat The value of the object
   */
  virtual hyFloat Value(void);

  /**
   * @brief Compute the value of a math object
   *
   * @return _MathObject* The computed object
   */
  virtual _MathObject *Compute(void) { return this; }
  /**
   * @brief Scan for variables in a math object
   *
   * @param dest The destination list
   * @param b Whether to scan for variables
   * @param avl The AVL list
   * @param l The long value
   */
  virtual void ScanForVariables(_AVLList &, bool = false, _AVLListX * = nil,
                                long = 0) const {}

  /**
   * @brief Make a dynamic copy of a math object
   *
   * @return BaseRef The dynamic copy of the object
   */
  virtual BaseRef makeDynamic(void) const;
  /**
   * @brief Duplicate a math object
   *
   * @param brc The object to duplicate
   */
  virtual void Duplicate(BaseRefConst);

  /**
   * @brief Check if a math object is a variable
   *
   * @return true if the object is a variable, false otherwise
   */
  virtual bool IsVariable(void) { return false; }
  /**
   * @brief Check if a math object is empty
   *
   * @return true if the object is empty, false otherwise
   */
  virtual bool IsObjectEmpty(void) { return true; }
  /**
   * @brief Check if a math object is printable
   *
   * @return true if the object is printable, false otherwise
   */
  virtual bool IsPrintable(void) { return false; }

  /**
   * @brief Check if a math object is independent
   *
   * @return true if the object is independent, false otherwise
   */
  virtual bool IsIndependent(void) { return true; }
  /**
   * @brief Get the object class of a math object
   *
   * @return unsigned long The object class
   */
  virtual unsigned long ObjectClass(void) const { return HY_UNDEFINED; }
  // returns a unique ID for this object
  // 0 - undefined
  // 1 - number
  // 4 - matrix

  /**
   * @brief Execute a single operation
   *
   * @param opCode The opcode of the operation
   * @param arguments The arguments of the operation
   * @param context The execution context
   * @param cache The cache to use
   * @return _MathObject* The result of the operation
   */
  virtual _MathObject *
  ExecuteSingleOp(long opCode, _List *arguments = nil,
                  _hyExecutionContext *context = _hyDefaultExecutionContext,
                  _MathObject *cache = nil);
  // execute this operation with the list of Args

  /**
   * @brief Check if a math object has changed
   *
   * @param b Whether the object has changed
   * @param avl The AVL list
   * @return true if the object has changed, false otherwise
   */
  virtual bool HasChanged(bool = false, _AVLListX * = nil) { return false; }

  /**
   * @brief Check if a math object is a constant
   *
   * @return true if the object is a constant, false otherwise
   */
  virtual bool IsConstant(void) { return true; }

private:
  _MathObject *_null_handler();
};

// pointer to a math object

/** @brief A pointer to a math object */
typedef _MathObject *HBLObjectRef;
/** @brief A pointer to a constant math object */
typedef _MathObject const *HBLObjectRefConst;

/**
 * @brief Return a constant or use the cache
 *
 * @param value The value to return
 * @param cache The cache to use
 * @return HBLObjectRef The constant
 */
HBLObjectRef _returnConstantOrUseCache(hyFloat value, HBLObjectRef cache);
/**
 * @brief Return null or use the cache
 *
 * @param cache The cache to use
 * @return HBLObjectRef Null or the cached object
 */
HBLObjectRef _returnNullOrUseCache(HBLObjectRef cache);

#endif
