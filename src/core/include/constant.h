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

#ifndef __CONSTANT__
#define __CONSTANT__

#include "global_things.h"
#include "mathobj.h"

#define _HY_CONSTANT_PREALLOCATE_SLOTS 16384

/**
 * @brief A constant value class
 *
 */
class _Constant : public _MathObject {

private:
  /**
   * @brief A template function to check the type of an operand and compute the result of a binary operation.
   *
   * @tparam T The type of the functor to apply.
   * @param operand The operand to compute.
   * @param functor The functor to apply.
   * @param cache The cache to use for the result.
   * @return The result of the computation.
   */
  template <class T>
  HBLObjectRef _check_type_and_compute(HBLObjectRef operand, T functor,
                                       HBLObjectRef cache) {
    if (operand) {
      if (operand->ObjectClass() == NUMBER) {
        return _returnConstantOrUseCache(
            functor(Value(), ((_Constant *)operand)->Value()), cache);
      } else {
        _String error_message =
            _String("<'constant' operation 'X'>, where 'X' is not a number. "
                    "\nconstant = ") &
            (_String((_String *)toStr())) & "\n'X' = " &
            (_String((_String *)operand->toStr()));

        throw error_message;
      }
    } else {
      _String error_message =
          _String("<'constant' operation 'null'>, where constant = ") &
          (_String((_String *)toStr()));

      throw error_message;
    }
    return new _MathObject;
  }

  /**
   * @brief A template function to check the type of two operands and compute the result of a ternary operation.
   *
   * @tparam T The type of the functor to apply.
   * @param operand The first operand to compute.
   * @param operand2 The second operand to compute.
   * @param functor The functor to apply.
   * @param cache The cache to use for the result.
   * @return The result of the computation.
   */
  template <class T>
  HBLObjectRef _check_type_and_compute_3(HBLObjectRef operand,
                                         HBLObjectRef operand2, T functor,
                                         HBLObjectRef cache) {
    if (operand && operand2 && operand->ObjectClass() == NUMBER &&
        operand2->ObjectClass() == NUMBER) {
      return _returnConstantOrUseCache(
          functor(Value(), ((_Constant *)operand)->Value(),
                  ((_Constant *)operand2)->Value()),
          cache);
    }

    _String err_msg(
        "Not a numeric 'X' type in a <'constant' operation 'X'> call");
    throw(err_msg);
    return new _MathObject;
  }

public:
  /**
   * @brief Construct a new _Constant object from a floating point value.
   *
   * @param v The value of the constant.
   */
  _Constant(hyFloat);
  /**
   * @brief Construct a new _Constant object from a string representation.
   *
   * @param s The string representation of the constant.
   */
  _Constant(_String &);
  /**
   * @brief Construct a new _Constant object with a default value of 0.0.
   */
  _Constant(void);
  /**
   * @brief Destroy the _Constant object
   *
   */
  ~_Constant(void) {}

  /**
   * @brief Add a value to this constant.
   *
   * @param o The value to add. Must be a `_Constant` object.
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the sum.
   */
  virtual HBLObjectRef Add(HBLObjectRef, HBLObjectRef cache = nil);
  /**
   * @brief Subtract a value from this constant.
   *
   * @param o The value to subtract. Must be a `_Constant` object.
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the difference.
   */
  virtual HBLObjectRef Sub(HBLObjectRef, HBLObjectRef cache = nil);
  /**
   * @brief Negate this constant.
   *
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the negated value.
   */
  virtual HBLObjectRef Minus(HBLObjectRef cache = nil);
  /**
   * @brief Sum this constant with another value.
   *
   * @param cache An optional cache to store the result.
   * @return This object.
   */
  virtual HBLObjectRef Sum(HBLObjectRef cache = nil);
  /**
   * @brief Multiply this constant by a value.
   *
   * @param o The value to multiply by. Must be a `_Constant` object.
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the product.
   */
  virtual HBLObjectRef Mult(HBLObjectRef, HBLObjectRef cache = nil);
  /**
   * @brief Divide this constant by a value.
   *
   * @param o The value to divide by. Must be a `_Constant` object.
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the quotient.
   */
  virtual HBLObjectRef Div(HBLObjectRef, HBLObjectRef cache = nil);
  /**
   * @brief Left divide this constant by a value.
   *
   * @param o The value to divide by. Must be a `_Constant` object.
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the quotient.
   */
  virtual HBLObjectRef lDiv(HBLObjectRef, HBLObjectRef cache = nil);
  /**
   * @brief Long divide this constant by a value.
   *
   * @param o The value to divide by. Must be a `_Constant` object.
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the quotient.
   */
  virtual HBLObjectRef longDiv(HBLObjectRef, HBLObjectRef cache = nil);
  /**
   * @brief Raise this constant to a power.
   *
   * @param o The power to raise to. Must be a `_Constant` object.
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the result of the exponentiation.
   */
  virtual HBLObjectRef Raise(HBLObjectRef, HBLObjectRef cache = nil);
  /**
   * @brief Check if this constant is equal to another value.
   *
   * @param o The value to compare with. Must be a `_Constant` object.
   * @return true if the values are equal, false otherwise.
   */
  virtual bool Equal(HBLObjectRef);
  /**
   * @brief Get the absolute value of this constant.
   *
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the absolute value.
   */
  virtual HBLObjectRef Abs(HBLObjectRef cache = nil);
  /**
   * @brief Get the sine of this constant.
   *
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the sine.
   */
  virtual HBLObjectRef Sin(HBLObjectRef cache = nil);
  /**
   * @brief Get the cosine of this constant.
   *
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the cosine.
   */
  virtual HBLObjectRef Cos(HBLObjectRef cache = nil);
  /**
   * @brief Get the tangent of this constant.
   *
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the tangent.
   */
  virtual HBLObjectRef Tan(HBLObjectRef cache = nil);
  /**
   * @brief Get the exponential of this constant.
   *
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the exponential.
   */
  virtual HBLObjectRef Exp(HBLObjectRef cache = nil);
  /**
   * @brief Get the natural logarithm of this constant.
   *
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the natural logarithm.
   */
  virtual HBLObjectRef Log(HBLObjectRef cache = nil);
  /**
   * @brief Get the square root of this constant.
   *
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the square root.
   */
  virtual HBLObjectRef Sqrt(HBLObjectRef cache = nil);
  /**
   * @brief Get the current time.
   *
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the current time.
   */
  virtual HBLObjectRef Time(HBLObjectRef cache = nil);
  /**
   * @brief Get the arctangent of this constant.
   *
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the arctangent.
   */
  virtual HBLObjectRef Arctan(HBLObjectRef cache = nil);
  /**
   * @brief Get the gamma function of this constant.
   *
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the gamma function.
   */
  virtual HBLObjectRef Gamma(HBLObjectRef cache = nil);
  /**
   * @brief Get the natural logarithm of the gamma function of this constant.
   *
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the natural logarithm of the gamma function.
   */
  virtual HBLObjectRef
  LnGamma(HBLObjectRef cache = nil); /* <- added by afyp, February 8, 2007 */
  /**
   * @brief Convert this constant to a polynomial.
   *
   * @param cache An optional cache to store the result.
   * @return A new `_Polynomial` object representing this constant.
   */
  virtual HBLObjectRef ConvertToPolynomial(HBLObjectRef cache = nil);
  /**
   * @brief Get the beta function of this constant and another value.
   *
   * @param o The second parameter of the beta function. Must be a `_Constant` object.
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the beta function.
   */
  virtual HBLObjectRef Beta(HBLObjectRef, HBLObjectRef cache = nil);
  /**
   * @brief Get the minimum of this constant and another value.
   *
   * @param o The value to compare with. Must be a `_Constant` object.
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the minimum of the two values.
   */
  virtual HBLObjectRef Min(HBLObjectRef, HBLObjectRef cache = nil);
  /**
   * @brief Get the maximum of this constant and another value.
   *
   * @param o The value to compare with. Must be a `_Constant` object.
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the maximum of the two values.
   */
  virtual HBLObjectRef Max(HBLObjectRef, HBLObjectRef cache = nil);
  /**
   * @brief Get the gamma distribution of this constant.
   *
   * @param o The shape parameter of the gamma distribution. Must be a `_Constant` object.
   * @param o2 The scale parameter of the gamma distribution. Must be a `_Constant` object.
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the gamma distribution.
   */
  virtual HBLObjectRef GammaDist(HBLObjectRef, HBLObjectRef,
                                 HBLObjectRef cache = nil);
  /**
   * @brief Get the cumulative gamma distribution of this constant.
   *
   * @param o The shape parameter of the gamma distribution. Must be a `_Constant` object.
   * @param o2 The scale parameter of the gamma distribution. Must be a `_Constant` object.
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the cumulative gamma distribution.
   */
  virtual HBLObjectRef CGammaDist(HBLObjectRef, HBLObjectRef,
                                  HBLObjectRef cache = nil);
  /**
   * @brief Get the incomplete beta function of this constant.
   *
   * @param o The second parameter of the incomplete beta function. Must be a `_Constant` object.
   * @param o2 The third parameter of the incomplete beta function. Must be a `_Constant` object.
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the incomplete beta function.
   */
  virtual HBLObjectRef IBeta(HBLObjectRef, HBLObjectRef,
                             HBLObjectRef cache = nil);
  /**
   * @brief Get the incomplete gamma function of this constant.
   *
   * @param o The second parameter of the incomplete gamma function. Must be a `_Constant` object.
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the incomplete gamma function.
   */
  virtual HBLObjectRef IGamma(HBLObjectRef, HBLObjectRef cache = nil);
  /**
   * @brief Get the cumulative chi-squared distribution of this constant.
   *
   * @param o The degrees of freedom. Must be a `_Constant` object.
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the cumulative chi-squared distribution.
   */
  virtual HBLObjectRef CChi2(HBLObjectRef, HBLObjectRef cache = nil);
  /**
   * @brief Get the inverse chi-squared distribution of this constant.
   *
   * @param o The degrees of freedom. Must be a `_Constant` object.
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the inverse chi-squared distribution.
   */
  virtual HBLObjectRef InvChi2(HBLObjectRef, HBLObjectRef cache = nil);
  /**
   * @brief Get the error function of this constant.
   *
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the error function.
   */
  virtual HBLObjectRef Erf(HBLObjectRef cache = nil);
  /**
   * @brief Get the standard normal cumulative distribution function of this constant.
   *
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the standard normal cumulative distribution function.
   */
  virtual HBLObjectRef ZCDF(HBLObjectRef cache = nil);
  /**
   * @brief Check if this constant is less than another value.
   *
   * @param o The value to compare with. Must be a `_Constant` object.
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing 1 if this constant is less than the other value, 0 otherwise.
   */
  virtual HBLObjectRef Less(HBLObjectRef, HBLObjectRef cache = nil);
  /**
   * @brief Check if this constant is greater than another value.
   *
   * @param o The value to compare with. Must be a `_Constant` object.
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing 1 if this constant is greater than the other value, 0 otherwise.
   */
  virtual HBLObjectRef Greater(HBLObjectRef, HBLObjectRef cache = nil);
  /**
   * @brief Check if this constant is less than or equal to another value.
   *
   * @param o The value to compare with. Must be a `_Constant` object.
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing 1 if this constant is less than or equal to the other value, 0 otherwise.
   */
  virtual HBLObjectRef LessEq(HBLObjectRef, HBLObjectRef cache = nil);
  /**
   * @brief Check if this constant is greater than or equal to another value.
   *
   * @param o The value to compare with. Must be a `_Constant` object.
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing 1 if this constant is greater than or equal to the other value, 0 otherwise.
   */
  virtual HBLObjectRef GreaterEq(HBLObjectRef, HBLObjectRef cache = nil);
  /**
   * @brief Check if this constant is equal to another value.
   *
   * @param o The value to compare with. Must be a `_Constant` object.
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing 1 if this constant is equal to the other value, 0 otherwise.
   */
  virtual HBLObjectRef AreEqual(HBLObjectRef, HBLObjectRef cache = nil);
  /**
   * @brief Check if this constant is not equal to another value.
   *
   * @param o The value to compare with. Must be a `_Constant` object.
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing 1 if this constant is not equal to the other value, 0 otherwise.
   */
  virtual HBLObjectRef NotEqual(HBLObjectRef, HBLObjectRef cache = nil);
  /**
   * @brief Logical AND of this constant and another value.
   *
   * @param o The value to AND with. Must be a `_Constant` object.
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the result of the logical AND.
   */
  virtual HBLObjectRef LAnd(HBLObjectRef, HBLObjectRef cache = nil);
  /**
   * @brief Logical OR of this constant and another value.
   *
   * @param o The value to OR with. Must be a `_Constant` object.
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the result of the logical OR.
   */
  virtual HBLObjectRef LOr(HBLObjectRef, HBLObjectRef cache = nil);
  /**
   * @brief Logical NOT of this constant.
   *
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the result of the logical NOT.
   */
  virtual HBLObjectRef LNot(HBLObjectRef cache = nil);
  /**
   * @brief Generate a random number.
   *
   * @param o The type of random number to generate.
   * @param cache An optional cache to store the result.
   * @return A new `_Constant` object representing the random number.
   */
  virtual HBLObjectRef Random(HBLObjectRef, HBLObjectRef cache = nil);
  /**
   * @brief Get the value of this constant.
   *
   * @return The value of this constant.
   */
  virtual hyFloat Value(void);
  /**
   * @brief Format this constant as a string.
   *
   * @param o The format string.
   * @param o2 The second format string.
   * @param cache An optional cache to store the result.
   * @return A new `_String` object representing the formatted constant.
   */
  virtual HBLObjectRef FormatNumberString(HBLObjectRef, HBLObjectRef,
                                          HBLObjectRef cache = nil);
  /**
   * @brief Compute the value of this constant.
   *
   * @return This object.
   */
  virtual HBLObjectRef Compute(void) { return this; };

  /**
   * @brief Initialize this constant.
   *
   * @param b Whether to initialize this constant.
   */
  virtual void Initialize(bool = true);
  /**
   * @brief Duplicate this constant.
   *
   * @param brc The constant to duplicate.
   */
  virtual void Duplicate(BaseRefConst);
  /**
   * @brief Make this constant dynamic.
   *
   * @return A new `_Constant` object representing the dynamic constant.
   */
  virtual BaseRef makeDynamic(void) const;
  /**
   * @brief Convert this constant to a string.
   *
   * @param ul The format to use.
   * @return A new `_String` object representing the constant.
   */
  virtual BaseRef toStr(unsigned long = 0UL);

  /**
   * @brief Get the object class of this constant.
   *
   * @return The object class, which is `NUMBER`.
   */
  virtual unsigned long ObjectClass(void) const { return NUMBER; }

  /**
   * @brief Set the value of this constant.
   *
   * @param pl The value to set.
   */
  virtual void SetValue(hyFloat pl) { theValue = pl; }

  /**
   * @brief The new operator for the `_Constant` class.
   *
   * @param size The size of the object to allocate.
   * @return A pointer to the allocated memory.
   */
  void *operator new(size_t size);
  /**
   * @brief The delete operator for the `_Constant` class.
   *
   * @param p A pointer to the memory to deallocate.
   */
  void operator delete(void *p);

#ifndef _USE_EMSCRIPTEN_
  static _SimpleList free_slots;
  static unsigned char preallocated_buffer[];
#endif

public:
  hyFloat theValue;
};

/**
 * @brief Compute the natural logarithm of the gamma function.
 *
 * @param alpha The parameter of the gamma function.
 * @return The natural logarithm of the gamma function.
 */
hyFloat _ln_gamma(hyFloat alpha);
/**
 * @brief Generate a normally distributed random number.
 *
 * @return A normally distributed random number.
 */
hyFloat gaussDeviate(void);
/**
 * @brief Generate an exponentially distributed random number.
 *
 * @return An exponentially distributed random number.
 */
hyFloat exponDeviate(void);
/**
 * @brief Generate a gamma distributed random number.
 *
 * @param a The shape parameter of the gamma distribution.
 * @param scale The scale parameter of the gamma distribution.
 * @return A gamma distributed random number.
 */
hyFloat gammaDeviate(hyFloat a, hyFloat scale = 1.);
/**
 * @brief Generate a chi-squared distributed random number.
 *
 * @param df The degrees of freedom.
 * @return A chi-squared distributed random number.
 */
hyFloat chisqDeviate(double df);

#endif
