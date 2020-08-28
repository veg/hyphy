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

#include "parser.h"
#include "mathobj.h"
#include "global_things.h"

using namespace hy_global;


_MathObject  default_null_argument;

_MathObject* _MathObject:: _null_handler() {
    HandleApplicationError (kErrorStringNullOperand);
    this->AddAReference();
    return this;
}

_MathObject* _MathObject:: Add        (_MathObject*, _MathObject*)     {
    return _null_handler();
}
_MathObject* _MathObject:: Sub        (_MathObject*, _MathObject*)     {
    return _null_handler();
}
_MathObject* _MathObject:: Minus      (_MathObject*)             {
    return _null_handler();
}
_MathObject* _MathObject:: Sum        (_MathObject*)             {
    return _null_handler();
}
_MathObject* _MathObject:: Mult       (_MathObject*,_MathObject*)     {
    return _null_handler();
}
_MathObject* _MathObject:: Div        (_MathObject*,_MathObject*)     {
    return _null_handler();
}
_MathObject* _MathObject:: lDiv       (_MathObject*,_MathObject*)     {
    return _null_handler();
}
_MathObject* _MathObject:: longDiv    (_MathObject*,_MathObject*)     {
    return _null_handler();
}
_MathObject* _MathObject:: Raise      (_MathObject*,_MathObject*)     {
    return _null_handler();
}

bool _MathObject::         Equal      (_MathObject* rhs)     {
    if (rhs->ObjectClass() == HY_UNDEFINED) {
        return true;
    }
    return false;
    // null is equal to null, otherwise false\
    return false;
}
_MathObject* _MathObject:: Abs        (_MathObject*)             {
    return _null_handler();
}
_MathObject* _MathObject:: Sin        (_MathObject*)             {
    return _null_handler();
}
_MathObject* _MathObject:: Cos        (_MathObject*)             {
    return _null_handler();
}
_MathObject* _MathObject:: Tan        (_MathObject*)             {
    return _null_handler();
}
_MathObject* _MathObject:: Exp        (_MathObject*)             {
    return _null_handler();
}
_MathObject* _MathObject:: Log        (_MathObject*)             {
    return _null_handler();
}
_MathObject* _MathObject:: Sqrt       (_MathObject*)             {
    return _null_handler();
}
_MathObject* _MathObject:: Gamma      (_MathObject*)             {
    return _null_handler();
}
_MathObject* _MathObject:: Erf        (_MathObject*)             {
    return _null_handler();
}
_MathObject* _MathObject:: LnGamma    (_MathObject*)             {
    return _null_handler();
}

_MathObject* _MathObject:: Beta       (_MathObject*,_MathObject*)     {
    return _null_handler();
}
_MathObject* _MathObject:: IGamma     (_MathObject*,_MathObject*)     {
    return _null_handler();
}
_MathObject* _MathObject:: CChi2      (_MathObject*,_MathObject*)     {
    return _null_handler();
}
_MathObject* _MathObject:: IBeta      (_MathObject*,_MathObject*,_MathObject*) {
    return _null_handler();
}
_MathObject* _MathObject:: Simplex    (_MathObject*)             {
    return _null_handler();
}

_MathObject* _MathObject:: Simplify    (_MathObject*)             {
    return _null_handler();
}

_MathObject* _MathObject:: Min        (_MathObject*,_MathObject*)     {
    return _null_handler();
}
_MathObject* _MathObject:: Max        (_MathObject*,_MathObject*)     {
    return _null_handler();
}
_MathObject* _MathObject:: InvChi2    (_MathObject*,_MathObject*)     {
    return _null_handler();
}
_MathObject* _MathObject:: ZCDF       (_MathObject*)             {
    return _null_handler();
}
_MathObject* _MathObject:: Time       (_MathObject*)             {
    return _null_handler();
}
_MathObject* _MathObject:: Arctan     (_MathObject*)             {
    return _null_handler();
}
_MathObject* _MathObject:: Less       (_MathObject*,_MathObject*)     {
    return _null_handler();
}
_MathObject* _MathObject:: Random     (_MathObject*,_MathObject*)     {
    return _null_handler();
}
_MathObject* _MathObject:: Greater    (_MathObject*,_MathObject*)     {
    return _null_handler();
}
_MathObject* _MathObject:: LessEq     (_MathObject*,_MathObject*)     {
    return _null_handler();
}
_MathObject* _MathObject:: GreaterEq  (_MathObject*,_MathObject*)     {
    return _null_handler();
}
_MathObject* _MathObject:: AreEqual   (_MathObject*,_MathObject*)     {
    return _null_handler();
}
_MathObject* _MathObject:: NotEqual   (_MathObject*,_MathObject*)     {
    return _null_handler();
}
_MathObject* _MathObject:: LAnd       (_MathObject*,_MathObject*)     {
    return _null_handler();
}
_MathObject* _MathObject:: LOr        (_MathObject*,_MathObject*)     {
    return _null_handler();
}
_MathObject* _MathObject:: GammaDist  (_MathObject*,_MathObject*,_MathObject*) {
    return _null_handler();
}
_MathObject* _MathObject:: CGammaDist (_MathObject*,_MathObject*,_MathObject*) {
    return _null_handler();
}
_MathObject* _MathObject:: LNot       (_MathObject*)             {
    return _null_handler();
}
_MathObject* _MathObject:: TipCount   (_MathObject*)             {
    return _null_handler();
}
_MathObject* _MathObject:: BranchCount (_MathObject*)            {
    return _null_handler();
}
_MathObject* _MathObject:: TipName     (_MathObject*,_MathObject*)    {
    return _null_handler();
}
_MathObject* _MathObject:: BranchName  (_MathObject*,_MathObject*)    {
    return _null_handler();
}
_MathObject* _MathObject:: BranchLength(_MathObject*,_MathObject*)    {
    return _null_handler();
}
_MathObject* _MathObject:: RerootTree  (_MathObject*,_MathObject*)    {
    return _null_handler();
}
_MathObject* _MathObject:: TEXTreeString(_MathObject*,_MathObject*) const {
    HandleApplicationError (kErrorStringNullOperand);
    return new _MathObject;
}

_MathObject* _MathObject:: PlainTreeString(_MathObject*,_MathObject*,_MathObject*) {
    return _null_handler();
}
_MathObject* _MathObject:: FormatNumberString (_MathObject*,_MathObject*,_MathObject*) {
    return _null_handler();
}
hyFloat _MathObject::   Value (void)              {
    HandleApplicationError (kErrorStringNullOperand);
    return 0.0;
}


_MathObject* _MathObject:: _extract_argument (_List * arguments, unsigned long index, bool fill_in) const {
  if (arguments && index < arguments->lLength) {
    return (_MathObject*)arguments->GetItem(index);
  }
  return fill_in ? &default_null_argument : nil;
}


HBLObjectRef _returnConstantOrUseCache (hyFloat value, HBLObjectRef cache) {
    if (cache && cache->ObjectClass() == NUMBER) {
        ((_Constant*)cache)->theValue = value;
        return cache;
    }
    return new _Constant (value);
}

//SW: This calls the function with the opcode after it's been parsed
HBLObjectRef _MathObject::ExecuteSingleOp (long opCode, _List* arguments, _hyExecutionContext* context, HBLObjectRef cache) {

    switch (opCode) { // first check operations without arguments
      case HY_OP_CODE_NOT: // !
        return LNot(cache);
      case HY_OP_CODE_ABS: // Abs
        return Abs(cache);
      case HY_OP_CODE_ARCTAN: // Arctan
        return Arctan(cache);
      case HY_OP_CODE_COS: // Cos
        return Cos(cache);
      case HY_OP_CODE_COLUMNS: // Columns
      case HY_OP_CODE_ROWS: // Rows
        return _returnConstantOrUseCache (0., cache);
      case HY_OP_CODE_ERF: // Erf
        return Erf(cache);
      case HY_OP_CODE_EVAL:
        return (HBLObjectRef)Compute()->makeDynamic();
      case HY_OP_CODE_EXP: // Exp
        return Exp(cache);
      case HY_OP_CODE_GAMMA: // Gamma
        return Gamma(cache);
      case HY_OP_CODE_LNGAMMA: // LnGamma
        return LnGamma(cache);
      case HY_OP_CODE_LOG: // Log
        return Log(cache);
      case HY_OP_CODE_MACCESS:
        return new _MathObject; // indexing None returns None
      case HY_OP_CODE_SIMPLEX: // Simplex
        return Simplex(cache);
      case HY_OP_CODE_SIN: // Sin
        return Sin(cache);
      case HY_OP_CODE_SQRT: // Sqrt
        return Sqrt(cache);
      case HY_OP_CODE_TAN: // Tan
        return Tan(cache);
      case HY_OP_CODE_TIME: // Time
        return Time(cache);
      case HY_OP_CODE_TYPE: // Type
        return Type(cache);
      case HY_OP_CODE_ZCDF: // ZCDF
        return ZCDF(cache);
    }
  
    _MathObject * arg0 = _extract_argument (arguments, 0UL, false);
  
    switch (opCode) { // next check operations without arguments or with one argument
      case HY_OP_CODE_ADD: // +
        if (arg0)
          return Add(arg0,cache);
        return Sum (cache);
      case HY_OP_CODE_SUB: // -
        if (arg0)
          return Sub(arg0,cache);
        return Minus(cache);
        break;
    }
  
    if (arg0) {
      switch (opCode) { // operations that require exactly one argument
          case HY_OP_CODE_NEQ: // !=
            if (ObjectClass () == HY_UNDEFINED) {
              if (arg0->ObjectClass () == HY_UNDEFINED)
                return _returnConstantOrUseCache (0., cache);
                //return new HY_CONSTANT_FALSE;
              else
                return _returnConstantOrUseCache (1., cache);
                //return new HY_CONSTANT_TRUE;
            }
            if (arg0->ObjectClass() == NUMBER)
              return NotEqual(arg0,cache);
            return new HY_CONSTANT_TRUE;
          case HY_OP_CODE_IDIV: // $
            return longDiv(arg0,cache);
          case HY_OP_CODE_MOD: // %
            return lDiv(arg0,cache);
          case HY_OP_CODE_AND: // &&
            return LAnd(arg0,cache);
          case HY_OP_CODE_MUL: // *
            return Mult(arg0,cache);
          case HY_OP_CODE_DIV: // /
            return Div(arg0,cache);
          case HY_OP_CODE_LESS: // <
            return Less(arg0,cache);
          case HY_OP_CODE_LEQ: // <=
            return LessEq(arg0,cache);
          case HY_OP_CODE_EQ: // ==
          {
            if (ObjectClass () == HY_UNDEFINED) {
              if (arg0->ObjectClass () == HY_UNDEFINED)
                return _returnConstantOrUseCache (1., cache);
                //return new HY_CONSTANT_TRUE;
              else
                return _returnConstantOrUseCache (0., cache);
                //return new HY_CONSTANT_FALSE;
            }
            if (arg0->ObjectClass() == NUMBER)
              return AreEqual(arg0,cache);
            return new HY_CONSTANT_FALSE;
          }
            break;
          case HY_OP_CODE_GREATER: // >
            return Greater(arg0,cache);
          case HY_OP_CODE_GEQ: // >=
            return GreaterEq(arg0,cache);
          case HY_OP_CODE_BETA: // Beta
            return Beta(arg0,cache);
          case HY_OP_CODE_CCHI2: // CChi2
            return CChi2(arg0,cache);
          case HY_OP_CODE_IGAMMA: // IGamma
            return IGamma(arg0,cache);
          case HY_OP_CODE_INVCHI2: // InvChi2
            return InvChi2(arg0,cache);
          case HY_OP_CODE_MAX: // Max
            return Max(arg0,cache);
          case HY_OP_CODE_MIN: // Min
            return Min(arg0,cache);
          case HY_OP_CODE_RANDOM: // Random
            return Random(arg0,cache);
          case HY_OP_CODE_POWER: // ^
            return Raise(arg0,cache);
          case HY_OP_CODE_OR: // ||
            return LOr(arg0,cache);
      }
      
      _MathObject * arg1 = _extract_argument (arguments, 1UL, false);
      
      if (arg1) {
        /** ops that require exactly TWO arguments */
        
        switch (opCode) {
          case HY_OP_CODE_CGAMMADIST: // CGammaDist
            return CGammaDist(arg0,arg1,cache);
          case HY_OP_CODE_FORMAT: // Format
            return FormatNumberString(arg0,arg1,cache);
          case HY_OP_CODE_GAMMADIST: // GammaDist
            return GammaDist(arg0,arg1,cache);
          case HY_OP_CODE_IBETA: // IBeta
            return IBeta(arg0,arg1,cache);
        }
      }
 
    }
  
    switch (opCode) {
      case HY_OP_CODE_CGAMMADIST: // CGammaDist
      case HY_OP_CODE_FORMAT: // Format
      case HY_OP_CODE_GAMMADIST: // GammaDist
      case HY_OP_CODE_IBETA: // IBeta
      case HY_OP_CODE_NEQ: // !=
      case HY_OP_CODE_IDIV: // $
      case HY_OP_CODE_MOD: // %
      case HY_OP_CODE_AND: // &&
      case HY_OP_CODE_MUL: // *
      case HY_OP_CODE_DIV: // /
      case HY_OP_CODE_LESS: // <
      case HY_OP_CODE_LEQ: // <=
      case HY_OP_CODE_EQ: // ==
      case HY_OP_CODE_GREATER: // >
      case HY_OP_CODE_GEQ: // >=
      case HY_OP_CODE_BETA: // Beta
      case HY_OP_CODE_CCHI2: // CChi2
      case HY_OP_CODE_IGAMMA: // IGamma
      case HY_OP_CODE_INVCHI2: // InvChi2
      case HY_OP_CODE_MAX: // Max
      case HY_OP_CODE_MIN: // Min
      case HY_OP_CODE_RANDOM: // Random
      case HY_OP_CODE_POWER: // ^
      case HY_OP_CODE_OR: // ||
        WarnWrongNumberOfArguments (this, opCode,context, arguments);
        break;
      default:
        WarnNotDefined (this, opCode,context);
    }
  
    return new _MathObject;
}

//__________________________________________________________________________________

BaseRef _MathObject::makeDynamic (void) const {
    return new _MathObject;
}

//__________________________________________________________________________________

void _MathObject::Duplicate (BaseRefConst) {
}



//__________________________________________________________________________________

  //SW: Why do we need a string for the type?
HBLObjectRef _MathObject::Type (HBLObjectRef cache) {
  
  static const _FString kNumber       ("Number", false);
  static const _FString kMatrix       ("Matrix", false);
  static const _FString kContainer    ("Container", false);
  static const _FString kTreeNode     ("TreeNode", false);
  static const _FString kTree         ("Tree", false);
  static const _FString kString       ("String", false);
  static const _FString kAssociativeList
                                      ("AssociativeList", false);
  static const _FString kTopology     ("Topology", false);
  static const _FString kPolynomial   ("Polynomial", false);
  static const _FString kUnknown      ("Unknown", false);
  
  switch (ObjectClass()) {
      
    case NUMBER:
      return new _FString (kNumber);
    case MATRIX:
      return new _FString (kMatrix);
    case CONTAINER:
      return new _FString (kContainer);
    case TREE_NODE:
      return new _FString (kTreeNode);
    case TREE:
      return new _FString (kTree);
    case STRING:
      return new _FString (kString);
    case ASSOCIATIVE_LIST:
      return new _FString (kAssociativeList);
    case TOPOLOGY:
      return new _FString (kTopology);
    case POLYNOMIAL:
      return new _FString (kPolynomial);
      
  }
  
  return new _FString (kUnknown);
}
