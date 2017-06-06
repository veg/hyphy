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

_MathObject* _MathObject:: Add        (_MathObject*)     {
    HandleApplicationError (kErrorStringNullOperand); 
    return this;
}
_MathObject* _MathObject:: Sub        (_MathObject*)     {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: Minus      (void)             {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: Sum        (void)             {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: Mult       (_MathObject*)     {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: Div        (_MathObject*)     {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: lDiv       (_MathObject*)     {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: longDiv    (_MathObject*)     {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: Raise      (_MathObject*)     {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
void _MathObject::         Assign     (_MathObject*)     {
    HandleApplicationError (kErrorStringNullOperand);
}
bool _MathObject::         Equal      (_MathObject*)     {
    HandleApplicationError (kErrorStringNullOperand);
    return false;
}
_MathObject* _MathObject:: Abs        (void)             {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: Sin        (void)             {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: Cos        (void)             {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: Tan        (void)             {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: Exp        (void)             {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: Log        (void)             {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: Sqrt       (void)             {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: Gamma      (void)             {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: Erf        (void)             {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: LnGamma    (void)             {
    HandleApplicationError (kErrorStringNullOperand);    // <-- added by afyp, February 7, 2007
    return this;
}
_MathObject* _MathObject:: Beta       (_MathObject*)     {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: IGamma     (_MathObject*)     {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: CChi2      (_MathObject*)     {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: IBeta      (_MathObject*,_MathObject*) {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: Simplex    (void)             {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}

_MathObject* _MathObject:: Simplify    (void)             {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}

_MathObject* _MathObject:: Min        (_MathObject*)     {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: Max        (_MathObject*)     {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: InvChi2    (_MathObject*)     {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: ZCDF       (void)             {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: Time       (void)             {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: Arctan     (void)             {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: Less       (_MathObject*)     {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: Random     (_MathObject*)     {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: Greater    (_MathObject*)     {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: LessEq     (_MathObject*)     {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: GreaterEq  (_MathObject*)     {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: AreEqual   (_MathObject*)     {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: NotEqual   (_MathObject*)     {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: LAnd       (_MathObject*)     {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: LOr        (_MathObject*)     {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: GammaDist  (_MathObject*,_MathObject*) {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: CGammaDist (_MathObject*,_MathObject*) {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: LNot       (void)             {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: TipCount   (void)             {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: BranchCount (void)            {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: TipName     (_MathObject*)    {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: BranchName  (_MathObject*)    {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: BranchLength(_MathObject*)    {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: RerootTree  (_MathObject*)    {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: TEXTreeString(_MathObject*) const {
    HandleApplicationError (kErrorStringNullOperand);
    return new _MathObject;
}

_MathObject* _MathObject:: PlainTreeString(_MathObject*,_MathObject*) {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
_MathObject* _MathObject:: FormatNumberString (_MathObject*,_MathObject*) {
    HandleApplicationError (kErrorStringNullOperand);
    return this;
}
hy_float _MathObject::   Value (void)              {
    HandleApplicationError (kErrorStringNullOperand);
    return 0.0;
}


_MathObject* _MathObject:: _extract_argument (_List * arguments, unsigned long index, bool fill_in) const {
  if (arguments && index < arguments->lLength) {
    return (_MathObject*)arguments->GetItem(index);
  }
  return fill_in ? &default_null_argument : nil;
}


//SW: This calls the function with the opcode after it's been parsed
_PMathObj _MathObject::ExecuteSingleOp (long opCode, _List* arguments, _hyExecutionContext* context) {

    switch (opCode) { // first check operations without arguments
      case HY_OP_CODE_NOT: // !
        return LNot();
      case HY_OP_CODE_ABS: // Abs
        return Abs();
      case HY_OP_CODE_ARCTAN: // Arctan
        return Arctan();
      case HY_OP_CODE_COS: // Cos
        return Cos();
      case HY_OP_CODE_COLUMNS: // Columns
      case HY_OP_CODE_ROWS: // Rows
        return new _Constant (0.0);
      case HY_OP_CODE_ERF: // Erf
        return Erf();
      case HY_OP_CODE_EVAL:
        return (_PMathObj)Compute()->makeDynamic();
      case HY_OP_CODE_EXP: // Exp
        return Exp();
      case HY_OP_CODE_GAMMA: // Gamma
        return Gamma();
      case HY_OP_CODE_LNGAMMA: // LnGamma
        return LnGamma();
      case HY_OP_CODE_LOG: // Log
        return Log();
      case HY_OP_CODE_MACCESS:
        return new _MathObject; // indexing None returns None
      case HY_OP_CODE_SIMPLEX: // Simplex
        return Simplex();
      case HY_OP_CODE_SIN: // Sin
        return Sin();
      case HY_OP_CODE_SQRT: // Sqrt
        return Sqrt();
      case HY_OP_CODE_TAN: // Tan
        return Tan();
      case HY_OP_CODE_TIME: // Time
        return Time();
      case HY_OP_CODE_TYPE: // Type
        return Type();
      case HY_OP_CODE_ZCDF: // ZCDF
        return ZCDF();
    }
  
    _MathObject * arg0 = _extract_argument (arguments, 0UL, false);
  
    switch (opCode) { // next check operations without arguments or with one argument
      case HY_OP_CODE_ADD: // +
        if (arg0)
          return Add(arg0);
        return Sum ();
      case HY_OP_CODE_SUB: // -
        if (arg0)
          return Sub(arg0);
        return Minus();
        break;
    }
  
    if (arg0) {
      switch (opCode) { // operations that require exactly one argument
          case HY_OP_CODE_NEQ: // !=
            if (ObjectClass () == HY_UNDEFINED) {
              if (arg0->ObjectClass () == HY_UNDEFINED)
                return new HY_CONSTANT_FALSE;
              else
                return new HY_CONSTANT_TRUE;
            }
            if (arg0->ObjectClass() == NUMBER)
              return NotEqual(arg0);
            return new HY_CONSTANT_TRUE;
          case HY_OP_CODE_IDIV: // $
            return longDiv(arg0);
          case HY_OP_CODE_MOD: // %
            return lDiv(arg0);
          case HY_OP_CODE_AND: // &&
            return LAnd(arg0);
          case HY_OP_CODE_MUL: // *
            return Mult(arg0);
          case HY_OP_CODE_DIV: // /
            return Div(arg0);
          case HY_OP_CODE_LESS: // <
            return Less(arg0);
          case HY_OP_CODE_LEQ: // <=
            return LessEq(arg0);
          case HY_OP_CODE_EQ: // ==
          {
            if (ObjectClass () == HY_UNDEFINED) {
              if (arg0->ObjectClass () == HY_UNDEFINED)
                return new HY_CONSTANT_TRUE;
              else
                return new HY_CONSTANT_FALSE;
            }
            if (arg0->ObjectClass() == NUMBER)
              return AreEqual(arg0);
            return new HY_CONSTANT_FALSE;
          }
            break;
          case HY_OP_CODE_GREATER: // >
            return Greater(arg0);
          case HY_OP_CODE_GEQ: // >=
            return GreaterEq(arg0);
          case HY_OP_CODE_BETA: // Beta
            return Beta(arg0);
          case HY_OP_CODE_CCHI2: // CChi2
            return CChi2(arg0);
          case HY_OP_CODE_IGAMMA: // IGamma
            return IGamma(arg0);
          case HY_OP_CODE_INVCHI2: // InvChi2
            return InvChi2(arg0);
          case HY_OP_CODE_MAX: // Max
            return Max(arg0);
          case HY_OP_CODE_MIN: // Min
            return Min(arg0);
          case HY_OP_CODE_RANDOM: // Random
            return Random(arg0);
          case HY_OP_CODE_POWER: // ^
            return Raise(arg0);
          case HY_OP_CODE_OR: // ||
            return LOr(arg0);
      }
      
      _MathObject * arg1 = _extract_argument (arguments, 1UL, false);
      
      if (arg1) {
        /** ops that require exactly TWO arguments */
        
        switch (opCode) {
          case HY_OP_CODE_CGAMMADIST: // CGammaDist
            return CGammaDist(arg0,arg1);
          case HY_OP_CODE_FORMAT: // Format
            return FormatNumberString(arg0,arg1);
          case HY_OP_CODE_GAMMADIST: // GammaDist
            return GammaDist(arg0,arg1);
          case HY_OP_CODE_IBETA: // IBeta
            return IBeta(arg0,arg1);
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
_PMathObj _MathObject::Type (void) {
    _FString * ts = new _FString();
    switch (ObjectClass()) {

    case NUMBER:
        *(ts->theString)="Number";
        break;
    case MATRIX:
        *(ts->theString)="Matrix";
        break;
    case CONTAINER:
        *(ts->theString)="Container";
        break;
    case TREE_NODE:
        *(ts->theString)="TreeNode";
        break;
    case TREE:
        *(ts->theString)="Tree";
        break;
    case STRING:
        *(ts->theString)="String";
        break;
    case ASSOCIATIVE_LIST:
        *(ts->theString)="AssociativeList";
        break;
    case TOPOLOGY:
        *(ts->theString)="Topology";
        break;
    case POLYNOMIAL:
        *(ts->theString)="Polynomial";
        break;
    default:
        *(ts->theString) = "Unknown";

    }

    return ts;
}
