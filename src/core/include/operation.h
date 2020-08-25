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

#ifndef     __OPERATION__
#define     __OPERATION__

#include "baseobj.h"
#include "list.h"
#include "trie.h"
#include "hy_strings.h"
#include "mathobj.h"
#include "global_things.h"

extern  _List BuiltInFunctions;


class _Stack;
class _VariableContainer;
class _Variable;
class _Formula;

_Variable * FetchVar (long, unsigned long = HY_ANY_OBJECT);

//__________________________________________________________________________________
class   _Operation : public BaseObj
{

    friend class    _Formula;
    friend class    _Variable;
    friend class    _VariableContainer;
    
protected:
    long           opCode;         // internal operation code
    long           numberOfTerms,  // 1 - unary, 2 - binary, etc
                   theData;
    
    HBLObjectRef   theNumber,
                   cachedResult;


public:
    _Operation  (void);
    _Operation  (_String const&, const long);
    // construct the operation by its symbol and, if relevant -
    // number of operands
    _Operation  (const long,const long);
    _Operation  (const _Operation&);

    _Operation  (bool, _String&, bool isG = false, _VariableContainer const*  = nil, bool take_a_reference = false);
    // store a variable or a constant
    _Operation  (HBLObjectRef);
    // store a non-numeric constant
    _Operation  (_Variable const&);
    // create an operation that references a variable


    virtual ~_Operation (void);

    virtual   BaseObj*      makeDynamic         (void) const;
    virtual   void          Duplicate           (BaseRefConst);
    void      operator = (_Operation const&);

    bool            Execute             (_Stack&, _VariableContainer const* = nil, _String* errMsg = nil, bool canCache = false); //execute this operation
    // see the commend for _Formula::ExecuteFormula for the second argument
    virtual   void          StackDepth          (long&);

    bool            ExecutePolynomial   (_Stack&,_VariableContainer* nameSpace = nil, _String* errMsg = nil);
    virtual   BaseObj*      toStr               (unsigned long = 0UL);    //convert the op to string

    virtual   void          Initialize          (bool = true);
    const _String&    GetCode             (void) {
        return (opCode>-1)&&(numberOfTerms>=0)?*(_String*)BuiltInFunctions(opCode): hy_global::kEmptyString;
    }
    long&       TheCode             (void) {
        return opCode;
    }
    virtual  bool           IsAVariable         (bool = true) ; // is this object a variable or not?
    virtual  bool           IsConstant          (bool strict = false);         // does this object depend on any independent variables or not?
    virtual  bool           IsHBLFunctionCall   (void) const;
    virtual  long           GetHBLFunctionID    (void) const;

    virtual  long           UserFunctionID      (void) {
        return numberOfTerms < 0 ? opCode : -1;
    };
    // return a non-neg number (function index) if this is a user function,
    // otherwise, return -1

    virtual  long           GetAVariable        (void) const{
        if (theData >= -1) {
            return theData;
        }
        if (theData < -2) {
            return -theData - 3;
        }
        return -numberOfTerms-1;
        // return the index of the variable
    }

    virtual  long           IsValueSubstitution        (void) const{
        return theData == -2;
    }

    virtual  void           SetAVariable        (long d) {  // return the index of the variable
        theData=d;
    }

    _Variable *             RetrieveVar         (void) const {
      if (theData != -1) {
        long var_idx = GetAVariable();
        if (var_idx >= 0) {
            return FetchVar(var_idx);
        } else {
            if (var_idx == -2) {
                return FetchVar (-numberOfTerms-1);
            }
        }
      }
      return nil;
    }

    virtual  bool           AssignmentVariable  (void) {
        return theData<-2;
    }

    virtual  bool           HasChanged          (void);

    virtual  void           SetTerms            (long d) {
        numberOfTerms=d;
    }

    long                    StackDepth          (void) const;

    virtual  HBLObjectRef      GetANumber          (void) {
        return theNumber;
    }

    virtual  void           SetNumber           (HBLObjectRef d) {
        theNumber=d;
    }

    long            GetNoTerms          (void) {
        return numberOfTerms;
    }
    long            PrecedenceLevel     (void);

    bool            CanResultsBeCached (_Operation *, bool exp_only = false);

    bool            IsConstantOfType   (const long type) const;


    virtual bool            EqualOp             (_Operation*);

    static  long            BinOpCode           (_String const &, long = -1);
    
    static bool             AreOpsInverse       (long op_code1, long op_code2) {
        long compose = op_code1 < op_code2 ? (op_code1 << 16) + (op_code2) : (op_code2 << 16) + (op_code1);
        return ListOfInverseOps.Find (compose) != kNotFound;
    }
    
    static      _SimpleList         ListOfInverseOps;
protected:


    bool        ReportOperationExecutionError ( _String, _String*);

};

#endif
