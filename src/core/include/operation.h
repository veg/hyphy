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

#ifndef     __OPERATION__
#define     __OPERATION__

#include "baseobj.h"
#include "list.h"
#include "hy_strings.h"
#include "mathobj.h"

extern  _List BuiltInFunctions;


class _Stack;
class _VariableContainer;
class _Formula;
//__________________________________________________________________________________
class   _Operation : public BaseObj
{

    friend class _Formula;
    friend class _Variable;
    friend class _VariableContainer;

public:
    _Operation  (void);
    _Operation  (_String&, const long);
    // construct the operation by its symbol and, if relevant -
    // number of operands
    _Operation  (const long,const long);

    _Operation  (bool, _String&, bool isG = false, _VariableContainer*  = nil, bool take_a_reference = false);
    // store a variable or a constant
    _Operation  (_PMathObj);
    // store a non-numeric constant

    virtual ~_Operation (void);

    virtual   BaseObj*      makeDynamic         (void);

    bool            Execute             (_Stack&, _VariableContainer* = nil, _String* errMsg = nil); //execute this operation
    // see the commend for _Formula::ExecuteFormula for the second argument
    virtual   void          StackDepth          (long&);

    bool            ExecutePolynomial   (_Stack&,_VariableContainer* nameSpace = nil, _String* errMsg = nil);
    virtual   BaseObj*      toStr               (void);    //convert the op to string

    virtual   void          Initialize          (void);
    virtual   void          Duplicate           (BaseRef);
    _String&    GetCode             (void) {
        return (opCode>-1)&&(numberOfTerms>=0)?*(_String*)BuiltInFunctions(opCode):empty;
    }
    long&       TheCode             (void) {
        return opCode;
    }
    virtual  bool           IsAVariable         (bool = true) ; // is this object a variable or not?
    virtual  bool           IsConstant          (void);         // does this object depend on any independent variables or not?
    virtual  bool           IsAFunctionCall          (void);       

    virtual  long           UserFunctionID      (void) {
        return numberOfTerms < 0 ? -numberOfTerms-1 : -1;
    };
    // return a non-neg number (function index) if this is a user function,
    // otherwise, return -1

    virtual  long           GetAVariable        (void) {    // return the index of the variable
        return theData>=-2?theData:-theData-3;
    }

    virtual  void           SetAVariable        (long d) {  // return the index of the variable
        theData=d;
    }

    virtual  bool           AssignmentVariable  (void) {
        return theData<-2;
    }

    virtual  bool           HasChanged          (void);

    virtual  void           SetTerms            (long d) {
        numberOfTerms=d;
    }

    virtual  _PMathObj      GetANumber          (void) {
        return theNumber;
    }

    virtual  void           SetNumber           (_PMathObj d) {
        theNumber=d;
    }

    long            GetNoTerms          (void) {
        return numberOfTerms;
    }
    long            PrecedenceLevel     (void);

    bool            CanResultsBeCached (_Operation *, bool exp_only = false);


    virtual bool            EqualOp             (_Operation*);

protected:

    bool        ReportOperationExecutionError ( _String, _String*);

    long        opCode;         // internal operation code
    long        numberOfTerms,  // 1 - unary, 2 - binary, etc
                theData;
    _PMathObj   theNumber;
};

#endif
