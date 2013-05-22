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

#ifndef __OPERATION__
#define __OPERATION__

#include "baseobj.h"
#include "list.h"
#include "hy_strings.h"
#include "mathobj.h"


#define _HY_OPERATION_INVALID_REFERENCE (-1L)


#define _HY_OPERATION_NOOP   0L  
// this operation does nothing
#define _HY_OPERATION_VALUE  1L  
// this operation contains an object to be pushed on the stack
#define _HY_OPERATION_VAR    2L  
// this operation contains a reference to a variable whose value will be pushed on the stack
#define _HY_OPERATION_REF    3L
// this operation contains a reference to a string variable whose value will
// used to look up another variable whose value will be pushed on the stack
#define _HY_OPERATION_BUILTIN 4L
// this operation refers to a built-in function or operation
#define _HY_OPERATION_FUNCTION_CALL 5L
// this operation will call an HBL function
#define _HY_OPERATION_DEFERRED_FUNCTION_CALL 6L
// this operation contains a reference to an HBL function ID
// whose name will be looked up and bound at the time of first call


extern _List BuiltInFunctions;

class _Stack;
class _VariableContainer;
class _Formula;
//__________________________________________________________________________________
class _Operation : public BaseObj {

  friend class _Formula;
  friend class _Variable;
  friend class _VariableContainer;

protected:
    
    bool ReportOperationExecutionError(_String, _String *);
    
    long          operationKind;
    // what KIND of an operation is this operation
    // one of the _HY_OPERATION_ #defines
    
    
    long          reference;
    long          attribute;
    _PMathObj     payload;
    
/*
 
operationKind               |  reference                        | attribute                         | payload 
----------------------------+-----------------------------------+-----------------------------------+---------+
_HY_OPERATION_NOOP          |  _HY_OPERATION_INVALID_REFERENCE  | _HY_OPERATION_INVALID_REFERENCE   | NULL        
_HY_OPERATION_VALUE         |  _HY_OPERATION_INVALID_REFERENCE  | _HY_OPERATION_INVALID_REFERENCE   | object to push on stack
_HY_OPERATION_VAR           |  index of the variable            | _HY_OPERATION_INVALID_REFERENCE   | NULL
_HY_OPERATION_REF           |  index of the variable            | _HY_OPERATION_INVALID_REFERENCE   | NULL
_HY_OPERATION_BUILTIN       |  opCode (e.g. HY_OP_CODE_ADD)     | number of terms to consume from   | NULL
                            |                                   | the stack                         |
_HY_OPERATION_FUNCTION_CALL |  index of the function to call    | number of terms to consume from   | NULL OR
                            |                                   | the stack                         | named argument list
_HY_OPERATION_DEFERRED_     |  _HY_OPERATION_INVALID_REFERENCE  | number of terms to consume from   | function id
FUNCTION_CALL               |                                   | the stack                         | 

 
*/
   
public:
  _Operation(void);
  _Operation(_String &, const long);
  // construct the operation by its symbol and, if relevant -
  // number of operands
  _Operation(const long, const long);

  _Operation(bool, _String &, bool isG = false, _VariableContainer * = nil,
             bool take_a_reference = false);
  // store a variable or a constant
    
  _Operation(_PMathObj);
  virtual void Duplicate(BaseRef);
  
  #ifdef __NEW_GRAMMAR__
     bool ResolveDeferredAction (void);
  #endif    
  
  // store a non-numeric constant

  virtual ~_Operation(void);

  virtual void    Initialize(void);
  virtual BaseObj *makeDynamic(void);

  bool Execute(_Stack &, _VariableContainer * = nil,
               _String *errMsg = nil); //execute this operation
  // see the commend for _Formula::ExecuteFormula for the second argument
  bool ExecutePolynomial(_Stack &, _VariableContainer *nameSpace = nil,
                           _String *errMsg = nil);
  virtual void StackDepth(long &);

  virtual BaseObj *toStr(void); //convert the op to string

  _String &GetCode(void) {
    return (opCode > -1) && (numberOfTerms >= 0)
               ? *(_String *)BuiltInFunctions(opCode)
               : empty;
  }
  long &TheCode(void) { return opCode; }
  virtual bool IsAVariable(bool = true); // is this object a variable or not?
  virtual bool IsConstant(
      void); // does this object depend on any independent variables or not?
  virtual bool IsAFunctionCall(void);

  virtual long UserFunctionID(void) {
    return numberOfTerms < 0 ? -numberOfTerms - 1 : -1;
  }
  ;

  // return a non-neg number (function index) if this is a user function,
  // otherwise, return HY_NOT_FOUND

  virtual long GetAVariable(void) { // return the index of the variable
    return theData >= -2 ? theData : -theData - 3;
  }

  virtual void SetAVariable(long d) { // return the index of the variable
    theData = d;
  }

  virtual bool AssignmentVariable(void) { return theData < -2; }

  virtual bool HasChanged(void);

  virtual void SetTerms(long d) { numberOfTerms = d; }

  virtual _PMathObj GetANumber(void) { return theNumber; }

  virtual void SetNumber(_PMathObj d) { theNumber = d; }

  long GetNoTerms(void) { return numberOfTerms; }
  long PrecedenceLevel(void);

  bool CanResultsBeCached(_Operation *, bool exp_only = false);

  virtual bool EqualOp(_Operation *);

};

#endif
