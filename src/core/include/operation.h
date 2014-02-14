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
#include "trie.h"
#include "mathobj.h"


#define _HY_OPERATION_INVALID_REFERENCE (-1L)
#define _HY_OPERATION_TOGGLE   1L
#define _HY_OPERATION_MIN_PRECEDENCE (-1L)
#define _HY_OPERATION_MAX_PRECEDENCE 0x07FFFFFFF

#define _HY_OPERATION_DUMMY_ARGUMENT_PLACEHOLDER empty,


#define _HY_OPERATION_NOOP                0x000000L  
// this operation does nothing
#define _HY_OPERATION_VALUE               0x000001L  
// this operation contains an object to be pushed on the stack
#define _HY_OPERATION_VAR                 0x000002L  
// this operation contains a reference to a variable whose computed value will be pushed on the stack
#define _HY_OPERATION_VAR_OBJ             0x000004L  
// this operation contains a reference to a variable whose value will be pushed on the stack 
// e.g. for matrix/avl assignments
#define _HY_OPERATION_REF                 0x000008L
// this operation contains a reference to a string variable whose value will
// used to look up another variable whose computed value will be pushed on the stack
#define _HY_OPERATION_BUILTIN             0x000010L
// this operation refers to a built-in function or operation
#define _HY_OPERATION_FUNCTION_CALL          0x000020L
// this operation will call an HBL function
#define _HY_OPERATION_DEFERRED_FUNCTION_CALL 0x000040L
// this operation contains a reference to an HBL function ID
// whose name will be looked up and bound at the time of first call
#define _HY_OPERATION_DEFERRED_INLINE        0x000080L
// this operation represents an ident__ call, where the value of ident at
// the time of execute/deferral resolution is substituted as _HY_OPERATION_VALUE

#define _HY_OPERATION_SPARSE_MATRIX  0x000100L
// this operation constructs a sparse matrix at run time
// this for example, can handle matrices of variable dimensions
// in the new check syntax now/execute later framework

#define _HY_OPERATION_DICTIONARY      0x000200L
// this operation constructs a dictionary at run time

#define _HY_OPERATION_ASSIGNMENT_VALUE       0x000800L

// assign a value of the form LHS or RHS
// handles x = y
// x += y
// x -= y
// x *= y
// x /= y etc


#define _HY_OPERATION_ASSIGNMENT_EXPRESSION  0x001000L
#define _HY_OPERATION_ASSIGNMENT_BOUND       0x002000L

#define _HY_OPERATION_ASSIGNMENT_BOUND_LOWER 0x000000L
#define _HY_OPERATION_ASSIGNMENT_BOUND_UPPER 0x000001L


#define _HY_OPERATION_FAST_EXEC_VALUE        0x010000L
// the analog of _HY_OPERATION_VALUE for 'Simple' formulas
#define _HY_OPERATION_FAST_EXEC_VAR          0x020000L
// the analog of _HY_OPERATION_VAR for 'Simple' formulas
#define _HY_OPERATION_FAST_EXEC_VAR_OBJ      0x040000L
// the analog of _HY_OPERATION_VAR_OBJ for 'Simple' formulas
#define _HY_OPERATION_FAST_EXEC_BUILTIN      0x080000L
// the analog of _HY_OPERATION_BUILTIN for 'Simple' formulas
#define _HY_OPERATION_FAST_EXEC_BUILTIN_REF  0x100000L
// the analog of _HY_OPERATION_BUILTIN for 'Simple' formulas
// but using 'references' instead of 'values' to access arguments
// used by matrix[]

#define _HY_OPERATION_FAST_EXEC              (_HY_OPERATION_FAST_EXEC_VALUE | _HY_OPERATION_FAST_EXEC_VAR | _HY_OPERATION_FAST_EXEC_VAR_OBJ | _HY_OPERATION_FAST_EXEC_BUILTIN |_HY_OPERATION_FAST_EXEC_BUILTIN_REF)

#define _HY_OPERATION_OP_CLASS               (_HY_OPERATION_BUILTIN | _HY_OPERATION_FAST_EXEC_BUILTIN | _HY_OPERATION_FAST_EXEC_BUILTIN_REF)

extern _Trie BuiltInFunctions;

class _Stack;
class _VariableContainer;
class _Formula;
union _SimpleFormulaDatum;

#define _HY2OPERATION(X) (dynamic_cast<_Operation*>(X))

//__________________________________________________________________________________
class _Operation : public BaseObj {

  friend class _Formula;
  friend class _Variable;
  friend class _VariableContainer;

protected:
    
    bool ReportOperationExecutionError(_String, _String *) const;
    
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
                                                                  or _HY_OPERATION_TOGGLE
                                                                  to place the INDEX of the 
                                                                  variable on the stack 
                                                                  (for assignments)
                                                                  
_HY_OPERATION_VAR_OBJ       |  index of the variable            | same as _HY_OPERATION_VAR         | NULL

_HY_OPERATION_REF           |  index of the variable            | _HY_OPERATION_INVALID_REFERENCE   | NULL

_HY_OPERATION_BUILTIN       |  opCode (e.g. HY_OP_CODE_ADD)     | number of terms to consume from   | NULL
                            |                                   | the stack                         |

_HY_OPERATION_FUNCTION_CALL |  index of the function to call    | number of terms to consume from   | NULL OR
                            |                                   | the stack                         | named argument list

_HY_OPERATION_DEFERRED_     |  _HY_OPERATION_INVALID_REFERENCE  | number of terms to consume from   | function id
FUNCTION_CALL               |                                   | the stack                         | 

_HY_OPERATION_DEFERRED_     |  index of the variable            | _HY_OPERATION_INVALID_REFERENCE   | NULL
INLINE                      |                                   |                                   | 

_HY_OPERATION_SPARSE_       |  _HY_OPERATION_INVALID_REFERENCE  | _HY_OPERATION_INVALID_REFERENCE   | matrix spec object
MATRIX                      |                                   |                                   | see _Matrix::_Matrix  (_PMathObj)

_HY_OPERATION_DICTIONARY    |  _HY_OPERATION_INVALID_REFERENCE  | _HY_OPERATION_INVALID_REFERENCE   | matrix spec object
                            |                                   |                                   | see ::_AssociativeArray (_PMathObj)

_HY_OPERATION_ASSIGNMENT    |  _HY_OPERATION_INVALID_REFERENCE  | opCode (which operation to        | NULL
_VALUE                      |  or the number of arguments to    | perform on the LHS)               | 
                            |  consume for x[][] = y assignment |                                   |

_HY_OPERATION_ASSIGNMENT    |  _HY_OPERATION_INVALID_REFERENCE  |_HY_OPERATION_..._BOUND_LOWER or   | NULL
_BOUND                      |                                   |_HY_OPERATION_..._BOUND_UPPER      |

_HY_OPERATION_ASSIGNMENT    |  _HY_OPERATION_INVALID_REFERENCE  | _HY_OPERATION_INVALID_REFERENCE   | *_Formula to assign
_EXPRESSION                 |                                   |                                   | to the LHS

 

_HY_OPERATION_FAST_EXEC     |  _HY_OPERATION_INVALID_REFERENCE  | _HY_OPERATION_INVALID_REFERENCE   | object to push on stack
_VALUE                      |                                   |                                   | (must be a scalar)

_HY_OPERATION_FAST_EXEC     |  index of the variable            | index of the variable in the      | NULL
_VAR                        |                                   | 'compiled' _SimpleFormulaDatum*   | 

_HY_OPERATION_FAST_EXEC     |  index of the variable            | index of the variable in the      | NULL
_VAR_OBJ                    |                                   | 'compiled' _SimpleFormulaDatum*   | 

_HY_OPERATION_FAST_EXEC     |  opCode (e.g. HY_OP_CODE_ADD)     | number of terms to consume from   | the 'shortcut' function
BUILTIN                     |                                   | the stack                         | to call

_HY_OPERATION_FAST_EXEC     |  opCode (e.g. HY_OP_CODE_ADD)     | number of terms to consume from   | the 'shortcut' function
BUILTIN_REF                 |                                   | the stack                         | to call

 
*/
   
public:
  _Operation(void);
  
  _Operation(bool isUserFunction, _String &, const long);
  // construct the operation by its symbol and, if relevant -
  // number of operands
  _Operation(const long, const long);
  _Operation(const long, const long, const long, _PMathObj);

  // store a variable or a constant
  _Operation(_String &, bool, bool = false, _VariableContainer * = nil,
             bool take_a_reference = false, bool deferred_inline = false);
    
  _Operation(_PMathObj);
  _Operation(_Operation&);
  virtual void Duplicate(BaseRef);
  

  // store a non-numeric constant

  virtual ~_Operation(void);

  virtual void    Initialize(const long kind = _HY_OPERATION_NOOP, 
                             const long ref  = _HY_OPERATION_INVALID_REFERENCE,
                             const long attr = _HY_OPERATION_INVALID_REFERENCE,
                             const _PMathObj  load = NULL,
                             bool clean_load = false);
                             
  virtual BaseObj *makeDynamic(void);

  bool ExecuteBuiltIn(_Stack &, _hyExecutionContext* = _hyDefaultExecutionContext); //execute a built-in function

  bool ExecuteFunctionCall (_Stack &, _hyExecutionContext* = _hyDefaultExecutionContext); //execute an HBL function

  bool ExecuteBounds (_Stack &, _hyExecutionContext* = _hyDefaultExecutionContext); //execute an assignment 

  bool ExecuteAssignment (_Stack &, _hyExecutionContext* = _hyDefaultExecutionContext); //execute an assignment 
               
  bool Execute(_Stack &, _hyExecutionContext* = _hyDefaultExecutionContext); //execute this operation
  // see the commend for _Formula::ExecuteFormula for the second argument
  
  bool      ExecuteFast (_SimpleFormulaDatum *stack,
                    _SimpleFormulaDatum *varValues,
                    long& stackTop,
                    _String *errMsg = nil); //execute this operation

  void ResolveDeferredAction (_hyExecutionContext*);
  
  bool ExecutePolynomial(_Stack &, _VariableContainer *nameSpace = nil,
                           _String *errMsg = nil);
  virtual void StackDepth(long &) const;

  virtual BaseObj *toStr(void); //convert the op to string

  _String GetCode(void);
  
  long &TheCode(void) { return operationKind == _HY_OPERATION_BUILTIN ? reference : attribute; }
  virtual bool IsAVariable(bool = true) const; // is this object a variable or not?
  virtual bool IsConstant(
      void) const; // does this object depend on any independent variables or not?
  virtual bool IsAFunctionCall(void) const;

  virtual long UserFunctionID(void) const {
    return operationKind == _HY_OPERATION_FUNCTION_CALL ? reference : _HY_OPERATION_INVALID_REFERENCE;
  }
  

  // return a non-neg number (function index) if this is a user function,
  // otherwise, return HY_NOT_FOUND

  virtual long GetAVariable(void) const { // return the index of the variable
    return reference;
  }

  virtual void SetReference(long d, long op_kind = _HY_OPERATION_VAR) { // return the index of the variable
    Initialize (op_kind, d, _HY_OPERATION_INVALID_REFERENCE, NULL, true);
  }

  virtual bool AssignmentVariable(void) const { return operationKind == _HY_OPERATION_VAR_OBJ; }

  virtual bool HasChanged(bool ignore_cats, const _SimpleList * variable_index = nil) const;

  virtual void SetAttribute(long d) { attribute  = d; }

  virtual _PMathObj GetPayload(void) const { return payload; }

  virtual void SetPayload(_PMathObj d) { payload = d; }
  virtual bool IsDeferred (void) const {
    return operationKind == _HY_OPERATION_DEFERRED_FUNCTION_CALL || 
           operationKind == _HY_OPERATION_DEFERRED_INLINE; }
           
  virtual bool IsComputed (void) const {
    return operationKind == _HY_OPERATION_BUILTIN || 
           operationKind == _HY_OPERATION_FUNCTION_CALL; }
           
  long GetOperationPrecedence (void) const;
    

  long GetAttribute (void) const { return attribute; }
  long GetOpKind    (void) const { return operationKind;}
  long GetReference  (void) const { return reference;}
  bool IsAssociativeOp (void) const;
  bool IsVolatileOp    (void) const;
  bool HasSparseMatrixChanged (void) const;
  void UpdateLValue   (_SimpleList&, _SimpleList&, const long) const;
  long PrepareLHS     (void);
  long OperandCount   (void) const;

  bool CanResultsBeCached(const _Operation *, bool exp_only = false) const;
  void ToggleVarRef    (bool);

  virtual bool EqualOp  (const _Operation *) const;
  _PMathObj    FetchReferencedObject (void) const;
  bool         ToggleFastExec (bool on_off, const _SimpleList* variable_index);
  // returns TRUE if the operation is volatile (e.g Time or Random)

};

#endif
