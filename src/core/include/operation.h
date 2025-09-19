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
    /**
     * @brief Construct a new _Operation object
     */
    _Operation  (void);
    /**
     * @brief Construct a new _Operation object
     *
     * @param s The string representation of the operation
     * @param l The number of operands
     */
    _Operation  (_String const&, const long);
    /**
     * @brief Construct a new _Operation object
     *
     * @param op The opcode
     * @param terms The number of terms
     */
    _Operation  (const long,const long);
    /**
     * @brief Construct a new _Operation object
     *
     * @param op The operation to copy
     */
    _Operation  (const _Operation&);

    /**
     * @brief Construct a new _Operation object
     *
     * @param is_var Whether the operation is a variable
     * @param name The name of the variable
     * @param isG Whether the variable is global
     * @param vc The variable container
     * @param take_a_reference Whether to take a reference to the variable
     */
    _Operation  (bool, _String&, bool isG = false, _VariableContainer const*  = nil, bool take_a_reference = false);
    /**
     * @brief Construct a new _Operation object
     *
     * @param o The object to store
     */
    _Operation  (HBLObjectRef);
    /**
     * @brief Construct a new _Operation object
     *
     * @param v The variable to reference
     */
    _Operation  (_Variable const&);


    /**
     * @brief Destroy the _Operation object
     */
    virtual ~_Operation (void);

    /**
     * @brief Create a dynamic copy of the object
     *
     * @return BaseObj* The dynamic copy
     */
    virtual   BaseObj*      makeDynamic         (void) const;
    /**
     * @brief Duplicate the object
     *
     * @param brc The object to duplicate
     */
    virtual   void          Duplicate           (BaseRefConst);
    /**
     * @brief The assignment operator
     *
     * @param op The operation to assign from
     */
    void      operator = (_Operation const&);

    /**
     * @brief Push the value of the operation onto the stack
     *
     * @param s The stack to push to
     * @param vc The variable container
     * @param errMsg The error message
     * @return true if successful, false otherwise
     */
    bool            PushValue             (_Stack&, _VariableContainer const* = nil, _String* errMsg = nil); 
    /**
     * @brief Execute the operation
     *
     * @param s The stack to use
     * @param vc The variable container
     * @param errMsg The error message
     * @param canCache Whether the result can be cached
     * @return true if successful, false otherwise
     */
    bool            Execute             (_Stack&, _VariableContainer const* = nil, _String* errMsg = nil, bool canCache = false); //execute this operation
    /**
     * @brief Get the stack depth of the operation
     *
     * @param l The stack depth
     */
    virtual   void          StackDepth          (long&);

    /**
     * @brief Execute the operation as a polynomial
     *
     * @param s The stack to use
     * @param nameSpace The namespace
     * @param errMsg The error message
     * @return true if successful, false otherwise
     */
    bool            ExecutePolynomial   (_Stack&,_VariableContainer* nameSpace = nil, _String* errMsg = nil);
    /**
     * @brief Convert the operation to a string
     *
     * @param ul The format to use
     * @return BaseObj* The string representation of the operation
     */
    virtual   BaseObj*      toStr               (unsigned long = 0UL);    //convert the op to string

    /**
     * @brief Initialize the operation
     *
     * @param b Whether to initialize the operation
     */
    virtual   void          Initialize          (bool = true);
    /**
     * @brief Get the code of the operation
     *
     * @return const _String& The code of the operation
     */
    const _String&    GetCode             (void) {
        return (opCode>-1)&&(numberOfTerms>=0)?*(_String*)BuiltInFunctions(opCode): hy_global::kEmptyString;
    }
    /**
     * @brief Get the opcode of the operation
     *
     * @return long& The opcode of the operation
     */
    long&       TheCode             (void) {
        return opCode;
    }
    /**
     * @brief Check if the operation is a variable
     *
     * @param b Whether to check if the operation is a variable
     * @return true if the operation is a variable, false otherwise
     */
    virtual  bool           IsAVariable         (bool = true) ; // is this object a variable or not?
    /**
     * @brief Check if the operation is a constant
     *
     * @param strict Whether to perform a strict check
     * @return true if the operation is a constant, false otherwise
     */
    virtual  bool           IsConstant          (bool strict = false);         // does this object depend on any independent variables or not?
    /**
     * @brief Check if the operation is a HBL function call
     *
     * @return true if the operation is a HBL function call, false otherwise
     */
    virtual  bool           IsHBLFunctionCall   (void) const;
    /**
     * @brief Get the HBL function ID
     *
     * @return long The HBL function ID
     */
    virtual  long           GetHBLFunctionID    (void) const;

    /**
     * @brief Get the user function ID
     *
     * @return long The user function ID
     */
    virtual  long           UserFunctionID      (void) {
        return numberOfTerms < 0 ? opCode : -1;
    };
    // return a non-neg number (function index) if this is a user function,
    // otherwise, return -1

    /**
     * @brief Get the variable index
     *
     * @return long The variable index
     */
    long           GetAVariable        (void) const{
        if (theData >= -1) {
            return theData;
        }
        if (theData < -2) {
            return -theData - 3;
        }
        return -numberOfTerms-1;
        // return the index of the variable
    }

    /**
     * @brief Check if the operation is a value substitution
     *
     * @return true if the operation is a value substitution, false otherwise
     */
    virtual  long           IsValueSubstitution        (void) const{
        return theData == -2;
    }

    /**
     * @brief Set the variable index
     *
     * @param d The variable index
     */
    virtual  void           SetAVariable        (long d) {  // return the index of the variable
        theData=d;
    }

    /**
     * @brief Retrieve the variable
     *
     * @return _Variable* The variable
     */
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

    /**
     * @brief Check if the operation is an assignment variable
     *
     * @return true if the operation is an assignment variable, false otherwise
     */
    virtual  bool           AssignmentVariable  (void) {
        return theData<-2;
    }

    /**
     * @brief Check if the operation has changed
     *
     * @return true if the operation has changed, false otherwise
     */
    virtual  bool           HasChanged          (void);

    /**
     * @brief Set the number of terms
     *
     * @param d The number of terms
     */
    virtual  void           SetTerms            (long d) {
        numberOfTerms=d;
    }

    /**
     * @brief Get the stack depth of the operation
     *
     * @return long The stack depth
     */
    long                    StackDepth          (void) const;

    /**
     * @brief Get the number of the operation
     *
     * @return HBLObjectRef The number
     */
    virtual  HBLObjectRef      GetANumber          (void) {
        return theNumber;
    }

    /**
     * @brief Set the number of the operation
     *
     * @param d The number
     */
    virtual  void           SetNumber           (HBLObjectRef d) {
        theNumber=d;
    }

    /**
     * @brief Get the number of terms
     *
     * @return long The number of terms
     */
    long            GetNoTerms          (void) {
        return numberOfTerms;
    }
    /**
     * @brief Get the precedence level of the operation
     *
     * @return long The precedence level
     */
    long            PrecedenceLevel     (void);

    /**
     * @brief Check if the results can be cached
     *
     * @param op The operation to check
     * @param exp_only Whether to only check for expressions
     * @return true if the results can be cached, false otherwise
     */
    bool            CanResultsBeCached (_Operation *, bool exp_only = false);

    /**
     * @brief Check if the operation is a constant of a given type
     *
     * @param type The type to check
     * @return true if the operation is a constant of the given type, false otherwise
     */
    bool            IsConstantOfType   (const long type) const;


    /**
     * @brief Check if two operations are equal
     *
     * @param op The operation to compare with
     * @return true if the operations are equal, false otherwise
     */
    virtual bool            EqualOp             (_Operation*);

    /**
     * @brief Get the binary opcode of a string
     *
     * @param s The string
     * @param l The long value
     * @return long The binary opcode
     */
    static  long            BinOpCode           (_String const &, long = -1);
    
    /**
     * @brief Check if two operations are inverse
     *
     * @param op_code1 The first opcode
     * @param op_code2 The second opcode
     * @return true if the operations are inverse, false otherwise
     */
    static bool             AreOpsInverse       (long op_code1, long op_code2) {
        long compose = op_code1 < op_code2 ? (op_code1 << 16) + (op_code2) : (op_code2 << 16) + (op_code1);
        return ListOfInverseOps.Find (compose) != kNotFound;
    }
    
    static      _SimpleList         ListOfInverseOps;
protected:


    bool        ReportOperationExecutionError ( _String, _String*);

};

#endif
