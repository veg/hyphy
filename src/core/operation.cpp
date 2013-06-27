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

#include "defines.h"
#include "variable.h"
#include "operation.h"
#include "legacy_parser.h"
#include "polynoml.h"
#include "batchlan.h"
#include "hy_globals.h"

#include "executionlist.h"

extern _SimpleList simpleOperationFunctions, simpleOperationCodes;

//______________________________________________________________________________
_Operation::_Operation(void) {
  Initialize();
}

//______________________________________________________________________________
void _Operation::Initialize(const long kind, const long ref,const long attr,const _PMathObj  load, bool clean) {
  operationKind = kind;
  reference     = ref;
  attribute     = attr;
  if (clean) 
    DeleteObject(payload);
  payload       = load;
}

//______________________________________________________________________________
BaseRef _Operation::makeDynamic(void) {
  _Operation *res = new _Operation;
  res->Duplicate(this);
  return res;
}

//______________________________________________________________________________
_Operation::_Operation(_Operation& copyFrom) {
    Duplicate (&copyFrom);
}

//______________________________________________________________________________
void _Operation::Duplicate(BaseRef r) {
  _Operation *o = (_Operation *)r;
  operationKind = o->operationKind;
  reference = o->reference;
  attribute = o->attribute;
  payload   = o->payload;
  if (payload) {
    payload->nInstances++;
  }
}

//______________________________________________________________________________
BaseRef _Operation::toStr(void) {
  _String * res = new _String (128L, true);
  
  (*res) << "<_Operation: ";
  switch (operationKind) {
    case _HY_OPERATION_NOOP:
      (*res) << "no-op";
      break;
    case _HY_OPERATION_VALUE:
      (*res) << "push value ";
      res->AppendNewInstance((_String*)payload->toStr());
      break;
    case _HY_OPERATION_VAR:
      if (attribute == _HY_OPERATION_TOGGLE) {
        (*res) << "push the index of variable ";      
      } else {
        (*res) << "push computed value of variable ";
      }
      (*res) << *LocateVar(reference)->GetName();
      break;
    case _HY_OPERATION_VAR_OBJ:
      if (attribute == _HY_OPERATION_TOGGLE) {
        (*res) << "push the index of variable ";      
      } else {
        (*res) << "push object value of variable ";
      }
      (*res) << *LocateVar(reference)->GetName();
      break;
    case _HY_OPERATION_DEFERRED_INLINE:
      (*res) << "deferred substitute value of variable ";
      (*res) << *LocateVar(reference)->GetName();
      break;
    case _HY_OPERATION_REF:
      (*res) << "push reference in variable ";
      (*res) << *LocateVar(reference)->GetName();
      break;
    case _HY_OPERATION_BUILTIN:
      (*res) << "execute a built-in operation ";
      (*res) << GetCode();
      (*res) << " (consume ";
      (*res) << _String(attribute);
      (*res) << " values off the stack)";
      break;
    case _HY_OPERATION_FUNCTION_CALL:
      (*res) << "call an HBL function ";
      (*res) << *_HBLObjectNameByType(HY_BL_HBL_FUNCTION, reference);
      break;
    case _HY_OPERATION_DEFERRED_FUNCTION_CALL:
      (*res) << "call an HBL function (deferred)";
      (*res) << (_String*)payload;
      break;
    case _HY_OPERATION_DICTIONARY:
      (*res) << "Create a dictionary";
       break;
    case _HY_OPERATION_SPARSE_MATRIX:
      (*res) << "Create a sparse matrix";
       break;
    case _HY_OPERATION_ASSIGNMENT_VALUE:
      (*res) << "LHS = RHS";
       break;
    case _HY_OPERATION_ASSIGNMENT_EXPRESSION:
      (*res) << "LHS := RHS";
      break;
  }
  (*res) << ">";
  res->Finalize();
  return res;
  
}

//______________________________________________________________________________
_String _Operation::GetCode(void) {
    if (operationKind == _HY_OPERATION_BUILTIN) {
      return BuiltInFunctions.RetrieveKeyByPayload (reference);
    }
    return empty;
}


//______________________________________________________________________________
_Operation::_Operation(_PMathObj theObj) {
  Initialize(_HY_OPERATION_VALUE,_HY_OPERATION_INVALID_REFERENCE,_HY_OPERATION_INVALID_REFERENCE,theObj);
}

//______________________________________________________________________________
_Operation::~_Operation(void) {
  if ((operationKind & _HY_OPERATION_FAST_EXEC) == 0L) {
    if (operationKind == _HY_OPERATION_SPARSE_MATRIX ||
       operationKind == _HY_OPERATION_DICTIONARY) {
      ((_SimpleList*)payload)->ClearFormulasInList();
    } else {
      if (operationKind == _HY_OPERATION_ASSIGNMENT_EXPRESSION) {
        delete ((_Formula*) payload);
        payload = nil;
      }
    }
    DeleteObject(payload);
  }
}

//______________________________________________________________________________
_Operation::_Operation(const long opk, const long ref, const long attr, _PMathObj data) {
    Initialize(opk, ref, attr, data);
}


//______________________________________________________________________________
  // a built-in
_Operation::_Operation(const long theCode, const long opNo) {
  Initialize (_HY_OPERATION_BUILTIN,theCode,opNo);
}

//______________________________________________________________________________
  // construct the operation by built-in

_Operation::_Operation(bool isUserFunction, _String &opc, const long opNo = 2) {
  
  long opCode = isUserFunction?FindBFFunctionName(opc):BuiltInFunctions.GetValueFromString(opc);
  
  if (opCode != HY_NOT_FOUND) {
    Initialize (isUserFunction?_HY_OPERATION_FUNCTION_CALL:_HY_OPERATION_BUILTIN,opCode,opNo);
  } else {
    if (isUserFunction) {
      Initialize (_HY_OPERATION_DEFERRED_FUNCTION_CALL,_HY_OPERATION_INVALID_REFERENCE,opNo,(_PMathObj)opc.makeDynamic());
    } else {
      WarnError(_String("'") & opc & "' is not a defined built-in function or operation.");
      Initialize ();
    }
  }
}


//______________________________________________________________________________
_Operation::_Operation(_String &stuff, bool defineAVar, bool global_flag, 
  _VariableContainer *theParent, bool take_a_reference, bool deferred_inline) {
  
    // creating a variable
  if (defineAVar) {
    long f;
    _String variable_name (stuff);
    if (theParent /*&&(!isG)*/) {
        // 20070620: SLKP the commenting may break
        // default behavior!
      
      f = LocateVarByName(variable_name);
      
      if (f != HY_NOT_FOUND && !FetchVar(f)->IsGlobal()) {
        f = HY_NOT_FOUND;
      }
      
      if (f == HY_NOT_FOUND) {
        variable_name = (*theParent->theName) & '.' & variable_name;
      }
    }
    
    _Variable * v = CheckReceptacle (&variable_name, empty, false, global_flag);
    Initialize (deferred_inline ? _HY_OPERATION_DEFERRED_INLINE : 
      (take_a_reference ? _HY_OPERATION_REF : _HY_OPERATION_VAR), v->GetAVariable());
    
  } else {
    Initialize(_HY_OPERATION_VALUE, _HY_OPERATION_INVALID_REFERENCE, _HY_OPERATION_INVALID_REFERENCE, 
               stuff.Equal(&noneToken) ? new _MathObject : new _Constant(stuff));
  }
  
}

//______________________________________________________________________________

bool _Operation::ExecuteBuiltIn (_Stack& theScrap, _hyExecutionContext* context) {
  
  if (theScrap.theStack.lLength < (unsigned long)attribute) {
    return ReportOperationExecutionError(
                                         _String((_String *)toStr()) & " needs " & _String(attribute) &
                                         " arguments; " & _String(theScrap.StackDepth()) & " were supplied",
                                         context->GetErrorBuffer());
    
  }
  
  _PMathObj terms[3] = {NULL, NULL, NULL},
  op_result;
  
  long sL = theScrap.theStack.lLength;
  for (long k = attribute-1; k>=0; k--) {
    terms[k] = (_PMathObj) theScrap.theStack.lData[--sL];
  }
  theScrap.theStack.lLength = sL;
  op_result = terms[0]->Execute(reference, terms[1], terms[2], context);
  for (long k = 0; k < attribute; k++) {
    DeleteObject (terms[k]);
  }
  
  if (op_result) {
    theScrap.theStack.Place(op_result);
    return true;
  } 
  
  return false;  
  
}

//______________________________________________________________________________

long  _Operation::PrepareLHS(void) {
  switch (operationKind) {
    case _HY_OPERATION_VAR:
    case _HY_OPERATION_VAR_OBJ:
      attribute = _HY_OPERATION_TOGGLE;
      break;
    case _HY_OPERATION_BUILTIN:
      if (reference == HY_OP_CODE_MACCESS) {
        return attribute;
      }
      break;
  }
  return _HY_OPERATION_INVALID_REFERENCE;
}

//______________________________________________________________________________

bool _Operation::ExecuteAssignment (_Stack& theScrap, _hyExecutionContext* context) {
  
  if (operationKind == _HY_OPERATION_ASSIGNMENT_VALUE
      || _HY_OPERATION_ASSIGNMENT_EXPRESSION) {
      
      _PMathObj rhs = operationKind == _HY_OPERATION_ASSIGNMENT_VALUE ? theScrap.Pop() : NULL;
      if (reference == _HY_OPERATION_INVALID_REFERENCE) {
        // this handles x ?= y assignments
        _PMathObj lhs = theScrap.Pop();
                  
        _Variable *lhs_var = LocateVar((long)lhs->Compute()->Value()); 
        
        if (operationKind == _HY_OPERATION_ASSIGNMENT_VALUE) {
          if (attribute != HY_OP_CODE_NONE) {
             _PMathObj op_result = lhs_var->Compute()->Execute(attribute, rhs);
             if (!op_result) {
                return false;
             }
             DeleteObject (rhs);
             rhs = op_result;
          }
          
          lhs_var -> SetValue (rhs, false);
          rhs->AddAReference();
          theScrap.theStack.Place (rhs);
        } else {
          lhs_var->SetFormula(*(_Formula*)payload);
          theScrap.theStack.Place (new _MathObject);
         }
        
       DeleteObject (lhs);
       return true;
      } else {
        // this handles x[][] ?= y assignments
        
        if (reference > 3) {
          return false;
        }
        _PMathObj lhs_accessors [3] = {nil, nil, nil};
        
                  
        for (long k = 0; k < reference; k++) {
          lhs_accessors[reference-k-1] = theScrap.Pop();
        }
       
        _PMathObj res = nil; 
        
        
        if (attribute != HY_OP_CODE_NONE) {
          _PMathObj  base_value = lhs_accessors[0] -> Execute(HY_OP_CODE_MACCESS, lhs_accessors[1], lhs_accessors[2], context),
          op_result = nil;
          
          if (base_value) {
            op_result = base_value->Execute(attribute, rhs);
            if (op_result) {
              DeleteObject (rhs);
            }
            DeleteObject(base_value);
          }
          
          rhs = op_result;
          
        }
                          
        if (rhs) {
          res = lhs_accessors[0] -> Execute(HY_OP_CODE_MSTORE, lhs_accessors[1], lhs_accessors[2], context, rhs);
        } else {
          res = lhs_accessors[0] -> Execute(HY_OP_CODE_FSTORE, lhs_accessors[1], lhs_accessors[2], context, (_PMathObj)payload);
        }
        
        for (long k = 0; k < reference; k++) {
          DeleteObject (lhs_accessors[reference-k-1]);
        }
        DeleteObject(rhs);
        
        if (res) {
          theScrap.theStack.Place (res);        
          return true;
        }
        return false;
          
      }
  }
  return false;
}

//______________________________________________________________________________

bool _Operation::ExecuteFunctionCall (_Stack& theScrap, _hyExecutionContext* context) {
  
  
  _List * func_parameters = _HYFetchFunctionParameters (reference);
  
  if (!func_parameters) {
    return ReportOperationExecutionError(
                                         "Attempted to call an undefined HBL function ", context->GetErrorBuffer());
  }
  
  long expected_arguments = func_parameters->lLength;
  
  if (theScrap.StackDepth() < expected_arguments || attribute != expected_arguments) {
    return ReportOperationExecutionError(
                                         _String("User-defined function '") &
                                         &*_HBLObjectNameByType(HY_BL_HBL_FUNCTION, reference) &
                                         "' needs " & _String(expected_arguments) & " parameters: " &
                                         _String(theScrap.StackDepth()) & " were supplied ",
                                         context->GetErrorBuffer());
  }
  
  _List displacedVars,
  displacedValues, 
  referenceArgs;
  
  _SimpleList existingIVars, 
  existingDVars, 
  displacedReferences;
  
  _String     *argNameString;
  
  bool purgeFlas = false;
  
  for (long k = attribute - 1; k >= 0; k--) {
    bool isRefVar = false;
    argNameString = (_String *)(*func_parameters)(k);
    
    if (argNameString->sData[argNameString->sLength - 1] == '&') {
      argNameString->Trim(0, argNameString->sLength - 2);
      isRefVar  = true;
      purgeFlas = true;
    }
    
    
    _PMathObj nthterm = theScrap.Pop();
    
    if (isRefVar && nthterm->ObjectClass() != STRING) {
      _FString *type = (_FString *)nthterm->Type();
      _String errText =
      _String("User-defined function '") &
      *_HBLObjectNameByType(HY_BL_HBL_FUNCTION, reference) &
      "' expected a string for the reference variable '" &
      *argNameString & "' but was passed a " & *type->theString &
      " with the value of " & _String((_String *)nthterm->toStr());
      
      DeleteObject(type);
      return ReportOperationExecutionError(errText, context->GetErrorBuffer());
    }
    
    _Variable *theV = CheckReceptacle (argNameString, empty, false, false);
    
    if (!isRefVar) {
        // if the variable exists and is independent then
        // simply swap the value of the var, otherwise
        // duplicate the entire variable
      if (theV->IsIndependent()) {
        theV->varFlags &= HY_VARIABLE_SET;
        if (!theV->varValue) {
          theV->Compute();
        }
        displacedValues << theV->varValue;
        theV->varFlags |= HY_VARIABLE_CHANGED;                
        theV->varValue = nthterm;
        existingIVars << theV->GetAVariable();
      } else {
        _Variable newV(*argNameString);
        newV.SetValue(nthterm);
        existingDVars << theV->GetAVariable();
        displacedVars << theV;
        variablePtrs.Replace(theV->GetAVariable(),
                             (_PMathObj) newV.makeDynamic());
        DeleteObject(nthterm);
      }
    } else {
      referenceArgs << theV->GetName();
      displacedReferences << theV->GetAVariable();
      
      _String *refArgName = ((_FString *)nthterm)->theString;
      
      if (context->GetContext()) {
        *refArgName = AppendContainerName(*refArgName, context->GetContext());
      }
      
      _Variable * newRefArgument = CheckReceptacle (refArgName, empty, false, false);
      
      variableNames.SetXtra(theV->GetAVariable(), variableNames.GetXtra(newRefArgument->GetAVariable()));
      *argNameString = *argNameString & '&';
      DeleteObject(nthterm);
    }
  }
  
  _ExecutionList *functionBody = _HYFetchFunctionBody(reference);
  
  if (purgeFlas) {
    functionBody->ResetFormulae();
  }
  
  
  if (currentExecutionList && currentExecutionList->stdinRedirect) {
    functionBody->stdinRedirect = currentExecutionList->stdinRedirect;
    functionBody->stdinRedirectAux = currentExecutionList->stdinRedirectAux;
  }
  
  _PMathObj ret = functionBody->Execute();
  
  functionBody->stdinRedirect    = nil;
  functionBody->stdinRedirectAux = nil;
  
  if (terminateExecution) {
    theScrap.Push(new _Constant(0.0));
    return true;
  }
  
  if (ret) {
    theScrap.Push(ret);
  }
  
  for (unsigned long di = 0; di < referenceArgs.lLength; di++) {
    variableNames.SetXtra(LocateVarByName(*(_String *)referenceArgs(di)),
                          displacedReferences.lData[di]);
  }
  
  
  for (unsigned long dv = 0; dv < displacedVars.lLength; dv++) {
    variablePtrs.Replace(existingDVars.lData[dv],
                         (_PMathObj) displacedVars(dv));
  }
  
  
  for (unsigned long dv2 = 0; dv2 < displacedValues.lLength; dv2++) {
    _Variable *theV = LocateVar(existingIVars.lData[dv2]);
    DeleteObject(theV->varValue);
    theV->varValue = ((_PMathObj) displacedValues(dv2));
  }
  
  
  return true;
}  


//______________________________________________________________________________

void _Operation::UpdateLValue (_SimpleList & stack, _SimpleList & indices, const long index) const {
  switch (operationKind) {
    case _HY_OPERATION_VAR:
    case _HY_OPERATION_REF:
    case _HY_OPERATION_VAR_OBJ:
      stack   << reference;
      indices << index;
      break;
    case _HY_OPERATION_BUILTIN:
    case _HY_OPERATION_FUNCTION_CALL:
    case _HY_OPERATION_DEFERRED_FUNCTION_CALL: {
        bool is_index = reference == HY_OP_CODE_MACCESS && operationKind == _HY_OPERATION_BUILTIN;
        long up_to = attribute - is_index;
        for (long k = 0; k < up_to; k++) {
          stack.Pop();
          indices.Pop();
        }
        if (! is_index) {
          stack << HY_NOT_FOUND;
          indices << index;
        } else {
          stack << stack.GetElement(-1L);
          indices << index;
        }
    }
    break;
    default:
      stack << HY_NOT_FOUND;
      indices << index;
      break;
  }
}
//______________________________________________________________________________
bool _Operation::Execute(_Stack &theScrap, _hyExecutionContext* context) {
  
  switch (operationKind) {
    case _HY_OPERATION_VALUE:
      theScrap.Push(payload);
      return true;
    case _HY_OPERATION_VAR:
      if (attribute == _HY_OPERATION_INVALID_REFERENCE) {
        theScrap.Push(((_Variable *)((BaseRef *)variablePtrs.lData)[reference])->Compute());
      } else {
        theScrap.Push(new _Constant (reference), false);
      }
      return true;
    case _HY_OPERATION_REF:
      theScrap.Push(((_Variable *)((BaseRef *)variablePtrs.lData)[reference])
                    ->ComputeReference(nil), false);
      return true;
    case _HY_OPERATION_VAR_OBJ:
      if (attribute == _HY_OPERATION_INVALID_REFERENCE) {
        theScrap.Push(((_Variable *)((BaseRef *)variablePtrs.lData)[reference])
                      ->GetValue());
      } else {
        theScrap.Push(new _Constant (reference), false);      
      }
      
      return true;

    case _HY_OPERATION_SPARSE_MATRIX:
      theScrap.Push(new _Matrix (payload, false), false);
      return true;

    case _HY_OPERATION_DICTIONARY:
      theScrap.Push(new _AssociativeList (payload), false);
      return true;

    case _HY_OPERATION_BUILTIN:  
      return ExecuteBuiltIn (theScrap, context);
    case _HY_OPERATION_FUNCTION_CALL:
      return ExecuteFunctionCall (theScrap, context);
    case _HY_OPERATION_DEFERRED_FUNCTION_CALL:
    case _HY_OPERATION_DEFERRED_INLINE:
      ResolveDeferredAction (context);
      return Execute (theScrap, context);
      
    case _HY_OPERATION_ASSIGNMENT_VALUE: 
    case _HY_OPERATION_ASSIGNMENT_EXPRESSION: 
      return ExecuteAssignment (theScrap, context);
    
      
    case _HY_OPERATION_NOOP:
      return true;
      
      
  }
  return false;
  
  
}

//______________________________________________________________________________
bool _Operation::ExecuteFast (_SimpleFormulaDatum *stack,
                    _SimpleFormulaDatum *varValues,
                    long& stackTop,
                    _String *errMsg) {
  
  
  switch (operationKind) {
    case _HY_OPERATION_FAST_EXEC_VALUE:
      stack[stackTop++].value = payload->Value();
      return true;

    case _HY_OPERATION_FAST_EXEC_VAR:
    case _HY_OPERATION_FAST_EXEC_VAR_OBJ:
      stack[stackTop++] = varValues [attribute];
      return true;
      
    case _HY_OPERATION_FAST_EXEC_BUILTIN:
    case _HY_OPERATION_FAST_EXEC_BUILTIN_REF: {
        if (attribute==2) {
          stackTop--;
          if (stackTop < 1) {
            ReportOperationExecutionError("Internal error in _Formula::ComputeSimple - stack underflow.)", errMsg);
            return false;
          }
          if (operationKind == _HY_OPERATION_FAST_EXEC_BUILTIN) {
            _Parameter(*theFunc)(_Parameter, _Parameter);
            theFunc = (_Parameter(*)(_Parameter, _Parameter)) payload;
            stack[stackTop - 1].value =
              (*theFunc)(stack[stackTop - 1].value, stack[stackTop].value);
          } else {
            _Parameter(*theFunc)(Ptr, _Parameter);
            theFunc = (_Parameter(*)(Ptr, _Parameter)) payload;
            stack[stackTop - 1].value = (*theFunc)(
                stack[stackTop - 1].reference, stack[stackTop].value);
          }
        } else {
          if (attribute == 1) {
            if (operationKind == _HY_OPERATION_FAST_EXEC_BUILTIN) {
              _Parameter(*theFunc)(_Parameter);
              theFunc = (_Parameter(*)(_Parameter)) payload;
              stack[stackTop++].value = (*theFunc)(stack[stackTop].value);
            } else {
              _Parameter(*theFunc)(Ptr);
              theFunc = (_Parameter(*)(Ptr)) payload;
              stack[stackTop++].value = (*theFunc)(stack[stackTop].reference);
            }
          } else {
            return false;
          }
        }
      }
      return true;
      

    case _HY_OPERATION_NOOP:
      return true;

}
  
  return false;
  
  
}

//______________________________________________________________________________
bool _Operation::HasChanged(bool ignore_cats, const _SimpleList * variable_index) const{
  switch (operationKind) {
    case _HY_OPERATION_VALUE:
      return payload->HasChanged();
      
    case _HY_OPERATION_FAST_EXEC_VALUE:
        return false;
        
    case _HY_OPERATION_VAR:
    case _HY_OPERATION_REF:
    case _HY_OPERATION_VAR_OBJ:
    case _HY_OPERATION_DEFERRED_INLINE:
      return LocateVar(reference)->HasChanged(ignore_cats);
      
    case _HY_OPERATION_FAST_EXEC_VAR:
    case _HY_OPERATION_FAST_EXEC_VAR_OBJ:
      return LocateVar(variable_index->GetElement (attribute)) -> HasChanged (false);
      
    case _HY_OPERATION_BUILTIN:
    case _HY_OPERATION_FAST_EXEC_BUILTIN:
    case _HY_OPERATION_FAST_EXEC_BUILTIN_REF:
      return IsVolatileOp();
        
    case _HY_OPERATION_SPARSE_MATRIX:
    case _HY_OPERATION_DICTIONARY:
      return HasSparseMatrixChanged ();

    case _HY_OPERATION_DEFERRED_FUNCTION_CALL:
      return true;
      
    case _HY_OPERATION_FUNCTION_CALL:
      return _HYGetHBLCallType (reference) != BL_FUNCTION_NORMAL_UPDATE;
  }
  
  return false;
}

//______________________________________________________________________________
bool _Operation::IsAVariable(bool deep) const{    
  switch (operationKind) {
    case _HY_OPERATION_VALUE:
      if (deep) {
        return payload->IsVariable();
      }
      return false;            
      
    case _HY_OPERATION_SPARSE_MATRIX:
    case _HY_OPERATION_DICTIONARY:
    case _HY_OPERATION_VAR:
    case _HY_OPERATION_REF:
    case _HY_OPERATION_VAR_OBJ:
    case _HY_OPERATION_DEFERRED_INLINE:
      return true;
      
  }
  
  return false;
}

//______________________________________________________________________________
_PMathObj    _Operation::FetchReferencedObject (void) const {
  switch (operationKind) {
    case _HY_OPERATION_VALUE:
      return payload;
    case _HY_OPERATION_VAR:
    case _HY_OPERATION_VAR_OBJ: {
      _Variable * r = LocateVar(reference);
      if (r) {
        return r->GetValue();
      }
    }
  
  }
  return nil;
}

//______________________________________________________________________________
bool _Operation::IsAFunctionCall(void) const{
  return operationKind == _HY_OPERATION_FUNCTION_CALL || operationKind == _HY_OPERATION_DEFERRED_FUNCTION_CALL ;
}

//______________________________________________________________________________
bool _Operation::IsConstant(void) const {
  
  switch (operationKind) {
    case _HY_OPERATION_VALUE:
      return payload->IsConstant();
      
    case _HY_OPERATION_VAR:
    case _HY_OPERATION_REF:
    case _HY_OPERATION_VAR_OBJ:
    case _HY_OPERATION_DEFERRED_INLINE:
      return LocateVar(reference)->IsConstant();
      
    case _HY_OPERATION_BUILTIN:
      return !IsVolatileOp();
      
    case _HY_OPERATION_NOOP:
      return true;
      
    case _HY_OPERATION_ASSIGNMENT_VALUE:
      return true;
      
  }
  return false;
}

//______________________________________________________________________________
bool _Operation::ReportOperationExecutionError(_String text, _String *errMsg) const {
  
  _String theError = text & ". ";
  
  if (errMsg) {
    *errMsg = theError;
  } else {
    WarnError(theError);
  }
  
  return false;
}

//______________________________________________________________________________
void _Operation::StackDepth(long &depth) const {
  
  switch (operationKind) {
    case _HY_OPERATION_VALUE:
    case _HY_OPERATION_VAR:
    case _HY_OPERATION_REF:
    case _HY_OPERATION_VAR_OBJ:
    case _HY_OPERATION_DEFERRED_INLINE:
    case _HY_OPERATION_SPARSE_MATRIX:
    case _HY_OPERATION_DICTIONARY:
    case _HY_OPERATION_ASSIGNMENT_VALUE:
      depth++;
      break;
      
    case _HY_OPERATION_BUILTIN:
    case _HY_OPERATION_FUNCTION_CALL:
    case _HY_OPERATION_DEFERRED_FUNCTION_CALL:
      depth += attribute-1;  
      
    case _HY_OPERATION_ASSIGNMENT_EXPRESSION:
    case _HY_OPERATION_ASSIGNMENT_BOUND:
      depth --;
      
      
  }
  
}


//______________________________________________________________________________
bool _Operation::CanResultsBeCached(const _Operation *prev, bool exp_only) const {

  if (operationKind == _HY_OPERATION_BUILTIN && attribute == 1) {
     if ((prev->operationKind == _HY_OPERATION_VALUE && prev->GetPayload()->ObjectClass() == MATRIX) ||
        ((prev->operationKind == _HY_OPERATION_VAR_OBJ || prev->operationKind == _HY_OPERATION_VAR) &&
         LocateVar(prev->GetAVariable())->ObjectClass() == MATRIX)) {
          if (!exp_only || reference == HY_OP_CODE_EXP)
            return true;
        }
  }
  return false;
}

//______________________________________________________________________________

void _Operation::ResolveDeferredAction (_hyExecutionContext* context ){
    if (operationKind == _HY_OPERATION_DEFERRED_FUNCTION_CALL) {
        long hbl_type = HY_BL_HBL_FUNCTION;
        _PMathObj func = (_PMathObj)_HYRetrieveBLObjectByName(*(_String*)payload,hbl_type , &reference);
        if (func) {
          Initialize (_HY_OPERATION_FUNCTION_CALL,  reference, attribute, NULL, true);
        } 
        else {
          ReportOperationExecutionError(
                                         _String("Failed to bind an HBL function '") &
                                         &*(_String*)payload &
                                         "' at run time ",
                                         context->GetErrorBuffer());
          Initialize(_HY_OPERATION_NOOP, _HY_OPERATION_INVALID_REFERENCE, _HY_OPERATION_INVALID_REFERENCE,
                NULL, true);
        
        }
    } 
    else {
      if (operationKind == _HY_OPERATION_DEFERRED_INLINE) {
        Initialize (_HY_OPERATION_VALUE, _HY_OPERATION_INVALID_REFERENCE, _HY_OPERATION_INVALID_REFERENCE,
                    (_PMathObj) LocateVar(reference)->Compute()->makeDynamic());
      }
    }
}


//______________________________________________________________________________

long _Operation::GetOperationPrecedence (void) const {
  switch (operationKind) {
    case _HY_OPERATION_FUNCTION_CALL:
    case _HY_OPERATION_DEFERRED_FUNCTION_CALL:
      return _HY_OPERATION_MAX_PRECEDENCE;
      
    case _HY_OPERATION_BUILTIN: {
      if (attribute == 2) {
        switch (reference) {
            case HY_OP_CODE_OR:
              return 1L;
            
            case HY_OP_CODE_AND:
              return 2L;
            
            case HY_OP_CODE_EQ:
            case HY_OP_CODE_NEQ:
            case HY_OP_CODE_LESS:
            case HY_OP_CODE_GREATER:
            case HY_OP_CODE_LEQ:
            case HY_OP_CODE_GEQ:
              return 3L;
            
            case HY_OP_CODE_ADD:
            case HY_OP_CODE_SUB:
              return 4L;

            case HY_OP_CODE_MUL:
            case HY_OP_CODE_DIV:
            case HY_OP_CODE_IDIV:
            case HY_OP_CODE_MOD:
              return 5L;
              
            case HY_OP_CODE_POWER:
              return 6L;

            default:
              return _HY_OPERATION_MAX_PRECEDENCE;
              
        }
      } else {
        return _HY_OPERATION_MAX_PRECEDENCE;
      }
    }
    
  }
  
  return _HY_OPERATION_MIN_PRECEDENCE;
}

//______________________________________________________________________________

bool _Operation::IsAssociativeOp (void) const {
  
  if (operationKind & _HY_OPERATION_OP_CLASS) {
    return (reference == HY_OP_CODE_ADD || reference == HY_OP_CODE_MUL) && attribute == 2;
  }
  
  return false;
}

//______________________________________________________________________________

bool _Operation::IsVolatileOp (void) const {
  
  if (operationKind & _HY_OPERATION_OP_CLASS) {
    if (reference == HY_OP_CODE_RANDOM || reference == HY_OP_CODE_TIME) {
      return true;
    }
  }
  
  return false;
}

//______________________________________________________________________________
bool _Operation::EqualOp(const _Operation *otherOp) const{
  
  if (otherOp->operationKind == operationKind) {
    
    switch (operationKind) {
      case _HY_OPERATION_VALUE: {
        
        unsigned long oc = payload->ObjectClass();
        
        if (oc == NUMBER && oc == otherOp->payload->ObjectClass()) {
          return CheckEqual(payload->Value(), otherOp->payload->Value());
        }      
        
        return false;
      }
        
      case _HY_OPERATION_BUILTIN:
      case _HY_OPERATION_VAR:
      case _HY_OPERATION_REF:
      case _HY_OPERATION_VAR_OBJ:
      case _HY_OPERATION_FUNCTION_CALL:
      case _HY_OPERATION_DEFERRED_INLINE:
        return reference == otherOp->reference && attribute == otherOp->attribute;
        
      case _HY_OPERATION_DEFERRED_FUNCTION_CALL:
        return attribute == otherOp->attribute && ((_String*)payload)->Equal((_String*)otherOp->payload);
        
      case _HY_OPERATION_NOOP:
        return true;
    }
        
  }
    
  return false;
}


//______________________________________________________________________________
bool _Operation::ExecutePolynomial(_Stack &theScrap,
                                   _VariableContainer *nameSpace,
                                   _String *errMsg) {
  
  switch (operationKind) {
    case _HY_OPERATION_VALUE:
      theScrap.Push(new _Polynomial(payload->Value()), false);
      return true;
    case _HY_OPERATION_VAR:
    case _HY_OPERATION_VAR_OBJ:
      theScrap.Push(new _Polynomial(*LocateVar(reference)), false);
      return true;
    case _HY_OPERATION_BUILTIN:  {
        if (theScrap.StackDepth() < attribute) {
            _String s((_String *)toStr());
            return ReportOperationExecutionError(
                                                 s & " needs " & _String(attribute) & " arguments. Only " &
                                                 _String(theScrap.StackDepth()) & " were given",
                                                 errMsg);
          }
          
          _PMathObj term1, term2 = nil, temp;
          
          bool opResult = true;
          
          if (attribute == 2) {
            term2 = theScrap.Pop();
          }
          
          _hyExecutionContext localContext(nameSpace, errMsg);
          term1 = theScrap.Pop();
          temp = term1->Execute(reference, term2, nil, &localContext);
          DeleteObject(term1);
          
          if (temp) {
            theScrap.Push(temp, false);
          } else {
            opResult = false;
          }
          
          if (term2) {
            DeleteObject(term2);
          }
          
          return opResult;
      
      
    }
    case _HY_OPERATION_NOOP:
      return true;

      
  }
  return false;

}

//______________________________________________________________________________
bool _Operation::HasSparseMatrixChanged(void) const {
  if (operationKind == _HY_OPERATION_SPARSE_MATRIX
     || operationKind == _HY_OPERATION_DICTIONARY) {
    _SimpleList * listSpec = (_SimpleList *) payload;
    for (unsigned long k = 0; k < listSpec->lLength; k++) {
       if (((_Formula*)listSpec->Element(k))->HasChanged()) {
        return true;
       }
    } 
  }
  return false;
}
//______________________________________________________________________________
void _Operation::ToggleVarRef(bool on_off) {
  if (on_off) {
    if (operationKind == _HY_OPERATION_VAR) {
      operationKind = _HY_OPERATION_VAR_OBJ;
    }
  } else {
    if (operationKind == _HY_OPERATION_VAR_OBJ) {
      operationKind = _HY_OPERATION_VAR;
    }
  }
}

//______________________________________________________________________________
bool _Operation::ToggleFastExec(bool on_off, const _SimpleList* variable_index) {

  if (on_off) {
    switch (operationKind) {
        case _HY_OPERATION_VALUE:
          operationKind = _HY_OPERATION_FAST_EXEC_VALUE;
          break;
        case _HY_OPERATION_VAR:
          operationKind = _HY_OPERATION_FAST_EXEC_VAR;
          attribute = variable_index->Find (reference);
          break;
        case _HY_OPERATION_VAR_OBJ:
          operationKind = _HY_OPERATION_FAST_EXEC_VAR_OBJ;
          attribute = variable_index->Find (reference);
          break;
        case _HY_OPERATION_BUILTIN: {
            operationKind = (reference == HY_OP_CODE_MACCESS) ? _HY_OPERATION_FAST_EXEC_BUILTIN_REF :  _HY_OPERATION_FAST_EXEC_BUILTIN; 
            if (reference == HY_OP_CODE_SUB && attribute == 1L) {
                payload = (_PMathObj)MinusNumber;
            } else {
                payload = (_PMathObj)simpleOperationFunctions(simpleOperationCodes.Find(reference)); 
            }
          }
          break;
        case _HY_OPERATION_FAST_EXEC_VALUE:
        case _HY_OPERATION_FAST_EXEC_BUILTIN_REF:
        case _HY_OPERATION_FAST_EXEC_BUILTIN:
        case _HY_OPERATION_FAST_EXEC_VAR:
        case _HY_OPERATION_FAST_EXEC_VAR_OBJ:
          break;
          
        default:
          WarnError ("Internal error in _Operation::ToggleFastExec: invalid operation type");
          break;
        
        
    }
    return IsVolatileOp();
  }


    switch (operationKind) {
        case _HY_OPERATION_FAST_EXEC_VALUE:
          operationKind = _HY_OPERATION_VALUE;
          break;
        case _HY_OPERATION_FAST_EXEC_VAR:
          operationKind = _HY_OPERATION_VAR;
          attribute = _HY_OPERATION_INVALID_REFERENCE;
          break;
        case _HY_OPERATION_FAST_EXEC_VAR_OBJ:
          operationKind = _HY_OPERATION_VAR_OBJ;
          attribute = _HY_OPERATION_INVALID_REFERENCE;
          break;
        case _HY_OPERATION_FAST_EXEC_BUILTIN_REF: 
        case _HY_OPERATION_FAST_EXEC_BUILTIN: {
            operationKind = _HY_OPERATION_BUILTIN; 
            payload = NULL;
          }
          break;        
    }

  return false;
}

