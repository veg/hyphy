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

#include "parser.h"
#include "polynoml.h"
#include "batchlan.h"
#include "parser.h"

extern _SimpleList BinOps,
       opPrecedence,
       FunctionArgumentCount,
       associativeOps;

//__________________________________________________________________________________

_Operation::_Operation  (void)
{
    numberOfTerms = 0;
    theData = -1;
    theNumber = nil;
}

//__________________________________________________________________________________
void    _Operation::Initialize(bool)
{
    numberOfTerms = 0;
    theData = -1;
    theNumber = nil;
}

//__________________________________________________________________________________
BaseRef _Operation::makeDynamic (void)
{
    _Operation * res = new _Operation;
    checkPointer(res);
    //memcpy ((char*)res, (char*)this, sizeof (_Operation));
    res->Duplicate(this);
    return res;
}

//__________________________________________________________________________________
void    _Operation::Duplicate(BaseRef r)
{
    _Operation * o = (_Operation*)r;
    numberOfTerms  = o->numberOfTerms;
    theData        = o->theData;
    theNumber      = o->theNumber;
    opCode         = o->opCode;
    if (theNumber) {
        theNumber->nInstances++;
    }
}


//__________________________________________________________________________________
BaseRef _Operation::toStr(void)
{
    _String * res = new _String;
  
    if (theData != -1) {
        *res = _String("Variable ")& *LocateVar(GetAVariable())->GetName();
    } else if (theNumber) {
        _FString * type = (_FString*)theNumber->Type();
       *res = _String("Constant (")& *type->theString & ")" & _String((_String*)theNumber->toStr());
        DeleteObject(type);
    } else {
        if (IsAFunctionCall())
          *res = GetBFFunctionNameByIndex(UserFunctionID());
        else
          *res = _String("Operation ")&*(_String*)BuiltInFunctions(opCode);
    }

  
    return res;

}

//__________________________________________________________________________________
_Operation::_Operation  (const long theCode, const long opNo = 2)
// by opcode
{
    opCode = theCode;
    numberOfTerms = opNo;
    theData       = -1;
    theNumber     = nil;
}

//__________________________________________________________________________________
_Operation::_Operation  (_String& opc, const long opNo = 2)
// construct the operation by its symbol and, if relevant -
// number of operands
{
    if(opNo>=0) {
        opCode = BuiltInFunctions.BinaryFindObject (&opc);
    } else {
        opCode = -opNo-1;
    }

    if (opCode<0) {
        WarnError(_String ("Operation: '") & opc &"' is not defined." );
        opCode = 0;
    }

    numberOfTerms = opNo;
    theData       = -1;
    theNumber     = nil;
}
//__________________________________________________________________________________
_Operation::_Operation  (_PMathObj theObj)
// construct the operation by its symbol and, if relevant -
// number of operands
{
    numberOfTerms = 0;
    theData       = -1;
    opCode        = -1;
    theNumber     = theObj;
}

//__________________________________________________________________________________
bool _Operation::CanResultsBeCached (_Operation* prev, bool exp_only)
{
    if (theNumber == nil && theData == -1 && numberOfTerms == 1) {
        if ((prev->theNumber && prev->theNumber->ObjectClass() == MATRIX)
           || (prev->theData >= 0 && LocateVar (prev->theData)->ObjectClass () == MATRIX)) {
            if (! exp_only || opCode == HY_OP_CODE_EXP)
                return true;
        }
    }
    return false;
}

//__________________________________________________________________________________

bool _Operation::HasChanged (void)
{
    if (theNumber) {
        return theNumber->HasChanged();
    }
    if (theData >= 0) {
        return LocateVar (GetAVariable())->HasChanged();
    }

    return false;
}

//__________________________________________________________________________________
_Operation::_Operation  (bool isVar, _String& stuff, bool isG, _VariableContainer* theParent, bool take_a_reference)
{
    if (isVar) { // creating a variable
        long f;
        _String theS (stuff);
        if (theParent/*&&(!isG)*/) { // 20070620: SLKP the commenting may break default behavior!
            f = LocateVarByName(theS);

            if (f>=0 && !FetchVar(f)->IsGlobal()) {
                f = -1;
            }

            if (f<0) {
                theS = (*theParent->theName)&"."&theS;
            }
        }

        f = LocateVarByName(theS);

        if (f<0) {
            _Variable v (theS, isG);
            f = v.theIndex;
        } else {
            f = variableNames.GetXtra(f);
        }

        theData       = f;
        theNumber     = nil;
        numberOfTerms = take_a_reference?(1):0;
        
    } else {
        numberOfTerms = 0;
        if (stuff.Equal (&noneToken))
            theNumber = new _MathObject;            
        else
            theNumber = new _Constant (stuff);
        theData = -1;
    }
    opCode = -1;

}


//__________________________________________________________________________________
_Operation::~_Operation (void)
{

    if (theNumber) {
        DeleteObject (theNumber);
    }
}

//__________________________________________________________________________________
bool _Operation::IsAVariable(bool deep)
{
    if (theData==-1) {
        if (deep&&theNumber) {
            return theNumber->IsVariable();
        }
        return false;
    }
    return true;

}

//__________________________________________________________________________________
bool _Operation::IsAFunctionCall (void)
{
    return theData == -1 && numberOfTerms < 0;

}

//__________________________________________________________________________________
bool _Operation::IsConstant (void)
{
    if (theData==-1) {
        if (theNumber) {
            return theNumber->IsConstant();
        }

        return !(opCode == HY_OP_CODE_BRANCHLENGTH||opCode == HY_OP_CODE_RANDOM||opCode == HY_OP_CODE_TIME);
    }
    return LocateVar(GetAVariable())->IsConstant();

}

//__________________________________________________________________________________
bool        _Operation::EqualOp (_Operation* otherOp)
{
    if (theNumber) {
        if (otherOp->theNumber) {
            long oc = theNumber->ObjectClass();

            if ((oc == NUMBER) && (oc==otherOp->theNumber->ObjectClass())) {
                return CheckEqual (theNumber->Value(), otherOp->theNumber->Value());
            }
        }
        return false;
    } else if (theData==-1) {
        if (numberOfTerms<0) {
            return numberOfTerms == otherOp->numberOfTerms;
        } else {
            return opCode == otherOp->opCode;
        }
    } else {
        return theData == otherOp->theData;
    }

    return false;
}

//__________________________________________________________________________________


long     _Operation::BinOpCode           (_String const & op_token, long index) {
  if (index >= 0)
    return BinOps.Find (op_token.getChar(index-1L) * 256L + op_token.getChar (index));
  
  return BinOps.Find (op_token.sLength == 2 ? op_token.getChar(0) * 256L + op_token.getChar (1) : op_token.getChar (0));
}


//__________________________________________________________________________________

bool        _Operation::ReportOperationExecutionError(_String text, _String * errMsg) {
    _String theError = text & ". "; 
    
    if (errMsg) {
        *errMsg = theError;
    } else {
        WarnError (theError);
    }
    
    return false;
}

//__________________________________________________________________________________
bool        _Operation::Execute (_Stack& theScrap, _VariableContainer* nameSpace, _String* errMsg) {
  if (theNumber) { // push value
    theScrap.Push(theNumber);
    return true;
  }
  
  if (theData >= 0L) { // variable reference
    if (numberOfTerms <= 0L) { // compute and push value
      theScrap.Push(((_Variable*)((BaseRef*)variablePtrs.lData)[theData])->Compute());
    } else { // compute reference and push value
      theScrap.Push (((_Variable*)((BaseRef*)variablePtrs.lData)[theData])->ComputeReference(nameSpace),false);
    }
    return true;
  }
  if (theData < -2L) { // WTF is is this?
    theScrap.Push(((_Variable*)((BaseRef*)variablePtrs.lData)[-theData-3])->GetValue());
    return true;
  }
  if (numberOfTerms<0) { // execute a user-defined function
      
      if (!IsBFFunctionIndexValid(opCode)) {
        return ReportOperationExecutionError ("Attempted to call an undefined/deleted user function ", errMsg);
      }
      
      long arguments = GetBFFunctionArgumentCount (opCode);
      
      if (theScrap.StackDepth()<arguments) {
        return ReportOperationExecutionError (_String("User-defined function:") &
                                              & GetBFFunctionNameByIndex (opCode)
                                              &" needs "&_String(long(arguments))& " parameters: "&_String(theScrap.StackDepth())&" were supplied ", errMsg);
      }
      
      _List       displacedVars,
      *funcVarList = &GetBFFunctionArgumentList(opCode),
      displacedValues,
      referenceArgs;
      
      _SimpleList existingIVars,
      existingDVars,
      displacedReferences,
      *funcVarTypes = &GetBFFunctionArgumentTypes (opCode);
      
      bool        need_to_purge = false;
    
      //printf ("***** Calling %s\n", GetBFFunctionNameByIndex (opCode).sData);
    
      for (long k = arguments-1L; k >= 0; k--) {
        bool            isRefVar = (funcVarTypes->Element (k) == BL_FUNCTION_ARGUMENT_REFERENCE);
        
        _String         *argument_k = (_String*)funcVarList->Element(k);
        _PMathObj       nthterm = theScrap.Pop();
        
        //printf ("\tArgument %d : %s (%s)\n", k, argument_k->sData, _String((_String*)nthterm->toStr()).sData);

        if (isRefVar) {
          need_to_purge = true;
          if (nthterm->ObjectClass()!=STRING) {
            _FString * type = (_FString*)nthterm->Type();
            _String errText = _String ("User-defined function '")
            & GetBFFunctionNameByIndex(opCode)
            &"' expected a string for the reference variable '"
            & *argument_k
            &"' but was passed a " & *type->theString
            & " with the value of " & _String((_String*)nthterm->toStr());
            
            DeleteObject (type);
            return ReportOperationExecutionError (errText, errMsg);
          }
        }
        
        _Variable* argument_var = CheckReceptacle (argument_k, empty, false, false);
        
        if (!isRefVar) {
          if (argument_var->IsIndependent()) {
            // if the variable exists and is independent then
            // simply swap the value of the var, otherwise
            // duplicate the entire variable
            argument_var->varFlags &= HY_VARIABLE_SET;
            if (!argument_var->varValue) {
              argument_var->Compute();
            }
            displacedValues<<argument_var->varValue;
            argument_var->varFlags |= HY_VARIABLE_CHANGED;
            argument_var->varValue = nthterm;
            existingIVars<<argument_var->GetAVariable();
          } else {
            _Variable *newV = new _Variable (*argument_k);
            newV->SetValue(nthterm,false);
            nthterm->AddAReference();
            existingDVars<<argument_var->GetAVariable();
            displacedVars<<argument_k;
            variablePtrs.Replace (argument_var->GetAVariable(),newV,false);
          }
        } else {
          
          long new_index = argument_var->GetAVariable();
          
          referenceArgs << argument_var->GetName();
          displacedReferences<<new_index;
          
          _String * refArgName = ((_FString*)nthterm)->theString;
          
          if (nameSpace) {
            *refArgName = AppendContainerName (*refArgName, nameSpace);
          }
          
          _Variable* reference_var = CheckReceptacle (refArgName, empty, false, false);
          
          variableNames.SetXtra (new_index, reference_var->GetAVariable());
          
          DeleteObject (nthterm);
        }
      }
      
      _ExecutionList * function_body = &GetBFFunctionBody(opCode);
      
      if (need_to_purge) {
        function_body->ResetFormulae();
      }
      
      if (currentExecutionList && currentExecutionList->stdinRedirect) {
        function_body -> stdinRedirect    = currentExecutionList->stdinRedirect;
        function_body -> stdinRedirectAux = currentExecutionList->stdinRedirectAux;
      }
      
      _PMathObj ret = function_body->Execute();
      
      function_body -> stdinRedirect    = nil;
      function_body -> stdinRedirectAux = nil;
      
      if (terminateExecution) {
        theScrap.Push (new _Constant (0.0));
        return true;
      }
      
      if (ret) {
        theScrap.Push (ret);
      }
    
    //printf ("\nFunction result = %s\n", _String ((_String*)theScrap.Pop (false)->toStr()).getStr());
    
      for (unsigned long di = 0UL; di < referenceArgs.lLength; di++) {
        variableNames.SetXtra(LocateVarByName (*(_String*)referenceArgs(di)),displacedReferences.lData[di]);
      }
      
      for (unsigned long dv = 0UL; dv < displacedVars.lLength; dv++) {
        variablePtrs.Replace (existingDVars.lData[dv],(_PMathObj)displacedVars(dv));
      }
      
      
      for (unsigned long dv2 = 0; dv2 < displacedValues.lLength; dv2++) {
        _Variable* theV = LocateVar (existingIVars.lData[dv2]);
        DeleteObject(theV->varValue);
        theV->varValue = ((_PMathObj)displacedValues(dv2));
      }
    
      return true;
    
  }
  
  if (theScrap.theStack.lLength<numberOfTerms) {
    return ReportOperationExecutionError (_String((_String*)toStr())&
                                          " needs "&_String(numberOfTerms)& " arguments; "&_String(theScrap.StackDepth())&" were supplied", errMsg);
    
  }
  
  _PMathObj term1, term2 = nil, term3 = nil, temp;
  
  _hyExecutionContext localContext (nameSpace, errMsg);
  
  
  if (numberOfTerms >= 3) {
    long sL = theScrap.theStack.lLength-1;
    term3 = (_PMathObj)theScrap.theStack.lData[sL--];
    term2 = (_PMathObj)theScrap.theStack.lData[sL--];
    term1 = (_PMathObj)theScrap.theStack.lData[sL];
    theScrap.theStack.lLength = sL;
    temp = term1->Execute (opCode, term2, term3, &localContext);
    DeleteObject (term1);
    DeleteObject (term2);
    DeleteObject (term3);
  } else if (numberOfTerms == 2) {
    long sL = theScrap.theStack.lLength-1;
    term2 = (_PMathObj)theScrap.theStack.lData[sL--];
    term1 = (_PMathObj)theScrap.theStack.lData[sL];
    theScrap.theStack.lLength = sL;
    temp = term1->Execute (opCode, term2, nil, &localContext);
    DeleteObject (term1);
    DeleteObject (term2);
  } else {
    term1 = (_PMathObj)theScrap.theStack.lData[--theScrap.theStack.lLength];
    temp = term1->Execute (opCode, nil, nil, &localContext);
    DeleteObject (term1);
  }
  
  
  if (temp) {
    theScrap.theStack.Place(temp);
    return true;
  } else {
    return false;
  }
  
}

//__________________________________________________________________________________
void        _Operation::StackDepth (long& depth)
{
    if (theNumber || theData > -1 || theData < -2) {
        depth++;
        return;
    }

    if (numberOfTerms<0) { // execute a user-defined function
        depth -= GetBFFunctionArgumentCount(opCode) - 1L;
        return;
    }
    depth -= numberOfTerms - 1;
}


//__________________________________________________________________________________
bool        _Operation::ExecutePolynomial (_Stack& theScrap, _VariableContainer* nameSpace, _String* errMsg)
{
    if (theData<=-2 || numberOfTerms < 0) {
        return false;
    }

    _Polynomial*p = nil;
    if (theNumber) {
        p= (_Polynomial*)checkPointer(new _Polynomial(theNumber->Value()));
    }

    if (theData>-1) {
        p= (_Polynomial*)checkPointer(new _Polynomial(*LocateVar(theData>-1?theData:-theData-2)));
    }

    if (p) {
        theScrap.Push(p);
        DeleteObject(p);
        return true;
    }

    if (theScrap.StackDepth()<numberOfTerms) {
        _String s((_String*)toStr());
        return ReportOperationExecutionError (s&
                   " needs "&_String(numberOfTerms)& " arguments. Only "&_String(theScrap.StackDepth())&" were given", errMsg);
    }

    _PMathObj term1,
              term2 = nil,
              temp;

    bool      opResult = true;

    if (numberOfTerms == 2) {
        term2 = theScrap.Pop();
    }

    _hyExecutionContext localContext (nameSpace, errMsg);
    term1 = theScrap.Pop();
    temp  = term1->Execute (opCode, term2, nil, &localContext);
    DeleteObject (term1);

    if (temp) {
        theScrap.Push (temp, false);
    } else {
        opResult = false;
    }

    if (term2) {
        DeleteObject (term2);
    }

    return opResult;
}

