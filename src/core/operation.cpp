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

#include "defines.h"
#include "variable.h"
#include "operation.h"

#include "parser.h"
#include "polynoml.h"
#include "batchlan.h"
#include "parser.h"
#include "global_things.h"

using namespace hy_global;


extern _SimpleList BinOps,
       opPrecedence,
       FunctionArgumentCount,
       associativeOps;

//__________________________________________________________________________________

_Operation::_Operation  (void) {
    Initialize();
 }

//__________________________________________________________________________________
void    _Operation::Initialize(bool)
{
    numberOfTerms = 0;
    theData = -1;
    theNumber = nil;
    opCode = -1;
}

//__________________________________________________________________________________
BaseRef _Operation::makeDynamic (void) const {
    _Operation * res = new _Operation;
    res->Duplicate(this);
    return res;
}

//__________________________________________________________________________________

_Operation::_Operation (_Operation const& rhs) {
    numberOfTerms  = rhs.numberOfTerms;
    theData        = rhs.theData;
    theNumber      = rhs.theNumber;
    opCode         = rhs.opCode;
    if (theNumber) {
        theNumber->AddAReference();
    }
}


//__________________________________________________________________________________
void    _Operation::Duplicate(BaseRefConst r) {
    _Operation const * o = (_Operation const*)r;
    numberOfTerms  = o->numberOfTerms;
    theData        = o->theData;
    theNumber      = o->theNumber;
    opCode         = o->opCode;
    if (theNumber) {
        theNumber->AddAReference();
    }
}


//__________________________________________________________________________________
BaseRef _Operation::toStr (unsigned long) {


    if (theData != -1) {
        return new _String (_String("Variable ")& *LocateVar(GetAVariable())->GetName());
    } else if (theNumber) {
        _FString * type = (_FString*)theNumber->Type();
        _String res = _String("Constant (")& *type->theString & ")" & _String((_String*)theNumber->toStr());
        DeleteObject(type);
        return new _String (res);
    } else {
        if (IsHBLFunctionCall())
          return new _String (_String ("Call HBL function ") & GetBFFunctionNameByIndex(UserFunctionID()));
        else
          return new _String (_String("Operation ") & *(_String*)BuiltInFunctions(opCode) & " with " & _String((long)numberOfTerms) & " arguments");
    }


    return new _String ("This should not happen");

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
_Operation::_Operation  (_String const& opc, const long opNo = 2) {
// construct the operation by its symbol and, if relevant -
// number of operands
    if(opNo>=0) {
        opCode = BuiltInFunctions.BinaryFindObject (&opc);
    } else {
        opCode = -opNo-1;
    }

    if (opCode<0L) {
        HandleApplicationError (_String ("Operation: ") & opc.Enquote() &" is not defined." );
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

bool   _Operation::IsConstantOfType   (const long type) const {
  if (theNumber && (theNumber->ObjectClass() & type)) {
    return true;
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
_Operation::_Operation  (bool isVar, _String& stuff, bool isG, _VariableContainer const* theParent, bool take_a_reference)
{
    if (isVar) { // creating a variable
        long f;
        _String theS (stuff);
        if (theParent) {
            f = LocateVarByName(theS);

            if (f>=0L && !FetchVar(f)->IsGlobal()) {
                f = -1L;
            }

            if (f<0L) {
                theS = (*theParent->theName)&"."&theS;
            }
        }

        f = LocateVarByName(theS);

        if (f<0L) {
            _Variable v (theS, isG);
            f = v.theIndex;
        } else {
            f = variableNames.GetXtra(f);
        }

        theData       = f;
        theNumber     = nil;
        numberOfTerms = take_a_reference?(1):0;

    } else {
        numberOfTerms = 0L;
        if (stuff.Equal (&noneToken))
            theNumber = new _MathObject;
        else
            theNumber = new _Constant (stuff);
        theData = -1L;
    }
    opCode = -1L;

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
    if (theData==-1 || theData == -2) {
        if (deep&&theNumber) {
            return theNumber->IsVariable();
        }
        return false;
    }
    return true;

}

//__________________________________________________________________________________
bool _Operation::IsHBLFunctionCall (void) const {
    return theData == -1 && numberOfTerms < 0;
}

//__________________________________________________________________________________
long _Operation::GetHBLFunctionID (void) const {
  return opCode;
}


//__________________________________________________________________________________
bool _Operation::IsConstant (bool strict)
{
    if (theData==-1) {
        if (theNumber) {
            return theNumber->IsConstant();
        }

        return !(opCode == HY_OP_CODE_BRANCHLENGTH||opCode == HY_OP_CODE_RANDOM||opCode == HY_OP_CODE_TIME);
    }
    if (strict) {
      return false;
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
    return BinOps.Find (op_token.get_char(index-1L) * 256L + op_token.get_char (index));
  
  return BinOps.Find (op_token.length () == 2 ? op_token.get_char(0) * 256L + op_token.get_char (1) : op_token.get_char (0));
}


//__________________________________________________________________________________

bool        _Operation::ReportOperationExecutionError(_String text, _String * errMsg) {
    _String theError = text & ". ";

    if (errMsg) {
        *errMsg = theError;
    } else {
        HandleApplicationError (theError);
    }

    return false;
}

//__________________________________________________________________________________
long        _Operation::StackDepth (void) const {
  if (theNumber || theData != -1L) {
    return 1L;
  }

  if (numberOfTerms<0L) { // execute a user-defined function
    return (numberOfTerms + 2); // assumes a return value
  }

  return (1-numberOfTerms);
}


//__________________________________________________________________________________
bool        _Operation::Execute (_Stack& theScrap, _VariableContainer const* nameSpace, _String* errMsg) {
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
  if (theData < -2L) { // place variable value (no compute, i.e. pass by reference)
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
                                              GetBFFunctionNameByIndex (opCode)
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
        bool            isRefVar = (funcVarTypes->Element (k) == kBLFunctionArgumentReference);
        
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
        
        _Variable* argument_var = CheckReceptacle (argument_k, kEmptyString, false, false);
        
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
            displacedVars<<argument_var; // 2 references
            argument_var->AddAReference(); // 3 references
            variablePtrs.Replace (argument_var->GetAVariable(),newV,false); // 2 references
          }
        } else {

          long new_index = argument_var->GetAVariable();

          referenceArgs << argument_var->GetName();
          displacedReferences<<new_index;

          _String * refArgName = ((_FString*)nthterm)->theString;

          if (nameSpace) {
            *refArgName = AppendContainerName (*refArgName, nameSpace);
          }
          
          _Variable* reference_var = CheckReceptacle (refArgName, kEmptyString, false, false);
          
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
      
      if (terminate_execution) {
        theScrap.Push (new _Constant (0.0));
        return true;
      }

      if (ret) {
        theScrap.Push (ret);
      } else {
        theScrap.Push (new _MathObject);
      }

    //printf ("\nFunction result = %s\n", _String ((_String*)theScrap.Pop (false)->toStr()).getStr());

      for (unsigned long di = 0UL; di < referenceArgs.lLength; di++) {
        variableNames.SetXtra(LocateVarByName (*(_String*)referenceArgs(di)),displacedReferences.lData[di]);
      }

      for (unsigned long dv = 0UL; dv < displacedVars.lLength; dv++) {
        variablePtrs.Replace (existingDVars.lData[dv],(_PMathObj)displacedVars(dv), false);
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

  _PMathObj arg0 = ((_PMathObj)theScrap.theStack.lData[theScrap.theStack.lLength-numberOfTerms]),
            temp;

  _hyExecutionContext localContext (nameSpace, errMsg);

  if (numberOfTerms > 1) {
    _List arguments;

    for (long k = numberOfTerms-1; k >= 1; k --) {
      arguments.AppendNewInstance ((_PMathObj)theScrap.theStack.lData[theScrap.theStack.lLength-k]);
    }

    temp =  arg0->ExecuteSingleOp(opCode, &arguments, &localContext);
    theScrap.theStack.lLength -= (numberOfTerms);

  } else {
    temp =  ((_PMathObj)theScrap.theStack.lData[theScrap.theStack.lLength-1])->ExecuteSingleOp(opCode, nil, &localContext);
    theScrap.theStack.lLength--;
  }

  DeleteObject (arg0);


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
        p= new _Polynomial(theNumber->Value());
    }

    if (theData>-1) {
        p=  new _Polynomial(*LocateVar(theData>-1?theData:-theData-2));
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


    _PMathObj arg0 = ((_PMathObj)theScrap.theStack.lData[theScrap.theStack.lLength-numberOfTerms]),
    temp;

    _hyExecutionContext localContext (nameSpace, errMsg);

    if (numberOfTerms > 1) {
      _List arguments;

      for (long k = numberOfTerms-1; k >= 1; k --) {
        arguments.AppendNewInstance ((_PMathObj)theScrap.theStack.lData[theScrap.theStack.lLength-k]);
      }

      temp =  arg0->ExecuteSingleOp(opCode, &arguments, &localContext);
      theScrap.theStack.lLength -= (numberOfTerms);

    } else {
      temp =  ((_PMathObj)theScrap.theStack.lData[theScrap.theStack.lLength-1])->ExecuteSingleOp(opCode, nil, &localContext);
      theScrap.theStack.lLength--;
    }

    DeleteObject (arg0);


    if (temp && temp->ObjectClass() != HY_UNDEFINED ) {
      theScrap.theStack.Place(temp);
      return true;
    } else {
      return false;
    }


}

