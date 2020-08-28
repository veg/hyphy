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

_SimpleList     _Operation::ListOfInverseOps;


//__________________________________________________________________________________

_Operation::_Operation  (void) {
    Initialize();
 }

//__________________________________________________________________________________
void    _Operation::Initialize(bool) {
    numberOfTerms = 0;
    theData = -1;
    theNumber = nil;
    cachedResult = nil;
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
    cachedResult   = rhs.cachedResult;
    if (theNumber) {
        theNumber->AddAReference();
    }
    if (cachedResult) {
        cachedResult->AddAReference();
    }
}



//__________________________________________________________________________________
void    _Operation::Duplicate(BaseRefConst r) {
    _Operation const * o = (_Operation const*)r;
    numberOfTerms  = o->numberOfTerms;
    theData        = o->theData;
    theNumber      = o->theNumber;
    opCode         = o->opCode;
    cachedResult   = o->cachedResult;
    if (theNumber) {
        theNumber->AddAReference();
    }
    if (cachedResult) {
        cachedResult->AddAReference();
    }
}

//__________________________________________________________________________________
void    _Operation::operator = (_Operation const & rhs) {
    Duplicate (&rhs);
}


//__________________________________________________________________________________
BaseRef _Operation::toStr (unsigned long) {
    
    if (theData != -1) {
        return new _String (_String("Variable ")& *LocateVar(GetAVariable())->GetName());
    } else if (theNumber) {
        _FString * type = (_FString*)theNumber->Type();
        _String res = _String("Constant (")& type->get_str() & ")" & _String((_String*)theNumber->toStr());
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
_Operation::_Operation  (const long theCode, const long opNo = 2) {
// by opcode
    opCode = theCode;
    numberOfTerms = opNo;
    theData       = -1;
    theNumber     = nil;
    cachedResult  = nil;
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
    cachedResult  = nil;
}
//__________________________________________________________________________________
_Operation::_Operation  (HBLObjectRef theObj) {
    numberOfTerms = 0;
    theData       = -1;
    opCode        = -1;
    theNumber     = theObj;
    cachedResult = nil;
}

//__________________________________________________________________________________
_Operation::_Operation  (_Variable const & v) {
    numberOfTerms = 0;
    theData       = v.get_index();
    opCode        = -1;
    theNumber     = nil;
    cachedResult  = nil;
}

//__________________________________________________________________________________
bool _Operation::CanResultsBeCached (_Operation* prev, bool exp_only) {
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

bool _Operation::HasChanged (void) {
    if (theNumber) {
        return theNumber->HasChanged();
    }
    if (theData >= 0) {
        return LocateVar (GetAVariable())->HasChanged();
    }

    return false;
}

//__________________________________________________________________________________
_Operation::_Operation  (bool isVar, _String& stuff, bool isG, _VariableContainer const* theParent, bool take_a_reference) {
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
        if (stuff == kNoneToken || stuff == kNullToken)
            theNumber = new _MathObject;
        else
            theNumber = new _Constant (stuff);
        theData = -1L;
    }
    opCode = -1L;
    cachedResult = nil;
}


//__________________________________________________________________________________
_Operation::~_Operation (void) {
   DeleteAndZeroObject (theNumber);
   DeleteAndZeroObject (cachedResult);
}

//__________________________________________________________________________________
bool _Operation::IsAVariable(bool deep) {
    if (theData == -1 || theData == -2) {
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
    return BinOps.Find (op_token.char_at (index-1L) * 256L + op_token.char_at (index));
  
  return BinOps.Find (op_token.length () == 2 ? op_token.char_at(0) * 256L + op_token.char_at (1) : op_token.char_at (0));
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
bool        _Operation::Execute (_Stack& theScrap, _VariableContainer const* nameSpace, _String* errMsg, bool canCache) {
  if (theNumber) { // push value
    theScrap.Push(theNumber);
    return true;
  }

  if (theData >= 0L) { // variable reference
    if (numberOfTerms <= 0L) { // compute and push value
      theScrap.Push(((_Variable*)((BaseRef*)variablePtrs.list_data)[theData])->Compute());
    } else { // compute reference and push value
      theScrap.Push (((_Variable*)((BaseRef*)variablePtrs.list_data)[theData])->ComputeReference(nameSpace),false);
    }
    return true;
  }
  if (theData < -2L) { // place variable value (no compute, i.e. pass by reference)
    theScrap.Push(((_Variable*)((BaseRef*)variablePtrs.list_data)[-theData-3])->GetValue());
    return true;
  }
  if (numberOfTerms<0) { // execute a user-defined function

      if (!IsBFFunctionIndexValid(opCode)) {
        return ReportOperationExecutionError ("Attempted to call an undefined/deleted user function ", errMsg);
      }

      long arguments = GetBFFunctionArgumentCount (opCode);

      if (theScrap.StackDepth()<arguments) {
        return ReportOperationExecutionError (_String("User-defined function ") &
                                              GetBFFunctionNameByIndex (opCode).Enquote()
                                              &" needs "&_String(long(arguments))& " parameters, but "&_String(theScrap.StackDepth())&" were supplied ", errMsg);
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
        HBLObjectRef       nthterm = theScrap.Pop();

        //printf ("\tArgument %d : %s (%s)\n", k, argument_k->get_str(), _String((_String*)nthterm->toStr()).get_str());

        if (isRefVar) {
          need_to_purge = true;
          if (nthterm->ObjectClass()!=STRING) {
            _FString * type = (_FString*)nthterm->Type();
            _String errText = _String ("User-defined function '")
            & GetBFFunctionNameByIndex(opCode)
            &"' expected a string for the reference variable '"
            & *argument_k
            &"' but was passed a " & type->get_str()
            & " with the value of " & _String((_String*)nthterm->toStr());

            DeleteObject (type);
            return ReportOperationExecutionError (errText, errMsg);
          }
        }
        
        _Variable* argument_var = CheckReceptacle (argument_k, kEmptyString, false, false, false);
        
        if (!isRefVar) {
          if (argument_var->IsIndependent() && (argument_var->ObjectClass() & (TREE|TOPOLOGY)) == 0) {
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
            existingIVars<<argument_var->get_index();
          } else {
            _Variable *newV = new _Variable (*argument_k);
            newV->SetValue(nthterm,false,true, NULL);
            existingDVars<<argument_var->get_index();
            displacedVars<<argument_var; // 2 references
            argument_var->AddAReference(); // 3 references
            variablePtrs.Replace (argument_var->get_index(),newV,false); // 2 references
          }
        } else {

          long new_index = argument_var->get_index();

          referenceArgs << argument_var->GetName();
          displacedReferences<<new_index;

          _String const * refArgName = &((_FString*)nthterm)->get_str();
            
          _Variable* reference_var;
          // TODO SLKP 20170928 Check that this doesn't break reference variable calls
          // used to write namespaced variable to nthterm
          if (nameSpace) {
              reference_var =  CheckReceptacle (&AppendContainerName (*refArgName, nameSpace), kEmptyString, false, false, false);
          } else {
              reference_var =  CheckReceptacle (refArgName, kEmptyString, false, false, false);
          }
          
          variableNames.SetXtra (new_index, reference_var->get_index());
          DeleteObject (nthterm);
        }
      }

      _ExecutionList * function_body = &GetBFFunctionBody(opCode);

      if (need_to_purge) {
        function_body->ResetFormulae();
      }

      HBLObjectRef ret;
      
      bool update_kw = false;
      if (currentExecutionList && (currentExecutionList->has_stdin_redirect() || currentExecutionList->has_keyword_arguments())) {
          // 20180620: SLKP, need to split this off because if Execute fails
          // then there will be a double free on stdinRedirect

        auto stash1 = currentExecutionList->stdinRedirect;
        auto stash2 = currentExecutionList->stdinRedirectAux;
        auto stash_kw_tags = currentExecutionList->kwarg_tags;
        auto stash_kw = currentExecutionList->kwargs;
          
        if (currentExecutionList->has_stdin_redirect()) {
            function_body -> stdinRedirect    = currentExecutionList->stdinRedirect;
            function_body -> stdinRedirectAux = currentExecutionList->stdinRedirectAux;
          
            currentExecutionList->stdinRedirect->AddAReference();
            currentExecutionList->stdinRedirectAux->AddAReference();
        }
          
          
        if (currentExecutionList->has_keyword_arguments()) {
            function_body -> kwarg_tags    = currentExecutionList->kwarg_tags;
            function_body -> kwargs = currentExecutionList->kwargs;
            function_body -> currentKwarg = currentExecutionList->currentKwarg;
            
            if (stash_kw_tags) currentExecutionList->kwarg_tags->AddAReference();
            if (stash_kw) currentExecutionList->kwargs->AddAReference();
            update_kw = true;
            
        }

        ret = function_body->Execute();
          
        if (stash1) stash1 -> RemoveAReference();
        if (stash2) stash2 -> RemoveAReference();
        if (stash_kw_tags) stash_kw_tags->RemoveAReference();
        if (stash_kw) stash_kw->RemoveAReference();

      } else {
          ret = function_body->Execute();
      }

      function_body -> stdinRedirect    = nil;
      function_body -> stdinRedirectAux = nil;
      if (update_kw) {
          function_body -> kwarg_tags    = nil;
          function_body -> kwargs = nil;
          currentExecutionList->currentKwarg = function_body->currentKwarg;
      }

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
        variableNames.SetXtra(LocateVarByName (*(_String*)referenceArgs(di)),displacedReferences.list_data[di]);
      }

      for (unsigned long dv = 0UL; dv < displacedVars.lLength; dv++) {
        variablePtrs.Replace (existingDVars.list_data[dv],(HBLObjectRef)displacedVars(dv), false);
      }


      for (unsigned long dv2 = 0; dv2 < displacedValues.lLength; dv2++) {
        _Variable* theV = LocateVar (existingIVars.list_data[dv2]);
        DeleteObject(theV->varValue);
        theV->varValue = ((HBLObjectRef)displacedValues(dv2));
      }

      return true;

  }

  if (theScrap.theStack.lLength<numberOfTerms) {
    return ReportOperationExecutionError (_String((_String*)toStr())&
                                          " needs "&_String(numberOfTerms)& " arguments; "&_String(theScrap.StackDepth())&" were supplied", errMsg);

  }

  HBLObjectRef arg0 = ((HBLObjectRef)theScrap.theStack.list_data[theScrap.theStack.lLength-numberOfTerms]),
            temp;

  _hyExecutionContext localContext (nameSpace, errMsg);

  if (numberOfTerms > 1) {
    _List arguments;

    if (numberOfTerms == 2) {
        arguments.AppendNewInstance ((HBLObjectRef)theScrap.theStack.list_data[theScrap.theStack.lLength - 1]);
    } else {
        for (long k = 1; k < numberOfTerms; k ++) {
          arguments.AppendNewInstance ((HBLObjectRef)theScrap.theStack.list_data[theScrap.theStack.lLength - numberOfTerms + k]);
        }
    }
    
    if (canCache) {
        temp =  arg0->ExecuteSingleOp(opCode, &arguments, &localContext, cachedResult);
        if (temp != cachedResult) {
            DeleteObject(cachedResult);
            cachedResult = temp;
        }
        temp->AddAReference();
    } else {
        temp =  arg0->ExecuteSingleOp(opCode, &arguments, &localContext);
    }
    theScrap.theStack.lLength -= (numberOfTerms);

  } else {
    temp =  ((HBLObjectRef)theScrap.theStack.list_data[theScrap.theStack.lLength-1])->ExecuteSingleOp(opCode, nil, &localContext);
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
bool        _Operation::ExecutePolynomial (_Stack& theScrap, _VariableContainer* nameSpace, _String* errMsg) {
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


    HBLObjectRef arg0 = ((HBLObjectRef)theScrap.theStack.list_data[theScrap.theStack.lLength-numberOfTerms]),
    temp;

    _hyExecutionContext localContext (nameSpace, errMsg);

    if (numberOfTerms > 1) {
      _List arguments;

      for (long k = numberOfTerms-1; k >= 1; k --) {
        arguments.AppendNewInstance ((HBLObjectRef)theScrap.theStack.list_data[theScrap.theStack.lLength-k]);
      }

      temp =  arg0->ExecuteSingleOp(opCode, &arguments, &localContext);
      theScrap.theStack.lLength -= (numberOfTerms);

    } else {
      temp =  ((HBLObjectRef)theScrap.theStack.list_data[theScrap.theStack.lLength-1])->ExecuteSingleOp(opCode, nil, &localContext);
      theScrap.theStack.lLength--;
    }

    DeleteObject (arg0);


    if (temp && temp->ObjectClass() != HY_UNDEFINED ) {
      theScrap.theStack.Place(temp);
      return true;
    } else {
     DeleteObject (temp);
     return false;
    }


}

