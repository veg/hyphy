//__________________________________________________________________________________

_Operation::_Operation  (void)
{
    numberOfTerms = 0;
    theData = -1;
    theNumber = nil;
}

//__________________________________________________________________________________

void    _Operation::Initialize(void)
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
    _String res, *dump = nil;
    if (theData!=-1) {
        dump = (_String*)((_Variable*)LocateVar(theData))->toStr();
        res = _String("Variable ")& *dump;
    } else if (theNumber) {
        dump = (_String*)theNumber->toStr();
        res = _String("Constant ")& *dump;
    } else {
        res = _String("Operation ")&*(_String*)BuiltInFunctions(opCode);
    }

    if(dump) {
        DeleteObject (dump);
    }
    return res.makeDynamic();

}
//__________________________________________________________________________________
_Operation::_Operation  (_String& opc, long opNo = 2)
// construct the operation by its symbol and, if relevant -
// number of operands
{
    if(opNo>=0) {
        opCode = BuiltInFunctions.BinaryFind(&opc);
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

bool _Operation::CanResultsBeCached (_Operation* prev)
{
    if (theNumber == nil && theData == -1 && numberOfTerms == 1) {
        if (prev->theNumber && prev->theNumber->ObjectClass() == MATRIX || prev->theData >= 0 && LocateVar (prev->theData)->ObjectClass () == MATRIX) {
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
_Operation::_Operation  (bool isVar, _String& stuff, bool isG, _VariableContainer* theParent)
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
                /*f = stuff.Find('.');
                do
                {
                    if (f!=-1)
                    {
                        f = stuff.Find(".",f+1,-1); // skip scope levels as needed
                        continue;
                    }
                    else
                        break;
                }
                while ((theParent = theParent->theParent));
                if (theParent)*/
                theS = (*theParent->theName)&"."&theS;
            }
            /*StringToConsole(stuff);
            NLToConsole();
            StringToConsole(theS);
            NLToConsole();*/
        }

        f = LocateVarByName(theS);

        if (f<0) {
            _Variable v (theS, isG);
            f = v.theIndex;
        } else {
            f = variableNames.GetXtra(f);
        }

        theData = f;
        theNumber = nil;
        numberOfTerms = 0;
    } else {
        numberOfTerms = 0;
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

bool        _Operation::Execute (_Stack& theScrap, _VariableContainer* nameSpace)
{
    if (theNumber) {
        theScrap.Push(theNumber);
        return true;
    }
    if (theData>-1) {
        theScrap.Push(((_Variable*)((BaseRef*)variablePtrs.lData)[theData])->Compute());
        return true;
    }
    if (theData<-2) {
        theScrap.Push(((_Variable*)((BaseRef*)variablePtrs.lData)[-theData-3])->GetValue());
        return true;
    }
    if (numberOfTerms<0) { // execute a user-defined function

        {
            long functionID = -numberOfTerms-1;
            if (functionID >= batchLanguageFunctionParameters.lLength) {
                WarnError ("Attempted to call an undefined user function.");
                return false;
            }
            long numb = batchLanguageFunctionParameters(functionID);
            if (theScrap.StackDepth()<numb) {
                _String errMsg ("User-defined function:");
                errMsg = errMsg&*(_String*)batchLanguageFunctionNames(functionID)
                         &" needs "&_String(numb)& " parameters, but only "&_String(theScrap.StackDepth())&" were supplied.";
                WarnError (errMsg);
                return false;
            }

            _List       displacedVars,
                        *funcVarList = (_List*)batchLanguageFunctionParameterLists(functionID),
                         displacedValues,
                         referenceArgs;

            _SimpleList existingIVars,
                        existingDVars,
                        displacedReferences;

            _String     *argNameString;
            long        i;
            bool        purgeFlas = false;

            for (long k=numb-1; k>=0; k--) {
                bool            isRefVar = false;
                argNameString = (_String*)(*funcVarList)(k);

                if (argNameString->sData[argNameString->sLength-1]=='&') {
                    argNameString->Trim(0,argNameString->sLength-2);
                    isRefVar  = true;
                    purgeFlas = true;
                }

                long      f = LocateVarByName (*argNameString);

                _PMathObj nthterm = theScrap.Pop();

                if (isRefVar && nthterm->ObjectClass()!=STRING) {
                    _FString * type = (_FString*)nthterm->Type();
                    _String errMsg = _String ("User-defined function '")
                                     &*(_String*)batchLanguageFunctionNames(-numberOfTerms-1)
                                     &"' expected a string for the reference variable '"
                                     & *argNameString
                                     &"' but was passed a " & *type->theString
                                     & " with the value of " & _String((_String*)nthterm->toStr());

                    DeleteObject (type);
                    WarnError    (errMsg);
                    return false;
                }

                if (f<0) { // not an existing var
                    _Variable newV (*argNameString);
                    f =  LocateVarByName (*argNameString);
                }
                _Variable* theV = FetchVar (f);

                if (!isRefVar) {
                    if (theV->IsIndependent())
                        // if the variable exists and is independent then
                        // simply swap the value of the var, otherwise
                        // duplicate the entire variable
                    {
                        theV->varFlags &= HY_VARIABLE_SET;
                        if (!theV->varValue) {
                            theV->Compute();
                        }
                        displacedValues<<theV->varValue;
                        theV->varFlags |= HY_VARIABLE_CHANGED;
                        /*if (*theV->GetName() == _String("gbdd1"))
                        {
                             printf ("Setting argument value of %s to %s\n",
                             theV->GetName()->sData,
                             _String((_String*)nthterm->toStr()).sData
                         );
                        }*/

                        theV->varValue = nthterm;
                        existingIVars<<theV->GetAVariable();
                    } else {
                        _Variable newV (*argNameString);
                        newV.SetValue(nthterm);
                        existingDVars<<theV->GetAVariable();
                        displacedVars<<theV;
                        variablePtrs.Replace (theV->GetAVariable(),(_PMathObj)newV.makeDynamic());
                        DeleteObject (nthterm);
                    }
                } else {
                    referenceArgs<< variableNames.Retrieve(f);
                    displacedReferences<<theV->GetAVariable();

                    _String * refArgName = ((_FString*)nthterm)->theString;

                    if (nameSpace) {
                        *refArgName = AppendContainerName (*refArgName, nameSpace);
                    }

                    i = LocateVarByName (*refArgName);

                    if (i<0) {
                        _Variable newAV (*refArgName);
                        i =  LocateVarByName (*refArgName);
                    }
                    variableNames.SetXtra (f, variableNames.GetXtra (i));
                    *argNameString = *argNameString&'&';
                    DeleteObject (nthterm);
                }
            }

            if (purgeFlas) {
                ((_ExecutionList*)batchLanguageFunctions(opCode))->ResetFormulae();
            }

            _ExecutionList * functionBody = ((_ExecutionList*)batchLanguageFunctions(opCode));
            if (currentExecutionList && currentExecutionList->stdinRedirect) {
                functionBody -> stdinRedirect    = currentExecutionList->stdinRedirect;
                functionBody -> stdinRedirectAux = currentExecutionList->stdinRedirectAux;
            }
            _PMathObj ret = functionBody->Execute();

            functionBody -> stdinRedirect    = nil;
            functionBody -> stdinRedirectAux = nil;

            if (terminateExecution) {
                theScrap.Push (new _Constant (0.0));
                return true;
            }

            if (ret) {
                theScrap.Push (ret);
            }

            for (long di = 0; di < referenceArgs.lLength; di++) {
                variableNames.SetXtra(LocateVarByName (*(_String*)referenceArgs(di)),displacedReferences.lData[di]);
            }

            //for (i= referenceArgs.lLength-1; i>=0; i--)
            //  variableNames.SetXtra(LocateVarByName (*(_String*)referenceArgs(i)),displacedReferences.lData[i]);

            for (long dv = 0; dv < displacedVars.lLength; dv++) {
                variablePtrs.Replace (existingDVars.lData[dv],(_PMathObj)displacedVars(dv));
            }

            //for (i= displacedVars.lLength-1; i>=0; i--)
            //  variablePtrs.Replace (existingDVars(i),(_PMathObj)displacedVars(i));

            for (long dv2 = 0; dv2 < displacedValues.lLength; dv2++) {
                _Variable* theV = LocateVar (existingIVars.lData[dv2]);
                DeleteObject(theV->varValue);
                theV->varValue = ((_PMathObj)displacedValues(dv2));
            }

            //for (i= displacedValues.lLength-1; i>=0; i--)
            //{
            //  _Variable* theV = LocateVar (existingIVars(i));
            //  DeleteObject(theV->varValue);
            //  theV->varValue = ((_PMathObj)displacedValues(i));
            //}


        }
        return true;

    }

    if (theScrap.theStack.lLength<numberOfTerms) {
        _String *s,errMsg;
        s = (_String*)toStr();
        errMsg = *s;
        DeleteObject(s);
        errMsg = errMsg&
                 " needs "&_String(numberOfTerms)& " arguments. Only "&_String(theScrap.StackDepth())&" were given.";
        WarnError (errMsg);
        return nil;
    }

    _PMathObj term1, term2 = nil, term3 = nil, temp;


    if (numberOfTerms >= 3) {
        long sL = theScrap.theStack.lLength-1;
        term3 = (_PMathObj)theScrap.theStack.lData[sL--];
        term2 = (_PMathObj)theScrap.theStack.lData[sL--];
        term1 = (_PMathObj)theScrap.theStack.lData[sL];
        theScrap.theStack.lLength = sL;
        temp = term1->Execute (opCode, term2, term3);
        DeleteObject (term1);
        DeleteObject (term2);
        DeleteObject (term3);
    } else if (numberOfTerms == 2) {
        long sL = theScrap.theStack.lLength-1;
        term2 = (_PMathObj)theScrap.theStack.lData[sL--];
        term1 = (_PMathObj)theScrap.theStack.lData[sL];
        theScrap.theStack.lLength = sL;
        temp = term1->Execute (opCode, term2, nil);
        DeleteObject (term1);
        DeleteObject (term2);
    } else {
        term1 = (_PMathObj)theScrap.theStack.lData[--theScrap.theStack.lLength];
        temp = term1->Execute (opCode, nil, nil);
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
        depth -= batchLanguageFunctionParameters(-numberOfTerms-1) - 1;
        return;
    }
    depth -= numberOfTerms - 1;
}


//__________________________________________________________________________________

bool        _Operation::ExecutePolynomial (_Stack& theScrap)
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
        WarnError (s&
                   " needs "&_String(numberOfTerms)& " arguments. Only "&_String(theScrap.StackDepth())&" were given.");
        return    false;
    }

    _PMathObj term1,
              term2 = nil,
              temp;

    bool      opResult = true;

    if (numberOfTerms == 2) {
        term2 = theScrap.Pop();
    }

    term1 = theScrap.Pop();
    temp  = term1->Execute (opCode, term2, nil);
    DeleteObject (term1);

    if (temp) {
        theScrap.Push (temp);
        DeleteObject (temp);
    } else {
        opResult = false;
    }

    if (term2) {
        DeleteObject (term2);
    }

    return opResult;
}

