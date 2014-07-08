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

#include "likefunc.h"
#include "parser.h"
#include "polynoml.h"
#include "batchlan.h"

extern _SimpleList freeSlots;
extern _SimpleList deferIsConstant;

//__________________________________________________________________________________

_Variable::_Variable (void)
{
    varFormula = nil;
    varFlags   = HY_VARIABLE_NOTSET;
    theName    = nil;
    varValue   = nil;
    theIndex   = -1;
    SetBounds (DEFAULTLOWERBOUND, DEFAULTUPPERBOUND);
}

//__________________________________________________________________________________
void _Variable::Initialize (void)
{
    //_Formula::Initialize();
    _Constant::Initialize();
    theName = (_String*)checkPointer(new _String());
    varValue = nil;
    theIndex = -1;
    varFlags = HY_VARIABLE_NOTSET;
    SetBounds (DEFAULTLOWERBOUND, DEFAULTUPPERBOUND);
}

//__________________________________________________________________________________
void _Variable::Duplicate (BaseRef r)
{
    _Variable *v = (_Variable*)r;
    //theFormula.Duplicate (&(v->theFormula));
    if (v->varFormula) {
        varFormula = new _Formula();
        varFormula->theFormula.Duplicate (&v->varFormula->theFormula);
    } else {
        varFormula = nil;
    }

    //theStack.theStack.Duplicate (&(v->theStack.theStack));
    //theStack.theStack.Clear();
    theValue = v->theValue;
    varValue = v->varValue;
    if (varValue) {
        varValue->nInstances++;
    }
    theIndex = v->theIndex;
    theName = v->theName;
    theName->nInstances++;
    lowerBound = v->lowerBound;
    upperBound = v->upperBound;
    //hasBeenChanged = v->hasBeenChanged;
    varFlags = v->varFlags;
}

//__________________________________________________________________________________
BaseRef _Variable::makeDynamic (void)
{
    _Variable * res = new _Variable;
    if (!res) {
        isError(0);
        return nil;
    }
    //memcpy ((char*)res, (char*)this, sizeof (_Variable));
    res->Duplicate(this);
    return res;
}

//__________________________________________________________________________________
bool _Variable::CheckFForDependence (long idx, bool opt)
{
    if (varFormula) {
        return varFormula->CheckFForDependence (idx, opt);
    }

    return false;
}

//__________________________________________________________________________________
BaseRef _Variable::toStr(void)
{
    if (varValue&&varValue->IsPrintable()) {
        return varValue->toStr();
    }
    _PMathObj vv = Compute();
    if (!vv) {
        return new _String("NAN");
    }
    return new _String((_String*)vv->toStr());
}

//__________________________________________________________________________________
void _Variable::toFileStr(FILE* f)
{
    if (varValue&&varValue->IsPrintable()) {
        varValue->toFileStr(f);
    } else {
        _PMathObj vv = Compute();
        if (!vv) {
            fprintf(f,"NAN");
        } else {
            vv->toFileStr(f);
        }
    }

}
//__________________________________________________________________________________

_Variable::_Variable (_String&s, bool isG)
{
    theName         = (_String*)checkPointer(new _String(s));
    varFlags        = HY_VARIABLE_NOTSET|(isG?HY_VARIABLE_GLOBAL:0);
    varValue        = nil;
    varFormula      = nil;
    SetBounds       (DEFAULTLOWERBOUND, DEFAULTUPPERBOUND);
    InsertVar       (this);
}

//__________________________________________________________________________________

_Variable::_Variable (_String&s, _String&f, bool isG)//:  _Formula (f)
{
    //hasBeenChanged = false;
    //isGlobal = isG;
    theName     = (_String*)checkPointer(new _String(s));
    varFlags    = isG?HY_VARIABLE_GLOBAL:0;
    varValue    = nil;
    SetBounds   (DEFAULTLOWERBOUND, DEFAULTUPPERBOUND);
    InsertVar   (this);
    varFormula = new _Formula (f);
    if (varFormula->IsAConstant()) {
        _PMathObj theP = varFormula->Compute();
        if (theP) {
            SetValue (theP);
            delete   (varFormula);
            varFormula = nil;
        } else {
            return;
        }
    }
}



//__________________________________________________________________________________

_Variable::~_Variable (void)
{
    nInstances++;
    if (varValue) {
        DeleteObject (varValue);
    }
    if (theName) {
        DeleteObject (theName);
    }
    if (varFormula) {
        delete (varFormula);
    }
}

//__________________________________________________________________________________
bool    _Variable::IsVariable (void)
{
    return true;
}

//__________________________________________________________________________________

void        _Variable::ScanForVariables (_AVLList& l, bool globals, _AVLListX* tagger, long weight) {
    if (varValue) {
        varValue->ScanForVariables (l, globals,tagger, weight);
    }
    if (varFormula && varFormula->theFormula.lLength) {
        varFormula->ScanFForVariables(l,globals, false, true, false,tagger, weight);
    }
}


//__________________________________________________________________________________

_PMathObj  _Variable::Compute (void) // compute or return the value
{
    if (varFormula == nil) { // no formula, just return the value
        if (varValue) {
            return varValue->Compute();
        }

        if (varFlags & HY_VARIABLE_NOTSET) {
            ReportWarning (_String ("Variable '") & *GetName() & "' was not initialized prior to being used");
        }

        varValue =  new _Constant(theValue);
    } else {
        //printf ("Recomputing value of %s\n", theName->sData);
        if (useGlobalUpdateFlag) {
            if ((varFlags & HY_DEP_V_COMPUTED) && varValue) {
                return varValue;
            } else if (varFormula->HasChanged()||!varValue) {
                DeleteObject (varValue);
                varValue = (_PMathObj)varFormula->Compute()->makeDynamic();
            }

            varFlags |= HY_DEP_V_COMPUTED;
        } else if (varFormula->HasChanged()||!varValue) {
            DeleteObject (varValue);
            varValue = (_PMathObj)varFormula->Compute()->makeDynamic();
        }

    }

    return varValue;
}

//__________________________________________________________________________________

void  _Variable::CompileListOfDependents (_SimpleList& rec)
{
    _SimpleList tcache;
    long        iv,
                i = variableNames.Traverser (tcache,iv,variableNames.GetRoot());

    for (; i >= 0; i = variableNames.Traverser (tcache,iv)) {
        _Variable* thisVar = FetchVar (i);
        if (!thisVar->IsIndependent()) {
            if (thisVar->CheckFForDependence (theIndex)) {
                long f = thisVar->GetAVariable();
                if (rec.Find(f)<0) {
                    rec<<f;
                }
            }
        }
    }
}

//__________________________________________________________________________________
void  _Variable::SetValue (_PMathObj theP, bool dup) // set the value of the var
{
    //hasBeenChanged = true;
    varFlags &= HY_VARIABLE_SET;
    varFlags |= HY_VARIABLE_CHANGED;

    long     valueClass = theP->ObjectClass();

    if (valueClass==NUMBER) {
        if (varFormula) {

            // also update the fact that this variable is no longer dependent in all declared
            // variable containers which contain references to this variable
            for (unsigned long i = 0UL; i<variablePtrs.lLength; i++) {
                if (freeSlots.Find(i)>=0) {
                    continue;
                }
                _Variable* theV = (_Variable*)variablePtrs(i);
                if (theV->IsContainer()) {
                    _VariableContainer* theVC = (_VariableContainer*)theV;
                    if (!theVC->RemoveDependance (theIndex)) {
                        ReportWarning ((_String("Can't make variable ")&*GetName()&" independent in the context of "&*theVC->GetName()&" because its template variable is not independent."));
                        continue;
                    }
                }
            }
            for (unsigned long i = 0UL; i<likeFuncList.lLength; i++)
                if (((_String*)likeFuncNamesList(i))->sLength) {
                    ((_LikelihoodFunction*)likeFuncList(i))->UpdateDependent(theIndex);
                }

            //_Formula::Clear();
            delete (varFormula);
            varFormula = nil;
        }
        if (varValue) {
            DeleteObject (varValue);
            varValue=nil;
        }

        theValue = theP->Value();

        if (!dup) {
            DeleteObject (theP);
        }

        if (theValue<lowerBound || theValue>upperBound) {
            if (theValue <= lowerBound+1e-50) {
                theValue = lowerBound;
            } else {
                theValue = upperBound;
            }
        }
    } else {
        if (varFormula) {
            delete (varFormula);
            varFormula = nil;
            //theFormula.Clear();
        }
        if (varValue) {
            DeleteObject (varValue);
            varValue=nil;
        }
        if (valueClass==TREE) {
            variablePtrs.lData[theIndex] = (long)(((_TheTree*)theP)->makeDynamicCopy(GetName()));
            DeleteObject(this);
        } else {
            if (dup) {
                varValue = (_PMathObj)theP->makeDynamic();
            } else {
                varValue = theP;
            }
        }
    }
}

//__________________________________________________________________________________
void  _Variable::SetNumericValue (_Parameter v) // set the value of the var to a number
{
    //hasBeenChanged = true;
    varFlags &= HY_VARIABLE_SET;
    varFlags |= HY_VARIABLE_CHANGED;
    theValue = v;

    if (theValue<lowerBound || theValue>upperBound) {
        if (theValue<=lowerBound+1e-50) {
            theValue = lowerBound;
        } else {
            theValue = upperBound;
        }
    }
}

//__________________________________________________________________________________

void  _Variable::CheckAndSet (_Parameter c, bool oob) // set the value of the var
{
    //hasBeenChanged = true;
    varFlags &= HY_VARIABLE_SET;
    varFlags |= HY_VARIABLE_CHANGED;
    _Parameter l = lowerBound+1.0e-30,
               u = upperBound-1.0e-30;
    if (c<l || c>u ) {
        if (oob) {
            return;
        }

        if (c<l) {
            theValue = l;
        } else {
            theValue = u;
        }
    } else {
        theValue = c;
    }

    if (varValue) {
        DeleteObject (varValue);
    }

    varValue =  new _Constant(theValue);
}

//__________________________________________________________________________________
void    _Variable::SetBounds (_Parameter lb, _Parameter ub)
{
    lowerBound = lb;
    upperBound = ub;

    /*_String * myName = GetName();
    if (myName)
        ReportWarning (_String ("Set variable bounds for '") & *myName & "' to [" & lb & ',' & ub & "].");
     */
}

//__________________________________________________________________________________
void    _Variable::EnsureTheValueIsInBounds (void)
{
    if (ObjectClass () == NUMBER && IsIndependent()) {
        _Constant*   myValue = (_Constant*) Compute();
        if (myValue->Value() < lowerBound) {
            SetValue (new _Constant (lowerBound),false);
        } else if (myValue->Value() > upperBound) {
            SetValue (new _Constant (upperBound),false);
        }
    }
}


//__________________________________________________________________________________
void    _Variable::ClearConstraints (void)
{
    if (IsCategory ()) {
        _Variable newVar (*GetName(), IsGlobal());
        newVar.SetValue ((_PMathObj)Compute()->makeDynamic(),false);
        ReplaceVar ( &newVar);
        /*_Matrix * modelMatrix = (_Matrix*)LocateVar(modelMatrixIndices.lData[1])->GetValue();
        for (long k=0; k<4; k++)
            for (long k2 = 0; k2<4; k2++)
                if (k!=k2)
                {
                    StringToConsole (*(_String*)modelMatrix->GetFormula(k,k2)->toStr());
                    BufferToConsole ("\n");
                }
        */
    } else {
        if (!IsIndependent()) {
            SetValue ((_PMathObj)Compute()->makeDynamic(),false);
        }
        SetBounds (DEFAULTLOWERBOUND,DEFAULTUPPERBOUND);
    }
}

//__________________________________________________________________________________
bool _Variable::IsConstant (void)
{
    if (varFormula && varFormula->theFormula.lLength) {
        return varFormula->IsConstant();
    }

    if (varValue && varValue->ObjectClass () != NUMBER) {
        return varValue->IsConstant();
    }

    return false;
}

//__________________________________________________________________________________

void  _Variable::SetFormula (_Formula& theF) // set the value of the var to a formula
{
    bool changeMe    = false,
         isAConstant = theF.IsAConstant();

    _Formula* myF = &theF;

    if (isAConstant) {
        _PMathObj theP = theF.Compute();
        if (theP) {
            myF = new _Formula ((_PMathObj)theP->makeDynamic(),false);
            checkPointer (myF);
        } else {
            return;
        }
    }

    _SimpleList vars;
    {
        _AVLList vA (&vars);
        theF.ScanFForVariables (vA,true);
        vA.ReorderList();
    }

    if (vars.BinaryFind(theIndex)>=0) {
        _String * sf = (_String*)theF.toStr();
        WarnError ((_String("Can't set variable ")&*GetName()&" to "&*sf&" because it would create a circular dependance."));
        DeleteObject(sf);
        if (&theF!=myF) {
            delete myF;
        }
        return;
    }

    varFlags &= HY_VARIABLE_SET;

    if (varFlags & HY_VARIABLE_CHANGED) {
        varFlags -= HY_VARIABLE_CHANGED;
    }


    if (varFormula) {
        delete (varFormula);
        varFormula = nil;
    } else {
        changeMe = true;
    }

    if (varValue) {
        DeleteObject (varValue);
        varValue=nil;
    }

    //_Formula::Duplicate ((BaseRef)myF);
    varFormula = new _Formula;
    varFormula->Duplicate ((BaseRef)myF);

    // mod 20060125 added a call to simplify constants
    varFormula->SimplifyConstants ();

    // also update the fact that this variable is no longer independent in all declared
    // variable containers which contain references to this variable
    if (changeMe)
        if (deferSetFormula) {
            *deferSetFormula << theIndex;
            deferIsConstant  << isAConstant;
        } else {
            long i;
            _SimpleList tcache;
            long        iv;

            i = variableNames.Traverser (tcache,iv,variableNames.GetRoot());

            for (; i >= 0; i = variableNames.Traverser (tcache,iv)) {
                _Variable* theV = FetchVar(i);
                if (theV->IsContainer()) {
                    _VariableContainer* theVC = (_VariableContainer*)theV;
                    if (theVC->SetDependance(theIndex) == -2) {
                        ReportWarning ((_String("Can't make variable ")&*GetName()&" dependent in the context of "&*theVC->GetName()&" because its template variable is bound by another relation in the global context."));
                        continue;
                    }
                }
            }
            {
                for (long i = 0; i<likeFuncList.lLength; i++)
                    if (((_String*)likeFuncNamesList(i))->sLength) {
                        ((_LikelihoodFunction*)likeFuncList(i))->UpdateIndependent(theIndex,isAConstant);
                    }
            }
        }

    if (&theF!=myF) {
        delete myF;
    }
}

//__________________________________________________________________________________
void  _Variable::PreMarkChanged (void)
{
    if (varFormula) {
        varFlags &= HY_DEP_V_INSPECTED_CLR;

        if (HasChanged(false)) {
            varFlags |= HY_DEP_V_MODIFIED;
        }
        if (HasChanged(true)) {
            varFlags |= HY_DEP_V_MODIFIED_CATS;
        }

        varFlags |= HY_DEP_V_INSPECTED;
    }
}

//__________________________________________________________________________________
void  _Variable::PostMarkChanged (void)
{
    varFlags &= HY_DEP_CLEAR_MASK;
}


//__________________________________________________________________________________
bool  _Variable::HasChanged (bool ignoreCats) // does this variable need recomputing
{
    if (varFormula) {
        if (useGlobalUpdateFlag && (varFlags&HY_DEP_V_COMPUTED)) {
            return false;
        }

        if (varFlags&HY_DEP_V_INSPECTED) {
            return ignoreCats?(varFlags&HY_DEP_V_MODIFIED_CATS):(varFlags&HY_DEP_V_MODIFIED);
        }

        return  varFormula->HasChanged(ignoreCats);

    } else {
        if (varValue&&(varValue->IsVariable())) {
            return varValue->HasChanged();
        }
        if (ignoreCats && IsCategory()) {
            return false;
        }
        return varFlags & HY_VARIABLE_CHANGED;
    }

}

//__________________________________________________________________________________

void _Variable::MarkDone (void)
{
    if (!varFormula && (varFlags & HY_VARIABLE_CHANGED) && !(varValue && varValue->IsVariable())) {
        varFlags -= HY_VARIABLE_CHANGED;
    }
}

//__________________________________________________________________________________
_PMathObj    _Variable::ComputeReference (_PMathObj context)
{
    _String reference_string (*GetName());
    reference_string = AppendContainerName(reference_string, (_VariableContainer*)context);
    
    return new _FString (reference_string, false);
}

//__________________________________________________________________________________
_String    _Variable::ContextFreeName(void) {
    long location = theName->FindBackwards (".", 0, -1);
    if (location > 0) {
       return theName->Cut (location+1,-1); 
    }  
    return *theName;
}

//__________________________________________________________________________________
_String    _Variable::ParentObjectName(void) {
    long location = theName->FindBackwards (".", 0, -1);
    if (location > 0) {
       return theName->Cut (0,location-1); 
    }  
    return empty;
}

//__________________________________________________________________________________
long    DereferenceString (_PMathObj v, _PMathObj context, char reference_type){
    if (v && v->ObjectClass () == STRING) {
        _FString * value = (_FString*)v;
        _String referencedVariable = *value->theString;
        if (reference_type == HY_STRING_LOCAL_DEREFERENCE && context) {
            referencedVariable = AppendContainerName(referencedVariable, (_VariableContainer*)context);
        }
        return LocateVarByName(referencedVariable);
    }
    return -1;
}

//__________________________________________________________________________________
long    DereferenceVariable (long index, _PMathObj context, char reference_type){
    if (reference_type == HY_STRING_DIRECT_REFERENCE) {
        return index;
    }
    
    return  DereferenceString (FetchObjectFromVariableByTypeIndex(index, STRING), context, reference_type);
}

