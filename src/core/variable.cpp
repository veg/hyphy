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

#include "likefunc.h"
#include "parser.h"
#include "polynoml.h"
#include "batchlan.h"

#include "function_templates.h"
#include "global_things.h"

using namespace hy_global;

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
void _Variable::Initialize (bool)
{
    //_Formula::Initialize();
    _Constant::Initialize();
    theName = new _String();
    varValue = nil;
    theIndex = -1;
    varFlags = HY_VARIABLE_NOTSET;
    SetBounds (DEFAULTLOWERBOUND, DEFAULTUPPERBOUND);
}

//__________________________________________________________________________________
void _Variable::Duplicate (BaseRefConst r) {
    _Variable const *v = (_Variable const*)r;
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
        varValue->AddAReference();
    }
    theIndex = v->theIndex;
    theName = v->theName;
    theName->AddAReference();
    lowerBound = v->lowerBound;
    upperBound = v->upperBound;
    //hasBeenChanged = v->hasBeenChanged;
    varFlags = v->varFlags;
}

//__________________________________________________________________________________
BaseRef _Variable::makeDynamic (void) const{
    _Variable * res = new _Variable;
    res->Duplicate(this);
    return res;
}

//__________________________________________________________________________________
void * _Variable::operator new (size_t size) {
     return MemAllocate (size);
}

//__________________________________________________________________________________
void  _Variable::operator delete (void * p) {
    free (p);
}

//__________________________________________________________________________________

unsigned long        _Variable::ObjectClass (void) const {
    
    if (varValue) {
        return varValue->ObjectClass();
    }
    
    if (varFormula && !varFormula->IsEmpty()) {
        return varFormula->ObjectClass();
    }
    
    return NUMBER;
}

//__________________________________________________________________________________
bool _Variable::CheckFForDependence (long idx, bool opt) {
    if (varFormula) {
        return varFormula->CheckFForDependence (idx, opt);
    }
    return false;
}

//__________________________________________________________________________________
bool _Variable::CheckFForDependence (_AVLList const& indices, bool opt) {
    if (varFormula) {
        return varFormula->CheckFForDependence (indices, opt);
    }
    return false;
}

//__________________________________________________________________________________
BaseRef _Variable::toStr(unsigned long padding)
{
    if (varValue&&varValue->IsPrintable()) {
        return varValue->toStr(padding);
    }
    HBLObjectRef vv = Compute();
    if (!vv) {
        return new _String("NAN");
    }
    return new _String((_String*)vv->toStr(padding));
}

//__________________________________________________________________________________
void _Variable::toFileStr(FILE* f, unsigned long padding)
{
    if (varValue&&varValue->IsPrintable()) {
        varValue->toFileStr(f, padding);
    } else {
        HBLObjectRef vv = Compute();
        if (!vv) {
            fprintf(f,"NAN");
        } else {
            vv->toFileStr(f, padding);
        }
    }

}
//__________________________________________________________________________________

_Variable::_Variable (_String const&s, bool isG) {
    theName         = new _String(s);
    varFlags        = HY_VARIABLE_NOTSET|(isG?HY_VARIABLE_GLOBAL:0);
    varValue        = nil;
    varFormula      = nil;
    SetBounds       (DEFAULTLOWERBOUND, DEFAULTUPPERBOUND);
    InsertVar       (this);
}

//__________________________________________________________________________________

_Variable::_Variable (_String const&s, _String const &f, bool isG) {
    //hasBeenChanged = false;
    //isGlobal = isG;
    theName     = new _String(s);
    varFlags    = isG?HY_VARIABLE_GLOBAL:0;
    varValue    = nil;
    SetBounds   (DEFAULTLOWERBOUND, DEFAULTUPPERBOUND);
    InsertVar   (this);
    varFormula = new _Formula (f);
    if (varFormula->IsAConstant()) {
        HBLObjectRef theP = varFormula->Compute();
        if (theP) {
            SetValue (theP,true,true,NULL);
            delete   (varFormula);
            varFormula = nil;
        } else {
            return;
        }
    }
}



//__________________________________________________________________________________

_Variable::~_Variable (void) {
    DeleteObject (varValue);
    DeleteObject (theName);
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

void        _Variable::ScanForVariables (_AVLList& l, bool globals, _AVLListX* tagger, long weight) const {
    if (varValue) {
        varValue->ScanForVariables (l, globals,tagger, weight);
    }
    if (varFormula && varFormula->theFormula.lLength) {
        varFormula->ScanFForVariables(l,globals, false, true, false,tagger, weight);
    }
}
  // long call_count = 0L;

  //__________________________________________________________________________________

HBLObjectRef  _Variable::ComputeMatchingType(long type) {
  HBLObjectRef computed_value = Compute();
  if (computed_value && (computed_value->ObjectClass() & type) > 0L) {
    return computed_value;
  }
  return nil;
}


//__________________________________________________________________________________

HBLObjectRef  _Variable::Compute (void) {
    // compute or return the value
    // call_count++;
    
    auto update_var_value = [this] () -> void {
        if (!varValue || varFormula->HasChanged()) {
            HBLObjectRef new_value = (HBLObjectRef)varFormula->Compute()->makeDynamic();
            DeleteObject (varValue);
            varValue = new_value;
            //DeleteObject (varValue);
            //(varValue = varFormula->Compute())->AddAReference();
        }
    };
  
    if (varFlags & HY_VARIABLE_COMPUTING) {
      HandleApplicationError (_String ("A recursive dependency error in _Variable::Compute; this is an HBL implementation bug; offending variable is '") & *GetName() & "'");
      return new _MathObject;
    }
  
    varFlags |= HY_VARIABLE_COMPUTING;
  
    if (varFormula == nil) { // no formula, just return the value
        if (varValue) {
            varFlags &= HY_VARIABLE_COMPUTING_CLR;
            return varValue->Compute();
        }

        if (varFlags & HY_VARIABLE_NOTSET) {
            ReportWarning (_String ("Variable '") & *GetName() & "' was not initialized prior to being used");
        }

        varValue =  new _Constant(theValue);
        
        //printf ("Recomputing value of %s => %g\n", theName->sData, varValue->Value());

    } else {
        if (useGlobalUpdateFlag) {
            if ((varFlags & HY_DEP_V_COMPUTED) && varValue) {
                varFlags &= HY_VARIABLE_COMPUTING_CLR;
                return varValue;
            } else {
                update_var_value ();
            }
                
            varFlags |= HY_DEP_V_COMPUTED;
          
        } else {
            update_var_value ();
        }

    }

    varFlags &= HY_VARIABLE_COMPUTING_CLR;
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
                long f = thisVar->get_index();
                if (rec.Find(f)<0) {
                    rec<<f;
                }
            }
        }
    }
}

//__________________________________________________________________________________
void  _Variable::SetValue (hyFloat new_value) {
// set the value of the var
  this->SetValue (new _Constant (new_value), false,true,NULL);
}

//__________________________________________________________________________________
void        SetVariablesToOwnValues (_AVLList const & indices) {
    /* compute current variable values*/
    _SimpleList  hist;
    long         ls = -1L, cn;
    
    cn = indices.Traverser (hist,ls,indices.GetRoot());
    
    _SimpleList _index2value;
    _AVLList    index2value (&_index2value);
    
    while (cn>=0) {
        _Variable* this_var = LocateVar (indices.RetrieveLong(cn));
        if (this_var->IsIndependent()) {
            // do nothing
        } else {
            HBLObjectRef curValue = this_var->Compute();
            curValue->AddAReference();
            index2value.InsertNumber(this_var->get_index());
            this_var->SetValue(curValue, false, false, NULL);
            //printf ("SetVariablesToOwnValues : Set own value for %s\n", this_var->GetName()->get_str());
        }
        cn = indices.Traverser (hist,ls);
    }
    
    // now clean up variable containers that may have had these as dependent template variables
    
    DoForEachVariable([&index2value] (_Variable* v, long v_idx) -> void {
        if (v->IsContainer()) {
            ((_VariableContainer*)v)->RemoveDependance (index2value);
        }
    });
    
    DoForEachLikelihoodFunction ([&index2value] (_LikelihoodFunction *lf, long idx) -> void {
        lf->UpdateDependent(index2value);
    });
    
    // update likelihood functions that may have had these as dependent variables
    
}

//__________________________________________________________________________________
void  _Variable::SetValue (HBLObjectRef theP, bool dup, bool do_checks, _AVLList* keep_track_of_changes) { // set the value of the var
    //hasBeenChanged = true;
    if (varFlags & HY_VARIABLE_COMPUTING) {
        HandleApplicationError (_String ("A recursive dependency error in _Variable::SetValue; this is an HBL implementation bug; offending variable is ") & GetName()->Enquote());
        return ;
    }
  
    varFlags &= HY_VARIABLE_SET;
    varFlags |= HY_VARIABLE_CHANGED;
    
    if (keep_track_of_changes) {
        keep_track_of_changes->InsertNumber (get_index());
    }

    long     valueClass = theP->ObjectClass();
  
    if (valueClass==NUMBER) {
      
        if (varFormula) {

            if (do_checks) {
                // also update the fact that this variable is no longer dependent in all declared
                // variable containers which contain references to this variable
                
                DoForEachVariable([this] (_Variable* v, long v_idx) -> void {
                    if (v->IsContainer()) {
                        if (!((_VariableContainer*)v)->RemoveDependance (this->theIndex)) {
                            ReportWarning ((_String("Can't make variable ")&GetName()->Enquote()&" independent in the context of "&*v->GetName()&" because its template variable is not independent."));
                        }
                    }
                });
                
                for (unsigned long i = 0UL; i<likeFuncList.lLength; i++)
                    if (((_String*)likeFuncNamesList(i))->nonempty()) {
                        ((_LikelihoodFunction*)likeFuncList(i))->UpdateDependent(theIndex);
                    }
            }

            //_Formula::Clear();
            delete varFormula;
            varFormula = nil;
        }
        
        theValue = theP->Value();
        
        if (theValue<lowerBound || theValue>upperBound) {
            /*if (verbosity_level >= 100) {
                printf ("\n<========== OUT OF BOUNDS FOR VARIABLE %s : value %.16g, range [%.16g, %.16g] =========>\n" , theName->get_str(), theValue, lowerBound, upperBound);
            }*/
            if (theValue <= lowerBound+1e-50) {
                theValue = lowerBound;
            } else {
                theValue = upperBound;
            }
        }
        
        if (varValue) {
            if (varValue->ObjectClass() == NUMBER && varValue->SingleReference()) {
                ((_Constant*)varValue)->SetValue(theValue);
            } else {
                DeleteAndZeroObject (varValue);
            }
        }


        if (!dup) {
            DeleteObject (theP);
        }

        
    } else {
        if (varFormula) {
            /*if (doPrint) {
              printf ("Resetting constraint on %s which was %s\n", theName->sData,_String((_String*)varFormula->toStr()).sData);
            }*/
            delete (varFormula);
            varFormula = nil;
            //theFormula.Clear();
        }
        
        
        if (varValue) {
            if (dup && varValue->SingleReference() && varValue->ObjectClass() == theP->ObjectClass()) {
                switch (varValue->ObjectClass()) {
                    case STRING: {
                        _FString *existing   = (_FString*)varValue,
                                 *new_string = ((_FString*)theP);
                        
                        if (new_string->empty()) {
                            existing->SetData (kEmptyString);
                        } else {
                            existing->SetData (((_FString*)theP)->get_str());
                        }
                        return;
                    }
                        
                }
            }
                        
            DeleteObject (varValue);
            varValue=nil;
        }
      
        
        if (dup) {
            varValue = (HBLObjectRef)theP->makeDynamic();
        } else {
            varValue = theP;
        }
    }
}

//__________________________________________________________________________________
void  _Variable::SetNumericValue (hyFloat v) // set the value of the var to a number
{
    //hasBeenChanged = true;
    if (varFlags & HY_VARIABLE_COMPUTING) {
      HandleApplicationError (_String ("A recursive dependency error in _Variable::SetNumericValue; this is an HBL implementation bug; offending variable is '") & *GetName() & "'");
      return;
    }

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

void  _Variable::CheckAndSet (hyFloat c, bool oob, _AVLList * keep_track_of_changes) // set the value of the var
{
    //hasBeenChanged = true;
    hyFloat l = lowerBound+1.0e-30,
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
    
    varFlags &= HY_VARIABLE_SET;
    varFlags |= HY_VARIABLE_CHANGED;

    if (keep_track_of_changes) {
        keep_track_of_changes->InsertNumber (get_index());
    }

    if (varValue) {
        if (varValue->ObjectClass() == NUMBER && varValue->SingleReference()) {
            ((_Constant*)varValue)->theValue = theValue;
            return;
        }
        DeleteObject (varValue);
    }
    
    varValue =  new _Constant(theValue);
    
}

//__________________________________________________________________________________
void    _Variable::SetBounds (hyFloat lb, hyFloat ub)
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
            SetValue (new _Constant (lowerBound),false,true,NULL);
        } else if (myValue->Value() > upperBound) {
            SetValue (new _Constant (upperBound),false,true,NULL);
        }
    }
}


//__________________________________________________________________________________
void    _Variable::ClearConstraints (void)
{
    if (IsCategory ()) {
        _Variable newVar (*GetName(), IsGlobal());
        newVar.SetValue ((HBLObjectRef)Compute()->makeDynamic(),false,true,NULL);
        ReplaceVar ( &newVar);
        /*_Matrix * modelMatrix = (_Matrix*)LocateVar(modelMatrixIndices.list_data[1])->GetValue();
        for (long k=0; k<4; k++)
            for (long k2 = 0; k2<4; k2++)
                if (k!=k2)
                {
                    StringToConsole (*(_String*)modelMatrix->GetFormula(k,k2)->toStr());
                    BufferToConsole ("\n");
                }
        */
    } else {
        //printf ("ClearConstraints %s %x\n", GetName()->get_str(), varFormula);

        if (!IsIndependent()) {
            SetValue ((HBLObjectRef)Compute()->makeDynamic(),false,true,NULL);
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

    if (varValue) {
        long varType = varValue->ObjectClass ();
        if (!(varType & HY_MUTABLE_OBJECT)) {
            return varValue->IsConstant();
        }
    }

    return false;
}

  //long wp_count = 0;

//__________________________________________________________________________________

void  _Variable::SetFormula (_Formula& theF) {
//  bind the variable to an expression
    if (varFlags & HY_VARIABLE_COMPUTING) {
      HandleApplicationError (_String ("A recursive dependency error in _Variable::SetFormula; this is an HBL implementation bug; offending variable name is '") & *GetName() & "'");
      return ;
    }

    bool changeMe    = false,
         isAConstant = theF.IsAConstant();

    //if (doPrint) {
    //printf ("Constraining %s to %s\n", theName->sData, _String((_String*)theF.toStr()).sData);
    //}

    _Formula* right_hand_side = &theF;

    if (isAConstant) {
        HBLObjectRef theP = theF.Compute();
        if (theP) {
            right_hand_side = new _Formula ((HBLObjectRef)theP->makeDynamic(),false);
        } else {
            return;
        }
    }

    _SimpleList vars;
    _AVLList vA (&vars);
    theF.ScanFForVariables (vA,true);

    if (vA.Find((BaseRefConst)theIndex)>=0) {
        HandleApplicationError ((_String("Can't set variable ")&*GetName()&" to "&*((_String*)theF.toStr(kFormulaStringConversionNormal))&" because it would create a circular dependance."));
        if (&theF!=right_hand_side) {
            delete right_hand_side;
        }
        return;
    }

    varFlags &= HY_VARIABLE_SET;

    if (varFlags & HY_VARIABLE_CHANGED) {
        varFlags -= HY_VARIABLE_CHANGED;
    }

    if (varFlags & HY_DEP_V_COMPUTED) {
        varFlags -= HY_DEP_V_COMPUTED;
    }
    
    if (varFormula) {
        delete (varFormula);
        varFormula = nil;
    } else {
        changeMe = true;
    }

    if (varValue) {
        DeleteAndZeroObject (varValue);
    }

    //_Formula::Duplicate ((BaseRef)myF);
    varFormula = new _Formula;
    varFormula->Duplicate (right_hand_side);

    // mod 20060125 added a call to simplify constants
    varFormula->SimplifyConstants ();

    // also update the fact that this variable is no longer independent in all declared
    // variable containers which hold references to this variable
    if (changeMe) {
          if (deferSetFormula) {
              *deferSetFormula << theIndex;
              deferIsConstant  << isAConstant;
          } else {
              
              for (AVLListXIteratorKeyValue variable_record : AVLListXIterator (&variableNames)) {
                  _Variable * theV = LocateVar(variable_record.get_value());
                  //printf ("%s\n", theV->GetName()->get_str());
                  if (theV->IsContainer()) {
                      _VariableContainer* theVC = (_VariableContainer*)theV;
                      if (theVC->SetDependance(theIndex) == -2) {
                          ReportWarning ((_String("Can't make variable ")&*GetName()&" dependent in the context of "&*theVC->GetName()&" because its template variable is bound by another relation in the global context."));
                          continue;
                      }
                  }
              }
              
              likeFuncNamesList.ForEach ([this, isAConstant] (BaseRef lf_name, unsigned long idx) -> void {
                  if (((_String*)lf_name)->nonempty()) {
                      ((_LikelihoodFunction*)likeFuncList(idx))->UpdateIndependent(theIndex,isAConstant);
                  }
              });
          }
    }

    if (&theF!=right_hand_side) {
        delete right_hand_side;
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
void  _Variable::PostMarkChanged (void) {
    varFlags &= HY_DEP_CLEAR_MASK;
}


//__________________________________________________________________________________
bool  _Variable::HasChanged (bool ignoreCats) {
    // does this variable need recomputing 
    if (varFormula) {
        if (useGlobalUpdateFlag && (varFlags&HY_DEP_V_COMPUTED)) {
            return false;
        }

       return  !varValue || varFormula->HasChanged(ignoreCats) ;

    } else {
        if (varValue && varValue->IsVariable()) {
            return varValue->HasChanged();
        }
        if (ignoreCats && IsCategory()) {
            return false;
        }
        return varFlags & HY_VARIABLE_CHANGED;
    }

}

//__________________________________________________________________________________

void _Variable::MarkDone (void) {
    if (!varFormula && (varFlags & HY_VARIABLE_CHANGED) && !(varValue && varValue->IsVariable())) {
        varFlags -= HY_VARIABLE_CHANGED;
    }
}

//__________________________________________________________________________________
HBLObjectRef    _Variable::ComputeReference (_MathObject const * context) const {
    _String reference_string (*GetName());
    reference_string = AppendContainerName(reference_string, (_VariableContainer const*)context);
    
    return new _FString (reference_string, false);
}

//__________________________________________________________________________________
_String const    _Variable::ContextFreeName(void) const {
    static const _String kDot (".");
    
    long location = theName->FindBackwards (kDot, 0L, kStringEnd);
    if (location > 0L) {
       return theName->Cut (location+1L,kStringEnd);
    }  
    return *theName;
}

//__________________________________________________________________________________
_StringBuffer&     _Variable::ContextFreeName(_StringBuffer & storage) const {
    static const _String kDot (".");
    
    long location = theName->FindBackwards (kDot, 0L, kStringEnd);
    if (location > 0L) {
       return storage.AppendSubstring(*theName,location+1L,kStringEnd);
    } else {
        storage << *theName;
    }
    return storage;
}

//__________________________________________________________________________________
_String const   _Variable::ParentObjectName(void) const {
    static const _String kDot (".");
    
    long location = theName->FindBackwards (kDot, 0L, kStringEnd);
    if (location > 0L) {
        return theName->Cut (0, location-1L);
    }
    return kEmptyString;
}

//__________________________________________________________________________________

_String const WrapInNamespace (_String const& name, _String const* context) {
  if (context) {
    return *context & '.' & name;
  }
  return name;
}

//__________________________________________________________________________________
long    DereferenceString (HBLObjectRef v, _MathObject const * context, char reference_type){
    if (v && v->ObjectClass () == STRING) {
        _FString * value = (_FString*)v;
        _String referencedVariable = value->get_str();
        if (reference_type == kStringLocalDeference && context) {
            referencedVariable = AppendContainerName(referencedVariable, (_VariableContainer*)context);
        }
        return LocateVarByName(referencedVariable);
    }
    return -1;
}

//__________________________________________________________________________________
long    DereferenceVariable (long index, _MathObject const * context, char reference_type){
    if (reference_type == kStringDirectReference) {
        return index;
    }
    
    return  DereferenceString (FetchObjectFromVariableByTypeIndex(index, STRING), context, reference_type);
}

