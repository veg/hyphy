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
#include "variablecontainer.h"
#include "operation.h"

#include "likefunc.h"
#include "parser.h"
#include "polynoml.h"
#include "batchlan.h"

#include "global_object_lists.h"
#include "global_things.h"

#include "function_templates.h"


using namespace hy_global;
using namespace hyphy_global_objects;


//__________________________________________________________________________________

_VariableContainer::_VariableContainer (void) : theParent(nil), theModel (-1L), iVariables(nil), dVariables(nil), gVariables(nil){
 }

//__________________________________________________________________________________

void    _VariableContainer::Duplicate (BaseRefConst theO) {
    _Variable::Duplicate (theO);
    _VariableContainer const *theVC = (_VariableContainer const*)theO;
    theParent= theVC->theParent;
    theModel = theVC->theModel;
    if (theVC->iVariables) {
        if (iVariables) {
            iVariables->Clear();
        } else {
            iVariables = new _SimpleList;
        }
        iVariables->Duplicate (theVC->iVariables);
    } else {
        if (iVariables) {
            delete (iVariables);
            iVariables = nil;
        }
    }
    if (theVC->dVariables) {
        if (dVariables) {
            dVariables->Clear();
        } else {
            dVariables = new _SimpleList;
        }
        dVariables->Duplicate (theVC->dVariables);
    } else {
        if (dVariables) {
            delete (dVariables);
            dVariables = nil;
        }
    }
    if (theVC->gVariables) {
        if (gVariables) {
            gVariables->Clear();
        } else {
            gVariables = new _SimpleList;
        }
        gVariables->Duplicate (theVC->gVariables);
    } else {
        if (gVariables) {
            delete (gVariables);
            gVariables = nil;
        }
    }
}


//__________________________________________________________________________________

void    _VariableContainer::TrimMemory () {
    _SimpleList* var_arrays [3] = {iVariables,dVariables,gVariables};
    for (_SimpleList* an_array : var_arrays) {
        if (an_array) {
            an_array->TrimMemory();
        }
    }
}

//__________________________________________________________________________________

BaseRef _VariableContainer::makeDynamic (void) const {
    _VariableContainer * res = new _VariableContainer;
    res->Duplicate(this);
    return res;
}

//__________________________________________________________________________________

BaseRef _VariableContainer::toStr (unsigned long) {
    _StringBuffer * res = new _StringBuffer (128L);

    (*res) << "Container Class:" << theName << ":{ Independent Variables:";

    if (iVariables)
        for (long i = 0L; i<iVariables->lLength; i+=2L) {
            res->AppendNewInstance ((_String*)variablePtrs(iVariables->lData[i])->toStr());

            if (i<iVariables->lLength-2) {
                (*res) << ',';
            }
            
        }

    (*res) << "; Dependent Variables:";

    if (dVariables)
        for (long i2 = 0L; i2<dVariables->lLength; i2+=2L) {
            res->AppendNewInstance ((_String*)variablePtrs(dVariables->lData[i2])->toStr());
            if (i2<dVariables->lLength-2) {
                (*res) << ',';
            }
        }

    (*res) << '}';
    res->TrimSpace();
    return res;
}

//__________________________________________________________________________________

_VariableContainer::_VariableContainer (_String theName, _String theTmplt, _VariableContainer* theP) : iVariables(nil), dVariables(nil), gVariables(nil) {
    InitializeVarCont (theName, theTmplt,theP);
}

//__________________________________________________________________________________

bool _VariableContainer::HasExplicitFormModel (void) const {
    if (theModel == -1L) {
        return false;
    }
    return (modelTypeList.lData[theModel]);
}

//__________________________________________________________________________________

_Formula* _VariableContainer::GetExplicitFormModel (void) const {
    if (theModel < 0L) {
        return nil;
    }
    if (modelTypeList.lData[theModel]) { // an explicit formula based matrix
        return (_Formula*)modelMatrixIndices.lData[theModel];
    }
    return nil;
}

//__________________________________________________________________________________

_String const* _VariableContainer::GetModelName (void)  const{
    _String const * res = GetObjectNameByType (HY_BL_MODEL, theModel, false);
    if (res) {
        return res;
    }
    return &kEmptyString;
}

//__________________________________________________________________________________

_Matrix* _VariableContainer::GetModelMatrix (_List* queue, _SimpleList* tags) const {
    if (theModel < 0L) {
        return nil;
    }

    if (modelTypeList.lData[theModel]) { // an explicit formula based matrix
        if (queue && tags) {
            long currentQueueLength = ((_Formula*)modelMatrixIndices.lData[theModel])->ExtractMatrixExpArguments (queue);
            if (currentQueueLength) {
                for (unsigned long k = 0; k < currentQueueLength; k++)
                  (*tags) << currentQueueLength;
                return nil;
            }
        }
        _Matrix* result = (_Matrix *)((_Formula *)modelMatrixIndices.lData[theModel])->Compute();
        result->CheckIfSparseEnough(true);
        return result;
    }

    return (_Matrix*) (LocateVar(modelMatrixIndices.lData[theModel])->GetValue());
}

//__________________________________________________________________________________

long _VariableContainer::GetModelDimension (void) {
    long matrixDim = 0L;
    if (theModel >= 0L) {
        matrixDim = modelTypeList.lData[theModel];
        if (matrixDim == 0L) {
            return GetModelMatrix()->GetHDim();
        }
    }
    return matrixDim;
}

//__________________________________________________________________________________

_Matrix* _VariableContainer::GetFreqMatrix (void) const  {
    if (theModel>=0) {
        long freqID = modelFrequenciesIndices.lData[theModel];
        if (freqID>=0) {
            return (_Matrix*) (LocateVar(freqID)->GetValue());
        } else {
            return (_Matrix*) (LocateVar(-freqID-1)->GetValue());
        }
    }
    return nil;
}

//__________________________________________________________________________________
void    _VariableContainer::ScanModelBasedVariables (_String& fullName, _AVLListXL* varCache) {
    if (theModel!= HY_NO_MODEL) { // build the matrix variables
        _SimpleList       mVars;
        _String           varName;
        
        {
            
            long cachedID = -1;
            bool doScan   = !varCache || (cachedID = varCache->Find ((BaseRef) theModel)) < 0L ;

            if (doScan) {

                _AVLList                ma (&mVars);
                ScanModelForVariables   (GetModelIndex(), ma,true,theModel,false);

                long freqID     = modelFrequenciesIndices.lData[theModel];
                if (freqID>=0) {
                    ((_Matrix*) (LocateVar(freqID)->GetValue()))->ScanForVariables2(ma,true,-1,false);
                }

                ma.ReorderList();

                if (varCache) {
                    varCache->Insert ((BaseRef)theModel, (long)mVars.makeDynamic(),false);
                }
            } else if (varCache) {
                mVars.Duplicate (varCache->GetXtra (cachedID));
            }

        }

        for (long i=0L; i<mVars.lLength; i++) {
            _Variable * aVar = (_Variable*)variablePtrs (mVars.lData[i]);
            if (aVar->IsGlobal()) {
                PushGlobalVariable(aVar->GetIndex());
            } else {
                varName = fullName&'.'&ContextFreeName();
                _Variable * spawnedVar = CheckReceptacle(&varName, kEmptyString, false, false);
                spawnedVar->SetBounds (aVar->GetLowerBound(), aVar->GetUpperBound());

                if (aVar->IsIndependent()) {
                    PushIndVariable(spawnedVar->GetIndex(), mVars.get(i));
                } else {
                    PushDepVariable(spawnedVar->GetIndex(), mVars.get(i));
                }
            }
        }
    }
}

//__________________________________________________________________________________
void    _VariableContainer::SetModel (long modelID, _AVLListXL* varCache) {
    theModel = modelID;
    ScanModelBasedVariables (*theName, varCache);
    SortVars();
}

//__________________________________________________________________________________
void    _VariableContainer::InitializeVarCont (_String const& aName, _String& theTmplt, _VariableContainer* theP, _AVLListXL* varCache) {
    
    theParent = theP;

    if (aName.nonempty()) {
        /*
            SLKP
            this entire section may be deprecated, and may actually
        */
        theName = new _String (aName);

        long f = aName.Find('.');

        while (theP) {
            if (f != -1L) {
                f = aName.Find('.',f+1L,-1L);
            } else {
                break;
            }
            theP = theP->theParent;
        }

        if (theP) {
            *theName = (*(theP->theName))&'.'&*theName;
        }

        InsertVar (this);
    }
    SetModel (FindModelName(theTmplt), varCache);
}

//__________________________________________________________________________________

void _VariableContainer::ScanAndAttachVariables (void) {
    _Variable* var;
    _SimpleList travcache;

    long f = variableNames.Find (theName,travcache);
    
    if (f >= 0L) {
        _String theNameAndADot = *theName & '.';

        for (f = variableNames.Next (f, travcache); f>=0; f = variableNames.Next (f, travcache)) {
            var = FetchVar (f);

            if (var->theName->BeginsWith(theNameAndADot)) {
                if (!var->IsContainer()) {
                    long   vix = variableNames.GetXtra (f);
                    
                    if (var->IsIndependent()) {
                        if ( ! (HasIndVariable(vix) || HasDepVariable(vix))) {
                            PushIndVariable(vix, -1);
                        }
                    } else {
                        if ( ! (HasIndVariable(vix) || HasDepVariable(vix))) {
                            PushDepVariable(vix, -1);
                       }
                    }
                }
            } else {
                break;
            }
        }
    }

}
//__________________________________________________________________________________

_VariableContainer::~_VariableContainer(void) {
    if (iVariables) {
        delete iVariables;
    }
    if (dVariables) {
        delete dVariables;
    }
    if (gVariables) {
        delete gVariables;
    }
}

//__________________________________________________________________________________

bool _VariableContainer::HasChanged (bool) {
 
    auto has_changed = [] (long var_index, long, unsigned long) -> bool {
        return LocateVar (var_index) -> HasChanged ();
    };
    auto has_changed_global = [=] (long var_index, unsigned long) -> bool {
        return LocateVar (var_index) -> HasChanged ();
    };

    return AnyLocalVariable (iVariables, has_changed) ||
           gVariables && gVariables->Any(has_changed_global) ||
           AnyLocalVariable (dVariables, has_changed);
}

//__________________________________________________________________________________

_Variable* _VariableContainer::GetIthIndependent (long index) const {
    if (iVariables && (index = index << 1)<iVariables->countitems()) {
        return LocateVar (iVariables->get(index));
    } else {
        return nil;
    }
}

//__________________________________________________________________________________

_Variable* _VariableContainer::GetIthDependent (long index) const {
    if (dVariables && (index = index << 1) < dVariables->countitems()) {
        return LocateVar (dVariables->get(index));
    } else {
        return nil;
    }
}

//__________________________________________________________________________________

_Variable* _VariableContainer::GetIthParameter (long index) const {
    if (iVariables) {
        if ( (index = index << 1 ) <iVariables->countitems()) {
            return LocateVar (iVariables->get(index));
        } else {
            if (dVariables) {
                index-=iVariables->countitems();
                if (index<dVariables->countitems()) {
                    return LocateVar (dVariables->get(index));
                }
            }
        }
    } else {
        if (dVariables && (index = index << 1) <dVariables->countitems()) {
            return LocateVar (dVariables->get(index));
        }
    }
    return nil;
}

//__________________________________________________________________________________

bool _VariableContainer::NeedToExponentiate (bool ignoreCats) const {
    if (HY_VC_NO_CHECK&varFlags) {
        return false;
    }

    auto has_changed = [=] (long var_index, long ref_index, unsigned long) -> bool {
        if (ref_index >= 0L) {
            return LocateVar (var_index) -> HasChanged (ignoreCats);
        }
        return false;
    };
    auto has_changed_global = [=] (long var_index, unsigned long) -> bool {
        return LocateVar (var_index) -> HasChanged (ignoreCats);
    };

    
    return AnyLocalVariable (iVariables, has_changed) ||
           gVariables && gVariables->Any(has_changed_global) ||
           AnyLocalVariable (dVariables, has_changed);
    
}

//__________________________________________________________________________________
void      _VariableContainer::SortVars(void) {
    // sort independents 1st
    // use dumb bubble sort
    
    auto bubble_sort = [] (_SimpleList * array) -> void {
        if (array && array->countitems ()>2) {
            bool        done = false;
            
            _String     *s1,
            *s2;
            while (!done) {
                done = true;
                s1 = LocateVar(array->lData[0])->GetName();
                for (long index = 2L; index<array->countitems(); index+=2L) {
                    s2 = LocateVar(array->lData[index])->GetName();
                    if (s2->Compare(*s1) == kCompareLess) {
                        done = false;
                        array->Swap (index, index-2);
                        array->Swap (index+1, index-1);
                    }
                }
            }
        }
    };
    
    bubble_sort (iVariables);
    bubble_sort (dVariables);
}
//__________________________________________________________________________________

void     _VariableContainer::PushGlobalVariable (long var_ref) {
    if (gVariables) {
        *gVariables << var_ref;
    } else {
        gVariables = new _SimpleList;
        *gVariables << var_ref;
    }
}

//__________________________________________________________________________________
void      _VariableContainer::PushIndVariable (long var_ref, long local_ref) {
    if (iVariables) {
        *iVariables << var_ref << local_ref;
    } else {
        iVariables = new _SimpleList;
        *iVariables << var_ref << local_ref;
    }
}
//__________________________________________________________________________________

void    _VariableContainer::PushDepVariable (long var_ref, long local_ref) {
    if (dVariables) {
        *dVariables << var_ref << local_ref;
    } else {
        dVariables = new _SimpleList;
        *dVariables << var_ref << local_ref;
    }
}
//__________________________________________________________________________________

bool    _VariableContainer::HasIndVariable  (long var_ref) const {
    return iVariables && iVariables->FindStepping(var_ref, 2L);
}
//__________________________________________________________________________________

bool    _VariableContainer::HasDepVariable  (long var_ref) const {
    return dVariables && dVariables->FindStepping(var_ref, 2L);
}
//__________________________________________________________________________________

void    _VariableContainer:: RemoveLocalVariable (_SimpleList*& array, long array_index) {
    if (array->countitems() >2UL) {
        array->Delete(array_index);
        array->Delete(array_index);
        array->TrimMemory();
    } else {
        delete array;
        array = nil;
    }
}


//__________________________________________________________________________________
bool      _VariableContainer::RemoveDependence (long varIndex) {
    if (dVariables) {
        long array_index = dVariables->FindStepping(varIndex,2L);

        if (array_index >= 0) {

            InsertVariableInSortedList(iVariables,
                                       *LocateVar (dVariables->lData[array_index])->GetName(),
                                       varIndex,
                                       dVariables->get(array_index+1));
            RemoveLocalVariable (dVariables, array_index);
        }
    }
    return true;
}

//__________________________________________________________________________________
long      _VariableContainer::CheckAndAddUserExpression (_String& parameter_name, long start_with) {
    _String localized_name = WrapInNamespace (parameter_name, theName),
            unused_name (localized_name);
    
    long    k = MAX (start_with, 2L);
    if (start_with>=2L) {
        unused_name = localized_name&start_with;
    }

    while (LocateVarByName(unused_name)>=0L) {
        unused_name = localized_name & _String (k++);
    }

    if (start_with<0L) {
        return k>2?k-1:0;
    }

    if (start_with<2) {
        if (k>2) {
            parameter_name = parameter_name&_String(k-1L);
        }
    } else {
        if (k>start_with) {
            parameter_name = parameter_name & _String (k-1L);
        } else {
            parameter_name = parameter_name & _String (start_with);
        }
    }

    _Variable newVar (unused_name);
    k =  newVar.GetIndex();

    PushDepVariable(k, -1);
    return k;
}

//__________________________________________________________________________________
void      _VariableContainer::CopyMatrixParameters (_VariableContainer* source, bool match_by_name) {
    if (iVariables && (source->iVariables || source->dVariables)) {
        if (match_by_name) {
            _List source_vars, target_vars;
            
            _SimpleList model_vars_in_source, model_vars_in_target;
            
            ForEachLocalVariable(source->iVariables, [&] (long var_idx, long ref_idx, long array_index) {
                if (ref_idx >= 0L) {
                    source_vars << LocateVar (ref_idx)->GetName();
                    model_vars_in_source << array_index;
                }
            });
            ForEachLocalVariable(source->dVariables, [&] (long var_idx, long ref_idx, long array_index) {
                if (ref_idx >= 0L) {
                    source_vars << LocateVar (ref_idx)->GetName();
                    model_vars_in_source << (-2L-array_index);
                }
            });
            ForEachLocalVariable(iVariables, [&] (long var_idx, long ref_idx, long array_index) {
                if (ref_idx >= 0L) {
                    target_vars << LocateVar (ref_idx)->GetName();
                    model_vars_in_target << array_index;
                }
            });
            
            _SimpleList the_mapping;
            target_vars.Map (source_vars, the_mapping);
            the_mapping.Each ([=] (long source_var, unsigned long index) -> void {
                if (source_var >= 0UL) {
                    long which_idx = model_vars_in_source.lData[source_var];
                    which_idx = which_idx >= 0 ? source->iVariables->get (which_idx) : source->dVariables->get (-which_idx-2L);
                    LocateVar (iVariables->get (model_vars_in_target.get(index)))->SetValue (LocateVar (which_idx)->Compute());
                }
            });
            
        } else {
            if (source->iVariables) {
                for (unsigned long i=0UL; i<iVariables->lLength && i< source->iVariables->lLength; i+=2UL) {
                    LocateVar (iVariables->get(i))->SetValue(LocateVar (source->iVariables->get(i))->Compute());
                }
            }
        }
    }
    SetValue (source->Compute());
}

//__________________________________________________________________________________
void      _VariableContainer::KillUserExpression (long varID) {
    if (dVariables) {
        long f = dVariables->FindStepping(varID,2);
        if (f>=0) {
            DeleteVariable (*LocateVar(varID)->GetName(),true);
            RemoveLocalVariable(dVariables, f);
         }
    }
}

//__________________________________________________________________________________

long    _VariableContainer::InsertVariableInSortedList (_SimpleList * & list, _String const & var_name, long var_idx, long ref_idx) {
    
    long    insert_here = 0L;

    if (!list) {
        list = new _SimpleList;
    } else {
        unsigned long array_l = list->countitems();
        while (insert_here< array_l) {
            _Variable *existing_var = LocateVar (list->get(insert_here));
            if (!existing_var) {
                HandleApplicationError ("Internal error in InsertVariableInSortedList()", false);
                return -1;
            }
            if (var_name.Compare (*existing_var->GetName()) != kCompareGreater) {
                break;
            }
            insert_here+=2;
        }
    }

    list->InsertElement ((BaseRef)var_idx, insert_here, false, false);
    list->InsertElement ((BaseRef)ref_idx, insert_here+1, false, false);
    
    return insert_here;
}


//__________________________________________________________________________________
long      _VariableContainer::SetDependence (long varIndex) {
    if (iVariables) {
        long f;

        if (varIndex>=0) {
            f = iVariables->FindStepping(varIndex,2);
            if (f<0) {
                return -1;
            }
        } else {
            f = -varIndex-1;
            varIndex = iVariables->lData[f];
        }


        //printf ("Moving ind->dep for %s from %s\n", LocateVar (varIndex)->GetName()->sData,
        //      GetName()->sData);

        if (iVariables->lData[f+1]>=0) {
            //printf ("Local variable %s\n", LocateVar (iVariables->lData[f+1])->GetName()->sData);
            if (!LocateVar(iVariables->lData[f+1])->IsIndependent()) {
                return -2;
            }
        }
        
        InsertVariableInSortedList (dVariables, *LocateVar (iVariables->get(f))->GetName(), varIndex,iVariables->get(f+1));
        RemoveLocalVariable(iVariables,f);
        return varIndex;
    }
    return -1;
}

//__________________________________________________________________________________
bool      _VariableContainer::SetMDependence (_SimpleList& mDep)
{
  if (iVariables) {
    if (mDep.lLength*2 > iVariables->lLength)
      for (long k=iVariables->lLength-2; k>=0; k-=2) {
        long f = mDep.BinaryFind (iVariables->lData[k]);
        if (f>=0) {
          SetDependence (-k-1);
        }
      }
    else
      for (unsigned long k=0UL; iVariables && k<mDep.lLength; k++) {
        SetDependence (mDep.lData[k]);
      }
  }
  
  return true;
}


//__________________________________________________________________________________
void      _VariableContainer::Clear(void) {
    theModel = HY_NO_MODEL;
    if (iVariables) {
        delete iVariables;
        iVariables = nil;
    }
    if (dVariables) {
        delete dVariables;
        dVariables = nil;
    }
    if (gVariables) {
        delete gVariables;
        gVariables = nil;
    }
}

//__________________________________________________________________________________
long      _VariableContainer::CountAll(void) const {
    return (iVariables? (iVariables->countitems() >> 1) :0L)+(dVariables?(dVariables->countitems() >> 1):0L);
}

//__________________________________________________________________________________
long      _VariableContainer::CountIndependents(void) const {
    return  (iVariables? (iVariables->countitems() >> 1) :0L);
}

//__________________________________________________________________________________
bool      _VariableContainer::HasLocals  (void) {
    return  iVariables && iVariables->countitems() > 0UL || dVariables && dVariables->countitems() > 0UL;
}

//__________________________________________________________________________________
bool      _VariableContainer::is_model_var  (long i) const {
    return dVariables->lData[2*i+1]>=0;
}

//__________________________________________________________________________________

_String*    _VariableContainer::GetSaveableListOfUserParameters (void) {
    _StringBuffer * result = new _StringBuffer (64L);
    
    ForEachLocalVariable(dVariables, [&] (long var_index, long ref_index, unsigned long array_index) -> void {
        if (ref_index < 0) {
            _Variable * userParm  = (_Variable*) LocateVar (var_index);
            result->AppendAnAssignmentToBuffer(userParm->GetName(),
                                               (_String*)userParm->GetFormulaString(kFormulaStringConversionNormal),
                                               kAppendAnAssignmentToBufferFree | kAppendAnAssignmentToBufferAssignment);
        }
    });

    result->TrimSpace ();
    return result;
}

//__________________________________________________________________________________
void      _VariableContainer::ClearConstraints(void) {
    while (dVariables) {
        LocateVar(dVariables->lData[0])->ClearConstraints();
    }
}

//__________________________________________________________________________________

void  _VariableContainer::CompileListOfDependents (_SimpleList& rec) {
    
    auto push_var = [&] (long var_idx, long ref_idx, unsigned long index) -> void {
        LocateVar (var_idx)->CompileListOfDependents(rec);
    };
    
    ForEachLocalVariable(iVariables, push_var);
    if (gVariables) {
        gVariables->Each ([&] (long var_idx, unsigned long index) -> void {
            LocateVar (var_idx)->CompileListOfDependents(rec);
        });
    }
    ForEachLocalVariable(dVariables, push_var);
    if (dVariables) {
        for (unsigned long i=0UL; i<dVariables->countitems(); i+=2UL) {
            long f = rec.Find (dVariables->get (i));
            if (f>=0L) {
                rec.Delete (f);
            }
        }
    }
}


//__________________________________________________________________________________

void _VariableContainer::MarkDone (void) {
    
    ForEachLocalVariable(iVariables, [] (long var_idx, long ref_idx, unsigned long index) -> void {
        LocateVar (var_idx)->MarkDone();
    });
    if (gVariables) {
        gVariables->Each ([&] (long var_idx, unsigned long index) -> void {
            LocateVar (var_idx)->MarkDone();
        });
    }
}

//__________________________________________________________________________________

void _VariableContainer::MatchParametersToList (_List& suffixes, bool doAll, bool indOnly) {
    /** TODO SLKP 20171130; what is this for?? Likely can be deprecated */
    
    if (doAll) {
        for (long i=suffixes.lLength-1; i>=0; i--) {
            long j;
            if (!indOnly) {
                if (dVariables) {
                    for (j=0; j<dVariables->lLength; j+=2)
                        if (LocateVar(dVariables->lData[j])->GetName()->EndsWith (*(_String*)suffixes.lData[i])) {
                            break;
                        }

                    if (j<dVariables->lLength) {
                        continue;
                    }
                }
            }
            if (iVariables) {
                for (j=0; j<iVariables->lLength; j+=2) {
                    if (LocateVar(iVariables->lData[j])->GetName()->EndsWith (*(_String*)suffixes.lData[i])) {
                        break;
                    }
                }
                if (j==iVariables->lLength) {
                    suffixes.Delete (i);
                }
            } else {
                suffixes.Delete (i);
            }
        }
    } else {
        for (long i=suffixes.lLength-1; i>=0; i--) {
            long j;
            if (dVariables) {
                for (j=0; j<dVariables->lLength; j+=2) {
                    if (dVariables->lData[j+1]<0) {
                        if (LocateVar(dVariables->lData[j])->GetName()->EndsWith (*(_String*)suffixes.lData[i])) {
                            break;
                        }
                    }
                }
                if (j==dVariables->lLength) {
                    suffixes.Delete (i);
                }
            } else {
                suffixes.Delete(i);
            }
        }
    }
}

//__________________________________________________________________________________

bool _VariableContainer::IsConstant (void) const {
    if (iVariables) {
        return false;
    }
    
    return ! (AnyLocalVariable(dVariables, [] (long var_idx, long ref_idx, unsigned long) -> bool {
        return !LocateVar(var_idx)->IsConstant();
    }) ||
        gVariables && gVariables->Any ([] (long var_idx,  unsigned long) -> bool {
        return !LocateVar(var_idx)->IsConstant();
    }));

}

//__________________________________________________________________________________

void _VariableContainer::ScanContainerForVariables (_AVLList& l,_AVLList& l2, _AVLListX * tagger, long weight) {
    
    ForEachLocalVariable(iVariables, [&] (long var_idx, long ref_idx, unsigned long) -> bool {
        l.Insert((BaseRef)var_idx);
        if (tagger) {
            tagger->UpdateValue ((BaseRef)var_idx, weight, 0);
        }
    });
    
    if (dVariables) {
        _SimpleList temp;
        _AVLList  ta (&temp);

        ForEachLocalVariable(dVariables, [&] (long var_idx, long ref_idx, unsigned long) -> bool {
            l2.Insert((BaseRef)var_idx);
            LocateVar (var_idx)->ScanForVariables(ta, true, tagger, weight);
        });

        ta.ReorderList();
        temp.Each([&] (long var_index, unsigned long) -> void {
            _Variable * v = LocateVar(var_index);
            if (!v->IsGlobal() && v->IsIndependent()) {
                l.Insert ((BaseRef)var_index);
                if (tagger) {
                    tagger->UpdateValue ((BaseRef)var_index, weight, 0);
                }
            }
        });
    }
}

//__________________________________________________________________________________

void _VariableContainer::ScanForDVariables (_AVLList& l,_AVLList&) const {
    ForEachLocalVariable(dVariables, [&] (long var_idx, long ref_idx, unsigned long) -> bool {
        l.Insert((BaseRef)var_idx);
    });
}

//__________________________________________________________________________________

void _VariableContainer::GetListOfModelParameters (_List& rec) {
    ForEachLocalVariable(iVariables, [&] (long var_idx, long ref_idx, unsigned long) -> bool {
        if (ref_idx >= 0) {
            rec << LocateVar(ref_idx)->GetName();
        }
    });
}

//__________________________________________________________________________________

void _VariableContainer::ScanForGVariables (_AVLList& independent,_AVLList& dependent, _AVLListX* tagger, long weight) const {
    
    
    auto insert_g_var = [&] (_Variable *v, long var_idx) -> void {
        if (v->IsIndependent()) {
            independent.Insert ((BaseRef)var_idx);
            if (tagger) {
                tagger->UpdateValue((BaseRef)var_idx, weight, 0);
            }
        } else {
            dependent.Insert ((BaseRef)var_idx);
        }
    };
    
    if (gVariables) {
        gVariables->Each([&] (long var_idx, unsigned long) -> void {
            insert_g_var (LocateVar (var_idx), var_idx);
        });
    }
    
    
    // additionally, check to see if there is any implicit dependence on the global variables yet unseen
    if (dVariables) {
        _SimpleList var_list;
        _AVLList  al (&var_list);
        ForEachLocalVariable(dVariables, [&] (long var_idx, long ref_idx, unsigned long) -> void {
            LocateVar (var_idx)->ScanForVariables(al, true);
        });
        al.ReorderList();
        var_list.Each ([&] (long var_idx, unsigned long) -> void {
            _Variable * v = LocateVar(var_idx);
            if (v->IsGlobal()) {
                insert_g_var (LocateVar (var_idx), var_idx);
            }
        });
     }
}
