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

#include <math.h>
#include <float.h>
#include <limits.h>
#include "string.h"
#include <stdio.h>
#include <time.h>
#include <ctype.h>

#include "likefunc.h"
#include "parser.h"
#include "matrix.h"
#include "polynoml.h"
#include "batchlan.h"
#include "global_things.h"

using namespace hy_global;




extern
_SimpleList     simpleOperationCodes,
                simpleOperationFunctions;

_List           hyReservedWords,
                varNamesSupportList,
                variablePtrs;   // stores all the variables declared so far

_AVLList        *lookAside = nil;
_AVLListX       variableNames (&varNamesSupportList);




// indices of all independent variables

_Trie           FunctionNameList;
_List           BuiltInFunctions;



_SimpleList     freeSlots,
                deferIsConstant,
                *deferSetFormula = nil;

bool            useGlobalUpdateFlag = false;


_String         HalfOps (":<>=!&|");

_Trie           UnOps;

_SimpleList opPrecedence,
            BinOps,
            associativeOps;

hyFloat pi_const = 3.141592653589793,
long_max = (hyFloat)LONG_MAX;

/**********************************/
/* Defining Globals here for now */
/********************************/
 
//Used in formula, and constant

hyFloat  tolerance  = DBL_EPSILON;



/*********************************/
/*          End Globals         */
/*******************************/

//__________________________________________________________________________________
//SW: Helper functions




//__________________________________________________________________________________
void     parameterToCharBuffer (hyFloat value, char* dump, long length, bool json)
{
    if (json) {
      if (isnan (value)) {
        snprintf (dump, length, "null");
        return;
      }
      if (isinf(value)) {
        snprintf (dump, length, value < 0 ? "-1e9999" : "1e9999");
        return;
      }
    }
  
    long digs = print_digit_specification;
    if (digs<=0 || digs>15) {
        if (round(value) == value && fabs (value) < long_max) {
            snprintf (dump,length, "%ld",lrint (value));
        } else {
            snprintf (dump,length, PRINTF_FORMAT_STRING,value);
        }
    } else {
        _String format("%-");
#ifdef __USE_LONG_DOUBLE__
        format = format&_String(digs)&"Lg";
#else
        format = format&_String(digs)&'g';
#endif
        snprintf (dump,length,(const char*)format,value);
    }
}


//__________________________________________________________________________________
BaseRef     parameterToString (hyFloat value)
{
    char dump [256];
    parameterToCharBuffer (value, dump, 256);
    return new _String (dump);
}

//__________________________________________________________________________________
void SplitVariableIDsIntoLocalAndGlobal (const _SimpleList& theList, _List& splitStorage)
{
    splitStorage.Clear();
    splitStorage.AppendNewInstance(new _SimpleList);
    splitStorage.AppendNewInstance(new _SimpleList);

    for (unsigned long k=0; k<theList.lLength; k++) {
        long varID = theList.list_data[k];
        (*(_SimpleList*)splitStorage(1-LocateVar (varID)->IsGlobal())) << varID;
    }
}

//__________________________________________________________________________________
_String FetchObjectNameFromType (const unsigned long objectClass) {
    
    
    if (objectClass == HY_ANY_OBJECT) {
        return "Any object";
    }
    _StringBuffer result;
    
    auto push_with_spacer = [&result] (const char * v) -> void {
        if (result.nonempty()) {
            result << " or ";
        }
        result << v;
    };
    
    if (objectClass & HY_UNDEFINED) {
        result << "Undefined";
    }
    if (objectClass & NUMBER) {
        push_with_spacer ("Number");
    }
    if (objectClass & MATRIX) {
        push_with_spacer ("Matrix");
    }
    if (objectClass & CONTAINER) {
        push_with_spacer ("Variable container");
    }
    if (objectClass & TREE_NODE) {
        push_with_spacer ("Tree Node");
    }
    if (objectClass & TREE) {
        push_with_spacer ("Tree");
    }
    if (objectClass & STRING) {
        push_with_spacer ("String");
    }
    if (objectClass & ASSOCIATIVE_LIST) {
        push_with_spacer ("Associative Array");
    }
    if (objectClass & TOPOLOGY) {
        push_with_spacer ("Topology");
    }
    if (objectClass & POLYNOMIAL) {
        push_with_spacer ("Polynomial");
    }
    
    return result;
}

//__________________________________________________________________________________
_String*   FetchMathObjectNameOfTypeByIndex (const unsigned long objectClass, const long objectIndex)
{
    if (objectIndex >=0 && objectIndex < variableNames.countitems()) {
            long tc = 0;
            _SimpleList nts;
            long        rt,
                        vi = variableNames.Traverser (nts, rt, variableNames.GetRoot());

            for (; vi >= 0; vi = variableNames.Traverser (nts, rt))
                if (FetchVar(variableNames.GetXtra (vi))->ObjectClass () == objectClass) {
                    if (tc==objectIndex) {
                        return (_String*)variableNames.Retrieve(vi);
                        break;
                    } else {
                        tc++;
                    }
                }    
    }
    return nil;
}

//__________________________________________________________________________________
HBLObjectRef   FetchObjectFromFormulaByType (_Formula& f, const unsigned long objectClass, long command_id, _String *errMsg) {
  HBLObjectRef v = f.Compute();
  if (v) {
    if (objectClass == HY_ANY_OBJECT || v->ObjectClass () == objectClass) {
      return v;
    }
    if (command_id >= 0 || errMsg) {
      if (command_id >= 0) {
        HandleApplicationError (_String ((_String*)f.toStr(kFormulaStringConversionNormal)).Enquote() & (" must evaluate to a ") & FetchObjectNameFromType (objectClass) & " in call to "
                   &_HY_ValidHBLExpressions.RetrieveKeyByPayload(command_id) & '.');
      }
    }
  }
  return nil;
}

//__________________________________________________________________________________
HBLObjectRef   FetchObjectFromVariableByType (_String const* id, const unsigned long objectClass, long command_id, _String *errMsg) {
    if (id) {
        _Variable * v = FetchVar (LocateVarByName (*id));
        if (v && (objectClass == HY_ANY_OBJECT || (v->ObjectClass () & objectClass))) {
            return v->Compute();
        }
        if (command_id >= 0 || errMsg) {
            if (command_id >= 0) {
                HandleApplicationError (id->Enquote() & (" must refer to a ") & FetchObjectNameFromType (objectClass) & " in call to "
                                         &_HY_ValidHBLExpressions.RetrieveKeyByPayload(command_id) & '.');
            } else {
                HandleApplicationError (errMsg->Replace ("_VAR_NAME_ID_", *id, true));
            }
        }
    }
    return nil;
}


//__________________________________________________________________________________
HBLObjectRef   FetchObjectFromVariableByTypeIndex (long idx, const unsigned long objectClass, long command_id, _String *errMsg) {
    _Variable * v = FetchVar (idx);
    if (v) {
        if (objectClass == HY_ANY_OBJECT || v->ObjectClass () == objectClass) {
            return v->Compute ();
        }
        if (command_id >= 0 || errMsg) {
            if (command_id >= 0) {
                HandleApplicationError (v->GetName()->Enquote() & (" must refer to a ") & FetchObjectNameFromType (objectClass) & " in call to "
                                         &_HY_ValidHBLExpressions.RetrieveKeyByPayload(command_id) & '.');
            } else {
                HandleApplicationError (errMsg->Replace ("_VAR_NAME_ID_", *v->GetName(), true));
            }
        }
    }
    return nil;
}

//__________________________________________________________________________________
long LocateVarByName (_String const& name) {
    return variableNames.Find (&name);
}

//__________________________________________________________________________________
_Variable* FetchVar (long index, unsigned long type_check) {
    if (index >= 0) {
        _Variable * var = (_Variable *)variablePtrs.GetItemRangeCheck(variableNames.GetXtra(index));
        if (var) {
            if (type_check == HY_ANY_OBJECT || (var->ObjectClass() & type_check)) {
                return var;
            }
        }
    }
    return nil;
}

//__________________________________________________________________________________
void       UpdateChangingFlas (long vN) {
    // check to see if formulae contain a reference to this var
    // if so: "decompile" them

    long topLimit = compiledFormulaeParameters.countitems();

    _SimpleList * toDelete = nil;

    for (long k = 0L; k<topLimit; k++) {
        long g = ((_SimpleList*)compiledFormulaeParameters.list_data[k])->BinaryFind (vN,0);

        if (g>=0) {
            ((_ElementaryCommand*)listOfCompiledFormulae.list_data[k])->DecompileFormulae();
          

            if (!toDelete) {
                toDelete = new _SimpleList;
            }

            *toDelete << k;
        }
    }

    if (toDelete) {
        listOfCompiledFormulae.DeleteList       (*toDelete);
        compiledFormulaeParameters.DeleteList   (*toDelete);
        DeleteObject (toDelete);
    }
}

//__________________________________________________________________________________
void       UpdateChangingFlas (_SimpleList const & involvedVariables) {
    long          topLimit         = compiledFormulaeParameters.lLength;
    _SimpleList * toDelete         = nil;

    for (long k = 0; k<topLimit; k++) {
        long g = ((_SimpleList*)compiledFormulaeParameters.list_data[k])->CountCommonElements (involvedVariables,true);

        if (g>0) {
            ((_ElementaryCommand*)listOfCompiledFormulae.list_data[k])->DecompileFormulae();

            if (!toDelete) {
                toDelete = new _SimpleList;
            }

            *toDelete << k;
        }
    }

    if (toDelete) {
        listOfCompiledFormulae.DeleteList       (*toDelete);
        compiledFormulaeParameters.DeleteList   (*toDelete);
        DeleteObject (toDelete);
    }
}

//__________________________________________________________________________________
void DeleteVariable (long dv, bool deleteself, bool do_checks) {
    if (dv>=0L) {

        _String *name   = (_String*)variableNames.Retrieve (dv);
        _String my_name = *name&'.';
        long    vidx    = variableNames.GetXtra (dv);

        UpdateChangingFlas (vidx);

        _SimpleList recCache;
        variableNames.Find (name,recCache);
        _String     nextVarID;// = *(_String*)variableNames.Retrieve(variableNames.Next (dv,recCache));
        long        nvid;
        if ((nvid = variableNames.Next (dv,recCache))>=0L) {
            nextVarID = *(_String*)variableNames.Retrieve(nvid);
        }

        if (deleteself) {
            
            
            _Variable * self_variable = FetchVar(dv);
            
            if (do_checks) {
                for (AVLListXIteratorKeyValue variable_iterator : AVLListXIterator (&variableNames)) {
                    _Variable * check_variable = FetchVar(variable_iterator.get_value());
                    if (self_variable != check_variable && check_variable->CheckFForDependence (vidx,false)) {
                        HBLObjectRef current_variable_value = check_variable->Compute();
                        current_variable_value->AddAReference(); // if this isn't done; the object will be deleted when the formula is cleared in SetValue
                        check_variable->SetValue (current_variable_value,true,true,NULL);
                        DeleteObject (current_variable_value);
                    }
                }
            }
        

            variableNames.Delete (variableNames.Retrieve(dv),true);
            variablePtrs[vidx] = nil;
            DeleteObject (self_variable);
            freeSlots<<vidx;
        } else {
            _Variable* delvar = (FetchVar(dv));
            if (delvar->IsContainer()) {
                _VariableContainer* dc = (_VariableContainer*)delvar;
                dc->Clear();
            }
        }
        
        _List       toDelete;

        recCache.Clear();
        long nextVar = variableNames.Find (&nextVarID,recCache);

        for (; nextVar>=0; nextVar = variableNames.Next (nextVar, recCache)) {
            _String dependent = *(_String*)variableNames.Retrieve (nextVar);
            if (dependent.BeginsWith(my_name)) {
                toDelete && & dependent;
            } else {
                break;
            }
        }

        for (unsigned long k=0; k< toDelete.lLength; k++) {
            DeleteVariable (*(_String*)toDelete(k));
        }
    }
}

//__________________________________________________________________________________
void DeleteTreeVariable (long tree_variable_index, _SimpleList & parms, _String * prefix, bool leftovers) {
    if (tree_variable_index>=0) {
        
        /*
            1. Indentify all the variables that are in the tree.* and tree.*.* domain
            2. Scan all the variables that are dependend on the list here and "batch" unconstrain them
            3. Delete the tree and tree nodes; store the remainign variables in parms
        */
        
        
        _String const *tree_name  = prefix ? prefix : (_String const*)variableNames.Retrieve (tree_variable_index),
                       kDot = ".",
                       tree_prefix = *tree_name&kDot;
        
        _SimpleList traversal_stack,
                    indices_to_delete,
                    all_subdomain_indices;
        
        _AVLList    indices_to_delete_avl (&indices_to_delete),
                    all_subdomain_indices_avl (&all_subdomain_indices),
                    will_keep_these (&parms);
                    
        
        long tree_array_index = tree_variable_index;
        variableNames.Find (tree_name,traversal_stack);
        
        while (long next_var = variableNames.Next (tree_array_index,traversal_stack)) {
            if (next_var < 0) {
                break;
            }
            _String const * next_var_name = (_String const *)variableNames.Retrieve(next_var);
            if (next_var_name->BeginsWith(tree_prefix)) {
                all_subdomain_indices_avl.InsertNumber(next_var);
                
                _Variable * secondary_variable = FetchVar (next_var);
                                
                if (next_var_name->FindBackwards(kDot, tree_prefix.length() + 1, kStringEnd) != kNotFound && !secondary_variable->IsContainer()) {
                    // 20200220: deleting all containers here; only "plain" variables will be kept
                    // subnode level
                    will_keep_these.InsertNumber(next_var);
                 } else {
                    //printf ("Deleting variable %s\n" , next_var_name->get_str());
                    indices_to_delete_avl.InsertNumber (next_var);
                }
            } else {
                break;
            }
            tree_array_index = next_var;
        }
        
        //if ((nvid = variableNames.Next (dv,recCache))>=0)
        
        all_subdomain_indices_avl.ReorderList();
        indices_to_delete_avl.ReorderList();
        will_keep_these.ReorderList();
        
        _SimpleList _touched_dependent_variables,
                    touched_containers;
        
        _AVLList    touched_dependent_variables (&_touched_dependent_variables);
        
        DoForEachVariable([&touched_dependent_variables, &indices_to_delete_avl, &all_subdomain_indices_avl] (_Variable* v, long idx) -> void {
            // only care about variables that are not going in the subdomain
            if (all_subdomain_indices_avl.FindLong(idx) < 0) {
                if (v->IsContainer()) {
                    //_VariableContainer * vc = (_VariableContainer*)vc;
                    /* TODO, it is theoretically possible that a node / tree variable is a template variable
                       in another container, but there are no use cases of this;
                       need to check, but later
                     */
                } else {
                    if (v->CheckFForDependence(indices_to_delete_avl)) {
                        //printf ("Touched dependent variable %s\n" , v->GetName()->get_str());
                        touched_dependent_variables.InsertNumber(v->get_index());
                    }
                }
            }
        });
        

        
        //SetVariableToOwn
        
        UpdateChangingFlas(all_subdomain_indices);
        // reset compiled formulas that may contain variable indices being deleted
        
        if (leftovers) {
            SetVariablesToOwnValues (will_keep_these);
        } else {
            indices_to_delete << parms;
        }
        
        if (touched_dependent_variables.countitems()) {
            SetVariablesToOwnValues (touched_dependent_variables);
        }
        
        /*if (prefix) {
            printf ("CLEARING %d variables\n", indices_to_delete.countitems() );
        }*/
        
        indices_to_delete.Each ([] (long var_idx, unsigned long) -> void {
            _Variable * delvar = LocateVar (var_idx);
            //printf ("Deleting variable %s\n" , delvar->GetName()->get_str());
            if (delvar->ObjectClass() != TREE) {
                variableNames.Delete (delvar->GetName(),true);
                (*((_SimpleList*)&variablePtrs))[delvar->get_index()]=0;
                freeSlots<<delvar->get_index();
                DeleteObject (delvar);
            } else {
                ((_VariableContainer*)delvar)->Clear();
            }
        });
        

        return;
         
        /*_String const *tree_name  = (_String const*)variableNames.Retrieve (dv),
                       tree_prefix = *tree_name&".";
        
        long    variable_index   = variableNames.GetXtra (dv);

        UpdateChangingFlas (variable_index);

        _SimpleList recCache;
        variableNames.Find (tree_name,recCache);
        _String     nextVarID;
        long        nvid;
        if ((nvid = variableNames.Next (dv,recCache))>=0) {
            nextVarID = *(_String*)variableNames.Retrieve(nvid);
        }


        
        {
            _SimpleList tcache;
            long        iv,
                        k = variableNames.Traverser (tcache, iv, variableNames.GetRoot());

            for (; k>=0; k = variableNames.Traverser (tcache, iv)) {
                _Variable * thisVar = FetchVar(k);

                if (thisVar->CheckFForDependence (variable_index,false)) {
                    HBLObjectRef curValue = thisVar->Compute();
                    curValue->AddAReference();
                    thisVar->SetValue (curValue, false);
                    //DeleteObject (curValue);
                }
            }
        }

        _Variable* delvar = (FetchVar(dv));
        if (delvar->ObjectClass() != TREE) {
            variableNames.Delete (variableNames.Retrieve(dv),true);
            (*((_SimpleList*)&variablePtrs))[variable_index]=0;
            freeSlots<<variable_index;
            DeleteObject (delvar);
        } else {
            ((_VariableContainer*)delvar)->Clear();
        }
        if (doDeps) {
            _List toDelete;
            recCache.Clear();
            long nextVar = variableNames.Find (&nextVarID,recCache);
            for (; nextVar>=0; nextVar = variableNames.Next (nextVar, recCache)) {
                _String dependent = *(_String*)variableNames.Retrieve (nextVar);
                if (dependent.BeginsWith(tree_prefix)) {
                    if (dependent.Find ('.', tree_prefix.length()+1, -1)>=0) {
                        _Variable * checkDep = FetchVar (nextVar);
                        if (!checkDep->IsIndependent()) {
                            HBLObjectRef curValue = checkDep->Compute();
                            curValue->AddAReference();
                            checkDep->SetValue (curValue);
                            DeleteObject (curValue);
                        }
                        parms << variableNames.GetXtra (nextVar);
                    } else {
                        toDelete && & dependent;
                    }
                } else {
                    break;
                }
            }

            for (unsigned long k=0; k<toDelete.lLength; k++) {
                //StringToConsole (*(_String*)toDelete(k));
                //BufferToConsole ("\n");
                DeleteTreeVariable (LocateVarByName(*(_String*)toDelete(k)),parms,false);
            }
        }*/
    }
}
//__________________________________________________________________________________
void DeleteVariable (_String const &name, bool deleteself) {
    DeleteVariable(LocateVarByName (name), deleteself);
}


//__________________________________________________________________________________
_Variable* CheckReceptacle (_String const * name, _String const & fID, bool checkValid, bool isGlobal, bool clear_trees) {
    if (checkValid && (!name->IsValidIdentifier(fIDAllowCompound))) {
        HandleApplicationError(name->Enquote() & " is not a valid variable identifier");
        return nil;
    }

    long    f = LocateVarByName (*name);

    if (f>=0L) {
      _Variable * existing = FetchVar (f);
    
      if (clear_trees && (existing->ObjectClass() == TREE || existing->ObjectClass() == TOPOLOGY)) {
        DeleteVariable (*existing->GetName());
        f = -1L;
      }
      else {
        return existing;
      }
    }

    {
        _Variable dummy (*name, isGlobal);
        f = LocateVarByName (*name);
    }

    return FetchVar(f);
}

//__________________________________________________________________________________
_Variable* CheckReceptacleCommandIDException (_String const* name, const long id, bool checkValid, bool isGlobal, _ExecutionList* context) {
    // TODO: allow ^name and such to constitute valid run-time references
   if (checkValid && (!name->IsValidIdentifier(fIDAllowCompound))) {
    throw  (name->Enquote('\'') & " is not a valid variable identifier");
   }
  
  long    f = LocateVarByName (*name);
  
  if (f>=0L) {
    _Variable * existing = FetchVar (f);
    if (existing->ObjectClass() == TREE) {
      DeleteVariable (*existing->GetName());
      f = -1L;
    }
    else {
      return existing;
    }
  }
  
  {
    _Variable (*name, isGlobal);
    f = LocateVarByName (*name);
  }
  
  return FetchVar(f);
}


//__________________________________________________________________________________
_Variable* CheckReceptacleCommandID (_String const* name, const long id, bool checkValid, bool isGlobal, _ExecutionList* context) {
  try {
    return CheckReceptacleCommandIDException (name, id, checkValid, isGlobal, context);
  } catch (_String const& err_msg) {
    if (context) {
      context->ReportAnExecutionError(err_msg);
    } else {
      HandleApplicationError (err_msg);
    }
  }
  return nil;
}

//__________________________________________________________________________________
bool CheckReceptacleCommandIDAndStore (_String const* name, const long id, bool checkValid, HBLObjectRef v, bool dup, bool isGlobal)
{
    _Variable *theV = CheckReceptacleCommandID (name, id, checkValid, isGlobal);
    if (theV) {
        theV->SetValue (v, dup, true, NULL);
        return true;
    }
    if (!dup) {
        DeleteObject (v);
    }
    return false;
}


//__________________________________________________________________________________
bool CheckReceptacleAndStore (_String const* name, _String fID, bool checkValid, HBLObjectRef v, bool dup) {
    _Variable * theV = CheckReceptacle(name, fID, checkValid);
    if (theV) {
        theV->SetValue (v, dup, true, NULL);
        return true;
    }
    if (!dup) {
        DeleteObject (v);
    }
    return false;
}

//__________________________________________________________________________________
bool CheckReceptacleAndStore (_String name, _String fID, bool checkValid, HBLObjectRef v, bool dup)
{
    return CheckReceptacleAndStore (&name, fID, checkValid, v, dup);
}

//__________________________________________________________________________________
void  InsertVar (_Variable* theV) {
        
    long pos = variableNames.Insert (theV->theName);

    if (pos < 0 && isDefiningATree == kTreeNodeBeingCreated)
        // automatically fix duplicate autogenerated tree node name
    {
        long trySuffix  = 1;
        _String * tryName = new _String;
        do {
            *tryName = *theV->theName & "_" & trySuffix;
            pos      = variableNames.Insert (tryName);
            trySuffix ++;
        } while (pos < 0);
        DeleteObject(theV->theName);
        theV->theName = tryName;
    }

    if (pos < 0) {
        if (isDefiningATree == kTreeIsBeingParsed) {
            HandleApplicationError(_String("Error while creating a tree: duplicate node name ") & *theV->GetName()->Enquote());
            return;
        }
        theV->theIndex = variableNames.GetXtra(-pos-1);
        return;
    } else {
        theV->theName->AddAReference();
    }

    if (freeSlots.lLength) {
        theV->theIndex = freeSlots.list_data[freeSlots.lLength-1];
        variablePtrs[theV->theIndex]=theV->makeDynamic();
        freeSlots.Delete(freeSlots.lLength-1);
    } else {
        theV->theIndex = variablePtrs.lLength;
        variablePtrs&&theV;
    }
    variableNames.SetXtra (pos, theV->theIndex);
}

//__________________________________________________________________________________
_String const&  AppendContainerName (_String const& inString, _VariableContainer const* theP) {
    return AppendContainerName (inString, theP?theP->GetName():nil);
}

//__________________________________________________________________________________
_String const&  AppendContainerName (_String const& inString, _String const* namescp) {
    static _String returnMe;
    
    if (_hy_application_globals.Find (&inString) >= 0) {
        return inString;
    }
    
    hy_reference_type reference_type = inString.ProcessVariableReferenceCases (returnMe, namescp && !namescp -> empty() ? namescp : nil);
    
    if (reference_type != kStringInvalidReference) {
        return returnMe;
    }
    return inString;
}


//__________________________________________________________________________________
void  RenameVariable (_String* oldName, _String* newName)
{
    _String     oldNamePrefix (*oldName&'.'),
                newNamePrefix (*newName&'.');

    _List           toRename;
    _SimpleList     xtras,
                    traverser;

    long f = variableNames.Find (oldName, traverser);
    if (f>=0) {
        toRename << oldName;
        xtras    << variableNames.GetXtra (f);
        f = variableNames.Next (f, traverser);

        for  (; f>=0 && ((_String*)variableNames.Retrieve (f))->BeginsWith (oldNamePrefix); f = variableNames.Next (f, traverser)) {
            toRename << variableNames.Retrieve (f);
            xtras << variableNames.GetXtra (f);
        }
    }

    for (unsigned long k = 0; k < toRename.lLength; k++) {
        _Variable * thisVar = FetchVar (xtras.list_data[k]);
        thisVar->GetName()->RemoveAReference();
        if (k) {
            thisVar->theName = new _String(thisVar->GetName()->Replace(oldNamePrefix,newNamePrefix,true));
        } else {
            thisVar->theName = new _String(*newName);
        }

        variableNames.Delete (toRename (k), true);
        variableNames.Insert (thisVar->GetName(),xtras.list_data[k]);
        thisVar->GetName()->AddAReference();
    }
}

//__________________________________________________________________________________
void  ReplaceVar (_Variable* theV) {
    long pos = variableNames.Find (theV->theName);
    if (pos>=0) {
        pos = variableNames.GetXtra(pos);
        UpdateChangingFlas   (pos);
        variablePtrs.Replace (pos,theV,true);
    } else {
        InsertVar (theV);
    }
}


//__________________________________________________________________________________
void    SetupOperationLists (void) {

    auto package_ops = [] (long op1, long op2) -> long {
        if (op1 < op2) {
            return (op1 << 16) + op2;
        }
        return (op2 << 16) + op1;
    };
  
    UnOps  < "-" <
             "!" <
             "+" <
             "*" <
             "^" <
             "&" <
             "Abs" <
             "Sin" <
             "Cos" <
             "Tan" <
             "Exp" <
             "Log" <
             "Arctan" <
             "Time" <
             "Gamma" <
             "Transpose" <
             "Sqrt" <
             "Erf" <
             "Rows" <
             "Columns" <
             "LUDecompose" <
             "Inverse" <
             "BranchCount" <
             "TipCount" <
             "ZCDF" <
             "Eigensystem" <
             "Simplex" <
             "Type" <
             "Eval" <
             "LnGamma";
 

    BinOps<<'|'*256+'|';
    opPrecedence<<1;
    associativeOps << opPrecedence.lLength;

    BinOps<<'&'*256+'&';
    opPrecedence<<2;
    associativeOps << opPrecedence.lLength;

    BinOps<<'='*256+'=';
    opPrecedence<<3;
    BinOps<<'!'*256+'=';
    opPrecedence<<3;
  
    BinOps<<'<';
    opPrecedence<<4;
    BinOps<<'>';
    opPrecedence<<4;
    BinOps<<'<'*256+'=';
    opPrecedence<<4;
    BinOps<<'>'*256+'=';
    opPrecedence<<4;
  
    BinOps<<'+';
    associativeOps << opPrecedence.lLength;
    opPrecedence<<5;
    BinOps<<'-';
    opPrecedence<<5;
    BinOps<<'*';
    associativeOps << opPrecedence.lLength;
  
    opPrecedence<<6;
    BinOps<<'/';
    opPrecedence<<6;
    BinOps<<'%';
    opPrecedence<<6;
    BinOps<<'$';
    opPrecedence<<6;
  
    BinOps<<'^';
    opPrecedence<<7;
  
    BinOps<<'+'*256+'=';
    opPrecedence<<8;

    _Operation::ListOfInverseOps
    << package_ops(HY_OP_CODE_NOT, HY_OP_CODE_NOT)
    << package_ops (HY_OP_CODE_EXP, HY_OP_CODE_LOG)
    << package_ops (HY_OP_CODE_TAN, HY_OP_CODE_ARCTAN);

    if (BuiltInFunctions.empty()) {
        // construct a list of operations
        // don't forget to update SimplifyConstants, simpleOperationCodes, InternalDifferentiate, InternalSimplify, Formula::HasChanged and all Execute commands
        // also MAccess and MCoord codes are used in Parse to merge multiple matrix access operations

        //HY_OP_CODE_NOT
        BuiltInFunctions.AppendNewInstance (new _String ('!'));

        //HY_OP_CODE_NEQ
        BuiltInFunctions.AppendNewInstance (new _String ("!="));

        //HY_OP_CODE_IDIV
        BuiltInFunctions.AppendNewInstance (new _String ('$'));

        //HY_OP_CODE_MOD
        BuiltInFunctions.AppendNewInstance (new _String ('%'));

        //HY_OP_CODE_REF
        BuiltInFunctions.AppendNewInstance (new _String ('&'));
 
        //HY_OP_CODE_AND
        BuiltInFunctions.AppendNewInstance (new _String ("&&"));
        simpleOperationCodes    << HY_OP_CODE_AND;
        simpleOperationFunctions<< (long)AndNumbers;

        //HY_OP_CODE_MUL
        BuiltInFunctions.AppendNewInstance (new _String ('*'));
        simpleOperationCodes    << HY_OP_CODE_MUL;
        simpleOperationFunctions<< (long)MultNumbers;

        //HY_OP_CODE_ADD
        BuiltInFunctions.AppendNewInstance (new _String ('+'));
        simpleOperationCodes<< HY_OP_CODE_ADD;
        simpleOperationFunctions<<(long)AddNumbers;

        //HY_OP_CODE_SUB
        BuiltInFunctions.AppendNewInstance (new _String ('-'));
        simpleOperationCodes<<HY_OP_CODE_SUB;
        simpleOperationFunctions<<(long)SubNumbers;

        //HY_OP_CODE_DIV
        BuiltInFunctions.AppendNewInstance (new _String ('/'));
        simpleOperationCodes<<HY_OP_CODE_DIV;
        simpleOperationFunctions<<(long)DivNumbers;

        //HY_OP_CODE_LESS
        BuiltInFunctions.AppendNewInstance (new _String ('<'));
        simpleOperationCodes<<HY_OP_CODE_LESS;
        simpleOperationFunctions<<(long)LessThan;

        //HY_OP_CODE_LEQ
        BuiltInFunctions.AppendNewInstance (new _String ("<="));
        simpleOperationCodes<<HY_OP_CODE_LEQ;
        simpleOperationFunctions<<(long)LessThanE;

        //HY_OP_CODE_EQ
        BuiltInFunctions.AppendNewInstance (new _String ("=="));
        simpleOperationCodes<<HY_OP_CODE_EQ;
        simpleOperationFunctions<<(long)EqualNumbers;

        //HY_OP_CODE_GREATER
        BuiltInFunctions.AppendNewInstance (new _String ('>'));
        simpleOperationCodes<<HY_OP_CODE_GREATER;
        simpleOperationFunctions<<(long)GreaterThan;

        //HY_OP_CODE_GEQ
        BuiltInFunctions.AppendNewInstance (new _String (">="));
        simpleOperationCodes<<HY_OP_CODE_GEQ;
        simpleOperationFunctions<<(long)GreaterThanE;

        //HY_OP_CODE_ABS
        BuiltInFunctions.AppendNewInstance (new _String ("Abs"));
        simpleOperationCodes<<HY_OP_CODE_ABS;
        simpleOperationFunctions<<(long)AbsNumber;

        //HY_OP_CODE_ARCTAN
        BuiltInFunctions.AppendNewInstance (new _String ("Arctan"));

        //HY_OP_CODE_BETA
        BuiltInFunctions.AppendNewInstance (new _String ("Beta"));
        FunctionNameList.Insert (*(_String*)BuiltInFunctions (HY_OP_CODE_BETA), 2L);

        //HY_OP_CODE_BRANCHCOUNT
        BuiltInFunctions.AppendNewInstance (new _String ("BranchCount"));
        FunctionNameList.Insert (*(_String*)BuiltInFunctions (HY_OP_CODE_BRANCHCOUNT), 2L);

        //HY_OP_CODE_BRANCHLENGTH
        BuiltInFunctions.AppendNewInstance (new _String ("BranchLength"));
        FunctionNameList.Insert (*(_String*)BuiltInFunctions (HY_OP_CODE_BRANCHLENGTH), 2L);

        //HY_OP_CODE_BRANCHNAME
        BuiltInFunctions.AppendNewInstance (new _String ("BranchName"));
        FunctionNameList.Insert (*(_String*)BuiltInFunctions (HY_OP_CODE_BRANCHNAME), 2L);
      

        //HY_OP_CODE_CCHI2
        BuiltInFunctions.AppendNewInstance (new _String ("CChi2"));
        FunctionNameList.Insert (*(_String*)BuiltInFunctions (HY_OP_CODE_CCHI2), 2L);

        //HY_OP_CODE_CGAMMADIST
        BuiltInFunctions.AppendNewInstance (new _String ("CGammaDist"));
        FunctionNameList.Insert (*(_String*)BuiltInFunctions (HY_OP_CODE_CGAMMADIST), 3L);
 
        //HY_OP_CODE_CALL
        BuiltInFunctions.AppendNewInstance (new _String ("Call"));
        FunctionNameList.Insert (*(_String*)BuiltInFunctions (HY_OP_CODE_CALL), 1L + (0x7fff << 16));
            // mimimum of 1 argument

        //HY_OP_CODE_COLUMNS
        BuiltInFunctions.AppendNewInstance (new _String ("Columns"));

        //HY_OP_CODE_COS
        BuiltInFunctions.AppendNewInstance (new _String ("Cos"));

        //HY_OP_CODE_DIFF
        BuiltInFunctions.AppendNewInstance (new _String ("Differentiate"));
        FunctionNameList.Insert (*(_String*)BuiltInFunctions (HY_OP_CODE_DIFF), 2L);

        //HY_OP_CODE_EIGENSYSTEM
        BuiltInFunctions.AppendNewInstance (new _String ("Eigensystem"));

        //HY_OP_CODE_ERF
        BuiltInFunctions.AppendNewInstance (new _String ("Erf"));

        //HY_OP_CODE_EVAL
        BuiltInFunctions.AppendNewInstance (new _String ("Eval"));

        //HY_OP_CODE_EXP
        BuiltInFunctions.AppendNewInstance (new _String ("Exp"));
        simpleOperationCodes<<HY_OP_CODE_EXP;
        simpleOperationFunctions<<(long)ExpNumbers;

        //HY_OP_CODE_FORMAT
        BuiltInFunctions.AppendNewInstance (new _String ("Format"));
        FunctionNameList.Insert (*(_String*)BuiltInFunctions (HY_OP_CODE_FORMAT), 3L);

        //HY_OP_CODE_GAMMA
        BuiltInFunctions.AppendNewInstance (new _String ("Gamma"));

        //HY_OP_CODE_GAMMADIST
        BuiltInFunctions.AppendNewInstance (new _String ("GammaDist"));
        FunctionNameList.Insert (*(_String*)BuiltInFunctions (HY_OP_CODE_GAMMADIST), 3L);

        //HY_OP_CODE_IBETA
        BuiltInFunctions.AppendNewInstance (new _String ("IBeta"));
        FunctionNameList.Insert (*(_String*)BuiltInFunctions (HY_OP_CODE_IBETA), 3L);

        //HY_OP_CODE_IGAMMA
        BuiltInFunctions.AppendNewInstance (new _String ("IGamma"));
        FunctionNameList.Insert (*(_String*)BuiltInFunctions (HY_OP_CODE_IGAMMA), 2L);

        //HY_OP_CODE_INVCHI2
        BuiltInFunctions.AppendNewInstance (new _String ("InvChi2"));
        FunctionNameList.Insert (*(_String*)BuiltInFunctions (HY_OP_CODE_INVCHI2), 2L);

        //HY_OP_CODE_INVERSE
        BuiltInFunctions.AppendNewInstance (new _String ("Inverse"));

        //HY_OP_CODE_JOIN
        BuiltInFunctions.AppendNewInstance (new _String ("Join"));
        FunctionNameList.Insert (*(_String*)BuiltInFunctions (HY_OP_CODE_JOIN), 2L);

        //HY_OP_CODE_LUDECOMPOSE
        BuiltInFunctions.AppendNewInstance (new _String ("LUDecompose"));

        //HY_OP_CODE_LUSOLVE
        BuiltInFunctions.AppendNewInstance (new _String ("LUSolve"));
        FunctionNameList.Insert (*(_String*)BuiltInFunctions (HY_OP_CODE_LUSOLVE), 2L);

        //HY_OP_CODE_LOG_GAMMA
        BuiltInFunctions.AppendNewInstance (new _String ("LnGamma"));


        //HY_OP_CODE_LOG
        BuiltInFunctions.AppendNewInstance (new _String ("Log"));
        simpleOperationCodes<<HY_OP_CODE_LOG;
        simpleOperationFunctions<<(long)LogNumbers;

        //HY_OP_CODE_MACCESS
        BuiltInFunctions.AppendNewInstance (new _String ("MAccess"));
        simpleOperationCodes<<HY_OP_CODE_MACCESS;
        simpleOperationFunctions<<(long)FastMxAccess;

        //HY_OP_CODE_MCOORD
        BuiltInFunctions.AppendNewInstance (new _String ("MCoord"));
        simpleOperationCodes<<HY_OP_CODE_MCOORD;
        simpleOperationFunctions<<(long)FastMxWrite;

        //HY_OP_CODE_MAX
        BuiltInFunctions.AppendNewInstance (new _String ("Max"));
        FunctionNameList.Insert (*(_String*)BuiltInFunctions (HY_OP_CODE_MAX), 1L + ((2L) << 16)); // (1 or 2 arguments)
        
        simpleOperationCodes<<HY_OP_CODE_MAX;
        simpleOperationFunctions<<(long)MaxNumbers;

        //HY_OP_CODE_MIN
        BuiltInFunctions.AppendNewInstance (new _String ("Min"));
        FunctionNameList.Insert (*(_String*)BuiltInFunctions (HY_OP_CODE_MIN), 1L + ((2L) << 16)); // (1 or 2 arguments)

        simpleOperationCodes<<HY_OP_CODE_MIN;
        simpleOperationFunctions<<(long)MinNumbers;

        //HY_OP_CODE_PSTREESTRING
        BuiltInFunctions.AppendNewInstance (new _String ("PSTreeString"));
        FunctionNameList.Insert (*(_String*)BuiltInFunctions (HY_OP_CODE_PSTREESTRING), 3L);


        //HY_OP_CODE_RANDOM
        BuiltInFunctions.AppendNewInstance (new _String ("Random"));
        FunctionNameList.Insert (*(_String*)BuiltInFunctions (HY_OP_CODE_RANDOM), 2L);

        simpleOperationCodes<<HY_OP_CODE_RANDOM;
        simpleOperationFunctions<<(long)RandomNumber;

        //HY_OP_CODE_REROOTTREE
        BuiltInFunctions.AppendNewInstance (new _String ("RerootTree"));
        FunctionNameList.Insert (*(_String*)BuiltInFunctions (HY_OP_CODE_REROOTTREE), 2L);

        //HY_OP_CODE_ROWS
        BuiltInFunctions.AppendNewInstance (new _String ("Rows"));

        //HY_OP_CODE_SIMPLEX
        BuiltInFunctions.AppendNewInstance (new _String ("Simplex"));

        //HY_OP_CODE_EXPRESSION
        BuiltInFunctions.AppendNewInstance (new _String ("Simplify"));
        FunctionNameList.Insert (*(_String*)BuiltInFunctions (HY_OP_CODE_SIMPLIFY), 2L);

        //HY_OP_CODE_SIN
        BuiltInFunctions.AppendNewInstance (new _String ("Sin"));

        //HY_OP_CODE_SQRT
        BuiltInFunctions.AppendNewInstance (new _String ("Sqrt"));

        //HY_OP_CODE_TEXTREESTRING
        BuiltInFunctions.AppendNewInstance (new _String ("TEXTreeString"));
        FunctionNameList.Insert (*(_String*)BuiltInFunctions (HY_OP_CODE_TEXTREESTRING), 2L);

        //HY_OP_CODE_TAN
        BuiltInFunctions.AppendNewInstance (new _String ("Tan"));

        //HY_OP_CODE_TIME
        BuiltInFunctions.AppendNewInstance (new _String ("Time"));

        //HY_OP_CODE_TIPCOUNT
        BuiltInFunctions.AppendNewInstance (new _String ("TipCount"));

        //HY_OP_CODE_TIPNAME
        BuiltInFunctions.AppendNewInstance (new _String ("TipName"));
        FunctionNameList.Insert (*(_String*)BuiltInFunctions (HY_OP_CODE_TIPNAME), 2L);

        //HY_OP_CODE_TRANSPOSE
        BuiltInFunctions.AppendNewInstance (new _String ("Transpose"));

        //HY_OP_CODE_TYPE
        BuiltInFunctions.AppendNewInstance (new _String ("Type"));

        //HY_OP_CODE_ZCDF
        BuiltInFunctions.AppendNewInstance (new _String ("ZCDF"));

        //HY_OP_CODE_POWER
        BuiltInFunctions.AppendNewInstance (new _String ('^'));
        simpleOperationCodes<<HY_OP_CODE_POWER;
        simpleOperationFunctions<<(long)Power;

        //HY_OP_CODE_OR
        BuiltInFunctions.AppendNewInstance (new _String ("||"));

        hyReservedWords << BuiltInFunctions;
        hyReservedWords.AppendNewInstance (new _String("global"));
      
      
        hyReservedWords.Sort();
        /*hyReservedWords.ForEach ([] (BaseRef item, unsigned long index) -> void {
            printf ("%d %s\n", index, ((_String*)item)->get_str());
        });*/
    }



}



//__________________________________________________________________________________

void  FinishDeferredSF (void) {
    if (deferSetFormula->nonempty()) {
        SortLists (deferSetFormula, &deferIsConstant);
        
        for (AVLListXIteratorKeyValue variable_record : AVLListXIterator (&variableNames)) {
            _Variable * theV = LocateVar(variable_record.get_value());
            if (theV->IsContainer()) {
                ((_VariableContainer*)theV)->SetMDependance (*deferSetFormula);
            }
        }
        
        likeFuncList.ForEach([] (BaseRef lf_object, unsigned long idx) -> void {
            if (((_String*)likeFuncNamesList(idx))->nonempty()) {
                _LikelihoodFunction * lf = (_LikelihoodFunction*)lf_object;
                for (long k = 0L; k < deferSetFormula->countitems(); k++) {
                    lf->UpdateIndependent(deferSetFormula->get(k),deferIsConstant.get(k));
                }
            }
        });
        

       
    }
    DeleteObject (deferSetFormula);
    deferSetFormula = nil;
    deferIsConstant.Clear();
}
