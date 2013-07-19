/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2009
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon              (apoon@cfenet.ubc.ca)

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

#include "likefunc.h"
#include "legacy_parser.h"
#include "matrix.h"
#include "stdlib.h"
#include "string.h"

#include "trieiterator.h"

#include <stdio.h>
#include "time.h"

#include "ctype.h"
#include "polynoml.h"
#include "batchlan.h"
#include "hy_globals.h"
#include "executionlist.h"

extern _SimpleList simpleOperationCodes, simpleOperationFunctions;

// stores all the variables declared so far
_List globalNamesSupportList, hyReservedWords, varNamesSupportList,
    variablePtrs; 
    
_SimpleList   _builtInArgumentCountsAux;


_AVLListX variableNames   (&varNamesSupportList),
    _hyApplicationGlobals (&globalNamesSupportList),
    _builtInArgumentCounts (&_builtInArgumentCountsAux);       

_Parameter printDigits;

// indices of all independent variables
_SimpleList deferIsConstant, 
            *deferSetFormula = nil;
            
bool useGlobalUpdateFlag = false;
_String HalfOps(":<>=!&|");
_Trie UnOps, BuiltInFunctions;
_SimpleList opPrecedence, BinOps, associativeOps;
_Parameter pi_const = 3.141592653589793, long_max = (_Parameter) LONG_MAX;

/**********************************/
/* Defining Globals here for now */
/********************************/

//Used in formula, and constant

#ifndef __HYALTIVEC__
_Parameter machineEps = 2. * DBL_EPSILON, tolerance = DBL_EPSILON;
#else
_Parameter machineEps = 1e-7, tolerance = FLT_EPSILON;
#endif

//Used in formula
_String intPrecFact("INTEGRATION_PRECISION_FACTOR"),
    intMaxIter("INTEGRATION_MAX_ITERATES");

//Used in parser2 and formula
_Parameter sqrtPi = 1.77245385090551603;
_Parameter twoOverSqrtPi = 2. / sqrtPi;

/*********************************/
/*          End Globals         */
/*******************************/

//______________________________________________________________________________
//Helper functions
//______________________________________________________________________________
void DeleteTreeVariable(long, _SimpleList &, bool);

//______________________________________________________________________________
_Variable *LocateVar(long index) {
  return (_Variable *)(((BaseRef *)variablePtrs.lData)[index]);
}

//______________________________________________________________________________
void parameterToCharBuffer(_Parameter value, char *dump, long length) {
  long digs = printDigits;
  if (digs <= 0 || digs > 15) {
    if (round(value) == value && fabs(value) < long_max) {
      snprintf(dump, length, "%ld", lrint(value));
    } else {
      snprintf(dump, length, PRINTF_FORMAT_STRING, value);
    }
  } else {
    _String format("%-");
#ifdef __USE_LONG_DOUBLE__
    format = format & _String(digs) & "Lg";
#else
    format = format & _String(digs) & 'g';
#endif
    snprintf(dump, length, (const char *)format.sData, value);
  }
}

//______________________________________________________________________________
BaseRef parameterToString(_Parameter value) {
  char dump[256];
  parameterToCharBuffer(value, dump, 256);
  return new _String(dump);
}

//______________________________________________________________________________
void SplitVariableIDsIntoLocalAndGlobal(const _SimpleList &theList,
                                        _List &splitStorage) {
  splitStorage.Clear();
  splitStorage.AppendNewInstance(new _SimpleList);
  splitStorage.AppendNewInstance(new _SimpleList);

  for (unsigned long k = 0; k < theList.lLength; k++) {
    long varID = theList.lData[k];
    (*(_SimpleList *)splitStorage(1 - LocateVar(varID)->IsGlobal())) << varID;
  }
}

//______________________________________________________________________________
_String FetchObjectNameFromType(const unsigned long objectClass) {
  switch (objectClass) {
  case HY_UNDEFINED:
    return "Undefined";
  case NUMBER:
    return "Number";
  case MATRIX:
    return "Container variable";
  case TREE_NODE:
    return "Tree node";
  case TREE:
    return "Tree";
  case STRING:
    return "String";
  case ASSOCIATIVE_LIST:
    return "Associative Array";
  case TOPOLOGY:
    return "Topology";
  case POLYNOMIAL:
    return "Polynomial";
  case HY_ANY_OBJECT:
    return "Any HyPhy object";
  }
  return empty;
}

//______________________________________________________________________________
_String *FetchMathObjectNameOfTypeByIndex(const unsigned long objectClass,
                                          const long objectIndex) {

  if (objectIndex >= 0 && objectIndex < variableNames.countitems()) {
    long tc = 0;
    _SimpleList nts;
    long rt, vi = variableNames.Traverser(nts, rt, variableNames.GetRoot());

    for (; vi >= 0; vi = variableNames.Traverser(nts, rt))
      if (FetchVar(variableNames.GetXtra(vi))->ObjectClass() == objectClass) {
        if (tc == objectIndex) {
          return (_String *)variableNames.Retrieve(vi);
          break;
        } else {
          tc++;
        }
      }
  }
  return nil;
}

//______________________________________________________________________________
_PMathObj FetchObjectFromVariableByType(_String *id,
                                        const unsigned long objectClass,
                                        long command_id, _String *errMsg) {
  if (id) {
    _Variable *v = FetchVar(LocateVarByName(*id));
    if (v &&
        (objectClass == HY_ANY_OBJECT || v->ObjectClass() == objectClass)) {
      return v->Compute();
    }
    if (command_id >= 0 || errMsg) {
      if (command_id >= 0) {
        WarnError(_String("'") & *id & ("' must refer to a ") &
                  FetchObjectNameFromType(objectClass) & " in call to " &
                  _HY_ValidHBLExpressions.RetrieveKeyByPayload(command_id) &
                  '.');
      } else {
        WarnError(errMsg->Replace("_VAR_NAME_ID_", *id, true));
      }
    }
  }
  return nil;
}

//______________________________________________________________________________
_PMathObj FetchObjectFromVariableByTypeIndex(long idx,
                                             const unsigned long objectClass,
                                             long command_id, _String *errMsg) {
  _Variable *v = FetchVar(idx);

  if (v && (objectClass == HY_ANY_OBJECT || v->ObjectClass() == objectClass)) {
    return v->GetValue();
  }

  if (command_id >= 0 || errMsg) {
    if (command_id >= 0) {
      WarnError(_String("'") & *v->GetName() & ("' must refer to a ") &
                FetchObjectNameFromType(objectClass) & " in call to " &
                _HY_ValidHBLExpressions.RetrieveKeyByPayload(command_id) & '.');
    } else {
      WarnError(errMsg->Replace("_VAR_NAME_ID_", *v->GetName(), true));
    }
  }
  return nil;
}

//______________________________________________________________________________
long LocateVarByName(_String &name) { return variableNames.Find(&name); }

//______________________________________________________________________________
_Variable *FetchVar(long index) {
  return index >= 0 ? (_Variable *)variablePtrs(variableNames.GetXtra(index))
                    : nil;
}

//______________________________________________________________________________
void UpdateChangingFlas(long vN) {
  // check to see if formulae contain a reference to this var
  // if so: "decompile" them

  long topLimit = compiledFormulaeParameters.lLength;

  _SimpleList *toDelete = nil;

  for (long k = 0; k < topLimit; k++) {
    long g =
        ((_SimpleList *)compiledFormulaeParameters.lData[k])->BinaryFind(vN, 0);

    if (g >= 0) {
      _ElementaryCommand *thisCommand =
          (_ElementaryCommand *)listOfCompiledFormulae.lData[k];
      _Formula *f = (_Formula *)(thisCommand->simpleParameters.lData[1]),
               *f2 = (_Formula *)(thisCommand->simpleParameters.lData[2]);

      delete f;
      delete f2;

      thisCommand->simpleParameters.Clear();

      // MOD 10/21/2005
      if (!toDelete) {
        checkPointer(toDelete = new _SimpleList);
      }

      //listOfCompiledFormulae.    Delete(k);
      //compiledFormulaeParameters.Delete(k);

      //k--;
      //topLimit--;

      *toDelete << k;
    }
  }

  if (toDelete) {
    listOfCompiledFormulae.DeleteList(*toDelete);
    compiledFormulaeParameters.DeleteList(*toDelete);
    DeleteObject(toDelete);
  }
}

//______________________________________________________________________________
void UpdateChangingFlas(_SimpleList &involvedVariables) {

  long topLimit = compiledFormulaeParameters.lLength;
  _SimpleList *toDelete = nil;

  for (long k = 0; k < topLimit; k++) {
    long g = ((_SimpleList *)compiledFormulaeParameters.lData[k])
        ->CountCommonElements(involvedVariables, true);

    if (g > 0) {
      _ElementaryCommand *thisCommand =
          (_ElementaryCommand *)listOfCompiledFormulae.lData[k];

      _Formula *f = (_Formula *)(thisCommand->simpleParameters.lData[1]),
               *f2 = (_Formula *)(thisCommand->simpleParameters.lData[2]);

      delete f;
      delete f2;

      thisCommand->simpleParameters.Clear();

      if (!toDelete) {
        checkPointer(toDelete = new _SimpleList);
      }

      *toDelete << k;
    }
  }

  if (toDelete) {
    listOfCompiledFormulae.DeleteList(*toDelete);
    compiledFormulaeParameters.DeleteList(*toDelete);
    DeleteObject(toDelete);
  }
}

//______________________________________________________________________________
void DeleteVariable(long dv, bool deleteself) {
  if (dv >= 0) {

    _String *name = (_String *)variableNames.Retrieve(dv);
    _String myName = *name & '.';
    long vidx = variableNames.GetXtra(dv);

    UpdateChangingFlas(vidx);

    _SimpleList recCache;
    variableNames.Find(name, recCache);
    _String
    nextVarID; // = *(_String*)variableNames.Retrieve(variableNames.Next
               // (dv,recCache));
    long nvid;
    if ((nvid = variableNames.Next(dv, recCache)) >= 0) {
      nextVarID = *(_String *)variableNames.Retrieve(nvid);
    }

    if (deleteself) {
      _SimpleList tcache;
      long iv, k = variableNames.Traverser(tcache, iv, variableNames.GetRoot());

      for (; k >= 0; k = variableNames.Traverser(tcache, iv)) {
        _Variable *thisVar = FetchVar(k);

        if (thisVar->CheckFForDependence(vidx, false)) {
          _PMathObj curValue = thisVar->Compute();
          curValue->nInstances++; // this could be a leak 01/05/2004.
          thisVar->SetValue(curValue);
          DeleteObject(curValue);
        }
      }

      _Variable *delvar = (FetchVar(dv));
      DeleteObject(delvar);

      variableNames.Delete(variableNames.Retrieve(dv), true);
      (*((_SimpleList *)&variablePtrs))[vidx] = 0;
      freeSlots << vidx;
    } else {
      _Variable *delvar = (FetchVar(dv));
      if (delvar->IsContainer()) {
        _VariableContainer *dc = (_VariableContainer *)delvar;
        dc->Clear();
      }
    }

    _List toDelete;

    recCache.Clear();
    long nextVar = variableNames.Find(&nextVarID, recCache);

    for (; nextVar >= 0; nextVar = variableNames.Next(nextVar, recCache)) {
      _String dependent = *(_String *)variableNames.Retrieve(nextVar);
      if (dependent.startswith(myName)) {
        toDelete &&&dependent;
      } else {
        break;
      }
    }

    for (unsigned long k = 0; k < toDelete.lLength; k++) {
      DeleteVariable(*(_String *)toDelete(k));
    }
  }
}

//______________________________________________________________________________
void DeleteTreeVariable(long dv, _SimpleList &parms, bool doDeps) {
  if (dv >= 0) {
    _String *name = (_String *)variableNames.Retrieve(dv);
    _String myName = *name & ".";
    long vidx = variableNames.GetXtra(dv);

    UpdateChangingFlas(vidx);

    _SimpleList recCache;
    variableNames.Find(name, recCache);
    _String nextVarID;
    long nvid;
    if ((nvid = variableNames.Next(dv, recCache)) >= 0) {
      nextVarID = *(_String *)variableNames.Retrieve(nvid);
    }

    {
      _SimpleList tcache;
      long iv, k = variableNames.Traverser(tcache, iv, variableNames.GetRoot());

      for (; k >= 0; k = variableNames.Traverser(tcache, iv)) {
        _Variable *thisVar = FetchVar(k);

        if (thisVar->CheckFForDependence(vidx, false)) {
          _PMathObj curValue = thisVar->Compute();
          curValue->nInstances++;
          thisVar->SetValue(curValue);
          DeleteObject(curValue);
        }
      }
    }

    _Variable *delvar = (FetchVar(dv));
    if (delvar->ObjectClass() != TREE) {
      variableNames.Delete(variableNames.Retrieve(dv), true);
      (*((_SimpleList *)&variablePtrs))[vidx] = 0;
      freeSlots << vidx;
      DeleteObject(delvar);
    } else {
      ((_VariableContainer *)delvar)->Clear();
    }
    if (doDeps) {
      _List toDelete;
      recCache.Clear();
      long nextVar = variableNames.Find(&nextVarID, recCache);
      for (; nextVar >= 0; nextVar = variableNames.Next(nextVar, recCache)) {
        _String dependent = *(_String *)variableNames.Retrieve(nextVar);
        if (dependent.startswith(myName)) {
          if (dependent.Find('.', myName.sLength + 1, -1) >= 0) {
            _Variable *checkDep = FetchVar(nextVar);
            if (!checkDep->IsIndependent()) {
              _PMathObj curValue = checkDep->Compute();
              curValue->nInstances++;
              checkDep->SetValue(curValue);
              DeleteObject(curValue);
            }
            parms << variableNames.GetXtra(nextVar);
          } else {
            toDelete &&&dependent;
          }
        } else {
          break;
        }
      }

      for (unsigned long k = 0; k < toDelete.lLength; k++) {
        //StringToConsole (*(_String*)toDelete(k));
        //BufferToConsole ("\n");
        DeleteTreeVariable(*(_String *)toDelete(k), parms, false);
      }
    }
  }
}

//______________________________________________________________________________
void DeleteVariable(_String &name, bool deleteself) {
  DeleteVariable(LocateVarByName(name), deleteself);
}

//______________________________________________________________________________
void DeleteTreeVariable(_String &name, _SimpleList &parms, bool doDeps) {
  DeleteTreeVariable(LocateVarByName(name), parms, doDeps);
}

//______________________________________________________________________________
_Variable *CheckReceptacle(_String *name, _String fID, bool checkValid,
                           bool isGlobal) {
  if (checkValid && (!name->IsValidIdentifier())) {
    _String errMsg =
        *name & " is not a valid variable identifier in call to " & fID;
    WarnError(errMsg);
    return nil;
  }

  long f = LocateVarByName(*name);
  if (f < 0) {
    _Variable dummy(*name, isGlobal);
    f = LocateVarByName(*name);
  }

  return FetchVar(f);
}

//______________________________________________________________________________
_Variable *CheckReceptacleCommandID(_String *name, const long id,
                                    bool checkValid, bool isGlobal,
                                    _ExecutionList *context) {
  if (checkValid && (!name->IsValidIdentifier())) {
    _String errMsg = _String("'") & *name &
                     "' is not a valid variable identifier in call to " &
                     _HY_ValidHBLExpressions.RetrieveKeyByPayload(id) & '.';
    if (context) {
      context->ReportAnExecutionError(errMsg);
    } else {
      WarnError(errMsg);
    }
    return nil;
  }

  long f = LocateVarByName(*name);
  if (f < 0) {
    _Variable dummy(*name, isGlobal);
    f = LocateVarByName(*name);
  }

  return FetchVar(f);
}

//______________________________________________________________________________
bool CheckReceptacleCommandIDAndStore(_String *name, const long id,
                                      bool checkValid, _PMathObj v, bool dup,
                                      bool isGlobal) {

  _Variable *theV = CheckReceptacleCommandID(name, id, checkValid, isGlobal);
  if (theV) {
    theV->SetValue(v, dup);
    return true;
  }
  if (!dup) {
    DeleteObject(v);
  }
  return false;

}

//______________________________________________________________________________
bool CheckReceptacleAndStore(_String *name, _String fID, bool checkValid,
                             _PMathObj v, bool dup) {

  _Variable *theV = CheckReceptacle(name, fID, checkValid);
  if (theV) {
    theV->SetValue(v, dup);
    return true;
  }
  if (!dup) {
    DeleteObject(v);
  }
  return false;

}

//______________________________________________________________________________
bool CheckReceptacleAndStore(_String name, _String fID, bool checkValid,
                             _PMathObj v, bool dup) {

  return CheckReceptacleAndStore(&name, fID, checkValid, v, dup);

}

//______________________________________________________________________________
void InsertVar(_Variable *theV) {
  long pos = variableNames.Insert(theV->theName);

  /*if (theV->GetName()->Equal (&_String("PS_2")))
  {
      printf ("Making...\n");
  }*/

  if (pos < 0 && isDefiningATree > 1)
      // automatically fix duplicate autogenerated tree node name
      {
    long trySuffix = 1;
    _String *tryName = new _String;
    do {
      *tryName = *theV->theName & "_" & trySuffix;
      pos = variableNames.Insert(tryName);
      trySuffix++;
    } while (pos < 0);
    DeleteObject(theV->theName);
    theV->theName = tryName;
  }

  if (pos < 0) {
    if (isDefiningATree == 1) {
      WarnError(_String("Error while creating a tree: duplicate node name '") &
                *theV->GetName() & "'");
      return;
    }

    theV->theIndex = variableNames.GetXtra(-pos - 1);
    return;
  } else {
    theV->theName->nInstances++;
  }

  if (freeSlots.lLength) {
    theV->theIndex = freeSlots.lData[freeSlots.lLength - 1];
    variablePtrs[theV->theIndex] = theV->makeDynamic();
    freeSlots.Delete(freeSlots.lLength - 1);
  } else {
    theV->theIndex = variablePtrs.lLength;
    variablePtrs &&theV;
  }
  variableNames.SetXtra(pos, theV->theIndex);
}

//______________________________________________________________________________
_String &AppendContainerName(_String &inString, _VariableContainer *theP) {
  return AppendContainerName(inString, theP ? theP->GetName() : nil);
}

//______________________________________________________________________________
_String &AppendContainerName(_String &inString, _String *namescp) {
  static _String returnMe;

  if (_hyApplicationGlobals.Find(&inString) >= 0) {
    return inString;
  }

  unsigned char reference_type = inString.ProcessVariableReferenceCases(
      returnMe, namescp && namescp->sLength ? namescp : nil);

  if (reference_type != HY_STRING_INVALID_REFERENCE) {
    return returnMe;
  }
  return inString;
}

//______________________________________________________________________________
void RenameVariable(_String *oldName, _String *newName) {
  _String oldNamePrefix(*oldName & '.'), newNamePrefix(*newName & '.');

  _List toRename;
  _SimpleList xtras, traverser;

  long f = variableNames.Find(oldName, traverser);
  if (f >= 0) {
    toRename << oldName;
    xtras << variableNames.GetXtra(f);
    f = variableNames.Next(f, traverser);

    for (; f >= 0 && ((_String *)variableNames.Retrieve(f))
                         ->startswith(oldNamePrefix);
         f = variableNames.Next(f, traverser)) {
      toRename << variableNames.Retrieve(f);
      xtras << variableNames.GetXtra(f);
    }
  }

  for (unsigned long k = 0; k < toRename.lLength; k++) {
    _Variable *thisVar = FetchVar(xtras.lData[k]);
    thisVar->GetName()->RemoveAReference();
    if (k) {
      thisVar->theName = new _String(
          thisVar->GetName()->Replace(oldNamePrefix, newNamePrefix, true));
    } else {
      thisVar->theName = new _String(*newName);
    }

    variableNames.Delete(toRename(k), true);
    variableNames.Insert(thisVar->GetName(), xtras.lData[k]);
    thisVar->GetName()->nInstances++;
  }
}

//______________________________________________________________________________
void ReplaceVar(_Variable *theV) {
  long pos = variableNames.Find(theV->theName);
  if (pos >= 0) {
    pos = variableNames.GetXtra(pos);
    UpdateChangingFlas(pos);
    variablePtrs.Replace(pos, theV, true);
  } else {
    InsertVar(theV);
  }
}

//______________________________________________________________________________

void  InsertOpDescription (const char* str, const long opcode, void* simpleOp, long arguments) {
  BuiltInFunctions.Insert(str, opcode);
  if (simpleOp) {
    simpleOperationCodes << opcode;
    simpleOperationFunctions << (long)simpleOp;
  }
  hyReservedWords.AppendNewInstance(new _String (str));
  if (arguments >= 0) {
    _builtInArgumentCounts.Insert((BaseRef)opcode, arguments);
  }
}

//______________________________________________________________________________
void SetupOperationLists(void) {

  _List all_unary_ops(
      "-", 29, "!", "+", "*", "^", "&", "Abs", "Sin", "Cos", "Tan", "Exp",
      "Log", "Arctan", "Time", "Gamma", "Transpose", "Sqrt", "Erf", "Rows",
      "Columns", "LUDecompose", "Inverse", "BranchCount", "TipCount", "ZCDF",
      "Eigensystem", "Simplex", "Type", "Eval", "LnGamma");
      //_all_binary_ops (
      
      //);

  UnOps.Insert(all_unary_ops);

  //BinOps<<'=';
  //opPrecedence<<1;

  BinOps << '|' * 256 + '|';
  opPrecedence << 2;

  BinOps << '&' * 256 + '&';
  opPrecedence << 3;

  BinOps << '=' * 256 + '=';
  opPrecedence << 4;
  BinOps << '!' * 256 + '=';
  opPrecedence << 4;
  BinOps << '<';
  opPrecedence << 4;
  BinOps << '>';
  opPrecedence << 4;
  BinOps << '<' * 256 + '=';
  opPrecedence << 4;
  BinOps << '>' * 256 + '=';
  opPrecedence << 4;

  BinOps << '+';
  associativeOps << opPrecedence.lLength;
  opPrecedence << 5;
  BinOps << '-';
  opPrecedence << 5;
  BinOps << '*';
  associativeOps << opPrecedence.lLength;

  opPrecedence << 6;
  BinOps << '/';
  opPrecedence << 6;
  BinOps << '%';
  opPrecedence << 6;
  BinOps << '$';
  opPrecedence << 6;

  BinOps << '^';
  opPrecedence << 7;

  BinOps << '+' * 256 + '=';
  opPrecedence << 8;

  if (BuiltInFunctions.countitems() == 0)
      // construct a list of operations
      // don't forget to update SimplifyConstants, simpleOperationCodes,
      // InternalDifferentiate, InternalSimplify, Formula::HasChanged and all
      // Execute commands
      // also MAccess and MCoord codes are used in Parse to merge multiple
      // matrix access operations
      {
    

    for (long k = 0; k < HY_OP_COUNT; k++) {
      InsertOpDescription (_hyBuiltInOperationProperties[k].op_string,
                           _hyBuiltInOperationProperties[k].op_code,
                           _hyBuiltInOperationProperties[k].fast_exec,
                           _hyBuiltInOperationProperties[k].argument_count);
    }

    hyReservedWords.AppendNewInstance(new _String("global"));
    hyReservedWords.Sort();
  }
  
  _TrieIterator ti (BuiltInFunctions);
  _String *key = ti.Last();
  while (key) {
    printf ("%s (%ld)\n", key->sData, BuiltInFunctions.GetValue (ti.CurrentIndex()));
    DeleteObject(key);
    key = ti.Previous ();
  }
  DeleteObject(key);
}

//______________________________________________________________________________
void CompileListOfUserExpressions(_SimpleList &varRefs, _List &rec,
                                  bool doAll) {
  rec.Clear();
  if (varRefs.lLength == 0) {
    return;
  }

  long i;
  _SimpleList startVars;
  _VariableContainer *firstVar =
      (_VariableContainer *)LocateVar(varRefs.lData[0]);

  firstVar->ScanAndAttachVariables();

  {
    _AVLList sA(&startVars);
    if (doAll) {

      firstVar->ScanForVariables(sA, sA);
      firstVar->ScanForGVariables(sA, sA);
    }

    firstVar->ScanForDVariables(sA, sA);
    sA.ReorderList();
  }

  if (!doAll) {
    for (i = startVars.lLength - 1; i >= 0; i--) {
      if (firstVar->IsModelVar(i)) {
        startVars.Delete(i);
      }
    }
  }

  for (i = 0; i < startVars.lLength; i++) {
    _String thisName(LocateVar(startVars.lData[i])->GetName()->Cut(
        LocateVar(startVars.lData[i])->GetName()->FindBackwards('.', 0, -1),
        -1));
    rec &&&thisName;
  }

  for (i = varRefs.lLength - 1; i >= 1; i--) {
    firstVar = (_VariableContainer *)LocateVar(varRefs.lData[i]);
    firstVar->ScanAndAttachVariables();
    firstVar->MatchParametersToList(rec, doAll);
  }

  for (i = rec.lLength - 1; i >= 0; i--) {
    _String *thisLine = ((_String *)rec(i));
    thisLine->Trim(1, -1);
    if (doAll)
      if (LocateVarByName(*thisLine) < 0) {
        *thisLine = _String('!') & *thisLine;
      }
  }

}

//______________________________________________________________________________
void FindUnusedObjectName(_String &prefix, _String &partName, _List &names,
                          bool sorted) {
  if (partName.sLength == 0) {
    partName = prefix;
  }

  _String tryName(partName);
  long k = 1;

  if (sorted)
    while (names.BinaryFind(&tryName) >= 0) {
      k++;
      tryName = partName & k;
    }
  else
    while (names.Find(&tryName) >= 0) {
      k++;
      tryName = partName & k;
    }

  partName = tryName;
}

//______________________________________________________________________________
void FindUnusedObjectName(_String &prefix, _String &partName, _AVLListX &names,
                          bool) {
  if (partName.sLength == 0) {
    partName = prefix;
  }

  _String tryName(partName);
  long k = 1;

  while (names.Find(&tryName) >= 0) {
    k++;
    tryName = partName & k;
  }

  partName = tryName;
}


//______________________________________________________________________________
void FinishDeferredSF(void) {
  if (deferSetFormula->lLength) {
    SortLists(deferSetFormula, &deferIsConstant);
    _SimpleList tcache;
    long iv, i = variableNames.Traverser(tcache, iv, variableNames.GetRoot());

    for (; i >= 0; i = variableNames.Traverser(tcache, iv)) {
      _Variable *theV = FetchVar(i);
      if (theV->IsContainer()) {
        ((_VariableContainer *)theV)->SetMDependance(*deferSetFormula);
      }
    }

    for (long j = 0; j < likeFuncList.lLength; j++)
      if (((_String *)likeFuncNamesList(j))->sLength) {
        _LikelihoodFunction *lf = (_LikelihoodFunction *)likeFuncList(j);
        for (long k = 0; k < deferSetFormula->lLength; k++) {
          lf->UpdateIndependent(deferSetFormula->lData[k],
                                deferIsConstant.lData[k]);
        }
      }
  }
  DeleteObject(deferSetFormula);
  deferSetFormula = nil;
  deferIsConstant.Clear();
}
