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

#include <math.h>
#include <float.h>
#include <limits.h>

#include "likefunc.h"
#include "parser.h"
#include "matrix.h"
#include "stdlib.h"
#include "string.h"

#include <stdio.h>
#include "time.h"

#include "ctype.h"
#include "polynoml.h"
#include "batchlan.h"


#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

extern
_SimpleList     simpleOperationCodes,
                simpleOperationFunctions;

_List           globalNamesSupportList,
                hyReservedWords,
                varNamesSupportList,
                variablePtrs;   // stores all the variables declared so far

_AVLList        *lookAside = nil;

_AVLListX       variableNames (&varNamesSupportList),
                _hyApplicationGlobals (&globalNamesSupportList);

_Parameter printDigits;



// indices of all independent variables


_List           FunctionNameList,
                BuiltInFunctions;

_SimpleList     FunctionArgumentCount,
                freeSlots,
                deferIsConstant,
                *deferSetFormula = nil;

bool            useGlobalUpdateFlag = false;


_String         HalfOps (":<>=!&|");

_Trie           UnOps;

_SimpleList opPrecedence,
            BinOps,
            associativeOps;

_Parameter pi_const = 3.141592653589793,
long_max = (_Parameter)LONG_MAX;

/**********************************/
/* Defining Globals here for now */
/********************************/
 
//Used in formula, and constant

#ifndef  __HYALTIVEC__
_Parameter  machineEps = 2.*DBL_EPSILON,
tolerance  = DBL_EPSILON;
#else
_Parameter  machineEps = 1e-7,
tolerance  = FLT_EPSILON;
#endif

//Used in formula
_String         intPrecFact ("INTEGRATION_PRECISION_FACTOR"),
                intMaxIter  ("INTEGRATION_MAX_ITERATES");


//Used in parser2 and formula
_Parameter sqrtPi = 1.77245385090551603;
_Parameter twoOverSqrtPi = 2./sqrtPi;

/*********************************/
/*          End Globals         */
/*******************************/

//__________________________________________________________________________________
//SW: Helper functions

//__________________________________________________________________________________
void            DeleteTreeVariable      (long, _SimpleList &, bool);

//__________________________________________________________________________________
_Variable * LocateVar (long index)
{
    return (_Variable *)(((BaseRef*)variablePtrs.lData)[index]);
}

//__________________________________________________________________________________
void     parameterToCharBuffer (_Parameter value, char* dump, long length, bool json)
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
  
    long digs = printDigits;
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
        snprintf (dump,length,(const char*)format.sData,value);
    }
}


//__________________________________________________________________________________
BaseRef     parameterToString (_Parameter value)
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
        long varID = theList.lData[k];
        (*(_SimpleList*)splitStorage(1-LocateVar (varID)->IsGlobal())) << varID;
    }
}

//__________________________________________________________________________________
_String FetchObjectNameFromType (const unsigned long objectClass) {
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
_PMathObj   FetchObjectFromVariableByType (_String* id, const unsigned long objectClass, long command_id, _String *errMsg)
{
    if (id) {
        _Variable * v = FetchVar (LocateVarByName (*id));
        if (v && (objectClass == HY_ANY_OBJECT || v->ObjectClass () == objectClass)) {
            return v->Compute();
        }
        if (command_id >= 0 || errMsg) {
            if (command_id >= 0) {
                WarnError (_String ("'") & *id & ("' must refer to a ") & FetchObjectNameFromType (objectClass) & " in call to " 
                                         &_HY_ValidHBLExpressions.RetrieveKeyByPayload(command_id) & '.');
            } else {
                WarnError (errMsg->Replace ("_VAR_NAME_ID_", *id, true));
            }
        }
    }
    return nil;
}


//__________________________________________________________________________________
_PMathObj   FetchObjectFromVariableByTypeIndex (long idx, const unsigned long objectClass, long command_id, _String *errMsg)
{
    _Variable * v = FetchVar (idx);
    if (v && (objectClass == HY_ANY_OBJECT || v->ObjectClass () == objectClass)) {
        return v->GetValue();
    }
    if (command_id >= 0 || errMsg) {
        if (command_id >= 0) {
            WarnError (_String ("'") & *v->GetName() & ("' must refer to a ") & FetchObjectNameFromType (objectClass) & " in call to " 
                                     &_HY_ValidHBLExpressions.RetrieveKeyByPayload(command_id) & '.');
        } else {
            WarnError (errMsg->Replace ("_VAR_NAME_ID_", *v->GetName(), true));
        }
    }    
    return nil;
}

//__________________________________________________________________________________
long LocateVarByName (_String& name)
{
    return variableNames.Find (&name);
}

//__________________________________________________________________________________
_Variable* FetchVar (long index)
{
    return index>=0?(_Variable *)variablePtrs(variableNames.GetXtra(index)):nil;
}

//__________________________________________________________________________________
void       UpdateChangingFlas (long vN)
{
    // check to see if formulae contain a reference to this var
    // if so: "decompile" them

    long topLimit = compiledFormulaeParameters.lLength;

    _SimpleList * toDelete = nil;

    for (long k = 0; k<topLimit; k++) {
        long g = ((_SimpleList*)compiledFormulaeParameters.lData[k])->BinaryFind (vN,0);

        if (g>=0) {
            _ElementaryCommand* thisCommand = (_ElementaryCommand*)listOfCompiledFormulae.lData[k];
            _Formula  *f  = (_Formula*)(thisCommand->simpleParameters.lData[1]),
                       *f2 = (_Formula*)(thisCommand->simpleParameters.lData[2]);

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
        listOfCompiledFormulae.DeleteList       (*toDelete);
        compiledFormulaeParameters.DeleteList   (*toDelete);
        DeleteObject (toDelete);
    }
}

//__________________________________________________________________________________
void       UpdateChangingFlas (_SimpleList & involvedVariables)
{

    long          topLimit         = compiledFormulaeParameters.lLength;
    _SimpleList * toDelete         = nil;

    for (long k = 0; k<topLimit; k++) {
        long g = ((_SimpleList*)compiledFormulaeParameters.lData[k])->CountCommonElements (involvedVariables,true);

        if (g>0) {
            _ElementaryCommand* thisCommand = (_ElementaryCommand*)listOfCompiledFormulae.lData[k];

            _Formula  *f  = (_Formula*)(thisCommand->simpleParameters.lData[1]),
                       *f2 = (_Formula*)(thisCommand->simpleParameters.lData[2]);

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
        listOfCompiledFormulae.DeleteList       (*toDelete);
        compiledFormulaeParameters.DeleteList   (*toDelete);
        DeleteObject (toDelete);
    }
}

//__________________________________________________________________________________
void DeleteVariable (long dv, bool deleteself)
{
    if (dv>=0) {

        _String *name  = (_String*)variableNames.Retrieve (dv);
        _String myName = *name&'.';
        long    vidx   = variableNames.GetXtra (dv);

        UpdateChangingFlas (vidx);

        _SimpleList recCache;
        variableNames.Find (name,recCache);
        _String     nextVarID;// = *(_String*)variableNames.Retrieve(variableNames.Next (dv,recCache));
        long        nvid;
        if ((nvid = variableNames.Next (dv,recCache))>=0) {
            nextVarID = *(_String*)variableNames.Retrieve(nvid);
        }

        if (deleteself) {
            _SimpleList tcache;
            long        iv,
                        k = variableNames.Traverser (tcache, iv, variableNames.GetRoot());

            for (; k>=0; k = variableNames.Traverser (tcache, iv)) {
                _Variable * thisVar = FetchVar(k);

                if (thisVar->CheckFForDependence (vidx,false)) {
                    _PMathObj curValue = thisVar->Compute();
                    curValue->nInstances++; // this could be a leak 01/05/2004.
                    thisVar->SetValue (curValue);
                    DeleteObject (curValue);
                }
            }

            _Variable* delvar = (FetchVar(dv));
            DeleteObject (delvar);

            variableNames.Delete (variableNames.Retrieve(dv),true);
            (*((_SimpleList*)&variablePtrs))[vidx]=0;
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
            if (dependent.startswith(myName)) {
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
void DeleteTreeVariable (long dv, _SimpleList & parms, bool doDeps)
{
    if (dv>=0) {
        _String *name  = (_String*)variableNames.Retrieve (dv);
        _String myName = *name&".";
        long    vidx   = variableNames.GetXtra (dv);

        UpdateChangingFlas (vidx);

        _SimpleList recCache;
        variableNames.Find (name,recCache);
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

                if (thisVar->CheckFForDependence (vidx,false)) {
                    _PMathObj curValue = thisVar->Compute();
                    curValue->nInstances++;
                    thisVar->SetValue (curValue);
                    DeleteObject (curValue);
                }
            }
        }

        _Variable* delvar = (FetchVar(dv));
        if (delvar->ObjectClass() != TREE) {
            variableNames.Delete (variableNames.Retrieve(dv),true);
            (*((_SimpleList*)&variablePtrs))[vidx]=0;
            freeSlots<<vidx;
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
                if (dependent.startswith(myName)) {
                    if (dependent.Find ('.', myName.sLength+1, -1)>=0) {
                        _Variable * checkDep = FetchVar (nextVar);
                        if (!checkDep->IsIndependent()) {
                            _PMathObj curValue = checkDep->Compute();
                            curValue->nInstances++;
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
                DeleteTreeVariable (*(_String*)toDelete(k),parms,false);
            }
        }
    }
}
//__________________________________________________________________________________
void DeleteVariable (_String&name, bool deleteself)
{
    DeleteVariable(LocateVarByName (name), deleteself);
}

//__________________________________________________________________________________
void DeleteTreeVariable (_String&name, _SimpleList& parms, bool doDeps)
{
    DeleteTreeVariable(LocateVarByName (name), parms,doDeps);
}

//__________________________________________________________________________________
_Variable* CheckReceptacle (_String* name, _String fID, bool checkValid, bool isGlobal)
{
    if (checkValid && (!name->IsValidIdentifier())) {
        _String errMsg = *name & " is not a valid variable identifier in call to " & fID;
        WarnError (errMsg);
        return nil;
    }

    long    f = LocateVarByName (*name);
    if (f<0) {
        _Variable dummy (*name, isGlobal);
        f = LocateVarByName (*name);
    }

    return FetchVar(f);
}
//__________________________________________________________________________________
_Variable* CheckReceptacleCommandID (_String* name, const long id, bool checkValid, bool isGlobal, _ExecutionList* context)
{
    if (checkValid && (!name->IsValidIdentifier())) {
        _String errMsg = _String ("'") & *name & "' is not a valid variable identifier in call to " & _HY_ValidHBLExpressions.RetrieveKeyByPayload(id) & '.';
        if (context) {
            context->ReportAnExecutionError(errMsg);
        } else {
            WarnError (errMsg);
        }
        return nil;
    }
    
    long    f = LocateVarByName (*name);
    if (f<0) {
        _Variable dummy (*name, isGlobal);
        f = LocateVarByName (*name);
    }
    
    return FetchVar(f);
}

//__________________________________________________________________________________
bool CheckReceptacleCommandIDAndStore (_String* name, const long id, bool checkValid, _PMathObj v, bool dup, bool isGlobal)
{
    _Variable *theV = CheckReceptacleCommandID (name, id, checkValid, isGlobal);
    if (theV) {
        theV->SetValue (v, dup);
        return true;
    }
    if (!dup) {
        DeleteObject (v);
    }
    return false;
}


//__________________________________________________________________________________
bool CheckReceptacleAndStore (_String* name, _String fID, bool checkValid, _PMathObj v, bool dup)
{
    _Variable * theV = CheckReceptacle(name, fID, checkValid);
    if (theV) {
        theV->SetValue (v, dup);
        return true;
    }
    if (!dup) {
        DeleteObject (v);
    }
    return false;
}

//__________________________________________________________________________________
bool CheckReceptacleAndStore (_String name, _String fID, bool checkValid, _PMathObj v, bool dup)
{
    return CheckReceptacleAndStore (&name, fID, checkValid, v, dup);
}

//__________________________________________________________________________________
void  InsertVar (_Variable* theV)
{
    long pos = variableNames.Insert (theV->theName);

    /*if (theV->GetName()->Equal (&_String("PS_2")))
    {
        printf ("Making...\n");
    }*/


    if (pos < 0 && isDefiningATree > 1)
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
        if (isDefiningATree == 1) {
            WarnError(_String("Error while creating a tree: duplicate node name '") & *theV->GetName() & "'");
            return;
        }

        theV->theIndex = variableNames.GetXtra(-pos-1);
        return;
    } else {
        theV->theName->nInstances++;
    }

    if (freeSlots.lLength) {
        theV->theIndex = freeSlots.lData[freeSlots.lLength-1];
        variablePtrs[theV->theIndex]=theV->makeDynamic();
        freeSlots.Delete(freeSlots.lLength-1);
    } else {
        theV->theIndex = variablePtrs.lLength;
        variablePtrs&&theV;
    }
    variableNames.SetXtra (pos, theV->theIndex);
}

//__________________________________________________________________________________
_String&  AppendContainerName (_String& inString, _VariableContainer* theP)
{
    return AppendContainerName (inString, theP?theP->GetName():nil);
}

//__________________________________________________________________________________
_String&  AppendContainerName (_String& inString, _String* namescp)
{
    static _String returnMe;

    if (_hyApplicationGlobals.Find (&inString) >= 0) {
        return inString;
    }
    
    unsigned char reference_type = inString.ProcessVariableReferenceCases (returnMe, namescp && namescp -> sLength? namescp : nil);
    

    if (reference_type != HY_STRING_INVALID_REFERENCE) {
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

        for  (; f>=0 && ((_String*)variableNames.Retrieve (f))->startswith (oldNamePrefix); f = variableNames.Next (f, traverser)) {
            toRename << variableNames.Retrieve (f);
            xtras << variableNames.GetXtra (f);
        }
    }

    for (unsigned long k = 0; k < toRename.lLength; k++) {
        _Variable * thisVar = FetchVar (xtras.lData[k]);
        thisVar->GetName()->RemoveAReference();
        if (k) {
            thisVar->theName = new _String(thisVar->GetName()->Replace(oldNamePrefix,newNamePrefix,true));
        } else {
            thisVar->theName = new _String(*newName);
        }

        variableNames.Delete (toRename (k), true);
        variableNames.Insert (thisVar->GetName(),xtras.lData[k]);
        thisVar->GetName()->nInstances++;
    }
}

//__________________________________________________________________________________
void  ReplaceVar (_Variable* theV)
{
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
void    SetupOperationLists (void)
{

    _List all_unary_ops ("-",29,
                         "!",
                         "+",
                         "*",
                         "^",
                         "&",
                         "Abs",
                         "Sin",
                         "Cos",
                         "Tan",
                         "Exp",
                         "Log",
                         "Arctan",
                         "Time",
                         "Gamma",
                         "Transpose",
                         "Sqrt",
                         "Erf",
                         "Rows",
                         "Columns",
                         "LUDecompose",
                         "Inverse",
                         "BranchCount",
                         "TipCount",
                         "ZCDF",
                         "Eigensystem",
                         "Simplex",
                         "Type",
                         "Eval",
                         "LnGamma"
                         );
 

    UnOps.Insert (all_unary_ops);
        
    BinOps<<'|'*256+'|';
    opPrecedence<<1;
    BinOps<<'&'*256+'&';
    opPrecedence<<2;
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

    if (BuiltInFunctions.lLength==0)
        // construct a list of operations
        // don't forget to update SimplifyConstants, simpleOperationCodes, InternalDifferentiate, InternalSimplify, Formula::HasChanged and all Execute commands
        // also MAccess and MCoord codes are used in Parse to merge multiple matrix access operations
    {
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
        FunctionNameList << BuiltInFunctions (HY_OP_CODE_BETA);
        FunctionArgumentCount << 2;

        //HY_OP_CODE_BRANCHCOUNT
        BuiltInFunctions.AppendNewInstance (new _String ("BranchCount"));
        FunctionNameList << BuiltInFunctions (HY_OP_CODE_BRANCHCOUNT);
        FunctionArgumentCount << 2;

        //HY_OP_CODE_BRANCHLENGTH
        BuiltInFunctions.AppendNewInstance (new _String ("BranchLength"));
        FunctionNameList << BuiltInFunctions (HY_OP_CODE_BRANCHLENGTH);
        FunctionArgumentCount << 2;

        //HY_OP_CODE_BRANCHNAME
        BuiltInFunctions.AppendNewInstance (new _String ("BranchName"));
        FunctionNameList << BuiltInFunctions (HY_OP_CODE_BRANCHNAME);
        FunctionArgumentCount << 2;

        //HY_OP_CODE_CCHI2
        BuiltInFunctions.AppendNewInstance (new _String ("CChi2"));
        FunctionNameList << BuiltInFunctions (HY_OP_CODE_CCHI2);
        FunctionArgumentCount << 2;

        //HY_OP_CODE_CGAMMADIST
        BuiltInFunctions.AppendNewInstance (new _String ("CGammaDist"));
        FunctionNameList << BuiltInFunctions (HY_OP_CODE_CGAMMADIST);
        FunctionArgumentCount << 3;

        //HY_OP_CODE_COLUMNS
        BuiltInFunctions.AppendNewInstance (new _String ("Columns"));

        //HY_OP_CODE_COS
        BuiltInFunctions.AppendNewInstance (new _String ("Cos"));

        //HY_OP_CODE_DIFF
        BuiltInFunctions.AppendNewInstance (new _String ("Differentiate"));
        FunctionNameList << BuiltInFunctions (HY_OP_CODE_DIFF);
        FunctionArgumentCount << 2;

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
        FunctionNameList << BuiltInFunctions (HY_OP_CODE_FORMAT);
        FunctionArgumentCount << 3;

        //HY_OP_CODE_GAMMA
        BuiltInFunctions.AppendNewInstance (new _String ("Gamma"));

        //HY_OP_CODE_GAMMADIST
        BuiltInFunctions.AppendNewInstance (new _String ("GammaDist"));
        FunctionNameList << BuiltInFunctions (HY_OP_CODE_GAMMADIST);
        FunctionArgumentCount << 3;

        //HY_OP_CODE_IBETA
        BuiltInFunctions.AppendNewInstance (new _String ("IBeta"));
        FunctionNameList << BuiltInFunctions (HY_OP_CODE_IBETA);
        FunctionArgumentCount << 3;

        //HY_OP_CODE_IGAMMA
        BuiltInFunctions.AppendNewInstance (new _String ("IGamma"));
        FunctionNameList << BuiltInFunctions (HY_OP_CODE_IGAMMA);
        FunctionArgumentCount << 2;

        //HY_OP_CODE_INVCHI2
        BuiltInFunctions.AppendNewInstance (new _String ("InvChi2"));
        FunctionNameList << BuiltInFunctions (HY_OP_CODE_INVCHI2);
        FunctionArgumentCount << 2;

        //HY_OP_CODE_INVERSE
        BuiltInFunctions.AppendNewInstance (new _String ("Inverse"));

        //HY_OP_CODE_JOIN
        BuiltInFunctions.AppendNewInstance (new _String ("Join"));
        FunctionNameList << BuiltInFunctions (HY_OP_CODE_JOIN);
        FunctionArgumentCount << 2;

        //HY_OP_CODE_LUDECOMPOSE
        BuiltInFunctions.AppendNewInstance (new _String ("LUDecompose"));

        //HY_OP_CODE_LUSOLVE
        BuiltInFunctions.AppendNewInstance (new _String ("LUSolve"));
        FunctionNameList << BuiltInFunctions (HY_OP_CODE_LUSOLVE);
        FunctionArgumentCount << 2;

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

        //HY_OP_CODE_MAX
        BuiltInFunctions.AppendNewInstance (new _String ("Max"));
        FunctionNameList << BuiltInFunctions (HY_OP_CODE_MAX);
        FunctionArgumentCount << 2;
        simpleOperationCodes<<HY_OP_CODE_MAX;
        simpleOperationFunctions<<(long)MaxNumbers;

        //HY_OP_CODE_MIN
        BuiltInFunctions.AppendNewInstance (new _String ("Min"));
        FunctionNameList << BuiltInFunctions (HY_OP_CODE_MIN);
        FunctionArgumentCount << 2;
        simpleOperationCodes<<HY_OP_CODE_MIN;
        simpleOperationFunctions<<(long)MinNumbers;

        //HY_OP_CODE_PSTREESTRING
        BuiltInFunctions.AppendNewInstance (new _String ("PSTreeString"));
        FunctionNameList << BuiltInFunctions (HY_OP_CODE_PSTREESTRING);
        FunctionArgumentCount << 3;

        //HY_OP_CODE_RANDOM
        BuiltInFunctions.AppendNewInstance (new _String ("Random"));
        FunctionNameList << BuiltInFunctions (HY_OP_CODE_RANDOM);
        FunctionArgumentCount << 2;
        simpleOperationCodes<<HY_OP_CODE_RANDOM;
        simpleOperationFunctions<<(long)RandomNumber;

        //HY_OP_CODE_REROOTTREE
        BuiltInFunctions.AppendNewInstance (new _String ("RerootTree"));
        FunctionNameList << BuiltInFunctions (HY_OP_CODE_REROOTTREE);
        FunctionArgumentCount << 2;

        //HY_OP_CODE_ROWS
        BuiltInFunctions.AppendNewInstance (new _String ("Rows"));

        //HY_OP_CODE_SIMPLEX
        BuiltInFunctions.AppendNewInstance (new _String ("Simplex"));

        //HY_OP_CODE_SIN
        BuiltInFunctions.AppendNewInstance (new _String ("Sin"));

        //HY_OP_CODE_SQRT
        BuiltInFunctions.AppendNewInstance (new _String ("Sqrt"));

        //HY_OP_CODE_TEXTREESTRING
        BuiltInFunctions.AppendNewInstance (new _String ("TEXTreeString"));
        FunctionNameList << BuiltInFunctions (HY_OP_CODE_TEXTREESTRING);
        FunctionArgumentCount << 2;

        //HY_OP_CODE_TAN
        BuiltInFunctions.AppendNewInstance (new _String ("Tan"));

        //HY_OP_CODE_TIME
        BuiltInFunctions.AppendNewInstance (new _String ("Time"));

        //HY_OP_CODE_TIPCOUNT
        BuiltInFunctions.AppendNewInstance (new _String ("TipCount"));

        //HY_OP_CODE_TIPNAME
        BuiltInFunctions.AppendNewInstance (new _String ("TipName"));
        FunctionNameList << BuiltInFunctions (HY_OP_CODE_TIPNAME);
        FunctionArgumentCount << 2;

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
    }



}

//__________________________________________________________________________________
void    CompileListOfUserExpressions (_SimpleList& varRefs,_List& rec, bool doAll)
{
    rec.Clear();
    if (varRefs.lLength == 0) {
        return;
    }

    long i;
    _SimpleList startVars;
    _VariableContainer*  firstVar = (_VariableContainer*)LocateVar(varRefs.lData[0]);

    firstVar->ScanAndAttachVariables();

    {
        _AVLList sA (&startVars);
        if (doAll) {

            firstVar->ScanForVariables (sA,sA);
            firstVar->ScanForGVariables (sA,sA);
        }

        firstVar->ScanForDVariables (sA,sA);
        sA.ReorderList ();
    }

    if (!doAll) {
        for (i=startVars.lLength-1; i>=0; i--) {
            if (firstVar->IsModelVar(i)) {
                startVars.Delete(i);
            }
        }
    }

    for (i=0; i<startVars.lLength; i++) {
        _String thisName (LocateVar(startVars.lData[i])->GetName()->Cut
                          (LocateVar(startVars.lData[i])->GetName()->FindBackwards('.',0,-1),-1));
        rec && &thisName;
    }

    for (i=varRefs.lLength-1; i>=1; i--) {
        firstVar = (_VariableContainer*)LocateVar(varRefs.lData[i]);
        firstVar->ScanAndAttachVariables();
        firstVar->MatchParametersToList (rec,doAll);
    }

    for (i=rec.lLength-1; i>=0; i--) {
        _String* thisLine = ((_String*)rec(i));
        thisLine->Trim(1,-1);
        if (doAll)
            if (LocateVarByName(*thisLine)<0) {
                *thisLine = _String('!')&*thisLine;
            }
    }

}

//__________________________________________________________________________________

void  FindUnusedObjectName (_String& prefix, _String& partName, _List& names, bool sorted)
{
    if (partName.sLength==0) {
        partName = prefix;
    }

    _String tryName (partName);
    long    k = 1;

    if (sorted)
        while (names.BinaryFind(&tryName)>=0) {
            k++;
            tryName = partName&k;
        }
    else
        while (names.Find(&tryName)>=0) {
            k++;
            tryName = partName&k;
        }

    partName = tryName;
}
//__________________________________________________________________________________

void  FindUnusedObjectName (_String& prefix, _String& partName, _AVLListX& names, bool)
{
    if (partName.sLength==0) {
        partName = prefix;
    }

    _String tryName (partName);
    long    k = 1;

    while (names.Find(&tryName)>=0) {
        k++;
        tryName = partName&k;
    }

    partName = tryName;
}

//__________________________________________________________________________________

void  FinishDeferredSF (void)
{
    if (deferSetFormula->lLength) {
        SortLists (deferSetFormula, &deferIsConstant);
        _SimpleList tcache;
        long        iv,
                    i = variableNames.Traverser (tcache,iv,variableNames.GetRoot());

        for (; i >= 0; i = variableNames.Traverser (tcache,iv)) {
            _Variable* theV = FetchVar(i);
            if (theV->IsContainer()) {
                ((_VariableContainer*)theV)->SetMDependance (*deferSetFormula);
            }
        }

        for (long j = 0; j<likeFuncList.lLength; j++)
            if (((_String*)likeFuncNamesList(j))->sLength) {
                _LikelihoodFunction * lf = (_LikelihoodFunction*)likeFuncList(j);
                for (long k = 0; k < deferSetFormula->lLength; k++) {
                    lf->UpdateIndependent(deferSetFormula->lData[k],deferIsConstant.lData[k]);
                }
            }
    }
    DeleteObject (deferSetFormula);
    deferSetFormula = nil;
    deferIsConstant.Clear();
}
