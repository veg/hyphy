/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-now
Core Developers:
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon    (apoon42@uwo.ca)
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

#ifndef     __PARSER__
#define     __PARSER__

#include "baseobj.h"
#include "list.h"
#include "avllistx.h"
#include "hy_strings.h"
#include "hy_string_buffer.h"
#include "stdio.h"
#include "classes.h"
#include "associative_list.h"

//Parser specific includes
#include "constant.h"
#include "defines.h"
#include "formula.h"
#include "fstring.h"
#include "mathobj.h"
#include "operation.h"
#include "stack.h"
#include "variable.h"
#include "variablecontainer.h"
#include "trie.h"
#include "hbl_env.h"
#include "global_things.h"


//End parser specific includes

#include "matrix.h"

extern  _List BuiltInFunctions;

class   _ExecutionList; 
class   _hyExecutionContext;


//__________________________________________________________________________________

#define  USE_POINTER_VC

extern      _List           variablePtrs,
            hyReservedWords;

extern      _SimpleList     BuiltInFunctionParameterCount,
            *deferSetFormula,
            opPrecedence,
            BinOps,
            associativeOps,
            simpleOperationCodes,
            simpleOperationFunctions;


extern      _AVLListX       variableNames;
extern      _AVLList        *deferClearConstraint;

extern      _String         HalfOps;

extern      _Trie           UnOps,
                            FunctionNameList;


extern      long            lastMatrixDeclared;

long        LocateVarByName (_String const&);
inline _Variable*  LocateVar       (long index) {
        return (_Variable *)(((BaseRef*)variablePtrs.list_data)[index]);
}
HBLObjectRef   FetchObjectFromVariableByType       (_String const*, const unsigned long, long = -1, _String* = nil);
HBLObjectRef   FetchObjectFromVariableByTypeIndex  (long, const unsigned long, long = -1, _String* = nil);
HBLObjectRef   FetchObjectFromVariableByTypeIndexWithoutCompute  (long, const unsigned long, long = -1, _String* = nil);
HBLObjectRef   FetchObjectFromFormulaByType         (_Formula&, const unsigned long, long = -1, _String* = nil);
_String     FetchObjectNameFromType (const unsigned long);
_String const&    AppendContainerName     (_String const&, _VariableContainer const*);
_String const&    AppendContainerName     (_String const&, _String const*);
_String*    FetchMathObjectNameOfTypeByIndex (const unsigned long objectClass, const long objectIndex);

void        DeleteVariable  (_String const&, bool deleteself = true, bool do_checks = true);
void        DeleteVariable  (long, bool deleteself, bool do_checks = true);

void        SetVariablesToOwnValues (_AVLList const & indices);

void        DeleteTreeVariable (long, _SimpleList&, _String * prefix = nil, bool return_leftovers = true);

template <typename ACTION>
void DoForEachVariable(ACTION do_this) {
    for (unsigned long i = 0; i < variablePtrs.lLength; i++) {
        _Variable * ith_var = (_Variable *)variablePtrs.GetItem (i);
        if (ith_var) {
            do_this (ith_var, i);
        }
    }
}


void        stashParameter  (_String const& name, hyFloat  newVal, bool);
void        setParameter    (_String const& name, hyFloat def, _String* = nil);
void        setParameter    (_String const& name, HBLObjectRef  def, _String* = nil, bool = true);

long        VerbosityLevel (void);
void        ReplaceVar      (_Variable*);
void        RenameVariable  (_String*,_String*);
void        CompileListOfUserExpressions
(_SimpleList&,_List&, bool doAll = false);

bool        ExpressionCalculator(void);
bool        ExpressionCalculator(_String data);

_Variable*  CheckReceptacle
(_String const*, _String const & , bool = true, bool = false, bool clear_trees = true);

bool        CheckReceptacleAndStore
(_String const*,_String, bool, HBLObjectRef, bool = true);

bool        CheckReceptacleAndStore
(_String,_String, bool, HBLObjectRef, bool = true);

_Variable*  CheckReceptacleCommandIDException
(_String const* name, const long id, bool checkValid, bool isGlobal = false, _ExecutionList* context = nil);

_Variable*  CheckReceptacleCommandID
(_String const* name, const long id, bool checkValid, bool isGlobal = false, _ExecutionList* context = nil);

bool        CheckReceptacleCommandIDAndStore
(_String const* name, const long id, bool checkValid, HBLObjectRef v, bool dup = true, bool isGlobal = false);

void        FinishDeferredSF(void);

void        SetupOperationLists (void);
void        ExportIndVariables
(_StringBuffer&, _StringBuffer&, _SimpleList*, _AssociativeList* global_remap = nil, _AssociativeList * subs = nil);
void        ExportDepVariables
(_StringBuffer&, _StringBuffer&, _SimpleList*, _AssociativeList* global_remap = nil, _AssociativeList * subs = nil);
void        ExportCatVariables
(_StringBuffer&, _SimpleList*);

void        SplitVariablesIntoClasses
(_SimpleList&, _SimpleList&, _SimpleList&, _SimpleList&);

bool        CheckEqual      (hyFloat,hyFloat,hyFloat = hy_global::kMachineEpsilon);
bool        CheckRange      (hyFloat value,hyFloat lb, hyFloat ub, bool exeptions = false);
bool        CheckArgumentType  (HBLObjectRef object, long type, bool exceptions = false);

extern      _AVLListX       _hyApplicationGlobals;

/**
    A utility function to take a list of variable indices and split them into local and global

    @param   inList supplies the input list (all variable indices are assumed to be valid!)
    @param   outList a list of two _SimpleLists: index 0 for global variables, and index 1 for local variables
    @author  SLKP
    @version 20110608
*/
void        SplitVariableIDsIntoLocalAndGlobal (const _SimpleList& inList, _List& outList);

hyFloat  AddNumbers  (hyFloat, hyFloat);
hyFloat  SubNumbers  (hyFloat, hyFloat);
hyFloat  MultNumbers (hyFloat, hyFloat);
hyFloat  AndNumbers  (hyFloat, hyFloat);
hyFloat  DivNumbers  (hyFloat, hyFloat);
hyFloat  EqualNumbers(hyFloat, hyFloat);
hyFloat  LessThan    (hyFloat, hyFloat);
hyFloat  GreaterThan (hyFloat, hyFloat);
hyFloat  LessThanE   (hyFloat, hyFloat);
hyFloat  GreaterThanE(hyFloat, hyFloat);
hyFloat  Power       (hyFloat, hyFloat);
hyFloat  RandomNumber(hyFloat, hyFloat);
hyFloat  ExpNumbers  (hyFloat);
hyFloat  LogNumbers  (hyFloat);
hyFloat  AbsNumber   (hyFloat);
hyFloat  MinusNumber (hyFloat);
hyFloat  MaxNumbers  (hyFloat, hyFloat);
hyFloat  MinNumbers  (hyFloat, hyFloat);
hyFloat  FastMxAccess(hyPointer, hyFloat);
void        FastMxWrite (hyPointer, hyFloat, hyFloat);

BaseRef parameterToString       (hyFloat);
void    parameterToCharBuffer   (hyFloat, char*, long, bool json = false);

hyFloat  InterpolateValue        (hyFloat*, hyFloat*, long, hyFloat*, hyFloat*, hyFloat, hyFloat&);
hyFloat  TrapezoidLevelK         (_Formula&, _Variable*, hyFloat, hyFloat, long);
hyFloat  TrapezoidLevelKSimple   (_Formula&, _Variable*, hyFloat, hyFloat, long, _SimpleFormulaDatum*, _SimpleFormulaDatum*, _SimpleList&, _SimpleList&);


void        PopulateArraysForASimpleFormula
(_SimpleList&, _SimpleFormulaDatum*);

void        WarnNotDefined (HBLObjectRef, long, _hyExecutionContext* );
void        WarnWrongNumberOfArguments (HBLObjectRef, long, _hyExecutionContext*, _List *);
  
extern      hyFloat  pi_const, tolerance;


extern      bool        useGlobalUpdateFlag;
extern      _AVLList    *_keepTrackOfDepVars;

#endif
