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
#include "stdio.h"
#include "classes.h"

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

extern      _String         HalfOps;

extern      _Trie           UnOps;

extern      long            printDigits,
                            verbosityLevel;

extern      long            lastMatrixDeclared,
            subNumericValues;

long        LocateVarByName (_String const&);
_Variable*  LocateVar       (long index);
_PMathObj   FetchObjectFromVariableByType       (_String const*, const unsigned long, long = -1, _String* = nil);
_PMathObj   FetchObjectFromVariableByTypeIndex  (long, const unsigned long, long = -1, _String* = nil);
_PMathObj   FetchObjectFromFormulaByType         (_Formula&, const unsigned long, long = -1, _String* = nil);
_String     FetchObjectNameFromType (const unsigned long);
_String const&    AppendContainerName     (_String const&, _VariableContainer const*);
_String const&    AppendContainerName     (_String const&, _String const*);
_String*    FetchMathObjectNameOfTypeByIndex (const unsigned long objectClass, const long objectIndex);

void        DeleteVariable  (_String const&, bool deleteself = true);
void        DeleteVariable  (long, bool deleteself);

void        DeleteTreeVariable
(_String&, _SimpleList&,bool);



void        stashParameter  (_String const& name, hyFloat  newVal, bool);
void        setParameter    (_String const& name, hyFloat def, _String* = nil);
void        setParameter    (_String const& name, _PMathObj  def, _String* = nil, bool = true);

long        VerbosityLevel (void);
void        ReplaceVar      (_Variable*);
void        RenameVariable  (_String*,_String*);
void        CompileListOfUserExpressions
(_SimpleList&,_List&, bool doAll = false);

bool        ExpressionCalculator(void);
bool        ExpressionCalculator(_String data);

_Variable*  CheckReceptacle
(_String const*, _String const & , bool = true, bool = false);

bool        CheckReceptacleAndStore
(_String const*,_String, bool, _PMathObj, bool = true);

bool        CheckReceptacleAndStore
(_String,_String, bool, _PMathObj, bool = true);

_Variable*  CheckReceptacleCommandID
(_String const* name, const long id, bool checkValid, bool isGlobal = false, _ExecutionList* context = nil);

bool        CheckReceptacleCommandIDAndStore
(_String const* name, const long id, bool checkValid, _PMathObj v, bool dup = true, bool isGlobal = false);

void        FinishDeferredSF(void);

void        SetupOperationLists (void);
void        ExportIndVariables
(_String&, _String&, _SimpleList*);
void        ExportDepVariables
(_String&, _String&, _SimpleList*);
void        ExportCatVariables
(_String&, _SimpleList*);

void        SplitVariablesIntoClasses
(_SimpleList&, _SimpleList&, _SimpleList&, _SimpleList&);
bool        CheckEqual      (hyFloat,hyFloat);

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

void        WarnNotDefined (_PMathObj, long, _hyExecutionContext* );
void        WarnWrongNumberOfArguments (_PMathObj, long, _hyExecutionContext*, _List *);
  
extern      hyFloat  pi_const;
extern      bool        useGlobalUpdateFlag;
extern      _String     noneToken;

#endif
