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

#ifndef     __PARSER__
#define     __PARSER__

#include "baseobj.h"
#include "list.h"
#include "avllistx.h"
#include "hy_strings.h"
#include "errorfns.h"
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
            *deferSetFormula;

extern      _AVLListX       variableNames;

extern      _String         HalfOps;

extern      _Trie           UnOps;

extern      _Parameter      printDigits,
            verbosityLevel;

extern      long            lastMatrixDeclared,
            subNumericValues;

long        LocateVarByName (_String&);
_Variable*  LocateVar       (long index);
_Variable*  FetchVar        (long index);
_PMathObj   FetchObjectFromVariableByType       (_String*, const unsigned long, long = -1, _String* = nil);
_PMathObj   FetchObjectFromVariableByTypeIndex  (long, const unsigned long, long = -1, _String* = nil);
_String     FetchObjectNameFromType (const unsigned long);
_String&    AppendContainerName     (_String&, _VariableContainer*);
_String&    AppendContainerName     (_String&, _String*);
_String*    FetchMathObjectNameOfTypeByIndex (const unsigned long objectClass, const long objectIndex);

void        DeleteVariable  (_String&, bool deleteself = true);
void        DeleteVariable  (long, bool deleteself);

void        DeleteTreeVariable
(_String&, _SimpleList&,bool);
void        checkParameter  (_String& name, _Parameter& dest, _Parameter def, _VariableContainer* = nil);
void        stashParameter  (_String& name, _Parameter  newVal, bool);
void        setParameter    (_String& name, _Parameter def, _String* = nil);
void        setParameter    (_String& name, _PMathObj  def, bool = true, _String* = nil);

long        VerbosityLevel (void);
void        ReplaceVar      (_Variable*);
void        RenameVariable  (_String*,_String*);
void        CompileListOfUserExpressions
(_SimpleList&,_List&, bool doAll = false);

void        FindUnusedObjectName
(_String&, _String&, _List&, bool = false);

void        FindUnusedObjectName
(_String&, _String&, _AVLListX&, bool = false);

bool        ExpressionCalculator(void);
bool        ExpressionCalculator(_String data);

_Variable*  CheckReceptacle
(_String*,_String, bool = true, bool = false);

bool        CheckReceptacleAndStore
(_String*,_String, bool, _PMathObj, bool = true);

bool        CheckReceptacleAndStore
(_String,_String, bool, _PMathObj, bool = true);

_Variable*  CheckReceptacleCommandID
(_String* name, const long id, bool checkValid, bool isGlobal = false, _ExecutionList* context = nil);

bool        CheckReceptacleCommandIDAndStore
(_String* name, const long id, bool checkValid, _PMathObj v, bool dup = true, bool isGlobal = false);

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
bool        CheckEqual      (_Parameter,_Parameter);

extern      _AVLListX       _hyApplicationGlobals;

/**
    A utility function to take a list of variable indices and split them into local and global

    @param   inList supplies the input list (all variable indices are assumed to be valid!)
    @param   outList a list of two _SimpleLists: index 0 for global variables, and index 1 for local variables
    @author  SLKP
    @version 20110608
*/
void        SplitVariableIDsIntoLocalAndGlobal (const _SimpleList& inList, _List& outList);

_Parameter  AddNumbers  (_Parameter, _Parameter);
_Parameter  SubNumbers  (_Parameter, _Parameter);
_Parameter  MultNumbers (_Parameter, _Parameter);
_Parameter  AndNumbers  (_Parameter, _Parameter);
_Parameter  DivNumbers  (_Parameter, _Parameter);
_Parameter  EqualNumbers(_Parameter, _Parameter);
_Parameter  LessThan    (_Parameter, _Parameter);
_Parameter  GreaterThan (_Parameter, _Parameter);
_Parameter  LessThanE   (_Parameter, _Parameter);
_Parameter  GreaterThanE(_Parameter, _Parameter);
_Parameter  Power       (_Parameter, _Parameter);
_Parameter  RandomNumber(_Parameter, _Parameter);
_Parameter  ExpNumbers  (_Parameter);
_Parameter  LogNumbers  (_Parameter);
_Parameter  AbsNumber   (_Parameter);
_Parameter  MinusNumber (_Parameter);
_Parameter  MaxNumbers  (_Parameter, _Parameter);
_Parameter  MinNumbers  (_Parameter, _Parameter);
_Parameter  FastMxAccess(Ptr, _Parameter);

BaseRef parameterToString       (_Parameter);
void    parameterToCharBuffer   (_Parameter, char*, long, bool json = false);

_Parameter  InterpolateValue        (_Parameter*, _Parameter*, long, _Parameter*, _Parameter*, _Parameter, _Parameter&);
_Parameter  TrapezoidLevelK         (_Formula&, _Variable*, _Parameter, _Parameter, long);
_Parameter  TrapezoidLevelKSimple   (_Formula&, _Variable*, _Parameter, _Parameter, long, _SimpleFormulaDatum*, _SimpleFormulaDatum*, _SimpleList&, _SimpleList&);


void        PopulateArraysForASimpleFormula
(_SimpleList&, _SimpleFormulaDatum*);

void        WarnNotDefined (_PMathObj, long, _hyExecutionContext* );

extern      _Parameter  pi_const;
extern      bool        useGlobalUpdateFlag;
extern      _String     noneToken;

#endif
