/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2006
Primary Development:
  Sergei L Kosakovsky Pond (sergeilkp@mac.com)
Significant contributions from:
  Spencer V Muse (muse@stat.ncsu.edu)
  Simon DW Frost (sdfrost@ucsd.edu)
  Art FY Poon    (apoon@biomail.ucsd.edu)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

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
//End parser specific includes

#include "matrix.h"

//__________________________________________________________________________________
extern  _List BuiltInFunctions;
//__________________________________________________________________________________

#define  USE_POINTER_VC

extern      _List           variablePtrs,
            hyReservedWords;

extern      _SimpleList     BuiltInFunctionParameterCount,
            *deferSetFormula;

extern      _AVLListX       variableNames;

extern      _String         UnOps,
            HalfOps;

extern      _Parameter      printDigits,
            verbosityLevel;

extern      long            lastMatrixDeclared,
            subNumericValues;

long        LocateVarByName (_String&);
_Variable*  LocateVar       (long index);
_Variable*  FetchVar        (long index);
_PMathObj   FetchObjectFromVariableByType
(_String*, int);
_PMathObj   FetchObjectFromVariableByTypeIndex
(long, int);
_String&    AppendContainerName
(_String&, _VariableContainer*);

void        DeleteVariable  (_String&, bool deleteself = true);
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

bool        ExpressionCalculator
(void);

_Variable*  CheckReceptacle
(_String*,_String, bool = true, bool = false);

bool        CheckReceptacleAndStore
(_String*,_String, bool, _PMathObj, bool = true);

bool        CheckReceptacleAndStore
(_String,_String, bool, _PMathObj, bool = true);

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
_Parameter  MinusNumber (_Parameter);
_Parameter  MaxNumbers  (_Parameter, _Parameter);
_Parameter  MinNumbers  (_Parameter, _Parameter);
_Parameter  FastMxAccess(Ptr, _Parameter);

BaseRef parameterToString (_Parameter);
_Parameter  InterpolateValue        (_Parameter*, _Parameter*, long, _Parameter*, _Parameter*, _Parameter, _Parameter&);
_Parameter  TrapezoidLevelK         (_Formula&, _Variable*, _Parameter, _Parameter, long);
_Parameter  TrapezoidLevelKSimple   (_Formula&, _Variable*, _Parameter, _Parameter, long, _SimpleFormulaDatum*, _SimpleFormulaDatum*, _SimpleList&, _SimpleList&);


void        PopulateArraysForASimpleFormula
(_SimpleList&, _SimpleFormulaDatum*);

void        WarnNotDefined (_PMathObj, long);

extern      _Parameter  pi_const;
extern      bool        useGlobalUpdateFlag;

#endif
