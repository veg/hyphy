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

#include "avllistx.h"
#include "avllistxl.h"
#include "baseobj.h"
#include "hy_strings.h"
#include "errorfns.h"
#include "stdio.h"
#include "classes.h"


#define  HY_UNDEFINED       0x00
#define  NUMBER             0x01
#define  MATRIX             0x04
#define  CONTAINER          0x08
#define  TREE_NODE          0x10
#define  TREE               0x20
#define  STRING             0x40
#define  ASSOCIATIVE_LIST   0x80
#define  TOPOLOGY           0x100
#define  POLYNOMIAL         0x200
#define  HY_ANY_OBJECT      0xFFFF

#define  DEFAULTLOWERBOUND  -1e26
#define  DEFAULTUPPERBOUND  1e26

#define  GLOBAL_VARIABLE   1
#define  CATEGORY_VARIABLE 2
#define  RANDOM_VARIABLE   3

#define  HY_VARIABLE_GLOBAL     0x01
#define  HY_VARIABLE_CHANGED    0x02
#define  HY_DEP_V_COMPUTED      0x04
#define  HY_DEP_V_INSPECTED     0x08
#define  HY_DEP_V_INSPECTED_CLR 0xF7
#define  HY_DEP_V_MODIFIED      0x10
#define  HY_DEP_V_MODIFIED_CATS 0x20
#define  HY_VC_NO_CHECK         0x40 // do not check this variable container in 
// NeedToExponentiate
#define  HY_VARIABLE_NOTSET     0x80
#define  HY_VARIABLE_SET        0x7F

#define  HY_VC_CLR_NO_CHECK     0xBF

#define  HY_DEP_CLEAR_MASK      0xC7

#define  HY_NO_MODEL            (-1)


// START OPCODES

#define  HY_OP_CODE_NONE            -1                              // noop

#define  HY_OP_CODE_NOT             0                               // !
#define  HY_OP_CODE_NEQ             (1+HY_OP_CODE_NOT)              // !=
#define  HY_OP_CODE_IDIV            (1+HY_OP_CODE_NEQ)              // $
#define  HY_OP_CODE_MOD             (1+HY_OP_CODE_IDIV)             // %
#define  HY_OP_CODE_AND             (1+HY_OP_CODE_MOD)              // &&
#define  HY_OP_CODE_MUL             (1+HY_OP_CODE_AND)              // *
#define  HY_OP_CODE_ADD             (1+HY_OP_CODE_MUL)              // +
#define  HY_OP_CODE_SUB             (1+HY_OP_CODE_ADD)              // -
#define  HY_OP_CODE_DIV             (1+HY_OP_CODE_SUB)              // /
#define  HY_OP_CODE_LESS            (1+HY_OP_CODE_DIV)              // <

//Check for correctness, then move on
#define  HY_OP_CODE_LEQ             (1+HY_OP_CODE_LESS)             // <=
#define  HY_OP_CODE_EQ              (1+HY_OP_CODE_LEQ)              // ==
#define  HY_OP_CODE_GREATER         (1+HY_OP_CODE_EQ)               // >
#define  HY_OP_CODE_GEQ             (1+HY_OP_CODE_GREATER)          // >=

#define  HY_OP_CODE_ABS             (1+HY_OP_CODE_GEQ)              // Abs
#define  HY_OP_CODE_ARCTAN          (1+HY_OP_CODE_ABS)              // Arctan
#define  HY_OP_CODE_BETA            (1+HY_OP_CODE_ARCTAN)           // Beta
#define  HY_OP_CODE_BRANCHCOUNT     (1+HY_OP_CODE_BETA)             // BranchCount
#define  HY_OP_CODE_BRANCHLENGTH    (1+HY_OP_CODE_BRANCHCOUNT)      // BranchLength
#define  HY_OP_CODE_BRANCHNAME      (1+HY_OP_CODE_BRANCHLENGTH)     // BranchName

#define  HY_OP_CODE_CCHI2           (1+HY_OP_CODE_BRANCHNAME)       // CChi2
#define  HY_OP_CODE_CGAMMADIST      (1+HY_OP_CODE_CCHI2)            // CGammaDist
#define  HY_OP_CODE_COLUMNS         (1+HY_OP_CODE_CGAMMADIST)       // Columns
#define  HY_OP_CODE_COS             (1+HY_OP_CODE_COLUMNS)          // Cos
#define  HY_OP_CODE_DIFF            (1+HY_OP_CODE_COS)              // D
#define  HY_OP_CODE_EIGENSYSTEM     (1+HY_OP_CODE_DIFF)             // Eigensystem
#define  HY_OP_CODE_ERF             (1+HY_OP_CODE_EIGENSYSTEM)      // Erf
#define  HY_OP_CODE_EVAL            (1+HY_OP_CODE_ERF)              // Eval
#define  HY_OP_CODE_EXP             (1+HY_OP_CODE_EVAL)             // Exp
#define  HY_OP_CODE_FORMAT          (1+HY_OP_CODE_EXP)              // Format
#define  HY_OP_CODE_GAMMA           (1+HY_OP_CODE_FORMAT)           // Gamma

#define  HY_OP_CODE_GAMMADIST       (1+HY_OP_CODE_GAMMA)            // GammaDist
#define  HY_OP_CODE_IBETA           (1+HY_OP_CODE_GAMMADIST)        // IBeta
#define  HY_OP_CODE_IGAMMA          (1+HY_OP_CODE_IBETA)            // IGamma
#define  HY_OP_CODE_INVCHI2         (1+HY_OP_CODE_IGAMMA)           // InvChi2
#define  HY_OP_CODE_INVERSE         (1+HY_OP_CODE_INVCHI2)          // Inverse
#define  HY_OP_CODE_JOIN            (1+HY_OP_CODE_INVERSE)          // Join
#define  HY_OP_CODE_LUDECOMPOSE     (1+HY_OP_CODE_JOIN)             // LUDecompose
#define  HY_OP_CODE_LUSOLVE         (1+HY_OP_CODE_LUDECOMPOSE)      // LUSolve
#define  HY_OP_CODE_LNGAMMA         (1+HY_OP_CODE_LUSOLVE)          // LnGamma
#define  HY_OP_CODE_LOG             (1+HY_OP_CODE_LNGAMMA)          // Log
#define  HY_OP_CODE_MACCESS         (1+HY_OP_CODE_LOG)              // MAccess
#define  HY_OP_CODE_MCOORD          (1+HY_OP_CODE_MACCESS)          // MCoord

#define  HY_OP_CODE_MAX             (1+HY_OP_CODE_MCOORD)           // Max
#define  HY_OP_CODE_MIN             (1+HY_OP_CODE_MAX)              // Min
#define  HY_OP_CODE_PSTREESTRING    (1+HY_OP_CODE_MIN)              // PSTreeString
#define  HY_OP_CODE_RANDOM          (1+HY_OP_CODE_PSTREESTRING)     // Random
#define  HY_OP_CODE_REROOTTREE      (1+HY_OP_CODE_RANDOM)           // RerootTree
#define  HY_OP_CODE_ROWS            (1+HY_OP_CODE_REROOTTREE)       // Rows
#define  HY_OP_CODE_SIMPLEX         (1+HY_OP_CODE_ROWS)             // Simplex
#define  HY_OP_CODE_SIN             (1+HY_OP_CODE_SIMPLEX)          // Sin
#define  HY_OP_CODE_SQRT            (1+HY_OP_CODE_SIN)              // Sqrt
#define  HY_OP_CODE_TEXTREESTRING   (1+HY_OP_CODE_SQRT)             // TEXTreeString

#define  HY_OP_CODE_TAN             (1+HY_OP_CODE_TEXTREESTRING)    // Tan
#define  HY_OP_CODE_TIME            (1+HY_OP_CODE_TAN)              // Time
#define  HY_OP_CODE_TIPCOUNT        (1+HY_OP_CODE_TIME)             // TipCount
#define  HY_OP_CODE_TIPNAME         (1+HY_OP_CODE_TIPCOUNT)         // TipName
#define  HY_OP_CODE_TRANSPOSE       (1+HY_OP_CODE_TIPNAME)          // Transpose
#define  HY_OP_CODE_TYPE            (1+HY_OP_CODE_TRANSPOSE)        // Type
#define  HY_OP_CODE_ZCDF            (1+HY_OP_CODE_TYPE)             // ZCDF
#define  HY_OP_CODE_POWER           (1+HY_OP_CODE_ZCDF)             // ^
#define  HY_OP_CODE_OR              (1+HY_OP_CODE_POWER)            // ||

// END OPCODES

// START FORMULA RETURN CODES

#define  HY_FORMULA_EXPRESSION                          0
#define  HY_FORMULA_FAILED                              (-1)
#define  HY_FORMULA_VARIABLE_VALUE_ASSIGNMENT           1
#define  HY_FORMULA_VARIABLE_FORMULA_ASSIGNMENT         2
#define  HY_FORMULA_FORMULA_FORMULA_ASSIGNMENT          3
#define  HY_FORMULA_FORMULA_VALUE_ASSIGNMENT            4
#define  HY_FORMULA_VARIABLE_LOWER_BOUND_ASSIGNMENT     5
#define  HY_FORMULA_VARIABLE_UPPER_BOUND_ASSIGNMENT     6
#define  HY_FORMULA_FORMULA_VALUE_INCREMENT             7
// END FORMULA RETURN CODES

#define  HY_INVALID_RETURN_VALUE                        NAN

// pointer to a math object
typedef _MathObject* _PMathObj ;

//__________________________________________________________________________________

extern  _List BuiltInFunctions;



//__________________________________________________________________________________

union       _SimpleFormulaDatum {
    _Parameter value;
    Ptr        reference;
};


//__________________________________________________________________________________
#include "matrix.h"

//__________________________________________________________________________________

#define  USE_POINTER_VC

//__________________________________________________________________________________


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

void        PopulateArraysForASimpleFormula
(_SimpleList&, _SimpleFormulaDatum*);

void        WarnNotDefined (_PMathObj, long);

extern      _Parameter  pi_const;
extern      bool        useGlobalUpdateFlag;

#endif
