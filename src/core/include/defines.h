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

#ifndef     __DEFINES__
#define     __DEFINES__

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

#define  GOLDEN_RATIO 1.618034
#define MAX_BRENT_ITERATES 100L

// START OPCODES

#define  HY_OP_CODE_NONE            -1                              // noop

#define  HY_OP_CODE_NOT             0                               // !
#define  HY_OP_CODE_NEQ             (1+HY_OP_CODE_NOT)              // !=
#define  HY_OP_CODE_IDIV            (1+HY_OP_CODE_NEQ)              // $
#define  HY_OP_CODE_MOD             (1+HY_OP_CODE_IDIV)             // %
#define  HY_OP_CODE_AND             (1+HY_OP_CODE_MOD)              // &&
#define  HY_OP_CODE_MUL             (1+HY_OP_CODE_AND)              // *
#define  HY_OP_CODE_ADD             (1+HY_OP_CODE_MUL)              // +
#define  HY_OP_CODE_SUB             (1+HY_OP_CODE_ADD)          // -
#define  HY_OP_CODE_DIV             (1+HY_OP_CODE_SUB)              // /
#define  HY_OP_CODE_LESS            (1+HY_OP_CODE_DIV)  // <

#define  HY_OP_CODE_LEQ             (1+HY_OP_CODE_LESS) // <=
#define  HY_OP_CODE_EQ              (1+HY_OP_CODE_LEQ) // ==
#define  HY_OP_CODE_GREATER         (1+HY_OP_CODE_EQ) // >
#define  HY_OP_CODE_GEQ             (1+HY_OP_CODE_GREATER) // >=
#define  HY_OP_CODE_ABS             (1+HY_OP_CODE_GEQ) // Abs
#define  HY_OP_CODE_ARCTAN          (1+HY_OP_CODE_ABS) // Arctan
#define  HY_OP_CODE_BETA            (1+HY_OP_CODE_ARCTAN) // Beta
#define  HY_OP_CODE_BRANCHCOUNT     (1+HY_OP_CODE_BETA) // BranchCount
#define  HY_OP_CODE_BRANCHLENGTH    (1+HY_OP_CODE_BRANCHCOUNT) // BranchLength
#define  HY_OP_CODE_BRANCHNAME      (1+HY_OP_CODE_BRANCHLENGTH) // BranchName

#define  HY_OP_CODE_CCHI2           (1+HY_OP_CODE_BRANCHNAME) // CChi2
#define  HY_OP_CODE_CGAMMADIST      (1+HY_OP_CODE_CCHI2) // CGammaDist
#define  HY_OP_CODE_COLUMNS         (1+HY_OP_CODE_CGAMMADIST) // Columns
#define  HY_OP_CODE_COS             (1+HY_OP_CODE_COLUMNS) // Cos
#define  HY_OP_CODE_DIFF            (1+HY_OP_CODE_COS) // D
#define  HY_OP_CODE_EIGENSYSTEM     (1+HY_OP_CODE_DIFF) // Eigensystem
#define  HY_OP_CODE_ERF             (1+HY_OP_CODE_EIGENSYSTEM) // Erf
#define  HY_OP_CODE_EVAL            (1+HY_OP_CODE_ERF) // Eval
#define  HY_OP_CODE_EXP             (1+HY_OP_CODE_EVAL) // Exp
#define  HY_OP_CODE_FORMAT          (1+HY_OP_CODE_EXP) // Format
#define  HY_OP_CODE_GAMMA           (1+HY_OP_CODE_FORMAT) // Gamma

#define  HY_OP_CODE_GAMMADIST       (1+HY_OP_CODE_GAMMA) // GammaDist
#define  HY_OP_CODE_IBETA           (1+HY_OP_CODE_GAMMADIST) // IBeta
#define  HY_OP_CODE_IGAMMA          (1+HY_OP_CODE_IBETA) // IGamma
#define  HY_OP_CODE_INVCHI2         (1+HY_OP_CODE_IGAMMA) // InvChi2
#define  HY_OP_CODE_INVERSE         (1+HY_OP_CODE_INVCHI2) // Inverse
#define  HY_OP_CODE_JOIN            (1+HY_OP_CODE_INVERSE) // Join
#define  HY_OP_CODE_LUDECOMPOSE     (1+HY_OP_CODE_JOIN) // LUDecompose
#define  HY_OP_CODE_LUSOLVE         (1+HY_OP_CODE_LUDECOMPOSE) // LUSolve
#define  HY_OP_CODE_LNGAMMA         (1+HY_OP_CODE_LUSOLVE) // LnGamma
#define  HY_OP_CODE_LOG             (1+HY_OP_CODE_LNGAMMA) // Log
#define  HY_OP_CODE_MACCESS         (1+HY_OP_CODE_LOG) // MAccess
#define  HY_OP_CODE_MCOORD          (1+HY_OP_CODE_MACCESS) // MCoord

#define  HY_OP_CODE_MAX             (1+HY_OP_CODE_MCOORD) // Max
#define  HY_OP_CODE_MIN             (1+HY_OP_CODE_MAX) // Min
#define  HY_OP_CODE_PSTREESTRING    (1+HY_OP_CODE_MIN) // PSTreeString
#define  HY_OP_CODE_RANDOM          (1+HY_OP_CODE_PSTREESTRING) // Random
#define  HY_OP_CODE_REROOTTREE      (1+HY_OP_CODE_RANDOM) // RerootTree
#define  HY_OP_CODE_ROWS            (1+HY_OP_CODE_REROOTTREE) // Rows
#define  HY_OP_CODE_SIMPLEX         (1+HY_OP_CODE_ROWS) // Simplex
#define  HY_OP_CODE_SIN             (1+HY_OP_CODE_SIMPLEX) // Sin
#define  HY_OP_CODE_SQRT            (1+HY_OP_CODE_SIN) // Sqrt
#define  HY_OP_CODE_TEXTREESTRING   (1+HY_OP_CODE_SQRT) // TEXTreeString

#define  HY_OP_CODE_TAN             (1+HY_OP_CODE_TEXTREESTRING) // Tan
#define  HY_OP_CODE_TIME            (1+HY_OP_CODE_TAN) // Time
#define  HY_OP_CODE_TIPCOUNT        (1+HY_OP_CODE_TIME) // TipCount
#define  HY_OP_CODE_TIPNAME         (1+HY_OP_CODE_TIPCOUNT) // TipName
#define  HY_OP_CODE_TRANSPOSE       (1+HY_OP_CODE_TIPNAME) // Transpose
#define  HY_OP_CODE_TYPE            (1+HY_OP_CODE_TRANSPOSE) // Type
#define  HY_OP_CODE_ZCDF            (1+HY_OP_CODE_TYPE) // ZCDF
#define  HY_OP_CODE_POWER           (1+HY_OP_CODE_ZCDF) // ^
#define  HY_OP_CODE_OR              (1+HY_OP_CODE_POWER) // ||

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

//Some useful defines
#define  HY_INVALID_RETURN_VALUE                        NAN
#define  HY_CONSTANT_FALSE                              _Constant (0.0);
#define  HY_CONSTANT_TRUE                               _Constant (1.0);

#define   BL_FUNCTION_ALWAYS_UPDATE     0
#define   BL_FUNCTION_NORMAL_UPDATE     1

#endif
