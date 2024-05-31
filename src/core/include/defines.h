/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-now
Core Developers:
  Sergei L Kosakovsky Pond (sergeilkp@mac.com)
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


#ifndef     __HY_DEFINES__
#define     __HY_DEFINES__

#include    <math.h>

#ifndef FALSE
    #define FALSE 0
#endif

#ifndef TRUE
    #define TRUE  1
#endif

#ifndef nil
#define nil   NULL
#endif

#ifdef __GNUC__
#define _hprestrict_ __restrict
#else
#define _hprestrict_
#endif

#define  HY_WIDTH_OF_LONG          ((long)(sizeof(long)*8))

#define  PRINTF_FORMAT_STRING    "%.16g"

#define  fNone              0L      // default 0 bit flags

#define  kNotFound          (-1L)

#define  kMaxDialogPrompts  10L

#define  HY_UNDEFINED       0x400UL
#define  NUMBER             0x001UL
#define  MATRIX             0x004UL
#define  CONTAINER          0x008UL
#define  TREE_NODE          0x010UL
#define  TREE               0x020UL
#define  STRING             0x040UL
#define  ASSOCIATIVE_LIST   0x080UL
#define  TOPOLOGY           0x100UL
#define  POLYNOMIAL         0x200UL

#define  HY_ANY_OBJECT      0xFFFFUL
#define  HY_MUTABLE_OBJECT  (NUMBER | STRING | HY_UNDEFINED)

#define  DEFAULTLOWERBOUND  -1e26
#define  DEFAULTUPPERBOUND  1e26

#define  GLOBAL_VARIABLE   1
#define  CATEGORY_VARIABLE 2
#define  RANDOM_VARIABLE   3

#define  HY_VARIABLE_GLOBAL                         0x0001
#define  HY_VARIABLE_CHANGED                        0x0002
#define  HY_VARIABLE_CHANGED_CLEAR                  0xFFFD
#define  HY_DEP_V_COMPUTED                          0x0004
#define  HY_DEP_V_COMPUTED_CLEAR                    0xFFFB
#define  HY_DEP_V_INSPECTED                         0x0008
#define  HY_DEP_V_MODIFIED                          0x0010
#define  HY_DEP_V_MODIFIED_CLEAR                    0xFFEF
#define  HY_DEP_V_MODIFIED_CATS                     0x0020
#define  HY_VC_NO_CHECK                             0x0040 // do not check this variable container in
// NeedToExponentiate
#define  HY_VARIABLE_NOTSET                         0x0080
#define  HY_VARIABLE_COMPUTING                      0x0100
#define  HY_VARIABLE_SET_EXTERNALLY                 0x0200

#define  HY_VARIABLE_SET                            0xFF7F
#define  HY_DEP_V_INSPECTED_CLR                     0xFFF7
#define  HY_VARIABLE_COMPUTING_CLR                  0xFEFF
#define  HY_VARIABLE_SET_EXTERNALLY_CLR             0xFDFF
#define  HY_VC_CLR_NO_CHECK                         0xFFBF
#define  HY_DEP_CLEAR_MASK                          0xFFC7
#define  HY_HY_VARIABLE_CHANGED_CLEAR               0xFFFC

#define  HY_VARIABLE_PARTITION                      0x1000

#define  HY_NO_MODEL            (-1)

#define  GOLDEN_RATIO 1.618034
#define MAX_BRENT_ITERATES 100L

// START OPCODES

#define  HY_OP_CODE_NONE            -1                              // noop

#define  HY_OP_CODE_NOT             0                               // !
#define  HY_OP_CODE_NEQ             (1+HY_OP_CODE_NOT)              // !=
#define  HY_OP_CODE_IDIV            (1+HY_OP_CODE_NEQ)              // $
#define  HY_OP_CODE_MOD             (1+HY_OP_CODE_IDIV)             // %
#define  HY_OP_CODE_REF             (1+HY_OP_CODE_MOD)              // &
#define  HY_OP_CODE_AND             (1+HY_OP_CODE_REF)              // &&
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
#define  HY_OP_CODE_CALL            (1+HY_OP_CODE_CGAMMADIST) // Call
#define  HY_OP_CODE_COLUMNS         (1+HY_OP_CODE_CALL) // Columns
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
#define  HY_OP_CODE_SIMPLIFY        (1+HY_OP_CODE_SIMPLEX) // Simplify
#define  HY_OP_CODE_SIN             (1+HY_OP_CODE_SIMPLIFY) // Sin
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
#define  HY_FORMULA_REFERENCE_VALUE_ASSIGNMENT          20
#define  HY_FORMULA_REFERENCE_FORMULA_ASSIGNMENT        30

#define  HY_FORMULA_REFERENCE_LOWER_BOUND_ASSIGNMENT    40
#define  HY_FORMULA_REFERENCE_UPPER_BOUND_ASSIGNMENT    41

// END FORMULA RETURN CODES

//Some useful defines
#define  HY_INVALID_RETURN_VALUE                        NAN
#define  HY_CONSTANT_FALSE                              _Constant (0.0)
#define  HY_CONSTANT_TRUE                               _Constant (1.0)
#define  HY_NULL_RETURN                                 _MathObject ()

//!  Batch Language 'Object' type codes
/*!
     20110608 SLKP introduced.

     Bit flag style type tags and masks
     This is primarily to be used for object retrieval
     using the _HYRetrieveBLObjectByName function
     and in the _HY_GetStringGlobalTypes list
*/

#define   HY_BL_NOT_DEFINED             0
#define   HY_BL_DATASET                 1
#define   HY_BL_DATASET_FILTER          2
#define   HY_BL_LIKELIHOOD_FUNCTION     4
#define   HY_BL_SCFG                    8
#define   HY_BL_BGM                     16
#define   HY_BL_MODEL                   32
#define   HY_BL_HBL_FUNCTION            64
#define   HY_BL_TREE                    128
#define   HY_BL_VARIABLE                256

#define   HY_BL_ANY                     65535

//!  Batch Lanuage Command Codes 
/*!
     20111222 SLKP introduced.
     Define tags to replace ExecuteCaseWhatever

*/

#define   HY_HBL_COMMAND_FORMULA                                        0L

#define   HY_HBL_COMMAND_FOR                                            500L
#define   HY_HBL_COMMAND_WHILE                                          501L
#define   HY_HBL_COMMAND_FUNCTION                                       502L
#define   HY_HBL_COMMAND_FFUNCTION                                      503L
#define   HY_HBL_COMMAND_RETURNSPACE                                    504L
#define   HY_HBL_COMMAND_RETURNPAREN                                    505L
#define   HY_HBL_COMMAND_IF                                             506L
#define   HY_HBL_COMMAND_ELSE                                           507L
#define   HY_HBL_COMMAND_DO                                             508L
#define   HY_HBL_COMMAND_BREAK                                          509L
#define   HY_HBL_COMMAND_CONTINUE                                       510L
#define   HY_HBL_COMMAND_INCLUDE                                        511L
#define   HY_HBL_COMMAND_DATA_SET                                       512L
#define   HY_HBL_COMMAND_DATA_SET_FILTER                                513L
#define   HY_HBL_COMMAND_HARVEST_FREQUENCIES                            514L
#define   HY_HBL_COMMAND_CONSTRUCT_CATEGORY_MATRIX                      515L
#define   HY_HBL_COMMAND_TREE                                           516L
#define   HY_HBL_COMMAND_LIKELIHOOD_FUNCTION                            517L
#define   HY_HBL_COMMAND_LIKELIHOOD_FUNCTION_3                          518L
#define   HY_HBL_COMMAND_OPTIMIZE                                       519L
#define   HY_HBL_COMMAND_COVARIANCE_MATRIX                              520L
#define   HY_HBL_COMMAND_MOLECULAR_CLOCK                                521L
#define   HY_HBL_COMMAND_FPRINTF                                        522L
#define   HY_HBL_COMMAND_FSCANF                                         523L
#define   HY_HBL_COMMAND_SSCANF                                         524L
#define   HY_HBL_COMMAND_GET_STRING                                     525L
#define   HY_HBL_COMMAND_EXPORT                                         526L
#define   HY_HBL_COMMAND_IMPORT                                         527L
#define   HY_HBL_COMMAND_CATEGORY                                       528L
#define   HY_HBL_COMMAND_CLEAR_CONSTRAINTS                              529L
#define   HY_HBL_COMMAND_SET_DIALOG_PROMPT                              530L
#define   HY_HBL_COMMAND_SELECT_TEMPLATE_MODEL                          531L
#define   HY_HBL_COMMAND_USE_MODEL                                      532L
#define   HY_HBL_COMMAND_MODEL                                          533L
#define   HY_HBL_COMMAND_SET_PARAMETER                                  534L
#define   HY_HBL_COMMAND_CHOICE_LIST                                    535L
#define   HY_HBL_COMMAND_GET_INFORMATION                                537L
#define   HY_HBL_COMMAND_EXECUTE_COMMANDS                               538L
#define   HY_HBL_COMMAND_EXECUTE_A_FILE                                 539L
#define   HY_HBL_COMMAND_LOAD_FUNCTION_LIBRARY                          541L
#define   HY_HBL_COMMAND_DIFFERENTIATE                                  544L
#define   HY_HBL_COMMAND_FIND_ROOT                                      545L
#define   HY_HBL_COMMAND_MPI_RECEIVE                                    546L
#define   HY_HBL_COMMAND_MPI_SEND                                       547L
#define   HY_HBL_COMMAND_GET_DATA_INFO                                  548L
#define   HY_HBL_COMMAND_STATE_COUNTER                                  549L
#define   HY_HBL_COMMAND_INTEGRATE                                      550L
#define   HY_HBL_COMMAND_LFCOMPUTE                                      551L
#define   HY_HBL_COMMAND_GET_URL                                        552L
#define   HY_HBL_COMMAND_DO_SQL                                         553L
#define   HY_HBL_COMMAND_TOPOLOGY                                       554L
#define   HY_HBL_COMMAND_ALIGN_SEQUENCES                                555L
#define   HY_HBL_COMMAND_GET_NEUTRAL_NULL                               556L
#define   HY_HBL_COMMAND_PROFILE                                        557L
#define   HY_HBL_COMMAND_DELETE_OBJECT                                  558L
#define   HY_HBL_COMMAND_REQUIRE_VERSION                                559L
#define   HY_HBL_COMMAND_SCFG                                           560L
#define   HY_HBL_COMMAND_NEURAL_NET                                     561L
#define   HY_HBL_COMMAND_BGM                                            562L
#define   HY_HBL_COMMAND_SIMULATE_DATA_SET                              563L
#define   HY_HBL_COMMAND_ASSERT                                         564L
#define   HY_HBL_COMMAND_REPLICATE_CONSTRAINT                           565L
#define   HY_HBL_COMMAND_NESTED_LIST                                    566L
#define   HY_HBL_COMMAND_KEYWORD_ARGUMENT                               567L
#define   HY_HBL_COMMAND_INIT_ITERATOR                                  568L
#define   HY_HBL_COMMAND_ADVANCE_ITERATOR                               569L
#define   HY_HBL_COMMAND_CONVERT_BRANCH_LENGTH                          570L


//!  HyPhy standard directory locations 
/*!
    retrieve using _HYStandardDirectory

*/

#define   HY_HBL_DIRECTORY_TEMPLATE_MODELS                              1000L

#define   HY_MAX_LONG_VALUE                                             0xffffffffL


#define   HY_HBL_GET_STRING_BGM_SCORE                                   0L
#define   HY_HBL_GET_STRING_BGM_SERIALIZE                               1L


//TODO 20170413: not sure where to put the conditional includes below
#ifdef _SLKP_USE_ARM_NEON
    #include <arm_neon.h>
#endif

#ifdef _SLKP_USE_SSE_INTRINSICS
    //#include "sse2neon.h"
    //#else
    #include <pmmintrin.h>
    //#endif
#endif

#ifdef _SLKP_USE_AVX_INTRINSICS
#include <immintrin.h>
#endif

#endif
