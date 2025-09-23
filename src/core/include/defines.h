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

#ifndef __HY_DEFINES__
#define __HY_DEFINES__

/**
 * @file defines.h
 * @brief Global definitions for the HyPhy project.
 *
 * This file contains a large number of #define macros that are used
 * for constants, flags, and opcodes throughout the HyPhy codebase.
 */

#include <math.h>

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef nil
#define nil NULL
#endif

/** @name Compiler-specific macros */
///@{
#ifdef __GNUC__
/** @brief Restrict keyword for GCC */
#define _hprestrict_ __restrict
#else
/** @brief Restrict keyword for other compilers */
#define _hprestrict_
#endif
///@}

/** @name General constants */
///@{
/** @brief Width of a long in bits */
#define HY_WIDTH_OF_LONG ((long)(sizeof(long) * 8))

/** @brief Format string for printing floats */
#define PRINTF_FORMAT_STRING "%.16g"

/** @brief Default 0 bit flags */
#define fNone 0L

/** @brief Not found value */
#define kNotFound (-1L)

/** @brief Maximum number of dialog prompts */
#define kMaxDialogPrompts 10L
///@}

/** @name Object type codes */
///@{
/** @brief Undefined object type */
#define HY_UNDEFINED 0x400UL
/** @brief Number object type */
#define NUMBER 0x001UL
/** @brief Matrix object type */
#define MATRIX 0x004UL
/** @brief Container object type */
#define CONTAINER 0x008UL
/** @brief Tree node object type */
#define TREE_NODE 0x010UL
/** @brief Tree object type */
#define TREE 0x020UL
/** @brief String object type */
#define STRING 0x040UL
/** @brief Associative list object type */
#define ASSOCIATIVE_LIST 0x080UL
/** @brief Topology object type */
#define TOPOLOGY 0x100UL
/** @brief Polynomial object type */
#define POLYNOMIAL 0x200UL

/** @brief Any object type */
#define HY_ANY_OBJECT 0xFFFFUL
/** @brief Mutable object types */
#define HY_MUTABLE_OBJECT (NUMBER | STRING | HY_UNDEFINED)
///@}

/** @name Default bounds */
///@{
/** @brief Default lower bound */
#define DEFAULTLOWERBOUND -1e26
/** @brief Default upper bound */
#define DEFAULTUPPERBOUND 1e26
///@}

/** @name Variable types */
///@{
/** @brief Global variable type */
#define GLOBAL_VARIABLE 1
/** @brief Category variable type */
#define CATEGORY_VARIABLE 2
/** @brief Random variable type */
#define RANDOM_VARIABLE 3
///@}

/** @name Variable flags */
///@{
/** @brief Global variable flag */
#define HY_VARIABLE_GLOBAL 0x0001
/** @brief Variable changed flag */
#define HY_VARIABLE_CHANGED 0x0002
/** @brief Clear variable changed flag */
#define HY_VARIABLE_CHANGED_CLEAR 0xFFFD
/** @brief Dependent variable computed flag */
#define HY_DEP_V_COMPUTED 0x0004
/** @brief Clear dependent variable computed flag */
#define HY_DEP_V_COMPUTED_CLEAR 0xFFFB
/** @brief Dependent variable inspected flag */
#define HY_DEP_V_INSPECTED 0x0008
/** @brief Dependent variable modified flag */
#define HY_DEP_V_MODIFIED 0x0010
/** @brief Clear dependent variable modified flag */
#define HY_DEP_V_MODIFIED_CLEAR 0xFFEF
/** @brief Dependent variable modified categories flag */
#define HY_DEP_V_MODIFIED_CATS 0x0020
/** @brief Do not check this variable container in NeedToExponentiate */
#define HY_VC_NO_CHECK 0x0040
/** @brief Variable not set flag */
#define HY_VARIABLE_NOTSET 0x0080
/** @brief Variable computing flag */
#define HY_VARIABLE_COMPUTING 0x0100
/** @brief Variable set externally flag */
#define HY_VARIABLE_SET_EXTERNALLY 0x0200

/** @brief Set variable flag */
#define HY_VARIABLE_SET 0xFF7F
/** @brief Clear dependent variable inspected flag */
#define HY_DEP_V_INSPECTED_CLR 0xFFF7
/** @brief Clear variable computing flag */
#define HY_VARIABLE_COMPUTING_CLR 0xFEFF
/** @brief Clear variable set externally flag */
#define HY_VARIABLE_SET_EXTERNALLY_CLR 0xFDFF
/** @brief Clear no check variable container flag */
#define HY_VC_CLR_NO_CHECK 0xFFBF
/** @brief Clear dependent variable mask */
#define HY_DEP_CLEAR_MASK 0xFFC7
/** @brief Clear variable changed flag */
#define HY_HY_VARIABLE_CHANGED_CLEAR 0xFFFC

/** @brief Variable partition flag */
#define HY_VARIABLE_PARTITION 0x1000
///@}

/** @brief No model constant */
#define HY_NO_MODEL (-1)

/** @name Mathematical constants */
///@{
/** @brief Golden ratio */
#define GOLDEN_RATIO 1.618034
/** @brief Maximum number of Brent iterates */
#define MAX_BRENT_ITERATES 100L
///@}

/** @name Opcodes */
///@{
/** @brief No operation */
#define HY_OP_CODE_NONE -1

/** @brief Logical NOT */
#define HY_OP_CODE_NOT 0
/** @brief Not equal */
#define HY_OP_CODE_NEQ (1 + HY_OP_CODE_NOT)
/** @brief Integer division */
#define HY_OP_CODE_IDIV (1 + HY_OP_CODE_NEQ)
/** @brief Modulo */
#define HY_OP_CODE_MOD (1 + HY_OP_CODE_IDIV)
/** @brief Reference */
#define HY_OP_CODE_REF (1 + HY_OP_CODE_MOD)
/** @brief Logical AND */
#define HY_OP_CODE_AND (1 + HY_OP_CODE_REF)
/** @brief Multiplication */
#define HY_OP_CODE_MUL (1 + HY_OP_CODE_AND)
/** @brief Addition */
#define HY_OP_CODE_ADD (1 + HY_OP_CODE_MUL)
/** @brief Subtraction */
#define HY_OP_CODE_SUB (1 + HY_OP_CODE_ADD)
/** @brief Division */
#define HY_OP_CODE_DIV (1 + HY_OP_CODE_SUB)
/** @brief Less than */
#define HY_OP_CODE_LESS (1 + HY_OP_CODE_DIV)

/** @brief Less than or equal */
#define HY_OP_CODE_LEQ (1 + HY_OP_CODE_LESS)
/** @brief Equal */
#define HY_OP_CODE_EQ (1 + HY_OP_CODE_LEQ)
/** @brief Greater than */
#define HY_OP_CODE_GREATER (1 + HY_OP_CODE_EQ)
/** @brief Greater than or equal */
#define HY_OP_CODE_GEQ (1 + HY_OP_CODE_GREATER)
/** @brief Absolute value */
#define HY_OP_CODE_ABS (1 + HY_OP_CODE_GEQ)
/** @brief Arctangent */
#define HY_OP_CODE_ARCTAN (1 + HY_OP_CODE_ABS)
/** @brief Beta function */
#define HY_OP_CODE_BETA (1 + HY_OP_CODE_ARCTAN)
/** @brief Branch count */
#define HY_OP_CODE_BRANCHCOUNT (1 + HY_OP_CODE_BETA)
/** @brief Branch length */
#define HY_OP_CODE_BRANCHLENGTH (1 + HY_OP_CODE_BRANCHCOUNT)
/** @brief Branch name */
#define HY_OP_CODE_BRANCHNAME (1 + HY_OP_CODE_BRANCHLENGTH)

/** @brief Cumulative chi-squared distribution */
#define HY_OP_CODE_CCHI2 (1 + HY_OP_CODE_BRANCHNAME)
/** @brief Cumulative gamma distribution */
#define HY_OP_CODE_CGAMMADIST (1 + HY_OP_CODE_CCHI2)
/** @brief Call a function */
#define HY_OP_CODE_CALL (1 + HY_OP_CODE_CGAMMADIST)
/** @brief Number of columns */
#define HY_OP_CODE_COLUMNS (1 + HY_OP_CODE_CALL)
/** @brief Cosine */
#define HY_OP_CODE_COS (1 + HY_OP_CODE_COLUMNS)
/** @brief Differentiation */
#define HY_OP_CODE_DIFF (1 + HY_OP_CODE_COS)
/** @brief Eigensystem */
#define HY_OP_CODE_EIGENSYSTEM (1 + HY_OP_CODE_DIFF)
/** @brief Error function */
#define HY_OP_CODE_ERF (1 + HY_OP_CODE_EIGENSYSTEM)
/** @brief Evaluate an expression */
#define HY_OP_CODE_EVAL (1 + HY_OP_CODE_ERF)
/** @brief Exponential */
#define HY_OP_CODE_EXP (1 + HY_OP_CODE_EVAL)
/** @brief Format a string */
#define HY_OP_CODE_FORMAT (1 + HY_OP_CODE_EXP)
/** @brief Gamma function */
#define HY_OP_CODE_GAMMA (1 + HY_OP_CODE_FORMAT)

/** @brief Gamma distribution */
#define HY_OP_CODE_GAMMADIST (1 + HY_OP_CODE_GAMMA)
/** @brief Incomplete beta function */
#define HY_OP_CODE_IBETA (1 + HY_OP_CODE_GAMMADIST)
/** @brief Incomplete gamma function */
#define HY_OP_CODE_IGAMMA (1 + HY_OP_CODE_IBETA)
/** @brief Inverse chi-squared distribution */
#define HY_OP_CODE_INVCHI2 (1 + HY_OP_CODE_IGAMMA)
/** @brief Inverse of a matrix */
#define HY_OP_CODE_INVERSE (1 + HY_OP_CODE_INVCHI2)
/** @brief Join strings */
#define HY_OP_CODE_JOIN (1 + HY_OP_CODE_INVERSE)
/** @brief LU decomposition */
#define HY_OP_CODE_LUDECOMPOSE (1 + HY_OP_CODE_JOIN)
/** @brief LU solve */
#define HY_OP_CODE_LUSOLVE (1 + HY_OP_CODE_LUDECOMPOSE)
/** @brief Natural logarithm of the gamma function */
#define HY_OP_CODE_LNGAMMA (1 + HY_OP_CODE_LUSOLVE)
/** @brief Natural logarithm */
#define HY_OP_CODE_LOG (1 + HY_OP_CODE_LNGAMMA)
/** @brief Matrix access */
#define HY_OP_CODE_MACCESS (1 + HY_OP_CODE_LOG)
/** @brief Matrix coordinate */
#define HY_OP_CODE_MCOORD (1 + HY_OP_CODE_MACCESS)

/** @brief Maximum */
#define HY_OP_CODE_MAX (1 + HY_OP_CODE_MCOORD)
/** @brief Minimum */
#define HY_OP_CODE_MIN (1 + HY_OP_CODE_MAX)
/** @brief PostScript tree string */
#define HY_OP_CODE_PSTREESTRING (1 + HY_OP_CODE_MIN)
/** @brief Polynomial */
#define HY_OP_CODE_POLYNOMIAL (1 + HY_OP_CODE_PSTREESTRING)
/** @brief Random number */
#define HY_OP_CODE_RANDOM (1 + HY_OP_CODE_POLYNOMIAL)
/** @brief Reroot a tree */
#define HY_OP_CODE_REROOTTREE (1 + HY_OP_CODE_RANDOM)
/** @brief Number of rows */
#define HY_OP_CODE_ROWS (1 + HY_OP_CODE_REROOTTREE)
/** @brief Simplex optimization */
#define HY_OP_CODE_SIMPLEX (1 + HY_OP_CODE_ROWS)
/** @brief Simplify an expression */
#define HY_OP_CODE_SIMPLIFY (1 + HY_OP_CODE_SIMPLEX)
/** @brief Sine */
#define HY_OP_CODE_SIN (1 + HY_OP_CODE_SIMPLIFY)
/** @brief Square root */
#define HY_OP_CODE_SQRT (1 + HY_OP_CODE_SIN)
/** @brief Text tree string */
#define HY_OP_CODE_TEXTREESTRING (1 + HY_OP_CODE_SQRT)

/** @brief Tangent */
#define HY_OP_CODE_TAN (1 + HY_OP_CODE_TEXTREESTRING)
/** @brief Time */
#define HY_OP_CODE_TIME (1 + HY_OP_CODE_TAN)
/** @brief Tip count */
#define HY_OP_CODE_TIPCOUNT (1 + HY_OP_CODE_TIME)
/** @brief Tip name */
#define HY_OP_CODE_TIPNAME (1 + HY_OP_CODE_TIPCOUNT)
/** @brief Transpose a matrix */
#define HY_OP_CODE_TRANSPOSE (1 + HY_OP_CODE_TIPNAME)
/** @brief Type of an object */
#define HY_OP_CODE_TYPE (1 + HY_OP_CODE_TRANSPOSE)
/** @brief Z-score cumulative distribution function */
#define HY_OP_CODE_ZCDF (1 + HY_OP_CODE_TYPE)
/** @brief Power */
#define HY_OP_CODE_POWER (1 + HY_OP_CODE_ZCDF)
/** @brief Logical OR */
#define HY_OP_CODE_OR (1 + HY_OP_CODE_POWER)

///@}

/** @name Formula return codes */
///@{
/** @brief Expression return code */
#define HY_FORMULA_EXPRESSION 0
/** @brief Failed return code */
#define HY_FORMULA_FAILED (-1)
/** @brief Variable value assignment return code */
#define HY_FORMULA_VARIABLE_VALUE_ASSIGNMENT 1
/** @brief Variable formula assignment return code */
#define HY_FORMULA_VARIABLE_FORMULA_ASSIGNMENT 2
/** @brief Formula formula assignment return code */
#define HY_FORMULA_FORMULA_FORMULA_ASSIGNMENT 3
/** @brief Formula value assignment return code */
#define HY_FORMULA_FORMULA_VALUE_ASSIGNMENT 4
/** @brief Variable lower bound assignment return code */
#define HY_FORMULA_VARIABLE_LOWER_BOUND_ASSIGNMENT 5
/** @brief Variable upper bound assignment return code */
#define HY_FORMULA_VARIABLE_UPPER_BOUND_ASSIGNMENT 6
/** @brief Formula value increment return code */
#define HY_FORMULA_FORMULA_VALUE_INCREMENT 7
/** @brief Reference value assignment return code */
#define HY_FORMULA_REFERENCE_VALUE_ASSIGNMENT 20
/** @brief Reference formula assignment return code */
#define HY_FORMULA_REFERENCE_FORMULA_ASSIGNMENT 30

/** @brief Reference lower bound assignment return code */
#define HY_FORMULA_REFERENCE_LOWER_BOUND_ASSIGNMENT 40
/** @brief Reference upper bound assignment return code */
#define HY_FORMULA_REFERENCE_UPPER_BOUND_ASSIGNMENT 41
///@}

/** @name Useful defines */
///@{
/** @brief Invalid return value */
#define HY_INVALID_RETURN_VALUE NAN
/** @brief Constant false value */
#define HY_CONSTANT_FALSE _Constant(0.0)
/** @brief Constant true value */
#define HY_CONSTANT_TRUE _Constant(1.0)
/** @brief Null return value */
#define HY_NULL_RETURN _MathObject()
///@}

/** @name Batch Language 'Object' type codes
    @{
*/

/** @brief Not defined */
#define HY_BL_NOT_DEFINED 0
/** @brief Dataset object */
#define HY_BL_DATASET 1
/** @brief Dataset filter object */
#define HY_BL_DATASET_FILTER 2
/** @brief Likelihood function object */
#define HY_BL_LIKELIHOOD_FUNCTION 4
/** @brief SCFG object */
#define HY_BL_SCFG 8
/** @brief BGM object */
#define HY_BL_BGM 16
/** @brief Model object */
#define HY_BL_MODEL 32
/** @brief HBL function object */
#define HY_BL_HBL_FUNCTION 64
/** @brief Tree object */
#define HY_BL_TREE 128
/** @brief Variable object */
#define HY_BL_VARIABLE 256

/** @brief Any object */
#define HY_BL_ANY 65535
///@}

/** @name Batch Language Command Codes
    @{
*/

/** @brief Formula command */
#define HY_HBL_COMMAND_FORMULA 0L
/** @brief Conditional command */
#define HY_HBL_COMMAND_CONDITIONAL 4L

/** @brief For loop command */
#define HY_HBL_COMMAND_FOR 500L
/** @brief While loop command */
#define HY_HBL_COMMAND_WHILE 501L
/** @brief Function definition command */
#define HY_HBL_COMMAND_FUNCTION 502L
/** @brief Function definition command */
#define HY_HBL_COMMAND_FFUNCTION 503L
/** @brief Return from function command */
#define HY_HBL_COMMAND_RETURNSPACE 504L
/** @brief Return from function command */
#define HY_HBL_COMMAND_RETURNPAREN 505L
/** @brief If command */
#define HY_HBL_COMMAND_IF 506L
/** @brief Else command */
#define HY_HBL_COMMAND_ELSE 507L
/** @brief Do-while loop command */
#define HY_HBL_COMMAND_DO 508L
/** @brief Break command */
#define HY_HBL_COMMAND_BREAK 509L
/** @brief Continue command */
#define HY_HBL_COMMAND_CONTINUE 510L
/** @brief Include file command */
#define HY_HBL_COMMAND_INCLUDE 511L
/** @brief Dataset command */
#define HY_HBL_COMMAND_DATA_SET 512L
/** @brief Dataset filter command */
#define HY_HBL_COMMAND_DATA_SET_FILTER 513L
/** @brief Harvest frequencies command */
#define HY_HBL_COMMAND_HARVEST_FREQUENCIES 514L
/** @brief Construct category matrix command */
#define HY_HBL_COMMAND_CONSTRUCT_CATEGORY_MATRIX 515L
/** @brief Tree command */
#define HY_HBL_COMMAND_TREE 516L
/** @brief Likelihood function command */
#define HY_HBL_COMMAND_LIKELIHOOD_FUNCTION 517L
/** @brief Likelihood function 3 command */
#define HY_HBL_COMMAND_LIKELIHOOD_FUNCTION_3 518L
/** @brief Optimize command */
#define HY_HBL_COMMAND_OPTIMIZE 519L
/** @brief Covariance matrix command */
#define HY_HBL_COMMAND_COVARIANCE_MATRIX 520L
/** @brief Molecular clock command */
#define HY_HBL_COMMAND_MOLECULAR_CLOCK 521L
/** @brief Fprintf command */
#define HY_HBL_COMMAND_FPRINTF 522L
/** @brief Fscanf command */
#define HY_HBL_COMMAND_FSCANF 523L
/** @brief Sscanf command */
#define HY_HBL_COMMAND_SSCANF 524L
/** @brief Get string command */
#define HY_HBL_COMMAND_GET_STRING 525L
/** @brief Export command */
#define HY_HBL_COMMAND_EXPORT 526L
/** @brief Import command */
#define HY_HBL_COMMAND_IMPORT 527L
/** @brief Category command */
#define HY_HBL_COMMAND_CATEGORY 528L
/** @brief Clear constraints command */
#define HY_HBL_COMMAND_CLEAR_CONSTRAINTS 529L
/** @brief Set dialog prompt command */
#define HY_HBL_COMMAND_SET_DIALOG_PROMPT 530L
/** @brief Select template model command */
#define HY_HBL_COMMAND_SELECT_TEMPLATE_MODEL 531L
/** @brief Use model command */
#define HY_HBL_COMMAND_USE_MODEL 532L
/** @brief Model command */
#define HY_HBL_COMMAND_MODEL 533L
/** @brief Set parameter command */
#define HY_HBL_COMMAND_SET_PARAMETER 534L
/** @brief Choice list command */
#define HY_HBL_COMMAND_CHOICE_LIST 535L
/** @brief Get information command */
#define HY_HBL_COMMAND_GET_INFORMATION 537L
/** @brief Execute commands command */
#define HY_HBL_COMMAND_EXECUTE_COMMANDS 538L
/** @brief Execute a file command */
#define HY_HBL_COMMAND_EXECUTE_A_FILE 539L
/** @brief Load function library command */
#define HY_HBL_COMMAND_LOAD_FUNCTION_LIBRARY 541L
/** @brief Differentiate command */
#define HY_HBL_COMMAND_DIFFERENTIATE 544L
/** @brief Find root command */
#define HY_HBL_COMMAND_FIND_ROOT 545L
/** @brief MPI receive command */
#define HY_HBL_COMMAND_MPI_RECEIVE 546L
/** @brief MPI send command */
#define HY_HBL_COMMAND_MPI_SEND 547L
/** @brief Get data info command */
#define HY_HBL_COMMAND_GET_DATA_INFO 548L
/** @brief State counter command */
#define HY_HBL_COMMAND_STATE_COUNTER 549L
/** @brief Integrate command */
#define HY_HBL_COMMAND_INTEGRATE 550L
/** @brief LFCompute command */
#define HY_HBL_COMMAND_LFCOMPUTE 551L
/** @brief Get URL command */
#define HY_HBL_COMMAND_GET_URL 552L
/** @brief Do SQL command */
#define HY_HBL_COMMAND_DO_SQL 553L
/** @brief Topology command */
#define HY_HBL_COMMAND_TOPOLOGY 554L
/** @brief Align sequences command */
#define HY_HBL_COMMAND_ALIGN_SEQUENCES 555L
/** @brief Get neutral null command */
#define HY_HBL_COMMAND_GET_NEUTRAL_NULL 556L
/** @brief Profile command */
#define HY_HBL_COMMAND_PROFILE 557L
/** @brief Delete object command */
#define HY_HBL_COMMAND_DELETE_OBJECT 558L
/** @brief Require version command */
#define HY_HBL_COMMAND_REQUIRE_VERSION 559L
/** @brief SCFG command */
#define HY_HBL_COMMAND_SCFG 560L
/** @brief Neural net command */
#define HY_HBL_COMMAND_NEURAL_NET 561L
/** @brief BGM command */
#define HY_HBL_COMMAND_BGM 562L
/** @brief Simulate dataset command */
#define HY_HBL_COMMAND_SIMULATE_DATA_SET 563L
/** @brief Assert command */
#define HY_HBL_COMMAND_ASSERT 564L
/** @brief Replicate constraint command */
#define HY_HBL_COMMAND_REPLICATE_CONSTRAINT 565L
/** @brief Nested list command */
#define HY_HBL_COMMAND_NESTED_LIST 566L
/** @brief Keyword argument command */
#define HY_HBL_COMMAND_KEYWORD_ARGUMENT 567L
/** @brief Init iterator command */
#define HY_HBL_COMMAND_INIT_ITERATOR 568L
/** @brief Advance iterator command */
#define HY_HBL_COMMAND_ADVANCE_ITERATOR 569L
/** @brief Convert branch length command */
#define HY_HBL_COMMAND_CONVERT_BRANCH_LENGTH 570L
///@}

/** @name HyPhy standard directory locations
    @{
*/

/** @brief Template models directory */
#define HY_HBL_DIRECTORY_TEMPLATE_MODELS 1000L
///@}

/** @brief Maximum long value */
#define HY_MAX_LONG_VALUE 0xffffffffL

/** @name BGM Get String codes */
///@{
/** @brief BGM score */
#define HY_HBL_GET_STRING_BGM_SCORE 0L
/** @brief BGM serialize */
#define HY_HBL_GET_STRING_BGM_SERIALIZE 1L
///@}

// TODO 20170413: not sure where to put the conditional includes below
#ifdef _SLKP_USE_ARM_NEON
#include <arm_neon.h>
#endif

#ifdef _SLKP_USE_SSE_INTRINSICS
// #include "sse2neon.h"
// #else
#include <pmmintrin.h>
// #endif
#endif

#ifdef _SLKP_USE_AVX_INTRINSICS
#include <immintrin.h>
#endif

#endif
