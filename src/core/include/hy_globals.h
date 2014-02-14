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

#ifndef __HY_GLOBALS__
#define __HY_GLOBALS__

#include "batchlan.h"
#include "legacy_parser.h"

extern bool dropIntoDebugMode, isInFunction, terminateExecution,
    skipWarningMessages;

;

extern char isDefiningATree;

extern FILE *globalErrorFile, *globalMessageFile;

extern _String blockWiseMatrix, errorFileName, messageFileName, randomSeed,
    scanfLastFilePath, siteWiseMatrix;

extern _SimpleList freeSlots;

extern long globalRandSeed, systemCPUCount;

extern _Parameter dFPrintFormat, dFDefaultWidth;

extern _Variable *_n_, *_x_;

extern _AVLListX openFileHandles, _builtInArgumentCounts;

#ifdef __HYPHYMPI__
extern int _hy_mpi_node_rank;
#endif

/* operation list */

struct _hy_op_definition {
  long  op_code;
  const char * op_string;
  long  argument_count; 
  void  * fast_exec;
} const _hyBuiltInOperationProperties [HY_OP_COUNT] = { 
  {HY_OP_CODE_NOT,"!",-1L, nil},
  {HY_OP_CODE_NEQ,"!=",-1L, nil},
  {HY_OP_CODE_IDIV,"$",-1L, nil},
  {HY_OP_CODE_MOD,"%",-1L, nil},
  {HY_OP_CODE_REF,"&",-1L, nil},
  {HY_OP_CODE_AND,"&&",-1L,(Ptr)AndNumbers},
  {HY_OP_CODE_MUL,"*",-1L, (Ptr)MultNumbers},
  {HY_OP_CODE_ADD,"+",-1L, (Ptr)AddNumbers},
  {HY_OP_CODE_SUB,"-",-1L, (Ptr)SubNumbers},
  {HY_OP_CODE_DIV,"/",-1L, (Ptr)DivNumbers},
  {HY_OP_CODE_LESS,"<",-1L,(Ptr)LessThan},
  {HY_OP_CODE_LEQ,"<=",-1L, (Ptr)LessThanE},
  {HY_OP_CODE_EQ,"==",-1L, (Ptr)EqualNumbers},
  {HY_OP_CODE_GREATER,">",-1L, (Ptr)GreaterThan},
  {HY_OP_CODE_GEQ,">=",-1L, (Ptr)GreaterThanE},
  {HY_OP_CODE_ABS,"Abs",1L, (Ptr)AbsNumber},
  {HY_OP_CODE_ARCTAN,"Arctan",1L, nil},
  {HY_OP_CODE_BETA,"Beta",2L, nil},
  {HY_OP_CODE_BRANCHCOUNT,"BranchCount",2L, nil},
  {HY_OP_CODE_BRANCHLENGTH,"BranchLength",2L, nil},
  {HY_OP_CODE_BRANCHNAME,"BranchName",2L, nil},
  {HY_OP_CODE_CCHI2,"CChi2",2L, nil},
  {HY_OP_CODE_CGAMMADIST,"CGammaDist",3L, nil},
  {HY_OP_CODE_COLUMNS,"Columns",1L, nil},
  {HY_OP_CODE_COS,"Cos",1L, nil},
  {HY_OP_CODE_DIFF,"Differentiate",2L, nil},
  {HY_OP_CODE_EIGENSYSTEM,"Eigensystem",1L, nil},
  {HY_OP_CODE_ERF,"Erf",1L, nil},
  {HY_OP_CODE_EVAL,"Eval",1L, nil},
  {HY_OP_CODE_EXP,"Exp",1L, (Ptr)ExpNumbers},
  {HY_OP_CODE_FORMAT,"Format",3L, nil},
  {HY_OP_CODE_GAMMA,"Gamma",1L, nil},
  {HY_OP_CODE_GAMMADIST,"GammaDist",3L, nil},
  {HY_OP_CODE_IBETA,"IBeta",3L, nil},
  {HY_OP_CODE_IGAMMA,"IGamma",2L, nil},
  {HY_OP_CODE_INVCHI2,"InvChi2",2L, nil},
  {HY_OP_CODE_INVERSE,"Inverse",1L, nil},
  {HY_OP_CODE_JOIN,"Join",2L, nil},
  {HY_OP_CODE_LUDECOMPOSE,"LUDecompose",1L, nil},
  {HY_OP_CODE_LUSOLVE,"LUSolve",2L, nil},
  {HY_OP_CODE_LNGAMMA,"LnGamma",1L, nil},
  {HY_OP_CODE_LOG,"Log",1L, (Ptr)LogNumbers},
  {HY_OP_CODE_MACCESS,"MAccess",3L, (Ptr)FastMxAccess},
  {HY_OP_CODE_MCOORD,"MCoord",3L, nil},
  {HY_OP_CODE_MSTORE,"MStore",3L, nil},
  {HY_OP_CODE_MAX,"Max",2L, (Ptr)MaxNumbers},
  {HY_OP_CODE_MIN,"Min",2L, (Ptr)MinNumbers},
  {HY_OP_CODE_PSTREESTRING,"PSTreeString",3L, nil},
  {HY_OP_CODE_RANDOM,"Random",2L, (Ptr) RandomNumber},
  {HY_OP_CODE_REROOTTREE,"RerootTree",2L, nil},
  {HY_OP_CODE_ROWS,"Rows",1L, nil},
  {HY_OP_CODE_SIMPLEX,"Simplex",1L, nil},
  {HY_OP_CODE_SIN,"Sin",1L, nil},
  {HY_OP_CODE_SQRT,"Sqrt",1L, nil},
  {HY_OP_CODE_TEXTREESTRING,"TEXTreeString",2L, nil},
  {HY_OP_CODE_TAN,"Tan",1L, nil},
  {HY_OP_CODE_TIME,"Time",1L, nil},
  {HY_OP_CODE_TIPCOUNT,"TipCount",1L, nil},
  {HY_OP_CODE_TIPNAME,"TipName",2L, nil},
  {HY_OP_CODE_TRANSPOSE,"Transpose",1L, nil},
  {HY_OP_CODE_TYPE,"Type",1L, nil},
  {HY_OP_CODE_ZCDF,"ZCDF",1L, nil},
  {HY_OP_CODE_POWER,"^",-1L, (Ptr) Power},
  {HY_OP_CODE_OR,"||",-1L, (Ptr)OrNumbers},
  {HY_OP_CODE_FSTORE,"FStore",3L, nil}
};

#endif
