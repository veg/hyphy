/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2009
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon              (apoon@cfenet.ubc.ca)

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

#ifndef __CALCNODEGLOBAL__
#define __CALCNODEGLOBAL__

#include "math.h"
#include "float.h"
#include "ctype.h"

extern long *nonZeroNodes, nonZeroNodesDim;

#ifdef __MP__

#include <pthread.h>
struct ThreadMatrixTask {
  long cID, tcat, startAt, endAt;
  _SimpleList *updateCN;

};

#include <pthread.h>
//struct   ThreadMatrixTask {
//    long   cID,
//           tcat,
//           startAt,
//           endAt;
//    _SimpleList* updateCN;

//};
extern pthread_t *matrixThreads;
extern ThreadMatrixTask *matrixTasks;
extern pthread_mutex_t matrixMutex;

void *MatrixUpdateFunction(void *arg);
#endif

extern char takeBranchLengths, autoSolveBranchLengths;

extern _Parameter ignoringInternalNames;

extern _Parameter treeLayoutVert, scalingLogConstant;

extern _SimpleList convertedMatrixExpressionsL;
extern _AVLListX convertedMatrixExpressions;

extern _String expectedNumberOfSubs, stringSuppliedLengths, noInternalLabels,
    includeModelSpecs, acceptRootedTrees, acceptBranchLengths, autoConvertBL,
    splitNodeNames, internalNodePrefix, ignoreUserINames, cotNode, cotSplit,
    cotBranchLength, cotDistance, cotCDF, cotSamples, cotSampler, cotToNode,
    treeOutputAVL, treeOutputBackground, treeOutputRightMargin, treeOutputEmbed,
    treeOutputXtraMargin, treeOutputSplit, treeOutputNotchesColor,
    treeOutputNotches, treeOutputLabel, treeOutputTLabel, treeOutputColor,
    treeOutputThickness, treeOutputLinecap, treeOutputDash, treeOutputOLabel,
    treeOutputSymbols, treeOutputSymbolSize, treeOutputExtraPS,
    treeOutputPrefixPS, treeOutputLayout, treeOutputNNPlaceH,
    treeOutputFSPlaceH, largeMatrixBranchLengthDimension,
    largeMatrixBranchLength, newNodeGraftName, newNodeGraftWhere,
    newNodeGraftParent, eqWithReroot, eqWithoutReroot, iNodePrefix;

extern _Parameter _timesCharWidths[256], _maxTimesCharWidth;


inline void _handle4x4_pruning_case(double *childVector, double *tMatrix,
                                    double *parentConditionals) {
#ifdef _SLKP_USE_SSE_INTRINSICS
  double tv[4]
      __attribute__((aligned(16))) = { childVector[0], childVector[1],
                                       childVector[2], childVector[3] };

  __m128d buffer0 = _mm_loadu_pd(tv), buffer1 = _mm_loadu_pd(tv + 2),
          matrix01 = _mm_loadu_pd(tMatrix),
          matrix12 = _mm_loadu_pd(tMatrix + 2),
          matrix34 = _mm_loadu_pd(tMatrix + 4),
          matrix56 = _mm_loadu_pd(tMatrix + 6),
          reg_storage = _mm_mul_pd(buffer0, matrix01),
          reg_storage2 = _mm_mul_pd(buffer0, matrix34);

  matrix34 = _mm_mul_pd(buffer1, matrix12);
  matrix56 = _mm_mul_pd(buffer1, matrix56);
  reg_storage = _mm_add_pd(reg_storage, matrix34);
  reg_storage2 = _mm_add_pd(reg_storage2, matrix56);
  reg_storage = _mm_hadd_pd(reg_storage, reg_storage2);
  matrix01 = _mm_loadu_pd(parentConditionals);
  matrix01 = _mm_mul_pd(reg_storage, matrix01);
  _mm_storeu_pd(parentConditionals, matrix01);

  matrix01 = _mm_loadu_pd(tMatrix + 8);
  matrix12 = _mm_loadu_pd(tMatrix + 10);
  matrix34 = _mm_loadu_pd(tMatrix + 12);
  matrix56 = _mm_loadu_pd(tMatrix + 14);
  reg_storage = _mm_mul_pd(buffer0, matrix01);
  reg_storage2 = _mm_mul_pd(buffer0, matrix34);

  matrix34 = _mm_mul_pd(buffer1, matrix12);
  matrix56 = _mm_mul_pd(buffer1, matrix56);
  reg_storage = _mm_add_pd(reg_storage, matrix34);
  reg_storage2 = _mm_add_pd(reg_storage2, matrix56);
  reg_storage = _mm_hadd_pd(reg_storage, reg_storage2);

  matrix01 = _mm_loadu_pd(parentConditionals + 2);
  matrix01 = _mm_mul_pd(reg_storage, matrix01);
  _mm_storeu_pd(parentConditionals + 2, matrix01);

#else
  _Parameter t1 = childVector[0] - childVector[3],
             t2 = childVector[1] - childVector[3],
             t3 = childVector[2] - childVector[3], t4 = childVector[3];

  parentConditionals[0] *=
      tMatrix[0] * t1 + tMatrix[1] * t2 + tMatrix[2] * t3 + t4;
  parentConditionals[1] *=
      tMatrix[4] * t1 + tMatrix[5] * t2 + tMatrix[6] * t3 + t4;
  parentConditionals[2] *=
      tMatrix[8] * t1 + tMatrix[9] * t2 + tMatrix[10] * t3 + t4;
  parentConditionals[3] *=
      tMatrix[12] * t1 + tMatrix[13] * t2 + tMatrix[14] * t3 + t4;
#endif

}

#endif
