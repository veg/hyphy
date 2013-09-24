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

#ifndef __CALCNODE__
#define __CALCNODE__

#include "legacy_parser.h"
#include "classes.h"
#include "site.h"

#define UNROOTED 0
#define ROOTED_LEFT 1
#define ROOTED_RIGHT 2

#define _HY2CALCNODE(X) (dynamic_cast<_CalcNode*>(X))

//_______________________________________________________________________________________________

class _CalcNode : public _VariableContainer {

public:

  // constructors

  _CalcNode(void); // default constructor, doesn't do much
  _CalcNode(_String, _String, int = 4, _VariableContainer * = nil,
            _AVLListXL * = nil);
  // construct a node from a string of the form
  // codeBase specifies the number of distinct states (4 for nucleotides, 61 for
  // codons etc)
  // matrix name, <optional comma separated variable declarations,
  // inititalizations>
  // also should be passed the pointer to a container tree

  _CalcNode(_CalcNode *source, _VariableContainer *parentTree);
  _CalcNode(_CalcNode&);

  virtual ~_CalcNode(void);
  virtual _PMathObj Compute(void);

  virtual unsigned long ObjectClass(void) { return TREE_NODE; }

  virtual void Duplicate(BaseRef);

  virtual long FreeUpMemory(long);

  void InitializeCN(_String &, int, _VariableContainer *, _AVLListXL * = nil);

  virtual BaseRef makeDynamic(void);
  // creates a dynamic copy of this object

  virtual BaseRef toStr(void);
  // converts this object to string

  _Parameter &operator[](unsigned long);
  // access the i-th element of the
  // probabilities (i = 0..codeBase-1)

  void SetCodeBase(int);
  // change the codeBase value for this node
  // this will resize the vector used to handle frequencies

  bool RecomputeMatrix(long = 0, long = 1, _Matrix * = nil, _List * = nil,
                       _SimpleList * = nil, _List * = nil);
  // reexponentiate the transition matrix and
  // store it in compExp.
  // return TRUE if the matrix is an explicit exponential form

  virtual bool HasChanged(void);
  virtual bool NeedToExponentiate(long = -1);
  virtual void SetModel(long, _AVLListXL *);

  bool IsFlagged(void) { return theProbs[0] == -3.1415296; }
  void SetFlag(void) { theProbs[0] = -3.1415296; }

  void SetSummedFlag(void) {
    if (theProbs[0] >= 0) {
      theProbs[0] -= 2.0;
    }
  }
  bool IsSummedFlagged(void) { return theProbs[0] < 0.0; }
  void RemoveSummedFlag(void) {
    if (theProbs[0] < 0) {
      theProbs[0] += 2.0;
    }
  }

  _Parameter GetProbs(long k) { return theProbs[k]; }
  _Parameter *GetProbs(void) { return theProbs; }

  void SetCompExp(_Matrix *, long = -1);
  void SetCompMatrix(long);
  _Matrix *GetCompExp(long catID = -1);

  _Formula *RecurseMC(long, node<long> *, bool first = false,
                      char rooted = UNROOTED);

  long GetCodeBase(void) { return cBase; }

  _Parameter BranchLength(void);
  virtual long SetDependance(long);

  node<long> *LocateMeInTree(void);
  // return the tree structure node corresponing to this one...
  long ConvertToSimpleMatrix(void);
  void ConvertFromSimpleMatrix(void);
  _Matrix *ComputeModelMatrix(bool expMe = false);
  long GetTheModelID(void) { return theModel; }
  bool MatchSubtree(_CalcNode *);
  virtual void RemoveModel(void);
  virtual void ReplaceModel(_String &modelName, _VariableContainer *parentTree);

  virtual long CheckForReferenceNode(void);

  void SetRefNode(long rn) {
    referenceNode = rn;
    slaveNodes = 0;
  }
  void AddRefNode(void) { referenceNode--; }

  virtual void ClearCategoryMap(void) { remapMyCategories.Clear(); }

  _VariableContainer *ParentTree(void);

  virtual void SetupCategoryMap(_List &, _SimpleList &, _SimpleList &);
  /* 20090324: SLKP
          This function will take a list of category variables (assumed to be a
          superset of categoryVariables)
          and the number of categories in each (second argument)
          and the multipliers for each category in the composite category
          and populate remapMyCategories
          (see comments)
   */

  friend class _TheTree;

public:
  _Parameter *theProbs; // list of transitional probabilities
  long lastState;

protected:

  _SimpleList categoryVariables, categoryIndexVars, remapMyCategories;

  /*
        20090324: SLKP
            because this calcnode may be a part of a likelihood
            function partition with more category variables,
            this mapper object takes the composite category ID
            from the container (likelihood function) and maps it
            to a category ID understood by this node; a composite
            ID is mapped to an N+1 tuplet (N is the number of
            category variables that this node depends on):

            composite ID inside this _CalcNode
            classes of each category variable in the same order
            as they appear in categoryVariables

            Hence the list will have (N+1)*(Container classes)
            entries.

            This list is populated using the SetupCategoryMap
            and cleared using ClearCategoryMap
     */

  _Matrix *compExp;      // matrix exponential computed previously
  _Matrix **matrixCache; // only meaningful for category computations

  long cBase, // dimension of theProbs
      nodeIndex, referenceNode, slaveNodes;
};

//_______________________________________________________________________________________________

#define HY_BRANCH_SELECT 0x01
#define HY_BRANCH_DESELECT 0xFFFFFFFE

struct nodeCoord {

  _Parameter h, v, auxD, bL, label1, label2;

  long varRef, auxL, textWidth, color, labelColor, flags;

  _String branchName, branchTag;

}; // used for tree imaging

//_______________________________________________________________________________________________

typedef bool _HYTopologyTraversalFunction(node<long> *, Ptr);

//_______________________________________________________________________________________________

class _TheTree; // forward declaration for xlc

//_______________________________________________________________________________________________
#if USE_SCALING_TO_FIX_UNDERFLOW
extern _Parameter scalingLogConstant;
#endif

//_______________________________________________________________________________________________

extern char isDefiningATree;
extern _String expectedNumberOfSubs, stringSuppliedLengths, includeModelSpecs,
    treeOutputAVL, treeOutputLayout;
#ifdef _OPENMP
#include "omp.h"
#endif
#endif
