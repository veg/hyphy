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

#ifndef __CALCNODE__
#define __CALCNODE__

#include "category.h"
#include "hy_strings.h"
#include "parser.h"
#include "variablecontainer.h"

#define UNROOTED 0
#define ROOTED_LEFT 1
#define ROOTED_RIGHT 2

#define HY_REPLACE_BAD_BRANCH_LENGTH_WITH_THIS 0.000000001
#define fReusePreviouslyAllocatedMatrixExponentials 0x01

//_______________________________________________________________________________________________

class _CalcNode : public _VariableContainer {

public:
  // constructors

  /**
   * @brief Construct a new _CalcNode object
   *
   */
  _CalcNode(void); // default constructor, doesn't do much
  /**
   * @brief Construct a new _CalcNode object
   *
   * @param s1
   * @param s2
   * @param i
   * @param vc
   * @param avl
   */
  _CalcNode(_String, _String, int = 4, _VariableContainer * = nil,
            _AVLListXL * = nil);
  // construct a node from a string of the form
  // codeBase specifies the number of distinct states (4 for nucleotides, 61 for
  // codons etc) matrix name, <optional comma separated variable declarations,
  // inititalizations> also should be passed the pointer to a container tree

  /**
   * @brief Construct a new _CalcNode object
   *
   * @param source
   * @param parentTree
   */
  _CalcNode(_CalcNode *source, _VariableContainer *parentTree);

  /**
   * @brief Destroy the _CalcNode object
   *
   */
  virtual ~_CalcNode(void);
  /**
   * @brief Compute the object
   *
   * @return HBLObjectRef
   */
  virtual HBLObjectRef Compute(void) { return this; }

  /**
   * @brief Get the Object Class
   *
   * @return unsigned
   */
  virtual unsigned long ObjectClass(void) const { return TREE_NODE; }

  /**
   * @brief Duplicate the object
   *
   * @param brc
   */
  virtual void Duplicate(BaseRefConst);

  /**
   * @brief Free up memory
   *
   * @param l
   * @return long
   */
  virtual long FreeUpMemory(long);

  /**
   * @brief Initialize the object
   *
   * @param sc
   * @param i
   * @param vc
   * @param avl
   */
  void InitializeCN(_String const &, int, _VariableContainer *,
                    _AVLListXL * = nil);

  /**
   * @brief Make a dynamic copy of the object
   *
   * @return BaseRef
   */
  virtual BaseRef makeDynamic(void) const;
  // creates a dynamic copy of this object

  /**
   * @brief Convert the object to a string
   *
   * @param ul
   * @return BaseRef
   */
  virtual BaseRef toStr(unsigned long = 0UL);
  // converts this object to string

  /**
   * @brief Access the i-th element of the probabilities (i = 0..codeBase-1)
   *
   * @param ul
   * @return hyFloat&
   */
  hyFloat &operator[](unsigned long);
  // access the i-th element of the
  // probabilities (i = 0..codeBase-1)

  /**
   * @brief Set the Code Base
   *
   * @param l
   */
  void SetCodeBase(long);
  // change the codeBase value for this node
  // this will resize the vector used to handle frequencies

  /**
   * @brief Clear the object
   *
   */
  virtual void Clear(void);

  /**
   * @brief Recompute the matrix
   *
   * @param l1
   * @param l2
   * @param m
   * @param l3
   * @param sl
   * @param l4
   * @param direct_copy
   * @return true
   * @return false
   */
  bool RecomputeMatrix(long = 0, long = 1, _Matrix * = nil, _List * = nil,
                       _SimpleList * = nil, _List * = nil,
                       bool direct_copy = false);
  // reexponentiate the transition matrix and
  // store it in compExp.
  // return TRUE if the matrix is an explicit exponential form

  /**
   * @brief Check if the object has changed
   *
   * @param b
   * @param cache
   * @return true
   * @return false
   */
  virtual bool HasChanged(bool = false, _AVLListX *cache = nil);
  /**
   * @brief Check if a new category exponential is needed
   *
   * @param l
   * @param cache
   * @return true
   * @return false
   */
  virtual bool NeedNewCategoryExponential(long = -1L,
                                          _AVLListX *cache = nil) const;
  /**
   * @brief Set the Model
   *
   * @param l
   * @param avl
   */
  virtual void SetModel(long, _AVLListXL *);

  /**
   * @brief Get the probability at a given index.
   *
   * @param k The index.
   * @return The probability.
   */
  hyFloat GetProbs(long k) { return theProbs[k]; }
  /**
   * @brief Get the array of probabilities.
   *
   * @return A pointer to the array of probabilities.
   */
  hyFloat *GetProbs(void) { return theProbs; }

  /**
   * @brief Check if the exponentials should be cleared
   *
   * @return true
   * @return false
   */
  bool clear_exponentials(void) const {
    return (flags & fReusePreviouslyAllocatedMatrixExponentials) == 0;
  }
  /**
   * @brief Reuse the exponentials
   *
   */
  void reuse_exponentials(void) {
    flags = flags | fReusePreviouslyAllocatedMatrixExponentials;
  }

  /**
   * @brief Set the CompExp
   *
   * @param m
   * @param l
   * @param do_exponentiation
   */
  void SetCompExp(_Matrix *, long = -1, bool do_exponentiation = false);
  /**
   * @brief Set the CompMatrix
   *
   * @param l
   */
  void SetCompMatrix(long);
  /**
   * @brief Get the CompExp
   *
   * @param catID
   * @param b
   * @return _Matrix*
   */
  _Matrix *GetCompExp(long catID = -1, bool = false) const;

  /**
   * @brief Recurse the Monte Carlo simulation
   *
   * @param l
   * @param n
   * @param first
   * @param rooted
   * @return _Formula*
   */
  _Formula *RecurseMC(long, node<long> *, bool first = false,
                      char rooted = UNROOTED);

  /**
   * @brief Get the Code Base
   *
   * @return long
   */
  long GetCodeBase(void) { return cBase; }

  /**
   * @brief Compute the branch length
   *
   * @return hyFloat
   */
  hyFloat ComputeBranchLength(void);
  /**
   * @brief Set the Dependance
   *
   * @param l
   * @return long
   */
  virtual long SetDependance(long);

  /**
   * @brief Locate the node in the tree
   *
   * @return node<long>*
   */
  node<long> *LocateMeInTree(void) const;
  // return the tree structure node corresponing to this one...
  /**
   * @brief Convert to a simple matrix
   *
   * @param category_count
   */
  void ConvertToSimpleMatrix(unsigned long category_count);
  /**
   * @brief Convert from a simple matrix
   *
   * @param category_count
   */
  void ConvertFromSimpleMatrix(unsigned long category_count);
  /**
   * @brief Compute the model matrix
   *
   * @param expMe
   * @return _Matrix*
   */
  _Matrix *ComputeModelMatrix(bool expMe = false);
  /**
   * @brief Get the The Model ID
   *
   * @return long
   */
  long GetTheModelID(void) { return theModel; }
  /**
   * @brief Match a subtree
   *
   * @param cn
   * @return true
   * @return false
   */
  bool MatchSubtree(_CalcNode *);
  /**
   * @brief Remove the model
   *
   */
  virtual void RemoveModel(void);
  /**
   * @brief Replace the model
   *
   * @param modelName
   * @param parentTree
   * @param aCache
   * @param clean_locals
   */
  virtual void ReplaceModel(_String &modelName, _VariableContainer *parentTree,
                            _AVLListXL *aCache = nil, bool clean_locals = true);

  /**
   * @brief Clear the category map
   *
   */
  virtual void ClearCategoryMap(void) { remapMyCategories.Clear(); }

  /**
   * @brief Get the Branch Spec
   *
   * @return _String*
   */
  _String *GetBranchSpec(void);

  /**
   * @brief Get the ith category
   *
   * @param i
   * @return _CategoryVariable*
   */
  _CategoryVariable *get_ith_category(long i) const {
    return (_CategoryVariable *)LocateVar(categoryVariables.get(i));
  }

  /**
   * @brief Map a global category to a local category
   *
   * @param l
   * @return long
   */
  long map_global_to_local_category(long) const;

  /**
   * @brief Get the Parent Tree
   *
   * @return _VariableContainer*
   */
  _VariableContainer *ParentTree(void);

  /**
   * @brief This function will take a list of category variables (assumed to be a superset of categoryVariables) and the number of categories in each (second argument) and the multipliers for each category in the composite category and populate remapMyCategories
   *
   * @param l
   * @param sl1
   * @param sl2
   */
  virtual void SetupCategoryMap(_List &, _SimpleList &, _SimpleList &);
  /* 20090324: SLKP
   This function will take a list of category variables (assumed to be a
   superset of categoryVariables) and the number of categories in each (second
   argument) and the multipliers for each category in the composite category and
   populate remapMyCategories (see comments)
   */

  friend class _TheTree;

  /**
   * @brief Converts a string of form ":[\d\.]\+" into a double
   *
   * @param branch_length the branch length specifier
   * @return the parsed branch length or default "small" value (or -1 when an empty string is passed)
   */
  static hyFloat ProcessTreeBranchLength(_String const &branch_length);

public:
  hyFloat *theProbs; // list of transitional probabilities
                     // deprecate??
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
      flags;

  /*nodeIndex,
  referenceNode,
  slaveNodes;*/

  // deprecate??
};

/*----------------------------------------------------------------------------------------------------------*/

/**
 * @brief Map a node to a _CalcNode.
 *
 * @tparam data_type The type of data in the node.
 * @param n The node to map.
 * @return A pointer to the _CalcNode.
 */
template <class data_type> _CalcNode *map_node_to_calcnode(node<data_type> *n) {
  if (n) {
    return (_CalcNode *)(LocateVar(n->in_object));
  }
  return nil;
}

/*----------------------------------------------------------------------------------------------------------*/

#endif
