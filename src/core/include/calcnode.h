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

#ifndef     __CALCNODE__
#define     __CALCNODE__

#include "variablecontainer.h"
#include "hy_strings.h"
#include "parser.h"
#include "category.h"

#define UNROOTED                        0
#define ROOTED_LEFT                     1
#define ROOTED_RIGHT                    2

#define HY_REPLACE_BAD_BRANCH_LENGTH_WITH_THIS  0.000000001
#define fReusePreviouslyAllocatedMatrixExponentials 0x01

//_______________________________________________________________________________________________

class       _CalcNode: public _VariableContainer {
    
public:
    
    // constructors
    
    _CalcNode           (void);                 // default constructor, doesn't do much
    _CalcNode           (_String, _String, int  = 4, _VariableContainer* = nil, _AVLListXL * = nil);
    // construct a node from a string of the form
    // codeBase specifies the number of distinct states (4 for nucleotides, 61 for codons etc)
    // matrix name, <optional comma separated variable declarations, inititalizations>
    // also should be passed the pointer to a container tree
    
    _CalcNode           (_CalcNode* source, _VariableContainer* parentTree);
    
    virtual             ~_CalcNode      (void);
    virtual             HBLObjectRef      Compute                             (void) {
        return this;
    }
    
    virtual unsigned long        ObjectClass     (void) const {
        return TREE_NODE;
    }
    
    virtual void        Duplicate       (BaseRefConst);
    
    virtual long        FreeUpMemory    (long);
    
    void                InitializeCN    ( _String const&, int, _VariableContainer*, _AVLListXL * = nil);
    
    virtual BaseRef     makeDynamic     (void) const;
    // creates a dynamic copy of this object
    
    virtual BaseRef     toStr           (unsigned long = 0UL);
    // converts this object to string
    
    hyFloat&         operator[]      (unsigned long);
    // access the i-th element of the
    // probabilities (i = 0..codeBase-1)
    
    void                SetCodeBase      (int);
    // change the codeBase value for this node
    // this will resize the vector used to handle frequencies
    
    virtual     void        Clear                       (void);

    bool                RecomputeMatrix  (long = 0, long = 1,_Matrix* = nil, _List* = nil, _SimpleList* = nil, _List* = nil, bool direct_copy = false);
    // reexponentiate the transition matrix and
    // store it in compExp.
    // return TRUE if the matrix is an explicit exponential form
    
    virtual bool        HasChanged       (bool = false, _AVLListX * cache = nil);
    virtual bool        NeedNewCategoryExponential (long = -1L, _AVLListX * cache = nil) const;
    virtual void        SetModel         (long, _AVLListXL*);
    
    hyFloat          GetProbs        (long k) {
        return theProbs[k];
    }
    hyFloat*         GetProbs        (void) {
        return theProbs;
    }
    
    bool                clear_exponentials (void) const { return (flags & fReusePreviouslyAllocatedMatrixExponentials) == 0;}
    void                reuse_exponentials (void) { flags = flags | fReusePreviouslyAllocatedMatrixExponentials;}

    void                SetCompExp      (_Matrix*, long = -1, bool do_exponentiation = false);
    void                SetCompMatrix   (long);
    _Matrix*            GetCompExp      (long catID = -1, bool = false) const;
    
    _Formula*           RecurseMC       (long , node<long>* , bool first = false, char rooted = UNROOTED);
    
    long                GetCodeBase     (void) {
        return cBase;
    }
    
    hyFloat          ComputeBranchLength    (void);
    virtual long        SetDependance   (long);
    
    node<long>*         LocateMeInTree  (void) const;
    // return the tree structure node corresponing to this one...
    void                ConvertToSimpleMatrix (unsigned long category_count) ;
    void                ConvertFromSimpleMatrix (unsigned long category_count);
    _Matrix*            ComputeModelMatrix(bool expMe=false);
    long                GetTheModelID   (void) {
        return theModel;
    }
    bool                MatchSubtree    (_CalcNode*);
    virtual void        RemoveModel     (void);
    virtual void        ReplaceModel    (_String & modelName, _VariableContainer* parentTree, _AVLListXL* aCache = nil, bool clean_locals = true);
    
    virtual void        ClearCategoryMap(void) {
        remapMyCategories.Clear();
    }
    
    _String*            GetBranchSpec (void);
  
    _CategoryVariable*  get_ith_category (long i) const {
      return (_CategoryVariable*)LocateVar(categoryVariables.get (i));
    }
    
    long                map_global_to_local_category (long) const;
    
    _VariableContainer*           ParentTree      (void);
    
    virtual void        SetupCategoryMap(_List&, _SimpleList&, _SimpleList&);
    /* 20090324: SLKP
     This function will take a list of category variables (assumed to be a superset of categoryVariables)
     and the number of categories in each (second argument)
     and the multipliers for each category in the composite category
     and populate remapMyCategories
     (see comments)
     */
    
    friend  class       _TheTree;
    
    /**
     * Converts a string of form ":[\d\.]\+" into a double
     * \n SLKP 20100831: a utility function to handle the
     * conversion of branch length strings to parameters
     * \n\n \b Example: string = ":3.14" returns 3.14
     * \n If it is not of correct form, it will return 1e-10
     * @param branch_length the branch length specifier
     * @return the parsed branch length or default "small" value (or -1 when an empty string is passed)
     * Revision history
     - SLKP 20170616 moved from _String v2.3, changed from being a member of _String
     */
    static hyFloat    ProcessTreeBranchLength (_String const& branch_length);
    
    
public:
    hyFloat*     theProbs;       // list of transitional probabilities
    // deprecate??
protected:
    
    _SimpleList     categoryVariables,
    categoryIndexVars,
    remapMyCategories;
    
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
    
    
    _Matrix   *     compExp;        // matrix exponential computed previously
    _Matrix   **    matrixCache;    // only meaningful for category computations
    
    long            cBase,          // dimension of theProbs
                    flags;

    /*nodeIndex,
    referenceNode,
    slaveNodes;*/
    
    // deprecate??
};

/*----------------------------------------------------------------------------------------------------------*/

template <class data_type> _CalcNode* map_node_to_calcnode (node<data_type>* n) {
  if (n) {
    return (_CalcNode*)(LocateVar(n->in_object));
  }
  return nil;
}

/*----------------------------------------------------------------------------------------------------------*/

#endif
