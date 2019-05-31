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

#include "parser.h"
#include "classes.h"
#include "dataset_filter.h"
#include "dataset_filter_numeric.h"
#include "associative_list.h"
#include "vector.h"
#include "site.h"


#define UNROOTED                        0
#define ROOTED_LEFT                     1
#define ROOTED_RIGHT                    2

#define HY_REPLACE_BAD_BRANCH_LENGTH_WITH_THIS  0.000000001

#ifdef MDSOCL

class _OCLEvaluator 
{

private:
	// OpenCL Vars

	// Forward Declarations
	// *********************************************************************
	void Cleanup (int iExitCode);
	unsigned int roundUpToNextPowerOfTwo(unsigned int x);
	double roundDoubleUpToNextPowerOfTwo(double x);
	// So the only thing that needs to be passed as an update for each LF is flatTree and flatCLeaves
	// as those are what goes into the new transition matrix stuff. 
	// So I could have a launchmdsocl that takes everything and if stuff is not NULL than update, if it is
	// null use the existing values. How about that?
	// The problem is that I need to have essentially the first LF's information to properly set everything up. 
	// And how do I have subsequent LF's not pass stuff. 
	// Oi, how do I not have to convert this into an object...
	// alright, I can probably keep all of this in likefunc. That is because the LFEvaluation is done in the
	// calc node
	double oclmain(void);
	bool contextSet;
	int setupContext(void);


public:
	void init(		long esiteCount,
						long ealphabetDimension,
						hyFloat* eiNodeCache);


	double launchmdsocl(	_SimpleList& updateNodes,
							_SimpleList& flatParents,
							_SimpleList& flatNodes,
							_SimpleList& flatCLeaves,
							_SimpleList& flatLeaves,
							_SimpleList& flatTree,
							hyFloat* theProbs,
							_SimpleList& theFrequencies,
							long* lNodeFlags,
							_SimpleList& taggedInternals,
							_Vector* lNodeResolutions);
	~_OCLEvaluator()
	{
		Cleanup(EXIT_SUCCESS);
	}


};
#endif




//_______________________________________________________________________________________________

class       _CalcNode: public _VariableContainer
{

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
    virtual             _PMathObj      Compute                             (void);

    virtual unsigned long        ObjectClass     (void) {
        return TREE_NODE;
    }

    virtual void        Duplicate       (BaseRefConst);

    virtual long        FreeUpMemory    (long);

    void                InitializeCN    ( _String&, int, _VariableContainer*, _AVLListXL * = nil);

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

    bool                RecomputeMatrix  (long = 0, long = 1,_Matrix* = nil, _List* = nil, _SimpleList* = nil, _List* = nil);
    // reexponentiate the transition matrix and
    // store it in compExp.
    // return TRUE if the matrix is an explicit exponential form

    virtual bool        HasChanged       (bool = false);
    virtual bool        NeedNewCategoryExponential (long = -1L) const;
    virtual void        SetModel         (long, _AVLListXL*);

    bool                IsFlagged        (void) {
        return theProbs[0]==-3.1415296;
    }
    void                SetFlag          (void) {
        theProbs[0]=-3.1415296;
    }

    void                SetSummedFlag    (void) {
        if (theProbs[0]>=0) {
            theProbs[0] -= 2.0;
        }
    }
    bool                IsSummedFlagged (void) {
        return theProbs[0]<0.0;
    }
    void                RemoveSummedFlag (void) {
        if (theProbs[0]<0) {
            theProbs[0]+= 2.0;
        }
    }

    hyFloat          GetProbs        (long k) {
        return theProbs[k];
    }
    hyFloat*         GetProbs        (void) {
        return theProbs;
    }

    void                SetCompExp      (_Matrix*, long = -1);
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
    void                ConvertToSimpleMatrix (void) const;
    void                ConvertFromSimpleMatrix (void);
    _Matrix*            ComputeModelMatrix(bool expMe=false);
    long                GetTheModelID   (void) {
        return theModel;
    }
    bool                MatchSubtree    (_CalcNode*);
    virtual void        RemoveModel     (void);
    virtual void        ReplaceModel    (_String & modelName, _VariableContainer* parentTree);
    
    virtual long        CheckForReferenceNode
    (void);

    void        SetRefNode      (long rn) {
        referenceNode = rn;
        slaveNodes = 0;
    }
    void        AddRefNode      (void) {
        referenceNode --;
    }

    virtual void        ClearCategoryMap(void) {
        remapMyCategories.Clear();
    }
  
    _String*            GetBranchSpec (void);
  
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
    long            lastState;

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
                    nodeIndex,
                    referenceNode,
                    slaveNodes;
};

//_______________________________________________________________________________________________

#define     HY_BRANCH_SELECT        0x01
#define     HY_BRANCH_DESELECT      0xFFFFFFFE

class      nodeCoord {

public:
    hyFloat  h,
                v,
                auxD,
                bL,
                label1,
                label2;

    long        varRef,
                auxL,
                textWidth,
                color,
                labelColor,
                flags;
  
    operator long (void) {return varRef;}

    _String     branchName,
                branchTag;

}; // used for tree imaging

//_______________________________________________________________________________________________

typedef bool _HYTopologyTraversalFunction (node<long>*, hyPointer);

//_______________________________________________________________________________________________

class _TheTree; // forward declaration for xlc


#define kGetNodeStringForTreeName   0x01
#define kGetNodeStringForTreeModel  0x02

//_______________________________________________________________________________________________

class _TreeTopology: public _CalcNode {

protected:

    virtual void            PreTreeConstructor                  (bool);
    virtual bool            MainTreeConstructor                 (_String const&,bool = true, _AssociativeList* mapping = nil);
    virtual void            PostTreeConstructor                 (bool);
    node<long>*     prepTree4Comparison                 (_List&, _SimpleList&, node<long>* = nil) const;
    void            destroyCompTree                     (node<long>*) const;
    _List*          SplitTreeIntoClustersInt            (node<long>*, _List*, _AVLListX&, long, long) const;
    char            internalTreeCompare                 (node<long>*, node<long>*, _SimpleList*, char, long, node<long>*, _TreeTopology const*, bool = false) const;
    char            internalNodeCompare                 (node<long>*, node<long>*, _SimpleList&, _SimpleList*, bool, long, node<long>*, _TreeTopology const*, bool = false) const;
    virtual _PMathObj       FlatRepresentation                  (void);
    void            FindCOTHelper                       (node<long>*, long, _Matrix&, _Matrix&, _Matrix&, _List&, _AVLListX&, hyFloat);
    void            FindCOTHelper2                      (node<long>*, _Matrix&, _Matrix&, _AVLListX&, node<long>*, hyFloat);
    void            AddANode                            (_PMathObj);
    /*

     20091006: SLKP

     given an AVL with at least three key -> string pairs
        "NAME" -> string (not a currently used node name)
        "WHERE" -> string (an existing node)
        "PARENT" -> string (not a currently used node name)



     */

    void            RemoveANode                          (_PMathObj);
    
    /*
     
        Delete a node from the tree by name
     
     */
  

public:

    node<long>      *theRoot;
  
    _List           flatTree,
                    flatCLeaves;

    char            rooted;

    virtual void            toFileStr                           (FILE*, unsigned long);
    virtual BaseRef         toStr                               (unsigned long = 0UL);
    void            RerootTreeInternalTraverser         (node<long>* iterator, long, bool,_StringBuffer&, long  = -1, bool = false) const;

    _TreeTopology                       (void);
    _TreeTopology                       (_String const, _String const&, bool = true, _AssociativeList* mapping = nil);
    _TreeTopology                       (_String*);
    _TreeTopology                       (_TheTree*);

    virtual                 ~_TreeTopology                      (void);

    virtual  _FString*      Compare                             (_PMathObj);
    virtual  BaseRef        makeDynamic                         (void) const;
    node<long>* CopyTreeStructure                   (node<long>*, bool) const;
    virtual  bool           FinalizeNode                        (node<long>*, long, _String, _String const&, _String&, _String* = NULL);


    virtual _PMathObj       ExecuteSingleOp                     (long, _List* = nil, _hyExecutionContext* context = _hyDefaultExecutionContext);
    virtual void            EdgeCount                           (long&, long&);
    // SLKP 20100827: a utility function to count edges in a tree
    //              : note that the root node WILL be counted as an internal node
    //              : writes [leaf count, internal node count] into the arguments

    virtual _PMathObj       TipCount                            (void);
    virtual _PMathObj       BranchCount                         (void);
    virtual _PMathObj       AVLRepresentation                   (_PMathObj);
    virtual unsigned long   ObjectClass                         (void) {
        return TOPOLOGY;
    }
    virtual bool IsDegenerate(void) { return theRoot && theRoot->get_num_nodes() == 1L; }
  
 
    virtual _AssociativeList*
    FindCOT                             (_PMathObj);

    node<long>      *FindNodeByName                     (_String const*) const;
    /*

     20091006: SLKP

     return the node with the name supplied by the argument
     or nil if no such node exists

     */

    /*
    void            DepthWiseT                          (bool = false, _HYTopologyTraversalFunction* = nil, hyPointer = nil);
    void            DepthWiseTRight                     (bool = false);
    void            DepthWiseTLevel                     (long& level, bool = false);
    void            StepWiseT                           (bool = false, _HYTopologyTraversalFunction* = nil, hyPointer = nil);
    void            StepWiseTLevel                      (long&, bool = false);
    void            LeafWiseT                           (bool = false);
    */
    
    _List*          MapNodesToModels                    (void);

    virtual _String const  GetNodeName                         (node<long> *, bool = false) const;
    virtual const _String*        GetNodeModel                        (node<long> *) const;
    virtual void            GetBranchLength                     (node<long> *, _String&, bool = false) const;
    // SLKP 20100901:
    //               added a boolean flag to ask to return branch length expression (if true) (returns "" for topologies)
    //               just the numeric value (if false)


    virtual hyFloat         GetBranchLength                     (node<long> *) const;
    virtual void            GetBranchValue                      (node<long> *, _String&) const;
    virtual void            GetBranchVarValue                   (node<long> *, _String&, long) const;
    virtual _String const  GetNodeStringForTree                (node<long> *, int flags) const;
    virtual void            PasteBranchLength                   (node<long> *, _StringBuffer &, long, hyFloat factor = 1.) const;

    node<long>&     GetRoot                             (void) const {
      return  *theRoot;
    }
    void            SetRoot                             (node<long>* r) {
        theRoot = r;
    }
  
    /*node<long>&     GetCurrentNode                      (void) {
        return *currentNode;
    }*/
  
    const _List     RetrieveNodeNames                   (bool doTips, bool doInternals, int travseralType) const;
    void            SubTreeString                       (node<long>* root, _StringBuffer &, bool all_names = false, hyB = -1, _AVLListXL* = nil) const;

    _String         CompareTrees                        (_TreeTopology*) const;
    const _String         MatchTreePattern                    (_TreeTopology const*) const;
    virtual _PMathObj       TipName                             (_PMathObj);
    _PMathObj       TreeBranchName                          (_PMathObj, bool = false, _PMathObj = nil);
    virtual _PMathObj       BranchLength                        (_PMathObj);
    virtual _PMathObj       RerootTree                          (_PMathObj);
    _List*          SplitTreeIntoClusters               (unsigned long, unsigned long) const;
    void            SetLeafName                         (long, _String*);
    _String         DetermineBranchLengthMappingMode    (_String*, char&);
    _AssociativeList*
    SplitsIdentity                      (_PMathObj);
    /* 20090609: SLKP
        given a tree agrument (p), the function returns an AVL with a 2x1 matrix (key "CLUSTERS")
        and a string (key "CONSENSUS");
        The first cell contains the number of splits in *this
        The second cell contains the number of splits in the argument that are present in *this

        This entry will contain -1 if the argument is invalid (nil or not a tree)
        and if the set of leaves differs between two trees

        The string will be empty (incomparable trees or another exception) or the Newick String
        with the consensus string

    */
    bool            ConvertToPSW                        (_AVLListX&,_List*, _SimpleList&,bool = false);
    /* 20090612: SLKP
       20100510: Modified the function to also return internal node names in the second AVL

        covert the topology into the post-order with weights representation
        The first argument maps node names to their internal indices in the traversal order
        (note that leaves are numbered 1..leaves-1 and internal indices as leaves-1...leaves+inodes-1)
        The second argument stores doubles for each node; the index in the array itself (mod 2),
        corresponds to the appropriate step in the post-order traversal;

        <node index, number of nodes in the subtree below>

        (((A,B)N1,C)N2,(D,E)N3)N0 will result in

        A->0; B->1; N1->5; C->2; N2->6; D->3; E->4; N3->7; N0->8
        <0,0>,<1,0>,<5,2>,<2,0>,<6,4>,<3,0>,<4,0>,<7,2>,<8,8>, <5,4>

        The last two entries store the number of leaves and internal nodes


        if the bool argument is TRUE, then each LEAF in the tree must be found in the reference
        dictionary supplied by the FIRST argument; false will be returned if this is not the case.
    */

    _String*        ConvertFromPSW                      (_AVLListX&,_SimpleList&);
    /* 20090612: SLKP
            given a PSW tree traversal order and a labeling legend,
            return the Newick string for the tree
    */

    void            ComputeClusterTable                 (_SimpleList&, _SimpleList&);
    /* given the PSW traversal representation (arg 2)
       compute the cluster table (as defined in William HE Day "Optimal Algorithms for Comparing Trees With Labels",
       Page 16) and store in arg 1
            the list a is a flat representation for an Nx3 (N = number of leaves) table
            clusters spanning L<->R leaves will be stored in either row L or row R
            i-th entry
            <L = leftmost cluster leaf (in traversal order),
             R = rightmost cluster leaf (in traversal order),
             F = a binary toggle (set to 0 by this procedure)
     */





};


#if USE_SCALING_TO_FIX_UNDERFLOW
extern hyFloat scalingLogConstant;
#endif

//_______________________________________________________________________________________________

class _TheTree: public _TreeTopology
{

// theModel matrix of _TheTree contains the column matrix of probabilities, which is computed
// based on the DataSetFilter passed on to the tree at the initialization stage

public:

    _TheTree ();                                                // default constructor - doesn't do much
    _TheTree (_String name, _String const& parms, bool = true);       // builds a tree from a string
    _TheTree (_String name, _TreeTopology*);                    // builds a tree from a tree topology
    _TheTree (_String name, _TheTree*);                    // builds a tree from another tree


    virtual                 ~_TheTree                   (void);
    virtual bool            HasChanged                  (bool = false);
    virtual void            MarkDone                    (void);
    bool            HasChanged2                 (void);

    /*
    _CalcNode*      DepthWiseTraversal          (bool = false);
    //performs a post-order traversal
    _CalcNode*      DepthWiseTraversalRight     (bool = false);
    //performs a post-order tree traversal going right first
    _CalcNode*      DepthWiseTraversalLevel     (long&, bool = false);
    //performs a post-order tree traversal
    //storing current node depth
    _CalcNode*      StepWiseTraversal           (bool = false);
    //performs a pre-order  tree traversal
    _CalcNode*      StepWiseTraversalLevel      (long&, bool = false);
    //performs a pre-order wise tree traversal
    //storing current node depth

    _CalcNode*      LeafWiseTraversal           (bool = false);
    //iterate through the leaves (left-to-right)
     */

    virtual  bool           FinalizeNode                (node<long>*, long, _String, _String const&, _String&, _String* = NULL);
    virtual  BaseRef        makeDynamic                 (void) const;

    virtual  BaseRef        makeDynamicCopy             (_String*);
    node<long>* DuplicateTreeStructure      (node<long>*, _String*, bool);
    virtual  BaseRef        toStr                       (unsigned long = 0UL);
    virtual unsigned long           ObjectClass                 (void) {
        return TREE;
    }

    virtual _PMathObj       ExecuteSingleOp                     (long, _List* = nil, _hyExecutionContext* context = _hyDefaultExecutionContext);
    virtual  _PMathObj      TEXTreeString               (_PMathObj) const;
    virtual  _PMathObj      PlainTreeString             (_PMathObj,_PMathObj);

    virtual _String const   GetNodeName                 (node<long> *, bool = false) const;
    virtual  hyFloat        GetBranchLength             (node<long> *) const ;
    virtual  void           GetBranchValue              (node<long> *, _String&) const ;
    virtual  void           GetBranchVarValue           (node<long> *, _String&, long) const ;
    virtual _String const*        GetNodeModel                (node<long> *) const;
    
    void            InitializeTreeFrequencies   (_Matrix *, bool = false);


    hyFloat      Process3TaxonNumericFilter  (_DataSetFilterNumeric*, long = 0);



    _List*      RecoverAncestralSequences       (_DataSetFilter const*, _SimpleList const&, _List const&, hyFloat *, hyFloat const*, long, long*, _Vector*, bool = false);
    void        RecoverNodeSupportStates        (_DataSetFilter const*, long, _Matrix&);
    void        RecoverNodeSupportStates2       (node<long>*,hyFloat*,hyFloat*,long);
    _List*      SampleAncestors                 (_DataSetFilter*, node<long>*);
    void        PurgeTree                       (void);

    long        ComputeReleafingCost            (_DataSetFilter const*, long, long, _SimpleList* = nil, long = 0) const;
    long        ComputeReleafingCostChar        (_DataSetFilter const*, long, long) const;
    void        DumpingOrder                    (_DataSetFilter*, _SimpleList&);
    void        SetTreeCodeBase                 (long);
    long        IsLinkedToALF                   (long&) const;

    bool        HasCache                        (void) {
        return topLevelNodes.lLength>0;
    }

    long        GetLeafCount                    (void) {
        return flatLeaves.lLength;
    }

    long        GetINodeCount                   (void) {
        return flatNodes.lLength    ;
    }

    void        ScanAndAttachVariables          (void) const;
    void        ScanContainerForVariables       (_AVLList& l, _AVLList& l2, _AVLListX* tagger = nil, long weight = 0) const;
    void        ScanForDVariables               (_AVLList& l, _AVLList& l2) const;
    void        ScanForGVariables               (_AVLList&, _AVLList&, _AVLListX* tagger = nil, long weight = 0) const;
    void        ScanForCVariables               (_AVLList&) const;
    void        MolecularClock                  (_String const&, _List&) const;

    void        SetUp                           (void);
    void        SetUpMatrices                   (long);
    void        CleanUpMatrices                 (void);
    //void        BuildTopLevelCache              (void);
    void        KillTopLevelCache               (void);
    void        SetCompMatrices                 (long) const;

    virtual void        ClearConstraints                (void);

    bool        FindScalingVariables            (_SimpleList&) const;
    bool        HaveStringBranchLengths         (void) const;
    void        AssignLabelsToBranches          (node<nodeCoord>*, _String*, bool);

    node<nodeCoord>*
    AlignedTipsMapping                          (node<long>*, bool first = false, bool respectRoot = true) const;

    void        AlignNodes                      (node<nodeCoord>*) const;

    node<nodeCoord>*
    ScaledBranchMapping             (node<nodeCoord>* , _String*, long, long&, char) const;

    node<nodeCoord>*
    RadialBranchMapping             (node<long>* , node<nodeCoord>*, _String*, hyFloat, long&, hyFloat&, char);

    void        ScaledBranchReMapping           (node<nodeCoord>*, hyFloat) const;
    char&       RootedFlag                      (void) {
        return rooted;
    }

    nodeCoord   TreeTEXRecurse                  (node<nodeCoord>*,_StringBuffer&,hyFloat,hyFloat,long,long) const;
    void        TreePSRecurse                   (node<nodeCoord>*,_StringBuffer&,hyFloat,hyFloat,long,long,long,long,_AssociativeList* = nil, char = 0, hyFloat* = nil) const;

    bool        AllBranchesHaveModels           (long) const;
    void        ScanSubtreeVars                 (_List&, char, _CalcNode*) const;
    void        BuildINodeDependancies          (void);
    void        AllocateResultsCache            (long);
    long        CountTreeCategories             (void);
    void        CompileListOfModels             (_SimpleList&);

    void        MarkMatches                     (_DataSetFilter*,long,long) const;
    long        GetLowerBoundOnCost             (_DataSetFilter*, _SimpleList* = nil) const;
    _SimpleList&GetLeftINodes                   (void) {
        return leftiNodes;
    }
    bool        MatchLeavesToDF                 (_SimpleList&, _DataSetFilter*, bool) const;
    virtual void
    RemoveModel                     (void);
    _String*    TreeUserParams                  (void) const;


    const _String     CompareSubTrees                 (_TheTree*, node<long>*);
    const _String     FindMaxCommonSubTree            (_TheTree const*, long&, _List*) const;
  
  
    void        AddNodeNamesToDS                (_DataSet*, bool, bool, char) const;
    // if the
    hyFloat  PSStringWidth                   (_String const&);

    hyFloat  DetermineBranchLengthGivenScalingParameter (long, _String&, char) const;

    _AVLListX*  ConstructNodeToIndexMap         (bool) const;
    // 20090206: SLKP
    // makes an AVL of with keys storing memory addresses of node<long> tree nodes
    // and values showing the order in either flatLeaves (bool = false) or flatNodes (bool = true)

    void        MapPostOrderToInOderTraversal   (_SimpleList&, bool = true) const;
    // 20090306: SLKP
    // 20100511: SLKP
    // construct a post-order -> in-order traveral map for internal nodes
    // bool = true (internal nodes), bool = false (leaf nodes)

    void        AddBranchToForcedRecomputeList  (long idx)      {
        forceRecalculationOnTheseBranches << idx;
    }
    void        ClearForcedRecomputeList        (void)          {
        forceRecalculationOnTheseBranches.Clear();
    }
    bool        HasForcedRecomputeList          (void)          {
        return forceRecalculationOnTheseBranches.lLength;
    }
    // 20090306: SLKP
    // the previous two functions are used to manipulate the list of
    // branches that are marked as 'dirty' for LF computation purposes
    // because of external manipulations of the cache (e.g. computing the LF with one of the interior
    // nodes set to

    void        SetupCategoryMapsForNodes       (_List& , _SimpleList&, _SimpleList& );
    /* 20090325: SLKP
       a wrapper function to set up category variable maps for nodes;
       see comments for remapMyCategories in _CalcNode
     */

    const _CalcNode * GetNodeFromFlatIndex (long index) const;
  
    hyFloat  VerySimpleLikelihoodEvaluator   (_SimpleList&            updateNodes,
            _DataSetFilter*      theFilter,
            hyFloat*          iNodeCache,
            long       *             lNodeFlags,
            _Vector*      lNodeResolutions);

#ifdef MDSOCL
			hyFloat OCLLikelihoodEvaluator (			_SimpleList&	     updateNodes, 
                                                        _DataSetFilter*		 theFilter,
                                                        hyFloat*			 iNodeCache,
                                                         long	   *		 lNodeFlags,
                                                        _Vector*		 lNodeResolutions,
														_OCLEvaluator& OCLEval);
#endif

#ifdef  _SLKP_LFENGINE_REWRITE_
    void            SampleAncestorsBySequence       (_DataSetFilter const*, _SimpleList const&, node<long>*, _AVLListX const*, hyFloat const*, _List&, _SimpleList*, _List&, hyFloat const*, long);

    hyFloat      ComputeTreeBlockByBranch        (_SimpleList&, _SimpleList&, _SimpleList*, _DataSetFilter const*, hyFloat*, long*, hyFloat*, _Vector*, long&, long, long, long = -1, hyFloat* = nil, long* = nil, long = -1, long * = nil);
    long            DetermineNodesForUpdate         (_SimpleList&,  _List* = nil, long = -1, long = -1, bool = true);
    void            ExponentiateMatrices            (_List&, long, long = -1);
    void            FillInConditionals              (_DataSetFilter const*, hyFloat*,  _SimpleList*);

    void            ComputeBranchCache              ( _SimpleList&,
            long nodeID,
            hyFloat*         cache,
            hyFloat*         iNodeCache,
            _DataSetFilter const*     theFilter,
            long           *        lNodeFlags,
            hyFloat*         scalingAdjustments,
            long*                   siteCorrectionCounts,
            _Vector*     lNodeResolutions,
            long&                   overallScaler,
            long                    siteFrom,
            long                    siteTo,
            long                    catID,
            _SimpleList*            = nil,
            hyFloat*         = nil
                                                    );

    hyFloat          ComputeLLWithBranchCache         (
        _SimpleList&            siteOrdering,
        long                    brID,
        hyFloat*         cache,
        _DataSetFilter const*     theFilter,
        long                    siteFrom,
        long                    siteTo,
        long                    catID,
        hyFloat*         storageVec = nil
    );

    hyFloat          ComputeTwoSequenceLikelihood    (
        _SimpleList&            siteOrdering,
        _DataSetFilter const*     theFilter,
        long           *        lNodeFlags,
        _Vector*     lNodeResolutions,
        long siteFrom,
        long siteTo,
        long catID,
        hyFloat* storageVec = nil);
#endif

    // --------------------------


    long      * nodeStates;
    char      * nodeMarkers;

    hyFloat* rootIChildrenCache,
                * marginalLikelihoodCache;

    _AVLListXL* aCache;

    long        categoryCount;

protected:


    bool        IntPopulateLeaves   (_DataSetFilter const*, long) const;

    virtual     void                PreTreeConstructor                  (bool);
    virtual     void                PostTreeConstructor                 (bool);


    // all of the following members exist to speed-up the pruning algorithm
    // the are created by calling the function set up
    _SimpleList flatLeaves,
                flatNodes,
                leftiNodes,
                flatParents,
                topLevelNodes,
                topLevelLeftL,
                topLevelRightL,
                forceRecalculationOnTheseBranches,
                nodesToUpdate;

};

#define       _HY_TREE_TRAVERSAL_MASK      0x07 // 111 binary
#define       _HY_TREE_TRAVERSAL_SKIP_ROOT 0x10 //
#define       _HY_TREE_TRAVERSAL_LEAVES    0x20

class _TreeIterator {
private:
  
  node_iterator<long> iterator;
  int                 flags;
  _SimpleList         history;
  _TheTree   const    *source_tree;
  
  
public:
  _TreeIterator (_TheTree const * source, int traversal_type);
  ~_TreeIterator (void);
  void                    Reset (void);
  
  _CalcNode *             Next (void);
  _CalcNode *             Current (void) const;
  long                    Depth (void) const;
  const _SimpleList&      History (void) const;
  bool  IsAtLeaf          (void) const;
  bool  IsAtRoot          (void) const;
  node <long>* GetNode    (void) const { return iterator.Current(); }
  
};

/*----------------------------------------------------------------------------------------------------------*/

template <class data_type> _CalcNode* map_node_to_calcnode (node<data_type>* n) {
  if (n) {
    return (_CalcNode*)(variablePtrs.GetItem(n->in_object));
  }
  return nil;
}

/*----------------------------------------------------------------------------------------------------------*/

extern hyTreeDefinitionPhase     isDefiningATree;

extern _String  expectedNumberOfSubs,
       stringSuppliedLengths,
       includeModelSpecs,
       treeOutputAVL,
       treeOutputLayout;
#ifdef _OPENMP
#include "omp.h"
#endif

#endif
