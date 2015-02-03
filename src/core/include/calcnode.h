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

#ifndef     __CALCNODE__
#define     __CALCNODE__

#include "parser.h"
#include "classes.h"
#include "site.h"


#define UNROOTED                        0
#define ROOTED_LEFT                     1
#define ROOTED_RIGHT                    2


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
						_Parameter* eiNodeCache);


	double launchmdsocl(	_SimpleList& updateNodes,
							_SimpleList& flatParents,
							_SimpleList& flatNodes,
							_SimpleList& flatCLeaves,
							_SimpleList& flatLeaves,
							_SimpleList& flatTree,
							_Parameter* theProbs,
							_SimpleList& theFrequencies,
							long* lNodeFlags,
							_SimpleList& taggedInternals,
							_GrowingVector* lNodeResolutions);
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

    virtual void        Duplicate       (BaseRef);

    virtual long        FreeUpMemory    (long);

    void                InitializeCN    ( _String&, int, _VariableContainer*, _AVLListXL * = nil);

    virtual BaseRef     makeDynamic     (void);
    // creates a dynamic copy of this object

    virtual BaseRef     toStr           (void);
    // converts this object to string

    _Parameter&         operator[]      (unsigned long);
    // access the i-th element of the
    // probabilities (i = 0..codeBase-1)

    void                SetCodeBase      (int);
    // change the codeBase value for this node
    // this will resize the vector used to handle frequencies

    bool                RecomputeMatrix  (long = 0, long = 1,_Matrix* = nil, _List* = nil, _SimpleList* = nil, _List* = nil);
    // reexponentiate the transition matrix and
    // store it in compExp.
    // return TRUE if the matrix is an explicit exponential form

    virtual bool        HasChanged       (void);
    virtual bool        NeedToExponentiate(long = -1L);
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

    _Parameter          GetProbs        (long k) {
        return theProbs[k];
    }
    _Parameter*         GetProbs        (void) {
        return theProbs;
    }

    void                SetCompExp      (_Matrix*, long = -1);
    void                SetCompMatrix   (long);
    _Matrix*            GetCompExp      (long catID = -1, bool = false);

    _Formula*           RecurseMC       (long , node<long>* , bool first = false, char rooted = UNROOTED);

    long                GetCodeBase     (void) {
        return cBase;
    }

    _Parameter          BranchLength    (void);
    virtual long        SetDependance   (long);

    node<long>*         LocateMeInTree  (void);
    // return the tree structure node corresponing to this one...
    long                ConvertToSimpleMatrix (void);
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

public:
    _Parameter*     theProbs;       // list of transitional probabilities
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

struct      nodeCoord {

    _Parameter  h,
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

    _String     branchName,
                branchTag;

}; // used for tree imaging

//_______________________________________________________________________________________________

typedef bool _HYTopologyTraversalFunction (node<long>*, Ptr);

//_______________________________________________________________________________________________

class _TheTree; // forward declaration for xlc


//_______________________________________________________________________________________________

class _TreeTopology: public _CalcNode
{

protected:

    virtual void            PreTreeConstructor                  (bool);
    virtual bool            MainTreeConstructor                 (_String&,bool = true);
    virtual void            PostTreeConstructor                 (bool);
    node<long>*     prepTree4Comparison                 (_List&, _SimpleList&, node<long>* = nil);
    void            destroyCompTree                     (node<long>*);
    _List*          SplitTreeIntoClustersInt            (node<long>*, _List*, _AVLListX&, long, long);
    char            internalTreeCompare                 (node<long>*, node<long>*, _SimpleList*, char, long, node<long>*, _TreeTopology*, bool = false);
    char            internalNodeCompare                 (node<long>*, node<long>*, _SimpleList&, _SimpleList*, bool, long, node<long>*, _TreeTopology*, bool = false);
    virtual _PMathObj       FlatRepresentation                  (void);
    void            FindCOTHelper                       (node<long>*, long, _Matrix&, _Matrix&, _Matrix&, _List&, _AVLListX&, _Parameter);
    void            FindCOTHelper2                      (node<long>*, _Matrix&, _Matrix&, _AVLListX&, node<long>*, _Parameter);
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

    node<long>      *theRoot,
         *currentNode;

    _List           flatTree,
                    flatCLeaves;

    char            rooted;

    virtual void            toFileStr                           (FILE*);
    virtual BaseRef         toStr                               (void);
    void            RerootTreeInternalTraverser         (long, bool,_String&, long  = -1, bool = false);

    _TreeTopology                       (void);
    _TreeTopology                       (_String, _String&, bool = true);
    _TreeTopology                       (_String*);
    _TreeTopology                       (_TheTree*);

    virtual                 ~_TreeTopology                      (void);

    virtual  _FString*      Compare                             (_PMathObj);
    virtual  BaseRef        makeDynamic                         (void);
    node<long>* CopyTreeStructure                   (node<long>*, bool);
    virtual  bool           FinalizeNode                        (node<long>*, long, _String&, _String&, _String&, _String* = NULL);


    bool            IsCurrentNodeATip                   (void);
    bool            IsCurrentNodeTheRoot                (void);
    bool            IsDegenerate                        (void);

    virtual _PMathObj       Execute                             (long, _PMathObj = nil , _PMathObj = nil, _hyExecutionContext* context = _hyDefaultExecutionContext);
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
    virtual _AssociativeList*
    FindCOT                             (_PMathObj);

    node<long>      *FindNodeByName                     (_String*);
    /*

     20091006: SLKP

     return the node with the name supplied by the argument
     or nil if no such node exists

     */

    void            DepthWiseT                          (bool = false, _HYTopologyTraversalFunction* = nil, Ptr = nil);
    void            DepthWiseTRight                     (bool = false);
    void            DepthWiseTLevel                     (long& level, bool = false);
    void            StepWiseT                           (bool = false, _HYTopologyTraversalFunction* = nil, Ptr = nil);
    void            StepWiseTLevel                      (long&, bool = false);
    void            LeafWiseT                           (bool = false);
    
    _List*          MapNodesToModels                    (void);

    virtual void            GetNodeName                         (node<long> *, _String&, bool = false);
    virtual _String*        GetNodeModel                        (node<long> *);
    virtual void            GetBranchLength                     (node<long> *, _String&, bool = false);
    // SLKP 20100901:
    //               added a boolean flag to ask to return branch length expression (if true) (returns "" for topologies)
    //               just the numeric value (if false)


    virtual void            GetBranchLength                     (node<long> *, _Parameter&);
    virtual void            GetBranchValue                      (node<long> *, _String&);
    virtual void            GetBranchVarValue                   (node<long> *, _String&, long);
    virtual void            PasteBranchLength                   (node<long> *, _String&, long, _Parameter factor = 1.);

    node<long>&     GetRoot                             (void) {
        return *theRoot;
    }
    void            SetRoot                             (node<long>* r) {
        theRoot = r;
    }
    node<long>&     GetCurrentNode                      (void) {
        return *currentNode;
    }
    void            SubTreeString                       (_String&, bool = false, long = -1, _AVLListXL* = nil);

    _String         CompareTrees                        (_TreeTopology*);
    _String         MatchTreePattern                    (_TreeTopology*);
    virtual _PMathObj       TipName                             (_PMathObj);
    virtual _PMathObj       BranchName                          (_PMathObj, bool = false, _PMathObj = nil);
    virtual _PMathObj       BranchLength                        (_PMathObj);
    virtual _PMathObj       RerootTree                          (_PMathObj);
    _List*          SplitTreeIntoClusters               (unsigned long, unsigned long);
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
extern _Parameter scalingLogConstant;
#endif

//_______________________________________________________________________________________________

class _TheTree: public _TreeTopology
{

// theModel matrix of _TheTree contains the column matrix of probabilities, which is computed
// based on the DataSetFilter passed on to the tree at the initialization stage

public:

    _TheTree ();                                                // default constructor - doesn't do much
    _TheTree (_String name, _String& parms, bool = true);       // builds a tree from a string
    _TheTree (_String name, _TreeTopology*);                    // builds a tree from a tree topology
    _TheTree (_String name, _TheTree*);                    // builds a tree from another tree


    virtual                 ~_TheTree                   (void);
    virtual bool            HasChanged                  (void);
    virtual void            MarkDone                    (void);
    bool            HasChanged2                 (void);

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

    virtual  bool           FinalizeNode                (node<long>*, long, _String&, _String&, _String&, _String* = NULL);
    virtual  BaseRef        makeDynamic                 (void);

    virtual  BaseRef        makeDynamicCopy             (_String*);
    node<long>* DuplicateTreeStructure      (node<long>*, _String*, bool);
    virtual  BaseRef        toStr                       (void);
    virtual unsigned long           ObjectClass                 (void) {
        return TREE;
    }

    virtual  _PMathObj      Execute                     (long, _PMathObj = nil , _PMathObj = nil, _hyExecutionContext* context = _hyDefaultExecutionContext);
    virtual  _PMathObj      TEXTreeString               (_PMathObj);
    virtual  _PMathObj      PlainTreeString             (_PMathObj,_PMathObj);

    virtual  void           GetNodeName                 (node<long> *, _String&, bool = false);
    virtual  void           GetBranchLength             (node<long> *, _String&, bool = false);
    virtual  void           GetBranchLength             (node<long> *, _Parameter&);
    virtual  void           GetBranchValue              (node<long> *, _String&);
    virtual  _String*       GetBranchSpec               (node<long> *);
    virtual  void           GetBranchVarValue           (node<long> *, _String&, long);
    virtual _String*        GetNodeModel                (node<long> *);
    
    void            InitializeTreeFrequencies   (_Matrix *, bool = false);

    _Parameter      ReleafTreeAndCheck          (_DataSetFilter*, long, bool, long categID = -1);
    _Parameter      ReleafTreeAndCheckChar4     (_DataSetFilter*, long, bool, long categID = -1);

    _Parameter      ReleafTree                  (_DataSetFilter*,long,long,long,long);
    _Parameter      ReleafTreeDegenerate        (_DataSetFilter*,long);

    _Parameter      ReleafTreeCache             (_DataSetFilter*,long,long,long,long,long);
#if USE_SCALING_TO_FIX_UNDERFLOW
    _Parameter      ThreadReleafTreeCache       (_DataSetFilter*,long,long,long,long,long,long offset = 0,long fixAttempt = 0, _Parameter = 690.);
    _Parameter      doScaling                   (_DataSetFilter*,long,long,long,long,_Parameter, bool, bool);
#else
    _Parameter      ThreadReleafTreeCache       (_DataSetFilter*,long,long,long,long,long,long offset = 0);
#endif

    _Parameter      Process3TaxonNumericFilter  (_DataSetFilterNumeric*, long = 0);


    void            ThreadMatrixUpdate          (long, bool);
    void            SerialMatrixUpdate          (long, bool);
    void            MatrixCacheUpdate           (void);

    _Parameter      ThreadReleafTreeCharCache   (_DataSetFilter*,long,long,long,long,long,long offset = 0);
    _Parameter      ReleafTreeCharDegenerate    (_DataSetFilter*,long);
    _Parameter      ReleafTreeChar4             (_DataSetFilter*,long,long,long,long,long);
    _Parameter      ReleafTreeChar4Degenerate   (_DataSetFilter*,long);

    _Parameter      ThreadReleafTreeChar4       (_DataSetFilter*,long,long,long,long,long,long offset = 0);
    _Parameter      ReleafTreeChar4             (_DataSetFilter*,long,long,long,long);

    _Parameter  Probij                          (long, long, _CalcNode*);
    _Parameter  PruneTree                       (long categID = -1);
    _Parameter  PruneTreeChar                   (long categID = -1);
    _Parameter  PruneTreeCharCache              (long categID = -1);
    _Parameter  PruneTreeChar4                  (long categID = -1);
    _Parameter  PruneTreeChar4Cache             (long categID = -1);

    _List*      RecoverAncestralSequences       (_DataSetFilter*, _SimpleList&, _List&, _Parameter*, _Parameter*, long, long*, _GrowingVector*, bool = false);
    void        RecoverNodeSupportStates        (_DataSetFilter*, long, long, _Matrix&);
    void        RecoverNodeSupportStates2       (node<long>*,_Parameter*,_Parameter*,long);
    _List*      SampleAncestors                 (_DataSetFilter*, node<long>*);
    void        PurgeTree                       (void);

    long        ComputeReleafingCost            (_DataSetFilter*, long, long, _SimpleList* = nil, long = 0);
    long        ComputeReleafingCostChar        (_DataSetFilter*, long, long);
    void        DumpingOrder                    (_DataSetFilter*, _SimpleList&);
    void        SetTreeCodeBase                 (long);
    long        IsLinkedToALF                   (long&);

    bool        HasCache                        (void) {
        return topLevelNodes.lLength>0;
    }

    long        GetLeafCount                    (void) {
        return flatLeaves.lLength;
    }

    long        GetINodeCount                   (void) {
        return flatNodes.lLength    ;
    }

    void        ScanAndAttachVariables          (void);
    void        ScanForVariables                (_AVLList& l, _AVLList& l2, _AVLListX* tagger = nil, long weight = 0);
    void        ScanForDVariables               (_AVLList& l, _AVLList& l2);
    void        ScanForGVariables               (_AVLList&, _AVLList&, _AVLListX* tagger = nil, long weight = 0);
    void        ScanForCVariables               (_AVLList&);
    void        MolecularClock                  (_String&, _List&);

    void        SetUp                           (void);
    void        SetUpMatrices                   (long);
    void        CleanUpMatrices                 (void);
    void        BuildTopLevelCache              (void);
    void        KillTopLevelCache               (void);
    void        SetCompMatrices                 (long);

    virtual void        ClearConstraints                (void);

    bool        FindScalingVariables            (_SimpleList&);
    bool        HaveStringBranchLengths         (void);
    void        AssignLabelsToBranches          (node<nodeCoord>*, _String*, bool);

    node<nodeCoord>*
    AlignedTipsMapping              (bool first = false, bool respectRoot = true);

    void        AlignNodes                      (node<nodeCoord>*);

    node<nodeCoord>*
    ScaledBranchMapping             (node<nodeCoord>* , _String*, long, long&, char);

    node<nodeCoord>*
    RadialBranchMapping             (node<long>* , node<nodeCoord>*, _String*, _Parameter, long&, _Parameter&, char);

    void        ScaledBranchReMapping           (node<nodeCoord>*, _Parameter);
    char&       RootedFlag                      (void) {
        return rooted;
    }

    nodeCoord   TreeTEXRecurse                  (node<nodeCoord>*,_String&,_Parameter,_Parameter,long,long);
    void        TreePSRecurse                   (node<nodeCoord>*,_String&,_Parameter,_Parameter,long,long,long,long,_AssociativeList* = nil, char = 0, _Parameter* = nil);

    bool        AllBranchesHaveModels           (long);
    void        ScanSubtreeVars                 (_List&, char, _CalcNode*);
    void        BuildINodeDependancies          (void);
    void        AllocateResultsCache            (long);
    long        CountTreeCategories             (void);
    void        CompileListOfModels             (_SimpleList&);

    void        MarkMatches                     (_DataSetFilter*,long,long);
    long        GetLowerBoundOnCost             (_DataSetFilter*);
    long        GetLowerBoundOnCostWithOrder    (_DataSetFilter*,_SimpleList*);
    _SimpleList&GetLeftINodes                   (void) {
        return leftiNodes;
    }
    bool        MatchLeavesToDF                 (_SimpleList&, _DataSetFilter*, bool);
    virtual void
    RemoveModel                     (void);
    _String*    TreeUserParams                  (void);


    _String     CompareSubTrees                 (_TheTree*, node<long>*);
    _String     FindMaxCommonSubTree            (_TheTree*, long&, _List*);
    void        WeightedCharacterDifferences    (_Parameter, _Matrix*, _Matrix*, long = -1);
    void        AddNodeNamesToDS                (_DataSet*, bool, bool, char);
    // if the
    _Parameter  PSStringWidth                   (_String&);

    _Parameter  DetermineBranchLengthGivenScalingParameter
    (long, _String&, char);

    _AVLListX*  ConstructNodeToIndexMap         (bool);
    // 20090206: SLKP
    // makes an AVL of with keys storing memory addresses of node<long> tree nodes
    // and values showing the order in either flatLeaves (bool = false) or flatNodes (bool = true)

    void        MapPostOrderToInOderTraversal   (_SimpleList&, bool = true);
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

    _Parameter  VerySimpleLikelihoodEvaluator   (_SimpleList&            updateNodes,
            _DataSetFilter*      theFilter,
            _Parameter*          iNodeCache,
            long       *             lNodeFlags,
            _GrowingVector*      lNodeResolutions);

#ifdef MDSOCL
			_Parameter OCLLikelihoodEvaluator (			_SimpleList&	     updateNodes, 
                                                        _DataSetFilter*		 theFilter,
                                                        _Parameter*			 iNodeCache,
                                                         long	   *		 lNodeFlags,
                                                        _GrowingVector*		 lNodeResolutions,
														_OCLEvaluator& OCLEval);
#endif

#ifdef  _SLKP_LFENGINE_REWRITE_
    void            SampleAncestorsBySequence       (_DataSetFilter*, _SimpleList&, node<long>*, _AVLListX*, _Parameter*, _List&, _SimpleList*, _List&, _Parameter*, long);

    _Parameter      ComputeTreeBlockByBranch        (_SimpleList&, _SimpleList&, _SimpleList*, _DataSetFilter*, _Parameter*, long*, _Parameter*, _GrowingVector*, long&, long, long, long = -1, _Parameter* = nil, long* = nil, long = -1, long * = nil);
    long            DetermineNodesForUpdate         (_SimpleList&,  _List* = nil, long = -1, long = -1, bool = true);
    void            ExponentiateMatrices            (_List&, long, long = -1);
    void            FillInConditionals              (_DataSetFilter*, _Parameter*,  _SimpleList*);

    void            ComputeBranchCache              ( _SimpleList&,
            long nodeID,
            _Parameter*         cache,
            _Parameter*         iNodeCache,
            _DataSetFilter*     theFilter,
            long           *        lNodeFlags,
            _Parameter*         scalingAdjustments,
            long*                   siteCorrectionCounts,
            _GrowingVector*     lNodeResolutions,
            long&                   overallScaler,
            long                    siteFrom,
            long                    siteTo,
            long                    catID,
            _SimpleList*            = nil,
            _Parameter*         = nil
                                                    );

    _Parameter          ComputeLLWithBranchCache         (
        _SimpleList&            siteOrdering,
        long                    brID,
        _Parameter*         cache,
        _DataSetFilter*     theFilter,
        long                    siteFrom,
        long                    siteTo,
        long                    catID,
        _Parameter*         storageVec = nil
    );

    _Parameter          ComputeTwoSequenceLikelihood    (
        _SimpleList&            siteOrdering,
        _DataSetFilter*     theFilter,
        long           *        lNodeFlags,
        _GrowingVector*     lNodeResolutions,
        long siteFrom,
        long siteTo,
        long catID,
        _Parameter* storageVec = nil);
#endif

    // --------------------------


    long      * nodeStates;
    char      * nodeMarkers;

    _Parameter* rootIChildrenCache,
                * marginalLikelihoodCache;

    _AVLListXL* aCache;

    long        categoryCount;

protected:


    bool        IntPopulateLeaves   (_DataSetFilter*, long, long);

    _Parameter  ConditionalBranchLikelihood
    (node<long>* , node<long>* , _Parameter* , _Parameter* , long, long);

    _Parameter  ConditionalNodeLikelihood
    (node<long>* , node<long>* , _Parameter* , _Parameter* , long ,long);

    _List*      MapCBaseToCharacters(_DataSetFilter*, bool = true);


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

extern char     isDefiningATree;
extern _String  expectedNumberOfSubs,
       stringSuppliedLengths,
       includeModelSpecs,
       treeOutputAVL,
       treeOutputLayout;
#ifdef _OPENMP
#include "omp.h"
#endif

#endif
