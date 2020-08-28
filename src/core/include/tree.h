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

#ifndef     __TREE__
#define     __TREE__

#include    "topology.h"
#include    "dataset_filter.h"
#include    "dataset_filter_numeric.h"
#include    "vector.h"

#define     HY_REPLACE_BAD_BRANCH_LENGTH_WITH_THIS  0.000000001
#define     HY_BRANCH_SELECT                    0x01
#define     HY_BRANCH_DESELECT                  0xFFFFFFFE

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

class _TheTree: public _TreeTopology
{

// theModel matrix of _TheTree contains the column matrix of probabilities, which is computed
// based on the DataSetFilter passed on to the tree at the initialization stage

public:

    _TheTree ();                                                // default constructor - doesn't do much
    _TheTree (_String const &name, _String const& parms, bool = true);       // builds a tree from a string
    _TheTree (_String const &name, _TreeTopology*);                    // builds a tree from a tree topology
    _TheTree (_String const& name, _TheTree*);                    // builds a tree from another tree


    virtual                 ~_TheTree                   (void);
    virtual bool            HasChanged                  (bool = false);
    virtual void            MarkDone                    (void);
    bool            HasChanged2                 (void);

    virtual  _String        FinalizeNode                (node<long>*, long, _String, _String const&, _String&, _TreeTopologyParseSettings const&);
    virtual  BaseRef        makeDynamic                 (void) const;

    virtual  BaseRef        makeDynamicCopy             (_String const*) const;
    node<long>* DuplicateTreeStructure      (node<long>*, _String const*, bool) const;
    virtual  BaseRef        toStr                       (unsigned long = 0UL);
    virtual unsigned long           ObjectClass                 (void) const {
        return TREE;
    }

    virtual HBLObjectRef       ExecuteSingleOp            (long, _List* = nil, _hyExecutionContext* context = _hyDefaultExecutionContext, HBLObjectRef cache = nil);
    virtual  HBLObjectRef      PlainTreeString            (HBLObjectRef,HBLObjectRef, HBLObjectRef cache);

    virtual _String const   GetNodeName                 (node<long> *, bool = false) const;
    virtual _String const   GetNodeName                 (_CalcNode const *cn, bool = false) const;

    virtual _String const   GetBranchLengthString       (node<long> *, bool get_expression = false) const;
    virtual  hyFloat        GetBranchLength             (node<long> *) const ;
    virtual  _String const  GetBranchValue              (node<long> *) const ;
    virtual  _String const  GetBranchVarValue           (node<long> *, long) const ;
    virtual _String const*  GetNodeModel                (node<long> *) const;
    
    void            InitializeTreeFrequencies   (_Matrix *, bool = false);


    hyFloat      Process3TaxonNumericFilter  (_DataSetFilterNumeric*, long = 0);



    _List*      RecoverAncestralSequences       (_DataSetFilter const*, _SimpleList const&, _List const&, hyFloat *, hyFloat const*, long, long*, _Vector*, bool = false);
    void        RecoverNodeSupportStates        (_DataSetFilter const*, long, _Matrix&);
    void        RecoverNodeSupportStates2       (node<long>*,hyFloat*,hyFloat*,long, _AVLListX &);
    _List*      SampleAncestors                 (_DataSetFilter*, node<long>*);
    void        PurgeTree                       (void);

    long        ComputeReleafingCost            (_DataSetFilter const*, long, long, _SimpleList* = nil, long = 0, _SimpleList const* = nil) const;
    long        ComputeReleafingCostChar        (_DataSetFilter const*, long, long, _SimpleList const*) const;
    void        DumpingOrder                    (_DataSetFilter*, _SimpleList&);
    void        SetTreeCodeBase                 (long);
    long        IsLinkedToALF                   (long&) const;

    bool        HasCache                        (void) const {
        return topLevelNodes.lLength>0;
    }

    long        GetLeafCount                    (void) const {
        return flatLeaves.lLength;
    }

    long        GetINodeCount                   (void) const {
        return flatNodes.lLength    ;
    }
    
    const _SimpleList& get_flat_nodes (void) {return flatNodes;}

    void        ScanAndAttachVariables          (void) const;
    void        ScanContainerForVariables       (_AVLList& l, _AVLList& l2, _AVLListX* tagger = nil, long weight = 0, _AVLListX* map_variables_to_nodes = nil) const;
    void        ScanForDVariables               (_AVLList& l, _AVLList& l2) const;
    void        ScanForGVariables               (_AVLList&, _AVLList&, _AVLListX* tagger = nil, long weight = 0) const;
    void        ScanForCVariables               (_AVLList&) const;
    void        MolecularClock                  (_String const&, _List&) const;

    void        SetUp                           (void);
    void        SetUpMatrices                   (long);
    void        CleanUpMatrices                 (void);
    //void        BuildTopLevelCache              (void);
    void        SetCompMatrices                 (long) const;

    virtual void        ClearConstraints                (void);

    void        AssignLabelsToBranches          (node<nodeCoord>*, _String*, bool);

    node<nodeCoord>*
    AlignedTipsMapping                          (node<long>*, hyFloat& current_offset, bool first = false, bool respectRoot = true) const;

    void        AlignNodes                      (node<nodeCoord>*) const;

    node<nodeCoord>*
    ScaledBranchMapping             (node<nodeCoord>* , _String*, long, long&, hyTopologyBranchLengthMode) const;

    node<nodeCoord>*
    RadialBranchMapping             (node<long>* , node<nodeCoord>*, _String*, hyFloat, long&, hyFloat&, hyTopologyBranchLengthMode);

    void        ScaledBranchReMapping           (node<nodeCoord>*, hyFloat) const;
    char&       RootedFlag                      (void) {
        return rooted;
    }

     void        TreePSRecurse                   (node<nodeCoord>*,_StringBuffer&,hyFloat,hyFloat,long,long,long,long,const _TreeTopologyParseSettings&,_AssociativeList* = nil, char = 0, hyFloat* = nil) const;

    bool        AllBranchesHaveModels           (long) const;
    void        AllocateResultsCache            (long);
    long        CountTreeCategories             (void);
    void        CompileListOfModels             (_SimpleList&);

    _SimpleList&GetLeftINodes                   (void) {
        return leftiNodes;
    }
    virtual void    RemoveModel                     (void);

  
  
    void        AddNodeNamesToDS                (_DataSet*, bool, bool, char) const;
    // if the
    static hyFloat  PSStringWidth                   (_String const&);

    hyFloat  DetermineBranchLengthGivenScalingParameter (long, _String&, hyTopologyBranchLengthMode) const;

    _AVLListX*  ConstructNodeToIndexMap         (bool) const;
    // 20090206: SLKP
    // makes an AVL of with keys storing memory addresses of node<long> tree nodes
    // and values showing the order in either flatLeaves (bool = false) or flatNodes (bool = true)

    void        MapPostOrderToInOrderTraversal   (_SimpleList&, bool = true) const;
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
  

#ifdef  _SLKP_LFENGINE_REWRITE_
    void            SampleAncestorsBySequence       (_DataSetFilter const*, _SimpleList const&, node<long>*, _AVLListX const*, hyFloat const*, _List&, _SimpleList*, _List&, hyFloat const*, long);

    hyFloat      ComputeTreeBlockByBranch        (_SimpleList&, _SimpleList&, _SimpleList*, _DataSetFilter const*, hyFloat*, long*, hyFloat*, _Vector*, long&, long, long, long = -1, hyFloat* = nil, long* = nil, long = -1, long * = nil);
    long            DetermineNodesForUpdate         (_SimpleList&,  _List* = nil, long = -1, long = -1, bool = true, _AVLListX * var_mapping = nil, _AVLList * changed_variables = nil);
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
            _Vector const*     lNodeResolutions,
            long&                   overallScaler,
            long const                    siteFrom,
            long                    siteTo,
            long const                   catID,
            _SimpleList const *            = nil,
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


    _AVLListXL* aCache;

    long        categoryCount;

protected:
  
    void        delete_associated_calcnode (node<long>*) const;

    virtual void _RemoveNodeList (_SimpleList const& list);

    bool        IntPopulateLeaves   (_DataSetFilter const*, long) const;

    virtual     void                PreTreeConstructor                  (bool);
    virtual     void                PostTreeConstructor                 (bool, _AssociativeList*);


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
    
    static      hyFloat _timesCharWidths[256],
                         _maxTimesCharWidth;
    
    static      _String const kTreeOutputLabel,
                              kTreeOutputTLabel,
                              kTreeOutputFSPlaceH;

};


#define    _HY_BITMASK_WIDTH_ (8*sizeof (unsigned long))

struct      bitMasks {
    unsigned long masks[_HY_BITMASK_WIDTH_];
    bitMasks (void) {
        unsigned long aBit = 1UL;
        for (unsigned long k=0UL; k<_HY_BITMASK_WIDTH_; k++) {
            masks[k] = aBit;
            aBit = aBit << 1;
        }
    }
};

extern bitMasks bitMaskArray;

#ifdef _OPENMP
#include "omp.h"
#endif

#endif
