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

#ifndef     __TOPOLOGY__
#define     __TOPOLOGY__

#include "calcnode.h"
#include "associative_list.h"
#include "fstring.h"

#define UNROOTED                        0
#define ROOTED_LEFT                     1
#define ROOTED_RIGHT                    2

#define HY_REPLACE_BAD_BRANCH_LENGTH_WITH_THIS  0.000000001


#define fGetNodeStringForTreeName   0x01
#define fGetNodeStringForTreeModel  0x02

class _TheTree;

struct _TreeTopologyParseSettings {
    
    _TreeTopologyParseSettings (void) {
        inode_prefix = "Node";
        auto_convert_lengths = false;
        accept_user_lengths = true;
        ingore_user_inode_names = false;
        parser_cache = nil;
    }
  
    ~_TreeTopologyParseSettings () {
      if (parser_cache) {
        parser_cache->ClearFormulasInList();
        DeleteObject (parser_cache->dataList);
        DeleteObject (parser_cache);
      }
    }
  
    void AllocateCache (void) {
      DeleteObject (parser_cache);
      parser_cache = new _AVLListX (new _SimpleList);
    }
  
    _String inode_prefix;
    bool    auto_convert_lengths,
            accept_user_lengths,
            ingore_user_inode_names;
  
    _AVLListX * parser_cache;
};

enum   hyTopologyBranchLengthMode {
    kTopologyBranchLengthLocalParameter = 0L,
    kTopologyBranchLengthNone         = -1L,
    kTopologyBranchLengthExpectedSubs = -3L,
    kTopologyBranchLengthUserLengths  = -2L
};

//_______________________________________________________________________________________________

class _TreeTopology: public _CalcNode {

protected:
    
    
    template <class CallBack> node<long>* ConditionalTraverser (CallBack && cb, bool do_root) {
        node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
        while (node<long>* iterator = ni.Next()) {
            if (!do_root && iterator->is_root() ) {
                return nil;
            }
            if (cb (iterator, ni)) {
                return iterator;
            }
        }
        return nil;
    }

    virtual void                PreTreeConstructor                 (bool);
    virtual _AssociativeList*  MainTreeConstructor                 (_String const&, _TreeTopologyParseSettings & settings, bool = true, _AssociativeList* mapping = nil);
    virtual void            PostTreeConstructor                    (bool, _AssociativeList* );
    node<long>*     prepTree4Comparison                            (_List&, _SimpleList&, node<long>* = nil) const;
    void            destroyCompTree                                (node<long>*) const;
    _List*          SplitTreeIntoClustersInt                       (node<long>*, _List*, _AVLListX&, long, long) const;
    char            internalTreeCompare                            (node<long>*, node<long>*, _SimpleList*, char, long, node<long>*, _TreeTopology const*, bool = false) const;
    
    /**
        Compare two nodes on the same (subset) of leaves and report if they create the same subpartition
     
        @param n1 the first node to compare
        @param n2 the second node to compare
        @param subTreeMap ?
        @param reindexer ?
        @param cangoup ?
        @param totalSize ?
        @param n22 ?
        @param tree2 ?
        @param isPattern ?
     
        @return whether or not the nodes are equal
    */
    bool            internalNodeCompare                 (node<long>* n1, node<long>*, _SimpleList& subTreeMap, _SimpleList* reindexer, bool cangoup, long totalSize, node<long>* n22, _TreeTopology const* tree2, bool isPattern = false) const;
    virtual HBLObjectRef       FlatRepresentation                  (HBLObjectRef cache);
    void            FindCOTHelper                       (node<long>*, long, _Matrix&, _Matrix&, _Matrix&, _List&, _AVLListX&, hyFloat);
    void            FindCOTHelper2                      (node<long>*, _Matrix&, _Matrix&, _AVLListX&, node<long>*, hyFloat);
    static          const   _TreeTopologyParseSettings  CollectParseSettings (void);
    void            AddANode                            (HBLObjectRef);
    /*

     20091006: SLKP

     given an AVL with at least three key -> string pairs
        "NAME" -> string (not a currently used node name)
        "WHERE" -> string (an existing node)
        "PARENT" -> string (not a currently used node name)



     */

    void            RemoveANode                          (HBLObjectRef);
    
    /*
     
        Delete a node from the tree by name
     
     */
    
  
    virtual void _RemoveNodeList (_SimpleList const& list);

public:
    // class constants
    
    static const _String kCompareEqualWithReroot,
                         kCompareEqualWithoutReroot,
                         kCompareUnequalToplogies,
                         kCompareUnequalLabelSets,
                         kMeta;

    node<long>      *theRoot;
  
    _List           flatTree,
                    flatCLeaves;

    char            rooted;

    virtual void            toFileStr                           (FILE*, unsigned long);
    virtual BaseRef         toStr                               (unsigned long = 0UL);
    void            RerootTreeInternalTraverser         (node<long>* iterator, long, bool,_StringBuffer&, _TreeTopologyParseSettings const& settings,  hyTopologyBranchLengthMode branch_length_mode, long variable_ref  = -1L, bool = false) const;

    _TreeTopology                       (void);
    _TreeTopology                       (_String const&, _String const&, bool = true, _AssociativeList* mapping = nil);
    _TreeTopology                       (_String const*);
    _TreeTopology                       (_TheTree const*);
    _TreeTopology                       (const _TreeTopology& );
    _TreeTopology const&                operator = (const _TreeTopology& );

    virtual                 ~_TreeTopology                      (void);

    virtual  _FString*      Compare                             (HBLObjectRefConst, HBLObjectRef cache) const;
    virtual  BaseRef        makeDynamic                         (void) const;
    node<long>* CopyTreeStructure                   (node<long>*, bool) const;
    virtual  _String           FinalizeNode                        (node<long>*, long, _String, _String const&, _String&, _TreeTopologyParseSettings const & settings);


    virtual HBLObjectRef       ExecuteSingleOp                     (long, _List* = nil, _hyExecutionContext* context = _hyDefaultExecutionContext,HBLObjectRef cache = nil);
    virtual void            EdgeCount                           (long&, long&) const;
    // SLKP 20100827: a utility function to count edges in a tree
    //              : note that the root node WILL be counted as an internal node
    //              : writes [leaf count, internal node count] into the arguments

    virtual HBLObjectRef       TipCount                            (HBLObjectRef cache);
    virtual HBLObjectRef       BranchCount                         (HBLObjectRef cache);
    virtual HBLObjectRef       AVLRepresentation                   (HBLObjectRef, HBLObjectRef cache);
    virtual unsigned long      ObjectClass                         (void) const {
        return TOPOLOGY;
    }
    virtual bool IsDegenerate(void) const { return theRoot && theRoot->get_num_nodes() == 1L && theRoot->go_down(1)->get_num_nodes() == 0L; }
  
 
    virtual _AssociativeList*           FindCOT                             (HBLObjectRef,HBLObjectRef cache);

    virtual HBLObjectRef                MaximumParsimony                    (HBLObjectRef,HBLObjectRef cache);

    node<long>      *FindNodeByName                     (_String const*) const;
    /*

     20091006: SLKP

     return the node with the name supplied by the argument
     or nil if no such node exists

     */

    
    _List*          MapNodesToModels                    (void);

    virtual _String const GetNodeName                         (node<long> *, bool = false) const;
    virtual _String const *GetNodeModel                        (node<long> *) const;
    virtual _String const GetBranchLengthString                (node<long> *, bool get_expression = false) const;
    // SLKP 20100901:
    //               added a boolean flag to ask to return branch length expression (if true) (returns "" for topologies)
    //               just the numeric value (if false)


    virtual hyFloat         GetBranchLength                     (node<long> *) const;
    virtual _String const   GetBranchValue                      (node<long> *) const;
    virtual _String const   GetBranchVarValue                   (node<long> *, long) const;
    virtual _String const   GetNodeStringForTree                (node<long> *, int flags) const;
    virtual void            PasteBranchLength                   (node<long> * node, _StringBuffer & result , hyTopologyBranchLengthMode const mode, long variable_reference , hyFloat factor = 1.) const;

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
    void            SubTreeString                       (node<long>* root, _StringBuffer & result, _TreeTopologyParseSettings const& settings, bool all_names = false, hyTopologyBranchLengthMode mode = kTopologyBranchLengthNone, long branch_length_variable = -1, _AVLListXL * substitutions  = nil) const;

    virtual HBLObjectRef       RandomizeTips            (HBLObjectRef, HBLObjectRef cache);
    /**
        Shuffle the order of tips in the tree by permuting the order of children of
        each node with certain probability
     
        The argument controls the rate at which this shuffling occurs
        if it's a number (between 0 and 1), then each internal node will be shuffled
        with this probability
     
        Returns a Topology object representing the reshuffled tree
     */

    _String                     CompareTrees                        (_TreeTopology*) const;
    const _String               MatchTreePattern                    (_TreeTopology const*) const;
    virtual HBLObjectRef        TipName                             (HBLObjectRef,HBLObjectRef cache);
    //HBLObjectRef                TreeBranchName                      (HBLObjectRef node_ref, bool get_subtree = false, HBLObjectRef mapping_mode = nil);
    HBLObjectRef                TreeBranchName                      (HBLObjectRef node_ref, bool get_subtree, HBLObjectRef mapping_mode, HBLObjectRef cache);
    virtual HBLObjectRef        BranchLength                        (HBLObjectRef,HBLObjectRef cache);
    virtual HBLObjectRef        RerootTree                          (HBLObjectRef, HBLObjectRef cache);
    _List*                      SplitTreeIntoClusters               (unsigned long, unsigned long) const;
    _String  const              DetermineBranchLengthMappingMode    (_String const*, hyTopologyBranchLengthMode&) const;
    _AssociativeList*           SplitsIdentity                      (HBLObjectRef, HBLObjectRef cache) const;
    
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
    bool            ConvertToPSW                        (_AVLListX&,_List*, _SimpleList&,bool = false) const;
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

    _String*        ConvertFromPSW                      (_AVLListX&,_SimpleList&) const;
    /* 20090612: SLKP
            given a PSW tree traversal order and a labeling legend,
            return the Newick string for the tree
    */

    void            ComputeClusterTable                 (_SimpleList&, _SimpleList&) const;
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
    
    
    /**
        return the default prefix for generating internal node names (e.g. 'Node')
    */

    const _String     CompareSubTrees                 (_TreeTopology const*, node<long>*) const;
    const _String     FindMaxCommonSubTree            (_TreeTopology const*, long&, _List*) const;




};



#endif
