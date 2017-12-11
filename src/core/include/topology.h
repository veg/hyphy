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








#define kGetNodeStringForTreeName   0x01
#define kGetNodeStringForTreeModel  0x02

class _TheTree;

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

    
    _List*          MapNodesToModels                    (void);

    virtual _String const  GetNodeName                         (node<long> *, bool = false) const;
    virtual const _String*        GetNodeModel                        (node<long> *) const;
    virtual void            GetBranchLength                     (node<long> *, _String&, bool = false) const;
    // SLKP 20100901:
    //               added a boolean flag to ask to return branch length expression (if true) (returns "" for topologies)
    //               just the numeric value (if false)


    virtual hyFloat      GetBranchLength                     (node<long> *) const;
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
    void            SubTreeString                       (node<long>* root, _StringBuffer &, bool = false, long = -1, _AVLListXL* = nil) const;

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
        given a tree argument (p), the function returns an AVL with a 2x1 matrix (key "CLUSTERS")
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



#endif
