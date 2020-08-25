/*

 HyPhy - Hypothesis Testing Using Phylogenies.

 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
 Art FY Poon    (apoon42@uwo.ca)
 Steven Weaver (sweaver@temple.edu)

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

#include "function_templates.h"

#include "topology.h"
#include "tree.h"
#include "vector.h"

#include "global_things.h"
#include "hbl_env.h"

#include <ctype.h>

using namespace hy_global;
using namespace hy_env;

const _String _TreeTopology::kCompareEqualWithReroot          = "Equal with reroot at ",
              _TreeTopology::kCompareEqualWithoutReroot       = "Equal without rerooting",
              _TreeTopology::kCompareUnequalToplogies         = "Unequal topologies (matching label sets).",
              _TreeTopology::kCompareUnequalLabelSets         = "Unequal label sets.",
              _TreeTopology::kMeta                            =  "__meta";

//_______________________________________________________________________________________________

_TreeTopology::_TreeTopology () {
    rooted = UNROOTED;
}

//_______________________________________________________________________________________________
_TreeTopology::_TreeTopology (_TheTree const *top):_CalcNode (*top->GetName(), kEmptyString) {
    _TreeTopology::PreTreeConstructor   (false);
    if (top->theRoot) {
        isDefiningATree         = kTreeIsBeingParsed;
        theRoot                 = top->theRoot->duplicate_tree ();
        
        _TreeTopologyParseSettings parse_settings = CollectParseSettings ();
        
        ConditionalTraverser (
              [&] (node<long>* iterator, node_iterator<long> const& ni) -> bool {
                  _String              nodeVS   (top->GetBranchValue (iterator)),
                                       *nodeSpec = map_node_to_calcnode(iterator)->GetBranchSpec();
                  
                  FinalizeNode (iterator, 0, top->GetNodeName    (iterator), *nodeSpec, nodeVS, parse_settings);
                  DeleteObject (nodeSpec);
                  return false;
              },
        false);
        isDefiningATree         = kTreeNotBeingDefined;
    } else {
        HandleApplicationError ("Can't create an empty tree");
        return;
    }
}

//_______________________________________________________________________________________________
_TreeTopology::_TreeTopology (_TreeTopology const &top) {
    //PreTreeConstructor   (false);
    *this = top;
}

//_______________________________________________________________________________________________
_TreeTopology const &  _TreeTopology::operator = (_TreeTopology const &top) {
    if (top.theRoot) {
        theRoot                 = top.theRoot->duplicate_tree ();
        flatTree.Duplicate(&top.flatTree),
        flatCLeaves.Duplicate(&top.flatCLeaves);
        rooted = top.rooted;
    } else {
        HandleApplicationError ("Can't create an empty tree");
    }
    return *this;
}


//_______________________________________________________________________________________________
_TreeTopology::_TreeTopology    (_String const & name, _String const & parms, bool dupMe, _AssociativeList* mapping):_CalcNode (name,kEmptyString)
// builds a tree from a string
{
    _TreeTopology::PreTreeConstructor   (dupMe);
    _TreeTopologyParseSettings parse_settings = CollectParseSettings();
  
    if (_AssociativeList * meta = MainTreeConstructor  (parms, parse_settings, false, mapping)) {
        _TreeTopology::PostTreeConstructor  (dupMe, meta);
    } else {
        DeleteObject     (compExp);
        compExp = nil;
    }
}

//_______________________________________________________________________________________________
_TreeTopology::_TreeTopology    (_String const * name):_CalcNode (*name,kEmptyString) {}

//_______________________________________________________________________________________________

_TreeTopology::~_TreeTopology (void) {
    if (theRoot) {
        theRoot->delete_tree();
        delete (theRoot);
        theRoot = nil;
    }
    if (compExp) {
        DeleteAndZeroObject(compExp);
    }
}

//_______________________________________________________________________________________________

void    _TreeTopology::PreTreeConstructor (bool) {
    rooted                  = UNROOTED;
    compExp                 = new _Vector;
}

//_______________________________________________________________________________________________

void    _TreeTopology::PostTreeConstructor (bool make_copy, _AssociativeList * meta) {
    
    auto variable_handler = [&] (void) -> void {
        /** TODO SLKP 20171211, make sure the semantics are unchanged */
        variablePtrs.Replace (get_index(), make_copy ? this->makeDynamic() : this, false);
        setParameter (WrapInNamespace (_TreeTopology::kMeta, GetName()), meta ? meta : new _MathObject, nil, false);

      
        /*
         BaseRef temp =  variablePtrs(theIndex);
         variablePtrs[theIndex]=make_copy ? this->makeDynamic() : this;
         DeleteObject(temp);
         */
    };
    
    bool accept_rooted = EnvVariableTrue(accept_rooted_trees);
    
    if (theRoot->get_num_nodes() <= 2) { // rooted tree - check
        if (accept_rooted == false) {
            
            long node_index = theRoot->get_data();
            
            bool recurse = false;
            
            if (theRoot->get_num_nodes() == 2) {
                for (int i = 1; i<=2; i++) {
                    node<long> *node_temp = theRoot->go_down(i);
                    if (node_temp->get_num_nodes()) { // an internal node - make it a root
                        node_temp->detach_parent();
                        node_temp->add_node(*theRoot->go_down(3-i));
                        delete theRoot;
                        theRoot = node_temp;
                        rooted = i == 1 ? ROOTED_LEFT : ROOTED_RIGHT;
                        ReportWarning (_String("Rooted topology. Removing one branch - the ") & (i==1 ? "left" : "right") & " root child has been promoted to be the new root");
                        break;
                    }
                }
                
                if (rooted==UNROOTED) {
                    ReportWarning ("One branch topology supplied - hopefully this IS what you meant to do.");
                    node<long> *node_temp = theRoot->go_down(1);
                    node_temp->detach_parent();
                    node_temp->add_node(*theRoot->go_down(2));
                    delete theRoot;
                    theRoot = node_temp;
                    rooted = ROOTED_LEFT;
                }
            } else {
                if (theRoot->get_num_nodes() == 0) {
                    ReportWarning ("An empty topology has been constructed");
                    variable_handler ();
                    return;
                }
                node<long> *node_temp = theRoot->go_down(1);
                node_temp->detach_parent();
                delete theRoot;
                theRoot = node_temp;
                ReportWarning ("The root has a single child, which is be promoted to the root");
                recurse = false;
            }
            
            flatTree.Delete (node_index);
            flatCLeaves.Delete (node_index);
            ((_Vector*)compExp)->Delete (node_index);
            
            node_iterator<long>  tree_iterator (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
            while (node<long>*topTraverser = tree_iterator.Next()) {
                if (topTraverser->get_data () > node_index) {
                    topTraverser->init (topTraverser->get_data () - 1);
                }
            }
          
            if (recurse) {
              PostTreeConstructor (make_copy, meta);
              return;
            }
        }
    }
    variable_handler ();
    
}

//_______________________________________________________________________________________________

struct _any_char_in_set {
    
    bool valid [256];
    _any_char_in_set (char * accept) {
        InitializeArray(valid, 256, false);
        long c = strlen (accept);
        for (unsigned long i = 0UL; i < c; i++) {
            valid [accept[i]] = true;
        }
    }
    
    
    void invert () {
        ArrayForEach(valid, 256, [] (bool v, unsigned long) -> bool {return !v;});
    }
    
    bool operator == (char c) {
        return valid [(unsigned char)c];
    }
    
};

//_______________________________________________________________________________________________

const _TreeTopologyParseSettings    _TreeTopology::CollectParseSettings (void) {
    _String const kInternalNodePrefix("INTERNAL_NODE_PREFIX"),
                  kIgnoreUserINames  ("IGNORE_INTERNAL_NODE_LABELS");

    _TreeTopologyParseSettings parse_settings;
    
    parse_settings.auto_convert_lengths = EnvVariableTrue(automatically_convert_branch_lengths);
    parse_settings.accept_user_lengths  = EnvVariableTrue(accept_branch_lengths);
    parse_settings.ingore_user_inode_names  = EnvVariableTrue(kIgnoreUserINames);
    
    HBLObjectRef user_node_name            = EnvVariableGet(kInternalNodePrefix, STRING);
    if (user_node_name) {
        parse_settings.inode_prefix = ((_FString*)user_node_name)->get_str();
    }
    
    return parse_settings;

}

//_______________________________________________________________________________________________
_AssociativeList*    _TreeTopology::MainTreeConstructor  (_String const& parms, _TreeTopologyParseSettings & parse_settings, bool checkNames, _AssociativeList* mapping) {
    
    /** TODO SLKP 20171211 this parser needs to be checked
     It admits invalid strings like
     ( a (1,2) , d (4,5),   6)
     strings like ( (,), (,))
     are parsed incorrectly
     */
    
    
    long i,
    nodeCount=0,
    lastNode;
    
    static _String       kBootstrap ("bootstrap"),
                         kComment   ("comment");
    
    static long          kBLLookahead = 128;
  
    _SimpleList nodeStack,
    nodeNumbers;
    
    _String     nodeParameters,
                nodeValue,
                nodeComment,
                nodeBootstrap;
    
    _StringBuffer      nodeName;
    _AssociativeList * node_comments = new _AssociativeList;
    
    _List              local_var_manager;
    local_var_manager < node_comments;
    
    _any_char_in_set non_space        (" \t\r\n"),
                     newick_delimiter (" \t\r\n(),{}[];:");
    non_space.invert();
    newick_delimiter.valid [0] = true;
    
    
    char        lastChar    = '\0';
    bool        isInLiteral = false;
    
    node<long>* currentNode = theRoot = nil,
    * newNode     = nil,
    * parentNode  = nil;
    
    isDefiningATree         = kTreeIsBeingParsed;
    
    
    try {
        auto       push_new_node = [&] (bool make_new) -> void {
            currentNode = newNode;
            nodeStack << (long) currentNode;
            nodeNumbers << nodeCount;
            nodeCount ++;
            if (make_new) {
                newNode = new node <long>;
            }
        };
        
        auto pop_node = [&] () -> void {
            nodeStack.Pop();
            nodeNumbers.Pop();
            //nodeName = kEmptyString;
            nodeName.Clear();
        };
        
        for (i=0; i<parms.length(); i++) {
            
            char   look_at_me = parms.char_at(i);
            
            if (isspace (look_at_me)) {
                continue;
            }
            
            switch (look_at_me) {
                case '(': { // creating a new internal node one level down
                            // a new node
                  
                    newNode = new node<long>;
                    
                    if (lastChar == '(' || lastChar == ',') {
                        if (!currentNode) {
                            delete newNode;
                            throw _String ("Unexpected '") & lastChar & "'";
                        }
                        currentNode->add_node (*newNode);
                    } else {
                        if (theRoot) {
                            parentNode = currentNode->get_parent();
                            if (!parentNode) {
                                throw _String ("Unexpected '") & lastChar & "'";
                            } else {
                                parentNode->add_node (*newNode);
                            }
                        } else {
                            theRoot = newNode;
                        }
                        push_new_node (true);
                        currentNode->add_node (*newNode);
                    }
                    push_new_node (false);
                    break;
                }
                    
                case ',':
                case ')': { // creating a new node on the same level and finishes updating the list of parameters
                    if (nodeStack.empty()) {
                        throw _String ("Unexpected '") & look_at_me & "'";
                    }
                    lastNode = nodeStack.countitems ()-1;
                    parentNode = (node<long>*)nodeStack.get (lastNode);
                    if (mapping) {
                        _FString * mapped_name = (_FString*)mapping->GetByKey (nodeName, STRING);
                        if (!mapped_name) {
                            mapped_name = (_FString*)mapping->GetByKey (nodeName.ChangeCase(kStringUpperCase), STRING);
                        }
                        if (mapped_name) {
                            //printf ("%s => %s\n", nodeName.get_str(), mapped_name->get_str().get_str());
                            nodeName = _String (mapped_name->get_str());
                        }
                    }
                    
                    // handle bootstrap support for this node
                    
                    if (parentNode->get_num_nodes()) {
                        // internal node, check to see if the node name is a bootstrap value
                        _SimpleList numerical_match (nodeName.RegExpMatch(hy_float_regex, 0));
                        if (numerical_match.nonempty() && numerical_match(0) == 0) {
                            nodeBootstrap = nodeName.Cut (0, numerical_match(1));
                            nodeName.Clear();
                        }
                    }
                    
                    nodeName = FinalizeNode (parentNode, nodeNumbers.get (lastNode), nodeName, nodeParameters, nodeValue, parse_settings);
                    
                    if (nodeComment.nonempty() || nodeBootstrap.nonempty()) {
                        _AssociativeList * node_info = new _AssociativeList;
                        if (nodeComment.nonempty()) {
                            if (node_info->MStore(&kComment, new _FString (nodeComment), false)) {
                                kComment.AddAReference();
                                // SLKP 20190507: this is not thread safe
                            }
                        }
                        if (nodeBootstrap.nonempty()) {
                            if (node_info->MStore(&kBootstrap, new _Constant (nodeBootstrap.to_float()), false)) {
                                kBootstrap.AddAReference();
                                // SLKP 20190507: this is not thread safe
                            }
                        }
                        _String * dict_key = new _String (&nodeName, false);
                        if (!node_comments->MStore (dict_key, node_info, false)) {
                            delete dict_key;
                        }
                    }
                    
                    nodeParameters.Clear();
                    nodeComment.Clear();
                    nodeBootstrap.Clear();
                    pop_node ();
                    
                    if (look_at_me == ',') { // also create a new node on the same level
                        newNode = new node<long>;
                        if (!(parentNode = parentNode->get_parent())) {
                            throw _String ("Unexpected '") & look_at_me & "'";
                        }
                        push_new_node (false);
                        parentNode->add_node(*currentNode);
                    }
                    break;
                }
                    
                case '{' : { // parameter list definition
                    long   closing_brace = parms.FindTerminator(i+1, '}');
                    
                    if (closing_brace<0) {
                        throw _String ("'{' has no matching '}'");
                    }
                    nodeParameters = parms.Cut (i+1,closing_brace-1);
                    i = closing_brace;
                    break;
                }
                    
                case '[' : { // hackish Newick annotation
                    long   closing_bracket = parms.FindTerminator(i+1, ']');
                    if (closing_bracket<0) {
                        throw _String ("'[' has no matching ']'");
                    }
                    nodeComment = parms.Cut (i+1,closing_bracket-1);
                    i = closing_bracket;
                    break;
                }
                    
                case ':' : { // tree branch definition
                    
                    // TODO SLKP 20171211: confirm that this will handle expressions with spaces like '..):  0.1'
                    /*_SimpleList numerical_match (parms.RegExpMatch(hy_float_regex, i+1));
                    if (numerical_match.empty()) {
                        throw _String("Failed to read a number for the branch length following ':'");
                    }
                    nodeValue = parms.Cut (numerical_match(0), numerical_match(1));
                    i = numerical_match(1);*/
                    
                    /** SLKP 20200813
                        sscanf is going to call strlen every time, which on long strings could be QUITE expensive!
                        to limit the cost of this, we are going to assume that the float, if it's there is present within kBLLookahead characters of the current position, copy out the next kBLLookahead chars, and read from there.
                    */
                    
                    int end_at;
                    hyFloat length = 0.;

                    if (i + kBLLookahead < parms.length()) {
                        char buffer [kBLLookahead+1];
                        memcpy (buffer, parms.get_str() + i + 1, kBLLookahead);
                        buffer [kBLLookahead] = 0;
                        if (sscanf (buffer, "%lf%n", &length, &end_at) != 1) {
                            throw _String("Failed to read a number for the branch length following ':'");
                        }
                    } else {
                        if (sscanf (parms.get_str()+i+1, "%lf%n", &length, &end_at) != 1) {
                            throw _String("Failed to read a number for the branch length following ':'");
                        }
                    }
                    nodeValue = parms.Cut (i+1, i+end_at);
                    i += end_at;
                    break;
                }
                    
                case ';' : { // done with tree definition
                    break;
                }
                    
                default: { // node name
                    
                    if (nodeName.nonempty()) {
                        throw _String ("Unexpected node name");
                    }
                    
                    long end_of_id = parms.ExtractEnclosedExpression(i, non_space, newick_delimiter, fExtractRespectQuote | fExtractRespectEscape | fExtractOneLevelOnly);
                    if (end_of_id == kNotFound) {
                        throw _String ("Unexpected end of tree string while searching for a node name");
                    } else {
                        nodeName = parms.Cut (i, end_of_id-1);
                        nodeName.StripQuotes("'\"","'\"");
                    }
                    
                    i = end_of_id - 1;
                    break;
                }
            }
          
            lastChar = look_at_me;
        }
    } catch (const _String& error) {
        isDefiningATree = kTreeNotBeingDefined;
        HandleApplicationError (   error & ", in the following string context " &
                                parms.Cut(i>31L?i-32L:0L,i)&
                                " <ERROR HERE> "&
                                parms.Cut(i+1L,parms.length()-i>32L?i+32L:-1L)
                                );
        
        return nil;
    }
    
    if (!theRoot) {
        isDefiningATree = kTreeNotBeingDefined;
        HandleApplicationError ("Can't create empty trees.");
        return nil;
    }
    
    if (nodeStack.countitems() > 1) {
        HandleApplicationError ("Unbalanced () in tree string");
    }
    if (nodeStack.countitems() == 1) {
        parentNode = (node<long>*)nodeStack(0);
        if (mapping) {
            _FString * mapped_name = (_FString*)mapping->GetByKey (nodeName, STRING);
            if (mapped_name) {
                nodeName = _String(mapped_name->get_str());
            }
        }
        FinalizeNode (parentNode, nodeNumbers.get (lastNode), nodeName, nodeParameters, nodeValue, parse_settings);
    }
    
    isDefiningATree = kTreeNotBeingDefined;
    
    //printf ("%s\n", _String ((_String*)node_comments->toStr()).get_str());
    node_comments->AddAReference();
    return node_comments;
    
}

//_______________________________________________________________________________________________

_String    _TreeTopology::FinalizeNode (node<long>* nodie, long number , _String nodeName, _String const& nodeParameters, _String& nodeValue, _TreeTopologyParseSettings const& settings) {
    
    if (nodeName.empty() || (settings.ingore_user_inode_names && nodie->get_num_nodes()>0)) {
        nodeName = settings.inode_prefix & number;
    }

    _String node_parameters = nodeParameters;

    if (nodie==theRoot) {
        node_parameters = kEmptyString;
        nodeValue = kEmptyString;
    }

    nodie->in_object = flatTree.countitems();
    flatTree          && & nodeName;
    flatCLeaves.AppendNewInstance (new _String (&node_parameters, false));

    ((_Vector*)compExp)->Store (ProcessTreeBranchLength(nodeValue));

    node_parameters.Clear();// = kEmptyString;
    nodeValue.Clear(); //      = kEmptyString;

    return nodeName;
}

//_______________________________________________________________________________________________

node<long>* _TreeTopology::FindNodeByName (_String const* match) const {
    node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
    while (node<long>* iterator = ni.Next()) {
        if (*match == GetNodeName  (iterator)) {
            return iterator;
        }
    }
    return nil;

}

//_______________________________________________________________________________________________

void    _TreeTopology::_RemoveNodeList (_SimpleList const& clean_indices) {
    flatTree.DeleteList(clean_indices);
    flatCLeaves.DeleteList(clean_indices);
    
    
    node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
    
    while (node<long>* iterator = ni.Next()) {
        if ((iterator->in_object = clean_indices.CorrectForExclusions(iterator->in_object, -1)) < 0L) {
            throw _String ("Internal error: index remap failed");
        }
    }
}

//_______________________________________________________________________________________________

void    _TreeTopology::RemoveANode (HBLObjectRef nodeName) {
    
    // TODO SLKP 20171213; why did the original implementation delete the entire path
    // up to the root?
    
    try {
        _SimpleList clean_indices;

        auto RemoveNodeByName = [this, &clean_indices] (_FString* remove_me) -> void {

            node<long>* remove_this_node = FindNodeByName (&remove_me->get_str()),
            * parent_of_removed_node;
            
            if (!remove_this_node || ( parent_of_removed_node = remove_this_node->get_parent()) == nil) {
                throw _String ("Node ") & remove_me->get_str().Enquote() & " not found in the tree or is the root node";
            }
            
            
            while (parent_of_removed_node) {
                clean_indices << remove_this_node->in_object;
                //  printf ("Removing %s\n", GetNodeName(remove_this_node).get_str());
                parent_of_removed_node->detach_child(remove_this_node->get_child_num());
                
                
                for (int orphans = 1; orphans <= remove_this_node->get_num_nodes(); orphans++) {
                    parent_of_removed_node->add_node(*remove_this_node->go_down(orphans));
                }
                
                delete remove_this_node;
                remove_this_node = parent_of_removed_node;
                parent_of_removed_node = parent_of_removed_node->get_parent();
                
                // SLKP 20171213: only delete nodes up the chain if they have a single child
                if (remove_this_node->get_num_nodes() == 1) {
                    if (parent_of_removed_node == nil) {
                        /* we are promoting the single remaining child of the current root to be the root */
                        theRoot = remove_this_node->go_down (1);
                        theRoot->detach_parent();
                        clean_indices << remove_this_node->in_object;
                        delete remove_this_node;

                        
                        //parent_of_removed_node = remove_this_node->get_parent();
                        break;
                    }
                } else {
                    break;
                }
            }
        };
        
        if (nodeName->ObjectClass () == STRING) {
            RemoveNodeByName ((_FString*)nodeName);
            clean_indices.Sort();
            _RemoveNodeList (clean_indices);

        } else if (nodeName->ObjectClass () == MATRIX) {
            _Matrix * names = (_Matrix* ) nodeName;
            if (names->IsAStringMatrix()) {
                names->ForEach([RemoveNodeByName] (HBLObjectRef n, unsigned long, unsigned long) -> void {if (n) RemoveNodeByName ((_FString*)n);},
                               [names] (unsigned long i) -> HBLObjectRef {return names->GetFormula(i, -1)->Compute();});
            } else {
                throw _String ("Matrix-valued argument was expected to contain strings");
                
            }
        } else {
           throw _String ("An invalid argument (not a string or a string matrix) supplied");
        }
    } catch (const _String& err) {
        HandleApplicationError (err & " in " & __PRETTY_FUNCTION__);
    }

}


//_______________________________________________________________________________________________

void    _TreeTopology::AddANode (HBLObjectRef newNode) {
    
    static const _String                kNewNodeGraftName  ("NAME"),
    kNewNodeGraftWhere      ("WHERE"),
    kNewNodeGraftParent   ("PARENT"),
    kNewNodeGraftLength   ("LENGTH"),
    kNewNodeGraftParentLength ("PARENT_LENGTH");
    
    try {
    
        if (newNode->ObjectClass () == ASSOCIATIVE_LIST) {
            
            _AssociativeList * newNodeSpec = (_AssociativeList*)newNode;
            _FString         * newName     = (_FString*)newNodeSpec->GetByKey (kNewNodeGraftName, STRING),
                             * newLocation = (_FString*)newNodeSpec->GetByKey (kNewNodeGraftWhere, STRING),
                             * newParent   = (_FString*)newNodeSpec->GetByKey (kNewNodeGraftParent, STRING);
            
            _Constant        * branchLengthSelf = (_Constant*)newNodeSpec->GetByKey (kNewNodeGraftLength, NUMBER),
                             * branchLengthParent = (_Constant*)newNodeSpec->GetByKey (kNewNodeGraftParentLength, NUMBER);
            
            
            if (!newLocation) {
                throw _String("Missing/invalid mandatory argument ") & kNewNodeGraftWhere.Enquote();
            }
            
            
            if (! (newName || newParent)) {
                throw _String("Either ") & kNewNodeGraftName.Enquote() & " or " &  kNewNodeGraftParent.Enquote() & " (or both) must be defined";
            }
            
            node<long>* graftAt = FindNodeByName (&newLocation->get_str());
            if (!graftAt || graftAt->get_parent() == nil) {
                throw _String ("Attachment node must be an exiting non-root node ");
            }
            
            node<long>* newp = newParent ? new node<long> : nil,
            * curp = graftAt->get_parent();
            
            if (newp) {
                newp->set_parent  (*curp);
                newp->add_node    (*graftAt);
                curp->replace_node(graftAt,newp);
            }
            
            if (newName && !newName->empty()) {
                node<long>* newt = new node<long>;
                if (newp) {
                    newp->add_node(*newt);
                } else {
                    graftAt->add_node (*newt);
                }
                _String bl (branchLengthSelf ? branchLengthSelf->Value() : kEmptyString);
                FinalizeNode (newt, 0, newName->get_str(),   kEmptyString, bl, CollectParseSettings());
            }
            
            if (newp && ! newParent->empty()) {
                _String bl (branchLengthParent ? branchLengthParent->Value() : kEmptyString);
                FinalizeNode (newp, 0, newParent->get_str(), kEmptyString, bl, CollectParseSettings());
            }
            
        } else {
            throw _String ("An invalid argument (not an associative array) supplied");
        }
    } catch (const _String& err) {
        HandleApplicationError (err & " in " & __PRETTY_FUNCTION__);
    }
    
}
//_______________________________________________________________________________________________

BaseRef  _TreeTopology::makeDynamic (void) const
{
    _TreeTopology* res = new _TreeTopology;
    res->_CalcNode::Duplicate (this);
    res->flatTree.Duplicate (&flatTree);
    res->flatCLeaves.Duplicate (&flatCLeaves);
    if (compExp) {
        res->compExp = (_Matrix*)compExp->makeDynamic();
    } else {
        res->compExp = nil;
    }

    res->rooted = rooted;
    res->theRoot = CopyTreeStructure (theRoot,true);
    return res;
}


//_______________________________________________________________________________________________

node<long>*  _TreeTopology::CopyTreeStructure (node <long>* theNode,  bool) const {
    node<long>* locNode = new node<long>;
    for (long i=1L; i<=theNode->get_num_nodes(); i++) {
        locNode->add_node(*CopyTreeStructure (theNode->go_down(i), false));
    }
    locNode->init (theNode->in_object);
    return locNode;
}

//_______________________________________________________________________________________________

const _String   _TreeTopology::GetNodeName (node<long>* n, bool fullName) const {
    if (fullName) {
        return *GetName()&'.'&*(_String*)(flatTree.GetItem(n->in_object));
    }
    return *(_String*)(flatTree.GetItem (n->in_object));
 }


//_______________________________________________________________________________________________

_String const*    _TreeTopology::GetNodeModel (node<long>* n) const {
  return (_String*)flatCLeaves.GetItem(n->in_object);
}

  //_______________________________________________________________________________________________

_List*     _TreeTopology::MapNodesToModels (void) {
    _List* map = new _List;
    
    node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
    
    while (node<long>* iterator = ni.Next()) {
        if (iterator->is_root()) {
            break;
        }
        _List * node_record = new _List;
        (*node_record) < new _String (GetNodeName(iterator)) < new _String(*GetNodeModel(iterator));
        (*map) < node_record;
    }
    return map;
}

  //_______________________________________________________________________________________________

_String const  _TreeTopology::GetNodeStringForTree  (node<long> * n , int flags) const {

  _String node_desc;

  if (flags & fGetNodeStringForTreeName) {
    node_desc = GetNodeName (n);
    // check to see if the node name has any chars that are a part of the Newick spec, and if so, enquote the string
    _any_char_in_set newick_delimiter (" (),{}[];:'\"");
    if (node_desc.Any ([&newick_delimiter](char c, unsigned long i) -> bool {
          return newick_delimiter==c;
      }
    ) != kNotFound) {
        node_desc = _StringBuffer ().SanitizeAndAppend(node_desc).Enquote('"');
    }
    
  }

  if (flags & fGetNodeStringForTreeModel) {
    _String const *mSpec = GetNodeModel (n);
    if (mSpec->nonempty()) {
        node_desc = node_desc & mSpec->Enquote('{','}');
    }
  }
  return node_desc;
}

  //_______________________________________________________________________________________________

BaseRef     _TreeTopology::toStr (unsigned long) {
    
  static const _String kNoInternalLabels ("NO_INTERNAL_LABELS");
    
  _StringBuffer     * res = new _StringBuffer((unsigned long)128),
  num;

  bool    skip_internal_labels = EnvVariableTrue(kNoInternalLabels),
           include_model_info = EnvVariableTrue(include_model_spec);


  int            leaf_flag      = fGetNodeStringForTreeName,
                 inode_flag     = skip_internal_labels ? 0 : fGetNodeStringForTreeName;

  if (include_model_info) {
    leaf_flag   |= fGetNodeStringForTreeModel;
    inode_flag  |= fGetNodeStringForTreeModel;
  }

  if (IsDegenerate ()) {
    node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
    (*res)<<'(';

    while (node<long>* iterator = ni.Next()) {
      (*res) << GetNodeStringForTree (iterator,  leaf_flag )
             << (iterator->is_root() ? ')' : ',');
    }
  } else {

    long        level       =   0L,
                lastLevel   =   0L;

    node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);

    node<long>*     curNode   =   ni.Next(),
              *     nextNode;

    level                     = ni.Level();

    nextNode = ni.Next();

    while (nextNode) {
      if (level>lastLevel) {
        if (lastLevel) {
          (*res)<<',';
        }
        res->AppendNCopies('(', level-lastLevel);
      } else if (level<lastLevel) {
        res->AppendNCopies(')', lastLevel-level);
      } else {
        (*res)<<',';
      }

      (*res) << GetNodeStringForTree (curNode,  curNode->is_leaf() ? leaf_flag : inode_flag);


      lastLevel = level;
      curNode   = nextNode;
      level     = ni.Level();
      nextNode  = ni.Next();
    }

    res->AppendNCopies(')', lastLevel-level);
  }
  (*res)<<';';
  res->TrimSpace();
  return res;
}

//__________________________________________________________________________________
void _TreeTopology::toFileStr(FILE* f, unsigned long padding) {
    _String * s = (_String*)toStr(padding);
    fprintf (f, "%s", s->get_str());
    DeleteObject(s);
}


//__________________________________________________________________________________

HBLObjectRef _TreeTopology::ExecuteSingleOp (long opCode, _List* arguments, _hyExecutionContext* context,HBLObjectRef cache) {
    const static _String kSplitNodeNames ("SPLIT_NODE_NAMES");
    
    switch (opCode) { // first check operations without arguments
        case HY_OP_CODE_ABS: // Abs
            return FlatRepresentation(cache);
        case HY_OP_CODE_BRANCHCOUNT: //BranchCount
            return BranchCount(cache);
        case HY_OP_CODE_TIPCOUNT: // TipCount
            return TipCount(cache);
        case HY_OP_CODE_TYPE: // Type
            return Type(cache);
    }
    
    _MathObject * arg0 = _extract_argument (arguments, 0UL, false);
    
    try {
        switch (opCode) { // next check operations without arguments or with one argument
            case HY_OP_CODE_ADD: // +
                if (!arg0) {
                    return Sum(cache);
                }
                AddANode (arg0);
                return _returnConstantOrUseCache(0.0, cache);
                
            case HY_OP_CODE_SUB:
                if (!arg0) {
                    return new _MathObject;
                }
                RemoveANode (arg0);
                return _returnConstantOrUseCache(0.0, cache);
        }
        
        if (arg0) {
            switch (opCode) { // operations that require exactly one argument
                case HY_OP_CODE_MUL: // compute the strict consensus between T1 and T2
                    return SplitsIdentity (arg0, cache);
                    
                case HY_OP_CODE_LEQ: { // MatchPattern (<=)
                    
                    if (arg0->ObjectClass()!=TREE && arg0->ObjectClass()!=TOPOLOGY) {
                        throw _String ("Invalid (not a tree/topology) 2nd argument is call to <= (MatchPattern) for trees/topologies.");
                    }
                    _String  res (((_TreeTopology*)arg0)->MatchTreePattern (this));
                    return _returnConstantOrUseCache (!res.BeginsWith("Unequal"), cache);
                }
                case HY_OP_CODE_EQ: // ==
                    return  Compare(arg0,cache);
                case HY_OP_CODE_BRANCHLENGTH: //BranchLength
                    return BranchLength(arg0,cache);
                case HY_OP_CODE_BRANCHNAME: //BranchName
                    return TreeBranchName(arg0,false,nil,cache);
                case HY_OP_CODE_TIPNAME: //BranchName
                    return TipName(arg0,cache);
                case HY_OP_CODE_MIN: // COT (Min)
                    return FindCOT (arg0,cache);
                case HY_OP_CODE_MAX: // Maximum parsimony
                    return MaximumParsimony (arg0,cache);
               case HY_OP_CODE_REROOTTREE: // RerootTree
                    return RerootTree(arg0, cache);
                case HY_OP_CODE_POWER: //^
                    return AVLRepresentation (arg0, cache);
                case HY_OP_CODE_RANDOM:
                    return RandomizeTips (arg0, cache);
                case HY_OP_CODE_IDIV: { // Split ($) - 2nd argument
                    if (arg0->ObjectClass()!=NUMBER) {
                        throw _String ("Invalid (not a number) 2nd argument is call to $ (split)for trees.");
                    }
                    
                    _List ref_manager;
                    
                    _Constant*  cc     = (_Constant*)TipCount(nil);
                    
                    ref_manager < cc;
                    
                    long        size   = cc->Value()/arg0->Value();
                    
                    if  (size<=4 || size * 2 >cc->Value()) {
                        throw _String("Poor choice of the 2nd numeric agrument in to $ (split) for tree. Either the resulting cluster size is too big(>half of the tree), or too small (<4)!");
                    }
                    
                    long        checkSize = 1L,
                    tol       = 0L;
                    
                    while (tol<size-2) {
                        _List*      resL   = SplitTreeIntoClusters (size,tol);
                        
                        ref_manager < resL;
                        
                        checkSize = cc->Value();
                        
                        if (resL->nonempty()) {
                            _Matrix*    mRes   = new _Matrix (resL->lLength, 2, false, true);
                            resL->ForEach(
                                  [&] (BaseRefConst cluster, unsigned long index) -> void {
                                      _List* cluster_list = (_List*)cluster;
                                      long node_count = ((_Constant*)cluster_list->GetItem(1))->Value();
                                      mRes->Store (index,0, node_count);
                                      mRes->Store (index,1, -2L + cluster_list->countitems());
                                  }
                            );

                        
                            
                            if (checkSize == 0UL) {
                                _List sorted_list;
                                resL->ForEach ([&] (BaseRefConst list, unsigned long) -> void {
                                    sorted_list << ((const _List*)list)->GetItem(0L);
                                });
                                sorted_list.Sort();
                                EnvVariableSet (kSplitNodeNames, new _Matrix (sorted_list, false), false);
                                return mRes;
                            }
                            
                            DeleteObject (mRes);
                        }
                        tol ++;
                    }
                    return new _Matrix (1,1,false, true);
                }
            }
            
            _MathObject * arg1 = _extract_argument (arguments, 1UL, false);
            
            switch (opCode) {
                case HY_OP_CODE_MACCESS: // MAccess
                    return TreeBranchName (arg0,true, arg1, cache);
            }
            
            
            if (arg1) {
                switch (opCode) {
                    case HY_OP_CODE_FORMAT: { // Format
                        _StringBuffer  *tStr = new _StringBuffer  (1024UL);
                        SubTreeString (theRoot, *tStr, CollectParseSettings(), arg0->Compute()->Value() > 0.1 , arg1->Compute()->Value() > 0.1 ? kTopologyBranchLengthExpectedSubs : kTopologyBranchLengthNone, -1, nil);
                        tStr->TrimSpace();
                        return _returnStringOrUseCache(tStr,cache);
                    }
                        
                }
            }
            
        }
        
        switch (opCode) {
            case HY_OP_CODE_MUL: // compute the strict consensus between T1 and T2
            case HY_OP_CODE_LEQ: // MatchPattern (<=)
            case HY_OP_CODE_EQ: // ==
            case HY_OP_CODE_BRANCHLENGTH: //BranchLength
            case HY_OP_CODE_BRANCHNAME: //BranchName
            case HY_OP_CODE_MIN: // COT (Min)
            case HY_OP_CODE_REROOTTREE: // RerootTree
            case HY_OP_CODE_POWER: //^
            case HY_OP_CODE_IDIV:  // Split ($) - 2nd argument
            case HY_OP_CODE_MACCESS: // MAccess
            case HY_OP_CODE_FORMAT:  // Format
                WarnWrongNumberOfArguments (this, opCode,context, arguments);
                break;
            default:
                WarnNotDefined (this, opCode,context);
        }
    } catch (const _String& err) {
        context->ReportError(err);
    }
    
    return new _MathObject;
}


  //_______________________________________________________________________________________________

const _List     _TreeTopology::RetrieveNodeNames (bool doTips, bool doInternals, int traversalType) const {
  _List result;
    
  node_iterator<long> ni (theRoot, traversalType);

  while (node<long> * iterator = ni.Next()) {
    bool is_leaf = iterator->is_leaf();
    if (is_leaf && doTips || doInternals && !is_leaf) {
      result < new _String (GetNodeName(iterator));
    }
  }

  return result;
}

//__________________________________________________________________________________

void _TreeTopology::FindCOTHelper (node<long>* aNode, long parentIndex, _Matrix& distances, _Matrix& rootDistances, _Matrix& branchLengths, _List& childLists, _AVLListX& addressToIndexMap2, hyFloat d)
{
    
    /** TODO: SLKP 20171218, the whole FindCOT stack  needs review **/
    
    long          myIndex     = addressToIndexMap2.GetDataByKey((BaseRef)aNode),
                  leafCount   = distances.GetVDim();
 
    
    _SimpleList * childLeaves = (_SimpleList*)childLists(myIndex);

    _Matrix     * lookup = parentIndex>=0?&distances:&rootDistances;

    if (parentIndex < 0L) {
        parentIndex = 0L;
    }

    long ci2 = 0L;

    hyFloat myLength = branchLengths.get_direct(myIndex);

    for (long ci = 0L; ci < leafCount; ci++) {
        if (ci == childLeaves->list_data[ci2]) {
            if (ci2 < childLeaves->lLength - 1) {
                ci2++;
            }
        } else {
            distances.Store (myIndex, ci, (*lookup)(parentIndex,ci) + myLength);
        }
    }


    for (long ci3 = aNode->get_num_nodes(); ci3; ci3--) {
        FindCOTHelper (aNode->go_down (ci3), myIndex, distances, rootDistances, branchLengths, childLists, addressToIndexMap2, myLength);
    }
}

//__________________________________________________________________________________

void _TreeTopology::FindCOTHelper2 (node<long>* aNode, _Matrix& branchSpans, _Matrix& branchLengths, _AVLListX& addressToIndexMap2, node<long>* referrer, hyFloat d)
{
    /** TODO: SLKP 20171218, the whole FindCOT stack  needs review **/
    long          myIndex     = aNode->parent?addressToIndexMap2.GetXtra(addressToIndexMap2.Find((BaseRef)aNode)):-1;
    hyFloat    myLength    = myIndex>=0?branchLengths.theData [myIndex]:0.0;

    for (long ci = aNode->get_num_nodes(); ci; ci--) {
        node <long>* daChild = aNode->go_down (ci);
        if (daChild != referrer) {
            FindCOTHelper2 (daChild,  branchSpans, branchLengths, addressToIndexMap2, nil, d+myLength);
        }
    }
    if (aNode->parent) { // not the root
        if (d>=0.0) {
            branchSpans.Store (myIndex,0,d);
        } else {
            branchSpans.Store (myIndex,0,0.);
        }

        branchSpans.Store (myIndex,1,d+myLength);
        if (referrer) {
            FindCOTHelper2 (aNode->parent, branchSpans, branchLengths, addressToIndexMap2, aNode, d+myLength);
        }
    }
}

//__________________________________________________________________________________

HBLObjectRef _TreeTopology::MaximumParsimony (HBLObjectRef parameters,HBLObjectRef cache) {
    static const _String
        kMPScore                ("score"),
        kMPLabels               ("labels"),
        kMPOutNodeScore         ("node-scores"),
        kMPOutSubstitutions     ("substitutions");
  
    try {
        CheckArgumentType(parameters, ASSOCIATIVE_LIST, true);
        _AssociativeList * arguments = (_AssociativeList *)parameters;
        
        _AssociativeList * labels    = (_AssociativeList *)arguments->GetByKeyException(kMPLabels, ASSOCIATIVE_LIST);
                         //* scores    = (_AssociativeList *)arguments->GetByKey(kMPScore, ASSOCIATIVE_LIST);
        
        _List           id2name; // integer label -> string label
        
        _Trie           unique_labels, // label -> unique integer ID
                        mapped_labels; // node_name -> integer label
        
        _Matrix         * score_matrix = nil;
        
        auto            score_pair = [&] (unsigned long from, unsigned long to) -> hyFloat {
            if (score_matrix) {
                return (*score_matrix)(from,to);
            }
            return from == to ? 0.0 : 1.0;
        };
        
        for (AVLListXLIteratorKeyValue key_value : labels->ListIterator()) {
            HBLObjectRef node_label = (HBLObjectRef)key_value.get_object();
            _String const * node_name = key_value.get_key();
            if (node_label->ObjectClass() != STRING) {
                throw _String ("Not a string-valued label for node ") & node_name->Enquote();
            }
            _String const * label_value = &((_FString*)node_label)->get_str();
            bool did_insert = false;
            long label_index  = unique_labels.GetValue (unique_labels.InsertExtended (*label_value, unique_labels.countitems(), false, &did_insert));
            if (did_insert) {
                id2name < new _String (*label_value);
            }
            mapped_labels.Insert (*node_name, label_index);
        }
        
        unsigned long unique_label_count = unique_labels.countitems();
        
        _Vector     conditional_scores;
        _SimpleList conditional_labels;
        
        /*
            for each node, in post-order traversal, these lists store K (== unique_label_count)
            in a serialized vector (for scores), and list (for integer labels)
        */

        _SimpleList has_labels;
        /*
            for each node, in post-order traversal
            stores whether or not the subtree starting at this node can be fully labeled
         (-1 : not labeled)
            only such trees are used for labeling
        */
        
        _SimpleList store_values;
        /*
            store node<long> values, because they will be overwritten
            temporarily during MP inference
        */

        node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
        
        // perform the first pass to
        
        long post_order_index = 0L;
        
        while (node<long>* iterator = ni.Next()) {
            
            if (iterator->is_leaf()) {
               
                long index = mapped_labels.FindKey (GetNodeName (iterator));
                if (index >= 0) {
                    index = mapped_labels.GetValue(index);
                    for (unsigned long parent_state = 0UL; parent_state < unique_label_count; parent_state++) {
                        conditional_scores << score_pair (parent_state, index);
                    }
                } else {
                    conditional_scores.AppendRange(unique_label_count, 0., 0.);
                }
                has_labels << index;
                conditional_labels.AppendRange(unique_label_count, index, 0);
            } else {
                conditional_labels.AppendRange(unique_label_count, 0, 0);
                conditional_scores.AppendRange(unique_label_count, 0., 0.);
                bool node_has_labels = true;
                for (unsigned long i = 1UL; i <= iterator->get_num_nodes(); i++) {
                    node_has_labels = node_has_labels && (has_labels.get (iterator->go_down(i)->in_object) >= 0);
                }
                has_labels << (node_has_labels ? 1 : -1);
            }
            store_values << iterator->in_object;
            iterator->in_object = post_order_index++;
        }

        ni.Reset(theRoot);
        
        while (node<long>* iterator = ni.Next()) {
            long my_index = iterator->in_object;
            long label = has_labels.get (my_index);
            if (label >= 0) {
                if (!iterator->is_leaf()) {
                    long offset = my_index * unique_label_count;
                    for (unsigned long parent_state = 0UL; parent_state < unique_label_count; parent_state++) {
                        hyFloat best_score = 1.e100;
                        long best_state = -1L;
                        
                        for (unsigned long my_state = 0UL; my_state < unique_label_count; my_state++) {
                            hyFloat running_score = 0.;
                            for (unsigned long i = 1UL; i <= iterator->get_num_nodes(); i++) {
                                long child_postorder_offset = iterator->go_down(i)->in_object * unique_label_count;
                                running_score += conditional_scores.get(0, child_postorder_offset + my_state);
                            }
                            if (StoreIfLess (best_score, running_score + score_pair (parent_state, my_state))) {
                                best_state = my_state;
                            }
                        }
                        conditional_labels [offset+parent_state] = best_state;
                        conditional_scores [offset+parent_state] = best_score;
                    }
                }
            }
        }
        
        node_iterator<long> ni_pre (theRoot, _HY_TREE_TRAVERSAL_PREORDER);
        
        hyFloat            total_score = 0.;
        _AssociativeList * inferred_labels          = new _AssociativeList;
        _AssociativeList * inferred_substitutions   = new _AssociativeList;
        _AssociativeList * node_scores              = new _AssociativeList;

        
        
        while (node<long>* iterator = ni_pre.Next()) {
            long my_index = iterator->in_object;
            long label = has_labels.get (my_index);
          
          
            if (!iterator->is_leaf () && label >= 0) { // labeled node
                long my_label;
                hyFloat best_score = 1.e100;
                if (!iterator->parent || has_labels.get (iterator->parent->in_object) < 0) {
                    // either a root or not a part of the labeled tree
                    
                    long offset = my_index * unique_label_count;
                    for (unsigned long k = 0UL; k < unique_label_count; k++) {
                        if (StoreIfLess (best_score, conditional_scores.get (0,offset + k))) {
                            my_label = k;
                        }
                    }
                    
                    total_score +=best_score;
                
                } else {
                    // read off the label for this node based on the state of the parent
                    my_label     = conditional_labels.get (my_index * unique_label_count + has_labels.get (iterator->parent->in_object));
                    best_score   = conditional_scores.get (my_index * unique_label_count + has_labels.get (iterator->parent->in_object),0);
                    
                }
                has_labels[my_index] = my_label;
                iterator->in_object = store_values.get (iterator->in_object);
              
                inferred_labels->MStore (GetNodeName(iterator), (*(_String*)id2name.GetItem(my_label)));
                iterator->in_object = my_index;
                node_scores->MStore (GetNodeName (iterator), new _Constant (best_score), false);
            }
            
            if (label >= 0 && iterator->parent) {
                long parent_label = has_labels.get (iterator->parent->in_object);
                if (parent_label >= 0) {
                    // not root, has labeled parent
                    long my_label = has_labels.get (my_index);
                    if (my_label != parent_label) {
                        iterator->in_object = store_values.get (iterator->in_object);
                        _AssociativeList * substitution_record = new _AssociativeList;
                        *substitution_record < (_associative_list_key_value){"from", new _FString  ((*(_String*)id2name.GetItem(parent_label)))}
                                             < (_associative_list_key_value){"to", new _FString  ((*(_String*)id2name.GetItem(my_label)))};
                        inferred_substitutions->MStore (GetNodeName(iterator), substitution_record, false);
                        iterator->in_object = my_index;
                    }
                }
            }
          
            
            //iterator->in_object = store_values.get (iterator->in_object);
        }
        
        _AssociativeList * result = new _AssociativeList;
        result->MStore (kMPScore, new _Constant (total_score), false);
        result->MStore (kMPLabels, inferred_labels, false);
        result->MStore (kMPOutSubstitutions, inferred_substitutions, false);
        result->MStore (kMPOutNodeScore, node_scores, false);

        //printf ("%s\n", _String ((_String*)conditional_scores.toStr(0)).get_str());
        
        ni.Reset(theRoot);
        
        while (node<long>* iterator = ni.Next()) {
            iterator->in_object = store_values.get (iterator->in_object);
            // restore node values
        }
        
        return result;
        
    } catch (const _String& err) {
        HandleApplicationError(err);
    }
    return new  _MathObject;
    
    
}
//__________________________________________________________________________________

_AssociativeList* _TreeTopology::FindCOT (HBLObjectRef p, HBLObjectRef cache) {
    // Find the Center of the Tree (COT) location
    // using an L_p metric (L_2 works well)

    /** TODO: SLKP 20171218, the whole FindCOT stack  needs review **/
    
    static const _String
    kCOTNode               ( "COT_NODE"),
    kCOTSplit             ("COT_SPLIT"),
    kCOTBranchLength         ( "COT_BRANCH_LENGTH"),
    kCOTDistance             ("COT_DISTANCE"),
    kCOTCDF                 ("COT_CDF"),
    kCOTSamples           ("COT_SAMPLES"),
    kCOTSampler            ("COT_SAMPLER"),
    kCOTToNode              ("COT_TO_NODE");

    
    hyFloat         power           = p->Compute()->Value(),
                       totalTreeLength = 0.0;

    _AssociativeList * resList = new _AssociativeList;

    if (power<=0.) {
        HandleApplicationError (_String("Invalid power argument in call to COT finder (Min on trees). Must be positive, had :") & power);
        return resList;
    }

    _SimpleList        avlSL,
                       avlSL2,
                       listOfNodes;

    _List              avlSL3,
                       childLists;

    _AVLListX          addressToIndexMap  (&avlSL),         // leaves only
                       addressToIndexMap2 (&avlSL2),        // complete index
                       lengthToIndexMap   (&avlSL3);

    long               leafCount   = 0,
                       branchCount = 0,
                       tIndex      = 0;

    node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);

    while (node<long>* iterator = ni.Next()) {
        if(iterator->is_root()) {
            break;
        } else if (iterator->is_leaf()) {
            addressToIndexMap.Insert ((BaseRef)iterator, leafCount++);
        } else {
            branchCount++;
        }

        addressToIndexMap2.Insert ((BaseRef)iterator, tIndex++);
    }

    // allocate the matrix of path lengths with hardwired (traversal order) indices
    // also allocate a list of sorted lists to store children nodes
    // and a map of (longed) node addresses to post order traversal indices

    _Matrix      distances          (branchCount+leafCount, leafCount, false, true),
                 rootDistances        (1, leafCount,  false, true),
                 branchLengths        (1, branchCount+leafCount,  false, true),
                 branchSpans      (branchCount+leafCount+1,2,false,true);

     // pass 1: fill up the nodes up to the root (i.e. below any internal node)

    ni.Reset (theRoot);
    tIndex            = 0;

    while (node<long>* iterator = ni.Next()) {
      if (iterator->is_root()) {
        break;
      }
      hyFloat          myLength = GetBranchLength     (iterator);
      lengthToIndexMap.Insert (new _String(totalTreeLength), tIndex, false, true);
      totalTreeLength      += myLength;

      branchLengths.Store (0, tIndex++, myLength);
      listOfNodes << (long)iterator;
      _SimpleList         *childIndices = new _SimpleList;

      if (iterator->is_leaf()) {
          (*childIndices) << addressToIndexMap.GetXtra(addressToIndexMap.Find((BaseRef)iterator));
      } else {
          long           myIndex = addressToIndexMap2.GetXtra(addressToIndexMap2.Find((BaseRef)iterator));
          _SimpleList    mappedLeaves (leafCount,0,0);

          for (long ci = iterator->get_num_nodes(); ci; ci--) {
              long          childIndex = addressToIndexMap2.GetXtra(addressToIndexMap2.Find((BaseRef)iterator->go_down (ci)));
              _SimpleList * childLeaves = (_SimpleList*)childLists(childIndex);

              myLength = branchLengths.theData[childIndex];

              for (long ci2 = 0; ci2 < childLeaves->lLength; ci2++) {
                  long ttIndex = childLeaves->list_data[ci2];
                  mappedLeaves.list_data[ttIndex] = 1;
                  distances.Store (myIndex, ttIndex, distances (childIndex, ttIndex) + myLength);
              }
          }

          for (long ci2 = 0; ci2 < leafCount; ci2++) {
                  if (mappedLeaves.list_data[ci2]) {
                      (*childIndices) << ci2;
                  }
          }
      }

        childLists.AppendNewInstance(childIndices);
    }

    // pass 2: fill the root vector

    //nodeName = "COT_DM1";
    //setParameter (nodeName, &distances);

    for (long ci = theRoot->get_num_nodes(); ci; ci--) {
        long          childIndex = addressToIndexMap2.GetXtra(addressToIndexMap2.Find((BaseRef)theRoot->go_down (ci)));
        _SimpleList * childLeaves = (_SimpleList*)childLists(childIndex);
        hyFloat       myLength = branchLengths.theData[childIndex];
        for (long ci2 = 0; ci2 < childLeaves->lLength; ci2++) {
            tIndex = childLeaves->list_data[ci2];
            rootDistances.Store (0, tIndex, distances (childIndex, tIndex) + myLength);
            //printf ("root->%s = %g\n", ((_String*)leafNames(tIndex))->sData, distances (childIndex, tIndex) + myLength);
        }
    }


    // pass 3: fill in the "other site" branch lengths

    for (long ci3 = theRoot->get_num_nodes(); ci3; ci3--) {
        FindCOTHelper (theRoot->go_down (ci3), -1, distances, rootDistances, branchLengths, childLists, addressToIndexMap2, 0);
    }

    //nodeName = "COT_DM2";
    //setParameter (nodeName, &distances);



    // now traverse the tree and look for the maximum

    tIndex                         = 0;         // stores the index of current min

    hyFloat  currentMin         = 1e100,
                currentBranchSplit = 0;


    for (long ci = distances.GetHDim()-1; ci>=0; ci--) {
        hyFloat    T           = branchLengths.theData[ci];
        _SimpleList * childLeaves = (_SimpleList*)childLists(ci);
        long          ci2         = 0;

        if (CheckEqual (power,2.0)) {
            hyFloat    sumbT  = 0.,
                          sumbT2 = 0.,
                          suma   = 0.,
                          suma2  = 0.;



            for (long ci3 = 0; ci3 < leafCount; ci3++) {
                hyFloat tt = distances(ci,ci3);

                /*
                printf ("%s->%s = %g\n", ttt2.sData, ((_String*)leafNames(ci3))->sData, tt);
                */

                if (ci3 == childLeaves->list_data[ci2]) {
                    if (ci2 < childLeaves->lLength-1) {
                        ci2++;
                    }

                    suma   += tt;
                    suma2  += tt*tt;
                } else {
                    sumbT  += tt;
                    sumbT2 += tt*tt;
                }
            }


            hyFloat tt = (sumbT-suma)/leafCount;/*(sumbT-suma)/leafCount*/;
            if (tt < 0.0) {
                tt = 0.;
            } else if (tt > T) {
                tt = T;
            }

            sumbT = tt*tt*leafCount + 2*tt*(suma-sumbT) + suma2 + sumbT2;

            if (sumbT < currentMin) {
                tIndex             = ci;
                currentBranchSplit = tt;
                currentMin         = sumbT;
            }
        } else {
            hyFloat  step        = T>0.0?T*0.0001:0.1,
                        currentT    = 0.;

            while (currentT<T) {
                hyFloat dTT = 0.0;

                ci2 = 0;

                for (long ci3 = 0; ci3 < leafCount; ci3++) {
                    hyFloat tt = distances(ci,ci3);
                    if (ci3 == childLeaves->get (ci2)) {
                        if (ci2 < childLeaves->lLength-1) {
                            ci2++;
                        }

                        dTT   += pow(tt+currentT,power);
                    } else {
                        dTT   += pow(tt-currentT,power);
                    }
                }

                if (dTT < currentMin) {
                    tIndex             = ci;
                    currentBranchSplit = currentT;
                    currentMin         = dTT;
                }
                currentT += step;
            }
        }
    }

    node <long>*    cotBranch = (node<long>*)listOfNodes.list_data[tIndex];
    

    resList->MStore (kCOTNode,  new _FString( GetNodeName     (cotBranch),false),false);
    resList->MStore (kCOTSplit, new _Constant (currentBranchSplit), false);
    resList->MStore (kCOTDistance, new _Constant (currentMin), false);
    resList->MStore (kCOTBranchLength, new _Constant (branchLengths.directIndex(tIndex)), false);

    //  compute the distribution of lengths away from COT

    FindCOTHelper2  (cotBranch, branchSpans, branchLengths, addressToIndexMap2, nil, currentBranchSplit-branchLengths.theData [tIndex]);
    if (cotBranch->parent) {
        hyFloat adjuster = branchLengths.theData [tIndex]-currentBranchSplit;
        branchSpans.Store (branchCount+leafCount,1,adjuster);
        node <long>* cotParent = cotBranch->parent;
        if (cotParent->parent) {
            adjuster -= branchLengths.theData[addressToIndexMap2.GetXtra(addressToIndexMap2.Find((BaseRef)cotParent))];
        }
        FindCOTHelper2  (cotParent, branchSpans, branchLengths, addressToIndexMap2, cotBranch, adjuster);
    }

    _List        timeSplits;
    _AVLListX    timeSplitsAVL (&timeSplits);

    for (long sc=0; sc <= branchCount+leafCount; sc++) {
        char       buffer[256];
        snprintf  (buffer, 255, "%.15f", branchSpans(sc,0));
        _String nodeName = buffer;
        branchSpans.Store(sc,0,nodeName.to_float());
        timeSplitsAVL.Insert (nodeName.makeDynamic(),0,false,true);
        snprintf  (buffer, 255, "%.15f", branchSpans(sc,1));
        nodeName = buffer;
        branchSpans.Store(sc,1,nodeName.to_float());
        timeSplitsAVL.Insert (nodeName.makeDynamic(),0,false,true);
    }

    _Matrix       cotCDFPoints (timeSplitsAVL.countitems(),3,false,true);

    _AssociativeList  * ctl = new _AssociativeList ();

    _SimpleList tcache;

    long        iv,
                k = timeSplitsAVL.Traverser (tcache, iv, timeSplitsAVL.GetRoot());

    for (long pc = 0; k>=0; k = timeSplitsAVL.Traverser (tcache, iv), pc++) {
        timeSplitsAVL.SetXtra (k, pc);
        cotCDFPoints.Store (pc,0,((_String*)(*((_List*)timeSplitsAVL.dataList))(k))->to_float());
    }


    for (long mxc=0; mxc <= branchCount+leafCount; mxc++) {
        hyFloat T0 =  branchSpans(mxc,0),
                   T1 =  branchSpans(mxc,1);

        if (mxc<branchCount+leafCount) {
            _String nodeName = GetNodeName     ((node<long>*)listOfNodes(mxc));
            ctl->MStore (nodeName, new _Constant (T1), false);
        }

        char       buffer[256];
        snprintf  (buffer, 256, "%.15f", T0);
        _String nn = buffer;
        tcache.Clear();
        long       startingPos =  timeSplitsAVL.Find (&nn,tcache),
                   k           = timeSplitsAVL.Next (startingPos,tcache);

        for (long pc = timeSplitsAVL.GetXtra (k); k>=0; k = timeSplitsAVL.Next (k,tcache), pc++) {
            hyFloat ub = ((_String*)(*((_List*)timeSplitsAVL.dataList))(k))->to_float();
            if (ub < T1 || CheckEqual (T1, ub)) {
                cotCDFPoints.Store (pc,1,cotCDFPoints(pc,1)+1);
            } else {
                break;
            }
        }
    }

    for (long mxc=1; mxc<cotCDFPoints.GetHDim() ; mxc++) {
        cotCDFPoints.Store (mxc,1,cotCDFPoints(mxc,1)*(cotCDFPoints(mxc,0)-cotCDFPoints(mxc-1,0))/totalTreeLength);
        cotCDFPoints.Store (mxc,2,cotCDFPoints(mxc,1)+cotCDFPoints(mxc-1,2));
    }
    timeSplitsAVL.Clear(true);


    resList->MStore (kCOTCDF, &cotCDFPoints, true);
    resList->MStore (kCOTToNode, ctl, false);

    //  sample  random branch placement
    if (totalTreeLength > 0.0) {
        hyFloat     sampler = EnvVariableGetNumber(kCOTSamples,0.0);

        tIndex = sampler;
        if (tIndex >= 1) {
            _Matrix sampledDs  (tIndex, 1, false, true);

            for (long its = 0; its < tIndex; its ++) {
                hyFloat tSample = (genrand_real2 () * totalTreeLength);
                _String nn = tSample;
                long branchIndex = 0;
                if (lengthToIndexMap.FindBest (&nn,branchIndex)<=0) {
                    branchIndex --;
                }

                hyFloat    T         = branchLengths.theData[branchIndex];
                _SimpleList * childLeaves = (_SimpleList*)childLists(branchIndex);

                hyFloat dTT = 0.0;
                tSample -= ((_String*)(*(_List*)lengthToIndexMap.dataList)(branchIndex))->to_float();

                long          ci2 = 0;
                for (long ci3 = 0; ci3 < leafCount; ci3++) {
                    hyFloat tt = distances(branchIndex,ci3);
                    if (ci3 == childLeaves->list_data[ci2]) {
                        if (ci2 < childLeaves->lLength-1) {
                            ci2++;
                        }

                        dTT   += pow(tt+tSample,power);
                    } else {
                        dTT   += pow(tt+T-tSample,power);
                    }
                }

                sampledDs.Store (its,0,dTT);
            }
            EnvVariableSet(kCOTSampler, &sampledDs, true);
        }

    }
    return resList;
}

//__________________________________________________________________________________

_FString*    _TreeTopology::Compare (HBLObjectRefConst p, HBLObjectRef cache) const {
// compare tree topologies
    long objClass = p->ObjectClass();

    if (objClass==TREE || objClass==TOPOLOGY) {
        _String cmp = CompareTrees ((_TreeTopology*)p);
        if (cmp.BeginsWith(kCompareEqualWithReroot)) {
            return (_FString*)_returnStringOrUseCache(cmp.Cut(kCompareEqualWithReroot.length() + ((_TreeTopology*)p)->GetName()->length() + 1, cmp.length()-2), cache);
        } else if (cmp == kCompareEqualWithoutReroot) {
            return (_FString*)_returnStringOrUseCache (_String (' '), cache);
        }
    }
    return (_FString*)_returnStringOrUseCache (kEmptyString, cache);
}


//__________________________________________________________________________________

bool     _TreeTopology::internalNodeCompare (node<long>* n1, node<long>* n2, _SimpleList& subTreeMap, _SimpleList* reindexer, bool cangoup, long totalSize, node<long>* n22, _TreeTopology const* tree2, bool isPattern) const {
    // compares whether the nodes create the same subpartition
    
    // TODO SLKP 20171224: this needs review
    
    // count the number of children
    long  nc1 = n1->get_num_nodes(),
          nc2 = n2->get_num_nodes() + (cangoup&&n2->parent) - (n22&&(!n2->parent));
    
  
    //printf ("Comparing nodes %ld %ld\n", nc1, nc2);
  
    if ( nc1 == nc2 || isPattern && nc2<nc1 ) {
        // now see if the descendants in all directions can be matched
        // first prepare the list of subtrees in the 2nd tree
        
        _List           nodeMap;
        
        _SimpleList*    complement = nil;
        _SimpleList     stSizes;
        
        if ((cangoup||n22)&&n2->parent) {
            complement = new _SimpleList ((unsigned long)totalSize, 1L, 0L);
        }
        
        nc2 = n2->get_num_nodes();
        
        long  skipped_child = -1,
        complementCorrection = 0;
        
        
        for (long k=1; k <= nc2; k++) {
            node<long>* kth_child = n2->go_down(k);
            _SimpleList * childLeaves   = (_SimpleList*)kth_child->in_object;
            if (kth_child == n22) {
                if (complement) {
                    if (reindexer) {
                         childLeaves->Each ([&] (long value, unsigned long index) -> void {
                            long lookup = reindexer->get (value);
                            if (lookup >= 0) {
                                (*complement) [lookup] = 0L;
                            } else {
                                complementCorrection++;
                            }
                        });
                   }
                   else {
                        childLeaves->Each ([&] (long value, unsigned long index) -> void {
                            (*complement) [value] = 0L;
                        });
                   }
                }
                
                skipped_child = k-1;
            } else {
                _SimpleList * all_descendants = new _SimpleList ((unsigned long)totalSize, 0L, 0L);
                
                try {
                    childLeaves->Each ([&] (long value, unsigned long) -> void {
                        long lookup = reindexer ? reindexer->get (value) : value;
                        if (lookup >= 0L) {
                            (*all_descendants) [lookup] = 1L;
                            if (complement) {
                                (*complement) [lookup] = 0L;
                            }
                        } else {
                            throw 0;
                        }
                    });
                } catch (int e) {
                    if (complement) delete complement;
                    delete all_descendants;
                    return 0;
                }
                
                
                stSizes << childLeaves->countitems();
                nodeMap < all_descendants;
            }
        }
        
        if (complement) {
            nodeMap < complement;
            stSizes << totalSize;
            
            for (long k4 = 0; k4 < stSizes.lLength-1; k4++) {
                stSizes.list_data[stSizes.lLength-1] -= stSizes.list_data[k4];
            }
            
            if (n22) {
                stSizes.list_data[stSizes.lLength-1] -= ((_SimpleList*)n22->in_object)->lLength-complementCorrection;
            }
        }
        
        // now go through the subtrees of the first tree and match up (if we can)
        // with the 2nd 1
        
        _SimpleList     unmatchedPatterns;
        
        for (long k5=1; k5 <= nc1; k5++) {
            _SimpleList * childNodes = (_SimpleList*) n1->go_down(k5)->in_object;
            long k6;
            for (k6=0; k6 < stSizes.lLength; k6++)
                if (stSizes.list_data[k6] == childNodes->lLength)
                    // potential subtree match
                {
                    _SimpleList* potMap = (_SimpleList*)nodeMap(k6);
                    long k7;
                    for (k7=0; k7<childNodes->lLength; k7++)
                        if (potMap->list_data[childNodes->list_data[k7]] == 0) {
                            break;
                        }
                    
                    if (k7 == childNodes->lLength) {
                        stSizes.list_data[k6] = -1;
                        
                        if (complement && (k6 == stSizes.lLength-1)) {
                            subTreeMap << -1;
                        } else {
                            if (n22&&(k6>=skipped_child)) {
                                subTreeMap << k6+1;
                            } else {
                                subTreeMap << k6;
                            }
                        }
                        break;
                    }
                }
            
            if (k6 == stSizes.lLength) {
                if (isPattern) {
                    unmatchedPatterns << k5;
                    subTreeMap << -2;
                } else {
                    return 0;
                }
            }
        }
        
        if (isPattern&&unmatchedPatterns.lLength) {
            _SimpleList rematchedPatterns;
            for (long k7 = 0; k7 < unmatchedPatterns.lLength; k7++) {
                rematchedPatterns << -1;
            }
            
            for (long k8 = 0; k8 < stSizes.lLength; k8++) {
                if (stSizes.list_data[k8]>0) {
                    _SimpleList* potMap = (_SimpleList*)nodeMap(k8);
                    for (long k9 = 0; k9 < unmatchedPatterns.lLength; k9++) {
                        if (rematchedPatterns.list_data[k9] < 0) {
                            _SimpleList * childNodes = (_SimpleList*) n1->go_down(unmatchedPatterns.list_data[k9])->in_object;
                            long k10 = 0;
                            for (; k10 < childNodes->lLength; k10++)
                                if (potMap->list_data[childNodes->list_data[k10]] == 0) {
                                    break;
                                }
                            if (k10 == childNodes->lLength) {
                                rematchedPatterns.list_data[k9] = k8;
                                if (complement && (k8 == stSizes.lLength-1)) {
                                    subTreeMap.list_data [unmatchedPatterns.list_data[k9]-1] = -(0xfffffff);
                                } else {
                                    if (n22&&(k8>=skipped_child)) {
                                        subTreeMap.list_data [unmatchedPatterns.list_data[k9]-1] = -k8-3;
                                    } else {
                                        subTreeMap.list_data [unmatchedPatterns.list_data[k9]-1] = -k8-2;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    stSizes.list_data[k8] = 0;
                }
            }
            
            if (rematchedPatterns.Find(-1)>=0) {
                return 0;
            }
            
            for (long k11 = 0; k11 < rematchedPatterns.lLength; k11++) {
                stSizes.list_data[rematchedPatterns.list_data[k11]] -= ((_SimpleList*) n1->go_down(unmatchedPatterns.list_data[k11])->in_object)->lLength;
            }
            
            for (long k12 = 0; k12 < stSizes.lLength; k12++)
                if (stSizes.list_data[k12]) {
                    return 0;
                }
        }
        return 1;
    }
    
    return 0;
}


//long   itcCount = 0;
//__________________________________________________________________________________

char     _TreeTopology::internalTreeCompare (node<long>* n1, node<long>* n2, _SimpleList* reindexer, char compMode, long totalSize, node<long>* n22, _TreeTopology const* tree2, bool isPattern) const
// compare tree topologies
// return values mean
{
    if (n1->get_num_nodes() == 0) {
        return 1;
    } else {
        _SimpleList mapper;
        if (internalNodeCompare (n1,n2, mapper, reindexer, compMode, totalSize, n22, tree2, isPattern)) {
            long nc1 = n1->get_num_nodes();

            _SimpleList     furtherMatchedPatterns;
            _List           patternList;

            for (long k1=0; k1 < nc1; k1++) {
                if (mapper.list_data[k1]>=0) {
                    if (internalTreeCompare (n1->go_down(k1+1),n2->go_down(mapper.list_data[k1]+1), reindexer, 0, totalSize, nil, tree2, isPattern)<1) {
                        return -1;
                    }
                } else if (mapper.list_data[k1] == -1) {
                    if (internalTreeCompare (n1->go_down(k1+1),n2->parent, reindexer, 0, totalSize, n2, tree2, isPattern)<1) {
                        return -1;
                    }
                } else {
                    long idx = -mapper.list_data[k1]-2,
                         k;

                    _SimpleList * patched;
                    if ((k=furtherMatchedPatterns.Find(idx))<0) {
                        k = furtherMatchedPatterns.lLength;
                        furtherMatchedPatterns << idx;
                        patternList < new _SimpleList;
                    }

                    patched = (_SimpleList*)patternList(k);
                    (*patched) << k1;
                }
            }
            for (long k2=0; k2 < furtherMatchedPatterns.lLength; k2++) {
                node <long>* dummy = new node<long>;
                dummy->parent = n1->parent;
                _SimpleList * children = (_SimpleList*)patternList (k2),
                              * newLeaves = new _SimpleList;

                dummy->in_object = (long)newLeaves;


                for (long k3 = 0; k3 < children->lLength; k3++) {
                    node<long>* aChild = n1->go_down(children->list_data[k3]+1);
                    dummy->add_node(*aChild);
                    _SimpleList  t;
                    t.Union (*newLeaves, *(_SimpleList*)aChild->in_object);
                    newLeaves->Clear();
                    newLeaves->Duplicate(&t);
                }


                long k4 = furtherMatchedPatterns.list_data[k2];
                char res = 1;

                if (newLeaves->lLength > 1) {
                    if (k4<n2->get_num_nodes()) {
                        res = internalTreeCompare(dummy, n2->go_down(k4+1),  reindexer, 0, totalSize, nil, tree2, isPattern);
                    } else {
                        res = internalTreeCompare (dummy,n2->parent, reindexer, 0, totalSize, n2, tree2, isPattern);
                    }
                }

                DeleteObject (newLeaves);
                delete (dummy);

                if (res<1) {
                    return -1;
                }

            }
        } else {
            return 0;
        }
    }
    return 1;
}


//__________________________________________________________________________________
void _TreeTopology::EdgeCount (long& leaves, long& internals) const {
    leaves    = 0L;
    internals = 0L;

    node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
    while (node<long>* iterator = ni.Next()) {
        if (iterator->is_leaf()) {
            leaves ++;
        } else {
            internals ++;
        }
    }

}


//__________________________________________________________________________________
HBLObjectRef _TreeTopology::TipCount (HBLObjectRef cache) {
    long leaves, ints;
    EdgeCount (leaves, ints);
    return _returnConstantOrUseCache(leaves, cache);
}

//__________________________________________________________________________________
HBLObjectRef _TreeTopology::BranchCount (HBLObjectRef cache) {
    long leaves, ints;
    EdgeCount (leaves, ints);
    return _returnConstantOrUseCache(ints-1., cache);
}

//__________________________________________________________________________________
HBLObjectRef _TreeTopology::FlatRepresentation (HBLObjectRef cache) {
    _SimpleList     flatTree;

    unsigned long      count = 0UL;
    node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);

    while     (node<long>* iterator = ni.Next()) {
        flatTree << iterator->in_object;
        iterator->in_object = count++;
    }


    _Matrix * res = (_Matrix*)_returnMatrixOrUseCache(1,count, _NUMERICAL_TYPE, false, cache);

    ni.Reset (theRoot);
    count = 0UL;

    while     (node<long>* iterator = ni.Next()) {
        if (iterator->is_root()) {
           res->theData[count] = -1;
        } else {
           res->theData[count] = iterator->parent->in_object;
        }

        iterator->in_object = flatTree.list_data[count++];
    }
    return res;
}

//__________________________________________________________________________________
HBLObjectRef _TreeTopology::AVLRepresentation (HBLObjectRef layoutOption, HBLObjectRef cache) {

    if (layoutOption->ObjectClass () == NUMBER) {
        bool               preOrder = layoutOption->Compute()->Value()>0.5;

        _AssociativeList * masterList = new _AssociativeList ();
        //             arrayKey;


        long         rootIndex = 0;

        _SimpleList  nodeList;
        _AVLListX    nodeIndexList (&nodeList);

        node_iterator<long> ni (theRoot, preOrder ? _HY_TREE_TRAVERSAL_PREORDER : _HY_TREE_TRAVERSAL_POSTORDER);

        while     (node<long>* iterator = ni.Next()) {
            nodeIndexList.Insert ((BaseObj*)iterator, nodeIndexList.countitems()+1L);

            if (iterator->is_root()) {
                rootIndex = nodeIndexList.countitems();
            }
         }

        ni.Reset (theRoot);

        while     (node<long>* iterator = ni.Next()) {
            _AssociativeList * nodeList = new _AssociativeList ();
            nodeList->MStore ("Name", new _FString (GetNodeName (iterator)), false);
            nodeList->MStore ("Length", new _Constant (GetBranchLength (iterator)), false);
            nodeList->MStore ("Depth", new _Constant (ni.Level()), false);
            if (! iterator->is_root()) {
                nodeList->MStore ("Parent", new _Constant(nodeIndexList.GetXtra (nodeIndexList.Find((BaseObj*)iterator->parent))), false);
            }

            long nCount = iterator->get_num_nodes();
            if (nCount) {
                _AssociativeList * childList = new _AssociativeList ();
                for (long k = 1; k<=nCount; k++ ) {
                    childList->MStore (_String((long)(k-1)),new _Constant(nodeIndexList.GetXtra (nodeIndexList.Find((BaseObj*)iterator->go_down(k)))) , false);
                }
                nodeList->MStore ("Children", childList, false);
            }
            masterList->MStore (_String((long)nodeIndexList.GetXtra (nodeIndexList.Find((BaseObj*)iterator))), nodeList, false);
        }

        _AssociativeList * headerList = new _AssociativeList ();

        headerList->MStore ("Name", new _FString (*GetName()), false);
        headerList->MStore ("Root", new _Constant(rootIndex), false);
        masterList->MStore ("0", headerList, false);

        return masterList;
    }
    return new _MathObject;
}

//__________________________________________________________________________________
HBLObjectRef _TreeTopology::TipName (HBLObjectRef p, HBLObjectRef cache) {

    if (p&& p->ObjectClass()==NUMBER) {
        long tip_index        = p->Value(),
             count            = -1L;

        if (tip_index < 0L) {
          return new _Matrix (RetrieveNodeNames(true, false, _HY_TREE_TRAVERSAL_POSTORDER));
        }

        node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
        while (node<long>* iterator = ni.Next()) {
          if (iterator->is_leaf()) {
            count++;
            if (count == tip_index) {
              return _returnStringOrUseCache (GetNodeName (iterator), cache);
            }
          }
        }
    }
    return _returnStringOrUseCache(kEmptyString, cache);
}

//__________________________________________________________________________________

bool _recurse_and_reshuffle (node<long>* root, long& from, long &to, long &leaf_index, _SimpleList& leaf_order, _TreeTopology const& T, hyFloat def_rate, _AssociativeList* node_rates) {
    int children = root->get_num_nodes();
    if (children) {
        
        hyFloat shuffle_rate = def_rate;
        
        if (node_rates) {
            try {
                shuffle_rate = node_rates->GetNumberByKey(T.GetNodeName(root));
            } catch (const _String& e) {
            }
        }
        
        long my_start = leaf_index;
        for (int k = 1; k <= root->get_num_nodes(); k ++) {
            long start, end;
            node<long>* kth_child = root->go_down (k);
            if (_recurse_and_reshuffle (kth_child, start, end, leaf_index, leaf_order,T, def_rate, node_rates)) {
                // reshuffled this child, lock the leaves
                for (long element = start; element <= end; element ++) {
                    if (leaf_order.get (element) >= 0) {
                        leaf_order[element] = -leaf_order.get (element) - 1;
                    }
                }
            }
        }
        from = my_start;
        to = leaf_index-1;

        if (to-from >= 2 && genrand_real1() <= shuffle_rate) { // not a cherry and selected
            //printf ("Shuffling at %s with rate %g\n", T.GetNodeName(root).get_str(), shuffle_rate);
            _SimpleList reshuffled_subset;
            for (long element = from; element <= to; element ++) {
                if (leaf_order.get (element) >= 0L) {
                    reshuffled_subset << leaf_order.get (element);
                }
            }
            reshuffled_subset.Permute(1L);
            long hits = 0;
            for (long element = from; element <= to; element ++) {
                if (leaf_order.get (element) >= 0L) {
                    leaf_order[element] = reshuffled_subset.get (hits++);
                }
            }
            return true;
        }
        return false;
    }
    leaf_index ++;
    return false;
}
//__________________________________________________________________________________
HBLObjectRef _TreeTopology::RandomizeTips (HBLObjectRef rate, HBLObjectRef cache) {
    _TreeTopology * reshuffled = nil;
    
    try {
        if (CheckArgumentType (rate, NUMBER | ASSOCIATIVE_LIST, true)) {
            
            hyFloat default_shuffle_rate = 0.1;
            _AssociativeList* node_level_shuffle_rates = nil;
            if (rate->ObjectClass() == NUMBER) {
                default_shuffle_rate = rate->Compute()->Value();
            } else {
                node_level_shuffle_rates = (_AssociativeList*)rate;
                try {
                    default_shuffle_rate = node_level_shuffle_rates->GetNumberByKey("default");
                } catch (const _String& err) { // no default shuffle rate
                    
                }
            }
            
            
            
            
            if (CheckRange (default_shuffle_rate, 0, 1, true)) {
                reshuffled = new _TreeTopology (*this);
                
                long leaves, ints;
                EdgeCount (leaves, ints);
                
                _SimpleList leaf_order;
                
                /* this will store the reshuffled order of leaves
                 leaf_order[i] = the index of the original leaf i (in post-order traversal) in the reshuffled tree
                 a negative value indicates that the particular leaf has been permuted already, and is locked for subsequent permutations
                 */
                
                node_iterator<long> ni (reshuffled->theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
                _Vector node_shuffle_rates;
                
                while (node<long>* iterator = ni.Next()) {
                    if (iterator->is_leaf()) {
                        leaf_order << iterator->in_object;
                    }
                }
                
                //node_shuffle_rates.Trim();
                //printf ("%s\n", _String((_String*)node_shuffle_rates.toStr()).get_str());
                
                long leaf_index = 0;
                long f,t;
                _recurse_and_reshuffle (&reshuffled->GetRoot(),f,t,leaf_index,leaf_order,*this,default_shuffle_rate,node_level_shuffle_rates);
                ni.Reset(reshuffled->theRoot);
                f = 0;
                //node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
                while (node<long>* iterator = ni.Next()) {
                    if (iterator->is_leaf()){
                        iterator->in_object = leaf_order.get (f++);
                        if (iterator->in_object < 0) {
                            iterator->in_object = -iterator->in_object - 1;
                        }
                    }
                }
                
                
                
                return reshuffled;
            }
        }

    } catch (const _String& e) {
        HandleApplicationError(e);
    }
    return new _MathObject;
}


//__________________________________________________________________________________
HBLObjectRef _TreeTopology::BranchLength (HBLObjectRef p, HBLObjectRef cache) {
  hyFloat branch_length = HY_INVALID_RETURN_VALUE;

  if (p) {
    node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
    if (p->ObjectClass()==NUMBER) {
      long res        = p->Value();

      if (res < 0L) {
          // get ALL branch lengths
        _Vector * branch_lengths = new _Vector;

        while (node<long>* iterator = ni.Next()) {
          if (!iterator->is_root()){
            branch_lengths->Store(GetBranchLength (iterator));
          }
        }
        branch_lengths->Store (0.0); // backward compatibility with storing 0 at the root
        branch_lengths->Trim();
        branch_lengths->Transpose();
        return branch_lengths;
      } else {
          // get a specific branch length
        long count = -1L;
        while (node<long>* iterator = ni.Next()) {
          if (!iterator->is_root()) {
            if (++count == res) {
              branch_length = GetBranchLength(iterator);
              break;
            }
          }
        }
      }
    } else {
      if (p->ObjectClass()==STRING) {
        _List twoIDs = ((_FString*)p->Compute())->get_str().Tokenize(";");

        if (twoIDs.countitems() == 2UL || twoIDs.countitems() == 1UL) {

          _String * nodes[2] = {(_String*)twoIDs.GetItem (0L),
                                 twoIDs.countitems()>1UL?(_String*)twoIDs.GetItem (1L):nil};


          node<long>* node_objects[2] = {nil,nil};
          long levels[2] = {0L,0L};


          while (node<long>* iterator = ni.Next()) {
            _String nn = GetNodeName(iterator);
            for (unsigned long i = 0L; i < 2UL; i++) {
              if (nodes[i] && nn == *nodes[i]) {
                node_objects[i] = iterator;
                levels[i] = ni.Level();
                if (node_objects[1-i]) {
                  break;
                }
              }
            }
          }

          if (node_objects[0] && node_objects[1]) {
            branch_length = 0.;

            for (long i = 0L; i < 2L; i++) { // walk up to the same depth
              while (levels[1-i] < levels[i]) {
                branch_length      += GetBranchLength(node_objects[i]);
                node_objects[i] = node_objects[i]->get_parent();
                levels[i]--;
              }
            }

            while (node_objects[0] != node_objects[1]) {
              branch_length += GetBranchLength (node_objects[0]) + GetBranchLength (node_objects[1]);
              node_objects[0] = node_objects[0]->parent;
              node_objects[1] = node_objects[1]->parent;
            }
          } else if (node_objects[0]) {
            if (nodes[1]) {
              if (*nodes[0] == *nodes[1]) {
                branch_length = 0.0;
              } else if (*nodes[1] == kExpectedNumberOfSubstitutions) {
                _String bl (GetBranchLengthString (node_objects[0], true));
                if (bl.nonempty()) {
                  return new _FString (bl);
                }
              }
            } else {
              branch_length = GetBranchLength(node_objects[0]);
            }
          }
        }
      }
    }
  }

  if (isnan (branch_length)) {
    return new _MathObject ();
  }

  return _returnConstantOrUseCache(branch_length, cache);
}

//__________________________________________________________________________________

HBLObjectRef _TreeTopology::TreeBranchName (HBLObjectRef node_ref, bool get_subtree, HBLObjectRef mapping_mode, HBLObjectRef cache) {
  _StringBuffer branch_name;
    

  if (node_ref) {
    if (node_ref->ObjectClass()==NUMBER) { // get by index
 
      long argument      = node_ref->Value(),
           count         = -1L;

      if (argument>=0L) { // get a specific internal node name/subtree
          
        node<long>* ith_internal_node = ConditionalTraverser (
              [&] (node <long>* iterator, node_iterator<long>& ni) -> bool {
                  if (iterator->is_leaf () == false) {
                      return (++count == argument);
                  }
                  return false;
              },
        false);
          
        if (!ith_internal_node) {
          ith_internal_node = theRoot;
        }
          
        if (ith_internal_node) {
            if (get_subtree) {
              hyTopologyBranchLengthMode  mapMode  = kTopologyBranchLengthNone;
              if (mapping_mode) {
                _String * t = (_String*)mapping_mode->Compute()->toStr();
                DetermineBranchLengthMappingMode (t,mapMode);
                DeleteObject (t);
              }
              SubTreeString         (ith_internal_node, branch_name,CollectParseSettings(), true,mapMode);
            } else {
              branch_name = GetNodeName (ith_internal_node);
            }
        }
      } else {
        _List branch_lengths;
         node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
         while (node<long>* iterator = ni.Next()) {
           branch_lengths.AppendNewInstance(new _String (GetNodeName (iterator)));
         }
        return new _Matrix (branch_lengths);
      }
    } else {
      if (node_ref->ObjectClass()==STRING) {
        _List twoIDs = ((_FString*)node_ref->Compute())->get_str().Tokenize(";");

        if (twoIDs.countitems() == 2UL || twoIDs.countitems() == 1UL) {

          _String * nodes[2] = {(_String*) twoIDs.GetItem(0),
                                (_String*)(twoIDs.countitems() > 1L?twoIDs.GetItem(1):nil)};



          if (twoIDs.countitems() == 1UL) {
            _AssociativeList * resList = new _AssociativeList;
            long              masterLevel = 0L;

            node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_PREORDER);
            while (node<long>* iterator = ni.Next()) {
              _String node_name = GetNodeName   (iterator);
              if (node_name == *nodes[0]) {
                masterLevel = ni.Level();
                  //resList->MStore(node_name,new _Constant (iterator->get_num_nodes()));
                do {
                  // 20151203: SLKP this will store the root name twice; commented the line above
                  resList->MStore(GetNodeName   (iterator),new _Constant (1L+(iterator->get_num_nodes()>0L)), false);
                  iterator = ni.Next();
                } while (iterator && ni.Level() > masterLevel);
                break;
              }
            }
            if (resList->countitems() == 0L) {
              // did not find the target node
              DeleteObject (resList);
              return new _MathObject;
            }
            return resList;
          } else {
              // this returns the sequence of nodes between node 1 and node 2
            node<long>* node_objects[2]= {nil, nil};
            long levels[2] = {0L, 0L};

            node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);

            while (node<long>* iterator = ni.Next()) {
              _String nn = GetNodeName(iterator);
              for (long i = 0L; i < 2L; i++) {
                if (nodes[i] && nn == *nodes[i]) {
                  node_objects[i] = iterator;
                  levels[i] = ni.Level();
                  if (node_objects[1-i]) {
                    break;
                  }
                }
              }
            }

            if (node_objects [0] && node_objects [1]) {
              _List partial_paths [2];

              for (long i = 0L; i < 2L; i++) { // walk up to the same depth
                while (levels[1-i] < levels[i]) {
                  partial_paths[i] < new _String (GetNodeName(node_objects[i]));
                  node_objects[i] = node_objects[i]->get_parent();
                  levels[i]--;
                }
              }

              while (node_objects[0] != node_objects[1]) {
                for (long i = 0L; i < 2L; i++) {
                  partial_paths[i] < new _String (GetNodeName(node_objects[i]));
                  node_objects[i] = node_objects[i]->parent;
               }
             }
              partial_paths[1].Flip();
              partial_paths[0] << partial_paths[1];
              return new _Matrix(partial_paths[0]);
            } else if (node_objects[0]) {
              return new _Matrix();
            }
            return new _MathObject();
          }
        }
      }
    }
  }
  //return new _FString (branch_name, false);
  return _returnStringOrUseCache(branch_name, cache);
}

  //__________________________________________________________________________________
void _TreeTopology::SubTreeString (node<long>* root, _StringBuffer &result, _TreeTopologyParseSettings const& settings, bool all_names, hyTopologyBranchLengthMode mode, long branch_length_variable, _AVLListXL* substitutions) const {
    
    
  auto handle_node = [&] (node<long>*  iterator, bool is_degenerate) -> void {
    _String node_name = GetNodeName (iterator);
      if (substitutions) {
            if (BaseRefConst replacement = substitutions->GetDataByKey(&node_name)) {
                node_name = *(_String*)replacement;
            }
        }

        if (is_degenerate || !iterator->is_root()) {
          if (all_names || !node_name.BeginsWith(settings.inode_prefix)) {
            result << node_name;
          }
          PasteBranchLength (iterator,result,mode,branch_length_variable);
        }
  };
                          
  if (root->parent == NULL && IsDegenerate()) {
    node_iterator<long> ni (root, _HY_TREE_TRAVERSAL_POSTORDER);
    result << '(';

    while (node<long>* iterator = ni.Next()) {
        handle_node (iterator, true);
        result << (iterator->is_root() ? ')' : ',');
    }
  } else {
      long    last_level        = 0L;

      node_iterator<long> ni (root, _HY_TREE_TRAVERSAL_POSTORDER);

      while (node <long>* iterator = ni.Next()) {

        if (ni.Level() > last_level) {
          if (last_level) {
            result<<',';
          }
          result.AppendNCopies('(', ni.Level()-last_level);
        } else if (ni.Level() < last_level) {
          result.AppendNCopies(')', last_level-ni.Level());
        } else if (last_level) {
          result<<',';
        }

        last_level = ni.Level();
        handle_node (iterator, false);
      }
  }
}


//__________________________________________________________________________________
void _TreeTopology::RerootTreeInternalTraverser (node<long>* iterator, long originator, bool passedRoot, _StringBuffer &res, _TreeTopologyParseSettings const& settings, hyTopologyBranchLengthMode branch_length_mode, long variable_ref, bool first_time) const {

  //printf ("[RerootTreeInternalTraverser]%s %ld %ld\n", GetNodeName(iterator).sData, originator, passedRoot);

    if (passedRoot) {
        SubTreeString (iterator, res, settings, false, branch_length_mode, variable_ref);
    } else {
        // move to parent now
        node<long>*     iterator_parent = iterator->get_parent();
        
        /*
        StringToConsole(GetNodeName(iterator)); NLToConsole();
        if (iterator_parent) {
            StringToConsole(GetNodeName(iterator_parent)); NLToConsole();
        }
        */

        if (iterator_parent) { // not root yet
            bool is_root_next = iterator_parent->get_parent() == NULL;
            if (!is_root_next) {
                res<<'(';
            }
            long the_index_of_this_child = iterator->get_child_num();
            RerootTreeInternalTraverser (iterator_parent, the_index_of_this_child ,false,res,settings,branch_length_mode,variable_ref,first_time);

            if (iterator_parent->get_parent()) {

              for (long i = 1L; i<=iterator_parent->get_num_nodes(); i++) {
                  if (i!=the_index_of_this_child) {
                    res<<',';
                    SubTreeString (iterator_parent->go_down(i),res, settings, false, branch_length_mode,variable_ref);
                  }
              }
             }
            
            if (!is_root_next) {
                res<<')';
            
                if (!first_time) {
                  _String node_name = GetNodeName (iterator);
                  if (!node_name.BeginsWith(settings.inode_prefix)) {
                    res<<node_name;
                  }
                }
                PasteBranchLength (iterator,res,branch_length_mode, variable_ref);
            }
        } else {
            /* passing old root
               create a new root with >=2 children nodes - this node,
               and one more containing all other children (>=2 of them)
            */
            long count               = 0L,
                 root_children_count = theRoot->get_num_nodes();

            if (root_children_count > 2) {
                res << '(';
            }

            node<long>* stash_originator = nil;

            for (long k = 1; k<=theRoot->get_num_nodes(); k++) {
                if (k==originator) {
                    stash_originator = theRoot->go_down(k);
                    continue;
                }
                if (count) {
                    res<<',';
                }
                count++;
                SubTreeString (theRoot->go_down(k), res,settings, false, branch_length_mode,variable_ref);
            }

            if (!stash_originator) {
              HandleApplicationError ("Internal error in RerootTreeInternalTraverser");
              return;
            }

            if (root_children_count > 2) {
                res<<')';
            }

            PasteBranchLength (stash_originator,res, branch_length_mode, variable_ref);
        }
    }
}


//__________________________________________________________________________________
void            _TreeTopology::PasteBranchLength (node<long>* node, _StringBuffer& result,
                                                  hyTopologyBranchLengthMode const mode, long variable_reference , hyFloat factor) const {
    
    auto decorate_string = [&] (_String const& value) -> void {
        if (value.nonempty()) {
            result << ':' << _String (value.to_float()*factor);
        }
    };
    
    switch (mode) {
        case kTopologyBranchLengthExpectedSubs: {
            result << ':' << _String(GetBranchLength (node)*factor);
            break;
        }
        case kTopologyBranchLengthUserLengths: {
            decorate_string (GetBranchValue (node));
            break;
        }
        case kTopologyBranchLengthLocalParameter: {
            decorate_string (GetBranchVarValue (node, variable_reference));
        }
    }
}

//__________________________________________________________________________________
const _String            _TreeTopology::GetBranchLengthString (node<long> * n, bool get_expression) const {
    if (get_expression) {
        return kEmptyString;
    } else {
        return compExp->theData[n->in_object];
    }
}

//__________________________________________________________________________________
hyFloat            _TreeTopology::GetBranchLength (node<long> * n) const {
    return compExp->theData[n->in_object];
}


//__________________________________________________________________________________
const _String            _TreeTopology::GetBranchValue     (node<long> *) const{
    return kEmptyString;
}
//__________________________________________________________________________________
const _String            _TreeTopology::GetBranchVarValue  (node<long> *, long) const {
    return kEmptyString;
}


//__________________________________________________________________________________
HBLObjectRef _TreeTopology::RerootTree (HBLObjectRef new_root, HBLObjectRef cache) {
    _StringBuffer * res = new _StringBuffer (256UL);

    _TreeTopologyParseSettings settings = CollectParseSettings();

    if (new_root && new_root->ObjectClass()==STRING) {
        if (rooted == UNROOTED) {
            ReportWarning ("Reroot was called with an unrooted tree. Rerooting was still performed.");
        }

        _String * new_root_name = (_String*)new_root->toStr();
        node<long>* reroot_at = FindNodeByName (new_root_name);
        DeleteObject (new_root_name);
        
        if (reroot_at) { // good node name, can reroot
            if (reroot_at->is_root()) {
                SubTreeString (theRoot, *res,settings,kTopologyBranchLengthUserLengths);
            } else {
              (*res)<<'('; // opening tree (
              RerootTreeInternalTraverser (reroot_at, reroot_at->get_child_num(),false,*res,settings,kTopologyBranchLengthUserLengths,-1,true);
              (*res)<<',';
              SubTreeString (reroot_at, *res,settings,kTopologyBranchLengthUserLengths);
              (*res)<<')';
            }
        }
    } else {
        HandleApplicationError ("Reroot Tree was passed an invalid branch argument.");
    }

    res->TrimSpace();
    return _returnStringOrUseCache(res, cache);
}


//__________________________________________________________________________________

_String  const    _TreeTopology::DetermineBranchLengthMappingMode (_String const * param, hyTopologyBranchLengthMode& map_mode) const {
    map_mode = kTopologyBranchLengthNone;
    if (param) {
        if (*param == hy_env::kExpectedNumberOfSubstitutions) {
            map_mode = kTopologyBranchLengthExpectedSubs;
        } else if (*param == hy_env::kStringSuppliedLengths) {
            map_mode = kTopologyBranchLengthUserLengths;
        } else {
            map_mode = kTopologyBranchLengthLocalParameter;
            return _String('.') & *param;
        }
    }
    return kEmptyString;
}

//_______________________________________________________________________________________________

node<long>* _TreeTopology::prepTree4Comparison (_List& leafNames, _SimpleList& mapping, node<long>* topNode) const {
  
    /**
        Creates a tree structure which facilitates comparison to other trees
     
        The in_object of each node is a numeric list of indices of all the leaves that are at or below this node
        The indices are in post-order traversal order
     
        `leafNames` stores lexicographically sorted leaf names
        `mapping` stores the mapping between the original post-order index of a leaf and its location in `leafNames`
     
     */
  
    node<long>* res     = topNode?topNode->duplicate_tree():theRoot->duplicate_tree();
    
    node_iterator<long> ni (res, _HY_TREE_TRAVERSAL_POSTORDER);
    
    _SimpleList     indexer;
    
    while (node<long> * iterator = ni.Next()) {
        
        _SimpleList * descendants = new _SimpleList;
        
        if (!iterator->is_leaf()) {
            for (unsigned long k = 1UL; k <= iterator->get_num_nodes(); k++) {
                (*descendants) << *(_SimpleList*)iterator->go_down(k)->in_object;
            }
        } else {
            (*descendants) << leafNames.countitems();
            indexer        << leafNames.countitems();
            leafNames < new _String (GetNodeName (iterator));
        }
        
        iterator->in_object = (long)descendants;
    }
    
    // now sort leaf names and build the indexer
    
    mapping.Clear();
    mapping.Duplicate (&indexer);
    
    SortLists (&leafNames, &indexer);
    SortLists (&indexer, &mapping);
    
    return res;
}

//_______________________________________________________________________________________________

void    _TreeTopology::destroyCompTree (node<long>* compRoot) const {
    for (unsigned long i=1UL; i<=compRoot->get_num_nodes(); i++) {
        destroyCompTree (compRoot->go_down(i));
    }
    DeleteObject ((BaseRef)compRoot->in_object);
    delete (compRoot);
}

//_______________________________________________________________________________________________

_String             _TreeTopology::CompareTrees      (_TreeTopology* compareTo) const {

    _List           myLeaves,
                    otherLeaves;

    _SimpleList     indexer,
                    otherIndexer;

    node<long>      *myCT,
                    *otherCT;

    _String         rerootAt;

    myCT            = prepTree4Comparison(myLeaves, indexer);
    otherCT         = compareTo->prepTree4Comparison(otherLeaves, otherIndexer);

    // first compare the set of leaf labels

    if (myLeaves.Equal (otherLeaves)) {
        _SimpleList * reindexer = nil;

        if (indexer.Equal (otherIndexer)) {
            otherIndexer.Populate (_SimpleList ((unsigned long)myLeaves.countitems()).Populate (indexer, _SimpleList::action_invert),
                                   _SimpleList::action_compose);
            reindexer = &otherIndexer;
        }

        char compRes;

        if (internalTreeCompare (myCT, otherCT, reindexer, 1, myLeaves.lLength, nil, compareTo)>0) {
            rerootAt = kCompareEqualWithoutReroot;
        } else {
            long   tCount = 0;

            node_iterator<long> ni (otherCT, _HY_TREE_TRAVERSAL_POSTORDER);
            node<long> * meNode;

            while (meNode = ni.Next()) {
              if (meNode == otherCT) {
                break;
              }
              if (meNode->get_num_nodes()) {
                  compRes = internalTreeCompare (myCT, meNode, reindexer, 1, myLeaves.lLength, nil, compareTo);
                  if (compRes>0) {
                      break;
                  } else if (compRes) {
                      meNode = otherCT;
                      break;
                  }
              }

              tCount ++;
            }

            if (meNode!=otherCT) {
                node_iterator<long> ni (compareTo->theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
                while (node<long> * meNode = ni.Next()) {
                    if (meNode == theRoot) {
                      break;
                    }
                    if (tCount==1) {
                        rerootAt = kCompareEqualWithReroot & *LocateVar (meNode->in_object)->GetName() & '.';
                        break;
                    } else {
                        tCount --;
                    }
                }
            }
        }
        if (rerootAt.empty()) {
            rerootAt = kCompareUnequalToplogies;
        }
    } else {
        rerootAt = kCompareUnequalLabelSets;
    }

    destroyCompTree (myCT);
    destroyCompTree (otherCT);

    return          rerootAt;
}




//_______________________________________________________________________________________________

_List*  _TreeTopology::SplitTreeIntoClustersInt (node<long>* , _List* , _AVLListX& , long , long ) const
// TODO SLKP 20180202: clearly this needs work (or needs to be deprecated)
// returns the list of nodes included in the splits BELOW trav1
{
    /*_List * dependancies = nil;

    for (long k = trav1->get_num_nodes(); k; k--)
    {
        _List * me = SplitTreeIntoClustersInt(trav1->go_down(k), res, counts, size, tol);
        if (me)
        {
            if (!dependancies)
                checkPointer (dependancies = new _List);

            (*dependancies) << (*me);
            DeleteObject (me);
        }
    }

    long distinctElements = counts.GetXtra (counts.Find((BaseRef)trav1->in_object)),
         dLength = dependancies?dependancies->lLength:0;

    if (((trav1->parent==nil)||(distinctElements+dLength>=size-tol))&&(distinctElements+dLength<=size+tol))
    {
        _List   entry;
        _String nName;
        GetNodeName (trav1, nName);
        entry && & nName;
        entry && & _Constant (distinctElements);
        if (dependancies)
            entry << (*dependancies);
        (*res) && & entry;

        trav1 = trav1->parent;
        while (trav1)
        {
            long aidx = counts.Find((BaseRef)trav1->in_object);
            counts.SetXtra(aidx,counts.GetXtra (aidx) - distinctElements);
            trav1 = trav1->parent;
        }

        if (dependancies)
            dependancies->Clear();
        else
            checkPointer (dependancies = new _List);

        (*dependancies) && & nName;
    }

    return dependancies;*/

    // decide if child or i-node

    /*if (trav1->get_num_nodes())
    {

    }
    else
    {
        node <long>* pn = trav1->parent;

        while (pn)
        {
            long pidx = counts.Find((BaseRef)pn->in_object);

        }
    }*/
    // TBI
    return new _List;
}

//_______________________________________________________________________________________________

_List*  _TreeTopology::SplitTreeIntoClusters (unsigned long size, unsigned long tol) const {
// TODO SLKP 20180202: clearly this needs work (or needs to be deprecated)
// assume fixed rooting for now
// returns the list of list speccing pointers to calcnodes rooting the splits
// the member list contains 2 or more entries for each node:
// itself, the number of (independent)leaves encompassed by that node,
// and an optional references to the node whose
// results it depends on.

// assumes that size-tol>=2

    _SimpleList     counts;
    _AVLListX       cavl (&counts);

    node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
    while (node<long> * iterator = ni.Next()) {
        long nC = iterator->get_num_nodes();
        if (nC) {
            long c = 0;
            for (long k=1; k<=nC; k++) {
                c += counts.list_data[iterator->go_down(k)->in_object];
            }
            cavl.Insert((BaseRef)iterator->in_object,c);
        } else {
            cavl.Insert((BaseRef)iterator->in_object,1);
        }
    }

    _List*                    result = new _List;
    DeleteObject(SplitTreeIntoClustersInt (theRoot, result, cavl, size, tol));
    return           result;
}

  //_______________________________________________________________________________________________

const _String _TreeTopology::MatchTreePattern (_TreeTopology const* compareTo) const {
    // TODO: 20180202, SLKP, this needs to be confirmed as working
    // the pattern is this
    // compare to is the potential match
  _List           myLeaves,
                  otherLeaves,
                  overlappedLeaves;

  _SimpleList     indexer,
                  otherIndexer;

  node<long>      *myCT,
                  *otherCT;

  _String         rerootAt;

  myCT            = prepTree4Comparison(myLeaves, indexer);
  otherCT         = compareTo->prepTree4Comparison(otherLeaves, otherIndexer);


  _SimpleList      matchedLeaves;
  overlappedLeaves.Intersect (myLeaves, otherLeaves,&matchedLeaves);
    
  unsigned long my_leaf_count = myLeaves.countitems();

  if ((my_leaf_count >=otherLeaves.countitems())&&(overlappedLeaves.countitems() == otherLeaves.countitems())) {
    if (my_leaf_count >otherLeaves.countitems()) {
        // map leaves back to ordering space

        // trim the pattern tree
      _SimpleList     allowedLeaves  (my_leaf_count, 0, 0),
                      recordTransfer (my_leaf_count),
                      inverted (my_leaf_count);
        
      inverted.Populate (indexer, _SimpleList::action_invert);
      
      for (long padresSuck = 0; padresSuck < matchedLeaves.lLength; padresSuck ++) {
        allowedLeaves.list_data[inverted.get(matchedLeaves.get(padresSuck))] = 1;
      }

      long slider = 0L;

      for (long dodgersSuck = 0L; dodgersSuck < my_leaf_count; dodgersSuck ++)
        if (allowedLeaves.get (dodgersSuck)) {
          recordTransfer.list_data[dodgersSuck] = slider++;
        }

        // pass 1 - delete all the superfluous leaves

      node_iterator<long> ni (myCT, _HY_TREE_TRAVERSAL_POSTORDER);
      while (node<long>* iterator = ni.Next()) {
        if (iterator != myCT && !iterator->is_leaf()) {
          _SimpleList*  descendants = ((_SimpleList*)iterator->in_object);

          if ((descendants&&(allowedLeaves.list_data[descendants->list_data[0]] == 0))||(!descendants)) {
              // mark this node for deletion
            if (descendants) {
              DeleteObject(descendants);
            }

            node<long>* sacLamb = iterator;
            ni.Next();

            if (sacLamb->parent->get_num_nodes()==1) {
              DeleteObject((BaseRef)sacLamb->parent->in_object);
              sacLamb->parent->in_object = 0L;
            }

            sacLamb->parent->detach_child (sacLamb->get_child_num());
            delete (sacLamb);
            continue;
          }
        }
      }

        // pass 2 - prune internal nodes with exactly one child

        // >O
      ni.Reset (myCT);

      while (node<long>* iterator = ni.Next()) {
        long nn = iterator->get_num_nodes();
        if (nn == 1) {
          node<long>* loneChild = iterator->go_down(1);
          DeleteObject((_SimpleList*)iterator->in_object);
          iterator->in_object = loneChild->in_object;
          iterator->detach_child (1);
          delete (loneChild);
        } else {
          if (nn > 1) {
            _SimpleList * myDescs = (_SimpleList*)iterator->in_object;
            myDescs->Clear();

            for (long cc = 1; cc <= nn; cc++) {
              _SimpleList temp;

              temp.Union (*myDescs, *(_SimpleList*)iterator->go_down(cc)->in_object);

              myDescs->Clear();
              myDescs->Duplicate (&temp);
            }
          } else {
            ((_SimpleList*)iterator->in_object)->list_data[0] = recordTransfer.list_data[((_SimpleList*)iterator->in_object)->list_data[0]];
          }

        }
      }

      _List   newLeaves;
      recordTransfer.Clear();
      inverted.Clear();

      for (long lc = 0; lc < allowedLeaves.lLength; lc++)
        if (allowedLeaves.get(lc)) {
          recordTransfer << newLeaves.lLength;
          inverted << newLeaves.lLength;
          newLeaves << myLeaves (indexer.get(lc));
        }

      SortLists (&newLeaves, &recordTransfer);
      SortLists (&recordTransfer,&inverted);
      indexer.Clear();
      indexer.Duplicate (&inverted);
      myLeaves.Clear();
      myLeaves.Duplicate (&newLeaves);

        // finally check whether the root still has 3 or more children

      if (myCT->get_num_nodes()==2) {
        node <long>* promoteMe = nil;

        if (myCT->go_down(1)->get_num_nodes ()) {
          promoteMe = myCT->go_down(1);
        } else {
          promoteMe = myCT->go_down(2);
        }

        long nn = promoteMe->get_num_nodes();
        if (nn) {
          for (long cc = 1; cc <= nn; cc++) {
            myCT->add_node (*promoteMe->go_down(cc));
          }

          myCT->detach_child (promoteMe->get_child_num());
          DeleteObject ((BaseRef)promoteMe->in_object);
          delete promoteMe;
        } else {
          HandleApplicationError ("Internal tree pattern error in MatchTreePattern");
          return   "Unequal: Error";
        }
      }
    }

    _SimpleList * reindexer = nil;

    if (!indexer.Equal (otherIndexer)) {
      otherIndexer.Populate (_SimpleList ((unsigned long)myLeaves.countitems()).Populate (indexer, _SimpleList::action_invert),
                             _SimpleList::action_compose);
      reindexer = &otherIndexer;
    }

    char compRes;

    if (internalTreeCompare (myCT, otherCT, reindexer, 1, myLeaves.lLength, nil, compareTo, true)>0) {
      rerootAt = kCompareEqualWithoutReroot;
    } else {
      long   tCount = 0;


      node_iterator <long> ni (otherCT, _HY_TREE_TRAVERSAL_POSTORDER);
      node<long> * iterator = ni.Next();

      while (iterator!=otherCT) {
        if (!iterator->is_leaf()) {
          compRes = internalTreeCompare (myCT, iterator, reindexer, 1, myLeaves.lLength, nil, compareTo, true);
          if (compRes>0) {
            break;
          } else if (compRes) {
            iterator = otherCT;
            break;
          }
        }
        iterator = ni.Next();

        tCount ++;
      }

      if (iterator!=otherCT) {
        ni.Reset (compareTo->theRoot);
        iterator = ni.Next();
        while (iterator!=theRoot) {
          if (tCount==1) {
            rerootAt = kCompareEqualWithReroot & *map_node_to_calcnode (iterator)->GetName() & '.';
            break;
          } else {
            tCount --;
          }

          iterator = ni.Next();
        }
      }
    }
    if (rerootAt.empty()) {
      rerootAt = kCompareUnequalToplogies;
    }
  } else {
    rerootAt = kCompareUnequalLabelSets;
  }

  destroyCompTree (myCT);
  destroyCompTree (otherCT);
  return          rerootAt;
}

//__________________________________________________________________________________

void        _TreeTopology::ComputeClusterTable (_SimpleList& result, _SimpleList& pswRepresentation) const {
    
    long            leaf_count = pswRepresentation.Element(-2L),
    upper_bound = pswRepresentation.countitems()-2L,
    leafCode  = 0L,
    L = 0L,
    R = 0L;
    
    result.Clear    ();
    result.Populate (3*leaf_count,-1,0);
    
    
    for (long k = 0; k < upper_bound; k+=2L) {
        if (pswRepresentation.get(k) < leaf_count) { // is a leaf
            R = leafCode++;
        } else {
            long row;
            L = pswRepresentation.get (k-2*pswRepresentation.get(k+1));
            if (k + 2L == upper_bound /* root */
                || pswRepresentation.get(k+3L) == 0L) {
                row = R;
            } else {
                row = L;
            }
            
            result[row*3]   = L;
            result[row*3+1L] = R;
        }
    }
}
//__________________________________________________________________________________
_String*        _TreeTopology::ConvertFromPSW                       (_AVLListX& nodeMap,_SimpleList& pswRepresentation) const {
    _StringBuffer * result = new _StringBuffer (128L);
    
    if (pswRepresentation.countitems() > 4L) {
        long leafCount = pswRepresentation.Element (-2);
        // traverse backwards
        bool lastLeaf = false;
        _SimpleList     bounds;
        
        for (long k = pswRepresentation.countitems()-4L; k>=0L; k-=2L) {
            if (lastLeaf) {
                (*result) << ',';
            }
            if (pswRepresentation.list_data[k] >= leafCount) { //
                (*result) << ')';
                lastLeaf = false;
                bounds   << k-2*pswRepresentation.get(k+1);
            } else {
                _String nodeName (*(_String*)nodeMap.Retrieve(pswRepresentation.get(k)));
                nodeName.Flip();
                (*result) << nodeName;
                while (bounds.Element(-1) == k && bounds.nonempty()) {
                    (*result) << '(';
                    bounds.Pop();
                }
                lastLeaf = true;
            }
        }
    }
    result->TrimSpace();
    result->Flip();
    return result;
}
/*
 if (reference == false) {
 nodeMap.Clear();
 }
 
 pswRepresentation.Clear();
 
 long    leafIndex  = 0,
 iNodeCount = -1;
 
 _SimpleList levelBuffer;
 
 node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
 
 while (node<long> * currentNode = ni.Next (&levelBuffer)) {
 _String nodeName = GetNodeName (currentNode);
 
 while (levelBuffer.countitems() <= ni.Level()) {
 levelBuffer << 0;
 }
 
 if (currentNode->is_leaf()) {
 pswRepresentation << leafIndex;
 pswRepresentation << 0;
 if (reference) {
 long remapped = nodeMap.Find(&nodeName);
 if (remapped < 0) {
 return false;
 } else {
 remapped = nodeMap.GetXtra (remapped);
 if (remapped >= 0) {
 pswRepresentation << remapped;
 } else {
 return false;
 }
 }
 
 leafIndex++;
 } else {
 nodeMap.Insert(nodeName.makeDynamic(), leafIndex++, false);
 }
 } else {
 pswRepresentation << iNodeCount;
 pswRepresentation << levelBuffer.list_data[ni.Level()];
 if (reference) {
 pswRepresentation << 0;
 } else {
 (*inames) && &nodeName;
 }
 
 iNodeCount--;
 }
 if (ni.Level()) {
 levelBuffer.list_data[ni.Level()-1] += levelBuffer.list_data[ni.Level()]+1;
 }
 levelBuffer.list_data[ni.Level()]   = 0;
 }
 
 for (long k = 0; k < pswRepresentation.lLength; k+=(reference?3:2))
 if (pswRepresentation.list_data[k] < 0) {
 pswRepresentation.list_data[k] = leafIndex-pswRepresentation.list_data[k]-1;
 }
 
 pswRepresentation << leafIndex;
 pswRepresentation << (-iNodeCount-1);
 
 return true;
 
*/
//__________________________________________________________________________________
bool        _TreeTopology::ConvertToPSW (_AVLListX& nodeMap, _List* inames, _SimpleList& pswRepresentation, bool reference) const {
    
    if (reference == false) {
        nodeMap.Clear();
    }
    
    pswRepresentation.Clear();
    
    long    leafIndex  = 0L,
            iNodeCount = -1L;
    
    _SimpleList levelBuffer;
    
    node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
    
    while (node<long> * currentNode = ni.Next ()) {
        _String nodeName = GetNodeName (currentNode);
        
        
        if (levelBuffer.countitems() <= ni.Level()) {
            levelBuffer.AppendRange(ni.Level() - levelBuffer.countitems() + 1 , 0, 0);
        }
        
        if (currentNode->is_leaf()) {
            pswRepresentation << leafIndex << 0;
            if (reference) {
                long remapped = nodeMap.Find(&nodeName);
                if (remapped == kNotFound) {
                    return false;
                } else {
                    remapped = nodeMap.GetXtra (remapped);
                    if (remapped >= 0) {
                        pswRepresentation << remapped;
                    } else {
                        return false;
                    }
                }
                
                leafIndex++;
            } else {
                nodeMap.Insert(nodeName.makeDynamic(), leafIndex++, false);
            }
        } else {
            pswRepresentation << iNodeCount << levelBuffer.get(ni.Level());
            if (reference) {
                pswRepresentation << 0L;
            } else {
                (*inames) && &nodeName;
            }
            iNodeCount--;
        }
        if (ni.Level()) {
            levelBuffer[ni.Level()-1] += levelBuffer.get(ni.Level())+1;
        }
        levelBuffer[ni.Level()]   = 0;
    }
    
    for (long k = 0; k < pswRepresentation.lLength; k+=(reference?3:2)) {
        if (pswRepresentation.get(k) < 0L) {
            pswRepresentation[k] = leafIndex-pswRepresentation.get(k)-1;
        }
    }
    
    pswRepresentation << leafIndex << (-iNodeCount-1);
    
    return true;
}


//__________________________________________________________________________________

_AssociativeList *   _TreeTopology::SplitsIdentity (HBLObjectRef p, HBLObjectRef cache)  const {
    // compare tree topologies
    _Matrix     * result  = new _Matrix (2,1,false,true),
                * result2 = nil;
    
    _String    * tree_repr  = nil;
    
    long leaves, internals;
    EdgeCount (leaves, internals);
    
    result->theData[0] = internals - 1.;
    result->theData[1] = -1.;
    
    if (p && (p->ObjectClass() == TOPOLOGY || p->ObjectClass() == TREE)) {
        _List           avlSupport,
        iNames;
        
        _AVLListX       nameMap  (&avlSupport);
        
        _SimpleList     psw, psw2, clusters, inodeList;
        
        
        ConvertToPSW            (nameMap, &iNames, psw);
        ComputeClusterTable     (clusters, psw);
        
        if (((_TreeTopology*)p)->ConvertToPSW           (nameMap, nil, psw2, true)) {
            _SimpleList           workSpace;
            long leafCount      = psw.Element (-2),
                 upper_bound    = psw2.countitems() - 3L;
            
            for (long k = 0L; k < upper_bound; k+=3L) {
                
                if (psw2.list_data[k] < leafCount) {
                    workSpace << 1L
                              << 1L
                              << psw2.get(k+2L)
                              << psw2.get(k+2L);
                } else {
                    _SimpleList quad;
                    
                    quad << leafCount+1 << 0 << 0 << 1;
                    
                    long  w = psw2.get (k+1);
                    while (w > 0) {
                        _SimpleList quad2;
                        quad2 << workSpace.Pop() << workSpace.Pop() << workSpace.Pop() << workSpace.Pop();
                        w -= quad2.get (3);
                        quad.list_data[0] = Minimum (quad2.get(0),quad.get(0));
                        quad.list_data[1] = Maximum (quad2.get(1),quad.get(1));
                        quad.list_data[2] += quad2.get(2);
                        quad.list_data[3] += quad2.get(3);
                    }
                    
                    if (quad.list_data[2] == quad.list_data[1] - quad.list_data[0] + 1) {
                        if (clusters.list_data[3*quad.list_data[0]] == quad.list_data[0] && clusters.list_data[3*quad.list_data[0]+1] == quad.list_data[1]) {
                            clusters.list_data[3*quad.list_data[0]+2] = 1;
                        } else if (clusters.list_data[3*quad.list_data[1]] == quad.list_data[0] && clusters.list_data[3*quad.list_data[1]+1] == quad.list_data[1]) {
                            clusters.list_data[3*quad.list_data[1]+2] = 1;
                        }
                    }
                    quad.Flip();
                    workSpace << quad;
                }
            }
            
            psw2.Clear();
            long matchCount = 0,
            iNodeCount = 0;
            
            long L, R = -1L;
            
            _SimpleList leafSpans (leafCount,0,0),
            iNodesTouched;
            
            for (unsigned long k = 0; k < psw.lLength-2; k+=2) {
                if (psw.list_data[k] < leafCount) {
                    R = psw.list_data[k];
                    psw2 << R;
                    psw2 << 0;
                    leafSpans.list_data[R] = (psw2.lLength>>1);
                } else {
                    long ll = k-2*psw.list_data[k+1];
                    L       = psw.list_data[ll];
                    if ((clusters.list_data[3*L] == L && clusters.list_data[3*L+1] == R && clusters.list_data[3*L+2] > 0)
                        || (clusters.list_data[3*R] == L && clusters.list_data[3*R+1] == R && clusters.list_data[3*R+2] > 0)) {
                        L = (psw2.lLength>>1) - leafSpans.list_data[L] + 1;
                        psw2 << leafCount+iNodeCount++;
                        psw2 << L;
                        
                        iNodesTouched << psw.list_data[k];
                    }
                }
            }
            
            for (unsigned long k = 0; k < psw2.lLength; k+=2)
                if (psw2.list_data[k] < leafCount) {
                    psw2.list_data[k+1] = 0;
                } else {
                    matchCount++;
                }
            
            psw2 << leafCount;
            psw2 << iNodeCount;
            
            result->theData[0] = psw.Element (-1);
            result->theData[1] = matchCount;
            
            
            tree_repr  = ConvertFromPSW (nameMap, psw2);
            
            _List sharedNames;
            for (long k = 0; k < iNodesTouched.lLength; k++) {
                sharedNames << iNames (iNodesTouched(k)-leafCount);
            }
            
            result2 = new _Matrix (sharedNames);
        }
        
    }
    
    //DeleteObject (bc);
    
    _AssociativeList * resultList = new _AssociativeList;
    resultList->MStore ("CLUSTERS", result, false);
    if (result2) {
        resultList->MStore ("NODES",    result2, false);
    }
    resultList->MStore ("CONSENSUS", tree_repr ? new _FString (tree_repr) : new _FString(), false);
    return resultList;
}

//__________________________________________________________________________________

const _String  _TreeTopology::FindMaxCommonSubTree (_TreeTopology const*  compare_to, long& size_tracker, _List* forest) const{
    _List           myLeaves,
                    otherLeaves,
                    sharedLeaves;
    
    _SimpleList     indexer,
                    otherIndexer,
                    sharedLeavesIDin1,
                    sharedLeavesIDin2;
    
    
    _String         rerootAt;
    
    node<long>      *myCT    = prepTree4Comparison(myLeaves, indexer),
                    *otherCT = compare_to->prepTree4Comparison(otherLeaves, otherIndexer);
    
    sharedLeaves.Intersect (otherLeaves,myLeaves,&sharedLeavesIDin1,&sharedLeavesIDin2);
    
    if (sharedLeaves.countitems() > 1UL) { // more than one common leaf
                                           // now we need to map shared leaves to a common indexing space
        
        _SimpleList     my_inverted_index, // was lidx1
                        other_inverted_index; // was lidx2
        
        my_inverted_index.Populate(indexer, _SimpleList::action_invert);
        other_inverted_index.Populate(otherIndexer, _SimpleList::action_invert);
        
        _SimpleList     reindexer (otherLeaves.countitems(),-1L,0L),
                        ldx1,
                        ldx2;
        
        
        for (long k2=0L; k2<sharedLeaves.countitems(); k2++) {
            reindexer [other_inverted_index.get (sharedLeavesIDin2.get (k2))] = my_inverted_index.get (sharedLeavesIDin1.get (k2));
            //          the original index of shared leaf k2  -> the original index of shared lead k1 (in post-order traveral)
        }
        
        // now we map actual leaf structures to their respective leaf indices
        
        node_iterator<long> ni (myCT, _HY_TREE_TRAVERSAL_POSTORDER);
        
        while (node<long>* iterator = ni.Next()) {
            if (iterator->is_leaf()) {
                ldx1 << (long)iterator;
            }
        }
        
        ni.Reset (otherCT);
        while (node<long>* iterator = ni.Next()) {
            if (iterator->is_leaf()) {
                ldx2 << (long)iterator;
            }
        }
        
        // now we loop through the list of leaves and try to match them all up
        
        _SimpleList     matchedTops,
                        matchedSize;
        
        for (long k3=0UL; k3<sharedLeaves.countitems()-1L; k3++) {
            long try_leaf_index = other_inverted_index.get (sharedLeavesIDin2.get (k3));
            
            if (reindexer.get (try_leaf_index) >= 0L) {
                // leaf still available
                node<long>*      ln1 = (node<long>*)ldx1.get (my_inverted_index.get (sharedLeavesIDin1.get (k3))),
                         *       ln2 = (node<long>*)ldx2.get (try_leaf_index),
                         *       p1  = ln1->parent,
                         *       p2  = ln2->parent;
                
                bool             comparison_result = false;
                
                while ((internalTreeCompare (p1,p2,&reindexer,0,myLeaves.lLength,p2->parent?nil:ln2,compare_to) == 1)&&p1&&p2) {
                    ln1 = p1;
                    ln2 = p2;
                    p1=p1->parent;
                    p2=p2->parent;
                    
                    comparison_result = true;
                }
                
                if (comparison_result) {
                    _SimpleList* matchedLeaves = (_SimpleList*)ln2->in_object;
                    
                    matchedTops << (long)ln1;
                    matchedSize << matchedLeaves->countitems();
                    
                    matchedLeaves->Each ([&] (long value, unsigned long) -> void {
                        reindexer[value] = -1;
                    });
                }
            }
        }
        
        if (matchedSize.nonempty()) {
            
            auto handle_subtree = [&] (node<long> const * top_node) -> _String * {
                ni.Reset(myCT);

                long internal_node_count = 0L;
                
                while (node<long>* iterator = ni.Next()) {
                    if (!iterator->is_leaf()) {
                        if (iterator == top_node) {
                            break;
                        }
                        internal_node_count ++;
                    }
                }
                ni.Reset(theRoot);
                while (node<long>* iterator = ni.Next()) {
                    if (!iterator->is_leaf()) {
                        if (internal_node_count == 0) {
                            return LocateVar(iterator->in_object)->GetName();
                        }
                        internal_node_count--;
                    }
                }
                return nil;
            };
            
            if (forest) {
                size_tracker = 0L;
                for (long k6=0L; k6<matchedSize.lLength; k6++) {
                    (*forest) << handle_subtree ((node<long>*)matchedTops.get (k6));
                    size_tracker += matchedSize.get (k6);
                }
            } else {
                long maxSz   = -1L,
                     maxIdx  =  0L;
                
                matchedSize.Each ([&] (long value, unsigned long index) -> void {
                    if (StoreIfGreater(maxSz, value)) {
                        maxIdx = index;
                    }
                });
                
                size_tracker = maxSz;
                _String const * res = handle_subtree ((node<long>*)matchedTops.list_data[maxIdx]);
                destroyCompTree(myCT);
                destroyCompTree(otherCT);
                return *res;
            }
        }
    }
    destroyCompTree(myCT);
    destroyCompTree(otherCT);
    return kEmptyString;
}


//__________________________________________________________________________________

const _String  _TreeTopology::CompareSubTrees (_TreeTopology const* compareTo, node<long>* topNode) const {
    // compare tree topologies
    
    _List           myLeaves,
    otherLeaves,
    sharedLeaves;
    
    _SimpleList     indexer,
    otherIndexer,
    sharedLeavesID;
    
    node<long>      *myCT,
    *otherCT;
    
    _String         rerootAt;
    
    myCT            = prepTree4Comparison(myLeaves, indexer);
    otherCT         = compareTo->prepTree4Comparison(otherLeaves, otherIndexer, topNode);
    
    sharedLeaves.Intersect (myLeaves, otherLeaves,&sharedLeavesID);
    
    // first compare the inclusion for the set of leaf labels
    
    if (sharedLeavesID.countitems() == otherLeaves.countitems()) {
        _SimpleList ilist;
        
        ilist.Populate(myLeaves.countitems(), -1L, 0L);
    
        for (long k2 = 0L; k2 < otherIndexer.countitems(); k2++) {
            ilist [sharedLeavesID.get(otherIndexer.get(k2))] = k2;
        }
        
        for (long k2 = 0L; k2<indexer.countitems(); k2++) {
            long         lidx      = ilist.list_data[indexer.get(k2)];
            indexer[k2] = lidx >=0L ? lidx : -1L;
        }
        
        _SimpleList *reindexer = &indexer;
        
        // now compare explore possible subtree matchings
        // for all internal nodes of this tree except the root
        
        char   compRes = 0;
        
        long   tCount = 1L,
        nc2 = topNode->get_num_nodes();
        
        node_iterator<long> ni (myCT, _HY_TREE_TRAVERSAL_POSTORDER);
        ni.Next();
        node<long>* iterator = ni.Next();
        
        while (iterator!=myCT) {
            long nc = iterator->get_num_nodes();
            if (nc == nc2) {
                long kk;
                for (kk = 1; kk <= nc; kk++) {
                    compRes = internalTreeCompare (otherCT,iterator, reindexer, 0, otherLeaves.lLength, iterator->go_down(kk),compareTo);
                    if (compRes) {
                        if (compRes == -1) {
                            iterator = myCT;
                        }
                        break;
                    }
                }
                if (kk>nc) {
                    compRes = internalTreeCompare (otherCT,iterator, reindexer, 0, otherLeaves.lLength, nil, compareTo);
                    if (compRes) {
                        if (compRes == -1) {
                            iterator = myCT;
                        }
                        break;
                    }
                } else {
                    break;
                }
            }
            tCount ++;
            iterator = ni.Next();
        }
        
        if (iterator != myCT) {
            ni.Reset (theRoot);
            iterator = ni.Next ();
            while (iterator != theRoot) {
                if (tCount==0) {
                    rerootAt = _String("Matched at the ") & *map_node_to_calcnode (iterator)->GetName() & '.';
                    break;
                } else {
                    tCount --;
                }
                
                iterator = ni.Next();
            }
        } else {
            long nc = myCT->get_num_nodes();
            
            if (nc == nc2+1)
                for (long kk = 1; kk <= nc; kk++) {
                    compRes = internalTreeCompare (otherCT,myCT, reindexer, 0, otherLeaves.lLength, myCT->go_down(kk),compareTo);
                    if (compRes==1) {
                        break;
                    }
                }
            
            if (compRes == 1) {
                rerootAt = "Matched at the root.";
            }
        }
        
        if (rerootAt.empty()) {
            rerootAt = "No match: Different topologies (matching label sets).";
        }
    } else {
        rerootAt = "No match: Unequal label sets.";
    }
    
    destroyCompTree (myCT);
    destroyCompTree (otherCT);
    
    return          rerootAt;
}






