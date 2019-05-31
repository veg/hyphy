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


#include "global_things.h"
using namespace hy_global;

#include "tree_iterator.h"

//_______________________________________________________________________________________________


_TreeIterator::_TreeIterator (_TheTree const* source, int traversal_type): iterator(source->theRoot, traversal_type & fTreeIteratorTraversalMask) {
    root_n = source->theRoot;
    root_node =  map_node_to_calcnode (source->theRoot);
    flags  = traversal_type;
}

_TreeIterator::_TreeIterator (_CalcNode const* root, node<long>* nroot, int traversal_type): iterator(nroot, traversal_type & fTreeIteratorTraversalMask) {
    root_node = root;
    root_n = nroot;
    flags  = traversal_type;
}

_TreeIterator::~_TreeIterator (void) {};

void     _TreeIterator:: Reset (void) {
    iterator.Reset (root_n);
}

_CalcNode *     _TreeIterator:: Next (void) {
    
    node<long> * nn = iterator.Next(&history);
    
    if (nn) {
        if (nn->is_root() && (flags & fTreeIteratorTraversalSkipRoot)) {
            return Next();
        }
        if (!nn->is_leaf() && (flags & fTreeIteratorTraversalLeaves)) {
            return Next();
        }
        
        return map_node_to_calcnode(nn);
    }
    return nil;
}

_CalcNode *     _TreeIterator::Current (void) const {
    return map_node_to_calcnode(iterator.Current());
}

long     _TreeIterator::Depth (void) const {
    return iterator.Level();
}

const _SimpleList&  _TreeIterator::History (void) const {
    return history;
}

bool  _TreeIterator::IsAtLeaf          (void) const {
    node <long>* const current = iterator.Current();
    if (current) {
        return current->get_num_nodes() == 0L;
    }
    return false;
}

bool  _TreeIterator::IsAtRoot          (void) const {
    node <long>* const current = iterator.Current();
    if (current) {
        return current == root_n;
    }
    return false;
}




