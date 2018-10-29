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

#ifndef     __TREE_ITERATOR__
#define     __TREE_ITERATOR__

#include "tree.h"

#define       fTreeIteratorTraversalMask      0x07 // 111 binary
#define       fTreeIteratorTraversalSkipRoot  0x10 //
#define       fTreeIteratorTraversalLeaves    0x20

class _TreeIterator {
private:
  
  node_iterator<long> iterator;
  int                 flags;
  _SimpleList         history;
  _CalcNode const     *root_node;
  node<long>          *root_n;
  
  
public:
  _TreeIterator (_TheTree const * source, int traversal_type);
  _TreeIterator (_CalcNode const * root,  node<long>* root_node, int traversal_type);
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


#endif
