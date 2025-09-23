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
  /**
   * @brief Construct a new _TreeIterator object
   *
   * @param source The tree to iterate over
   * @param traversal_type The type of traversal
   */
  _TreeIterator (_TheTree const * source, int traversal_type);
  /**
   * @brief Construct a new _TreeIterator object
   *
   * @param root The root of the tree
   * @param root_node The root node
   * @param traversal_type The type of traversal
   */
  _TreeIterator (_CalcNode const * root,  node<long>* root_node, int traversal_type);
  /**
   * @brief Destroy the _TreeIterator object
   */
  ~_TreeIterator (void);
  /**
   * @brief Reset the iterator
   */
  void                    Reset (void);
  
  /**
   * @brief Get the next node in the traversal
   *
   * @return _CalcNode* The next node
   */
  _CalcNode *             Next (void);
  /**
   * @brief Get the current node in the traversal
   *
   * @return _CalcNode* The current node
   */
  _CalcNode *             Current (void) const;
  /**
   * @brief Get the depth of the current node
   *
   * @return long The depth
   */
  long                    Depth (void) const;
  /**
   * @brief Get the history of the traversal
   *
   * @return const _SimpleList& The history
   */
  const _SimpleList&      History (void) const;
  /**
   * @brief Check if the iterator is at a leaf
   *
   * @return true if the iterator is at a leaf, false otherwise
   */
  bool  IsAtLeaf          (void) const;
  /**
   * @brief Check if the iterator is at the root
   *
   * @return true if the iterator is at the root, false otherwise
   */
  bool  IsAtRoot          (void) const;
  /**
   * @brief Get the current node
   *
   * @return node<long>* The current node
   */
  node <long>* GetNode    (void) const { return iterator.Current(); }
  
};


#endif
