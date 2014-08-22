/*
 
 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (spond@ucsd.edu)
 Steven Weaver (sweaver@ucsd.edu)
 Martin Smith (martin.audacis@gmail.com)
 
 Module Developers:
 Art FY Poon    (apoon@cfenet.ubc.ca)
 Lance Hepler (nlhepler@gmail.com)
 
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

#ifndef _AVLLIST_
#define _AVLLIST_
//#pragma once
#include "hy_list_orderable.h"
#include "hy_list_numeric.h"
#include "helperfunctions.h"
#include "hy_stack.h"
#include "defines.h"

#define HY_AVL_LEAF (-1L)
#define HY_AVL_MOVE_LEFT (-1L)
#define HY_AVL_MOVE_RIGHT (1L)

/**
 * A helper data object representing a tree node;
 * leftChild and rightChild are 0-based indices into the avl_structure array or HY_AVL_LEAF (-1L)
 * if there are no children
 */

struct _AVLListHelper {
  long left_child,
  right_child,
  balance_factor;
  
  _AVLListHelper (void) {
    left_child = HY_AVL_LEAF;
    right_child = HY_AVL_LEAF;
    balance_factor = 0L;
  }
  
  bool operator == (const _AVLListHelper&  rhs) {
    return left_child == rhs.left_child && right_child == rhs.right_child && balance_factor == rhs.balance_factor;
  }
};

typedef _hyStack <long> _AVLListTraversalHistory;


/**
 * This object implements a key:value dictionary functionality using
 * a balanced binary search tree.
 * @param KEYTYPE is the type used to represent keys; it must support stack copy, all pairwise comparison operations, i.e. ==, <= etc
 * @param PAYLOAD is the type which represents dictionary values; it must support stack copy and assignment (very little else is assumed)
 */


//_____________________________________________________________________________
template <typename KEYTYPE>
class _AVLListBase : public virtual BaseObj{
  
protected:
  
  /**
   * A list representing the flattened tree;
   * @sa root
   * @sa empty_slots
   */
  
  _hyList             <_AVLListHelper> avl_structure;
  
  /**
   * This list stores the keys present in the list
   * they are tied to tre tree nodes by index, i.e. the key stored in the i-th index
   */
  
  // _hyListOrderable    <KEYTYPE>        keys;
  
  /**
   * This list stores the payloads associated with each key
   * they are tied by index, i.e. the key stored in the i-th index
   * of avl_structure is bound to the i-th value in this list
   */
  
  //_hyList             <PAYLOAD>        values;
  
  
  /**
   * Stores indices of deleted nodes in avl_structure,
   * so that they can be used by subsequent insert operations
   */
  
  _SimpleList         empty_slots;
  
  
  /**
   * Stores the index of the root; if HY_AVL_LEAF, then the list is empty
   */
  
  long root;
  
  /**
   * this is a convenience function which moves from the current node in the tree
   * based on the result of a key comparison
   * @param node the index in the tree; no range checking is performed
   * @param direction move left if < 0 (less than), right if > 0 (greater then)
   * @return the index of the node moved to
   */
  
  inline long _MoveInTree (long node, long direction) const ;
  
  /**
   * Find the index of the smallest/largest key in the subtree rooted at node
   * @param node start traversal here (if HY_AVL_LEAF, the return value is HY_NOT_FOUND)
   * @param direction go left (small values) or right (large values)
   * @param history if not null, store the path leading to the extereme element here
   * @sa First
   * @sa Last
   * @return the >=0 index of smallest/largest key in this tree or HY_NOT_FOUND if not found
   */
  
  long   _DescendToTerminal (long node, long direction, _AVLListTraversalHistory * history = NULL) const;
  
  /**
   * Compares the key at a given index to a fixed value
   * Must be implemented by a derived class
   * @param node The index of the item to compare
   * @param key The value to compare the index to
   * @return -1 if keys[node] > key , 0 if keys[node] == key, or 1 if keys[node] < key
   */
  
  virtual long _CompareIndexToValue (long node, KEYTYPE const & key) const = 0;
  
  /**
   * Stores a key in the list at a given index
   * Must be implemented by a derived class
   * @param key The key to store (assumed to NOT already be in the tree, otherwise things will break)
   * @param index The index to store at (pass HY_LIST_INSERT_AT_END to find an unused index to store)
   * @return the index where key was stored
   */
  
  virtual long _StoreKey   (const KEYTYPE& key, long index = HY_LIST_INSERT_AT_END);
  
  /**
   * Remove the key from 'index' (note that the tree structure is handled by the ::Delete function)
   * Must be implemented by a derived class
   * @param index The index to delete from at
   */
  virtual void _RemoveKey  (long index);
  
  /**
   * Remove a kay/value pair from 'index' (note that the tree structure is handled by the ::Delete function)
   * Must be implemented by a derived class
   * @param index The index to delete from
   */
  
  inline long& RightChild (const long index) const;
  inline long& LeftChild (const long index) const ;
  inline long& BalanceFactor (const long index) const;
  
  void   _DeleteHelper (const long index, const long new_node, _SimpleList const &directions, _SimpleList const &nodes, bool update_root = true);
  
  virtual KEYTYPE const * _AtIndex (const unsigned long) const = 0;
  
public:
  
  
  /**
   * Construct an empty AVLList
   */
  
  _AVLListBase     (void);
  
  /**
   * Stack copy constructor
   */
  
  _AVLListBase     (_AVLListBase <KEYTYPE> const & );
  
 
  /**
   * Object destuctor
   */
  
  virtual ~_AVLListBase(void);
  
  /**
   * Initialize a freshly created object
   */
  virtual void Initialize (bool = false);

  
  /**
   * Clone an existing object into this object, which has been previously initilized
   * @param source the object to clone from
   */
  void Clone (_AVLListBase <KEYTYPE> const & source);
  
  /**
   * Length of the list
   * @return The number of elements in the list
   */
  inline unsigned long Length (void) const;
  
  /**
   * The basic search function
   * @param key The value to find
   * @param hist if not null, then the indices traversed on the way to this key will be stored
   * @return the >=0 index of the key in this list or HY_NOT_FOUND if not found
   */
  
  
  virtual long   Find(const KEYTYPE& key, _AVLListTraversalHistory * history = NULL) const;
  
  
  /**
   * Find the key in the tree; if not found, store the last key traversed
   * @param key The value to find
   * @param last_index store the index of the last traversed key here
   * @return the result of comparing the search key with the key stored at last_index
   */
  
  
  virtual long   FindBest (const KEYTYPE& key, long& last_index) const;
  
  // tree traversal functions
  // 20140808: SLKP TODO Move to an iterator class
  
  /**
   * Find the index of the smallest key in the tree
   * @param history if not null, store the path leading to the first element here
   * @return the >=0 index of smallest key in this tree or HY_NOT_FOUND if not found
   */
  
  long First(_AVLListTraversalHistory * history = NULL) const;
  
  /**
   * Find the index of the largest key in the tree
   * @param history if not null, store the path leading to the last element here
   * @return the >=0 index of largest key in this tree or HY_NOT_FOUND if not found
   */
  
  long Last(_AVLListTraversalHistory * history = NULL) const;
  
  /**
   * Find the index of the next largest key in the tree given traversal history
   * @param history The traversal history up to now (empty
   * @return the >=0 index of the next key in this tree or HY_NOT_FOUND if not found
   */
  
  long Next (_AVLListTraversalHistory& history) const;
  
  /**
   * Find the index of the next smallest key in the tree given traversal history
   * @param history The traversal history up to now (empty history returns HY_NOT_FOUND)
   * @return the >=0 index of the previous key in this tree or HY_NOT_FOUND if not found
   */
  
  long Prev (_AVLListTraversalHistory& history) const;
  
  /**
   * Traverse the tree and return key indices in sequences
   * @param history The traversal history up to now (empty history returns initializes the traversal)
   * @param reverse go in reverse (largest to smallest)
   * @return the >=0 index of the next key in sequence this tree or HY_NOT_FOUND if not found
   */
  
  virtual long Traverser(_AVLListTraversalHistory &, bool reverse = false) const;
  
  // insertion and deletion functions
  
  
  /**
   * Traverse the tree and insert the key-value pair (if not found); update and rebalance the tree
   * @param key The key to insert
   * @return the index>=0 if the key/value was inserted at position index, returns -index-1 if the key
   * found in position index
   */
  
  virtual long Insert  (KEYTYPE const& key);
  
  /**
   * Traverse the tree and delete the key (if found); update and rebalance the tree
   * @param key The key to delete
   * @return the index of the key deleted (>=0) of HY_NOT_FOUND the key was not in the dict
   */
  
  virtual long Delete  (KEYTYPE const& key);
  
  
  /**
   * An operator version of Find
   * @param key The key to find
   * @return the index of the key if found (>=0) of HY_NOT_FOUND the key was not in the dict
   */
  virtual long operator [] (const KEYTYPE& key) const;
  // virtual PAYLOAD* operator () (long index) = 0;
  

  /**
   * Compare two lists for equality (i.e. do they represent the same set of keys)
   * @param rhs The _AVLList to compare this one to
   * @return are the lists equal
   */
  virtual bool operator == (_AVLListBase <KEYTYPE> const & rhs) const;

  /**
   * Compare two lists for inequality (i.e. do they represent the same set of keys)
   * @param rhs The _AVLList to compare this one to
   * @return are the lists UNequal
   */
  virtual bool operator != (_AVLListBase <KEYTYPE> const & rhs) const;
  
  /**
   * Return the key at a given index
   * @param key the index
   * @return a pointer to the key, or NULL if the index is invalid of is poitning at a vacated
   * space
   */
  virtual KEYTYPE const * AtIndex (const unsigned long) const;
  
  /*virtual void Clear(bool = false);
   virtual bool HasData(long);
   
   virtual void ReorderList(_SimpleList * = nil);
   virtual long InsertData(BaseRef, long, bool);
   virtual BaseRef toStr(void);
   virtual void DeleteXtra(long); // {}
   
   virtual void DeleteAll(bool cL);
   //{
   //Clear(cL);
   //DeleteObject(dataList);
   //}
   
   
   // 20100623: a shortcut function to look for integers only
   // avoids calling ::Compare
   
   
   long GetByIndex(const long);
   
   // the 1st bool flag is to say whether to dup the object being inserted
   // the 2nd bool flag (if the first flag is false) if set to true,
   // will cause failed inserts (key already exists) to delete the key
   long Insert(BaseRef, long = 0, bool = true, bool = false);
   
   BaseRef Retrieve(long);
   
   void Delete(BaseRef, bool = false); */
  
   const char* ConsistencyCheck(void) const;
   void EchoList        (void) const;
   void _EchoListHelper  (long, long) const;
  
};

//*************** CONSTRUCTORS ***************//


template <typename KEYTYPE>
_AVLListBase<KEYTYPE>::_AVLListBase(void) {
  this->Initialize();
}

template <typename KEYTYPE>
_AVLListBase<KEYTYPE>::_AVLListBase(_AVLListBase<KEYTYPE> const & source) {
  //keys.Clone (&source.keys);
  //values.Clone (&source.values);
  this->Initialize();
  this->Clone (source);
}



//*************** INITIALIZER and CLONER ***************//
template <typename KEYTYPE>
void _AVLListBase<KEYTYPE>::Clone(_AVLListBase<KEYTYPE> const & source) {
  this->avl_structure.Clone (&source.avl_structure);
  this->empty_slots.Clone (&source.empty_slots);
  this->root = source.root;
}

template <typename KEYTYPE>
void _AVLListBase<KEYTYPE>::Initialize(bool) {
  this->root = HY_AVL_LEAF;
}


//*************** DESTRUCTOR ***************//

template <typename KEYTYPE>
_AVLListBase<KEYTYPE>::~_AVLListBase(void) {
  
}


//*************** ATTRIBUTE/ACCESSOR ACCESSORS ***************//

template <typename KEYTYPE>
inline unsigned long _AVLListBase<KEYTYPE>::Length (void) const {
  return this->avl_structure.Length() - this->empty_slots.Length();
}

template <typename KEYTYPE>
inline long& _AVLListBase<KEYTYPE>::LeftChild(const long index) const {
  return this->avl_structure.AtIndex (index).left_child;
}

template <typename KEYTYPE>
inline long& _AVLListBase<KEYTYPE>::RightChild(const long index) const {
  return this->avl_structure.AtIndex (index).right_child;
}

template <typename KEYTYPE>
inline long& _AVLListBase<KEYTYPE>::BalanceFactor(const long index) const {
  return this->avl_structure.AtIndex (index).balance_factor;
}

//*************** SEARCH FUNCTIONS  ***************//

template <typename KEYTYPE>
long _AVLListBase<KEYTYPE>::operator [] (const KEYTYPE& key) const{
  return this->Find (key);
}

template <typename KEYTYPE>
long _AVLListBase<KEYTYPE>::Find (const KEYTYPE& key, _AVLListTraversalHistory * history) const{
  long current_node = this->root;
  
  while (current_node != HY_AVL_LEAF) {
    long comp = this->_CompareIndexToValue (current_node, key);
    if (comp) {
      if (history) {
        history->push (current_node);
      }
      current_node = this->_MoveInTree (current_node, comp);
    } else {
      return current_node;
    }
  }
  
  return HY_NOT_FOUND;
}

template <typename KEYTYPE>
long _AVLListBase<KEYTYPE>::FindBest (const KEYTYPE& key, long &last_index ) const{
  long current_node = root,
  comp = HY_AVL_MOVE_RIGHT;
  
  while (current_node != HY_AVL_LEAF) {
    comp = this->_CompareIndexToValue (current_node, key);
    last_index = current_node;
    
    if (comp == 0L) {
      return 0L;
    }
    
    current_node = this->_MoveInTree (current_node, comp);
    
  }
  
  return comp;
}


//*************** HELPER FUNCTIONS ***************//

template <typename KEYTYPE>
long _AVLListBase<KEYTYPE>::_MoveInTree (long node, long direction) const{
  if (direction < 0L) {
    return this->avl_structure.AtIndex(node).left_child;
  } else {
    if (direction > 0L) {
      return this->avl_structure.AtIndex(node).right_child;
    }
  }
  
  return node;
}

template <typename KEYTYPE>
KEYTYPE const* _AVLListBase<KEYTYPE>::AtIndex (unsigned long index) const{
  if (index < this->avl_structure.Length()) {
    if (empty_slots.Length()) {
      if (empty_slots.Find (index) != HY_NOT_FOUND) {
        return NULL;
      }
    }
    return this->_AtIndex (index);
  }
  return NULL;
}


//*************** TREE TRAVERSAL FUNCTIONS ***************//

//______________________________________________________________________________
template <typename KEYTYPE>
long _AVLListBase<KEYTYPE>::_DescendToTerminal(long current_index, long direction, _AVLListTraversalHistory * history) const {
  
  long the_index = HY_NOT_FOUND;
  
  while (current_index != HY_AVL_LEAF) {
    if (history) {
      history->push (current_index);
    }
    the_index = current_index;
    current_index = this->_MoveInTree (current_index, direction);
  }
  
  return the_index;
}


//______________________________________________________________________________
template <typename KEYTYPE>
long _AVLListBase<KEYTYPE>::First(_AVLListTraversalHistory * history) const {
  return this->_DescendToTerminal (this->root, HY_AVL_MOVE_LEFT, history);
}

//______________________________________________________________________________
template <typename KEYTYPE>
long _AVLListBase<KEYTYPE>::Last(_AVLListTraversalHistory * history) const {
  return this->_DescendToTerminal (this->root, HY_AVL_MOVE_RIGHT, history);
}

//______________________________________________________________________________
template <typename KEYTYPE>
long _AVLListBase<KEYTYPE>::Next(_AVLListTraversalHistory& history) const {
  
  if (history.Length()) {
    long current_node = history.pop(),
    try_node = this->_MoveInTree (current_node, HY_AVL_MOVE_RIGHT);
    if (try_node != HY_AVL_LEAF) {
      return this->_DescendToTerminal (try_node, HY_AVL_MOVE_LEFT, &history);
    } else {
      while (history.Length()) {
        long parent_node = history.pop();
        try_node = this->_MoveInTree (parent_node, HY_AVL_MOVE_RIGHT);
        if (try_node == current_node) {
          current_node = parent_node;
        } else {
          history.push (parent_node);
          return parent_node;
        }
      }
    }
  }
  
  return HY_NOT_FOUND;
}

//______________________________________________________________________________
template <typename KEYTYPE>
long _AVLListBase<KEYTYPE>::Prev(_AVLListTraversalHistory& history) const {
  
  if (history.Length()) {
    long current_node = history.pop(),
    try_node = this->_MoveInTree (current_node, HY_AVL_MOVE_LEFT);
    if (try_node != HY_AVL_LEAF) {
      return this->_DescendToTerminal (try_node, HY_AVL_MOVE_RIGHT, &history);
    } else {
      while (history.Length()) {
        long parent_node = history.pop();
        try_node = this->_MoveInTree (parent_node, HY_AVL_MOVE_LEFT);
        if (try_node == current_node) {
          current_node = parent_node;
        } else {
          history.push (parent_node);
          return parent_node;
        }
      }
    }
  }
  
  return HY_NOT_FOUND;
}

//______________________________________________________________________________
template <typename KEYTYPE>
long _AVLListBase<KEYTYPE>::Traverser(_AVLListTraversalHistory &history, bool reverse) const {
  if (reverse) {
    if (history.Length()) {
      return this->Prev (history);
    } else {
      return this->Last (&history);
    }
    return HY_NOT_FOUND;
    
  }
  if (history.Length()) {
    return this->Next (history);
  } else {
    return this->First (&history);
  }
  return HY_NOT_FOUND;
}


//*************** LIST COMPARISON ***************//

//______________________________________________________________________________
template <typename KEYTYPE>
bool _AVLListBase<KEYTYPE>::operator == (_AVLListBase <KEYTYPE> const& rhs) const {
  if (Length () == rhs.Length()) {
    _AVLListTraversalHistory historyLHS,
                             historyRHS;
    
    long indexLHS = Traverser (historyLHS),
         indexRHS = rhs.Traverser (historyRHS);
    
    while (indexLHS != HY_NOT_FOUND && indexRHS != HY_NOT_FOUND) {
      
      if (this->_CompareIndexToValue (indexLHS, *rhs.AtIndex (indexRHS)) != 0L) {
        return false;
      }
      indexLHS = Next (historyLHS);
      indexRHS = rhs.Next (historyRHS);
    }
    
    return (indexLHS == HY_NOT_FOUND && indexRHS == HY_NOT_FOUND);
  }
  return false;
}

//______________________________________________________________________________
template <typename KEYTYPE>
bool _AVLListBase<KEYTYPE>::operator != (_AVLListBase <KEYTYPE> const& rhs) const {
  return ! (*this == rhs);
}

//*************** INSERTION AND DELETION FUNCTIONS ***************//

//______________________________________________________________________________
template <typename KEYTYPE>
long _AVLListBase<KEYTYPE>::_StoreKey   (const KEYTYPE&, long index) {
  if (index == HY_LIST_INSERT_AT_END) {
    if (empty_slots.Length()) {
      index = empty_slots.Pop(true);
    } else {
      index = avl_structure.Length();
      avl_structure.append (_AVLListHelper ());
      return index;
    }
  }
  this->LeftChild (index) = HY_AVL_LEAF;
  this->RightChild (index) = HY_AVL_LEAF;
  this->BalanceFactor (index) = 0L;
  return index;
}

//______________________________________________________________________________
template <typename KEYTYPE>
void _AVLListBase<KEYTYPE>::_RemoveKey   (long index) {
  empty_slots.append (index);
}


//______________________________________________________________________________
template <typename KEYTYPE>
long _AVLListBase<KEYTYPE>::Insert(KEYTYPE const& key) {
  if (this->Length() > 0UL) { // something to do
    
    
    _SimpleList directions;
    
    long y = root,
    z = HY_AVL_LEAF,
    p,
    q,
    n,
    w;
    
    //bool go_right = false;
    long move_direction = HY_AVL_MOVE_LEFT;
    
    for (q = z, p = y; p != HY_AVL_LEAF;
         q = p, p = this->_MoveInTree (p, move_direction)) {
      
      //long comp = dataList->Compare(b, p);
      move_direction = this->_CompareIndexToValue (p, key);
      
      if (move_direction == 0L) {
        return -p - 1L;
      }
      
      if (this->BalanceFactor (p) != 0L) {
        z = q;
        y = p;
        directions.Clear();
      }
      directions.append (move_direction);
    }
    
    
    // insert new node
    
    n = this->_StoreKey (key);
    
    if (move_direction == HY_AVL_MOVE_RIGHT) {
      this->RightChild (q) = n;
    } else {
      this->LeftChild (q)  = n;
    }
    
    // update balance factors
    
    p = y;
    
    for (long k = 0L; p != n;
         p =  this->_MoveInTree (p, directions.AtIndex(k)), k++) {
      this->BalanceFactor (p) += (directions.AtIndex(k) == HY_AVL_MOVE_RIGHT ? 1L : -1L);
    }
    
    if ( this->BalanceFactor (y)  == -2L) {
      long x = this->LeftChild (y);
      if (this->BalanceFactor(x) == -1L) {
        w = x;
        LEFT_SHIFT (this->LeftChild (y),this->RightChild (x),y);
        this->BalanceFactor(x) = 0L;
        this->BalanceFactor(y) = 0L;
      } else {
        w = this->RightChild (x);
        LEFT_SHIFT (this->RightChild(x),this->LeftChild (w),x);
        LEFT_SHIFT (this->LeftChild(y),this->RightChild(w),y);
        
        switch (this->BalanceFactor(w)) {
          case -1L:
            this->BalanceFactor(x) = 0L;
            this->BalanceFactor(y) = 1L;
            this->BalanceFactor(w) = 0L;
            break;
          case 0L:
            this->BalanceFactor(x) = 0L;
            this->BalanceFactor(y) = 0L;
            break;
          default:
            this->BalanceFactor(x) = -1L;
            this->BalanceFactor(y) = 0L;
            this->BalanceFactor(w) = 0L;
            break;
        }
        
      }
    } else if (this->BalanceFactor (y) == 2L) {
      long x = this->RightChild (y);
      if (this->BalanceFactor(x) == 1L) {
        w = x;
        LEFT_SHIFT (this->RightChild (y), this->LeftChild(x),y);
        this->BalanceFactor(x) = 0L;
        this->BalanceFactor(y) = 0L;
      } else {
        w = this->LeftChild (x);
        LEFT_SHIFT (this->LeftChild(x),this->RightChild (w),x);
        LEFT_SHIFT (this->RightChild(y),this->LeftChild (w),y);
        switch (this->BalanceFactor(w)) {
          case 1L:
            this->BalanceFactor(x) = 0L;
            this->BalanceFactor(y) = -1L;
            this->BalanceFactor(w) = 0L;
            break;
          case 0L:
            this->BalanceFactor(x) = 0L;
            this->BalanceFactor(y) = 0L;
            break;
          default:
            this->BalanceFactor(x) = 1L;
            this->BalanceFactor(y) = 0L;
            this->BalanceFactor(w) = 0L;
            break;
        }
        
      }
    } else {
      return n;
    }
    
    if (z != HY_AVL_LEAF) {
      if (y == this->LeftChild (z)) {
        this->LeftChild (z) = w;
      } else {
        this->RightChild (z) = w;
      }
    }
    
    if (y == this->root) {
      this->root = w;
    }
    return n;
    // WAS return p;
  }
  
  
  return (this->root = this->_StoreKey (key, HY_LIST_INSERT_AT_END));
  
}


//______________________________________________________________________________
template <typename KEYTYPE>
void _AVLListBase<KEYTYPE>::_DeleteHelper (const long index, const long new_node, _SimpleList const &directions, _SimpleList const& nodes, bool update_root) {
  if (index > 1L) {
    (directions.AtIndex(index-1L) > 0L ?
      this->RightChild (nodes.AtIndex (index-1L)):
      this->LeftChild (nodes.AtIndex (index-1L))) = new_node;
  } else {
    if (update_root) {
      this->root = new_node;
    }
  }
}
//______________________________________________________________________________

template <typename KEYTYPE>
long _AVLListBase<KEYTYPE>::Delete(const KEYTYPE  &key) {
  
  if (this->Length() > 0UL) { // something to do
    
    _SimpleList directions,
                nodes;
    
    long p = this->root,
         cmp = this->_CompareIndexToValue(p, key),
         k = 1L;
    
    nodes.append (HY_AVL_LEAF);
    directions.append (HY_AVL_MOVE_RIGHT);
    
    for (; cmp != 0L; cmp = this->_CompareIndexToValue (p, key), k++) {
      
      nodes.append (p);
      directions.append (cmp);
      p = this->_MoveInTree (p, cmp);
      
      if (p == HY_AVL_LEAF) {
        return HY_NOT_FOUND;
      }
    }
    
    _RemoveKey(p);
    
    long r = this->RightChild (p); //rightChild.lData[p];
    
    if (r == HY_AVL_LEAF) {
      this->_DeleteHelper (k, this->LeftChild(p), directions, nodes, false);
      if (p == this->root) {
        this->root = this->LeftChild (this->root); //leftChild.lData[root];
      }
    } else {
      if (this->LeftChild (r)  == HY_AVL_LEAF) {
        this->LeftChild (r)     = this->LeftChild (p);      // leftChild.lData[r] = leftChild.lData[p];
        this->BalanceFactor (r) = this->BalanceFactor (p);   // balanceFactor.lData[r] = balanceFactor.lData[p];
        this->_DeleteHelper (k, r, directions, nodes);
        
        directions.append (HY_AVL_MOVE_RIGHT);
        nodes.append (r);
        k++;
      } else {
        
        long s;
        int j = k++;
        
        nodes.append (0L);
        directions.append (HY_AVL_MOVE_RIGHT);
        
        
        for (;;) {
          directions.append (HY_AVL_MOVE_LEFT);
          nodes.append (r);
          k++;
          s = this->LeftChild (r); //this->avl_structure.AtIndex(r).left_child;
          if (this->LeftChild (s) == HY_AVL_LEAF) {
            break;
          }
          r = s;
        }
 
        this->LeftChild (s) = this->LeftChild (p);
        LEFT_SHIFT (this->LeftChild (r), this->RightChild (s), this->RightChild (p));
        this->BalanceFactor (s) = this->BalanceFactor (p);
        this->_DeleteHelper (j, s, directions, nodes, false);
        
        nodes.SetItem(j,s);
        
        if (p == this->root) {
          this->root = s;
        }
      }
    }
    
    while (--k > 0L) {
      long y = nodes.AtIndex(k);
      if (directions.AtIndex(k) == HY_AVL_MOVE_LEFT) {
        switch ( (++ this->BalanceFactor (y)) ) {
          case 1L: {
            k = 0L; // break out of the main loop
            break;
          }
          case 2L: {
            long x = this->RightChild (y);
            if (this->BalanceFactor (x) == -1L) {
              long w = this->LeftChild (x);
              LEFT_SHIFT (this->LeftChild (x), this->RightChild (w), x);
              LEFT_SHIFT (this->RightChild (y), this->LeftChild (w), y);
              switch (this->BalanceFactor (w)) {
                case 1L:
                  this->BalanceFactor (x) = 0L;
                  this->BalanceFactor (y) = -1L;
                  this->BalanceFactor (w) = 0L;
                  break;
                case 0L:
                  this->BalanceFactor (x) = 0L;
                  this->BalanceFactor (y) = 0L;
                  break;
                default:
                  this->BalanceFactor (x) = 1L;
                  this->BalanceFactor (y) = 0L;
                  this->BalanceFactor (w) = 0L;
              }
              this->_DeleteHelper (k, w, directions, nodes);
            } else {
              LEFT_SHIFT (this->RightChild (y),this->LeftChild (x),y);
              this->_DeleteHelper (k, x, directions, nodes);
              
              if ( this->BalanceFactor(x) == 0L) {
                this->BalanceFactor(x) = -1L;
                this->BalanceFactor(y) = 1L;
                k = 0L; // break out of the main loop
                break;
              } else {
                this->BalanceFactor(x) = 0L;
                this->BalanceFactor(y) = 0L;
              }
            }
            break; // switch
          }
        } // end switch
      } else {
        switch ( (-- this->BalanceFactor (y)) ) {
          case -1L : {
            k = 0L;
            break;
          }
          case -2L : {
            long x = this->LeftChild (y);
            if (this->BalanceFactor (x) == 1L) {
              long w = this->RightChild (x);
              LEFT_SHIFT (this->RightChild (x), this->LeftChild (w), x);
              LEFT_SHIFT (this->LeftChild (y), this->RightChild (w), y);
              switch (this->BalanceFactor (w)) {
                case -1L:
                  this->BalanceFactor (x) = 0L;
                  this->BalanceFactor (y) = 1L;
                  this->BalanceFactor (w) = 0L;
                  break;
                case 0L:
                  this->BalanceFactor (x) = 0L;
                  this->BalanceFactor (y) = 0L;
                  break;
                default:
                  this->BalanceFactor (x) = -1L;
                  this->BalanceFactor (y) = 0L;
                  this->BalanceFactor (w) = 0L;
              }
              this->_DeleteHelper (k, w, directions, nodes);
            } else {
              LEFT_SHIFT (this->LeftChild (y), this->RightChild (x), y);
              this->_DeleteHelper (k, x, directions, nodes);
              
              if ( this->BalanceFactor(x) == 0L) {
                this->BalanceFactor(x) = 1L;
                this->BalanceFactor(y) = -1L;
                k = 0L; // break out of the main loop
                break;
              } else {
                this->BalanceFactor(x) = 0L;
                this->BalanceFactor(y) = 0L;
              }
            }
            break;
          }
        }
      }
    }
    
    return p;
  }
  return HY_NOT_FOUND;
  
}


// DEBUG FUNCTIONS

//______________________________________________________________

template <typename KEYTYPE>
void _AVLListBase<KEYTYPE>::EchoList(void) const {
  if (root == HY_AVL_LEAF) {
    printf ("Empty list\n");
  } else {
    _EchoListHelper (root, 0);
  }
}


template <typename KEYTYPE>
void _AVLListBase<KEYTYPE>::_EchoListHelper(long node, long offset) const {
  if (node != HY_AVL_LEAF) {
    for (long i = 0; i < offset; i++) {
      printf ("-");
    }
    printf ("[%ld]=>%ld (balance = %ld)\n", node, (long)*AtIndex (node), BalanceFactor (node));
    _EchoListHelper (LeftChild(node), offset+1);
    _EchoListHelper (RightChild(node), offset+1);
  }
}
  
  

#include <math.h>

//______________________________________________________________

template <typename KEYTYPE>
const char* _AVLListBase<KEYTYPE>::ConsistencyCheck(void) const {
  _SimpleList nodeStack ;
  
  long        curNode  = this->root,
  lastNode = HY_AVL_LEAF;
  
  while (1) {
    
    while (curNode != HY_AVL_LEAF) {
      
      nodeStack << curNode;
      
      curNode = this->_MoveInTree (curNode, HY_AVL_MOVE_LEFT);
      if (curNode >= (long)this->avl_structure.Length()) {
        return "Index out of bounds";
      }
      
    }
    
    if (long h = nodeStack.Length()) {
      if (h>3*log (1.+this->Length())) {
        return "Unbalanced tree";
      }
      
      h--;
      
      curNode = nodeStack.AtIndex (h);
      
      if (lastNode != HY_AVL_LEAF && curNode != HY_AVL_LEAF) {
        
        if (_CompareIndexToValue (lastNode, *this->AtIndex(curNode)) <= 0) {
          return "Key comparisons contradtict tree structure";
        }
      }
      if (this->BalanceFactor (curNode) < -1 || this->BalanceFactor (curNode) > 1) {
        return "Balance factor is out of bounds";
      }
      
      lastNode = curNode;
      
      curNode = this->_MoveInTree (curNode, HY_AVL_MOVE_RIGHT);
      if (curNode >= (long)this->avl_structure.Length()) {
        return "Index out of bounds";
      }
      
      nodeStack.Pop();
    } else {
      break;
    }
  }
  
  return NULL;
}


#include "hy_avllistbase.cpp"


#endif
