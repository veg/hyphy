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

#define HY_AVL_LEAF (-1L)

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
};


/**
 * This object implements a key:value dictionary functionality using 
 * a balanced binary search tree.
 * @param KEYTYPE is the type used to represent keys; it must support stack copy, all pairwise comparison operations, i.e. ==, <= etc
 * @param PAYLOAD is the type which represents dictionary values; it must support stack copy and assignment (very little else is assumed)
 */


//_____________________________________________________________________________
template <typename KEYTYPE, typename PAYLOAD>
class _AVLList : public virtual BaseObj{
    
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
    
    _hyListOrderable    <KEYTYPE>        keys;
    
    /**
     * This list stores the payloads associated with each key
     * they are tied by index, i.e. the key stored in the i-th index
     * of avl_structure is bound to the i-th value in this list
     */

    _hyList             <PAYLOAD>        values;
    
    
    /**
     * Stores indices of deleted nodes in avl_structure,
     * so that they can be used by subsequent insert operations
     */

    _SimpleList         empty_slots;

    
    /**
     * Stores the index of the root; if HY_AVL_LEAF, then the list is empty
    */

    long root;

    
    //_SimpleList *dataList, emptySlots;
    
    /**
     * this is a convenience function which moves from the current node in the tree
     * based on the result of a key comparison
     * @param node the index in the tree; no range checking is performed
     * @param direction move left if < 0 (less than), right if > 0 (greater then)
     * @return the index of the node moved to
     */

    inline long _MoveInTree (long node, long direction) const ;
   
    /**
     * Compares the key at a given index to a fixed value
     * @param node The index of the item to compare
     * @param key The value to compare the index to
     * @return -1 if keys[node] < key , 0 if keys[node] == key, or 1 if keys[node] > key
     */
    
    virtual long _CompareIndexToValue (long node, KEYTYPE const & key) const ;
   
    
public:


    /**
     * Construct an empty AVLList
     */
    
  _AVLList     (void);

    /**
     * Stack copy constructor
     */
    
  _AVLList     (_AVLList <KEYTYPE, PAYLOAD> const & );

    /**
     * Object destuctor
     */
    
  virtual ~_AVLList(void);

    /**
     * The basic search function
     * @param key The value to find
     * @param hist if not null, then the indices traversed on the way to this key will be stored
     * @return the >=0 index of the key in this list or HY_NOT_FOUND if not found
     */
  
  long   Find(const KEYTYPE key, _SimpleList* history = NULL) const;

    /**
     * Length of the list
     * @return The number of elements in the list
     */
    
  unsigned long Length (void) const;
    
  /*virtual void Clear(bool = false);
  virtual bool HasData(long);

  virtual void ReorderList(_SimpleList * = nil);
  virtual long InsertData(BaseRef, long, bool);
  virtual BaseRef toStr(void);
  virtual long Traverser(_SimpleList &, long &, long = -1);
  virtual long GetRoot(void) { return root; }
  virtual void DeleteXtra(long); // {}

 virtual void DeleteAll(bool cL);
    //{
    //Clear(cL);
    //DeleteObject(dataList);
    //}


  // 20100623: a shortcut function to look for integers only
  // avoids calling ::Compare

  char FindBest(BaseRef, long &);

  long Next(long, _SimpleList &);
  long Prev(long, _SimpleList &);
  long First(void);
  long Last(void);

  long GetByIndex(const long);

  // the 1st bool flag is to say whether to dup the object being inserted
  // the 2nd bool flag (if the first flag is false) if set to true,
  // will cause failed inserts (key already exists) to delete the key
  long Insert(BaseRef, long = 0, bool = true, bool = false);

  BaseRef Retrieve(long);

  void Delete(BaseRef, bool = false);
  void ConsistencyCheck(void);*/

};

#include "avllist.cpp"


#endif
