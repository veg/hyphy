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

#ifndef _AVLLISTBASE_
#define _AVLLISTBASE_
//#pragma once
#include <hy_avllistbase.h>


/**
 * This object implements a key:value dictionary functionality using
 * a balanced binary search tree.
 * @param KEYTYPE is the type used to represent keys; it must support stack copy, all pairwise comparison operations, i.e. ==, <= etc
 * @param PAYLOAD is the type which represents dictionary values; it must support stack copy and assignment (very little else is assumed)
 */


//_____________________________________________________________________________
template <typename KEYTYPE>
class _AVLList : public virtual _AVLListBase<KEYTYPE>{
  
protected:
  
  /**
   * This list stores the keys present in the list
   * they are tied to tre tree nodes by index, i.e. the key stored in the i-th index
   */
  
   _hyListOrderable    <KEYTYPE>        keys;
    virtual long _CompareIndexToValue (long node, KEYTYPE const & key) const ;
  
  
  /**
   * Stores a key in the list at a given index
   * Must be implemented by a derived class
   * @param key The key to store (assumed to NOT already be in the tree, otherwise things will break)
   * @param index The index to store at (pass HY_LIST_INSERT_AT_END to find an unused index to store)
   * @return -1 if keys[node] < key , 0 if keys[node] == key, or 1 if keys[node] > key
   */
  
   virtual long _StoreKey   (const KEYTYPE& key, long index = HY_LIST_INSERT_AT_END);
  
  /**
   * Remove the key from 'index' (note that the tree structure is handled by the ::Delete function)
   * @param index The index to delete from (assumed to be valid)
   */
  virtual void _RemoveKey   (long index);
  virtual KEYTYPE const * _AtIndex (const unsigned long) const;
  
  
public:
  
  
  /**
   * Construct an empty AVLList
   */
  
  _AVLList     (void);
  
  /**
   * Stack copy constructor
   */
  
  _AVLList     (_AVLList <KEYTYPE> const & );

  /**
   * Single item constructor
   */
  
  _AVLList     (KEYTYPE const & );
  
  /**
   * List of items constructor
   */
  
  _AVLList     (_hyList <KEYTYPE> const &);

  
  /**
   * Object destuctor
   */
  
  virtual ~_AVLList(void);
  
  /**
   * Initialize a freshly created object
   */
  void Initialize (bool = false);
  
  /**
   * Clone an existing object into this object, which has been previously initilized
   * @param source the object to clone from
   */
  void Clone (_AVLList <KEYTYPE> const & source);

  /* Return the key at a given index
   * @param key the index
   * @return a pointer to the key, or NULL if the index is invalid of is poitning at a vacated
   * space
   */
  
  virtual BaseObj *makeDynamic(void) const;
  virtual void Duplicate(BaseObj const * ref);
  
  
  
};

//#include "hy_avllist.cpp"

//*************** CONSTRUCTORS ***************//


template <typename KEYTYPE>
_AVLList<KEYTYPE>::_AVLList(void):_AVLListBase<KEYTYPE> () {
  this->Initialize();
}

template <typename KEYTYPE>
_AVLList<KEYTYPE>::_AVLList(_AVLList<KEYTYPE> const & source):_AVLListBase<KEYTYPE> (source) {
  this->Clone (source);
}

template <typename KEYTYPE>
_AVLList<KEYTYPE>::_AVLList(KEYTYPE const & item) {
  this->Initialize();
  this->Insert (item);
}

template <typename KEYTYPE>
_AVLList<KEYTYPE>::_AVLList(_hyList<KEYTYPE> const & items) {
  this->Initialize();
  for (unsigned long k = 0UL; k < items.Length(); k+=1) {
    this->Insert (items.AtIndex (k));
  }
}


//*************** INITIALIZER and CLONER ***************//
template <typename KEYTYPE>
void _AVLList<KEYTYPE>::Clone(_AVLList<KEYTYPE> const & source) {
  this->_AVLListBase<KEYTYPE>::Clone (source);
  this->keys.Clone (&source.keys);
}

template <typename KEYTYPE>
void _AVLList<KEYTYPE>::Initialize(bool) {
}


//*************** DESTRUCTOR ***************//

template <typename KEYTYPE>
_AVLList<KEYTYPE>::~_AVLList(void) {
}

//*************** REQUIRED FUNCTION DEFINITIONS ***************//

template <typename KEYTYPE>
long _AVLList<KEYTYPE>::_CompareIndexToValue(long node, KEYTYPE const &key ) const {
  return -this->keys.CompareToValue (node, key);
}

template <typename KEYTYPE>
long _AVLList<KEYTYPE>::_StoreKey(KEYTYPE const &key, long index) {
  index = this->_AVLListBase <KEYTYPE>::_StoreKey (key, index);
  this->keys.append_or_insert (key, index);
  return index;
}

template <typename KEYTYPE>
void _AVLList<KEYTYPE>::_RemoveKey(long index) {
  this->_AVLListBase <KEYTYPE>::_RemoveKey (index);
  //this->keys.Delete (index);
}

template <typename KEYTYPE>
KEYTYPE const * _AVLList<KEYTYPE>::_AtIndex(unsigned long index) const {
  return &this->keys.AtIndex (index);
}

template <typename KEYTYPE>
BaseObj* _AVLList<KEYTYPE>::makeDynamic(void) const {
  return new _AVLList <KEYTYPE> (*this);
}

template <typename KEYTYPE>
void _AVLList<KEYTYPE>::Duplicate (BaseObj const * ref) {
  this->Clone (*dynamic_cast<_AVLList<KEYTYPE> const *> (ref));
}

#endif
