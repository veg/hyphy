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
