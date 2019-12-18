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

#include "avllistx_iterator.h"

#define  AVL_LISTX_ITERATOR_ENDINDEX (-1L)

AVLListXIterator::AVLListXIterator (_AVLListX const * the_list):
  current_location (AVL_LISTX_ITERATOR_ENDINDEX),
  traversal_lookup (AVL_LISTX_ITERATOR_ENDINDEX),
  traversal_history () {
  this->the_list = the_list;
}

AVLListXIterator::AVLListXIterator (AVLListXIterator const& copy_from) :
  traversal_history(copy_from.traversal_history) {
  this->the_list = copy_from.the_list;
  this->current_location = copy_from.current_location;
  this->traversal_lookup = copy_from.traversal_lookup;
}

AVLListXIterator& AVLListXIterator::begin (void) {
  this->current_location = the_list->Traverser (traversal_history, traversal_lookup, the_list->GetRoot());
  return *this;
}

AVLListXIterator& AVLListXIterator::end (void) {
  this->current_location = AVL_LISTX_ITERATOR_ENDINDEX;
  return *this;
}

AVLListXIterator& AVLListXIterator::operator++ (void) {
  this->current_location = the_list->Traverser (traversal_history, traversal_lookup);
  return *this;
}

AVLListXIterator& AVLListXIterator::operator =  (const AVLListXIterator& assign_from) {
  if (this != &assign_from) {
    this->the_list = assign_from.the_list;
    this->current_location = assign_from.current_location;
    this->traversal_lookup = assign_from.traversal_lookup;
    this->traversal_history = assign_from.traversal_history;
  }
  return *this;
}

AVLListXIteratorKeyValue const  AVLListXIterator::operator* (void) {
  if (this->current_location >= 0) {
    return AVLListXIteratorKeyValue (this->current_location, the_list->GetXtra(this->current_location));
  }
  return AVLListXIteratorKeyValue (AVL_LISTX_ITERATOR_ENDINDEX,0);
}

bool AVLListXIterator::operator == (AVLListXIterator const & compare) {
  return this->the_list == compare.the_list && this->current_location == compare.current_location;
  
}

bool AVLListXIterator::operator != (AVLListXIterator const & compare) {
  return this->the_list != compare.the_list || this->current_location != compare.current_location;
  
}
