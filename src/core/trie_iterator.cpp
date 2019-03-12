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

#include "trie_iterator.h"

//#define  TRIE_ITERATOR_ENDINDEX (-1L)

TrieIterator::TrieIterator (_Trie const * trie): traversal_history (), root_list(NULL) {
  this->the_trie = trie;
}

TrieIterator::TrieIterator (TrieIterator const& copy_from) : traversal_history(copy_from.traversal_history), root_list(copy_from.root_list) {
  this->the_trie = copy_from.the_trie;
}

TrieIterator& TrieIterator::begin (void) {
    traversal_history.Clear();
    this->root_list = (_SimpleList*)the_trie->GetItem (0L);
    traversal_history << 0L << -2L << 0L << 0L;
    if (!Next()) {
        end();
    }
    return *this;
}

TrieIterator& TrieIterator::end (void) {
  traversal_history.Clear();
  this->root_list = NULL;
  return *this;
}

bool TrieIterator::Next (void) {

    traversal_history.Pop();
    traversal_history.Pop();
    if (traversal_history.countitems() >= 2) {
        traversal_history[traversal_history.lLength-1] += 2; // advance the counter in the parent
    }
    
    while (!(traversal_history.countitems() == 2 && traversal_history.get(1) == root_list->countitems())) {
        _SimpleList* current_list = (_SimpleList*)the_trie->GetItem(traversal_history.Element(-2L));
        long current_position = traversal_history.Element (-1L);
        // if current list is empty, then generate a string based on the path, and advance up the chain
        if (current_list && current_list->nonempty()) {
            if (current_position < current_list->countitems()) {
                traversal_history << current_list->get(current_position+1);
                traversal_history << 0;
            } else {
                traversal_history.Pop();
                traversal_history.Pop();
                traversal_history[traversal_history.lLength-1] += 2; // advance the counter in the parent
            }
        } else {
            return true;
            //traversal_history.Pop();
            //traversal_history.Pop();
            //traversal_history[traversal_history.lLength-1] += 2; // advance the counter in the parent
        }
    }
    return false;
}

TrieIterator& TrieIterator::operator++ (void) {
  if (!Next()) {
    this->end();
  }
  return *this;
}

TrieIterator& TrieIterator::operator =  (const TrieIterator& assign_from) {
  if (this != &assign_from) {
    this->the_trie = assign_from.the_trie;
    this->root_list = assign_from.root_list;
  }
  return *this;
}

TrieIteratorKeyValue const  TrieIterator::operator* (void) {
  if (this->root_list) {
    _String alph (the_trie->Alphabet());
    return TrieIteratorKeyValue (the_trie->RetrieveStringFromPath(traversal_history, &alph),traversal_history.Element (-2));
  }
  return TrieIteratorKeyValue (nil,0);
}

bool TrieIterator::operator == (TrieIterator const & compare) {
  return this->the_trie == compare.the_trie && this->traversal_history.Equal (compare.traversal_history);
  
}

bool TrieIterator::operator != (TrieIterator const & compare) {
    return this->the_trie != compare.the_trie || !this->traversal_history.Equal (compare.traversal_history);
}
