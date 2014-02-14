/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-now
Core Developers:
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon    (apoon@cfenet.ubc.ca)
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

#include "trieiterator.h"

_TrieIterator::_TrieIterator (const _Trie& s) {
  data_source = &s;
  alphabet = s.Alphabet();
  Init ();
}

void _TrieIterator::Init (bool last) {
  traversal_history.Clear();
  traversal_history << 0L;
  if (last) {
    _SimpleList *current_traversal_list = (_SimpleList*)data_source->GetElement(0L);
    traversal_history << current_traversal_list->lLength-2;
  } else {
    traversal_history << 0L;
  }
  last_retrieved_index = HY_NOT_FOUND;
}

_String* _TrieIterator::First (void) {
  Init();
  if (MoveDown() != HY_NOT_FOUND) {
    return BuildPath();
  }
  return NULL;
}

_String* _TrieIterator::Last (void) {
  Init(true);
  if (MoveDown(true) != HY_NOT_FOUND) {
    return BuildPath();
  }
  return NULL;
}

_String* _TrieIterator::Next (void) {
  if (RollUp (false) != HY_NOT_FOUND && MoveDown(true) != HY_NOT_FOUND) {
    return BuildPath();
  }
  return NULL;
}

_String* _TrieIterator::Previous (void) {
  if (RollUp (true) != HY_NOT_FOUND && MoveDown(false) != HY_NOT_FOUND) {
    return BuildPath();
  }
  return NULL;
}

long     _TrieIterator::MoveDown (bool move_left) {
  if (traversal_history.lLength) {
    
    _SimpleList *current_traversal_list = (_SimpleList*)data_source->GetElement(traversal_history.GetElement(-2L));
    long         current_child_index    = traversal_history.GetElement(-1L);
    
    while (current_traversal_list && current_child_index < current_traversal_list->lLength) {
      long next_index = current_traversal_list->GetElement(current_child_index+1);
      current_traversal_list = (_SimpleList*)data_source->GetElement(next_index);
      if (current_traversal_list && current_traversal_list->lLength) {
        traversal_history << next_index;
        current_child_index = move_left? 0L: current_traversal_list->lLength - 2;
        traversal_history << current_child_index;
      } else {
        traversal_history << next_index;
        traversal_history << 0L;
        return next_index;
      }
    }
  }
  return HY_NOT_FOUND;
}

long     _TrieIterator::RollUp (bool move_left) {
  if (traversal_history.lLength > 2L) {
    _SimpleList *current_traversal_list = (_SimpleList*)data_source->GetElement(traversal_history.GetElement(-2L));
    if (current_traversal_list) {
      return HY_NOT_FOUND;
    }
    while (traversal_history.lLength > 2L) {
      
      traversal_history.Pop();
      traversal_history.Pop();
      current_traversal_list = (_SimpleList*)data_source->GetElement(traversal_history.GetElement(-2L));
      
      if (move_left) {
        if (traversal_history.GetElement (-1L) >= 2L) {
          traversal_history.lData[traversal_history.lLength-1] -= 2L;
          return traversal_history.GetElement (-2L);
        }
      } else {
        if (traversal_history.GetElement (-1L) < current_traversal_list->lLength - 2L) {
          traversal_history.lData[traversal_history.lLength-1] += 2L;
          return traversal_history.GetElement (-2L);
        }
      }
    }
    
  }
  return HY_NOT_FOUND;
}

_String*     _TrieIterator::BuildPath (void){
  unsigned long pathL = traversal_history.lLength/2;
  _String * res = new _String (pathL, true);
  
  for (unsigned long k = 0L; k < traversal_history.lLength-2L; k+=2L) {
    _SimpleList *current_list = (_SimpleList*)data_source->GetElement(traversal_history.lData[k]);
    long current_position = traversal_history.lData[k + 1];
    (*res) << alphabet.sData[current_list->lData[current_position]];
  }
  
  last_retrieved_index = traversal_history.GetElement(-2L);
  return res;
}