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

_TrieIterator::_TrieIterator(const _Trie &s) {

  this->data_source = &s;
  this->alphabet = s.alphabet();
  this->Init();

}

void _TrieIterator::Init(bool last) {

  this->traversal_history.Clear();
  this->traversal_history << 0L;

  if (last) {
    _hyListNumeric<long> *current_traversal_list =
        (_hyListNumeric<long> *)this->data_source->linear_list.Element(0L);

    this->traversal_history << current_traversal_list->Length() - 2;

  } else {
    this->traversal_history << 0L;
  }

  last_retrieved_index = HY_NOT_FOUND;

}

_String *_TrieIterator::First(void) {
  this->Init();
  if (MoveDown() != HY_NOT_FOUND) {
    return this->BuildPath();
  }
  return NULL;
}

_String *_TrieIterator::Last(void) {
  this->Init(true);
  if (MoveDown(true) != HY_NOT_FOUND) {
    return this->BuildPath();
  }
  return NULL;
}

_String *_TrieIterator::Next(void) {
  if (this->RollUp(false) != HY_NOT_FOUND && this->MoveDown(true) != HY_NOT_FOUND) {
    return this->BuildPath();
  }
  return NULL;
}

_String *_TrieIterator::Previous(void) {
  if (this->RollUp(true) != HY_NOT_FOUND && this->MoveDown(false) != HY_NOT_FOUND) {
    return this->BuildPath();
  }
  return NULL;
}

long _TrieIterator::MoveDown(bool move_left) {

  if ( this->traversal_history.Length() ) {

    _hyListNumeric<long> *current_traversal_list =
        (_hyListNumeric<long> *)this->data_source->linear_list.Element(
            this->traversal_history.Element(-2L));

    long current_child_index = this->traversal_history.Element(-1L);

    while (current_traversal_list &&
           current_child_index < current_traversal_list->Length()) {

      long next_index =
          current_traversal_list->Element(current_child_index + 1);

      current_traversal_list =
          (_hyListNumeric<long> *)this->data_source->linear_list.Element(next_index);

      if (current_traversal_list && current_traversal_list->Length()) {
        this->traversal_history << next_index;
        current_child_index =
            move_left ? 0L : current_traversal_list->Length() - 2;
        this->traversal_history << current_child_index;

      } else {
        this->traversal_history << next_index;
        this->traversal_history << 0L;
        return next_index;
      }
    }
  }

  return HY_NOT_FOUND;

}

long _TrieIterator::RollUp(bool move_left) {

  if (this->traversal_history.Length() > 2L) {
    _hyListNumeric<long> *current_traversal_list =
        (_hyListNumeric<long> *)this->data_source->linear_list.Element(
            this->traversal_history.Element(-2L));

    if (current_traversal_list->Length()) {
      return HY_NOT_FOUND;
    }

    while (this->traversal_history.Length() > 2L) {

      this->traversal_history.Pop();
      this->traversal_history.Pop();

      current_traversal_list = (_hyListNumeric<long> *)this->data_source->linear_list.Element(
          this->traversal_history.Element(-2L));

      if (move_left) {
        if (this->traversal_history.Element(-1L) >= 2L) {
          this->traversal_history[this->traversal_history.Length() - 1] -= 2L;
          return this->traversal_history.Element(-2L);
        }

      } else {
        if (this->traversal_history.Element(-1L) <
            current_traversal_list->Length() - 2L) {

          this->traversal_history[this->traversal_history.Length() - 1] += 2L;
          return this->traversal_history.Element(-2L);

        }
      }
    }
  }

  return HY_NOT_FOUND;

}

_StringBuffer* _TrieIterator::BuildPath(void) {

  unsigned long pathL = this->traversal_history.Length() / 2;
  _StringBuffer *res = new _StringBuffer(pathL);

  for (unsigned long k = 0L; k < this->traversal_history.Length() - 2L; k += 2L) {
    _hyListNumeric<long> *current_list =
        (_hyListNumeric<long> *)this->data_source->linear_list.Element(this->traversal_history[k]);
    long current_position = this->traversal_history[k + 1];
    (*res) << (char)this->alphabet[current_list->Element(current_position)];
  }

  last_retrieved_index = this->traversal_history.Element(-2L);
  return res;

}


