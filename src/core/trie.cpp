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

#include "trie.h"
#include "hy_string_buffer.h"
#include "helperfunctions.h"

//______________________________________________________________________________
_Trie::_Trie(_String* alph) {
  this->linear_list.Initialize();
  this->setAlphabet(alph, true);
  this->linear_list.AppendNewInstance(new _hyListNumeric<long>);
  this->payload.append(0L);
  this->parents.append(-1L);
}

//______________________________________________________________________________
_Trie::~_Trie(void) {}


//______________________________________________________________________________
void _Trie::clear(bool all) {

  this->linear_list.Clear(all);
  payload.Clear(all);
  this->linear_list.AppendNewInstance(new _hyListNumeric<long>);
  payload << 0L;
  parents << -1L;

}

//______________________________________________________________________________
_hyList<char> _Trie::alphabet(void) const {

  _hyList<char> result;

  for (unsigned long charIndex = 0; charIndex < 256; charIndex++) {
    if (char_map.Element(charIndex) >= 0) {
      result.append((char)charIndex);
    }
  }

  return result;

}

//______________________________________________________________________________
void _Trie::setAlphabet(const _String* alph, bool do_clear) {

  if (do_clear) {
    clear(true);
    char_map.Clear();
  }

  if (alph) {
    this->char_map.Populate(257, -1, 0);
    unsigned long char_counter = 0;

    // always allow the '\0' character
    this->char_map[0] = 1; 

    for (unsigned long charIndex = 0; charIndex < alph->Length();
         charIndex++) {
      this->char_map[(unsigned char) alph->getChar(charIndex)] = 1;
    }

    // now sort alphically
    for (unsigned long charIndex = 0; charIndex < 256; charIndex++) {
      if (this->char_map[charIndex] == 2)
        this->char_map[charIndex] = char_counter++;
    }

  } else {
    this->char_map.Populate(256, 0, 1);
  }

}

////______________________________________________________________________________
//_Trie::_Trie(_Trie & source) {
//  Duplicate(&source);
//}

//______________________________________________________________________________
unsigned long _Trie::countitems(void) const {
  return this->Length();
}

//______________________________________________________________________________
unsigned long _Trie::Length(void) const {
  return payload.countitems() - 1L;
}

//______________________________________________________________________________
long _Trie::FindNextLetter(const char letter,
                           const unsigned long current_index) const {

  long letterKey = this->char_map.Element((const unsigned char) letter);

  if (letterKey >= 0) {

    _hyListNumeric<long> *thisList = this->linear_list.Element(current_index);
    letterKey = thisList->FindStepping(letterKey, 2, 0);

    if (letterKey < 0) {
      return HY_TRIE_NOTFOUND;
    }

    return thisList->Element(letterKey + 1);

  }

  return HY_TRIE_INVALID_LETTER;

}

//______________________________________________________________________________
long _Trie::FindNextUnusedIndex(bool alloc) {

  this->payload << 0;
  this->parents << 0;
  this->linear_list.AppendNewInstance(new _hyListNumeric<long>);
  return this->linear_list.Length() - 1;

}

//______________________________________________________________________________
long _Trie::InsertNextLetter(const char letter,
                             const unsigned long current_index) {

  long letter_key = char_map.Element((const unsigned char) letter);

  if (letter_key >= 0) {
    long next_index = FindNextUnusedIndex(letter != 0);
    _hyListNumeric<long> *currentList = this->linear_list.Element(current_index);
    (*currentList) << letter_key;
    (*currentList) << next_index;
    parents[next_index] = current_index;
    return next_index;
  }

  return HY_TRIE_INVALID_LETTER;
}

//______________________________________________________________________________
long _Trie::Find(const _String key, _hyListNumeric<long> *path, bool prefixOK) const {

  long current_index = 0, next_index = 0;

  for (long k = 0; k <= key.Length() && current_index >= 0; k++) {

    next_index = FindNextLetter(key[k], current_index);

    if (path) {
      path->append(next_index);
    }

    if (next_index < 0 && prefixOK) {
      next_index = FindNextLetter(0, current_index);
      current_index = next_index;
      break;
    }

    current_index = next_index;

  }

  return current_index;

}

//______________________________________________________________________________
long _Trie::Find(const char key, bool prefixOK) const {

  long current_index = 0, next_index = FindNextLetter(key, current_index);
  if (next_index < 0 && prefixOK) {
    next_index = FindNextLetter(0, current_index);
  }
  current_index = next_index;
  return current_index;
}

//______________________________________________________________________________
long _Trie::GetValueFromString(const _String &key) {
  long keyIndex = Find(key);
  if (keyIndex != HY_TRIE_NOTFOUND) {
    return GetValue(keyIndex);
  }
  return HY_TRIE_NOTFOUND;
}

//______________________________________________________________________________
void _Trie::UpdateValue(const long key, const long value) {
  if (key >= 0 && key < this->payload.Length())
    this->payload.SetItem(key, value);
}

//______________________________________________________________________________

long _Trie::GetValue(const long key) const {
  if (key >= 0 && key < payload.Length())
    return payload.Element(key);
  return 0L;
}

//______________________________________________________________________________

long _Trie::Insert(const char *key, const long value, bool return_index) {
  _String key_string(key);
  long ret_value = this->Insert(key_string, value);
  if (ret_value >= 0 && return_index == false) {
    return key_string.Length();
  }
  return ret_value;
}

//______________________________________________________________________________

long _Trie::Insert(const _String &key, const long value) {

  // the root is always at index 0
  long current_index = 0, current_char = 0,
       next_index = this->FindNextLetter(key[current_char++], current_index);

  while (next_index >= 0 && current_char <= key.Length()) {
    current_index = next_index;
    next_index = this->FindNextLetter(key[current_char++], current_index);
  }

  if (next_index == HY_TRIE_INVALID_LETTER) {
    return HY_TRIE_INVALID_LETTER;
  }

  if (current_char == key.Length() && next_index >= 0) {
    return next_index;
  }

  current_char--;

  // validate the rest of the string

  for (long k = current_char; k <= key.Length(); k++) {

    if (char_map[key[k]] < 0) {
      return HY_TRIE_INVALID_LETTER;
    }

  }

  // insert the rest of the string
  for (; current_char <= key.Length(); current_char++) {
    //printf ("\nInserting %c\n", key[current_char]);
    current_index = this->InsertNextLetter(key[current_char], current_index);
  }

  this->UpdateValue(current_index, value);

  return current_index;
}

//______________________________________________________________________________

unsigned long _Trie::Insert(const _hyListReference<_String> key, const _hyListNumeric<long> values) {

  unsigned long how_many = 0;

  for (long k = 0; k < key.Length(); k++) {

    _String serialized_key(*key.Element(k));
    long val =  values.Length() ? values.Element(k) : 0;

    long this_index = this->Insert(serialized_key, val);

    if (this_index >= 0) {
      how_many++;
    }

  }

  return how_many;

}

//______________________________________________________________________________
bool _Trie::Delete(_String key) {

  _hyListNumeric<long> history;

  long found_key = this->Find(key, &history);
  
  if (found_key >= 0L) {
    // now traverse the history list backwards and delete all keys that have no
    // children

    for (long k = history.Length() - 1; k >= 0; k--) {
      _hyListNumeric<long> *current_list = this->linear_list.Element(history[k]);

      if (current_list == nil || current_list->Length() <= 1L) {
        this->payload[history[k]] = 0L;
        this->parents[history[k]] = -1L;

        _hyListNumeric<long> *parentList = this->linear_list.Element(history[k - 1]);

        unsigned long parentNode =
            parentList->FindStepping(history[k], 2, 1) - 1;

        parentList->Delete(parentNode);
        parentList->Delete(parentNode);

      }

    }

    // SW 20140811: Breaks deconstructor
    //for (long k = history.Length() - 1; k >= 0; k--) {
    //  std::cout << "Setting null " << k << std::endl;
    //  this->linear_list.SetItem(history[k], NULL);
    //}

    return true;

  }

  return false;

}

//______________________________________________________________________________
bool _Trie::Delete(const char *key) {
  _String string_key(key);
  return Delete(string_key);
}

//______________________________________________________________________________
unsigned long _Trie::Delete(const _hyListReference<_String> &key) {

  unsigned long how_many = 0;
  for (long k = 0; k < key.Length(); k++) {
    _String serializedKey(*key.Element(k));

    long this_index = Delete(serializedKey);

    if (this_index >= 0) {
      how_many++;
    }
  }

  return how_many;

}


//______________________________________________________________________________

_StringBuffer _Trie::RetrieveStringFromPath(_hyListNumeric<long> path,
                                            _hyList<char> alphabet) {

  _StringBuffer this_string(256UL);
  _hyList<char> my_alph = alphabet.Length() ? alphabet : this->alphabet();

  for (long k = 0; k < path.Length() - 4; k += 2) {
    _hyListNumeric<long> current_list = *this->linear_list.Element(path.Element(k));
    long current_position = path.Element(k + 1);
    this_string << my_alph[current_list[current_position]];
  }

  //if (!alphabet) {
  //  this->DeleteObject(my_alph);
  //}

  return this_string;

}


//______________________________________________________________________________
void _Trie::DumpRaw() {

  for (long k = 0; k < this->linear_list.Length(); k++) {

      printf("Position %ld:\n", k);
      _hyListNumeric<long> *this_list = this->linear_list.Element(k);
      for (long m = 0; m < this_list->Length(); m += 2) {
        printf("'%c'(%ld) -> %ld\n", (char) this_list->Element(m),
               this_list->Element(m), this_list->Element(m + 1));
      }

  }

}

//______________________________________________________________________________

_StringBuffer _Trie::toStr() {

  _StringBuffer serialized = _StringBuffer(128UL); 
  _hyList<char> alph = this->alphabet();

  // 2 indices per entry: node and current position (in multiples of 2)
  _hyListNumeric<long> traversal_history,
                      root_list = *this->linear_list.Element(0L);

  traversal_history.append(0L);
  traversal_history.append(0L);
  bool doComma = false;

  serialized << '{';

  while (!(traversal_history.Length() == 2 &&
           traversal_history[1] == root_list.Length())) {

    _hyListNumeric<long> current_list = _hyListNumeric<long>(*this->linear_list[
        traversal_history[traversal_history.Length() - 2]]);

    long current_position =
        traversal_history[traversal_history.Length() - 1];


    // if current list is empty, then generate a string based on the path, and
    // advance up the chain
    if (current_list.Length()) {

      if (current_position < current_list.Length()) {

        traversal_history.append(current_list[current_position + 1]);
        traversal_history.append(0L);

      } else {

        traversal_history.Pop();
        traversal_history.Pop();

        // advance the counter in the parent
        long to_set = traversal_history[traversal_history.Length() - 1] + 2; 
        traversal_history.SetItem(traversal_history.Length() - 1, to_set); 

      }

    } else {

      _String this_string = this->RetrieveStringFromPath(traversal_history, alph);
      serialized << '"';
      serialized << this_string;
      serialized << "\":";

      serialized << _String(this->GetValue(
                   traversal_history[
                   traversal_history.Length() - 2]));

      if (doComma) {
        serialized << ',';
      } else {
        doComma = true;
      }

      serialized << '\n';
      traversal_history.Pop();
      traversal_history.Pop();

      // advance the counter in the parent
      long to_set = traversal_history[traversal_history.Length() - 1] + 2; 
      traversal_history.SetItem(traversal_history.Length() - 1, to_set); 

    }
  }

  serialized << '}';
  return serialized;

}

//______________________________________________________________________________

_String _Trie::RetrieveKeyByPayload(const long key) {

  long key_index = payload.Find(key);
  if (key_index >= 0) {

    _hyListNumeric<long> parent_indices, traversal_history;
    long keyer = key_index;

    do {
      parent_indices << keyer;
      keyer = parents[keyer];
    } while (keyer > 0);

    parent_indices << 0;
    parent_indices.Flip();

    for (long i = 0; i < parent_indices.Length() - 1; i++) {
      traversal_history << parent_indices[i];
      traversal_history << this->linear_list.Element(parent_indices[i])->FindStepping(parent_indices[i + 1], 2, 1) - 1;
    }

    traversal_history << key_index;
    traversal_history << 0L;
    _hyList<char> alph = this->alphabet();
    return _String(RetrieveStringFromPath(traversal_history, alph));

  }

  return nil;

}
