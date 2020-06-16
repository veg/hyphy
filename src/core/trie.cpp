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

#include "trie.h"
#include "global_things.h"
#include "hy_string_buffer.h"

//----------------------------------------------------------------------------------------------------------------------

_Trie::_Trie (const _String* alphabet) {
    InitializeTrie(alphabet);
}

//----------------------------------------------------------------------------------------------------------------------

void _Trie::InitializeTrie (const _String* alphabet) {
    SetAlphabet (alphabet, false);
    AppendNewInstance(new _SimpleList);
    payload << 0L;
    parents <<-1L;
    inserted_values = 0UL;
}


//----------------------------------------------------------------------------------------------------------------------

_Trie::~_Trie (void){
    
}
       
//----------------------------------------------------------------------------------------------------------------------

void _Trie::Clear (bool all){
    _List::Clear (all);
    payload.Clear(all);
    emptySlots.Clear(all);
    AppendNewInstance(new _SimpleList);
    inserted_values = 0UL;
    payload << 0L;
    parents <<-1L;
}

//----------------------------------------------------------------------------------------------------------------------

_String _Trie::Alphabet (void) const {
    _StringBuffer result (256UL);
    for (unsigned long charIndex = 0UL; charIndex < 256UL; charIndex++) {
        if (charMap.list_data[charIndex] >= 0)
            result << char (charIndex);
    }
    return result;
}

//----------------------------------------------------------------------------------------------------------------------

void _Trie::SetAlphabet (const _String* alphabet, bool doClear)
{
    if (doClear) {
        Clear (true);
        charMap.Clear();
    }

    if (alphabet) {
        charMap.Populate (256,-1,0);
        unsigned long charCounter = 0;
        charMap.list_data[0] = 1; // always allow the '\0' character
        for (unsigned long charIndex = 0; charIndex < alphabet->length(); charIndex++) {
            charMap.list_data [(unsigned char)alphabet->char_at(charIndex)] = 1;
        }
        // now sort alphabetically
        for (unsigned long charIndex = 0; charIndex < 256; charIndex++) {
            if (charMap.list_data[charIndex] == 1)
                charMap.list_data[charIndex] = charCounter++;
        }
    } else {
        charMap.Populate (256,0,1);
    }
}

//----------------------------------------------------------------------------------------------------------------------

BaseRef _Trie::makeDynamic (void) const {
    _Trie *newTrie = new _Trie ();
    newTrie->Duplicate (this);
    return newTrie;
          
}

//----------------------------------------------------------------------------------------------------------------------

void _Trie::Duplicate (BaseRefConst t) {
    _Trie* source = (_Trie*)t;
    _String a (source->Alphabet());
    SetAlphabet(&a, true);
    _List::Duplicate ((_List*)source);
    charMap.Duplicate (&source->charMap);
    emptySlots.Duplicate (&source->emptySlots);
    payload.Duplicate(&source->payload);
    parents.Duplicate(&source->parents);
    inserted_values = source->inserted_values;
          
}        
    
//----------------------------------------------------------------------------------------------------------------------

long    _Trie::FindNextLetter (const char letter, const unsigned long current_index) const {
    long letterKey = charMap.list_data[(const unsigned char)letter];
    if (letterKey >= 0) {
        _SimpleList* thisList = ((_SimpleList**)list_data)[current_index];
        letterKey = thisList->FindStepping (letterKey, 2, 0);
        if (letterKey < 0)
            return kNotFound;
        return thisList->list_data[letterKey+1];
    }
    return kTrieInvalidLetter;
}

//----------------------------------------------------------------------------------------------------------------------
long  _Trie::FindNextUnusedIndex (bool alloc){
    
    if (emptySlots.lLength) {
        long newIndex = emptySlots.Pop();
        if (alloc)
            ((_SimpleList**)list_data)[newIndex] = new _SimpleList;
        return newIndex;
    }  
    payload << 0;
    parents << 0;
    if (alloc)
        AppendNewInstance(new _SimpleList);
    else
        *((_SimpleList*)this)<<0L;

    return lLength - 1;
}

//----------------------------------------------------------------------------------------------------------------------

long    _Trie::InsertNextLetter (const char letter, const unsigned long current_index) {
    long letter_key = charMap.list_data[(const unsigned char)letter];
    if (letter_key >= 0) {
        long next_index = FindNextUnusedIndex (letter != 0);
        _SimpleList * currentList = ((_SimpleList**)list_data)[current_index];
        (*currentList) << letter_key;
        (*currentList) << next_index;
        parents.list_data[next_index] = current_index;
        return next_index;
    }
    return kTrieInvalidLetter;
}


//----------------------------------------------------------------------------------------------------------------------
long     _Trie::FindKey (const _String& key, _SimpleList* path, bool prefixOK, unsigned long * start_index) const{
    long current_index = 0L,
         next_index    = 0L;
    for (unsigned long k = start_index ? *start_index : 0UL; k <= key.length() && current_index >= 0L; k++){
       next_index = FindNextLetter (key.char_at(k), current_index);
       if (path)
            (*path) << next_index;
       if (next_index < 0 && prefixOK) {
           next_index = FindNextLetter (0, current_index);
           current_index = next_index;
           if (start_index) {
             *start_index = k;
           }
           return current_index;
       }
       current_index = next_index;
    }
    if (start_index) {
      *start_index = key.length();
    }
    return current_index;
}


//----------------------------------------------------------------------------------------------------------------------
long     _Trie::FindKey (const char key, bool prefixOK) const {
    long current_index = 0,
    next_index    = FindNextLetter (key, current_index);
    if (next_index < 0 && prefixOK) {
        next_index = FindNextLetter (0, current_index);
    }
    current_index = next_index;
    return current_index;
}

//----------------------------------------------------------------------------------------------------------------------
long     _Trie::GetValueFromString (const _String& key){
    long keyIndex = FindKey (key);
    if (keyIndex != kNotFound) {
        return GetValue (keyIndex);
    }
    return kNotFound;
}


//----------------------------------------------------------------------------------------------------------------------
void     _Trie::UpdateValue(const long key, const long value) {
    if (key >= 0 && key < payload.lLength)
        payload.list_data[key] = value;
}

//----------------------------------------------------------------------------------------------------------------------
long     _Trie::GetValue(const long key) const{
    if (key >= 0 && key < payload.lLength)
        return payload.list_data[key];
    
    return 0L;
}

//----------------------------------------------------------------------------------------------------------------------

long    _Trie::Insert (const char* key, const long value, bool return_index) {
    _String key_string(key);
    long ret_value = Insert  (key_string, value);
    if (ret_value >= 0 && return_index == false) {
        return key_string.length();
    }
    return ret_value;
}

//----------------------------------------------------------------------------------------------------------------------

_Trie& _Trie::operator < (const char* key) {
  Insert (_String (key), 0L);
  return *this;
}

//----------------------------------------------------------------------------------------------------------------------

long    _Trie::Insert (const _String& key, const long value) {
    return InsertExtended (key, value);
}
    
//----------------------------------------------------------------------------------------------------------------------

long    _Trie::InsertExtended (const _String& key, const long value, bool update_value, bool * did_insert) {
    // the root is always at index 0
    long current_index  = 0L,
         current_char   = 0L,
         next_index     = FindNextLetter(key.char_at (current_char++), current_index);
    
    while (next_index >= 0L && current_char <= key.length()) {
        current_index = next_index;
        next_index     = FindNextLetter(key.char_at (current_char++), current_index);
    }
    
    if (next_index == kTrieInvalidLetter)
        return kTrieInvalidLetter;
    
    if (current_char > key.length() ) { // key already present
        if (update_value) {
            UpdateValue (next_index, value);
        }
        if (did_insert) *did_insert = false;
        return next_index;
    }
    
    current_char --;
    
    // validate the rest of the string
    
    for (long k = current_char; k <= key.length(); k++) {
        if (charMap[key.char_at(k)] < 0)
            return kTrieInvalidLetter;
    }
    
     // insert the rest of the string
    for (; current_char <= key.length(); current_char++) {
        //printf ("\nInserting %c\n", key.sData[current_char]);
        current_index = InsertNextLetter (key.char_at (current_char), current_index);
        //DumpRaw ();
    }
    inserted_values ++;
    UpdateValue (current_index, value);
    if (did_insert) *did_insert = true;

    return current_index;
}

//----------------------------------------------------------------------------------------------------------------------

unsigned long    _Trie::Insert (const _List& key, const _SimpleList* values) {
    unsigned long how_many = 0;
    for (long k = 0; k < key.lLength; k++) {
        _String serializedKey ((_String*)((BaseRef*)key.list_data)[k]->toStr());
        bool did_insert = false;
        
        InsertExtended (serializedKey, values?values->list_data[k]:0, false, &did_insert);
        if (did_insert) {
            how_many ++;
        }
    }
    return how_many;
}

//----------------------------------------------------------------------------------------------------------------------
bool    _Trie::Delete (const _String& key){
    _SimpleList history;
    long found_key = FindKey (key, &history);
    if (found_key >= 0) {
        // now traverse the history list backwards and delete all keys that have no children
        for (long k = history.lLength-1; k>=0; k--) {
            _SimpleList * current_list = ((_SimpleList**)list_data)[history.list_data[k]];
            if (current_list == nil || current_list->lLength <= 1){
                emptySlots << history.list_data[k];
                payload.list_data[history.list_data[k]] = 0L;
                parents.list_data[history.list_data[k]] = -1L;
                _SimpleList * parentList = ((_SimpleList**)list_data)[history.list_data[k-1]];
                unsigned long parentNode = parentList->FindStepping (history.list_data[k],2, 1) - 1;
                parentList->Delete (parentNode);
                parentList->Delete (parentNode);
                DeleteObject (current_list);
                ((_SimpleList**)list_data)[history.list_data[k]] = nil;
                inserted_values--;
            }
        }
        return true;
    }
    return false;
}

//----------------------------------------------------------------------------------------------------------------------
bool    _Trie::Delete (const char* key){
    _String string_key (key);
    return Delete (string_key);
}

//----------------------------------------------------------------------------------------------------------------------
unsigned long     _Trie::Delete (const _List& key){
     unsigned long how_many = 0;
    for (long k = 0; k < key.lLength; k++) {
        _String serializedKey ((_String*)((BaseRef*)key.list_data)[k]->toStr());
        if (Delete (serializedKey)) {
            how_many ++;
        }
    }
    return how_many;
}

 //----------------------------------------------------------------------------------------------------------------------
_String*         _Trie::RetrieveStringFromPath (const _SimpleList& path, _String const * alphabet) const {
  
    _StringBuffer * this_string = new _StringBuffer (128UL);
  
  
    _String     * local_alph      = nil;
    _String  const * read_from = alphabet;
  
    if (!read_from) {
      local_alph = new _String (Alphabet());
      read_from = local_alph;
    }
  
    
    for (long k = 0; k < path.lLength - 4; k+=2) {
         _SimpleList* current_list     = ((_SimpleList**)list_data)[path.list_data[k]];
         long         current_position = path.list_data[k+1];
         (*this_string) << read_from->char_at (current_list->list_data[current_position]);
    }
    
  
    DeleteObject(local_alph);
    return this_string;
}
  
//----------------------------------------------------------------------------------------------------------------------
void    _Trie::DumpRaw() {
    _String alph       = Alphabet(); 
    for (long k = 0; k < lLength; k++) {
        if (emptySlots.Find(k) < 0) {
            printf ("Position %ld:\n", k);
            _SimpleList * this_list = ((_SimpleList**)list_data)[k];
            for (long m = 0; m < this_list->lLength; m+=2) {
                printf ("'%c'(%ld) -> %ld\n", (char)this_list->list_data[m], this_list->list_data[m], this_list->list_data[m+1]);
            }
            
        } else {
            printf ("Position %ld is unused\n", k);
        }
    }
}
                                                            
//----------------------------------------------------------------------------------------------------------------------

BaseRef     _Trie::toStr(unsigned long) {

    _StringBuffer * serialized = new _StringBuffer (128UL);
    _String const   alph       = Alphabet();
  
    _SimpleList       traversal_history, 
                        // 2 indices per entry: node and current position (in multiples of 2)
                      *root_list = ((_SimpleList**)list_data)[0];
                      
    traversal_history << 0; traversal_history << 0;
    
    bool doComma = false;
    
    (*serialized) << '{';
    while (!(traversal_history.lLength == 2 && traversal_history.list_data[1] == root_list->lLength)) {
        _SimpleList* current_list = ((_SimpleList**)list_data)[traversal_history.list_data[traversal_history.lLength-2]];
        long current_position = traversal_history.list_data[traversal_history.lLength-1];
        // if current list is empty, then generate a string based on the path, and advance up the chain
        if (current_list && current_list->lLength) {
            if (current_position < current_list->lLength) {
                traversal_history << current_list->list_data[current_position+1];
                traversal_history << 0;
            } else {
                traversal_history.Pop();
                traversal_history.Pop();
                traversal_history.list_data[traversal_history.lLength-1] += 2; // advance the counter in the parent
           }
        } else {
            _String * this_string = RetrieveStringFromPath(traversal_history, &alph);
            (*serialized) << '"' << this_string << "\":" << _String (GetValue (traversal_history.list_data[traversal_history.lLength-2]));
            if (doComma) {
                (*serialized) << ',';
            } else {
                doComma = true;
            }
            (*serialized) << '\n';
            traversal_history.Pop();
            traversal_history.Pop();            
            traversal_history.list_data[traversal_history.lLength-1] += 2; // advance the counter in the parent
       }
    }
    
    (*serialized) << '}';
    return serialized;
}
 
//----------------------------------------------------------------------------------------------------------------------

_String  const _Trie::RetrieveKeyByPayload (const long key){
    long key_index = payload.Find (key);
    if (key_index >= 0) {
        _SimpleList parent_indices,
                    traversal_history;
        long keyer = key_index;
        do{
            parent_indices << keyer;
            keyer = parents.list_data[keyer];
        }
        while (keyer > 0);
        parent_indices << 0;
        parent_indices.Flip();
        
        
        for (long i = 0; i < parent_indices.lLength-1; i++) {
            traversal_history << parent_indices.list_data[i];
            traversal_history << (((_SimpleList**)list_data)[parent_indices.list_data[i]])->FindStepping (parent_indices.list_data[i+1],2,1)-1;
        }
        traversal_history << key_index;
        traversal_history << 0L;
        _String alph = Alphabet();
        return _String(RetrieveStringFromPath(traversal_history, &alph));
        
    }
    return hy_global::kEmptyString;
}


