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

//----------------------------------------------------------------------------------------------------------------------

_Trie::_Trie (const _String* alphabet) {
    SetAlphabet (alphabet, false);
    AppendNewInstance(new _SimpleList);
    payload << 0L;
    parents <<-1L;
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
    payload << 0L;
    parents <<-1L;
}

//----------------------------------------------------------------------------------------------------------------------

_String _Trie::Alphabet (void)
{
    _String result (256L, true);
    for (unsigned long charIndex = 0; charIndex < 256; charIndex++) {
        if (charMap.lData[charIndex] >= 0)
            result << char (charIndex);
    }
    result.Finalize();
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
        charMap.lData[0] = 1; // always allow the '\0' character
        for (unsigned long charIndex = 0; charIndex < alphabet->sLength; charIndex++) {
            charMap.lData [(unsigned char)alphabet->sData[charIndex]] = 1;
        }
        // now sort alphabetically
        for (unsigned long charIndex = 0; charIndex < 256; charIndex++) {
            if (charMap.lData[charIndex] == 1)
                charMap.lData[charIndex] = charCounter++;
        }
    } else {
        charMap.Populate (256,0,1);
    }
}

//----------------------------------------------------------------------------------------------------------------------

BaseRef _Trie::makeDynamic (void) {
    _Trie *newTrie = new _Trie ();
    newTrie->Duplicate (newTrie);
    return newTrie;
          
}

//----------------------------------------------------------------------------------------------------------------------

void _Trie::Duplicate (BaseRef storage) {
    _Trie* newTrie = (_Trie*)storage;
    _String myAlphabet = Alphabet();
    newTrie->SetAlphabet (&myAlphabet, true);
    newTrie->_List::Duplicate ((_List*)this);
    newTrie->charMap.Duplicate (&charMap);
    newTrie->emptySlots.Duplicate (&emptySlots);
    newTrie->payload.Duplicate(&payload);
    newTrie->parents.Duplicate(&parents);
          
}        
    
//----------------------------------------------------------------------------------------------------------------------

long    _Trie::FindNextLetter (const char letter, const unsigned long current_index) {
    long letterKey = charMap.lData[(const unsigned char)letter];
    if (letterKey >= 0) {
        _SimpleList* thisList = ((_SimpleList**)lData)[current_index];
        letterKey = thisList->FindStepping (letterKey, 2, 0);
        if (letterKey < 0)
            return HY_TRIE_NOTFOUND;
        return thisList->lData[letterKey+1];
    }
    return HY_TRIE_INVALID_LETTER;
}

//----------------------------------------------------------------------------------------------------------------------
long  _Trie::FindNextUnusedIndex (bool alloc){
    
    if (emptySlots.lLength) {
        long newIndex = emptySlots.Pop();
        if (alloc)
            ((_SimpleList**)lData)[newIndex] = new _SimpleList;
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
    long letter_key = charMap.lData[(const unsigned char)letter];
    if (letter_key >= 0) {
        long next_index = FindNextUnusedIndex (letter != 0);
        _SimpleList * currentList = ((_SimpleList**)lData)[current_index];
        (*currentList) << letter_key;
        (*currentList) << next_index;
        parents.lData[next_index] = current_index;
        return next_index;
    }
    return HY_TRIE_INVALID_LETTER;
}


//----------------------------------------------------------------------------------------------------------------------
long     _Trie::Find (const _String& key, _SimpleList* path, bool prefixOK){
    long current_index = 0,
         next_index    = 0;
    for (long k = 0; k <= key.sLength && current_index >= 0; k++){
       next_index = FindNextLetter (key.sData[k], current_index);
       if (path)
            (*path) << next_index;
       if (next_index < 0 && prefixOK) {
           next_index = FindNextLetter (0, current_index);
           current_index = next_index;
           break;
       }
       current_index = next_index;
    }
    return current_index;
}

//----------------------------------------------------------------------------------------------------------------------
long     _Trie::Find (const char key, bool prefixOK){
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
    long keyIndex = Find(key);
    if (keyIndex != HY_TRIE_NOTFOUND) {
        return GetValue (keyIndex);
    }
    return HY_TRIE_NOTFOUND;
}


//----------------------------------------------------------------------------------------------------------------------
void     _Trie::UpdateValue(const long key, const long value) {
    if (key >= 0 && key < payload.lLength)
        payload.lData[key] = value;
}

//----------------------------------------------------------------------------------------------------------------------
long     _Trie::GetValue(const long key) {
    if (key >= 0 && key < payload.lLength)
        return payload.lData[key];
    
    return 0L;
}

//----------------------------------------------------------------------------------------------------------------------

long    _Trie::Insert (const char* key, const long value, bool return_index) {
    _String key_string(key);
    long ret_value = Insert  (key_string, value);
    if (ret_value >= 0 && return_index == false) {
        return key_string.sLength;
    }
    return ret_value;
}


//----------------------------------------------------------------------------------------------------------------------

long    _Trie::Insert (const _String& key, const long value) {
    // the root is always at index 0
    long current_index  = 0, 
         current_char   = 0,
         next_index     = FindNextLetter(key.sData[current_char++], current_index);
    
    while (next_index >= 0 && current_char <= key.sLength) {
        current_index = next_index;
        next_index     = FindNextLetter(key.sData[current_char++], current_index);
    }
    
    if (next_index == HY_TRIE_INVALID_LETTER)
        return HY_TRIE_INVALID_LETTER;
    
    if (current_char == key.sLength && next_index >= 0)
        return next_index;
    
    current_char --;
    
    // validate the rest of the string
    
    for (long k = current_char; k <= key.sLength; k++) {
        if (charMap[key.sData[k]] < 0)
            return HY_TRIE_INVALID_LETTER;
    }
    
     // insert the rest of the string
    for (; current_char <= key.sLength; current_char++) {
        //printf ("\nInserting %c\n", key.sData[current_char]);
        current_index = InsertNextLetter (key.sData[current_char], current_index);      
        //DumpRaw ();
   }
    
    UpdateValue (current_index, value);

    return current_index;
}

//----------------------------------------------------------------------------------------------------------------------

unsigned long    _Trie::Insert (const _List& key, const _SimpleList* values) {
    unsigned long how_many = 0;
    for (long k = 0; k < key.lLength; k++) {
        _String serializedKey ((_String*)((BaseRef*)key.lData)[k]->toStr());
        
        long this_index = Insert (serializedKey, values?values->lData[k]:0);
        if (this_index >= 0) {
            how_many ++;
        }
    }
    return how_many;
}

//----------------------------------------------------------------------------------------------------------------------
bool    _Trie::Delete (const _String& key){
    _SimpleList history;
    long found_key = Find (key, &history);
    if (found_key >= 0) {
        // now traverse the history list backwards and delete all keys that have no children
        for (long k = history.lLength-1; k>=0; k--) {
            _SimpleList * current_list = ((_SimpleList**)lData)[history.lData[k]];
            if (current_list == nil || current_list->lLength <= 1){
                emptySlots << history.lData[k];
                payload.lData[history.lData[k]] = 0L;
                parents.lData[history.lData[k]] = -1L;
                _SimpleList * parentList = ((_SimpleList**)lData)[history.lData[k-1]];
                unsigned long parentNode = parentList->FindStepping (history.lData[k],2, 1) - 1;
                parentList->Delete (parentNode);
                parentList->Delete (parentNode);
                DeleteObject (current_list);
                ((_SimpleList**)lData)[history.lData[k]] = nil;
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
        _String serializedKey ((_String*)((BaseRef*)key.lData)[k]->toStr());
        
        long this_index = Delete (serializedKey);
        
        if (this_index >= 0) {
            how_many ++;
        }
    }
    return how_many;
}

 //----------------------------------------------------------------------------------------------------------------------
_String*         _Trie::RetrieveStringFromPath (const _SimpleList& path, _String* alphabet) {
    _String* this_string = new _String (128L,true),
           * my_alph      = alphabet? alphabet : new _String (Alphabet());
    
    
    for (long k = 0; k < path.lLength - 4; k+=2) {
         _SimpleList* current_list     = ((_SimpleList**)lData)[path.lData[k]];
         long         current_position = path.lData[k+1];
         (*this_string) << my_alph->sData[current_list->lData[current_position]];
    }
    
    this_string->Finalize();
    
    if (!alphabet)
        DeleteObject(my_alph);
    
    return this_string;
}
  
//----------------------------------------------------------------------------------------------------------------------
void    _Trie::DumpRaw() {
    _String alph       = Alphabet(); 
    for (long k = 0; k < lLength; k++) {
        if (emptySlots.Find(k) < 0) {
            printf ("Position %ld:\n", k);
            _SimpleList * this_list = ((_SimpleList**)lData)[k];
            for (long m = 0; m < this_list->lLength; m+=2) {
                printf ("'%c'(%ld) -> %ld\n", (char)this_list->lData[m], this_list->lData[m], this_list->lData[m+1]);
            }
            
        } else {
            printf ("Position %ld is unused\n", k);
        }
    }
}
                                                            
//----------------------------------------------------------------------------------------------------------------------

BaseRef     _Trie::toStr() {
    _String         * serialized = new _String (128L, true),
                      alph       = Alphabet();
                      
    _SimpleList       traversal_history, 
                        // 2 indices per entry: node and current position (in multiples of 2)
                      *root_list = ((_SimpleList**)lData)[0];
                      
    traversal_history << 0; traversal_history << 0;
    
    bool doComma = false;
    
    (*serialized) << '{';
    while (!(traversal_history.lLength == 2 && traversal_history.lData[1] == root_list->lLength)) {
        _SimpleList* current_list = ((_SimpleList**)lData)[traversal_history.lData[traversal_history.lLength-2]];
        long current_position = traversal_history.lData[traversal_history.lLength-1];
        // if current list is empty, then generate a string based on the path, and advance up the chain
        if (current_list && current_list->lLength) {
            if (current_position < current_list->lLength) {
                traversal_history << current_list->lData[current_position+1];
                traversal_history << 0;
            } else {
                traversal_history.Pop();
                traversal_history.Pop();
                traversal_history.lData[traversal_history.lLength-1] += 2; // advance the counter in the parent
           }
        } else {
            _String * this_string = RetrieveStringFromPath(traversal_history, &alph);
            (*serialized) << '"';
            (*serialized) << this_string;
            (*serialized) << "\":";
            (*serialized) << _String (GetValue (traversal_history.lData[traversal_history.lLength-2]));
            if (doComma) {
                (*serialized) << ',';
            } else {
                doComma = true;
            }
            (*serialized) << '\n';
            traversal_history.Pop();
            traversal_history.Pop();            
            traversal_history.lData[traversal_history.lLength-1] += 2; // advance the counter in the parent
       }
    }
    
    (*serialized) << '}';
    serialized->Finalize();
    return serialized;
}
 
//----------------------------------------------------------------------------------------------------------------------

_String  _Trie::RetrieveKeyByPayload (const long key){
    long key_index = payload.Find (key);
    if (key_index >= 0) {
        _SimpleList parent_indices,
                    traversal_history;
        long keyer = key_index;
        do{
            parent_indices << keyer;
            keyer = parents.lData[keyer];
        }
        while (keyer > 0);
        parent_indices << 0;
        parent_indices.Flip();
        
        
        for (long i = 0; i < parent_indices.lLength-1; i++) {
            traversal_history << parent_indices.lData[i];
            traversal_history << (((_SimpleList**)lData)[parent_indices.lData[i]])->FindStepping (parent_indices.lData[i+1],2,1)-1;
        }
        traversal_history << key_index;
        traversal_history << 0L;
        _String alph = Alphabet();
        return _String(RetrieveStringFromPath(traversal_history, &alph));
        
    }
    return empty;
}


