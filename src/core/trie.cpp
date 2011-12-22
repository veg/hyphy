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
        charMap.Populate (0,256,-1);
        unsigned long charCounter = 0;
        charMap.lData[0] == 1; // always allow the '\0' character
        for (unsigned long charIndex = 0; charIndex < alphabet->sLength; charIndex++) {
            charMap.lData [(unsigned char)alphabet->sData[charIndex]] = 1;
        }
        // now sort alphabetically
        for (unsigned long charIndex = 0; charIndex < 256; charIndex++) {
            if (charMap.lData[charIndex] == 1)
                charMap.lData[charIndex] = charCounter++;
        }
    } else {
        charMap.Populate (0,256,1);
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
          
}        
    
//----------------------------------------------------------------------------------------------------------------------

long    _Trie::FindNextLetter (const char letter, const unsigned long current_index) {
    long letterKey = charMap.lData[(const unsigned char)letter];
    if (letterKey >= 0) {
        _SimpleList* thisList = ((_SimpleList**)lData)[current_index];
        letterKey = thisList->FindStepping (letterKey, 0, 2);
        if (letterKey < 0)
            return HY_TRIE_NOTFOUND;
        return letterKey;
    }
    return HY_TRIE_INVALID_LETTER;
}

//----------------------------------------------------------------------------------------------------------------------
long  _Trie::FindNextUnusedIndex (void){
    
    if (emptySlots.lLength) {
        long newIndex = emptySlots.Pop();
        ((_SimpleList**)lData)[newIndex] = new _SimpleList;
        return newIndex;
    }  
    payload << 0;
    AppendNewInstance(new _SimpleList);
    return lLength - 1;
}

//----------------------------------------------------------------------------------------------------------------------

long    _Trie::InsertNextLetter (const char letter, const unsigned long current_index) {
    long letter_key = charMap.lData[(const unsigned char)letter];
    if (letter_key >= 0) {
        long next_index = FindNextUnusedIndex ();
        _SimpleList * currentList = ((_SimpleList**)lData)[current_index];
        (*currentList) << letter_key;
        (*currentList) << next_index;
        return next_index;
    }
    return HY_TRIE_INVALID_LETTER;
}


//----------------------------------------------------------------------------------------------------------------------
long     _Trie::Find (const _String& key, _SimpleList* path){
    long current_index = 0;
    for (long k = 0; k <= key.sLength; k++){
       current_index = FindNextLetter (key.sData[k], current_index);
       if (path)
            (*path) << current_index;
       if (current_index < 0) 
          break;
    }
    return current_index;
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

long    _Trie::Insert (const _String& key, const long value) {
    // the root is always at index 0
    long current_index = 0, 
         current_char = 0,
         next_index     = FindNextLetter(key.sData[current_char++], current_index);
    
    while (next_index >= 0 && current_char <= key.sLength) {
        current_index = next_index;
        next_index     = FindNextLetter(key.sData[current_char++], current_index);
    }
    
    if (next_index == HY_TRIE_INVALID_LETTER)
        return HY_TRIE_INVALID_LETTER;
    
    if (current_char == key.sLength)
        return next_index;
    
    current_char --;
    
    // validate the rest of the string
    
    for (long k = current_char; k <= key.sLength; k++) {
        if (charMap[key.sData[k]] < 0)
            return HY_TRIE_INVALID_LETTER;
    }
    
    next_index = current_index;
    // insert the rest of the string
    for (; current_char < key.sLength; current_char++) {
        current_index = next_index;
        next_index = InsertNextLetter (key.sData[current_char++], current_index);      
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
        // now traverse the history list backwards and delete 
    }
    return false;
}
                      
//----------------------------------------------------------------------------------------------------------------------
 
        virtual BaseRef toStr(void);
        /**
         * Return a string representation of this object
         * @return A _String reference to the list of strings (one per line) currently stored in the trie 
         */

       
        

        bool    Delete (const _String& key);
        /**
         * Delete the key from the trie
         * @param key -- the string to delete
         * @return True if the delete was successful (also returned if key is not in this trie), otherwise False 
         */
    
        unsigned long    Delete (const _List& key);
        /**
         * Delete all keys in the list from the trie
         * @param key -- the list of strings (non string objects will be cast to strings) to delete
         * @return the number of elements successfully deleted (including those not present)
         */


