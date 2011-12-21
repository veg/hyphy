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
}
       
//----------------------------------------------------------------------------------------------------------------------

_Trie::~_Trie (void){
    
}
       
//----------------------------------------------------------------------------------------------------------------------

void _Trie::Clear (bool all){
    _List::Clear (all);
    payload.Clear(all);
    emptySlots.Clear(all);
}

//----------------------------------------------------------------------------------------------------------------------

_String _Trie::Alphabet (void)
{
    _String result (256L, true);
    for (unsigned long charIndex = 1; charIndex < 255; charIndex++) {
        if (charMap.lData[charIndex] >= 0)
            result << char (charIndex+1);
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
        for (unsigned long charIndex = 0; charIndex < alphabet->sLength; charIndex++) {
            charMap.lData [(unsigned char)alphabet->sData[charIndex]] = 1;
        }
        // now sort alphabetically
        for (unsigned long charIndex = 0; charIndex < 256; charIndex++) {
            if (charMap.lData[charIndex] == 1)
                charMap.lData[charIndex] = charCounter++;
        }
    } else {
        charMap.Populate (-1,256,1);
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
 
                        virtual BaseRef toStr(void);
        /**
         * Return a string representation of this object
         * @return A _String reference to the list of strings (one per line) currently stored in the trie 
         */

       
        virtual void    Duplicate(BaseRef storage);
        /**
         * Perform a deep copy of this object into storage (an allocated empty trie)
         * @param storage -- the _Trie object to copy this one into
         * @return Nothing. 
         */
       

        long     Find (const _String& key);
        /**
         * Determine if 'key' is in the trie
         * @param  key -- the string to search for
         * @return the index of the key in 'nodes' if found, -1 otherwise  
         */
        
        long    Insert (const _String& key, const long value);
        /**
         * Insert the key into the trie
         * @param key -- the string to insert
         * @param value -- the value to associate with the key
         * @return non-negative index if the insert was successful (also returned if key is already in this trie), otherwise -1 
         */
    
        void     UpdateValue (const long key, const long value);
        /**
         * Update the value associated with the key _index_
         * @param  key -- the index of the key (returned by Find for example); if key < 0 or key >= nodes.lLength, nothing is done
         * @param  value -- new value to associate with the key
         * @return None
         */
        
        unsigned long    Insert (const _List& key, const _SimpleList* values = nil);
        /**
         * Insert all keys in the list into the trie
         * @param key -- the list of strings (non string objects will be cast to strings) to insert
         * @param values -- the list of values to associate with the keys (the index of the key in the list by default)
         * @return the number of elements successfully inserted (including those already present)
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


