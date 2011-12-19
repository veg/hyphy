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

#ifndef _HY_TRIE_
#define _HY_TRIE_

#include "simplelist.h"
#include "list.h"

/*_____________________________________________________________________________
    This is a simple class for representing prefix tries with integer values 
    attached to each string key.
*/

//_____________________________________________________________________________
class _Trie: public BaseObj
{
    protected:
        // data members
        _SimpleList charMap,
            /** charMap[c] maps a valid character to the internal index (0..validChars.sLength)
             invalid characters are mapped to -1
             e.g. if the alphabet is "CGTA", then charMap ['A']  = 3, and charMap['z'] = -1 
             */
            emptySlots,
            /** allocated entries in the 'nodes' list that can be reused (e.g. those created by delete operations)
             */
            payload;
            /** the values associated with each key in 'nodes' 
             */
            
    
        _List       nodes;
            /** a linear representation of this trie    
             each node is a _SimpleList that contains pairs of entries 
                (character index, integer index of the child node in 'nodes')
                for each extension of the prefix encoded by the node
             
             for example, if the alphabet is "ABC" and if nodes[1] has children nodes[5] for "B" and nodes[7] for "A", then the _SimpleList for nodes[1] will be 1 (index of B),5 ('B' child index),0 (index of A),7 ('A' child index).
             
             */
        

    public:
        _Trie (const _String* alphabet = nil);
        /**
         * Construct an empty trie over a given alphabet
         * @param alphabet -- a string listing all valid characters (e.g. "ACGT"). By default (or if an empty string is passed), all ASCII characters are allowed
         * @return Nothing. 
         */
       
        virtual BaseRef toStr(void);
        /**
         * Return a string representation of this object
         * @return A _String reference to the list of strings (one per line) currently stored in the trie 
         */
        
        virtual BaseRef makeDynamic(void);
        /**
         * Return a dynamic representation of this object
         * @return A _Trie reference created on the heap which is 'deep-copied' (i.e. all dynamic objects have a reference count of 1)
         */
       
        virtual void    Duplicate(BaseRef storage);
        /**
         * Perform a deep copy of this object into storage (an allocated empty trie)
         * @param storage -- the _Trie object to copy this one into
         * @return Nothing. 
         */
       
    
        virtual void    Clear(void);
        /**
         * Clear this trie; everything except the alphabet will be deleted
         * @return Nothing
         */
   
        virtual ~_Trie ();
        /**
         * The destructor 
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

        
        _String  Alphabet (void);
        /**
         * Return the valid alphabet for this Trie
         * @return The string containing all the letters allowed for strings in this trie. The ordering of the letters is ASCII-alphabetical. 
         */
    
    
           

     
};

#endif
