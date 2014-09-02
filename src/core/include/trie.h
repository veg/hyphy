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

#include "hy_strings.h"
#include "hy_list.h"
#include "hy_list_numeric.h"
#include "hy_list_reference.h"

#define HY_TRIE_NOTFOUND HY_NOT_FOUND
#define HY_TRIE_INVALID_LETTER -2L

/*
  This is a simple class for representing prefix tries with integer values
  attached to each string key.

   The index of the parent nodes (in 'nodes') for each node   
   base class will store the linear representation of this trie 

   A linear representation of this trie:
     each node is a _hyListNumeric<long> that contains N pairs of entries (character
     index, integer index of the child node in 'nodes') for each extension of
     the prefix encoded by the node
   
   For example, if the alphabet is "ABC" and if nodes[1] has children nodes[5]
   for "B" and nodes[7] for "A", then the _hyListNumeric<long> for nodes[1] will be 1
   (index of B),5 ('B' child index),0 (index of A),7 ('A' child index).

*/

class _Trie {

private:

protected:

  // data members
  _hyListNumeric<char> char_map,          /** 
                                  char_map[c] maps a valid character to the internal
                                  index (0..validChars.sLength)
                                  invalid characters are mapped to -1
                                  e.g. if the alphabet is "CGTA", then char_map ['A']  = 3,
                                  and char_map['z'] = -1
                                  */
                payload,          /** the values associated with each key in 'nodes' */
                parents;



public:
  _hyListReference< _hyListNumeric<long> > linear_list;

  /**
   * Construct an empty trie over a given alphabet
   * @param alphabet -- a string listing all valid characters (e.g. "ACGT"). By
   * default (or if an empty string is passed), all ASCII characters are allowed
   * @return Nothing.
   */
   _Trie(_String* alphabet = nil);
   _Trie(_Trie&);

  /**
   * Return a string representation of this object
   * @return An dictionary representation of the trie, i.e. {"key1":"value1",
   * "key2":"value2", ...} pairs
   */
  virtual _StringBuffer toStr(void);

  /**
   * Clear this trie; everything except the alphabet will be deleted
   * @param all -- whether or not to clear out the lists completely
   * @return Nothing
   */

  virtual void clear(bool all = TRUE);

  /**
   * The destructor
   * @return Nothing.
   */

  virtual ~_Trie(void);
   
  /**
   * Counts the number of strings in this Trie
   * @return the number of strings.
   */
  virtual unsigned long countitems (void) const;

  /**
   * Counts the number of strings in this Trie
   * @return the number of strings.
   */
  virtual unsigned long Length     (void) const;
  
  

  /**
   * Determine if 'key' is in the trie
   * @param  key      -- the string to search for
   * @param  path     -- store the indices for the trie traversal history (if
   * supplied)
   * @param  prefixOK -- returns a match if a prefix of 'key' in the trie
   * @return the index of the key in 'nodes' if found,
   * HY_TRIE_NOTFOUND/HY_TRIE_INVALID_LETTER otherwise
   */
  long Find(const _String key, _hyListNumeric<long> *path = nil, bool prefixOK = false) const;

  /**
   * Determine if 'key' is in the trie
   * @param  key      -- the character to search for
   * @param  prefixOK -- returns a match if a prefix of 'key' in the trie
   * @return the index of the key in 'nodes' if found,
   * HY_TRIE_NOTFOUND/HY_TRIE_INVALID_LETTER otherwise
   */
  long Find(const char key, bool prefixOK = false) const;

  /**
   * A convenience function which calls Find and then GetValue if the key is
   * found
   * @param  key      -- the string to search for
   * @return the value associated with the key if found, HY_TRIE_NOTFOUND
   * otherwise
   */
  long GetValueFromString(const _String &key);

  /**
   * Insert the key into the trie
   * @param key -- the string to insert
   * @param value -- the value to associate with the key
   * @return non-negative index if the insert was successful (also returned if
   * key is already in this trie), otherwise
   * HY_TRIE_NOTFOUND/HY_TRIE_INVALID_LETTER
   */
  long Insert(const _String &key, const long value);

  /**
   * Insert the key into the trie
   * @param key -- the string to insert
   * @param value -- the value to associate with the key
   * @param return_index - whether or not to return the index of the string in
   * the trie (if true) or the length of the key (if false)
   * @return non-negative index if the insert was successful (also returned if
   * key is already in this trie), otherwise
   * HY_TRIE_NOTFOUND/HY_TRIE_INVALID_LETTER; if return_index == false, return
   * strlen (key) if insert was successful
   */
  long Insert(const char *key, const long value, bool return_index = true);

  /**
   * Update the value associated with the key _index_
   * @param  key -- the index of the key (returned by Find for example); if key
   * < 0 or key >= nodes.lLength, nothing is done
   * @param  value -- new value to associate with the key
   * @return None
   */
  void UpdateValue(const long key, const long value);

  /**
     * Retrieve the value associated with the key _index_
     * @param  key -- the index of the key (returned by Find for
     * example); if key < 0 or key >= nodes.lLength, nothing is done
     * @return the value associated with the key
     */
  long GetValue(const long key) const;

  /**
   * Insert all keys in the list into the trie
   * @param key -- the list of strings (non string objects will be cast to
   * strings) to insert
   * @param values -- the list of values to associate with the keys (the index
   * of the key in the list by default)
   * @return the number of elements successfully inserted (including those
   * already present)
   */
  unsigned long Insert(const _hyListReference<_String> key, const _hyListNumeric<long> values = nil);

  /**
   * Delete the key from the trie
   * @param key -- the string to delete
   * @return True if the delete was successful (also returned if key is not in
   * this trie), otherwise False
   */
  bool Delete(_String key);

  /**
   * Delete the key from the trie
   * @param key -- the string to delete
   * @return True if the delete was successful (also returned if key is not in
   * this trie), otherwise False
   */
  bool Delete(const char *key);

  /**
   * Delete all keys in the list from the trie
   * @param key -- the list of strings (non string objects will be cast to
   * strings) to delete
   * @return the number of elements successfully deleted (including those not
   * present)
   */
  unsigned long Delete(const _hyListReference<_String> &key);

  /**
  * Given a traversal path of the trie (and an optional cached alphabet),
  * retrive the _String object spelling the path
  * @param path -- the traversal path (pairs of node index, character index)
  * @param
  * @return the string spelling the path
  */
  _StringBuffer RetrieveStringFromPath(_hyListNumeric<long> path, _hyList<char> alphabet);

  /**
   * Return the valid alphabet for this Trie
   * @return The string containing all the letters allowed for strings in this
   * trie. The ordering of the letters is ASCII-alphabetical.
   */
  _hyList<char> alphabet(void) const;

  /**
   * Return the string spelling the pay to the 'key'
   * @param key -- the key for which we fish to retrieve the path
   * @return The string spelling the path from the root to the node tagged with
   * value 'key'. Empty string is returned if 'key' is not in this trie
   */
  _String RetrieveKeyByPayload(const long key);

private:

  void setAlphabet(const _String*, bool);

  /**
   * Given a current position in the trie (current_index), try to walk down the
   * next character
   * @param  letter -- the next letter
   * @param  current_index -- where in the trie are we currently located
   * @return A non-negative index (next position) in the trie; HY_TRIE_NOTFOUND/
   * if the letter were valid but no extension could be found, and
   * HY_TRIE_INVALID_LETTER if the letter were invalid
   */
  long FindNextLetter(const char letter, const unsigned long currentIndex) const;

  /**
   * Given a current position in the trie (current_index), insert the character
   * (this assumes that the character is NOT present)
   * @param  letter -- the next letter
   * @param  current_index -- where in the trie are we currently located
   * @return A non-negative index (next position) in the trie; HY_TRIE_NOTFOUND/
   * if the letter were valid but no extension could be found, and
   * HY_TRIE_INVALID_LETTER if the letter were invalid
   */
  long InsertNextLetter(const char letter, const unsigned long currentIndex);

  /**
      Find the next index to store something to: place an empty _hyListNumeric<long>
      there
      This will either go to the end the list or to the last freed block (in
      empty_slots)
      @param alloc -- allocate a new _hyListNumeric<long> storage object
      @return the index of the new empty _hyListNumeric<long> in (this) List
    */
  long FindNextUnusedIndex(bool alloc = TRUE);

  void DumpRaw(void);

};

#endif
