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

#ifndef _HY_TRIE_ITERATOR_
#define _HY_TRIE_ITERATOR_

#include "trie.h"

/*_____________________________________________________________________________
  This is a simple class for representing prefix tries with integer values
  attached to each string key.
*/

//_____________________________________________________________________________
class _TrieIterator {
  
public:
  _TrieIterator (const _Trie&);
  ~_TrieIterator (void) {}

  _String* First(void);
  _String* Last(void);
  _String* Next(void);
  _String* Previous(void);
  
  //long CurrentIndex(void) const { return last_retrieved_index; }
  
private:
  
  void Init(bool = false);
  long MoveDown(bool move_left = true);

  /*
   * Starting at a current (incomplete) traversal position, move down to the nearest tip,
   * taking either left-most or right-most branches every time
   * @param move_left -- take left-most (false) or right-most (true) branches
   * @return the index of the terminal leaf found or HY_NOT_FOUND if already at a tip,
   or if the initial traversal position is invalid
  */
  long RollUp(bool move_left = false);

   /*
   * Starting at a current traversal position (tip), move up until there is,
   * taking either left-most or right-most branches every time
   * @param move_left -- move left whenever possible, otherwise more right
   * @return the index of the parent where the move was made,
   or HY_NOT_FOUND if the move could not be made
   */
  _StringBuffer* BuildPath(void);
  
  long last_retrieved_index;
  _hyListNumeric<long> traversal_history;
  _hyList<char> alphabet;
  const _Trie* data_source;

};

#endif
