/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-now
Core Developers:
  Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
  Art FY Poon    (apoon42@uwo.ca)
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

#ifndef _AVLLISTX_
#define _AVLLISTX_
// #pragma once
#include "avllist.h"
#include "avllistx_iterator.h"
#include "list.h"

#define MEMORYSTEP 8

//_____________________________________________________________________________

/**
 * @brief An AVL list with an extra data member for each node.
 */
class _AVLListX : public _AVLList {

public:
  /* SLKP: 20090817
      add key: index values from the list of strings
   */

  // Data Members
  _SimpleList xtraD;

  // Methods
  /**
   * @brief Construct a new _AVLListX object.
   *
   * @param sl A _SimpleList to construct from. The new list will contain
   * pointers to the data in the _SimpleList, but will not own the data. The
   * _SimpleList will be cleared after the new _AVLListX is created.
   */
  _AVLListX(_SimpleList *);

  /**
   * @brief Destroy the _AVLListX object.
   */
  virtual ~_AVLListX(void) {}
  /**
   * @brief Convert the list to a string representation.
   *
   * @param padding The number of spaces to use for padding.
   * @return A _String object representing the list.
   */
  virtual BaseRef toStr(unsigned long = 0UL);

  /**
   * @brief Clear the list.
   * If shallow is true, the data in the list will not be deleted.
   *
   * @param shallow Whether to perform a shallow clear.
   */
  virtual void Clear(bool = false);
  /**
   * @brief This is a virtual function that does nothing in this class.
   * It is meant to be overridden in derived classes to delete extra data
   * associated with a node. In this class, the extra data is a long integer
   * stored in the `xtraD` list.
   *
   * @param l The index of the node.
   */
  virtual void DeleteXtra(long);
  /**
   * @brief Populate the list from a _List.
   * The new list will contain pointers to the data in the _List, but will not
   * own the data. The _List will be cleared after the new _AVLListX is created.
   * The extra data will be the index of the item in the original list.
   *
   * @param l The _List to populate from.
   *
   * @example
   * \code{.cpp}
   * _List list;
   * list << (BaseRef)new _String("A") << (BaseRef)new _String("B");
   * _AVLListX avl(nil);
   * avl.PopulateFromList(list);
   * // avl will contain "A" and "B", with extra data 0 and 1 respectively.
   * \endcode
   */
  virtual void PopulateFromList(_List &);
  /**
   * @brief Get the extra data associated with a key.
   *
   * @param brc The key to search for.
   * @return The extra data associated with the key, or kNotFound if the key is
   * not found.
   */
  long GetDataByKey(BaseRefConst) const;
  /**
   * @brief Get the extra data associated with a key.
   *
   * @param l The key to search for.
   * @return The extra data associated with the key, or kNotFound if the key is
   * not found.
   */
  long GetDataByKey(long) const;
  /**
   * @brief Find a key and get the extra data associated with it.
   *
   * @param brc The key to search for.
   * @param not_found_value The value to return if the key is not found.
   * @return The extra data associated with the key, or not_found_value if the
   * key is not found.
   */
  long FindAndGetXtra(BaseRefConst, long not_found_value = kNotFound) const;

  /**
   * @brief Insert data into the list.
   *
   * @param br The data to insert.
   * @param l The extra data to associate with the key.
   * @param b If true, the data will not be copied.
   * @return The index of the inserted item.
   */
  virtual long InsertData(BaseRef, long, bool);
  /**
   * @brief Update the value of an item in the list.
   *
   * @param br The new value.
   * @param l1 The index of the item to update.
   * @param l2 The new extra data.
   * @return The index of the updated item.
   */
  virtual long UpdateValue(BaseRef, long, long);
  virtual long SetOrIncrement(BaseRef, long, long = 1);

  /**
   * @brief Set the extra data of an item in the list.
   *
   * @param l1 The index of the item.
   * @param l2 The new extra data.
   */
  void SetXtra(long, long);
  /**
   * @brief Get the extra data of an item in the list.
   *
   * @param l The index of the item.
   * @return The extra data.
   */
  long GetXtra(long) const;

  /**
   * @brief Clear the formulas in the list.
   * This is a hack to automate clearing lists that have pointers to formulas in
   * them.
   */
  void ClearFormulasInList(void) { xtraD.ClearFormulasInList(); }
};

#endif
