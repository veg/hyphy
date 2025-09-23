/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-now
Core Developers:
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
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

#ifndef _AVLLIST_
#define _AVLLIST_
// #pragma once
#include "list.h"

#define MEMORYSTEP 8

//_____________________________________________________________________________
/**
 * @brief An AVL tree implementation that stores pointers to BaseObj objects.
 * This class is used to store sorted lists of objects, and provides fast insertion, deletion, and searching.
 */
class _AVLList : public BaseObj {

public:
  // Data Members
  _SimpleList *dataList, leftChild, rightChild, balanceFactor, emptySlots;

  long root;

  // Methods
  /**
   * @brief Construct a new _AVLList object from a _SimpleList.
   * The new list will contain pointers to the data in the _SimpleList, but will not own the data.
   * The _SimpleList will be cleared after the new _AVLList is created.
   *
   * @param sl The _SimpleList to construct from.
   * @return A new _AVLList object.
   *
   * @example _SimpleList sl; sl << 1 << 2 << 3; _AVLList avl(&sl); // avl now contains 1, 2, 3
   */
  /**
   * @brief Construct a new _AVLList object
   *
   */
  _AVLList(void);

  _AVLList(_SimpleList *);
  /**
   * @brief Construct a new _AVLList object from another _AVLList.
   * This is a copy constructor.
   *
   * @param src The _AVLList to copy from.
   * @return A new _AVLList object.
   *
   * @example
   * \code{.cpp}
   * _AVLList avl1;
   * avl1.Insert((BaseRef)new _String("A"));
   * _AVLList avl2(avl1); // avl2 is a copy of avl1
   * \endcode
   */
  _AVLList(_AVLList const &src);

  /**
   * @brief Destroy the _AVLList object
   *
   */
  virtual ~_AVLList(void) {}
  /**
   * @brief Clear the list.
   * If shallow is true, the data in the list will not be deleted.
   *
   * @param shallow Whether to perform a shallow clear.
   */
  virtual void Clear(bool = false);
  /**
   * @brief Check if the list has data at a given index.
   *
   * @param index The index to check.
   * @return true if the list has data at the given index, false otherwise.
   */
  virtual bool HasData(long);
  /**
   * @brief Make a dynamic copy of the list.
   *
   * @return A new _AVLList object that is a copy of this one.
   */
  virtual BaseRef makeDynamic(void) const;
  /**
   * @brief Duplicate the list from a reference.
   * This is a deep copy.
   *
   * @param ref The reference to duplicate from.
   */
  virtual void Duplicate(BaseRefConst);
  /**
   * @brief Assignment operator.
   * This is a deep copy.
   *
   * @param rhs The _AVLList to assign from.
   */
  void operator=(_AVLList const &rhs);

  /**
   * @brief Reorder the list.
   * This function is used to re-balance the AVL tree after a series of insertions or deletions.
   *
   * @param sl A _SimpleList to store the reordered indices. If nil, a new list will be created.
   */
  virtual void ReorderList(_SimpleList * = nil);
  /**
   * @brief Insert data into the list.
   *
   * @param data The data to insert.
   * @param index The index to insert at. If -1, the data will be inserted at the end.
   * @param shallow If true, the data will not be copied.
   * @return The index where the data was inserted.
   */
  virtual long InsertData(BaseRef, long, bool);
  /**
   * @brief Convert the list to a string representation.
   *
   * @param padding The number of spaces to use for padding.
   * @return A _String object representing the list.
   */
  virtual BaseRef toStr(unsigned long = 0UL);
  /**
   * @brief Traverse the list in order and return the indices of the nodes.
   *
   * @param sl The _SimpleList to store the indices in.
   * @param index The starting index for the traversal.
   * @param node The node to start the traversal from. If -1, the traversal starts from the root.
   * @return The number of nodes traversed.
   */
  virtual long Traverser(_SimpleList &, long &, long = -1) const;
  /**
   * @brief Get the root of the AVL tree.
   *
   * @return The index of the root node.
   */
  virtual long GetRoot(void) const { return root; }
  /**
   * @brief This is a virtual function that does nothing in this class.
   * It is meant to be overridden in derived classes to delete extra data associated with a node.
   * @sergeilkp
   *
   * @param l The index of the node.
   */
  virtual void DeleteXtra(long) {};
  /**
   * @brief Delete all items in the list.
   *
   * @param cL If true, the data in the list will be deleted.
   */
  virtual void DeleteAll(bool cL) {
    Clear(cL);
    DeleteObject(dataList);
  }

  /**
   * @brief Count the number of items in the list.
   *
   * @return The number of items in the list.
   */
  unsigned long countitems(void) const;

  /**
   * @brief Find an item in the list.
   *
   * @param brc The item to find.
   * @return The index of the item if found, otherwise -1.
   */
  long Find(BaseRefConst) const;
  /**
   * @brief Find an item in the list and return the path to it.
   *
   * @param brc The item to find.
   * @param sl A _SimpleList to store the path to the item.
   * @return The index of the item if found, otherwise -1.
   */
  long Find(BaseRefConst, _SimpleList &) const;

  /**
   * @brief Find a long integer in the list.
   * This is a shortcut function that avoids calling ::Compare.
   *
   * @param l The long integer to find.
   * @return The index of the long integer if found, otherwise -1.
   */
  long FindLong(long) const;
  /**
   * @brief Find the best match for an item in the list.
   * "Best" is defined as the largest element smaller than or equal to the search key.
   *
   * @param brc The item to find the best match for.
   * @param l The index of the best match will be stored here.
   * @return The comparison result between the best match and the search key.
   */
  hyComparisonType FindBest(BaseRefConst, long &) const;

  /**
   * @brief Get the next item in the list in sorted order.
   *
   * @param l The index of the current item.
   * @param sl A _SimpleList to store the path to the next item.
   * @return The index of the next item.
   */
  long Next(long, _SimpleList &) const;
  /**
   * @brief Get the previous item in the list in sorted order.
   *
   * @param l The index of the current item.
   * @param sl A _SimpleList to store the path to the previous item.
   * @return The index of the previous item.
   */
  long Prev(long, _SimpleList &) const;
  /**
   * @brief Get the first item in the list in sorted order.
   *
   * @return The index of the first item.
   */
  long First(void) const;
  /**
   * @brief Get the last item in the list in sorted order.
   *
   * @return The index of the last item.
   */
  long Last(void) const;

  /**
   * @brief Check if an index is valid.
   *
   * @param l The index to check.
   * @return true if the index is valid, false otherwise.
   */
  bool IsValidIndex(long) const;

  /**
   * @brief Get an item by its index in the sorted list.
   *
   * @param l The index of the item to get.
   * @return The index of the item in the data list.
   */
  long GetByIndex(const long);

  /**
   * @brief Insert an item into the list.
   *
   * @param br The item to insert.
   * @param l The index to insert at. Not used.
   * @param b1 If true, the item will be duplicated before insertion.
   * @param b2 If b1 is false and the item already exists, the item will be deleted.
   * @return The index of the inserted item.
   */
  long Insert(BaseRef, long = 0, bool = true, bool = false);

  /**
   * @brief Insert a number into the list.
   *
   * @param v The number to insert.
   * @return The index of the inserted number.
   */
  long InsertNumber(long v) { return Insert((BaseRef)v); }

  /**
   * @brief Retrieve an item from the list by its index in the data list.
   *
   * @param l The index of the item to retrieve.
   * @return The item at the given index.
   */
  BaseRef Retrieve(long) const;
  /**
   * @brief Retrieve a long integer from the list by its index in the data list.
   *
   * @param l The index of the item to retrieve.
   * @return The long integer at the given index.
   */
  long RetrieveLong(long) const;

  /**
   * @brief Delete an item from the list.
   *
   * @param brc The item to delete.
   * @param b If true, the data will be deleted.
   */
  void Delete(BaseRefConst, bool = false);
  /**
   * @brief Check the consistency of the AVL tree.
   * This is a debugging function.
   */
  void ConsistencyCheck(void);

  /**
   * @brief Get the keys of the list.
   * The keys are the items in the list.
   *
   * @return A _List object containing the keys.
   */
  const _List Keys(void) const;
};

/**
 * @brief A utility function to populate and sort a list.
 *
 * @tparam AGGREGARTOR A function object that takes an _AVLList& as an argument and populates it.
 * @param agg The aggregator function.
 * @return A _SimpleList containing the sorted indices of the populated list.
 *
 * @example
 * \code{.cpp}
 * _SimpleList sorted_indices = PopulateAndSort([](_AVLList& list) {
 *   list.Insert(new _String("C"));
 *   list.Insert(new _String("A"));
 *   list.Insert(new _String("B"));
 * });
 * // sorted_indices will contain the indices of the sorted list, e.g., {1, 2, 0}
 * \endcode
 */
template <typename AGGREGARTOR>
_SimpleList const PopulateAndSort(AGGREGARTOR agg) {
  _SimpleList indexer;
  _AVLList avl(&indexer);
  agg(avl);
  avl.ReorderList();
  return indexer;
}

#endif
