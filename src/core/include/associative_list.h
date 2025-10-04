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

#ifndef __ASSOCIATIVE_LIST__
#define __ASSOCIATIVE_LIST__

#include "avllistx.h"
#include "avllistxl_iterator.h"
#include "hy_strings.h"
#include "trie.h"
#include "variablecontainer.h"

/**
 * @struct _associative_list_key_value
 * @brief A simple struct to hold a key-value pair for an associative list.
 */
struct _associative_list_key_value {
  const char *key;      /**< The key of the pair. */
  HBLObjectRef payload; /**< The value (payload) of the pair. */
};

/**
 * @class _AssociativeList
 * @brief An associative list (dictionary/map) implementation, mapping string
 * keys to HBLObjectRef values.
 *
 * This class provides a dictionary-like data structure with string keys. It is
 * implemented using an AVL tree for efficient lookups, insertions, and
 * deletions.
 */
class _AssociativeList : public _MathObject {
public:
  /**
   * @brief Default constructor.
   */
  _AssociativeList(void);
  /**
   * @brief Virtual destructor.
   */
  virtual ~_AssociativeList(void) {}

  /**
   * @brief Parse the list represented as `{"key": value, ...}`.
   * @param s The string representation to parse.
   * @param context The formula parsing context.
   * @return Returns `true` on successful parse, `false` otherwise.
   */
  bool ParseStringRepresentation(_String &, _FormulaParsingContext &);

  /**
   * @brief Converts the associative list to its string representation.
   * @param ulong An optional parameter for formatting.
   * @return A `BaseRef` containing the string representation.
   */
  virtual BaseRef toStr(unsigned long = 0UL);
  /**
   * @brief Execute a single operation on the list.
   * @param opCode The operation code to execute.
   * @param arguments A list of arguments for the operation.
   * @param context The execution context.
   * @param cache An optional cache for the operation.
   * @return An `HBLObjectRef` representing the result of the operation.
   */
  virtual HBLObjectRef
  ExecuteSingleOp(long opCode, _List *arguments = nil,
                  _hyExecutionContext *context = _hyDefaultExecutionContext,
                  HBLObjectRef cache = nil);
  /**
   * @brief Creates a dynamic copy of the object.
   * @return A `BaseRef` pointing to the new dynamic copy.
   */
  virtual BaseRef makeDynamic(void) const;
  /**
   * @brief Computes the object. For an associative list, this typically returns
   * a reference to itself.
   * @return An `HBLObjectRef` to the computed object.
   */
  virtual HBLObjectRef Compute(void);
  /**
   * @brief Clears all elements from the associative list.
   */
  void Clear(void);
  /**
   * @brief Merges another associative list into this one.
   * The combined list will have a key set equal to the union of the two input
   * key sets. If there are conflicting values for a given key, an undefined
   * value will be stored for that key.
   * @param The `HBLObjectRef` to merge with.
   */
  virtual void Merge(HBLObjectRef);

  /**
   * @brief Checks for equality with another HBLObjectRef.
   * @param The object to compare with.
   * @return `true` if the objects are equal, `false` otherwise.
   */
  virtual bool Equal(HBLObjectRef);

  /**
   * @brief Duplicates the content of another object into this one.
   * @param The object to duplicate.
   */
  virtual void Duplicate(BaseRefConst);
  /**
   * @brief Randomize key-value assignments, sampling values with or without
   * replacement.
   * @param with_replacement If truthy, sample with replacement.
   * @param cache An optional cache for the operation.
   * @return A new `_AssociativeList` with randomized assignments.
   */
  HBLObjectRef Random(HBLObjectRef, HBLObjectRef cache);
  /**
   * @brief Access an element by key.
   * @param key The key of the element to access.
   * @param cache An optional cache for the operation.
   * @return The value associated with the key.
   */
  HBLObjectRef MAccess(HBLObjectRef, HBLObjectRef cache);

  /**
   * @brief Iterates over the list and performs a function call on each item.
   *
   * This method performs a function call (ID stored in the first argument) on
   * each key-value pair, optionally after a conditional check on the key
   * (function ID in the second argument, or empty for no-op). The main function
   * must take two arguments (key, value), and the conditional function must
   * take one argument (key).
   *
   * @param function_id The ID of the function to call on each item.
   * @param condition_id The ID of the conditional function to check against the
   * key.
   * @param cache An optional cache for the operation.
   * @return The number of items processed.
   */
  HBLObjectRef MIterator(HBLObjectRef, HBLObjectRef, HBLObjectRef cache);

  /**
   * @brief Get a value by its key.
   * @param key The key to look up.
   * @param index An additional index parameter (usage may vary).
   * @return The value associated with the key, or `nil` if not found.
   */
  HBLObjectRef GetByKey(_String const &, long) const;
  /**
   * @brief Get a value by its key, throwing an exception if not found.
   * @param key The key to look up.
   * @param index An additional index parameter.
   * @return The value associated with the key.
   */
  HBLObjectRef GetByKeyException(_String const &, long) const;
  /**
   * @brief Get a value by its key.
   * @param key The key to look up.
   * @return The value associated with the key, or `nil` if not found.
   */
  HBLObjectRef GetByKey(_String const &) const;
  /**
   * @brief Get a value by its index in the underlying storage.
   * @param index The index of the key.
   * @param index2 An additional index parameter for the value.
   * @return The value at the specified location.
   */
  HBLObjectRef GetByKey(long, long) const;
  /**
   * @brief Delete a key-value pair by key.
   * @param key The key of the pair to delete, passed as an `HBLObjectRef`.
   */
  void DeleteByKey(HBLObjectRef);
  /**
   * @brief Delete a key-value pair by key.
   * @param key The key of the pair to delete.
   */
  void DeleteByKey(_String const &);
  /**
   * @brief Get the coordinate (index) of a key.
   * @param key The key to find.
   * @param cache An optional cache.
   * @return The index of the key, or -1 if not found.
   */
  HBLObjectRef MCoord(HBLObjectRef, HBLObjectRef);
  /**
   * @brief Store a key-value pair.
   * @param key The key.
   * @param value The value.
   * @param ref_count Whether to increment the reference count of the value.
   * @param op_code An optional operation code for assignment (e.g., `+=`).
   * @return `true` if the store was successful.
   */
  bool MStore(_String *, HBLObjectRef, bool = true, long = HY_OP_CODE_NONE);
  /**
   * @brief Store a key-value pair.
   * @param key The key, passed as an `HBLObjectRef`.
   * @param value The value.
   * @param ref_count Whether to increment the reference count of the value.
   * @param op_code An optional operation code for assignment.
   */
  void MStore(HBLObjectRef, HBLObjectRef, bool = true, long = HY_OP_CODE_NONE);

  /**
   * @brief Store a key-value pair.
   * @param key The key.
   * @param value The value.
   * @param ref_count Whether to increment the reference count of the value.
   */
  void MStore(const _String &, HBLObjectRef, bool = true);

  /**
   * @brief A convenience build-out function to push key-value pairs.
   * This operator adds a reference count to the payload.
   * @param pair The key-value pair to add.
   * @return A reference to this `_AssociativeList`.
   */
  _AssociativeList &operator<<(_associative_list_key_value pair);
  /**
   * @brief A convenience build-out function to push key-value pairs.
   * This operator does NOT add a reference count to the payload.
   * @param pair The key-value pair to add.
   * @return A reference to this `_AssociativeList`.
   */
  _AssociativeList &operator<(_associative_list_key_value pair);

  /**
   * @brief Store a key-value pair where the value is a string.
   * @param key The key.
   * @param value The string value.
   */
  void MStore(const _String &, const _String &);
  /**
   * @brief Returns the object's class type.
   * @return The class type identifier `ASSOCIATIVE_LIST`.
   */
  virtual unsigned long ObjectClass(void) const { return ASSOCIATIVE_LIST; }
  /**
   * @brief Get all keys as a `_List`.
   * @return A `_List` containing all the keys.
   */
  _List *GetKeys(void) const;
  /**
   * @brief Get the smallest key that can be parsed as a number.
   * @return The smallest numerical key as a `_String`, or `nil` if none exist.
   */
  _String *GetSmallestNumericalKey(void) const;
  /**
   * @brief Fill a given `_List` with the values from this associative list.
   * @param list_to_fill The list to be filled with values.
   */
  void FillInList(_List &);
  /**
   * @brief Get the number of items in the list.
   * @return The number of key-value pairs.
   */
  unsigned long Length(void) const { return avl.countitems(); }
  /**
   * @brief Serialize the object to a `_StringBuffer`.
   * @param ulong An optional parameter for serialization options.
   * @return A `_StringBuffer` containing the serialized data.
   */
  _StringBuffer *Serialize(unsigned long) const;
  /**
   * @brief Get the number of items in the list.
   * @return The number of key-value pairs.
   */
  unsigned long countitems(void) const { return avl.countitems(); }

  /**
   * @brief Traverse the dictionary, and store keys and corresponding values as
   * two lists: key[x] and value[x], where x is 0..|Dict|-1.
   * @param list_to_fill The list to contain key and value lists (at indices 0
   * and 1, respectively).
   */
  void KeysValuesAsLists(_List &);

  /**
   * @brief Obtain an iterator over list elements.
   * @return An `AVLListXLIterator` for this list.
   */
  AVLListXLIterator ListIterator(void) { return AVLListXLIterator(&avl); }

  /**
   * @brief Traverse the dictionary, cast each value into a float, and return
   * their sum. Note that matrices and dictionary values will be processed
   * recursively (i.e., `Sum` will be called on them). All values that cannot be
   * cast to a float will be treated as 0.
   * @param cache An optional cache for the operation.
   * @return The sum of all dictionary elements as an `HBLObjectRef`.
   */
  HBLObjectRef Sum(HBLObjectRef cache);
  /**
   * @brief Traverse the dictionary and find the extreme (min/max) numeric
   * value. Returns `{"key": key, "value": min/max}`. All values that cannot be
   * cast to a float will be IGNORED. If no valid numbers are found, "key" will
   * be `None`, and min/max will be +/-Inf.
   * @param do_minimum If `true`, find the minimum value; otherwise, find the
   * maximum.
   * @param cache An optional cache for the operation.
   * @return An `_AssociativeList` containing the key and value of the extreme
   * element.
   */
  HBLObjectRef ExtremeValue(bool do_mimimum, HBLObjectRef cache) const;

  /**
   * @brief A convenience function to get a numeric value by key.
   * @param key The key to look up.
   * @return The numeric value.
   * @throws A `const _String` error if the key is not found or the value is not
   * numeric.
   */
  hyFloat GetNumberByKey(const _String &key) const;

  /**
   * @brief A convenience function to get a numeric value by key, with a default
   * value.
   * @param key The key to look up.
   * @param default_value The value to return if the key is not found or the
   * value is not numeric.
   * @return The numeric value or the default.
   */
  hyFloat GetNumberByKeyDefault(const _String &key, hyFloat = 0.) const;

private:
  _AVLListXL avl; /**< The underlying AVL tree for storing key indices. */
  _List theData;  /**< The list holding the actual value objects. */
};

/**
 * @brief Insert a list of strings into an associative list under a given key.
 * @param target_list The associative list to modify.
 * @param key The key under which to insert the list.
 * @param string_indices A list of indices for the strings.
 * @param string_values The list of actual string values.
 */
void InsertStringListIntoAVL(_AssociativeList *, _String const &,
                             _SimpleList const &, _List const &);
/**
 * @brief Insert variable IDs into an associative list.
 * @param target_list The associative list to modify.
 * @param key The key for the insertion.
 * @param var_indices A list of variable indices.
 */
void InsertVarIDsInList(_AssociativeList *, _String const &,
                        _SimpleList const &);

#endif
