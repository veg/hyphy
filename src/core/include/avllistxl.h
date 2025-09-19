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

#ifndef _AVLLISTXL_
#define _AVLLISTXL_

//#pragma once
#include "list.h"
#include "avllist.h"

#define  MEMORYSTEP 8

//_____________________________________________________________________________

/**
 * @brief An AVL list with an extra data member (_List) for each node.
 */
class _AVLListXL: public _AVLList {

public:

    /**
     * @brief Construct a new _AVLListXL object.
     *
     * @param sl A _SimpleList to construct from. The new list will contain pointers to the data in the _SimpleList, but will not own the data. The _SimpleList will be cleared after the new _AVLListXL is created.
     */
    _AVLListXL(_SimpleList*);
    /**
     * @brief Construct a new _AVLListXL object.
     *
     * @param sl A _SimpleList to construct from. The new list will contain pointers to the data in the _SimpleList, but will not own the data. The _SimpleList will be cleared after the new _AVLListXL is created.
     * @param l A _List to use as the extra data.
     */
    _AVLListXL(_SimpleList*, _List&&);
    /**
     * @brief Get the extra data associated with a key.
     *
     * @param l The index of the key.
     * @return The extra data.
     */
    BaseRef GetXtra(long) const;

    /**
     * @brief Set the extra data associated with a key.
     *
     * @param l The index of the key.
     * @param br The extra data.
     * @param b If true, the extra data will be duplicated.
     */
    void SetXtra(long,BaseRef,bool);

    /**
     * @brief Destroy the _AVLListXL object.
     */
    virtual ~_AVLListXL(void){}
    /**
     * @brief Convert the list to a string representation.
     *
     * @param ul The number of spaces to use for padding.
     * @return A _String object representing the list.
     */
    virtual BaseRef toStr(unsigned long = 0UL);
    /**
     * @brief Get the extra data associated with a key.
     *
     * @param brc The key to search for.
     * @return The extra data associated with the key, or nil if the key is not found.
     */
    virtual BaseRef GetDataByKey (BaseRefConst) const;

    /**
     * @brief Insert data into the list.
     *
     * @param br The data to insert.
     * @param l The extra data to associate with the key.
     * @param b If true, the data will not be copied.
     * @return The index of the inserted item.
     */
    virtual long InsertData(BaseRef, long,bool);
    /**
     * @brief Clear the list.
     * If shallow is true, the data in the list will not be deleted.
     *
     * @param b Whether to perform a shallow clear.
     */
    virtual void Clear(bool = false);
    /**
     * @brief This is a virtual function that does nothing in this class.
     * It is meant to be overridden in derived classes to delete extra data associated with a node.
     * @sergeilkp What is this extra data?
     *
     * @param l The index of the node.
     */
    virtual void DeleteXtra(long);
    /**
     * @brief Update the value of an item in the list.
     *
     * @param br1 The new key.
     * @param br2 The new extra data.
     * @param b1 If true, the new key will be duplicated.
     * @param b2 If true, the new extra data will be duplicated.
     * @return The index of the updated item.
     */
    virtual long UpdateValue (BaseRef, BaseRef, bool = false, bool = true);
    
    /**
     * @brief Push a key-value pair to the list. The key is copied.
     *
     * @param key The key to push.
     * @param br The value to push.
     * @return A reference to this list.
     */
    _AVLListXL&  PushPairCopyKey  (_String const key, BaseRef);
    /**
     * @brief Push a key-value pair to the list. The key is not copied.
     *
     * @param key The key to push.
     * @param br The value to push.
     * @return A reference to this list.
     */
    _AVLListXL&  PushPair         (_String* key, BaseRef);
    

    _List xtraD;

};


//_____________________________________________________________________________

#endif
