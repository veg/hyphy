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

#ifndef     __ASSOCIATIVE_LIST__
#define     __ASSOCIATIVE_LIST__

#include "hy_strings.h"
#include "avllistx.h"
#include "variablecontainer.h"
#include "trie.h"



/*__________________________________________________________________________________________________________________________________________ */

struct _associative_list_key_value {
  const char * key;
  _PMathObj  payload;
};

/*__________________________________________________________________________________________________________________________________________ */

class           _AssociativeList: public _MathObject {
public:
    _AssociativeList                    (void);
    virtual ~_AssociativeList           (void) {}

    bool    ParseStringRepresentation   (_String&, _FormulaParsingContext&);
    /* SLKP 20090803

        Parse the list represented as
            {"key": value, ...}

        the boolean argument is supplied to request reporting/suppression of error messages

        returns true on successful parse

     */

    virtual BaseRef     toStr           (unsigned long = 0UL);
    virtual _PMathObj   ExecuteSingleOp (long opCode, _List* arguments = nil, _hyExecutionContext* context = _hyDefaultExecutionContext);
    // execute this operation with the list of Args
    virtual BaseRef     makeDynamic     (void) const;
    virtual _PMathObj   Compute         (void);
    void                Clear           (void);
    virtual void        Merge           (_PMathObj);
    /* 20100907: SLKP
            A simple function to merge two lists;
            the combined list will have the key set equal to the union of the two input key sets
            if there are conflicting values for a given key, an undefined value will be stored in
            for the corresponding key


     */

    virtual void        Duplicate       (BaseRefConst);
    _PMathObj           MAccess         (_PMathObj);

    _PMathObj           MIterator       (_PMathObj, _PMathObj);
    /* perform a function call (ID stored in the first argument)
       having performed [an optional] conditional check on the associated key (either empty for noop or a function ID)
       Both functional IDs MUST be defined and take TWO and ONE arguments respectively

       returns the number of items processed
    */

    _PMathObj           GetByKey        (_String const&, long) const;
    _PMathObj           GetByKey        (_String const&) const;
    _PMathObj           GetByKey        (long, long) const;
    void                DeleteByKey     (_PMathObj);
    _PMathObj           MCoord          (_PMathObj);
    void                MStore          (_PMathObj, _PMathObj, bool = true, long = HY_OP_CODE_NONE);
    // SLKP 20100811: see the comment for _Matrix::MStore

    void                MStore          (const _String&  , _PMathObj, bool = true);
  
    /* a convenience build-out function to push key-value pairs
       << adds a reference count to the payload
    */
  
    _AssociativeList &  operator <<     (_associative_list_key_value pair);
    _AssociativeList &  operator <     (_associative_list_key_value pair);
  
    void                MStore          (const _String&  , const _String&);
    virtual unsigned long        ObjectClass     (void)      {
        return ASSOCIATIVE_LIST;
    }
    _List*              GetKeys         (void) const;
    void                FillInList      (_List&);
    unsigned long       Length          (void) const {
      return avl.countitems();
    }
    _StringBuffer *            Serialize       (unsigned long) const;
    unsigned   long     countitems      (void) const {
        return avl.countitems();
    }
    
    /**
     * Traverse the dictionary, cast each value into a float and return their sum.
     * Note that matrices and dictionary values will be processed recursively, i.e. "Sum" will be called on them.
     * All values that cannot be cast to a float will be treated as 0.
     * @return The sum of all dictionary elements.
     */
    _PMathObj           Sum             (void);
    /**
     * Traverse the dictionary, and return { "key" : key, "value" : min / max over the list}
     * All values that cannot be cast to a float will be IGNORED.
     * If no valid numbers could be found, "key" will be None, and min/max will be an +/-Inf
     * @return The minimum or maximum numeric value and corresponding key
     */
    _PMathObj           ExtremeValue    (bool do_minimum) const;


private:
    _AVLListXL          avl;
    _List           theData;
};

void       InsertStringListIntoAVL  (_AssociativeList* , _String const&, _SimpleList const&, _List const&);
void       InsertVarIDsInList       (_AssociativeList* , _String const&, _SimpleList const&);


#endif
