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

#ifndef __ASSOCIATIVELIST__
#define __ASSOCIATIVELIST__

#include "mathobj.h"
#include "avllistxl.h"
#include "elementarycommand.h"

class _AssociativeList : public _MathObject {

public:
  _AssociativeList(void);
  _AssociativeList(_PMathObj); 

  /* construct an object from a list of key : value pairs 
    supplied as formulas in the _SimpleList argument */

  virtual ~_AssociativeList(void) {}

  //SLKP 20090803
  // Parse the list represented as {"key": value, ...}
  // the boolean argument is supplied to request reporting/suppression
  // of error messages returns true on successful parse
  bool ParseStringRepresentation(_String &, bool = true,
                                 _VariableContainer * = nil);

  virtual BaseRef toStr(void);
  virtual _PMathObj
  Execute(long opCode, _PMathObj = nil, _PMathObj = nil,
          _hyExecutionContext *context = _hyDefaultExecutionContext);
  virtual BaseRef makeDynamic(void);
  virtual _PMathObj Compute(void);

  // 20100907: SLKP
  // A simple function to merge two lists;the combined list will have the key
  // set equal to the union of the two input key sets if there are
  // conflicting values for a given key, an undefined value will be stored
  // in for the corresponding key
  virtual void Merge(_PMathObj);
  virtual void Duplicate(BaseRef);
  _PMathObj MAccess(_PMathObj);

  //perform a function call (ID stored in the first argument) having performed
  //[an optional] conditional check on the associated key (either empty for noop
  //or a function ID) Both functional IDs MUST be defined and take TWO and ONE
  //argumens respectively returns the number of items processed
  _PMathObj MIterator(_PMathObj, _PMathObj);
  _PMathObj GetByKey(_String &, long);
  _PMathObj GetByKey(_String &);
  _PMathObj GetByKey(long, long);
  void DeleteByKey(_PMathObj);
  _PMathObj MCoord(_PMathObj);

  // SLKP 20100811: see the comment for _Matrix::MStore
  void MStore(_PMathObj, _PMathObj, bool = true, long = HY_OP_CODE_NONE);
  void MStore(_String, _PMathObj, bool = true);
  void MStore(_String, _String);
  virtual unsigned long ObjectClass(void) { return ASSOCIATIVE_LIST; }
  _List *GetKeys(void);
  void FillInList(_List &);
  _String *Serialize(_String &);

  //Traverse the dictionary, cast each value into a float and return their sum.
  //Note that matrices and dictionary values will be processed recursively,
  //i.e. "Sum" will be called on them.
  //All values that cannot be cast to a float will be treated as 0.
  //@return The sum of all dictionary elements.
  _PMathObj Sum(void);

  _AVLListXL avl;

private:
  _List theData;
};

void InsertStringListIntoAVL(_AssociativeList *, _String, _SimpleList &,
                             _List &);
void InsertVarIDsInList(_AssociativeList *, _String, _SimpleList &);

#endif
