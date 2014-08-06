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

#ifndef _GENSITE_
#define _GENSITE_
//#pragma once
#include "sequence.h"
#include "legacy_parser.h"
#include "simplelist.h"
#include "list.h"
#include "avllist.h"
#include "avllistx.h"
#include "avllistxl.h"
#include <stdlib.h>

#include "dataset.h"
#include "translationtable.h"

#define NUCLEOTIDEDATA 0
#define CODONDATA 1

/*
class _DataSetFilter;

class _Site : public _CString // compressible string
              {

public:
  _Site(void);
  //does nothing
  _Site(_String &);
  // data constructor
  _Site(char);
  // data constructor
  _Site(long);
  // reference constructor

  virtual ~_Site(void);
  //destructor

  void Complete(void); // mark this site as complete and compress it

  virtual BaseRef makeDynamic(void);
  virtual void Duplicate(BaseRef);
  virtual void Clear(void);

  void PrepareToUse(void); // decompress the site preparing for intensive use
  void Archive(void);      // archive the site for later use

  long GetRefNo(void) { return refNo < 0 ? -refNo - 2 : refNo - 2; }

  bool IsComplete(void) { return refNo < 0; }

  void SetRefNo(long r) { refNo = -r - 2; }

private:

  long refNo; // if this site contains a reference to another one
  // if refNo is negative, then shows whether the definition of this datatype
  // has been completed
};
*/

extern _TranslationTable defaultTranslationTable;

void ReadNextLine(FILE *fp, _String *s, FileState *fs, bool append = false,
                  bool upCase = true);
_DataSet *ReadDataSetFile(FILE *, char = 0, _String * = nil, _String * = nil,
                          _String * = nil,
                          _TranslationTable * = &defaultTranslationTable);
void fillDefaultCharTable(void);
void printFileResults(_DataSet *);
void printDSFilter(_DataSetFilter *d);

bool StoreADataSet(_DataSet *, _String *);

extern _String dataFileTree, dataFileTreeString, nexusFileTreeMatrix,
    dataFilePartitionMatrix, defaultLargeFileCutoff, nexusBFBody;

extern _DataSet *lastNexusDataMatrix;

#endif
