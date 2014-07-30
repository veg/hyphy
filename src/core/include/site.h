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
#include "hy_string_buffer.h"
#include <stdlib.h>

class _DataSetFilter;

class _Site : public _StringBuffer {

public:
  //does nothing
  _Site(void);
  // data constructor
  _Site(_String &);
  // data constructor
  _Site(char);
  // reference to another site constructor
  _Site(long);

  //destructor
  virtual ~_Site(void);

  // mark this site as complete
  void Complete(void);

  //virtual BaseRef makeDynamic(void);
  void duplicate(BaseRef s) {
    _StringBuffer::duplicate(s);
    refNo = -1;
  }
  virtual void Clear(void);

  // XXX maybe remove these
  //void PrepareToUse(void); // decompress the site preparing for intensive use
  //void Archive(void);      // archive the site for later use

  long GetRefNo(void) { return refNo < 0 ? -refNo - 2 : refNo - 2; }

  // Complete and IsComplete are never used, but they don't make sense as
  // implemented before refactoring, where complete is when refNo is < 0. I
  // believe this because complete makes refNo positive always.
  bool IsComplete(void) { return refNo > 0; }

  void SetRefNo(long r) { refNo = -r - 2; }

private:

  long refNo; // if this site contains a reference to another one
  // if refNo is negative, then shows whether the definition of this datatype
  // has been completed
};

/*
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
*/

#endif
