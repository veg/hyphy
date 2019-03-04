/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-now
Core Developers:
   Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
   Art FY Poon    (apoon42@uwo.ca)
   Steven Weaver (sweaver@temple.edu)

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

#pragma once

#include "dataset_filter.h"
#include "matrix.h"

class _DataSetFilterNumeric : public _DataSetFilter {

public:
  _DataSetFilterNumeric(void) {}
  _DataSetFilterNumeric(_Matrix *, _List &, _DataSet *, long);
  virtual ~_DataSetFilterNumeric(void);

  virtual bool IsNormalFilter(void) const { return false; }

  virtual BaseRef makeDynamic(void) const;
  virtual unsigned long GetDimension(bool) const { return dimension; }

  hyFloat *getProbabilityVector(long, long, long = 0);
  virtual bool CompareTwoSites(unsigned long, unsigned long,
                               unsigned long) const;

  long shifter, categoryShifter, categoryCount;
  _Matrix probabilityVectors;
  // N x M dense matrix
  // N = spec count
  // M = unique sites * dimension
};
