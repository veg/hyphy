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

#include <site.h>
#include "helperfunctions.h"
#include "errorfns.h"

#include "hy_string_buffer.h"
#include <string.h> // for memcpy

#if !defined __UNIX__ || defined __HEADLESS__
#include "preferences.h"
#endif

//______________________________________________________________________________
_Site::_Site(void) : _StringBuffer(16UL) { ref_no = -1; }

//______________________________________________________________________________
_Site::_Site(_String &s) : _StringBuffer(s.s_length) {
  ref_no = -1;
  (*this) << &s;
}

//______________________________________________________________________________
_Site::_Site(char s) : _StringBuffer(16UL) {
  ref_no = -1;
  (*this) << s;
}

//______________________________________________________________________________
_Site::_Site(long s) {
  this->setRefNo(s);
}

//______________________________________________________________________________
_Site::~_Site(void) {}

//______________________________________________________________________________
void _Site::complete(void) {
  ref_no = ref_no < 0 ? -ref_no : ref_no;
}

void _Site::duplicate(BaseRef s) {
  _StringBuffer::duplicate(s);
  ref_no = -1;
}

//______________________________________________________________________________
// It isn't immediately obvious that resetting the ref_no is necessary, but
// it seems like it would be the expected behaviour
void _Site::clear(void) {
  if (s_data) {
    free(s_data);
  }
  this->initialize();
  ref_no = -1;
}

