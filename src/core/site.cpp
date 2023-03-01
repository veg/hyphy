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

#include <ctype.h>
#include <math.h>
#include <string.h>

#include "avllistxl_iterator.h"
#include "function_templates.h"
#include "hbl_env.h"
#include "likefunc.h"
#include "site.h"

#include "string_file_wrapper.h"

#include "global_object_lists.h"
#include "global_things.h"
#include "hbl_env.h"

using namespace hyphy_global_objects;
using namespace hy_global;

#define DEFAULT_SITE_ALLOCATION 16L



//_________________________________________________________

_Site::_Site(void) : _StringBuffer (DEFAULT_SITE_ALLOCATION) { refNo = -1; }

//_________________________________________________________

_Site::_Site(unsigned long space, long ref_no) : _StringBuffer (space) { refNo = ref_no; }


//_________________________________________________________
_Site::_Site(_String const &s) : _StringBuffer (s.length()) {
  refNo = -1;
  (*this) << &s;
}

//_________________________________________________________
_Site::_Site(char s) : _StringBuffer (DEFAULT_SITE_ALLOCATION) {
  refNo = -1;
  (*this) << s;
}

//_________________________________________________________
_Site::_Site(long s) {
  SetRefNo(s);
}

//_________________________________________________________
_Site::~_Site(void) {}

//_________________________________________________________
void _Site::Complete(void) {
  refNo = refNo < 0L ? -refNo : refNo;
}

//__________________________________________________________________________________
void * _Site::operator new (size_t size) {
    return MemAllocate (size);
}

//__________________________________________________________________________________
void  _Site::operator delete (void * p) {
    free (p);
}
