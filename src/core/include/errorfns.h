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

#ifndef _HERROR_
#define _HERROR_
//#pragma once
#include "baseobj.h"

#define  _HYNOERROR 1

//successful execution marker

extern  int      gError;
extern  bool  isFixable;


bool gStatus(void);
// global Error Status

bool isError (long);
// checks ToolBox function return, and if it is abnormal,
// sets the error status to true;

//bool acknError (bool);
// function acknoledges Error with either fixable (True) or terminal (False)
// flag, the latter causes acknError to terminate the program
// writes error number to stderror

void acknError (const char*);
// function acknoledges Error with either fixable (True) or terminal (False)
// flag, the latter causes acknError to terminate the program
// writes error number to stderror

void    warnError (const char*);
// warns user of current error with a supplied string
// writes to stderror

void    warnError (long);
// warns user of current error with a built in string
// writes to stderror

void    flagError (long);
// reports the text of current error to the user with a built in string
// writes to stderror
// terminates execution of current BF

void*   checkPointer  (void*);

#if !defined __UNIX__  && !defined __HEADLESS__
extern bool skipWarningMessages;
#endif

#endif


