/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2006
Primary Development:
  Sergei L Kosakovsky Pond (sergeilkp@mac.com)
Significant contributions from:
  Spencer V Muse (muse@stat.ncsu.edu)
  Simon DW Frost (sdfrost@ucsd.edu)
  Art FY Poon    (apoon@biomail.ucsd.edu)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

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


