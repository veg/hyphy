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

#pragma once

#include <QtCore>
#include <QtGui>
#include "defines.h"
#include "batchlan.h"
#include "ui_hyphymain.h"
#include "hyphymain.h"
#include "hy_strings.h"

extern HyphyMain* _hyPrimaryConsoleWindow;

_String _hyQTFileDialog (_String caption, _String defaultFileName, bool isWrite);

/**
 * Display a platform appropriate file selection dialog
 * @param caption: dialog caption string
 * @param defaultFileName (for write dialogs): the default file name to display
 * @param isWrite: whether or not to display a read/write dialog
 * @return full path to the selected file name; empty if canceled.
 */


long  HandleListSelection (_List& data, _SimpleList& choices, _SimpleList& validChoices, _String titleInfo, _SimpleList& selections, long fixedLength, Ptr prt);
long  HandleListSelection (_List& data, _String titleInfo, Ptr prt);




void    SetStatusLineUser(_String s);
void    SetStatusBarValue (long l, _Parameter max, _Parameter rate);
bool    Get_a_URL (_String& urls, _String* fileName);

void ClearStatusLine();

void    DoApplicationSettings (void);
