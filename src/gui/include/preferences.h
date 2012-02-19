/*
HyPhy - Hypothesis Testing Using Phylogenies.
This file implements shared and platform specific
(via ifdefs) functions for reading/writing and
setting preferences.

Written by SL Kosakovsky Pond
June 8, 2007
Copyright (C) 1997-2006

Primary Development:
  Sergei L Kosakovsky Pond (sergeilkp@mac.com)
Significant contributions from:
  Spencer V Muse (muse@stat.ncsu.edu)
  Simon DW Frost (sdfrost@ucsd.edu)
  Art FY Poon    (apoon@biomail.ucsd.edu)

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



#ifndef __HYPHY_PREFERENCES__

#define  PREFITEM_TEXTBOX 0
#define  PREFITEM_POPUP   1
#define  RECENT_FILE_ITEMS 10

#include    "list.h"
#include    "hy_strings.h"
#ifndef     __HEADLESS__
#include    "HYBaseGUI.h"
#endif

/* read HyPhy preferences from a platform - dependent location */
void      ReadPreferences       (void);

/* write HyPhy preferences to a platform - dependent location */
void      WritePreferences      (void);

/* apply relevant settings before executing a batch file */
void      ApplyPreferences      (void);

/* immediately apply those preferences which can be set
 interactively; e.g. console font */
void      SetPreferences        (void);

extern    bool showDialogAtStartup,
          doAutoConsoleMove ;

extern    _List                globalPreferencesList,
          recentPaths,
          recentFiles;



#ifndef     __HEADLESS__

extern    _HYRect              consolePositionRectangle;

#endif



extern    _String              recentFilesList;



void      AddStringToRecentMenu (_String&, _String&);

/* add a string/path to the recent analyses menu */



void      AddItemToPreferences   (long,long,_String,_String,_String,_List*,_List&,bool);

/* a function used to append items to a preferences list */



char      AutoOpenTreeWindow     (void);

void      SetShowDialogAtStartup (bool);



#endif





//EOF

