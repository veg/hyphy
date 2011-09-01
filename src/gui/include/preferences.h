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

#ifndef __HYPHY_PREFERENCES__

#define  PREFITEM_TEXTBOX 0
#define  PREFITEM_POPUP   1

#define  RECENT_FILE_ITEMS 10

#include    "hy_lists.h"
#include    "hy_strings.h"

#ifndef     __HEADLESS__
#include    "HYBaseGUI.h"
#endif

void      ReadPreferences       (void);
/* read HyPhy preferences from a platform - dependent location */

void      WritePreferences      (void);
/* write HyPhy preferences to a platform - dependent location */

void      ApplyPreferences      (void);
/* apply relevant settings before executing a batch file */

void      SetPreferences        (void);
/* immediately apply those preferences which can be set
     interactively; e.g. console font */

extern    bool                 showDialogAtStartup,
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
