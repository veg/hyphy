/*



HyPhy - Hypothesis Testing Using Phylogenies.



This file defined shared function used by 'mains' in the GUI



Written by SL Kosakovsky Pond

June 11, 2007



Copyright (C) 1997-2007

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



#ifndef _HY_SHARED_MAINS_
#define _HY_SHARED_MAINS_



#include    "likefunc.h"
#include    "preferences.h"
#include    <time.h>



void        PrepareToExecuteBatchFile    (void);

// set up menus and such prior to executing a batch file

bool        ExecuteBatchFile             (void);

// execute the file with the path name stored in argFileName

void        DoneWithExecutionOfBatchFile (bool = true);

// set up menus and such following the execution of a batch file



void        updateTimerF                 (_String&, time_t);

// converts a long seconds field into hrs:mins:secs format



bool        OpenBatchFile                (bool = true, _String* = nil);

// propmt for a file (if flag = true), read it and set paths



void        ReadInTemplateFiles          (void);

// scan and read in TemplateBatchFiles and result processors



long        SelectATemplate              (void);

// choose a standard analysis to run



void        RunStandardAnalyses          (void);

// Choose and run a standard analysis



_String     MatrixExpCounter             (void);

// run a very simple benchmark (shown in the About Box)



void        SpoolFile                    (void);

// display a text file (path in argFileName) in the console window



void        RunTemplate                  (long);



void        ExecuteAPostProcessor        (_String);



// execute a post-processing module



void        ReportAnalysisAsFinished     (_String, bool = false);

// reflect the fact that an analysis has finished

// if appropriate for the platform



/* GLOBALS */



extern      _String*                     
                        argFileName,
                        errorFileName,
                        messageFileName,
                        menuSeparator;



extern      bool                         isSuspended,
                                         hasTemplates,
                                         highLevelQuit,
                                         isRerunAvailable,
                                         updateTimer,

            addToRecent,

            echoPaused,

            calculatorMode;



extern      _ExecutionList               ex;



extern      _SimpleList                  windowPtrs,

            windowObjects,

            treeIDReferences,

            windowObjectRefs;



extern      _List                        availableTemplateFiles,

            availablePostProcessors;



#ifdef      __WINDOZE__

#define          UPDATE_TIMER         WM_USER + 4

#endif

#endif



//EOF

