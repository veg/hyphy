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


#include "HYDialogs.h"
#include "hy_strings.h"
#include "HYTreePanel.h"
#include "HYDataPanel.h"
#include "HYObjectInspector.h"
#include "HYConsoleWindow.h"
#include "HYEventTypes.h"

#include <gtk/gtk.h>
#include <gdk/gdk.h>
#include <unistd.h>
#include <time.h>
#include <stdio.h>
#include "batchlan.h"
#include "calcnode.h"
#include <unistd.h>
//#include <curses.h>
#include <termios.h>
#ifdef   __MP2__
#include <pthread.h>
#endif

#include "likefunc.h"
#include "preferences.h"
#include "HYSharedMain.h"

#ifndef __HYPHY_NO_CURL__
#define __HYPHYCURL__
#endif

#ifdef  __HYPHYCURL__
#include <curl/curl.h>
#endif

#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif


#include "HYWindow.h"
#include "HYUtils.h"


void    ProcessConfigStr (_String&);
long    DisplayListOfPostChoices (void);


bool    usePostProcessors = false,
        updateMode        = false,
        pipeMode          = false,
        logInputMode    = false;

#ifdef  __HYPHYMPI__
bool mpiParallelOptimizer   = false,
     mpiPartitionOptimizer     = false;
void mpiNormalLoop    (int, int, _String &);
void mpiOptimizerLoop (int, int);


#endif


/* GLOBALS */

double                  fontConversionFactor = 1.0;

PangoContext*           screenPContext;
GdkPixbuf*              tablePDMenuIcon,
                        *redButtonIcon,
                        *yellowButtonIcon,
                        *orangeButtonIcon,
                        *greenButtonIcon;

GdkCursor               *hSizeCursor,
                        *pickUpCursor,
                        *dropOffCursor;

clock_t                 timerStart = 0,
                        lastTimer  = 0;

//__________________________________________________________________________________
void    ProcessConfigStr (_String& conf)
{
    _String errMsg;
    for (long i=1; i<conf.sLength; i++) {
        switch (conf.sData[i]) {
        case 'p':
        case 'P': {
            usePostProcessors = true;
            break;
        }
        case 'c':
        case 'C': {
            calculatorMode = true;
            break;
        }
        case 'u':
        case 'U': {
            updateMode = true;
            break;
        }
        case 'l':
        case 'L': {
            logInputMode = true;
            break;
        }
        //case 'i':
        //case 'I':
        //{
        //pipeMode = true;
        //break;
        //}
        default: {
            errMsg = "Option ";
            errMsg = errMsg & conf.sData[i] & " is not valid and is ignored";
            ReportWarning (errMsg);
        }
        }
    }
}

//__________________________________________________________________

gboolean progressTimerFunction (Ptr* userData)
{
    if (updateTimer) {
        clock_t time_diff = clock()-lastTimer;

        if (time_diff>CLOCKS_PER_SEC) { // update dispay
            lastTimer+=time_diff;
            _String r;
            updateTimerF(r,clock()-timerStart);
            SetStatusLine (empty, empty, r, 0, HY_SL_TIMER);
        }
    }
    return true;
}

//__________________________________________________________________

gboolean GlobalQueueTimer (Ptr* userData)
{
    if (GlobalGUIEventQueue.lLength) {
        HandleGlobalQueueEvent ();
    }
    //if (updateTimer)
    //  SendMessage ((HWND)hyphyConsoleWindow->GetOSWindowData(), UPDATE_TIMER, 0, 0);

    return true;
}

//__________________________________________________________________

void SetUpStatusBarStuff (GtkWidget* aWindow)
{
    _String            fName = libDirectory & "GTKResources/striped.xpm";
    statusBarLayout          = pango_layout_new (screenPContext);
    statusBarFontDesc        = pango_font_description_new ();
    stripedFill              = gdk_pixmap_create_from_xpm (GDK_DRAWABLE(aWindow->window), NULL, NULL, fName.sData);
    stripedFillGC            = gdk_gc_new (GDK_DRAWABLE(aWindow->window));
    if (stripedFill) {
        gdk_gc_set_fill (stripedFillGC,GDK_TILED);
        gdk_gc_set_tile (stripedFillGC,stripedFill);
    } else {
        printf ("Failed to load a status bar .xpm from %s\n", fName.sData);
    }

    gdk_gc_set_line_attributes        (stripedFillGC, 1, GDK_LINE_SOLID, GDK_CAP_NOT_LAST, GDK_JOIN_MITER);
    GdkColor saveFG = {0,0,0,0};
    gdk_gc_set_foreground             (stripedFillGC, &saveFG);

    pango_font_description_set_family (statusBarFontDesc, statusBarFont.face.sData);
    pango_font_description_set_style  (statusBarFontDesc, (statusBarFont.style & HY_FONT_ITALIC) ? PANGO_STYLE_ITALIC : PANGO_STYLE_NORMAL);
    pango_font_description_set_weight (statusBarFontDesc, (statusBarFont.style & HY_FONT_BOLD) ? PANGO_WEIGHT_BOLD : PANGO_WEIGHT_NORMAL);
    pango_font_description_set_size   (statusBarFontDesc, statusBarFont.size*PANGO_SCALE);
    pango_layout_set_font_description (statusBarLayout, statusBarFontDesc ); // ref ?
    pango_layout_set_width            (statusBarLayout, -1);

    redButtonIcon = (GdkPixbuf*)ProcureIconResource(4000);
    yellowButtonIcon = (GdkPixbuf*)ProcureIconResource(4001);
    greenButtonIcon = (GdkPixbuf*)ProcureIconResource(4002);
    orangeButtonIcon = (GdkPixbuf*)ProcureIconResource(4003);

}
//_________________________________________________________________________

void    ShowObjectInspector (void)
{
    long f = FindWindowByName (objectInspectorTitle);
    if (f>=0) {
        gtk_window_present(GTK_WINDOW(windowPtrs(f)));
    } else {
        _HYObjectInspector* newOI = new _HYObjectInspector ();
        newOI->BringToFront       ( );
    }
}

//_________________________________________________________________________

void displayAbout (bool splash)
{
    // TBI
}

//__________________________________________________________________

int main( int   argc, char *argv[] )
{

#ifdef  __HYPHYMPI__
    int            rank,
                   size;

    MPI_Init       (&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    _hy_mpi_node_rank = rank;

    setParameter  (mpiNodeID, (_Parameter)rank);
    setParameter    (mpiNodeCount, (_Parameter)size);

    if (rank == 0)
#endif

        gtk_init (&argc, &argv);

    /* set up globals */

    char curWd[4096];
    getcwd (curWd,4096);

    _String baseDir (curWd);
    baseDir=baseDir & '/';

#ifdef _HYPHY_LIBDIRECTORY_
    _String libDir (_HYPHY_LIBDIRECTORY_);
    pathNames&& &libDir;
#else
    _String libDir (baseDir);
#endif

    if (libDir.sData[libDir.sLength-1] != '/') {
        libDir = libDir & "/";
    }

    pathNames&& &baseDir;
    baseDirectory = baseDir;
    libDirectory = libDir;
    for (long i=1; i<argc; i++) {
        _String thisArg (argv[i]);
        if (thisArg.beginswith ("BASEPATH=")) {
            baseDirectory = thisArg.Cut(9,-1);
            if (baseDirectory.sLength) {
                if (baseDirectory.sData[baseDirectory.sLength-1]!='/') {
                    baseDirectory = baseDirectory & "/";
                }
            } else {
                baseDirectory = baseDir;
            }
        } else if (thisArg.beginswith ("LIBPATH=")) {
            libDirectory = thisArg.Cut(8,-1);
            if (libDirectory.sLength) {
                if (libDirectory.sData[baseDirectory.sLength-1] != '/') {
                    libDirectory = libDirectory & "/";
                }
            } else {
                libDirectory = libDir;
            }
        } else if (thisArg.beginswith ("USEPATH=")) {
            _String     baseArgDir          (thisArg,8,-1);
            errorFileName                   = baseArgDir & errorFileName;
            messageFileName                 = baseArgDir & messageFileName;
            pathNames.Delete                (0);
            pathNames&&                     &baseDir;
        } else if (thisArg.beginswith ("CPU=")) {
#ifdef __MP__
            _String cpus = thisArg.Cut(4,-1);
            systemCPUCount = cpus.toNum();
            if (systemCPUCount<1) {
                systemCPUCount = 1;
            }
#ifdef __MP2__
            pthread_setconcurrency (systemCPUCount+1);
#endif
#endif
        }
#ifdef  __HYPHYMPI__
        else if (thisArg == _String("MPIOPTIMIZER")) {
            mpiParallelOptimizer = true;
            setParameter    (mpiNodeCount, 0.0);
        } else if (thisArg == _String("MPIPARTITIONS")) {
            mpiPartitionOptimizer = true;
            setParameter    (mpiNodeCount, 0.0);
        }
#endif
    }

#ifdef  __HYPHYMPI__
    if (rank == 0)
#endif
    {
        baseDir = libDirectory & "GTKResources";
        _List scanRes;
        ScanDirectoryForFileNames(baseDir,scanRes,false);
        if (scanRes.lLength == 0) {
            GtkWidget * noRez = gtk_message_dialog_new (NULL, GTK_DIALOG_MODAL, GTK_MESSAGE_ERROR, GTK_BUTTONS_OK, "HYPHY_GTK was unable to find a required GTKResources directory in %s. Please use BASEPATH= command line option to specify where the installation directory of HyPhy can be found.", baseDirectory.sData);
            gtk_dialog_run (GTK_DIALOG (noRez));
            gtk_widget_destroy (noRez);
            return 1;
        }
        _String rcPath = baseDir & "/theme/theme.rc";
        //printf ("Loading res files from %s\n", rcPath.sData);
        gtk_rc_parse (rcPath.sData);
    }

    GlobalStartup();

#ifdef  __HYPHYMPI__
    if (rank == 0) {
#endif
        GdkDisplay * defDisplay = gdk_screen_get_display (gdk_screen_get_default());
        hSizeCursor = gdk_cursor_new_for_display (defDisplay,GDK_SB_H_DOUBLE_ARROW);
        pickUpCursor = gdk_cursor_new_for_display (defDisplay,GDK_TARGET);
        dropOffCursor = gdk_cursor_new_for_display (defDisplay,GDK_TCROSS);

        screenPContext = gdk_pango_context_get_for_screen (gdk_screen_get_default());
        tablePDMenuIcon = (GdkPixbuf*)ProcureIconResource(4020);

        /*{
            GdkScreen * defD = gdk_screen_get_default();
            fontConversionFactor = 72.27 / (gdk_screen_get_height (defD) *25.4 / gdk_screen_get_height_mm(defD));
            printf ("Pango conversion factor computed at: %g\n", fontConversionFactor);
        }*/


        ReadInTemplateFiles ();

        hyphyConsoleWindow = new _HYConsoleWindow ("HYPHY Console");
        ReadPreferences     ();
        SetStatusLine ("None","Idle","00:00:00");
        while (gtk_events_pending()) {
            gtk_main_iteration();
        }

        SetPreferences      ();
        ReadGeneticCodes    ();
        ReadModelTemplates  ();
        ReadTreeProcessors ();
        MoveConsoleWindow  (consolePositionRectangle);
        StringToConsole (hyphyCiteString);
        hyphyConsoleWindow->BringToFront();

#ifdef __HYPHYMPI__
        {
            char statBuffer[1024];
            snprintf (statBuffer, sizeof(statBuffer),"MPI version of HyPhy running on %d nodes (a master and %d compute nodes) in %s mode\n",
                     size,
                     size-1,
                     mpiPartitionOptimizer?"partition":(mpiParallelOptimizer?"rate heterogeneity":"normal"));
            BufferToConsole (statBuffer);
        }
#endif

        g_timeout_add  (100,(GSourceFunc)GlobalQueueTimer,nil);
        g_timeout_add  (1000,(GSourceFunc)progressTimerFunction,nil);
        gtk_main ();

        WritePreferences();
#ifdef  __HYPHYMPI__
    } else { // slave node
        if (mpiParallelOptimizer || mpiPartitionOptimizer) {
            mpiOptimizerLoop (rank, size);
        } else {
            mpiNormalLoop (rank, size, baseDir);
        }
    }
#endif

    GlobalShutdown();
    return 0;
}


