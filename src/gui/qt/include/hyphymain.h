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

#include <QtGui>
#include <QProgressBar>
#include <QtNetwork>
#include <QElapsedTimer>

#include "ui_hyphymain.h"
#include "qterminal.h"
#include "hy_strings.h"

#define  HY_CONSOLE_CAN_COPY    0x01
#define  HY_CONSOLE_CAN_PASTE   0x02
#define  HY_CONSOLE_CAN_UNDO    0x04
#define  HY_CONSOLE_CAN_REDO    0x08

#define  HY_SL_FILE             0x01
#define  HY_SL_TASK             0x02
#define  HY_SL_TIMER            0x04
#define  HY_SL_PERCENT          0x08
#define  HY_SL_SUSPEND          0x10
#define  HY_SL_RESUME           0x20
#define  HY_SL_FORCE            0x40
#define  HY_SL_DONE             0x80

class HyphyMain : public QMainWindow, private Ui::MainWindow
{

    Q_OBJECT

public:
    HyphyMain(QMainWindow *parent = NULL);
    void initialText();
    void StartBarTimer();
    void StopBarTimer();
    void SetStatusLine     (_String);
    void SetStatusLine     (_String, _String, _String, long l);
    void SetStatusLine     (_String, _String, _String);
    void SetStatusLine     (_String, _String, _String, long, char);
    void SetStatusBarValue (long, _Parameter, _Parameter);
    void ClearStatusLine     ();
    void AddStringToRecentMenu (const _String, const _String);
    void AppendTextToEnd   (const QString&, bool isHTML = false, _SimpleList* color = nil);
    void DisplayPrompt     (void);
    void setWaitingOnStringFromConsole (bool);
    const _String & getUserData (void) { return userData;}

protected:
    bool eventFilter(QObject *obj, QEvent *ev);
    virtual void closeEvent(QCloseEvent *);

signals:
    void handled_user_input ();
    
private slots:
    //File Menu
    void hy_open();
    void hy_save();
    void quit();
    void handle_user_input (const QString);

    //Analysis Menu
    void hy_cancelexecution();
    void hy_suspendexecution();
    void hy_viewlog();
    void hy_standardanalysis();
    void hy_results();
    void hy_rerunlastanalysis();

    //Window Menu
    void hy_minimize();
    void hy_consolewindow();
    void hy_objectinspector();
    void hy_cyclethroughwindows();

    //Status Timer
    void    update_timer_display();

private:
    QTimer *_timer;
    QElapsedTimer  *_elapsed_timer;

    //File Actions
    QAction *_hyConsoleOpenAction;
    QAction *_hyConsoleSaveAction;
    QAction *_hyConsoleExitAction;

    //Edit Actions
    QAction *_hyConsoleUndoAction;
    QAction *_hyConsoleRedoAction;
    QAction *_hyConsoleCutAction;
    QAction *_hyConsoleCopyAction;
    QAction *_hyConsolePasteAction;
    QAction *_hyConsoleFindAction;
    QAction *_hyConsoleSelectAllAction;
    QAction *_hyConsoleClearWindowAction;

    //View Actions
    QAction *_hyConsoleHideStatusBarAction;

    //Analysis Actions
    QAction *_hyConsoleCancelExecutionAction;
    QAction *_hyConsoleSuspendExecutionAction;
    QAction *_hyConsoleViewLogAction;
    QAction *_hyConsoleStandardAnalysisAction;
    QAction *_hyConsoleResultsAction;
    QAction *_hyConsoleRerunLastAnalysisAction;


    //Window Actions
    QAction *_hyConsoleMinimizeAction;
    QAction *_hyConsoleConsoleWindowAction;
    QAction *_hyConsoleObjectInspectorAction;
    QAction *_hyConsoleCycleThroughWindowsAction;

    QMenu   *_hyConsoleMenu;
    void    initializeMenuBar();
    void    initializeStatusBar();
    QString  lastAnalysisFilePath;
    
    // preferences options
    void     ReadSettings  (void);
    void     WriteSettings (void);
    
    bool     waitingOnStringFromConsole;
    _String  userData;

    //Status Bar
    QStatusBar *status;
    QLabel *filename_status;
    QLabel *updated_status;
    QProgressBar *progress_bar;

    void printHelp();
    };

extern HyphyMain* _hyPrimaryConsoleWindow;
