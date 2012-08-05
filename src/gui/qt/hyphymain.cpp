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

#include <QtGui>
#include <QDebug>
#include "hyphymain.h"

//For OSX number of cpus
#include <sys/sysctl.h>

#include "hy_strings.h"
#include "HYUtils.h"

#include "qterminal.h"
#include "hyphyhierarchicalselector.h"
#include "hyphyevents.h"
#include "hyphy_qt_helpers.h"

#include "HYSharedMain.h"

HyphyMain::HyphyMain(QMainWindow *parent) : QMainWindow(parent) {
    setupUi(this);
    waitingOnStringFromConsole = false;
    initialText();
    initializeMenuBar();
    initializeStatusBar();
    installEventFilter(this);

    this->_statusBarVisible = true;
}

void HyphyMain::initialText() {
    //Options
    console->setLineWrapMode(QTextEdit::FixedColumnWidth);
    console->setWordWrapMode(QTextOption::WordWrap);
    console->setLineWrapColumnOrWidth(80);
    console->setLineWrapMode(QTextEdit::WidgetWidth);

    //HyPhy version
    _String version = GetVersionString();
    console->insertHtml((QString)version + ". ");
 
    // SW20120702
    // MPProcessors is deprecated as of OSX10.7, using sysctl
    size_t size;
#ifdef _OPENMP
    systemCPUCount = omp_get_max_threads();
#else
    systemCPUCount = 1;
#endif
    
    console->insertHtml("Up to " + QString::number(systemCPUCount) + " threads can be used for analyses.<p>");

    //The HyPhy Citation request
    const char* qtHyphyCiteString = "<p>If you use HyPhy in a publication, please cite <b>SL Kosakovsky Pond, SDW Frost and SV Muse. (2005)</b> HyPhy: hypothesis testing using phylogenies. <i>Bioinformatics 21: 676-679</i></p><p>If you are a new HyPhy user, the tutorial located at <a href='http://www.hyphy.org/docs/HyphyDocs.pdf'>http://www.hyphy.org/docs/HyphyDocs.pdf</a> may be a good starting point.</p><br>";

    console->insertHtml((QString)qtHyphyCiteString);

    //Begin prompting for user input
    console->prompt();
}

void HyphyMain::initializeMenuBar() {

    //Create File Menu Options
    _hyConsoleOpenAction = new QAction(tr("&Open Batch File"), this);
    _hyConsoleOpenAction->setShortcuts(QKeySequence::Open);
    _hyConsoleOpenAction->setStatusTip(tr("Open a HyPhy Batch Language File"));
    _hyConsoleSaveAction = new QAction(tr("&Save Console"), this);
    _hyConsoleSaveAction->setShortcuts(QKeySequence::Save);
    _hyConsoleSaveAction->setStatusTip(tr("Save HyPhy console content"));
    _hyConsoleExitAction = new QAction(tr("&Quit"), this);

    //Connect File Menu Events to appropriate slots
    connect(_hyConsoleOpenAction, SIGNAL(triggered()), this, SLOT(hy_open()));
    connect(_hyConsoleSaveAction, SIGNAL(triggered()), this, SLOT(hy_save()));
    connect(_hyConsoleExitAction, SIGNAL(triggered()), qApp, SLOT(quit()));

    //Add the File Menu to the Menu Bar
    _hyConsoleMenu = menuBar()->addMenu(tr("&File"));
    _hyConsoleMenu->addAction(_hyConsoleOpenAction);
    _hyConsoleMenu->addAction(_hyConsoleSaveAction);
    _hyConsoleMenu->addSeparator();
    _hyConsoleMenu->addAction(_hyConsoleExitAction);

    //Create Edit Menu Options
    _hyConsoleUndoAction = new QAction(tr("&Undo"),this);
    _hyConsoleUndoAction->setShortcuts(QKeySequence::Undo);
    _hyConsoleRedoAction = new QAction(tr("&Redo"),this);
    _hyConsoleRedoAction->setShortcuts(QKeySequence::Redo);

    //No Cutting allowed for right now
    _hyConsoleCopyAction = new QAction(tr("&Copy"),this);
    _hyConsoleCopyAction->setShortcuts(QKeySequence::Copy);
    _hyConsolePasteAction = new QAction(tr("&Paste"),this);
    _hyConsolePasteAction->setShortcuts(QKeySequence::Paste);
    _hyConsoleFindAction = new QAction(tr("&Find"),this);
    _hyConsoleFindAction->setShortcuts(QKeySequence::Find);
    _hyConsoleSelectAllAction = new QAction(tr("&Select All"),this);
    _hyConsoleSelectAllAction->setShortcuts(QKeySequence::SelectAll);
    _hyConsoleClearWindowAction = new QAction(tr("&Clear Window"),this);
    _hyConsoleClearWindowAction->setStatusTip("Clears the Console Window");

    //Connect Edit Menu Events to appropriate slots
    connect(_hyConsoleUndoAction, SIGNAL(triggered()), console, SLOT(undo()));
    connect(_hyConsoleRedoAction, SIGNAL(triggered()), console, SLOT(redo()));
    //connect(_hyConsoleCutAction, SIGNAL(triggered()), console, SLOT(cut()));
    connect(_hyConsoleCopyAction, SIGNAL(triggered()), console, SLOT(copy()));
    connect(_hyConsolePasteAction, SIGNAL(triggered()), console, SLOT(paste()));
    //connect(_hyConsoleFindAction, SIGNAL(triggered()), console, SLOT(find()));
    connect(_hyConsoleSelectAllAction, SIGNAL(triggered()), console, SLOT(selectAll()));
    //connect(_hyConsoleClearWindowAction, SIGNAL(triggered()), console, SLOT(clearwindow()));

    connect(console, SIGNAL(userEnteredString(const QString)), this, SLOT(handle_user_input(const QString)));

    //Add the Edit Menu to the Menu Bar
    _hyConsoleMenu = menuBar()->addMenu(tr("&Edit"));
    _hyConsoleMenu->addAction(_hyConsoleUndoAction);
    _hyConsoleMenu->addAction(_hyConsoleRedoAction);
    _hyConsoleMenu->addSeparator();
    //_hyConsoleMenu->addAction(_hyConsoleCutAction);
    _hyConsoleMenu->addAction(_hyConsoleCopyAction);
    _hyConsoleMenu->addAction(_hyConsolePasteAction);
    _hyConsoleMenu->addSeparator();
    _hyConsoleMenu->addAction(_hyConsoleFindAction);
    _hyConsoleMenu->addSeparator();
    _hyConsoleMenu->addAction(_hyConsoleSelectAllAction);
    _hyConsoleMenu->addAction(_hyConsoleClearWindowAction);


    //Create View Menu Options
    _hyConsoleHideStatusBarAction = new QAction(trUtf8("&\U00002713 Status Bar"),this);

    //Connect View Menu Events to appropriate slots
    connect(_hyConsoleHideStatusBarAction, SIGNAL(triggered()), this, SLOT(toggle_status_bar()));

    //Add the View Menu to the Menu Bar
    _hyConsoleMenu = menuBar()->addMenu(tr("&View"));
    _hyConsoleMenu->addAction(_hyConsoleHideStatusBarAction);

    //Create Analysis Menu Options
    _hyConsoleCancelExecutionAction = new QAction(tr("&Cancel Execution"),this);
    _hyConsoleCancelExecutionAction->setStatusTip("Cancels Execution of Analysis Currently Running");
    _hyConsoleSuspendExecutionAction = new QAction(tr("&Suspend Execution"),this);
    _hyConsoleViewLogAction = new QAction(tr("&View Log"),this);
    _hyConsoleStandardAnalysisAction = new QAction(tr("&Standard Analysis"),this);
    QList <QKeySequence>keyList;
    keyList << QKeySequence ("Ctrl+E");
    _hyConsoleStandardAnalysisAction->setShortcuts(keyList);
    _hyConsoleStandardAnalysisAction->setStatusTip("Standard Analysis");
    _hyConsoleResultsAction = new QAction(tr("&Results"),this);
    _hyConsoleRerunLastAnalysisAction = new QAction(tr("&Rerun Last Analysis"),this);

    //Connect Analysis Menu Events to appropriate slots
    connect(_hyConsoleCancelExecutionAction, SIGNAL(triggered()), this, SLOT(hy_cancelexecution()));
    connect(_hyConsoleSuspendExecutionAction, SIGNAL(triggered()), this, SLOT(hy_suspendexecution()));
    connect(_hyConsoleViewLogAction, SIGNAL(triggered()), this, SLOT(hy_viewlog()));
    connect(_hyConsoleStandardAnalysisAction, SIGNAL(triggered()), this, SLOT(hy_standardanalysis()));
    connect(_hyConsoleResultsAction, SIGNAL(triggered()), this, SLOT(hy_results()));
    connect(_hyConsoleRerunLastAnalysisAction, SIGNAL(triggered()), this, SLOT(hy_rerunlastanalysis()));

    //Add the Analysis Menu to the Menu Bar
    _hyConsoleMenu = menuBar()->addMenu(tr("&Analysis"));
    _hyConsoleMenu->addAction(_hyConsoleCancelExecutionAction);
    _hyConsoleMenu->addAction(_hyConsoleSuspendExecutionAction);
    _hyConsoleMenu->addSeparator();
    _hyConsoleMenu->addAction(_hyConsoleViewLogAction);
    _hyConsoleMenu->addSeparator();
    _hyConsoleMenu->addAction(_hyConsoleStandardAnalysisAction);
    _hyConsoleMenu->addAction(_hyConsoleResultsAction);
    _hyConsoleMenu->addAction(_hyConsoleRerunLastAnalysisAction);

    //Create Windows Menu Options
    _hyConsoleMinimizeAction = new QAction(tr("&Minimize"),this);
    _hyConsoleConsoleWindowAction = new QAction(tr("&Console Window"),this);
    _hyConsoleObjectInspectorAction = new QAction(tr("&Object Inspector"),this);
    _hyConsoleCycleThroughWindowsAction = new QAction(tr("&Cycle Through Windows"),this);

    //Connect Windows Menu Events to appropriate slots
    connect(_hyConsoleMinimizeAction, SIGNAL(triggered()), this, SLOT(hy_minimize()));
    connect(_hyConsoleConsoleWindowAction, SIGNAL(triggered()), this, SLOT(hy_consolewindow()));
    connect(_hyConsoleObjectInspectorAction, SIGNAL(triggered()), this, SLOT(hy_objectinspector()));
    connect(_hyConsoleCycleThroughWindowsAction, SIGNAL(triggered()), this, SLOT(hy_cyclethroughwindows()));

    //Add the Windows Menu to the Menu Bar
    _hyConsoleMenu = menuBar()->addMenu(tr("&Window"));
    _hyConsoleMenu->addAction(_hyConsoleMinimizeAction);
    _hyConsoleUndoAction->setShortcuts(QKeySequence::Undo);
    _hyConsoleMenu->addAction(_hyConsoleConsoleWindowAction);
    _hyConsoleMenu->addAction(_hyConsoleObjectInspectorAction);
    _hyConsoleMenu->addSeparator();
    _hyConsoleMenu->addAction(_hyConsoleCycleThroughWindowsAction);
}

void HyphyMain::initializeStatusBar() {
    //Status elements
    //progressBar;
    //statusLine;
    //timerDisplay;

    //Initialize progressBar to 0
    progressBar->setValue(0);

    //Timer stuff
    this->_elapsed_timer = new QElapsedTimer();
    this->_timer = new QTimer(this);
    connect(this->_timer, SIGNAL(timeout()), this, SLOT(update_timer_display()));
    //this->StartBarTimer();
}

//File Menu Options
void HyphyMain::hy_open() {
    if (OpenBatchFile (true)) {
         ExecuteBatchFile();
    }
}

void HyphyMain::hy_save() {
    console->insertPlainText((QString)_hyQTFileDialog ("Save console contents to", "HyPhy Console.txt", true));
}

void HyphyMain::quit() {
}

//View Menu Options
void HyphyMain::toggle_status_bar() {

    //Use of bool just in case a widget gets out of sync
    if(this->_statusBarVisible) {
        this->_hyConsoleHideStatusBarAction->setText("&Status Bar");
        this->progressBar->hide();
        this->statusLine->hide();
        this->timerDisplay->hide();
    }

    //Display Check Mark next to Status Bar menu option if visible
    else {
        this->_hyConsoleHideStatusBarAction->setText(trUtf8("&\U00002713 Status Bar"));
        this->progressBar->show();
        this->statusLine->show();
        this->timerDisplay->show();
    }

    this->_statusBarVisible = !this->_statusBarVisible;
}

//Analysis Menu
void HyphyMain::hy_cancelexecution(){}
void HyphyMain::hy_suspendexecution(){}
void HyphyMain::hy_viewlog(){}

void HyphyMain::hy_standardanalysis() {

    _SimpleList std,
                vc (availableTemplateFiles.lLength,0,1),
                selections;
                
    std<<2;
    std<<1;

    qDebug() << availableTemplateFiles.lLength;
    _HY_HierarchicalSelector *hs = new _HY_HierarchicalSelector(this, availableTemplateFiles, std, vc, "Standard Analyses", &selections, 1, true);
    hs->setWindowModality(Qt::WindowModal);
    hs->exec(); 
    if (selections.lLength == 1) {
        RunTemplate (selections.lData[0]);
    }
}

void HyphyMain::hy_results(){}
void HyphyMain::hy_rerunlastanalysis(){}

//Window Menu
void HyphyMain::hy_minimize(){}
void HyphyMain::hy_consolewindow(){}
void HyphyMain::hy_objectinspector(){}
void HyphyMain::hy_cyclethroughwindows(){}

void HyphyMain::AppendTextToEnd(const QString& text, bool isHTML, _SimpleList* color) {
    console->moveCursor(QTextCursor::End);
    if (isHTML) 
        console->insertHtml(text);
    else {
        QColor currentColor = console->textColor();
        if (color && color->lLength >= 3) {
            console->setTextColor (QColor (color->lData[0], color->lData[1], color->lData[2]));
        } else {
            console->setTextColor (QColor (0,0,212));            
        }
        console->insertPlainText(text);
        console->setTextColor (currentColor);
    }
}

void HyphyMain::DisplayPrompt     (void) {
    console->moveCursor(QTextCursor::End);
    console->prompt();
}

bool HyphyMain::eventFilter(QObject *obj, QEvent *event) {
    if (event->type() == BufferToStringType) {
        AppendTextToEnd (((QBufferToConsoleEvent*)event)->buffer(), false, &((QBufferToConsoleEvent*)event)->color());
        return true;
    }

    return QMainWindow::eventFilter(obj,event);
}

//Status line specific
void HyphyMain::StartBarTimer() {
      this->_timer->start(1000);
      this->_elapsed_timer->start();
}

void HyphyMain:: StopBarTimer() {
      this->_timer->stop();
}

void HyphyMain::update_timer_display() {
    //Need to format the elapsed time
    this->timerDisplay->display((int)(_elapsed_timer->elapsed()/1000));
}

void HyphyMain::SetStatusLine(_String) {

}

void HyphyMain::SetStatusLine     (_String, _String, _String, long l){

}

void HyphyMain::SetStatusLine     (_String, _String, _String){

}

void HyphyMain::SetStatusLine     (_String, _String, _String, long, char){
}

void HyphyMain::SetStatusBarValue (long, _Parameter, _Parameter){
}

void HyphyMain::AddStringToRecentMenu (const _String, const _String) {

}

void HyphyMain::setWaitingOnStringFromConsole (bool value) {
    if (value) {
        console->newline(true);
        console->prompt(true);
    }
    waitingOnStringFromConsole = value;
}

void HyphyMain::handle_user_input(const QString data) {
    if (waitingOnStringFromConsole) {
        setWaitingOnStringFromConsole (false);
        userData = (_String)(char *)data.toAscii().data();
        emit handled_user_input();
    } else {
        ExpressionCalculator((_String)(char *)data.toAscii().data());
        console->newline();
        console->prompt();
    }
}
