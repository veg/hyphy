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

//For OSX number of cpus
#include <sys/sysctl.h>
#include "hyphy_main.h"
#include "hy_strings.h"
#include "qterminal.h"
#include "hyphyevents.h"

HyphyMain::HyphyMain(QMainWindow *parent) : QMainWindow(parent) {
    setupUi(this);
    this->initialText();
    textEdit->installEventFilter(this);
}

void HyphyMain::initialText() {
    //Options
    textEdit->setLineWrapMode(QTextEdit::FixedColumnWidth);
    textEdit->setWordWrapMode(QTextOption::WordWrap);
    textEdit->setLineWrapColumnOrWidth(80);

    //HyPhy version

    //Number of cpus

    //SW20120702:MPProcessors is deprecated as of OSX10.7, using sysctl
    //hw.physicalcpu
    size_t size;
    int systemCPUCount;
    size = sizeof systemCPUCount;

    sysctlbyname("hw.physicalcpu", &systemCPUCount, &size, NULL, 0);

    if (systemCPUCount == 1) {
        textEdit->insertHtml("One processor detected.\n");
    } else {
        textEdit->insertPlainText(QString::number(systemCPUCount) + " processors detected.\n\n");
    }

    //Model Templates

    //The HyPhy Citation request
    const char* qtHyphyCiteString = "<p>If you use HyPhy in a publication, please cite:<br />S.L. Kosakovsky Pond, S. D. W. Frost"
                                "and S.V. Muse. (2005) HyPhy: hypothesis testing using phylogenies. Bioinformatics 21: 676-679</p>"
                                "<p>If you are a new HyPhy user:"
                                "<br />The tutorial located at <a href='http://www.hyphy.org/docs/HyphyDocs.pdf'>http://www.hyphy.org/docs/HyphyDocs.pdf</a> may be a good starting point.</p><br />";

    textEdit->insertHtml((QString)qtHyphyCiteString);

    //Begin prompting for user input
    textEdit->prompt();
}


bool HyphyMain::eventFilter(QObject *obj, QEvent *event)
{
    
    if (event->type() == BufferToStringType) 
    {
        QBufferToConsoleEvent* e = (QBufferToConsoleEvent*)event;
        textEdit->insertPlainText(e->buffer());
        return true;
    }

    else 
    {
        return false;
    }

}
