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

#include "std_analysis.h"
#include <QList>
#include <QFile>
#include <QTextStream>
#include <QMessageBox>
#include <QRegExp>
#include <QDebug>

//HyPhy specific
#include "batchlan.h"

SelectAnalysisDialog::SelectAnalysisDialog(QWidget *parent) : QDialog(parent) {
    setupUi(this);

    //Tree Widget setup
    treeWidget->setHeaderLabels(QStringList("Standard Analyses"));
    this->loadTree();

    //Connect buttons
    connect(buttonBox, SIGNAL(accepted()), this, SLOT(ok()));
    connect(buttonBox, SIGNAL(rejected()), this, SLOT(cancel()));
}

void SelectAnalysisDialog::loadTree() {
    //Regex for parsing files.lst
    //"","Various tools for estimating the distribution of substitution rates.","!Substitution Rates";
    QRegExp rx("\"(.*)\",\"(.*)\",\"(.*)\";");

    //"SR","DNARates like site specific rate estimation. Other model parameters are inferred from the entire alignment.","SiteRates.bf";

    //Read files.lst
    QFile file("res/TemplateBatchFiles/files.lst");
    if(!file.open(QIODevice::ReadOnly)) {
        //QMessageBox::information(0, "error", file.errorString());
        QMessageBox::information(0, "error", "error");
    }

    QTextStream in(&file);

    //Before looping, let's set the first top level item
    QString line = in.readLine();    

    //Have to implement RegEx
    rx.indexIn(line);
    QStringList fields = rx.capturedTexts();

    //treeWidget->setColumnCount(1);
    QList<QTreeWidgetItem *> items;

    //First item must always be a top level item
    if(!fields[3].startsWith("!")) {
        qDebug() << fields;
        QMessageBox::information(0, "error", "Invalid File, first line is not a top level item");
    }

    //Set the first item. The remove is for the beginning "!"
    QTreeWidgetItem* item = new QTreeWidgetItem((QTreeWidget*)0, QStringList(fields[3].remove(0,1)));
    //item.setWhatsThis(fields[0]);

    while(!in.atEnd()) {

        //Before looping, let's set the first top level item
        line = in.readLine();    
        rx.indexIn(line);
        fields = rx.capturedTexts();

        //Need to do a while here, and then check if the label starts with a bang
        if(fields[3].startsWith("!")) {
            items.append(item);
            item = new QTreeWidgetItem((QTreeWidget*)0, QStringList(fields[3].remove(0,1)));
        }

        //not a bang, add to the current item until we find the next one
        else {
            item->addChild(new QTreeWidgetItem((QTreeWidget*)0,QStringList(fields[3])));
        }
    }

    //Append last item
    items.append(item);

    treeWidget->insertTopLevelItems(0, items);
    file.close();
    return;
}

void SelectAnalysisDialog::ok() {
    //Get currently selected treewidget item
    QTreeWidgetItem* selected =  treeWidget->currentItem();
    _ExecutionList ex;

    //Get filename
    QString batchfile_selected = selected->text(0); 
    _String whole_path = _String("res/TemplateBatchFiles/")&(_String)batchfile_selected.toLocal8Bit().data();
    ReadBatchFile(whole_path,ex);
    this->accept();
}

void SelectAnalysisDialog::cancel() {
    this->reject();
}
