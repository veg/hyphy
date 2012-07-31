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

#include "hyphyhierarchicalselector.h"
#include <QList>
#include <QFile>
#include <QTextStream>
#include <QMessageBox>
#include <QRegExp>
#include <QDebug>
#include <QPushButton>
#include <QKeyEvent>

//HyPhy specific
#include "batchlan.h"

_HY_HierarchicalSelector::_HY_HierarchicalSelector(QWidget *parent, _List& definition, _SimpleList& c, _SimpleList& vc, _String n, _SimpleList* s, long r, bool is_modal) : QDialog(parent) {
    setupUi(this);

    //Tree Widget setup
    setWindowTitle(QString(n.sData));
    
    if (r >= 1) {
        n = n & " (select " & r &")";
    } else {
        n = n & " (select 1 or more)";
    }
    
    
    
    itemList->setHeaderLabels(QStringList(n.sData));
    if (r != 1) {
        itemList->setSelectionMode (QAbstractItemView::ExtendedSelection);
    }
    //this->loadTree();
    if (is_modal) {
        if (parent) {
            setWindowFlags (Qt::Sheet);
        } else {
            this->setModal (true);
        }
    } else {
        this->setModal (false);
    }
    
    validChoices.Duplicate(&vc);

    unsigned long  counter;
    
    QTreeWidgetItem       *  current_item = NULL;
    
    for (counter = 0; counter< validChoices.lLength; counter++) {
        _String* tStr = (_String*)(*(_List*)(definition(validChoices.lData[counter])))(c.lData[0]);
        bool is_a_header = tStr->sData[0] == '!';
        if (is_a_header) {
            iData.AppendNewInstance (new _String(*tStr,1,-1));
            offsets << counter;
        } else {
            iData  << tStr;
        }
        dData << (*(_List*)(definition(validChoices.lData[counter])))(c.lData[1]);
        QTreeWidgetItem * this_item = new QTreeWidgetItem((QTreeWidget*)NULL, QStringList(((_String*)iData(iData.lLength-1))->sData));
        if (is_a_header) {
            current_item = this_item;
            //this_item->setFlags (Qt::ItemIsEnabled);
        }  else {
            if (current_item) {
                current_item -> addChild (this_item);
                continue;
            }
        }
        if (!is_a_header)
            offsets << counter;
        itemList->addTopLevelItem(this_item);
    }
    
    selections = s;
    result     = r;
    validSelection = false;
    toggleAcceptStatus (0);
    
    buttonBox->installEventFilter(this);
    itemList->installEventFilter(this);
    explanationLabel->installEventFilter (this);
    


   //Connect buttons
    connect(buttonBox, SIGNAL(accepted()), this, SLOT(ok()));
    connect(buttonBox, SIGNAL(rejected()), this, SLOT(cancel()));
    connect(itemList, SIGNAL(itemSelectionChanged()), this, SLOT(handle_selection_change()));
}

/*
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

        line = in.readLine();    
        rx.indexIn(line);
        fields = rx.capturedTexts();

        //Check if the label starts with a bang
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
*/

void _HY_HierarchicalSelector::ok() {
    //Get currently selected treewidget item
    this->accept();
}


void _HY_HierarchicalSelector::cancel() {
    selections->Clear();
    this->reject();
}

void _HY_HierarchicalSelector::toggleAcceptStatus (long selected) {
    validSelection = (result>0 && selected == result) || (selected > 0);
    QPushButton* OK = buttonBox->button(QDialogButtonBox::Ok);
    if (OK) {
        OK->setEnabled(validSelection);
    }   
}

void _HY_HierarchicalSelector::handle_selection_change() {
    QList <QTreeWidgetItem*> selection = itemList->selectedItems();
    
    _String aboutText (128L, true);
    long    totalSelected = 0,
            rawSelected = selection.count();
            
    selections->Clear();
    
    for (long k = 0; k < rawSelected; k ++ )    {
        QTreeWidgetItem* theSelection = selection.at(k),
                       * parent = theSelection->parent();
        long           item_offset = 0;
        if (parent) {
            item_offset = offsets.lData[itemList->indexOfTopLevelItem(parent)] + parent->indexOfChild(theSelection)+1;
        } else {
            item_offset = offsets.lData[itemList->indexOfTopLevelItem(theSelection)];   
        }
        if (theSelection->childCount() == 0) {
            (*selections) << item_offset;
            totalSelected++;
        }
        aboutText << *(_String*)dData(item_offset);
        if (k) {
            aboutText << '\n';
        }
    }
    aboutText.Finalize();    
    
    explanationLabel->document()->setPlainText (QString(aboutText.sData));
    toggleAcceptStatus (totalSelected);
}

void _HY_HierarchicalSelector::closeEvent(QCloseEvent *event) {
    if (!validSelection) {
        selections->Clear();
    }
    QDialog::closeEvent (event);
}

bool _HY_HierarchicalSelector::eventFilter ( QObject * watched, QEvent * event ) {
    if (event->type () == QEvent::KeyPress) {
        QKeyEvent *keyEvent = (QKeyEvent*)(event);
        if (keyEvent->key () == Qt::Key_Escape) {
            selections->Clear();
            return false;
        }
     }
    return QDialog::eventFilter(watched, event);
}

