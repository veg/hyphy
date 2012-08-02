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
#include "avllistxl.h"
#include "hyphyhierarchicalselector.h"

_List      _HY_HierarchicalSelectorPreviousSelections_Aux;
_AVLListXL _HY_HierarchicalSelectorPreviousSelections(&_HY_HierarchicalSelectorPreviousSelections_Aux);

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
    connect(itemList, SIGNAL(itemDoubleClicked ( QTreeWidgetItem*, int ) ), this, SLOT(handle_double_click()));
    SetInitialSelection ();
}

//______________________________________________________________


void _HY_HierarchicalSelector::ok() {
    //Get currently selected treewidget item
    QList <QTreeWidgetItem*> selection = itemList->selectedItems();
    if (selection.count()) {
        _String * selectedValue = new _String(selection.at(0)->text(0).toAscii().data()),
                wTitle (windowTitle().toAscii().data());
        
         _HY_HierarchicalSelectorPreviousSelections.UpdateValue(&wTitle, selectedValue, false, true);
    }
    this->accept();
}

//______________________________________________________________

void _HY_HierarchicalSelector::SetInitialSelection () {
    _String wTitle (windowTitle().toAscii().data()),
            *value = (_String*)_HY_HierarchicalSelectorPreviousSelections.GetDataByKey(&wTitle);
    QTreeWidgetItem * theItem;        
    if (value) {
        QString theValue (value->sData);
        QList<QTreeWidgetItem *> matches = itemList->findItems(theValue, Qt::MatchStartsWith);
        if (matches.count() > 0) {
            theItem = matches.at(0);
        }
    } else {
        theItem = itemList->topLevelItem(0);
    }
    itemList->setCurrentItem(theItem);
}


//______________________________________________________________

void _HY_HierarchicalSelector::cancel() {
    selections->Clear();
    this->reject();
}

//______________________________________________________________

void _HY_HierarchicalSelector::toggleAcceptStatus (long selected) {
    validSelection = (result>0 && selected == result) || (selected > 0);
    QPushButton* OK = buttonBox->button(QDialogButtonBox::Ok);
    if (OK) {
        OK->setEnabled(validSelection);
    }   
}

//______________________________________________________________

void _HY_HierarchicalSelector::handle_double_click() {
    if (validSelection) {
        accept();
    }
}

//______________________________________________________________

void _HY_HierarchicalSelector::handle_selection_change() {
    QList <QTreeWidgetItem*> selection = itemList->selectedItems();
    
    _String aboutText (128L, true);
    long    rawSelected = selection.count();
            
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
        }
        aboutText << *(_String*)dData(item_offset);
        if (k) {
            aboutText << '\n';
        }
    }
    aboutText.Finalize();    
    
    explanationLabel->document()->setPlainText (QString(aboutText.sData));
    toggleAcceptStatus (selections->lLength);
}

//______________________________________________________________

void _HY_HierarchicalSelector::closeEvent(QCloseEvent *event) {
    if (!validSelection) {
        selections->Clear();
    }
    QDialog::closeEvent (event);
}

//______________________________________________________________
 
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

