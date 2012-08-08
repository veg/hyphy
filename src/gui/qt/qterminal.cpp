//TODO: Add better command line editing

#include <QtGui>
#include <QFont>
#include "qterminal.h"
#include "hy_strings.h"
#include "batchlan.h"

QTerminal::QTerminal(QWidget *parent, Qt::WindowFlags f) : QTextBrowser(parent) {
    setWindowFlags(f);
    cmdStr = "";
    curCursorLoc = textCursor();
    inputCharCount = 0;
    histLocation = -1;
    tempCmd = "";

    //Font Settings
    //QFont font("Monaco Helvetica");
    QFont font("Monaco Helvetica");
    //setTextColor(Qt::darkCyan);
    setCurrentFont(font);

    //Make links clickable
    setTextInteractionFlags(Qt::TextEditorInteraction|Qt::LinksAccessibleByMouse);

    //Allow links to open browser
    setOpenExternalLinks(true);
}

QTerminal::~QTerminal() {

}

void QTerminal::changeDir(const QString & dir) {
    /* Change directories */
    QString theDir = QString(dir);
    theDir = theDir.replace(QChar('/'), "\\");
}

void QTerminal::handleUserLineEntry (void) {
    cmdHistory.push_back(cmdStr);
    histLocation = -1;
    emit userEnteredString(cmdStr);
    cmdStr = "";
    tempCmd = "";

}

void QTerminal::keyPressEvent(QKeyEvent * event) {
    int key = event->key();

    if(!event->modifiers()) {
        setTextCursor(curCursorLoc);
    }

    if (key != Qt::Key_Backspace) {

        if (key == Qt::Key_Return || key == Qt::Key_Enter) {
            inputCharCount = 0;
        } 

        //Command History
        else if (key == Qt::Key_Up) {
            if (cmdHistory.size()) {
                if (histLocation == -1) {
                    histLocation = cmdHistory.size() - 1;
                    tempCmd = cmdStr;
                } 

                else if (histLocation == 0) {
                    QApplication::beep();
                    event->ignore();
                    return;
                } 

                else {
                    --histLocation;
                }

                for (int i = 0; i < inputCharCount; ++i) {
                    QTextEdit::keyPressEvent(new QKeyEvent(QEvent::KeyPress, Qt::Key_Backspace, Qt::NoModifier));
                }

                inputCharCount = cmdHistory.at(histLocation).length();
                insertPlainText(cmdHistory.at(histLocation));
                cmdStr = cmdHistory.at(histLocation);
            }

            event->ignore();
            return;
        } 

        else if (key == Qt::Key_Down) {
            QString str = "";

            if (histLocation == -1) {
                QApplication::beep();
                event->ignore();
                return;
            } 

            else if (histLocation == cmdHistory.size() - 1) {
                histLocation = -1;
                str = tempCmd;
            } 

            else {
                ++histLocation;
                str = cmdHistory.at(histLocation);
            }

            for (int i = 0; i < inputCharCount; ++i) {
                QTextEdit::keyPressEvent(new QKeyEvent(QEvent::KeyPress, Qt::Key_Backspace, Qt::NoModifier));
            }

            inputCharCount = str.length();
            insertPlainText(str);
            cmdStr = str;
        } 

        else if (key == Qt::Key_Left) {
            if (inputCharCount) {
                --inputCharCount;
                QTextEdit::keyPressEvent(event);
            } 

            else {
                QApplication::beep();
            }
        } 

        else if (key == Qt::Key_Right) {
            QTextCursor cursor = textCursor();

            if (cursor.movePosition(QTextCursor::Right)) {
                ++inputCharCount;
                setTextCursor(cursor);
            } 

            else {
                QApplication::beep();
            }
        } 

        else if (key == Qt::Key_Tab) {
            // TODO: Tab Completion
        } 

        else {
            QString text = event->text();
            for (int i = 0; i < text.length(); ++i) {
                if (text.at(i).isPrint()) {
                    //only increase input counter for printable characters
                    ++inputCharCount;
                }
            }

            QTextEdit::keyPressEvent(event);
        }
    } 

    else {
        if (inputCharCount) {
            --inputCharCount;
            QTextEdit::keyPressEvent(event);
            cmdStr.remove(inputCharCount, 1);
        } 
        else {
            QApplication::beep();
        }
    }

    // now pass a char* copy of the input to the shell process
    if (key == Qt::Key_Return || key == Qt::Key_Enter) {
        moveCursor(QTextCursor::End);
        QTextEdit::keyPressEvent(event);
        handleUserLineEntry();
    } 

    else {
        QString input = event->text();
        for (int i = 0; i < input.length(); ++i) {
            if (input.at(i).isPrint()) {
                cmdStr.insert(inputCharCount - 1, input.at(i));
            }
        }
    }

    if(!event->modifiers()) {
        moveCursor(QTextCursor::End);
    }
}

void QTerminal::newline(bool onlyIfLineNonEmpty) {
    moveCursor(QTextCursor::End);
    if (onlyIfLineNonEmpty) {
        if (textCursor().columnNumber() == 0) {
            return;
        }
    }
    insertPlainText ("\n");
}

void QTerminal::prompt(bool hbl) {
    moveCursor(QTextCursor::End);
    if (hbl) {
        _String tag = currentExecutionList?currentExecutionList->GetFileName():empty;
        if (tag.sLength == 0)
                tag = "Analysis waiting for input";
        insertHtml("<font color='#A60000'>["+QString(tag.sData)+"]&gt;</font> ");
    } else {
        insertHtml("<font color='#A60000'>&gt;</font> ");
    }
    moveCursor(QTextCursor::End);

    //Clear undo stack trick
    /*QString text = toHtml();
    setHtml(text);
    moveCursor(QTextCursor::End);*/
}

void QTerminal::insertFromMimeData(const QMimeData * source) {
    //setTextCursor(curCursorLoc);
    moveCursor(QTextCursor::End);

    QString     pastedText = source->text();
    bool        endsWithNL = pastedText.endsWith ('\n');
    QStringList commands = pastedText.split('\n',QString::SkipEmptyParts);

    for (long item = 0; item < commands.count(); item ++) {
        QString str = commands.at (item);
        insertPlainText(str);
        cmdStr += str;
        if (endsWithNL || item < commands.count()-1) {
            newline (true);
            handleUserLineEntry();
            cmdStr = "";
        } else {
            inputCharCount = str.length();
        }
    }

     moveCursor(QTextCursor::End);
}

void QTerminal::contextMenuEvent(QContextMenuEvent *event) {
    QMenu *menu = new QMenu();
    menu->addAction(tr("Undo"),this,SLOT(undo()),QKeySequence::Undo);
    menu->addAction(tr("Redo"),this,SLOT(redo()),QKeySequence::Redo);
    menu->addAction(tr("Copy"),this,SLOT(copy()),QKeySequence::Copy);
    menu->addAction(tr("Paste"),this,SLOT(paste()),QKeySequence::Paste);
    menu->addAction(tr("Find"),this,SLOT(find()),QKeySequence::Find);
    menu->addAction(tr("Select All"),this,SLOT(selectAll()),QKeySequence::SelectAll);
    //menu->addAction(tr("Clear Window"),this,SLOT(clearwindow),QKeySequence::SelectAll);
    menu->exec(event->globalPos());
    delete menu;
}

//Edit Menu Options
//void QTerminal::find(){}
//void QTerminal::clearwindow(){}
