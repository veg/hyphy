//TODO: Add better command line editing

#include <QtGui>
#include <QFont>
#include "qterminal.h"
#include "hy_strings.h"
#include "parser.h"

QTerminal::QTerminal(QWidget *parent, Qt::WindowFlags f) : QTextEdit(parent) {
    setWindowFlags(f);
    cmdStr = "";
    curCursorLoc = this->textCursor();
    inputCharCount = 0;
    histLocation = -1;
    tempCmd = "";
    //setTextColor(Qt::darkCyan);

    //QFont font("Monaco Helvetica");
    QFont font("Monaco Helvetica");
    this->setCurrentFont(font);
    //this->setFontFamily("Monaco Helvetica");
}

QTerminal::~QTerminal() {

}

void QTerminal::readBufferOut() {
    this->moveCursor(QTextCursor::End, QTextCursor::KeepAnchor);
}

void QTerminal::changeDir(const QString & dir) {
    /* Change directories */
    QString theDir = QString(dir);
    theDir = theDir.replace(QChar('/'), "\\");
}

void QTerminal::keyPressEvent(QKeyEvent * event) {

    int key = event->key();

    this->setTextCursor(curCursorLoc);

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
                this->insertPlainText(cmdHistory.at(histLocation));
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
            this->insertPlainText(str);
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
            QTextCursor cursor = this->textCursor();

            if (cursor.movePosition(QTextCursor::Right)) {
                ++inputCharCount;
                this->setTextCursor(cursor);
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
        this->moveCursor(QTextCursor::End);
        this->insertPlainText("\n");

        //Execute command string
        ExpressionCalculator((_String)(char *)cmdStr.toAscii().data());

        QTextEdit::keyPressEvent(event);
        cmdHistory.push_back(cmdStr);
        histLocation = -1;
        cmdStr = "";
        tempCmd = "";
        this->prompt();
    } 

    else {
        QString input = event->text();
        for (int i = 0; i < input.length(); ++i) {
            if (input.at(i).isPrint()) {
                cmdStr.insert(inputCharCount - 1, input.at(i));
            }
        }
    }

    curCursorLoc = this->textCursor();

}

void QTerminal::insertFromMimeData(const QMimeData * source) {
    //this->setTextCursor(curCursorLoc);
    this->moveCursor(QTextCursor::End);

    QString pastedText = source->text();
    QStringList commands = pastedText.split("\n");

    foreach (const QString &str, commands) {
        this->insertPlainText(str);
        this->insertPlainText("\n");
        //ExpressionCalculator((_String)(char *)str.toAscii().data());
        cmdHistory.push_back(str);
        histLocation = -1;
        cmdStr = "";
        tempCmd = "";
        this->insertPlainText("\n");
        this->prompt();
     }

    curCursorLoc = this->textCursor();
}

/*
 *void QTerminal::undo() {
 *    qDebug() << "undo was called";
 *}
 */

void QTerminal::prompt() {
    this->insertHtml("<font color=\"#A60000\">></font> ");
}
