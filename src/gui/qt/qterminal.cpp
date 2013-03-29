//TODO: Add better command line editing

#include <QtGui>
#include <QFont>
#include "qterminal.h"
#include "hy_strings.h"
#include "batchlan.h"

//#define __QTERMINAL_DEBBUGGING

QTerminal::QTerminal(QWidget *parent, Qt::WindowFlags f) : QTextBrowser(parent) {
    setWindowFlags(f);
    
    current_command                   = "";
    location_of_last_propmpt          = textCursor();
    location_of_insertion_point       = textCursor();
    connect(this, SIGNAL(cursorPositionChanged()), this, SLOT(handleCursorMove()));
    
    length_of_input_buffer            = 0;

    command_history_location           = -1;
 
    //Font Settings
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
    // SLKP 20121205: seems unused
    QString theDir = QString(dir);
    theDir = theDir.replace(QChar('/'), "\\");
}


void QTerminal::handleCursorMove (void) {
    int position = textCursor().position();
    if (location_of_insertion_point.position() != position && position > location_of_last_propmpt.position()) {
        location_of_insertion_point = textCursor();
    }
}

void QTerminal::handleUserLineEntry (void) {
    
    command_history.push_back(current_command);
    command_history_location = -1;
    
#ifdef __QTERMINAL_DEBBUGGING
    printf ("handleUserLineEntry: %s\n", current_command.toAscii().data());
#endif
    
    emit userEnteredString(current_command);
    
    current_command         = "";
    stash_current_command   = "";
    length_of_input_buffer  = 0;

}

void        QTerminal::beepAndIgnoreEvent (QEvent* event, bool do_ignore) {
    QApplication::beep();
    if (do_ignore)
        event->ignore();
}

void        QTerminal::replaceCurrentCommand   (const QString* replace_with, bool stash_current){
    if (stash_current) {
        stash_current_command = current_command;
    }
    current_command = *replace_with;
    for (int i = 0; i < length_of_input_buffer; ++i) {
        QTextEdit::keyPressEvent(new QKeyEvent(QEvent::KeyPress, Qt::Key_Backspace, Qt::NoModifier));
    }
    insertPlainText(*replace_with);
    length_of_input_buffer      = replace_with->length();
    moveCursor(QTextCursor::End);
    location_of_insertion_point = textCursor();
}


void QTerminal::keyPressEvent(QKeyEvent * event) {
    
    int     key        = event->key();
 
    /*if(event->modifiers()) {
        QTextEdit::keyPressEvent(event);
        return;
    }*/
    
    switch (key) {
        case Qt::Key_Up: { // fetch an older command from history
            if (command_history.size()) {
                if (command_history_location == -1) {
                    replaceCurrentCommand (&command_history.at(command_history_location = command_history.size()-1),true);
                    event->ignore();
                    return;
                } else {
                    if (command_history_location > 0) {
                       replaceCurrentCommand (&command_history.at(--command_history_location),false);
                       event->ignore();
                       return;
                    }
                }
            }
            beepAndIgnoreEvent(event);
            return;
        }
        case Qt::Key_Down: { // fetch a newer command from history
            if (command_history_location >= 0) {
                if (command_history_location == command_history.size()-1) {
                    replaceCurrentCommand (&stash_current_command,false);
                    command_history_location = -1;
                } else {
                    replaceCurrentCommand (&command_history.at(++command_history_location),false);
                }
                event->ignore();
            } else {
                beepAndIgnoreEvent(event);
            }
            return;
        }
        case Qt::Key_Left: { // move left in the command string
            setTextCursor(location_of_insertion_point);
            if (location_of_insertion_point.position() > location_of_last_propmpt.position() + 1) {
                location_of_insertion_point.movePosition (QTextCursor::Left);
                setTextCursor (location_of_insertion_point);
            } else {
                beepAndIgnoreEvent(event);
            }
            return;
        }
        case Qt::Key_Right: { // move right in the command string
            setTextCursor(location_of_insertion_point);
            if (location_of_insertion_point.movePosition (QTextCursor::Right)) {
                 setTextCursor (location_of_insertion_point);
            } else {
                beepAndIgnoreEvent(event);
            }
            return;
        }
        /*
        case Qt::Key_Tab: {
            // auto-complete later
            break;
        }
        */
        case Qt::Key_Return:
        case Qt::Key_Enter: {
            moveCursor(QTextCursor::End);
            QTextEdit::keyPressEvent(event);
            location_of_insertion_point = textCursor();
            handleUserLineEntry();
#ifdef __QTERMINAL_DEBBUGGING
            printf ("new line start = %d insert = %d current = %d\n", location_of_last_propmpt.position(), location_of_insertion_point.position(),textCursor().position());
#endif
            return;
        }
            
        case Qt::Key_Backspace: {
            
            setTextCursor(location_of_insertion_point);
        
            long diff = location_of_insertion_point.position() - location_of_last_propmpt.position();
            if (diff>1) {
               current_command.remove(diff-2, 1);
               QTextEdit::keyPressEvent(event);
               location_of_insertion_point = textCursor();
               length_of_input_buffer --;
            }
            else {
                beepAndIgnoreEvent(event);
            }
            return;
        }
            
        default:{
             QString text = event->text();
            int  diff = location_of_insertion_point.position() - location_of_last_propmpt.position(),
                 inserted = 0;
            
            for (int i = 0; i < text.length(); ++i) {
                if (text.at(i).isPrint()) {
                    current_command.insert(diff + inserted++ - 1, text.at(i));
                }
            }
#ifdef __QTERMINAL_DEBBUGGING
            printf ("keyPressEvent enter text: start = %d insert = %d printable %d buffer %s\n", location_of_last_propmpt.position(), location_of_insertion_point.position(), inserted, current_command.toAscii().data());
#endif
            if (inserted) {
                setTextCursor(location_of_insertion_point);
                length_of_input_buffer+=inserted;
            }
            QTextEdit::keyPressEvent(event);
        }
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
    location_of_last_propmpt = textCursor();
    location_of_last_propmpt.movePosition(QTextCursor::Left);
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
        current_command += str;
        if (endsWithNL || item < commands.count()-1) {
            newline (true);
            handleUserLineEntry();
            current_command = "";
        } else {
            length_of_input_buffer = str.length();
        }
    }

    moveCursor(QTextCursor::End);
    location_of_insertion_point = textCursor();
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
