#pragma once

#include <QTextBrowser>
#include <QStringList>

class QTerminal : public QTextBrowser {
    Q_OBJECT

public:
    QTerminal(QWidget * parent = 0, Qt::WindowFlags f = 0);
    ~QTerminal();
    
    void keyPressEvent          (QKeyEvent * event);
    void contextMenuEvent       (QContextMenuEvent* event);
    void insertFromMimeData     (const QMimeData * source);
    void changeDir              (const QString & dir);
        // SLKP 20121205: seems unused
    
    void prompt                 (bool = false);
    void newline                (bool = false);
    
signals:
    void userEnteredString (const QString);

private slots:
    
    void handleCursorMove  ();
    //void find();

private:
    QStringList command_history;
    
    int         command_history_location,
                length_of_input_buffer;
 
    QTextCursor location_of_last_propmpt,
                location_of_insertion_point;
    
    QString     current_command,
                stash_current_command;
    
    void        handleUserLineEntry     (void);
    void        beepAndIgnoreEvent      (QEvent* event, bool do_ignore = true);
    void        replaceCurrentCommand   (const QString* replace_with, bool stash_current);
};
