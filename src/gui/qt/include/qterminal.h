#pragma once

#include <QTextBrowser>
#include <QStringList>

class QTerminal : public QTextBrowser {
    Q_OBJECT

public:
    QTerminal(QWidget * parent = 0, Qt::WindowFlags f = 0);
    ~QTerminal();
    void keyPressEvent(QKeyEvent * event);
    void contextMenuEvent(QContextMenuEvent* event);
    void insertFromMimeData(const QMimeData * source);
    void changeDir(const QString & dir);
    void prompt();

//private slots:
    //void find();

private:
    int histLocation;
    int inputCharCount;
    QTextCursor curCursorLoc;
    QString cmdStr;
    QString tempCmd;
    QStringList cmdHistory;
    void printPrompt();

};
