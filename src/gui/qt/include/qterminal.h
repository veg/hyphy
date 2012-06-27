#pragma once

#include <QTextEdit>
#include <QStringList>

class QTerminal : public QTextEdit {
    Q_OBJECT

public:
    QTerminal(QWidget * parent = 0, Qt::WindowFlags f = 0);
    ~QTerminal();
    void keyPressEvent(QKeyEvent * event);
    void changeDir(const QString & dir);

private slots:
    void readBufferOut();
    void readStandardErr();

private:
    int inputCharCount;
    QTextCursor curCursorLoc;
    QString cmdStr;
    void printPrompt();
    QStringList cmdHistory;
    int histLocation;
    QString tempCmd;
};
