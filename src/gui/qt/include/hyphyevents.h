#pragma once

#include <QEvent>
#include <QString>

class Q_GUI_EXPORT QBufferToConsoleEvent : public QEvent
{
public:
    QBufferToConsoleEvent(Type type, QString bufferStr);
    //~QBufferToConsoleEvent();
    inline QString buffer() const { return bufferStr; }
protected:
    QString bufferStr;
};

