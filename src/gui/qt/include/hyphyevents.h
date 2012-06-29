#pragma once

#include <QEvent>
#include <QString>

const static QEvent::Type BufferToStringType = QEvent::Type(QEvent::User+1);

class Q_GUI_EXPORT QBufferToConsoleEvent : public QEvent
{
public:
    QBufferToConsoleEvent(QString bufferStr);
    //~QBufferToConsoleEvent();
    inline QString buffer() const { return bufferStr; }
protected:
    QString bufferStr;
};

