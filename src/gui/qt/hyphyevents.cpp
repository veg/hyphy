#include "hyphyevents.h"

QBufferToConsoleEvent::QBufferToConsoleEvent(Type type, QString bufferStr) : QEvent(QEvent::User) {
    this->bufferStr = bufferStr;
}
