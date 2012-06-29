#include "hyphyevents.h"

QBufferToConsoleEvent::QBufferToConsoleEvent(QString bufferStr) : QEvent(BufferToStringType) {
    this->bufferStr = bufferStr;
}
