#include <QApplication>
#include <QEvent>
#include "baseobj.h"
#include "hyphymain.h"
#include "hyphyevents.h"

int main(int argc, char *argv[])
{
    GlobalStartup();

    QApplication app(argc, argv);
    HyphyMain mainWindow;
    mainWindow.show();
    return app.exec();
}
