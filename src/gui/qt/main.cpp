#include <QApplication>
#include "hyphy_main.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    HyphyMain mainWindow;
    mainWindow.show();
    return app.exec();
}
