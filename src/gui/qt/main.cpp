#include <QApplication>
#include <QEvent>
#include <unistd.h>
#include "baseobj.h"
#include "hyphymain.h"
#include "hyphyevents.h"
#include "hyphy_qt_helpers.h"

#include "HYSharedMain.h"

_String baseArgDir,
        libArgDir;
        
HyphyMain * _hyPrimaryConsoleWindow;

int main(int argc, char *argv[])
{
    GlobalStartup();
    DoApplicationSettings ();
    
    char    curWd[4096],
            dirSlash = GetPlatformDirectoryChar ();
    getcwd (curWd,4096);

    _String baseDir (curWd);

    if (baseDir.getChar (baseDir.sLength-1) != dirSlash) {
        baseDir=baseDir & dirSlash;
    }

    _String libDir (_HYPHY_LIBDIRECTORY_);

    if (libDir.getChar (libDir.sLength-1) != dirSlash) {
        libDir=libDir & dirSlash;
    }

    pathNames&& &libDir;

    libDirectory  = libDir;
    libArgDir     = libDirectory;
    baseDirectory = baseDir;
    baseArgDir    = baseDirectory;


    QApplication app(argc, argv);
    HyphyMain mainWindow;
    ReadInTemplateFiles     ();

    _hyPrimaryConsoleWindow = & mainWindow;
 
    /*mainWindow.setGeometry(
    QStyle::alignedRect(
        Qt::LeftToRight,
        Qt::AlignCenter,
        mainWindow.size(),
        qApp->desktop()->availableGeometry()
    ));*/
    
    mainWindow.show();
    mainWindow.raise();
    return app.exec();
}
