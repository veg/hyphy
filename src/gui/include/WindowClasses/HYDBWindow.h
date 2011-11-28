/*
    SQLite DB Viewer/Editor Widget

    Sergei L. Kosakovsky Pond, October 2005.
*/

#ifndef _HYDBWINDOW_
#define _HYDBWINDOW_
//#pragma once

#include "HYTableWindow.h"
#include "HYTableComponent.h"
#include "sqlite3.h"


//__________________________________________________________________

class _HYDBWindow: public _HYTWindow
{

public:

    _HYDBWindow             (_String,_String* = nil);
    virtual             ~_HYDBWindow            (void);

    virtual bool        ProcessEvent            (_HYEvent*);

    /*virtual void      _SetMenuBar             (void);
    virtual void        _UnsetMenuBar           (void);
    virtual bool        _ProcessMenuSelection   (long);*/

    virtual bool        IsSaveEnabled           (void) {
        return false;
    }
    virtual bool        IsPrintEnabled          (void) {
        return true;
    }

    void                SetDB                   (_String);
    void                ScanTables              (bool = false);
    void                CloseCurrentDB          (void);
    void                UpdateStatusBar         (void);
    _List*              ExecuteSQLBlurb         (_String&);
    bool                LoadTable               (long,_String* = nil);
    void                EditLongEntry           (void);
    _List*              RetrieveRecord          (long, _String&, bool, _String* = nil);
    void                RunSQLQuery             (void);
    virtual bool        _ProcessMenuSelection   (long);

    long                dbID,
                        currentTable,
                        lengthLimit,
                        recordLimit;

    _HYColor            backColor1,
                        backColor2;

    _HYFont             labelFont1,
                        labelFont2,
                        labelFont3,
                        headerFont;

    _List               tableList,
                        sqlHistory;
};

//__________________________________________________________________

void        ReadSQLPlugins                      (void);
bool        ExecuteSomeCodeAndCheck             (_String&);
int         _HYDBWCallBack                      (void*,int, char**, char**);
int         _HYDBWTablePopulatorCallBack        (void*,int, char**, char**);
int         _HYDBWTableCellCounter              (void*,int, char**, char**);

#endif

//EOF
