/*
    Console Window

    Sergei L. Kosakovsky Pond, December 2003.
*/

#ifndef _HYCONSOLEWINDOW_
#define _HYCONSOLEWINDOW_
//#pragma once
#include "HYTableWindow.h"

#define  HY_CONSOLE_CAN_COPY    0x01
#define  HY_CONSOLE_CAN_PASTE   0x02
#define  HY_CONSOLE_CAN_UNDO    0x04
#define  HY_CONSOLE_CAN_REDO    0x08

#define  HY_SL_FILE             0x01
#define  HY_SL_TASK             0x02
#define  HY_SL_TIMER            0x04
#define  HY_SL_PERCENT          0x08
#define  HY_SL_SUSPEND          0x10
#define  HY_SL_RESUME           0x20
#define  HY_SL_FORCE            0x40
#define  HY_SL_DONE             0x80

//__________________________________________________________________

class _HYConsoleWindow: public _HYTWindow
{

public:

    _HYConsoleWindow        (_String);
    virtual             ~_HYConsoleWindow       ();

    virtual bool        ProcessEvent            (_HYEvent*);
    virtual bool        ConfirmClose            (void);

    virtual void        _SetMenuBar             (void);
    virtual void        _UnsetMenuBar           (void);
    virtual bool        _ProcessMenuSelection   (long);
    virtual bool        _ProcessOSEvent         (Ptr);
    virtual void         _PaintStatusBar        (Ptr = nil, bool = false);
    void                _UpdateEditMenu         (void);

    virtual bool        IsSaveEnabled           (void) {
        return true;
    }
    virtual bool        IsPrintEnabled          (void) {
        return true;
    }
    void                PrintString             (_String&);
    _String*            ReadString              (void);
    void                SetFont                 (_HYFont);
    void                DoFind                  (bool);

    virtual void        DoSave                  (void);


    void                _DoFind                 (_String&);
    void                _DoSave                 (void);
    void                _DoPrint                (void);
    void                _DoPageSetup            (void);

    void                DoEcho                  (void);
    virtual bool        _Close                  (Ptr);

    FILE*               echoFileRef;
    char                echoStatus;

    _List               recentInputs,
                        recentSearches;

    _String             fileName,
                        action,
                        timer,
                        searchTerm;

    long                percentDone,
                        inputLocation;

    int                 editOptions;

    char                inputStatus;

};

extern  _String         cState,
        cTask,
        cInput;

extern  _HYConsoleWindow* hyphyConsoleWindow;

void    SetStatusLine               (_String);
void    SetStatusLine               (_String, _String, _String, long l);
void    SetStatusLine               (_String, _String, _String);
void    SetStatusLine               (_String, _String, _String, long, char);
void    SetStatusBarValue           (long, _Parameter, _Parameter);
void    SaveConsole                 (void);
#endif

//EOF
