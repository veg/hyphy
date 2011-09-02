/*
    Object Inspector Panel

    Sergei L. Kosakovsky Pond, December 2000.
*/

#ifndef _HYOBJECTINSPECTOR_
#define _HYOBJECTINSPECTOR_
//#pragma once
#include "HYTableWindow.h"
#define  HY_OBJECT_INSPECTOR_TABLE_ROW 2

//__________________________________________________________________

class _HYObjectInspector: public _HYTWindow
{

public:

    _HYObjectInspector(void);

    virtual             ~_HYObjectInspector();

    virtual bool        ProcessEvent            (_HYEvent*);
    virtual bool        ProcessGEvent           (_HYEvent*);
    virtual void        Update                  (Ptr);
    virtual void        Paint                   (Ptr);
    virtual void        Activate                (void);
    void                SortObjectsByName       (long);
    virtual bool        ConfirmClose            (void);
//virtual bool      _ProcessOSEvent         (Ptr);

private:

    void                BuildListOfObjects      (long);
    void                UpdateButtonsAndInfo    (void);
    void                OpenObjectWindow        (void);
    void                KillObject              (void);
    void                OpenObject              (void);
    void                NewObject               (void);

    long                lastKillID;
    bool                firstTime;

};

extern  _String     objectInspectorTitle;
#endif

//EOF
