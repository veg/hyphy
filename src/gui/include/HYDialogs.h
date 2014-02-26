/*
    A collection of orphan dialogs

    Sergei L. Kosakovsky Pond, October 2001-May 2002.
*/

#ifndef _HYDIALOGS_
#define _HYDIALOGS_

#include "HYComponent.h"
#include "HYWindow.h"
#include "HYTableWindow.h"
#include "preferences.h"

//__________________________________________________________________

class _HYFontDialog: public _HYTWindow
{

public:

    _HYFontDialog    (_HYFont&,_HYWindow*);
    virtual     ~_HYFontDialog   (void) {}

    virtual bool ProcessEvent    (_HYEvent*);


    _HYFont     myFont,
                firstFont;

    _HYWindow*  mr;

};

//__________________________________________________________________

class _HYListSelectDialog: public _HYTWindow
{

public:

    _HYListSelectDialog  (_List*, _SimpleList*, _SimpleList*, _String, _SimpleList*, long*);
    virtual     ~_HYListSelectDialog     (void) {}

    virtual bool ProcessEvent            (_HYEvent*);
    virtual void SetInitialSelection     (void);


    _List       * data,
                dData;
    _SimpleList * choices,
                * validChoices,
                * selections;

    long        * result;
    _String     lastString;

};

//__________________________________________________________________

class _HYSimpleListSelectDialog: public _HYTWindow
{

public:

    _HYSimpleListSelectDialog    (_List*, _SimpleList*, _SimpleList*, _String, _SimpleList*, long, long*,Ptr = nil);
    virtual     ~_HYSimpleListSelectDialog   (void) {}

    virtual bool ProcessEvent                (_HYEvent*);
    virtual void SetInitialSelection         (void);


    _List       dData;
    _SimpleList * choices,
                * validChoices,
                * selections;

    long        * result,
                reqSel;

    _String     lastString;

};

//__________________________________________________________________

class _HYPreferencesDialog: public _HYTWindow
{

public:

    _HYPreferencesDialog         (_List&, _String, bool, bool*);
    virtual     ~_HYPreferencesDialog        (void)
    {}

    virtual bool ProcessEvent                (_HYEvent*);
    virtual void SetInitialSelection         (void);


    _List    dData,
             *prefList,
             backupValues;

    bool     setFont,
             *result;
};

//__________________________________________________________________

class _HYCombDialog: public _HYTWindow
{

public:

    _HYCombDialog        (_SimpleList*, bool*, Ptr = nil);
    virtual     ~_HYCombDialog       (void)
    {}

    virtual bool ProcessEvent                (_HYEvent*);


    _SimpleList*
    combResult;

    bool    *result;
};


//__________________________________________________________________

class _HYPartitionDialog: public _HYTWindow
{

public:

    _HYPartitionDialog       (_String&, _String*, _HYColor*, bool*, char*, long*, long,  bool*, char, Ptr);
    virtual     ~_HYPartitionDialog      (void)
    {}

    virtual bool ProcessEvent            (_HYEvent*);


    _String*  partName;
    _HYColor* color;

    bool*     direction,
              *   result;

    long*     codeRef,
              dfID;

    char*     readFrame;
};


//__________________________________________________________________

class _HYTextDialog: public _HYTWindow
{

public:

    _HYTextDialog        (_String*,_String&, bool*, Ptr = nil, _hyStringValidatorType = nil);
    virtual     ~_HYTextDialog       (void)  {}

    virtual bool ProcessEvent                (_HYEvent*);

    //**********************************************//

    bool    * result,
            lastValidationState;

    _hyStringValidatorType  validator;
    _String * textOut;
};

//__________________________________________________________________

class _HYTextDialog2: public _HYTextDialog
{

public:

    _HYTextDialog2       (_String*,_String*,_String&,_String&, bool*, Ptr = nil);
    virtual     ~_HYTextDialog2      (void)  {}

    virtual bool ProcessEvent                (_HYEvent*);

    _String * textOut2;
};

//__________________________________________________________________

class _HYTextDialogWithCheckbox: public _HYTextDialog
{

public:

    _HYTextDialogWithCheckbox        (_String*,_String&, _String& ,bool*, bool*, Ptr = nil);
    virtual     ~_HYTextDialogWithCheckbox       (void)
    {}

    virtual bool ProcessEvent                    (_HYEvent*);

    bool    * checkState;
};

//__________________________________________________________________

class _HYTextDialogWithPulldown: public _HYTextDialogWithCheckbox
{

public:

    _HYTextDialogWithPulldown        (_String*,_String&, _String& , _List&, _List*, bool*, bool*, long*, long, Ptr = nil);
    virtual     ~_HYTextDialogWithPulldown       (void)
    {}

    virtual bool ProcessEvent                    (_HYEvent*);

    long      * menuSelection;
    _List     * autoFill;
};

//__________________________________________________________________

class _HYProceedPromptBox: public _HYTWindow
{

public:

    _HYProceedPromptBox          (_String&,_String&, _String&,  _HYFont&, long, bool*, Ptr);
    virtual     ~_HYProceedPromptBox             (void)
    {}

    virtual bool ProcessEvent                    (_HYEvent*);

    bool    * result;
};

//__________________________________________________________________

class _HYProceedPromptBoxWCheck: public _HYProceedPromptBox
{

public:

    _HYProceedPromptBoxWCheck            (_String&, _String&, _String&, _String&, _HYFont&, long, bool*, bool*, Ptr);
    virtual     ~_HYProceedPromptBoxWCheck           (void) {}

    virtual bool ProcessEvent                    (_HYEvent*);

    bool    * checkState;
};

//__________________________________________________________________

class _HYInferenceConstraints: public _HYTextDialog
{

public:

    _HYInferenceConstraints          (_String&, _List&, _List&, _List&, _List&,
                                      _List&, _List&, _List&, _List&,
                                      _HYFont&, bool*, _String*, Ptr);
    virtual     ~_HYInferenceConstraints             (void)
    {}

    virtual bool ProcessEvent                        (_HYEvent*);

    bool    * checkState;
    _String * topC;

    _List   * il1,
            * il2,
            * il3,
            * il4;
};


//__________________________________________________________________

class _HYCIDialog: public _HYTWindow
{

public:

    _HYCIDialog          (bool*, char*, _Parameter*,_Parameter*, Ptr);
    virtual     ~_HYCIDialog         (void) {}

    virtual bool ProcessEvent        (_HYEvent*);


    bool        *result;

    char        *option,
                lastState;

    _Parameter  *sigValue,
                *sigValue2;
};

//__________________________________________________________________

class _HYLProfDialog: public _HYTWindow
{

public:

    _HYLProfDialog           (bool*,_Parameter*,_Parameter*,long*,bool*,_Parameter,_Parameter,_Parameter,Ptr);
    virtual     ~_HYLProfDialog          (void) {}

    virtual bool ProcessEvent            (_HYEvent*);


    bool        *result,
                *doNorm;

    char        lastState[3];

    _Parameter  *left,
                *right,
                mle,
                vlb,
                vub;

    long*       intervals;

};

//__________________________________________________________________

class _HYNewChart: public _HYTWindow
{

public:

    _HYNewChart  (long*,long*,bool*,_String*);
    virtual     ~_HYNewChart     (void)
    {}

    virtual bool ProcessEvent    (_HYEvent*);

    bool    * result;
    long    * resR,
            * resC;

    _String * title;
};

//__________________________________________________________________

class _HYTreePrintPrefs: public _HYTWindow
{

public:

    _HYTreePrintPrefs    (long*,long*,bool*,Ptr);
    virtual     ~_HYTreePrintPrefs   (void) {}

    virtual bool ProcessEvent    (_HYEvent*);

    bool    * result;
    long    * resW,
            * resH;
};


//__________________________________________________________________

long     HandleHierListSelection (_List&, _SimpleList&, _SimpleList&, _String, _SimpleList&, long);
long     HandleListSelection     (_List&, _SimpleList&, _SimpleList&, _String, _SimpleList&, long, Ptr = nil);
long     HandleListSelection     (_List&, _String, Ptr = nil);
bool     HandlePreferences       (_List&, _String , bool = true);

_String  NewTreeWindow           (long);
void     BuildBalancedTree       (long, long , _String& , _List&);

extern  _List                   cachedDialogTitles,
        cachedDialogSelections;

_String*RetrieveCachedSelection (_String*);
void    StoreDialogSelection    (_String*,_String*);
bool    EnterStringDialog       (_String& res, _String& prompt, Ptr = nil, _hyStringValidatorType = nil);
bool    EnterString2Dialog      (_String& , _String&, _String&, _String& , Ptr = nil);
bool    EnterStringDialogWithCheckbox
(_String&, _String&, _String&, bool&, Ptr = nil);
bool    EnterStringDialogWithPulldown
(_String&, _String&, _String&, long&, _List&, _List *, bool&, long,  Ptr);

bool    ProceedPrompt           (_String&, Ptr = nil);
bool    ProceedPromptWithCheck  (_String&, _String&, bool&, Ptr = nil);
void    ProblemReport           (_String&, Ptr = nil);
bool    TreePrintSetup          (long&, long&, Ptr);
bool    HandleCIDialog          (char&, _Parameter&, _Parameter&, Ptr);
bool    HandlePLDialog          (_Parameter&, _Parameter&, long&, bool&, _Parameter, Ptr);
void    ShowMessagesLog         (void);
void    ExecuteSelection        (void);
bool    OpenDatabaseFile        (_String*);
bool    NewDatabaseFile         (_String*);


#endif

//EOF
