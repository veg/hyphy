/*
 
 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (spond@ucsd.edu)
 Art FY Poon    (apoon@cfenet.ubc.ca)
 Steven Weaver (sweaver@ucsd.edu)
 
 Module Developers:
 Lance Hepler (nlhepler@gmail.com)
 Martin Smith (martin.audacis@gmail.com)
 
 Significant contributions from:
 Spencer V Muse (muse@stat.ncsu.edu)
 Simon DW Frost (sdf22@cam.ac.uk)
 
 Permission is hereby granted, free of charge, to any person obtaining a
 copy of this software and associated documentation files (the
 "Software"), to deal in the Software without restriction, including
 without limitation the rights to use, copy, modify, merge, publish,
 distribute, sublicense, and/or sell copies of the Software, and to
 permit persons to whom the Software is furnished to do so, subject to
 the following conditions:
 
 The above copyright notice and this permission notice shall be included
 in all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 
 */

#include "batchlan.h"
#include "HYDBWindow.h"

#include "HYLabel.h"
#include "HYPullDown.h"
#include "HYButtonBar.h"
#include "HYTextBox.h"
#include "HYGraphicPane.h"
#include "HYModelWindow.h"
#include "HYUtils.h"
#include "HYDialogs.h"
#include "HYEventTypes.h"

#include "string.h"

#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif


extern      _SimpleList                         windowObjects,
            sqlDatabases;

extern      _HYColor                            tableDefaultBk2,
            labelBColor,
            labelFColor,
            black;

extern      _String                             donotWarnAgain,
            windowTypeTable,
            windowTypeDistribTable,
            none;


#define     HY_LARGE_COMPONENT_SIZE             10000L
#define     HY_DBW_DEFAULT_TRUNCATE_TEXT        50
#define     HY_DBW_DEFAULT_RECORD_LIMIT         1000

_String     internalSQLiteIndex                 ("__HYPHY_INTERNAL_SQLITE_INDEX"),
            _HY_DBW_longFieldEdit               ("Edit a long database field"),
            _HY_DBW_TemplateFileName            ("HY_DBW_TemplateList"),
            _HY_DBW_TablePlaceholder            ("_HY_DBW_TABLE_NAME_"),
            _HY_DBW_ColumnPlaceholder           ("_HY_DBW_COLUMN_NAMES_"),
            _HY_DBW_RowValues                   ("_HY_DBW_ROW_VALUES_"),
            _HY_DBW_OutputProcessorFunction     ("_HY_DBW_OUTPUT_PROCESSOR_FUNCTION_"),
            _HY_DBW_OutputProcessorName         ("_HY_DBW_OUTPUT_PROCESSOR_NAME_"),
            _HY_DBW_OutputDBID                  ("_HY_DBW_OUTPUT_DB_ID_"),
            _HY_DBW_SqlQuery                    ("_HY_DBW_SQL_QUERY"),
            _HY_DBW_RunMe                       ("_HY_DBW_OUTPUT_RUN_ME_");


_List       _HYDBW_Templates,
            _HYDBW_OutputProcessors,
            _HYDBW_FilePaths;


//____________________________________________________________________________________

int  _HYDBWCallBack (void* aL,int cc, char** rd, char** cn)
{
    _List * aList = (_List*)aL;

    if (aList && cc && aList) {
        if (aList->lLength == 0) { // first call; populate column headers
            for (long cnt = 0; cnt < cc; cnt++) {
                _List * newList     = new _List (2L),
                * emptyList = new _List ();

                checkPointer (newList);
                checkPointer (emptyList);

                if (cn[cnt]) {
                    _String * colData = new _String (cn[cnt]);
                    (*newList) << colData;
                    DeleteObject (colData);
                } else {
                    (*newList) && & empty;
                }

                (*newList) << emptyList;
                DeleteObject (emptyList);
                (*aList) << newList;
                DeleteObject (newList);
            }

        }
        for (long cnt = 0; cnt < cc; cnt++) {
            _List * storeIn = (_List *)(*((_List*)(*aList) (cnt)))(1);

            if (rd[cnt]) {
                _String * colData = new _String (rd[cnt]);
                (*storeIn) << colData;
                DeleteObject (colData);
            } else {
                (*storeIn) && & empty;
            }
        }

    }
    return 0;
}

//____________________________________________________________________________________

int  _HYDBWTablePopulatorCallBack (void* aL,int cc, char** rd, char** cn)
{
    _HYDBWindow * aWindow = (_HYDBWindow*)aL;

    if (cc && aWindow) {
        _HYTable * ttop   = (_HYTable*)aWindow->GetObject (1),
                   * tleft  = (_HYTable*)aWindow->GetObject (2),
                     * tmain  = (_HYTable*)aWindow->GetObject (0);


        if (ttop->horizontalSpaces.lLength < cc) {
            long defWidth  = (aWindow->componentR[1]-aWindow->componentL[1]+1),
                 defCWidth = aWindow->componentR[2]-aWindow->componentL[2]+1;

            for (long cnt = 0; cnt < cc; cnt++) {
                tmain->AddColumn (-1,defWidth/cc,HY_TABLE_STATIC_TEXT);
                ttop->AddColumn  (-1,defWidth/cc,HY_TABLE_STATIC_TEXT);
            }


            ttop-> AddRow (-1,18,HY_TABLE_STATIC_TEXT|HY_TABLE_BEVELED);
            tleft->AddColumn (-1,defCWidth,HY_TABLE_STATIC_TEXT);


            for (long cnt = 0; cnt < cc; cnt++) {
                if (cn[cnt]) {
                    _String * colData = new _String (cn[cnt]);
                    ttop->SetCellData (colData,0,cnt,HY_TABLE_STATIC_TEXT|HY_TABLE_BEVELED,false);
                } else {
                    ttop->SetCellData (&empty,0,cnt,HY_TABLE_STATIC_TEXT|HY_TABLE_BEVELED,true);
                }

                ttop->AutoFitColumn (cnt,true);
            }

            ttop  -> EnforceWidth   (defWidth,-1,true);

        }

        if (tleft->verticalSpaces.lLength < aWindow->recordLimit) {
            tleft-> AddRow (-1,tleft->textFont.size+6,HY_TABLE_STATIC_TEXT|HY_TABLE_BEVELED);
            tmain-> AddRow (-1,tleft->textFont.size+6,HY_TABLE_STATIC_TEXT);

            _String rowIndex ((long)tleft->verticalSpaces.lLength);
            tleft->SetCellData (&rowIndex,tleft->verticalSpaces.lLength-1,0,HY_TABLE_STATIC_TEXT|HY_TABLE_BEVELED,true);

            for (long cnt = 0; cnt < cc; cnt++) {
                if (rd[cnt]) {
                    long ll = strlen (rd[cnt]);

                    if (ll>HY_DBW_DEFAULT_TRUNCATE_TEXT) {
                        _String * cellData = new _String(HY_DBW_DEFAULT_TRUNCATE_TEXT,true);
                        checkPointer (cellData);
                        for (long k=0; k<HY_DBW_DEFAULT_TRUNCATE_TEXT; k++) {
                            (*cellData) << rd[cnt][k];
                        }
                        (*cellData) << "...";
                        cellData->Finalize();
                        tmain->SetCellData (cellData,tleft->verticalSpaces.lLength-1,cnt,HY_TABLE_STATIC_TEXT,false);
                    } else {
                        _String * cellData = new _String(rd[cnt]);
                        checkPointer (cellData);
                        tmain->SetCellData (cellData,tleft->verticalSpaces.lLength-1,cnt,HY_TABLE_STATIC_TEXT,false);
                    }
                }
            }
        }
    }
    return 0;
}

//____________________________________________________________________________________

int  _HYDBWTableCellCounter (void* aL,int cc, char**, char**)
{
    long * aCounter = (long*)aL;

    if (aCounter) {
        aCounter[0] ++;
        aCounter[1] = cc;
    }

    return 0;
}

//__________________________________________________________

_HYDBWindow::_HYDBWindow (_String name, _String *filePath):_HYTWindow (name, true)
{

    if (!_HYDBW_Templates.lLength && !_HYDBW_OutputProcessors.lLength) {
        ReadSQLPlugins();
    }

    _HYRect         canvasSettings = {30,50,30,50,HY_COMPONENT_NO_SCROLL};

    _HYLabel*       l1      = new _HYLabel (canvasSettings, GetOSWindowData());
    _HYLabel*       l2      = new _HYLabel (canvasSettings, GetOSWindowData());

    canvasSettings.left     = 100;
    canvasSettings.right    = 150;

    _HYPullDown*    p1      = new _HYPullDown (canvasSettings,GetOSWindowData());
    _HYPullDown*    p2      = new _HYPullDown (canvasSettings,GetOSWindowData());

    canvasSettings.left     = 100;
    canvasSettings.right    = 100;

    canvasSettings.top      = 200;
    canvasSettings.bottom   = HY_LARGE_COMPONENT_SIZE;
    canvasSettings.right    = HY_LARGE_COMPONENT_SIZE;

    canvasSettings.width    = HY_COMPONENT_V_SCROLL|HY_COMPONENT_H_SCROLL|HY_TABLE_DONT_GROW_HORIZ;
    _HYTable*   table       = new _HYTable (canvasSettings,GetOSWindowData(),1,1,100,16,HY_TABLE_STATIC_TEXT);

    canvasSettings.top      = canvasSettings.bottom     = 20;
    canvasSettings.width    = HY_COMPONENT_NO_SCROLL|HY_TABLE_DONT_GROW_HORIZ;

    _HYTable*   tableHead   = new _HYTable (canvasSettings,GetOSWindowData(),1,1,100,18,HY_TABLE_STATIC_TEXT);
    canvasSettings.left     = canvasSettings.right = 40;
    canvasSettings.top      = 200;
    canvasSettings.bottom   = HY_LARGE_COMPONENT_SIZE;

    canvasSettings.width    = HY_COMPONENT_NO_SCROLL|HY_COMPONENT_BORDER_R|HY_TABLE_DONT_GROW_HORIZ|HY_TABLE_DONT_GROW_VERT;
    _HYTable*   tableLeft   = new _HYTable (canvasSettings,GetOSWindowData(),1,1,40,16,HY_TABLE_STATIC_TEXT);

    table->selectionType    = HY_TABLE_FOCUSABLE|HY_TABLE_HORZ_STRETCH|HY_TABLE_VERT_STRETCH;
    tableHead->selectionType= HY_TABLE_HORZ_STRETCH;
    tableLeft->selectionType= HY_TABLE_DONT_GROW_VERT|HY_TABLE_VERT_STRETCH;

    canvasSettings.top      = canvasSettings.bottom     = 18;
    canvasSettings.width    = HY_COMPONENT_NO_SCROLL|HY_COMPONENT_BORDER_T|HY_COMPONENT_BORDER_R;
    _HYLabel*       l3      = new _HYLabel (canvasSettings, GetOSWindowData());

    canvasSettings.top      = canvasSettings.bottom     = HY_SCROLLER_WIDTH+1;
    canvasSettings.width    = HY_COMPONENT_NO_SCROLL|HY_COMPONENT_BORDER_T|HY_COMPONENT_BORDER_R;
    _HYLabel*       l4      = new _HYLabel (canvasSettings, GetOSWindowData());

    canvasSettings.left     = 100;
    canvasSettings.right    = HY_LARGE_COMPONENT_SIZE;
    canvasSettings.bottom   = 61;
    canvasSettings.top      = 61;

    canvasSettings.width    = HY_COMPONENT_NO_SCROLL;

    _HYTextBox  * sqlView   = new _HYTextBox (canvasSettings,GetOSWindowData(),false);

    canvasSettings.left = canvasSettings.right = 50;
    _HYButtonBar*   bb2     = new _HYButtonBar (canvasSettings,GetOSWindowData());

    _String btnTT ("Run typed SQL command");
    bb2->SetButtonDim     (16);
    bb2->SetButtonLayoutW (2);
    bb2->AddButton (ProcureIconResource(HY_MODELWINDOW_ICON_BASE+20),&btnTT);
    btnTT = "Paste a predefined SQL command";
    bb2->AddButton (ProcureIconResource(6042),&btnTT);
    btnTT = "SQL Command History";
    bb2->AddButton (ProcureIconResource(6006),&btnTT);
    bb2->MarkAsPullDown (1,true);
    bb2->MarkAsPullDown (2,true);
    bb2->EnableButton  (0,false);
    bb2->EnableButton  (1,false);
    bb2->EnableButton  (2,false);

    table->SetMessageRecipient      (this);
    tableHead->SetMessageRecipient  (this);
    tableLeft->SetMessageRecipient  (this);
    p1->SetMessageRecipient         (this);
    p2->SetMessageRecipient         (this);
    sqlView->SetMessageRecipient    (this);
    bb2->SetMessageRecipient        (this);

    AddObject (table);     // 0
    AddObject (tableHead); // 1
    AddObject (tableLeft); // 2
    AddObject (p1);        // 3
    AddObject (p2);        // 4
    AddObject (l1);        // 5
    AddObject (l2);        // 6
    AddObject (bb2);       // 7
    AddObject (sqlView);   // 8
    AddObject (l3);        // 9
    AddObject (l4);        // 10

    SetTableDimensions (5,5);

    SetCell   (0,0,l1);
    SetCell   (0,1,p1);
    SetCell   (0,2,sqlView);
    SetCell   (0,3,sqlView);
    SetCell   (0,4,bb2);

    SetCell   (1,0,l2);
    SetCell   (1,1,p2);
    SetCell   (1,2,sqlView);
    SetCell   (1,3,sqlView);
    SetCell   (1,4,bb2);

    SetCell   (2,0,l3);
    SetCell   (2,1,tableHead);
    SetCell   (2,2,tableHead);
    SetCell   (2,3,tableHead);
    SetCell   (2,4,tableHead);

    SetCell   (3,0,tableLeft);
    SetCell   (3,1,table);
    SetCell   (3,2,table);
    SetCell   (3,3,table);
    SetCell   (3,4,table);

    SetCell   (4,0,l4);
    SetCell   (4,1,table);
    SetCell   (4,2,table);
    SetCell   (4,3,table);
    SetCell   (4,4,table);

    p1->SetBackColor (labelBColor);
    p2->SetBackColor (labelBColor);

    l1->SetForeColor (labelFColor);
    l2->SetForeColor (labelFColor);

    l1->SetBackColor (labelBColor);
    l2->SetBackColor (labelBColor);
    l3->SetBackColor (tableHead->backColor2);
    l4->SetBackColor (tableHead->backColor2);

    sqlView->SetBackColor (labelBColor);

    bb2->SetBackColor (labelBColor);

    _HYFont  labelFont;

#ifdef __WINDOZE__
    labelFont.face = "Arial";
    labelFont.size = 14;
#else
    labelFont.face = "Geneva";
    labelFont.size = 10;
#endif

    labelFont.style = HY_FONT_PLAIN;

    table->SetFont(labelFont);
    tableHead->SetFont(labelFont);
    tableLeft->SetFont(labelFont);

#ifdef __WINDOZE__
    labelFont.size = 14;
#else
    labelFont.size = 12;
#endif

    p1->SetAlignFlags (HY_ALIGN_LEFT);
    p2->SetAlignFlags (HY_ALIGN_LEFT);

    p1->AddMenuItem (none,-1);
    p2->AddMenuItem ("Table",-1);

    if (_HYDBW_OutputProcessors.lLength) {
        p2->AddMenuItem (menuSeparator,-1);
        for (long pi = 0; pi < _HYDBW_OutputProcessors.lLength; pi++) {
            p2->AddMenuItem (*(_String*)_HYDBW_OutputProcessors(pi),-1);
        }
    }

    l1->SetAlignFlags (HY_ALIGN_LEFT);
    l2->SetAlignFlags (HY_ALIGN_LEFT);
    l3->SetAlignFlags (HY_ALIGN_LEFT);

    l1->SetShadow (true);
    l2->SetShadow (true);

    l1->SetFont (labelFont);
    l2->SetFont (labelFont);
    l3->SetFont (tableHead->textFont);

    l1->SetText ("Table:");
    l2->SetText ("Output:");
    l3->SetText ("Index");

    backColor1.R = backColor1.G = backColor1.B =
                                      backColor2.R = backColor2.G = backColor2.B = 255;

    labelFont1.face  = "Times";
    labelFont1.style = HY_FONT_PLAIN;
#ifdef __WINDOZE__
    labelFont1.size = 14;
#else
    labelFont1.size = 12;
#endif

    labelFont2.face  = "Times";
    labelFont2.style = HY_FONT_PLAIN;
#ifdef __WINDOZE__
    labelFont2.size = 12;
#else
    labelFont2.size = 10;
#endif

    labelFont3.face  = "Times";
    labelFont3.style = HY_FONT_ITALIC;
#ifdef __WINDOZE__
    labelFont3.size = 14;
#else
    labelFont3.size = 12;
#endif

    headerFont.face  = "Times";
    headerFont.style = HY_FONT_ITALIC;
#ifdef __WINDOZE__
    labelFont.size = 20;
#else
    labelFont.size = 18;
#endif

    dbID         = -1;
    currentTable = -1;
    lengthLimit  = HY_DBW_DEFAULT_TRUNCATE_TEXT;
    recordLimit  = HY_DBW_DEFAULT_RECORD_LIMIT;

    SetPosition        (70,70);
    SetWindowRectangle (0,0,400,400);

    if (filePath) {
        SetDB (*filePath);
    }

    //LoadTable (0);

    sqlView->SetText ("SQL Query");

    DeleteObject (table);
    DeleteObject (tableHead);
    DeleteObject (tableLeft);
    DeleteObject (l1);
    DeleteObject (l2);
    DeleteObject (l3);
    DeleteObject (l4);
    DeleteObject (bb2);
    DeleteObject (p1);
    DeleteObject (p2);
    DeleteObject (sqlView);
}

//__________________________________________________________
_HYDBWindow::~_HYDBWindow()
{
    CloseCurrentDB ();
}

//__________________________________________________________
void    _HYDBWindow::SetDB (_String filePath)
{
    CloseCurrentDB ();

    _String             openCommands = _String("DoSQL(SQL_OPEN,\"") &filePath & "\"," & internalSQLiteIndex & ");";

    if (ExecuteSomeCodeAndCheck (openCommands)) {
        _Parameter db;
        checkParameter  (internalSQLiteIndex, db, -1.);
        dbID = db;
        ScanTables (true);
    }
}

//__________________________________________________________
void    _HYDBWindow::CloseCurrentDB (void)
{
    if (dbID >= 0) {
        _String             closeCommands = _String("DoSQL(SQL_CLOSE,\"\",") & internalSQLiteIndex & ");";
        ExecuteSomeCodeAndCheck (closeCommands);
        tableList.Clear();
        currentTable    = -1;
        dbID            = -1;
    }
}

//__________________________________________________________
void    _HYDBWindow::UpdateStatusBar (void)
{
    _String statusLine;

    if (dbID >= 0) {
        _HYTable        * tmain  = (_HYTable*)GetObject (0);

        statusLine = _String ("Tables: ") & (long)tableList.lLength & ". Displaying " & (long)tmain->horizontalSpaces.lLength & " fields for " & (long)(tmain->verticalSpaces.lLength-1) & " records.";
    } else {
        statusLine = "No database is currently loaded";
    }

    SetStatusBar (statusLine);
}

//__________________________________________________________
void    _HYDBWindow::ScanTables (bool forceCreate)
{
    _String getTables ("SELECT name FROM sqlite_master\nWHERE type='table'\nORDER BY name"),
            curTableName;

    _List * checkTables = ExecuteSQLBlurb (getTables);

    if (currentTable >= 0 && currentTable < tableList.lLength) {
        curTableName = *(_String*)tableList (currentTable);
    }

    tableList.Clear();

    _HYPullDown * p1 = (_HYPullDown*) GetObject (3);
    p1->DeleteAllItems();
    if (checkTables->lLength) {
        tableList.Duplicate((*((_List*)(*checkTables) (0)))(1));
        if (tableList.lLength) {
            for (long mi = 0; mi < tableList.lLength; mi++) {
                p1->AddMenuItem (*(_String*)tableList(mi), -1);
            }

            currentTable = tableList.Find (&curTableName);
            if (currentTable >= 0) {
                p1->ChangeSelection (currentTable,false);
            } else {
                p1->ChangeSelection (0,true);
            }
        }
    } else {
        if (forceCreate) {
            _String createTable ("CREATE TABLE TABLE_ID (\na INTEGER PRIMARY KEY,\nb INTEGER);"),
                    createTableP ("Enter SQL instructions to make a table");

            /*BringToFront();*/
            if (EnterStringDialog (createTable,createTableP,nil)) {
                DeleteObject (ExecuteSQLBlurb (createTable));
                ScanTables (true);
                return;
            }
        }

        getTables = "Could not load any SQLite tables";
        ProblemReport (getTables,(Ptr)this);
        p1->AddMenuItem     (none,-1);
        p1->ChangeSelection (0,false);

    }

    DeleteObject (checkTables);
    UpdateStatusBar ();
}

//__________________________________________________________
bool    _HYDBWindow::LoadTable (long tableIndex, _String* code)
{
    bool result = false;
    if (dbID >= 0 && tableIndex < tableList.lLength && tableIndex >= 0) {
        _HYTable        * ttop   = (_HYTable*)GetObject (1),
                          * tleft  = (_HYTable*)GetObject (2),
                            * tmain  = (_HYTable*)GetObject (0);

        _HYButtonBar    * bb2    = (_HYButtonBar*)GetObject (7);

        bb2->EnableButton(0,false);
        bb2->EnableButton(1,false);

        ttop->ClearTable  (true);
        tleft->ClearTable (true);
        tmain->ClearTable (true);

        _String getTables    = code ? *code:(_String("SELECT * FROM ") & *(_String*)tableList (tableIndex));

        char  * errMsg         = nil;
        long    recordCount[2] = {0,0};

        sqlite3_exec((sqlite3*)sqlDatabases.lData[dbID], getTables.sData, _HYDBWTableCellCounter, (Ptr)recordCount, &errMsg);

        if (recordCount [0]) {
            if (recordCount [0] > recordLimit) {
                _String     res (recordLimit),
                            prompt = _String("Large table warning: ") & recordCount [0] & " rows. Display:";

                tmain->AddColumn (-1,10,HY_TABLE_STATIC_TEXT|HY_TABLE_CANTSELECT);
                ttop ->AddColumn (-1,10,HY_TABLE_STATIC_TEXT|HY_TABLE_CANTSELECT);
                tleft->AddColumn (-1,tleft->settings.left,HY_TABLE_STATIC_TEXT|HY_TABLE_CANTSELECT);
                ttop ->AddRow    (-1,1,HY_TABLE_STATIC_TEXT|HY_TABLE_CANTSELECT);
                tmain -> AddRow (-1,1,HY_TABLE_STATIC_TEXT|HY_TABLE_CANTSELECT);
                tleft -> AddRow (-1,1,HY_TABLE_STATIC_TEXT|HY_TABLE_CANTSELECT);

                if (EnterStringDialog (res,prompt,nil)) {
                    recordLimit = res.toNum();
                    if (recordLimit < 1) {
                        recordLimit = HY_DBW_DEFAULT_RECORD_LIMIT;
                    }
                }
                ttop->ClearTable  (true);
                tleft->ClearTable (true);
                tmain->ClearTable (true);
            }
            recordCount[0] = MIN (recordCount[0],recordLimit);

            tmain -> RequestSpace (recordCount[0], recordCount[1]);
            tleft -> RequestSpace (recordCount[0], 1);
            sqlite3_exec((sqlite3*)sqlDatabases.lData[dbID], getTables.sData, _HYDBWTablePopulatorCallBack, (Ptr)this, &errMsg);

        } else if (tmain->horizontalSpaces.lLength == 0 && tmain->verticalSpaces.lLength == 0) {
            tmain->AddColumn (-1,10,HY_TABLE_STATIC_TEXT|HY_TABLE_CANTSELECT);
            ttop ->AddColumn (-1,10,HY_TABLE_STATIC_TEXT|HY_TABLE_CANTSELECT);
            tleft->AddColumn (-1,tleft->settings.left,HY_TABLE_STATIC_TEXT|HY_TABLE_CANTSELECT);
            ttop ->AddRow    (-1,1,HY_TABLE_STATIC_TEXT|HY_TABLE_CANTSELECT);
        }

        tmain -> AddRow (-1,1,HY_TABLE_STATIC_TEXT|HY_TABLE_CANTSELECT);
        tleft -> AddRow (-1,1,HY_TABLE_STATIC_TEXT|HY_TABLE_CANTSELECT);

        tmain -> AutoFitWidth   (*ttop);
        ttop  -> SetVisibleSize (ttop->rel);
        tleft -> SetVisibleSize (tleft->rel);
        tmain -> SetVisibleSize (tmain->rel);
        ttop  -> _MarkForUpdate ();
        tleft -> _MarkForUpdate ();
        tmain -> _MarkForUpdate ();

        bb2->EnableButton(0,true);
        bb2->EnableButton(1,true);


        if (errMsg) {
            _String         errStr (errMsg);
            sqlite3_free    (errMsg);
            ProblemReport   (errStr, (Ptr)this);
        } else {
            result = true;
        }

        currentTable = tableIndex;
        UpdateStatusBar ();
    }
    return result;
}

//__________________________________________________________
_List*  _HYDBWindow::ExecuteSQLBlurb (_String& theBlurb)
{
    _List      * outList = new _List;
    checkPointer (outList);
    if (dbID >= 0) {
        char  * errMsg = nil;
        if (sqlite3_exec((sqlite3*)sqlDatabases.lData[dbID], theBlurb.sData, _HYDBWCallBack, (Ptr)outList, &errMsg) != SQLITE_OK) {
            _String errStr = _String("SQL Error:") & errMsg;
            ProblemReport (errStr, (Ptr)this);
            sqlite3_free (errMsg);
            outList->Clear();
            return outList;
        }
    }
    return outList;
}



//__________________________________________________________

bool ExecuteSomeCodeAndCheck (_String & code)
{
    _ExecutionList      exl  (code);
    exl.Execute ();
    if (terminateExecution) {
        terminateExecution = false;
        return false;
    }
    return true;
}
//__________________________________________________________

_List*  _HYDBWindow::RetrieveRecord (long row, _String& columns, bool report, _String * newValue)
{
    if (dbID >= 0 && currentTable >= 0) {
        _String    queryString (128L, true),
                   whereString (128L, true);

        _HYTable * tmain = (_HYTable*)GetObject (0),
                   * ttop    = (_HYTable*)GetObject (1);

        queryString << "SELECT ALL ";
        queryString << columns;
        queryString << " FROM ";
        queryString << (_String*)tableList(currentTable);
        whereString << " WHERE (";

        bool        addAnd = false;

        for (long k=0; k<tmain->horizontalSpaces.lLength; k++) {
            _String * cellData = (_String*)tmain->GetCellData (k,row);
            if (cellData->sLength) {
                if (addAnd) {
                    whereString << " AND ";
                }

                whereString << (_String*)ttop->GetCellData (k,0);
                if (cellData->endswith ("...")) {
                    whereString << " LIKE '";
                    whereString << cellData->Cut (0, cellData->sLength-4);
                    whereString << "%";
                } else {
                    whereString << " = '";
                    whereString << cellData;
                }
                whereString << '\'';

                addAnd = true;

            }
        }

        whereString << ")";
        whereString.Finalize();
        queryString << whereString;
        queryString.Finalize();

        _List * ml = ExecuteSQLBlurb (queryString);

        if (ml->lLength) {
            if (((_List*)(*ml)(0))->lLength > 2) {
                if (report) {
                    queryString = "Internal Error: could not extract a record from the database, because multiple records matched selected column values";
                    ProblemReport (queryString,(Ptr)this);
                }
                ml->Clear();
                return ml;
            } else {
                if (newValue) {
                    ml->Clear();

                    _String replaceString (128L,true);

                    replaceString << "UPDATE OR REPLACE  ";
                    replaceString << (_String*)tableList(currentTable);
                    replaceString << " SET ";
                    replaceString << columns;
                    replaceString << " = '";
                    replaceString << newValue;
                    replaceString << "' ";
                    replaceString << whereString;
                    replaceString.Finalize();

                    char  * errMsg = nil;
                    if (sqlite3_exec((sqlite3*)sqlDatabases.lData[dbID], replaceString.sData, NULL, NULL, &errMsg) != SQLITE_OK) {
                        _String errStr = _String("SQL Error:") & errMsg;
                        ProblemReport (errStr, (Ptr)this);
                        sqlite3_free (errMsg);
                        return ml;
                    }
                    (*ml) && & empty;
                    return ml;
                } else {
                    _List * records = new _List;
                    checkPointer (records);

                    for (long k=0; k<ml->lLength; k++) {
                        (*records) << (*((_List*)(*((_List*)(*ml)(k)))(1)))(0);
                    }

                    DeleteObject (ml);
                    return records;
                }
            }
        } else {
            if (report) {
                queryString = "Internal Error: could not extract a record from the database, because no records matched column values";
                ProblemReport (queryString,(Ptr)this);
            }
        }

        return ml;
    }

    return new _List;

}

//__________________________________________________________

void    _HYDBWindow::EditLongEntry (void)
{
    _HYTable * tmain = (_HYTable*)GetObject (0),
               * ttop  = (_HYTable*)GetObject (1);

    _SimpleList tsel;
    tmain->GetSelection (tsel);

    if (tsel.lLength == 1) {
        _List* outRec = RetrieveRecord (tsel.lData[0]/tmain->horizontalSpaces.lLength,
                                        *(_String*) ttop->GetCellData(tsel.lData[0]%tmain->horizontalSpaces.lLength,0),
                                        true
                                       );

        if (outRec->lLength) {
            _String * newValue = new _String (*(_String*)(*outRec)(0));
            checkPointer (newValue);

            if (!EnterStringDialog (*newValue, _HY_DBW_longFieldEdit, (Ptr)this)) {
                DeleteObject (newValue);
                DeleteObject (outRec);
                return;
            }


            _List * repResult = RetrieveRecord (tsel.lData[0]/tmain->horizontalSpaces.lLength,
                                                *(_String*) ttop->GetCellData(tsel.lData[0]%tmain->horizontalSpaces.lLength,0),
                                                true,
                                                newValue
                                               );

            if (repResult->lLength) {
                if (newValue->sLength > lengthLimit) {
                    newValue->Trim(0,lengthLimit-1);
                    *newValue = *newValue & "...";
                }
                tmain->SetCellData (newValue,tsel.lData[0]/tmain->horizontalSpaces.lLength,tsel.lData[0]%tmain->horizontalSpaces.lLength,HY_TABLE_STATIC_TEXT,false);
                tmain->_MarkCellsForUpdate (tsel);
            } else {
                DeleteObject (newValue);
            }


            DeleteObject (repResult);

        }
        DeleteObject (outRec);
    }
}


//__________________________________________________________

void    _HYDBWindow::RunSQLQuery (void)
{
    if (dbID >= 0 && currentTable >=0) {
        _HYTextBox * sqlCommand = (_HYTextBox*)GetObject (8);
        _String      sqlText    (sqlCommand->GetText());

        bool         addToHistory = true;
        long         curSel       = ((_HYPullDown*)GetObject (4))->GetSelection();

        if (curSel == 0) {
            addToHistory = LoadTable (currentTable, &sqlText);
        } else {
            FILE * thisFile = doFileOpen (((_String*)_HYDBW_FilePaths(curSel-2))->sData,"rb");
            _String processorCode (thisFile);
            fclose (thisFile);
            processorCode = _HY_DBW_RunMe & "=1;\n" & _HY_DBW_OutputDBID & "=" & (long)dbID & ";\n" & _HY_DBW_SqlQuery & " = \"" & sqlText.Replace ("\"","\\\"",true) & "\";\n" & processorCode;
            addToHistory = ExecuteSomeCodeAndCheck (processorCode);
            DeleteVariable (_HY_DBW_RunMe,true);
            DeleteVariable (_HY_DBW_OutputDBID,true);
            DeleteVariable (_HY_DBW_SqlQuery,true);
        }

        if (addToHistory) {
            if (sqlHistory.lLength > 100) {
                sqlHistory.Delete (0);
            }

            sqlHistory && & sqlText;
        }
        ((_HYButtonBar*)GetObject(7))->EnableButton (2,sqlHistory.lLength);
        ScanTables ();
    }
}

//__________________________________________________________

bool    _HYDBWindow::ProcessEvent (_HYEvent* e)
{
    _String firstArg,
            secondArg,
            thirdArg;

    long    k,i,f;

    bool    done = false;

    if (e->EventClass()==_hyScrollingEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        long g = e->EventCode().Find(',',f+1);
        if (g>=0) {
            i = MatchComponentID (firstArg);
            _HYTable* theTable = nil;
            if ( i==0 || i==2 ) {
                firstArg = e->EventCode().Cut (g+1,-1);
                k = firstArg.toNum();
                if (k) {
                    theTable    = (_HYTable*) ((i==0)?GetObject(2):
                                               GetObject(0));
                    theTable->SetMessageRecipient (nil);
                    theTable->ProcessEvent (generateScrollEvent (0,k));
                } else if (i==0) {
                    firstArg = e->EventCode().Cut (f+1,g-1);
                    k = firstArg.toNum();
                    if (k) {
                        theTable    = (_HYTable*) GetObject(1);
                        theTable->SetMessageRecipient (nil);
                        theTable->ProcessEvent (generateScrollEvent (k,0));
                    }
                }
                done = true;
            } else if (i==1) {
                firstArg = e->EventCode().Cut (f+1,g-1);
                k = firstArg.toNum();
                if (k) {
                    theTable    = (_HYTable*)     GetObject  (0);
                    theTable->ProcessEvent (generateScrollEvent (k,0));
                    theTable->SetMessageRecipient (nil);
                }
                done = true;
            }
            if (theTable) {
                theTable->SetMessageRecipient (this);
            }
        }
    } else if (e->EventClass()==_hyTableResizeCEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        i = MatchComponentID (firstArg);
        if (i==0 || i==1) {
            firstArg = e->EventCode().Cut (f+1,-1);
            f = firstArg.Find(',');
            k = firstArg.Cut(f+1,-1).toNum(); // shift
            f = firstArg.Cut(0,f-1).toNum();  // column
            _HYTable*     table = (_HYTable*)     GetObject  (1-i);
            table->SetColumnSpacing (f,k,true);
            dim = MinMaxWindowDimensions ();
            done = true;
        }
    } else if (e->EventClass()==_hyTableChangeSelEvent) {
        i = MatchComponentID (e->EventCode());
        done = true;

        _HYTable* tmain = (_HYTable*)GetObject (0),
                  * ttop  = (_HYTable*)GetObject (1),
                    * tleft = (_HYTable*)GetObject (2);

        if (i==1 || i==2) {
            _SimpleList      tSel,
                             tSel2;

            ttop  -> GetSelection (tSel);
            tleft -> GetSelection (tSel2);

            if (tSel.lLength||tSel2.lLength) {
                if (i==1) {
                    tleft->ClearSelection();
                    tmain->SetColumnSelection(tSel);
                } else {
                    ttop->ClearSelection();
                    tmain->SetRowSelection(tSel2);
                }
            }
        }
    } else {
        if (e->EventClass() == _hyTableDblClickEvent) {
            i = MatchComponentID (e->EventCode());
            if (i==0) {
                done = true;
                EditLongEntry ();
            }
        } else if (e->EventClass()==_hyButtonPushEvent) {
            firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
            i = MatchComponentID (firstArg);
            if (i==7) {
                firstArg        = e->EventCode().Cut (f+1,-1);
                k               = firstArg.toNum();
                switch  (k) {
                case 0:
                    RunSQLQuery ();
                    break;
                case 1:
                case 2: {
                    int h,v;

                    _HYButtonBar * bb2 = (_HYButtonBar*)GetObject(7);

                    _HYTable     * tmain = (_HYTable*)GetObject (0),
                                   * ttop  = (_HYTable*)GetObject (1);

                    bb2->GetButtonLoc(k,h,v,true);
                    _List  menuOptions;

                    if (k==1)
                        for (f=0; f<_HYDBW_Templates.lLength; f++) {
                            if (((_String*)(*((_List*)_HYDBW_Templates(f)))(1))->sLength) {
                                menuOptions << (*((_List*)_HYDBW_Templates(f)))(0);
                            } else {
                                menuOptions << & menuSeparator;
                            }
                        }
                    else
                        for (f=0; f<sqlHistory.lLength; f++) {
                            _String option ((long)(1+f)),
                                    *command = (_String*)sqlHistory(f);

                            option = option & ". ";
                            i = command->Find (' ',command->FirstNonSpace(0,-1,1),-1);
                            long f2 = command->Find ('(');
                            if (i<0) {
                                i = 50;
                            }
                            if (f2<0) {
                                f2 = 50;
                            }
                            i = MIN(50,MIN(i,f2))-1;
                            if (i<command->sLength) {
                                option = option & command->Cut (0,i) & "...";
                            } else {
                                option = option & *command;
                            }

                            menuOptions && &option;
                        }

                    firstArg = HandlePullDown (menuOptions,h,v,0);
                    bb2->_UnpushButton();

                    if (firstArg.sLength) {
                        f = menuOptions.Find (&firstArg);
                        if (f>=0) {
                            if (k==1) {
                                _String command (*(_String*)(*((_List*)_HYDBW_Templates(f)))(1));
                                if (dbID >= 0 && currentTable >= 0) {
                                    command = command.Replace (_HY_DBW_TablePlaceholder,*(_String*)tableList(currentTable),true);

                                    _SimpleList colSelection;

                                    tmain->ScanRowSelection (colSelection);
                                    if (colSelection.lLength == 1) {
                                        colSelection.Clear();
                                        tmain->GetSelection (colSelection);

                                        _String rowSel (128L,true);
                                        for (long i2=0; i2<colSelection.lLength; i2++) {
                                            long i2c = colSelection.lData[i2]%tmain->horizontalSpaces.lLength,
                                                 i2r = colSelection.lData[i2]/tmain->horizontalSpaces.lLength;

                                            if (i2) {
                                                rowSel << " AND ";
                                            }
                                            rowSel << (_String*)ttop->GetCellData (i2c,0);
                                            rowSel << " = '";
                                            rowSel << (_String*)tmain->GetCellData (i2c,i2r);
                                            rowSel << "' ";
                                        }
                                        rowSel.Finalize();
                                        command = command.Replace (_HY_DBW_RowValues,rowSel,true);
                                        colSelection.Clear();
                                        tmain->ScanColumnSelection (colSelection);
                                    } else {
                                        colSelection.Clear();
                                        ttop->GetColumnSelection (colSelection);
                                    }

                                    if (colSelection.lLength) {
                                        _String colNames (128L,true);
                                        for (long i2=0; i2<colSelection.lLength; i2++) {
                                            if (i2) {
                                                colNames << ',';
                                            }
                                            colNames << (_String*)ttop->GetCellData (colSelection.lData[i2],0);
                                        }
                                        colNames.Finalize();
                                        command = command.Replace (_HY_DBW_ColumnPlaceholder,colNames,true);
                                    }
                                }

                                ((_HYTextBox*)GetObject(8))->SetText (command,true);
                            } else {
                                ((_HYTextBox*)GetObject(8))->SetText (*(_String*)sqlHistory(f),true);
                            }
                        }
                    }
                }

                }
                done = true;
            }
        } else if (e->EventClass() == _hyMenuSelChangeEvent) {
            firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
            i =  MatchComponentID (firstArg);
            firstArg = e->EventCode().Cut (f+1,-1);
            k = firstArg.toNum();
            if (i == 3) {
                LoadTable (k);
            }
        }
    }

    if (done) {
        DeleteObject (e);
        return true;
    }
    return _HYTWindow::ProcessEvent(e);
}

//__________________________________________________________

void     ReadSQLPlugins (void)
{
    _String     pathToModelTemplates;
    _List       receptacle;

    pathToModelTemplates = baseDirectory&"ChartAddIns"&baseDirectory.sData[baseDirectory.sLength-1]&"DBAddIns";
    ScanDirectoryForFileNames (pathToModelTemplates,receptacle,false);

    for (long k=0; k<receptacle.lLength; k++) {
        FILE * thisFile = doFileOpen (((_String*)receptacle(k))->sData,"rb");
        if (thisFile) {
            _String buffer (thisFile);
            fclose (thisFile);
            if (buffer.sLength) {
                _ExecutionList ex (buffer);
                long f = ex.ExecuteAndClean (batchLanguageFunctionNames.lLength,&_HY_DBW_OutputProcessorFunction);

                if ((((_String*)receptacle(k)))->endswith(_HY_DBW_TemplateFileName)) {
                    _Variable * opList = CheckReceptacle (&_HY_DBW_TemplateFileName,empty,false);
                    if (opList && opList->ObjectClass () == ASSOCIATIVE_LIST) {
                        _AssociativeList * aList = (_AssociativeList*)opList->GetValue ();
                        _List * meKeys = aList->GetKeys ();

                        for (long k=0; k<meKeys->lLength; k++) {
                            _String*  aKey   = (_String*)(*meKeys)(k);
                            _PMathObj aValue = aList->GetByKey (*aKey,STRING);
                            if (aValue) {
                                _List * aPair = new _List;
                                checkPointer (aPair);
                                (*aPair) << aKey;
                                (*aPair) << ((_FString*)aValue)->theString;
                                _HYDBW_Templates << aPair;
                            }
                        }

                    }
                } else {
                    _Variable * postProcName = CheckReceptacle (&_HY_DBW_OutputProcessorName,empty,false);
                    if (f >= 0 && postProcName && postProcName->ObjectClass () == STRING) {
                        _HYDBW_FilePaths << receptacle (k);
                        _HYDBW_OutputProcessors << ((_FString*)postProcName->GetValue())->theString;
                    }
                    DeleteVariable (_HY_DBW_OutputProcessorName,true);
                }

            }
        }
    }
}


// EOF