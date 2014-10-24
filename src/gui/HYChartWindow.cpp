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

/*
    Chart Window

    Sergei L. Kosakovsky Pond, December 2001
    Rev                        April 2002,
                                 Added contrast plot
                                 Added 3D Plot
                                 Added Legends and Labels
                               January 2003
                                 Added more options to saving
                                 Added 3D scaled plots
*/


#include "HYChartWindow.h"
#include "HYLabel.h"
#include "HYButtonBar.h"
#include "HYTextBox.h"
#include "HYUtils.h"
#include "HYEventTypes.h"
#include "HYCanvas.h"
#include "HYDialogs.h"
#include "parser.h"
#include "likefunc.h"
#include "batchlan.h"
#include "math.h"
#include "HYButton.h"

#include "HYDataPanel.h"

#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

#define     HY_CHART_COLOR_AXIS_LABEL           0
#define     HY_CHART_COLOR_COORD_LABEL          1
#define     HY_CHART_COLOR_BACKGROUND           2
#define     HY_CHART_COLOR_VALUELINE            3
#define     HY_CHART_COLOR_BARBORDER            4
#define     HY_CHART_COLOR_LEGEND               5
#define     HY_CHART_COLOR_OVERPLOT             6

#define     HY_CHART_COLOR_3D_FLOOR             7
#define     HY_CHART_COLOR_3D_LEFT              8
#define     HY_CHART_COLOR_3D_BACK              9
#define     HY_CHART_COLOR_3D_VALUELINE         10
#define     HY_CHART_COLOR_3D_HISTOGRAM         11
#define     HY_CHART_COLOR_WINDOW_BACKGROUND    12

#define     HY_CHART_COLOR_OFFSET               13

extern      _SimpleList                         windowObjects;

extern      _HYColor                            tableDefaultBk2,
            labelBColor,
            labelFColor,
            black;

extern      _String                             donotWarnAgain,
            windowTypeTable,
            windowTypeDistribTable;

_HYColor            chartColors [HY_CHART_COLOR_COUNT] = {
    {255*.94, 255*.12, 255*.11 },//(Red)
    {255*.41, 255*.46, 255*.91 },//(Evening Blue)
    {255    , 255*.91, 255*.34 },//(Banana)
    {255*.18, 255*.55, 255*.13 },//(Clover)
    {255*.55, 255*.38, 255*.21 },//(Dirt)
    {255*.42, 255*.09, 255*.69 },//(Royal Violet)
    {255*.09, 255*.29, 255*.51 },//(Sea Blue)
    {255   ,  255*.57, 255*.09 },//(Orange)
    {255*.67, 255*.67, 255*.67 },//(Concrete)
    {255*.85, 255*.27, 255*.42 } //(Carnation)
};

extern      _Parameter  pi_const;

_List       chartProcessors,
            distribProcessors;

_String     chartProcessorDataMatrix ("SELECTED_CHART_DATA"),
            chartProcessorRowMatrix  ("SELECTED_CHART_ROWS"),
            chartProcessorColMatrix  ("SELECTED_CHART_COLS"),
            updateCellData           ("UPDATE_CELL_DATA"),
            updateTableSelection     ("UPDATE_CHART_SELECTION"),
            tableDimensions          ("CHART_DIMENSIONS"),
            chartColumnHeaders       ("CHART_COLUMN_HEADERS"),
            nonEmptyChartSelection   ("NON_EMPTY_SELECTION"),
            chartColorChange         ("Change chart color"),
            chartCustomColumnColors  ("_CHART_CUSTOM_COLUMN_COLORS_"),

chartColorLabels        [HY_CHART_COLOR_OFFSET] = {
    "Axis label",
    "Coordinate labels",
    "Chart background",
    "Value line",
    "Bar border",
    "Legend text and border",
    "Overplot",
    "3D Floor",
    "3D Left Wall",
    "3D Back Wall",
    "3D Value line",
    "3D Histogram",
    "Window Background"
};

#define     HY_LARGE_COMPONENT_SIZE  10000L



//__________________________________________________________
_HYChartWindow::_HYChartWindow (_String name, _List& columns, _Matrix& data, _HYWindow*):_HYTWindow (name, true)
{

    long            hDim = data.GetHDim();

    _HYRect         canvasSettings =
    {25,40,25,40,HY_COMPONENT_NO_SCROLL|HY_COMPONENT_BORDER_B};

    _HYLabel*       l1      = new _HYLabel (canvasSettings, GetOSWindowData());

    canvasSettings.left     = 100;
    canvasSettings.right    = 100;
    _HYPullDown*    p1      = new _HYPullDown (canvasSettings,GetOSWindowData());

    canvasSettings.left     = 20;
    canvasSettings.right    = 20;
    _HYLabel*       l2      = new _HYLabel (canvasSettings,GetOSWindowData());
    _HYLabel*       l3      = new _HYLabel (canvasSettings,GetOSWindowData());

    canvasSettings.left     = 100;
    canvasSettings.right    = 100;
    _HYPullDown*    p2      = new _HYPullDown (canvasSettings,GetOSWindowData());
    canvasSettings.right    = HY_LARGE_COMPONENT_SIZE;
    _HYPullDown*    p3      = new _HYPullDown (canvasSettings,GetOSWindowData());

    canvasSettings.top      = HY_CHART_TABLE_HEIGHT+HY_SCROLLER_WIDTH;
    canvasSettings.bottom   = HY_CHART_TABLE_HEIGHT+HY_SCROLLER_WIDTH;

    canvasSettings.width    = HY_COMPONENT_V_SCROLL|HY_COMPONENT_H_SCROLL;
    _HYTable*   table       = new _HYTable (canvasSettings,GetOSWindowData(),1,1,100,16,HY_TABLE_STATIC_TEXT);

    canvasSettings.top      = canvasSettings.bottom     = 18;
    canvasSettings.width    = HY_COMPONENT_NO_SCROLL | HY_COMPONENT_BORDER_T | HY_COMPONENT_BORDER_B;

    _HYTable*   tableHead   = new _HYTable (canvasSettings,GetOSWindowData(),1,1,100,18,HY_TABLE_STATIC_TEXT);
    canvasSettings.left     = canvasSettings.right = 40;
    canvasSettings.top      = HY_CHART_TABLE_HEIGHT+HY_SCROLLER_WIDTH;
    canvasSettings.bottom   = HY_CHART_TABLE_HEIGHT+HY_SCROLLER_WIDTH;

    canvasSettings.width    = HY_COMPONENT_NO_SCROLL|HY_COMPONENT_BORDER_R;
    _HYTable*   tableLeft   = new _HYTable (canvasSettings,GetOSWindowData(),1,1,40,16,HY_TABLE_STATIC_TEXT);

    table->selectionType    = HY_TABLE_DONT_SIZE|HY_TABLE_DONT_GROW_VERT|HY_TABLE_FOCUSABLE|HY_TABLE_HORZ_STRETCH;
    tableHead->selectionType= HY_TABLE_DONT_GROW_VERT|HY_TABLE_HORZ_STRETCH;
    tableLeft->selectionType= HY_TABLE_DONT_SIZE|HY_TABLE_DONT_GROW_VERT;

    canvasSettings.left     = canvasSettings.right  = 40;
    canvasSettings.top      = canvasSettings.bottom = 18;

    canvasSettings.width    = HY_COMPONENT_NO_SCROLL | HY_COMPONENT_BORDER_T | HY_COMPONENT_BORDER_R | HY_COMPONENT_BORDER_B;
    _HYLabel*   l4          = new _HYLabel (canvasSettings, GetOSWindowData());

    canvasSettings.left     = 380;
    canvasSettings.right    = HY_LARGE_COMPONENT_SIZE;
    canvasSettings.bottom   = HY_LARGE_COMPONENT_SIZE;
    canvasSettings.top      = HY_CHART_TABLE_HEIGHT;

    canvasSettings.width    = HY_COMPONENT_NO_SCROLL;

    _HYStretchCanvas* cn    = new _HYStretchCanvas (canvasSettings,GetOSWindowData(),200,380,32,HY_SCANVAS_HORIZONTAL|HY_SCANVAS_VERTICAL);

    table->SetMessageRecipient (this);
    tableHead->SetMessageRecipient (this);
    tableLeft->SetMessageRecipient (this);
    p1->SetMessageRecipient (this);
    p2->SetMessageRecipient (this);
    p3->SetMessageRecipient (this);
    cn->SetMessageRecipient (this);

    AddObject (cn);        // 0
    AddObject (table);     // 1
    AddObject (tableHead); // 2
    AddObject (tableLeft); // 3
    AddObject (p1);        // 4
    AddObject (p2);        // 5
    AddObject (p3);        // 6
    AddObject (l1);        // 7
    AddObject (l2);        // 8
    AddObject (l3);        // 9
    AddObject (l4);        // 10

    SetTableDimensions (4,6);

    SetCell   (HY_CHART_MENU_ROW,0,l1);
    SetCell   (HY_CHART_MENU_ROW,1,p1);
    SetCell   (HY_CHART_MENU_ROW,2,l2);
    SetCell   (HY_CHART_MENU_ROW,3,p2);
    SetCell   (HY_CHART_MENU_ROW,4,l3);
    SetCell   (HY_CHART_MENU_ROW,5,p3);

    SetCell   (HY_CHART_CHART_ROW,0,cn);
    SetCell   (HY_CHART_CHART_ROW,1,cn);
    SetCell   (HY_CHART_CHART_ROW,2,cn);
    SetCell   (HY_CHART_CHART_ROW,3,cn);
    SetCell   (HY_CHART_CHART_ROW,4,cn);
    SetCell   (HY_CHART_CHART_ROW,5,cn);

    SetCell   (HY_CHART_HEADER_ROW,0,l4);
    SetCell   (HY_CHART_HEADER_ROW,1,tableHead);
    SetCell   (HY_CHART_HEADER_ROW,2,tableHead);
    SetCell   (HY_CHART_HEADER_ROW,3,tableHead);
    SetCell   (HY_CHART_HEADER_ROW,4,tableHead);
    SetCell   (HY_CHART_HEADER_ROW,5,tableHead);

    SetCell   (HY_CHART_TABLE_ROW,0,tableLeft);
    SetCell   (HY_CHART_TABLE_ROW,1,table);
    SetCell   (HY_CHART_TABLE_ROW,2,table);
    SetCell   (HY_CHART_TABLE_ROW,3,table);
    SetCell   (HY_CHART_TABLE_ROW,4,table);
    SetCell   (HY_CHART_TABLE_ROW,5,table);


    p1->SetBackColor (labelBColor);
    p2->SetBackColor (labelBColor);
    p3->SetBackColor (labelBColor);

    l1->SetForeColor (labelFColor);
    l2->SetForeColor (labelFColor);
    l3->SetForeColor (labelFColor);

    l1->SetBackColor (labelBColor);
    l2->SetBackColor (labelBColor);
    l3->SetBackColor (labelBColor);

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
    p3->SetAlignFlags (HY_ALIGN_LEFT);

    p1->AddMenuItem ("None",-1);
    p1->AddMenuItem ("Bar Chart",-1);
    p1->AddMenuItem ("Pie Chart",-1);
    p1->AddMenuItem ("Scatterplot",-1);
    p1->AddMenuItem ("Line Plot",-1);
    p1->AddMenuItem ("Step Plot",-1);
    p1->AddMenuItem ("Stacked Bars",-1);
    p1->AddMenuItem ("Contrast Bars",-1);
    p1->AddMenuItem ("Error Bars",-1);
    p1->AddMenuItem ("3D Bar Chart",-1);
    p1->AddMenuItem ("3D Histogram",-1);
    p1->AddMenuItem ("3D Scatterplot",-1);


    p2->EnableMenu  (false);
    p3->EnableMenu  (false);

    l1->SetAlignFlags (HY_ALIGN_LEFT);
    l2->SetAlignFlags (HY_ALIGN_LEFT);
    l3->SetAlignFlags (HY_ALIGN_LEFT);

    l1->SetShadow (true);
    l2->SetShadow (true);
    l3->SetShadow (true);

    l1->SetFont (labelFont);
    l2->SetFont (labelFont);
    l3->SetFont (labelFont);

    l1->SetText ("Type:");
    l2->SetText ("X:");
    l3->SetText ("Y:");

    l4->SetBackColor (tableDefaultBk2);

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

    // define default colors

    theColors << HYColorToLong (black);
    theColors << HYColorToLong (black);
    theColors << HYColorToLong ((_HYColor) {
        210,210,210
    });
    theColors << HYColorToLong ((_HYColor) {
        255,255,255
    });
    theColors << HYColorToLong (black);
    theColors << HYColorToLong (black);
    theColors << HYColorToLong ((_HYColor) {
        100,100,100
    });

    theColors << HYColorToLong ((_HYColor) {
        180,180,180
    });
    theColors << HYColorToLong ((_HYColor) {
        200,200,200
    });
    theColors << HYColorToLong ((_HYColor) {
        220,220,220
    });
    theColors << HYColorToLong (black);
    theColors << HYColorToLong ((_HYColor) {
        60,60,60
    });
    theColors << HYColorToLong ((_HYColor) {
        255,255,255
    });

    for (hDim = 0; hDim < HY_CHART_COLOR_COUNT; hDim ++) {
        theColors << HYColorToLong (chartColors[hDim]);
    }

    SetTable           (columns,data);
    SetPosition        (70,70);
    SetWindowRectangle (0,0,400,400);

    userMin = 0.;
    userMax = 0.;

    projectionMatrix  = nil;

    showLegend = HY_CHART_LEGEND_NONE;

    xyAngle = pi_const/2.4;
    zAngle  = pi_const/4;
    oR = 10;
    ComputeProjectionSettings ();
    xScale = 1.;
    yScale = 1.;
    xShift = 0.;
    yShift = 0.;

    xAxis3DScale = -1;
    yAxis3DScale = -1;
    surfaceDivs  = 16;

    projectionMatrix = ComputeProjectionMatrix();
    suspendDraw      = false;


    if (chartProcessors.lLength==0) {
        ReadChartProcessors ();
    }

    DeleteObject (table);
    DeleteObject (tableHead);
    DeleteObject (tableLeft);
    DeleteObject (l1);
    DeleteObject (l2);
    DeleteObject (l3);
    DeleteObject (l4);
    DeleteObject (p1);
    DeleteObject (p2);
    DeleteObject (p3);
    DeleteObject (cn);
}

//__________________________________________________________
_HYChartWindow::~_HYChartWindow()
{
    if (projectionMatrix) {
        DeleteObject (projectionMatrix);
    }
}

//__________________________________________________________

void    _HYChartWindow::ComputeProjectionSettings (void)
{
    projectionSettings[0] = oR*sin(zAngle)*cos(xyAngle);
    projectionSettings[1] = oR*sin(zAngle)*sin(xyAngle);
    projectionSettings[2] = -oR*cos(zAngle);
    projectionSettings[3] = -projectionSettings[0]/oR;
    projectionSettings[4] = -projectionSettings[1]/oR;
    projectionSettings[5] = -projectionSettings[2]/oR;
    projectionSettings[6] = oR;
}


//__________________________________________________________

void    _HYChartWindow::SetColumnHeaders (_List& h)
{
    _HYTable* th = (_HYTable*)GetObject (2),
              * t  = (_HYTable*)GetObject (1);

    if (h.lLength >= th->horizontalSpaces.lLength) {
        for (long k=0; k<th->horizontalSpaces.lLength; k++) {
            th->SetCellData (h(k),0,k,HY_TABLE_EDIT_TEXT|HY_TABLE_BEVELED|HY_TABLE_BOLD,true);
        }

        th->AutoFitWidth (*t);

        if (h.lLength ==  th->horizontalSpaces.lLength+1) {
            _String * rowHeaders = (_String*)h(th->horizontalSpaces.lLength);
            _List   * split      = rowHeaders->Tokenize (";");
            SetRowHeaders   (*split);
            DeleteObject (split);
        }
    }

    t->_MarkForUpdate();
    th->_MarkForUpdate();
    SetPullDowns (h);
}

//__________________________________________________________

void    _HYChartWindow::SetDataMatrix (_Matrix& h)
{
    _HYTable* th = (_HYTable*)GetObject (2),
              * t  = (_HYTable*)GetObject (1);

    if ((h.GetVDim() >= th->horizontalSpaces.lLength)&&(h.GetHDim() >= t->verticalSpaces.lLength)) {
        _Matrix * numberMatrix = (_Matrix*)h.ComputeNumeric();
        for (long k=0; k<th->horizontalSpaces.lLength; k++)
            for (long j=0; j<t->verticalSpaces.lLength; j++) {
                _String cellData ((*numberMatrix)(j,k));
                th->SetCellData (&cellData,j,k,HY_TABLE_EDIT_TEXT,true);
            }
        t->AutoFitWidth (*th);
    }
    t->_MarkForUpdate();
    th->_MarkForUpdate();
}

//__________________________________________________________

bool    _HYChartWindow::SetUserBounds (_Parameter minY, _Parameter maxY, bool redraw)
{
    bool changed = false;
    if (!CheckEqual (minY, userMin)) {
        changed = true;
        userMin = minY;
    }
    if (!CheckEqual (maxY, userMax)) {
        changed = true;
        userMax = maxY;
    }

    if (!CheckEqual (userMin, userMax)) {
        if (userMax < userMin) {
            changed = true;
            userMax = userMin + 1.;
        }
    }

    if (changed && redraw) {
        DrawChart ();
    }

    return changed;
}

//__________________________________________________________

void    _HYChartWindow::SetRowHeaders (_List& l)
{
    _HYTable* tl = (_HYTable*)GetObject (3);

    _String   space (" ");

    if (l.lLength == tl->verticalSpaces.lLength) {
        for (long k=1; k<l.lLength; k++) {
            _String* rowHeader = (_String*)l(k);
            if (rowHeader->sLength==0) {
                rowHeader = &space;
            }

            tl->SetCellData (rowHeader,k-1,0,tl->cellTypes[k-1],true);
        }

        _HYLabel* l4 = (_HYLabel*)GetObject   (10);
        l4->SetAlignFlags (HY_ALIGN_LEFT);
        _HYFont hFont = tl->textFont;
        hFont.style = HY_FONT_BOLD;
        l4->SetFont (hFont);
        _String* rowHeader = (_String*)l(0);
        l4->SetText (*rowHeader);
#ifdef __WINDOZE__
        l4->SetAlignFlags (HY_ALIGN_BOTTOM);
#endif
        // now compute new dimensions

        _HYRect     labelSize = l4->_SuggestDimensions();
        tl->AutoFitColumn (0, false, false);

        long        newWidth = tl->GetColumnSpacing(0);

        if (labelSize.right > newWidth) {
            tl->SetColumnSpacing (0,labelSize.right-newWidth,false);
            newWidth = labelSize.right;
        }

        l4->settings.left = l4->settings.right = newWidth;
        tl->settings.left = tl->settings.right = newWidth;

        _HYRect         windowRect = GetWindowRect();

        SetWindowRectangle (windowRect.top, windowRect.left, windowRect.bottom, windowRect.right);
    }
}

//__________________________________________________________

void    _HYChartWindow::SetTable (_List& l, _Matrix& h)
{
    _HYTable* th = (_HYTable*)GetObject (2),
              * t  = (_HYTable*)GetObject (1),
                * tl = (_HYTable*)GetObject (3);

    SetUserBounds ();

    if (((h.GetVDim() == l.lLength)||(h.GetVDim()+1 == l.lLength))&&(h.GetHDim() >= 1)&&(l.lLength>=1)) {
        long      cC = h.GetVDim();
        _Matrix * n = (_Matrix*)h.ComputeNumeric();

        if (terminateExecution) {
            terminateExecution = false;
            _String errMsg ("Matrix evaluation failed in 'SetTable'.");
            ProblemReport (errMsg);
            return;
        }

        th->ClearTable(true);
        tl->ClearTable();
        t->ClearTable (true);

        long k;

        for (k=0; k<cC; k++) {
            th->AddColumn   (-1,40,HY_TABLE_STATIC_TEXT);
            t->AddColumn   (-1,40,HY_TABLE_STATIC_TEXT);
        }

        th->AddRow (-1,18,HY_TABLE_STATIC_TEXT);
        for (k=0; k<cC; k++) {
            th->SetCellData (l(k),0,k,HY_TABLE_EDIT_TEXT|HY_TABLE_BEVELED|HY_TABLE_BOLD,true);
        }

        th->AutoFitWidth();
        //th->SetColumnSpacing (cC-1,HY_LARGE_COMPONENT_SIZE-th->GetColumnSpacing (cC-1),false);

        for (k=0; k<h.GetHDim(); k++) {
            _String cellData (k+1);
            tl->AddRow (-1,16,HY_TABLE_STATIC_TEXT);
            tl->SetCellData (&cellData,k,0,HY_TABLE_STATIC_TEXT|HY_TABLE_BEVELED|HY_TABLE_BOLD,true);
            /*t->AddRow (-1,16,HY_TABLE_STATIC_TEXT);
            for (j=0; j<cC; j++)
            {
                cellData = (*n)(k,j);
                t->SetCellData (&cellData,k,j,HY_TABLE_EDIT_TEXT,true);
            }*/
        }
        t->SetTableFromMx (n,16,40,HY_TABLE_EDIT_TEXT);
        t->AutoFitWidth (*th);

        //t->SetColumnSpacing (cC-1,HY_LARGE_COMPONENT_SIZE-t->GetColumnSpacing(cC-1),false);
        //th->SetColumnSpacing (cC-1,HY_LARGE_COMPONENT_SIZE-th->GetColumnSpacing(cC-1),false);

        t->EnforceWidth  (componentR.lData[1]-componentL.lData[1]+1-HY_SCROLLER_WIDTH, 0x7fffffff,true);
        th->EnforceWidth (componentR.lData[2]-componentL.lData[2]+1, 0x7fffffff,true);

        t->EnforceHeight(HY_CHART_TABLE_HEIGHT,0x7fffffff);
        tl->EnforceHeight(HY_CHART_TABLE_HEIGHT,0x7fffffff);
        tl->AddRow (-1,HY_SCROLLER_WIDTH,HY_TABLE_STATIC_TEXT);
        t->_MarkForUpdate();
        tl->_MarkForUpdate();
        th->_MarkForUpdate();
        SetPullDowns (l);

        if (l.lLength > cC) {
            _String * rowHeaders = (_String*)l(cC);
            _List   * split      = rowHeaders->Tokenize (";");
            SetRowHeaders   (*split);
            DeleteObject (split);
        }
    }
}


//__________________________________________________________

void    _HYChartWindow::SetPullDowns (_List& l)
{
    _HYTable* t  = (_HYTable*)GetObject (1);

    _HYPullDown* p2 = (_HYPullDown*)GetObject (5),
                 * p3 = (_HYPullDown*)GetObject (6);

    p2->DeleteAllItems();
    p3->DeleteAllItems();

    p2->AddMenuItem     ("Index",-1);
    p3->AddMenuItem     ("None",-1);

    long        upto = t->horizontalSpaces.lLength;

    for (long k=0; k<upto; k++) {
        p2->AddMenuItem (*((_String*)l(k)),-1);
        p3->AddMenuItem (*((_String*)l(k)),-1);
    }
}
//__________________________________________________________

void    _HYChartWindow::SetYPullDown (void)
{
    _HYPullDown* p3 = (_HYPullDown*)GetObject (6);

    long f = p3->MenuItemCount(),
         k;

    for (k=1; k<f; k++) {
        p3->_SetMenuItemTextStyle (k,HY_FONT_PLAIN);
        _String itemText = * p3->GetMenuItem (k);
        long g = itemText.FindBackwards (" [",0,-1);
        if (g>=0) {
            itemText.Trim (0,g-1);
            p3->SetMenuItem (itemText,k);
        }
    }

    if (ySeries.lLength) {
        for (k=0; k<ySeries.lLength; k++) {
            p3->_SetMenuItemTextStyle (ySeries.lData[k]+1,HY_FONT_ITALIC);
            _String itemText = * p3->GetMenuItem (ySeries.lData[k]+1);

            itemText = itemText & " [" & (long)(k+1) & ']';
            p3->SetMenuItem (itemText,ySeries.lData[k]+1);
        }
        p3->ChangeSelection (ySeries.lData[0]+1,false);
    } else {
        p3->ChangeSelection (0,false);
    }
}

//__________________________________________________________

void    _HYChartWindow::SetChartType (_String t1, _String t2, _String t3, bool doDraw)
{
    _HYPullDown* p1 = (_HYPullDown*)GetObject (4);
    _HYPullDown* p2 = (_HYPullDown*)GetObject (5);
    _HYPullDown* p3 = (_HYPullDown*)GetObject (6);

    _List   *ySer = t3.Tokenize(";");

    _String* ySerV = (_String*)(*ySer)(0);

    long     f3 = p3->FindMenuItem (*ySerV),
             f2 = p2->FindMenuItem (t2),
             f1 = p1->FindMenuItem (t1);

    suspendDraw = true;
    if ((f1>=0)&&(f2>=0)&&(f3>=0)) {
        ySeries.Clear();
        p3->ChangeSelection (f3,true);
        p2->ChangeSelection (f2,true);
        p1->ChangeSelection (f1,true);
    }

    for (long k=1; k<ySer->lLength; k++) {
        ySerV = (_String*)(*ySer)(k);
        f3 = p3->FindMenuItem (*ySerV);
        if (f3>=0) {
            p3->ChangeSelection (f3,true);
        }
    }

    DeleteObject (ySer);
    suspendDraw = false;
    if (doDraw) {
        DrawChart    ();
    }
}


//__________________________________________________________

bool    _HYChartWindow::ProcessEvent (_HYEvent* e)
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
            k = firstArg.toNum();
            for (i=0; i<components.lLength; i++) {
                if (((_HYGuiObject*)components(i))->MatchID(k)) {
                    break;
                }
            }
            _HYTable* theTable = nil;
            if ((i==1)||(i==3)) {
                firstArg = e->EventCode().Cut (g+1,-1);
                k = firstArg.toNum();
                if (k) {
                    theTable    = (_HYTable*) ((i==1)?GetObject(3):
                                               GetObject(1));
                    theTable->SetMessageRecipient (nil);
                    theTable->ProcessEvent (generateScrollEvent (0,k));
                } else if (i==1) {
                    firstArg = e->EventCode().Cut (f+1,g-1);
                    k = firstArg.toNum();
                    if (k) {
                        theTable    = (_HYTable*) GetObject(2);
                        theTable->SetMessageRecipient (nil);
                        theTable->ProcessEvent (generateScrollEvent (k,0));
                    }
                }
                done = true;
            } else if (i==2) {
                firstArg = e->EventCode().Cut (f+1,g-1);
                k = firstArg.toNum();
                if (k) {
                    theTable    = (_HYTable*)     GetObject  (1);
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
        k = firstArg.toNum();
        for (i=0; i<components.lLength; i++) {
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }
        }
        if (i==2) {
            firstArg = e->EventCode().Cut (f+1,-1);
            f = firstArg.Find(',');
            k = firstArg.Cut(f+1,-1).toNum(); // shift
            f = firstArg.Cut(0,f-1).toNum();  // column
            _HYTable*     table = (_HYTable*)     GetObject  (1);
            table->SetColumnSpacing (f,k,true);
            dim = MinMaxWindowDimensions ();
            done = true;
        }
    } else if (e->EventClass()==_hyTableChangeSelEvent) {
        k = e->EventCode().toNum();
        for (i=0; i<components.lLength; i++) {
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }
        }

        done = true;

        _HYTable* t1 = (_HYTable*)GetObject (1),
                  * t2 = (_HYTable*)GetObject (3),
                    * t3 = (_HYTable*)GetObject (2);

        if ((i==2)||(i==3)) {
            _SimpleList tSel, tSel2;
            t2->GetSelection (tSel);
            t3->GetSelection (tSel2);
            if (tSel.lLength||tSel2.lLength) {
                if (i==3) {
                    t3->ClearSelection();
                    t2->SetMessageRecipient (nil);
                    t1->SetRowSelection(tSel);
                    t2->SetMessageRecipient (this);
                } else {
                    t2->ClearSelection();
                    t2->SetMessageRecipient (nil);
                    t1->SetColumnSelection(tSel2);
                    t2->SetMessageRecipient (this);
                }
            }
        } else {
            if (i==1) {
                if (t2->messageRecipient) {
                    t2->ClearSelection();
                    t3->ClearSelection();
                }
                _SimpleList newSelection;
                t1->GetSelection (newSelection);
            }
        }
    } else if (e->EventClass()==_hyMenuSelChangeEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        for (i=0; i<components.lLength; i++) {
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }
        }
        k = e->EventCode().Cut (f+1,-1).toNum();
        if (i==4) {
            HandleChartTypeChange (k);
            done = true;
        } else {
            if (i==6) {
                _HYPullDown* p3 = (_HYPullDown*)GetObject (6);
                f = p3->MenuItemCount();
                if (k) {
                    i = ySeries.Find (k-1);
                    if (i<0) {
                        ySeries << k-1;
                    } else {
                        ySeries.Delete (i);
                    }
                } else {
                    ySeries.Clear();
                }
                SetYPullDown();
            }
            DrawChart ();
            done = true;
        }
    } else if (e->EventClass()==_hyRebuildSCanvasEvent) {
        k = e->EventCode().toNum();
        for (i=0; i<components.lLength; i++) {
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }
        }
        if (i==0) {
            DrawChart();
        }

        done = true;
    } else if (e->EventClass()==_hyTableEditCellEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        for (i=0; i<components.lLength; i++) {
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }
        }
        if (i==1) {
            firstArg        = e->EventCode().Cut (f+1,-1);
            HandleCellEditEvent (firstArg.toNum());
            done = true;
        } else if (i==2) {
            // update the pulldown menu as well
            _HYTable    *header = (_HYTable*)GetObject    (2);
            _HYPullDown *p2     = (_HYPullDown*)GetObject (5);
            _HYPullDown *p3     = (_HYPullDown*)GetObject (6);

            k = e->EventCode().Cut (f+1,-1).toNum();

            _String*    newName  = (_String*)header->GetCellData(k,0),
                        menuItem = *p3->GetMenuItem(k+1),
                        appendix;

            if (menuItem.sLength>3) {
                if (menuItem.sData[menuItem.sLength-1] == ']') {
                    f = menuItem.FindBackwards ("[",0,menuItem.sLength-2);
                    if (f>=0) {
                        appendix = _String(" ") & menuItem.Cut (f,-1);
                    }
                }
            }

            p2->SetMenuItem (*newName,k+1);
            p2->_MarkForUpdate();
            menuItem = *newName & appendix;
            p3->SetMenuItem (menuItem,k+1);
            p3->_MarkForUpdate();
            if (showLegend!=HY_CHART_LEGEND_NONE) {
                DrawChart();
            }

            done = true;
        }
    }

    /*else
        if (e->EventClass()==_hyMenuOpenEvent)
        {
            k = e->EventCode().toNum();
            for (i=0;i<components.lLength;i++)
            {
                if (((_HYGuiObject*)components(i))->MatchID(k))
                    break;
            }
            if (i==3)
            {
                PrepareLFMenu();
                done = true;
            }
        }
        */
    if (done) {
        DeleteObject (e);
        return true;
    }
    return _HYTWindow::ProcessEvent(e);
}

//__________________________________________________________

void    _HYChartWindow::HandleChartTypeChange (long type)
{
    _HYPullDown         *p2 = (_HYPullDown*)GetObject (5),
                         *p3 = (_HYPullDown*)GetObject (6);

    p2->EnableMenu      (type&&(type<=8));
    p3->EnableMenu      (type&&(type!=2));

    DrawChart           ();
}

//__________________________________________________________

void    _HYChartWindow::DrawLegend (_HYRect & plotRect, long chartType, bool isPie)
{
    if ((showLegend != HY_CHART_LEGEND_NONE) && chartType&& (chartType!=2 || isPie ) && (ySeries.lLength || isPie)) {
        // compute the size of the legend box

        _HYStretchCanvas    *sc = (_HYStretchCanvas*)GetObject (0);
        sc->SetFont (labelFont2);

        long  legendWidth  = (labelFont2.size*5)/2,
              legendHeight = (3*ySeries.lLength*labelFont2.size)/2,
              k,
              m,
              s;

        _HYTable            *th  = (_HYTable*)GetObject (isPie?3:2);
        _List seriesLabels;

        if (isPie) {
            legendHeight = (3*(th->verticalSpaces.lLength-1)*labelFont2.size)/2;
            for (k=0; k<th->verticalSpaces.lLength-1; k++) {
                seriesLabels << th->GetCellData (0,k);
            }
        } else
            for (k=0; k<ySeries.lLength; k++) {
                seriesLabels << th->GetCellData (ySeries.lData[k],0);
            }

        m = 0;
        for (k=0; k<seriesLabels.lLength; k++) {
            s = GetVisibleStringWidth (*(_String*)seriesLabels(k), labelFont2);
            if (s>m) {
                m = s;
            }
        }

        legendWidth += m;

        _HYRect     legendBox,
                    shrunkPlot = plotRect;

        switch (showLegend) {
        case HY_CHART_LEGEND_TOP_LEFT:
            legendBox.left   = plotRect.left+5;
            legendBox.top    = 5;
            shrunkPlot.left += legendWidth+5;
            break;

        case HY_CHART_LEGEND_TOP_MID:
            legendBox.left   = (plotRect.right+plotRect.left-legendWidth)/2;
            legendBox.top    = 5;
            shrunkPlot.top  += legendHeight+5;
            break;

        case HY_CHART_LEGEND_TOP_RIGHT:
            legendBox.left   = plotRect.right-5-legendWidth;
            legendBox.top    = 5;
            shrunkPlot.right-= legendWidth+5;
            break;

        case HY_CHART_LEGEND_BOT_LEFT:
            legendBox.left   = plotRect.left+5;
            legendBox.top    = plotRect.bottom-5-legendWidth;
            shrunkPlot.left += legendWidth+5;
            break;

        case HY_CHART_LEGEND_BOT_MID:
            legendBox.left   = (plotRect.right+plotRect.left-legendWidth)/2;
            legendBox.top    = plotRect.bottom-5-legendWidth;
            shrunkPlot.bottom -= legendHeight+5;
            break;

        case HY_CHART_LEGEND_BOT_RIGHT:
            legendBox.left   = plotRect.right-5-legendWidth;
            legendBox.top    = plotRect.bottom-5-legendWidth;
            shrunkPlot.right-= legendWidth+5;
            break;

        case HY_CHART_LEGEND_MID_LEFT:
            legendBox.left   = plotRect.left+5;
            legendBox.top    = (plotRect.bottom+plotRect.top-legendHeight)/2;
            shrunkPlot.left += legendWidth+5;
            break;

        case HY_CHART_LEGEND_MID_RIGHT:
            legendBox.left   = plotRect.right-5-legendWidth;
            legendBox.top    = (plotRect.bottom+plotRect.top-legendHeight)/2;
            shrunkPlot.right-= legendWidth+5;
            break;
        }

        legendBox.right  = legendBox.left +  legendWidth;
        legendBox.bottom = legendBox.top  +  legendHeight;
        legendBox.width = 1;

        if ((shrunkPlot.right-shrunkPlot.left<100)||(shrunkPlot.bottom-shrunkPlot.top<100)) {
            sc->DisplayText ("Can't display legend - the chart is too small",labelFont2.size+2,2,true);
            return;
        }

        plotRect = shrunkPlot;

        _HYColor blk =  LongToHYColor(theColors[HY_CHART_COLOR_LEGEND]);

        sc->SetColor (blk);
        sc->DrawRect (legendBox);

        m = legendBox.left + labelFont2.size/4;
        s = legendBox.top  + (labelFont2.size*5)/4;

        for (k=0; k<seriesLabels.lLength; k++) {
            legendBox.left      = m;
            legendBox.bottom    = s;
            legendBox.right     = m+labelFont2.size;
            legendBox.top       = s-labelFont2.size;

            sc->SetColor (LongToHYColor(theColors[HY_CHART_COLOR_OFFSET+k%HY_CHART_COLOR_COUNT]));
            if (chartType != 3) {
                sc->FillRect (legendBox);
                sc->SetColor (blk);
                sc->DrawRect (legendBox);
            } else {
                switch (k) {
                case 1:
                case 5: {
                    _SimpleList diamond;
                    diamond << legendBox.left;
                    diamond << legendBox.top+labelFont2.size/2;
                    diamond << legendBox.left+labelFont2.size/2;
                    diamond << legendBox.top;
                    diamond << legendBox.right;
                    diamond << legendBox.top+labelFont2.size/2;
                    diamond << legendBox.left+labelFont2.size/2;
                    diamond << legendBox.bottom;
                    if (k==5) {
                        sc->FillPolygon (diamond);
                    } else {
                        sc->DrawPolygon (diamond);
                    }
                    break;
                }
                case 2: {
                    sc->DrawRect (legendBox);
                    break;
                }
                case 3: {
                    sc->FillOval (legendBox);
                    break;
                }
                case 4: {
                    sc->FillRect (legendBox);
                    break;
                }
                default: {
                    sc->DrawOval (legendBox);
                    break;
                }
                }

            }

            sc->DisplayText (*(_String*)seriesLabels(k),s-labelFont2.size/4,m+2*labelFont2.size,true);

            s += (labelFont2.size*3)/2;
        }

    }
}

//__________________________________________________________

void    _HYChartWindow::DrawChart (_HYRect * printRect)
{
    if (suspendDraw) {
        return;
    }

    _HYPullDown         *p2 = (_HYPullDown*)GetObject (5),
                         *p1 = (_HYPullDown*)GetObject (4);
    _HYStretchCanvas    *sc = (_HYStretchCanvas*)GetObject (0);
    _HYRect             plotRect;

    if (!printRect) {
        sc->StartDraw();
        sc->EraseAll();
        plotRect = sc->GetCanvasSize();
    } else {
        plotRect = *printRect;
    }


    long    mc = p1->GetSelection();

    sc->SetColor (LongToHYColor (theColors[HY_CHART_COLOR_WINDOW_BACKGROUND]));
    sc->FillRect (plotRect);

    DrawLegend   (plotRect,mc);
    sc->SetColor ((_HYColor) {
        0,0,0
    });
    _HYFont      bailOut;

    bailOut.face  = "System Font";
    bailOut.size  = 10;
    bailOut.style = HY_FONT_PLAIN;

    sc->SetFont (bailOut);

    switch  (mc) {
    case 1:
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
    case 8: {
        _SimpleList      choices;
        choices << p2->GetSelection()-1;
        choices << ySeries;
        GenerateBarChart        (choices,plotRect,(mc>=3)?mc-2:0);
        break;
    }
    case 2:
        GeneratePieChart (p2->GetSelection(),plotRect);
        break;

    case 9:
    case 10:
    case 11: {
        _SimpleList      choices;
        choices << ySeries;
        Generate3DChart         (choices,plotRect,mc-9);
        break;
    }

    }

    if (!printRect) {
        sc->EndDraw ();
        sc->_MarkForUpdate();
    }
}

//__________________________________________________________

void    _HYChartWindow::GeneratePieChart (long source, _HYRect plotRect)
{
    _HYTable            *t  = (_HYTable*)GetObject (1);
    _HYStretchCanvas    *sc = (_HYStretchCanvas*)GetObject (0);

    DrawLegend   (plotRect,2,true);

    long        k, dividerR = 1, dividerG = 1, dividerB = 1;
    _Matrix     numbers (1,t->verticalSpaces.lLength,false,true);
    _Parameter   startAngle = 0.,
                 endAngle   = 0.;

    if  (source)
        for (k=0; k<numbers.GetVDim(); k++) {
            numbers.theData [k] = ((_String*)t->GetCellData (source-1,k))->toNum();
            endAngle += numbers.theData [k];
        }
    else
        for (k=0; k<numbers.GetVDim(); k++) {
            numbers.theData [k] = k+1;
            endAngle += numbers.theData [k];
        }

    numbers *= (360./endAngle);

    _HYRect      chartRect  = {plotRect.top+20,plotRect.left+20,plotRect.bottom-20,plotRect.right-20,3};

    sc->SetColor (LongToHYColor(theColors[HY_CHART_COLOR_BARBORDER]));
    sc->DrawOval (chartRect);

    chartRect.width = 2;

    for    (k=0; k<numbers.GetVDim(); k++) {
        if (k&&(k%HY_CHART_COLOR_COUNT==0))
            switch ((k/HY_CHART_COLOR_COUNT)%3) {
            case 0:
                dividerR *= 2;
                break;
            case 1:
                dividerB *=2;
                break;
            case 2:
                dividerG *=2;
                break;
            }

        _HYColor    currentColor = LongToHYColor(theColors[HY_CHART_COLOR_OFFSET+k%HY_CHART_COLOR_COUNT]);

        currentColor.R /= dividerR;
        currentColor.B /= dividerB;
        currentColor.G /= dividerG;

        sc->SetColor (currentColor);
        sc->FillArc  (chartRect, startAngle, ceil(numbers.theData[k]));
        startAngle += numbers.theData[k];
    }

    sc->SetColor (LongToHYColor(theColors[HY_CHART_COLOR_BARBORDER]));
    sc->DrawOval (chartRect);
}

//__________________________________________________________

void    _HYChartWindow::Generate3DChart (_SimpleList& columns, _HYRect plotRect, char options)
{
    _HYStretchCanvas    *sc = (_HYStretchCanvas*)GetObject (0);

    _HYTable            *t  = (_HYTable*)GetObject (1);

    long                k,
                        j,
                        m = t->verticalSpaces.lLength,
                        counter,
                        lS;

    _Parameter          minY,
                        maxY,
                        x,
                        yTick,
                        ySc,
                        tp,
                        x3Dtick = -1.,
                        y3Dtick = -1.,
                        x3DMax  = 0.,
                        y3DMax  = 0.;

    _String             label;

    _SimpleList         cache;

    _Matrix             xSpacing (1,options<2?columns.lLength:2,false,true),
                        ySpacing (1,options<2?m:2,false,true);

    if (options == 2 && columns.lLength != 3) {
        label = "Choose exactly 3 columns (X,Y,Value) to produce a 3D histogram.";
        sc->DisplayText (label,15,10,true);
        return;
    } else if (columns.lLength == 0) {
        label = "Choose some data from the 'Y' menu to plot.";
        sc->DisplayText (label,15,10,true);
        return;
    }

    _Matrix             numbers (columns.lLength,m,false,true),
                        unitCube (8,4,false,true),
                        *projP;

    for (j=0; j<columns.lLength; j++)
        for (k=0; k<m; k++) {
            numbers.theData [j*m+k] = ((_String*)t->GetCellData (columns.lData[j],k))->toNum();
        }


    // determine the range of the plot
    minY        = 1e100;
    maxY        = -1e100;
    if (options == 2) {
        xSpacing.theData[0] = 1.e100;
        xSpacing.theData[1] =-1.e100;

        for (j=0; j<m; j++) {
            x = numbers.theData [j];

            if (x>xSpacing.theData[1]) {
                xSpacing.theData[1] = x;
            }
            if (x<xSpacing.theData[0]) {
                xSpacing.theData[0] = x;
            }

            x = numbers.theData [m+j];

            if (x>xSpacing.theData[1]) {
                xSpacing.theData[1] = x;
            }
            if (x<xSpacing.theData[0]) {
                xSpacing.theData[0] = x;
            }

            x = numbers.theData [2*m+j];

            if (x>maxY) {
                maxY = x;
            }
            if (x<minY) {
                minY = x;
            }
        }

        if (xSpacing.theData[0] == xSpacing.theData[1]) {
            xSpacing.theData[1] = xSpacing.theData[0]+0.001;
        }

        ySpacing.theData[0] = xSpacing.theData[0];
        ySpacing.theData[1] = xSpacing.theData[1];
    } else
        for         (k=0; k<m; k++) {
            if (options==0)
                for (j=0; j<columns.lLength; j++) {
                    x = numbers.theData [j*m+k];

                    if (x>maxY) {
                        maxY = x;
                    }
                    if (x<minY) {
                        minY = x;
                    }
                }
            else {
                x = 0;
                for (j=0; j<columns.lLength; j++) {
                    x += numbers.theData [j*m+k];
                }
                if (x>maxY) {
                    maxY = x;
                }
                if (x<minY) {
                    minY = x;
                }
            }
        }

    if (options==1)
        for         (k=0; k<columns.lLength; k++) {
            x = 0;
            for (j=0; j<m; j++) {
                x += numbers.theData [k*m+j];
            }

            if (x>maxY) {
                maxY = x;
            }
            if (x<minY) {
                minY = x;
            }
        }

    if ((minY<0)||(maxY<0)) {
        label = "Negative valued data can't be properly represented by this charting option.";
        sc->DisplayText (label,15,10,true);
        return;
    }


    if ((minY>0)&&(maxY>0)) {
        minY = 0.0;
    }

    if  (maxY==minY) {
        maxY = minY+1.;
    }

    yTick = (maxY-minY);

    if (yTick>0.0) {
        yTick = pow (10,floor(log (yTick)/log (10.)));
        while ((maxY-minY)/yTick<2) {
            yTick /= 2.;
        }
    }



    plotRect.top    += 17;
    plotRect.bottom -= (3*labelFont1.size)/2;
    plotRect.right  -= 10;

    if (minY<0.0) {
        x = - floor (-minY/yTick) * yTick;
    } else {
        x = minY;
    }

    tp = x;
    sc->SetFont (labelFont1);
    counter = 0;
    while (tp<=maxY) {
        label = tp;
        k = GetVisibleStringWidth (label, labelFont1);
        cache << k;
        if (k>counter) {
            counter = k;
        }
        tp += yTick;
    }

    plotRect.left    += counter+7;
    plotRect.right   -= counter+7;
    plotRect.bottom  -= 2;
    plotRect.top     += 2;

    lS = counter+7;

    if (maxY-minY) {
        ySc = 1./(maxY-minY);
    }

    /* calculate scaling constants */


    // (0,0,0,1) [0-3]
    unitCube.theData[3] = 1;
    // (0,0,1,1) [4-7]
    unitCube.theData[6] = 1;
    unitCube.theData[7] = 1;
    // (0,1,0,1) [8-11]
    unitCube.theData[9]  = 1;
    unitCube.theData[11] = 1;
    // (0,1,1,1) [12-15]
    unitCube.theData[13] = 1;
    unitCube.theData[14] = 1;
    unitCube.theData[15] = 1;
    // (1,0,0,1) [16-19]
    unitCube.theData[16] = 1;
    unitCube.theData[19] = 1;
    // (1,0,1,1) [20-23]
    unitCube.theData[20] = 1;
    unitCube.theData[22] = 1;
    unitCube.theData[23] = 1;
    // (1,0,1,1) [24-27]
    unitCube.theData[24] = 1;
    unitCube.theData[26] = 1;
    unitCube.theData[27] = 1;
    // (1,1,1,1) [28-31]
    unitCube.theData[28] = 1;
    unitCube.theData[29] = 1;
    unitCube.theData[30] = 1;
    unitCube.theData[31] = 1;

    xShift = 1e100;
    yShift = 1e100;
    xScale  = 0.0;
    yScale  = 0.0;

    projP = (_Matrix*)unitCube.MultObj (projectionMatrix);

    for (counter = 0; counter < 32; counter+=4) {
        _Parameter scaling = 1./projP->theData[counter+3],
                   a,
                   b;

        a = (projP->theData[counter]     *= scaling);
        b = (projP->theData[counter+1]   *= scaling);

        if (a<xShift) {
            xShift = a;
        }
        if (b<yShift) {
            yShift = b;
        }
    }

    for (counter = 0; counter < 32; counter+=4) {
        _Parameter a = projP->theData[counter],
                   b = projP->theData[counter+1];

        if (a-xShift>xScale) {
            xScale = a-xShift;
        }
        if (b-yShift>yScale) {
            yScale = b-yShift;
        }
    }


    xScale = (plotRect.right-plotRect.left)/xScale;
    yScale = (plotRect.bottom-plotRect.top)/yScale;

    /*if (yScale>xScale)
        yScale = xScale;
    if (xScale>yScale)
        xScale = yScale;*/

    xShift *= xScale;
    yShift *= yScale;
    DeleteObject (projP);

    /* plot the coordinate walls */

    // back wall
    _Matrix      planeP (4,3,false, true);
    planeP.theData[0] = 0;
    planeP.theData[1] = 1;
    planeP.theData[2] = 1;
    planeP.theData[3] = 1;
    planeP.theData[4] = 1;
    planeP.theData[5] = 1;
    planeP.theData[6] = 0;
    planeP.theData[7] = 0;
    planeP.theData[8] = 1;
    planeP.theData[9] = 1;
    planeP.theData[10] = 0;
    planeP.theData[11] = 1;

    sc->SetColor (LongToHYColor (theColors(HY_CHART_COLOR_3D_BACK)));

    DrawAPlane   (planeP,0,1,plotRect);

    // left wall
    planeP.theData[2]  = 0;
    planeP.theData[3]  = 0;
    planeP.theData[8]  = 0;
    planeP.theData[9]  = 0;
    planeP.theData[10] = 0;
    planeP.theData[11] = 1;

    sc->SetColor (LongToHYColor (theColors(HY_CHART_COLOR_3D_LEFT)));

    DrawAPlane   (planeP,0,1,plotRect);
    // ground wall

    planeP.theData[0] = 1;
    planeP.theData[1] = 0;
    planeP.theData[4] = 0;
    planeP.theData[5] = 0;
    planeP.theData[6] = 1;
    planeP.theData[8] = 1;

    sc->SetColor (LongToHYColor (theColors(HY_CHART_COLOR_3D_FLOOR)));
    DrawAPlane   (planeP,0,1,plotRect);

    tp = x;
    k  = 0;

    // plot z lines and ticks

    sc->SetColor (LongToHYColor (theColors(HY_CHART_COLOR_3D_VALUELINE)));
    while (tp<=maxY) {
        label = tp;
        _Matrix coordLine (2,3,false,true);
        coordLine.theData[0] = 1;
        coordLine.theData[1] = coordLine.theData[4] = (tp-minY)*ySc;
        coordLine.theData[2] = coordLine.theData[5] = 1;
        coordLine.theData[3] = 0;
        DrawALine   (coordLine,0,1.,plotRect,1,&label,4,0);
        coordLine.theData[0] = coordLine.theData[2] = 0;
        DrawALine   (coordLine,0,1.,plotRect,1,&label,-4-cache.lData[k++],0);
        tp += yTick;
    }


    _Parameter      stashMY   = maxY,
                    stashMinY = minY;


    if (options == 2)
        // draw a 3D Historgram
    {
        minY = xSpacing.theData[1]-xSpacing.theData[0];

        x =  xSpacing.theData[0]/minY;

        x3Dtick = pow (10,floor(log (minY)/log (10.)));

        while (minY/x3Dtick<2) {
            x3Dtick /= 2.;
        }

        if (x>maxY) {
            x=maxY;
        } else {
            maxY=x;
        }

        tp     = 1./minY;
        yTick  = 1./minY;

        x3DMax  = xSpacing.theData[1];
        y3DMax  = x3DMax;
        y3Dtick = x3Dtick;


        for (j=0; j<m; j++) {
            numbers.theData [j]   = (numbers.theData [j]-xSpacing.theData[0])/minY;
            numbers.theData [m+j] = (numbers.theData [m+j]-xSpacing.theData[0])/minY;
        }

        /* set those for coordinate gridlines */

        j = 0;
        k = 0;

        minY = 0.0;
        maxY = 0.0;
    } else {
        // now draw the bars, starting with the back row, left to right

        minY = 1./(2*m);               // spacing per column width-wise
        maxY = 1./(2*columns.lLength); // spacing per column depth wise

        j = -1;
        k = -1;

        if ((xAxis3DScale>=0)&&(xAxis3DScale<t->horizontalSpaces.lLength)&&(columns.lLength>1)) {
            if (t->verticalSpaces.lLength>=columns.lLength) {
                for (counter = 0; counter < columns.lLength; counter++) {
                    xSpacing.theData[counter] = ((_String*)t->GetCellData (xAxis3DScale,counter))->toNum();
                }


                x = 1e100;

                for (counter = 1; counter < columns.lLength; counter++) {
                    yTick = (xSpacing.theData[counter]-xSpacing.theData[counter-1])/2;
                    if (yTick<=0) {
                        break;
                    }
                    if (yTick<x) {
                        x = yTick;
                    }
                }

                if (counter == columns.lLength) { // valid coordinates
                    x /= (xSpacing.theData[columns.lLength-1]-xSpacing.theData[0]);

                    x3Dtick = (xSpacing.theData[columns.lLength-1]-xSpacing.theData[0]);

                    x3Dtick = pow (10,floor(log (x3Dtick)/log (10.)));
                    while ((xSpacing.theData[columns.lLength-1]-xSpacing.theData[0])/x3Dtick<2) {
                        x3Dtick /= 2.;
                    }


                    if (x>maxY) {
                        x=maxY;
                    } else {
                        maxY=x;
                    }

                    tp = (1-2*maxY)/(xSpacing.theData[columns.lLength-1]-xSpacing.theData[0]);
                    x3DMax = xSpacing.theData[columns.lLength-1];
                    j = 0;
                }
            }
        }

        if ((yAxis3DScale>=0)&&(yAxis3DScale<t->horizontalSpaces.lLength)&&(m>1)) {
            if (t->verticalSpaces.lLength>=m) {
                for (counter = 0; counter < m; counter++) {
                    ySpacing.theData[counter] = ((_String*)t->GetCellData (yAxis3DScale,counter))->toNum();
                }


                x = 1e100;

                for (counter = 1; counter < m; counter++) {
                    yTick = (ySpacing.theData[counter]-ySpacing.theData[counter-1])/2;
                    if (yTick<=0) {
                        break;
                    }
                    if (yTick<x) {
                        x = yTick;
                    }
                }

                if (counter == m) {
                    x /= (ySpacing.theData[m-1]-ySpacing.theData[0]);

                    y3Dtick = (ySpacing.theData[m-1]-ySpacing.theData[0]);

                    y3Dtick = pow (10,floor(log (y3Dtick)/log (10.)));
                    while ((ySpacing.theData[m-1]-ySpacing.theData[0])/y3Dtick<2) {
                        y3Dtick /= 2.;
                    }

                    if (x>minY) {
                        x=minY;
                    } else {
                        minY=x;
                    }

                    y3DMax = ySpacing.theData[m-1];
                    yTick = (1-2*minY)/(ySpacing.theData[m-1]-ySpacing.theData[0]);
                    k = 0;
                }
            }
        }
    }


    if ((j>=0)&&(k>=0)) {
        if (yTick>tp) {
            minY *= tp/yTick;
            yTick = tp;
        } else {
            maxY *= yTick/tp;
            tp = yTick;
        }

        if (maxY>minY) {
            maxY = minY;
        } else {
            minY = maxY;
        }

        if (x3Dtick>y3Dtick) {
            y3Dtick = x3Dtick;
        } else {
            x3Dtick = y3Dtick;
        }

        if (y3DMax-ySpacing.theData[0]>x3DMax-xSpacing.theData[0]) {
            x3DMax = y3DMax-ySpacing.theData[0]+xSpacing.theData[0];
        } else {
            y3DMax = x3DMax-xSpacing.theData[0]+ySpacing.theData[0];
        }
    }

    if (j<0) {
        for (counter = 0; counter < columns.lLength; counter ++) {
            xSpacing.theData [columns.lLength-counter-1] = 1- (2*counter+1)*maxY;
        }
    } else {
        if (x3Dtick>0.0) {
            x = ceil (xSpacing.theData[0]/x3Dtick) * x3Dtick;

            while (x<=x3DMax) {
                label = x;
                _Matrix coordLine (2,3,false,true);

                coordLine.theData[0] = 1;
                coordLine.theData[1] = coordLine.theData[4] = 0; // z-coordinate
                coordLine.theData[2] = coordLine.theData[5] = (x-xSpacing.theData[0])*tp + maxY;
                coordLine.theData[3] = 0;

                DrawALine   (coordLine,0,1.,plotRect,1,&label,labelFont1.size,labelFont1.size);
                x += x3Dtick;
            }
        }

        if (options < 2) {
            for (counter = 1; counter < columns.lLength; counter++) {
                xSpacing.theData[counter] = (xSpacing.theData[counter]-xSpacing.theData[0]) * tp + maxY;
            }

            xSpacing.theData[0] = maxY;
        }
    }

    if (k<0) {
        for (counter = 0; counter < m; counter ++) {
            ySpacing.theData [counter] = minY + 2*counter*minY;
        }
    } else {
        if (y3Dtick>0.0) {

            x = ceil (ySpacing.theData[0]/y3Dtick) * y3Dtick;

            while (x<=y3DMax) {
                label = x;
                _Matrix coordLine (2,3,false,true);

                coordLine.theData[0] = coordLine.theData[3] = (x-ySpacing.theData[0])*yTick + minY;
                coordLine.theData[1] = coordLine.theData[4] = 0; // z-coordinate

                coordLine.theData[2] = 0;
                coordLine.theData[5] = 1;

                DrawALine   (coordLine,0,1.,plotRect,1,&label,-GetVisibleStringWidth (label, labelFont1) - labelFont1.size/2,labelFont1.size);
                x += y3Dtick;
            }
        }

        if (options < 2) {
            for (counter = 1; counter < m; counter++) {
                ySpacing.theData[counter] = (ySpacing.theData[counter]-ySpacing.theData[0]) * yTick + minY;
            }

            ySpacing.theData[0] = minY;
        }
    }

    maxY/=2;
    minY/=2;


    if (options==1)
        // generate coordinate histogram
    {
        sc->SetColor (LongToHYColor (theColors(HY_CHART_COLOR_3D_HISTOGRAM)));
        // back wall first

        planeP.theData[1]  =
            planeP.theData[4]  = 0;

        planeP.theData[2]  =
            planeP.theData[5]  =
                planeP.theData[8]  =
                    planeP.theData[11] = 1;

        for (j=0; j<m; j++) {
            planeP.theData[0]  =
                planeP.theData[6]  = ySpacing.theData [j] - minY;
            planeP.theData[3]  =
                planeP.theData[9]  = ySpacing.theData [j] + minY;

            _Parameter sum = 0;
            for (counter = 0; counter < columns.lLength; counter++) {
                sum += numbers.theData [j+counter*m];
            }

            planeP.theData[7]  =
                planeP.theData[10] = sum*ySc;

            DrawAPlane (planeP,0,1.,plotRect);

        }

        // now - the left wall
        planeP.theData[0] =
            planeP.theData[1] =
                planeP.theData[3] =
                    planeP.theData[4] =
                        planeP.theData[6] =
                            planeP.theData[9] = 0;


        for (j=columns.lLength-1; j>=0; j--) {
            planeP.theData[2] =
                planeP.theData[8] = xSpacing.theData [j] - maxY;

            planeP.theData[5]  =
                planeP.theData[11] = xSpacing.theData[j] + maxY;

            _Parameter sum = 0;
            for (counter = 0; counter < m; counter++) {
                sum += numbers.theData [m*j+counter];
            }

            planeP.theData[7] = planeP.theData[10] = sum*ySc;
            DrawAPlane (planeP,0,1.,plotRect);
        }
    }

    // handle possible spacing issues
    if (options < 2)
        for (counter = columns.lLength-1; counter >=0; counter --) {
            // current position witdh wise
            sc->SetColor (LongToHYColor(theColors[HY_CHART_COLOR_OFFSET+counter%HY_CHART_COLOR_COUNT]));

            unitCube.theData[2]  = unitCube.theData[5]  = unitCube.theData[14] = unitCube.theData[17] =
                                       xSpacing.theData[counter] - maxY;

            unitCube.theData[8]  = unitCube.theData[11]  = unitCube.theData[20] = unitCube.theData[23] =
                                       xSpacing.theData[counter] + maxY;

            unitCube.theData[0]  =
                unitCube.theData[6]  =
                    unitCube.theData[12] =
                        unitCube.theData[18] = ySpacing.theData[0] - minY;

            unitCube.theData[3]  =
                unitCube.theData[9]  =
                    unitCube.theData[15] =
                        unitCube.theData[21] = ySpacing.theData[0] + minY;

            unitCube.theData[13] = 0;
            unitCube.theData[16] = 0;
            unitCube.theData[19] = 0;
            unitCube.theData[22] = 0;

            for (long counter2 = 0; counter2 < m; counter2++) {
                if (counter2) {
                    for (j=0; j<24; j+=6) {
                        unitCube.theData[j] = ySpacing.theData[counter2] - minY;
                    }

                    for (j=3; j<24; j+=6) {
                        unitCube.theData[j] = ySpacing.theData[counter2] + minY;
                    }
                }

                unitCube.theData[1] = unitCube.theData[4] = unitCube.theData[7] = unitCube.theData[10] =
                                          numbers.theData [counter2+counter*m]*ySc;

                DrawAParallelogram (unitCube,0,1,plotRect);

            }
            //x -= 2*maxY;
            // current position depth wise
        }
    else {
        sc->SetColor (LongToHYColor(theColors[HY_CHART_COLOR_OFFSET]));

        for (counter = 0; counter<m; counter=counter+1) {

            unitCube.theData[2]  = unitCube.theData[5]  = unitCube.theData[14] = unitCube.theData[17] =
                                       numbers.theData [counter]-0.01;

            unitCube.theData[8]  = unitCube.theData[11]  = unitCube.theData[20] = unitCube.theData[23] =
                                       numbers.theData [counter]+0.01;

            unitCube.theData[0]  =
                unitCube.theData[6]  =
                    unitCube.theData[12] =
                        unitCube.theData[18] = numbers.theData [counter+m]-0.01;

            unitCube.theData[3]  =
                unitCube.theData[9]  =
                    unitCube.theData[15] =
                        unitCube.theData[21] = numbers.theData [counter+m]+0.01;

            unitCube.theData[13] =
                unitCube.theData[16] =
                    unitCube.theData[19] =
                        unitCube.theData[22] = numbers.theData [counter+2*m]*ySc-0.01;

            unitCube.theData[1] =
                unitCube.theData[4] =
                    unitCube.theData[7] =
                        unitCube.theData[10] = numbers.theData [counter+2*m]*ySc+0.01;;

            DrawAParallelogram (unitCube,0,1,plotRect);
        }

    }
    // produce an overplot

    if (!overlayPlot.IsEmpty()) {
        _Matrix         plotValues    (surfaceDivs,surfaceDivs,false,true),
                        plotBase      (2,surfaceDivs,false,true);

        for (j=0; j<surfaceDivs; j++) {
            plotBase.Store (0,j,.1*j);
            plotBase.Store (1,j,.1*j);
        }

        if (x3DMax!=0.0) {
            x3DMax = (xSpacing.theData[0]-maxY)/tp;

            for (j=0; j<surfaceDivs; j++) {
                plotBase.theData[j] = plotBase.theData[j]/tp + x3DMax;
            }
        }

        if (y3DMax!=0.0) {
            y3DMax = (xSpacing.theData[0]-minY)/yTick;

            for (j=surfaceDivs; j<2*surfaceDivs; j++) {
                plotBase.theData[j] = plotBase.theData[j]/yTick + y3DMax;
            }
        }

        label = "_x_";
        _Variable*  xVar = CheckReceptacle (&label, empty);
        label = "_y_";
        _Variable*  yVar = CheckReceptacle (&label, empty);

        for (j=0; j<surfaceDivs; j++) {
            xVar->SetValue (new _Constant(plotBase.theData[j]), false);
            for (k=0; k<surfaceDivs; k++) {
                yVar->SetValue (new _Constant (plotBase.theData[surfaceDivs+k]), false);
                y3DMax = overlayPlot.Compute()->Value();
                if (y3DMax < 0) {
                    y3DMax = 0.;
                } else if (y3DMax > stashMY) {
                    y3DMax = stashMY;
                }
                plotValues.Store(j,k,(y3DMax-stashMinY)*ySc);
            }
        }

        // now plot the thing

        x3DMax = (1-2*maxY)/(surfaceDivs-1.);
        y3DMax = (1-2*minY)/(surfaceDivs-1.);

        sc->SetColor (LongToHYColor(theColors[HY_CHART_COLOR_OVERPLOT]));

        for (j=0; j<surfaceDivs-1; j++) {

            planeP.theData[0] =
                planeP.theData[3] = maxY+j*x3DMax;

            planeP.theData[2] =
                planeP.theData[8] = 1-minY-y3DMax;

            planeP.theData[5] =
                planeP.theData[11] = 1-minY;

            planeP.theData[6] =
                planeP.theData[9] = planeP.theData[0]+x3DMax;


            for (k=surfaceDivs-2; k>=0; k--) {
                planeP.theData[1]  = plotValues (j,k);
                planeP.theData[4]  = plotValues (j,k+1);
                planeP.theData[7]  = plotValues (j+1,k);
                planeP.theData[10] = plotValues (j+1,k+1);

                DrawAPlane (planeP,2,1.,plotRect);

                planeP.theData[2] -= y3DMax;
                planeP.theData[5] -= y3DMax;
                planeP.theData[8] -= y3DMax;
                planeP.theData[11]-= y3DMax;
            }
        }

    }
}


//__________________________________________________________

void    _HYChartWindow::GenerateBarChart (_SimpleList& columns, _HYRect plotRect, char options)
{
    _HYTable            *t  = (_HYTable*)GetObject (1);
    _HYStretchCanvas    *sc = (_HYStretchCanvas*)GetObject (0);

    if ((xLabel.sLength+zLabel.sLength)||(yLabel.sLength&&(options==5))) {
        if ((xLabel.sLength)||(yLabel.sLength&&(options==5))) {
            plotRect.bottom -= labelFont3.size;
        }

        if (zLabel.sLength) {
            plotRect.top += (3*labelFont3.size)/2-8;
        }

        if (plotRect.bottom-plotRect.top<50) {
            sc->DisplayText (_String("Can't display x-axis label: chart is too short"),15,2,true);
            return;
        }
    }


    long        k,
                j,
                m = t->verticalSpaces.lLength,
                barWidth = 10,
                s,
                lS;

    _Matrix     numbers (columns.lLength,m,false,true);

    _Parameter  minX,
                maxX,
                minY,
                maxY,
                x,
                yTick,
                //xTick,
                yScale,
                xScale,
                tp;

    _String     label;

    _HYRect     lr;

    _SimpleList cache;

    if (columns.lLength < 2) {
        label = "Choose some data from the 'Y' menu to plot.";
        sc->DisplayText (label,15,10,true);
        return;
    }

    if ((options==5)&&(columns.lLength!=3)) {
        label = "Choose two columns of positive data from the 'Y' menu to plot.";
        sc->DisplayText (label,15,10,true);
        return;
    }

    if ((options==6)&&(columns.lLength!=4)&&(columns.lLength!=5)) {
        label = "Choose three or four columns (low, high, median, optional extra point) of data from the 'Y' menu to plot.";
        sc->DisplayText (label,15,10,true);
        return;
    }


    for (j=0; j<columns.lLength; j++) {
        if (columns.lData[j]==-1)
            for (k=0; k<m; k++) {
                numbers.theData [j*m+k] = k+1;
            }
        else
            for (k=0; k<m; k++) {
                numbers.theData [j*m+k] = ((_String*)t->GetCellData (columns.lData[j],k))->toNum();
            }
    }

    if (options==5)
        for (k=0; k<m; k++) {
            numbers.theData [2*m+k] *= -1;
        }

    _SimpleList  barPositions;


    // determine the range of the plot
    minX        = 1e100;
    maxX        = -1e100;
    minY        = 1e100;
    maxY        = -1e100;
    for         (k=0; k<m; k++) {
        x = numbers.theData[k];
        if (isfinite (x)) {
            if (x>maxX) {
                maxX = x;
            }
            if (x<minX) {
                minX = x;
            }
            if (options!=4)
                for (j=1; j<columns.lLength; j++) {
                    x = numbers.theData [j*m+k];
                    if (isfinite(x)) {
                        if (x>maxY) {
                            maxY = x;
                        }
                        if (x<minY) {
                            minY = x;
                        }
                    }
                }
            else {
                x = 0.;
                for (j=1; j<columns.lLength; j++) {
                    x += numbers.theData [j*m+k];
                }

                if (x>maxY) {
                    maxY = x;
                }
                if (x<minY) {
                    minY = x;
                }
            }
        }
    }

    if (minY>0 && maxY>0) {
        minY = 0.0;
    } else if (minY<0 && maxY<0) {
        maxY = 0.0;
    }

    if (options==5) {
        if (-minY>maxY) {
            maxY = -minY;
        } else {
            minY = -maxY;
        }
    }

    if (maxY==minY) {
        maxY = minY+1;
    }

    if (!CheckEqual (userMin, userMax)) {
        minY = userMin;
        maxY = userMax;
    }



    yTick = (maxY-minY);

    if (yTick>0.0) {
        yTick = pow (10,floor(log (yTick)/log (10.)));
        /*while ((maxY-minY)/yTick>5)
            yTick *= 2.;*/
        while ((maxY-minY)/yTick<2) {
            yTick /= 2.;
        }
    }



//  plotRect = sc->GetCanvasSize();

    plotRect.top    += 17;
    plotRect.bottom -= (labelFont1.size*3)/2;
    plotRect.right  -= 10;

    if (minY<0.0) {
        x = - floor (-minY/yTick) * yTick;
    } else {
        x = minY;
    }

    tp = x;
    sc->SetFont (labelFont1);
    while (tp<=maxY) {
        if (options==5 && tp<0) {
            label = -tp;
        } else {
            label = tp;
        }

        k = GetVisibleStringWidth (label, labelFont1);
        cache << k;
        if (k>plotRect.left) {
            plotRect.left = k;
            lS = k;
        }
        tp += yTick;
    }

    plotRect.left += 7;
    lS += 7;

    if (maxY-minY) {
        yScale = (plotRect.bottom-plotRect.top-4)/(maxY-minY);
    }




    sc->SetColor (LongToHYColor (theColors[HY_CHART_COLOR_BACKGROUND]));
    sc->FillRect (plotRect);

    plotRect.bottom -= 2;
    plotRect.top +=2;

    sc->SetColor (LongToHYColor (theColors[HY_CHART_COLOR_COORD_LABEL]));
    tp = x;
    k=0;
    while (tp<=maxY) {
        if ((options==5)&&(tp<0)) {
            label = -tp;
        } else {
            label = tp;
        }

        sc->DisplayText (label,plotRect.bottom - (tp-minY)*yScale,plotRect.left-4-cache.lData[k++],true);
        tp += yTick;
    }
    lr.left = plotRect.left;
    lr.right = plotRect.right;
    lr.width = 1;

    sc->SetColor(LongToHYColor (theColors[HY_CHART_COLOR_VALUELINE]));
    tp = x;
    while (tp<=maxY) {
        lr.top = lr.bottom = plotRect.bottom - (tp-minY)*yScale;
        sc->DrawLine (lr);
        tp += yTick;
    }

    if ((options == 0)||(options == 4)||(options == 5)) {
        if (options) {
            s = 1;
        } else {
            s = columns.lLength - 1;
        }

        barWidth = (plotRect.right-plotRect.left)/(s*m+.5*(s+1));

        if (barWidth <= 1) {
            barWidth = 1;
        }

        plotRect.left  += barWidth*s/2+5;
        plotRect.right -= barWidth*s/2+5;
    }

    if (maxX-minX) {
        xScale = (plotRect.right-plotRect.left)/(maxX-minX);
    } else {
        xScale = 0.0;
    }


    if (options&&(options!=4)&&(options!=5)) {
        _SimpleList sortedOrder;

        for (j=0; j<m; j++) {
            sortedOrder << j;
        }

        bool doneSorting = false;

        while (!doneSorting) {
            doneSorting = true;
            for (j=0; j<m-1; j++) {
                if (numbers.theData [sortedOrder.lData[j]] > numbers.theData [sortedOrder.lData[j+1]]) {
                    k = sortedOrder.lData[j];
                    sortedOrder.lData[j] = sortedOrder.lData[j+1];
                    sortedOrder.lData[j+1] = k;
                    doneSorting = false;
                }
            }
        }
        lr.width = 1;
        for (j=1; j<columns.lLength; j++) {
            sc->SetColor (LongToHYColor(theColors[HY_CHART_COLOR_OFFSET+(j-1)%HY_CHART_COLOR_COUNT]));
            for (k=0; k<m; k++) {
                x = (numbers.theData[sortedOrder.lData[k]]-minX)*xScale;
                lr.left   = plotRect.left + x - 1;
                lr.top    = plotRect.bottom-(numbers.theData [j*m+sortedOrder.lData[k]]-minY)*yScale;
                if (options==2) {
                    if (k!=m-1) {
                        lr.right = plotRect.left + (numbers.theData[sortedOrder.lData[k+1]]-minX)*xScale;
                        lr.bottom = plotRect.bottom-(numbers.theData [j*m+sortedOrder.lData[k+1]]-minY)*yScale;
                        lr.width = 2;
                        sc->DrawLine (lr);
                    }
                    lr.bottom = lr.top + 3;
                    lr.top   -= 3;
                    lr.right  = lr.left +3;
                    lr.left  -= 3;
                    lr.width = 1;
                    sc->FillRect (lr);
                } else {
                    if (options == 3) {
                        if (k!=m-1) {
                            lr.right  =  plotRect.left + (numbers.theData[sortedOrder.lData[k+1]]-minX)*xScale - 1;
                            lr.bottom = lr.top;
                            lr.width  = 2;
                            sc->DrawLine (lr);
                            lr.left   = lr.right;
                            lr.bottom = plotRect.bottom-(numbers.theData [j*m+sortedOrder.lData[k+1]]-minY)*yScale;
                            sc->DrawLine (lr);
                        }
                    } else {
                        lr.left = lr.right =  plotRect.left + (numbers.theData[sortedOrder.lData[k]]-minX)*xScale - 1;
                        lr.top = lr.bottom = plotRect.bottom-(numbers.theData [j*m+sortedOrder.lData[k]]-minY)*yScale;

                        if (options == 6) {
                            if (j==4) {
                                _SimpleList triag;
                                triag << lr.left-4;
                                triag << lr.bottom;
                                triag << lr.left;
                                triag << lr.bottom-4;
                                triag << lr.left+4;
                                triag << lr.bottom;

                                sc->FillPolygon (triag);
                            } else {
                                if (j==3) {
                                    lr.top      -=  3;
                                    lr.bottom   +=  4;
                                    lr.right    +=  4;
                                    lr.left     -=  3;

                                    sc->FillOval (lr);
                                } else {

                                    long tcoord = plotRect.bottom-(numbers.theData [m*2+sortedOrder.lData[k]]-minY)*yScale;

                                    lr.left     -=  3;
                                    lr.right    +=  3;

                                    sc->DrawLine (lr);
                                    lr.bottom ++;
                                    lr.top ++;
                                    sc->DrawLine (lr);
                                    lr.bottom--;
                                    lr.top++;

                                    lr.left     +=  3;
                                    lr.right    -=  3;

                                    if (tcoord < lr.bottom) {
                                        lr.top = tcoord;
                                    } else {
                                        lr.bottom = tcoord;
                                    }

                                    if (lr.bottom>=lr.top) {
                                        sc->DrawLine (lr);
                                    }
                                }
                            }
                        } else {
                            lr.top      -=3;
                            lr.bottom   +=3;
                            lr.right    +=3;
                            lr.left     -=3;
                            switch (j) {
                            case 2:
                            case 6: {
                                _SimpleList diamond;
                                diamond << lr.left;
                                diamond << lr.top+3;
                                diamond << lr.left+3;
                                diamond << lr.top;
                                diamond << lr.right;
                                diamond << lr.top+3;
                                diamond << lr.left+3;
                                diamond << lr.bottom;
                                if (j==6) {
                                    sc->FillPolygon (diamond);
                                } else {
                                    sc->DrawPolygon (diamond);
                                }
                                break;
                            }
                            case 3: {
                                sc->DrawRect (lr);
                                break;
                            }
                            case 4: {
                                sc->FillOval (lr);
                                break;
                            }
                            case 5: {
                                sc->FillRect (lr);
                                break;
                            }
                            default: {
                                sc->DrawOval (lr);
                                break;
                            }
                            }
                        }
                    }
                }
            }
        }
    } else {
        lr.width = 1;
        if (options==4) {
            _Matrix     sums    (1,m,false,true);

            yTick = -barWidth*.5;

            for (j=1; j<columns.lLength; j++) {
                sc->SetColor (LongToHYColor(theColors[HY_CHART_COLOR_OFFSET+(j-1)%HY_CHART_COLOR_COUNT]));
                for (k=0; k<m; k++) {
                    x = (numbers.theData[k]-minX)*xScale;
                    lr.left = plotRect.left+x+yTick+1;
                    lr.right = lr.left+barWidth;
                    tp = numbers.theData [j*m+k];

                    if (tp>=0.0) {
                        lr.bottom = plotRect.bottom+(minY-sums.theData[k])*yScale;
                        lr.top    = lr.bottom - tp*yScale;
                    } else {
                        lr.top    = plotRect.bottom+minY*yScale;
                        lr.bottom = plotRect.bottom-(tp-minY-sums.theData[k])*yScale;
                    }

                    sc->FillRect (lr);
                    sums.theData[k]+=tp;

                    if ((barWidth>2)&&(j==columns.lLength-1)) {
                        tp = sums.theData[k];
                        if (tp>=0.0) {
                            lr.bottom = plotRect.bottom+minY*yScale;
                            lr.top = lr.bottom - tp*yScale;
                        } else {
                            lr.top    = plotRect.bottom+minY*yScale;
                            lr.bottom = plotRect.bottom-(tp-minY)*yScale;
                        }
                        sc->SetColor (LongToHYColor(theColors[HY_CHART_COLOR_BARBORDER]));
                        sc->DrawRect (lr);
                        sc->SetColor (LongToHYColor(theColors[HY_CHART_COLOR_OFFSET+(j-1)%HY_CHART_COLOR_COUNT]));
                    }
                }
            }
        } else {

            // check for custom colors

            _Variable * ucs                 = FetchVar (LocateVarByName(chartCustomColumnColors));
            _Matrix   * customColColors     = nil;
            bool        tags = false;

            if (ucs && ucs->ObjectClass () == MATRIX) {
                customColColors = (_Matrix*)ucs->Compute();
                if (customColColors->GetHDim()!=m || customColColors->GetVDim()<columns.lLength-1 || customColColors->IsAStringMatrix()) {
                    customColColors = nil;
                } else {
                    tags = (options == 5 && customColColors->GetVDim() == 4);
                }
            }

            long  checkBarWidth = barWidth;

            for (k=0; k<m; k++) {
                x = (numbers.theData[k]-minX)*xScale;

                lr.left = x;

                if (k) {
                    if (lr.left<=lr.right) {
                        if (options==5) {
                            checkBarWidth -= lr.right-lr.left+1;
                        } else {
                            checkBarWidth -= ceil((lr.right-lr.left)/((_Parameter)columns.lLength-1)) + 1;
                        }

                        if (checkBarWidth<1) {
                            checkBarWidth = 1;
                        }

                    }
                }

                if (options==5) {
                    lr.right = lr.left+checkBarWidth;
                } else {
                    lr.right = lr.left+(columns.lLength-1)*checkBarWidth;
                }
            }

            for (j=1; j<columns.lLength; j++) {
                if (options==5) {
                    yTick = -.5*checkBarWidth;
                } else {
                    yTick = checkBarWidth*(j-(columns.lLength-1.)/2-1);
                }

                if (!customColColors) {
                    sc->SetColor (LongToHYColor(theColors[HY_CHART_COLOR_OFFSET+(j-1)%HY_CHART_COLOR_COUNT]));
                }

                for (k=0; k<m; k++) {
                    x = (numbers.theData[k]-minX)*xScale;

                    if (options==5) {
                        lr.left = plotRect.left + x + yTick + 1;
                    } else {
                        lr.left = plotRect.left + x + yTick-j+1;
                    }

                    lr.right = lr.left+checkBarWidth;
                    tp = numbers.theData [j*m+k];
                    if (tp>=0.0) {
                        lr.bottom = plotRect.bottom+minY*yScale;
                        lr.top = lr.bottom - tp*yScale;
                    } else {
                        lr.top    = plotRect.bottom+minY*yScale;
                        lr.bottom = plotRect.bottom-(tp-minY)*yScale;
                    }

                    if (customColColors) {
                        sc->SetColor (LongToHYColor((*customColColors)(k,j-1)));
                    }

                    sc->FillRect (lr);

                    if (checkBarWidth>2 && !tags) {
                        sc->SetColor (LongToHYColor(theColors[HY_CHART_COLOR_BARBORDER]));
                        sc->DrawRect (lr);
                        if (!customColColors) {
                            sc->SetColor (LongToHYColor(theColors[HY_CHART_COLOR_OFFSET+(j-1)%HY_CHART_COLOR_COUNT]));
                        }
                    }

                    if (tags && (*customColColors)(k,1+j) > 0.0) {
                        sc->SetColor (black);
                        long     diam = MIN(lr.bottom-lr.top,lr.right-lr.left);
                        _HYRect  circRect;
                        circRect.left  = (lr.right+lr.left-diam)/2;
                        circRect.right = circRect.left + diam;
                        if (j==1) {
                            circRect.top    = lr.top - diam*3/2;
                            circRect.bottom = circRect.top+diam;
                        } else {
                            circRect.bottom = lr.bottom + diam*3/2;
                            circRect.top    = circRect.bottom-diam;
                        }
                        sc->FillOval (circRect);
                    }
                }
            }
        }

        plotRect.left  -= barWidth*s/2+5;
        plotRect.right += barWidth*s/2+5;
    }

    lr.width = 2;

    if ((options!=5)&&(!overlayPlot.IsEmpty())) {
        label = "_x_";

        _Variable* thisVar = CheckReceptacle(&label,label,false);

        _Parameter   ox1 = minX,
                     ox2,
                     shifter = 0;

        ox2 = ox1+1./xScale;

        if ((columns.lData[0]==-1)&&((options==0)||(options==4)||(options==5))) {
            //maxX += .5;
            shifter = -.5;
        }


        setParameter (label,ox1);

        x = overlayPlot.Compute ()->Value();

        if (terminateExecution) {
            terminateExecution = false;
        }

        else {

            sc->SetColor (LongToHYColor(theColors[HY_CHART_COLOR_OVERPLOT]));

            while (1) {
                if (ox2>maxX) {
                    ox2 = maxX;
                }

                _Constant cv (ox2);

                thisVar->SetValue (&cv);

                yTick = overlayPlot.Compute ()->Value();

                if (terminateExecution) {
                    terminateExecution = false;
                    break;
                }

                lr.left  = (ox1-(minX+shifter))*xScale+plotRect.left;
                lr.right = (ox2-(minX+shifter))*xScale+plotRect.left;

                bool      plotThisPair = false;

                if (yTick>=minY) {
                    if (yTick<=maxY) {
                        plotThisPair = true;
                        lr.bottom = plotRect.bottom-(yTick-minY)*yScale;
                    } else {
                        lr.bottom = plotRect.bottom-(maxY-minY)*yScale;
                    }
                } else {
                    lr.bottom = plotRect.bottom;
                }

                if (x>=minY) {
                    if (x<=maxY) {
                        plotThisPair = true;
                        lr.top = plotRect.bottom-(x-minY)*yScale;
                    } else {
                        lr.top = plotRect.bottom-(maxY-minY)*yScale;
                    }
                } else {
                    lr.top = plotRect.bottom;
                }



                if (ox2==maxX) {
                    break;
                }

                x = yTick;

                if (plotThisPair) {
                    sc->DrawLine (lr);
                }

                ox1 = ox2;
                ox2+= 1./xScale;
            }
        }
        if ((columns.lData[0]==-1)&&((options==0)||(options==4)||(options==5))) {
            minX-=.5;
            maxX+=.5;
        }
    } else if ((columns.lData[0]==-1)&&((options==0)||(options==4)||(options==5))) {
        minX-=.5;
        maxX+=.5;
    }

    if (options != 6) {
        if (CheckEqual (userMin, userMax)) {
            lr.top = lr.bottom = plotRect.bottom + minY*yScale-1;
        } else {
            lr.top = lr.bottom = plotRect.bottom -1;
        }
    } else {
        lr.top = lr.bottom = plotRect.bottom -1;
    }

    lr.right = plotRect.right;
    lr.left  = plotRect.left;

    sc->SetColor (LongToHYColor(theColors[HY_CHART_COLOR_COORD_LABEL]));
    sc->DrawLine (lr);
    _List       labels;
    _SimpleList labelInfo;

    ComputeHashMarkPlacement (lr,minX,maxX,minY,maxY,labelInfo,labels,4,labelFont1);

    lr.width = 5;

    sc->DrawHashes (lr,labelInfo,labels,labelFont1.size/2+1,1+labelFont1.size/10);

    if (xLabel.sLength+zLabel.sLength+(yLabel.sLength&&(options==5))) {
        sc->SetColor (LongToHYColor(theColors[HY_CHART_COLOR_AXIS_LABEL]));
        sc->SetFont  (labelFont3);

        if (xLabel.sLength) {
            sc->DisplayText (xLabel,plotRect.bottom+labelFont3.size+labelFont1.size, (plotRect.left+plotRect.right-GetVisibleStringWidth(xLabel,labelFont3))/2,true);
        }
        if (zLabel.sLength) {
            sc->DisplayText (zLabel,(5*labelFont3.size)/4,plotRect.left-lS+labelFont3.size/4,true);
        }
        if (yLabel.sLength&&(options==5)) {
            sc->DisplayText (yLabel,plotRect.bottom+labelFont3.size+labelFont1.size,plotRect.left-lS+labelFont3.size/4,true);
        }
    }

}





//__________________________________________________________



void    _HYChartWindow::HandleCellEditEvent (long index)
{

    _HYTable*   table = (_HYTable*)GetObject (1);

    long        h = index%table->horizontalSpaces.lLength,
                v = index/table->horizontalSpaces.lLength;

    bool        good = true;

    _String     *currentValue = (_String*)table->GetCellData (h,v);

    _Formula    f (*currentValue);

    if          (f.IsEmpty()) {
        good = false;
    } else {
        _PMathObj p = (_PMathObj)f.Compute();
        if (p&&(p->ObjectClass()==NUMBER)) {
            currentValue = (_String*)p->toStr();
        } else {
            good = false;
        }
    }

    if (good) {
        table->SetCellData (currentValue,v,h,table->cellTypes.lData[index],false);
    } else {
        table->SetCellData (&empty,v,h,table->cellTypes.lData[index],true);
    }

    table->_MarkCellForUpdate (index);

    if (NeedToRedrawChart(h)) {
        DrawChart();
    }
}

//__________________________________________________________

bool    _HYChartWindow::NeedToRedrawChart (long index)
{
    bool res = false;
    _HYPullDown         *p1 = (_HYPullDown*)GetObject (4);

    long m = p1->GetSelection();

    if (m==2) {
        p1 = (_HYPullDown*)GetObject (5);
        res = (p1->GetSelection()==index);
    } else {
        p1 = (_HYPullDown*)GetObject (5);
        if (m<9) {
            res = (p1->GetSelection()==index);
        }
        if (!res) {
            res = (ySeries.Find (index)>=0)||(index==xAxis3DScale)||(index==yAxis3DScale);
        }
    }

    return res;
}

//__________________________________________________________

void    _HYChartWindow::ToggleSuspend (bool s)
{
    suspendDraw = s;
}

//__________________________________________________________

void    _HYChartWindow::SetProjection (_Parameter a, _Parameter b, _Parameter c)
{
    oR      = a;
    xyAngle = b;
    zAngle  = c;
    ComputeProjectionSettings ();
    projectionMatrix = ComputeProjectionMatrix   ();
}

//__________________________________________________________

void    _HYChartWindow::SetFonts      (_List * data)
{
    _HYFont * fonts2Set [4] = {&labelFont1, &labelFont2, &labelFont3, &headerFont};
    for (long k = 0; k < data->lLength; k++) {
        _List * fontSpec = ((_String*)(*data)(k))->Tokenize (":");
        if (fontSpec->lLength == 3) {
            fonts2Set [k]->face  = *((_String*)(*fontSpec)(0));
            fonts2Set [k]->size  = ((_String*)(*fontSpec)(1))->toNum();
            fonts2Set [k]->style = ((_String*)(*fontSpec)(2))->toNum();
        }
        DeleteObject (fontSpec);
        if (k==3) {
            break;
        }
    }
}

//__________________________________________________________

void    _HYChartWindow::SetColors      (_List * data)
{
    for (long k = 0; k < data->lLength; k++) {
        theColors [k] = ((_String*)(*data)(k))->toNum();

        if (k==theColors.lLength-1) {
            break;
        }
    }
}

//__________________________________________________________

bool    ReadDataFromFile (_String fileName, char delimiter, _Matrix& data, _List& names)
{
    _List readStrings;
    long  columns = 0,
          lastRead;

    FILE* f = doFileOpen (fileName.sData,"rb");
    if (!f) {
        _String errMsg;
        errMsg = _String ('"') & fileName & "\" could not be open for writing.";
        ProblemReport (errMsg);
        return false;
    }

    int  c = fgetc (f);
    while (!feof(f)) {
        lastRead = readStrings.lLength;

        _String  *currentTerm = new _String (16L, true);
        checkPointer (currentTerm);

        while ((c!='\n')&&(c!='\r')&&(c!=-1)) {
            if (c==delimiter) {
                currentTerm->Finalize();
                readStrings << currentTerm;
                currentTerm = new _String (16L, true);
                checkPointer (currentTerm);
            } else {
                (*currentTerm) << (char)c;
            }
            c = fgetc (f);
        }

        currentTerm->Finalize();
        if ((readStrings.lLength>lastRead)||(currentTerm->FirstNonSpaceIndex(0,-1,1)>=0)) {
            readStrings << currentTerm;
        } else {
            DeleteObject (currentTerm);
        }

        c = fgetc (f);
        while ((c=='\n')||(c=='\r')) {
            c = fgetc(f);
        }

        if (lastRead==0) {
            columns = readStrings.lLength;
        }
    }

    fclose (f);

    if (columns&&(readStrings.lLength>=2*columns)&&(readStrings.lLength%columns==0)) {
        data.Clear();
        CreateMatrix(&data,readStrings.lLength/columns-1,columns,false,true,false);

        _Constant h, v;

        for (lastRead = 1; lastRead < readStrings.lLength/columns; lastRead++) {
            h.SetValue (lastRead-1);
            for (long k=0; k<columns; k++) {
                _String * thisString = (_String*)readStrings (lastRead*columns+k);
                if ((thisString->sLength)&&(thisString->FirstNonSpaceIndex (0,-1)>=0)) {
                    _Formula f (*thisString,nil,false);
                    v.SetValue (k);
                    if (!f.IsEmpty()) {
                        data.MStore (&h,&v,f);
                    }
                }
            }
        }

        for (lastRead = 0; lastRead < columns; lastRead++) {
            names << readStrings (lastRead);
        }
    } else {
        _String errMsg ("HyPhy didn't find a well-defined '");
        errMsg = errMsg & delimiter & "'-separated data table in the input file.";
        if (columns == 0) {
            errMsg = errMsg & " Input file was empty...";
        } else if (readStrings.lLength<2*columns) {
            errMsg = errMsg & " Input file contained less than a single valid entry line";
        } else {
            errMsg = errMsg & " The number of fields read (" &(long)readStrings.lLength& ") was not divisible by the number of columns ("& columns &").";
        }

        ProblemReport (errMsg);
        return false;
    }

    return true;
}

//__________________________________________________________

void    _HYChartWindow::DoSave (long option, _String* distrib)
{

    if (option==2) {
        _HYTable*       table = (_HYTable*)    GetObject (1),
                        *       head  = (_HYTable*)    GetObject (2),
                                *       left  = (_HYTable*)    GetObject (3);

        _List           menuOptions;
        _String         menuChoice,
                        filePath,
                        str1 ("Save data table as"),
                        str2;

        _HYTable::GetTableFormats(menuOptions);

        str2 = GetTitle() & " Data Table";

        long menuSel = SaveFileWithPopUp (filePath, str1 ,str2 , empty, menuOptions);

        if (menuSel>=0) {
            FILE*   outFile = doFileOpen (filePath.sData,"w");
            if (!outFile) {
                menuChoice = filePath & " could not be opened for writing.";
                ProblemReport (menuChoice,(Ptr)this);
                return;
            }
            table->SaveTable (head,left,menuSel,outFile,GetTitle());
            fclose (outFile);
        }
    } else {
        if (option==1) {
            _HYCanvas* sc =  (_HYCanvas*)GetObject (0);
            sc->_SavePicture (_String("Chart ")&GetTitle());
        } else {
            _String filePath,
                    cPrompt  ("Save chart as:"),
                    cPrompt2 (GetTitle()),
                    cSaveAs  ("Format:");

            filePath = "HYPHY Batch File";

            _List   fFormats;

            fFormats && & filePath;

            option = SaveFileWithPopUp (filePath,cPrompt,cPrompt2,cSaveAs,fFormats);
            if (option!=-1) {
                FILE * f = doFileOpen (filePath.getStr(), "w");
                if (!f) {
                    filePath = _String ("I couldn't open \"") & filePath & "\" for writing.";
                    ProblemReport (filePath,(Ptr)this);
                }

                fprintf (f,"columnHeaders = {{");
                _HYTable*       table = (_HYTable*)    GetObject (1),
                                *       head  = (_HYTable*)    GetObject (2),
                                        *       tl    = (_HYTable*)    GetObject (3);

                for (option = 0; option < head->horizontalSpaces.lLength; option++) {
                    fprintf (f,"\"%s\",",((_String*)head->GetCellData(option,0))->getStr());
                }


                _HYLabel* l4 = (_HYLabel*)GetObject   (10);
                fprintf (f,"\"%s", l4->GetText().getStr());

                for (option = 0; option < tl->verticalSpaces.lLength-1; option++) {
                    fprintf (f,";%s",((_String*)tl->GetCellData(0,option))->getStr());
                }

                fprintf (f,"\"}};\ntableData = {\n");

                for (option = 0; option < table->verticalSpaces.lLength; option++) {
                    fprintf (f,"{");
                    fprintf (f,"%s",((_String*)table->GetCellData(0,option))->getStr());

                    for (long option2 = 1; option2 < table->horizontalSpaces.lLength; option2++) {
                        fprintf (f,",%s",((_String*)table->GetCellData(option2,option))->getStr());
                    }
                    fprintf (f,"}\n");
                }

                _HYRect wr = GetWindowRect ();

                _HYPullDown         *p1 = (_HYPullDown*)GetObject (4);
                _HYPullDown         *p2 = (_HYPullDown*)GetObject (5);

                _String             windowType = distrib?windowTypeDistribTable:windowTypeTable;

                fprintf (f,"};\nOpenWindow (%s,{{\"%s\"}\n\t\t{\"columnHeaders\"}\n\t\t{\"tableData\"}\n\t\t{\"%s\"}\n\t\t{\"%s\"}\n\t\t{\"",
                         windowType.getStr(),
                         GetTitle().getStr(),
                         p1->GetMenuItem(p1->GetSelection())->getStr(),
                         p2->GetMenuItem(p2->GetSelection())->getStr());

                if (ySeries.lLength==0) {
                    fprintf (f,"None");
                } else {
                    fprintf (f,"%s",((_String*)head->GetCellData(ySeries.lData[0],0))->getStr());
                    for (option = 1; option < ySeries.lLength; option++) {
                        fprintf (f,";%s",((_String*)head->GetCellData(ySeries.lData[option],0))->getStr());
                    }
                }

                _String * fla = (_String*)overlayPlot.toStr();

                fprintf (f,
                         "\"}\n\t\t{\"%s\"}\n\t\t{\"%s\"}\n\t\t{\"%s\"}\n\t\t{\"%d\"}\n\t\t{\"%s\"}\n\t\t{\"%ld;%ld\"}\n\t\t{\"%g;%g;%g\"}\n\t\t{\"%s:%d:%d;%s:%d:%d;%s:%d:%d\"}\n\t\t{\"",
                         xLabel.getStr(),
                         yLabel.getStr(),
                         zLabel.getStr(),
                         (int)showLegend,
                         fla->getStr(),
                         xAxis3DScale+1,
                         yAxis3DScale+1,
                         oR,
                         xyAngle,
                         zAngle,
                         labelFont1.face.sData,
                         labelFont1.size,
                         labelFont1.style,
                         labelFont2.face.sData,
                         labelFont2.size,
                         labelFont2.style,
                         labelFont3.face.sData,
                         labelFont3.size,
                         labelFont3.style);

                fprintf (f, "%ld", theColors.lData[0]);
                for (option = 1; option<theColors.lLength; option++) {
                    fprintf (f, ";%ld", theColors.lData[option]);
                }


                if (distrib)
                    fprintf (f,"\"}\n\t\t{\"%ld,%g,%g\"}\n\t\t{\"%s\"}\n\t\t}",
                             surfaceDivs,userMin,userMax,
                             distrib->getStr());

                else
                    fprintf (f,"\"}\n\t\t{\"%ld,%g,%g\"}\n\t\t}",
                             surfaceDivs,userMin,userMax);

                fprintf (f,",\n\t\t\"%ld;%ld;%ld;%ld\");",
                         wr.right-wr.left,
                         wr.bottom-wr.top,
                         wr.left,
                         wr.top);

                DeleteObject (fla);
                fclose (f);
            }
        }
    }
}

//__________________________________________________________

void    _HYChartWindow::DoPrint (long option)
{
    if (option) {
        _HYTable*       table = (_HYTable*)    GetObject (1),
                        *       head  = (_HYTable*)    GetObject (2);

        table->_PrintTable (head);
    } else {
        _PrintChart();
    }

}

//__________________________________________________________

_Matrix* _HYChartWindow::ComputeProjectionMatrix (void)
{
    if (projectionMatrix) {
        DeleteObject (projectionMatrix);
        projectionMatrix = nil;
    }

    _Matrix      *res = new _Matrix (4,4,false,true);
    checkPointer (res);

    _Parameter   den = 1+projectionSettings[5],
                 a = projectionSettings[0],
                 b = projectionSettings[1],
                 c = projectionSettings[2],
                 d = projectionSettings[3],
                 e = projectionSettings[4],
                 f = projectionSettings[5],
                 r = 1./projectionSettings[6];

    res->theData[0] = (e*e+f*(1+f))/den;
    res->theData[1] = res->theData[4] = -d*e/den;
    res->theData[2] = res->theData[6] = res->theData[10] = res->theData[14];
    res->theData[3] = d*r;
    res->theData[5] = (d*d+f*(1+f))/den;
    res->theData[7] = e*r;
    res->theData[8] = -d;
    res->theData[9] = -e;
    res->theData[11] = f*r;
    res->theData[12] = (c*d+b*d*e-a*e*e-a*f+c*d*f-a*f*f)/den;
    res->theData[13] = (-b*d*d+c*e+a*d*e-b*f+c*e*f-b*f*f)/den;
    res->theData[15] = -(a*d+b*e+c*f)*r;

    return res;
}

//__________________________________________________________

void    _HYChartWindow::DrawAParallelogram  (_Matrix& pData ,char  ,_Parameter ,_HYRect& cRect)

/* order of points:

    x direction -> left
    y direction -> up
    z direction -> away

    3------4 (farther)
    |      |
    |      | (top face) (x direction ---->)
    |      |
    1------2 (closer)

    7------8 (farther)
    |      |
    |      | (bottom face)
    |      |
    5------6 (closer)

    order of painting the faces:

        3487 (back)
        5687 (bottom)
        1375 (left)
        2486 (right)
        1243 (top)
        5621 (front)
*/

{

    _Matrix    *projP = ComputeProjection (pData);

    PaintAFace (*projP,2,3,7,6,cRect);
    PaintAFace (*projP,4,5,7,6,cRect);
    PaintAFace (*projP,0,2,6,4,cRect);
    PaintAFace (*projP,1,3,7,5,cRect,.7);
    PaintAFace (*projP,0,1,3,2,cRect,.8);
    PaintAFace (*projP,4,5,1,0,cRect);

    DeleteObject (projP);
}

//__________________________________________________________

_Matrix*    _HYChartWindow::ComputeProjection   (_Matrix& pData)
{
    _Matrix             hompc (pData.GetHDim(),4,false,true),
                        *projP;

    long                counter;

    for (counter = 0; counter < hompc.GetHDim(); counter++) {
        hompc.theData[counter*4]   = pData.theData[counter*3];
        hompc.theData[counter*4+1] = pData.theData[counter*3+1];
        hompc.theData[counter*4+2] = pData.theData[counter*3+2];
        hompc.theData[counter*4+3] = 1;
    }

    projP = (_Matrix*)hompc.MultObj (projectionMatrix);

    for (counter = 0; counter < hompc.GetHDim(); counter++) {
        _Parameter scaling = 1./projP->theData[counter*4+3];

        projP->theData[counter*4]   *= scaling;
        projP->theData[counter*4+1] *= scaling;
        projP->theData[counter*4+2] *= scaling;
    }

    return projP;
}

//__________________________________________________________

void    _HYChartWindow::DrawALine           (_Matrix& pData ,char  ,_Parameter ,_HYRect& cRect, long width,
        _String* label, long xSpace, long ySpace)

{
    _HYStretchCanvas    *sc = (_HYStretchCanvas*)GetObject (0);
    _Matrix             *projP = ComputeProjection (pData);
    _HYRect             lineRect;

    lineRect.left       = cRect.left   - xShift + xScale*(*projP)(0,0);
    lineRect.right      = cRect.left   - xShift + xScale*(*projP)(1,0);
    lineRect.top        = cRect.bottom + yShift - yScale*(*projP)(0,1);
    lineRect.bottom     = cRect.bottom + yShift - yScale*(*projP)(1,1);
    lineRect.width      = width;

    sc->DrawLine (lineRect);

    if (label) {
        _HYColor       bc = sc->GetColor();

        sc->SetColor    (LongToHYColor (theColors[HY_CHART_COLOR_COORD_LABEL]));
        sc->DisplayText (*label, lineRect.top+ySpace, lineRect.left+xSpace, true);
        sc->SetColor    (bc);
    }
}

//__________________________________________________________

void    _HYChartWindow::DrawAPlane          (_Matrix& pData ,char opt ,_Parameter ,_HYRect& cRect)
/* specify the plane (rectangle) as follows

    1----2
    |    |
    |    |
    3----4
*/

{
    _Matrix             *projP = ComputeProjection (pData);
    PaintAFace          (*projP,0,1,3,2,cRect,1., opt);
    DeleteObject        (projP);
}

//__________________________________________________________

void    _HYChartWindow::PaintAFace          (_Matrix& pData ,long v1, long v2, long v3, long v4,_HYRect& cRect,_Parameter shading, char opt)
{
    _HYStretchCanvas    *sc = (_HYStretchCanvas*)GetObject (0);

    _SimpleList         points;

    points << cRect.left   - xShift + xScale*pData(v1,0);
    points << cRect.bottom + yShift - yScale*pData(v1,1);
    points << cRect.left   - xShift + xScale* pData(v2,0);
    points << cRect.bottom + yShift - yScale* pData(v2,1);
    points << cRect.left   - xShift + xScale* pData(v3,0);
    points << cRect.bottom + yShift - yScale* pData(v3,1);
    points << cRect.left   - xShift + xScale* pData(v4,0);
    points << cRect.bottom + yShift - yScale* pData(v4,1);

    _HYColor       bc = sc->GetColor(),
                   nc = sc->GetColor();

    if (opt == 0) {
        if (shading!=1.0) {
            nc.R   *= shading;
            nc.B   *= shading;
            nc.G   *= shading;
            sc->SetColor (nc);
        }
        sc->FillPolygon  (points);
    }

    if (opt != 2) {
        sc->SetColor (LongToHYColor(theColors[HY_CHART_COLOR_BARBORDER]));
    } else {
        sc->SetColor (LongToHYColor(theColors[HY_CHART_COLOR_OVERPLOT]));
    }

    sc->DrawPolygon  (points,1);
    sc->SetColor (bc);

}

//__________________________________________________________

void    _HYChartWindow::SetFont     (_HYFont& newFont)
{
    if (whichFont == 0) {
        labelFont1 = newFont;
        DrawChart ();
    } else {
        if (whichFont == 1) {
            labelFont2 = newFont;
            if (showLegend != HY_CHART_LEGEND_NONE) {
                DrawChart();
            }
        } else {
            labelFont3 = newFont;
            if (xLabel.sLength+yLabel.sLength+zLabel.sLength) {
                DrawChart();
            }
        }
    }
}

//__________________________________________________________

void    _HYChartWindow::DoChangeFont        (long item)
{
    _HYFont *inFont;
    switch (item) {
    case 0:
        inFont = &labelFont1;
        break;
    case 1:
        inFont = &labelFont2;
        break;
    case 2:
        inFont = &labelFont3;
        break;
    }
    _HYFontDialog * fD = new _HYFontDialog (*inFont,this);
    checkPointer  (fD);

    fD->BringToFront();

    whichFont = item;

    while (windowObjectRefs.Find ((long)fD)>=0) {
        handleGUI();
    }
}

//__________________________________________________________

void    _HYChartWindow::HandleChartOptions      (void)
{
    _List       args;

    _String*    overlayString = (_String*)overlayPlot.toStr();

    args && & xLabel;
    args && & yLabel;
    args && & zLabel;
    args << overlayString;
    args.AppendNewInstance(new _String (xAxis3DScale));
    args.AppendNewInstance(new _String (yAxis3DScale));
    args.AppendNewInstance(new _String (surfaceDivs));
    if (!CheckEqual (userMin, userMax)) {
        args.AppendNewInstance(new _String (userMin));
        args.AppendNewInstance(new _String (userMax));
    } else {
        args && & empty;
        args && & empty;
    }

    long        lsel = showLegend;
    bool        res  = false;
    _SimpleList clst;

    clst.Duplicate (&theColors);

    _HYChartOptionDialog * fD = new _HYChartOptionDialog (&args,&lsel,&res,&clst, this);
    checkPointer  (fD);

    fD->BringToFront();

    while (windowObjectRefs.Find ((long)fD)>=0) {
        handleGUI();
    }

    if (res) {
        res = false;

        if (!((_String*)args(0))->Equal(&xLabel)) {
            res = true;
            xLabel = *(_String*)args(0);
        }
        if (!((_String*)args(1))->Equal(&yLabel)) {
            res = true;
            yLabel = *(_String*)args(1);
        }
        if (!((_String*)args(2))->Equal(&zLabel)) {
            res = true;
            zLabel = *(_String*)args(2);
        }
        if (lsel!=showLegend) {
            res = true;
            showLegend = lsel;
        }

        if (!((_String*)args(3))->Equal(overlayString)) {
            DeleteObject (overlayString);
            overlayString = (_String*)args(3);
            res = true;
            overlayPlot.Clear();



            if (overlayString->sLength) {
                long varRef = 0;
                _String errMsg;
              _FormulaParsingContext fpc(&errMsg, nil);
		if (Parse (&overlayPlot, *(_String*)args(3),fpc,nil)!= HY_FORMULA_EXPRESSION) {
                    overlayPlot.Clear();
                }
            }
        } else {
            DeleteObject (overlayString);
        }

        lsel = ((_String*)args(4))->toNum();

        if (xAxis3DScale!=lsel) {
            xAxis3DScale = lsel;
        }

        lsel = ((_String*)args(5))->toNum();

        if (yAxis3DScale!=lsel) {
            yAxis3DScale = lsel;
        }

        lsel = ((_String*)args(6))->toNum();

        if ((surfaceDivs != lsel)&&(lsel>1)) {
            surfaceDivs = lsel;
            res = true;
        }

        res |= SetUserBounds (((_String*)args(7))->toNum(), ((_String*)args(8))->toNum());

        for (lsel = 0; lsel < clst.lLength; lsel++)
            if (clst.lData[lsel] != theColors.lData[lsel]) {
                theColors.Duplicate (&clst);
                res = true;
                break;
            }

        if (res) {
            DrawChart();
        }
    }

}

//__________________________________________________________

void    _HYChartWindow::SetLabels       (_String xx, _String yy, _String zz, long vv, _String fla, long xs, long ys, long sDiv, _Parameter uMin, _Parameter uMax)
{
    xLabel      =   xx;
    yLabel      =   yy;
    zLabel      =   zz;
    showLegend  =   vv;
    overlayPlot.Clear();

    if (fla.sLength) {
        long varRef = 0;
        _String errMsg;
      _FormulaParsingContext fpc (&errMsg);
        if (Parse (&overlayPlot, fla,fpc, nil)!= HY_FORMULA_EXPRESSION) {
            overlayPlot.Clear();
        }
    }

    xAxis3DScale = xs;
    yAxis3DScale = ys;

    if (sDiv>1) {
        surfaceDivs = sDiv;
    }

    SetUserBounds (uMin, uMax);

    DrawChart   ();
}

//__________________________________________________________

long    _HYChartWindow::Get3DScaling        (bool isY)
{
    return isY?yAxis3DScale:xAxis3DScale;
}

//__________________________________________________________________
void        _HYChartWindow::HandleCopyPaste (bool paste)
{
    if (!paste) {
        _HYTable*               table = (_HYTable*)    GetObject (1);
        _SimpleList             tSel;

        table->GetSelection     (tSel);

        if (tSel.lLength == 0) {
            _CopyChart ();
            return;
        }
    }

    _HYTWindow::HandleCopyPaste (paste);
}

//__________________________________________________________________
void        _HYChartWindow::RenameChartWindow (void)
{
    _String currentName = GetTitle(),
            prompt      = "New Chart Name:",
            newName     = currentName;

    if (EnterStringDialog (newName, prompt, (Ptr)this)) {
        if (!newName.Equal (&currentName)) {
            SetTitle (newName);
        }
    }
}

//__________________________________________________________

void    _HYChartWindow::ExecuteProcessor (long cID)
{
    if (cID<chartProcessors.lLength) {
        _HYTable*               table   = (_HYTable*)      GetObject (1);
        _HYTable*               tableH = (_HYTable*)       GetObject (2);
        _SimpleList             tableSel;

        table->GetSelection     (tableSel);

        if (tableSel.lLength > 0) {
            //_String         eM ("Please select some entries in the data table before calling result processing modules.");
            //ProblemReport (eM);
            //return;
            //}
            //else
            //{
            _Matrix selData(1,tableSel.lLength,false,true),
                    selRows(1,tableSel.lLength,false,true),
                    selCols(1,tableSel.lLength,false,true);

            for (long k=0; k<tableSel.lLength; k++) {
                long r = tableSel.lData[k]/table->horizontalSpaces.lLength,
                     c = tableSel.lData[k]%table->horizontalSpaces.lLength;

                selData.theData[k] = ((_String*)table->GetCellData(c,r))->toNum();
                selRows.theData[k] = r;
                selCols.theData[k] = c;
            }

            setParameter (chartProcessorDataMatrix,&selData);
            setParameter (chartProcessorRowMatrix, &selRows);
            setParameter (chartProcessorColMatrix, &selCols);
        }

        _List            colNList;

        for (long cc = 0; cc < tableH->horizontalSpaces.lLength; cc++) {
            colNList << tableH->GetCellData (cc, 0);
        }

        _Matrix       colNames (colNList);
        setParameter (chartColumnHeaders, &colNames);
        setParameter (nonEmptyChartSelection, tableSel.lLength);

        _Matrix tableDims (1,2,false,true);
        tableDims.theData[0] = table->verticalSpaces.lLength;
        tableDims.theData[1] = table->horizontalSpaces.lLength;
        setParameter (tableDimensions, &tableDims);

        FILE * thisFile = doFileOpen (((_String*)chartProcessors(cID))->sData,"rb");
        if (thisFile) {
            _String buffer (thisFile);
            fclose (thisFile);

            long   g = batchLanguageFunctionNames.lLength;

            setParameter (updateCellData,0.0);
            setParameter (updateTableSelection,0.0);

            _String cpp = *((_String*)chartProcessors(cID));

            PushFilePath (cpp);
            _ExecutionList   thisList;
            terminateExecution = false;
            thisList.BuildList (buffer);
            thisList.ExecuteAndClean(g);
            terminateExecution = false;

            PopFilePath ();

            _Parameter setValues = 0.0;

            checkParameter (nonEmptyChartSelection, setValues, 0.0);

            if (setValues < -.1) {
                _String       eM ("Please select some entries in the data table before calling result processing modules.");
                ProblemReport (eM);
            } else {
                checkParameter (updateCellData, setValues, 0.0);

                if (setValues>0.1) {
                    _Variable * rM = FetchVar(LocateVarByName (chartProcessorRowMatrix)),
                                * cM = FetchVar(LocateVarByName (chartProcessorColMatrix)),
                                  * cV = FetchVar(LocateVarByName (chartProcessorDataMatrix));

                    if ((rM->ObjectClass()==MATRIX)&&(cM->ObjectClass()==MATRIX)&&(cV->ObjectClass()==MATRIX)) {
                        _Matrix   * rowIndices      = (_Matrix*)rM->Compute(),
                                    * colIndices     = (_Matrix*)cM->Compute(),
                                      * cellValues      = (_Matrix*)cV->Compute();

                        if ((rowIndices->GetHDim()==1)&&(colIndices->GetHDim()==1)&&(cellValues->GetHDim()==1)
                                &&(rowIndices->GetVDim()==colIndices->GetVDim())&&(cellValues->GetVDim()==colIndices->GetVDim())) {
                            _SimpleList   toUpdate;

                            for (long jj = 0; jj < cellValues->GetVDim(); jj++) {
                                _String newVal ((*cellValues)(0,jj));

                                long    cI = (*colIndices)(0,jj),
                                        rI = (*rowIndices)(0,jj);

                                if ((rI<table->verticalSpaces.lLength)&&(cI<table->horizontalSpaces.lLength)) {
                                    long idx = rI*table->horizontalSpaces.lLength + cI;
                                    table->SetCellData (&newVal,rI,cI,table->cellTypes.lData[idx],true);
                                    toUpdate << idx;
                                }
                            }

                            table->_MarkCellsForUpdate (toUpdate);
                            DrawChart ();

                        } else {
                            _String       eM = _String ("Can't update cell values because the three matrices: ") & chartProcessorRowMatrix &
                                               ',' & chartProcessorColMatrix & " or " & chartProcessorDataMatrix & " are not column vector of the same dimension.";

                            ProblemReport (eM);
                            return;
                        }
                    } else {
                        _String       eM = _String ("Can't update cell values because one of the three variables: ") & chartProcessorRowMatrix &
                                           ',' & chartProcessorColMatrix & " or " & chartProcessorDataMatrix & " is not a matrix, as expected";
                        ProblemReport (eM);
                        return;
                    }
                } else {
                    checkParameter (updateTableSelection, setValues, 0.0);
                    if (setValues > 0.1) {
                        _Variable * cM = FetchVar(LocateVarByName (chartProcessorColMatrix)),
                                    * rM = FetchVar(LocateVarByName (chartProcessorRowMatrix));

                        if ((rM->ObjectClass()==MATRIX)&&(cM->ObjectClass()==MATRIX)) {
                            _Matrix   * rowIndices      = (_Matrix*)rM->Compute(),
                                        * colIndices     = (_Matrix*)cM->Compute();

                            if ((rowIndices->GetHDim()==1)&&(colIndices->GetHDim()==1)&&(rowIndices->GetVDim()==colIndices->GetVDim())) {
                                _SimpleList   toUpdate;
                                for (long jj = 0; jj < rowIndices->GetVDim(); jj++) {
                                    toUpdate << rowIndices->theData[jj] * table->horizontalSpaces.lLength+
                                             colIndices->theData[jj];
                                }
                                table->SetSelection (toUpdate,true);
                                table->_MarkForUpdate();
                            } else {
                                _String       eM = _String ("Can't update cell values because the two matrices: ") & chartProcessorRowMatrix &
                                                   ',' & chartProcessorColMatrix & " are not column vector of the same dimension.";

                                ProblemReport (eM);
                                return;
                            }
                        } else {
                            _String       eM = _String ("Can't update cell values because one of the two variables: ") & chartProcessorRowMatrix &
                                               ',' & chartProcessorColMatrix & " is not a matrix, as expected";
                            ProblemReport (eM);
                            return;
                        }
                    }
                }
            }
        } else {
            _String eMsg = _String("Problem reading file:") & *((_String*)chartProcessors(cID));
            ProblemReport (eMsg);
        }

        DeleteVariable (chartProcessorRowMatrix);
        DeleteVariable (chartProcessorColMatrix);
        DeleteVariable (chartProcessorDataMatrix);
        DeleteVariable (tableDimensions);
        DeleteVariable (chartColumnHeaders);
    }
}

//__________________________________________________________

_HYChartOptionDialog::_HYChartOptionDialog (_List* argList, long* legOpt, bool* retVal, _SimpleList* clst,_HYChartWindow* pw):_HYTWindow ("HyPhy Chart Preferences", false, true, (Ptr)pw)
{
    _HYRect         canvasSettings = {30,100,30,100,HY_COMPONENT_NO_SCROLL|HY_COMPONENT_TRANSP_BG};

    args                = argList;
    msel                = legOpt;
    res                 = retVal;
    parWin              = pw;
    colorsPicked        = clst;

    _HYLabel*       l1      = new _HYLabel (canvasSettings, GetOSWindowData());
    _HYLabel*       l2      = new _HYLabel (canvasSettings, GetOSWindowData());
    _HYLabel*       l3      = new _HYLabel (canvasSettings, GetOSWindowData());
    _HYLabel*       l4      = new _HYLabel (canvasSettings, GetOSWindowData());
    _HYLabel*       l5      = new _HYLabel (canvasSettings, GetOSWindowData());
    _HYLabel*       l6      = new _HYLabel (canvasSettings, GetOSWindowData());
    _HYLabel*       l7      = new _HYLabel (canvasSettings, GetOSWindowData());
    _HYLabel*       l8      = new _HYLabel (canvasSettings, GetOSWindowData());
    _HYLabel*       l9      = new _HYLabel (canvasSettings, GetOSWindowData());
    _HYLabel*       l10     = new _HYLabel (canvasSettings, GetOSWindowData());

    checkPointer   (l1);
    checkPointer   (l2);
    checkPointer   (l3);
    checkPointer   (l4);
    checkPointer   (l5);
    checkPointer   (l6);
    checkPointer   (l7);
    checkPointer   (l8);
    checkPointer   (l9);
    checkPointer   (l10);

    _HYColor        bgc     = GetDialogBackgroundColor();

    canvasSettings.left = canvasSettings.right = 200;

    _HYPullDown*    p1      = new _HYPullDown (canvasSettings, GetOSWindowData());
    _HYPullDown*    p2      = new _HYPullDown (canvasSettings, GetOSWindowData());
    _HYPullDown*    p3      = new _HYPullDown (canvasSettings, GetOSWindowData());
    _HYPullDown*    p4      = new _HYPullDown (canvasSettings, GetOSWindowData());

    checkPointer   (p1);
    checkPointer   (p2);
    checkPointer   (p3);
    checkPointer   (p4);

    _HYTextBox *    t1      = new _HYTextBox (canvasSettings, GetOSWindowData());
    _HYTextBox *    t2      = new _HYTextBox (canvasSettings, GetOSWindowData());
    _HYTextBox *    t3      = new _HYTextBox (canvasSettings, GetOSWindowData());
    _HYTextBox *    t4      = new _HYTextBox (canvasSettings, GetOSWindowData());
    _HYTextBox *    t5      = new _HYTextBox (canvasSettings, GetOSWindowData());

    checkPointer    (t1);
    checkPointer    (t2);
    checkPointer    (t3);
    checkPointer    (t4);
    checkPointer    (t5);

    canvasSettings.left = canvasSettings.right = 100;

    _HYTextBox *    t6      = new _HYTextBox (canvasSettings, GetOSWindowData());
    _HYTextBox *    t7      = new _HYTextBox (canvasSettings, GetOSWindowData());

    checkPointer    (t6);
    checkPointer    (t7);

    canvasSettings.left = canvasSettings.right = 101;
    _HYButton*      b1      = new _HYButton (canvasSettings, GetOSWindowData());

    canvasSettings.left = canvasSettings.right = 101;
    _HYButton*      b2      = new _HYButton (canvasSettings, GetOSWindowData());

    canvasSettings.left = canvasSettings.right = 100;
    _HYButton*      b3      = new _HYButton (canvasSettings, GetOSWindowData());

    canvasSettings.left = canvasSettings.right = 300;
    _HYCanvas*      color   = new _HYCanvas (canvasSettings, GetOSWindowData(), canvasSettings.top, canvasSettings.left, 32);
    checkPointer    (color);

    b1->SetMessageRecipient (this);
    b2->SetMessageRecipient (this);
    b3->SetMessageRecipient (this);
    t1->SetMessageRecipient (this);
    t2->SetMessageRecipient (this);
    t3->SetMessageRecipient (this);
    t4->SetMessageRecipient (this);
    t5->SetMessageRecipient (this);
    t6->SetMessageRecipient (this);
    t7->SetMessageRecipient (this);

    p4->SetMessageRecipient (this);
    color->SetMessageRecipient (this);
    color->SetMouseClick       (true);

    AddObject (l1); // 0
    AddObject (l2); // 1
    AddObject (l3); // 2
    AddObject (l4); // 3
    AddObject (l5); // 4

    AddObject (p1); // 5

    AddObject (t1); // 6
    AddObject (t2); // 7
    AddObject (t3); // 8
    AddObject (t4); // 9

    AddObject (b1); // 10
    AddObject (b2); // 11
    AddObject (b3); // 12

    AddObject (l6); // 13
    AddObject (l7); // 14

    AddObject (p2); // 15
    AddObject (p3); // 16

    AddObject (l8); // 17
    AddObject (p4); // 18
    AddObject (color); // 19

    AddObject (l9); // 20
    AddObject (t5); // 21

    AddObject (l10); // 22
    AddObject (t6); // 23
    AddObject (t7); // 24




    keyboardFocusChain << 6;
    keyboardFocusChain << 7;
    keyboardFocusChain << 8;
    keyboardFocusChain << 9;
    keyboardFocusChain << 21;
    keyboardFocusChain << 23;
    keyboardFocusChain << 24;

    SetTableDimensions (12,4);

    SetCell   (0,0,l1);
    SetCell   (0,1,t1);
    SetCell   (0,2,t1);
    SetCell   (0,3,t1);

    SetCell   (1,0,l2);
    SetCell   (1,1,t2);
    SetCell   (1,2,t2);
    SetCell   (1,3,t2);

    SetCell   (2,0,l3);
    SetCell   (2,1,t3);
    SetCell   (2,2,t3);
    SetCell   (2,3,t3);

    SetCell   (3,0,l5);
    SetCell   (3,1,t4);
    SetCell   (3,2,t4);
    SetCell   (3,3,t4);

    SetCell   (4,0,l9);
    SetCell   (4,1,t5);
    SetCell   (4,2,t5);
    SetCell   (4,3,t5);

    SetCell   (5,0,l10);
    SetCell   (5,1,t6);
    SetCell   (5,2,t7);
    SetCell   (5,3,t7);

    SetCell   (6,0,l4);
    SetCell   (6,1,p1);
    SetCell   (6,2,p1);
    SetCell   (6,3,p1);

    SetCell   (7,0,l6);
    SetCell   (7,1,p2);
    SetCell   (7,2,p2);
    SetCell   (7,3,p2);

    SetCell   (8,0,l7);
    SetCell   (8,1,p3);
    SetCell   (8,2,p3);
    SetCell   (8,3,p3);

    SetCell   (9,0,l8);
    SetCell   (9,1,p4);
    SetCell   (9,2,p4);
    SetCell   (9,3,p4);

    SetCell   (10,0,color);
    SetCell   (10,1,color);
    SetCell   (10,2,color);
    SetCell   (10,3,color);

    SetCell   (11,0,b3);
    SetCell   (11,1,b1);
    SetCell   (11,2,b2);
    SetCell   (11,3,b2);


    _HYFont  labelFont;
    l1->SetBackColor (bgc);
    l2->SetBackColor (bgc);
    l3->SetBackColor (bgc);
    l4->SetBackColor (bgc);
    l5->SetBackColor (bgc);
    l6->SetBackColor (bgc);
    l7->SetBackColor (bgc);
    l8->SetBackColor (bgc);
    l9->SetBackColor (bgc);
    l10->SetBackColor (bgc);

    p1->SetBackColor (bgc);

    t1->SetBackColor (bgc);
    t2->SetBackColor (bgc);
    t3->SetBackColor (bgc);
    t4->SetBackColor (bgc);
    t5->SetBackColor (bgc);
    t6->SetBackColor (bgc);
    t7->SetBackColor (bgc);

    b1->SetBackColor (bgc);
    b2->SetBackColor (bgc);
    b3->SetBackColor (bgc);


    labelFont.face  = "System Font";
    labelFont.size  = 12;
    labelFont.style = HY_FONT_PLAIN;

    l1->SetFont (labelFont);
    l2->SetFont (labelFont);
    l3->SetFont (labelFont);
    l4->SetFont (labelFont);
    l5->SetFont (labelFont);
    l6->SetFont (labelFont);
    l7->SetFont (labelFont);
    l8->SetFont (labelFont);
    l9->SetFont (labelFont);
    l10->SetFont (labelFont);

    b1->SetFont (labelFont);
    b2->SetFont (labelFont);
    b3->SetFont (labelFont);

    labelFont.face  = "Times";

    t1->SetFont (labelFont);
    t2->SetFont (labelFont);
    t3->SetFont (labelFont);
    t4->SetFont (labelFont);
    t5->SetFont (labelFont);
    t6->SetFont (labelFont);
    t7->SetFont (labelFont);

    color->StartDraw();
    color->SetFont (labelFont);
    color->SetDialogBG ();
    color->EndDraw();

    l1->SetAlignFlags (HY_ALIGN_LEFT);
    l2->SetAlignFlags (HY_ALIGN_LEFT);
    l3->SetAlignFlags (HY_ALIGN_LEFT);
    l4->SetAlignFlags (HY_ALIGN_LEFT);
    l5->SetAlignFlags (HY_ALIGN_LEFT);
    l6->SetAlignFlags (HY_ALIGN_LEFT);
    l7->SetAlignFlags (HY_ALIGN_LEFT);
    l8->SetAlignFlags (HY_ALIGN_LEFT);
    l9->SetAlignFlags (HY_ALIGN_LEFT);
    l10->SetAlignFlags (HY_ALIGN_LEFT);

    t1->SetAlignFlags (HY_ALIGN_LEFT);
    t2->SetAlignFlags (HY_ALIGN_LEFT);
    t3->SetAlignFlags (HY_ALIGN_LEFT);
    t4->SetAlignFlags (HY_ALIGN_LEFT);
    t5->SetAlignFlags (HY_ALIGN_LEFT);
    t6->SetAlignFlags (HY_ALIGN_LEFT);
    t7->SetAlignFlags (HY_ALIGN_LEFT);

    //b3->SetAlignFlags (HY_ALIGN_LEFT);
    b1->SetAlignFlags (HY_ALIGN_RIGHT);

    p1->SetAlignFlags (HY_ALIGN_LEFT);
    p2->SetAlignFlags (HY_ALIGN_LEFT);
    p3->SetAlignFlags (HY_ALIGN_LEFT);
    p4->SetAlignFlags (HY_ALIGN_LEFT);

    l1->SetText       ("X Axis:");
    l2->SetText       ("Y Axis:");
    l3->SetText       ("Value Axis:");

    l4->SetText       ("Legend:");
    l5->SetText       ("Overplot f(_x_)=");

    l6->SetText       ("3D X scaling:");
    l7->SetText       ("3D Y scaling:");
    l8->SetText       ("Colors:");

    l9->SetText       ("Surface Divisions:");

    l10->SetText      ("Y Range:");

    b1->SetText       ("  OK  ");
    b1->SetButtonKind (HY_BUTTON_OK);
    b2->SetText       (" Cancel ");
    b2->SetButtonKind (HY_BUTTON_CANCEL);
    b3->SetText       (" Apply ");


    p1->AddMenuItem   ("None",-1);
    p1->AddMenuItem   ("Top Left",-1);
    p1->AddMenuItem   ("Top Middle",-1);
    p1->AddMenuItem   ("Top Right",-1);
    p1->AddMenuItem   ("Bottom Left",-1);
    p1->AddMenuItem   ("Bottom Middle",-1);
    p1->AddMenuItem   ("Bottom Right",-1);
    p1->AddMenuItem   ("Middle Left",-1);
    p1->AddMenuItem   ("Middle Right",-1);
    p1->ChangeSelection
    (*msel,false);

    p2->AddMenuItem   ("None",-1);
    p2->AddMenuItem   (menuSeparator,-1);
    p3->AddMenuItem   ("None",-1);
    p3->AddMenuItem   (menuSeparator,-1);

    _HYTable*         thead = (_HYTable*)pw->GetObject (2);

    for (long k=0; k<thead->horizontalSpaces.lLength; k++) {
        _String* tentry = (_String*)thead->GetCellData (k,0);
        p2->AddMenuItem (*tentry,-1);
        p3->AddMenuItem (*tentry,-1);
    }

    if (pw->Get3DScaling (false)>=0) {
        p2->ChangeSelection (pw->Get3DScaling (false)+2, false);
    }

    if (pw->Get3DScaling (true)>=0) {
        p3->ChangeSelection (pw->Get3DScaling (true)+2, false);
    }

    t1->SetText       (*(_String*)(*args)(0),false);
    t2->SetText       (*(_String*)(*args)(1),false);
    t3->SetText       (*(_String*)(*args)(2),false);
    t4->SetText       (*(_String*)(*args)(3),false);
    t5->SetText       (*(_String*)(*args)(6),false);
    t6->SetText       (*(_String*)(*args)(7),false);
    t7->SetText       (*(_String*)(*args)(8),false);

    for (long k=0; k<HY_CHART_COLOR_OFFSET; k++) {
        p4->AddMenuItem (chartColorLabels[k], -1);
    }

    for (long k=1; k<=HY_CHART_COLOR_COUNT; k++) {
        _String lbl = "Series ";
        lbl = lbl & k & " color";
        p4->AddMenuItem (lbl, -1);
    }

    DrawColorCanvas ();

    DeleteObject    (l1);
    DeleteObject    (l2);
    DeleteObject    (l3);
    DeleteObject    (l4);
    DeleteObject    (l5);
    DeleteObject    (l6);
    DeleteObject    (l7);
    DeleteObject    (l8);
    DeleteObject    (l9);
    DeleteObject    (l10);

    DeleteObject    (p1);
    DeleteObject    (p2);
    DeleteObject    (p3);
    DeleteObject    (p4);

    DeleteObject    (b1);
    DeleteObject    (b2);
    DeleteObject    (t1);
    DeleteObject    (t2);
    DeleteObject    (t3);
    DeleteObject    (t4);
    DeleteObject    (t5);
    DeleteObject    (t6);
    DeleteObject    (t7);

    DeleteObject    (color);

    _HYRect     dim = MinMaxWindowDimensions();

    SetWindowRectangle  (0,0,dim.bottom,dim.right);
    CenterWindow        (this);
}
//__________________________________________________________

void    _HYChartOptionDialog::DrawColorCanvas (void)
{
    _HYCanvas  * cnvs = (_HYCanvas*)   GetObject (19);
    _HYPullDown* pd   = (_HYPullDown*) GetObject (18);

    _HYColor    clr = LongToHYColor ((*colorsPicked)[pd->GetSelection()]);

    _HYRect     colorRect = {5,80,25,120,1};
    cnvs->StartDraw();
    cnvs->EraseAll ();
    cnvs->SetColor (clr);
    cnvs->FillRect (colorRect);
    cnvs->SetColor (black);
    cnvs->DrawRect (colorRect);

    cnvs->DisplayText ("Click to change", 5+cnvs->GetFont().size, 120+cnvs->GetFont().size, true);
    cnvs->EndDraw  ();
    cnvs->_MarkForUpdate ();
}

//__________________________________________________________

bool    _HYChartOptionDialog::ProcessEvent (_HYEvent* e)
{
    bool        done = false;
    _String     firstArg;
    long        i,f,k;

    if (e->EventClass()==_hyButtonPushEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        i = MatchComponentID (firstArg);

        if (i==10) { // OK
            *res = true;
            postWindowCloseEvent (GetID());

            for (k=6; k<=9; k++) {
                _String tbv (((_HYTextBox*)GetObject (k))->GetText());
                args->Replace (k-6,&tbv, true);
            }

            args->Replace (4,new _String (((_HYPullDown*)GetObject(15))->GetSelection()-2),false);
            args->Replace (5,new _String (((_HYPullDown*)GetObject(16))->GetSelection()-2),false);
            _String       *tb;
            ((_HYTextBox*)GetObject (21))->StoreText (tb);
            args->Replace (6,tb,false);
            ((_HYTextBox*)GetObject (23))->StoreText (tb);
            args->Replace (7,tb,false);
            ((_HYTextBox*)GetObject (24))->StoreText (tb);
            args->Replace (8,tb,false);

            *msel = ((_HYPullDown*)GetObject (5))->GetSelection();
            postWindowCloseEvent (GetID());
        } else if (i==11) { // Cancel
            *res = false;
            postWindowCloseEvent (GetID());
        } else if (i==12) {
            parWin->GetColors().Duplicate (colorsPicked);
            parWin->SetLabels (((_HYTextBox*)GetObject (6))->GetText(),
                               ((_HYTextBox*)GetObject (7))->GetText(),
                               ((_HYTextBox*)GetObject (8))->GetText(),
                               ((_HYPullDown*)GetObject(5))->GetSelection(),
                               ((_HYTextBox*)GetObject (9))->GetText(),
                               ((_HYPullDown*)GetObject(15))->GetSelection()-2,
                               ((_HYPullDown*)GetObject(16))->GetSelection()-2,
                               ((_HYTextBox*)GetObject (21))->GetText().toNum(),
                               ((_HYTextBox*)GetObject (23))->GetText().toNum(),
                               ((_HYTextBox*)GetObject (24))->GetText().toNum());
        }

        done = true;
    } else if (e->EventClass()==_hyMenuSelChangeEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        i = MatchComponentID (firstArg);
        if (i==18) {
            DrawColorCanvas ();
            done = true;
        }
    } else if (e->EventClass()==_hyContextPopUp) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        i = MatchComponentID (firstArg);
        if (i==19) {
            i = ((_HYPullDown*)GetObject (18))->GetSelection();
            _HYColor curColor = LongToHYColor (colorsPicked->lData[i]);
            colorsPicked->lData[i] = HYColorToLong(SelectAColor (curColor, chartColorChange));
            DrawColorCanvas ();
            done = true;
        }
    }

    if (done) {
        DeleteObject (e);
        return true;
    }
    return _HYTWindow::ProcessEvent(e);
}

//__________________________________________________________

void ReadChartProcessors (bool distribs)
{
    _String     pathToModelTemplates;
    _List       receptacle;

    pathToModelTemplates =libDirectory&"ChartAddIns";

    if (distribs) {
        pathToModelTemplates = pathToModelTemplates&GetPlatformDirectoryChar()&"DistributionAddIns";
    }

    ScanDirectoryForFileNames (pathToModelTemplates,receptacle,false);

    for (long k=0; k<receptacle.lLength; k++) {
        FILE * thisFile = doFileOpen (((_String*)receptacle(k))->sData,"rb");
        if (thisFile) {
            _String buffer (thisFile);
            fclose (thisFile);
            if (buffer.sLength) {
                if (distribs) {
                    distribProcessors << receptacle (k);
                } else {
                    chartProcessors << receptacle(k);
                }
            }
        }
    }
}

//__________________________________________________________
// _HYDistributionChartWindow
//__________________________________________________________

_String         catMarginalMatrix ("_MARGINAL_MATRIX_"),
                catVariableCDF    ("_CATEGORY_VARIABLE_CDF_"),
                catVariableEvent  ("_CATEGORY_VARIABLE_EVENT_"),
                catVariableID     ("_CATEGORY_VARIABLE_ID_");

//__________________________________________________________

_HYDistributionChartWindow::_HYDistributionChartWindow (_String name, _List& columns, _Matrix& data, _List& vars, _HYWindow* pw):_HYChartWindow (name, columns, data, pw)
{
    _HYRect         canvasSettings = {HY_CHART_TABLE_HEIGHT+HY_SCROLLER_WIDTH,260,HY_CHART_TABLE_HEIGHT+HY_SCROLLER_WIDTH,260,
                                      HY_COMPONENT_V_SCROLL|HY_COMPONENT_BORDER_R
                                     };

    _List*    cvi     = BuildAtoms (data, vars);
    _HYHList* catVars = new _HYHList (canvasSettings, GetOSWindowData(), *cvi, false);
    checkPointer (catVars);
    DeleteObject (cvi);

    canvasSettings.top = canvasSettings.bottom = 18;
    canvasSettings.width = HY_COMPONENT_BORDER_B|HY_COMPONENT_BORDER_T|HY_COMPONENT_BORDER_R;

    _HYLabel* l5     =  new _HYLabel (canvasSettings, GetOSWindowData());
    checkPointer (l5);

    l5->SetBackColor (tableDefaultBk2);
    l5->SetAlignFlags(HY_ALIGN_LEFT);
    l5->SetText      ("Category Variables");

    _HYFont tf = ((_HYTable*)GetObject (2))->textFont;
    tf.style = HY_FONT_BOLD;
    l5->SetFont (tf);

    AddObject (catVars);           // 11
    AddObject (l5);                // 12


    SetCell   (HY_CHART_HEADER_ROW,0,l5);
    SetCell   (HY_CHART_HEADER_ROW,1,l5);
    SetCell   (HY_CHART_HEADER_ROW,2,l5);
    SetCell   (HY_CHART_HEADER_ROW,3,l5);
    SetCell   (HY_CHART_HEADER_ROW,4,(_HYGuiObject*)components(10));

    SetCell   (HY_CHART_TABLE_ROW,0,catVars);
    SetCell   (HY_CHART_TABLE_ROW,1,catVars);
    SetCell   (HY_CHART_TABLE_ROW,2,catVars);
    SetCell   (HY_CHART_TABLE_ROW,3,catVars);
    SetCell   (HY_CHART_TABLE_ROW,4,(_HYGuiObject*)components(3));

    catVars->SetMessageRecipient (this);

    if (distribProcessors.lLength==0) {
        ReadChartProcessors (true);
    }

    SetWindowRectangle (0,0,600,600);
    DeleteObject (catVars);
    DeleteObject (l5);
}

//__________________________________________________________

bool    _HYDistributionChartWindow::ProcessEvent (_HYEvent* e)
{
    return _HYChartWindow::ProcessEvent (e);
}

//__________________________________________________________

_List*  _HYDistributionChartWindow::BuildAtoms (_Matrix& data, _List& vars)
{
    atoms.Clear();
    atomNames.Clear();
    derived.Clear();
    derivedDependencies.Clear();
    derivedNames.Clear();
    derivedRemaps.Clear();
    atomSizes.Clear();
    atomMultiples.Clear();


    _List* catVarHL = new   _List;
    checkPointer (catVarHL);

    marginals.Duplicate (&data);

    for (long k = 0; k < vars.lLength; k++) {
        _List               *thisList = (_List*)vars(k);

        _List               anEntry,
                            dist;

        anEntry << (*thisList)(0);

        _Matrix *wts    =  (_Matrix*)(*thisList)(1),//cv->GetWeights(),
                 *vals    =  (_Matrix*)(*thisList)(2);//cv->GetValues();


        long          k2 = wts->GetHDim() * wts->GetVDim();
        atomSizes     << k2;

        _List       atomInfo;

        atomInfo    && wts;
        atomInfo    && vals;

        atoms && & atomInfo;
        atomNames << (*thisList)(0);
        atomMultiples << 1;

        _Constant   f1 (8.), f2 (6.);

        for (long k3 = 0; k3 < k2; k3++) {
            _Constant f3 (vals->theData[k3]);
            _FString* fs = (_FString*)f3.FormatNumberString(&f1,&f2);

            _String dEntry = *fs->theString & " (pr =";
            f3.SetValue (wts->theData[k3]);
            DeleteObject (fs);
            fs = (_FString*)f3.FormatNumberString(&f1,&f2);
            dEntry = dEntry & *fs->theString & ')';
            DeleteObject (fs);
            dist && & dEntry;
        }

        anEntry && & dist;
        (*catVarHL) && & anEntry;
    }

    for (long k = vars.lLength-2; k >= 0; k--) {
        atomMultiples.lData[k] = atomMultiples.lData[k+1]*atomSizes.lData[k+1];
    }

    return catVarHL;
}
//__________________________________________________________

void    _HYDistributionChartWindow::SetAtoms (_Matrix& data, _List& vars)
{
    _List       *cinf = BuildAtoms (data, vars);

    _HYHList    *clst = (_HYHList*)GetObject (11);
    clst->ClearTable();

    for (long k=0; k<cinf->lLength; k++) {
        _List * anItem = (_List*)(*cinf)(k);
        clst->AddRubrik (*(_String*)(*anItem)(0),*(_List*)(*anItem)(1),-1);
    }

    DeleteObject (cinf);
}

//__________________________________________________________

void    _HYDistributionChartWindow::AddVariable (_String * expr)
{
    _String newExpression,
            aPrompt ("Expression for the new variable"),
            cPrompt ("Remember this expression");

    long    msel;

    bool    resC = savedExpression.sLength;

    if (resC) {
        newExpression = savedExpression;
    }

    if (expr || EnterStringDialogWithPulldown (newExpression, aPrompt, cPrompt, msel, atomNames, &atomNames, resC, 0, (Ptr)this)) {
        if (expr) {
            newExpression = *expr;
        }

        if (resC) {
            savedExpression = newExpression;
        } else {
            savedExpression = empty;
        }

        aPrompt = newExpression;

        _Formula f (newExpression, nil,false);

        if (f.IsEmpty()) {
            newExpression = _String("Failed to parse the expression :") & aPrompt;
            ProblemReport (newExpression, (Ptr)this);
            return;
        }

        // now check to see that all variables have been properly defined

        _SimpleList         scanVariables,
                            map;

        {
            _AVLList sva (&scanVariables);
            f.ScanFForVariables (sva, true, true, true);
            sva.ReorderList();
        }

        long k = 0;

        for (; k < scanVariables.countitems(); k++) {
            _Variable* thisV = LocateVar (scanVariables.lData[k]);
            long f = atomNames.Find (thisV->GetName());
            if (f == -1) {
                if (!thisV->IsCategory()) {
                    scanVariables.Delete (k);
                    k--;
                    continue;
                }
                break;
            }
            map << f;
        }

        if (k && k == scanVariables.lLength) {
            SortLists (&map, &scanVariables);
            _HYHList * hl = (_HYHList*)GetObject (11);

            _List      newDistribution,
                       rubrikItems;

            ProduceDistribution (scanVariables, map, newDistribution, f, 0, 1.);

            derivedDependencies && & map;
            derivedNames && & aPrompt;

            k = newDistribution.lLength/2;

            _SimpleList     derivedRemap (k, 0, 1);

            _Matrix         wts     (k,1,false,true),
                            probs   (k,1,false,true);

            _Constant   f1 (8.),
                        f2 (6.);

            bool        goOn = true;

            while (goOn) {
                goOn = false;

                for (long k3 = 1; k3 < k; k3++) {
                    _Constant   *f3 = (_Constant*)newDistribution(2*derivedRemap[k3-1]),
                                 *f4 = (_Constant*)newDistribution(2*derivedRemap[k3]);

                    if (f3->Value () > f4->Value()) {
                        goOn = true;
                        long   t = derivedRemap[k3];
                        derivedRemap[k3] = derivedRemap[k3-1];
                        derivedRemap[k3-1] = t;
                    }
                }
            }


            for (long k2 = 0; k2 < k; k2++) {
                _Constant   *f3 = (_Constant*)newDistribution(2*derivedRemap[k2]),
                             *f4 = (_Constant*)newDistribution(2*derivedRemap[k2]+1);

                _FString* fs = (_FString*)f3->FormatNumberString(&f1,&f2);
                _String dEntry = *fs->theString & " (pr =";
                DeleteObject (fs);
                fs = (_FString*)f4->FormatNumberString(&f1,&f2);
                dEntry = dEntry & *fs->theString & ')';
                DeleteObject (fs);
                rubrikItems && & dEntry;

                wts.theData[k2] = f4->Value();
                probs.theData[k2] = f3->Value();
            }

            newDistribution.Clear();
            newDistribution && & wts;
            newDistribution && & probs;

            derived && & newDistribution;
            derivedRemaps && & derivedRemap;

            hl->AddRubrik (aPrompt, rubrikItems, -1);
        } else {
            newExpression = aPrompt & " must depend on at least one model category variable and no other variables";
            ProblemReport (newExpression, (Ptr)this);
            return;
        }
    }
}
//__________________________________________________________

void        _HYDistributionChartWindow::ProduceDistribution (_SimpleList& vars, _SimpleList& map, _List& res, _Formula& f, long idx, _Parameter p)
{
    long vIdx = map.lData[idx];

    _Matrix * wts   = (_Matrix*) (*(_List*)atoms (vIdx))(0),
              * probs = (_Matrix*) (*(_List*)atoms (vIdx))(1);

    _Variable*v     = LocateVar (vars.lData[idx]);

    for (long k = 0; k < wts->GetHDim () * wts->GetVDim (); k++) {
        v->SetValue (new _Constant ((*probs)[k]),false);

        if (idx < vars.lLength-1) {
            ProduceDistribution (vars, map, res, f, idx+1, p*(*wts)[k]);
        } else {
            res.AppendNewInstance (new _Constant (f.Compute ()->Value()));
            res.AppendNewInstance (new _Constant (p*(*wts)[k]));
        }
    }
}

//__________________________________________________________

void        _HYDistributionChartWindow::ShowMarginals   (void)
{
    _List           choices;

    _SimpleList     all,
                    selectors;

    selectors << 0;
    selectors << 1;

    for (long vi = 0; vi < atomNames.lLength; vi++) {
        _List   aChoice;
        aChoice << atomNames (vi);
        aChoice << atomNames (vi);
        choices && & aChoice;
        all << all.lLength;
    }

    for (long vi = 0; vi < derivedNames.lLength; vi++) {
        _List   aChoice;
        aChoice << derivedNames (vi);
        aChoice << derivedNames (vi);
        choices && & aChoice;
        all << all.lLength;
    }

    long idx = HandleListSelection (choices,selectors, all, "Variable to Delete", all,1,(Ptr)this);
    if (idx >= 0) {

        _Matrix * res = ComputeConditionals (idx);
        _String aName = _String ("Marginal distributions in ") & *(_String*) (*(_List*)choices(idx))(0),
                testName = aName;

        _List     labels;
        for (idx = 1; idx <= marginals.GetHDim(); idx++) {
            _String aLabel = _String ("Site ") & idx;
            labels && & aLabel;
        }

        idx = 2;
        while (FindWindowByName (testName) >= 0) {
            testName = aName & " " & idx++;
        }

        _HYChartWindow * conditionals = new _HYChartWindow (testName,labels,*res);
        checkPointer (conditionals);
        conditionals->BringToFront ();
        DeleteObject (res);
    }
}

//__________________________________________________________

_Matrix*    _HYDistributionChartWindow::ComputeConditionals (long idx)
{
    _SimpleList varstc,
                *remap = nil;

    if (idx < atoms.lLength) {
        varstc << idx;
    } else {
        varstc << *(_SimpleList*)derivedDependencies (idx-atoms.lLength);
        varstc.Sort();
        remap = (_SimpleList*)derivedRemaps (idx-atoms.lLength);
    }

    return  ComputeConditionals (varstc, remap);

}

//__________________________________________________________

_Matrix*    _HYDistributionChartWindow::MakeCDF (long idx)
{
    _Matrix*    wts,
                *    prb;

    if (idx < atoms.lLength) {
        wts = (_Matrix*)(*(_List*)atoms(idx))(0);
        prb = (_Matrix*)(*(_List*)atoms(idx))(1);

    } else {
        idx -= atoms.lLength;
        wts = (_Matrix*)(*(_List*)derived(idx))(0);
        prb = (_Matrix*)(*(_List*)derived(idx))(1);
    }

    idx = wts->GetHDim() * wts->GetVDim();
    _Matrix * res = new _Matrix (2, idx , false, true);
    checkPointer (res);

    for (long k=0; k<idx; k++) {
        res->Store (0,k, wts->theData[k]);
        res->Store (1,k, prb->theData[k]);
    }

    return res;
}

//__________________________________________________________

void    _HYDistributionChartWindow::RemoveVariable (void)
{
    if (derived.lLength) {
        _List           choices;

        _SimpleList     all,
                        selectors;

        selectors << 0;
        selectors << 1;

        for (long vi = 0; vi < derivedNames.lLength; vi++) {
            _List   aChoice;
            aChoice << derivedNames (vi);
            aChoice << derivedNames (vi);
            choices && & aChoice;
            all << all.lLength;
        }

        long choice = HandleListSelection (choices,selectors, all, "Variable to Delete", all,1,(Ptr)this);
        if (choice >= 0) {
            _HYHList * hl = (_HYHList*)GetObject (11);
            hl->DeleteRubrik            (choice + atoms.lLength);
            derivedNames.Delete         (choice);
            derived.Delete              (choice);
            derivedDependencies.Delete  (choice);
            derivedRemaps.Delete        (choice);
        }
    } else {
        _String errMsg = "There are no user defined variables to delete";
        ProblemReport (errMsg, (Ptr)this);
    }
}

//__________________________________________________________

_Matrix*    _HYDistributionChartWindow::ComputeConditionals (_SimpleList& varstc, _SimpleList* remap)
{
    if (varstc.lLength == atoms.lLength) {
        if (remap) {
            _Matrix * cond = new _Matrix (marginals.GetVDim(), marginals.GetHDim(), false, true);
            checkPointer (cond);

            long    tc = marginals.GetHDim(),
                    tr = marginals.GetVDim();

            for (long c=0; c<tc; c++)
                for (long r=0; r<tr; r++) {
                    cond->Store (r,c,marginals(c,remap->lData[r]));
                }

            return      cond;
        } else {
            _Matrix* cnd = new _Matrix (marginals);
            checkPointer (cnd);
            cnd->Transpose();
            return cnd;
        }
    } else {
        long            ts = 1,
                        colCount = marginals.GetHDim();

        for (long k=0;  k < varstc.lLength; k++) {
            long       tsz = atomSizes.lData [varstc.lData[k]];
            ts *= tsz;
        }

        _Matrix * cond = new _Matrix (ts, colCount , false, true);
        checkPointer (cond);

        _SimpleList  toAdjust   (varstc.lLength,0,0),
                     expectOver    (atoms.lLength,0,1),
                     adjustable;

        adjustable.Subtract (expectOver, varstc);
        toAdjust.lData[toAdjust.lLength-1] = -1;

        for (long cs = 0; cs < ts; cs++) {
            long sp = 0;
            while (toAdjust.lData[sp] == atomSizes.lData[varstc.lData[sp]]-1) {
                toAdjust.lData[sp++] = 0;
            }

            toAdjust.lData[sp] ++;

            long   totalAdditive = 0;

            for (sp=0; sp<varstc.lLength; sp++) {
                totalAdditive += toAdjust.lData[sp] * atomMultiples.lData[varstc.lData[sp]];
            }

            _SimpleList adjOver (adjustable.lLength,0,0);

            long r = remap?remap->lData[cs]:cs;
            adjOver.lData[adjOver.lLength-1] = -1;

            while (1) {
                long sp2 = 0;
                while ((sp2<adjustable.lLength) && (adjOver.lData[sp2] == atomSizes.lData[adjustable.lData[sp2]]-1)) {
                    adjOver.lData[sp2++] = 0;
                }

                if (sp2 == adjustable.lLength) {
                    break;
                }

                adjOver.lData[sp2] ++;

                _Parameter  pr = 1.;

                long ta2  = totalAdditive;

                for  (sp2=0; sp2<adjustable.lLength; sp2++) {
                    ta2 += adjOver.lData[sp2] * atomMultiples.lData[adjustable.lData[sp2]];
                    pr *= ((_Matrix*)((*(_List*)atoms(adjustable.lData[sp2]))(0)))->theData[adjOver.lData[sp2]];
                }

                for (long cc = 0; cc < colCount; cc++) {
                    cond->Store (r,cc,(*cond)(r,cc) + pr*marginals(cc,ta2));
                }
            }
        }

        return cond;
    }
}

//__________________________________________________________

void    _HYDistributionChartWindow::HandleCatPostProcessor (long procIdx)
{
    if (procIdx<distribProcessors.lLength) {
        _HYHList*               hList = (_HYHList*)    GetObject (11);
        _SimpleList             tableSel;

        hList->GetRowSelection  (tableSel,1);

        if (tableSel.lLength > 0) {
            long ridx = hList->IsSingleRubrik (tableSel);
            if (ridx>=0) {
                FILE * thisFile = doFileOpen (((_String*)distribProcessors(procIdx))->sData,"rb");
                if (thisFile) {
                    _String buffer (thisFile);
                    fclose (thisFile);

                    tableSel.Offset (-1);
                    if (tableSel.lLength>1)
                        if (tableSel.lData[0] == -1) {
                            tableSel.Delete (0);
                        }

                    _Matrix * m = ComputeConditionals (ridx),
                              * c = MakeCDF           (ridx),
                                * e = new _Matrix         (tableSel);

                    _FString*   vn;

                    if (ridx < atomNames.lLength) {
                        vn = new _FString (*(_String*)atomNames(ridx));
                    } else {
                        vn = new _FString (*(_String*)derivedNames(ridx-atomNames.lLength));
                    }


                    checkPointer (c);
                    checkPointer (e);
                    checkPointer (m);
                    checkPointer (vn);


                    setParameter (catMarginalMatrix,m);
                    setParameter (catVariableCDF,   c);
                    setParameter (catVariableEvent, e);
                    setParameter (catVariableID, vn);

                    DeleteObject (c);
                    DeleteObject (e);
                    DeleteObject (m);
                    DeleteObject (vn);

                    long   g = batchLanguageFunctionNames.lLength;

                    _String dpp = *((_String*)distribProcessors(procIdx));
                    PushFilePath (dpp);

                    _ExecutionList   thisList;
                    terminateExecution = false;
                    thisList.BuildList (buffer);
                    thisList.ExecuteAndClean(g);

                    PopFilePath ();

                    DeleteVariable (catMarginalMatrix);
                    DeleteVariable (catVariableCDF);
                    DeleteVariable (catVariableEvent);
                } else {
                    _String eMsg = _String("Problem reading file:") & *((_String*)distribProcessors(procIdx));
                    ProblemReport (eMsg);
                }
            } else {
                _String wm2 ("All values (events) must come from a single category variable.");
                ProblemReport (wm2, (Ptr)this);
            }
        } else {
            _String wm    ("Please select some category variables to process");
            ProblemReport (wm, (Ptr)this);
        }
    }
}

//__________________________________________________________

void    _HYDistributionChartWindow::DoSave (long option, _String*)
{
    _String result (128L, true);

    for (long k=0; k<atoms.lLength; k++) {
        if (k) {
            result << ';';
        }
        result << (_String*)atomNames (k);
        for (long p = 0; p < 2; p++) {
            _Matrix * m = (_Matrix*)(*(_List*)atoms(k))(p);
            for (long n=0; n <atomSizes.lData[k]; n++) {
                result << ':';
                result << _String(m->theData[n]);
            }
        }
    }

    for (long k=0; k<derivedNames.lLength; k++) {
        result << ';';
        result << (_String*)derivedNames (k);
    }

    result.Finalize();
    _HYChartWindow::DoSave (option, &result);
}



// EOF
