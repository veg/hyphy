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

#include <QApplication>
#include <QWidget>
#include <QEvent>
#include <QMouseEvent>
#include <QDebug>
#include "hy_strings.h"
#include "hyphymain.h"
#include "hyphyevents.h"
#include "hyphyhierarchicalselector.h"

bool needExtraNL = true; 
bool dropIntoDebugMode=false;

void StringToConsole(_String&);
void BufferToConsole(const char* buffer);
_String* StringFromConsole (bool echo);
void SetStatusLine(_String arg);
void SetStatusLineUser(_String s);
void SetStatusBarValue (long l, _Parameter max, _Parameter rate);
bool Get_a_URL (_String& urls, _String* fileName);
void NLToConsole();

void StringToConsole(_String& str)
{
   BufferToConsole(str); 
}

void BufferToConsole(const char* buffer)
{
    
    //Get MainWindow object specifically
    QBufferToConsoleEvent event((QString)buffer);

    foreach (QWidget *widget, QApplication::allWidgets()) {
        //Write to main console
        QApplication::sendEvent(widget, &event);
    }

    return;
}

_String* StringFromConsole (bool echo)
{
    //Do nothing for right now
    return nil;
}

void SetStatusLine(_String arg)
{

}

void SetStatusLineUser(_String s)
{

}

void SetStatusBarValue (long l, _Parameter max, _Parameter rate)
{
}

bool Get_a_URL (_String& urls, _String* fileName)
{
    return false;
}

void NLToConsole()
{

}

//____________________________________________________________________________________________
long  HandleListSelection (_List& data, _SimpleList& choices, _SimpleList& validChoices, _String titleInfo, _SimpleList& selections, long fixedLength, Ptr prt)
{
    long res = -1;
    if (data.lLength < 1 || validChoices.lLength < 1) {
        _String errMsg ("An empty list of choices was passed to 'HandleListSelection'");
        WarnError (errMsg);
    } else {
        selections.Clear();
        _HY_HierarchicalSelector *hs = new _HY_HierarchicalSelector((QWidget*)prt, data, choices, validChoices, titleInfo, &selections, fixedLength, true);
        hs->setWindowModality(Qt::WindowModal);
        hs->exec();
        res = selections.lLength;
        
    }
    return res;
}

//____________________________________________________________________________________________
long  HandleListSelection (_List& data, _String titleInfo, Ptr prt)
{
    _SimpleList validChoices,
    choices,
    sels;
    
    _List       menuData;
    
    validChoices << 0;
    validChoices << 1;
    
    for (long k=0; k<data.lLength; k+=2) {
        _List aChoice;
        
        aChoice << data (k);
        aChoice << data (k+1);
        
        menuData && & aChoice;
        
        choices << choices.lLength;
    }
    
    return  HandleListSelection (menuData, validChoices, choices, titleInfo, sels, 1, prt);
}

