
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

#pragma once
#include <QDialog>
#include "ui_hyphyhierarchicalselector.h"
#include "list.h"
#include "hy_strings.h"


class _HY_HierarchicalSelector : public QDialog, private Ui::Dialog
{
    Q_OBJECT

public:
    _HY_HierarchicalSelector(QWidget *parent, _List& definition, _SimpleList& c, _SimpleList& vc, _String n, _SimpleList* s, long r, bool is_modal = true);
virtual
    ~_HY_HierarchicalSelector () {};
    void SetInitialSelection ();
    
protected:

    bool eventFilter ( QObject * watched, QEvent * event );

private:   
    void  toggleAcceptStatus(long);
    void  closeEvent(QCloseEvent *event);

    _List data,
          iData,
          dData;
          
    _SimpleList validChoices,
                *selections,
                offsets;
                
    long          result;
    _String       lastString;
    bool          validSelection;

private slots:
    void ok();
    void cancel();
    void handle_selection_change ();
    void handle_double_click ();
};
