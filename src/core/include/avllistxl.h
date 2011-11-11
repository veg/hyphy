/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2011
Primary Development:
  Sergei L Kosakovsky Pond (sergeilkp@mac.com)
Significant contributions from:
  Spencer V Muse (muse@stat.ncsu.edu)
  Simon DW Frost (sdfrost@ucsd.edu)
  Art FY Poon    (apoon@biomail.ucsd.edu)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/

#ifndef _AVLLISTXL_
#define _AVLLISTXL_

//#pragma once
#include "list.h"
#include "avllist.h"

#define  MEMORYSTEP 8

//_____________________________________________________________________________

class _AVLListXL: public _AVLList
{

public:

    _AVLListXL(_SimpleList*);
    BaseRef GetXtra(long);

    void SetXtra(long,BaseRef,bool);

    virtual ~_AVLListXL(void){}
    virtual BaseRef toStr(void);

    virtual long InsertData(BaseRef, long,bool);
    virtual void Clear(bool = false);
    virtual void DeleteXtra(long);

    _List xtraD;

};

//_____________________________________________________________________________

#endif
