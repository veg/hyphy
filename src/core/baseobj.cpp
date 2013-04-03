/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2007

Primary Development:
  Sergei L Kosakovsky Pond (sergeilkp@mac.com)
Significant contributions from:
  Spencer V Muse (muse@stat.ncsu.edu)
  Simon DW Frost (sdfrost@ucsd.edu)
  Art FY Poon    (apoon@biomail.ucsd.edu)

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

#include "baseobj.h"
#ifdef   __HYPHYMPI__
#include "likefunc.h"
extern int _hy_mpi_node_rank;
#endif

#if defined   __UNIX__ || defined __HYPHY_GTK__
#include <sys/time.h>
#include <unistd.h>
#endif

#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

#ifdef   __HYPHYXCODE__
#include "HYUtils.h"
#endif

#ifdef __WINDOZE__
    #include <Windows.h>
#endif




#include    "baseobj.h"
#include    "helperfunctions.h"
#include    "hy_globals.h"


//____________________________________________________________________________________

BaseObj::BaseObj()
{
    nInstances=1;
}


//____________________________________________________________________________________
BaseRef   BaseObj::toStr (void)
{
    return new _String ("<HyPhy Base Object>");
}

//____________________________________________________________________________________
BaseRef   BaseObj::toErrStr (void)
{
    return toStr();
}

//____________________________________________________________________________________
void     BaseObj::toFileStr (FILE* dest)
{
    _String* s = (_String*)toStr();
    fwrite(s->sData,1,s->Length(),dest);
    DeleteObject (s);
}

//____________________________________________________________________________________
BaseObj*  BaseObj::makeDynamic(void)
{
    warnError(-112);
    return nil;
}





//EOF
