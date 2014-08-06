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

#include "hy_globals.h"

bool terminateExecution = false, isInFunction = false;

char isDefiningATree = 0;


long globalRandSeed;

FILE *globalErrorFile = nil, *globalMessageFile = nil;



_String compileDate = __DATE__,
        __KERNEL__VERSION__ =
                        _String("3.00") & compileDate.Cut(7, 10) &
                        compileDate.Cut(0, 2).Replace("Jan", "01", true)
                        .Replace("Feb", "02", true).Replace("Mar", "03", true)
                        .Replace("Apr", "04", true).Replace("May", "05", true)
                        .Replace("Jun", "06", true).Replace("Jul", "07", true)
                        .Replace("Aug", "08", true).Replace("Sep", "09", true)
                        .Replace("Oct", "10", true).Replace("Nov", "11", true)
                        .Replace("Dec", "12", true) &
                        compileDate.Cut(4, 5).Replace(" ", "0", true) & "alpha",
        empty(""),
        emptyAssociativeList("{}"),
        hyphyCiteString(
                "\nPlease cite S.L. Kosakovsky Pond, S. D. W. Frost and S.V. Muse. "
                "(2005) HyPhy: hypothesis testing using phylogenies. Bioinformatics "
                "21: 676-679 if you use HyPhy in a publication\nIf you are a new HyPhy "
                "user, the tutorial located at http://www.hyphy.org/docs/HyphyDocs.pdf "
                "may be a good starting point.\n"),
        scanfLastFilePath,
        errorFileName("errors.log"),
        messageFileName("messages.log");



#ifndef HY_2014_REWRITE_MASK

_SimpleList freeSlots;
_List openFileHandlesBackend;

_AVLListX openFileHandles (&openFileHandlesBackend);

#endif