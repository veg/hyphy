/*
 
 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
 Art FY Poon    (apoon42@uwo.ca)
 Steven Weaver (sweaver@temple.edu)
 
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

#include "stack.h"

//__________________________________________________________________________________
_Stack::_Stack (void)
{
}

//__________________________________________________________________________________
void _Stack::Initialize (void)
{
    theStack.Initialize();
}

//__________________________________________________________________________________
void _Stack::Duplicate (BaseRefConst s) {
    theStack.Duplicate(&((_Stack*)s)->theStack);
}

//__________________________________________________________________________________
_Stack::~_Stack (void)
{
}

//__________________________________________________________________________________
bool _Stack::Push (HBLObjectRef newObj, bool dup) {   // push object onto the stack
    if (dup)
        theStack<<(newObj);
    else
        theStack.AppendNewInstance(newObj);
    return true;
}

//__________________________________________________________________________________
HBLObjectRef _Stack::Pop (bool del)   {      // pop object from the top of the stack
    HBLObjectRef r = (HBLObjectRef)theStack.list_data[theStack.lLength-1];
    if (del) {
        theStack.lLength--;
    }
    return r;
}

//__________________________________________________________________________________
HBLObjectRef _Stack::Peek (long offset)   {      // pop object from the top of the stack
    if (offset < theStack.lLength) {
        return (HBLObjectRef)theStack.list_data[theStack.lLength-1-offset];
    }
    return nil;
}

//__________________________________________________________________________________
long _Stack::StackDepth (void) const { // returns the depth of the stack
    return theStack.lLength;
}

//__________________________________________________________________________________
void _Stack::Reset (void) {  // clears the stack
    /*for (long i = 0; i < theStack.lLength; i++) {
        HBLObjectRef this_item = (HBLObjectRef)theStack.GetItem(i);
        if (this_item->ObjectClass() == MATRIX) {
            printf ("%ld=>%ld (%p)\n", i, this_item->CanFreeMe(), this_item);
        }
    }*/
    theStack.Clear();
}
