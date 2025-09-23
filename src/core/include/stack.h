/*
 
 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (spond@ucsd.edu)
 Art FY Poon    (apoon42@uwo.ca)
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

#ifndef     __STACK__
#define     __STACK__

#include "mathobj.h"

//__________________________________________________________________________________
class   _Stack   //computational stack
{

    friend class _Formula;
    friend class _Operation;
    friend class _Variable;

public:

    /**
     * @brief Construct a new _Stack object
     */
    _Stack (void);
    /**
     * @brief Destroy the _Stack object
     */
    virtual ~_Stack (void);

    /**
     * @brief Push an item onto the stack
     *
     * @param data The item to push
     * @param store Whether to store the item
     * @return true if successful, false otherwise
     */
    bool      Push (HBLObjectRef, bool = true);     // push object onto the stack
    /**
     * @brief Pop an item from the stack
     *
     * @param del Whether to delete the item
     * @return HBLObjectRef The popped item
     */
    HBLObjectRef Pop (bool del = true);            // pop object from the top of the stack
    /**
     * @brief Get the stack depth
     *
     * @return long The stack depth
     */
    long      StackDepth (void) const;    // returns the depth of the stack
    /**
     * @brief Reset the stack
     */
    void      Reset (void);         // clear the stack
    /**
     * @brief Peek at an item on the stack
     *
     * @param offset The offset from the top of the stack
     * @return HBLObjectRef The item
     */
    HBLObjectRef    Peek (long offset = 0L);
        // peek at the object 'offset'
        // units from the top of the stack
    
    /**
     * @brief Initialize the stack
     */
    virtual   void    Initialize (void);
    /**
     * @brief Duplicate the stack
     *
     * @param brc The stack to duplicate
     */
    virtual   void    Duplicate (BaseRefConst);

protected:

    _List  theStack;
};

#endif
