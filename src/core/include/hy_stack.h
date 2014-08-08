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

#ifndef _HY_STACK_
#define _HY_STACK_

#include "hy_list.h"

template <typename PAYLOAD, typename STACK_CONTAINER = _hyList<PAYLOAD> >
class _hyStack {

private:
  STACK_CONTAINER c;

public:
  //does nothing
  _hyStack();

  // data constructor (1 member list)
  _hyStack(const PAYLOAD);

  //destructor
  virtual ~_hyStack(void);


  /*
  ==============================================================
  Methods
  ==============================================================
  */

  /**
   * Push a value on the stack
   * Argument is a pointer to make it possible to overload Clone
   * @return None.
   */
  void push(const PAYLOAD);
  void push(PAYLOAD *);

  /**
   * Clear the current list and make a copy from the argment
   * Argument is a pointer to make it possible to overload Clone
   * @param clone_from the object to clone from
   * @return None.
   */
  PAYLOAD pop();    

  // returns the depth of the stack
  unsigned long Length () const;

  // clear the stack
  void reset();                  

};

#include "hy_stack.cpp"

#endif
