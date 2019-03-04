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

#ifndef _TimeDifference_
#define _TimeDifference_

#if defined   __UNIX__
  #include <sys/time.h>
  #include <unistd.h>
#endif

#if defined(__APPLE__) && defined(__MACH__)
  #include <mach/mach.h>
  #include <mach/mach_time.h>
#endif



class TimeDifference {
  /** 
    A platform wrapper for Time Difference Functions
   */
public:
  
  TimeDifference (void);
  /** Create an instance of the time difference object;
      starts the timer by default
   */
  
  void Start (void);
  /** Start the timer, or reset an existing timer */
  
  double TimeSinceStart (void) const;
  /** Return time, in seconds, since the last 'Start' call */
  
  
private:
  #if defined(__APPLE__) && defined(__MACH__)
    uint64_t base_time;
  #elif defined   __UNIX__
    timeval base_time;
  #endif

};

#endif
