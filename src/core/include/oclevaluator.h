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

//From calcnode.h
#ifndef __OCLEVALUATOR__
#define __OCLEVALUATOR__

#ifdef MDSOCL

class _OCLEvaluator {

private:
  // OpenCL Vars

  // Forward Declarations
  // *********************************************************************
  void Cleanup(int iExitCode);
  unsigned int roundUpToNextPowerOfTwo(unsigned int x);
  double roundDoubleUpToNextPowerOfTwo(double x);
  // So the only thing that needs to be passed as an update for each LF is
  // flatTree and flatCLeaves
  // as those are what goes into the new transition matrix stuff.
  // So I could have a launchmdsocl that takes everything and if stuff is not
  // NULL than update, if it is
  // null use the existing values. How about that?
  // The problem is that I need to have essentially the first LF's information
  // to properly set everything up.
  // And how do I have subsequent LF's not pass stuff.
  // Oi, how do I not have to convert this into an object...
  // alright, I can probably keep all of this in likefunc. That is because the
  // LFEvaluation is done in the
  // calc node
  double oclmain(void);
  bool contextSet;
  int setupContext(void);

public:
  void init(long esiteCount, long ealphabetDimension, _Parameter *eiNodeCache);

  double launchmdsocl(_SimpleList &updateNodes, _SimpleList &flatParents,
                      _SimpleList &flatNodes, _SimpleList &flatCLeaves,
                      _SimpleList &flatLeaves, _SimpleList &flatTree,
                      _Parameter *theProbs, _SimpleList &theFrequencies,
                      long *lNodeFlags, _SimpleList &taggedInternals,
                      _GrowingVector *lNodeResolutions);
  ~_OCLEvaluator() { Cleanup(EXIT_SUCCESS); }

};

#endif
#endif
