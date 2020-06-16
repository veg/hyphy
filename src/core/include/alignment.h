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

#ifndef __ALIGNMENT_HEADER_FILE__

#define __ALIGNMENT_HEADER_FILE__
#include "baseobj.h"
#include "likefunc.h"
#include "matrix.h"
#include "simplelist.h"

double AlignStrings( char const * r_str
                   , char const * q_str
                   , char * & r_res
                   , char * & q_res
                   , long * char_map
                   , double * cost_matrix
                   , const long cost_stride
                   , const char gap
                   , double open_insertion
                   , double extend_insertion
                   , double open_deletion
                   , double extend_deletion
                   , double miscall_cost
                   , const bool do_local
                   , const bool do_affine
                   , const bool do_codon
                   , const long char_count
                   , double * codon3x5
                   , double * codon3x4
                   , double * codon3x2
                   , double * codon3x1
                   , const bool do_true_local = false
                   );

hyFloat LinearSpaceAlign( _String const * s1           // first string
                           , _String const * s2           // second string
                           , long*  cmap     // char -> position in scoring matrix mapper
                           , _Matrix * ccost        // NxN matrix of edit distances on characters
                           , hyFloat gopen       // the cost of opening a gap in sequence 1
                           , hyFloat gextend     // the cost of extending a gap in sequence 1 (ignored unless doAffine == true)
                           , hyFloat gopen2      // the cost of opening a gap in sequence 2
                           , hyFloat gextend2    // the cost of opening a gap in sequence 2   (ignored unless doAffine == true)
                           , bool doLocal           // ignore prefix and suffix gaps
                           , bool doAffine          // use affine gap penalties
                           , _SimpleList & ops      // edit operations for the optimal alignment
                           , hyFloat scoreCheck  // check the score of the alignment
                           , long from1
                           , long to1
                           , long from2
                           , long to2
                           , _Matrix ** buffer      // matrix storage,
                           , char parentGapLink
                           , char * ha
                           );

#endif
