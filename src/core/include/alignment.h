
#include "baseobj.h"
#include "likefunc.h"
#include "matrix.h"
#include "simplelist.h"

_Parameter AlignStrings( char * r_str
                       , char * q_str
                       , long * char_map
                       , double * cost_matrix
                       , const long cost_stride
                       , const char gap
                       , _Parameter open_insertion
                       , _Parameter extend_insertion
                       , _Parameter open_deletion
                       , _Parameter extend_deletion
                       , _Parameter miscall_cost
                       , const bool do_local
                       , const bool do_affine
                       , const bool do_codon
                       , _List & store
                       , const long char_count
                       , double * codon3x5
                       , double * codon3x4
                       , double * codon3x2
                       , double * codon3x1
                       );

_Parameter LinearSpaceAlign( _String * s1           // first string
                           , _String * s2           // second string
                           , _SimpleList & cmap     // char -> position in scoring matrix mapper
                           , _Matrix * ccost        // NxN matrix of edit distances on characters
                           , _Parameter gopen       // the cost of opening a gap in sequence 1
                           , _Parameter gextend     // the cost of extending a gap in sequence 1 (ignored unless doAffine == true)
                           , _Parameter gopen2      // the cost of opening a gap in sequence 2
                           , _Parameter gextend2    // the cost of opening a gap in sequence 2   (ignored unless doAffine == true)
                           , bool doLocal           // ignore prefix and suffix gaps
                           , bool doAffine          // use affine gap penalties
                           , _SimpleList & ops      // edit operations for the optimal alignment
                           , _Parameter scoreCheck  // check the score of the alignment
                           , long from1
                           , long to1
                           , long from2
                           , long to2
                           , _Matrix ** buffer      // matrix storage,
                           , char parentGapLink
                           , char * ha
                           );
