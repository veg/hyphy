LoadFunctionLibrary("../IOFunctions.bf");
LoadFunctionLibrary("../all-terms.bf");


/**
 * Define scoring parameters for alignment tasks
 * @name mapping.DefineScoringMatrices
 * @param code - genetic code definition
 * @returns {Dictionary} a structure with various scoring options (opaque)
 */

namespace mapping {
 
    lfunction DefineScoringMatrices (code) {
    
            score_matrix = {
                {                 6, -3, -4, -4, -2, -2, -2, -1, -3, -3, -3, -2, -2, -4, -2,  0, -1, -5, -3, -1, -4, -2, -2, -7}
                {                -3,  8, -2, -4, -6,  0, -2, -5, -2, -6, -4,  1, -3, -5, -4, -2, -3, -5, -3, -5, -3, -1, -2, -7}
                {                -4, -2,  8,  0, -5, -1, -2, -2,  0, -6, -6, -1, -4, -5, -4,  0, -1, -7, -4, -5,  6, -2, -2, -7}
                {                -4, -4,  0,  8, -6, -2,  0, -3, -3, -5, -6, -2, -6, -6, -3, -1, -3, -7, -6, -6,  6,  0, -3, -7}
                {                -2, -6, -5, -6, 10, -5, -7, -5, -5, -3, -3, -6, -3, -5, -5, -2, -2, -4, -4, -2, -5, -6, -4, -7}
                {                -2,  0, -1, -2, -5,  8,  1, -4,  0, -6, -4,  0, -1, -6, -3, -1, -2, -3, -3, -4, -1,  6, -2, -7}
                {                -2, -2, -2,  0, -7,  1,  7, -4, -1, -6, -5,  0, -4, -6, -3, -1, -2, -5, -4, -4,  0,  6, -2, -7}
                {                -1, -5, -2, -3, -5, -4, -4,  7, -4, -7, -6, -3, -5, -5, -4, -2, -4, -4, -5, -6, -2, -4, -4, -7}
                {                -3, -2,  0, -3, -5,  0, -1, -4, 10, -6, -5, -2, -3, -3, -4, -2, -4, -5,  0, -6, -1, -1, -3, -7}
                {                -3, -6, -6, -5, -3, -6, -6, -7, -6,  6,  0, -5,  0, -1, -5, -5, -2, -5, -3,  2, -5, -6, -2, -7}
                {                -3, -4, -6, -6, -3, -4, -5, -6, -5,  0,  6, -5,  1, -1, -5, -5, -3, -3, -3,  0, -6, -5, -2, -7}
                {                -2,  1, -1, -2, -6,  0,  0, -3, -2, -5, -5,  7, -3, -6, -2, -1, -2, -5, -3, -4, -2,  0, -2, -7}
                {                -2, -3, -4, -6, -3, -1, -4, -5, -3,  0,  1, -3,  9, -1, -5, -3, -2, -3, -3,  0, -5, -2, -1, -7}
                {                -4, -5, -5, -6, -5, -6, -6, -5, -3, -1, -1, -6, -1,  8, -6, -4, -4,  0,  1, -3, -6, -6, -3, -7}
                {                -2, -4, -4, -3, -5, -3, -3, -4, -4, -5, -5, -2, -5, -6,  9, -2, -3, -6, -5, -4, -4, -3, -4, -7}
                {                 0, -2,  0, -1, -2, -1, -1, -2, -2, -5, -5, -1, -3, -4, -2,  7,  0, -5, -3, -4, -1, -1, -2, -7}
                {                -1, -3, -1, -3, -2, -2, -2, -4, -4, -2, -3, -2, -2, -4, -3,  0,  7, -4, -3, -1, -2, -2, -2, -7}
                {                -5, -5, -7, -7, -4, -3, -5, -4, -5, -5, -3, -5, -3,  0, -6, -5, -4, 12,  0, -6, -7, -4, -4, -7}
                {                -3, -3, -4, -6, -4, -3, -4, -5,  0, -3, -3, -3, -3,  1, -5, -3, -3,  0,  9, -3, -5, -3, -2, -7}
                {                -1, -5, -5, -6, -2, -4, -4, -6, -6,  2,  0, -4,  0, -3, -4, -4, -1, -6, -3,  6, -6, -4, -2, -7}
                {                -4, -3,  6,  6, -5, -1,  0, -2, -1, -5, -6, -2, -5, -6, -4, -1, -2, -7, -5, -6,  7, -1, -3, -7}
                {                -2, -1, -2,  0, -6,  6,  6, -4, -1, -6, -5,  0, -2, -6, -3, -1, -2, -4, -3, -4, -1,  7, -2, -7}
                {                -2, -2, -2, -3, -4, -2, -2, -4, -3, -2, -2, -2, -1, -3, -4, -2, -2, -4, -2, -2, -3, -2, -2, -7}
                {                -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7,  1}
            };

            max_score = Max (score_matrix,0);
            penalty = Max (max_score, -Min (score_matrix,0));


            options = {
                "SEQ_ALIGN_CHARACTER_MAP" : "ACGTN",
                "SEQ_ALIGN_GAP_OPEN"	 : 	1.5*penalty,
                "SEQ_ALIGN_GAP_OPEN2"	 :	1.5*penalty,
                "SEQ_ALIGN_GAP_EXTEND"   : 	0.1*penalty,
                "SEQ_ALIGN_GAP_EXTEND2"  : 	0.1*penalty,
                "SEQ_ALIGN_FRAMESHIFT"   :  2*penalty,
                "SEQ_ALIGN_NO_TP"        :  1,
                "SEQ_ALIGN_AFFINE"       :  1,
                "SEQ_ALIGN_CODON_ALIGN"  :  1
            };

           base_frequencies = {
                                    {   0.069352301}
                                    {   0.021000333}
                                    {   0.049862283}
                                    {   0.026563029}
                                    {   0.033026253}
                                    {     0.1054545}
                                    {   0.008970554}
                                    {   0.036648952}
                                    {   0.036329907}
                                    {   0.065164842}
                                    {   0.021805663}
                                    {   0.032992805}
                                    {     0.0348093}
                                    {   0.036818766}
                                    {   0.054168098}
                                    {    0.12149678}
                                    {   0.082464001}
                                    {   0.053564744}
                                    {   0.038383113}
                                    {    0.07112377}
                                };



            options["SEQ_ALIGN_SCORE_MATRIX"] = pSM2cSM(score_matrix, "ACDEFGHIKLMNPQRSTVWY", code["code"], code["ordering"]);
            shift_penalty = computeExpectedPerBaseScore (.5,score_matrix,base_frequencies);

            _cdnaln_partialScoreMatrices = cSM2partialSMs(options["SEQ_ALIGN_SCORE_MATRIX"],
                    {{shift_penalty__*1.5,shift_penalty__,shift_penalty__,shift_penalty*1.5}});


            options ["SEQ_ALIGN_PARTIAL_3x1_SCORES"] = _cdnaln_partialScoreMatrices["3x1"];
            options ["SEQ_ALIGN_PARTIAL_3x2_SCORES"] = _cdnaln_partialScoreMatrices["3x2"];
            options ["SEQ_ALIGN_PARTIAL_3x4_SCORES"] = _cdnaln_partialScoreMatrices["3x4"];
            options ["SEQ_ALIGN_PARTIAL_3x5_SCORES"] = _cdnaln_partialScoreMatrices["3x5"];

            options ["MATCH"] = computeExpectedPerBaseScore (1,score_matrix,base_frequencies);
            options ["E"] = Max (0.1, computeExpectedPerBaseScore (.4,score_matrix,base_frequencies));

            return options;
    }
    
    
    lfunction cSM2partialSMs(_cdnScoreMatrix, penalties) {

        m3x5  =  { 126, 1250 };
        m3x4  =  { 126, 500 };
        m3x2  =  { 126,  75 };
        m3x1  =  { 126,  15 };


        if (utility.Array1D (penalties) == 4) {
            p3x1 = penalties [0];
            p3x2 = penalties [1];
            p3x4 = penalties [2];
            p3x5 = penalties [3];
        } else {
            p3x5 = 0;
            p3x4 = 0;
            p3x2 = 0;
            p3x1 = 0;
        }

        for ( thisCodon = 0; thisCodon < 64; thisCodon += 1 ) {
            for ( d1 = 0; d1 < 5; d1 += 1 ) {
                max100 = -1e100;
                max010 = -1e100;
                max001 = -1e100;

                for ( d2 = 0; d2 < 5; d2 += 1 ) {
                    partialCodon = 5 * d1 + d2;
                    max110 = -1e100;
                    max101 = -1e100;
                    max011 = -1e100;

                    for ( d3 = 0; d3 < 5; d3 += 1 ) {
                        thisCodon2 = 5 * partialCodon + d3;
                        thisScore = _cdnScoreMatrix[ thisCodon ][ thisCodon2 ];

                        // this is the trivial and stupid way of doing it, but it should work
                        m3x5[ thisCodon ][ 10 * thisCodon2 + 0 ] = thisScore - p3x5;
                        m3x5[ thisCodon ][ 10 * thisCodon2 + 1 ] = thisScore - p3x5;
                        m3x5[ thisCodon ][ 10 * thisCodon2 + 2 ] = thisScore - p3x5;
                        m3x5[ thisCodon ][ 10 * thisCodon2 + 3 ] = thisScore - p3x5;
                        m3x5[ thisCodon ][ 10 * thisCodon2 + 4 ] = thisScore - p3x5;
                        m3x5[ thisCodon ][ 10 * thisCodon2 + 5 ] = thisScore - p3x5;
                        m3x5[ thisCodon ][ 10 * thisCodon2 + 6 ] = thisScore - p3x5;
                        m3x5[ thisCodon ][ 10 * thisCodon2 + 7 ] = thisScore - p3x5;
                        m3x5[ thisCodon ][ 10 * thisCodon2 + 8 ] = thisScore - p3x5;
                        m3x5[ thisCodon ][ 10 * thisCodon2 + 9 ] = thisScore - p3x5;

                        m3x4[ thisCodon ][ 4 * thisCodon2 + 0 ] = thisScore - p3x4;
                        m3x4[ thisCodon ][ 4 * thisCodon2 + 1 ] = thisScore - p3x4;
                        m3x4[ thisCodon ][ 4 * thisCodon2 + 2 ] = thisScore - p3x4;
                        m3x4[ thisCodon ][ 4 * thisCodon2 + 3 ] = thisScore - p3x4;

                        // d1 is 1
                        max100 = Max( max100, _cdnScoreMatrix[ thisCodon ][ 25 * d1 + 5 * d2 + d3 ] );
                        max010 = Max( max010, _cdnScoreMatrix[ thisCodon ][ 25 * d2 + 5 * d1 + d3 ] );
                        max001 = Max( max001, _cdnScoreMatrix[ thisCodon ][ 25 * d2 + 5 * d3 + d1 ] );

                        // d1 and d2 are 1
                        max110 = Max( max110, _cdnScoreMatrix[ thisCodon ][ 25 * d1 + 5 * d2 + d3 ] );
                        max101 = Max( max101, _cdnScoreMatrix[ thisCodon ][ 25 * d1 + 5 * d3 + d2 ] );
                        max011 = Max( max011, _cdnScoreMatrix[ thisCodon ][ 25 * d3 + 5 * d1 + d2 ] );
                    }

                    m3x2[ thisCodon ][ 3 * partialCodon + 0 ] = max110 - p3x2;
                    m3x2[ thisCodon ][ 3 * partialCodon + 1 ] = max101 - p3x2;
                    m3x2[ thisCodon ][ 3 * partialCodon + 2 ] = max011 - p3x2;
                }

                m3x1[ thisCodon ][ 3 * d1 + 0 ] = max100 - p3x1;
                m3x1[ thisCodon ][ 3 * d1 + 1 ] = max010 - p3x1;
                m3x1[ thisCodon ][ 3 * d1 + 2 ] = max001 - p3x1;
            }
        }


        return { "3x1": m3x1, "3x2": m3x2, "3x4": m3x4, "3x5": m3x5 };
    }



    lfunction _digits (index) {
        return {{index__$25, index__ % 25 $ 5, index__ % 5}};
    };

    lfunction _hazN (digits) {
        return digits [0] == 4 || digits [1] == 4 || digits [2] == 4;
    };

    lfunction _map_to_nuc (digits) {
        return digits [0] * 16 + digits [1] * 4 + digits[2];
    };

    lfunction _generate_resolutions (digits) {
        resolutions = {};
        for (k = 0; k < 4; k += 1) {
            try = digits;
            if (digits[0] == 4) {
                try [0] = k;
            }
            for (k2 = 0; k2 < 4; k2 += 1) {
                if (digits[1] == 4) {
                    try [1] = k;
                }
                for (k3 = 0; k3 < 4; k3 += 1) {
                    if (digits[2] == 4) {
                        try [2] = k;
                    }
                    resolutions[_map_to_nuc (try)] = 1;
                }
            }
        }
        return Rows(resolutions);
    };

    // -------------------------------------------------------------------------- //

    lfunction pSM2cSM (_scorematrix, _letters, code, ordering) {

        _cdnScoreMatrix  = { 126 , 126 };
        _mapping      = utility.MapStrings ( ordering, _letters );


        for ( _k = 0; _k < 125; _k += 1 ) {

            letters1 = _digits (_k);

            if (_hazN (letters1)) {
                 letters1 = _generate_resolutions (letters1);
                  for ( _k2 = _k; _k2 < 125 ; _k2 += 1 ) {
                    letters2 = _digits (_k2);
                    if (_hazN (letters2) == 0) {
                        codon2 = _map_to_nuc(letters2);
                        _mappedK2 = _mapping[ code[ codon2 ] ];
                        _aScore = -1e4;
                        if (_mappedK2 >= 0) {

                            res_count = utility.Array1D (letters1);
                            for (r = 0; r < res_count; r += 1) {
                                resolution_codon = 0 + letters1[r];
                                resolution_aa = _mapping[ code[ resolution_codon  ] ];
                                if (resolution_aa >= 0) {
                                    try_score = _scorematrix[ resolution_aa ][ _mappedK2 ] - 1;
                                    if (resolution_aa == _mappedK2 && codon2 != resolution_codon) {
                                        try_score = try_score - 1;
                                    }
                                    _aScore = Max (_aScore, try_score);
                                }
                            }
                        }


                        _cdnScoreMatrix[ _k ][ _k2 ] = _aScore;
                        _cdnScoreMatrix[ _k2 ][ _k ] = _aScore;
                    }
                 }
            } else {
                codon1 = _map_to_nuc(letters1);

                _mappedK = _mapping[ code[ codon1 ] ];

                if ( _mappedK >= 0) {
                    for ( _k2 = _k; _k2 < 125 ; _k2 += 1 ) {
                        letters2 = _digits (_k2);

                        if (_hazN (letters2)) {
                           _aScore = -1e4;
                           letters2 = _generate_resolutions (letters2);
                            res_count = utility.Array1D (letters2);
                            for (r = 0; r < res_count; r += 1) {
                                resolution_codon = 0 + letters2[r];
                                resolution_aa = _mapping[ code[ resolution_codon  ] ];
                                if (resolution_aa >= 0) {
                                    try_score = _scorematrix[ _mappedK ][ resolution_aa ] - 1;
                                    if (resolution_aa == _mappedK && codon1 != resolution_codon) {
                                        try_score = try_score - 1;
                                    }
                                    _aScore = Max (_aScore, try_score);
                                }
                            }
                            _cdnScoreMatrix[ _k ][ _k2 ] = _aScore;
                            _cdnScoreMatrix[ _k2 ][ _k ] = _aScore;
                            continue;
                        }

                        codon2 = _map_to_nuc(letters2);

                        _mappedK2 = _mapping[ code[ codon2 ] ];
                        if ( _mappedK2 >= 0 ) {
                            _aScore = _scorematrix[ _mappedK ][ _mappedK2 ];
                            if ( _mappedK == _mappedK2 && _k2 > _k ) {
                                _aScore = _aScore - 1; // synonymous match
                            }
                        } else {
                            // stop codons don't match anything
                            _aScore = -1e4;
                        }
                        _cdnScoreMatrix[ _k ][ _k2 ] = _aScore;
                        _cdnScoreMatrix[ _k2 ][ _k ] = _aScore;
                    }
                } else { // stop codons here
                    for ( _k2 = _k; _k2 < 125; _k2 += 1 ) {

                        letters2 = _digits (_k2);

                        if (_hazN (letters2)) {
                            continue;
                        }

                        codon2 = _map_to_nuc(letters2);

                        _mappedK2 = _mapping[ code[ codon2 ] ];

                        if ( _mappedK2 < 0 ) {
                            // don't penalize stop codons matching themselves
                            _cdnScoreMatrix[ _k ][ _k2 ] = 0;
                            _cdnScoreMatrix[ _k2 ][ _k ] = 0;
                        } else {
                            _cdnScoreMatrix[ _k ][ _k2 ] = -1e4;
                            _cdnScoreMatrix[ _k2 ][ _k ] = -1e4;
                        }
                    }
                }
            }
        }

        return _cdnScoreMatrix;
    }

    // -------------------------------------------------------------------------- //

    lfunction computeExpectedPerBaseScore( _expectedIdentity, _cdnaln_scorematrix, _cdnaln_base_freqs ) {
        meanScore = 0;

        for (_aa1 = 0; _aa1 < 20; _aa1 += 1) {
            for (_aa2 = 0; _aa2 < 20; _aa2 += 1) {
                if ( _aa1 != _aa2 ) {
                    meanScore += ( 1 - _expectedIdentity ) * _cdnaln_scorematrix[_aa1][_aa2] * _cdnaln_base_freqs[_aa1] * _cdnaln_base_freqs[_aa2];
                } else {
                    meanScore += _expectedIdentity * _cdnaln_scorematrix[_aa1][_aa1] * _cdnaln_base_freqs[_aa1] * _cdnaln_base_freqs[_aa1];
                }
            }
        }

        return meanScore;
    }
    
    // -------------------------------------------------------------------------- //

    lfunction	computeCorrection (str) {
        result = {1,2};

        result[0]	 = (str$"^\\-+")[1]+1;
        result[1]	 = (str$"\\-+$")[0];

        if (result[1] >= 0) {
            result[1] = Abs(str)-result[1];
        }
        else {
            result[1] = 0;
        }
        return result;
    }

    
    // -------------------------------------------------------------------------- //

    lfunction igg_alignment_cleanup (reference, query, offset_nuc, code) {

        too_short = 0;
        too_long  = 0;
        span      = 0; // how many nucleotides in the reference were covered by non-gaps
        _seqL     = Abs (reference);

   
        ref_cleaned = ""; ref_cleaned * 128;
        qry_cleaned = ""; qry_cleaned * 128;

        _codon_in_reference = 0;

        for ( _rcidx = 0; _rcidx < _seqL; _rcidx += 1 ) {
            _del1 = reference [_rcidx] != (reference [_rcidx]&&1);
            if (_del1) {
                too_short += 1;
                _codon_in_reference += 1;
                ref_cleaned * (reference [_rcidx]&&1);
                qry_cleaned * (query [_rcidx]&&1);
            } else {
                _del1 = query [_rcidx] != (query [_rcidx]&&1);
                if (_del1) {
                    if (_seqL-_rcidx < 3 && _codon_in_reference % 3 == 0) {
                        break;
                    }
                    too_long += 1;
                } else {
                    ref_cleaned * (reference [_rcidx]&&1);
                    qry_cleaned * (query [_rcidx]&&1);
                    span += 1;
                    _codon_in_reference +=1;
                }
            }
        }
        ref_cleaned * 0; qry_cleaned * 0;

        return {"REF": ref_cleaned, "QRY": qry_cleaned, "TOO_SHORT" : too_short, "TOO_LONG": too_long, "SPAN": span, "OFFSET_AA" :  offset_nuc$3 + (offset_nuc % 3 > 0),"OFFSET" :  offset_nuc, "AA" : alignments.TranslateCodonsToAminoAcids (qry_cleaned, (3-offset_nuc%3)%3, code), "AA_REF" : alignments.TranslateCodonsToAminoAcids (ref_cleaned, (3-offset_nuc%3)%3, code)};
    }
    
    // -------------------------------------------------------------------------- //

    lfunction correctReadUsingCodonAlignedData (ref, qry, code) {
        reference_shifts = computeCorrection(ref);

        /*reference_shifts is the starting,ending nucleotide on the reference relative to the read. if reference is longer than the read, then both are 0*/

        offsetFrom = (qry$"^\\-+")[1]+1;
        offsetTo   = (qry$"\\-+$")[0]-1;

        /* the $ looks for the regular expression in bestAl[2] and returns a 2x1 array with the starting and ending 0-based positions of the regular expression. in this case multiple indels, -. returns -1 for both if the regular expression is not found.
            i.e. 0-based index leading indels start at (bestAl[2]$"^\\-+")[0] and end at (bestAl[2]$"^\\-+")[1]; trailing indels start at (bestAl[2]$"\\-+$")[0] and end at (bestAl[2]$"\\-+$")[0];

            so offSetFrom to offSetTo will return the reference sequence co-ordinates overlapping with the read.
        */


        if (offsetTo < 0) {
            offsetTo = Abs(qry)-1; /*if no trailing indels then to end of read*/
        }

        // check to see if the prefix in REF has some out-of-frame indels so that we can start offsetFrom in the correct frame */

        if (offsetFrom > 0) {
            frame_skips = 0;
            for (i = 0; i < offsetFrom; i += 1) {
               if ((ref[i] && 1) != ref[i]) {
                    frame_skips += 1;
               }
            }
            offsetFrom += -frame_skips;
        }

        seqOffset  = offsetFrom;          /*set the offset of the read relative to the reference. ie the number of indels needed on the read to align to the reference */
        offsetFrom +=  reference_shifts[0];           /*if the read starts before the reference then shift to start of reference ie. by reference_shifts[0] */
        offsetTo    =  offsetTo	- reference_shifts[1];           /*if the read extends beyond the reference then shift to end of reference ie. by reference_shifts[1] */

        theSeq     = qry;
        theSeq	   = qry[reference_shifts[0]][Abs(qry)-reference_shifts[1]-1]; /*the nucleotide sequence of the read that overlaps with the reference sequence */

        nucSeq	   = qry[offsetFrom][offsetTo]; /*read sequence pruned to exactly overlapping region*/
        nucSeqRef  = ref[offsetFrom][offsetTo]; /*reference sequence pruned to exactly overlapping region*/

        extra_keys = {};

        if (reference_shifts[0] > 0) { // 'qry' has a prefix that's not aligned to the reference
            extra_keys['PREFIX'] = qry[0][reference_shifts[0]-1];
        }
        if (reference_shifts[1] > 0) {
            l = Abs (qry);
            extra_keys['SUFFIX'] = qry[l-reference_shifts[1]][l-1];
        }


        return utility.Extend (igg_alignment_cleanup (nucSeqRef, nucSeq,seqOffset, code), extra_keys);
    }


    // -------------------------------------------------------------------------- //

    
    lfunction align_sequence_to_reference_set (seq, references, alignment_settings, seq_name) {


        input     = {1,2};
        input [1] = seq;

        overall = {'SCORE' : -1e100, 'REFNAME' : ''};
        
        for (id, seq; in; references) {
            input [0] = seq;
            AlignSequences (result, input, alignment_settings);

            result = result[0];

            if (result [0] >= overall['SCORE']) {

                if (alignment_settings ["REPORT_VARIANTS"] && result [0] == overall['SCORE'] && Abs (overall['REFNAME'])) {
                    overall['REFNAME'] += "|" + id;
                } else {
                    overall['REFNAME'] = id;
                }
                overall['SCORE'] = result[0];
                overall['RAW-REF'] = result[1];
                overall['RAW-QRY'] = result[2];
            }
        }


        if (Type (overall['RAW-REF']) == "String") {

            computed_score = (overall["SCORE"] - 30 * alignment_settings["MATCH"] * Exp (-Abs(seq)/3) ) / Abs (seq) * 3 ;

            if (alignment_settings["E"] <= computed_score) {
                utility.Extend (overall, correctReadUsingCodonAlignedData (overall['RAW-REF'], overall['RAW-QRY'], alignment_settings["code"]));


                if (alignment_settings["SEQ_ALIGN_CODON_ALIGN"] == TRUE && overall["SPAN"] <= 3) {
                    //assert (0, "Internal error in align_sequence_to_reference_set" + overall + "\nComputed score `computed_score`; expected score " + alignment_settings["E"] + "; match score " +  alignment_settings["MATCH"] + "\nInput sequence: `seq`");
                    return None;
                }
                return overall;
            }
        }


        return None;
     }
}

