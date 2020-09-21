LoadFunctionLibrary("../IOFunctions.bf");
LoadFunctionLibrary("../all-terms.bf");
LoadFunctionLibrary("../models/frequencies.bf");


/**
 * Compute  all pairwise distances between sequences in a data set
 * @name    distances.nucleotide.tn93
 * @param   filter {String} - id of dataset
 * @param   freqs {Matrix/null} - if not null, use these nucleotide frequencies
 * @param   options{Dict/null} - options for ambiguity treatment
 * @returns {Matrix} r - pairwise TN93 distances
 */
lfunction distances.nucleotide.tn93 (filter, freqs, options) {
    if (null == freqs) {
        freqs = (frequencies._aux.empirical.singlechar ({}, null, filter))[^"terms.efv_estimate"];
    }

    fY = freqs[1] + freqs[3];
    fR = 1 - fY;

	if (Min (freqs, 0) == 0) {
		K2P = TRUE;
	}
	else {
	    K2P = FALSE;
		K1 = 2*freqs[0]*freqs[2]/fR;
		K2 = 2*freqs[1]*freqs[3]/fY;
		K3 = 2*(fR*fY-freqs[0]*freqs[2]*fY/fR-freqs[1]*freqs[3]*fR/fY);
	}

    sequence_count = ^(filter + ".species");
    distances      = {sequence_count,sequence_count};
    resolution     = "RESOLVE_AMBIGUITIES";
    if (utility.Has (options, "ambigs","String")) {
        resolution = options["ambigs"];
    }

    for (s1 = 0; s1 < sequence_count; s1 += 1) {
        for (s2 = s1 + 1; s2 < sequence_count; s2 += 1) {
             ExecuteCommands ("GetDataInfo (count, ^filter, s1, s2, `resolution`)");
             totalSitesCompared = +count;
             d = 1000.;
             if (totalSitesCompared > 0) {            
                 count = count * (1/totalSitesCompared);
                 AG                     = count[0][2] + count[2][0];
                 CT                     = count[1][3] + count[3][1];
                 transversions          = 1 - AG - CT - count[0][0] - count[1][1] - count[2][2] - count[3][3];
                 if (K2P) {
                    d1 = 1-2*(AG+CT)-transversions;
                    d2 = 1-2*transversions;

                    if (d1 >0 && d2>0) {
                        d = -(0.5*Log(d1)+.25*Log(d2));
                    }

                 } else {
                     d1 = 1-AG/K1-0.5*transversions/fR;
                     d2 = 1-CT/K2-0.5*transversions/fY;
                     d3 = 1-0.5*transversions/fR/fY;
                     if (d1>0 && d2 >0 && d3 >0) {
                        d = -K1*Log(d1)-K2*Log(d2)-K3*Log(d3);
                     }
                }
            }
            distances [s1][s2] = d;
            distances [s2][s1] = d;
        }
    }


    return distances;
}

/**
 * Compute  all pairwise percent distances between sequences in a data set
 * @name    distances.nucleotide.p_distance
 * @param   filter  {String} - id of dataset
 * @param   options {Dict/null} - options for ambiguity treatment
 * @returns {Matrix} r - pairwise TN93 distances
 */
lfunction distances.p_distance (filter, options) {
    sequence_count = ^(filter + ".species");
    distances      = {sequence_count,sequence_count};
    resolution     = "RESOLVE_AMBIGUITIES";
    if (utility.Has (options, "ambigs","String")) {
        resolution = options["ambigs"];
    }

    if (sequence_count) {
        ExecuteCommands ("GetDataInfo (count, ^filter, 0, 0, `resolution`)");
        off_diag = count["_MATRIX_ELEMENT_COLUMN_==_MATRIX_ELEMENT_ROW_"];
    }

    expr := "GetDataInfo (count, ^filter, s1, s2, `resolution`)";
    
    for (s1 = 0; s1 < sequence_count; s1 += 1) {
        for (s2 = s1 + 1; s2 < sequence_count; s2 += 1) {
             ExecuteCommands (expr);
             totalSitesCompared = +count;
             if (totalSitesCompared != 0) {
                identicalSites = +count[off_diag];
                d = (totalSitesCompared-identicalSites) / totalSitesCompared;
             } else {
                d = 1000;
             }
             distances [s1][s2] = d;
             distances [s2][s1] = d;
        }
    }

    return distances;
}

