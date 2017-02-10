LoadFunctionLibrary("libv3/convenience/matrix.bf");
LoadFunctionLibrary("libv3/UtilityFunctions.bf");

genetic_code.hyphyAAOrdering        = "FLIMVSPTAYXHQNKDECWRG";
genetic_code.alphabeticalAAOrdering = "ACDEFGHIKLMNPQRSTVWY";

genetic_code.singleAALetterToFullName = {
    "A": "Alanine",
    "C": "Cysteine",
    "D": "Aspartic Acid",
    "E": "Glutamic Acid",
    "F": "Phenylalanine",
    "G": "Glycine",
    "H": "Histidine",
    "I": "Isoleucine",
    "K": "Lysine",
    "L": "Leucine",
    "M": "Methionine",
    "N": "Aspargine",
    "P": "Proline",
    "Q": "Glutamine",
    "R": "Arginine",
    "S": "Serine",
    "T": "Theronine",
    "V": "Valine",
    "W": "Tryptophan",
    "Y": "Tyrosine",
    "X": "Stop Codon"
};

genetic_code.stop_code = 10;


/*----------------------------------------------------------------------------------------------------------*/

lfunction genetic_code.DefineCodonToAAGivenCode(code) {

    codonToAAMap = {};
    nucChars = "ACGT";

    for (p1 = 0; p1 < 64; p1 += 1) {
        codonToAAMap[nucChars[p1$16] + nucChars[p1 % 16 $4] + nucChars[p1 % 4]] = (^"genetic_code.hyphyAAOrdering")[code[p1]];
    }

    return codonToAAMap;
}

/*----------------------------------------------------------------------------------------------------------*/

lfunction genetic_code.ComputeCodonCodeToStringMap(genCode) {
    _codonMap = {};
    _nucLetters = "ACGT";
    for (_idx = 0; _idx < Columns(genCode); _idx += 1) {
        if (genCode[_idx] != ^ "genetic_code.stop_code") {
            _codonMap + (_nucLetters[_idx$16] + _nucLetters[(_idx % 16) $4] + _nucLetters[_idx % 4]);
        }
    }
    return _codonMap;
}

/*----------------------------------------------------------------------------------------------------------*/

lfunction genetic_code.ComputeCodonCodeToStringMapStop (genCode) {
	_codonMap = {};
	_nucLetters = "ACGT";
	for (_idx = 0; _idx < Columns(genCode); _idx += 1) {
		if (genCode[_idx] == ^ "genetic_code.stop_code") {
			_codonMap + (_nucLetters[_idx$16] + _nucLetters[(_idx%16)$4] + _nucLetters[_idx%4]);
		}
	}
	return _codonMap;
}

/*----------------------------------------------------------------------------------------------------------*/

lfunction genetic_code.partition_codon(codon) {
    return Eval({
        {
            codon$16, (codon % 16) $4, codon % 4
        }
    });
}

lfunction genetic_code.assemble_codon(positions) {
    return positions[0] * 16 + positions[1] * 4 + positions[2];
}

/*----------------------------------------------------------------------------------------------------------*/

lfunction genetic_code.ComputeBranchLengthStencils(genCode) {
 /*
            given a genetic code (`genCode`), computes a matrix of N x N entries (N = sense codons)
            where a value of 1 in cell (i,j) means that i <-> j is a substitution of a specific type
            (e.g. synonymous or non-synonymous), and a value of 0 is assigned to all other cells.

            returns a dictionary with 'synonymous' and 'non-synonymous' keys

            Also see inline comments

        */

    sense_codons = CountSenseCodons (genCode);
    SS = {sense_codons, sense_codons};
    NS = {sense_codons, sense_codons};

    stop_code = ^ "genetic_code.stop_code";
    codon_offset = 0;


    for (codon = 0; codon < 64; codon += 1) {
        if (genCode[codon] == stop_code) {
            codon_offset += 1;
        } else {
            aa1 = genCode [codon];
            codon_offset2 = codon_offset;
            for (codon2 = codon + 1; codon2 < 64; codon2 += 1) {
                if (genCode [codon2] == stop_code) {
                    codon_offset2 += 1;
                } else {
                    if (aa1 == genCode [codon2]) {
                        SS [codon-codon_offset][codon2-codon_offset2] = 1;
                    } else {
                        NS [codon-codon_offset][codon2-codon_offset2] = 1;
                    }
                }
            }
        }
    }

    matrix.Symmetrize(SS);
    matrix.Symmetrize(NS);

    return {"synonymous" : SS, "non-synonymous" : NS};

}

/*----------------------------------------------------------------------------------------------------------*/

lfunction genetic_code.DefineCodonToAAMapping (code) {
    codonToAAMap = {};
    nucChars = "ACGT";

    for (p = 0; p < 64; p += 1) {
        codonToAAMap[nucChars[p$16] + nucChars[p % 16 $4] + nucChars[p % 4]] = (^"genetic_code.hyphyAAOrdering")[code[p]];
    }

    return codonToAAMap;
}

/*----------------------------------------------------------------------------------------------------------*/

lfunction genetic_code.IsStop(codon, code)

/*
	given:
		 codon (a number between 0 and 63 in AAA...TTT encoding)
		 code (the genetic code)

	returns
		 whether or not the codon is a stop codon
*/

{
    return code[codon] == ^ "genetic_code.stop_code";
}

/*----------------------------------------------------------------------------------------------------------*/

lfunction genetic_code.DefineIntegerToAAMapping (code, only_sense) {
    codon_code_map = {};

    shift = 0;

    for (p = 0; p < 64; p += 1) {
        if (genetic_code.IsStop (p, code)) {
            if (only_sense) {
               shift += 1;
               continue;
            }
        }

        codon_code_map[p-shift] = (^"genetic_code.hyphyAAOrdering")[code[p]];
    }

    return codon_code_map;
}

/**
Given a genetic code (`genCode`), computes a number of per-codon (or per pair of codons)
quantities that relate to numbers of synonymous and non-synonymous sites or
substitutions.

`options` can be null, or have any of the following keys:

    `weighting-matrix` is expected to be a set of 3 4x4 matrices showing relative frequencies of
    various nucleotide->nucleotide substitutions stratified by codon position; by default they
    are all equal

    `count-stop-codons` treat mutations to stop codons as non-synonymous changes for counting purposes
    (by the default they are not counted at all)

Also see inline comments

*/
lfunction genetic_code.ComputePairwiseDifferencesAndExpectedSites(genCode, options) {

        SS = {
            64, 1
        }; // raw codon index -> # of synonymous sites     [0-3]
        NS = SS; // raw codon index -> # of non-synonymous sites [0-3]

        stop_code = ^ "genetic_code.stop_code";

        if (Type(options["weighting-matrix"]) == "AssociativeList") {
            weighting_matrix = options["weighting-matrix"];
        } else {
            equal = {
                4, 4
            }["1"];
            weighting_matrix = {};
            weighting_matrix + equal;
            weighting_matrix + equal;
            weighting_matrix + equal;
        }

        keep_stop_codons = FALSE;

        if (Type(options["count-stop-codons"]) == "Number") {
            keep_stop_codons = options["count-stop-codons"];
        }

        codon_offset = 0;

        for (codon = 0; codon < 64; codon += 1) {

            if (genCode[codon] == stop_code) {
                codon_offset += 1;
            } else {

                codon_info = genetic_code.partition_codon(codon);
                aa = genCode[codon];

                for (codon_position = 0; codon_position < 3; codon_position += 1) {
                    norm_factor = 0.0;
                    sSites = 0.0;
                    nsSites = 0.0;
                    // mutagenize 'codon' at 'codon_position'

                    copy_codon = codon_info;
                    for (new_nuc = 0; new_nuc < 4; new_nuc += 1) {
                        if (new_nuc != codon_info[codon_position]) {
                            copy_codon[codon_position] = new_nuc;
                            new_codon = genetic_code.assemble_codon(copy_codon);
                            w = (weighting_matrix[codon_position])[codon_info[codon_position]][new_nuc];
                            if (keep_stop_codons || stop_code == genCode[new_codon] == 0) {
                                if (genCode[new_codon] != aa) {
                                    nsSites += w;
                                } else {
                                    sSites += w;
                                }
                            }
                            norm_factor += w;
                        }
                    }

                    if (norm_factor > 0) {
                        SS[codon] += sSites / norm_factor;
                        NS[codon] += nsSites / norm_factor;
                    }

                }
            }
        }

        senseCodonCount = 64 - codon_offset;

        EPS = {
            senseCodonCount, senseCodonCount
        };
        EPN = EPS;
        OPS = EPS;
        OPN = EPS;
        NTP = EPS["-1"];

        empty_dict = {};

        codon_offset_1 = 0;

        all_permutations = {
            "0": {
                {
                    0, 1, 2
                }
            },
            "1": {
                {
                    0, 2, 1
                }
            },
            "2": {
                {
                    1, 0, 2
                }
            },
            "3": {
                {
                    1, 2, 0
                }
            },
            "4": {
                {
                    2, 0, 1
                }
            },
            "5": {
                {
                    2, 1, 0
                }
            }
        };

        ntp_matrix = {
            {
                0, 0, 1, 2
            } {
                0, 0, 3, 4
            } {
                0, 0, 0, 5
            } {
                0, 0, 0, 0
            }
        };

        matrix.Symmetrize(ntp_matrix);


        for (codon_1 = 0; codon_1 < 64; codon_1 += 1) {
            if (genCode[codon_1] == stop_code) {
                codon_offset_1 += 1;
                continue;
            }

            codon_info_1 = genetic_code.partition_codon(codon_1);
            aa_1 = genCode[codon_1];
            direct_index_1 = codon_1 - codon_offset_1;

            EPS[direct_index_1][direct_index_1] = SS[codon_1];
            EPN[direct_index_1][direct_index_1] = NS[codon_1];

            codon_offset_2 = codon_offset_1;


            for (codon_2 = codon_1 + 1; codon_2 < 64; codon_2 += 1) {
                if (genCode[codon_2] == stop_code) {
                    codon_offset_2 += 1;
                    continue;
                }


                codon_info_2 = genetic_code.partition_codon(codon_2);
                aa_2 = genCode[codon_2];
                direct_index_2 = codon_2 - codon_offset_2;


                path_count = 0;
                eps = 0;
                epn = 0;
                ops = 0;
                opn = 0;
                ntp = None;

                for (path = 0; path < 6; path += 1) {
                    current_codon = codon_info_1;
                    current_aa = aa_1;
                    codon_sequence = empty_dict;
                    codon_sequence + codon_1;

                    ps = 0;
                    pn = 0;

                    for (path_step = 0; path_step < 3; path_step += 1) {
                        change_index = (all_permutations[path])[path_step];
                        if (current_codon[change_index] != codon_info_2[change_index]) {
                            current_codon[change_index] = codon_info_2[change_index];
                            current_codon_index = genetic_code.assemble_codon(current_codon);
                            next_aa = genCode[current_codon_index];
                            if (next_aa == stop_code) {
                                break;
                            }
                            codon_sequence + current_codon_index;
                            if (current_aa == next_aa) {
                                ps += 1;
                            } else {
                                pn += 1;
                            }
                            current_aa = next_aa;
                        }
                    }

                    if (path_step == 3) {
                        path_count += 1;
                        path_length = Abs(codon_sequence);

                        if (path_length == 2 && ntp == None) {
                            for (position = 0; position < 3; position += 1) {
                                if (codon_info_1[position] != codon_info_2[position]) {
                                    ntp = ntp_matrix[codon_info_1[position]][codon_info_2[position]];
                                    break;
                                }
                            }
                        }

                        pes = 0;
                        pns = 0;
                        for (path_step = 0; path_step < path_length; path_step += 1) {
                            pes += SS[codon_sequence[path_step]];
                            pns += NS[codon_sequence[path_step]];
                        }
                        eps += pes / path_length;
                        epn += pns / path_length;
                        ops += ps;
                        opn += pn;
                    }
                }

                if (path_count > 0) {
                    EPS[direct_index_1][direct_index_2] = eps / path_count;
                    EPN[direct_index_1][direct_index_2] = epn / path_count;
                    OPS[direct_index_1][direct_index_2] = ops / path_count;
                    OPN[direct_index_1][direct_index_2] = opn / path_count;
                    if (None != ntp) {
                        NTP[direct_index_1][direct_index_2] = ntp;
                    }
                }

            }
        }

        matrix.Symmetrize(EPS);
        matrix.Symmetrize(EPN);
        matrix.Symmetrize(OPS);
        matrix.Symmetrize(OPN);
        matrix.Symmetrize(NTP);

        SS_sense = {senseCodonCount, 1};
        NS_sense = {senseCodonCount, 1};

        codon_offset_1 = 0;
        for (codon_1 = 0; codon_1 < 64; codon_1 += 1) {
            if (genCode[codon_1] == stop_code) {
                codon_offset_1 += 1;
            }
            SS_sense [codon_1 - codon_offset_1] = SS[codon_1];
            NS_sense [codon_1 - codon_offset_1] = NS[codon_1];
        }

        return {"EPS" : EPS, "EPN": EPN, "OPS" : OPS, "OPN" : OPN, "NTP" : NTP, "SS" : SS_sense, "NS": NS_sense};
}
