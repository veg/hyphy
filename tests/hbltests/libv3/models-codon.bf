LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/models/codon.bf");

function test_model_codon_mapcode() {

    genetic_code = {
        {14, 13, 14, 13, 7, 7, 7, 7, 19, 5, 19, 5, 2, 2, 3, 2, 12, 11, 12, 11, 6, 6, 6, 6, 19, 19, 19, 19, 1, 1, 1, 1, 16, 15, 16, 15, 8, 8, 8, 8, 20, 20, 20, 20, 4, 4, 4, 4, 10, 9, 10, 9, 5, 5, 5, 5, 10, 17, 18, 17, 1, 0, 1, 0}
    };

    mapped_code = models.codon.MapCode (genetic_code);
    
    expected_sense = { {"AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"} };
    expected_stop = {{"TAA", "TAG", "TGA"}};

    assert((mapped_code[terms.sense_codons])[0] == expected_sense[0], "wrong sense");
    assert((mapped_code[terms.stop_codons])[0] == expected_stop[0], "wrong sense");
    assert(utility.KeyExists(mapped_code, terms.translation_table), "translation-table key not found");

}

test_model_codon_mapcode();


