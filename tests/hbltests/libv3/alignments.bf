LoadFunctionLibrary("libv3/tasks/alignments.bf");

function test_read_nucleotide_alignment() {

    file_name = "`PATH_TO_CURRENT_BF`data/CD2.nex";
    hky85.nuc_data = {};
    hky85.nuc_filter = {};

    results = alignments.ReadNucleotideAlignment(file_name, "hky85.nuc_data", "hky85.nuc_filter");

    assert(results["sequences"] == 10, "parsed wrong number of sequences");
    assert(results["sites"] == 561, "parsed wrong number of sites");
    assert(utility.KeyExists(results, terms.json.partitions), "partitions key not found");

}

test_read_nucleotide_alignment();


