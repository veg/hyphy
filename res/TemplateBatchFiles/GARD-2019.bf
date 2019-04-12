LoadFunctionLibrary ("libv3/tasks/trees.bf");
LoadFunctionLibrary ("libv3/tasks/alignments.bf");
LoadFunctionLibrary ("libv3/convenience/regexp.bf");
LoadFunctionLibrary ("libv3/convenience/math.bf");
LoadFunctionLibrary ("libv3/IOFunctions.bf");
LoadFunctionLibrary ("libv3/UtilityFunctions.bf");

gard.analysis_description = {terms.io.info : "GARD : Genetic Algorithms for Recombination Detection. Implements a heuristic
approach to screening alignments of sequences for recombination, by using the CHC genetic algorithm to search for phylogenetic
incongruence among different partitions of the data. The number of partitions is determined using a step-up procedure, while the
placement of breakpoints is searched for with the GA. The best fitting model (based on c-AIC) is returned; and additional post-hoc
tests run to distinguish topological incongruence from rate-variation.",
                           terms.io.version : "0.1",
                           terms.io.reference : "**Automated Phylogenetic Detection of Recombination Using a Genetic Algorithm**, _Mol Biol Evol 23(10), 1891–1901",
                           terms.io.authors : "Sergei L Kosakovsky Pond",
                           terms.io.contact : "spond@temple.edu",
                           terms.io.requirements : "A sequence alignment."
                          };


namespace terms.gard {
    nucleotide = "Nucleotide";
    protein    = "Protein";
    codon      = "Codon";
};

namespace gard.parameters {

};

gard.master_list = {}; // "model string" -> "model fitness"

io.DisplayAnalysisBanner (gard.analysis_description);

KeywordArgument ("type",        "The type of data to perform screening on", "Nucleotide");
KeywordArgument ("code",        "Genetic code to use (for codon alignments)", "Universal", "Choose Genetic Code");
KeywordArgument ("alignment",   "Sequence alignment to screen for recombination");

gard.data_type = io.SelectAnOption  ({terms.gard.nucleotide : "A nucleotide (DNA/RNA) alignment",
                                      terms.gard.protein : "A protein alignment",
                                      terms.gard.codon : "An in-frame codon alignment"},
                                      "The type of data to perform screening on");


if (gard.data_type == terms.gard.nucleotide) {
    gard.alignment = alignments.ReadNucleotideDataSet ("gard.sequences", null);

    DataSetFilter gard.filter = CreateFilter (gard.sequences, 1);

    //console.log (distances.nucleotide.tn93 (&gard.filter, null, null));
} else {
    if (gard.data_type == terms.gard.protein) {
        gard.alignment = alignments.ReadProteinDataSet ("gard.sequences", null);
        DataSetFilter gard.filter = CreateFilter (gard.sequences, 1);
    } else {
        gard.alignment = alignments.LoadGeneticCodeAndAlignment ("gard.sequences", "gard.filter", null);
    }
}

console.log (tree.infer.NJ (&gard.filter, null));

