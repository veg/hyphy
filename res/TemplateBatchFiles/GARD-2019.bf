LoadFunctionLibrary ("libv3/tasks/trees.bf");
LoadFunctionLibrary ("libv3/tasks/alignments.bf");
LoadFunctionLibrary ("libv3/tasks/estimators.bf");
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
    LoadFunctionLibrary ("libv3/models/DNA/GTR.bf");
    gard.model.generator = "models.DNA.GTR.ModelDescription";
    gard.alignment = alignments.ReadNucleotideDataSet ("gard.sequences", null);
    DataSetFilter gard.filter = CreateFilter (gard.sequences, 1);
} else { 
    // TODO: implement these branches
    if (gard.data_type == terms.gard.protein) {
        gard.alignment = alignments.ReadProteinDataSet ("gard.sequences", null);
        DataSetFilter gard.filter = CreateFilter (gard.sequences, 1);
    } else {
        gard.alignment = alignments.LoadGeneticCodeAndAlignment ("gard.sequences", "gard.filter", null);
    }
}

gard.model = model.generic.DefineModel (gard.model.generator, "gard.overall_model", {"0" : "terms.global"}, "gard.filter", null); 
 
console.log (gard.fit_partitioned_model ({{100,200,300}}, gard.model, null)); 

//--------------------------------------------------------------------------------------------------------------------


/**
 * @name tgard.fit_partitioned_model
 * Given a list of partitions, specified as increasing breakpoint locations,
   fit the specified model to said partitions, using neighbor joining trees on each partition
   return LogL and IC values
   
 * @param {Matrix} parts : sorted, 0-based breakpoints, e.g.
    {{100,200}} -> 3 partitions : 0-100, 101-200, 201-end
 * @param {Dict} model : an instantiated model to be used for all partitions
 * @param {Dict/null} initial_values : if provided, use as initial values
    
 * @returns a {Dictionary} :
    terms.fit.log_likelihood -> log likelihood
    terms.fit.AICc -> small sample AIC

 */
 
lfunction gard.fit_partitioned_model (parts, model, initial_values) {

    current_index = 0;
    current_start = 0;
    part_count = utility.Array1D (parts);
    lf_components = {2 * (part_count + 1), 1};
    trees = {};
            
    for (p = 0; p < part_count; p += 1) {
        lf_components [2*p] = "gard.filter.part_" + p;
        lf_components [2*p+1] = "gard.tree.part_" + p;
        DataSetFilter ^(lf_components[2*p]) = CreateFilter (^"gard.filter", 1, "" + current_start + "-" + parts[p]);
        current_start = parts[p] + 1;
        trees[p] = trees.ExtractTreeInfo (tree.infer.NJ ( lf_components[2*p], null));
        model.ApplyModelToTree(lf_components[2 * p + 1], trees[p], {
            "default": model
        }, None);
        
    }

    lf_components [2*p] = "gard.filter.part_" + p;
    lf_components [2*p+1] = "gard.tree.part_" + p;
    DataSetFilter ^(lf_components[2*p]) = CreateFilter (^"gard.filter", 1, "" + current_start + "-" + (^"gard.filter.sites" - 1));
    trees[p] = trees.ExtractTreeInfo (tree.infer.NJ ( lf_components[2*p], null));
    model.ApplyModelToTree(lf_components[2 * p + 1], trees[p], {
        "default": model
    }, None);

    lf_id = &likelihoodFunction;
    utility.ExecuteInGlobalNamespace ("LikelihoodFunction `lf_id` = (`&lf_components`)");
    
    model_objects = {
        "gard.overall_model" : model
    };
    
    df = 0;
    if (Type(initial_values) == "AssociativeList") {
        utility.ToggleEnvVariable("USE_LAST_RESULTS", 1);
        df = estimators.ApplyExistingEstimates(&likelihoodFunction, model_objects, initial_values, None);
    }
    
    return estimators.FitExistingLF (&likelihoodFunction, model_objects);    
}


