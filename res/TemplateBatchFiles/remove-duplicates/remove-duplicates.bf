RequireVersion ("2.5.0");

LoadFunctionLibrary     ("libv3/tasks/alignments.bf");
LoadFunctionLibrary     ("libv3/tasks/trees.bf");
LoadFunctionLibrary     ("libv3/convenience/regexp.bf");

filter.analysis_description = {terms.io.info :
                            "
                            Read an alignment (and, optionally, a tree) remove duplicate sequences, and prune the tree accordingly.
                            ",
                            terms.io.version :          "0.2",
                            terms.io.reference :        "TBD",
                            terms.io.authors :          "Sergei L Kosakovsky Pond",
                            terms.io.contact :          "spond@temple.edu",
                            terms.io.requirements :     "An MSA and, optionally, a tree. v0.2 adds support for --preserve."
                          };


io.DisplayAnalysisBanner (filter.analysis_description);

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", FALSE);

KeywordArgument                     ("msa", "The MSA to remove duplicate sequences from");
SetDialogPrompt                     ("The MSA to remove duplicate sequences from");

DataSet filter.dataset              = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter filter.input          = CreateFilter (filter.dataset,1);

KeywordArgument     ("preserve", "If desired, specify a comma-separated list of sequences to use representatives of duplicate sequence clusters", "None");
filter.preserve = io.PromptUserForString ("A comma-separated list of sequences to use representatives of duplicate sequence clusters (or None)");

if (filter.preserve == "None") {
    filter.preserve = None;
} else {
    filter.temp = {};
    for (k; in; regexp.Split (filter.preserve, " *, *")) {
        filter.temp [k&&1] = 1;
    }
    filter.preserve = filter.temp;
    
}

console.log ("> Loaded an alignment with `filter.input.species` sequences and `filter.input.sites` sites from `LAST_FILE_PATH`");

if (None == filter.preserve) {
    filter.unique_count = alignments.CompressDuplicateSequences ("filter.input", "filter.datafilter.unique", FALSE);
} else {
    GetDataInfo (filter.duplicate_info, filter.input, -2);
    filter.unique_count = filter.duplicate_info["UNIQUE_SEQUENCES"];
    filter.uniques = {};
    filter.duplicate_clusters = {};
    filter.seq_names = alignments.GetSequenceNames ("filter.input");
    for (i,k; in; filter.duplicate_info["SEQUENCE_MAP"]) {
        if (filter.uniques / k == FALSE) {
            filter.uniques [k] = i;
        } else {
            if (filter.preserve [filter.seq_names[i]&&1]) {
                console.log ("\nSelecting **`filter.seq_names[i]`** to represent cluster `k`");
                filter.uniques [k] = i;
            }
        }
    }
    filter.uniques = Join (",",filter.uniques);  
    DataSetFilter filter.datafilter.unique = CreateFilter (filter.input, 1, "",filter.uniques);
}

console.log ("\nThere are **`filter.unique_count`** unique sequences in alignment ");

if (filter.unique_count == filter.input.species) {
    console.log ("### No duplicate sequences found");
}

KeywordArgument     ("tree", "An optional tree file to trim", "None");
filter.tree = io.PromptUserForString ("An optional tree file to trim");

if (filter.tree != "None") {
    fscanf (filter.tree, "Raw", filter.tree_string);
    filter.tree  = trees.LoadAnnotatedTopology (filter.tree_string);

    Topology T = filter.tree[terms.trees.newick_with_lengths];

    filter.valid_names = {};
    for (n; in; alignments.GetSequenceNames ("filter.datafilter.unique")) {
        filter.valid_names [n] = TRUE;
    }
    
    filter.existing = filter.tree [terms.trees.model_map];

    filter.delete_leaves = {};
    for (k, s; in; filter.tree[terms.trees.partitioned]) {
        if (s == terms.tree_attributes.leaf) {
            if ( filter.valid_names / k == FALSE) {
                filter.delete_leaves [k] = TRUE;
            }
        }
    }

    if (utility.Array1D (filter.delete_leaves)) {
        T - utility.Keys(filter.delete_leaves);
    }
    
    function relabel_and_annotate (node_name) {
        _label = "";
        if (Abs(filter.existing [node_name]) > 0 ) {
            _label = "{" + filter.existing [node_name] + "}";
        } 
        return node_name + _label;
    }
    
    filter.tree_string = tree.Annotate ("T", "relabel_and_annotate", "{}", TRUE);
    
    utility.SetEnvVariable ("DATAFILE_TREE", filter.tree_string);
    utility.SetEnvVariable ("IS_TREE_PRESENT_IN_DATA", TRUE);
}


KeywordArgument ("output", "Write de-duplicated MSA to");
filter.path = io.PromptUserForFilePath ("Write de-duplicated MSA to");
fprintf (filter.path, CLEAR_FILE, filter.datafilter.unique);
