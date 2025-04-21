RequireVersion ("2.5.15");

LoadFunctionLibrary ("libv3/tasks/trees.bf");
LoadFunctionLibrary ("libv3/tasks/alignments.bf");
LoadFunctionLibrary ("libv3/convenience/regexp.bf");
LoadFunctionLibrary ("libv3/IOFunctions.bf");

LoadFunctionLibrary     ("libv3/tasks/alignments.bf");
LoadFunctionLibrary     ("libv3/tasks/trees.bf");

labeler.analysis_description = {terms.io.info :
                            "
                            Read a tree and a sequence alignment and reconcile their content.
                            Sequences NOT in the tree will be trimmed from the alignment,
                            and tips NOT in the alignment will be trimmed from the tree.
                            Additionally, sequences that are too short (are mostly gaps), will be removed from the alignment/tree
                            ",
                            terms.io.version :          "0.1",
                            terms.io.reference :        "TBD",
                            terms.io.authors :          "Sergei L Kosakovsky Pond",
                            terms.io.contact :          "spond@temple.edu",
                            terms.io.requirements :     "An MSA and, optionally, a tree"
                          };


io.DisplayAnalysisBanner (labeler.analysis_description);

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", FALSE);
utility.SetEnvVariable ("COUNT_GAPS_IN_FREQUENCIES", FALSE);

KeywordArgument  ("tree", "The tree to annotate (Newick format)");
SetDialogPrompt  ("Load the tree to annotate (Newick format)");

SetDialogPrompt  ("Load an alignment to trim the tree to)");

labeler.tree = trees.LoadAnnotatedTopology (TRUE);
labeler.ts = labeler.tree [^"terms.trees.newick_with_lengths"];
KeywordArgument  ("exclude", "Remove these sequences from the tree (semicolon-separated; 'null' to NOT remove any sequences)", "null");
labeler.exclude = io.PromptUserForString("Remove these sequences from the tree (comma-separated; 'null' to NOT remove any sequences)");

Topology T = labeler.ts;

if (labeler.exclude != "null") {
    labeler.exclude  = utility.DictToArray(regexp.Split (labeler.exclude, ";"));
    T - labeler.exclude;
}


labeler.existing = labeler.tree [terms.trees.model_map];

KeywordArgument  ("regexp", "Use the following regular expression to select a subset of leaves", "()");
labeler.regexp = io.PromptUserForString ("Use the following regular expression to select a subset of leaves");

KeywordArgument  ("label", "Use the following label for annotation", "Foreground");
labeler.tag = io.PromptUserForString ("Use the following label for annotation");

KeywordArgument  ("reroot", "Reroot the tree on this node ('None' to skip rerooting)", "None");
labeler.reroot = io.PromptUserForString ("Reroot the tree on this node ('None' to skip rerooting)");

if (labeler.reroot != "None") {
    rerooted = RerootTree (T, labeler.reroot);
    Topology T = rerooted;
}
KeywordArgument  ("code", "Which genetic code should be used", "Universal");
KeywordArgument  ("msa", "Load a coding alignment to trim the tree to");
labeler.info = alignments.PromptForGeneticCodeAndAlignmentWarnStops ("labeler.file", "labeler.filter");

KeywordArgument ("threshold", "Trim sequences whose non-gap content fraction is less that this", .20);
labeler.threshold = io.PromptUser ("Trim sequences whose non-gap content fraction is less that this", 0.2, 0, 1, FALSE);


labeler.seq_names = alignments.GetSequenceNames ("labeler.file");
labeler.tree_names = TipName (T, -1);

console.log (">Loaded an alignment with " + utility.Array1D (labeler.seq_names) + " sequences");
console.log (">Loaded an tree with " + utility.Array1D (labeler.tree_names) + " tips");

labeler.seq_names_dict  = {};
labeler.rename_seqs     = {};
labeler.same_taxon      = {};
labeler.content_by_id   = {};

for (i = 0; i < labeler.filter.species; i+=1) {
    GetDataInfo (labeler.counter, labeler.filter, i, -1);
    GetString   (seq, labeler.filter, i);
    labeler.content_by_id [i] = +labeler.counter;
}


labeler.too_short = {};

console.log ("");

for (i,n; in; labeler.seq_names) {
    n_trimmed = n; 
    if (labeler.same_taxon / n_trimmed == FALSE) {
        labeler.same_taxon[n_trimmed] = {};
    } 
    labeler.same_taxon[n_trimmed] + i;
    labeler.rename_seqs[n] = n_trimmed;
    if (labeler.content_by_id [i] < labeler.filter.sites * labeler.threshold) {
        console.log ("- Trimming " + n + " as too short (" + labeler.content_by_id [i] + ", " + Format (100*labeler.content_by_id [i]/labeler.filter.sites, 5, 3) + "%)");      
        labeler.too_short [n] = 1;
        labeler.too_short [n && 1] = 1;
        continue;
    }
        
    labeler.seq_names_dict [n_trimmed] = 1;
    labeler.seq_names_dict [n_trimmed && 1] = 1;
}


labeler.exclude_duplicates = {};

alignments.GetSequenceByName ("labeler.filter", None);

for (k,v;in;labeler.same_taxon) {
    if (Abs (v) > 1) { 
        labeler.non_gaps = {};
        labeler.max_score = 0;
        for (seq;in;v) {
            if (labeler.content_by_id [seq] > labeler.max_score ) {
                labeler.max_score = labeler.content_by_id [seq];
                labeler.keep_me   = seq;
            }
        }
        for (seq;in;v) {
            if (seq != labeler.keep_me ) {
                labeler.exclude_duplicates [seq] = 1;
            }
        }
        
    }
}


labeler.tree_names_dict = {};
for (n; in; labeler.tree_names) {
    labeler.tree_names_dict [n] = 1;
}

labeler.include_species = {};

console.log ("");

for (k, v; in; labeler.tree_names_dict) {
    if (labeler.seq_names_dict[k]) {
        labeler.include_species[k] = 1;
        labeler.include_species[k&&1] = 1;
    } else {
        if (labeler.too_short / k == 0) {
            console.log ("- Label `k` has no matching sequence");
        }
    }
}



labeler.tree_names_dict - labeler.include_species;

if (Abs (labeler.tree_names_dict)) {
    T - Rows(labeler.tree_names_dict);
}


function labeler.includeSeq (name, seq) {
    return labeler.include_species[labeler.rename_seqs[name]];
}

DataSetFilter labeler.reduced_filter1 = CreateFilter (labeler.filter, 1, "", labeler.exclude_duplicates[speciesIndex]!=1);
DataSetFilter labeler.reduced_filter = CreateFilter (labeler.reduced_filter1, 1, "", "labeler.includeSeq");

console.log (">Writing out a reconciled filter and tree with `labeler.reduced_filter.species` sequences/tips");

for (i = 0; i < labeler.reduced_filter.species; i+=1) {
    GetString (labeler.name_i, labeler.reduced_filter, i);
    SetParameter (  labeler.reduced_filter, i, labeler.rename_seqs[labeler.name_i]);
}

//console.log ("" + labeler.reduced_filter.species + " " + TipCount (T) + " " +  Abs (labeler.include_species ) + " " + Abs (labeler.rename_seqs));


labeler.labels = {};

if (labeler.regexp == "()") {
  KeywordArgument  ("list", "Line list of sequences to include in the set (required if --regexp is not supplied)");
  SetDialogPrompt  ("Line list of sequences to include in the set");
  labeler.list_dict = io.ReadDelimitedFile (null,",",FALSE);
  labeler.list = {};
  for (k, v; in; labeler.list_dict["rows"]) {
    labeler.labels[regexp.Replace (v[0],"\\ +$", "")] = labeler.tag;
  }
} else {
  for (n; in; T) {
    if (regexp.Find (n,labeler.regexp)) {
      labeler.labels[n] = labeler.tag;
    }
  }
}

assert (utility.Array1D (labeler.labels) > 0, "A non-empty set of selected branches is required");

label.core = utility.Array1D (labeler.labels);

console.log ("\nSelected " + label.core + " branches to label\n");

KeywordArgument  ("internal-nodes", "Strategy for labeling internal nodes", "All descendants");
labeler.kind  = io.SelectAnOption (
  {
    {"None","Only assign labels to selected nodes"}
    {"All descendants","Only label an internal node if all its descendants are labeled"}
    {"Some descendants","Only label an internal node if some of its descendants are labeled"}
    {"Parsimony","Use maximum parsimony to label internal nodes"}
  },
  "Strategy for labeling internal nodes"
);

if (labeler.kind == "All descendants") {
  labeler.labels * ((trees.ConjunctionLabel ("T", labeler.labels))["labels"]);
}

if (labeler.kind == "Some descendants") {
  labeler.labels * ((trees.DisjunctionLabel ("T", labeler.labels))["labels"]);
}

if (labeler.kind == "Parsimony") {
  labeler.labels * ((trees.ParsimonyLabel ("T", labeler.labels))["labels"]);
}

console.log ("\nLabeled " + (utility.Array1D (labeler.labels) - label.core) + " additional branches\n");

KeywordArgument ("output", "Write labeled Newick tree to");
labeler.path = io.PromptUserForFilePath ("Write labeled Newick tree to");


utility.SetEnvVariable ("IS_TREE_PRESENT_IN_DATA", TRUE);
utility.SetEnvVariable ("DATAFILE_TREE", tree.Annotate ("T", "relabel_and_annotate", "{}", TRUE));

fprintf (labeler.path, labeler.reduced_filter);

function relabel_and_annotate (node_name) {
    _label = "";
    if (Abs(labeler.existing [node_name]) > 0 ) {
        _label = "{" + labeler.existing [node_name] + "}";
    } else {
        if (labeler.labels / node_name) {
            _label = "{" + labeler.labels[node_name] + "}";
        }
    }
    return node_name + _label;
}
