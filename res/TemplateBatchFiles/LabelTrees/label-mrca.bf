RequireVersion ("2.5.15");

LoadFunctionLibrary ("libv3/tasks/trees.bf");
LoadFunctionLibrary ("libv3/tasks/alignments.bf");
LoadFunctionLibrary ("libv3/convenience/regexp.bf");
LoadFunctionLibrary ("libv3/IOFunctions.bf");

LoadFunctionLibrary     ("libv3/tasks/alignments.bf");
LoadFunctionLibrary     ("libv3/tasks/trees.bf");

labeler.analysis_description = {terms.io.info :
                            "
                            Read a tree and a list of taxon pairs, identify their MRCA, and label them programmaticaly.
                            ",
                            terms.io.version :          "0.2",
                            terms.io.reference :        "TBD",
                            terms.io.authors :          "Sergei L Kosakovsky Pond",
                            terms.io.contact :          "spond@temple.edu",
                            terms.io.requirements :     "A tree to label. Existing tree labels will be respected (not overwritten)"
                          };


io.DisplayAnalysisBanner (labeler.analysis_description);

KeywordArgument  ("tree", "The tree to annotate (Newick format)");
SetDialogPrompt  ("Load the tree to annotate (Newick format)");

labeler.tree = trees.LoadAnnotatedTopology (TRUE);

labeler.ts = labeler.tree [^"terms.trees.newick_with_lengths"];
Topology T = labeler.ts;

labeler.existing = labeler.tree [terms.trees.model_map];


KeywordArgument  ("reroot", "Reroot the tree on this node ('None' to skip rerooting)", "None");
labeler.reroot = io.PromptUserForString ("Reroot the tree on this node ('None' to skip rerooting)");

if (labeler.reroot != "None") {
    rerooted = RerootTree (T, labeler.reroot);
    Topology T = rerooted;
}

// build an upper case list of taxa found in the tree

labeler.tree_info        = trees.ExtractTreeInfoFromTopology ("T");

labeler.available_leaves = {};

for (k,l; in; labeler.tree_info[terms.trees.partitioned]) {
    if (l == terms.tree_attributes.leaf) {
        labeler.available_leaves [k && 1]  = 1;
    }
}

console.log ("\nRead a tree on " + (utility.Array1D (labeler.available_leaves)) + " leaves.\n" + Join (", ", utility.Keys (labeler.available_leaves)) + "\n");

KeywordArgument  ("taxa", "Find and label MRCA of these taxa pairs (NOT case sensitive); for example to label the MRCAs of 'human' and 'chimp' and of 'mouse' and 'rat', specify 'human,chimp;mouse,rat");

labeler.mrca_pairs = io.PromptUserForString ("Find and label MRCA of these taxa pairs (NOT case sensitive); for example to label the MRCAs of 'human' and 'chimp' and of 'mouse' and 'rat', specify 'human,chimp;mouse,rat");

KeywordArgument  ("label", "Use the following label for annotation", "Foreground");
labeler.tag = io.PromptUserForString ("Use the following label for annotation");

KeywordArgument  ("label-tips", "If not 'None', use this tag to label selected tips", "None");
labeler.tip_tag = io.PromptUserForString ("If not 'None', use this tag to label selected tips");


labeler.mrca_pair_set = {};
labeler.mrca_taxa = {};
for (pair; in; regexp.Split (labeler.mrca_pairs, " *; *")) {
    labeler.pair = {};
    for (seq; in; regexp.Split (pair, " *, *") ) {
        if (labeler.available_leaves [seq&&1]) {
            labeler.pair  + (seq && 1);
        }
    }
    if (utility.Array1D ( labeler.pair ) == 2) {
        labeler.mrca_taxa [labeler.pair[0]] = 1;
        labeler.mrca_taxa [labeler.pair[1]] = 1;
        labeler.mrca_pair_set + labeler.pair ;
    }
}

assert (utility.Array1D (labeler.mrca_pair_set) > 0, "A non-empty set of pairs of terminal branches is required");

console.log ("\nWill label MRCAs for " + (utility.Array1D (labeler.mrca_pair_set)) + " pairs of taxa");

labeler.parents     = trees.ParentMap (&T);
labeler.descendants = {};
labeler.core        = {};

for (n; in; T) {   
    labeler.parent = labeler.parents[n];
    //if (labeler.parent) {
        if (labeler.descendants[labeler.parent] == 0) {
            labeler.descendants[labeler.parent] = {};
        }
        if (labeler.available_leaves [n && 1]) {
            (labeler.descendants[labeler.parent])[n && 1] = 1;
        } else {
            // check to see if any pairs of unlabeled MRCAs go into this node
            labeler.pass = {};
            
            if (utility.Array1D (labeler.descendants [n])) {
                for (i,v; in; labeler.mrca_pair_set) {
                     //console.log (labeler.descendants [n]);
                     if ((labeler.descendants [n])[v[0]] && (labeler.descendants [n])[v[1]] ) {
                        labeler.core [n] = 1;
                        labeler.pass + i;
                    }
                }
                        
                for (i; in ; labeler.pass) {
                    labeler.mrca_pair_set - i;
                }
            
                for (nd,v; in; labeler.descendants [n]) {
                    (labeler.descendants[labeler.parent])[nd] = 1;
                }
            
                labeler.descendants - n;
            }
        }
    //} 
}

console.log ("\nWill label " + (utility.Array1D (labeler.core)) + " internal branches");

if (labeler.tip_tag != "None") {
    console.log ("\nWill also label " + (utility.Array1D (labeler.mrca_taxa)) + " taxa (terminal branches)");
}


KeywordArgument ("output", "Write labeled Newick tree to");
labeler.path = io.PromptUserForFilePath ("Write labeled Newick tree to");
fprintf (labeler.path, CLEAR_FILE, tree.Annotate ("T", "relabel_and_annotate", "{}", TRUE));

function relabel_and_annotate (node_name) {
    _label = "";
    if (Abs(labeler.existing [node_name]) > 0 ) {
        _label = "{" + labeler.existing [node_name] + "}";
    } else {
        if (labeler.core  / node_name) {
            _label = "{" + labeler.tag  + "}";
        } 
        if (labeler.mrca_taxa / ( node_name && 1)) {
            if (labeler.tip_tag != "None") {
                _label = "{" + labeler.tip_tag  + "}";
            }
        }
    }
    return node_name + _label;
}
