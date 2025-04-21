RequireVersion ("2.5.15");

LoadFunctionLibrary ("libv3/tasks/trees.bf");
LoadFunctionLibrary ("libv3/tasks/alignments.bf");
LoadFunctionLibrary ("libv3/convenience/regexp.bf");
LoadFunctionLibrary ("libv3/IOFunctions.bf");

LoadFunctionLibrary     ("libv3/tasks/alignments.bf");
LoadFunctionLibrary     ("libv3/tasks/trees.bf");

labeler.analysis_description = {terms.io.info :
                            "
                            Read a tree and (optionally) a list of sequence names, and annotate branches for subsequent analyses that accept branch partitions.
                            ",
                            terms.io.version :          "0.2",
                            terms.io.reference :        "TBD",
                            terms.io.authors :          "Sergei L Kosakovsky Pond",
                            terms.io.contact :          "spond@temple.edu",
                            terms.io.requirements :     "A tree to label. Existing tree labels will be respected, allowing for multiple sets to be labeled"
                          };


io.DisplayAnalysisBanner (labeler.analysis_description);

KeywordArgument  ("tree", "The tree to annotate (Newick format)");
SetDialogPrompt  ("Load the tree to annotate (Newick format)");

labeler.tree = trees.LoadAnnotatedTopology (TRUE);


labeler.ts = labeler.tree [^"terms.trees.newick_with_lengths"];
Topology T = labeler.ts;

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

KeywordArgument  ("invert", "Invert selection", "No");

labeler.inverse  = io.SelectAnOption (
  {
    {"No","Matching list/regexp"}
    {"Yes","NOT matching list/regexp"}
    
  },
  "Invert selection"
) == "Yes";

labeler.labels = {};

if (labeler.regexp == "()") {
  KeywordArgument  ("list", "Line list of sequences to include in the set (required if --regexp is not supplied)");
  SetDialogPrompt  ("Line list of sequences to include in the set");
  labeler.list_dict = io.ReadDelimitedFile (null,",",FALSE);
  labeler.list = {};
  for (k, v; in; labeler.list_dict["rows"]) {
    labeler.list[regexp.Replace (v[0],"\\ +$", "")] = 1;
  }
  if (labeler.inverse) {
     for (n; in; T) {
        if (labeler.list[n] == 0) {
          labeler.labels[n] = labeler.tag;
        }
      }
  } else {
      for (n; in; T) {
        if (labeler.list[n]) {
          labeler.labels[n] = labeler.tag;
        }
      }
    }
} else {
  if (labeler.inverse) {
     for (n; in; T) {
        if (None == regexp.Find (n,labeler.regexp)) {
          labeler.labels[n] = labeler.tag;
        }
      }
  } else {
      for (n; in; T) {
        if (regexp.Find (n,labeler.regexp)) {
          labeler.labels[n] = labeler.tag;
        }
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
    {"All descendants, no MRCA","Only label an internal node if all its descendants are labeled"}
    {"Some descendants","Only label an internal node if some of its descendants are labeled"}
    {"Parsimony","Use maximum parsimony to label internal nodes"}
  },
  "Strategy for labeling internal nodes"
);

KeywordArgument  ("leaf-nodes", "Strategy for labeling selected leaves", "Label");

labeler.tips  = io.SelectAnOption (
  {
    {"Label","Include tips"}
    {"Skip","Only label internal nodes using the selected strategy"}
    
  },
  "Include tips in labeling"
) == "Label";

if (labeler.tips) {

    labeler.internal_node_labels = {};
    
    if (labeler.kind == "All descendants" || labeler.kind == "All descendants, no MRCA" ) {
    labeler.internal_node_labels  =  ((trees.ConjunctionLabel ("T", labeler.labels))["labels"]);
    }

    if (labeler.kind == "Some descendants") {
      labeler.internal_node_labels = ((trees.DisjunctionLabel ("T", labeler.labels))["labels"]);
    }

    if (labeler.kind == "Parsimony") {
       labeler.internal_node_labels = ((trees.ParsimonyLabel ("T", labeler.labels))["labels"]);
    }
    
    if (labeler.kind == "All descendants, no MRCA") {
        /* 
            here, we are going to perform a POST-order traversal of the tree, and remove labels for the internal nodes that are closest to the root,
            i.e. nodes that are roots of subtrees
        */
        
        labeler.remove = {};
        
        for (n,p; in; trees.ParentMap ("T")) {
            if (p) {
                if (labeler.internal_node_labels / n) {
                    if (labeler.internal_node_labels[n] != labeler.internal_node_labels[p]) {
                        labeler.remove + n;
                    }
                }
            }
        }   
        
        for (n; in;  labeler.remove) {
            labeler.internal_node_labels - n;
        }
    }
    labeler.labels * labeler.internal_node_labels;
    

    console.log ("\nLabeled " + (utility.Array1D (labeler.labels) - label.core) + " additional internal branches\n");
} else {
    if (labeler.kind == "All descendants") {
      labeler.labels = ((trees.ConjunctionLabel ("T", labeler.labels))["labels"]);
    }

    if (labeler.kind == "Some descendants") {
      labeler.labels = ((trees.DisjunctionLabel ("T", labeler.labels))["labels"]);
    }

    if (labeler.kind == "Parsimony") {
      labeler.labels = ((trees.ParsimonyLabel ("T", labeler.labels))["labels"]);
    }
    console.log ("\nLabeled " + (utility.Array1D (labeler.labels) ) + " internal branches\n");
}

KeywordArgument ("output", "Write labeled Newick tree to");
labeler.path = io.PromptUserForFilePath ("Write labeled Newick tree to");


fprintf (labeler.path, CLEAR_FILE, tree.Annotate ("T", "relabel_and_annotate", "{}", TRUE));

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
