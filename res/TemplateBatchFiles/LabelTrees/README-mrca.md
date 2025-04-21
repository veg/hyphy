## label-mrca -- label most recent common ancestors in a given tree.

Given an input Newick tree, and a list of **pairs** of taxa (tips), some of which may be absent from the tree, identify all most recent common ancestors (MRCAs) and label them. Existing labels will be **preserved**.

## Invokation

This analysis has three **required** arguments

- `--tree` the file containing the Newick-format tree to annotate
- `--output` the file where the labelled tree will be written
- `--taxa` a semi-colon separated list of pairs of taxa


### Complete options list



```
Available analysis command line options
---------------------------------------
Use --option VALUE syntax to invoke
If a [reqired] option is not provided on the command line, the analysis will prompt for its value
[conditionally required] options may or not be required based on the values of other options

tree [required]
	The tree to annotate (Newick format)

reroot
	Reroot the tree on this node ('None' to skip rerooting)
	default value: None

taxa [required]
	Find and label MRCA of these taxa pairs (NOT case sensitive); for example to label the MRCAs of 'human' and 'chimp' and of 'mouse' and 'rat', specify 'human,chimp;mouse,rat

label
	Use the following label for annotation
	default value: Foreground

label-tips
	If not 'None', use this tag to label selected tips
	default value: None

output [required]
	Write labeled Newick tree to

```


## Example runs

Results can be viewed using [Phylotree.js](phylotree.hyphy.org)

1. Label MRCA for two pairs of taxa

```
hyphy label-mrca.bf --tree data/unlabeled.nwk --taxa "Arctopsyche_sp, Psychomyia_flavida;Alloecella_grisea,Molanna_blenda" --output data/labeled-mrca-1.nwk
```

![MRCA-1](figures/mrca-1.png)

2. Label MRCA for two pairs of taxa and also label the original taxa

```
hyphy label-mrca.bf --tree data/unlabeled.nwk --label-tips Source --taxa "Arctopsyche_sp, Psychomyia_flavida;Alloecella_grisea,Molanna_blenda" --output data/labeled-mrca-2.nwk
```

![MRCA-2](figures/mrca-2.png)

