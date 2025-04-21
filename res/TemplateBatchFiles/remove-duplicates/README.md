## Remove duplicate sequences

This is a simple script to read in a multiple sequence alignment (and optionally a tree), remove all duplicate sequences (keep one copy), and trim the tree file accordingly. The resulting de-duplicated alignment (and trimmed tree, if provided) is saved to a new file.

### Example

```
hyphy remove-duplicates.bf --msa example.fas --tree example.nwk --output uniques.nxh
```

#### Output

Analysis Description
--------------------
 Read an alignment (and, optionally, a tree) remove duplicate sequences,
and prune the tree accordingly. 

- __Requirements__: An MSA and, optionally, a tree. v0.2 adds support for --preserve.

- __Citation__: TBD

- __Written by__: Sergei L Kosakovsky Pond

- __Contact Information__: spond@temple.edu

- __Analysis Version__: 0.2

A comma-separated list of sequences to use representatives of duplicate sequence clusters (or None) : preserve: None
> Loaded an alignment with 90 sequences and 1260 sites from /Users/sergei/Development/hyphy-analyses/remove-duplicates/example.fas

There are **11** unique sequences in alignment 
An optional tree file to trim : tree: example.nwk


### Analysis arguments

```
Available analysis command line options
---------------------------------------
Use --option VALUE syntax to invoke
If a [reqired] option is not provided on the command line, the analysis will prompt for its value
[conditionally required] options may or not be required based on the values of other options

msa [required]
	The MSA to remove duplicate sequences from

preserve
	If desired, specify a comma-separated list of sequences to use representatives of duplicate sequence clusters
	default value: None

tree
	An optional tree file to trim
	default value: None

output [required]
	Write de-duplicated MSA to

```

### Tips

To change the output format, add (in version 2.5.8 or later)

```
ENV="DATA_FILE_PRINT_FORMAT=N"
```

to the command line.

`N` is a number between 0 and 11 (4 is the default, see below). For example, to produce a FASTA output, use 

```
hyphy remove-duplicates.bf --msa example.fas --tree example.nwk 
--output uniques.fas ENV="DATA_FILE_PRINT_FORMAT=9"
```

```
kFormatMEGASequential             = 0,
kFormatMEGAInterleaved            = 1,
kFormatPHYLIPSequential           = 2,
kFormatPHYLIPInterleaved          = 3,
kFormatNEXUSLabelsSequential      = 4,
kFormatNEXUSLabelsInterleaved     = 5,
kFormatNEXUSSequential            = 6,
kFormatNEXUSInterleaved           = 7,
kFormatCharacterList              = 8,
kFormatFASTASequential            = 9,
kFormatFASTAInterleaved           = 10,
kFormatPAML                       = 11

```

### Selecting specific sequences to keep

If one wishes to use particular sequence(s) to retain as representatives of clusters of duplicate sequences, these sequences can be supplied as `--preserve seq1,seq2,...` command line arguments. For example:


```
$hyphy remove-duplicates.bf --msa example.fas --tree example.nwk --output uniques.fas ENV="DATA_FILE_PRINT_FORMAT=9"
...

$grep '>' uniques.fas


>seq_1
>seq_309
>seq_947
>seq_969
>seq_991
>seq_1343
>seq_1387
>seq_1541
>seq_1651
>seq_1783
>seq_1893

```

```
hyphy LIBPATH=/Users/sergei/Development/hyphy/res remove-duplicates.bf --msa example.fas --tree example.nwk --output uniques.fas ENV="DATA_FILE_PRINT_FORMAT=9" --preserve seq_1321

...

Selecting **seq_1321** to represent cluster 0

$grep '>' uniques.fas

>seq_1321
>seq_309
>seq_947
>seq_969
>seq_991
>seq_1343
>seq_1387
>seq_1541
>seq_1651
>seq_1783
>seq_1893
...
