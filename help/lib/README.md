### AnalyzeCodonData (ACD): standard codon model fits

This analysis fits standard codon model templates (such as Muse-Gaut 94 or Goldman-Yang 94) to a codon alignment. It provides an interface to estimate evolutionary rates, compare different branch length scaling options (unconstrained maximum likelihood estimates vs. branch lengths proportional to an input tree), and fit either global or local branch-specific parameters.

This analysis requires an in-frame **codon** alignment and a **phylogenetic tree** (either embedded in the alignment file or provided separately).

The output is a `JSON` file containing details of the analysis description, input metadata, optimized tree topology with branch lengths, and the model fit results (including log-likelihood, parameter counts, and parameter estimates).

---

### Model and Constraint Options

The `AnalyzeCodonData` analysis supports fitting models under different options and constraints, comparing branch length parameters or rate parameters.

#### Substitution Model Templates

1. **MG94** (Muse-Gaut 94): Substitution rates depend on nucleotide mutational biases and a codon-level selection parameter.
2. **GY94** (Goldman-Yang 94): Substitution rates include an explicit transition/transversion ratio parameter (`kappa`).

#### Model Options

1. **Global**: A single alignment-wide selection parameter dN/dS (`R`) is inferred across all branches.
2. **Local**: Each branch gets its own independent selection parameter dN/dS (represented by local synonymous and non-synonymous rates).

#### Branch Length Scaling Constraints

1. **Estimate**: All branch lengths are estimated independently by Maximum Likelihood.
2. **Proportional to input tree**: Branch lengths are constrained to be proportional to those in the input tree. A single global scaling parameter (`treeScaler`) is estimated to scale the entire tree.

#### Summary of Models

For a tree with B branches:

| Model | Branch Length Constraint | Selection Rate Parameter (dN/dS) | # Parameters |
| :--- | :--- | :--- | :--- |
| **MG94 Global (Proportional)** | treeScaler * ref | Single global `R` | 2 |
| **MG94 Global (Estimate)** | Independent (B parameters) | Single global `R` | B + 1 |
| **MG94 Local (Estimate)** | Independent (B parameters) | Independent per branch | 2B |

---

### Example 1: Standard Fit with Independent Branch Lengths (Estimate)

Below is an example invocation of the `AnalyzeCodonData` analysis using the default `MG94` model with `Global` dN/dS estimation, where branch lengths are fully estimated by Maximum Likelihood. We use the rhodopsin codon alignment test dataset:

```bash
hyphy res/TemplateBatchFiles/AnalyzeCodonData.bf \
  --alignment tests/data/yokoyama.rh1.cds.mod.1-990.nex \
  --model MG94 \
  --model-options Global \
  --use-tree Yes \
  --branch-lengths Estimate \
  --output results/yokoyama.ACD.json
```

#### Output Tables printed to stdout:

#### Model Fit Summary

|        Statistic        |          Value          |
|:-----------------------:|:-----------------------:|
|     Log Likelihood      |       -13667.3145       |
|  Estimated Parameters   |           71            |
|           AIC           |       27476.6289        |
|           BIC           |       28004.6331        |

#### Global Parameter Estimates

|   Parameter    |    Estimate    |
|:--------------:|:--------------:|
|       R        |    0.124496    |

*(Note: `R` is the global dN/dS selection ratio, indicating strong purifying selection across the alignment).*

#### Branch Lengths (Subset of 70 branches shown)

|     Branch     |     Length     |
|:--------------:|:--------------:|
|      EELA      |    0.096183    |
|    CONGERA     |    0.050266    |
|     Node2      |    0.064208    |
|      EELB      |    0.053518    |
|    CONGERB     |    0.068551    |
|     ...        |    ...         |
|     BOVINE     |    0.050857    |
|    DOLPHIN     |    0.052454    |

---

### Example 2: Fit with Proportional Branch Length Constraint

Here is the same analysis run on the same dataset, but constraining the branch lengths to be proportional to the reference tree:

```bash
hyphy res/TemplateBatchFiles/AnalyzeCodonData.bf \
  --alignment tests/data/yokoyama.rh1.cds.mod.1-990.nex \
  --model MG94 \
  --model-options Global \
  --use-tree Yes \
  --branch-lengths "Proportional to input tree" \
  --output results/yokoyama.ACD.proportional.json
```

#### Output Tables printed to stdout:

#### Model Fit Summary

|        Statistic        |          Value          |
|:-----------------------:|:-----------------------:|
|     Log Likelihood      |       -13713.5641       |
|  Estimated Parameters   |            2            |
|           AIC           |       27431.1283        |
|           BIC           |       27446.0016        |

#### Global Parameter Estimates

|   Parameter    |    Estimate    |
|:--------------:|:--------------:|
|   treeScaler   |    1.153235    |
|       R        |    0.123803    |

---

### Hypothesis Testing: Constrained vs. Unconstrained Branch Lengths

Because the **Proportional** branch length model is nested within the **Unconstrained (Estimate)** model, they can be compared using a Likelihood Ratio Test (LRT) to determine if the relative branch lengths in the input tree are a good fit for this specific alignment:

| Model | Log Likelihood | # Parameters | AIC |
| :--- | :--- | :--- | :--- |
| **Proportional** | -13713.5641 | 2 | 27431.1283 |
| **Unconstrained** | -13667.3145 | 71 | 27476.6289 |

#### Likelihood Ratio Test (LRT) Results

- **Null Hypothesis (H0)**: Relative branch lengths are proportional to the input tree.
- **Alternative Hypothesis (HA)**: Branch lengths are unconstrained.
- **Test Statistic (D)**: 2 * (ln L_HA - ln L_H0) = 2 * (-13667.3145 - (-13713.5641)) = 92.4992
- **Degrees of Freedom (df)**: 71 - 2 = 69
- **p-value**: `0.0311`

#### Interpretation

The Likelihood Ratio Test yields a p-value of `0.0311`. Since p < 0.05, we reject the null hypothesis at the 5% significance level, suggesting that the relative branch lengths of the rhodopsin gene tree deviate significantly from the input tree (which may reflect gene-specific rate variation or differences in evolutionary dynamics).

---

### Key Analysis Options

#### `--code`
The genetic code table to use for translation. Default is `Universal`. Options include:
- `Universal` (Standard nuclear code)
- `Vertebrate-mtDNA`
- `Yeast-mtDNA`
- `Mold-mtDNA`
- `Invertebrate-mtDNA`
- `Ciliate-Nuclear`
- `Echinoderm-mtDNA`
- `Euplotid-Nuclear`
- `Alternative-Flatworm-mtDNA`
- `Blepharisma-Nuclear`

#### `--alignment`
Path to the input file containing the in-frame codon sequence alignment. Supported formats include FASTA, PHYLIP, Nexus, etc.

#### `--use-tree`
Whether to use the tree definition embedded in the alignment file (if one exists).
- `Yes` (default)
- `No`

#### `--tree`
Path to a separate Newick tree file. This option is used when `--use-tree` is set to `No` or when the alignment file does not contain a tree. If you want to automatically build a Neighbor-Joining tree from the alignment data, you can pass:
`neighbor-joining`

#### `--model`
The standard substitution model template.
- `MG94` (default)
- `GY94`

#### `--model-options`
Determines if parameters are globally shared or branch-local.
- `Global` (default)
- `Local`

#### `--branch-lengths`
The scaling constraint style for branch lengths.
- `Estimate` (default, all branch lengths estimated independently by Maximum Likelihood)
- `Proportional to input tree` (constrains relative branch lengths to match the input tree, scaling them with a single `treeScaler` parameter)

#### `--output`
Path to write the output JSON results. If not specified, the output is saved to the same directory as the alignment file with a `.json` suffix.
