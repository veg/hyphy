# Standard Codon Substitution Models in HyPhy

This document provides a comprehensive overview of the 23 standard codon substitution models available in HyPhy, detailing their mathematical formulations, parameterizations, frequency options, and scientific references.

---

## 1. Mechanistic Models: MG94 vs. GY94

Codon substitution models describe the evolutionary process of protein-coding DNA sequences. They are typically defined by a 61x61 rate matrix Q (excluding stop codons under the selected genetic code). The two foundational formulations are **Muse-Gaut (1994)** and **Goldman-Yang (1994)**.

### Muse-Gaut 94 (MG94)
In the MG94 formulation, rate matrix entries are scaled by the frequency of the target nucleotide at the changing codon position. For a substitution from codon i to codon j:

```
q_ij = alpha * pi_nj,  if i -> j is a synonymous transition/transversion
q_ij = beta * pi_nj,   if i -> j is a non-synonymous transition/transversion
q_ij = 0,              if i -> j requires 2 or more nucleotide substitutions
```

where alpha is the synonymous rate, beta is the non-synonymous rate, and pi_nj is the frequency of the target nucleotide.

### Goldman-Yang 94 (GY94)
In the GY94 formulation, rate matrix entries are scaled by the frequency of the target codon (pi_j). Transitions are scaled by a transition/transversion ratio (kappa):

```
q_ij = mu * pi_j,                 synonymous transversion
q_ij = mu * kappa * pi_j,          synonymous transition
q_ij = mu * omega * pi_j,          non-synonymous transversion
q_ij = mu * kappa * omega * pi_j,  non-synonymous transition
q_ij = 0,                         2 or more nucleotide substitutions
```

where omega = dN/dS, kappa is the transition/transversion ratio, and pi_j is the target codon frequency.

---

## 2. Standard Codon Model Directory

Below is the complete list of the 23 standard codon models available in HyPhy, organized by their core methodology.

### 2.1 Standard & Custom Mechanistic Models

These models build on the nucleotide substitution background (e.g., GTR, HKY85) and extend it to codon space using either MG94 or GY94 structures.

| Model Name | Background Nucleotide Model | Frequency Parameterization | Description & Key Options |
|:---|:---|:---|:---|
| **MG94** | HKY85 / GTR | F3x4 / F1x4 | Standard Muse-Gaut 94 model with selection parameter R = dN/dS. |
| **GY94** | HKY85 / GTR | F61 | Standard Goldman-Yang 94 model with transition/transversion ratio kappa and dN/dS ratio omega. |
| **MG94W9** | 9-parameter transition matrix | F3x4 | MG94 with separate exchange rates for transition/transversion classes. |
| **GY94W9** | 9-parameter transition matrix | F61 | GY94 with separate exchange rates for transition/transversion classes. |
| **MG94CUSTOM** | Custom (User-specified) | F3x4 | MG94 with custom GTR-like rate-designation string (e.g., 010010 for HKY85). |
| **GY94CUSTOMF3X4** | Custom (User-specified) | F3x4 | GY94 using position-specific nucleotide frequencies (F3x4) instead of full F61 frequencies. |
| **MG94CUSTOMFREQS** | Custom (User-specified) | User-defined | MG94 with custom nucleotide rates and customizable frequency estimators. |
| **MG94CUSTOMF1X4** | Custom (User-specified) | F1x4 | MG94 with a single set of nucleotide frequencies across all codon positions. |
| **MG94CUSTOMCF3X4**| Custom (User-specified) | CF3x4 | MG94 with corrected position-specific nucleotide frequencies. |

### 2.2 Synonymous & Non-synonymous Rate Heterogeneity

| Model Name | Rate Variation Type | Description & Key Options |
|:---|:---|:---|
| **MG94X2** | Synonymous and Non-synonymous | Allows independent or correlated discrete/gamma distributions for both synonymous and non-synonymous rates across sites or branches. Options include:<br>- `rate-variation` (`Syn:Gamma, Non-syn:Gamma`, `Independent Discrete`, etc.)<br>- `syn-rate-classes` (number of synonymous categories)<br>- `non-syn-rate-classes` (number of non-synonymous categories)<br>- `rate-designation` (custom nucleotide bias matrix string) |

### 2.3 Empirical Codon Models (ECM)

Empirical codon models use substitution rates derived from large genomic databases rather than calculating them strictly from nucleotide biases.

| Model Name | Frequency Scaling | Description & Key Options |
|:---|:---|:---|
| **ECM** | Unscaled | Fits the raw Empirical Codon Model (based on Kosiol et al. 2007). Prompted for:<br>- `multiple-substitutions` (`Yes`/`No` to allow double/triple changes) |
| **ECM+F** | F61 | Empirical Codon Model scaled by target codon frequencies. |
| **ECM+F+OMEGA** | F61 + dN/dS | Scales empirical non-synonymous substitution rates by an estimated dN/dS parameter (omega). |
| **ECM+MLFREQS** | ML-estimated | ECM with frequencies estimated via Maximum Likelihood. |
| **ECM+F3X4** | F3x4 | ECM using position-specific nucleotide frequencies (F3x4). |

### 2.4 Biophysical Property & Mixed Models

| Model Name | Model Class | Description & Key Options |
|:---|:---|:---|
| **LCAP** | Biophysical | **Linear Combination of Amino Acid Properties**. Parameterizes non-synonymous rates using physical property differences between amino acids (Chemical Composition, Polarity, Volume, Iso-electric point, and Hydropathy) with weights (alpha_0 to alpha_4) estimated from the data.<br>Options:<br>- `freq-type` (`F1x4`, `F3x4`, `F61`, `ML3x4`) |
| **MEC** | Mixed | **Mixed Empirical-Mechanistic model**. Combines site-specific mechanistic rates with empirical exchangeabilities to model protein evolution under positive selection.<br>Options:<br>- `aa-model` (e.g., `WAG`, `JTT`, `LG`, `Dayhoff`) |

### 2.5 Mechanistic Amino Acid & Empirical Exchangeability Models

These models weigh codon substitution rates using empirical amino acid exchangeability matrices.

| Model Name | Exchangeability Background | Frequencies | Description & Key Options |
|:---|:---|:---|:---|
| **MG94WAA** | AA Matrices (JTT, WAG, etc.) | F3x4 | Muse-Gaut model weighted by amino acid exchangeabilities. |
| **MG94WAAUSERFREQS**| AA Matrices | User-defined | MG94 weighted by AA exchangeabilities with custom frequency inputs. |
| **MG94WAAFREQS** | AA Matrices | Estimated | MG94 weighted by AA exchangeabilities with ML-estimated frequencies. |
| **MG94WAAF61** | AA Matrices | F61 | MG94 weighted by AA exchangeabilities with F61 frequencies. |
| **MG94WAAF61MULTIPLE**| AA Matrices | F61 | MG94 weighted by AA exchangeabilities allowing multiple substitutions. |
| **MG94WEX** | Empirical AA Exchange | F3x4 | MG94 weighted by database-derived empirical AA exchange rates. |

---

## 3. Key Configuration Options

### Frequency Parameterizations
- **F1x4**: Uses nucleotide frequencies averaged across all three codon positions (4 parameters).
- **F3x4**: Uses separate nucleotide frequencies for each of the three codon positions (12 parameters).
- **CF3x4**: Position-specific nucleotide frequencies corrected for stop-codon exclusions.
- **F61**: Estimates the frequencies of all 61 sense codons directly (60 parameters).

### Nucleotide Rate Designation Strings
Custom GTR-like models prompt for a 6-character rate designation string defining which exchange rates are shared:
- `010010`: HKY85 model (equivalent to 2 parameters: transition/transversion).
- `012345`: Full GTR model (equivalent to 6 independent exchange parameters).

---

## 4. References

1. **MG94**: Muse, S. V., & Gaut, B. S. (1994). A likelihood approach for comparing synonymous and nonsynonymous nucleotide substitution rates, with application to the chloroplast genome. *Molecular Biology and Evolution*, 11(5), 715-724. [https://doi.org/10.1093/oxfordjournals.molbev.a040152](https://doi.org/10.1093/oxfordjournals.molbev.a040152)
2. **GY94**: Goldman, N., & Yang, Z. (1994). A codon-based model of molecular evolution and its application to estimating evolutionary rates in proteins. *Molecular Biology and Evolution*, 11(5), 725-736. [https://doi.org/10.1093/oxfordjournals.molbev.a040153](https://doi.org/10.1093/oxfordjournals.molbev.a040153)
3. **ECM**: Kosiol, C., Goldman, N., & Buttimore, N. H. (2007). A new empirical model of codon substitution. *Molecular Biology and Evolution*, 24(7), 1464-1479. [https://doi.org/10.1093/molbev/msm064](https://doi.org/10.1093/molbev/msm064)
4. **MEC**: Doron-Faigenboim, A., & Pupko, T. (2007). A combined empirical-mechanistic method for estimating positive selection. *Systematic Biology*, 56(4), 627-638. [https://doi.org/10.1080/10635150701546231](https://doi.org/10.1080/10635150701546231)
5. **LCAP**: Conant, G. C., Wagner, G. P., & Stadler, P. F. (2007). Modeling amino acid substitution patterns in orthologous and paralogous genes. *Molecular Phylogenetics and Evolution*, 42(2), 298-307. [https://doi.org/10.1016/j.ympev.2006.07.009](https://doi.org/10.1016/j.ympev.2006.07.009)
6. **Custom GTR Models & Site Selection**: Kosakovsky Pond, S. L., & Frost, S. D. (2005). Not so different after all: a comparison of methods for detecting amino acid sites under selection. *Molecular Biology and Evolution*, 22(5), 1208-1222. [https://doi.org/10.1093/molbev/msi105](https://doi.org/10.1093/molbev/msi105)
7. **Rate Variation**: Yang, Z. (2000). Complexity of the evolutionary process at the codon level. *Molecular Biology and Evolution*, 17(10), 1600-1609. [https://doi.org/10.1093/oxfordjournals.molbev.a026258](https://doi.org/10.1093/oxfordjournals.molbev.a026258)
