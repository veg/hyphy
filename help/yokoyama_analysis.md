# Codon Substitution Analysis: Yokoyama Rhodopsin Dataset

This report provides an interpretive comparison of eight different **global** codon substitution model configurations fitted to the Yokoyama vertebrate rhodopsin dataset (`tests/data/yokoyama.rh1.cds.mod.1-990.nex`). 

The dataset consists of 36 sequences and 71 branches (990 nucleotides/330 codons in length).

---

## 1. Executive Summary

- **Best Fit Model**: **MG94CUSTOM Global w/HM** (GTR HMM rate variation) provides the best overall fit with an AIC-c of **26315.54**, outperforming the biophysical properties model (**LCAP Global**, AIC-c = 26476.38) and the empirical codon model (**ECM+F+OMEGA Global**, AIC-c = 26676.22).
- **HMM Rate Variation**: Modelling site-to-site rate variation via a Hidden Markov Model with 4 Gamma-distributed rate categories yields a massive improvement (delta AIC-c = -160.84 compared to the best model without site-to-site variation, LCAP).
- **Autocorrelation**: The estimated HMM transition parameter lambda (0.284) indicates a strong positive correlation in evolutionary rates between adjacent codon sites, reflecting structural domains of rhodopsin.
- **Custom bias matrix**: Specifying custom GTR-like nucleotide biases (e.g., HKY85-like `010010` bias in `MG94CUSTOM`) improves the fit dramatically compared to baseline MG94.
- **Selection Pressure**: All global selection estimates (ranging from 0.08 to 0.29) indicate that the rhodopsin gene is overall under strong purifying (selective) constraint.

---

## 2. Model Fit Comparison

Below is the summary table of the log-likelihood (logL), parameter count, information criteria (AIC, AIC-c, BIC), estimated tree length, and global selection parameter (R or omega) for each model:

| Model Configuration | Log Likelihood | Parameters | AIC | AIC-c | BIC | Total Tree Length | dN/dS Estimate |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **MG94CUSTOM Global w/HM**| -13083.32 | 74 | 26314.65 | **26315.54** | 26864.96 | 6.192 | 0.087 |
| **LCAP Global** | -13160.71 | 77 | 26475.42 | **26476.38** | 27048.04 | 3.878 | *N/A (Biophysical)* |
| **ECM+F+OMEGA Global** | -13266.70 | 71 | 26675.40 | **26676.22** | 27203.40 | 4.378 | 0.292 |
| **MG94CUSTOM Global** | -13425.62 | 72 | 26995.25 | **26996.09** | 27530.69 | 3.845 | 0.110 |
| **GY94CUSTOMF3X4 Global**| -13500.36 | 72 | 27144.71 | **27145.55** | 27680.15 | 4.050 | 0.079 |
| **GY94 Global** | -13626.35 | 72 | 27396.69 | **27397.54** | 27932.13 | 3.902 | 0.132 |
| **MG94 Global** | -13667.31 | 71 | 27476.63 | **27477.45** | 28004.63 | 3.794 | 0.125 |
| **MEC Global** | -15179.54 | 76 | 30511.08 | **30512.02** | 31076.27 | 2.936 | *N/A (WAG empirical)* |

---

## 3. Visualizations

### 3.1 Model Information Criteria Comparison
The delta AIC-c score measures the information loss relative to the best-fitting model. Lower delta AIC-c indicates a better model:

![AICc Comparisons](yokoyama/aicc_comparison.png)

### 3.2 Global Selection Pressure Estimates
Comparison of the estimated global selection parameters (dN/dS or R) across the global models:

![dN/dS Estimates](yokoyama/selection_estimates.png)

### 3.3 Total Tree Length Estimates
Comparison of the total estimated tree lengths (expected substitutions per codon):

![Tree Lengths](yokoyama/tree_lengths.png)

---

## 4. Key Interpretations & Discussions

### 4.1 Why GTR-HMM and LCAP Fit Best
1. **MG94CUSTOM Global w/HM (HMM Rate Variation)**:
   This model achieves the highest likelihood and lowest AIC-c because it models site-to-site rate variation and rate autocorrelation. Codon sites inside transmembrane alpha-helices or binding pockets are highly constrained, whereas outer loop residues are variable. By using a Hidden Markov Model, rates at adjacent sites are correlated. The estimated autocorrelation transition parameter lambda is 0.284 and the Gamma shape parameter alpha is 0.882, indicating strong spatial rate autocorrelation along the rhodopsin sequence.
2. **LCAP (Linear Combination of Amino Acid Properties)**:
   This model also fits very well because it scales rates according to biophysical differences between amino acids (hydrophobicity, volume, polarity, charge). Because rhodopsin is a GPCR, transmembrane hydrophobic regions do not tolerate mutations that alter polarity.
3. **ECM (Empirical Codon Model)**:
   The empirical transition matrix accounts for realistic substitution patterns and multi-nucleotide changes (double and triple mutations), which mechanistic models assume are impossible.

### 4.2 dN/dS Estimates Divergence
The dN/dS values estimated under global models differ:
- **GY94 Global** estimates omega = 0.132, whereas **MG94 Global** estimates R = 0.125.
- **MG94CUSTOM w/HM** estimates a lower selection ratio of 0.087, showing that accounting for site-to-site rate variation prevents overestimation of selection by separating site-specific constraints.
- **ECM+F+OMEGA** estimates a higher omega of 0.292 because the baseline empirical matrix already accounts for amino acid similarity.

### 4.3 Tree Length Variations
The total tree length (sum of branch lengths) represents the total expected substitutions per codon:
- **MG94CUSTOM w/HM** estimates the longest tree (6.19), because accounting for rate variation increases tree length estimates due to high-rate variable sites.
- Standard mechanistic models are highly consistent, estimating a tree length of ~3.8 to 4.0 expected substitutions per site.
- **MEC Global** estimates a much shorter tree (2.94) because it restricts substitution patterns to WAG empirical rates.
