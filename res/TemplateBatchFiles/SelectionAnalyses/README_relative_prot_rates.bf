## Site-wise relative rate estimator for protein multiple sequence alignments

_Written by Sergei L Kosakovsky Pond [spond@temple.edu] and Stephanie J Spielman_

> 2017-01-25: v0.1. Initial release.


### Motivation

This analysis performs a "non-paramteric" estimation of site-level substitution rates in a multiple sequence alignment of protein sequences. This allows one to evaluate levels of substitutional rate heterogeneity and, by extension, conservation. This analysis is based on the ["Rate4Site" method](http://www.tau.ac.il/~itaymay/cp/rate4site.html).

### Analysis workflow

#### Input

1. A protein sequence alignment (**file**)
2. A phylogenetic tree

#### Output

1. **standard output**: a MarkDown file (see sample.md)
2. **file.json**: a JSON object (see sample.json), which contains the following keys
	* `Relative site rate estimates`: for each site, a record like <pre>"1":{
     "LB":0.9712850593725352,
     "MLE":1.343044028469595,
     "UB":1.821044718637831
    }</pre> is provided. **MLE** is the point estimate of the relative rate at the site, **UB** and **LB** are the upper and lower bounds, respectively, of the profile likelihood confidence interval.
    * `alignment`: file path for the alignment used to infer the rates, e.g. _/Users/sergei/Dropbox/Work/Collaborations/rates4sites/sim178.fasta_
    * `analysis`: an object describing the version of the analysis run

#### Procedure

1. Fit a protein model of sequence evolution to the entire alignment to obtain branch lengths.
2. For each site, fixing all the other model parameters, estimate site level scaler for branch lengths: **r<sub>i</sub>**, i.e. for site **i**, the following relationship holds, for each branch **b**: **length(b | data @ site i) = r<sub>i</sub> length (b | entire alignment)**
3. The MLE for **r<sub>i</sub>**, along with a profile likelihood confidence interval, is obtained.

#### Features

* MPI Enabled