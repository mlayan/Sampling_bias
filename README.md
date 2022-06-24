# Impact and mitigation of sampling bias to determine viral spread: evaluating discrete phylogeography through CTMC modeling and structured coalescent model approximations

Full text and supplementary materials are available [here]().

## Abstract

Bayesian phylogeographic inference is a powerful tool in molecular epidemiological studies that enables reconstructing the origin and subsequent geographic spread of pathogens. Such inference is, however, potentially affected by geographic sampling bias. Here, we investigated the impact of sampling bias on the spatiotemporal reconstruction of viral epidemics using Bayesian discrete phylogeographic models and explored different operational strategies to mitigate this impact. We considered the continuous-time Markov chain (CTMC) model and two structured coalescent approximations (BASTA and MASCOT). For each approach, we compared the estimated and simulated spatiotemporal histories in biased and unbiased conditions based on simulated epidemics of rabies virus (RABV) in dogs in Morocco. While the reconstructed spatiotemporal histories were impacted by sampling bias under CTMC, BASTA and MASCOT reconstructions were biased even when employing unbiased samples. This resulted from assuming incorrect population dynamics but can be overcome by accommodating time-varying population sizes. Increasing the number of analyzed genomes led to more robust estimates at low sampling bias, notably for CTMC. Alternative sampling strategies that maximize the spatiotemporal coverage greatly improved the inference at intermediate sampling bias for CTMC, and to a lesser extent, for BASTA and MASCOT. We further applied these approaches to two empirical datasets: a RABV dataset from the Philippines, and a SARS-CoV-2 dataset from its early spread across the world. In conclusion, sampling biases are ubiquitous in phylogeographic analyses but may be accommodated by increasing sample size, balancing spatial and temporal composition in the samples, and informing structured coalescent models with reliable case count data.

## Simulation study
### Overview
Passive surveillance systems are suspected to report infectious disease cases in a biased way by overrepresenting cases from urban and administratively important areas at the expense of rural and lowly populated areas, especially in low- and middle-income countries. Biased sampling of cases might have an impact on discrete phylogeographic predictions. To test the potential impact of spatial sampling bias on discrete phylogeographic inference, 
1. We simulated realistic RABV epidemics among domestic dog populations in Morocco using a stochastic metapopulation model
2. RABV whole-genomes associated with each case were simulated using a simple HKY model 
3. RABV genomes were sampled uniformly, in a biased way or in a way expected to mitigate the effects of bias    
4. The spatiotemporal history of rabies epidemics was infered using CTMC, BASTA, MASCOT, and MASCOT-GLM, all implemented in BEAST
5. Inferred and simulated migration histories were compared in biased, unbiased and mitigating scenarios

### Simulations
RABV epidemics were simulated using a discrete-time spatially-explicit stochastic model. Dog mobility was parametrized using a radiation model fitted to human population data in Morocco. Sequence evolution follows a simple HKY model without considering site heterogeneity, selection processes, varying evolutionary rate, or within-host evolution.

The human population data, the inferred mobility matrix, the Moroccan shapefiles, and the fasta file of the empirical RABV genome used for the index case are all located in `inputfiles`.

Two simulation frameworks were tested: 
- `3demes` that corresponds to epidemics across 3 Moroccan demes/regions,
- `7demes` that corresponds to epidemics across 7 Moroccan demes/regions. 

In both frameworks, approximatively 1% of infection events happen between demes/regions. 

Simulations can be launched using the `R` scripts in the `generate` subfolders of `3demes` and `7demes`. This script depends on the following R libraries: `dplyr`, `Rcpp`, `seqinr`, `doParallel`, `doRNG`, and `foreach`. Simulated transmission chains, fasta files and traits files (location of the sampled sequences) are generated and stored in the following subfolders: `3demes/simulation1/files` or `7demes/simulation1/files`.

### XML generation
XML files are generated using an in-house python module based on the `lxml` python module and tailored BEAST XML template files. The python module and templates files are available in `python` and one can use the Jupyter notebook in `generate_xml` to generate XML files from simulated data and template files. XML files are stored in `3demes/simulation1/"BEAST_model"` or `7demes/simulation1/"BEAST_model"`. `"BEAST_model"` subfolders corresponds to the BEAST models tested in the original analysis:
- `dta` (CTMC model)
- `mascot` 
- `basta`
- `glm` (MASCOT-GLM model)

An HKY prior was used to model sequence evolution. Population dynamics are described by a constant size coalescent model in CTMC and we assumed constant and equal deme sizes in the approximations of the structured coalescent model (BASTA and MASCOT). Concerning the migration rates, asymmetric migration matrices were considered in CTMC and we implemented Bayesian stochastic search variable selection (BSSVS) on the migration rates for all algorithms to avoid over-parametrization. Native priors and operators are used for each algorithm. 

### Analysis of BEAST outputs and simulated epidemics
Scripts to analyze the log and tree files generated by BEAST are located in `analyses`. They rely on multiple in-house python modules and the following Python libraries: `pandas`, `numpy`, `dendropy`, `Biopython`, `sklearn`, `scipy`, `io`, `concurrent`, and `itertools`.

BEAST log files generated by the same algorithm are summarized into tab-delimited files. These files store basic summary statistics of inferred parameters (mean, median, 95% CrI, 95% HPD, standard deviation, minimum, maximum) as well as the effective sample size (ESS) computed using BEAST2 function and converted in cython and the Bayes factor associated to each migration rate. For BASTA chains, backwards-in-time migration counts were converted into forwards-in-time migration counts. Root location probabilities are either retrieved from log files (CTMC and BASTA) or from the maximum clade credibility (MCC) trees (MASCOT).  

We also extracted the lineage migration counts and introduction dates from tree posterior distributions. 

In parallel, total migration counts, lineage migration counts, introduction dates, and root location of the most recent common ancestor (MRCA) of the sampled tips are calculated from simulated transmission chains and collated in tab-delimited files. Nexus trees for each sample are also extracted from simulated transmission chains and stored in the `files` sub-directory.

### Comparison between simulated and estimated values 
BEAST chains that displayed at least one continuous parameter with an ESS value below 200 were discarded. 

We compared the simulated and inferred phylogenies using linear regression.

We compared the inferred and simulated total migration counts, lineage migration counts, and introduction dates using:
- Kendall's tau correlation coefficient
- Probability that the true value of the parameter is in (referred to as calibration in the manuscript)
- Average relative 95% HPD width
- Average relative bias
- Weigthed interval score

All figures presented in `figures` can be reproduced using the following R scripts: `plot_phylogeny.R`, `plot_lineage_migrations.R`, `plot_genetic_params_and_total_migrations.R` and `plot_introduction_dates.R`. 

## Empirical data analyses
### RABV
To illustrate the differences between the discrete phylogeographic approaches, we reanalyzed a dataset of 233 RABV glycoprotein gene sequences isolated in the Philippines between 2004 and 2010 ([Saito et al., 2013](https://doi.org/10.1371/journal.pntd.0002144) and [Tohma et al., 2014](https://doi.org/10.1016/j.meegid.2014.01.026)).

XML, log files, and MCC trees are available in `rabies`. Log files and tree posterior distributions were analyzed using the Jupyter notebook and then the results were ploted using the R script in `rabies`.

### SARS-CoV-2
We also reanalyzed a dataset of SARS-CoV-2 genomes from its early spread in the world ([Lemey et al., 2020](https://doi.org/10.1038/s41467-020-18877-9)).

XML, log files, and MCC trees are available in `sars_cov_2`. Log files and tree posterior distributions were analyzed using the Jupyter notebook and then the results were ploted using the R script in `sars_cov_2`.
