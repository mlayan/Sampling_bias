# Impact and mitigation of sampling bias to determine viral spread: evaluating discrete phylogeography through CTMC modeling and structured coalescent model approximations

Full text and supplementary materials are available [here]().

## Abstract

Bayesian phylogeographic inference is a powerful tool in molecular epidemiological studies that enables reconstructing the origin and subsequent geographic spread of pathogens. Such inference is, however, potentially affected by geographic sampling bias. Here, we investigated the impact of sampling bias on the spatiotemporal reconstruction of viral epidemics using Bayesian discrete phylogeographic models and explored different operational strategies to mitigate this impact. We considered the continuous-time Markov chain (CTMC) model and two structured coalescent approximations (BASTA and MASCOT). For each approach, we compared the estimated and simulated spatiotemporal histories in biased and unbiased conditions based on simulated epidemics of rabies virus (RABV) in dogs in Morocco. While the reconstructed spatiotemporal histories were impacted by sampling bias under CTMC, BASTA and MASCOT reconstructions were biased even when employing unbiased samples. This resulted from assuming incorrect population dynamics but can be overcome by accommodating time-varying population sizes. Increasing the number of analyzed genomes led to more robust estimates at low sampling bias, notably for CTMC. Alternative sampling strategies that maximize the spatiotemporal coverage greatly improved the inference at intermediate sampling bias for CTMC, and to a lesser extent, for BASTA and MASCOT. We further applied these approaches to two empirical datasets: a RABV dataset from the Philippines, and a SARS-CoV-2 dataset from its early spread across the world. In conclusion, sampling biases are ubiquitous in phylogeographic analyses but may be accommodated by increasing sample size, balancing spatial and temporal composition in the samples, and informing structured coalescent models with reliable case count data.

## Study overview 
Passive surveillance systems are suspected to report infectious disease cases in a biased way by overrepresenting cases from urban and administratively important areas at the expense of rural and lowly populated areas, especially in low- and middle-income countries. Biased sampling of cases might have an impact on discrete phylogeographic predictions. To test the potential impact of spatial sampling bias on discrete phylogeographic inference, 
1. We simulated realistic RABV epidemics among domestic dog populations in Morocco using a stochastic metapopulation model
2. RABV whole-genomes associated with each case were simulated using a simple HKY model 
3. RABV genomes were sampled uniformly, in a biased way or in a way to mitigate the effects of bias    
4. The spatiotemporal history of rabies epidemics was infered by CTMC, BASTA, and MASCOT, all implemented in BEAST
5. Inferred and simulated migration histories were compared in biased, unbiased and mitigating scenarios

## Workflow
### Simulations
RABV epidemics were simulated using a discrete-time spatially-explicit stochastic model. Dog mobility was parametrized using a radiation model fitted to human population data in Morocco. Sequence evolution follows a simple HKY model without considering selection processes, varying evolutionary rate, within-host evolution, or site heterogeneity.
All the simulation steps are gathered in the `epidemics` directory. The human population data, the inferred mobility matrix, the Moroccan shapefiles and the fasta sequences of a real RABV are all located in `inputfiles`.
Two simulation frameworks were tested: `3demes` that corresponds to epidemics across 3 Moroccan demes/regions, and `7demes` that corresponds to epidemics across 7 Moroccan demes/regions. In both frameworks, approximatively 1% of infection events happen between demes/regions. 
Simulated transmission chains, fasta files and traits files (location of the sampled sequences) are generated and stored in the following subfolders: `3demes/simulation1/files` or `7demes/simulation1/files`.

### XML generation
XML files are generated using an in-house python module based on the `lxml` python module and tailored BEAST XML template files. The python module and templates files are available in the `python` directory.
XML files are stored in `3demes/simulation1/"BEAST_model"` or `7demes/simulation1/"BEAST_model"`. `"BEAST_model"` subfolders corresponds to specific BEAST models implemented with a certain set of priors:
- `dta` (CTMC model)
- `mascot` 
- `basta`
- `glm` (MASCOT-GLM model)
An HKY prior is used to model sequence evolution. Population dynamics are described by a constant size coalescent model in CTMC and contant and equal deme sizes in the approximations of the structured coalescent model (BASTA and MASCOT). Concerning the migration rates, asymmetric migration matrices are considered and BSSVS is implemented to avoid over-parametrization. Native priors and operators are used for each approach. 

### Analysis of BEAST outputs and simulated epidemics
Scripts to analyze the log and tree files generated by BEAST are located in `analyses`.
BEAST log files from the same BEAST model are summarized into a single tab-delimited file. This file stores basic summary statistics of inferred parameters (mean, median, 95% CrI, 95% HPD, standard deviation, minimum, maximum) as well as the effective sample size (ESS) computed using BEAST2 function and converted in cython. For BASTA files, backwards-in-time migration counts are converted into forwards-in-time migration counts. BSSVS is implemented on migration rates, consequently the Bayes factor is calculated using the corresponding indicator variables. Finally, root location probabilities are retrieved either from log files (CTMC and BASTA) or from the maximum clade credibility (MCC) trees (MASCOT).  
In parallel, migration rates and root location of the most recent common ancestor (MRCA) of the sampled tips are calculated from simulated transmission chains and collated in tab-delimited files. Nexus trees for each sample are also extracted from simulated transmission chains and stored in the `files` sub-directory.

All BEAST chains displaying at least one ESS value below 200 are discarded.

Inferred and simulated parameters are compared using:
- Kendall's tau correlation coefficient
- Probability that the true value of the parameter is in (referred to as calibration in the manuscript)
- Average relative 95% HPD width
- Average relative bias
- Weigthed interval score


Finally, the tree topology of the mcc trees and the inferred Newick trees are compared using linear regressions.

## Empirical data analyses
### RABV

### SARS-CoV-2



