# Impact and mitigation of reporting bias on discrete phylogeography inference 

## Study overview 
Passive surveillance systems are suspected to report infectious disease cases in a biased way by overrepresenting cases from urban and administratively important areas at the expense of rural and lowly populated areas, especially in low- and middle-income countries. Biased reporting might have an impact on modelling predictions as well as phylogeography inference. To test the potential impact of spatial bias on phylogeography inference, 
1. We simulated realistic rabies epidemics among domestic dog populations in Morocco using a stochastic metapopulation model
2. Whole-genome rabies virus sequences of the complete epidemics were simulated using a simple HKY model
3. Rabies virus sequences were sampled uniformly, in a biased way or in a way to mitigate the effects of bias    
4. The spatiotemporal history of rabies epidemics was infered by the DTA and MASCOT models implemented in BEAST based on sequence alignments
5. Inferred and simulated migration histories were compared in biased, unbiased and mitigating scenarios

## Structure of the project
### Simulations
Rabies epidemics were simulated using a discrete-time spatially-explicit stochastic model. Dog mobility was parametrized using a radiation model fitted to human population data in Morocco. Sequence evolution follows a simple HKY model without considering selection processes nor sequence heterogeneity.
All the simulation steps are gathered in the `epidemics` directory with the human population data, the inferred mobility matrix, the Moroccan shapefiles and the fasta sequences of a real RABV in `inputfiles`.
Simulated transmission chains, fasta files and traits files (location of the sampled sequences) are generated and stored in `mig1/simulation1/files`.

### XML generation
XML files are generated using an in-house python module which is based on the `lxml`python module and tailored BEAST XML template files.
XML files are stored in `mig1/simulation/"BEAST_model"`.
An HKY prior is used to model sequence evolution. Population dynamics are described by a constant size population in DTA and equal deme sizes in the structured coalescent models. Concerning the migration rates, asymmetric migration matrices are considered and BSSVS is implemented to avoid over-parametrization.  
Default operators and priors specific to each model are implemented.   

### Analysis of BEAST outputs and simulated epidemics
Scripts and intermediary files are located in `analyses`.
BEAST log files from the same BEAST model are summarized into a single tab-delimited file. This file stores basic summary statistics of inferred parameters (mean, median, 95%-CI, 95%-HPD, standard deviation, minimum, maximum) as well as the ESS computed using BEAST2 function. For MASCOT files, backwards-in-time migration rates are converted into forwards-in-time migration rates. Since BSSVS is implemented on migration rates, the Bayes factor is calculated using the corresponding indicator variables. Finally, root location probabilities is retrieved either from log files (DTA model) or from mcc trees (MASCOT model).  
In parallel, migration rates and root location of the MRCA of the sampled tips are calculated from simulated transmission chains and collated in tab-delimited files. Newick trees for each sample are also extracted from simulated transmission chains and stored in the `files` sub-directory.
BEAST chains with the ESS of the prior, posterior or likelihood lower than 200 are discarded.
Inferred and simulated parameters are compared using majorly linear regressions.
Finally, the tree topology of the mcc trees and the inferred Newick trees are compared using linear regressions.
