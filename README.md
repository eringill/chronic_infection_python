# chronic_infection_python
### If you use this application, please cite: 
**[Gill, E.E. et al. SMDP: SARS-CoV-2 Mutation Distribution Profiler for rapid estimation of mutational histories of unusual lineages. arXiv 2024; 2407.11201v3](https://doi.org/10.48550/arXiv.2407.11201v3)**
### Overview
This application was developed by the Computational Analysis, Modelling and Evolutionary Outcomes ([CAMEO](https://covarrnet.ca/computational-analysis-modelling-and-evolutionary-outcomes-cameo/)) pillar of Canada's Coronavirus Variants Rapid Response Network ([CoVaRR-Net](https://covarrnet.ca/)). Data analysis, code and maintenance of the application are conducted by Erin E. Gill, Fiona S.L. Brinkman, and Sarah Otto. 

Given a user-provided set of SARS-CoV-2 nucleotide mutations, this application compares the probability of generating this set from the following three distributions:
- Mutations observed during the first nine months of the pandemic (pre-VoC) (global pre-VoC distribution)
- Mutations observed during the Omicron era (global Omicron distribution)
- Mutations observed in chronic infections (chronic distribution)
- Mutations observed in zoonotic spillovers from humans to white-tailed deer (deer distribution)
In addition, the application will inform the user if the mutation pattern is:
- Consistent with molnupiravir use (via examination of the transition:transversion ratio)
- A mutator lineage (contains a mutation in nsp14 that is known to increase the mutation rate of the lineage)
See Application Notes tab for more information.
   	 
### Background
SARS-CoV-2 evolution exhibits a strong clock-like signature with mutational changes accumulating over time, but this pattern is punctuated by “saltational changes”, where lineages appear with a higher number of mutations than expected from their divergence time from other lineages ([Neher (2022)](https://academic.oup.com/ve/article/8/2/veac113/6887176)). Such unusual lineages are thought to reflect long passage times within immunocompromised individuals, sharing many of the same signatures seen in chronic infections ([Harari et al. (2022)](https://www.nature.com/articles/s41591-022-01882-4)). 

When unusual lineages arise, however, it is challenging to know the evolutionary history leading to the observed genomic changes.  Other processes, including passage through animals, ([Bashor et al. 2021](https://www.pnas.org/doi/full/10.1073/pnas.2105253118), [Naderi et al. (2023)](https://elifesciences.org/articles/83685)) mutator lineages with error-prone polymerases ([Takeda et al. (2023)](https://doi.org/10.1016/j.isci.2023.106210)), and exposure to mutagens such as molnupiravir ([Gruber et al. (2024)](https://onlinelibrary.wiley.com/doi/10.1002/jmv.29642)), can also leave unusual genomic signatures. 

Given a user-provided set of nucleotide mutations or genome consensus sequence defining an unusual lineage of SARS-CoV-2, this application compares the probability of generating this set from the following four distributions:
- The list of mutations observed during the first nine months of the pandemic, prior to the spread of VoC [Harari et al. (2022)](https://www.nature.com/articles/s41591-022-01882-4).
- The list of mutations observed in Omicron-era sequences by Harari et al., included submission dates only up to 25 May 2022.
- The list of mutations compiled from 27 chronic infections of immunocompromised individuals [Harari et al. (2022)](https://www.nature.com/articles/s41591-022-01882-4).
- The list of mutations inferred from 109 separate zoonotic spillovers from humans to white-tailed deer [Feng et al. (2023)](https://www.nature.com/articles/s41467-023-39782-x). 

In the first paper, the authors demonstrate that specific lineage-defining mutation patterns occur in SARS-CoV-2 genomes that are sequenced from chronic infections vs. mutations that occurred in SARS-CoV-2 genomes sequenced around the globe at the start of the pandemic (before the rise of Variants of Concern (VOCs)). They also analyzed lineage-defining mutation patterns in VOCs, and concluded that “mutations in chronic infections are predictive of lineage-defining mutations of VOCs”.

Feng et al. sequenced hundreds of SARS-CoV-2 samples obtained from white-tailed deer in the United States. They observed Alpha, Gamma, Delta and Omicron VOCs and determined that the deer infections arose from a minimum of 109 separate transmission events from humans. In addition, the deer were then able to transmit the virus to each other. Deer infections resulted in three documented human zoonoses. The SARS-CoV-2 virus displayed specific adaptation patterns in deer, which differ from adaptations seen in humans. 

In addition, the app informs the user whether the data contain signals consistent with:
- **Past molnupiravir Use:** The transition-to-transversion ratio of mutations is calculated in the focal lineage and compared to a background ratio of ~2:1 for SARS-CoV-2 and to case-control cohort studies indicate a ratio of ~14:1 under molnupiravir treatment ([Gruber et al. (2024)](https://onlinelibrary.wiley.com/doi/10.1002/jmv.29642)). A high ratio may thus suggest past exposure to molnupiravir or a similar factor inducing transitions.
- **Mutator lineages:** Mutator alleles may contribute to the unusual features of a lineage by increasing the rate and type of mutation. Known mutators have been observed in nsp14 within the ExoN proofreading domain of SARS-CoV-2.  P203L in nsp14 was shown to have an elevated substitution rate in phylogenetic analyses, which was confirmed to double the mutation rate when passaged through hamsters ([Takeda et al. (2023)](https://doi.org/10.1016/j.isci.2023.106210)). Sites F60S and C39F in nsp14 were associated with a 22-fold and 6-fold higher substitution rate in phylogenetic analyses ([Mack et al. (2023)](https://link.springer.com/article/10.1186/s12967-020-02344-6)). We considered mutations at sites 39, 60, and 203 in nsp14 to be known mutators and mutations in sites 90, 92, 191, 268, and 273, which fall within the ExoN proofreading domain of nsp14, to be potential mutators.

**Table 1: Mutator Sites.** Known and Potential mutator sites (denoted by “Confirmed” and “Potential” in the “Site Type” column, respectively) are listed in the table below. Known sites have been confirmed experimentally, and the specific amino acid / nucleotide changes leading to mutator phenotypes are shown. Potential sites lie within the ExoN proofreading domain of nsp14 (as shown in Mack et al. 2023). The wild type amino acids, their positions within the mature nsp14 protein, encoding nucleotides and genomic locations are shown for these sites, but changes that would lead to mutator phenotypes have not been confirmed.
| **Gene** 	| **Amino Acid Change** 	| **Nucleotide Change** 	| **Site Type** 	|     **Reference**    	|
|:--------:	|:---------------------:	|:---------------------:	|:-------------:	|:--------------------:	|
| nsp14    	|          C39F         	|        G18,155T       	|   Confirmed   	|  (Mack et al. 2023)  	|
| nsp14    	|          F60S         	|        T18,218C       	|   Confirmed   	| (Takada et al. 2023) 	|
| nsp14    	|         P203L         	|        C18,647T       	|   Confirmed   	|  (Mack et al. 2023)  	|
| nsp14    	|          D90          	|  18,307-18,309 (GAT)  	|   Potential   	|  (Mack et al. 2023)  	|
| nsp14    	|          E92          	|  18,313-18,315 (GAG)  	|   Potential   	|  (Mack et al. 2023)  	|
| nsp14    	|          E191         	|  18,610-18,612 (GAG)  	|   Potential   	|  (Mack et al. 2023)  	|
| nsp14    	|          H268         	|  18,841-18,843 (CAT)  	|   Potential   	|  (Mack et al. 2023)  	|
| nsp14    	|          D273         	|  18,856-18,858 (GAT)  	|   Potential   	|  (Mack et al. 2023)  	|


### Application Use
This application accepts a list of comma separated nucleotide positions in a SARS-CoV-2 genome where lineage-defining mutations occur. **Lineage-defining mutations are the subset of mutations in a lineage that have occurred since divergence from the larger SARS-CoV-2 tree.** A list of lineage-defining mutations (the “mutation set”) for [pangolin-designated SARS-CoV-2 lineages](https://www.pango.network/) can be found [here](https://github.com/cov-lineages/pango-designation?tab=readme-ov-file). The tool will also accept a FASTA file containing a **SINGLE** SARS-CoV-2 genome consensus sequence. In this case, the [NextClade CLI](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli/index.html) is used to determine lineage-defining mutations (called private mutations in NextClade).

The application determines the likelihood of observing the mutation set as a random draw from each distribution (chronic infection, deer-specific mutations, global (pre-VOC) and global (Omicron era)). The log likelihood of observing the mutation set from each distribution is displayed (in natural log units)12.

Because the mutational data sets are sparse, the method bins sites across the genome when calculating likelihoods. The user can define the bin of interest: genes, genes splitting the spike protein into regions of interest, genome split into 500 nucleotide windows, or genome split into 1000 nucleotide windows. For a given bin choice, the log-likelihood of drawing the user-defined mutation set from each distribution is calculated from the multinomial distribution as:
```
sum(log(((distribution bin counts + 1) / sum(distribution bin counts + 1))^user bin counts))
```
The addition of one to each bin ensures that there are no bins lacking data.

### CLI

A command line interface (CLI) is available for this application. The CLI is a Python script. You can install the necessary packages with conda using the following command:

```sh
conda env create -f environment.yaml
```

Here is an example of how to run the CLI with a list of mutations and the output you can expect:

```sh
$ python covid_mutation_distribution/cli.py "C241T, C3037T, A23403G, G28881A, G28882A, G28883C"
Number of mutations: 6
Transition/Transversion ratio: 5.00

Log Likelihoods:
  chronic: -11.04
  total chronic: -11.50
  deer: -12.92
  total deer: -12.11

Best fit distribution: (np.float64(-11.040868182380382), 'global_pre-VoC')
(1.58 times more likely than the global Omicron distribution)

Mutator lineage analysis:
  No mutator lineage detected
```

Full usage information can be found by running:

```txt
usage: cli.py [-h] [--bin-size {genes_split,gene,500,1000}] [--output {text,json}] [--plot] [--plot-output PLOT_OUTPUT] [--color-palette {plasma,viridis,inferno,seaborn}] [--verbose] mutations

SARS-CoV-2 Mutation Distribution Profiler (SMDP) CLI

positional arguments:
  mutations             Comma-separated list of mutations or path to a file containing mutations

options:
  -h, --help            show this help message and exit
  --bin-size {genes_split,gene,500,1000}
                        Bin size for analysis (default: gene)
  --output {text,json}  Output format (default: text)
  --plot                Generate a plot of mutation distribution
  --plot-output PLOT_OUTPUT
                        Output file for the plot (default: mutation_distribution.png)
  --color-palette {plasma,viridis,inferno,seaborn}
                        Color palette for the plot (default: plasma)
  --verbose             Print detailed information during analysis
```

Currently, the CLI only supports a single query at a time.

## Notes on Input
- Your list can be formatted **with** or **without** nucleotide abbreviations. e.g. `C897A, G3431T, A7842G, C8293T,...`  OR `897, 3431, 7842, 8293,...`
- These coordinates MUST be **genomic** coordinates, **not gene** coordinates like `S:G107Y`
- Indels should be reported by including the first position only e.g. `ins21608` **NOT** `ins21608TCATGCCGCTGT`
- If you would like to convert gene coordinates to nucleotide coordinates, try using Theo Sanderson’s [tool](https://codon2nucleotide.theo.io/).
-- FASTA files must contain a single sequence with a canonical header (e.g. `>genome_sequence`), have one of the following suffixes: `.FASTA`, `.fasta` or `.fa` and **ALL** U nucleotides must be converted to T before upload


### Feedback
We're pleased to accept any feedback you have. You can submit an issue in the GitHub repository [here](https://github.com/eringill/chronic_infection_python).
You can also email questions, comments or suggestions to erin.gill81(at)gmail.com. You can also leave comments in the Discussions tab.

 
