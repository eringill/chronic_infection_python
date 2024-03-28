# chronic_infection_python
### Overview
This application was developed by the Computational Analysis, Modelling and Evolutionary Outcomes ([CAMEO](https://covarrnet.ca/computational-analysis-modelling-and-evolutionary-outcomes-cameo/)) pillar of Canada's Coronavirus Variants Rapid Response Network ([CoVaRR-Net](https://covarrnet.ca/)). Data analysis, code and maintenance of the application are conducted by Erin E. Gill, Fiona S.L. Brinkman, and Sarah Otto.
   	 
### Background
This application draws from work conducted by [Harari et al. (2022)](https://www.nature.com/articles/s41591-022-01882-4) and [Feng et al. (2023)](https://www.nature.com/articles/s41467-023-39782-x). 


In the first paper, the authors demonstrate that specific lineage-defining mutation patterns occur in SARS-CoV-2 genomes that are sequenced from chronic infections vs. mutations that occurred in SARS-CoV-2 genomes sequenced around the globe at the start of the pandemic (before the rise of Variants of Concern (VOCs)). They also analyzed lineage-defining mutation patterns in VOCs, and concluded that “mutations in chronic infections are predictive of lineage-defining mutations of VOCs”.


Feng et al. sequenced hundreds of SARS-CoV-2 samples obtained from white-tailed deer in the United States. They observed Alpha, Gamma, Delta and Omicron VOCs and determined that the deer infections arose from a minimum of 109 separate transmission events from humans. In addition, the deer were then able to transmit the virus to each other. Deer infections resulted in three documented human zoonoses. The SARS-CoV-2 virus displayed specific adaptation patterns in deer, which differ from adaptations seen in humans. 

### Application Use
This application accepts a list of comma separated nucleotide positions in a SARS-CoV-2 genome where lineage-defining mutations occur. A list of lineage-defining mutations for [pangolin-designated SARS-CoV-2 lineages](https://www.pango.network/) can be found [here](https://github.com/cov-lineages/pango-designation?tab=readme-ov-file). 
The application determines which mutation distribution best fits your list of mutations (chronic, deer, global (pre-VOC)) via likelihood calculations. The log likelihood that your list fits each distribution is displayed.
Likelihoods are calculated based on user-defined bin size (genes with split spike protein, genes, genome split into 500nt windows or genome split into 1000nt windows) as follows:
```
sum(log(((distribution bin counts + 1) / sum(distribution bin counts + 1)) ^ user bin counts))
```

#### Notes on Input
* Your list can be formatted **with** or **without** nucleotide abbreviations. e.g. `C897A, G3431T, A7842G, C8293T,...`  OR `897, 3431, 7842, 8293,...`.
* Do **NOT** include insertions or deletions (indels) e.g. `ins21608TCATGCCGCTGT, ∆23009-23011`.
* If you have an unaligned SARS-CoV-2 genome sequence and would like to use this tool, you must first place it into a phylogeny so that you can detect lineage-defining mutations. To get started, you may wish to access the tools associated with the [UCSC SARS-CoV-2 Genome Browser](https://genome.ucsc.edu/goldenPath/help/covidBrowserIntro.html#data).

### Feedback
We're pleased to accept any feedback you have. You can submit an issue in the GitHub repository [here](https://github.com/eringill/chronic_infection_python).
You can also email questions, comments or suggestions to erin.gill81(at)gmail.com. You can also leave comments in the Discussions tab.

 
