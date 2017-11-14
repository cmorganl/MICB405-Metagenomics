# MICB405-Metagenomics
A metagenome annotation workflow for UBC's MICB405 Microbial Bioinformatics course 

## Outline
 
Data used for this project is from Saanich Inlet, a seasonally-anoxic fjord on the coast of Vancouver Island British Columbia that provides a model ecosystem for studying microbial commuity responses to ocean deoxygenation. These data are described further in the series of Scientific Data publications ([geochemical data](https://www.nature.com/articles/sdata2017159) and [multi-omic data](https://www.nature.com/articles/sdata2017160)). You can access these papers as well as additional background information on the study site and phenomenon of oxygen minimum zone expansion on Connect. 

The primary objective of this workflow is to annotate metagenome-assembled genomes (MAGs) and reconstruct the nitrogen cycle in relation to defined water column oxygen gradients in Saanich Inlet. 

To begin, metagenome sequencing reads were downloaded for samples sourced from different water column depths on August 2013 (cruise 72):

| Biosample    | SRA_accession | Paired_reads | Sample_ID         |
|--------------|---------------|--------------|-------------------|
| SAMN05224440 | SRR3719539    | 41332070     | SI072_LV_10m_DNA  |
| SAMN05224441 | SRR3719544    | 46927189     | SI072_LV_100m_DNA |
| SAMN05224512 | SRR3719545    | 43008747     | SI072_LV_120m_DNA |
| SAMN05224518 | SRR3719562    | 47580985     | SI072_LV_150m_DNA |
| SAMN05224513 | SRR3719563    | 42915270     | SI072_LV_135m_DNA |
| SAMN05224519 | SRR3719564    | 32004368     | SI072_LV_200m_DNA |
| SAMN05224523 | SRR3719654    | 64420137     | SI072_LV_165m_DNA |

SRA-tools' `fastq-dump --split-files --gzip` was used to download each dataset from the Sequence Read archive (SRA). Due to time and resource limitations each FASTQ file was subsetted using the `head` command to 32 million reads, leaving 64 million reads for each sample. The number 32 million was selected becasue this was the total number of reads in the sample with lowest coverage (SI072_LV_200m_DNA). Think of this subsetting as [rarefaction](https://en.wikipedia.org/wiki/Rarefaction_(ecology)). Generally speaking though we would try and work with all the data in a primary research project focused on producing a more complete representation of microbial community structure and function. 

Here is where the class takes over! The next step is to use [MEGAHIT](https://github.com/voutcn/megahit) to assemble the genomes, [MaxBin 2.0](https://sourceforge.net/projects/maxbin2/) to bin the genomes from metagenomes, [checkM](https://github.com/Ecogenomics/CheckM/wiki) to identify the best MAGs, [Mash](https://mash.readthedocs.io/en/latest/) to identify the closest genomic relative in RefSeq (i.e. assign taxonomy), and finally [Prokka](https://github.com/tseemann/prokka) to annotate the MAGs. These are just some of the tools that are available for metagenome analysis that were selected in part because they play nice on the server. 

Recommended usage for each software is provided at [workflow.md](https://github.com/cmorganl/MICB405-Metagenomics/blob/master/workflow.md).

The outline for the report is available at [Evaluation.md](https://github.com/cmorganl/MICB405-Metagenomics/blob/master/Evaluation.md)

## References

1. Torres-Beltrán, M., Hawley, A. K., Capelle, D., Zaikova, E., Walsh, D. A., Mueller, A., ... & Finke, J. (2017). A compendium of geochemical information from the Saanich Inlet water column. Scientific Data, 4, sdata2017159.
2. Hawley, A. K., Torres-Beltrán, M., Zaikova, E., Walsh, D. A., Mueller, A., Scofield, M., ... & Shevchuk, O. (2017). A compendium of multi-omic sequence information from the Saanich Inlet water column. Scientific Data, 4, sdata2017160.
3. Li, D., Liu, C. M., Luo, R., Sadakane, K., & Lam, T. W. (2015). MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics, 31(10), 1674-1676.
4. Wu, Y. W., Simmons, B. A., & Singer, S. W. (2015). MaxBin 2.0: an automated binning algorithm to recover genomes from multiple metagenomic datasets. Bioinformatics, 32(4), 605-607.
5. Parks, D. H., Imelfort, M., Skennerton, C. T., Hugenholtz, P., & Tyson, G. W. (2015). CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome research, 25(7), 1043-1055.
6. Ondov, B. D., Treangen, T. J., Melsted, P., Mallonee, A. B., Bergman, N. H., Koren, S., & Phillippy, A. M. (2016). Mash: fast genome and metagenome distance estimation using MinHash. Genome biology, 17(1), 132.
7. Seemann, T. (2014). Prokka: rapid prokaryotic genome annotation. Bioinformatics, 30(14), 2068-2069.

