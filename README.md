# MICB405-Metagenomics
A metagenome annotation workflow for UBC's Microbiology 405 course (Bioinformatics)

## Outline

Data used for this project is from Saanich Inlet, a seasonally-anoxic fjord. These data are described further in the series of Scientific Data publications ([geochemical data](https://www.nature.com/articles/sdata2017159) and [multi-omic data](https://www.nature.com/articles/sdata2017160)). The goal of this workflow is to annotate metagenome-assembled genomes (MAGs) and reconstruct the nitrogen cycle through the water column. 

To begin, metagenome sequencing reads were downloaded for samples from August 2013 (cruise 72):
```
SRR3719539
SRR3719544
SRR3719545
SRR3719563
SRR3719562
SRR3719654
SRR3719564
```

SRA-tools' `fastq-dump --split-files --gzip` was used to download each dataset from the Sequence Read archive (SRA). 

The next step is to use [MEGAHIT](https://github.com/voutcn/megahit) to assemble the genomes, [MaxBin 2.0](https://sourceforge.net/projects/maxbin2/) to bin the genomes from metagenomes, [checkM](https://github.com/Ecogenomics/CheckM/wiki) to identify the best MAGs, [MASH](https://mash.readthedocs.io/en/latest/) to identify the closest genomic relative in RefSeq (i.e. assign taxonomy), and finally [Prokka](https://github.com/tseemann/prokka) to annotate the MAGs.

Recommendations for each software is available at [workflow.md]().

## References

