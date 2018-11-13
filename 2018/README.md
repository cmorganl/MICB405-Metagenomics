# MICB405 Metagenomics Project
A metagenome annotation workflow for UBC's MICB405 Microbial Bioinformatics course 

## Outline
 
Data used for this project is from Saanich Inlet, a seasonally-anoxic fjord on the coast of Vancouver Island British Columbia that provides a model ecosystem for studying microbial commuity responses to ocean deoxygenation. These data are described further in the series of Scientific Data publications ([geochemical data](https://www.nature.com/articles/sdata2017159) and [multi-omic data](https://www.nature.com/articles/sdata2017160)). You can access these papers as well as additional background information on the study site and phenomenon of oxygen minimum zone expansion on Canvas. 

The primary objective of this workflow is to annotate metagenome-assembled genomes (MAGs) and reconstruct their metabolic potentials in relation to defined water column oxygen gradients in Saanich Inlet. 

Data you will be analyzing has been is derived from 10m, 100m, 120m, 135m, 150m and 200m depths of research cruises SI042, SI048, SI072, SI073, SI074, SI075. 

| Biosample    | SRA_accession | Paired_reads | Sample_ID         |
|--------------|---------------|--------------|-------------------|
| SAMN05224451 | NA            |              | SI042_LV_10m_DNA  |  
| SAMN05224447 | NA            |              | SI042_LV_100m_DNA |
| SAMN05224436 | NA            |              | SI042_LV_120m_DNA |
| SAMN05224437 | NA            |              | SI042_LV_135m_DNA |
| SAMN05224442 | NA            |              | SI042_LV_150m_DNA |
| SAMN05224443 | NA            |              | SI042_LV_200m_DNA |
| SAMN05224462 | SRR3724534    |              | SI048_LV_10m_DNA  |
| SAMN05224393 |               |              | SI048_LV_100m_DNA |
| SAMN05224394 |               |              | SI048_LV_120m_DNA |
| SAMN05224397 |               |              | SI048_LV_135m_DNA |
| SAMN05224398 |               |              | SI048_LV_150m_DNA |
| SAMN05224401 |               |              | SI048_LV_200m_DNA |
| SAMN05224440 | SRR3719539    | 41332070     | SI072_LV_10m_DNA  |
| SAMN05224441 | SRR3719544    | 46927189     | SI072_LV_100m_DNA |
| SAMN05224512 | SRR3719545    | 43008747     | SI072_LV_120m_DNA |
| SAMN05224513 | SRR3719563    | 42915270     | SI072_LV_135m_DNA |
| SAMN05224518 | SRR3719562    | 47580985     | SI072_LV_150m_DNA |
| SAMN05224519 | SRR3719564    | 32004368     | SI072_LV_200m_DNA |
| SAMN05224534 |               |              | SI073_LV_10m_DNA  |
| SAMN05224524 |               |              | SI073_LV_100m_DNA |
| SAMN05224525 |               |              | SI073_LV_120m_DNA |
| SAMN05224530 |               |              | SI073_LV_135m_DNA |
| SAMN05224531 |               |              | SI073_LV_150m_DNA |
| SAMN05224508 |               |              | SI073_LV_200m_DNA |
| SAMN05224529 |               |              | SI074_LV_10m_DNA  |
| SAMN05224509 |               |              | SI074_LV_100m_DNA |
| SAMN05224514 |               |              | SI074_LV_120m_DNA |
| SAMN05224515 |               |              | SI074_LV_135m_DNA |
| SAMN05224528 |               |              | SI074_LV_150m_DNA |
| SAMN05224520 |               |              | SI074_LV_200m_DNA |
| SAMN05224536 |               |              | SI075_LV_10m_DNA  |
| SAMN05224522 |               |              | SI075_LV_100m_DNA |
| SAMN05224493 |               |              | SI075_LV_120m_DNA |
| SAMN05224495 |               |              | SI075_LV_135m_DNA |
| SAMN05224497 |               |              | SI075_LV_150m_DNA |
| SAMN05224494 |               |              | SI075_LV_200m_DNA |

### What we have done

SRA-tools' `fastq-dump --split-files --gzip` could be used to download each dataset from the Sequence Read archive (SRA) using their respective SRA_accession identifiers, though we already had these data on our server.
The first processing step was [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) to automatically remove low-quality bases from ends of reads and discard reads with unacceptably low Phred quality scores. The remaining good-quality reads were assembled into metagenomic contigs using [MEGAHIT](https://github.com/voutcn/megahit). 
Binning was accomplished by combining the metagenome assemblies from the same depth, aligning each sample's trimmed reads to the concatenated metagenome to generate abundance profiles and providing these to [MetaBAT 2.11.2](https://bitbucket.org/berkeleylab/metabat) for binning the contigs into MAGs. [checkM](https://github.com/Ecogenomics/CheckM/wiki) was used to calculate the completeness and contamination of the resulting MAGs with a single-copy marker-gene analysis approach and the genome-taxonomy database toolkit (gtdbtk) was used to assign taxonomy.
The script used for the first two steps is called launch_MEGAHIT.sh and the script used for binning and taxonomic classification is called bin_multi.sh, both of which are on the GitHub page. 

Your group has been assigned a depth (one of 10, 100, 120,  135, 150, or 200 meters) to perform a time-series analysis. There are directories on the Orca sever containing data for each depth at `/projects/micb405/resources/project_2/2018/`. 
Abundance information can be found in the binned.rpkm.csv files. "Reads Per Kilobase per Million mappable reads" or *RPKM* is an abundance value normalized by the length of each sequence (contig) aligned to (we used BWA MEM with default parameters) and the number of reads that were generated for a sample during sequencing (number in a FASTQ file). 
Based on the `checkM` results we removed the low-quality MAGs (using the quality standard definition from Bowers *et al.*, 2017) leaving the medium- and high-quality MAGs in `SaanichInlet_$DEPTH/MetaBAT2_SaanichInlet_$DEPTH/MedQPlus_MAGs/` for you to analyze. 
The taxonomic lineages for the bacterial and archaeal MAGs are found in `SaanichInlet_$DEPTH/MetaBAT2_SaanichInlet_$DEPTH/gtdbtk_output/gtdbtk.bac120.classification_pplacer.tsv` and `SaanichInlet_$DEPTH/MetaBAT2_SaanichInlet_$DEPTH/gtdbtk_output/gtdbtk.ar122.classification_pplacer.tsv`, respectively. 

### What you need to do

Groups are expected to generate transcriptional abundance information to better understand the transcriptional activity of genes of interest for specific pathways. Metatranscriptomes have been copied to `/projects/micb405/resources/project_2/2018/MetaT`. Using BWA MEM as in tutorial and project 1, create SAM files for every cruise at your assigned depth. The FASTA file used for building a BWA index and the output SAM file can then be used to create a comma-separated value (csv) file with RPKM abundances for each sequence in the FASTA file.
 [Prokka](https://github.com/tseemann/prokka) should be used to predict open-reading frames (ORFs) in the MAGs you select. 
The next step is to annotate the sequences with KEGG Orthology (KO) numbers. These KO numbers correspond to gene sequences with rich, hierarchical attributions assigned to them thanks to the fine folks at the Kyoto Encyclopedia of Genes and Genomes (KEGG). This can be accomplished either by using a pairwise alignment tool (e.g. BLAST) to match homologous sequences in a KEGG database we have access to with your ORFs or use the [KEGG Automatic Annotation Server](https://www.genome.jp/kegg/kaas/) (KAAS).
Combining the RPKM abundance information and the KO numbers assigned to your sequences, visualise the metabolic map of your MAGs using [pathview](http://pathview.r-forge.r-project.org/) in R.

In summary, we have helped you identify who is in Saanich Inlet. It is your job to identify what they are doing and how are they responding to change!

Recommended usage for each software is provided at [workflow.md](https://github.com/cmorganl/MICB405-Metagenomics/blob/master/workflow.md).

The outline for the report is available at [Evaluation.md](https://github.com/cmorganl/MICB405-Metagenomics/blob/master/Evaluation.md).

## References

1. Torres-Beltrán, M., Hawley, A. K., Capelle, D., Zaikova, E., Walsh, D. A., Mueller, A., ... & Finke, J. (2017). A compendium of geochemical information from the Saanich Inlet water column. Scientific Data, 4, sdata2017159.
2. Hawley, A. K., Torres-Beltrán, M., Zaikova, E., Walsh, D. A., Mueller, A., Scofield, M., ... & Shevchuk, O. (2017). A compendium of multi-omic sequence information from the Saanich Inlet water column. Scientific Data, 4, sdata2017160.
3. Li, D., Liu, C. M., Luo, R., Sadakane, K., & Lam, T. W. (2015). MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics, 31(10), 1674-1676.
4. 
5. Parks, D. H., Imelfort, M., Skennerton, C. T., Hugenholtz, P., & Tyson, G. W. (2015). CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome research, 25(7), 1043-1055.
7. Seemann, T. (2014). Prokka: rapid prokaryotic genome annotation. Bioinformatics, 30(14), 2068-2069.

