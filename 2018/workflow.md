# MICB405-Metagenomics
A metagenome annotation workflow for UBC's MICB405 Microbial Bioinformatics course

## Data Provided

The data you will be working with on the orca server in `/projects/micb405/resources/project_2/2018/`. These can be divided into two 'omic types:
1. Metagenomic in all of the SaanichInlet\_\*m directories
2. Metatranscriptomic data in the `Metatranscriptomes` directory. These FASTQ files should __not be copied__ to your individual or group directories.

Metagenomic data are vast and varied! Generally, these are the result of a metagenomic binning experiment to generate metagenome-assembled genomes (MAGs) using the `bin_multi.sh` shell script. 
There are numerous genome assembly software tools available but just a handful that are commonly used to assemble metagenomes.
The one we used is MEGAHIT due to its memory-efficiency [succinct de Bruijn graph](https://link.springer.com/chapter/10.1007/978-3-642-33122-0_18)
(data structure used to store and traverse k-mers), rapid assembly algorithms, and high-quality results.
To construct representative bins for each of the 6 depths that are sampled in all the 6 cruises, the metagenome assemblies were concatenated into a single FASTA file, contigs shorter than 1500bp were removed and then the assembly was 'deduplicated'. This deduplication is necessary since a large protion of the genomic sequence across the several assemblies is identical or very similar. Specifically, sequences greater than 98% similar were identified and the longest representative sequence was retained. The result of these steps is in the `*_gt1500_dedup.fasta`.
Like assembly, there are a multiplicity of "binners" to choose from. We used MetaBAT 2 to perform this multi-sample binning experiment, leveraging both the tetranucleotide frequencies of each contig (this is basically a census of all the 4-mers) and differential abundance. This differential abundance is highly valuable information that can be used for teasing apart similar genomic content that actually originates from different organisms. The assumption is sequences that are derived from the same genome (or population) will increase and decrease in abundance __together__.
__Note:__ 'bin', 'MAG', and 'population genome' are all synonymous and are therefore used interchangably.
The MAGs you should use for further analyses are in the directory `SaanichInlet_$DEPTH/MetaBAT2_SaanichInlet_$DEPTH/MedQPlus_MAGs/`.

After binning the quality-controlled metagenomic reads from each sample of the corresponding depth were aligned to the binned contigs. The resulting SAM files were used for calculating RPKM of each contig for each sample. The csv file containing these data is in SaanichInlet\_$DEPTH/SaanichInlet\_$DEPTH_binned.rpkm.csv. These data can be used for determining the relative abundance of a contig and/or an entire MAG (since the contigs were renamed with the format SampleName_MAGNumber_ContigNumber you can sum by MAG to yield a total RPKM).

checkM is a popular software tool for determining the completeness and contamination of single-cell genome assemblies, MAGs, and reference genomes.
 Read the [paper](http://genome.cshlp.org/content/25/7/1043.full.pdf+html) to learn how it works... and what the hazards may be. A pseudo-command is:

```bash
checkm lineage_wf \
--tab_table \
-x .fasta \
--threads 4 \
--pplacer_threads 4 \
$BIN_DIR checkm_output/ >$sid\_checkM_stdout.tsv
```

These outputs are in the files: SaanichInlet\_$DEPTHm/MetaBAT2_SaanichInlet\_$DEPTHm/MetaBAT2_SaanichInlet\_$DEPTHm\_min1500_checkM_stdout.tsv
Information contained in these files are valuable for assessing which MAGs are worth investigating and those that are so incomplete or contaminated that pursuing these would be a waste of time. You will probably be swimming in chaff so think hard about what constitutes a good MAG and refer to this here [standards paper](http://www.nature.com/articles/nbt.3893). 
We have removed all low-quality MAGs and archived them in LowQ_MAGs.tar.gz. From here on you will only be analyzing the medium- and high-quality MAGs.

Finally, the taxonomic lineage was determined using the Genome-taxonomy Database's software toolkit, [gtdbtk](https://github.com/Ecogenomics/GTDBTk). 
This approach relies on a massive reference tree of tens-of-thousands of sequenced microbial genomes. The tree itself was inferred from concatenated gene alignments for >100 conserved single-copy marker genes. The first step is to find homologous genes in each MAG and adding the concatenated sequence (hereon referred to as query sequence) into the MSA. These query sequences are then mapped onto the reference tree. 
The position they are mapped to constrains their taxonomic lineage (say if the members of the clade it maps to are all *Lactobacillus* then the query sequence may be assigned to a taxonomy between Lacto and generally *Bacteria*) and several phylogenetic calculations are then used to confidently assign a taxonomy to each query sequence. 
The taxonomic assignments are in gtdbtk_output/gtdbtk.\*.classification_pplacer.tsv, one each for bacterial and archaeal classifications.

## Guiding analyses

To dig a little deeper into these data we offer 4 guiding questions, one of which your group is to select and attempt to answer:

1. [BIOCHEM] Identify 3-5 MAGs of interest. How do they "eat" and "breath"? 
The goal of this is to understand the lifestyle and metabolic niche of these organisms.
2. [GENOMICS] What is the genomic landscape of related MAGs?
For example, comparing the genomes of all the Epsilonproteobacteria MAGs of sufficient quality at their assigned depth.
3. ~~[STATS] How is the community transcriptionally responding to change? Mapping the transcriptional profile across time for all MAGs in the community, observing seasonal trends in differential expression for particular pathways (example finding could be transcriptional activity for the Nitrogen-cycling genes is greatest for the Proteobacteria in the summer months).~~
4. [BIOCHEM] How much do Saanich Inlet microbial communities work together?
Investigating the distributed (or not) metabolism of the carbon, nitrogen, and sulphur cycles by identifying the microbes participating at each step of relevant pathways.

The methods outlined below will help you answer the question your group is interested in. If you are *really* interested in some other aspect of the community that can be answered using data provided I encourage you to discuss it with the teaching team! We can probably find a way for it to be your group's guiding question.

## Next steps

### Genome annotation

[Prokka](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btu153)
is a popular microbial genome annotation tool.
It is able to rapidly annotate genomes (1-10 minutes) using custom databases and generates a variety of useful outputs.
To use it, I suggest providing it with the Kingdom/Domain that gtdbtk assigned to the MAG. An example command could be:

```
prokka \
--outdir Prokka_output/SaanichInlet_MAG3000/ \
--prefix MAG3000 \
--kingdom Bacteria \
--cpus 4 \
/projects/micb405/resources/project_2/2018/SaanichInlet\_10m/MetaBAT2\_SaanichInlet_10m/MedQPlus\_MAGs/MAG3000.fa
```

This should probably be wrapped inside a for-loop :smile:.

Among Prokka's various output files are predicted open-reading frames (ORFs) and common gene names where an ORF could be annotated. These are in the `.faa` and `.ffn` files. The .tsv file is also handy and provides a few more details in summarizing each ORF. Take a look a the other files to know what you have at your disposal but they likely will not be of much use for you. 

### Metabolic reconstruction

Using the amino acid ORFs in the .faa files we are able to determine the KEGG Orthology numbers to help uncover the metabolic potentials of these MAGs. You can upload your predicted ORFs to [KAAS](https://www.genome.jp/kegg/kaas/), submit jobs and download the text results files. You're probably fine sticking with the defaults here. 
KAAS uses one of three pairwise alignment methods (BLAST, GHOSTX or GHOSTZ) to perform homology searches between your ORFs (query sequences) and the KEGG database (reference sequences). Annotations are based on either the single best hit (highest scoring) or bi-directional (reciprocal) best hit. So nothing fancy here as they are using techniques you are mostly familiar with, all packaged up and hosted on a server for everyone to use!
The resulting files are tab-separated value (tsv) formatted files with two columns: the ORF name and the KO number. Unfortunately, we are not able to directly read these files into R since rows for which ORFs were not annotated only have a single column. However, it is easy to remove these rows and create a new file with only the annotated ORFs using `grep`:

```
grep '\sK' SaanichInlet_MAGx_ORFs_ko.txt >SaanichInlet_MAGx_ORFs_ko.cleaned.txt
```

These files are now ready to be loaded into R!

### Transcriptional activity

Groups are expected to generate transcriptional abundance information to better understand the transcriptional activity of genes of interest for specific pathways.
 This is necessary since presence of a gene does not always mean the gene is being transcribed and translated. So we are going to use RNA-Seq data to get one step closer. Although it would be ideal to sequence a metaproteome to actually determine what genes are expressed this is more difficult and expensive.
Metatranscriptome FASTQ files have been copied to `/projects/micb405/resources/project\_2/2018/Metatranscriptomes/`. Using `bwa mem` as in tutorial and project 1, create SAM files for every cruise at your assigned depth. The FASTA file used for building a BWA index and the output SAM file can then be used to create a comma-separated value (csv) file with RPKM abundances for each sequence in the FASTA file using the `rpkm` executable provided:

```
/projects/micb405/resources/project_2/2018/rpkm \
-c ~/ProcessedData/Prokka_output/SaanichInlet_MAG_ORFs.ffn \
-a ~/ProcessedData/MetaT_alignments/SI042_SaanichInlet_MAG_ORFs.sam
-o ~/ProcessedData/RPKM_outputs/SI042_SaanichInlet_MAG_ORFs_RPKM.csv
```

Normalization is required to account for variance in sequencing coverage and sequence lengths.
Very simply, longer sequences are aligned to more often than short sequences. To do this, we use [RPKM](https://www.nature.com/articles/nmeth.1226), a method that was first introduced in analyzing RNA-Seq data.

At this point it would be a good idea to ensure you understand why we can use the reads to estimate abundance of genes, contigs, or whole genomes.

### Figure generation

Alright, from these analyses you should be ready to load the resulting csv and tsv files into R for generating figures.

To get visualize the quality of your bins (from the perspective of completeness and contamination) you will need the checkM_stdout.tsv and the classification_pplacer.tsv files mentioned above. If you want to be all-stars and add the abundance layer, you will need the provided RPKM files that contain normalized genomic abundances for all binned contigs - this will probably comprise only but a small portion of the entire microbial community!

Pathway analyses with [pathview](http://bioconductor.org/packages/release/bioc/vignettes/pathview/inst/doc/pathview.pdf) will require the metatranscriptome RPKM files and the KAAS output(s), as well as a file mapping the Prokka sample IDs to MAGs that you are more familiar with... Prokka :pouting_cat:.
To create such a file I suggest using a for loop on either the .faa or .ffn files generated by Prokka like so:

```
for f in path/to/prokka_output/*faa
do
prokka_id=$( head -1 $f | awk -F_ '{ print $1 }' | sed 's/^>//g' )
mag_id=$( echo $f | sed 's/.faa//g')
echo $prokka_id,$mag_id
done >Prokka_MAG_map.csv
```

I have created a [short tutorial](https://htmlpreview.github.io/?https://github.com/cmorganl/MICB405-Metagenomics/blob/master/2018/pathview_trials.html) on how to view MAG nitrogen metabolism with pathview.

## Fin.

More analyses, tips and tricks may be added as we explore these samples together - stay tuned!

