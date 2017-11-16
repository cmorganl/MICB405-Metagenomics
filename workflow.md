# MICB405-Metagenomics
A metagenome annotation workflow for UBC's MICB405 Microbial Bioinformatics course

### Data

FASTQ files are located on the orca server in `/home/micb405/data/project_2/`. These are the subsetted FASTQ files (32 million reads per file) mentioned in the README. These data should __not be copied__ to your individual or group directories. 

### Metagenome assembly

There are numerous genome assembly software tools available but just a handful that are commonly used to assemble metagenomes. The one we will be using is MEGAHIT due to its memory-efficiency [succinct de Buijn graph](https://link.springer.com/chapter/10.1007/978-3-642-33122-0_18) (data structure used to store and tranverse k-mers), rapid assembly algorithms, and high-quality results. As we learned during the first project, the server has its RAM and CPU limitations so we will be using just 2 threads, and limiting the amount of RAM to 1/14th of the total memory:

```bash
megahit -1 $FORWARD_READS -2 $REVERSE_READS \ 
--k-min 27 --k-max 147 --k-step 20 \
--min-contig-len 1000 -m 0.07 -t 2 \
--out-dir /home/micb405/GroupX/project_2/MEGAHIT/SI072_LV_$DEPTH
```

An astute observer would also see the minimum contig length threshold is set to 1000bp. This is because metagenomic binning (and other software) perform best with long contigs so I chose the arbitrary cut-off of 1000bp. There are plenty of other options worth playing with and optimizing (hint hint `--min-count`) so if you end up testing more assemblies feel free to discuss this in your report! I just ask you to wait until all other groups have finished assembly and binning before you pin the server with many variations on a theme!

### Metagenomic binning

Like assembly, there are a multiplicity of "binners" to choose from. However, one that has been shown to perform consistently better than its competition is MaxBin 2.0. Unlike most other software tools used in the course, MaxBin 2.0 was not pre-installed as part of ORCA's docker container. I've installed this in `/home/micb405/resources/project_2/MaxBin-2.2.4/` for everyone to use by:

```bash
perl5.26.0 /home/micb405/resources/project_2/MaxBin-2.2.4/run_MaxBin.pl
```

__Note__: FragGeneScan, an open-reading frame (ORF) prediction tool, is not in your shell's path. You will need to include `/home/micb405/resources/project_2/FragGeneScan1.30` in your `$PATH`.

Using the default parameters will work, but again certain parameters can be tuned to change the quality and number of bins generated. For this step, please __limit the number of threads to 2__. For information on how to use the software and its outputs, refer to the [README](https://downloads.jbei.org/data/microbial_communities/MaxBin/README.txt) and the help statement.

__Note:__ 'bin', 'MAG', and 'population genome' are all synonymous and are therefore used interchangably.

### checkM

checkM is a popular software tool for determining the completeness and contamination of single-cell genome assemblies, MAGs, and reference genomes. Read the [paper](http://genome.cshlp.org/content/25/7/1043.full.pdf+html) to learn how it works... and what the hazards may be.

Unfortunately, we were unable to install checkM in time but this will be updated when it becomes available. Send me an email with the path to your MAGs when you have finished binning and __I will run checkM for your group__. The command I will be using is:

```bash
checkm lineage_wf --tab_table -x .fasta --threads 4 --pplacer_threads 4 $BIN_DIR \
/mnt/nfs/sharknado/Connor_MICB405_sandbox/ProcessedData/checkM/Reference/$sid\_checkm_output/ >/mnt/nfs/sharknado/Connor_MICB405_sandbox/ProcessedData/checkM/Reference/$sid\_checkM_stdout.tsv
```

By the end of this stage, your group will have MAGs also known as population genome bins for your metagenome and can easily separate the wheat from the chaff. You will probably be swimming in chaff so think hard about what contsitutes a good assemble [paper](http://www.nature.com/articles/nbt.3893).

### Taxonomic assignment

Great, so now you have genomes but you have no idea what organisms they represent! We will *try* to change that condition here with the caveat that the extent of completion and the availability of known reference genomes are needed to make the most accurate assignments. There are a few methods for doing this but one that is currently making positive waves is an alignment-free method called [Mash](http://mash.readthedocs.io/en/latest/). Briefly, it leverages a reduced data representation called a MinHash sketch that represents a genome with an array of hashed k-mers of length S. This is then used to compare against the sketches of other genomes using the [Jaccard Index](https://en.wikipedia.org/wiki/Jaccard_index). Genomes which are very similar have many shared k-mer hashes within their sketch and will therefore have a Jaccard Index approaching 1. The original publication has derived a novel distance using MinHash sketches which is specifically for estimating the mutation rate between genomes. By determining the Mash distance between two genomes, their evolutionary relatedness can theoretically be inferred. For more information, refer to the [Mash publication](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x).

Mash can be called using `mash`. There are several different functions within the master command. Your goal here is to compare your good bins to all genomes in RefSeq - a feat that would potentially require hours using alignment methods such as BLAST! The sketch of each genome is located in the file `/home/micb405/resources/project_2/refseq.genomes.k21s1000.msh`.

```
mash dist
```

Now, although RefSeq (a non-redundant database of all genomes in GenBank) is truly awesome with tens-of-thousands of genomes, we are still unable to assign relatives to many bins. Can you think of why? There are 91,282 sketches in the database we are using and that is *after* filtering out the redundant genomes! Yes, there was just shy of 10,000 *E. coli* genomes in there before... 

Therefore, we must turn to alignment-based methods for querying marker genes (16s rRNA) from our MAGs against a different reference database, SILVA. A LAST-formatted database is `/micb405/resources/project_2/db_SILVA_128_SSURef_tax_silva` (note: no extension!). LAST is similar to BLAST but orders of magnitude faster.

By the end of this stage, the goal is to have a tabular file mapping a bin to a taxonomy, though there will invariably be some "Unknown"s. For example:

Bin ID | Taxonomy
-------------------- | -------------------------------------------
SI072_LV_185m.Bin001 | Bacteria;Proteobacteria;Gammaproteobacteria
SI072_LV_185m.Bin002 | Unknown
. | .
. | .

### Genome annotation

[Prokka](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btu153) is a popular microbial genome annotation tool. It is able to rapidly annotate genomes (1-10 minutes) using custom databases and generates a variety of useful outputs. From these annotations, your group will reconstruct metabolism.

__To focus your enthusiasm we will be investigating the [nitrogen-cycle](http://www.genome.jp/kegg-bin/show_pathway?map=map00910&show_description=show) in Saanich Inlet.__

Genes that guide your analyses are:

Gene | Protein 
---- | ------------------------
napA | Periplasmic nitrate reductase NapA
narG | Nitrate reductase / nitrite oxidoreductase, alpha subunit 
nifD | Nitrogenase molybdenum-iron protein alpha chain
nifH | Nitrogenase iron protein NifH
nirA | Ferredoxin-nitrite reductase
nirB | Nitrite reductase (NADH) large subunit
nirK | Nitrite reductase (NO-forming)
nirS | Nitrite reductase (NO-forming) / hydroxylamine reductase 
norB | Nitric oxide reductase, subunit B
norC | Nitric oxide reductase, subunit C
nosZ | Nitrous-oxide reductase
nrfA | Nitrite reductase (Cytochrome c-552)
nxrB | Nitrite oxidoreductase, subunit B
amoA | Ammonia monooxygenase, subunit A
hzo  | Hydrazine oxidase
hao  | Hydroxylamine oxidoreductase

As a fail-safe, we have included the Kyoto Encyclopedia of Genes and Genomes (KEGG) database in `/home/micb405/resources/project_2` which can also be used to annotate your genomes using an alignment tool such as LAST.

### Abundance estimation

At this point, you know which bins are high quality (i.e. low contamination and high completeness) and you know which bins have which genes involved in the nitrogen cycle. From here it would be interesting to learn what the abundance of these genomes are in the environment. This informs the scale as to which each organism is potentially involved in cycling nitrogen in the Saanich Inlet water column. To estimate abundance, we map the reads back to the contigs using a short read aligner such as `bwa`. At this point it would be a good idea to ensure you understand why we can use the reads to estimate abundance of genes, contigs, or whole genomes.

Normalization is required at this point to account for variance in sequencing coverage (though, this is already accomplished by sub-sampling *in silico*) and sequence lengths. Very simply, longer sequences garner more read alignments. To do this, we use [RPKM](https://www.nature.com/articles/nmeth.1226), a method that was first introduced in analyzing RNA-Seq data. The executable for this is in `/home/micb405/resources/project_2/rpkm`.

### Figure generation and other analyses

Many of the suggested figures can be generated in R (using ggplot2!). If your group is interested in a specific visualization or analysis, for example using [Phylosift](https://phylosift.wordpress.com/), [EMIRGE](https://github.com/csmiller/EMIRGE), or [Metapathways](https://github.com/hallamlab/metapathways2), let Connor know.

### Timeline

1. MEGAHIT: 8-14 hours
2. MaxBin 2.0: 1-4 hours
3. checkM: 30-40 minutes
4. MASH: 20-45 minutes
5. Prokka: 1-2 hours
6. Figure generation: hours-days :smiley:

Let's head over to the [evaluation markdown](https://github.com/cmorganl/MICB405-Metagenomics/blob/master/Evaluation.md).
