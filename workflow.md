# MICB405-Metagenomics
A metagenome annotation workflow for UBC's Microbiology 405 course (Bioinformatics)

### Metagenome assembly

There are tens of genome assembly software available but just a handful that are commonly used to assemble metagenomes. The one we will be using is MEGAHIT due to its memory-efficient [succinct de Buijn graph](https://link.springer.com/chapter/10.1007/978-3-642-33122-0_18) (data structure used to store and tranverse k-mers), rapid assembly algorithms, and great results. As we learned during the first project, the server has its RAM and CPU limitations so we will be using just 2 threads, and limiting the amount of RAM to 1/14th of the total memory:

```bash
megahit -1 $FORWARD_READS -2 $REVERSE_READS \ 
--k-min 27 --k-max 147 --k-step 20 \
--min-contig-len 1000 -m 0.07 -t 2 \
--out-dir /home/micb405/GroupX/project_2/MEGAHIT/SI072_LV_$DEPTH
```

An astute observer would also see the minimum contig length threshold is set to 1000bp. This is because metagenomic binning (and other software) perform best with long contigs so I chose the arbitrary cut-off of 1000bp. There are plenty of other options worth playing with and optimizing (hint hint `--min-contig-len`) so if you end up testing more assemblies feel free to discuss this in your report! I just ask you to wait until all other groups have finished assembly and binning before you hammer the server!

### Metagenomic binning

Like assembly, there are plenty of "binners" to choose from but one that has been shown to perform consistently better than its competition is MaxBin 2.0. Unlile many other software we used in the course, this was not pre-installed as part of ORCA's docker container. I've installed this in my *bin* directory (`/home/cmlang/bin/MaxBin-2.2.4/run_MaxBin.pl`) for everyone to use:

```bash
perl5.26.0 /home/cmlang/bin/MaxBin-2.2.4/run_MaxBin.pl
```

Using the default parameters will work, but again there are some that can be tuned to change the quality and number of bins generated. For this step, please __limit the number of threads to 2__. For information on how to use the software and its outputs, refer to the [README](https://downloads.jbei.org/data/microbial_communities/MaxBin/README.txt) and the help statement.

### checkM

checkM is software for determining the completeness and contamination, among other things, of single-cell genome assemblies, MAGs, and well genomes in general. Read the [paper](http://genome.cshlp.org/content/25/7/1043.full.pdf+html) to learn how it works... and what the hazards may be.

Unfortunately, we were unable to install checkM in time but this will be updated when it becomes available. Send me an email with the path to your MAGs when you have finished binning and __I will run checkM for your group__. The command I will be using is:

```bash
checkm lineage_wf --tab_table -x .fasta --threads 4 --pplacer_threads 4 $BIN_DIR \
/mnt/nfs/sharknado/Connor_MICB405_sandbox/ProcessedData/checkM/Reference/$sid\_checkm_output/ >/mnt/nfs/sharknado/Connor_MICB405_sandbox/ProcessedData/checkM/Reference/$sid\_checkM_stdout.tsv
```

By the end of this stage, your group will have MAGs for your metagenome and can easily separate the wheat from the chaff. You will probably be swimming in chaff.

### Taxonomic assignment

Great, so now you have genomes but you have no idea what organisms they represent! We will *try* to change that here. There are a couple methods for doing this but one that is currently making waves is called [MASH](http://mash.readthedocs.io/en/latest/).

### Genome annotation

Prokka.

__To focus your enthusiasm we will be investigating the nitrogen-cycle in Saanich Inlet.__

### Timeline

__Gantt__
