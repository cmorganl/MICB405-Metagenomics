# Workflow for finding the relative abundances of all MAGs in your dataset

In this workflow, it is assumed that you have completed the assembly with MEGAHIT and binning portion with MaxBin 2.0.
For this step, you will align the FASTQ files you used for assembly to your contigs, calculate the RPKM for each contig,
then use this output to find the average RPKM of your Metagenome-Assembled Genomes (MAGs).

Throughout this example workflow, I will be using the SI_072_150m metagenome.

Estimated time for completion (1.5-2.5 hours)

## Sequence alignment (1-2 hours)

The first step for aligning your reads back to the contigs is indexing the reference, in this case your MEGAHIT assembly.

```bash
cmlang@cmlang-ORCA:~$ bwa index /home/cmlang/ProcessedData/MEGAHIT/SI072_LV_150m/SI072_LV_150m.contigs.fa
```

Now, align the reads using bwa mem. I'm writing to my home directory but please do this within your group
 directories under the shared resource.

```bash
nohup bwa mem -t 4 /home/cmlang/ProcessedData/MEGAHIT/SI072_LV_150m/SI072_LV_150m.contigs.fa \
/home/micb405/data/project_2/SI072_LV_150m_DNA_R1.fastq.gz /home/micb405/data/project_2/SI072_LV_150m_DNA_R2.fastq.gz \
1>/home/cmlang/ProcessedData/SI072_LV_150m_DNA.sam 2>/home/cmlang/ProcessedData/SI072_LV_150m_DNA.bwa.stderr &
```

## RPKM calculation (20-30 minutes)

Great! Now you have a SAM file. Let's calculate RPKM values ASAP so we can remove that large SAM file.
You can view the options and arguments of the script like so:

```bash
/home/micb405/resources/project_2/rpkm -h
```

As you can see, it requires three things: MEGAHIT assembly, SAM file, and the name of an output comma-separated value file.

```bash
/home/micb405/resources/project_2/rpkm -c /home/cmlang/ProcessedData/MEGAHIT/SI072_LV_150m/SI072_LV_150m.contigs.fa \
-a /home/cmlang/ProcessedData/SI072_LV_150m_DNA.sam -o /home/cmlang/ProcessedData/RPKM/SI072_LV_150m_DNA_RPKM.csv
```

In addition to the basics, I like to use the `--verbose` option so I can see what is going on.

## Mean, mean MAGs (1 minute)

The last step here is to calculate the mean RPKM of each MAG. This is accomplished by summing the RPKMs of each contig
belonging to a MAG and dividing it by the number of contigs comprising the MAG. This is performed by the script
`/home/micb405/resources/project_2/find_mag_rpkm_average.py`

It requires either a single MAG FASTA file or a text file with a list of FASTAs, the csv file made by `rpkm`, and
the name of an output file.

To create a list of all your bins you could use something like: `ls MaxBin/SI072_LV_150m/*fasta >mag_list.txt`.

And finally, to get the csv file you've been waiting for:

```bash
/home/micb405/resources/project_2/find_mag_rpkm_average.py -l mag_list.txt \
-r /home/cmlang/ProcessedData/RPKM/SI072_LV_150m_DNA_RPKM.csv \
-o /home/cmlang/ProcessedData/RPKM/SI072_LV_150m_MAG_RPKM.csv
```

Now, please remove that SAM file!

### Fin.