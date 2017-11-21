# Estimating N-cycling gene abundance within the microbial community

This workflow is meant to generate normalized abundances for specific genes of interest. The normalization statistic
we will be using here is again Reads Per Kilobase per Million mappable reads (RPKM). Our genes of interest are those
involved in the nitrogen cycle outlined in the [workflow markdown](https://github.com/cmorganl/MICB405-Metagenomics/blob/master/workflow.md).

The result of this workflow is a csv file with the Bin identifier, Gene, Prokka gene ID, and RPKM. This can then be easily
plotted in ggplot to represent the abundance of each gene from each bin that contains a gene of interest.

This workflow assumes you have functionally and taxonomically annotated your bins using Prokka,
and have all 16 gene names (e.g., nosZ, norC, etc.) in nitrogen_cyclers.txt. This set of gene names will be referred to
as "N-cycling genes". This term does not refer to all genes which are involved in metabolising Nitrogen!

Estimated time for completion (1 - 1.5 hours)


## Creating the reference set (5 minutes)

Step one is to find the names of the Prokka-annotated N-cycling genes. This information is contained within the
tab-separated value (tsv) file and can be extracted by pattern matching (with `grep`) for the name of each N-cycling gene.

```bash
while read line; do grep $line Prokka/Bins/SI072_LV_*/*tsv >>bin_nitrogen_cycler_genes.txt; done<nitrogen_cyclers.txt
```

The lines containing the gene names were written to bin_nitrogen_cycler_genes.txt. Here are the first two lines:

```
Prokka/Bins/SI072_LV_100m.041/SI072_LV_100m.041.tsv:OCIHAJHB_01798	CDS	napA	1.7.99.4	Periplasmic nitrate reductase
Prokka/Bins/SI072_LV_120m.013/SI072_LV_120m.013.tsv:NGLKKDPM_00530	CDS	napA	1.7.99.4	Periplasmic nitrate reductase
```

The reason we need to perform the above step is that the gene names are not present in Prokka's annotated fasta file (.ffn).
 The next step is to search for each contig annotated as an N-cycling gene by its internal Prokka identifier;
 this is the alphanumeric code following the file name, 'OCIHAJHB_01798' for the first sequence above.

```bash
while read line
do ffn=$( echo $line | awk -F':' '{ print $1 }' | sed 's/.tsv/.ffn/g' )
prefix=$( echo $line | awk '{ print $1 }' | awk -F':' '{ print $2 }' )
grep "$prefix" $ffn; done<bin_nitrogen_cycler_genes.txt >bin_nitrogen_cycler_headers.txt
```

bin_nitrogen_cycler_headers.txt contains the headers of all sequences that were annotated as an N-cycling gene. Concatenating all
Prokka-annotated gene sequences together will make our collective lives easier for the next step:

```bash
cat Prokka/Bins/SI072_LV_*/*ffn >tmp_All_bins.ffn
```

I've provided a python script for subsetting a FASTA file based on a list of headers of interest called `FastaSubsetter.py`.
With the list you just generated and the FASTA file containing all nucleotide sequences for all annotated genes, you
are now able to create a FASTA file containing the nitrogen-cycling genes annotated like so:

```bash
/home/micb405/resources/project_2/FastaSubsetter.py -l bin_nitrogen_cycler_headers.txt \
-i tmp_All_bins.ffn -o bin_nitrogen_cycler_genes.ffn -m 1 -v
```

I suggest "smell-testing" the output file. Use `less` to ensure it is not empty and count the number of sequences in it.
Does this match the number of headers in "bin_nitrogen_cycler_headers.txt"?

## Alignment (0.5 - 1 hour)

Back to `bwa` for aligning the reads to the nitrogen-cycling genes. I'm sure you're all too familiar with these commands
by now! Index first:

```bash
bwa index bin_nitrogen_cycler_genes.ffn
```

And now align. I'm using the 165m depth in this example:

```bash
nohup bwa mem -t 12 bin_nitrogen_cycler_genes.ffn \
/home/micb405/data/project_2/SI072_LV_165m_DNA_R1.fastq.gz \
/home/micb405/data/project_2/SI072_LV_165m_DNA_R2.fastq.gz \
1>bin_nitrogen_cycler_genes_165m.sam 2>bin_nitrogen_cycler_genes.bwa.stderr &
```

## Abundance R'Us (15 minutes)

Just run the `rpkm` script as before:

```bash
/home/micb405/resources/project_2/rpkm -c bin_nitrogen_cycler_genes.ffn \
-a bin_nitrogen_cycler_genes_165m.sam \
-o bin_nitrogen_cycler_genes_165m_RPKM.csv \
--verbose
```

Now that you have the abundance of each gene, it would be great to map these RPKM values back to the original bin names.
Currently, the `rpkm` output contains the Prokka gene ID (e.g. AKEAEFHP_00500). This can be accomplished using a bit of `grep` and `sed`.
The file `bin_nitrogen_cycler_genes.txt` we generated earlier will be useful here!

## Abundance of all N-cycling genes in the environment(!!)

Now that you know the abundance of the genes annotated within your bins specifically in the environment,
what about those genes that were not binned?

Stay tuned...