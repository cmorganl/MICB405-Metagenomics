# Taxonomic profiling using the 16S rRNA gene

This workflow will build a phylogenetic tree using the 16S genes annotated in your MAGs by Prokka.

This assumes you have binned your metagenome and analyzed the resulting MAGs with Prokka.

Estimated time to completion: 30 minutes

## Extracting 16S sequences from Prokka outputs (5 minutes)

In addition to annotating functional genes, Prokka also finds RNA genes (both ribosomal and transfer) in the input
genomic sequence. This first step pulls the 16S ribosomal RNA sequences from the Prokka .ffn files to create a single
 FASTA file with which we can build a multiple sequence alignment from.

The names of the open-reading frames (ORFs) will be in the .tsv files in the Prokka output directories.
 We need these annotations find the 16S rRNA genes in the .ffn files. You should not expect the Prokka .tsv file to always
 contain a 16S gene, especially when working with MAGs since they are often times pathetically incomplete.

```bash
grep "16S ribosomal RNA" Prokka/Bins/SI072_LV_*/*tsv >bin_16S_genes.txt
```

This file contains the tsv file name (thanks `grep`!), ORF name, and whether the 16S gene is predicted to be partial or
near full-length. The next step is to find all header names for all 16S rRNA within the .ffn files:

```bash
while read line
do ffn=$( echo $line | awk -F':' '{ print $1 }' | sed 's/.tsv/.ffn/g' )
prefix=$( echo $line | awk '{ print $1 }' | awk -F':' '{ print $2 }' )
grep "$prefix" $ffn; done<bin_16S_genes.txt >bin_16S_headers.txt
```

The next step uses a python script which, when provided the name of a FASTA file and a list of headers, will create a new
FASTA file which only contains those sequences listed (i.e. its a subset). To make our lives easier by only running the
command once, we can create a single .ffn file which will serve as an input for `FastaSubsetter.py`:

```bash
cat Prokka/Bins/SI072_LV_*/*ffn >tmp_All_bins.ffn
```

Then, launch `FastaSubsetter.py` accoring to the usage (invoked by `/home/micb405/resources/project_2/FastaSubsetter.py -h`)

```bash
/home/micb405/resources/project_2/FastaSubsetter.py -l bin_16S_headers.txt \
-i tmp_All_bins.ffn -o bin_16S_genes.ffn -m 1 -v
```

There! bin_16S_genes.ffn contains all 16S ribosomal RNA genes detected by Prokka.
We are now ready for multiple sequence alignment.

## Multiple alignment (10 minutes)

While there are an overwhelming number of options available to generate a multiple alignment, and perhaps some that are better,
we will use MAFFT.

```
mafft --maxiterate 1000 --localpair bin_16S_genes.ffn >bin_16S_genes.mfa
```

Now, this is not very useful since we don't immediately know which bin each sequence corresponds to each Prokka ORF name.
To replace the Prokka bin name with the MaxBin bin name:

```
while read line
do bin=$( echo $line | awk -F':' '{ print $1 }' | awk -F/ '{ print $NF }' | sed 's/.tsv//g' )
prefix=$( echo $line | awk '{ print $1 }' | awk -F':' '{ print $2 }' | awk -F_ '{ print $1 }' )
sed -i "s/$prefix/$bin/g" bin_16S_genes.mfa
done<bin_16S_genes.txt
```

## [OPTIONAL] Next level community analyses (4 hours)

[EMIRGE](https://github.com/csmiller/EMIRGE) is a tool for assembling 16S sequences from your metagenomic reads.
The first step is to perform the EMIRGE analysis on your raw reads. I used the shell script `run_EMIRGE.sh` in "Code/"
to generate the file `EMIRGE/SI072_LV_100m_DNA/SI072_LV_100m_DNA.emirge.fasta` which contains all 16S sequences for the
100m metagenome.

```bash
grep "^>" EMIRGE/SI072_LV_100m_DNA/SI072_LV_100m_DNA.emirge.fasta >EMIRGE/SI072_LV_100m_DNA/headers.txt
while read line
do sid=$( echo $line | awk -F'|' '{ print $2 }' | gawk -F. '{ print $1 }' )
tax=$( grep $sid SILVA_128_SSURef_taxa_headers.txt | sed 's/.* Bacteria;/Bacteria;/g' | sed 's/.* Archaea;/Archaea;/g' )
echo $line,$tax
done<EMIRGE/SI072_LV_100m_DNA/headers.txt >EMIRGE/SI072_LV_100m_DNA/header_map.csv
```

The next step will involve using sed to clean up the headers so there are no redundancies and the taxonomy will remain
after multiple alignment, which comes down to replacing spaces with underscores.

```bash
 while read line
  do
  sid=$( echo $line | awk -F'|' '{ print $2 }' | gawk -F. '{ print $1 }' )
  prior=$(echo $line | awk -F= '{ print $2 }' | gawk '{ print $1 }')
  tax=$( grep $sid SILVA_128_SSURef_taxa_headers.txt | sed 's/.* Bacteria;/Bacteria;/g' | sed 's/.* Archaea;/Archaea;/g' | sed 's/.* Eukaryota;/Eukaryota;/g')
  if [ -z "$tax" ]; then tax=Unknown; fi
  sed -i "s/$line/>$sid\_Prior=$prior\_$tax/g" EMIRGE/SI072_LV_165m_DNA/SI072_LV_165m_DNA.emirge_taxa.fasta
  done<EMIRGE/SI072_LV_165m_DNA/headers.txt
```

## Phylogenetic Tree construction (5 minutes)

Once you have your input multiple alignment (for example, bin_16S_genes.mfa), generating a phylogenetic tree is simple.
We can use FastTree since it is installed on the server:

```
FastTree -nt -gtr <bin_16S_genes.mfa >bin_16S_genes.tree
```

And that's it!
