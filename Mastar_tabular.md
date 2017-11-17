# Workflow generating summary figures from Mash and checkM outputs

## checkM filtering

This step is very simple, after you have decided on your Completeness/Contamination/Strain Heterogeneity thresholds! Honestly, this is the most difficult part and could suck hours of your precious life. Just a single command:

```bash
awk -F"\t" '{ if ($12>10 && $13<5) print $0 }' *checkM_stdout.tsv >GT10Complete_LT5Contam_MAGs_checkM.tsv
```

Let's dissect this: [awk](https://linux.die.net/man/1/awk) is a streaming pattern processing tool. I mostly use it for parsing data with a separator, which is specified by `-F` and tabs are "\t". So, this command parses lines (or rows) where the fields (or columns) are separated by tabs. Next we have the if statement. If the 13th column (completeness) is greater than 10 and the 12th column (contamination) is less than 5 the entire line is printed - $0 means the whole line in (g)awk. The checkM data for all of my "purified" MAGs are now in GT10Complete_LT5Contam_MAGs_checkM.tsv.

## Mash filtering

Like alignment tools, Mash will by default write multiple potential matches that meet whatever thresholds are specified. The most important thresholds affecting your output table (beside `-t` which I don't use) are `-v` and `-d`, determining the maximum p-value and maximum distance to be reported, respectively. Play around with them! I eventually settled on `-v 1E-8` which roughly translated into only matches exceeding 3 21-mers out of 1000 were reported. I ran my Mash command individually on both `/home/micb405/resources/project_2/refseq.genomes.k21s1000.msh` and `/home/micb405/resources/project_2/Saanich_QCd_SAGs_k21s1000.sig.msh` sketch databases. To select the best annotations you can again use awk. This command is adapted from a command on a great GitHub page that I frequent, [Stephen Turner's bash one-liners](https://github.com/stephenturner/oneliners). 

```bash
cat RefSeq_Mash_output.tsv Saanich_Mash_output.tsv | sort -t$'\t' -k2,2 | awk '{ if(!x[$2]++) {print $0; dist=($3-1)} else { if($3<dist) print $0} }' >Mash_classifications.BEST.tsv 
```

At this point there will have to be a lot of text editing to replace the genome identifiers with lineages. I'm sure there is a better way, I just haven't thought of it.

## Silva LAST-ing

Its almost guaranteed in the case of non-human associated environmental genomics, databases will not contain a representative genome of all organisms in your sample. Due to this, we must turn to using marker genes and the *de facto* standard of microbial ecology is the 16s ribosomal RNA gene. 

SILVA is one of many 16s databases available and we will use it to try to assign taxonomy to the remaining bins. 

```bash
while read line; do bin=$( echo $line | awk '{ print $1 }'); sid=$( echo $bin | awk -F. '{ print $1 }'); if [ -f MaxBin/$sid/$bin.fasta ]; then best_hit=$(lastal -f TAB -P 4 /home/micb405/resources/project_2/db_SILVA_128_SSURef_tax_silva MaxBin/$sid/$bin.fasta | grep -v "^#" | head -1); echo $bin,$sid,$best_hit | sed 's/,\| /\t/g'; fi; done<GT90Complete_LT5Contam_MAGs_checkM.tsv >LAST_SILVA_alignments.BEST.tsv
```

After looking in the output file you will see that there is NO TAXONOMY! What has happened is there are spaces in the header of each FASTA record in the database and LAST splits the header on the spaces and only uses the first as the reference name. This is unfortunate because these are totally useless to us humans, that is, without a mapping file! 

```bash
while read line; do accession=$( echo $line | awk '{ print $4 }'); bin=$( echo $line | awk '{ print $1 }' ); if [ ! -z $accession ]; then last_hit=$( grep "$accession" /home/micb405/resources/project_2/SILVA_128_SSURef_taxa_headers.txt | awk '{ $1=""; print $0 }'); echo $bin,$last_hit; fi; done<LAST_SILVA_alignments.BEST.tsv >LAST_SILVA_classifications.BEST.csv
```

