# Workflow generating summary figures from Mash and checkM outputs

## checkM filtering

This step is very simple, after you have decided on your Completeness/Contamination/Strain Heterogeneity thresholds! Honestly, choosing these thresholds is the most difficult part and could suck hours of your precious life. Once you have them, its only a single command:

```bash
awk -F"\t" '{ if ($12>10 && $13<5) print $0 }' *checkM_stdout.tsv >GT10Complete_LT5Contam_MAGs_checkM.tsv
```

Let's dissect this: [awk](https://linux.die.net/man/1/awk) is a streaming pattern processing tool. I mostly use it for parsing data with a separator, which is specified by `-F` and tabs are "\t". So, this command parses lines (or rows) where the fields (or columns) are separated by tabs. Next we have the if statement. If the 13th column (completeness) is greater than 10 and the 12th column (contamination) is less than 5 the entire line is printed - $0 means the whole line in (g)awk. The checkM data for all of my "purified" MAGs are now in GT10Complete_LT5Contam_MAGs_checkM.tsv.

## Mash filtering

Like alignment tools, Mash will by default write multiple potential matches that meet whatever thresholds are specified. The most important thresholds affecting your output table (beside `-t` which I don't use) are `-v` and `-d`, determining the maximum p-value and maximum distance to be reported, respectively. Play around with them! I eventually settled on `-v 1E-8` which roughly translated into only matches exceeding 3 21-mers out of 1000 were reported. I ran my Mash command individually on both `/home/micb405/resources/project_2/refseq.genomes.k21s1000.msh` and `/home/micb405/resources/project_2/Saanich_QCd_SAGs_k21s1000.sig.msh` sketch databases. Here is a bash while loop that will run `mash dist` with each of your MAGs in GT10Complete_LT5Contam_MAGs_checkM.tsv against the RefSeq genomes sketch (just the good stuff!).

```bash
while read line
do 
bin=$( echo $line | awk '{ print $1 }')
sid=$( echo $bin | awk -F. '{ print $1 }')
if [ -f MaxBin/$sid/$bin.fasta ]
    then
    mash dist -v 1E-8 /home/micb405/resources/project_2/refseq.genomes.k21s1000.msh MaxBin/$sid/$bin.fasta
fi
done<GT90Complete_LT5Contam_MAGs_checkM.tsv >RefSeq_Mash_output.tsv
```

Replace `/home/micb405/resources/project_2/refseq.genomes.k21s1000.msh` with `/home/micb405/resources/project_2/Saanich_QCd_SAGs_k21s1000.sig.msh` to find matches between your bins and the Saanich single-cell genomes.

As you can see, there are many reference matches for each query (MAG). To select the best annotations you can again use awk. This command is adapted from a command on a great GitHub page that I frequent, [Stephen Turner's bash one-liners](https://github.com/stephenturner/oneliners). 

```bash
cat RefSeq_Mash_output.tsv Saanich_Mash_output.tsv | sort -t$'\t' -k2,2 | awk '{ if(!x[$2]++) {print $0; dist=($3-1)} else { if($3<dist) print $0} }' >Mash_classifications.BEST.tsv 
```

At this point there will have to be a lot of text editing to replace the genome identifiers with lineages. I'm sure there is a better way, I just haven't thought of it. You may have noticed that the first column in the Mash output (reference) is a fasta file (e.g. GCF_000010405.1_ASM1040v1_genomic.fna.gz). This is neither the name of an organism nor lends any lineage-specific information. At this point, we do not have a map between GenBank assembly ID (such as GCF_000010405.1) and taxonomic lineage. SO! Head over to NCBI's website and search for each of your GenBank IDs then replace them with their taxonomy (e.g. Bacteria; Proteobacteria; Gammaproteobacteria; unclassified Gammaproteobacteria; sulfur-oxidizing symbionts; Calyptogena okutanii thioautotrophic gill symbiont) in a text editor or (SIGH) excel...

## Silva LAST-ing

Its almost guaranteed in the case of non-human associated environmental genomics, databases will not contain a representative genome of all organisms in your sample. Due to this, we must turn to using marker genes and the *de facto* standard of microbial ecology is the 16s ribosomal RNA gene. 

SILVA is one of many 16s databases available and we will use it to try to assign taxonomy to the remaining bins. 

```bash
while read line; do bin=$( echo $line | awk '{ print $1 }'); sid=$( echo $bin | awk -F. '{ print $1 }'); if [ -f MaxBin/$sid/$bin.fasta ]; then best_hit=$(lastal -f TAB -P 4 /home/micb405/resources/project_2/db_SILVA_128_SSURef_tax_silva MaxBin/$sid/$bin.fasta | grep -v "^#" | head -1); echo $bin,$sid,$best_hit | sed 's/,\| /\t/g'; fi; done<GT10Complete_LT5Contam_MAGs_checkM.tsv >LAST_SILVA_alignments.BEST.tsv
```

After looking in the output file you will see that there is NO TAXONOMY! What has happened is there are spaces in the header of each FASTA record in the database and LAST splits the header on the spaces and only uses the first as the reference name. This is unfortunate because these are totally useless to us humans, that is, without a mapping file! 

```bash
while read line; do accession=$( echo $line | awk '{ print $4 }'); bin=$( echo $line | awk '{ print $1 }' ); if [ ! -z $accession ]; then last_hit=$( grep "$accession" /home/micb405/resources/project_2/SILVA_128_SSURef_taxa_headers.txt | awk '{ $1=""; print $0 }'); echo $bin,$last_hit; fi; done<LAST_SILVA_alignments.BEST.tsv >LAST_SILVA_classifications.BEST.csv
```

## Finale

At the end of these stages you should have a Mash output file like this:

----- | ------------------------- | ---------- | -------- |
Bacteria.Proteobacteria.Gammaproteobacteria.Oceanospirillales.OM182_clade.uncultured_gamma_proteobacterium                                                 | SI072_LV_100m.041 | 0.0121679 | 0            | 632/1000 |
Bacteria.Proteobacteria.Gammaproteobacteria.unclassified Gammaproteobacteria.sulfur-oxidizing symbionts.Calyptogena okutanii thioautotrophic gill symbiont | SI072_LV_10m.005  | 0.243761  | 5.4271e-13   | 3/1000   |
Bacteria.FCB group.Bacteroidetes/Chlorobi group.Bacteroidetes.Flavobacteriia.unclassified Flavobacteriia                                                   | SI072_LV_10m.008  | 0.0166683 | 0            | 544/1000 |
Bacteria.Terrabacteria group.Actinobacteria.Actinobacteria.Micrococcales.Microbacteriaceae.Leifsonia.Leifsonia xyli.Leifsonia xyli subsp. xyli             | SI072_LV_10m.034  | 0.243761  | 1.43878e-12  | 3/1000   |
Bacteria.Terrabacteria group.Actinobacteria.Acidimicrobiia.Acidimicrobiales.Acidimicrobiaceae.Ilumatobacter.Ilumatobacter coccineus                        | SI072_LV_10m.041  | 0.243761  | 1.45335e-11  | 3/1000   |
Bacteria.Proteobacteria.Gammaproteobacteria.Legionellales.Legionellaceae.Legionella                                                                        | SI072_LV_120m.013 | 0.243761  | 4.70704e-12  | 3/1000   |
Bacteria.Proteobacteria.Gammaproteobacteria.Legionellales.Legionellaceae.Legionella                                                                        | SI072_LV_135m.005 | 0.243761  | 4.41006e-12  | 3/1000   |
Archaea.TACK group.Thaumarchaeota.Nitrosopumilales.Nitrosopumilaceae.Candidatus Nitrosoarchaeum.Candidatus Nitrosoarchaeum koreensis                       | SI072_LV_135m.004 | 0.23011   | 2.11252e-17  | 4/1000   |
Archaea.TACK group.Thaumarchaeota.Nitrosopumilales.Nitrosopumilaceae.Candidatus Nitrosoarchaeum.Candidatus Nitrosoarchaeum koreensis                       | SI072_LV_150m.007 | 0.243761  | 9.76439e-13  | 3/1000   |

Your LAST classifications table should look like this:

----- | -------------------------
SI072_LV_100m.026 | Bacteria;Actinobacteria;Acidimicrobiia;Acidimicrobiales;Sva0996 marine group;uncultured actinobacterium                                                                 |
SI072_LV_100m.041 | Bacteria;Proteobacteria;Gammaproteobacteria;Chromatiales;Chromatiaceae;Halochromatium;uncultured gamma proteobacterium                                                  |
SI072_LV_100m.059 | Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus durans                                                                            |
SI072_LV_10m.005  | Bacteria;Bacteroidetes;Flavobacteriia;Flavobacteriales;Cryomorphaceae;uncultured;uncultured bacterium                                                                   |
