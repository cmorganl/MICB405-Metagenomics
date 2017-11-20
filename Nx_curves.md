# Workflow to determine the most contiguous assembly

It is assumed at this point you have at least two different metagenome assemblies. Whether they are made
from trimmed reads, untrimmed reads, using different k-mers, or even different assemblers it doesn't matter - the
only constant must be they are in FASTA format. Everything moves very quickly once you have these!

Estimated time for completion (5 minutes)

## Generate Nx files (1 - 4 minutes)

The first step is to use a my script `getNx` in `/home/micb405/resources/project_2`. To automate things a bit
I'm going to create a list containing all the assembly files I want to build Nx curves for then run `getNx` on
all of these:

```bash
ls MEGAHIT/*/*fa >megahit_assemblies_list.txt
time while read line; do onx=$( echo $line | sed 's/.contigs.fa/_Nx.csv/g')
/home/micb405/resources/project_2/getNx -i $line -o $onx -v
done<megahit_assemblies_list.txt
```


Algorithm pseudocode:
```
    contig_lengths = list
    for contig in fasta:
        contig_lengths.append(length(contig))
    total_size = sum(contig_lengths)
    for i in range(0.00,1.00):
        lengths_position = -1
        contig_sum = 0
        while contig_sum < (i*total_size):
            lengths_position += 1
            contig_sum += contig_lengths[lengths_position]
        write(i,contig_lengths[lengths_position])
```

Output is a csv file with two fields:

| Genome_proportion | Contig_size |
|-------------------|-------------|
| 0                 | 372678      |
| 0.01              | 296671      |
| 0.02              | 194128      |
| 0.03              | 144621      |
| 0.04              | 120597      |
| 0.05              | 99507       |
| 0.06              | 83520       |
| 0.07              | 69846       |
| 0.08              | 60199       |
| 0.09              | 51418       |
| . | . |
| 1                 | 1000        |

## Plot the curves (10 seconds)

The R script to plot these data is `/home/micb405/resources/project_2/Nx_curve_generator.R`
and can be used like so:

```
Rscript /home/micb405/resources/project_2/Nx_curve_generator.R MEGAHIT/*/*_Nx.csv
```

Now, don't feel restricted to using the default ggplot2 rainbow colour palette! I encourage everyone to copy this
script and edit the ggplot2 command to create a custom plot. Starting with a new colour palette or editing the background
with `theme()` are good starting points.
