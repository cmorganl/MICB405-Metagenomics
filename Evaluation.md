# MICB405-Metagenomics: Evaluation
A metagenome annotation workflow for UBC's Microbiology 405 course (Bioinformatics)

## Outline
The project will be assessed by the following criteria:
1. Completion of the the suggested workflow. Expected outputs are available and formatted correctly. (20%)
2. A written report (3000 - 3500 words) conforming to the structure below (60-80%)
3. Oral presentation following the structure of the written report with emphasis on results (20%)

## Report Structure
__Introduction__ (750 words)
 - Introduce Saanich Inlet as a model system, its unique biogeochemistry, 
 - Explain metagenomics (how data are generated and its applications), pros, cons, and alternatives

__Methods__ (750 words)
 - Reference the bash script(s) containing all commands. These are to be sent to me and I will upload them to my MICB405-Metagenomics GitHub page under `student_scripts/Group*/`.
 - Explain why each step was used, what the outputs are, etc. For example, metagenome assembly was used to build larger sequences (contigs) from the sequencing reads to reduce data size and redundancy, as well as inputs for binning. __Don't be silly and plagiarise this__
 - Explain any potential deviations from the default parameters or the suggested workflow.

__Results__ (1000 words)
 - Provide a summary of the assembly, binning, checkM, and MASH (or other taxonomic assigner) outputs, as well as those from other software, if applicable.
 - Report on the Prokka annotations such as rRNA predicted, # genes annotated, nitrogen-cycling gene annotations, and potential metabolisms of the top 10 bins.
 - Report the abundance of genes belonging to specific taxa used in major nutrient cycling pathways (e.g. nitrogen)

__Discussion__ (1000 words)
 - Describe the pros/cons of using a MASH (or more generally, genome-based taxonomic assigners) compared to taxonomic assignment software using marker-genes.
 - State any challenges encountered (we know there will be some!) and troubleshooting (parameter or software changes) - briefly!
 - Relate the results back to the introduction and what is or is not sensical (i.e. discuss how the biology relates to the geochemical state of the environment). What are the ecosystem services and functions offered by the population genomes you have discovered?
 - Questions worth pursuing further.

__Figures and tables__ (>3 with captions)
Some potential figures to be included:
 - Geochemical gradients
 - RPKM bubble-plot of each N-cycling gene versus taxonomy
 - Barplot with the number of taxa assigned by each method, using *Kingdom* as a categorical variable
 - 16s (SSU) rRNA phylogenetic tree from those sequences identified (either by aligning to SILVA or using [EMIRGE](https://github.com/csmiller/EMIRGE))
 - Population structure diagram

__References__ (>20)
Necessary citations include:
1. Each software used
2. Data papers
Bonus 1: Meta-analysis of metagenome assembly and/or binning
Bonus 2: Pre-prints from any server you like!
