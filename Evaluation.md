# MICB405-Metagenomics: Evaluation
A metagenome annotation workflow for UBC's MICB405 Microbial Bioinformatics course

## Outline
The project will be assessed by the following criteria:
1. Completion of the proposed workflow. Expected outputs are available and formatted correctly. (20%)
2. A written report (3000 - 3500 words) adhering to the structure below (60-80%)
3. Oral presentation following the structure of the written report with emphasis on visual analysis of results (20%)

## Report Structure
__Introduction__ (750 words)
 - Introduce Saanich Inlet as a model ecosystem for studying microbial community responses to ocean deoxygenation e.g. seasonal cycles, relevant biogeochemistry, pervious studies, etc. 
 - Provide relevant information regarding metagenomic approaches (data generation and application), pros, cons, and alternatives

__Methods__ (750 words)
 - Reference the bash script(s) containing all commands. These are to be either hosted on a group member's GitHub page and/or emailed to Connor and I will upload them to the MICB405-Metagenomics GitHub page under `student_scripts/Group*/`.
 - Explain why each step was used, what the outputs are, etc. For example, metagenome assembly was used to build larger sequences (contigs) from the sequencing reads to reduce data size and redundancy, as well as inputs for binning. __Don't be silly and plagiarise this__.
 - Explain any potential deviations from the default parameters or the suggested workflow.

__Results__ (1000 words)
 - Provide a summary of the assembly, binning, checkM, and Mash (or other taxonomic assigner) outputs, as well as those from other software, if applicable.
 - Report on the Prokka annotations such as rRNA predicted, # genes annotated, nitrogen-cycling gene annotations, and potential metabolisms of the bins greater than 10% complete and less than 5% contaminated.
 - Report the abundance of genes belonging to specific taxa used in nitrogen cycling pathways

__Discussion__ (1000 words)
 - Describe the pros/cons of using a Mash (or more generally, genome-based taxonomic assigners) compared to taxonomic assignment software using marker-genes.
 - State any challenges encountered (we know there will be some!) and troubleshooting (parameter or software changes) - briefly!
 - Relate the results back to the introduction and what is or is not sensical. For example, discuss how the biology relates to the geochemical state of the environment). How do the different population genome bins in your data set contribute to the nitrgoen cycle? Do they encode partial or complete pathways? How modular is the nitrogen cycle based on your data?
 - Additional questions worth pursuing further.

__Figures and tables__ (≥4 with captions)
Some recommended figures to include:
 - Geochemical gradients (temperature, salinity, nutrients (phosphate, silicate), nitrogen compounds (nitrate, nitrite, ammonia, nitrous oxide), oxygen and sulfide.  
 - scatter plot comparing contamination % (Y-axis) versus completion % (X-axis) betweenpopulation genome bins. Consider color coding the bins by taxonomic rank e.g. phylum, order, etc. Indicate high, medium and low contamination bins (see Connect for an example).
 - RPKM bubble-plot of each N-cycling gene versus taxonomy
 - Barplot with the number of taxa assigned by each method, using *Phylum* as a categorical variable
 - Small subunit ribsomal RNA (SSU or 16S rRNA) gene phylogenetic tree from those sequences identified (either by aligning to SILVA or using [EMIRGE](https://github.com/csmiller/EMIRGE))
 - Population structure diagram

__References__ (>20)
Necessary citations include:
1. Each software used
2. Data papers
3. Bonus a) Meta-analysis of metagenome assembly and/or binning
4. Bonus b) Pre-prints from any server you like!

### This evaluation rubric should serve as a useful scaffold for your Group. How you expand on it will depend on your combined interests and analytic choices. Let the science take you somewhere!
