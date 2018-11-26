# MICB405-Metagenomics: Evaluation
A metagenome annotation workflow for UBC's MICB405 Microbial Bioinformatics course

## Outline
The project will be assessed by the following criteria:
1. Completion of the proposed workflow. Expected outputs are available and formatted correctly. (40%)
2. A written report (4000 - 4500 words) adhering to the structure below (40-80%)
3. Oral presentation following the structure of the written report with emphasis on visual analysis of results (20%)

## Report Structure

The word counts for each section are reccommended, can be adjusted by +/- 200 words, but the topics within must be discussed!

__Abstract__ (250 word max.)
 - Consisely summarize the major results and implications of the report

__Importance__ (120 word max.)
 - Lay explanation of the significance of the research performed to the field of environmental microbiology

__Introduction__ (1250 words)
 - Introduce Saanich Inlet as a model ecosystem for studying microbial community responses to ocean deoxygenation e.g. seasonal cycles, relevant biogeochemistry, pervious studies, etc. 
 - Provide relevant information regarding metagenomic approaches (data generation and application), pros, cons, and alternatives
 - Introduce one or all of: metagenomic binning, genome assembly, taxonomic profiling
 - Justify reasoning for using cultivation-independent methods

__Methods__ (750 words)
 - Explain the workflow used by Connor to generate inputs for the project to show you understand each step.
 - Explain the tools used at each step, and why they were chosen, to answer your group's specific question.
 - Reference the bash and R(md) scripts containing all commands.
 These are to be either hosted on a group member's GitHub page and/or emailed to Connor and I will upload them to the MICB405-Metagenomics GitHub page under `student_scripts/Group*/`.
 - Explain why each step was used, what the outputs are, etc.
 For example, metagenome assembly was used to build larger sequences (contigs) from the sequencing reads to reduce data size and redundancy, as well as inputs for binning. __Don't be silly and plagiarise this__.
 - Explain any potential deviations from the default parameters or the suggested workflow.

__Results__ (1000 words)
 - Provide a summary of the checkM, and gtdbtk outputs, as well as those from other software, if applicable.
 - Report on the Prokka annotations such as rRNA predicted, # genes annotated, and potential metabolisms of the Medium- and High-quality bins.
 - Report of the relevant results from your group's selected analysis

__Discussion__ (1000 words)
 - Relate the results back to the introduction and what is or is not sensible.
 For example, discuss how the biology relates to the geochemical state of the environment).
 How do the different population genome bins in your data set contribute to the nitrogen cycle?
  Do they encode partial or complete pathways? How modular is the nitrogen cycle based on your data?
 - Discuss your reasoning behind the steps you took to answer your group's specific question.
 - Additional questions worth pursuing further.
 - State any challenges encountered (we know there will be some!) and troubleshooting (parameter or software changes) - briefly!

__Figures and tables__ (≥4 with captions)
Some recommended (not required!) figures to include:
 - Geochemical gradients for nutrients (phosphate, silicate), nitrogen compounds (nitrate, nitrite, ammonia, nitrous oxide), oxygen and/or sulfur compounds.
These are depth-dependent.
 - scatter plot comparing contamination % (Y-axis) versus completion % (X-axis) between population genome bins.
 Consider colour coding the bins by taxonomic rank e.g. phylum, order, etc.
 Indicate high, medium and low contamination bins (see Canvas for an example).
 - RPKM bubble-plot of each Nitrogen/Sulphur-cycling gene versus taxonomy
 - Pathview pathway diagram 

__References__ (>20)
Necessary citations include:
1. Each software used (including those in the provided bash scripts (MEGAHIT, MetaBAT, checkM, etc.) 
2. Saanich Inlet data papers


## Fin.

__This evaluation rubric should serve as a useful scaffold for your group. How you expand on it will depend on your combined interests and analytic choices. Let the science take you somewhere!__

