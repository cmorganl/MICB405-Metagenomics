#!/usr/bin/env Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied!\n\n\tUSAGE: assembly1_Nx.txt [assembly2_Nx.txt] ...", call.=FALSE)
}

# vector_size <- length(args) * 100
lengths <- vector(mode = "double")
proportions <- vector(mode = "double")
Assembly <- vector(mode = "character")

for (n in 1:length(args)) {
  nx_file <- args[n]
  cat(paste0("Reading Nx stats for ", nx_file, "... "))
  
  Nx_stats <- read.table(file = nx_file,
                         header=T,
                         sep=',')
  
  proportions <- c(proportions, Nx_stats$Genome_proportion[2:101])
  lengths <- c(lengths, Nx_stats$Contig_size[2:101])
  Assembly <- c(Assembly, rep(basename(nx_file), 100))
  cat("done.\n")
}

nx_data <- data.frame(Assembly, proportions, lengths)

nx_curve <- ggplot(nx_data, aes(x=proportions, y=lengths, col=Assembly)) +
  geom_line()+
  xlab("Assembly Proportion") +
  ylab("Contig Length")

ggsave(filename = "Nx_plot.png", plot = nx_curve, dpi = 500, type = "cairo-png")
