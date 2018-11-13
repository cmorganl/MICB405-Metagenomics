#!/usr/bin/env Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

get_Length <- function(contigs, proportion) {
  i = 1
  genome_length <- sum(contigs)
  while (i <= length(contigs)) {
    total = sum(contigs[1:i])
    if (total >= genome_length * proportion) {
      return(contigs[i])
    }
    i = i + 1
  }
}

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied!\nUSAGE: (contig_lengths).txt [contig_lengths2.txt] ...", call.=FALSE)
}

vector_size <- length(args) * 100
Proportion <- vector(mode = "double", length = vector_size)
Length <- vector(mode = "double", length = vector_size)
Names <- vector(mode = "double", length = vector_size)

for (n in 1:length(args)) {
  cl_file <- args[n]
  cat(paste0("Reading contig lengths for ", cl_file, "... "))
  
  Contig_lengths <- as.vector(read.table(file = cl_file, header=T)$lgth)

  offset <- (n-1) * 100
  for (i in 1:100) {
    sort(Contig_lengths, decreasing = T)
    Proportion[i + offset] <- i
    # sort(Contig_lengths, decreasing = T) sorts the Contig_lengths vector in decreasing order
    Length[i + offset] <- get_Length(sort(Contig_lengths, decreasing = T), proportion = i/100)
    Names[i + offset] <- cl_file
  }
  cat("done.\n")
}

Length_data <- data.frame(Proportion, Length, Names)

Length_plot <- ggplot(Length_data, aes(x = Proportion, y = Length, col = Names)) + 
  geom_line()
ggsave(filename = "Nx_plot.png", plot = Length_plot, dpi = 500, type = "cairo-png")
