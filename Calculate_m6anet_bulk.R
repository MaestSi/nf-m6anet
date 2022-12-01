#!/usr/bin/env Rscript
### load input variables ###
args = commandArgs(trailingOnly=TRUE)

if (args[1] == "-h" | args[1] == "--help") {
  cat("", sep = "\n")
  cat(paste0("Usage: Rscript Calculate_m6anet_bulk.R m6anet_output_file=<data.result.csv> report_file=<report.txt>"), sep = "\n")
  cat(paste0("<data.result.csv>: m6Anet output file"), sep = "\n")
  cat(paste0("<report.txt>: report file with m6A bulk levels"), sep = "\n")
  stop(simpleError(sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "))))
}

for(v in args)
{
  vTmp <- strsplit(v,"=")[[1]]
  assign(vTmp[[1]],vTmp[[2]])
}

suppressMessages(library(Biostrings))
suppressMessages(library(GenomicRanges))

m6anet_output <- read.table(m6anet_output_file, sep = ",", header = TRUE)
sum_mod_ratio <- sum(m6anet_output$mod_ratio)
sum_prob_mod <- sum(m6anet_output$probability_modified)

num_DRACH <- dim(m6anet_output)[1]
bulk_m6A_mod_A_ratio <- sum_mod_ratio/num_DRACH
bulk_m6A_mod_A_prob <- sum_prob_mod/num_DRACH

cat(sprintf("Analysed %d DRACH sites\n", num_DRACH), file = report_file)
cat(sprintf("m6A bulk level estimate based on modification ratio: %.2f%%\n", 100*bulk_m6A_mod_A_ratio), file = report_file, append = TRUE)
cat(sprintf("m6A bulk level estimate based on modification probability: %.2f%%\n", 100*bulk_m6A_mod_A_prob), file = report_file, append = TRUE)
