#!/usr/bin/env Rscript
### load input variables ###
args = commandArgs(trailingOnly=TRUE)

if (args[1] == "-h" | args[1] == "--help") {
  cat("", sep = "\n")
  cat(paste0("Usage: Rscript Calculate_m6anet_bulk.R m6anet_output_file=<data.site_proba.csv> report_file=<report.txt> prob_mod_thr=<prob_mod_thr>"), sep = "\n")
  cat(paste0("<data.site_proba.csv>: m6Anet output file"), sep = "\n")
  cat(paste0("<report.txt>: report file with m6A bulk levels"), sep = "\n")
  cat(paste0("<prob_mod_thr>: probability modification threshold for calling a site as m6A+"), sep = "\n")
  prob_mod_thr
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
num_m6A_sites <- length(which(m6anet_output$probability_modified > as.numeric(prob_mod_thr)))

num_DRACH <- dim(m6anet_output)[1]
bulk_m6A_mod_A_ratio <- sum_mod_ratio/num_DRACH
bulk_m6A_mod_A_prob <- sum_prob_mod/num_DRACH
bulk_m6A_fraction <- num_m6A_sites/num_DRACH

cat(sprintf("Analysed %d DRACH sites\n", num_DRACH), file = report_file)
cat(sprintf("m6A bulk level estimate based on fraction of m6A+ sites: %.2f%%\n", 100*bulk_m6A_fraction), file = report_file, append = TRUE)
cat(sprintf("m6A bulk level estimate based on modification ratio: %.2f%%\n", 100*bulk_m6A_mod_A_ratio), file = report_file, append = TRUE)
cat(sprintf("m6A bulk level estimate based on modification probability: %.2f%%\n", 100*bulk_m6A_mod_A_prob), file = report_file, append = TRUE)