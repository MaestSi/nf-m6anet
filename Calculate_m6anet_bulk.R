#!/usr/bin/env Rscript
### load input variables ###
args = commandArgs(trailingOnly=TRUE)

if (args[1] == "-h" | args[1] == "--help") {
  cat("", sep = "\n")
  cat(paste0("Usage: Rscript Calculate_m6anet_bulk.R transcriptome_file=<transcriptome.fasta> m6anet_output_file=<data.result.csv> report_file=<report.txt>"), sep = "\n")
  cat(paste0("<transcriptome.fasta>: transcriptome fasta file"), sep = "\n")
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

chrs <- readDNAStringSet(transcriptome_file, format = "fasta")

RRACH <- c("AAACA","AAACT","AAACC","GAACA","GAACT","GAACC","GGACA","GGACT","GGACC","GAACA","GAACT","GAACC")
RRACH_bed <- data.frame()

for (x in RRACH) {
  match <- as.data.frame(unlist(vmatchPattern(x, chrs)))
  match$strand <- rep("+", nrow(match))
  match <-match[,c(4,1,2,5)]
  RRACH_bed <- rbind(RRACH_bed, match)
}
colnames(RRACH_bed) <- c("chr", "start", "end", "strand")
num_RRACH <- dim(RRACH_bed)[1]

m6anet_output <- read.table(m6anet_output_file, sep = ",", header = TRUE)
m6anet_output_RRACH <- m6anet_output[which(m6anet_output$kmer %in% RRACH), ]
sum_mod_ratio <- sum(m6anet_output_RRACH$mod_ratio)
sum_prob_mod <- sum(m6anet_output_RRACH$probability_modified)

bulk_m6A_mod_A_ratio <- sum_mod_ratio/num_RRACH
bulk_m6A_mod_A_prob <- sum_prob_mod/num_RRACH

cat(sprintf("Found %d RRACH motifs\n", num_RRACH), file = report_file)
cat(sprintf("m6A bulk level estimate based on modification ratio: %.2f%%\n", 100*bulk_m6A_mod_A_ratio), file = report_file, append = TRUE)
cat(sprintf("m6A bulk level estimate based on modification probability: %.2f%%\n", 100*bulk_m6A_mod_A_prob), file = report_file, append = TRUE)
