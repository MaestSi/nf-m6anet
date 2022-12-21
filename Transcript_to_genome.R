#!/usr/bin/env Rscript
### load input variables ###
args = commandArgs(trailingOnly=TRUE)

for(v in args)
{
  vTmp <- strsplit(v,"=")[[1]]
  assign(vTmp[[1]],vTmp[[2]])
}

library(IRanges)
library(ensembldb)
library(GenomicRanges)
library(stringr)
library(Biostrings)
library(parallel)

Lift_over_m6anet <- function(input_file, prob_mod_thr, genome_gtf, output_file, output_file_genome, mccores) {
  if (file.exists(input_file) && file.info(input_file)$size != 0) {
    data_m6anet <- read.table(input_file, header = TRUE, sep=",") 
    m6anet <- data.frame("TranscriptID" = data_m6anet[,1],
                         "Start" = data_m6anet[,2] - 1,
                         "End" = data_m6anet[,2],
                         "Status" = data_m6anet[,4],
                         "Prob_mod" = data_m6anet[,4]
    )
    m6anet$Status <- ifelse(!is.nan(m6anet$Status) & !is.na(m6anet$Status) & m6anet$Status > as.numeric(prob_mod_thr), "Mod", "Unmod")
    if (length(m6anet$Status > as.numeric(prob_mod_thr)) > 0) {
      m6anet <- m6anet[which(m6anet$Status == "Mod"), ]
      write.table(m6anet, file = output_file, quote = F, sep = "\t", row.names = F)
      # Creation Edb Database from genome GTF
      EnsDb <- suppressWarnings(suppressMessages(ensDbFromGtf(gtf = genome_gtf)))
      edb <- EnsDb(EnsDb)
      ## Lift-over + output bed
      test_m6anet <- IRanges(start = m6anet[,2], end = m6anet[,3], names = c(m6anet[,1]))
      num_rows_chunk <- 1
      mc.cores <- as.numeric(mccores)
      if (length(test_m6anet) < num_rows_chunk) {
        test_m6anet_split <- list(test_m6anet)
      } else {
        test_m6anet_split <- split(test_m6anet, rep(seq(from = 1, to = ceiling(length(test_m6anet)/num_rows_chunk)), each = num_rows_chunk)[1:length(test_m6anet)])
      }
      #remove tx version
      test_m6anet_split <- lapply(test_m6anet_split, function(x) { names(x) = gsub(x = names(x), pattern = "\\.\\d*", replacement = ""); return(x)})
      #remove tx not available in the gtf file
      txdb <- makeTxDbFromGFF(file = genome_gtf, format = "auto")
      tx_from_txdb <- mcols(GenomicFeatures::transcripts(txdb))$tx_name
      unannotated_tx <- setdiff(unlist(lapply(test_m6anet_split, names)), tx_from_txdb)
      if (length(unannotated_tx) > 0) {
        unannotated_tx_ind <- which(unlist(lapply(test_m6anet_split, function(x) names(x) %in% unannotated_tx)))
        test_m6anet_split <- test_m6anet_split[-unannotated_tx_ind]
        names(test_m6anet_split) <- 1:length(test_m6anet_split)
      }
      #perform lift-over
      tmp1 <- vector(mode = "list", length = length(test_m6anet_split))
      names(tmp1) <- 1:length(test_m6anet_split)
      tmp <- vector(mode = "list", length = length(test_m6anet_split))
      names(tmp) <- 1:length(test_m6anet_split)
      ind_retry <- 1:length(test_m6anet_split)
      while(any(unlist(lapply(tmp, is.null)))) {
        cat(sprintf("Starting new iteration for m6Anet; %d sites missing\n", length(which(unlist(lapply(tmp, is.null))))))
        tmp1 <- tmp1[ind_retry]
        tmp1 <- mclapply(test_m6anet_split[ind_retry], function(x) {
          tryCatch({
            coordinate_m6anet_unlisted <- unlist(transcriptToGenome(x, edb))
            return(coordinate_m6anet_unlisted)
          }, warning = function(w) {
            print("Warning")
            return(NULL)
          }, error = function(e) {
            print("Error")
            return(NULL)
          }
          )}, mc.cores = mc.cores)
        ind_retry_tmp <- unlist(lapply(tmp1, function(x) is.null(x)))
        if (length(ind_retry_tmp) > 0) {
          ind_retry <- names(which(ind_retry_tmp))
        } else {
          ind_retry <- c()  
        }
        ind_ok_tmp <- unlist(lapply(tmp1, function(x) !is.null(x)))
        if (length(ind_ok_tmp) > 0) {
          ind_ok <- names(which(ind_ok_tmp))
        } else {
          ind_ok <- c() 
        }
        tmp[ind_ok] <- tmp1[ind_ok]
        if (length(ind_retry) > 0) {
          tmp1 <- tmp1[ind_retry]
        }
      }
      
      coordinate_m6anet_unlisted <- unlist(as(tmp, "GRangesList"))
      df_m6anet <- as.data.frame(unname(coordinate_m6anet_unlisted[, c(0, 2, 4, 5)]))[, c(1:3, 5, 6, 7, 8)]
      
      tmp_rownames <- paste0(df_m6anet[, 6], "_", df_m6anet[, 7], "_", df_m6anet[, 5])
      dup_names <- names(which(table(tmp_rownames) > 1))
      ind_dup <- which(tmp_rownames %in% dup_names)
      ind_dup_rm <- ind_dup[which(duplicated(df_m6anet[ind_dup, c(5, 6, 7)]))]
      if (length(ind_dup_rm) > 0) {
        df_m6anet <- df_m6anet[-ind_dup_rm, ]
      }
      rownames(df_m6anet) <- paste0(df_m6anet[, 6], "_", df_m6anet[, 7], "_", df_m6anet[, 5])
      
      rownames(m6anet) <- paste0(m6anet[, 2], "_", m6anet[, 3], "_", m6anet[, 1])
      
      df_m6anet$Status <- m6anet[rownames(df_m6anet), 4]
      df_m6anet$Prob_Mod <- m6anet[rownames(df_m6anet), 5]
      df_m6anet_final <- df_m6anet[, c(1, 2, 3, 5, 9, 4)]
      colnames(df_m6anet_final) <- c("Chr", "Start", "End", "Transcript_id", "Prob_mod", "Strand")
      write.table(df_m6anet_final, file = output_file_genome, quote = F, sep = "\t", row.names = F)
      return(df_m6anet_final)
    } else {
      message("No m6A+ sites were detected")
    }
  }
  else {message("m6anet output file does not exist.")}
}

Filtered_m6anet <- Lift_over_m6anet(input_file, prob_mod_thr, genome_gtf, output_file, output_file_genome, mccores)
