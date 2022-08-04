#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(Biostrings))
f <- commandArgs(T)[1]
sin <- readDNAStringSet(f)
sinfo <- seqlengths(sin)
scon <- DNAStringSet(paste0(sin, collapse = ''))
names(scon) <- basename(f)
writeXStringSet(scon, filepath=paste0(f,".concatenated"))

write.table(as.data.frame(sinfo), file = paste0(f,".concatenated.info"),
            col.names = FALSE, row.names = TRUE, quote = FALSE)




