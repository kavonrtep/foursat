#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Biostrings)
  library(Rfast)
  library(parallel)
  library(optparse)
})

parser = OptionParser()
option_list = list(
  make_option(c("-i", "--input"), action = "store", type = "character",
              help = "help", default = NULL),
  make_option(c("-t", "--input_type"), action = "store" , type = "character",
              default = "FASTA", help = "input type is eithet FASTA of MATRIX"),
  make_option(c("-c", "--cpu"), type="integer", default=20, 
              help="Number of cpu [default %default]",
              metavar="number"),
  make_option(c("-w", "--window"), type="integer", default=10000, 
              help="window size (number of nt)",
              metavar="number"),

  make_option(c("-s", "--step"), type="integer", default=200, 
              help="Step size (number of nt) [default %default]",
              metavar="number"),
  make_option(c("-o", "--output", action='store', type='character', help=''))
)

description = paste(strwrap(""))
epilogue = ""
parser = OptionParser(option_list = option_list, epilogue = epilogue, description = description, 
                      usage = "usage: %prog COMMAND [OPTIONS]")
opt = parse_args(parser, args = commandArgs(TRUE))


convert.fft <- function(cs, sample.rate=1) {
  cs <- cs / length(cs) # normalize

  distance.center <- function(c)signif( Mod(c),        4)
  angle           <- function(c)signif( 180*Arg(c)/pi, 3)
  df <- data.frame(cycle    = 0:(length(cs)-1),
                   freq     = 0:(length(cs)-1) * sample.rate / length(cs),
                   strength = sapply(cs, distance.center),
                   delay    = sapply(cs, angle))
  df$monomer <- nrow(df)/df$cycle
  df
}

dna2vector <- function(s, method="basic"){
  ## assume single sequence!
  coding1 <- c(A = -1, C = 1, G = -2 , T = -2)
  coding2 <- c(A = -1, C = 0, G = 0 , T = 1)
  coding3 <- c(A = 0, C = -1, G = 1 , T = 0)

  if (method == "basic"){
    x <- as.numeric(as.factor(unlist(strsplit(as.character(s), split=""))))
  }
  if (method == "pupy"){
    x <- coding1[unlist(strsplit(as.character(s), split=""))]
  }
  if (method == "cumulative"){
    x <- cumsum(coding1_corrected[unlist(strsplit(as.character(s), split=""))])
    ## detrend:
    x <- x- lowess(x, f= 0.01)$y
  }
  if (method=="pupy-complement"){
    x1 <- coding2[unlist(strsplit(as.character(s), split=""))]
    x2 <- coding3[unlist(strsplit(as.character(s), split=""))]
    x <- cbind(x1, x2)
  }
  return(x)
}

get_periodicity <- function(x, window=5000, step=10, xrange=seq(1, 15000, 1), base = 1.2){
  L <- length(x)
  j <- 0
  y <- numeric(length(xrange))
  ## for (i in seq(0,window,step)){
  for (i in unique((round(base^(seq(0,log(window,base=base))))))){
    j <- j + 1
    df <- convert.fft(fft(x[1:(L-i)]))
    y <- y + approx(x = df$monomer[-1], y=df$strength[-1], xout = xrange)$y
  }
  y/j
}


ACF2 <- function(s, L){
  x <- unlist(strsplit(tolower(unname(as.character(s))), ""))
  S1_4 <- list()
  for (j in c("a","c","g","t")){
    A <- as.numeric(x == j)
    S1_4[[j]] <- acf(A, lag.max = L, plot=FALSE)$acf[-1, 1, 1]
  }
  ##do.call(cbind, S1_4)
  Reduce("+", S1_4)
}

## input fasta file
if (opt$input_type == "FASTA"){
  s_long=readDNAStringSet(opt$input)
  L=nchar(s_long)
  ## get periodicity using fft
  s = substring(s_long,
                first = seq(1, L - opt$window, opt$step),
                last = seq(1, L - opt$window, opt$step) + opt$window)
  xbig <- lapply(s, dna2vector, method = "basic")
  periodicity <- mclapply(xbig, FUN = get_periodicity, mc.cores = opt$cpu, window=200)
  periodicity_matrix <- do.call(cbind, periodicity )
  saveRDS(list(periodicity_matrix = periodicity_matrix, window=opt$window, step_size=opt$step, seq_length = L), file=opt$output)

}else{
  if (opt$input_type == "MATRIX"){
    xmat <- as.matrix(read.table(opt$input))
    xlist <- lapply(seq_len(ncol(xmat)), function(i)xmat[,i])
    periodicity <- mclapply(xlist, FUN = get_periodicity, mc.cores = opt$cpu, window=200, xrange=seq_len(length(xlist[[1]])))
    periodicity_matrix <- do.call(cbind, periodicity)
    ## normalization
    med=apply(periodicity_matrix, 2, median, na.rm=TRUE)
    periodicity_matrix_norm <- t(t(periodicity_matrix) - med)
    saveRDS(list(periodicity_matrix_norm = periodicity_matrix, window=NA, step_size=NA, seq_length = NA), file=opt$output)
  }
}






#periodicity = do.call(rbind,mclapply(s, FUN = ACF2, L=10000, mc.cores = 40))
