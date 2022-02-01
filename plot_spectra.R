#!/usr/bin/enr Rcsript

fin <-  commandArgs(T)[1]
window <- commandArgs(T)[2]
step_size <- commandArgs(T)[3]

if(FALSE){
  fin <-  "/mnt/ceph/454_data/MinION/analysis/centromere_assembly_Cameor/211209_CEN6_ANALYSIS/tmp/plots_sat_arrays_Petr/CEN6_ver_211209_79000-79700k_sliding_W50k_s0.1k.fasta_periodicity_fft.RDS"
  window <- 50000
  step_size <- 100
  fout <- "test.png"


  fin <-  "/mnt/ceph/454_data/MinION/analysis/centromere_assembly_Cameor/211209_CEN6_ANALYSIS/tmp/plots_sat_arrays_Petr/CEN6_ver_211209_79000-79700k_sliding_W30k_s2k.fasta_periodicity_fft.RDS"
  window <- 30000
  step_size <- 2000
  fout <- "test.png"

  fin <-  "/mnt/ceph/454_data/MinION/analysis/centromere_assembly_Cameor/211209_CEN6_ANALYSIS/tmp/plots_sat_arrays_Petr/CEN6_ver_211209_79000-79700k_sliding_W10k_s0.1k.fasta_periodicity_fft.RDS"
  fin <-  "/mnt/ceph/454_data/MinION/analysis/centromere_assembly_Cameor/211209_CEN6_ANALYSIS/tmp/plots_sat_arrays_Petr/CEN6_ver_211209_79000-79700k_sliding_W10k_s0.1k.fasta_periodicity_fft.RDS"

  window <- 10000
  step_size <- 100
  fout <- "test.png"




}

M <-  readRDS(fin)
R<- nrow(M)
Nintervals <- ncol(M)
xlog <-  sapply(1:Nintervals, function(i)approx(y = M[,i], x=1:R, xout=1.005^(1:log(R, 1.005)))$y)

xaxs <-  1.005^(1:log(R, 1.005))

Y <- c(1,2,3,4,5,7,10,20,30,40,50,70, 100,150, 200,300,400,500,700,1000,1500,2000,2500,3000,4000,5000,7000,10000, 12000,15000)
Yat <- log(Y,1.005)/log(R,1.005)

library(baseline)

substract_baseline <- function(S, method='lowpas'){
  S0 <-  S[!is.na(rowSums(S)),]
  bo <- baseline(t(S0), method = method)
  bline <-  getBaseline(bo)
  Sb <-  S
  Sb[!is.na(rowSums(S)),] <- Sb[!is.na(rowSums(S)),] - t(bline)
  return(Sb)
}

xlog_norm = substract_baseline(xlog, method = "irls")

X <-  pretty(step_size*(1:Nintervals) + window/2, 25)
i <- (X-window/2)/step_size
Xat <- i/Nintervals

png(fout, width = 5000, height = 2000, pointsize = 20)
#image(t(xlog)^(1/4), col = hcl.colors(25, "Grays", rev = TRUE))
par(mar=c(4,4,1,1), mfrow = c(1,1))
image(t(xlog)^(1/2), col = rainbow(25), axes = FALSE)
axis(2,at = Yat, labels = Y, las = 2)
axis(1,at = Xat, labels = X, las = 1, line=-1)
dev.off()

