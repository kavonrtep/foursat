#!/usr/bin/env Rscript

input_dir = commandArgs(T)[1]
four_in1 <-  paste0(input_dir,"/periodicity_fft.RDS")
four_in2 <-  paste0(input_dir,"/periodicity_blast_fft.RDS")

## density dotplot
params = read.table(paste0(input_dir,"/", pattern = "params.txt"), sep = ":", row.names=1)
seq_length = as.numeric(params['sequence_length',])
mfiles = dir(input_dir, pattern = "density.+[.]csv$", full.names = TRUE)
dname = dirname(mfiles)
fname = basename(mfiles)
wl = as.integer(gsub("_[FR].+", "", gsub("^.+_w","", fname)))
orientation = gsub(".csv", "", gsub("^.+[0-9]_","", fname))
pid = gsub("_.+", "", gsub("^.+pid","", fname))

ord = order(wl, as.numeric(pid), orientation)
coef = as.numeric(factor(wl[ord]))
alpha = c("00", "50", "AA", "FF");
alpha = c("00", "EE", "99", "50");
names(alpha) = c("0", 90, 95, 100)

print("loading data for dotplot")
j = 0
m=list()
m_alpha=list()
for (i in ord){
  j = j + 1
  m[[j]] = (as.matrix(read.table(paste0(dname[i],"/",fname[i]), sep="\t"))>0) * coef[j]
  m_alpha[[j]] = m[[j]]
  m_alpha[[j]][,] = 0
  m_alpha[[j]][m[[j]] > 0] = as.integer(pid[i])
}

Mfw = Reduce(pmax, m[orientation[ord]=="F"])
Mrc = Reduce(pmax, m[orientation[ord]=="R"])
MA = Reduce(pmax, m_alpha) 
L = ncol(Mfw)

brv = colorRampPalette(c("#000000","#330000", "#FF0000"))

Fcolors = substr(brv(max(coef) + 1), 2,3)


McolorsR = matrix(paste0("#", alpha[as.character(MA)], "00", "00", Fcolors[Mfw + 1]) ,ncol=L,nrow=L)
McolorsG = matrix(paste0("#", "00", alpha[as.character(MA)],  "00", Fcolors[Mrc + 1]) ,ncol=L,nrow=L)


## color legend
Ylegend = rbind(
  cbind(sort(unique(wl), decreasing=TRUE),"F"),
  c("-",""),
  cbind(sort(unique(wl), decreasing=FALSE),"R")
  )
Xlegend = sort(as.integer(unique(pid)))


R = c(rev(Fcolors),rep("00", length(unique(wl))))
G = c(rep("00", length(unique(wl))+ 1), Fcolors[-1])
lgR= sapply(Xlegend, function(x)paste0("#", alpha[as.character(x)], "00", "00",R))
lgG= sapply(Xlegend, function(x)paste0("#", "00", alpha[as.character(x)], "00",G))

NX=ncol(lgR)
NY=nrow(lgR)
x_label = pretty(c(1,seq_length), 20)
x_at = (pretty(c(1,seq_length), 20)/seq_length)*L
y_at = L-(pretty(c(1,seq_length), 20)/seq_length)*L


print("loading data for periodicity 1")
#############################################################
## fourier spectra
## from sequence
Minfo <-  readRDS(four_in1)
M = Minfo$periodicity_matrix
# normalize
med=apply(M, 2, median, na.rm=TRUE)
M = t(t(M) - med)
R<- nrow(M)
Nintervals <- ncol(M)
xlog1 <-  sapply(1:Nintervals, function(i)approx(y = M[,i], x=1:R, xout=1.005^(1:log(R, 1.005)))$y)
Y1 <- c(1,2,3,4,5,7,10,20,30,40,50,70, 100,150, 200,300,400,500,700,1000,1500,2000,2500,3000,4000,5000,7000,10000, 12000,15000)
Yat1 <- log(Y1,1.005)/log(R,1.005)

print("loading data for periodicity 2")
## from blast
Minfo <-  readRDS(four_in2)
M = Minfo$periodicity_matrix


max_thr=quantile(M, 0.999, na.rm=TRUE)
min_thr=quantile(M, 0.001, na.rm=TRUE)
M[M>max_thr] = max_thr
M[M<min_thr] = min_thr
min_val = apply(M, 2, min, na.rm=TRUE)
M = t(t(M) - min_val)


R<- nrow(M)
Nintervals <- ncol(M)
xlog2 <-  sapply(1:Nintervals, function(i)approx(y = M[,i], x=1:R, xout=1.005^(1:log(R, 1.005)))$y)
Y2 <- c(1,2,3,4,5,7,10,20,30,40,50,70, 100,150, 200,300,400,500,700,1000,1500,2000,2500,3000,4000,5000,7000,10000, 12000,15000)
Yat2 <- log(Y2,1.005)/log(R,1.005)

save.image("tmp.RData")





##############################################################
## plots
widths = c(1,0.25)
heights = c(1,.5)
graphics_width = sum(3000*widths)
graphics_height = sum(3000*heights)

print("ploting")
output_png1 = paste0(input_dir,"/density_dotplot_periodicity_v1.png")
png(output_png1, width = graphics_width, height = graphics_height, pointsize = 50)
layout(matrix(c(1,1,1,3,0,2,0,4), ncol=2), widths = widths, heights = heights)
par(xaxs = 'i', yaxs = 'i', lwd = 5, mar=c(0,5,5,0))
plot(c(1,L ), c(1, L), type = "n", xlab = "", ylab = "", axes=FALSE)
rasterImage(McolorsR, 1, 1, L, L, interpolate = FALSE)
rasterImage(McolorsG, 1, 1, L, L, interpolate = FALSE)
axis(3, at=x_at, labels = x_label)
axis(2, at=y_at, labels = x_label)
box()
## legend
print("ploting2")
par(mar=c(0,5,0,3))
plot(c(0,NX), c(0,NY), type='n', axes = FALSE, xlab="percentage of identity", ylab="word length")
rasterImage(lgR, 0,0,NX,NY, interpolate = FALSE)
rasterImage(lgG, 0,0,NX,NY, interpolate = FALSE)
axis(1,at = (0:(NX -1)) + 0.5, Xlegend, lty = 0 )
axis(2,at = (0:(NY -1)) + 0.5, Ylegend[,1], lty = 0 , las=2)
axis(4, at = c(NY/4, 3*NY/4), c("RC","F"), las=0, lty=0, cex.axis=1)
par(mar=c(0,5,0,0))
# image(t(xlog1)^(1/2), col = rainbow(25), axes = FALSE, ylim = c(0.05, 1))
image(t(xlog1)^(1/2), col = grey.colors(25, rev = TRUE, start=0, end=1, gamma = 1), axes = FALSE, ylim = c(0.05, 1), useRaster = TRUE)
axis(2,at = Yat1, labels = Y1, las = 2)
prof_x = rowSums(xlog1)
prof_y = seq_along(prof_x)/length(prof_x)
par(mar=c(0,0,0,5))
plot(prof_x, prof_y, ylim=c(0.05,1), type='l', axes = FALSE, xlab="", ylab="")
dev.off()


output_png2 = paste0(input_dir,"/density_dotplot_periodicity_v2.png")
png(output_png2, width = graphics_width, height = graphics_height, pointsize = 50)
layout(matrix(c(1,1,1,3,0,2,0,4), ncol=2), widths = widths, heights = heights)
par(xaxs = 'i', yaxs = 'i', lwd = 5, mar=c(0,5,5,0))
plot(c(1,L ), c(1, L), type = "n", xlab = "", ylab = "", axes=FALSE)
rasterImage(McolorsR, 1, 1, L, L, interpolate = FALSE)
rasterImage(McolorsG, 1, 1, L, L, interpolate = FALSE)
axis(3, at=x_at, labels = x_label)
axis(2, at=y_at, labels = x_label)
box()
## legend
par(mar=c(0,5,0,3))
plot(c(0,NX), c(0,NY), type='n', axes = FALSE, xlab="percentage of identity", ylab="word length")
rasterImage(lgR, 0,0,NX,NY, interpolate = FALSE)
rasterImage(lgG, 0,0,NX,NY, interpolate = FALSE)
axis(1,at = (0:(NX -1)) + 0.5, Xlegend, lty = 0 )
axis(2,at = (0:(NY -1)) + 0.5, Ylegend[,1], lty = 0 , las=2)
axis(4, at = c(NY/4, 3*NY/4), c("RC","F"), las=0, lty=0, cex.axis=1)
par(mar=c(0,5,0,0))
# image(t(xlog2)^(1/2), col = rainbow(25), axes = FALSE, ylim = c(0.05, 1))
image(t(xlog2), col = grey.colors(25, rev = TRUE, start=0, end=1, gamma = 1), axes = FALSE, ylim = c(0.05, 1), useRaster = TRUE)
axis(2,at = Yat2, labels = Y2, las = 2)
prof_x = rowSums(xlog2)
prof_y = seq_along(prof_x)/length(prof_x)
par(mar=c(0,0,0,5))
plot(prof_x, prof_y, ylim=c(0.05,1), type='l', axes = FALSE, xlab="", ylab="")

dev.off()



