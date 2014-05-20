#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data <- read.csv(args[1], sep = "\t", header=TRUE, check.names=FALSE)
data <- data[-1,]

svg(args[2])
par(mar = c(10, 5, 2, 2))  
barplot(as.matrix(data), beside=TRUE, col=rainbow(4), las=2, cex.names=0.5, ylab=args[4])
legend("topleft", c("RNAfold", "Deigan", "Zarringhalam", "Washietl"), fill=rainbow(4), bty="n", cex=0.5)

for (i in 1:dim(data)[2]) {
  for (j in 2:4) {
    data[j,i] <- data[j,i] - data[1,i]
  }
}

data <- data[-1,]

svg(args[3])
par(mar = c(10, 5, 2, 2))
barplot(as.matrix(data), beside=TRUE, col=rainbow(3), las=2, cex.names=0.5, ylab=args[4])
legend("topleft", c("Deigan", "Zarringhalam", "Washietl"), fill=rainbow(3), bty="n", cex=0.5)

