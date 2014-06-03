#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

data <- read.csv(args[1], sep = "\t", header=TRUE, check.names=FALSE)

data <- data[order(data$Length),]

values <- data[,c(3,4,5,6)]

colors = rainbow(3)

pdf(args[2])
plot(data[,2], data[,3], type = 'l', lwd=2, ylim=c(min(values), max(values)), ylab = "Runtime in seconds", xlab="Sequence length")
lines(data[,2], data[,4], col=colors[1])
lines(data[,2], data[,5], col=colors[2])
lines(data[,2], data[,6], col=colors[3])
legend("bottomright", c("RNAfold", "Deigan", "Zarringhalam", "Washietl"), lty=c(1,1,1,1), lwd=c(2,1,1,1), col=c("black", colors[1], colors[2], colors[3]), cex=.5)
grid()

values <- data[,c(6,7)]

pdf(args[3])
plot(data[,2], data[,7], type = 'l', ylim=c(min(values), max(values)), ylab = "Runtime in seconds", xlab="Sequence length")
#lines(data[,2], data[,7], lty=2)
#legend("bottomright", c("Folding recursion", "Perturbation vector calculation"), lty=c(1,2), lwd=c(1,1), cex=.5)
grid()
