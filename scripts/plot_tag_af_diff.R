# Visualizing the difference between taq SNP allele frequencies in Europeans and Egyptians.

# Getting input and output filenames
fname_in <- snakemake@input[[1]]
fname_out_hist <- snakemake@output[[1]]
fname_out_scatter <- snakemake@output[[2]]
fname_out_missing <- snakemake@output[[3]]

data <- read.table(fname_in, header=T)
dim(data)
data <- data[data$EGP_N_CHR>=200,]

# Plotting
pdf(fname_out_hist)
hist(data$EUR_ALT_AF-data$EGP_ALT_AF,breaks=100,xlab="European AF - Egyptian AF", main="")
dev.off()

pdf(fname_out_scatter)
plot(c(-0.1,1.1),c(-0.1,1.1),type="l",lty=1,col="gray",xlab="European AF",ylab="Egyptian AF", main="",xlim=c(0,1),ylim=c(0,1))
add=TRUE
abline(h=0.05,col="gray")
abline(v=0.05,col="gray")
abline(a=0.1,b=1,col="gray")
abline(a=-0.1,b=1,col="gray")
abline(a=0.2,b=1,col="gray")
abline(a=-0.2,b=1,col="gray")
abline(a=0.3,b=1,col="gray")
abline(a=-0.3,b=1,col="gray")
points(data$EUR_ALT_AF,data$EGP_ALT_AF,pch=20,cex=0.25)
dev.off()

pdf(fname_out_missing)
data <- data[data$EGP_ALT_AF==0,]
hist(data$EUR_ALT_AF,breaks=100,xlab="European AF for tag SNPs not occuring in Egyptians",main="",xlim=c(0,1))
dev.off()
