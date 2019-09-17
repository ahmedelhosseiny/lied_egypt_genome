# Visualizing the difference between taq SNP allele frequencies in Europeans and Egyptians.

# Getting input and output filenames
fname_in <- snakemake@input[[1]]
fname_out_hist <- snakemake@output[[1]]
fname_out_scatter <- snakemake@output[[2]]
fname_out_missing <- snakemake@output[[3]]
fname_out_boxplot_afdiff <- snakemake@output[[4]]

data <- read.table(fname_in, header=T)
dim(data)
data <- data[data$EGP_N_CHR>=200,]

# Plotting
pdf(fname_out_hist)
hist(data$EUR_ALT_AF-data$EGP_ALT_AF,breaks=100,xlab="European AF - Egyptian AF", main="")
dev.off()

pdf(fname_out_scatter)
plot(c(-0.1,1.1),c(-0.1,1.1),type="l",lty=1,col="gray",xlab="European AF",ylab="Egyptian AF", main="",xlim=c(0,1),ylim=c(0,1),lwd=5)
add=TRUE
blue <- rainbow(10, alpha=0.2)[7]
red <- rainbow(10, alpha=0.2)[1]
abline(h=0.05,col=red,lwd=5)
abline(v=0.05,col=blue,lwd=5)
abline(a=0.1,b=1,col=red,lwd=5)
abline(a=-0.1,b=1,col=blue,lwd=5)
abline(a=0.2,b=1,col=red,lwd=5)
abline(a=-0.2,b=1,col=blue,lwd=5)
abline(a=0.3,b=1,col=red,lwd=5)
abline(a=-0.3,b=1,col=blue,lwd=5)
points(data$EUR_ALT_AF,data$EGP_ALT_AF,pch=21,cex=0.5,col="grey40",bg="grey90")
dev.off()

pdf(fname_out_boxplot_afdiff,width=5, height=7)
data_wo_missing <- data[data$EGP_ALT_AF!=0,]
c1 <- rainbow(10)[4]
c2 <- rainbow(10, alpha=0.2)[4]
c3 <- rainbow(10, v=0.7)[4]
boxplot(data_wo_missing$EUR_ALT_AF-data_wo_missing$EGP_ALT_AF,ylab="European AF - Egyptian AF",cex.axis=0.9,outcex=0.5,col=c2,medcol=c3,whiskcol=c1,staplecol=c3,boxcol=c3,outcol=c2,pch=19)
dev.off()

pdf(fname_out_missing)
data <- data[data$EGP_ALT_AF==0,]
hist(data$EUR_ALT_AF,breaks=100,xlab="European AF for tag SNPs not occuring in Egyptians",main="",xlim=c(0,1))
dev.off()
