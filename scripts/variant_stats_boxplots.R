# Making boxplots for snp stat numbers

# Getting input and output filenames
fname_out_boxplots <- snakemake@output[[1]]
fname_out_het_missing <- snakemake@output[[2]]
fname_out_corplot <- snakemake@output[[3]]
fname_out_indelhist <- snakemake@output[[4]]
fname_out_varsperchrom <- snakemake@output[[5]]
fname_insnps <- snakemake@input[[1]]
fname_imiss <- snakemake@input[[2]]
fname_idepth <- snakemake@input[[3]]
fname_het <- snakemake@input[[4]]
fname_indelhist <- snakemake@input[[5]]
fname_varsperchrom <- snakemake@input[[6]]

# Plotting
pdf(fname_out_boxplots)

par(mfrow=c(2,2))
data <- read.table(fname_insnps, header=TRUE)
data <- data[order(data$INDV),] 
boxplot(data$N_SNPS,ylab="Number of SNPs",col=colors()[250])
stripchart(data$N_SNPS, vertical=TRUE, method="jitter", add=TRUE, pch=20, col=2)
tail(data)
num_snps <- data$N_SNPS

data <- read.table(fname_imiss, header=TRUE)
data <- data[order(data$INDV),] 
boxplot(100*data$F_MISS,ylab="Missing SNPs [%]",col=colors()[250])
stripchart(100*data$F_MISS, vertical=TRUE, method="jitter", add=TRUE, pch=20, col=3)
tail(data)
missing_percent <- 100*data$F_MISS

data <- read.table(fname_idepth, header=TRUE)
data <- data[order(data$INDV),] 
boxplot(data$MEAN_DEPTH,ylab="Mean Depth",col=colors()[250])
stripchart(data$MEAN_DEPTH, vertical=TRUE, method="jitter", add=TRUE, pch=20, col=4)
tail(data)
mean_depth <- data$MEAN_DEPTH

data <- read.table(fname_het, header=TRUE)
data <- data[order(data$INDV),] 
boxplot(data$F,ylab="Heterozygosity [F Coefficient]",col=colors()[250])
stripchart(data$F, vertical=TRUE, method="jitter", add=TRUE, pch=20, col=6)
tail(data)
het <- data$F

dev.off()


# Make a color vector
# Default color(for Pagani samples) is orange
cols = rep(rainbow(14)[2],length(num_snps))
# Delta Egypt samples are red
cols[102:107] <- rainbow(14)[1]
# Upper Egypt samples are yellow
cols[108:110] <- rainbow(14)[3]
# Egyptref black
cols[101] <- 1


pdf(fname_out_het_missing)

plot(missing_percent,het,xlab="Missing SNPs [%]",ylab="Heterozygosity [F Coefficient]",pch=21, bg=cols)
legend("topleft",c("EGD","EGP","EGU","EGYPTREF"),fill=c(rainbow(14)[1],rainbow(14)[2],rainbow(14)[3],1),cex=0.8)

dev.off()


pdf(fname_out_corplot)

# Correlation panel
panel.cor <- function(x, y){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- round(cor(x, y), digits=2)
    txt <- paste0("R = ", r)
    cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * abs(r))
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch = 19, col = cols)
}
# Create the plots
pairs(cbind(num_snps,missing_percent,mean_depth,het), 
      lower.panel = panel.cor,
      upper.panel = upper.panel)
      
dev.off()


pdf(fname_out_indelhist)

data <- read.table(fname_indelhist, header=TRUE)
plot(data$LENGTH,log10(data$COUNT),xlab="Length",ylab="log10(Number)",col=8,type="h")

dev.off()

pdf(fname_out_varsperchrom)

data <- read.table(fname_varsperchrom, header=FALSE)
barplot(data[,1], xlab="Chromosome", ylab="Number", col=8, las=2, names=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"))

dev.off()
