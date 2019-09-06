# Visualizations with respect to proxy SNPs of tag SNPs

# Getting input and output filenames
fname_in <- snakemake@input[[1]]
fname_out_hist <- snakemake@output[[1]]

data <- read.table(fname_in, header=T)
dim(data)
data <- data[data$EGP_N_CHR>=200,]

# Plotting
pdf(fname_out_hist)
hist(data$EUR_N_PROXY,breaks=100,xlab="number of proxies", main="")
dev.off()

