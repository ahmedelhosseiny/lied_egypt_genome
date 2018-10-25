# Visualizing the distribution of number of positions versus the corresponding 
# sequencing depth.

# Getting input and output filenames
fname_in <- snakemake@input[[1]]
fname_out <- snakemake@output[[1]]

# Plotting
pdf(fname_out)
data <- read.table(fname_in, header=F)
data <- data[order(data[,2]),]
plot(log10(data[,2]),log10(data[,1]),type="p",pch=20,xlab="log10(depth)",ylab="log10(#positions)")
dev.off()
