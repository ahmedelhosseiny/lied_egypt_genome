# Getting input and output filenames
fname_in <- snakemake@input[[1]]
fname_out <- snakemake@output[[1]]

read_numbers <- read.table(fname_in,header=TRUE)

# Plotting
pdf(fname_out)
par(mfrow=c(1,3))

boxplot(read_numbers$REFERENCE_UNMAPPED,ylab="Number reads unmapped w.r.t. reference and GATK bundle sequences",col=colors()[250],log="y")
stripchart(read_numbers$REFERENCE_UNMAPPED, vertical=TRUE, method="jitter", add=TRUE, pch=20, col=2,log="y")

boxplot(read_numbers$THEROF_ASSEMBLY_MAPPED,ylab="Number unmapped reads mapping to assembly",col=colors()[250],log="y")
stripchart(read_numbers$THEROF_ASSEMBLY_MAPPED, vertical=TRUE, method="jitter", add=TRUE, pch=20, col=3, log="y")

boxplot(read_numbers$PERCENT,ylab="Percent of previously unmapped reads mapping to assembly",col=colors()[250])
stripchart(read_numbers$PERCENT, vertical=TRUE, method="jitter", add=TRUE, pch=20, col=4)


dev.off()