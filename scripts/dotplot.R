# Making a dotplot from lastz dotplot output format

# Getting input and output filenames
fname_in <- snakemake@input[[1]]
fname_out <- snakemake@output[[1]]

# Plotting
pdf(fname_out)
dots <- read.table(fname_in, header=T)
plot(dots, type="l")
dev.off()
