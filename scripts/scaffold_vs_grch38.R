# Making all dot plots (combined in one figure) for a scaffold versus all chrom

# Getting the output filename
fname_out <- snakemake@output[[1]]

pdf(fname_out,width=25,height=25)

par(mfrow=c(5,5), oma = c(0, 0, 2, 0))

for (i in 1:length(snakemake@input))
{
    fname <- snakemake@input[[i]]
    # For some very small scaffolds, there are not alignments to some
    # chromosomes and the result is an empty rdot file, which we here skip
    scaffold <- ""
    if (file.size(fname) == 0)
    { 
        next
    }
    dotplot_data <- read.table(fname, header=TRUE)
    chrom <- colnames(dotplot_data)[1]
    split_chr_ids <- unlist(strsplit(chrom,"X"))
    # Only if the chromosome has a number, we have to remove the "X" R appends
    if(length(split_chr_ids) == 2)
    {
        chrom <- split_chr_ids[2]
        chrom <- paste("chr",chrom,sep="") 
    }  
    scaffold <- paste(colnames(dotplot_data)[2])   
    # Plotting
    plot(dotplot_data, type="l", xlab=chrom, ylab="", cex.lab=1.5, cex.axis=1.5)
}

# Adding overall title to plot
mtext(scaffold, outer=TRUE, cex=1)

dev.off()