# Making all dot plots (combined in one figure) for a scaffold versus all chrom

# Getting the output filename
fname_out <- snakemake@output[[1]]

# Geting the lengths of the GRCh38 chromosomes
numbases_grch38_file <- snakemake@input[[1]]
num_bases_grch38 <- read.table(numbases_grch38_file,header=FALSE,row.names=1,skip=1)
numbases_assembly_file <- snakemake@input[[2]]
num_bases_assembly <- read.table(numbases_assembly_file, header=TRUE)

pdf(fname_out,width=25,height=25)

par(mfrow=c(5,5), oma = c(0, 0, 2, 0))

for (i in 3:length(snakemake@input))
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
    chrom_wo_chr <- chrom
    split_chr_ids <- unlist(strsplit(chrom,"X"))
    # Only if the chromosome has a number, we have to remove the "X" R appends
    if(length(split_chr_ids) == 2)
    {
        chrom_wo_chr <- split_chr_ids[2]
        chrom <- paste("chr",chrom_wo_chr,sep="") 
    }  
    scaffold <- paste(colnames(dotplot_data)[2])   
    # Getting the lengths of the scaffold and corresponding chromosome in order
    # to scale axis accordingly
    len_scaff <- num_bases_assembly[scaffold,1]
    len_chr <- num_bases_grch38[chrom_wo_chr,1]
    # Plotting
    plot(dotplot_data, type="l", xlab=chrom, ylab="", cex.lab=1.5, cex.axis=1.5, xlim=c(0,len_chr), ylim=c(0,len_scaff))
}

# Adding overall title to plot
mtext(scaffold, outer=TRUE, cex=1)

dev.off()