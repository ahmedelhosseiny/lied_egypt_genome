# Writes the first 10 genotype principal components to file for later inclusion 
# into the linear model. Computes various visualizations, e.g. PCs, heatmaps and
# scree plot.

#library("lattice")

# Get the name of the input file
filename <- snakemake@input[[1]]
annotation_file <- snakemake@input[[2]]
path <- snakemake@params[[1]]


########## Get Eigenstrat PCs ##########

# Lese die Eigenvalues von der ersten Zeile (Kommentarzeile)
eigenvalues <- as.vector(read.table(filename,nrows=1,comment.char = "", row.names=1))
pcs <- read.table(filename,row.names=1,col.names=c("rowname","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","pheno"))


########## Sample annotation ##########

# Obtain the sample annotation
sample_annotation <- read.table(annotation_file, header=FALSE)

# Bring samples in genotype and in sample annotation in the same order
order <- sample_annotation[,1]
pcs<- pcs[order,]


########## Genotype principal components ##########

# Plot the PCs
colors <- rainbow(14) #c("#1F78B4","#A6CEE3","#FDBF6F","#E31A1C")
var_pc1 <- eigenvalues[1]
var_pc2 <- eigenvalues[2]
var_pc3 <- eigenvalues[3]
var_pc4 <- eigenvalues[4]
group <- sample_annotation[,2]
pdf(paste(path,"pca_1vs2.pdf",sep=""))
xyplot(PC2~PC1,groups=group, data=as.data.frame(pcs),pch=16,cex=1,col=colors,main=draw.key(key=list(rect=list(col=colors),text=list(levels(group)))),xlab=paste("PC1 (",var_pc1,")",sep=""),ylab=paste("PC2 (",var_pc2,")",sep=""))
dev.off()
pdf(paste(path,"pca_2vs3.pdf",sep=""))
xyplot(PC3~PC2,groups=group, data=as.data.frame(pcs),pch=16,cex=1,col=colors,main=draw.key(key=list(rect=list(col=colors),text=list(levels(group)))),xlab=paste("PC2 (",var_pc2,")",sep=""),ylab=paste("PC3 (",var_pc3,")",sep=""))
dev.off()
pdf(paste(path,"pca_1vs3.pdf",sep=""))
xyplot(PC3~PC1,groups=group, data=as.data.frame(pcs),pch=16,cex=1,col=colors,main=draw.key(key=list(rect=list(col=colors),text=list(levels(group)))),xlab=paste("PC1 (",var_pc1,")",sep=""),ylab=paste("PC3 (",var_pc3,")",sep=""))
dev.off()
pdf(paste(path,"pca_1vs4.pdf",sep=""))
xyplot(PC4~PC1,groups=group, data=as.data.frame(pcs),pch=16,cex=1,col=colors,main=draw.key(key=list(rect=list(col=colors),text=list(levels(group)))),xlab=paste("PC1 (",var_pc1,")",sep=""),ylab=paste("PC4 (",var_pc4,")",sep=""))
dev.off()
pdf(paste(path,"pca_2vs4.pdf",sep=""))
xyplot(PC4~PC2,groups=group, data=as.data.frame(pcs),pch=16,cex=1,col=colors,main=draw.key(key=list(rect=list(col=colors),text=list(levels(group)))),xlab=paste("PC2 (",var_pc2,")",sep=""),ylab=paste("PC4 (",var_pc4,")",sep=""))
dev.off()
pdf(paste(path,"pca_3vs4.pdf",sep=""))
xyplot(PC4~PC3,groups=group, data=as.data.frame(pcs),pch=16,cex=1,col=colors,main=draw.key(key=list(rect=list(col=colors),text=list(levels(group)))),xlab=paste("PC3 (",var_pc3,")",sep=""),ylab=paste("PC4 (",var_pc4,")",sep=""))
dev.off()
pdf(paste(path,"scree_plot.pdf",sep=""))
plot(1:10,eigenvalues,type="b",main="Scree plot genotype PCA")
dev.off()

