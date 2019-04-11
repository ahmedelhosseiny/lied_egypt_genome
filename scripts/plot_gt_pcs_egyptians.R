# Writes the first 10 genotype principal components to file for later inclusion 
# into the linear model. Computes various visualizations, e.g. PCs, heatmaps and
# scree plot.

library("latticeExtra")
library("scatterplot3d")

# Get the name of the input file
filename <- snakemake@input[[1]]
annotation_file <- snakemake@input[[2]]
path <- snakemake@params[[1]]
fname_pca_1vs2 <- snakemake@output[[1]]
fname_pca_1vs3 <- snakemake@output[[2]]
fname_pca_1vs4 <- snakemake@output[[3]]
fname_pca_2vs3 <- snakemake@output[[4]]
fname_pca_2vs4 <- snakemake@output[[5]]
fname_pca_3vs4 <- snakemake@output[[6]]
fname_scree <- snakemake@output[[7]]
fname_pca_3d <- snakemake@output[[8]]

########## Get Eigenstrat PCs ##########

# Lese die Eigenvalues von der ersten Zeile (Kommentarzeile)
eigenvalues <- as.vector(read.table(filename,nrows=1,comment.char = "", row.names=1))
pcs <- read.table(filename,row.names=1,col.names=c("rowname","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","pheno"))


########## Sample annotation ##########

# Obtain the sample annotation
sample_annotation <- read.table(annotation_file, header=FALSE)
egyptian_ids <- c("EGYPTREF","LU18","LU19","LU2","LU22","LU23","LU9","PD114","PD115","PD82","EGAN00001101667","EGAN00001101668","EGAN00001101669","EGAN00001101670","EGAN00001101671","EGAN00001101672","EGAN00001101676","EGAN00001101677","EGAN00001101678","EGAN00001101679","EGAN00001101680","EGAN00001101681","EGAN00001101682","EGAN00001101687","EGAN00001101688","EGAN00001101689","EGAN00001101690","EGAN00001101692","EGAN00001101694","EGAN00001101699","EGAN00001101700","EGAN00001101702","EGAN00001101705","EGAN00001101706","EGAN00001101711","EGAN00001101712","EGAN00001101713","EGAN00001101716","EGAN00001101717","EGAN00001101718","EGAN00001101719","EGAN00001101723","EGAN00001101724","EGAN00001101725","EGAN00001101732","EGAN00001101734","EGAN00001101735","EGAN00001101736","EGAN00001101737","EGAN00001101739","EGAN00001101742","EGAN00001101744","EGAN00001101748","EGAN00001101749","EGAN00001101750","EGAN00001101751","EGAN00001101752","EGAN00001101753","EGAN00001101754","EGAN00001101755","EGAN00001101756","EGAN00001101758","EGAN00001101759","EGAN00001101761","EGAN00001101767","EGAN00001101768","EGAN00001101769","EGAN00001101771","EGAN00001101772","EGAN00001101774","EGAN00001101776","EGAN00001101780","EGAN00001101781","EGAN00001101782","EGAN00001101783","EGAN00001101784","EGAN00001101786","EGAN00001101787","EGAN00001101788","EGAN00001101791","EGAN00001101792","EGAN00001101793","EGAN00001101794","EGAN00001101796","EGAN00001101797","EGAN00001101798","EGAN00001101799","EGAN00001101801","EGAN00001101802","EGAN00001101803","EGAN00001101804","EGAN00001101807","EGAN00001101808","EGAN00001101809","EGAN00001101813","EGAN00001101814","EGAN00001101816","EGAN00001101819","EGAN00001101820","EGAN00001101823","EGAN00001101824","EGAN00001101825","EGAN00001101827","EGAN00001101829","EGAN00001101830","EGAN00001101831","EGAN00001101835","EGAN00001101839","EGAN00001101840","EGAN00001101841")
sample_annotation <- sample_annotation[sample_annotation[,1] %in% egyptian_ids,]

# Get the PC values of the reference Egyptian
pc_egyptref <- pcs["EGYPTREF",]

# Get the PC values of "our" 6 upper Egyptians 
pc_upper_egyptians <- pcs[c("LU18","LU19","LU2","LU22","LU23","LU9","PD114","PD115","PD82"),]

# Get the PC values of "our" 3 lower Egyptians
pc_lower_egyptians <- pcs[c("PD114","PD115","PD82"),]

# Bring samples in genotype and in sample annotation in the same order
order <- sample_annotation[,1]
pcs <- pcs[as.character(order),]
# Make sure that rownames and sample annotation match
all(as.character(sample_annotation[,1]) == row.names(pcs))


########## Genotype principal components ##########

# Plot the PCs
pc_col <- rainbow(2)
pc_col[1] <- rainbow(14)[1]
pc_col[2] <- rainbow(14)[2]
pc_col[3] <- rainbow(14)[3]
eig_pc1 <- eigenvalues[1]
eig_pc2 <- eigenvalues[2]
eig_pc3 <- eigenvalues[3]
eig_pc4 <- eigenvalues[4]
group <- as.factor(as.vector(sample_annotation[,2]))

# PC1/PC2
pdf(fname_pca_1vs2)
# pch=16
fig <- xyplot(PC2~PC1,groups=group, data=as.data.frame(pcs),pch=16,cex=1,col=pc_col,main=draw.key(key=list(rect=list(col=pc_col),text=list(levels(group)))),xlab="PC1",ylab="PC2")#,panel=function(x,y,labels,...){panel.xyplot(pcs$PC1,pcs$PC2,...);ltext(x=pcs$PC1, y=pcs$PC2, labels=sample_annotation[,1], pos=1, offset=1, cex=0.4);})
## Overlay Egyptref samples
fig <- fig + layer(panel.points(pc_upper_egyptians$PC1, pc_upper_egyptians$PC2, pch=16, cex=1, col=pc_col[1]))
fig <- fig + layer(panel.points(pc_lower_egyptians$PC1, pc_lower_egyptians$PC2, pch=16, cex=1, col=pc_col[3]))
fig <- fig + layer(panel.points(pc_egyptref$PC1, pc_egyptref$PC2, pch=21, cex=1, col='black', fill=pc_col[1]))
fig
dev.off()

# PC1/PC3
pdf(fname_pca_1vs3)
fig <- xyplot(PC3~PC1,groups=group, data=as.data.frame(pcs),pch=16,cex=1,col=pc_col,main=draw.key(key=list(rect=list(col=pc_col),text=list(levels(group)))),xlab="PC1",ylab="PC3")#,panel=function(x,y,labels,...){panel.xyplot(pcs$PC1,pcs$PC3,...);ltext(x=pcs$PC1, y=pcs$PC3, labels=sample_annotation[,1], pos=1, offset=1, cex=0.4);})
## Overlay Egyptian samples
fig <- fig + layer(panel.points(pc_upper_egyptians$PC1, pc_upper_egyptians$PC3, pch=16, cex=1, col=pc_col[1]))
fig <- fig + layer(panel.points(pc_lower_egyptians$PC1, pc_lower_egyptians$PC3, pch=16, cex=1, col=pc_col[3]))
fig <- fig + layer(panel.points(pc_egyptref$PC1, pc_egyptref$PC3, pch=21, cex=1, col='black', fill=pc_col[1]))
fig
dev.off()

# PC1/PC4
pdf(fname_pca_1vs4)
fig <- xyplot(PC4~PC1,groups=group, data=as.data.frame(pcs),pch=16,cex=1,col=pc_col,main=draw.key(key=list(rect=list(col=pc_col),text=list(levels(group)))),xlab="PC1",ylab="PC4")#,panel=function(x,y,labels,...){panel.xyplot(pcs$PC1,pcs$PC4,...);ltext(x=pcs$PC1, y=pcs$PC4, labels=sample_annotation[,1], pos=1, offset=1, cex=0.4);})
## Overlay Egyptian samples
fig <- fig + layer(panel.points(pc_upper_egyptians$PC1, pc_upper_egyptians$PC4, pch=16, cex=1, col=pc_col[1]))
fig <- fig + layer(panel.points(pc_lower_egyptians$PC1, pc_lower_egyptians$PC4, pch=16, cex=1, col=pc_col[3]))
fig <- fig + layer(panel.points(pc_egyptref$PC1, pc_egyptref$PC4, pch=21, cex=1, col='black', fill=pc_col[1]))
fig
dev.off()

# PC2/PC3
pdf(fname_pca_2vs3)
fig <- xyplot(PC3~PC2,groups=group, data=as.data.frame(pcs),pch=16,cex=1,col=pc_col,main=draw.key(key=list(rect=list(col=pc_col),text=list(levels(group)))),xlab="PC2",ylab="PC3")#,panel=function(x,y,labels,...){panel.xyplot(pcs$PC2,pcs$PC3,...);ltext(x=pcs$PC2, y=pcs$PC3, labels=sample_annotation[,1], pos=1, offset=1, cex=0.4);})
## Overlay Egyptian samples
fig <- fig + layer(panel.points(pc_upper_egyptians$PC2, pc_upper_egyptians$PC3, pch=16, cex=1, col=pc_col[1]))
fig <- fig + layer(panel.points(pc_lower_egyptians$PC2, pc_lower_egyptians$PC3, pch=16, cex=1, col=pc_col[3]))
fig <- fig + layer(panel.points(pc_egyptref$PC2, pc_egyptref$PC3, pch=21, cex=1, col='black', fill=pc_col[1]))
fig
dev.off()

# PC2/PC4
pdf(fname_pca_2vs4)
fig <- xyplot(PC4~PC2,groups=group, data=as.data.frame(pcs),pch=16,cex=1,col=pc_col,main=draw.key(key=list(rect=list(col=pc_col),text=list(levels(group)))),xlab="PC2",ylab="PC4")#,panel=function(x,y,labels,...){panel.xyplot(pcs$PC2,pcs$PC4,...);ltext(x=pcs$PC2, y=pcs$PC4, labels=sample_annotation[,1], pos=1, offset=1, cex=0.4);})
## Overlay Egyptian samples
fig <- fig + layer(panel.points(pc_upper_egyptians$PC2, pc_upper_egyptians$PC4, pch=16, cex=1, col=pc_col[1]))
fig <- fig + layer(panel.points(pc_lower_egyptians$PC2, pc_lower_egyptians$PC4, pch=16, cex=1, col=pc_col[3]))
fig <- fig + layer(panel.points(pc_egyptref$PC2, pc_egyptref$PC4, pch=21, cex=1, col='black', fill=pc_col[1]))
fig
dev.off()

# PC3/PC4
pdf(fname_pca_3vs4)
fig <- xyplot(PC4~PC3,groups=group, data=as.data.frame(pcs),pch=16,cex=1,col=pc_col,main=draw.key(key=list(rect=list(col=pc_col),text=list(levels(group)))),xlab="PC3",ylab="PC4")#,panel=function(x,y,labels,...){panel.xyplot(pcs$PC3,pcs$PC4,...);ltext(x=pcs$PC3, y=pcs$PC4, labels=sample_annotation[,1], pos=1, offset=1, cex=0.4);})
## Overlay Egyptian samples
fig <- fig + layer(panel.points(pc_upper_egyptians$PC3, pc_upper_egyptians$PC4, pch=16, cex=1, col=pc_col[1]))
fig <- fig + layer(panel.points(pc_lower_egyptians$PC3, pc_lower_egyptians$PC4, pch=16, cex=1, col=pc_col[3]))
fig <- fig + layer(panel.points(pc_egyptref$PC3, pc_egyptref$PC4, pch=21, cex=1, col='black', fill=pc_col[1]))
fig
dev.off()

# Scree plot
pdf(fname_scree)
plot(1:10,eigenvalues,type="b",main="Scree plot genotype PCA")
dev.off()

# 3D Plot

# Make a vector of colors for the dots (samples) to be plotted
num_samples <- length(group)
sample_names <- sample_annotation[,1]
cols <- rep("",num_samples)
for (i in 1:num_samples){
  sample_class <- group[i]
  cols[i] <- pc_col[match(sample_class,as.vector(levels(group)))]
}

pdf(fname_pca_3d,useDingbats=FALSE)
s3d <- scatterplot3d(pcs$PC1, pcs$PC2, pcs$PC3,
              color='black', pch=21, box=FALSE, tick.marks=FALSE, bg=cols, angle=20,
              xlab="PC1",
              ylab="PC2",
              zlab="PC3")
s3d <- scatterplot3d(pc_upper_egyptians$PC1, pc_upper_egyptians$PC2, pc_upper_egyptians$PC3,
              color='black', pch=21, box=FALSE, tick.marks=FALSE, bg=pc_col[1], angle=20,
              xlab="PC1",
              ylab="PC2",
              zlab="PC3")
s3d <- scatterplot3d(pc_lower_egyptians$PC1, pc_lower_egyptians$PC2, pc_lower_egyptians$PC3,
              color='black', pch=21, box=FALSE, tick.marks=FALSE, bg=pc_col[3], angle=20,
              xlab="PC1",
              ylab="PC2",
              zlab="PC3")
s3d <- scatterplot3d(pc_egyptref$PC1, pc_egyptref$PC2, pc_egyptref$PC3,
              color='gray', pch=21, box=FALSE, tick.marks=FALSE, bg=pc_col[3], angle=20,
              xlab="PC1",
              ylab="PC2",
              zlab="PC3")
s3d.coords <- s3d$xyz.convert(pcs$PC1, pcs$PC2, pcs$PC3)
#text(s3d.coords$x, s3d.coords$y, labels=sample_names, pos=4, cex=.5)
legend("topleft",inset=.05,bty="n",title="",levels(group),fill=pc_col,cex=0.8)
dev.off()



