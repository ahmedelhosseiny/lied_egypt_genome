# Making boxplots for snp stat numbers
library("RColorBrewer")

# Getting input and output filenames
# ["insertions","deletions","inversions","duplications"],
# filter=["all","pass"])
fname_ins_all <- snakemake@input[[1]]
fname_ins_pass <- snakemake@input[[2]]
fname_del_all <- snakemake@input[[3]]
fname_del_pass <- snakemake@input[[4]]
fname_inv_all <- snakemake@input[[5]]
fname_inv_pass <- snakemake@input[[6]]
fname_dup_all <- snakemake@input[[7]]
fname_dup_pass <- snakemake@input[[8]]
fname_trans_all <- snakemake@input[[9]]
fname_trans_pass <- snakemake@input[[10]]
fname_out_hist_all <- snakemake@output[[1]]
fname_out_hist_pass <- snakemake@output[[2]]
fname_out_svs_all <- snakemake@output[[3]]
fname_out_svs_pass <- snakemake@output[[4]]

# Plotting
pdf(fname_out_hist_all)
par(mfrow=c(2,2))

hist_cols <- brewer.pal(n = 8, name = 'Set2')[3:6]

del_all <- read.table(fname_del_all, header=FALSE)
hist(log10(del_all$V1),xlab="log10(Length)",ylab="Number",main="Deletions",breaks=30,col=hist_cols[1])
num_del_all <- length(del_all$V1)

inv_all <- read.table(fname_inv_all, header=FALSE)
hist(log10(inv_all$V1),xlab="log10(Length)",ylab="Number",main="Inversions",breaks=30,col=hist_cols[2])
num_inv_all <- length(inv_all$V1)

dup_all <- read.table(fname_dup_all, header=FALSE)
hist(log10(dup_all$V1),xlab="log10(Length)",ylab="Number",main="Duplications",breaks=30,col=hist_cols[3])
num_dup_all <- length(dup_all$V1)

ins_all <- read.table(fname_ins_all, header=FALSE)
hist(ins_all$V1,xlab="Length",ylab="Number",main="Insertions",breaks=30,col=hist_cols[4])
num_ins_all <- length(ins_all$V1)

dev.off()

pdf(fname_out_hist_pass)
par(mfrow=c(2,2))

del_pass <- read.table(fname_del_pass, header=FALSE)
hist(log10(del_pass$V1),xlab="log10(Length)",ylab="Number",main="Deletions",breaks=30,col=hist_cols[1])
num_del_pass <- length(del_pass$V1)

inv_pass <- read.table(fname_inv_pass, header=FALSE)
hist(log10(inv_pass$V1),xlab="log10(Length)",ylab="Number",main="Inversions",breaks=30,col=hist_cols[2])
num_inv_pass <- length(inv_pass$V1)

dup_pass <- read.table(fname_dup_pass, header=FALSE)
hist(log10(dup_pass$V1),xlab="log10(Length)",ylab="Number",main="Duplications",breaks=30,col=hist_cols[3])
num_dup_pass <- length(dup_pass$V1)

ins_pass <- read.table(fname_ins_pass, header=FALSE)
hist(ins_pass$V1,xlab="Length",ylab="Number",main="Insertions",breaks=30,col=hist_cols[4])
num_ins_pass <- length(ins_pass$V1)

dev.off()

pdf(fname_out_svs_all,width=5,height=8)

num_trans_all <- read.table(fname_trans_all)$V1
num_trans_pass <- read.table(fname_trans_pass)$V1

hist_cols <- brewer.pal(n = 8, name = 'Set2')[3:7]

sv_nums <- c(num_del_all,num_inv_all,num_dup_all,num_ins_all,num_trans_all)
sv_names <- c("Deletions","Inversions","Duplications","Insertions","Translocations")
barplot(sv_nums,names=sv_names,ylab="Number",col=hist_cols,las=2,cex.names=0.77)

dev.off()


pdf(fname_out_svs_pass,width=5,height=8)

sv_nums <- c(num_del_pass,num_inv_pass,num_dup_pass,num_ins_pass,num_trans_pass)
sv_names <- c("Deletions","Inversions","Duplications","Insertions","Translocations")
barplot(sv_nums,names=sv_names,ylab="Number",col=hist_cols,las=2,cex.names=0.77)

dev.off()