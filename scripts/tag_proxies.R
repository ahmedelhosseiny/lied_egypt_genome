# Visualizations with respect to proxy SNPs of tag SNPs

# Getting input and output filenames
fname_in <- snakemake@input[[1]]
fname_out_proxynumber <- snakemake@output[[1]]
fname_out_proxysharing <- snakemake@output[[2]]
fname_out_af_difference <- snakemake@output[[3]]
fname_out_eurvsegp <- snakemake@output[[4]]
fname_out_eur_af_vs_proxynum <- snakemake@output[[5]]
fname_out_egp_af_vs_proxynum <- snakemake@output[[6]]

data <- read.table(fname_in,header=T,colClasses=c("numeric","numeric","factor","factor","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"))

# Data consistency check for European proxys, this should return nothing
which(!(data$PROXY_SHARED+data$PROXY_EUR_ONLY == data$EUR_N_PROXY))

# Same data consistency check, should also return nothing
which(!(data$PROXY_SHARED+data$PROXY_EGP_ONLY == data$EGP_N_PROXY))

# Plotting boxplot of number of proxies
c1 <- rainbow(10)[c(7,1)]
c2 <- rainbow(10, alpha=0.2)[c(7,1)]
c3 <- rainbow(10, v=0.7)[c(7,1)]
pdf(fname_out_proxynumber)
boxplot(na.omit(data[,9:10]+1),ylab="number of proxies",log="y",names=c("European","Egyptian"),outcex=0.5,col=c2,medcol=c3,whiskcol=c1,staplecol=c3,boxcol=c3,outcol=c2,pch=19,boxwex=0.5)
dev.off()

# Plotting boxplot of sharing of proxies
# In order not to divide by zero, we add +1, which means we also count the tag SNP here, which by definition is always shared
percent_shared_of_egyptian <- 100*(data$PROXY_SHARED+1)/(data$PROXY_EGP_ONLY+(data$PROXY_SHARED+1))
percent_shared_of_european <- 100*(data$PROXY_SHARED+1)/(data$PROXY_EUR_ONLY+(data$PROXY_SHARED+1))
percent_egp_only <- 100*data$PROXY_EGP_ONLY/(data$PROXY_EUR_ONLY+data$PROXY_EGP_ONLY+(data$PROXY_SHARED+1))
percent_eur_only <- 100*data$PROXY_EUR_ONLY/(data$PROXY_EUR_ONLY+data$PROXY_EGP_ONLY+(data$PROXY_SHARED+1))

c1 <- rainbow(10)[c(7,7,1,1)]
c2 <- rainbow(10, alpha=0.2)[c(7,7,1,1)]
c3 <- rainbow(10, v=0.7)[c(7,7,1,1)]
pdf(fname_out_proxysharing)
boxplot(cbind(percent_shared_of_european,percent_eur_only,percent_shared_of_egyptian,percent_egp_only),ylab="%",names=c("Shared European","European-only","Shared Egyptian","Egyptian-only"),cex.axis=0.9,outcex=0.5,col=c2,medcol=c3,whiskcol=c1,staplecol=c3,boxcol=c3,outcol=c2,pch=19,boxwex=0.5)
dev.off()

# Histogram of af differences
pdf(fname_out_af_difference)
hist(data$EUR_ALT_AF-data$EGP_ALT_AF,cex=0.3,breaks=100,xlab="European AF - Egyptian AF",col="gray",border="gray",main="")
abline(v=0,lty=2,lwd=3,col="black")
dev.off()

# European versus Egyptian proxy numbers
pdf(fname_out_eurvsegp)
plot(data$EUR_N_PROXY,data$EGP_N_PROXY,cex=0.3,xlab="European proxy number",ylab="Egyptian proxy number",pch=19)
dev.off()

# Plotting AF versus proxy number for European data
pdf(fname_out_eur_af_vs_proxynum)
plot(data$EUR_ALT_AF,data$EUR_N_PROXY,cex=0.3,xlab="European AF",ylab="proxy number")
dev.off()

# Plotting AF versus proxy number for Egyptian data
pdf(fname_out_egp_af_vs_proxynum)
plot(data$EGP_ALT_AF,data$EGP_N_PROXY,cex=0.3,xlab="Egyptian AF",ylab="proxy number")
dev.off()
