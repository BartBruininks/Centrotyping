
getAdjustedProbes <- function(x){
 dd <- read.table(paste("chr1/",dir("chr1/")[x],sep=""))
 rr <- NULL
 for(y in 1:nrow(dd)){
    pval <- anova(lm(t(dd[y,17:196]) ~ samplev2[,3]))[[5]][1]
	if(-log10(pval) > 3){
	  rr <- rbind(rr,c(dd[y,1:16],residuals(lm(t(dd[y,17:196]) ~ samplev2[,3]))))
	}else{
	  rr <- rbind(rr,c(dd[y,1:16],dd[y,17:196]))
	}
 }
 invisible(rr)
}

 apply(cor(t(apply(aa[,17:196],2,as.numeric))),2,function(x){rownames(bb)[which(x> 0.6)]})


boxplot(all)

## verwijder de verkeerde dag uit de data van het genotype en probematrix

newgenomatrix <- NULL
for(x in samplev2[9:172,2]){
  newgenomatrix <- rbind(newgenomatrix,genotypes[which(paste("RIL",genotypes[,1],sep="") == x),8:ncol(genotypes)])
}


for(x in 1:100){
	rawdata <- read.table(paste("chr1/",dir("chr1/")[x],sep=""))
	#rawdata <- getAdjustedProbes(x)
	nonfounders <- grep("RIL",samplev2[,2])
	probes <- rawdata[,17:196][,nonfounders]

	mm <- NULL
	for(t in 1:nrow(probes)){
		pvals <- NULL
		for(x in 1:ncol(newgenomatrix)){
		 pvals <- c(pvals, anova(lm(unlist(probes[t,]) ~ newgenomatrix[,x]))[[5]][1])
		}
		mm <- rbind(mm, pvals)
	}
	heatmap(apply(-log10(mm) > 3.5,2,as.numeric),col=c("white","black"))
}