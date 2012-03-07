
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

## verwijder de verkeerde dag uit de data van het genotype en de probematrix

rawfaulty <- c(128, 130, 131, 132, 133, 134, 135, 137, 138, 139, 141, 142, 143, 144, 145)
faulty <- strsplit("V128, V130, V131, V132, V133, V134, V135, V137, V138, V139, V141, V142, V143, V144, V145",", ")
pool <- (9:172)
samplev2 <- read.table("samplev2.txt")
genotypes <- read.table("genotypes.txt")

samples <- samplev2[grep("RIL",samplev2[,2]), ]
torem <- which(as.character(samples[,1]) %in% unlist(faulty))
samples <- samples[-torem,]

newgenomatrix <- NULL
for(x in samples[,2]){
  genrij <- which(paste("RIL",genotypes[,1],sep="") == x)
  cat(x, genrij,"\n")
  newgenomatrix <- rbind(newgenomatrix,genotypes[genrij,5:ncol(genotypes)])
}
rownames(newgenomatrix) <-  samples[,1]

Centrotype <- function(myrange = c(1 : length(dir("chr1/")))){
	res <- vector("list",length(myrange))
  cnt <- 1
  for(x in myrange){
		TijdA <-proc.time()
		rawdata <- read.table(paste("chr1/",dir("chr1/")[x],sep=""))
		#rawdata <- getAdjustedProbes(x)
		probes <- rawdata[,  rownames(newgenomatrix)]
    if(nrow(probes) >= 2){
      mm <- NULL
      for(p in 1:nrow(probes)){
        pvals <- NULL
        for(ele in 1:ncol(newgenomatrix)){
         pvals <- c(pvals, anova(lm(unlist(probes[p,]) ~ newgenomatrix[,ele]))[[5]][1])
        }
        mm <- rbind(mm, pvals)
      }
      name = strsplit(dir("chr1/")[x], "[.]")[[1]][1]
      #jpeg(filename = paste("regulatie", name, ".jpeg", sep=""), width = 1000, height = 1000, bg = "white", units = "px")
      heatmap(apply(-log10(mm) > 3.5,2,as.numeric),col=c("white","black"))
      #dev.off()
      cat(name, paste(" Took:", (proc.time()-TijdA)[3], "seconds", sep=" "), "\n")
      res[[cnt]]$name <- x
      res[[cnt]]$mmatrix <- mm
      res[[cnt]]$ok <- ok_function(mm)
    }
    cnt <- cnt + 1
	}
	invisible(res)
}

ok_function <- function(m){
  any(apply(apply(-log10(m) > 3.5,2,as.numeric),2,sum) >= 3) 
}

 unlist(lapply(aa,"[","name"))

# Krijg deze manier van een kolom uitrekenen en in een lijst zetten in de functie van Centrotype #
rij <- 1:20 
lijst <- NULL
matr <- matrix(rij, ncol=10, nrow=2)
for(x in 1:ncol(matr)){
  if(sum(matr[,x]) >= 15){
    lijst <- c(lijst, (matr[,x]))
    lijst
   }
}  




