
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
setwd("C:\\Users\\Bart\\Desktop\\Leren Programmeren\\Centrotype")

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

ok_function <- function(m){
  any(apply(apply(-log10(m) > 3.5,2,as.numeric),2,sum) >= 3) 
}

relatiefCentrotype <- function(mm){
	matrixCor<- matrix(apply(-log10(as.matrix(unlist(mm))) > 3.5,2, as.numeric), ncol=69)
	maxmarker <- which.max(apply(matrixCor, 2,  sum))
	tempprobes <- which(matrixCor[,maxmarker] == 1)
	output <- c("Marker:", maxmarker, "probes:", tempprobes)
}

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
      heatmap(apply(-log10(mm) > 3.5,2,as.numeric),col=c("white","black"), scale="none")
      #dev.off()
      cat(name, paste(" Took:", (proc.time()-TijdA)[3], "seconds", sep=" "), "\n")
      res[[cnt]]$Name <- name
      res[[cnt]]$Matrix <- mm
      res[[cnt]]$Ok <- ok_function(mm)
	  res[[cnt]]$Centrotype <- relatiefCentrotype(mm)
    }else{
	  res[[cnt]]$Name <- name
      res[[cnt]]$Matrix <- 0
	  res[[cnt]]$Ok <- FALSE
	  res[[cnt]]$Centrotype <- FALSE
	}
    cnt <- cnt + 1
	}
	invisible(res)
}

# Write the file 
save(myList, file="test1.bin") 

# Reload the data, under the same name, 'myList' 
load(file="test1.bin") 

#om een lijst uit te lezen op een naam:
# unlist(lapply(<lijst>,"[","<name>"))
lijst
for(ele in 1:69){
  if(ss[ele] == 1){ cat(ele, "juist", "\n")}
}

#Om je 0 en 1 waarden uit te lezen, dit kan ook in een matrix gezet worden met ncol= 69 en nrow is de hoeveelheid probes uit dd
 cc<-apply(-log10(as.matrix(unlist(dd), ncol=96, nrow= 20)) > 3.5,2, as.numeric)
 
 relatiefCentrotype <- function(pvalmatrix, x){
	matrixCor<- matrix(apply(-log10(as.matrix(unlist(pvalmatrix[[1]]$mmatrix))) > 3.5,2, as.numeric), ncol=69)
	maxmarker <- which.max(apply(matrixCor, 2,  sum))
	tempprobes <- which(matrixCor[,maxmarker] == 1)
	output <- c("Marker:", maxmarker, "probes:", tempprobes)
}
 
#Manier van lijst maken die nodig is om de data van een Centrotype in te stoppen
myrange <- 1:10
namen <- NULL
cnt <- 1
Basiclist <- vector("list", length(myrange))
for(x in myrange){
	output1 <- x^3
	output2 <- x+4
	output3 <- x+5
	output4 <- x+6
	output5 <- matrix((1:20)+x, ncol=5, nrow=4)
	Basiclist[[cnt]]$Macht <- output1
	Basiclist[[cnt]]$Lineair[1] <- output2
	Basiclist[[cnt]]$Lineair[2] <- output3
	Basiclist[[cnt]]$Lineair[3] <- output4
	Basiclist[[cnt]]$Matrix <- output5
	namen <- c(namen, paste("Gen", x, sep=""))
	cat(namen, "\n")
	cnt <- cnt+1
}
names(Basiclist) <- namen

#Gevonden Centrotypes
#> Centrotypes1...200 <- which(unlist(lapply(SUPERLIJST,"[","ok")))
#> Centrotypes1...200
# ok  ok  ok  ok  ok  ok  ok  ok  ok  ok  ok  ok  ok 
# 17  31  86 104 108 136 137 138 145 176 179 185 191 

print.hotjeknor <- function(x){
  cat("Results for ",length(x),"Genes")
  cat("Has ",length(which(unlist(lapply(x,"[","ok")))),"c")
}
