#(c) Bart Bruininks written in the [R] language
# 12-03-2012
# Centrotyping Probes

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

#apply(cor(t(apply(aa[,17:196],2,as.numeric))),2,function(x){rownames(bb)[which(x> 0.6)]})

 
#boxplot(all)
###############################################################################
## verwijder de verkeerde dag uit de data van het genotype en de probematrix ##
##                   CENTROTYPING FUNCTION WORKING!!!!                       ##
###############################################################################
setwd("C:\\Users\\Bart\\Desktop\\Leren Programmeren\\Centrotype")
#setwd("C:\\Users\\Students\\Desktop\\Hoofdmap")								# pathway voor de server!!!

rawfaulty <- c(128, 130, 131, 132, 133, 134, 135, 137, 138, 139, 141, 142, 143, 144, 145)
faulty <- unlist(strsplit("V128, V130, V131, V132, V133, V134, V135, V137, V138, V139, V141, V142, V143, V144, V145",", "))
ind <- paste("V", 1:180, sep="")
goodind <- ind[-rawfaulty]


pool <- (9:172)
samplev2 <- read.table("samplev2.txt")
genotypes <- read.table("genotypes.txt")


samples <- samplev2[grep("RIL",samplev2[,2]), ]
torem <- which(as.character(samples[,1]) %in% faulty)
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


Centrotype <- function(myrange = c(1 : length(dir("chr1/"))), enablematrix = FALSE){
  res <- vector("list",length(myrange))
  cnt <- 1
  for(x in myrange){
		TijdA <-proc.time()
		rawdata <- read.table(paste("chr1/",dir("chr1/")[x],sep=""))								#Waarden met omgevingsverschil
		#rawdata <- getAdjustedProbes(x)															#Waarden zonder omgevingsverschil
		probes <- rawdata[,  rownames(newgenomatrix)]
    if(nrow(probes) >= 2){
      #mm <- NULL
      #for(p in 1:nrow(probes)){
      #  pvals <- NULL
      #  for(ele in 1:ncol(newgenomatrix)){
      #   pvals <- c(pvals, anova(lm(unlist(probes[p,]) ~ newgenomatrix[,ele]))[[5]][1])
      #  }
      #  mm <- rbind(mm, pvals)
      #}
	  probes <- as.matrix(t(apply(probes,2,as.numeric)))
	  mm <- apply(newgenomatrix,2,function(m){
	    unlist(lapply(summary(aov(probes ~ as.factor(m))),"[",1,5),use.names=TRUE)
		}
	  )
      name = strsplit(dir("chr1/")[x], "[.]")[[1]][1]
      #jpeg(filename = paste("regulatie", name, ".jpeg", sep=""), width = 1000, height = 1000, bg = "white", units = "px")
      heatmap(apply(-log10(mm) > 3.5,2,as.numeric),col=c("white","black"), scale="none")
      #dev.off()
      cat(name, paste(" Took:", (proc.time()-TijdA)[3], "seconds", sep=" "), "\n")
		if(enablematrix == TRUE){
			res[[cnt]]$Matrix <- mm
			res[[cnt]]$Name <- name
			res[[cnt]]$Ok <- ok_function(mm)
			if(res[[cnt]]$Ok == FALSE){
			res[[cnt]]$Centrotype <- FALSE
			}else{
			res[[cnt]]$Centrotype <-relatiefCentrotype(mm)}
		}else{
			res[[cnt]]$Name <- name
			res[[cnt]]$Ok <- ok_function(mm)
			if(res[[cnt]]$Ok == FALSE){
			res[[cnt]]$Centrotype <- FALSE
			}else{
			res[[cnt]]$Centrotype <-relatiefCentrotype(mm)}
		}
    }else{
	  res[[cnt]]$Matrix <- "EMPTY"
	  res[[cnt]]$Name <- name
      res[[cnt]]$Ok <- FALSE
	  res[[cnt]]$Centrotype <- FALSE
	}
	cat(name, res[[cnt]]$Centrotype, "\n",sep="\t", file = "Centrotypes.txt", append = TRUE)
    cnt <- cnt + 1
	}
	invisible(res)
}

readCentrotypes <- function(x = "Centrotypes.txt"){
  filein <- file(x, "rt") #rt zorgt ervoor dat ik mijn file open laat staan terwijl ik hem lees
  nline <- readLines(filein, n=1,)
  s <- proc.time()
  splitit <- NULL
  namez <- NULL
  cnt <- 1
  while(length(nline)){
    elementen <- strsplit(nline, "\t")
    splitit <- c(splitit, elementen)
	splitit[[cnt]] <- splitit[[cnt]][-c(1)]
	namez <- c(namez, elementen[[1]][1])
	nline <- readLines(filein, n=1)
	cnt <- cnt+1
  }
  e <- proc.time()
  cat(e[3]-s[3], "seconds\n")
  close(filein)
  names(splitit) <- paste("", namez, sep="")
  invisible(splitit)
}




#write(x, file = "data",
#      ncolumns = if(is.character(x)) 1 else 5,
#      append = FALSE, sep = " ")


# Write the file 
#save(myList, file="test1.bin") 

# Reload the data, under the same name, 'myList' 
#load(file="test1.bin") 

#om een lijst uit te lezen op een naam:
#
#lijst
#for(ele in 1:69){
#  if(ss[ele] == 1){ cat(ele, "juist", "\n")}
#}

#Om je 0 en 1 waarden uit te lezen, dit kan ook in een matrix gezet worden met ncol= 69 en nrow is de hoeveelheid probes uit dd
# cc<-apply(-log10(as.matrix(unlist(dd), ncol=96, nrow= 20)) > 3.5,2, as.numeric)
 
# relatiefCentrotype <- function(pvalmatrix, x){
#	matrixCor<- matrix(apply(-log10(as.matrix(unlist(pvalmatrix[[1]]$mmatrix))) > 3.5,2, as.numeric), ncol=69)
#	maxmarker <- which.max(apply(matrixCor, 2,  sum))
#	tempprobes <- which(matrixCor[,maxmarker] == 1)
#	output <- c("Marker:", maxmarker, "probes:", tempprobes)
#}
 
#Manier van lijst maken die nodig is om de data van een Centrotype in te stoppen
#myrange <- 1:10
#namen <- NULL
#cnt <- 1
#Basiclist <- vector("list", length(myrange))
#for(x in myrange){
#	output1 <- x^3
#	output2 <- x+4
#	output3 <- x+5
#	output4 <- x+6
#	output5 <- matrix((1:20)+x, ncol=5, nrow=4)
#	Basiclist[[cnt]]$Macht <- output1
#	Basiclist[[cnt]]$Lineair[1] <- output2
#	Basiclist[[cnt]]$Lineair[2] <- output3
#	Basiclist[[cnt]]$Lineair[3] <- output4
#	Basiclist[[cnt]]$Matrix <- output5
#	namen <- c(namen, paste("Gen", x, sep=""))
#	cat(namen, "\n")
#	cnt <- cnt+1
#}
#names(Basiclist) <- namen

which(!unlist(lapply(Cent,"[",1))=="FALSE")
which(unlist(lapply(Cent,"[",1))=="FALSE")

which(lapply(Cent,length)==1)
which(lapply(Cent,length)!=1)


#Gevonden Centrotypes
#> Centrotypes1...200 <- which(unlist(lapply(SUPERLIJST,"[","ok")))
#> Centrotypes1...200
# ok  ok  ok  ok  ok  ok  ok  ok  ok  ok  ok  ok  ok 
# 17  31  86 104 108 136 137 138 145 176 179 185 191 

#print.hotjeknor <- function(x){
#  cat("Results for ",length(x),"Genes")
#  cat("Has ",length(which(unlist(lapply(x,"[","ok")))),"c")
#}

# Kijken welk type centrotyping meer oplevert

PowerLevel <- function(x){
	cat("Vegeta, what does the scouter say about his powerlevel?", "\n")
	Sys.sleep(2)
	Goten <- invisible(readCentrotypes(paste(x, ".txt", sep="")))
	Trunks <- which(lapply(Goten,length)!=1)
	Gotenks <- Goten[Trunks]
	SsjGotenks <- NULL
	for(Dragonball in 1:length(Gotenks)){
		SsjGotenks <- c(SsjGotenks, Gotenks[[Dragonball]][-c(1:3)])
	}
	powerlevel <- length(SsjGotenks)
	print(powerlevel)
	
	Sys.sleep(1)
	if(powerlevel > 9){cat("\n","\n", "OVER 9.000!?!?!?", "\n", "\n", "\n")}
	invisible(powerlevel)
}

# unadjusted   = pwrlvl 1914
# adjusted = pwrlvl 1834
# Conclusie: de unadjusted zal worden gebruikt, puur voor het hebben van een voorbeeld dat nog een beetje ergens op lijkt!!


OnlyCent <- function(x){
	Goten <- invisible(readCentrotypes(paste(x, ".txt", sep="")))
	Trunks <- which(lapply(Goten,length)!=1)
	Gotenks <- Goten[Trunks]
}

Centrotype_unadjusted <- OnlyCent("Centrotype_unadjusted")

####################################
preCentrotype <- function(x){
	varA <- Centrotype_unadjusted[[x]][-c(1:3)]				#inlezen van de probes
	varB <- read.table(paste("chr1/", x, ".txt", sep=""))	#inlezen van de dataset
	varC <- varB[as.numeric(varA),]							#alleen de data van de gewenste probes
	varD <- apply(varC[goodind], 2, mean)					#mean per ind van de probes
	output <- c(x, Centrotype_unadjusted[[x]][2], varD)		#output in: gennaam, marker, ind mean expression
}

preProbes <- function(x){
	varB <- read.table(paste("chr1/", x, ".txt", sep=""))	#inlezen van de dataset
	varD <- apply(varB[goodind], 2, mean)					#mean per ind van de probes
	output <- c(x, Centrotype_unadjusted[[x]][2], varD)		#output in: gennaam, marker, ind mean expression
}	


####################################

CentrotypeMatrix <- function(x = 1:length(Centrotype_unadjusted), y = preCentrotype){
	vectorcentrotype <- NULL
	colnamen <- NULL
	for(ele in x){
		gennaam <- names(Centrotype_unadjusted[ele])
		enkelcentrotype <- y(gennaam)
		vectorcentrotype <- c(vectorcentrotype, enkelcentrotype)
		cat(gennaam, "\n")
		colnamen <- c(colnamen, gennaam)
	}
	Centrotypematrix <- matrix(vectorcentrotype, nrow=167, ncol=length(x))
	colnames(Centrotypematrix) <- colnamen
	Centrotypematrix
}

centrotype_matrix 	<- CentrotypeMatrix()
probes_matrix		<- CentrotypeMatrix(, preProbes)
samplev2 <- read.table("samplev2.txt")
rawEnvirement <- as.numeric(samplev2[,3])		# 1 = 6H, 2 = Dry_AR, 3 = Dry_FRESH, 4 = RP 
Envirement <- rawEnvirement[-c(rawfaulty)]

##### ANOVA test #####
UltimateCentrotypeMatrix <- function(x = centrotype_matrix){
mm <- NULL
centroidprobes <- as.matrix(apply(x[11:159,],2,as.numeric))
tijdA <- proc.time()[3]
vals <- NULL
mm <- apply(newgenomatrix,2,function(marker){
	unlist(lapply(summary(aov(centroidprobes ~ Envirement[9:157]+as.factor(marker)+Envirement[9:157]:as.factor(marker))),"[",2,5))
})
tijdB <- proc.time()[3]
cat("Took:", tijdB-tijdA, "seconds", "\n", sep=" ")
colnames(mm) <- paste("Marker", 1:69, sep="")
rownames(mm) <- colnames(x) 
mm
}
#####################

ulTimatecentRotypEmatriX <- UltimateCentrotypeMatrix()
TREX <- -log10(ulTimatecentRotypEmatriX)

UltimateProbeMatrix <- UltimateCentrotypeMatrix(probes_matrix)
UPM <- -log10(UltimateProbeMatrix)

       
######################################################################
###################### Geniet en bewonder!!! #########################
######################################################################

for(ele in 1:length(TREX)){
	plot(TREX[ele,], type="l",col="green", main = paste("Centrotype",
	rownames(TREX)[ele], sep=" "), xlab= "Markers", ylab="-log10(pval)")
	points(UPM[ele,], type="l", col="red")
	legend("topright",c("Centrotype","Mean Probes"), col=c("green","red"),lwd=2)
	Sys.sleep(2)
}
























