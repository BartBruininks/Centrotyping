################################################
## Creëer een cetrotype m.b.v. de langste rij ##
################################################
source("http://bioconductor.org/biocLite.R")									# Inlezen van packages
biocLite("limma")
require("limma")

setwd("C:\\Users\\Bart\\Desktop\\Leren Programmeren\\Centrotype")				# Hoofdlocatie en mappen inlezen
all_genes <- dir("chr1/")
all_rawCentrotypes <- dir("cor0.6/")

# Inlezen van een rawCentrotype en omzetten naar een lijst
read <- function(x, y = "cor0.6/"){
  filein <- file(paste(y, all_rawCentrotypes[x], sep=""), "rt") #rt zorgt ervoor dat ik mijn file open laat staan terwijl ik hem lees
  nline <- readLines(filein, n=1,)
  s <- proc.time()
  splitit <- NULL
  namez <- NULL
  cnt <- 1
  while(length(nline)){
    elementen <- strsplit(nline, " ")
    splitit <- c(splitit, elementen)
	splitit[[cnt]] <- splitit[[cnt]][-c(1)]
	namez <- c(namez, elementen[[1]][1])
	nline <- readLines(filein, n=1)
	cnt <- cnt+1
  }
  e <- proc.time()
  cat(e[3]-s[3], "seconds\n")
  close(filein)
  names(splitit) <- paste("p", namez, sep="")
  invisible(splitit)
}
#lengte van de langste cor. probe met als aa de gewenste lijst met cor. probes
enkelemaxprobelength = function(x,aa){
	tempgene <- strsplit(all_rawCentrotypes[x], split="[.]")
	gene <- tempgene[[1]][1]
	ncorrelatedprobes <- max(unlist(lapply(aa, length )))
	mostcorr <- names(which.max(unlist(lapply(aa, length ))))
	return(c(gene,mostcorr,ncorrelatedprobes))
}

maxprobelength <- function(x){
	aa <- read(x)
	output <- enkelemaxprobelength(x,aa)
	output
}

# Voor alle genen de lengte van de langste cor. probe 
allmaxprobelength <- function(x = "cor0.6/"){
	all_rawCentrotypes <- dir(x)
	lijst <- NULL
	for(ele in 1:length(all_rawCentrotypes)){
		temp <- maxprobelength(ele)
		lijst <- c(lijst, temp)
	}	
	lijst
}

# Een matrix met de probes die met de meeste andere probes correleren binnen één gen,
# !!!Nog toe te voegen is de gemiddelde expressie van de individuën, iets met de baseparen en de richting!!! 

#### Tussenstation ####
matrixCentrotypes <- function(x = "cor0.6/"){
	matrix(allmaxprobelength(x), ncol=(3), nrow=length(all_rawCentrotypes), byrow=TRUE)
}
#### Tussenstation ####