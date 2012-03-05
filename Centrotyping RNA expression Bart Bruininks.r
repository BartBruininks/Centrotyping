# (c) Bart Bruininks written in the [R] language
# 02-03-2012 	10.15u-16.00u
# Computational biology 
# Research opdracht "centrotyping RNA expression"

setwd("C:/Users/Bart/Desktop/Leren Programmeren/Centrotype")

#Data analyse Chr 1
AT1G01010 = read.table("chr1/AT1G01010.txt")
dAT1G01010 = AT1G01010[,-c(1:16)] 					# alleen de meetwaarden van AT1G01010
dAT1G01010
image(t(dAT1G01010))

# Sorteren van de samplev2 op zaadconditie en image van de gesoorteeerde matrix van dAT1G01010
samplev2 = read.table("samplev2.txt")
samplev2[sort.list(samplev2[,3]), ]	
volgorde = sort.list(samplev2[,3])	
jpeg(filename = "zaadconditieAT1G01010.jpeg", width = 2000, height = 2000, bg = "white", units = "px")
image(t(dAT1G01010[,volgorde]))
dev.off()
jpeg(filename = "heatmapzaadconditieAT1G01010.jpeg", width = 2000, height = 2000, bg = "white", units = "px")
heatmap(t(dAT1G01010[,volgorde]), Rowv = NA)



# Heatmap van AT1G01010
jpeg(filename = "heatmapAT1G01010.jpeg", width = 2000, height = 2000, bg = "white", units = "px")
heatmap(t(dAT1G01010[,volgorde]))
dev.off()

# Kijken of er bepaalde plekken in het genoom is waar crossingover vaker voor komt
genotypes <- read.table("genotypes.txt",sep="\t")
dim(genotypes)		# controleren of de dimensies van je var goed zijn

paste("RIL",genotypes[,1],sep="")		# kijken of de volgorde van de columns van genotype en sample
as.character(samplev2[,2]) %in% paste("RIL",genotypes[,1],sep="")
which(as.character(samplev2[,2]) %in% paste("RIL",genotypes[,1],sep=""))
dgenotypes = genotypes[,-c(1:7)]
image(t(dgenotypes))

numericgeno <- matrix(unlist(apply(dgenotypes,2,function(x){
  lapply(x,function(e){
    which(LETTERS %in% e)
  })
})),nrow(dgenotypes),ncol(dgenotypes))
numericgeno
heatmap(numericgeno, Colv = NULL)

jpeg(filename = "Crossing-overhotspots.jpeg", width = 2000, height = 2000, bg = "white", units = "px")		#Crossingover hotspots
image(t(numericgeno), xlab = "genotypes", ylab = "individuals", main = "crossing-over hotspots")			#Hotspots gevonden, alleen nog een goede as nummering vinden
dev.off()

boxplot(dAT1G01010)
cordAT1G01010 = cor(dAT1G01010, dAT1G01010)
image(cordAT1G01010)
heatmap(cordAT1G01010)
heatmap(cordAT1G01010, scale = "none")
hm = heatmap(cordAT1G01010, scale = "none")

hm = heatmap(cordAT1G01010, scale = "none",keep.dendro=T)
plot(hm[[3]])

boxplot(cor(dAT1G01010[,1:8]))		#Ouder Bay-0 correlatie onderling
boxplot(cor(dAT1G01010[,173:180]))	#Ouder Sha
boxplot(cor(dAT1G01010[,9:172]))	#Nageslacht

# Normalizeren van dAT1G01010 en quality control
normdAT1G01010 = normalizeQuantiles(dAT1G01010)

jpeg(filename = "imagenormdAT1G01010.jpeg", width = 2000, height = 2000, bg = "white", units = "px")
image(cor(normdAT1G01010))
dev.off()			

jpeg(filename = "imagenormdAT1G01010_Bay-0.jpeg", width = 1000, height = 1000, bg = "white", units = "px")
image(cor(normdAT1G01010[,1:8]))
dev.off()

jpeg(filename = "imagenormdAT1G01010_Sha.jpeg", width = 1000, height = 1000, bg = "white", units = "px")
image(cor(normdAT1G01010[,173:180]))
dev.off()

jpeg(filename = "dendogramnormdAT1G01010.jpeg", width = 3000, height = 3000, bg = "white", units = "px")
hm = heatmap(cordAT1G01010, scale = "none",keep.dendro=T)
plot(hm[[3]])
dev.off()

jpeg(filename = "boxplotAT1G01010.jpeg", width = 1000, height = 1000, bg = "white", units = "px")
boxplot(dAT1G01010)
dev.off()


##############################################
##	Alle correlerende probes samenvoegen	##
##############################################

setwd("C:/Users/Bart/Desktop/Leren Programmeren/Centrotype")

# Simpel idee voor latere algemenisering
rawAT1G01010 = read.table("chr1/AT1G01020.txt")												#Data inlezen
dAT1G01010 = rawAT1G01010[,-c(1:16)]														#Bijgevoegde informatie weghalen
cleanAT1G01010 = dAT1G01010[,-c(142, 145, 135, 128, 133, 137, 143, 138, 132, 139)]			#Lijpe individuen weghalen
ncleanAT1G01010 = normalizeQuantiles(cleanAT1G01010)										#Genormalizeerde data

# Een probe is gecorreleerd als de gegeven correlatie hoger is dan de mean+(2*sd)			#Eerste these
meancor = (apply(abs(cor(t(ncleanAT1G01010))),2,mean))										#Vaststellen gemiddelde cor.
sdcor   = (apply(abs(cor(t(ncleanAT1G01010))),2,sd))										#Vaststellen sd cor.

# Een lijst aanmaken met de lengte van mijn hoeveelheid probes per gen
cor.probe.lijst = vector("list",length(meancor))
mynames <- NULL
x = 1
while(x <= length(meancor)){
	cor.probe = which((abs(cor(t(ncleanAT1G01010)))[x,]) >= (meancor+abs(sdcor)))			#Kijken met welke hij cor. 
	cor.probe.lijst[[x]] = names(cor.probe)
	mynames <- c(mynames, names(meancor[x]))
	x = x+1
}
names(cor.probe.lijst) <- mynames

# Lijst wegschrijven
cat("",file="t.txt")
for(mlist in cor.probe.lijst){
  cat(mlist,"\n",file="t.txt",append=TRUE)
}
#####################
# Lijst weer ophalen#
#####################
myfile <- file("t.txt")
mylist <- vector("list",0)
mline <- readLines(myfile,n=1)
while(mline){
  elem <- strsplit(mline,"\t")
  mylist <- c(mylist,elem)
  mline <- readLines(myfile,n=1)
}
close(myfile)
mylist

######################################################
##	Kijken of het te veralgemenizeren is	(GELUKT)##
######################################################
source("http://bioconductor.org/biocLite.R")
biocLite("limma")
require("limma")

setwd("C:\\Users\\Students\\Desktop\\Hoofdmap")
all_files <- dir("chr1/")

cleancor <- function(a){
	for(x in 1:nrow(a)){
		a[x,x] <- NA
    }
    return(a)
}

u=1
for(u in 1:length(all_files)){
# Simpel idee voor latere algemenisering
	rawX = read.table(paste("chr1/",all_files[u], sep = ""))										#Data inlezen
	if(nrow(rawX) >= 2){
		dX = rawX[,-c(1:16)]														#Bijgevoegde informatie weghalen
		cleanX = dX[,-c(128, 130, 131, 132, 133, 134, 135, 137, 138, 139, 141, 142, 143, 144, 145)]			#Lijpe individuen weghalen
		ncleanX = normalizeQuantiles(cleanX)										#Genormalizeerde data
		
		#Make sure there are no traits which show no variance !!!! this will propagate NA's
		variances <- apply(ncleanX,1,var)
		if(any(variances == 0)) ncleanX <- ncleanX[-which(variances==0),] 
		if(nrow(ncleanX) >= 2){
			abscormatrix <- abs(cor(t(ncleanX),use="pair"))
			(cat(u, "\n"))
			abscormatrix <- cleancor(abscormatrix)
# Een probe is gecorreleerd als de gegeven correlatie hoger is dan de mean+(2*sd)		#Eerste these
			meancor = mean(apply(abscormatrix,2,mean, na.rm = TRUE))										#Vaststellen gemiddelde cor.
			sdcor   = mean(apply(abscormatrix,2,sd, na.rm = TRUE))										#Vaststellen sd cor.

# Een lijst aanmaken met de lengte van mijn hoeveelheid probes per gen
			cor.probe.lijst = vector("list",ncol(abscormatrix))
			x = 1
			while(x <= ncol(abscormatrix)){
				cor.probe = which(abscormatrix[x,] >= (meancor+(2*sdcor)))				#Kijken met welke hij cor. 
				cor.probe.lijst[[x]] = c(rownames(ncleanX)[x],names(cor.probe))
				x = x+1
				names(cor.probe.lijst) <- rownames(ncleanX)

# Lijst wegschrijven
				cat("",file = all_files[u])
				for(mlist in cor.probe.lijst){
					cat(mlist,"\n", file = all_files[u], append=TRUE)
				}
			}	
		}
	}
	u <- u+1
}

#### Nomgaals quality control ####

# maak een matrix met 180 random rows uit 180 random genen op chr1
#colRemoval is specific for our case in which we dont want columns 1 to 16
BQC = function(x, directory = "chr1/", colRemoval = c(1:16)){
	if(missing(x)) stop("Please specify a number of rows!")
	setwd("C:/Users/Bart/Desktop/Leren Programmeren/Centrotype")
	pool <- dir(directory)
	tempmatrix1 = NULL
	for(n in 1:x){
		randomgenenaam = sample(pool, 1, replace = FALSE)
		randomgene = read.table(paste("chr1/",randomgenenaam, sep=""))
		randomrownumber = sample(nrow(randomgene), 1) 
		randomrow = randomgene[randomrownumber,]
		tempmatrix1 = rbind(tempmatrix1, randomrow)
		cat(nrow(tempmatrix1),"\n")
	}
	matrix1 = tempmatrix1[,-colRemoval]
	return(invisible(matrix1))
}
## Bereken het gemiddelde van de correlatie van deze waarden en de gemmiddelde sd
normStats = function(x){
	normMatrix = normalizeQuantiles(x)
	cornormMatrix = cor(normMatrix)
	meancorMatrix = mean(cornormMatrix)
	sdcorMatrix = apply(cornormMatrix, 2, sd)
	meansdcorMatrix = mean(sdcorMatrix)
	output = c(meancorMatrix, meansdcorMatrix)
	output
}
### In deze functie roep ik iets op wat ik al berekend heb ik de functie die er in zit
# Dit kost natuurlijk meer rekenkracht, kan dat ook anders? 
Significant = function(x){
	temp = normStats(x)
	normMatrix = normalizeQuantiles(x)
	cornormMatrix = cor(normMatrix)
	output = c(which(apply(cornormMatrix, 2, mean) < (temp[1]-temp[2])), temp)
	output
}

######## Run BQC en Significant met 2000 #########
aa = BQC(10000)
cc = Significant(aa)

# uitkomsten zijn: (10000 random rows)
#        V128         V133         V134         V135         V137         V139 
#128.00000000 133.00000000 134.00000000 135.00000000 137.00000000 139.00000000 
#        V141         V143         V144         V145   gem.         gem. sd.               
#141.00000000 143.00000000 144.00000000 145.00000000   0.74049332   0.08308256

#Deze individuen zitten allemaal in batch 4. Batch 4 heb ik helemaal verwijdert! 

#####################################################################################
## Probes controleren of ze niet vaker dan één keer binden op het DNA m.b.v. BLAST ##
#####################################################################################
	setwd("C:/Users/Bart/Desktop/Leren Programmeren/Centrotype")
	all_files <- dir("chr1/")
	u <- 1
	for(u in 1:length(all_files)){
		rawX <- read.table(paste("chr1/",all_files[u], sep = ""))
		probeX <- rawX[,-c(1:4, 7:196)]
		probeX
		write.table(probeX, file = paste("sequentie_probes", all_files[u], sep=""))
	}

	
################################################
## Creëer een cetrotype m.b.v. de langste rij ##
################################################
source("http://bioconductor.org/biocLite.R")			# Inlezen van packages
biocLite("limma")
require("limma")

setwd("C:\\Users\\Bart\\Desktop\\Leren Programmeren\\Centrotype")			# Hoofdlocatie en mappen inlezen
all_genes <- dir("chr1/")
all_rawCentrotypes <- dir("rawCentrotypes/")

# Inlezen van een rawCentrotype en omzetten naar een lijst
read <- function(x){
  filein <- file(paste(x, all_rawCentrotypes[x], sep= ""), "rt") #rt zorgt ervoor dat ik mijn file open laat staan terwijl ik hem lees
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
allmaxprobelength <- function(x = "rawCentrotypes/"){
	all_rawCentrotypes <- dir(x)
	lijst <- NULL
	for(ele in 1:length(all_rawCentrotypes)){
		temp <- maxprobelength(ele)
		lijst <- c(lijst, temp)
		ele <- ele+1
	}	
	lijst
}

# Een matrix met de probes die met de meeste andere probes correleren binnen één gen
#### Tussenstation ####
matrixCentrotypes <- function(x = "rawCentrotypes/"){
	matrix(allmaxprobelength(x), ncol=(3), nrow=length(all_rawCentrotypes), byrow=TRUE)
}
#### Tussenstation ####


 


 








