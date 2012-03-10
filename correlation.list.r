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

for(u in 1:length(all_files)){
	rawX = read.table(paste("chr1/",all_files[u], sep = ""))												#Data inlezen
	if(nrow(rawX) >= 2){
		dX = rawX[,-c(1:16)]																				#Bijgevoegde informatie weghalen
		cleanX = dX[,-c(128, 130, 131, 132, 133, 134, 135, 137, 138, 139, 141, 142, 143, 144, 145)]			#Lijpe individuen weghalen
		ncleanX = normalizeQuantiles(cleanX)																#Genormalizeerde data
		
		#Make sure there are no traits which show no variance !!!! this will propagate NA's
		variances <- apply(ncleanX,1,var)
		if(any(variances == 0)) ncleanX <- ncleanX[-which(variances==0),] 
		if(nrow(ncleanX) >= 2){
			abscormatrix <- abs(cor(t(ncleanX),use="pair"))
			(cat(u, "\n"))
			abscormatrix <- cleancor(abscormatrix)
# Een probe is gecorreleerd als de gegeven correlatie hoger is dan de mean+(2*sd)							#Eerste these
			meancor = mean(apply(abscormatrix,2,mean, na.rm = TRUE))										#Vaststellen gemiddelde cor.
			sdcor   = mean(apply(abscormatrix,2,sd, na.rm = TRUE))											#Vaststellen sd cor.

# Een lijst aanmaken met de lengte van mijn hoeveelheid probes per gen
			cor.probe.lijst = vector("list",ncol(abscormatrix))
			x = 1
			while(x <= ncol(abscormatrix)){
				cor.probe = which(abscormatrix[x,] >= (meancor+(2*sdcor)) & abscormatrix[x,] >= 0.6)		#Kijken met welke hij cor. met minimale cor. == 0.6  
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
}
