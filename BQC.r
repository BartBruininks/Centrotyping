# \file Basic quality control
#  Do basic quality control on RNA expression data
# (c) 2012 - Bart Bruininks, written in the [R] language

#maak een matrix met 180 random rows uit 180 random genen op chr1
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