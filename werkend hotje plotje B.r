#(c)Inne
#date created:8 maart

#een plotje
data <- read.csv("BayShatraitsAll.csv", sep=";")
chromos<- data[1,405:ncol(data)]
Morgan<- data[2,405:ncol(data)]
#plot grote

plotInne <- function(Morgan, chromos, yass1, yass2, gapsize=25,type='l',cuttoff=NULL,Title="Summarized QTL plot",Title_Y.as1="value of significant QTL"){
if(missing(Morgan)) stop(cat("Morgan afstanden van de Markers missin.")
if(missing(chromos)) stop("lijst met markers die gekoppeld zijn aan een chromosoom mist.")
if(missing(yass1)) stop("Vul de eerste te plotten lijst in.")
if(missing(yass2)) stop("Vul de tweede te plotten lijst in.")

#Getd bepaald de centimorgan afstand tussen de markers van één chromosoom
getD <- function(which.chr=1,distances,chr){
  distances[max(which(chr==which.chr))]
}

#GetCD bepaald de chromosoom distance en telt er een gap bij op
getCD <- function(which.chr=1, gapsize = 25,distances,chr){
  if(which.chr==0) return(0)
  D <- 0
  for(x in 1:which.chr){
    D <- D + getD(x,distances,chr) + gapsize
  }
  return(D)
}
#hier worden plot parameters meegegeven
  op <- par(las = 2)							#las is loodrechte x-aslabel notatie
  op <- par(cex.axis = 0.6)						#De grootte van de labels
  distances <- as.numeric(t(Morgan[1,]))		#is de afstand van het totaal
  chr <- as.numeric(unlist(t(chromos)))			#Chr is welk chromosoom en de lengte ervan
  nchr <-length(unique(chr))					#nchr is de lenget van alle chromosomen tot dusver
  plot(c(0, getCD(nchr, gapsize=gapsize,distances=distances,chr=chr)),c(0,max(max(yass1), max(yass2))),type="n",main=Title, ylab=Title_Y.as1, xlab="Markers",xaxt="n")	#De algemene plot
  axis(2, c(0,max(max(yass1),max(yass2))))		#De Y-as
  locs <- NULL									#Locaties op de x-as
  for(x in 1:nchr){								#De loop voor het plotten van je points
    locs <- c(locs,distances[which(chr==x)] + (getCD(x-1,gapsize,distances,chr)))
    points(x=distances[which(chr==x)] + (getCD(x-1,gapsize,distances,chr)),y=yass1[which(chr==x)],type=type,col="green",lwd=1)
	points(x=distances[which(chr==x)] + (getCD(x-1,gapsize,distances,chr)),y=yass2[which(chr==x)],type=type,col="red",lwd=1)
  }
  abline(h=cuttoff, lty="dashed")				#de cut-off
  axis(1,locs,labels=names(Morgan[1,]))			#de x-as notatie
}
plotInne(Morgan,chromos,kak, pop,cuttoff=3)