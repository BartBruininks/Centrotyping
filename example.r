#
# iqtlhd.utils.R
#
# copyright (c) 2012 Danny Arends
# last modified Mar, 2012
# first written Mar, 2012
# 

readThroughBigFile <- function(){
  filein <- file("t.txt")

  line <- readLines(filein, n=1)
  s <- proc.time()
  while(length(line)){
    #<DO STUFF>

    line <- readLines(filein, n=1)
  }
  e <- proc.time()
  cat(e[3]-s[3], "seconds\n")
  close(filein)
}
