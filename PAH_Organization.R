# Example for Jeff Hatten 
# Creator: N.L. Osborne
# Created: 3/16/2015
# Edited: 1/25/18 - modified for PAH quantification

# Objective: Transform series of .csv files into easier to handle format
#Set up the .csv files so that there is only one quantification (the Chemstation
#software will append quanitifications if you do multiple).
library(reshape2)
library(plyr)

#Set your working directory where the data .csv file(s) are located
################################
####   CHANGE DIRECTORY   ######
################################
setwd("~/Desktop/GCMS-LMM")



#Load the function
# Function processes chromatograph/mass spectrophotometer files into an easy to handle format
gcms.db <- function(filename){
	
  # load the csv file
  x = read.csv(filename, header=FALSE,sep=",",na.strings=c("", "NA"))
  # Convert all blanks to consistent value: "NA"
  x[x==""] <- NA
  x[x==" "] <- NA
  # Delete any columns or rows that are completely "NA". Column 4 is always useless so it can be specifically deleted
  x <- x[,colSums(is.na(x)) < nrow(x)]
  x[,4] <- NULL
  x <- x[rowSums(is.na(x)) != ncol(x),]
  
  #name the columns of data
  names(x) <- c('compound','area','rt')
  
  
  #determine where each sample (observation) begins (i.e. set the datum)
  #this will allow all quantification times to be included in the data set (and dropped or included or whatever)
  header.breaks = grep("Header Info", as.character(x[,1]))
  
  #define where each important piece of header information lies in relation to the datum
  name.y <- header.breaks + 1
  name.x <- 1
  locate.y <- header.breaks + 1
  locate.x <- 3
  rundate.y <- header.breaks + 2
  rundate.x <- 1
  gc.method.y <- header.breaks + 2
  gc.method.x <- 2
  sample.y <- header.breaks + 3
  sample.x <- 1
  quant.method.y <- header.breaks + 5
  quant.method.x <- 1
  
  #bind all of the header information togehter into a data frame
  header.complete = data.frame(
    sample = x[sample.y, sample.x],
    name = x[name.y,name.x],
    rundate = x[rundate.y, rundate.x],
    gc.method = x[gc.method.y, gc.method.x],
    quant.method = x[quant.method.y, quant.method.x],
    locate = x[locate.y, locate.x])
  
  # Make the formatting a little more friendly
  header.complete[] = lapply(header.complete,as.character)


  #find where each sample file begins, determine the length of each sample file
  z <- grep("Header Info", as.character(x[,1]))
  length <- z[2]-1
  
  #create column of data that labels each sample with a number
  sample.no <- rep(1:nrow(header.complete), each=length)
  
  #create a column of data that labels each row with each sample
  row.no <- rep(1:length, nrow(header.complete))

  #bind together the data, observation indices, and within observation indices and remove the first 11 rows of data from eahc observation
  x <- data.frame(x, sample.no, row.no)
  x <- x[ which(x$row.no > 6), ]
  #delete extra quantifications. Chemstation automaically adds most recent quantification to top of sample "section" in data file
  #currently assumes 132 rows (compounds per sample)
  #if the number of compounds changes this value will need to be adjusted
  x <- x[ which(x$row.no < 23), ]
  x[,1:ncol(x)] <- sapply(x[,1:ncol(x)],as.character)
  

  
  #create wide format to fit with header extraction
  melt.x <- melt(x, id=c("sample.no", "compound", "row.no"))
  x.cast <- dcast(melt.x, sample.no ~ compound + variable)
  
  #sort by sample number
  x.cast[,1:33] <- sapply(x.cast[,1:33],as.numeric)
  index <- order(x.cast$sample.no)
  x.cast <- x.cast[index,]
  
  #bind the header info and data
  out = data.frame(header.complete,x.cast)
  
  #delete unneeded columns and convert columns to appropriate format
  out$sample.no <- NULL
  out[,7:38] <- sapply(out[,7:38],as.numeric)
  out[,1:6] <- sapply(out[,1:6],as.character)
  
  #place the columns in desired order (see CuOox_col_order.csv)
  order <- read.csv("ColumnOrder/PAH_col_order.csv", header=FALSE)$V1
  order <- as.character(order)
  test.df1 <- data.frame(colnames(out),order)
  order
  out <-  out[ ,order]
  
        
				return(out)
			}
	
###Runs----
#################################
###   Change CSV filenames   ####
#################################
###Add as many runs as CSV files needed to organize
mydata <- gcms.db("quechers_triplicate_3.csv")


###Output----
###Change file name
write.csv(mydata, file = "quechers-triplicate-organized_3.csv")
