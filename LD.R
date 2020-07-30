
###############################################
#outfileABHTT <- outfileABHTT[-1]
#outfileABHT <- outfileABHTTTT
outfileABHTLast3 <- outfileABH3
##Extract ID column and also, Chromosome and Position rows from the data
library(dplyr)
library(tidyr)
outfileABH  <- X2018G2FPHGABH
outfileABH <- outfileABHT
outfileABHT <- outfileABH
outtest <- outfileABH[-c(1:1), -c(1:1)]
ABHdf <- outtest

##Create data frames that you will use as storage of your data

i=1
j = i+1
ABHdf2 <- ABHdf
Storageytoremove2 <- ABHdf
Newdataframe <- Storageytoremove2[,-c(1:ncol(Storageytoremove2))]
Newdataframe2 <- Newdataframe
Newdataframe3 <- outfileABH
Newdataframe3 <- Newdataframe3[,-c(1:ncol(outfileABH))]
Newdataframe4 <- Newdataframe
Newdataframe5 <- Newdataframe4
Newdataframe6 <- Newdataframe5
missingx <- 0
missingy <- 0
Storageytoremove <- Storageytoremove2
colnumx <- c(i:{ncol(ABHdf)-1})
colnumy <- c(j:ncol(ABHdf))

## Start the loop at this point.
for (i in 1:ncol(ABHdf)){
  missingx <- 0
  missingy <- 0
  x = colnumx[i]
  y = colnumy[i]
  namex <- colnames(ABHdf[x])
  namey <- colnames(ABHdf[y])
  colnames(ABHdf)[x] <- "x"
  colnames(ABHdf)[y] <- "y"
  columnmx <- ABHdf$x
  columnmx <- as.data.frame(columnmx)
  
  ## Take out the data which has >%20 of missing markers 
  missingx <- sum(is.na(columnmx))
  missingx <- missingx/306 ##!!Don´t forget to change this value if number of individuals change!
  missingy <- sum(is.na(columnmy))
  missingy <- missingy/306 ##!!Don´t forget to change this value if number of individuals change!
  ifelse(missingx>0.003, Newdataframe4 <-cbind(Newdataframe4,ABHdf2[x]),print("<0.0032"))
  ifelse(missingy>0.003, Newdataframe4 <-cbind(Newdataframe4,ABHdf2[x]),print("<0.0032"))
  #Newdataframe4 <-cbind(Newdataframe4,ABHdf2[x])
  #print(">0.05")

  columnmx <- drop_na(columnmx)
  columnmy <- drop_na(columnmy)
  mycount <- columnmx %>% count(columnmx, sort = FALSE)
  mycount [is.na(mycount)] = 0
  
  mycount2 <- columnmy %>% count(columnmy, sort = FALSE)
  mycount2 [is.na(mycount2)] = 0
  sum1 <- sum(mycount$n)
  sum2 <- sum(mycount2$n)
  
  ##Calculate allele frequencies
  p1 <- mycount$n[1]/sum1
  p2 <- mycount$n[2]/sum1
  
  q1 <- mycount2$n[1]/sum2
  q2 <- mycount2$n[2]/sum2
  
  ##In some cases data could have only "A" or "B" in that case your r2 will be NaN.
  ##For next run we have to extract those markers from the data
  if(is.na(p1)==TRUE){
    Newdataframe6 <-cbind(Newdataframe6,ABHdf2[x]) 
    p1=0
  }
  if(is.na(p2)==TRUE){
    Newdataframe6 <-cbind(Newdataframe6,ABHdf2[x]) 
    p2=0
  }
  if(is.na(q1)==TRUE){
    Newdataframe6 <-cbind(Newdataframe6,ABHdf2[y]) 
    q1=0
  }
  if(is.na(q2)==TRUE){
    Newdataframe6 <-cbind(Newdataframe6,ABHdf2[y]) 
    q2=0
  }
  
  ##Merge two markers and count the possible combination
  ##!!! If your data has NA you could end up having "NA A" or "NA B" 
  ##If the data isn´t DH lines you have to take into account of other possibilities and change the values below
  df <- data.frame(paste(ABHdf$x,ABHdf$y)) 
  names(df)[1] <- "xy"
  mycount3 <- df %>% count(df$xy, sort = FALSE)
  colnames(mycount3)[1] <- "Haplotypes"
  colnames(mycount3)[2] <- "number"
  Haplotypes <- c("A A", "A B", "B A", "B B")
  number <- c("0","0","0","0")
  numbertest <- mycount3
  mycount4 <- data.frame(Haplotypes,number)
  mycount5 <- merge(mycount4,numbertest, by = "Haplotypes", all.x = TRUE, all.y = TRUE)
  mycount5 <- mycount5[,-2]
  colnames(mycount5)[2] <- "number"
  
  ##Merge the frequencies before you calculated
  frequencies <- data.frame(p1,p2,q1,q2)
  mycount5 <- mycount5[c(1:4),c(1:2)]
  sum3 <- sum(mycount5$number[c(1:4)], na.rm = TRUE)
  fre <- mycount5$number[c(1:4)]/ sum3
  
  ##When a new data frame has been created there will be alway NA values if the data had 0
  ##Assign is.na 0 to prevent this
  mycount6 <- data.frame(mycount5, fre)
  mycount6 [is.na(mycount6)] = 0
  frequencies2 <- frequencies
  frequencies2[1] <- mycount6$fre[1]
  frequencies2[2] <- mycount6$fre[2]
  frequencies2[3] <- mycount6$fre[3]
  frequencies2[4] <- mycount6$fre[4]
  
  colnames(frequencies2)[1] <- "x11"
  colnames(frequencies2)[2] <- "x12"
  colnames(frequencies2)[3] <- "x21"
  colnames(frequencies2)[4] <- "x22"
  
  frequenciesx <- frequencies2
  
  ##Calculate frequencies again, taking into account of "x.."s
  frequenciesx[1] <- frequencies2[1]+frequencies2[2]
  frequenciesx[2] <- frequencies2[3]+frequencies2[4]
  frequenciesx[3] <- frequencies2[1]+frequencies2[3]
  frequenciesx[4] <- frequencies2[2]+frequencies2[4]
  
  colnames(frequenciesx)[1] <- "p1"
  colnames(frequenciesx)[2] <- "p2"
  colnames(frequenciesx)[3] <- "q1"
  colnames(frequenciesx)[4] <- "q2"
  
  ##Calculate LD value !(There is also another way to calculate LD)
  D = frequencies2[1] - (frequenciesx[1]*frequenciesx[3])
  D <- as.numeric(D)
  
  ##Calculate r2 !!!Don't forget r2 could be NA if a marker is monomorphic!
  ##(Not exactly monomorphic but because of NA values the value of "A" or "B" could be absent)
  r2 <- ((D^2)/(frequencies[1]*frequencies[2]*frequencies[3]*frequencies[4]))
  r2 <- as.numeric(r2)
  r2na = is.na(r2)
  ##!!!! Don't forget to roll the number because sometimes r2 could be 0.9999 instead of 1!
  r2 <- format(round(r2, 2), nsmall = 2)
  r2 <- as.numeric(r2)
  
  ##If r2 is lower than 1 keep those values into the data frame
  ##And also r2 could be bigger than 1 because of missing markers. Keep those markers into data frame also
  ##(Those markers could probably removed because of >20% missign)
  ##If r2 is equal to 1, storage those markers into different data frame
  ##If r2na = TRUE leave them without any change because they have already stored before.
  if(r2na == FALSE){
    if(r2 < 1){
      colnames(ABHdf)[x] <- namex
      colnames(ABHdf)[y] <- namey
      print("Pass")
    } else if(r2 > 1.00){
        colnames(ABHdf)[x] <- namex
        colnames(ABHdf)[y] <- namey
        print("fail")
    } else if(r2 == 1){
      colnames(ABHdf)[x] <- namex
      colnames(ABHdf)[y] <- namey
      Newdataframe3 <- cbind(Newdataframe3,outfileABH[y])
      Newdataframe <- ABHdf[y]
      Newdataframe2 <- cbind(Newdataframe2,Newdataframe)
      print("r2=1")
    } 
  }
  if(r2na == TRUE){
    if (r2 == "NaN") {
      colnames(ABHdf)[x] <- namex
      colnames(ABHdf)[y] <- namey
      print("r2 fail")
    }
  }
  
  ## In some cases there could be some mistakes which could miss if conditions 
  ## for those cases just rename the column back again
  colnames(ABHdf)[x] <- namex
  colnames(ABHdf)[y] <- namey   
  print(i)
  print(r2)
  
  ##!!Don't forget to assign 0 for those two value
  missingx <- 0
  missingy <- 0
}  

##Extract markers which has r2=1 from the original data 
outfileABHTT <- outfileABH
outfileABHTT <- outfileABHTT[,!(names(outfileABHTT) %in% colnames(Newdataframe2))]
outfileABHTT <- outfileABHTT[,!(names(outfileABHTT) %in% colnames(Newdataframe4))]
outfileABHTT <- outfileABHTT[,!(names(outfileABHTT) %in% colnames(Newdataframe5))]
outfileABHTT <- outfileABHTT[,!(names(outfileABHTT) %in% colnames(Newdataframe6))]
outfileABH <- outfileABHTT
outtest <- outfileABH[-c(1:1), -c(1:1)]
ABHdf <- outtest
ABHdf2 <- ABHdf
outfileABHT <- outfileABH

##Save your data
setwd("C://R(MoveLater)//2018GenoandPheno//Geno")
write.csv(outfileABHTT, file = "outfileABHT.csv", row.names = FALSE)
write.csv(Newdataframe3, file = "r2-2markers.csv")
write.csv(Newdataframe4, file = "missingmarkers3.csv")
write.csv(Newdataframe5, file = "20%markers2-3.csv")
write.csv(Newdataframe6, file = "Monomarkers.csv")
###############

marker.positions <- marker.positions[1,]/1000000
write.csv(marker.positions, file = "markerposs.csv")
