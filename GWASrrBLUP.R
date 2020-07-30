install.packages("rrBLUP")
library(rrBLUP)
###

outfileABHTwithpos2 <- t(outfileABHTwithpos)
write.table(outfileABHTwithpos2, row.names = F, col.names = F, file = "GWASgeno.txt")
GWASall <- GWAS(pheno = PhenoBLUP, geno = GWASgeno,plot = FALSE)
write.table(GWASall,row.names = F, col.names = T, file = "GWASall.txt")

#random population of 200 lines with 1000 markers
M <- matrix(rep(0,200*1000),1000,200)
for (i in 1:200) {
  M[,i] <- ifelse(runif(1000)<0.5,-1,1)
}
colnames(M) <- 1:200
genorandom <- data.frame(marker=1:1000,chrom=rep(1,1000),pos=1:1000,M,check.names=FALSE)

phenotessstt <- data.frame(line=1:200,y=y)
y <- g + rnorm(200,mean=0,sd=sqrt((1-h2)/h2*var(g)))
g <- as.vector(crossprod(M,u))
QTL <- 100*(1:5) #pick 5 QTL
u <- rep(0,1000) #marker effects
u[QTL] <- 1
h2 <- 0.5
