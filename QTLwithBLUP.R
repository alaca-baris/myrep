library(qtl)
data("fake.bc")
#################################
?read.cross
QTL2 <- read.cross(format = "csvs", dir = "C:/R(MoveLater)/", genfile = "outfileABHTLast3.csv",
                  phefile = "OutGoettingen.csv", crosstype = "dh",
                  map.function = "morgan", estimate.map = FALSE)

QTL <- QTL2

############################

summary(QTL)
plotMissing(QTL)
plotGeno(QTL)
plotMap(QTL, show.marker.names = FALSE)
QTL <- drop.nullmarkers(QTL)
totmar(QTL)
QTLnew <- est.rf(QTL)
plotPheno(QTL, pheno.col = 3)
fix(QTL)
QTL <- calc.errorlod(QTL, error.prob = 0.01)
top.errorlod(QTL)
summary(top.errorlod(QTL))
#################################################################################

result_hk <- scanone(QTL, pheno.col=2, method= "hk") #Scan the first column using Haley-Knott
result_em <- scanone(QTL, pheno.col=4, method= "em")     #Use expectation maximization
result_imp <- scanone(QTL, pheno.col=4, method= "imp")   #Use imputation

plot(result_hk, col = "red", lwd = "1")
legend

plot(result_hk,result_em,result_imp,lwd=c(1,1,1),col=c("red","green","blue"))
legend("topleft",c("hk","em","imp"),lwd=c(1,1,1),col=c("red","green","blue"))
save.image()

##################################################################################

perm_hk <- scanone(QTL, pheno.col=3, method= "hk",n.perm=1000)
plot(perm_hk)
summary(perm_hk)  #print the thresholds
plot(result_hk, ylim = c(0,5), lwd=c(1,1), col= c("red", "blue"))
abline(h=summary(perm_hk)[[1]],col="green",lty=2,lwd=2)
abline(h=summary(perm_hk)[[2]],col="orange",lty=2,lwd=2)
###legend("topleft", c("outcim","hk"), lwd = c(1,1), col = c("red", "blue"))


perm_imp <- scanone(QTL, pheno.col=13, method= "imp",n.perm=1000)
plot(perm_imp)
summary(perm_imp)  #print the thresholds
plot(result_imp, chr = 5)
abline(h=summary(perm_imp)[[1]],col="green",lty=2,lwd=2)
abline(h=summary(perm_imp)[[2]],col="orange",lty=2,lwd=2)

summary(result_imp, threshold = summary(perm_imp)[[1]])
summary(result_imp, threshold = summary(perm_imp)[[2]])

perm_em <- scanone(QTL, pheno.col=13, method= "em",n.perm=1000)
plot(perm_em)
summary(perm_em)  #print the thresholds
plot(result_em)
abline(h=summary(perm_em)[[1]],col="green",lty=2,lwd=2)
abline(h=summary(perm_em)[[2]],col="orange",lty=2,lwd=2)

##################################################################################

result_hk <- scanone(QTL, pheno.col=1, method= "hk")
plot(result_hk, col = "red", lwd = "1")

perm_hk <- scanone(QTL, pheno.col=1, method= "hk",n.perm=5000)
plot(perm_hk)
summary(perm_hk)  #print the thresholds
plot(result_hk)
abline(h=summary(perm_hk)[[1]],col="green",lty=2,lwd=2)
abline(h=summary(perm_hk)[[2]],col="orange",lty=2,lwd=2)

##################################################################################

data(QTL)
QTL <- calc.genoprob(QTL, step=2.5)

?cim
out <- scanone(hyper)
out.cim <- cim(QTL, n.marcovar=3, pheno.col = 11, method = "hk")
plot(out.cim)
abline(h=0.72, col="green", lty = 2, lwd=2)
plot(result_imp, out.cim ,col=c("blue", "red"))
legend("topleft",c("imp","CIM"),lwd=c(1,1),col=c("blue","red"))

perm_hk <- scanone(QTL, pheno.col=3, method= "hk",n.perm=1000)
plot(perm_hk)
summary(perm_hk)  #print the thresholds
plot(result_hk, ylim = c(0,5))
abline(h=summary(perm_hk)[[1]],col="green",lty=2,lwd=2)
abline(h=summary(perm_hk)[[2]],col="orange",lty=2,lwd=2)

add.cim.covar(out.cim)
#################################################################################
#################################################################################
#################################################################################
QTL2 <- read.cross(format = "csvs", dir = "C:/R(MoveLater)/", genfile = "outfileABHTLast3.csv",
                   phefile = "TestPhenoCovBLUP2.csv", crosstype = "dh",
                   map.function = "morgan", estimate.map = FALSE)


############################
est.rf(QTL2)
newmap <- est.map(QTL2, omit.noninformative = TRUE, map.function = "morgan")
plotMap(QTL2, newmap)

summary(QTL2)
plotMissing(QTL2)
plotGeno(QTL2)
plotMap(QTL2, show.marker.names = FALSE)
QTL2 <- drop.nullmarkers(QTL2)
QTL2 <- drop.dupmarkers(QTL2)

QTL2 <- sim.geno(QTL2, n.draws=16, step=2, off.end=0, error.prob=0.0001,map.function="morgan",stepwidth="fixed")
QTL2 <- calc.genoprob(QTL2, step=2)
QTL2sub <- reduce2grid(QTL)

totmar(QTL2sub)

totmar(QTL2)
QTLnew <- est.rf(QTL2)
plotPheno(QTL2, pheno.col = 3)
fix(QTL)

#################################################################################

result_hk2 <- scanone(QTL2, pheno.col=3, method= "hk") #Scan the first column using Haley-Knott
result_em2 <- scanone(QTL2, pheno.col=11, method= "em")     #Use expectation maximization
result_imp2 <- scanone(QTL2, pheno.col=11, method= "imp")   #Use imputation

plot(result_hk2, col = "red", lwd = "1")


plot(result_hk,result_em,result_imp,lwd=c(1,1,1),col=c("red","green","blue"))
legend("topleft",c("hk","em","imp"),lwd=c(1,1,1),col=c("red","green","blue"))

plot(result_hk,result_hk2, lwd=c(3,1),col=c("red","green"), ylim = c(0,5))
abline(h=summary(perm_hk)[[1]],col="green",lty=2,lwd=2)
abline(h=summary(perm_hk)[[2]],col="orange",lty=2,lwd=2)
legend("topleft",c("BLUPs","Goettingen"),lwd=c(1,1),col=c("red","green"), title = "Grain Yield")
save.image()

##################################################################################

?scantwo
QTL2 <- calc.genoprob(QTL2, step = 5, error.prob = 0.01)
out2.hk <- scantwo(QTL2, method = "hk", pheno.col = 3)
summary(out2.hk)
plot(out2.hk, lower = "fv1")
### We can also run perm test to see LOD scores (Time consuming)
operm2.hk <- scantwo(QTL, method="hk", n.perm=10, pheno.col = 3)
summary(out2.hk)
summary(operm2.hk)
summary(out2.hk, perms=operm2.hk, pvalues=TRUE, alphas=c(0.05, 0.05, 0, 0.05, 0.05))
setwd("C://R(MoveLater)//Genofile")
load("out2.hk.Rdata")
