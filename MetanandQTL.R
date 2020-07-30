library(metan)
library(qtl)
install.packages("lmem.qtler")
library(lmem.qtler)
###################
####Calculating the BLUP values######
?gamem_met
mixed_mod <- gamem_met(g2f_2018_hybrid_data_clean_tester,
                        env = Location,
                        gen = Pedigree,
                        rep = Replicate,
                       block = Tester,
                        resp = everything(),
                        verbose = FALSE)

mixed_mod1 <- gamem_met(TestPhenoCov,
                       env = Location,
                       gen = Pedigree,
                       rep = Replicate,
                       block = Tester,
                       resp = 'GrainYield[bu/A]',
                       verbose = FALSE)
#####Delete later#####
Positions <- read.csv("C:/R(MoveLater)/2018GenoandPheno/Geno/Positions.txt", row.names=1, sep="")
library(data.table)
Positions <- transpose(Positions)
write.table(Positions, file = "C:/R(MoveLater)/2018GenoandPheno/Geno/Positions2.txt", row.names = FALSE, col.names = FALSE)
subsetting <- outfileABHT[1]
subsetting <- subsetting[-1,]
sub2 <- data3[,!(rownames(subsetting) %in% rownames(data3))]
############

plot(mixed_mod1) # The S3 generic function plot() is used to generate diagnostic plots of residuals of the model.
plot(mixed_mod1, type = "re") # The normality of the random effects of genotype and interaction effects may be also obtained by using type = "re".
data <- get_model_data(mixed_mod, "lrt") # The output LRT contains the Likelihood Ratio Tests for genotype and genotype-vs-environment random effects.
write.table(data, file = "C://R(MoveLater)//Results//Pheno//Data.txt")
data2 <- get_model_data(mixed_mod) # Estimation of heritability, coefficient of determination of the interaction (GEIr2), 
                                   # Heritability of means (h2gm), Accuracy (Ac), genotype enviroment corelation (rge), CVg and CVr are the the genotypic coefficient of variation and the residual coefficient of variation estimated
write.table(data2, file = "C://R(MoveLater)//Results//Pheno//Data2.txt")
data3 <- get_model_data(mixed_mod, what = "blupg") # Get the blup values
write.table(data3, file = "C://R(MoveLater)//Results//Pheno//Data3.txt")
data4 <- get_model_data(mixed_mod, "details") # Get the usefull information from the model
write.table(data4, file = "C://R(MoveLater)//Results//Pheno//Data4.txt")



########
## Not run: 
data (SxM_geno)
data (SxM_map)
data (SxMxE_pheno)
data (DHpop_pheno)

P.data <- SxMxE_pheno
G.data <- SxM_geno
map.data <- SxM_map

cross.data <- qtl.cross (P.data, G.data, map.data, cross='dh',
                         heterozygotes=FALSE)

summary (cross.data)

## Pheno Quality
pq.diagnostics (crossobj=cross.data, boxplot =FALSE)

## Marker Quality
mq.diagnostics (crossobj=cross.data,I.threshold=0.1,
                p.val=0.01,na.cutoff=0.1)

# QTL_SIM
QTL.result <- qtl.memq (crossobj = cross.data, P.data = P.data,
                        env.label = c('ID91','ID92','MAN92','MTd91',
                                      'MTd92','MTi91','MTi92','SKs92','WA91','WA92'),
                        trait = 'yield', step = 10, method = 'SIM',
                        threshold = 'Li&Ji', distance = 50, cofactors = NULL,
                        window.size = 50)

## QTL_CIM
QTL.result <- qtl.memq (crossobj = cross.data, P.data = P.data,
                        env.label = c('ID91','ID92','MAN92','MTd91','MTd92',
                                      'MTi91','MTi92','SKs92','WA91','WA92'),
                        trait = 'yield', step = 10, method = 'CIM',
                        threshold = 'Li&Ji', distance = 50,
                        cofactors = QTL.result$selected$marker, window.size = 50)

## End(Not run)
######

data("data_ge")
### Assign Genotypes,Replications and Locations as factors.
TestPhenoCov$Genotypes <- as.factor(TestPhenoCov$Genotypes) 
TestPhenoCov$Replication<- as.factor(TestPhenoCov$Replication)
TestPhenoCov$Location <- as.factor(TestPhenoCov$Location)

inspect(TestPhenoCov, plot = TRUE)
TestPhenoCovA <- impute_missing_val(TestPhenoCov)

model <- performs_ammi(TestPhenoCov,
                       env = Location,
                       gen = Genotypes,
                       rep = Replication,
                       resp = everything(),
                       verbose = FALSE)


get_model_data(ge2, type = "GEN", verbose = TRUE)
predict(ge3)


?ge_reg()
ge_r <- ge_reg(TestPhenoCov, env =  Location, gen = Genotypes, rep = Replication,
               resp = everything(), verbose = FALSE)


desc_stat(TestPhenoCov)
get_model_data(model, "ASV")

bakbak <- model %>%
  AMMI_indexes() %>%
  get_model_data("ASV")

###########################################
means_ge <- ge_means(TestPhenoCov, env = Location, gen = Genotypes, resp = everything())

TestPhenoCov2<- remove_rows_na(TestPhenoCov2)
stats <- ge_stats(TestPhenoCov, env = Location, gen = Genotypes,verbose = FALSE,prob = 0.05)

#ge <- ge_effects(TestPhenoCov,type = "ge", env = Location, gen = Genotypes, resp = everything(),verbose = TRUE)## Get GE interactions
#ge2 <- ge_effects(TestPhenoCov,type = "gge", env = Location, gen = Genotypes, resp = everything(),verbose = TRUE) ## G+GE
?gge
ge3 <- gge(TestPhenoCov2,env = Location, gen = Genotypes, centering = "global") ## E+G+GE 
ge4 <- gge(TestPhenoCov2,env = Location, gen = Genotypes, centering = "environment") ## G+GE
ge5 <- gge(TestPhenoCov2,env = Location, gen = Genotypes, centering = "double") ## GE

## E+G+GE 
PollenDAPDaysEGGE <- ge3$`PollenDAP(Days)`[9]
SilkDAPDaysEGGE <- ge3$`SilkDAP(Days)`[9]
PHEGGE <- ge3$`PlantHigh(cm)`[9]
EHEGGE <- ge3$`EarHigh(cm)`[9]
StandCountnumEGGE <- ge3$`StandCount(number)`[9]
StandCountEGGE <- ge3$`StandCount%`[9]
RootLodgingEGGE <- ge3$`RootLodging(number)`[9]
StalkLodgingEGGE <- ge3$`StalkLodging(number)`[9]
MoistureEGGE <- ge3$`Moisture%`[9]
DrymatterEGGE <- ge3$`Dry Matter%`[9]
YieldEGGE <- ge3$`GrainYield[bu/A]`[9]

## G+GE
PollenDAPDaysGGE <- ge4$`PollenDAP(Days)`[9] 
SilkDAPDaysGGE <- ge4$`SilkDAP(Days)`[9] 
PHGGE <- ge4$`PlantHigh(cm)`[9] 
EHGGE <- ge4$`EarHigh(cm)`[9] 
StandCountnumGGE <- ge4$`StandCount(number)`[9] 
StandCountGGE <- ge4$`StandCount%`[9] 
RootLodgingGGE <- ge4$`RootLodging(number)`[9] 
StalkLodgingGGE <- ge4$`StalkLodging(number)`[9] 
MoistureGGE <- ge4$`Moisture%`[9] 
DrymatterGGE <- ge4$`Dry Matter%`[9] 
YieldGGE <- ge4$`GrainYield[bu/A]`[9] 

## GE
PollenDAPDaysGE <- ge5$`PollenDAP(Days)`[9]
SilkDAPDaysGE <- ge5$`SilkDAP(Days)`[9]
PHGE <- ge5$`PlantHigh(cm)`[9] 
EHGE <- ge5$`EarHigh(cm)`[9] 
StandCountnumGE <- ge5$`StandCount(number)`[9]
StandCountGE <- ge5$`StandCount%`[9] 
RootLodgingGE <- ge5$`RootLodging(number)`[9] 
StalkLodgingGE <- ge5$`StalkLodging(number)`[9]
MoistureGE <- ge5$`Moisture%`[9]
DrymatterGE <- ge5$`Dry Matter%`[9]
YieldGE <- ge5$`GrainYield[bu/A]`[9]

## E Effects
PollenDAPDaysE <-  PollenDAPDaysEGGE$ge_mat - PollenDAPDaysGGE$ge_mat
PollenDAPDaysE <- as.data.frame(PollenDAPDaysE)

## G Effects
PollenDAPDaysG <- ge2$`PollenDAP(Days)`[2] - ge$`PollenDAP(Days)`[2] 
SilkDAPDaysG <- ge2$`SilkDAP(Days)`[2] - ge$`SilkDAP(Days)`[2]
PHG <- ge2$`PlantHigh(cm)`[2] - ge$`PlantHigh(cm)`[2]
EHG <- ge2$`EarHigh(cm)`[2] - ge$`EarHigh(cm)`[2]
StandCountnumG <- ge2$`StandCount(number)`[2]-ge$`StandCount(number)`[2]
StandCountG <- ge2$`StandCount%`[2]-ge$`StandCount%`[2]
RootLodgingG <- ge2$`RootLodging(number)`[2]-ge$`RootLodging(number)`[2]
StalkLodgingG <- ge2$`StalkLodging(number)`[2]-ge$`StalkLodging(number)`[2]
MoistureG <- ge2$`Moisture%`[2]-ge$`Moisture%`[2]
DrymatterG <- ge2$`Dry Matter%`[2]-ge$`Dry Matter%`[2]
YieldG <- ge2$`GrainYield[bu/A]`[2]-ge$`GrainYield[bu/A]`[2]

write.table(PollenDAPDaysEGGE,file = "C://R(MoveLater)//Data//eggeYield(buA).txt" ,row.names = FALSE) ##Save E+G+GEI effects
write.table(PollenDAPDaysGGE,file = "C://R(MoveLater)//Data//ggeYield(buA).txt" ,row.names = FALSE) ##Save G+GEI effects
write.table(PollenDAPDaysG ,file = "C://R(MoveLater)//Data//GPollen.txt" ,row.names = FALSE) ##Save G effects seperately
write.table(PollenDAPDaysGE ,file = "C://R(MoveLater)//Data//GGEPollen.txt" ,row.names = FALSE) ##Save GEI effects seperately
write.table(PollenDAPDaysE ,file = "C://R(MoveLater)//Data//EPollen.txt" ,row.names = FALSE) ##Save E effects seperately

#####################################
TestPhenoCov2 <- TestPhenoCov

ge3 <- gge(TestPhenoCov,env = Location, gen = Genotypes, centering = "global")


#write.table(ge$`GrainYield[bu/A]`[2],file = "C://R(MoveLater)//Data//geYield(buA).txt" ,row.names = FALSE) ##save ge effects (It is not a data frame so spesify the ge)

##################################################
library(qtl)

QTL <- read.cross(format = "csvs", dir = "C:/R(MoveLater)/", genfile = "outfileABHTLast3.csv",
                   phefile = "OutGoettingen.csv", crosstype = "dh",
                   map.function = "morgan", estimate.map = FALSE)

#QTL <- calc.genoprob(QTL, step = 2.5)
outnocov <- scanone(QTL, pheno.col = 1, method = "hk") # Scan without covariates
plot(outnocov)
summary(outnocov)

Pollencov <- QTL2$pheno$PollenDapDayscov # read covariates
Silkcov <- QTL2$pheno$SilkDapDayscov # read covariates
PHcov <- QTL2$pheno$PlantHeightcov # read covariates
EHcov <- QTL2$pheno$EarHeightcov # read covariates
Standcountcov <- QTL2$pheno$StandCount.cov # read covariates
Rootlodgingcov <- QTL2$pheno$RootLodgingcov # read covariates
Stalklodgingcov <- QTL2$pheno$StalkLodgingcov # read covariates
Moisturecov <- QTL2$pheno$Moisturecov # read covariates
Yieldcov <- QTL2$pheno$Yieldcov # read covariates

out.acovar <- scanone(QTL, pheno.col=1, addcovar=Pollencov, method = "hk")
plot(out.acovar)
summary(outnocov, threshold=3, format="allpeaks")
summary(out.acovar, threshold=3, format="allpeaks")

plot(outnocov, out.acovar)
plot(outnocov, out.acovar, lodcolumn=1)

out.icov <- scanone(QTL, pheno.col = 1,addcovar = Pollencov ,intcovar = Pollencov, method = "hk")
summary(out.icov, threshold=3, format="allpeaks")
plot(out.acovar, out.icov, col=c("blue","red"))

out.Pollenint <- out.icov - out.acovar
plot(out.Pollenint,outnocov ,lodcolumn=1, col=c("green", "purple"))

#####
seed <- ceiling(runif(1, 0, 10^8))
set.seed(seed)
operm.acovar <- scanone(QTL, pheno.col=1, addcovar=Pollencov,method="hk", n.perm=100)
set.seed(seed)
operm.icovar <- scanone(QTL, pheno.col=1, addcovar=Pollencov,intcovar=Pollencov, method="hk", n.perm=100)

#Again, the differences concern the QTL×E interaction.
operm.Pollenint <- operm.icovar - operm.acovar

#We can use summary to get the genome-wide LOD thresholds.
summary(operm.Pollenint, alpha=c(0.05, 0.20))

#################
perm_hk <- scanone(QTL, pheno.col=1, method= "hk",n.perm=1000)
plot(perm_hk)
summary(perm_hk)  #print the thresholds
#################

plot(out.icov,outnocov, ylim = c(0,5), lwd=c(1,1), col= c("red", "blue"))
abline(h=summary(perm_hk)[[1]],col="green",lty=2,lwd=2)

#We can also use these results to look at evidence for QTL×E interaction in our initial scans.
summary(out.Pollenint, perms=operm.Pollenint, alpha=0.1,format="allpeaks", pvalues=TRUE)
##################

QTL2 <- calc.genoprob(QTL)
summary(outnocov)
qtl <- makeqtl(QTL2, chr=c(1,2,4,5,8,1,3,5,6,10), pos=c(289.2,18.1,243.6,214.4,143.7,115.4484,57.1927,28.0878,95.7534,148), what="prob")
plot(qtl)
summary(fitqtl(QTL2, qtl=qtl, method="hk", get.ests=TRUE, dropone=FALSE, pheno.col = 1))
addint(QTL2, qtl = qtl, method = "hk", pheno.col = 1) ## Look all the possible interactions between QTLs

## Create your model
out.fq <- fitqtl(QTL2, qtl = qtl, method = "hk", pheno.col = 1)
summary(out.fq)
out.fqi <- fitqtl(QTL2, qtl=qtl, method="hk", formula=y ~ Q1+Q2+Q3+Q4+Q5+Q6+Q7+Q8+Q9+Q10+Q1*Q7+Q1*Q8+Q1*Q9+Q1*Q10+Q2*Q7+Q2*Q8+Q2*Q9+Q2*Q10
                  +Q3*Q7+Q3*Q8+Q3*Q9+Q3*Q10+Q4*Q7+Q4*Q8+Q4*Q9+Q4*Q10+Q5*Q7+Q5*Q8+Q5*Q9+Q5*Q10+Q6*Q7+Q6*Q8+Q6*Q9+Q6*Q10, pheno.col = 1)
summary(out.fqi)
out.fqi2 <- fitqtl(QTL, qtl=qtl, method="hk", formula=y ~ Q1+Q2+Q3+Q4+Q6+Q7+Q8+Q9+Q1*Q3+Q2*Q3+Q2*Q5+Q3*Q9+Q4*Q8+Q7*Q8, pheno.col = 3)
summary(out.fqi2)
out.fqi3 <- fitqtl(QTL, qtl=qtl, method="hk", formula=y ~ Q1+Q2+Q3+Q4+Q6+Q7+Q8+Q9+Q2*Q3+Q2*Q5+Q3*Q9+Q4*Q8+Q7*Q8, pheno.col = 3)
summary(out.fqi3)
out.fqi4 <- fitqtl(QTL, qtl=qtl, method="hk", formula=y ~ Q1+Q2+Q3+Q4+Q6+Q7+Q8+Q9+Q2*Q3+Q3*Q9+Q4*Q8+Q7*Q8, pheno.col = 3)
summary(out.fqi4)
out.fqi5 <- fitqtl(QTL, qtl=qtl, method="hk", formula=y ~ Q1+Q2+Q3+Q4+Q6+Q7+Q8+Q9+Q2*Q3+Q3*Q9+Q4*Q8, pheno.col = 3)
summary(out.fqi5)
out.fqi6 <- fitqtl(QTL, qtl=qtl, method="hk", formula=y ~ Q1+Q2+Q3+Q4+Q6+Q7+Q8+Q9+Q2*Q3+Q3*Q9, pheno.col = 3)
summary(out.fqi6)
out.fqi7 <- fitqtl(QTL, qtl=qtl, method="hk", formula=y ~ Q1+Q2+Q3+Q4+Q6+Q7+Q8+Q9+Q3*Q9, pheno.col = 3)
summary(out.fqi7)
out.fqi8 <- fitqtl(QTL, qtl=qtl, method="hk", formula=y ~ Q1+Q2+Q3+Q4+Q6+Q7+Q9+Q3*Q9, pheno.col = 3)
summary(out.fqi8)

##Use the model you have created 

rqtl <- refineqtl(QTL2, qtl = qtl, method = "hk", pheno.col = 1, formula=y ~ Q1+Q2+Q3+Q4+Q5+Q6+Q7+Q8+Q9+Q10+Q1*Q7+Q1*Q8+Q1*Q9+Q1*Q10+Q2*Q7+Q2*Q8+Q2*Q9+Q2*Q10
                  +Q3*Q7+Q3*Q8+Q3*Q9+Q3*Q10+Q4*Q7+Q4*Q8+Q4*Q9+Q4*Q10+Q5*Q7+Q5*Q8+Q5*Q9+Q5*Q10+Q6*Q7+Q6*Q8+Q6*Q9+Q6*Q10)
summary(out.fqr <- fitqtl(QTL2, qtl=rqtl, method="hk", formula=y ~ Q1+Q2+Q3+Q4+Q5+Q6+Q7+Q8+Q9+Q10+Q1*Q7+Q1*Q8+Q1*Q9+Q1*Q10+Q2*Q7+Q2*Q8+Q2*Q9+Q2*Q10
                          +Q3*Q7+Q3*Q8+Q3*Q9+Q3*Q10+Q4*Q7+Q4*Q8+Q4*Q9+Q4*Q10+Q5*Q7+Q5*Q8+Q5*Q9+Q5*Q10+Q6*Q7+Q6*Q8+Q6*Q9+Q6*Q10, pheno.col = 1))
plotLodProfile(rqtl)
plot(outnocov, col="red", add = TRUE)
abline(h=summary(perm_hk)[[1]],col="green",lty=2,lwd=2)
abline(h=summary(perm_hk)[[2]],col="orange",lty=2,lwd=2)

##Look for additional qtls
out.aq <- addqtl(QTL2, qtl = rqtl, method = "hk", pheno.col = 1)
plot(out.aq)

#
saveRDS(QTL2, file = "QTL2.R")
out2 <- scantwo(QTL2, method = "hk", pheno.col = 1,addcovar = Pollencov ,intcovar = Pollencov)
saveRDS(pheno_data,file = "pheno_data.R")
print(pen <- calc.penalties(QTL2))
outpermpollen <- readRDS(file = "/home/uni08/balaca/R/Results/Multiple QTLs/out2permpollenDAPdays.R")

##save your file (I don't know why I saved additional QTLs. Check it out later!)
setwd("C://R(MoveLater)//Results//Multiple QTLs")
save(out.aq, file = "MQTLEarHeight.Rdata")

## Plot and compare other methods (You can skip that part)
#########################################################
plot(result_hk,result_em,result_imp,lwd=c(1,1,1),col=c("red","green","blue"))
legend("topleft",c("hk","em","imp"),lwd=c(1,1,1),col=c("red","green","blue"))

plot(result_hk,result_hk2, lwd=c(3,1),col=c("red","green"), ylim = c(0,5))
abline(h=summary(perm_hk)[[1]],col="green",lty=2,lwd=2)
abline(h=summary(perm_hk)[[2]],col="orange",lty=2,lwd=2)
legend("topleft",c("BLUPs","Goettingen"),lwd=c(1,1),col=c("red","green"), title = "Grain Yield")
save.image()
##########################################################


