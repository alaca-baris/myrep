install.packages("lme4")
install.packages("dplyr")
library(lme4)
library(dplyr)
fix(TestPhenoCov)
TestPhenoCov$Genotypes <- as.factor(TestPhenoCov$Genotypes)
TestPhenoCov$Replication <- as.factor(TestPhenoCov$Replication)
TestPhenoCov$Location <- as.factor(TestPhenoCov$Location)
TestPhenoCov$`Pollen DAP (Days)`<-as.factor(TestPhenoCov$`Pollen DAP (Days)`)


hist(TestPhenoCov$`Pollen DAP (Days)`)
hist(TestPhenoCov$`Plant High (cm)`)

Dataoutput <- data.frame(matrix(vector(), 316,1, dimnames = list(c(), c("Entry"))))

Dataoutput$Entry <- unique(TestPhenoCov[,1])
Dataoutput$Row <- c(1:316)

DataVarComp <- data.frame()
DataVarCompOutput <- data.frame()
Hdata <- data.frame()

drops <- c("var1", "var2", "sdcor")

colnames(TestPhenoCov)

colnum = c(4:ncol(TestPhenoCov))
str(TestPhenoCov[,4:ncol(TestPhenoCov)])

i = 1
for (i in 1:12) {
  x = colnum[i]
  trait = colnames(TestPhenoCov)[x]
  df1 <- TestPhenoCov
  colnames(df1)[x] <- "y"
  
  model <- lmer(y ~ (1|Genotypes)+(1|Location)+(1|Replication)+(1|Genotypes:Location),df1,na.action = "na.omit" )
  
  summary(model)
  varComp <- as.data.frame(VarCorr(model,comp="vcov"))
  blup <- coef(model)$Genotypes
  hist(blup[,1])
  colnames(blup) <- trait
  Dataoutput <- cbind(Dataoutput,blup)
  varComp <- varComp[ , !(names(varComp) %in% drops)]
  varComp$Trait <- trait
  DataVarComp <- rbind(DataVarComp,varComp)
  
}
DataVarCompOutput <- reshape(DataVarComp, idvar = "Trait", timevar = "grp", direction = "wide")

nloc = 316
Hdata <- ((DataVarCompOutput[,3]))/(((DataVarCompOutput[,3])) + (((DataVarCompOutput[,5]))/(nloc)))
