## Use Packages below....

library(snow)
library(doSNOW)
library(parallel)
detectCores()
cl <- makeCluster(6, type = "SOCK")
registerDoSNOW(cl)

Packages <- c("gdata", "readr", "tidyr", "ggthemes", "qqman", "dglm", "hglm",
              "CMplot", "dplyr", "ggplot2", "SNPRelate", "ggcorrplot", "statmod", "pheatmap", 
              "ggfortify")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SNPRelate")
lapply(Packages, library, character.only = TRUE)

##### Read PC values into R
setwd("C://R(MoveLater)//Results//vGWAS")
PCs <- read.table(file = "PCs.txt", header = TRUE)
info <- read.table(file = "Info.txt", header = TRUE)
pca_data_final <- merge(PCs, info, by = "sample.id")

pca <- pca_data_final[c(2:7)]
autoplot(prcomp(pca[, -6]))
autoplot(prcomp(pca[, -6]), data = pca_data_final, colour = "Type") + xlab("PC1") + 
  ylab("PC2") + theme(text = element_text(size = 20, family = "Times")) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) + 
  theme_bw()

## read your geno and positions file into r
geno <- read.csv("C:/R(MoveLater)/Results/vGWAS/outfileABHTLast5.csv", row.names=1)
#positionsss <- transpose(positionsss)
#write.table(positionsss,"Positions.txt", row.names = FALSE, col.names = F)
Positions <- read.csv("C:/R(MoveLater)/Results/vGWAS/Positions.txt", row.names=1, sep="")
map <- Positions
covar <- PCs

##Read your pheno data
pheno_data <- read.delim("C:/R(MoveLater)/Results/vGWAS/PhenoBLUP.txt", row.names=1)


i = i
y <- pheno_data[,11]
model <- dglm(y ~ geno[, i] + covar[, 2] + covar[, 3] + covar[, 4] + covar[, 5], ~geno[, i], family = gaussian(link = identity))
P.mean <- summary(model)$coef[2, 4]  # Extarct p values for mean part
P.disp <- anova(model)$Adj.P[2]  # Extract P values for dispersion part
s.model <- summary(model$dispersion.fit)
beta <- s.model$coef[2, 1]  # Extarct cofficients
se <- s.model$coef[2, 2]  # Extract standard errors
out <- data.frame(Beta = beta, SE = se, P.mean = P.mean, P.disp = P.disp, 
                  stringsAsFactors = FALSE)  # Save all the extracted variables in data frame out



##Create your function for dglm

my.pdglm <- function(cT, i, Phenos, geno, covar) {
  y <- Phenos[,cT]
  model <- dglm(y ~ geno[, i] + covar[, 2] + covar[, 3] + covar[, 4] + covar[, 5], ~geno[, i], family = gaussian(link = "identity"))
  P.mean <- summary(model)$coef[2, 4]  # Extarct p values for mean part
  P.disp <- anova(model)$Adj.P[2]  # Extract P values for dispersion part
  s.model <- summary(model$dispersion.fit)
  beta <- s.model$coef[2, 1]  # Extarct cofficients
  se <- s.model$coef[2, 2]  # Extract standard errors
  out <- data.frame(Beta = beta, SE = se, P.mean = P.mean, P.disp = P.disp, 
                    stringsAsFactors = FALSE)  # Save all the extracted variables in data frame out
  return(out)
  
  # print(i)
}

TF <- matrix(NA, nrow = dim(geno)[2], ncol = 4)
TF <- as.data.frame(TF)
row.names(TF) <- colnames(geno)
colnames(TF) <- c("Beta", "SE", "P.mean", "P.disp")
for (i in 1:dim(geno)[2]) {
  try({
    outm <- my.pdglm(cT = 11, i = i, Phenos = pheno_data, geno = geno, covar = covar)
    TF[i,] <- outm
    print(i)
  }, silent = FALSE)
}

saveRDS(TF, file = "TF")

#### Before applying this step save your progress!!!!

TF <- data.frame(merge(map, TF[, 3:4], by = "row.names", all.x = TRUE))  # add map info
colnames(TF) <- c("marker", "chrom", "pos", "p.mean", "p.disp")

write.csv(TF, file = "DGLMYield.csv", sep = ",", row.names = FALSE)

## Read your p values to R
gwas_cd <- read.delim2("C://R(MoveLater)//Results//vGWAS//DGLM&HGLMoutput//GwasYield.txt")
colnames(gwas_cd) <- c("marker", "chrom", "pos", "gwas_p")
dglm_cd <- read.csv("C://R(MoveLater)//Results//vGWAS//DGLM&HGLMoutput//DGLMYield.csv")
dglm_cd <- dglm_cd[, c(1, 4)]
colnames(dglm_cd) <- c("marker", "dglm_p.disp")
hglm_cd <- read.csv("C://R(MoveLater)//Results//vGWAS//DGLM&HGLMoutput//hglm_cdYield.csv")
hglm_cd <- hglm_cd[, c(1, 4)]
colnames(hglm_cd) <- c("marker", "hglm_p.disp")


all_p <- Reduce(merge, list(gwas_cd, dglm_cd, hglm_cd)) ## Merge them together


CMplot(all_p, plot.type = "c", r = 0.4, col = matrix(c("grey30", "slategray", 
          NA, "lightskyblue1", "lightseagreen", NA, "thistle", "pink", NA), 3, 3, 
          byrow = T), chr.labels = paste("", c("1","2","3","4","5","6","7", "8", "9","10"), sep = ""), cir.legend = TRUE, cir.legend.cex = 0.8, 
       threshold = alpha, outward = FALSE, amplify = TRUE, threshold.lty = c(1,2), signal.line = 1, signal.col = "orangered", threshold.col = "blue", 
       threshold.lwd = 0.5, signal.cex = 1.5, signal.pch = 19, cir.chr.h = 1.5, 
       chr.den.col = "black", file = "jpg", memo = "", dpi = 300)

?CMplot

## HGLM
genotest <- readRDS(file = "geno.rds")
# First upload the marker data
#geno <- readRDS(file = "~/Documents/GitHub/vGWAS/Data/Raw_data/geno.rds")
# Scale the marker data
# Calculate VanRaden
geno3 <- geno
geno2 <- geno
info2 <- info

compute_g_vanraden<-function(geno,lines){
  
  alelleFreq <- function(x, y) {
    (2 * length(which(x == y)) + length(which(x == 1)))/(2 * 
                                                           length(which(!is.na(x))))
  }
  Frequency <- cbind(apply(geno, 2, function(x) alelleFreq(x,0))
                     , apply(geno, 2, function(x) alelleFreq(x, 2)))
  FreqP <- matrix(rep(Frequency[, 2], each = nrow(geno)), 
                  ncol = ncol(geno))
  TwoPQ <- 2 * t(Frequency[, 1]) %*% (Frequency[, 2])
  geno <- geno- 2 * FreqP
  geno[is.na(geno)] <- 0
  Gmatrix <- (tcrossprod(as.matrix(geno)))/as.numeric(TwoPQ)
  
  Gmatrix2 <- cbind(lines,Gmatrix)
  return (Gmatrix2)
}

Gmatrix2 <- compute_g_vanraden(geno = geno2, lines = info2)
Gmatrix2 <- Gmatrix2[,-1]
Gmatrix2 <- as.matrix(Gmatrix2)

###########


Xs <- scale(geno[,-1], center = TRUE, scale = F)
WWl <- Xs %*% t(Xs)
Ga <- WWl/(sum(diag(WWl))/nrow(Xs)) + diag(1e-6, nrow(WWl))
dim(G)
chol.Ga <- chol(Ga)
Z0 <- diag(1, nrow = nrow(Ga), ncol = ncol(Ga))
Z <- Z0 %*% chol.Ga
# Construct G matrix
##Step 1: Second get Cholesky decomposition of G matrix
chol.G <- chol(Gmatrix2)
Z0 <- diag(1, nrow = nrow(Gmatrix2), ncol = ncol(Gmatrix2))
Z <- Z0 %*% chol.G
# Z0 is the identity matrix

?hglm
y <- pheno_data[,1]
i=i
hgl
geno <- geno2
## The basic syntax of HGLM model we are fitting is as follows:
## phenotype=marker_effect (fixed)+ Z (random) (modeling mean), marker_effect(modeling dispersion/variability)
# Run hglm model for all SNPs using for loop
my.hglm <- function(cT, i, Phenos, Z, X, X.disp) {
  y <- Phenos[,cT]
  y2 <- Phenos[,cT]
  outm <- hglm(y = y, X = as.matrix(geno[,i]), Z = Gmatrix2 , X.disp = as.matrix(geno[,i]), family = gaussian(link = "identity"))
  estimates_fix <- outm$fixef
  SE_Mean <- outm$SeFe
  DF <- outm$dfReFe
  DP_Mm <- outm$varFix
  DP_RM <- outm$varRanef
  estimates_rand <- outm$SummVC1[1]
  S.E_rand <- outm$SummVC1[2]
  out <- data.frame(estimates_fix = estimates_fix, SE_Mean = SE_Mean, DF = DF, 
                    DP_Mm = DP_Mm, DP_RM = DP_RM, estimates_rand = estimates_rand, S.E_rand = S.E_rand, 
                    stringsAsFactors = FALSE)
  return(out)
}

# Aanalysis for phenotypes
TF <- matrix(NA, nrow = dim(geno)[2], ncol = 7)

for (i in 1:dim(geno)[2]) {
  try({
    outm <- my.hglm(cT = 1, i = i, Phenos = pheno_data, Z = Gmatrix2, X = geno, X.disp = geno)
    TF[i, ] <- as.numeric(outm)
    print(i)
  }, silent = FALSE)
}

  # read the raw data file
my.hglm1 <- function(pheno, geno1, map1) {
  # add the column names
  colnames(pheno) <- c("estimates_fix", "SE_Mean", "DF", "DP_Mm", "DP_RM", "estimates_rand", 
                       "S.E_rand")
  # add the marker name and position
  markernames <- data.frame(colnames(geno1))
  # Now combine the output file
  pheno <- cbind(markernames, pheno)
  # Now estimate the p-values for mean and dispersion part using library dplyr
  pheno <- mutate(pheno, p.mean = 2 * pt(-abs(estimates_fix/SE_Mean), df = 1), 
                  p.disp = 2 * pt(-abs(estimates_rand/S.E_rand), df = 1))

  # match the markers between map file and outfile
  colnames(pheno) <- c("marker", "estimates_fix", "SE_Mean", "DF", "DP_Mm", "DP_RM", 
                       "estimates_rand", "S.E_rand", "p.mean", "p.disp")
  map <- map
  # now combine the mapfile and outputfile
  pheno <- cbind(map, pheno)
  # now remove the markers with NA values in the file
  pheno <- pheno %>% filter(!is.na(p.mean))
  # now select the appropriate columns for Manhattan plot
  pheno <- select(pheno, marker, Chromosome, Position, p.disp)
  #colnames(pheno) <- c(marker, Chromosome, Position, P)
}

# Now use the function above to extract hglm output and save it.  read the file
# obtained from hglm output
hglm_CD <- TF

# load marker data
geno1 <- geno

# now add the map info to the file for Manhattan plots
map_final <- map
pheno <- hglm_CD
# Now use the function
hglm_CD_final <- my.hglm1(pheno = hglm_CD, geno1 = geno1, map1 = map_final)

# save the file in folder
setwd("C://R(MoveLater)//Results//vGWAS")
write.csv(hglm_CD_final, file = "hglm_cdYield.csv", sep = ",", row.names = FALSE)


# Sample Codes for Epistasis analysis between markers
# Here a linear model is used to check pairwise interactions between markers. 
# Two markers are run at one time for interactions
# for loop is used to run all the pairwise interactions between marker subsets.

# Reading data files
map<-read.csv(file="../data/map.csv")
geno<-read.csv(file="../data/markers.csv")
pheno<-read.csv(file="../data/pheno.final.csv")
pca<-read.csv(file="../data/.csv")
geno2 <- geno[,-1]
geno2 <- geno
pheno <- pheno_data

# Create data.frames to save the output  
TF_M1_M <- matrix(NA,nrow=dim(geno2)[2],ncol=dim(geno2)[2])
rownames(TF_M1_M) <- colnames(TF_M1_M) <- colnames(geno2)
TF_M2_M <- TF_M1.M2_M <- TF_M1_SE <- TF_M2_SE <- TF_M1.M2_SE <- TF_M1_PV <- TF_M2_PV <- TF_M1.M2_PV <- TF_M1_M

# Use for loop to run interaction analysis
for(i in 1:dim(geno2)[2])
{
  for(j in 1:dim(geno2)[2])
  {
    # Run the linear model (lm) 
    model_1<-lm(pheno[,1]~geno2[,i]*geno2[,j]+pca[,2] + pca[, 3] + pca[, 4] + pca[, 5])
    summary(model_1)
    ts <- summary(model_1)$coefficients
    
    i1 <- rownames(ts) %in% 'geno2[, i]'
    i2 <- rownames(ts) %in% 'geno2[, j]'
    i3 <- rownames(ts) %in% 'geno2[, i]:geno2[, j]'
    
    # Coefficients for MAIN marker 1
    if(any(i1))
    {
      p1 = which(i1)
      TF_M1_M[i,j]     <- ts[p1,1] # Extract estimates for marker 1
      TF_M1_SE[i,j]    <- ts[p1,2] # Extract SE  for marker 1
      TF_M1_PV[i,j]    <- ts[p1,4] # Extract p-values  for marker 1
    }
    # Coefficients for MAIN marker 2
    if(any(i2))
    {
      p2 = which(i2)
      TF_M2_M[i,j]     <- ts[p2,1] # Extract estimates  for marker 2
      TF_M2_SE[i,j]    <- ts[p2,2] # Extract SE   for marker 2
      TF_M2_PV[i,j]    <- ts[p2,4] # Extract p-values   for marker 2
    }
    # Coefficients for INTERACTIONs between SNP1*SNP2
    if(any(i3))
    {
      p3 = which(i3)
      TF_M1.M2_M[i,j]  <- ts[p3,1] # Extract estimates  for interactions marker1*marker2
      TF_M1.M2_SE[i,j] <- ts[p3,2] # Extract SE   for interactions marker1*marker2
      TF_M1.M2_PV[i,j] <- ts[p3,4] # Extract p-values  for interactions marker1*marker2
    }
    print(c(i,j))
  }
}

# Save the files
write.csv(TF_M1_M,file='TF_M1_M.csv') # estimates for main marker 1
write.csv(TF_M2_M,file='TF_M2_M.csv') # estimates for for main marker 2
write.csv(TF_M1.M2_M,file='TF_M1.M2_M.csv') # estimates for interactions marker1*marker2

write.csv(TF_M1_SE,file='TF_M1_SE.csv') # standard error for main marker 1
write.csv(TF_M2_SE,file='TF_M2_SE.csv') # standard error for main marker 2
write.csv(TF_M1.M2_SE,file='TF_M1.M2_SE.csv') # standard error for interactions marker1*marker2

write.csv(TF_M1_PV,file='TF_M1_PV.csv') #  p-value  for main effect marker 1
write.csv(TF_M2_PV,file='TF_M2_PV.csv')  # p-value  for main effect marker 2
write.csv(TF_M1.M2_PV,file='TF_M1.M2_PV.csv') # p-values for interactions marker 1* marker 2

####################################################END#################################
