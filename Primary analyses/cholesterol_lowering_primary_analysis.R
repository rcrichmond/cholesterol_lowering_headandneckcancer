#####################
# PRELIMINARY STEPS #
#####################
install.packages("devtools")
devtools::install_github("MRCIEU/MRInstruments", force = TRUE)
devtools::install_github('MRCIEU/TwoSampleMR')
install.packages("MendelianRandomization", force = TRUE)
install.packages("LDlinkR")
install.packages("plyr")
install.packages("ggplot2")
install.packages("ggpubr")
install.packages("simex")
devtools::install_github("rondolab/MR-PRESSO")

library(devtools)
library(TwoSampleMR)
library(MRInstruments)
library(plyr) 
library(ggplot2)
library(MendelianRandomization)
library(gridExtra)
library(grid)
library(lattice)
library(LDlinkR)
library(ggpubr)
library(simex)
library(MRPRESSO)

# Set working directory 
setwd ("")

#1. SNP-EXPOSURE SUMMARY DATA#

drug_exp_data <- read_exposure_data("./drug_exposure_data.csv", sep=",")
head(drug_exp_data)
dim(drug_exp_data)

########################################################
#ORAL AND OROPHARYNGEAL CANCER COMBINED                #
########################################################

#2. SNP-OUTCOME SUMMARY DATA
drug_hnc <- read.csv("./drug_exp_out_hnc.csv", header=TRUE)
drug_hnc <- read_outcome_data("drug_exp_out_hnc.csv", drug_exp_data$SNP, sep=",")

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
dat_hnc <- harmonise_data(drug_exp_data, drug_hnc, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
mr_results_hnc <- mr(dat_hnc)
mr_results_hnc

#Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results_hnc, dat_hnc)
mr_scatter
pdf("drug_exp_hnc_scatter.pdf",  width = 15, height = 20)
ggarrange(mr_scatter$`4bR9XM.At7QtJ`,mr_scatter$`56iCmO.At7QtJ`,mr_scatter$`8WzrD5.At7QtJ`,mr_scatter$`A1eL1U.At7QtJ`,mr_scatter$`F5jola.At7QtJ`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#5. ESTIMATE AND PLOT THE CAUSAL EFFECTS OF THE TRAIT ON THE OUTCOME
# Flip the results (to protective effect of LDL-lowering)
mr_results_hnc$b <- -1*(mr_results_hnc$b)

# Estimate odds ratio and 95% confidence interval
mr_results_hnc$or <- exp(mr_results_hnc$b)
mr_results_hnc$cil <- exp(mr_results_hnc$b-1.96*mr_results_hnc$se)
mr_results_hnc$ciu <- exp(mr_results_hnc$b+1.96*mr_results_hnc$se)

results_hnc<-cbind.data.frame(mr_results_hnc$outcome,mr_results_hnc$exposure,mr_results_hnc$nsnp,mr_results_hnc$method,mr_results_hnc$b,mr_results_hnc$se,mr_results_hnc$pval,mr_results_hnc$or,mr_results_hnc$cil,mr_results_hnc$ciu)

#Export results
write.csv(results_hnc,"./gameon_hnc_results.csv")

#6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(dat_hnc)
res_single$b <- -1*(res_single$b)

mr_forest <- mr_forest_plot(res_single)
mr_forest

pdf("drug_exp_hnc_forest.pdf",  width = 15, height = 20)
ggarrange(mr_forest$`4bR9XM.At7QtJ`,mr_forest$`56iCmO.At7QtJ`,mr_forest$`8WzrD5.At7QtJ`,mr_forest$`A1eL1U.At7QtJ`,mr_forest$`F5jola.At7QtJ`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(dat_hnc)
mr_loo <- mr_leaveoneout_plot(res_loo)

pdf("drug_exp_hnc_loo.pdf",  width = 15, height = 20)
ggarrange(mr_loo$`4bR9XM.At7QtJ`,mr_loo$`56iCmO.At7QtJ`,mr_loo$`8WzrD5.At7QtJ`,mr_loo$`A1eL1U.At7QtJ`,mr_loo$`F5jola.At7QtJ`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#7. Calculate I-squared, r-squared and F-statistics
#I-squared function
Isq <- function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

#Calculate Isq weighted and unweighted
I2<-c()
dat <- harmonise_data(drug_exp_data, drug_hnc, action = 1)
str(dat)

#F-statistic
dat$samplesize.exposure <- 188578
dat$samplesize.outcome <- 6034
dat <- steiger_filtering(dat) 

N = dat$samplesize.exposure[1] #sample size
K = length(dat$SNP) #number of SNPs
total_r2 <- sum(dat$rsq.exposure) 
Fstat <- (N-K-1)/K * total_r2 / (1-total_r2)
total_r2
Fstat

#Rename required columns
dat $BetaXG<-dat $beta.exposure
dat $seBetaXG<-dat $se.exposure
BetaXG   = dat $BetaXG
seBetaXG = dat $seBetaXG 
seBetaYG<-dat $se.outcome

BXG = abs(BetaXG)         
# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F = BXG^2/seBetaXG^2
mF = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted

#Save results
output<-cbind (F, mF, Isq_unweighted, Isq_weighted)
I2<-rbind(I2, output)
colnames(I2) <- c("F", "Exposure", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="regression_dilution_isq_weighted_insertgenename.csv", row.names = FALSE)

#8. SIMEX correction
#Create empty dataframe to store output
simexegger<-c()

#Run simex 
#Rename required columns
dat$BetaXG<-dat$beta.exposure
dat$seBetaXG<-dat$se.exposure
dat$BetaYG<-dat$beta.outcome
dat$seBetaYG<-dat$se.outcome
BetaXG <- dat$BetaXG
BetaYG <- dat$BetaYG
seBetaXG <- dat$seBetaXG
seBetaYG <- dat$seBetaYG

BYG <- BetaYG*sign(BetaXG)#Pre-processing steps to ensure all gene--exposure estimates are positive
BXG <- abs(BetaXG)         

#MR-Egger regression (weighted) 
Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

#MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

#Simulation extrapolation 
mod.sim1 <- simex(Fit1,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod.sim2 <- simex(Fit2,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)
#Print results
mod1
mod2


#8. ASSESS HETEROGENEITY AND PLEIOTROPY 

#Hetereogenity test
mr_heterogeneity(dat_hnc)
write.csv(mr_heterogeneity(dat_hnc),"./heterogeneity_results_hnc.csv")

#Pleiotropy test
mr_pleiotropy_test(dat_hnc)
write.csv(mr_pleiotropy_test(dat_hnc),"./pleiotropy_results_hnc.csv")

# MR PRESSO
# Run MR-PRESSO global method
HMGCR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat_hnc[dat_hnc$exposure=="HMGCR"&dat_hnc$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
NPC1L1_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat_hnc[dat_hnc$exposure=="NPC1L1"&dat_hnc$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
PCSK9_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat_hnc[dat_hnc$exposure=="PCSK9"&dat_hnc$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LDLR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat_hnc[dat_hnc$exposure=="LDLR"&dat_hnc$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
CETP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat_hnc[dat_hnc$exposure=="CETP"&dat_hnc$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
#NB Not enough intrumental variables for LDLR 

headers<-c("outcome", "exposure", "RSSobs","Pvalue")
mr_presso_global <- as.data.frame(matrix(,ncol=4,nrow=5))
names(mr_presso_global)<-headers  
mr_presso_global$outcome <- "HNC" 
mr_presso_global$exposure <- c("HMGCR", "NPC1L1", "PCSK9", "LDLR", "CETP")

mr_presso_global$RSSobs[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$Pvalue

write.csv(mr_presso_global, "./mr_presso_hnc_global.csv")

#9. ACCOUNT FOR LD STRUCTURE

dat_hnc <- dat_hnc[dat_hnc$mr_keep=="TRUE",]

HMGCR <- dat_hnc[dat_hnc$exposure=="HMGCR",]
LDmatrix(HMGCR$SNP, pop="CEU", r2d = "r2", token = "b0a2f3f3c1f4", file = "HMGCR_cor.txt")
NPC1L1 <- dat_hnc[dat_hnc$exposure=="NPC1L1",]
LDmatrix(NPC1L1$SNP, pop="CEU", r2d = "r2", token = "b0a2f3f3c1f4", file = "NPC1L1_cor.txt")
PCSK9 <- dat_hnc[dat_hnc$exposure=="PCSK9",]
LDmatrix(PCSK9$SNP, pop="CEU", r2d = "r2", token = "b0a2f3f3c1f4", file = "PCSK9_cor.txt")
LDLR <- dat_hnc[dat_hnc$exposure=="LDLR",]
LDmatrix(LDLR$SNP, pop="CEU", r2d = "r2", token = "b0a2f3f3c1f4", file = "LDLR_cor.txt")
CETP <- dat_hnc[dat_hnc$exposure=="CETP",]
LDmatrix(CETP$SNP, pop="CEU", r2d = "r2", token = "b0a2f3f3c1f4", file = "CETP_cor.txt")

correl <- read.table("HMGCR_cor.txt", header=T)
correl <- correl[order(correl$RS_number),]
correl <- correl[,order(names(correl))]
HMGCR <- HMGCR[HMGCR$SNP %in% names(correl),]
correl <- within(correl, rm(RS_number))
colnames(correl) <- NULL
rownames(correl) <- NULL
correl <- data.matrix(correl)
mr_ivw_HMGCR <- mr_ivw(mr_input(bx = HMGCR$beta.exposure, bxse = HMGCR$se.exposure, by = HMGCR$beta.outcome, byse = HMGCR$se.outcome, correlation = correl))
mr_egger_HMGCR <- mr_egger(mr_input(bx = HMGCR$beta.exposure, bxse = HMGCR$se.exposure, by = HMGCR$beta.outcome, byse = HMGCR$se.outcome, correlation = correl))
mr_median_HMGCR <- mr_median(mr_input(bx = HMGCR$beta.exposure, bxse = HMGCR$se.exposure, by = HMGCR$beta.outcome, byse = HMGCR$se.outcome, correlation = correl))

correl <- read.table("NPC1L1_cor.txt", header=T)
correl <- correl[order(correl$RS_number),]
correl <- correl[,order(names(correl))]
NPC1L1 <- NPC1L1[NPC1L1$SNP %in% names(correl),]
correl <- within(correl, rm(RS_number))
colnames(correl) <- NULL
rownames(correl) <- NULL
correl <- data.matrix(correl)
mr_ivw_NPC1L1 <- mr_ivw(mr_input(bx = NPC1L1$beta.exposure, bxse = NPC1L1$se.exposure, by = NPC1L1$beta.outcome, byse = NPC1L1$se.outcome, correlation = correl))
mr_egger_NPC1L1 <- mr_egger(mr_input(bx = NPC1L1$beta.exposure, bxse = NPC1L1$se.exposure, by = NPC1L1$beta.outcome, byse = NPC1L1$se.outcome, correlation = correl))
mr_median_NPC1L1 <- mr_median(mr_input(bx = NPC1L1$beta.exposure, bxse = NPC1L1$se.exposure, by = NPC1L1$beta.outcome, byse = NPC1L1$se.outcome, correlation = correl))

correl <- read.table("PCSK9_cor.txt", header=T)
correl <- correl[order(correl$RS_number),]
correl <- correl[,order(names(correl))]
PCSK9 <- PCSK9[PCSK9$SNP %in% names(correl),]
correl <- within(correl, rm(RS_number))
colnames(correl) <- NULL
rownames(correl) <- NULL
correl <- data.matrix(correl)
mr_ivw_PCSK9 <- mr_ivw(mr_input(bx = PCSK9$beta.exposure, bxse = PCSK9$se.exposure, by = PCSK9$beta.outcome, byse = PCSK9$se.outcome, correlation = correl))
mr_egger_PCSK9 <- mr_egger(mr_input(bx = PCSK9$beta.exposure, bxse = PCSK9$se.exposure, by = PCSK9$beta.outcome, byse = PCSK9$se.outcome, correlation = correl))
mr_median_PCSK9 <- mr_median(mr_input(bx = PCSK9$beta.exposure, bxse = PCSK9$se.exposure, by = PCSK9$beta.outcome, byse = PCSK9$se.outcome, correlation = correl))

correl <- read.table("LDLR_cor.txt", header=T)
correl <- correl[order(correl$RS_number),]
correl <- correl[,order(names(correl))]
LDLR <- LDLR[LDLR$SNP %in% names(correl),]
correl <- within(correl, rm(RS_number))
colnames(correl) <- NULL
rownames(correl) <- NULL
correl <- data.matrix(correl)
mr_ivw_LDLR <- mr_ivw(mr_input(bx = LDLR$beta.exposure, bxse = LDLR$se.exposure, by = LDLR$beta.outcome, byse = LDLR$se.outcome, correlation = correl))
mr_egger_LDLR <- mr_egger(mr_input(bx = LDLR$beta.exposure, bxse = LDLR$se.exposure, by = LDLR$beta.outcome, byse = LDLR$se.outcome, correlation = correl))
mr_median_LDLR <- mr_median(mr_input(bx = LDLR$beta.exposure, bxse = LDLR$se.exposure, by = LDLR$beta.outcome, byse = LDLR$se.outcome, correlation = correl))

correl <- read.table("CETP_cor.txt", header=T)
correl <- correl[order(correl$RS_number),]
correl <- correl[,order(names(correl))]
CETP <- CETP[CETP$SNP %in% names(correl),]
correl <- within(correl, rm(RS_number))
colnames(correl) <- NULL
rownames(correl) <- NULL
correl <- data.matrix(correl)
mr_ivw_CETP <- mr_ivw(mr_input(bx = CETP$beta.exposure, bxse = CETP$se.exposure, by = CETP$beta.outcome, byse = CETP$se.outcome, correlation = correl))
mr_egger_CETP<- mr_egger(mr_input(bx = CETP$beta.exposure, bxse = CETP$se.exposure, by = CETP$beta.outcome, byse = CETP$se.outcome, correlation = correl))
mr_median_CETP <- mr_median(mr_input(bx = CETP$beta.exposure, bxse = CETP$se.exposure, by = CETP$beta.outcome, byse = CETP$se.outcome, correlation = correl))

# Flip the results (to protective effect of LDL-lowering)
headers<-c("outcome","exposure","method","b","se","pval")
mr_results_hnc_corr <- as.data.frame(matrix(,ncol=6,nrow=15))
names(mr_results_hnc_corr)<-headers

mr_results_hnc_corr$outcome <- "HNC" 
mr_results_hnc_corr$exposure <- c("HMGCR", "HMGCR", "HMGCR", "NPC1L1", "NPC1L1", "NPC1L1", "PCSK9", "PCSK9", "PCSK9", "LDLR", "LDLR", "LDLR", "CETP", "CETP", "CETP")
mr_results_hnc_corr$method <- c("ivw", "egger", "median", "ivw", "egger", "median", "ivw", "egger", "median", "ivw", "egger", "median","ivw", "egger", "median")

mr_results_hnc_corr$b[1] <- -1*(mr_ivw_HMGCR$Estimate)
mr_results_hnc_corr$b[2] <- -1*(mr_egger_HMGCR$Estimate)
mr_results_hnc_corr$b[3] <- -1*(mr_median_HMGCR$Estimate)
mr_results_hnc_corr$b[4] <- -1*(mr_ivw_NPC1L1$Estimate)
mr_results_hnc_corr$b[5] <- -1*(mr_egger_NPC1L1$Estimate)
mr_results_hnc_corr$b[6] <- -1*(mr_median_NPC1L1$Estimate)
mr_results_hnc_corr$b[7] <- -1*(mr_ivw_PCSK9$Estimate)
mr_results_hnc_corr$b[8] <- -1*(mr_egger_PCSK9$Estimate)
mr_results_hnc_corr$b[9] <- -1*(mr_median_PCSK9$Estimate)
mr_results_hnc_corr$b[10] <- -1*(mr_ivw_LDLR$Estimate)
mr_results_hnc_corr$b[11] <- -1*(mr_egger_LDLR$Estimate)
mr_results_hnc_corr$b[12] <- -1*(mr_median_LDLR$Estimate)
mr_results_hnc_corr$b[13] <- -1*(mr_ivw_CETP$Estimate)
mr_results_hnc_corr$b[14] <- -1*(mr_egger_CETP$Estimate)
mr_results_hnc_corr$b[15] <- -1*(mr_median_CETP$Estimate)

mr_results_hnc_corr$se[1] <- mr_ivw_HMGCR$StdError
mr_results_hnc_corr$se[2] <- mr_egger_HMGCR$StdError.Est
mr_results_hnc_corr$se[3] <- mr_median_HMGCR$StdError
mr_results_hnc_corr$se[4] <- mr_ivw_NPC1L1$StdError
mr_results_hnc_corr$se[5] <- mr_egger_NPC1L1$StdError.Est
mr_results_hnc_corr$se[6] <- mr_median_NPC1L1$StdError
mr_results_hnc_corr$se[7] <- mr_ivw_PCSK9$StdError
mr_results_hnc_corr$se[8] <- mr_egger_PCSK9$StdError.Est
mr_results_hnc_corr$se[9] <- mr_median_PCSK9$StdError
mr_results_hnc_corr$se[10] <- mr_ivw_LDLR$StdError
mr_results_hnc_corr$se[11] <- mr_egger_LDLR$StdError.Est
mr_results_hnc_corr$se[12] <- mr_median_LDLR$StdError
mr_results_hnc_corr$se[13] <- mr_ivw_CETP$StdError
mr_results_hnc_corr$se[14] <- mr_egger_CETP$StdError.Est
mr_results_hnc_corr$se[15] <- mr_median_CETP$StdError

mr_results_hnc_corr$pval[1] <- mr_ivw_HMGCR$Pvalue
mr_results_hnc_corr$pval[2] <- mr_egger_HMGCR$Pvalue.Est
mr_results_hnc_corr$pval[3] <- mr_median_HMGCR$Pvalue
mr_results_hnc_corr$pval[4] <- mr_ivw_NPC1L1$Pvalue
mr_results_hnc_corr$pval[5] <- mr_egger_NPC1L1$Pvalue.Est
mr_results_hnc_corr$pval[6] <- mr_median_NPC1L1$Pvalue
mr_results_hnc_corr$pval[7] <- mr_ivw_PCSK9$Pvalue
mr_results_hnc_corr$pval[8] <- mr_egger_PCSK9$Pvalue.Est
mr_results_hnc_corr$pval[9] <- mr_median_PCSK9$Pvalue
mr_results_hnc_corr$pval[10] <- mr_ivw_LDLR$Pvalue
mr_results_hnc_corr$pval[11] <- mr_egger_LDLR$Pvalue.Est
mr_results_hnc_corr$pval[12] <- mr_median_LDLR$Pvalue
mr_results_hnc_corr$pval[13] <- mr_ivw_CETP$Pvalue
mr_results_hnc_corr$pval[14] <- mr_egger_CETP$Pvalue.Est
mr_results_hnc_corr$pval[15] <- mr_median_CETP$Pvalue

# Estimate odds ratio and 95% confidence interval
mr_results_hnc_corr$or <- exp(mr_results_hnc_corr$b)
mr_results_hnc_corr$cil <- exp(mr_results_hnc_corr$b-1.96*mr_results_hnc_corr$se)
mr_results_hnc_corr$ciu <- exp(mr_results_hnc_corr$b+1.96*mr_results_hnc_corr$se)

write.csv(mr_results_hnc_corr, "./hnc_results_corr.csv")

#Hetereogenity test
headers<-c("outcome","exposure","method","Q","Q_df", "Q_pval")
mr_heterogenity_test_corr <- as.data.frame(matrix(,ncol=6,nrow=10))
names(mr_heterogenity_test_corr)<-headers
mr_heterogenity_test_corr$outcome <- "HNC" 
mr_heterogenity_test_corr$exposure <- c("HMGCR", "HMGCR", "NPC1L1", "NPC1L1", "PCSK9", "PCSK9", "LDLR", "LDLR", "CETP", "CETP")
mr_heterogenity_test_corr$method <- c("IVW", "MR Egger", "IVW", "MR Egger", "IVW", "MR Egger", "IVW", "MR Egger", "IVW", "MR Egger")

mr_heterogenity_test_corr$Q[1] <- mr_ivw_HMGCR$Heter.Stat[1]
mr_heterogenity_test_corr$Q[2] <- mr_egger_HMGCR$Heter.Stat[1]
mr_heterogenity_test_corr$Q[3] <- mr_ivw_NPC1L1$Heter.Stat[1]
mr_heterogenity_test_corr$Q[4] <- mr_egger_NPC1L1$Heter.Stat[1]
mr_heterogenity_test_corr$Q[5] <- mr_ivw_PCSK9$Heter.Stat[1]
mr_heterogenity_test_corr$Q[6] <- mr_egger_PCSK9$Heter.Stat[1]
mr_heterogenity_test_corr$Q[7] <- mr_ivw_LDLR$Heter.Stat[1]
mr_heterogenity_test_corr$Q[8] <- mr_egger_LDLR$Heter.Stat[1]
mr_heterogenity_test_corr$Q[9] <- mr_ivw_CETP$Heter.Stat[1]
mr_heterogenity_test_corr$Q[10] <- mr_egger_CETP$Heter.Stat[1]

mr_heterogenity_test_corr$Q_df[1] <- mr_ivw_HMGCR$SNPs - 1
mr_heterogenity_test_corr$Q_df[2] <- mr_egger_HMGCR$SNPs - 2
mr_heterogenity_test_corr$Q_df[3] <- mr_ivw_NPC1L1$SNPs - 1
mr_heterogenity_test_corr$Q_df[4] <- mr_egger_NPC1L1$SNPs - 2
mr_heterogenity_test_corr$Q_df[5] <- mr_ivw_PCSK9$SNPs - 1
mr_heterogenity_test_corr$Q_df[6] <- mr_egger_PCSK9$SNPs - 2
mr_heterogenity_test_corr$Q_df[7] <- mr_ivw_LDLR$SNPs - 1
mr_heterogenity_test_corr$Q_df[8] <- mr_egger_LDLR$SNPs - 2 
mr_heterogenity_test_corr$Q_df[9] <- mr_ivw_CETP$SNPs - 1
mr_heterogenity_test_corr$Q_df[10] <- mr_egger_CETP$SNPs - 2 

mr_heterogenity_test_corr$Q_pval[1] <- mr_ivw_HMGCR$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[2] <- mr_egger_HMGCR$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[3] <- mr_ivw_NPC1L1$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[4] <- mr_egger_NPC1L1$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[5] <- mr_ivw_PCSK9$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[6] <- mr_egger_PCSK9$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[7] <- mr_ivw_LDLR$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[8] <- mr_egger_LDLR$Heter.Stat[2] 
mr_heterogenity_test_corr$Q_pval[9] <- mr_ivw_CETP$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[10] <- mr_egger_CETP$Heter.Stat[2] 

write.csv(mr_heterogenity_test_corr,"./heterogeneity_results_hnc_corr.csv")

# Egger intercept 
headers<-c("outcome","exposure","egger_intercept","se","pval")
mr_pleiotropy_test_corr <- as.data.frame(matrix(,ncol=5,nrow=5))
names(mr_pleiotropy_test_corr)<-headers
mr_pleiotropy_test_corr$outcome <- "HNC" 
mr_pleiotropy_test_corr$exposure <- c("HMGCR", "NPC1L1","PCSK9","LDLR","CETP")
mr_pleiotropy_test_corr$egger_intercept[1] <- -1*(mr_egger_HMGCR$Intercept)
mr_pleiotropy_test_corr$egger_intercept[2] <- -1*(mr_egger_NPC1L1$Intercept)
mr_pleiotropy_test_corr$egger_intercept[3] <- -1*(mr_egger_PCSK9$Intercept)
mr_pleiotropy_test_corr$egger_intercept[4] <- -1*(mr_egger_LDLR$Intercept)
mr_pleiotropy_test_corr$egger_intercept[5] <- -1*(mr_egger_CETP$Intercept)
mr_pleiotropy_test_corr$se[1] <- mr_egger_HMGCR$StdError.Int
mr_pleiotropy_test_corr$se[2] <- mr_egger_NPC1L1$StdError.Int
mr_pleiotropy_test_corr$se[3] <- mr_egger_PCSK9$StdError.Int
mr_pleiotropy_test_corr$se[4] <- mr_egger_LDLR$StdError.Int
mr_pleiotropy_test_corr$se[5] <- mr_egger_CETP$StdError.Int
mr_pleiotropy_test_corr$pval[1] <- mr_egger_HMGCR$Pvalue.Int
mr_pleiotropy_test_corr$pval[2] <- mr_egger_NPC1L1$Pvalue.Int
mr_pleiotropy_test_corr$pval[3] <- mr_egger_PCSK9$Pvalue.Int
mr_pleiotropy_test_corr$pval[4] <- mr_egger_LDLR$Pvalue.Int
mr_pleiotropy_test_corr$se[5] <- mr_egger_CETP$StdError.Int

write.csv(mr_pleiotropy_test_corr,"./pleiotropy_results_hnc_corr.csv")

########################################################
# ORAL CANCER                                          #
########################################################

#2. SNP-OUTCOME SUMMARY DATA
drug_hnc <- read.csv("./drug_exp_out_oc.csv", header=TRUE)
drug_hnc <- read_outcome_data("drug_exp_out_oc.csv", drug_exp_data$SNP, sep=",")

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
dat_hnc <- harmonise_data(drug_exp_data, drug_hnc, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
mr_results_hnc <- mr(dat_hnc)
mr_results_hnc

#Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results_hnc, dat_hnc)
mr_scatter
pdf("drug_exp_oc_scatter.pdf",  width = 15, height = 20)
ggarrange(mr_scatter$`D1lrdN.vWEiKn`,mr_scatter$`j9cXa6.vWEiKn`,mr_scatter$`W6Ai9u.vWEiKn`,mr_scatter$`ZL81yn.vWEiKn`,mr_scatter$`zzRhsx.vWEiKn`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#5. ESTIMATE AND PLOT THE CAUSAL EFFECTS OF THE TRAIT ON THE OUTCOME
# Flip the results (to protective effect of LDL-lowering)
mr_results_hnc$b <- -1*(mr_results_hnc$b)

# Estimate odds ratio and 95% confidence interval
mr_results_hnc$or <- exp(mr_results_hnc$b)
mr_results_hnc$cil <- exp(mr_results_hnc$b-1.96*mr_results_hnc$se)
mr_results_hnc$ciu <- exp(mr_results_hnc$b+1.96*mr_results_hnc$se)

results_hnc<-cbind.data.frame(mr_results_hnc$outcome,mr_results_hnc$exposure,mr_results_hnc$nsnp,mr_results_hnc$method,mr_results_hnc$b,mr_results_hnc$se,mr_results_hnc$pval,mr_results_hnc$or,mr_results_hnc$cil,mr_results_hnc$ciu)

#Export results
write.csv(results_hnc,"./gameon_oc_results.csv")

#6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(dat_hnc)
res_single$b <- -1*(res_single$b)

mr_forest <- mr_forest_plot(res_single)
mr_forest

pdf("drug_exp_oc_forest.pdf",  width = 15, height = 20)
ggarrange(mr_forest$`4bR9XM.At7QtJ`,mr_forest$`56iCmO.At7QtJ`,mr_forest$`8WzrD5.At7QtJ`,mr_forest$`A1eL1U.At7QtJ`,mr_forest$`F5jola.At7QtJ`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(dat_hnc)
mr_loo <- mr_leaveoneout_plot(res_loo)

pdf("drug_exp_oc_loo.pdf",  width = 15, height = 20)
ggarrange(mr_loo$`4bR9XM.At7QtJ`,mr_loo$`56iCmO.At7QtJ`,mr_loo$`8WzrD5.At7QtJ`,mr_loo$`A1eL1U.At7QtJ`,mr_loo$`F5jola.At7QtJ`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#7. Calculate I-squared, r-squared and F-statistics
#I-squared function
Isq <- function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

#Calculate Isq weighted and unweighted
I2<-c()
dat <- harmonise_data(drug_exp_data, drug_hnc, action = 1)
str(dat)

#F-statistic
dat$samplesize.exposure <- 188578
dat$samplesize.outcome <- 6034
dat <- steiger_filtering(dat) 

N = dat$samplesize.exposure[1] #sample size
K = length(dat$SNP) #number of SNPs
total_r2 <- sum(dat$rsq.exposure) 
Fstat <- (N-K-1)/K * total_r2 / (1-total_r2)
total_r2
Fstat

#Rename required columns
dat $BetaXG<-dat $beta.exposure
dat $seBetaXG<-dat $se.exposure
BetaXG   = dat $BetaXG
seBetaXG = dat $seBetaXG 
seBetaYG<-dat $se.outcome

BXG = abs(BetaXG)         
# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F = BXG^2/seBetaXG^2
mF = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted

#Save results
output<-cbind (F, mF, Isq_unweighted, Isq_weighted)
I2<-rbind(I2, output)
colnames(I2) <- c("F", "Exposure", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="regression_dilution_isq_weighted_insertgenename.csv", row.names = FALSE)

#8. SIMEX correction
#Create empty dataframe to store output
simexegger<-c()

#Run simex 
#Rename required columns
dat$BetaXG<-dat$beta.exposure
dat$seBetaXG<-dat$se.exposure
dat$BetaYG<-dat$beta.outcome
dat$seBetaYG<-dat$se.outcome
BetaXG <- dat$BetaXG
BetaYG <- dat$BetaYG
seBetaXG <- dat$seBetaXG
seBetaYG <- dat$seBetaYG

BYG <- BetaYG*sign(BetaXG)#Pre-processing steps to ensure all gene--exposure estimates are positive
BXG <- abs(BetaXG)         

#MR-Egger regression (weighted) 
Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

#MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

#Simulation extrapolation 
mod.sim1 <- simex(Fit1,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod.sim2 <- simex(Fit2,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)
#Print results
mod1
mod2

#8. ASSESS HETEROGENEITY AND PLEIOTROPY 

#Hetereogenity test
mr_heterogeneity(dat_hnc)
write.csv(mr_heterogeneity(dat_hnc),"./heterogeneity_results_oc.csv")

#Hetereogenity test
mr_pleiotropy_test(dat_hnc)
write.csv(mr_pleiotropy_test(dat_hnc),"./pleiotropy_results_oc.csv")

# MR PRESSO
# Run MR-PRESSO global method
HMGCR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat_hnc[dat_hnc$exposure=="HMGCR"&dat_hnc$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
NPC1L1_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat_hnc[dat_hnc$exposure=="NPC1L1"&dat_hnc$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
PCSK9_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat_hnc[dat_hnc$exposure=="PCSK9"&dat_hnc$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LDLR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat_hnc[dat_hnc$exposure=="LDLR"&dat_hnc$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
CETP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat_hnc[dat_hnc$exposure=="CETP"&dat_hnc$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
#NB Not enough intrumental variables for LDLR 

headers<-c("outcome", "exposure", "RSSobs","Pvalue")
mr_presso_global <- as.data.frame(matrix(,ncol=4,nrow=5))
names(mr_presso_global)<-headers  
mr_presso_global$outcome <- "HNC" 
mr_presso_global$exposure <- c("HMGCR", "NPC1L1", "PCSK9", "LDLR", "CETP")

mr_presso_global$RSSobs[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$Pvalue

write.csv(mr_presso_global, "./mr_presso_oc_global.csv")

#9. ACCOUNT FOR LD STRUCTURE

dat_hnc <- dat_hnc[dat_hnc$mr_keep=="TRUE",]

HMGCR <- dat_hnc[dat_hnc$exposure=="HMGCR",]
LDmatrix(HMGCR$SNP, pop="CEU", r2d = "r2", token = "b0a2f3f3c1f4", file = "HMGCR_cor.txt")
NPC1L1 <- dat_hnc[dat_hnc$exposure=="NPC1L1",]
LDmatrix(NPC1L1$SNP, pop="CEU", r2d = "r2", token = "b0a2f3f3c1f4", file = "NPC1L1_cor.txt")
PCSK9 <- dat_hnc[dat_hnc$exposure=="PCSK9",]
LDmatrix(PCSK9$SNP, pop="CEU", r2d = "r2", token = "b0a2f3f3c1f4", file = "PCSK9_cor.txt")
LDLR <- dat_hnc[dat_hnc$exposure=="LDLR",]
LDmatrix(LDLR$SNP, pop="CEU", r2d = "r2", token = "b0a2f3f3c1f4", file = "LDLR_cor.txt")
CETP <- dat_hnc[dat_hnc$exposure=="CETP",]
LDmatrix(CETP$SNP, pop="CEU", r2d = "r2", token = "b0a2f3f3c1f4", file = "CETP_cor.txt")

correl <- read.table("HMGCR_cor.txt", header=T)
correl <- correl[order(correl$RS_number),]
correl <- correl[,order(names(correl))]
HMGCR <- HMGCR[HMGCR$SNP %in% names(correl),]
correl <- within(correl, rm(RS_number))
colnames(correl) <- NULL
rownames(correl) <- NULL
correl <- data.matrix(correl)
mr_ivw_HMGCR <- mr_ivw(mr_input(bx = HMGCR$beta.exposure, bxse = HMGCR$se.exposure, by = HMGCR$beta.outcome, byse = HMGCR$se.outcome, correlation = correl))
mr_egger_HMGCR <- mr_egger(mr_input(bx = HMGCR$beta.exposure, bxse = HMGCR$se.exposure, by = HMGCR$beta.outcome, byse = HMGCR$se.outcome, correlation = correl))
mr_median_HMGCR <- mr_median(mr_input(bx = HMGCR$beta.exposure, bxse = HMGCR$se.exposure, by = HMGCR$beta.outcome, byse = HMGCR$se.outcome, correlation = correl))

correl <- read.table("NPC1L1_cor.txt", header=T)
correl <- correl[order(correl$RS_number),]
correl <- correl[,order(names(correl))]
NPC1L1 <- NPC1L1[NPC1L1$SNP %in% names(correl),]
correl <- within(correl, rm(RS_number))
colnames(correl) <- NULL
rownames(correl) <- NULL
correl <- data.matrix(correl)
mr_ivw_NPC1L1 <- mr_ivw(mr_input(bx = NPC1L1$beta.exposure, bxse = NPC1L1$se.exposure, by = NPC1L1$beta.outcome, byse = NPC1L1$se.outcome, correlation = correl))
mr_egger_NPC1L1 <- mr_egger(mr_input(bx = NPC1L1$beta.exposure, bxse = NPC1L1$se.exposure, by = NPC1L1$beta.outcome, byse = NPC1L1$se.outcome, correlation = correl))
mr_median_NPC1L1 <- mr_median(mr_input(bx = NPC1L1$beta.exposure, bxse = NPC1L1$se.exposure, by = NPC1L1$beta.outcome, byse = NPC1L1$se.outcome, correlation = correl))

correl <- read.table("PCSK9_cor.txt", header=T)
correl <- correl[order(correl$RS_number),]
correl <- correl[,order(names(correl))]
PCSK9 <- PCSK9[PCSK9$SNP %in% names(correl),]
correl <- within(correl, rm(RS_number))
colnames(correl) <- NULL
rownames(correl) <- NULL
correl <- data.matrix(correl)
mr_ivw_PCSK9 <- mr_ivw(mr_input(bx = PCSK9$beta.exposure, bxse = PCSK9$se.exposure, by = PCSK9$beta.outcome, byse = PCSK9$se.outcome, correlation = correl))
mr_egger_PCSK9 <- mr_egger(mr_input(bx = PCSK9$beta.exposure, bxse = PCSK9$se.exposure, by = PCSK9$beta.outcome, byse = PCSK9$se.outcome, correlation = correl))
mr_median_PCSK9 <- mr_median(mr_input(bx = PCSK9$beta.exposure, bxse = PCSK9$se.exposure, by = PCSK9$beta.outcome, byse = PCSK9$se.outcome, correlation = correl))

correl <- read.table("LDLR_cor.txt", header=T)
correl <- correl[order(correl$RS_number),]
correl <- correl[,order(names(correl))]
LDLR <- LDLR[LDLR$SNP %in% names(correl),]
correl <- within(correl, rm(RS_number))
colnames(correl) <- NULL
rownames(correl) <- NULL
correl <- data.matrix(correl)
mr_ivw_LDLR <- mr_ivw(mr_input(bx = LDLR$beta.exposure, bxse = LDLR$se.exposure, by = LDLR$beta.outcome, byse = LDLR$se.outcome, correlation = correl))
mr_egger_LDLR <- mr_egger(mr_input(bx = LDLR$beta.exposure, bxse = LDLR$se.exposure, by = LDLR$beta.outcome, byse = LDLR$se.outcome, correlation = correl))
mr_median_LDLR <- mr_median(mr_input(bx = LDLR$beta.exposure, bxse = LDLR$se.exposure, by = LDLR$beta.outcome, byse = LDLR$se.outcome, correlation = correl))

correl <- read.table("CETP_cor.txt", header=T)
correl <- correl[order(correl$RS_number),]
correl <- correl[,order(names(correl))]
CETP <- CETP[CETP$SNP %in% names(correl),]
correl <- within(correl, rm(RS_number))
colnames(correl) <- NULL
rownames(correl) <- NULL
correl <- data.matrix(correl)
mr_ivw_CETP <- mr_ivw(mr_input(bx = CETP$beta.exposure, bxse = CETP$se.exposure, by = CETP$beta.outcome, byse = CETP$se.outcome, correlation = correl))
mr_egger_CETP<- mr_egger(mr_input(bx = CETP$beta.exposure, bxse = CETP$se.exposure, by = CETP$beta.outcome, byse = CETP$se.outcome, correlation = correl))
mr_median_CETP <- mr_median(mr_input(bx = CETP$beta.exposure, bxse = CETP$se.exposure, by = CETP$beta.outcome, byse = CETP$se.outcome, correlation = correl))

# Flip the results (to protective effect of LDL-lowering)
headers<-c("outcome","exposure","method","b","se","pval")
mr_results_hnc_corr <- as.data.frame(matrix(,ncol=6,nrow=15))
names(mr_results_hnc_corr)<-headers

mr_results_hnc_corr$outcome <- "HNC" 
mr_results_hnc_corr$exposure <- c("HMGCR", "HMGCR", "HMGCR", "NPC1L1", "NPC1L1", "NPC1L1", "PCSK9", "PCSK9", "PCSK9", "LDLR", "LDLR", "LDLR", "CETP", "CETP", "CETP")
mr_results_hnc_corr$method <- c("ivw", "egger", "median", "ivw", "egger", "median", "ivw", "egger", "median", "ivw", "egger", "median","ivw", "egger", "median")

mr_results_hnc_corr$b[1] <- -1*(mr_ivw_HMGCR$Estimate)
mr_results_hnc_corr$b[2] <- -1*(mr_egger_HMGCR$Estimate)
mr_results_hnc_corr$b[3] <- -1*(mr_median_HMGCR$Estimate)
mr_results_hnc_corr$b[4] <- -1*(mr_ivw_NPC1L1$Estimate)
mr_results_hnc_corr$b[5] <- -1*(mr_egger_NPC1L1$Estimate)
mr_results_hnc_corr$b[6] <- -1*(mr_median_NPC1L1$Estimate)
mr_results_hnc_corr$b[7] <- -1*(mr_ivw_PCSK9$Estimate)
mr_results_hnc_corr$b[8] <- -1*(mr_egger_PCSK9$Estimate)
mr_results_hnc_corr$b[9] <- -1*(mr_median_PCSK9$Estimate)
mr_results_hnc_corr$b[10] <- -1*(mr_ivw_LDLR$Estimate)
mr_results_hnc_corr$b[11] <- -1*(mr_egger_LDLR$Estimate)
mr_results_hnc_corr$b[12] <- -1*(mr_median_LDLR$Estimate)
mr_results_hnc_corr$b[13] <- -1*(mr_ivw_CETP$Estimate)
mr_results_hnc_corr$b[14] <- -1*(mr_egger_CETP$Estimate)
mr_results_hnc_corr$b[15] <- -1*(mr_median_CETP$Estimate)

mr_results_hnc_corr$se[1] <- mr_ivw_HMGCR$StdError
mr_results_hnc_corr$se[2] <- mr_egger_HMGCR$StdError.Est
mr_results_hnc_corr$se[3] <- mr_median_HMGCR$StdError
mr_results_hnc_corr$se[4] <- mr_ivw_NPC1L1$StdError
mr_results_hnc_corr$se[5] <- mr_egger_NPC1L1$StdError.Est
mr_results_hnc_corr$se[6] <- mr_median_NPC1L1$StdError
mr_results_hnc_corr$se[7] <- mr_ivw_PCSK9$StdError
mr_results_hnc_corr$se[8] <- mr_egger_PCSK9$StdError.Est
mr_results_hnc_corr$se[9] <- mr_median_PCSK9$StdError
mr_results_hnc_corr$se[10] <- mr_ivw_LDLR$StdError
mr_results_hnc_corr$se[11] <- mr_egger_LDLR$StdError.Est
mr_results_hnc_corr$se[12] <- mr_median_LDLR$StdError
mr_results_hnc_corr$se[13] <- mr_ivw_CETP$StdError
mr_results_hnc_corr$se[14] <- mr_egger_CETP$StdError.Est
mr_results_hnc_corr$se[15] <- mr_median_CETP$StdError

mr_results_hnc_corr$pval[1] <- mr_ivw_HMGCR$Pvalue
mr_results_hnc_corr$pval[2] <- mr_egger_HMGCR$Pvalue.Est
mr_results_hnc_corr$pval[3] <- mr_median_HMGCR$Pvalue
mr_results_hnc_corr$pval[4] <- mr_ivw_NPC1L1$Pvalue
mr_results_hnc_corr$pval[5] <- mr_egger_NPC1L1$Pvalue.Est
mr_results_hnc_corr$pval[6] <- mr_median_NPC1L1$Pvalue
mr_results_hnc_corr$pval[7] <- mr_ivw_PCSK9$Pvalue
mr_results_hnc_corr$pval[8] <- mr_egger_PCSK9$Pvalue.Est
mr_results_hnc_corr$pval[9] <- mr_median_PCSK9$Pvalue
mr_results_hnc_corr$pval[10] <- mr_ivw_LDLR$Pvalue
mr_results_hnc_corr$pval[11] <- mr_egger_LDLR$Pvalue.Est
mr_results_hnc_corr$pval[12] <- mr_median_LDLR$Pvalue
mr_results_hnc_corr$pval[13] <- mr_ivw_CETP$Pvalue
mr_results_hnc_corr$pval[14] <- mr_egger_CETP$Pvalue.Est
mr_results_hnc_corr$pval[15] <- mr_median_CETP$Pvalue

# Estimate odds ratio and 95% confidence interval
mr_results_hnc_corr$or <- exp(mr_results_hnc_corr$b)
mr_results_hnc_corr$cil <- exp(mr_results_hnc_corr$b-1.96*mr_results_hnc_corr$se)
mr_results_hnc_corr$ciu <- exp(mr_results_hnc_corr$b+1.96*mr_results_hnc_corr$se)

write.csv(mr_results_hnc_corr, "./oc_results_corr.csv")

#Hetereogenity test
headers<-c("outcome","exposure","method","Q","Q_df", "Q_pval")
mr_heterogenity_test_corr <- as.data.frame(matrix(,ncol=6,nrow=10))
names(mr_heterogenity_test_corr)<-headers
mr_heterogenity_test_corr$outcome <- "HNC" 
mr_heterogenity_test_corr$exposure <- c("HMGCR", "HMGCR", "NPC1L1", "NPC1L1", "PCSK9", "PCSK9", "LDLR", "LDLR", "CETP", "CETP")
mr_heterogenity_test_corr$method <- c("IVW", "MR Egger", "IVW", "MR Egger", "IVW", "MR Egger", "IVW", "MR Egger", "IVW", "MR Egger")

mr_heterogenity_test_corr$Q[1] <- mr_ivw_HMGCR$Heter.Stat[1]
mr_heterogenity_test_corr$Q[2] <- mr_egger_HMGCR$Heter.Stat[1]
mr_heterogenity_test_corr$Q[3] <- mr_ivw_NPC1L1$Heter.Stat[1]
mr_heterogenity_test_corr$Q[4] <- mr_egger_NPC1L1$Heter.Stat[1]
mr_heterogenity_test_corr$Q[5] <- mr_ivw_PCSK9$Heter.Stat[1]
mr_heterogenity_test_corr$Q[6] <- mr_egger_PCSK9$Heter.Stat[1]
mr_heterogenity_test_corr$Q[7] <- mr_ivw_LDLR$Heter.Stat[1]
mr_heterogenity_test_corr$Q[8] <- mr_egger_LDLR$Heter.Stat[1]
mr_heterogenity_test_corr$Q[9] <- mr_ivw_CETP$Heter.Stat[1]
mr_heterogenity_test_corr$Q[10] <- mr_egger_CETP$Heter.Stat[1]

mr_heterogenity_test_corr$Q_df[1] <- mr_ivw_HMGCR$SNPs - 1
mr_heterogenity_test_corr$Q_df[2] <- mr_egger_HMGCR$SNPs - 2
mr_heterogenity_test_corr$Q_df[3] <- mr_ivw_NPC1L1$SNPs - 1
mr_heterogenity_test_corr$Q_df[4] <- mr_egger_NPC1L1$SNPs - 2
mr_heterogenity_test_corr$Q_df[5] <- mr_ivw_PCSK9$SNPs - 1
mr_heterogenity_test_corr$Q_df[6] <- mr_egger_PCSK9$SNPs - 2
mr_heterogenity_test_corr$Q_df[7] <- mr_ivw_LDLR$SNPs - 1
mr_heterogenity_test_corr$Q_df[8] <- mr_egger_LDLR$SNPs - 2 
mr_heterogenity_test_corr$Q_df[9] <- mr_ivw_CETP$SNPs - 1
mr_heterogenity_test_corr$Q_df[10] <- mr_egger_CETP$SNPs - 2 

mr_heterogenity_test_corr$Q_pval[1] <- mr_ivw_HMGCR$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[2] <- mr_egger_HMGCR$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[3] <- mr_ivw_NPC1L1$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[4] <- mr_egger_NPC1L1$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[5] <- mr_ivw_PCSK9$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[6] <- mr_egger_PCSK9$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[7] <- mr_ivw_LDLR$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[8] <- mr_egger_LDLR$Heter.Stat[2] 
mr_heterogenity_test_corr$Q_pval[9] <- mr_ivw_CETP$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[10] <- mr_egger_CETP$Heter.Stat[2] 

write.csv(mr_heterogenity_test_corr,"./heterogeneity_results_oc_corr.csv")

# Egger intercept 
headers<-c("outcome","exposure","egger_intercept","se","pval")
mr_pleiotropy_test_corr <- as.data.frame(matrix(,ncol=5,nrow=5))
names(mr_pleiotropy_test_corr)<-headers
mr_pleiotropy_test_corr$outcome <- "HNC" 
mr_pleiotropy_test_corr$exposure <- c("HMGCR", "NPC1L1","PCSK9","LDLR","CETP")
mr_pleiotropy_test_corr$egger_intercept[1] <- -1*(mr_egger_HMGCR$Intercept)
mr_pleiotropy_test_corr$egger_intercept[2] <- -1*(mr_egger_NPC1L1$Intercept)
mr_pleiotropy_test_corr$egger_intercept[3] <- -1*(mr_egger_PCSK9$Intercept)
mr_pleiotropy_test_corr$egger_intercept[4] <- -1*(mr_egger_LDLR$Intercept)
mr_pleiotropy_test_corr$egger_intercept[5] <- -1*(mr_egger_CETP$Intercept)
mr_pleiotropy_test_corr$se[1] <- mr_egger_HMGCR$StdError.Int
mr_pleiotropy_test_corr$se[2] <- mr_egger_NPC1L1$StdError.Int
mr_pleiotropy_test_corr$se[3] <- mr_egger_PCSK9$StdError.Int
mr_pleiotropy_test_corr$se[4] <- mr_egger_LDLR$StdError.Int
mr_pleiotropy_test_corr$se[5] <- mr_egger_CETP$StdError.Int
mr_pleiotropy_test_corr$pval[1] <- mr_egger_HMGCR$Pvalue.Int
mr_pleiotropy_test_corr$pval[2] <- mr_egger_NPC1L1$Pvalue.Int
mr_pleiotropy_test_corr$pval[3] <- mr_egger_PCSK9$Pvalue.Int
mr_pleiotropy_test_corr$pval[4] <- mr_egger_LDLR$Pvalue.Int
mr_pleiotropy_test_corr$se[5] <- mr_egger_CETP$StdError.Int

write.csv(mr_pleiotropy_test_corr,"./pleiotropy_results_oc_corr.csv")

########################################################
# OROPHARYNGEAL CANCER                                 #
########################################################

#2. SNP-OUTCOME SUMMARY DATA
drug_hnc <- read.csv("./drug_exp_out_opc.csv", header=TRUE)
drug_hnc <- read_outcome_data("drug_exp_out_opc.csv", drug_exp_data$SNP, sep=",")

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
dat_hnc <- harmonise_data(drug_exp_data, drug_hnc, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
mr_results_hnc <- mr(dat_hnc)
mr_results_hnc

#Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results_hnc, dat_hnc)
mr_scatter
pdf("drug_exp_opc_scatter.pdf",  width = 15, height = 20)
ggarrange(mr_scatter$`4bR9XM.At7QtJ`,mr_scatter$`56iCmO.At7QtJ`,mr_scatter$`8WzrD5.At7QtJ`,mr_scatter$`A1eL1U.At7QtJ`,mr_scatter$`F5jola.At7QtJ`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#5. ESTIMATE AND PLOT THE CAUSAL EFFECTS OF THE TRAIT ON THE OUTCOME
# Flip the results (to protective effect of LDL-lowering)
mr_results_hnc$b <- -1*(mr_results_hnc$b)

# Estimate odds ratio and 95% confidence interval
mr_results_hnc$or <- exp(mr_results_hnc$b)
mr_results_hnc$cil <- exp(mr_results_hnc$b-1.96*mr_results_hnc$se)
mr_results_hnc$ciu <- exp(mr_results_hnc$b+1.96*mr_results_hnc$se)

results_hnc<-cbind.data.frame(mr_results_hnc$outcome,mr_results_hnc$exposure,mr_results_hnc$nsnp,mr_results_hnc$method,mr_results_hnc$b,mr_results_hnc$se,mr_results_hnc$pval,mr_results_hnc$or,mr_results_hnc$cil,mr_results_hnc$ciu)

#Export results
write.csv(results_hnc,"./gameon_opc_results.csv")

#6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(dat_hnc)
res_single$b <- -1*(res_single$b)

mr_forest <- mr_forest_plot(res_single)
mr_forest

pdf("drug_exp_opc_forest.pdf",  width = 15, height = 20)
ggarrange(mr_forest$`4bR9XM.At7QtJ`,mr_forest$`56iCmO.At7QtJ`,mr_forest$`8WzrD5.At7QtJ`,mr_forest$`A1eL1U.At7QtJ`,mr_forest$`F5jola.At7QtJ`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(dat_hnc)
mr_loo <- mr_leaveoneout_plot(res_loo)

pdf("drug_exp_opc_loo.pdf",  width = 15, height = 20)
ggarrange(mr_loo$`4bR9XM.At7QtJ`,mr_loo$`56iCmO.At7QtJ`,mr_loo$`8WzrD5.At7QtJ`,mr_loo$`A1eL1U.At7QtJ`,mr_loo$`F5jola.At7QtJ`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#7. Calculate I-squared, r-squared and F-statistics
#I-squared function
Isq <- function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

#Calculate Isq wieghted and unweighted
I2<-c()
dat <- harmonise_data(drug_exp_data, drug_hnc, action = 1)
str(dat)

#F-statistic
dat$samplesize.exposure <- 188578
dat$samplesize.outcome <- 6034
dat <- steiger_filtering(dat) 

N = dat$samplesize.exposure[1] #sample size
K = length(dat$SNP) #number of SNPs
total_r2 <- sum(dat$rsq.exposure) 
Fstat <- (N-K-1)/K * total_r2 / (1-total_r2)
total_r2
Fstat

#Rename required columns
dat $BetaXG<-dat $beta.exposure
dat $seBetaXG<-dat $se.exposure
BetaXG   = dat $BetaXG
seBetaXG = dat $seBetaXG 
seBetaYG<-dat $se.outcome

BXG = abs(BetaXG)         
# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F = BXG^2/seBetaXG^2
mF = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted

#Save results
output<-cbind (F, mF, Isq_unweighted, Isq_weighted)
I2<-rbind(I2, output)
colnames(I2) <- c("F", "Exposure", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="regression_dilution_isq_weighted_insertgenename.csv", row.names = FALSE)

#8. SIMEX correction
#Create empty dataframe to store output
simexegger<-c()

#Run simex 
#Rename required columns
dat$BetaXG<-dat$beta.exposure
dat$seBetaXG<-dat$se.exposure
dat$BetaYG<-dat$beta.outcome
dat$seBetaYG<-dat$se.outcome
BetaXG <- dat$BetaXG
BetaYG <- dat$BetaYG
seBetaXG <- dat$seBetaXG
seBetaYG <- dat$seBetaYG

BYG <- BetaYG*sign(BetaXG)#Pre-processing steps to ensure all gene--exposure estimates are positive
BXG <- abs(BetaXG)         

#MR-Egger regression (weighted) 
Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

#MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

#Simulation extrapolation 
mod.sim1 <- simex(Fit1,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod.sim2 <- simex(Fit2,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)
#Print results
mod1
mod2

#8. ASSESS HETEROGENEITY AND PLEIOTROPY 

#Hetereogenity test
mr_heterogeneity(dat_hnc)
write.csv(mr_heterogeneity(dat_hnc),"./heterogeneity_results_opc.csv")

#Hetereogenity test
mr_pleiotropy_test(dat_hnc)
write.csv(mr_pleiotropy_test(dat_hnc),"./pleiotropy_results_opc.csv")

# MR PRESSO
# Run MR-PRESSO global method
HMGCR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat_hnc[dat_hnc$exposure=="HMGCR"&dat_hnc$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
NPC1L1_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat_hnc[dat_hnc$exposure=="NPC1L1"&dat_hnc$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
PCSK9_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat_hnc[dat_hnc$exposure=="PCSK9"&dat_hnc$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LDLR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat_hnc[dat_hnc$exposure=="LDLR"&dat_hnc$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
CETP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat_hnc[dat_hnc$exposure=="CETP"&dat_hnc$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
#NB Not enough intrumental variables for LDLR 

headers<-c("outcome", "exposure", "RSSobs","Pvalue")
mr_presso_global <- as.data.frame(matrix(,ncol=4,nrow=5))
names(mr_presso_global)<-headers  
mr_presso_global$outcome <- "HNC" 
mr_presso_global$exposure <- c("HMGCR", "NPC1L1", "PCSK9", "LDLR", "CETP")

mr_presso_global$RSSobs[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$Pvalue

write.csv(mr_presso_global, "./mr_presso_opc_global.csv")

#9. ACCOUNT FOR LD STRUCTURE

dat_hnc <- dat_hnc[dat_hnc$mr_keep=="TRUE",]

HMGCR <- dat_hnc[dat_hnc$exposure=="HMGCR",]
LDmatrix(HMGCR$SNP, pop="CEU", r2d = "r2", token = "b0a2f3f3c1f4", file = "HMGCR_cor.txt")
NPC1L1 <- dat_hnc[dat_hnc$exposure=="NPC1L1",]
LDmatrix(NPC1L1$SNP, pop="CEU", r2d = "r2", token = "b0a2f3f3c1f4", file = "NPC1L1_cor.txt")
PCSK9 <- dat_hnc[dat_hnc$exposure=="PCSK9",]
LDmatrix(PCSK9$SNP, pop="CEU", r2d = "r2", token = "b0a2f3f3c1f4", file = "PCSK9_cor.txt")
LDLR <- dat_hnc[dat_hnc$exposure=="LDLR",]
LDmatrix(LDLR$SNP, pop="CEU", r2d = "r2", token = "b0a2f3f3c1f4", file = "LDLR_cor.txt")
CETP <- dat_hnc[dat_hnc$exposure=="CETP",]
LDmatrix(CETP$SNP, pop="CEU", r2d = "r2", token = "b0a2f3f3c1f4", file = "CETP_cor.txt")

correl <- read.table("HMGCR_cor.txt", header=T)
correl <- correl[order(correl$RS_number),]
correl <- correl[,order(names(correl))]
HMGCR <- HMGCR[HMGCR$SNP %in% names(correl),]
correl <- within(correl, rm(RS_number))
colnames(correl) <- NULL
rownames(correl) <- NULL
correl <- data.matrix(correl)
mr_ivw_HMGCR <- mr_ivw(mr_input(bx = HMGCR$beta.exposure, bxse = HMGCR$se.exposure, by = HMGCR$beta.outcome, byse = HMGCR$se.outcome, correlation = correl))
mr_egger_HMGCR <- mr_egger(mr_input(bx = HMGCR$beta.exposure, bxse = HMGCR$se.exposure, by = HMGCR$beta.outcome, byse = HMGCR$se.outcome, correlation = correl))
mr_median_HMGCR <- mr_median(mr_input(bx = HMGCR$beta.exposure, bxse = HMGCR$se.exposure, by = HMGCR$beta.outcome, byse = HMGCR$se.outcome, correlation = correl))

correl <- read.table("NPC1L1_cor.txt", header=T)
correl <- correl[order(correl$RS_number),]
correl <- correl[,order(names(correl))]
NPC1L1 <- NPC1L1[NPC1L1$SNP %in% names(correl),]
correl <- within(correl, rm(RS_number))
colnames(correl) <- NULL
rownames(correl) <- NULL
correl <- data.matrix(correl)
mr_ivw_NPC1L1 <- mr_ivw(mr_input(bx = NPC1L1$beta.exposure, bxse = NPC1L1$se.exposure, by = NPC1L1$beta.outcome, byse = NPC1L1$se.outcome, correlation = correl))
mr_egger_NPC1L1 <- mr_egger(mr_input(bx = NPC1L1$beta.exposure, bxse = NPC1L1$se.exposure, by = NPC1L1$beta.outcome, byse = NPC1L1$se.outcome, correlation = correl))
mr_median_NPC1L1 <- mr_median(mr_input(bx = NPC1L1$beta.exposure, bxse = NPC1L1$se.exposure, by = NPC1L1$beta.outcome, byse = NPC1L1$se.outcome, correlation = correl))

correl <- read.table("PCSK9_cor.txt", header=T)
correl <- correl[order(correl$RS_number),]
correl <- correl[,order(names(correl))]
PCSK9 <- PCSK9[PCSK9$SNP %in% names(correl),]
correl <- within(correl, rm(RS_number))
colnames(correl) <- NULL
rownames(correl) <- NULL
correl <- data.matrix(correl)
mr_ivw_PCSK9 <- mr_ivw(mr_input(bx = PCSK9$beta.exposure, bxse = PCSK9$se.exposure, by = PCSK9$beta.outcome, byse = PCSK9$se.outcome, correlation = correl))
mr_egger_PCSK9 <- mr_egger(mr_input(bx = PCSK9$beta.exposure, bxse = PCSK9$se.exposure, by = PCSK9$beta.outcome, byse = PCSK9$se.outcome, correlation = correl))
mr_median_PCSK9 <- mr_median(mr_input(bx = PCSK9$beta.exposure, bxse = PCSK9$se.exposure, by = PCSK9$beta.outcome, byse = PCSK9$se.outcome, correlation = correl))

correl <- read.table("LDLR_cor.txt", header=T)
correl <- correl[order(correl$RS_number),]
correl <- correl[,order(names(correl))]
LDLR <- LDLR[LDLR$SNP %in% names(correl),]
correl <- within(correl, rm(RS_number))
colnames(correl) <- NULL
rownames(correl) <- NULL
correl <- data.matrix(correl)
mr_ivw_LDLR <- mr_ivw(mr_input(bx = LDLR$beta.exposure, bxse = LDLR$se.exposure, by = LDLR$beta.outcome, byse = LDLR$se.outcome, correlation = correl))
mr_egger_LDLR <- mr_egger(mr_input(bx = LDLR$beta.exposure, bxse = LDLR$se.exposure, by = LDLR$beta.outcome, byse = LDLR$se.outcome, correlation = correl))
mr_median_LDLR <- mr_median(mr_input(bx = LDLR$beta.exposure, bxse = LDLR$se.exposure, by = LDLR$beta.outcome, byse = LDLR$se.outcome, correlation = correl))

correl <- read.table("CETP_cor.txt", header=T)
correl <- correl[order(correl$RS_number),]
correl <- correl[,order(names(correl))]
CETP <- CETP[CETP$SNP %in% names(correl),]
correl <- within(correl, rm(RS_number))
colnames(correl) <- NULL
rownames(correl) <- NULL
correl <- data.matrix(correl)
mr_ivw_CETP <- mr_ivw(mr_input(bx = CETP$beta.exposure, bxse = CETP$se.exposure, by = CETP$beta.outcome, byse = CETP$se.outcome, correlation = correl))
mr_egger_CETP<- mr_egger(mr_input(bx = CETP$beta.exposure, bxse = CETP$se.exposure, by = CETP$beta.outcome, byse = CETP$se.outcome, correlation = correl))
mr_median_CETP <- mr_median(mr_input(bx = CETP$beta.exposure, bxse = CETP$se.exposure, by = CETP$beta.outcome, byse = CETP$se.outcome, correlation = correl))

# Flip the results (to protective effect of LDL-lowering)
headers<-c("outcome","exposure","method","b","se","pval")
mr_results_hnc_corr <- as.data.frame(matrix(,ncol=6,nrow=15))
names(mr_results_hnc_corr)<-headers

mr_results_hnc_corr$outcome <- "HNC" 
mr_results_hnc_corr$exposure <- c("HMGCR", "HMGCR", "HMGCR", "NPC1L1", "NPC1L1", "NPC1L1", "PCSK9", "PCSK9", "PCSK9", "LDLR", "LDLR", "LDLR", "CETP", "CETP", "CETP")
mr_results_hnc_corr$method <- c("ivw", "egger", "median", "ivw", "egger", "median", "ivw", "egger", "median", "ivw", "egger", "median","ivw", "egger", "median")

mr_results_hnc_corr$b[1] <- -1*(mr_ivw_HMGCR$Estimate)
mr_results_hnc_corr$b[2] <- -1*(mr_egger_HMGCR$Estimate)
mr_results_hnc_corr$b[3] <- -1*(mr_median_HMGCR$Estimate)
mr_results_hnc_corr$b[4] <- -1*(mr_ivw_NPC1L1$Estimate)
mr_results_hnc_corr$b[5] <- -1*(mr_egger_NPC1L1$Estimate)
mr_results_hnc_corr$b[6] <- -1*(mr_median_NPC1L1$Estimate)
mr_results_hnc_corr$b[7] <- -1*(mr_ivw_PCSK9$Estimate)
mr_results_hnc_corr$b[8] <- -1*(mr_egger_PCSK9$Estimate)
mr_results_hnc_corr$b[9] <- -1*(mr_median_PCSK9$Estimate)
mr_results_hnc_corr$b[10] <- -1*(mr_ivw_LDLR$Estimate)
mr_results_hnc_corr$b[11] <- -1*(mr_egger_LDLR$Estimate)
mr_results_hnc_corr$b[12] <- -1*(mr_median_LDLR$Estimate)
mr_results_hnc_corr$b[13] <- -1*(mr_ivw_CETP$Estimate)
mr_results_hnc_corr$b[14] <- -1*(mr_egger_CETP$Estimate)
mr_results_hnc_corr$b[15] <- -1*(mr_median_CETP$Estimate)

mr_results_hnc_corr$se[1] <- mr_ivw_HMGCR$StdError
mr_results_hnc_corr$se[2] <- mr_egger_HMGCR$StdError.Est
mr_results_hnc_corr$se[3] <- mr_median_HMGCR$StdError
mr_results_hnc_corr$se[4] <- mr_ivw_NPC1L1$StdError
mr_results_hnc_corr$se[5] <- mr_egger_NPC1L1$StdError.Est
mr_results_hnc_corr$se[6] <- mr_median_NPC1L1$StdError
mr_results_hnc_corr$se[7] <- mr_ivw_PCSK9$StdError
mr_results_hnc_corr$se[8] <- mr_egger_PCSK9$StdError.Est
mr_results_hnc_corr$se[9] <- mr_median_PCSK9$StdError
mr_results_hnc_corr$se[10] <- mr_ivw_LDLR$StdError
mr_results_hnc_corr$se[11] <- mr_egger_LDLR$StdError.Est
mr_results_hnc_corr$se[12] <- mr_median_LDLR$StdError
mr_results_hnc_corr$se[13] <- mr_ivw_CETP$StdError
mr_results_hnc_corr$se[14] <- mr_egger_CETP$StdError.Est
mr_results_hnc_corr$se[15] <- mr_median_CETP$StdError

mr_results_hnc_corr$pval[1] <- mr_ivw_HMGCR$Pvalue
mr_results_hnc_corr$pval[2] <- mr_egger_HMGCR$Pvalue.Est
mr_results_hnc_corr$pval[3] <- mr_median_HMGCR$Pvalue
mr_results_hnc_corr$pval[4] <- mr_ivw_NPC1L1$Pvalue
mr_results_hnc_corr$pval[5] <- mr_egger_NPC1L1$Pvalue.Est
mr_results_hnc_corr$pval[6] <- mr_median_NPC1L1$Pvalue
mr_results_hnc_corr$pval[7] <- mr_ivw_PCSK9$Pvalue
mr_results_hnc_corr$pval[8] <- mr_egger_PCSK9$Pvalue.Est
mr_results_hnc_corr$pval[9] <- mr_median_PCSK9$Pvalue
mr_results_hnc_corr$pval[10] <- mr_ivw_LDLR$Pvalue
mr_results_hnc_corr$pval[11] <- mr_egger_LDLR$Pvalue.Est
mr_results_hnc_corr$pval[12] <- mr_median_LDLR$Pvalue
mr_results_hnc_corr$pval[13] <- mr_ivw_CETP$Pvalue
mr_results_hnc_corr$pval[14] <- mr_egger_CETP$Pvalue.Est
mr_results_hnc_corr$pval[15] <- mr_median_CETP$Pvalue

# Estimate odds ratio and 95% confidence interval
mr_results_hnc_corr$or <- exp(mr_results_hnc_corr$b)
mr_results_hnc_corr$cil <- exp(mr_results_hnc_corr$b-1.96*mr_results_hnc_corr$se)
mr_results_hnc_corr$ciu <- exp(mr_results_hnc_corr$b+1.96*mr_results_hnc_corr$se)

write.csv(mr_results_hnc_corr, "./opc_results_corr.csv")

#Hetereogenity test
headers<-c("outcome","exposure","method","Q","Q_df", "Q_pval")
mr_heterogenity_test_corr <- as.data.frame(matrix(,ncol=6,nrow=10))
names(mr_heterogenity_test_corr)<-headers
mr_heterogenity_test_corr$outcome <- "HNC" 
mr_heterogenity_test_corr$exposure <- c("HMGCR", "HMGCR", "NPC1L1", "NPC1L1", "PCSK9", "PCSK9", "LDLR", "LDLR", "CETP", "CETP")
mr_heterogenity_test_corr$method <- c("IVW", "MR Egger", "IVW", "MR Egger", "IVW", "MR Egger", "IVW", "MR Egger", "IVW", "MR Egger")

mr_heterogenity_test_corr$Q[1] <- mr_ivw_HMGCR$Heter.Stat[1]
mr_heterogenity_test_corr$Q[2] <- mr_egger_HMGCR$Heter.Stat[1]
mr_heterogenity_test_corr$Q[3] <- mr_ivw_NPC1L1$Heter.Stat[1]
mr_heterogenity_test_corr$Q[4] <- mr_egger_NPC1L1$Heter.Stat[1]
mr_heterogenity_test_corr$Q[5] <- mr_ivw_PCSK9$Heter.Stat[1]
mr_heterogenity_test_corr$Q[6] <- mr_egger_PCSK9$Heter.Stat[1]
mr_heterogenity_test_corr$Q[7] <- mr_ivw_LDLR$Heter.Stat[1]
mr_heterogenity_test_corr$Q[8] <- mr_egger_LDLR$Heter.Stat[1]
mr_heterogenity_test_corr$Q[9] <- mr_ivw_CETP$Heter.Stat[1]
mr_heterogenity_test_corr$Q[10] <- mr_egger_CETP$Heter.Stat[1]

mr_heterogenity_test_corr$Q_df[1] <- mr_ivw_HMGCR$SNPs - 1
mr_heterogenity_test_corr$Q_df[2] <- mr_egger_HMGCR$SNPs - 2
mr_heterogenity_test_corr$Q_df[3] <- mr_ivw_NPC1L1$SNPs - 1
mr_heterogenity_test_corr$Q_df[4] <- mr_egger_NPC1L1$SNPs - 2
mr_heterogenity_test_corr$Q_df[5] <- mr_ivw_PCSK9$SNPs - 1
mr_heterogenity_test_corr$Q_df[6] <- mr_egger_PCSK9$SNPs - 2
mr_heterogenity_test_corr$Q_df[7] <- mr_ivw_LDLR$SNPs - 1
mr_heterogenity_test_corr$Q_df[8] <- mr_egger_LDLR$SNPs - 2 
mr_heterogenity_test_corr$Q_df[9] <- mr_ivw_CETP$SNPs - 1
mr_heterogenity_test_corr$Q_df[10] <- mr_egger_CETP$SNPs - 2 

mr_heterogenity_test_corr$Q_pval[1] <- mr_ivw_HMGCR$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[2] <- mr_egger_HMGCR$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[3] <- mr_ivw_NPC1L1$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[4] <- mr_egger_NPC1L1$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[5] <- mr_ivw_PCSK9$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[6] <- mr_egger_PCSK9$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[7] <- mr_ivw_LDLR$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[8] <- mr_egger_LDLR$Heter.Stat[2] 
mr_heterogenity_test_corr$Q_pval[9] <- mr_ivw_CETP$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[10] <- mr_egger_CETP$Heter.Stat[2] 

write.csv(mr_heterogenity_test_corr,"./heterogeneity_results_opc_corr.csv")

# Egger intercept 
headers<-c("outcome","exposure","egger_intercept","se","pval")
mr_pleiotropy_test_corr <- as.data.frame(matrix(,ncol=5,nrow=5))
names(mr_pleiotropy_test_corr)<-headers
mr_pleiotropy_test_corr$outcome <- "HNC" 
mr_pleiotropy_test_corr$exposure <- c("HMGCR", "NPC1L1","PCSK9","LDLR","CETP")
mr_pleiotropy_test_corr$egger_intercept[1] <- -1*(mr_egger_HMGCR$Intercept)
mr_pleiotropy_test_corr$egger_intercept[2] <- -1*(mr_egger_NPC1L1$Intercept)
mr_pleiotropy_test_corr$egger_intercept[3] <- -1*(mr_egger_PCSK9$Intercept)
mr_pleiotropy_test_corr$egger_intercept[4] <- -1*(mr_egger_LDLR$Intercept)
mr_pleiotropy_test_corr$egger_intercept[5] <- -1*(mr_egger_CETP$Intercept)
mr_pleiotropy_test_corr$se[1] <- mr_egger_HMGCR$StdError.Int
mr_pleiotropy_test_corr$se[2] <- mr_egger_NPC1L1$StdError.Int
mr_pleiotropy_test_corr$se[3] <- mr_egger_PCSK9$StdError.Int
mr_pleiotropy_test_corr$se[4] <- mr_egger_LDLR$StdError.Int
mr_pleiotropy_test_corr$se[5] <- mr_egger_CETP$StdError.Int
mr_pleiotropy_test_corr$pval[1] <- mr_egger_HMGCR$Pvalue.Int
mr_pleiotropy_test_corr$pval[2] <- mr_egger_NPC1L1$Pvalue.Int
mr_pleiotropy_test_corr$pval[3] <- mr_egger_PCSK9$Pvalue.Int
mr_pleiotropy_test_corr$pval[4] <- mr_egger_LDLR$Pvalue.Int
mr_pleiotropy_test_corr$se[5] <- mr_egger_CETP$StdError.Int

write.csv(mr_pleiotropy_test_corr,"./pleiotropy_results_opc_corr.csv")

#####################################
# CARDIoGRAM data check
#####################################

#2. SELECT SNP-OUTCOME SUMMARY DATA                   
ao <- available_outcomes()
drug_cvd <- extract_outcome_data(drug_exp_data$SNP, c("ieu-a-7"))

#3. HARMONIZE DATASETS
dat_cvd <- harmonise_data(drug_exp_data, drug_cvd, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME e.g. HNC
mr_results_cvd <- mr(dat_cvd)
mr_results_cvd

#5. ESTIMATE AND PLOT THE CAUSAL EFFECTS OF THE TRAIT ON THE OUTCOME

# Flip the results (to protective effect of LDL-lowering)
mr_results_cvd$b <- -1*(mr_results_cvd$b)

# Estimate odds ratio and 95% confidence interval
mr_results_cvd$or <- exp(mr_results_cvd$b)
mr_results_cvd$cil <- exp(mr_results_cvd$b-1.96*mr_results_cvd$se)
mr_results_cvd$ciu <- exp(mr_results_cvd$b+1.96*mr_results_cvd$se)

results_cvd <-cbind.data.frame(mr_results_cvd$outcome,mr_results_cvd$exposure,mr_results_cvd$nsnp,mr_results_cvd$method,mr_results_cvd$b,mr_results_cvd$se,mr_results_cvd$pval,mr_results_cvd$or,mr_results_cvd$cil,mr_results_cvd$ciu)

#Export results
write.csv(results_cvd,"./cvd_hnc_results.csv")

########################################################
#COMBINED PLOTS                                        #
########################################################

#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer
res_single <- mr_singlesnp(dat_cvd)
res_single$b <- -1*(res_single$b)
mr_forest_plot(res_single)
mr_scatter_plot(mr_results_cvd, dat_cvd)

#by individual SNP Forest Plot
res_hmgcr <- mr_singlesnp(dat_cvd)
res_hmgcr$b <- -1*(res_hmgcr$b)
mr_forest_plot(res_hmgcr)
mr_scatter_plot(mr_results_cvd, dat_cvd)

res_pcsk9 <- mr_singlesnp(dat_cvd)
res_pcsk9 $b <- -1*(res_pcsk9 $b)
mr_forest_plot(res_pcsk9 )
mr_scatter_plot(mr_results_cvd, dat_cvd)

res_npc1l1<- mr_singlesnp(dat_cvd)
res_npc1l1 $b <- -1*(res_npc1l1 $b)
mr_forest_plot(res_npc1l1 )
mr_scatter_plot(mr_results_cvd, dat_cvd)

res_ldlr<- mr_singlesnp(dat_cvd)
res_ldlr $b <- -1*(res_ldlr $b)
mr_forest_plot(res_ldlr )
mr_scatter_plot(mr_results_cvd, dat_cvd)

#all MR estimates on one Forest plot 
dat <- read.table("HMGCR_cvd.txt", sep="\t", header=T)
dat <- read.table("LDLR_cvd.txt", sep="\t", header=T)
dat <- read.table("NPC1L1_cvd.txt", sep="\t", header=T)
dat <- read.table("PCSK9_cvd.txt", sep="\t", header=T)

dat <- dat[order(dat$method),]
row.names(dat) <- NULL

dat$xpos <- nrow(dat) - as.numeric(row.names(dat))
dat$cil <- dat$b - (1.96*dat$se)
dat$ciu <- dat$b + (1.96*dat$se)

dat_plot <- ggplot(data=dat, aes(x=xpos, y=b, ymin=cil, ymax=ciu, colour=as.factor(outcome))) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=0 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Method") + ylab("Effect (95% CI)") +
  theme_bw()  # use a white background + 

dat_plot <- dat_plot + scale_x_continuous(breaks = c(10, 7, 4, 1), labels = (levels(dat$method))) 
dat_plot <- dat_plot + theme(legend.position="none")

#dat_plot <- dat_plot + scale_y_continuous(limits = c(-2,7))
dat_plot1 <- dat_plot + ggtitle("HMGCR")
dat_plot2 <- dat_plot + ggtitle("LDLR")
dat_plot3 <- dat_plot + ggtitle("NPC1L1")
dat_plot4 <- dat_plot + ggtitle("PCSK9")

library(ggpubr)
ggarrange(dat_plot1, dat_plot2, ncol=2, nrow=1, common.legend=T)
ggarrange(dat_plot3, dat_plot4, ncol=2, nrow=1, common.legend=T)

# Run each plot separately
dat <- read.table("allsites_HMGCR.txt", sep="\t", header=T)
dat <- read.table("allsites_NPC1L1.txt", sep="\t", header=T)
dat <- read.table("allsites_CETP.txt", sep="\t", header=T)
dat <- read.table("allsites_PCSK9.txt", sep="\t", header=T)
dat <- read.table("allsites_LDLR.txt", sep="\t", header=T)

dat <- dat[order(dat$method),]
row.names(dat) <- NULL

dat$xpos <- nrow(dat) - as.numeric(row.names(dat))
dat$cil <- dat$or - (1.96*dat$se)
dat$ciu <- dat$or + (1.96*dat$se)

library(ggplot2)
dat_plot <- ggplot(data=dat, aes(x=xpos, y=or, ymin=cil, ymax=ciu, colour=as.factor(outcome))) +
  geom_pointrange(shape=18) + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=0 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Method") + ylab("Effect (95% CI)") +
  theme_bw()  # use a white background + 

dat_plot <- dat_plot + scale_color_grey(name = "Cancer", breaks=c("OPC+OC","OPC","OC")) +
  theme(legend.title = element_text(colour="black", size = 12, face = "bold") + theme_classic())

dat_plot <- dat_plot + scale_x_continuous(breaks = c(10, 7, 4, 1), labels = (levels(dat$method))) 

dat_plot <- dat_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"))

dat_plot <- dat_plot + geom_point() + 
  geom_rect(aes(xmin = 5.5, xmax = 2.5, ymin = -Inf, ymax = Inf),
            fill = "gray", alpha = 0.01, color="NA") +   
  geom_rect(aes(xmin = 11.5, xmax = 8.5, ymin = -Inf, ymax = Inf),
            fill = "gray", alpha = 0.01, color="NA")  

dat_plot1 <- dat_plot + ggtitle("HMGCR")
dat_plot2 <- dat_plot + ggtitle("NPC1L1")
dat_plot3 <- dat_plot + ggtitle("CETP")
dat_plot4 <- dat_plot + ggtitle("PCSK9")
dat_plot5 <- dat_plot + ggtitle("LDLR")

#install.packages("ggpubr")
library(ggpubr)
ggarrange(dat_plot1, ncol=1, nrow=1, common.legend=T)
ggarrange(dat_plot2, ncol=1, nrow=1, common.legend=T)
ggarrange(dat_plot3, ncol=1, nrow=1, common.legend=T)
ggarrange(dat_plot4, ncol=1, nrow=1, common.legend=T)
ggarrange(dat_plot5, ncol=1, nrow=1, common.legend=T)

ggarrange(dat_plot1, dat_plot2, dat_plot3, dat_plot4, dat_plot5,ncol=NULL, nrow=NULL, common.legend=T)

