####################
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
setwd("")

#1. SELECT SNP-EXPOSURE SUMMARY DATA#
# SNP Extraction from MRBase- for other circulating lipid SNPs of traits of interest
# 299- HDL, 300- LDL, 301- TC, 302- TG, 842- ApoA, 843- ApoB

ao <- available_outcomes()
exposure_dat <- extract_instruments(c("ieu-a-299", "ieu-a-300", "ieu-a-301", "ieu-a-302", "ieu-a-842", "ieu-a-843"))

write.table(exposure_dat$SNP, "snplist_allSNPs.txt", row.names=F, col.names=F)
write.table(exposure_dat, "exposure_data.txt", row.names=F, col.names=F)
# Convert file to UNIX format using dos2unix snplist_LDL.txt command in Putty

##############################################
# ORAL AND OROPHARYNGEAL CANCER              #  
##############################################

#2. SELECT SNP-OUTCOME SUMMARY DATA
exp_hnc <- read.csv("./snplist_out_allSNPs_hnc.csv", header=TRUE)
exp_hnc <- read_outcome_data("snplist_out_allSNPs_hnc.csv", exposure_dat$SNP, sep=",")


#3. HARMONIZE DATASETS
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
dat_hnc <- harmonise_data(exposure_dat, exp_hnc, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME e.g. HNC
mr_results_hnc <- mr(dat_hnc)
mr_results_hnc

#Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results_hnc, dat_hnc)
mr_scatter

pdf("allSNPs_hnc_scatter.pdf",  width = 15, height = 20)
ggarrange(mr_scatter$`ieu-a-299.D55bOj`,mr_scatter$`ieu-a-300.D55bOj`,mr_scatter$`ieu-a-301.D55bOj`,mr_scatter$`ieu-a-302.D55bOj`, ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

pdf("allSNPs_hnc_scatter_2.pdf",  width = 15, height = 20)
ggarrange(mr_scatter$`met-c-842.D55bOj`,mr_scatter$`met-c-843.D55bOj`, ncol=2, nrow=2, widths = 2, heights = 1)
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
write.csv(results_hnc,"./mrresults_allSNPs_hnc.csv")

#6. PLOT FOREST AND LEAVE-ONE-OUT 
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(dat_hnc)
res_single$b <- -1*(res_single$b)

mr_forest <- mr_forest_plot(res_single)
mr_forest

pdf("allSNPs_hnc_forest.pdf",  width = 15, height = 20)
ggarrange(mr_forest$`ieu-a-299.D55bOj`,mr_forest$`ieu-a-300.D55bOj`,mr_forest$`ieu-a-301.D55bOj`,mr_forest$`ieu-a-302.D55bOj`, ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

pdf("allSNPs_hnc_forest_2.pdf",  width = 15, height = 20)
ggarrange(mr_forest$`met-c-842.D55bOj`,mr_forest$`met-c-843.D55bOj`, ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(dat_hnc)
mr_loo <- mr_leaveoneout_plot(res_loo)
mr_loo

pdf("allSNPs_hnc_loo.pdf",  width = 15, height = 20)
ggarrange(mr_loo$`ieu-a-299.D55bOj`,mr_loo$`ieu-a-300.D55bOj`,mr_loo$`ieu-a-301.D55bOj`,mr_loo$`ieu-a-302.D55bOj`, ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

pdf("allSNPs_hnc_loo_2.pdf",  width = 15, height = 20)
ggarrange(mr_loo$`met-c-842.D55bOj`,mr_loo$`met-c-843.D55bOj`, ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#7. ASSESS HETEROGENEITY AND PLEIOTROPY 
# I-squared function
Isq <- function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

#calculate Isq wieghted and unweighted
I2<-c()
dat <- harmonise_data(exposure_dat, drug_hnc, action = 1)
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

BXG             = abs(BetaXG)         # gene--exposure estimates are positive  

# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F   = BXG^2/seBetaXG^2
mF  = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted

#Save results
output<-cbind (F, mF, Isq_unweighted, Isq_weighted)
I2<-rbind(I2, output)
colnames(I2) <- c("F", "Exposure", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="regression_dilution_isq_weighted_filename.csv", row.names = FALSE)

# SIMEX correction
#create empty dataframe to store output
simexegger<-c()

#Rename required columns
dat$BetaXG<-dat$beta.exposure
dat$seBetaXG<-dat$se.exposure
dat$BetaYG<-dat$beta.outcome
dat$seBetaYG<-dat$se.outcome
BetaXG <- dat$BetaXG
BetaYG <- dat$BetaYG
seBetaXG <- dat$seBetaXG
seBetaYG <- dat$seBetaYG

BYG <- BetaYG*sign(BetaXG)# Pre-processing steps to ensure all gene--exposure estimates are positive
BXG <- abs(BetaXG)         

# MR-Egger regression (weighted) 
Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

# MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

# Simulation extrapolation 
mod.sim1 <- simex(Fit1,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod.sim2 <- simex(Fit2,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)

mod1
mod2

#Hetereogenity test
mr_heterogeneity(dat_hnc)
write.csv(mr_heterogeneity(dat_hnc),"./heterogeneity_results_allSNPs_hnc.csv")

#Pleiotropy test
mr_pleiotropy_test(dat_hnc)
write.csv(mr_pleiotropy_test(dat_hnc),"./pleiotropy_results_allSNPs_hnc.csv")

# MR PRESSO
HDL_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat_hnc[dat_hnc$exposure=="HDL cholesterol || id:ieu-a-299"&dat_hnc$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LDL_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat_hnc[dat_hnc$exposure=="LDL cholesterol || id:ieu-a-300"&dat_hnc$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
TC_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat_hnc[dat_hnc$exposure=="Total cholesterol || id:ieu-a-301"&dat_hnc$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
TG_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat_hnc[dat_hnc$exposure=="Triglycerides || id:ieu-a-302"&dat_hnc$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
ApoA_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat_hnc[dat_hnc$exposure=="Apolipoprotein A-I || id:met-c-842"&dat_hnc$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
ApoB_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat_hnc[dat_hnc$exposure=="Apolipoprotein B || id:met-c-843"&dat_hnc$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)

headers<-c("outcome", "exposure", "RSSobs","Pvalue")
mr_presso_global <- as.data.frame(matrix(,ncol=4,nrow=6))
names(mr_presso_global)<-headers  
mr_presso_global$outcome <- "HNC" 
mr_presso_global$exposure <- c("HDL cholesterol", "LDL cholesterol", "Total cholesterol", "Triglycerides", "Apolipoprotein A", "Apolipoprotein B" )

mr_presso_global$RSSobs[1] <- HDL_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[1] <- HDL_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[2] <- LDL_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[2] <- LDL_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[3] <- TC_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[3] <- TC_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[4] <- TG_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[4] <- TG_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[5] <- ApoA_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[5] <- ApoA_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[6] <- ApoB_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[6] <- ApoB_presso$`MR-PRESSO results`$`Global Test`$Pvalue

write.csv(mr_presso_global, "./allSNPs_hnc_mr_presso_global.csv")

##############################################
# ORAL CANCER                                #  
##############################################
#2. SELECT SNP-OUTCOME SUMMARY DATA
exp_hnc <- read.csv("./snplist_out_allSNPs_oc.csv", header=TRUE)
exp_hnc <- read_outcome_data("snplist_out_allSNPs_oc.csv", exposure_dat$SNP, sep=",")


#3. HARMONIZE DATASETS
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
dat_hnc <- harmonise_data(exposure_dat, exp_hnc, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME e.g. HNC
mr_results_hnc <- mr(dat_hnc)
mr_results_hnc

#Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results_hnc, dat_hnc)
mr_scatter

pdf("allSNPs_oc_scatter.pdf",  width = 15, height = 20)
ggarrange(mr_scatter$`ieu-a-299.D55bOj`,mr_scatter$`ieu-a-300.D55bOj`,mr_scatter$`ieu-a-301.D55bOj`,mr_scatter$`ieu-a-302.D55bOj`, ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

pdf("allSNPs_oc_scatter_2.pdf",  width = 15, height = 20)
ggarrange(mr_scatter$`met-c-842.D55bOj`,mr_scatter$`met-c-843.D55bOj`, ncol=2, nrow=2, widths = 2, heights = 1)
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
write.csv(results_hnc,"./mrresults_allSNPs_oc.csv")

#6. PLOT FOREST AND LEAVE-ONE-OUT 
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(dat_hnc)
res_single$b <- -1*(res_single$b)

mr_forest <- mr_forest_plot(res_single)
mr_forest

pdf("allSNPs_oc_forest.pdf",  width = 15, height = 20)
ggarrange(mr_forest$`ieu-a-299.D55bOj`,mr_forest$`ieu-a-300.D55bOj`,mr_forest$`ieu-a-301.D55bOj`,mr_forest$`ieu-a-302.D55bOj`, ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

pdf("allSNPs_oc_forest_2.pdf",  width = 15, height = 20)
ggarrange(mr_forest$`met-c-842.D55bOj`,mr_forest$`met-c-843.D55bOj`, ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(dat_hnc)
mr_loo <- mr_leaveoneout_plot(res_loo)
mr_loo

pdf("allSNPs_oc_loo.pdf",  width = 15, height = 20)
ggarrange(mr_loo$`ieu-a-299.D55bOj`,mr_loo$`ieu-a-300.D55bOj`,mr_loo$`ieu-a-301.D55bOj`,mr_loo$`ieu-a-302.D55bOj`, ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

pdf("allSNPs_oc_loo_2.pdf",  width = 15, height = 20)
ggarrange(mr_loo$`met-c-842.D55bOj`,mr_loo$`met-c-843.D55bOj`, ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

##############################################
#OROPHARYNGEAL CANCER                        #  
##############################################
#2. SELECT SNP-OUTCOME SUMMARY DATA
exp_hnc <- read.csv("./snplist_out_allSNPs_opc.csv", header=TRUE)
exp_hnc <- read_outcome_data("snplist_out_allSNPs_opc.csv", exposure_dat$SNP, sep=",")


#3. HARMONIZE DATASETS
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
dat_hnc <- harmonise_data(exposure_dat, exp_hnc, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME e.g. HNC
mr_results_hnc <- mr(dat_hnc)
mr_results_hnc

#Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results_hnc, dat_hnc)
mr_scatter

pdf("allSNPs_opc_scatter.pdf",  width = 15, height = 20)
ggarrange(mr_scatter$`ieu-a-299.D55bOj`,mr_scatter$`ieu-a-300.D55bOj`,mr_scatter$`ieu-a-301.D55bOj`,mr_scatter$`ieu-a-302.D55bOj`, ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

pdf("allSNPs_opc_scatter_2.pdf",  width = 15, height = 20)
ggarrange(mr_scatter$`met-c-842.D55bOj`,mr_scatter$`met-c-843.D55bOj`, ncol=2, nrow=2, widths = 2, heights = 1)
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
write.csv(results_hnc,"./mrresults_allSNPs_opc.csv")

#6. PLOT FOREST AND LEAVE-ONE-OUT 
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(dat_hnc)
res_single$b <- -1*(res_single$b)

mr_forest <- mr_forest_plot(res_single)
mr_forest

pdf("allSNPs_opc_forest.pdf",  width = 15, height = 20)
ggarrange(mr_forest$`ieu-a-299.D55bOj`,mr_forest$`ieu-a-300.D55bOj`,mr_forest$`ieu-a-301.D55bOj`,mr_forest$`ieu-a-302.D55bOj`, ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

pdf("allSNPs_opc_forest_2.pdf",  width = 15, height = 20)
ggarrange(mr_forest$`met-c-842.D55bOj`,mr_forest$`met-c-843.D55bOj`, ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(dat_hnc)
mr_loo <- mr_leaveoneout_plot(res_loo)
mr_loo

pdf("allSNPs_opc_loo.pdf",  width = 15, height = 20)
ggarrange(mr_loo$`ieu-a-299.D55bOj`,mr_loo$`ieu-a-300.D55bOj`,mr_loo$`ieu-a-301.D55bOj`,mr_loo$`ieu-a-302.D55bOj`, ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

pdf("allSNPs_opc_loo_2.pdf",  width = 15, height = 20)
ggarrange(mr_loo$`met-c-842.D55bOj`,mr_loo$`met-c-843.D55bOj`, ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()