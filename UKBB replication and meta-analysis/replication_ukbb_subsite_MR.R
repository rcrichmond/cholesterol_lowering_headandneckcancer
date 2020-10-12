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
install.packages("meta")

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
library(meta)

#Set working directory to source file location
setwd("")

#1. SELECT SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("./exposure_data.csv", sep=",")
head(drug_exp_data)
dim(drug_exp_data)

########################################################
#UKBB ORAL AND OROPHARYNGEAL HNC COMBINED (NO LARYNX)  #
########################################################

#2. SELECT SNP-OUTCOME SUMMARY DATA FROM UKBB HNC GWAS
drug_hnc <- read.csv("./hnc_snps_out_ukbb.csv", header=TRUE)
drug_hnc <- read_outcome_data("hnc_snps_out_ukbb.csv", drug_exp_data$SNP, sep=",")

#3. HARMONIZE DATASETS
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
dat_hnc <- harmonise_data(drug_exp_data, drug_hnc, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME e.g. HNC
mr_results_hnc <- mr(dat_hnc)
mr_results_hnc

#Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results_hnc, dat_hnc)
mr_scatter

pdf("hnc_ukbb_scatter.pdf",  width = 15, height = 20)
ggarrange(mr_scatter$`1d3qRO.ojwwbg`,mr_scatter$`BXr0AZ.ojwwbg`,mr_scatter$`DHssei.ojwwbg`,mr_scatter$`fF1czi.ojwwbg`,mr_scatter$`Jm4sXl.ojwwbg`,ncol=2, nrow=2, widths = 2, heights = 1)
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
write.csv(results_hnc,"./hnc_ukbb_results.csv")

#6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(dat_hnc)
res_single$b <- -1*(res_single$b)
mr_forest <- mr_forest_plot(res_single)
mr_forest

pdf("hnc_ukbb_forest.pdf",  width = 15, height = 20)
ggarrange(mr_forest$`1d3qRO.ojwwbg`,mr_forest$`BXr0AZ.ojwwbg`,mr_forest$`DHssei.ojwwbg`,mr_forest$`fF1czi.ojwwbg`,mr_forest$`Jm4sXl.ojwwbg`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(dat_hnc)
mr_loo <- mr_leaveoneout_plot(res_loo)

pdf("hnc_ukbb_loo.pdf",  width = 15, height = 20)
ggarrange(mr_loo$`1d3qRO.ojwwbg`,mr_loo$`BXr0AZ.ojwwbg`,mr_loo$`DHssei.ojwwbg`,mr_loo$`fF1czi.ojwwbg`,mr_loo$`Jm4sXl.ojwwbg`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()


#7. ASSESS HETEROGENEITY AND PLEIOTROPY 

#Hetereogenity test
mr_heterogeneity(dat_hnc)
write.csv(mr_heterogeneity(dat_hnc),"./hnc_heterogeneity_results_ukbb.csv")

#Pleiotropy test
mr_pleiotropy_test(dat_hnc)
write.csv(mr_pleiotropy_test(dat_hnc),"./hnc_pleiotropy_results_ukbb.csv")

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

write.csv(mr_presso_global, "./hnc_ukbb_mr_presso_global.csv")

########################################################
#UKBB ORAL CANCER                                      #
########################################################
drug_hnc <- read.csv("./oc_snps_out_ukbb.csv", header=TRUE)
drug_hnc <- read_outcome_data("oc_snps_out_ukbb.csv", drug_exp_data$SNP, sep=",")

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
dat_hnc <- harmonise_data(drug_exp_data, drug_hnc, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME e.g. HNC
mr_results_hnc <- mr(dat_hnc)
mr_results_hnc

#Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results_hnc, dat_hnc)
mr_scatter

pdf("oc_ukbb_scatter.pdf",  width = 15, height = 20)
ggarrange(mr_scatter$`6rt0yh.ZGCsAt`,mr_scatter$`cx8FFm.ZGCsAt`,mr_scatter$`r6FySU.ZGCsAt`,mr_scatter$`wNpvB4.ZGCsAt`,mr_scatter$`ZbmyjK.ZGCsAt`,ncol=2, nrow=2, widths = 2, heights = 1)
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
write.csv(results_hnc,"./oc_ukbb_results.csv")

#6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(dat_hnc)
res_single$b <- -1*(res_single$b)

mr_forest <- mr_forest_plot(res_single)
mr_forest

pdf("oc_ukbb_forest.pdf",  width = 15, height = 20)
ggarrange(mr_forest$`6rt0yh.ZGCsAt`,mr_forest$`cx8FFm.ZGCsAt`,mr_forest$`r6FySU.ZGCsAt`,mr_forest$`wNpvB4.ZGCsAt`,mr_forest$`ZbmyjK.ZGCsAt`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(dat_hnc)
mr_loo <- mr_leaveoneout_plot(res_loo)

pdf("oc_ukbb_loo.pdf",  width = 15, height = 20)
ggarrange(mr_loo$`6rt0yh.ZGCsAt`,mr_loo$`cx8FFm.ZGCsAt`,mr_loo$`r6FySU.ZGCsAt`,mr_loo$`wNpvB4.ZGCsAt`,mr_loo$`ZbmyjK.ZGCsAt`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#7. ASSESS HETEROGENEITY AND PLEIOTROPY
#Hetereogenity test
mr_heterogeneity(dat_hnc)
write.csv(mr_heterogeneity(dat_hnc),"./oc_heterogeneity_results_ukbb.csv")

#Pleiotropy test
mr_pleiotropy_test(dat_hnc)
write.csv(mr_pleiotropy_test(dat_hnc),"./oc_pleiotropy_results_ukbb.csv")

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

write.csv(mr_presso_global, "./oc_ukbb_mr_presso_global.csv")

########################################################
#UKBB OROPHARYNGEAL CANCER                             #
########################################################
drug_hnc <- read.csv("./opc_snps_out_ukbb.csv", header=TRUE)
drug_hnc <- read_outcome_data("opc_snps_out_ukbb.csv", drug_exp_data$SNP, sep=",")

#3. HARMONIZE DATASETS                                
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
dat_hnc <- harmonise_data(drug_exp_data, drug_hnc, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME e.g. HNC
mr_results_hnc <- mr(dat_hnc)
mr_results_hnc

#Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results_hnc, dat_hnc)
mr_scatter

pdf("opc_ukbb_scatter.pdf",  width = 15, height = 20)
ggarrange(mr_scatter$`FVzsrn.S6Qog1`,mr_scatter$`IerR7p.S6Qog1`,mr_scatter$`SsKOnD.S6Qog1`,mr_scatter$`Zb9afg.S6Qog1`,mr_scatter$`ZX4fmm.S6Qog1`,ncol=2, nrow=2, widths = 2, heights = 1)
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
write.csv(results_hnc,"./op _ukbb_results.csv")

#6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(dat_hnc)
res_single$b <- -1*(res_single$b)

mr_forest <- mr_forest_plot(res_single)
mr_forest

pdf("opc_ukbb_forest.pdf",  width = 15, height = 20)
ggarrange(mr_forest$`FVzsrn.S6Qog1`,mr_forest$`IerR7p.S6Qog1`,mr_forest$`SsKOnD.S6Qog1`,mr_forest$`Zb9afg.S6Qog1`,mr_forest$`ZX4fmm.S6Qog1`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(dat_hnc)
mr_loo <- mr_leaveoneout_plot(res_loo)

pdf("opc_ukbb_loo.pdf",  width = 15, height = 20)
ggarrange(mr_loo$`FVzsrn.S6Qog1`,mr_loo$`IerR7p.S6Qog1`,mr_loo$`SsKOnD.S6Qog1`,mr_loo$`Zb9afg.S6Qog1`,mr_loo$`ZX4fmm.S6Qog1`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#7. ASSESS HETEROGENEITY AND PLEIOTROPY
#Hetereogenity test
mr_heterogeneity(dat_hnc)
write.csv(mr_heterogeneity(dat_hnc),"./opc_heterogeneity_results_ukbb.csv")

#Pleiotropy test
mr_pleiotropy_test(dat_hnc)
write.csv(mr_pleiotropy_test(dat_hnc),"./opc_pleiotropy_results_ukbb.csv")

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

write.csv(mr_presso_global, "./opc_ukbb_mr_presso_global.csv")

########################################################
#UKBB LARYNGEAL CANCER                                 #
########################################################
drug_hnc <- read.csv("./larynx_snps_out_ukbb.csv", header=TRUE)
drug_hnc <- read_outcome_data("larynx_snps_out_ukbb.csv", drug_exp_data$SNP, sep=",")

#3. HARMONIZE DATASETS
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
dat_hnc <- harmonise_data(drug_exp_data, drug_hnc, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME e.g. HNC
mr_results_hnc <- mr(dat_hnc)
mr_results_hnc

#Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results_hnc, dat_hnc)
mr_scatter

pdf("larynx_ukbb_scatter.pdf",  width = 15, height = 20)
ggarrange(mr_scatter$`JrPx2h.m615jM`,mr_scatter$`nu1lX1.m615jM`,mr_scatter$`RV8V41.m615jM`,mr_scatter$`tZq2Rc.m615jM`,mr_scatter$`wSAeB2.m615jM`,ncol=2, nrow=2, widths = 2, heights = 1)
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
write.csv(results_hnc,"./larynx_ukbb_results.csv")

#6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(dat_hnc)
res_single$b <- -1*(res_single$b)

mr_forest <- mr_forest_plot(res_single)
mr_forest

pdf("larynx_ukbb_forest.pdf",  width = 15, height = 20)
ggarrange(mr_forest$`JrPx2h.m615jM`,mr_forest$`nu1lX1.m615jM`,mr_forest$`RV8V41.m615jM`,mr_forest$`tZq2Rc.m615jM`,mr_forest$`wSAeB2.m615jM`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(dat_hnc)
mr_loo <- mr_leaveoneout_plot(res_loo)

pdf("larynx_ukbb_loo.pdf",  width = 15, height = 20)
ggarrange(mr_loo$`JrPx2h.m615jM`,mr_loo$`nu1lX1.m615jM`,mr_loo$`RV8V41.m615jM`,mr_loo$`tZq2Rc.m615jM`,mr_loo$`wSAeB2.m615jM`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#7. ASSESS HETEROGENEITY AND PLEIOTROPY
#Hetereogenity test
mr_heterogeneity(dat_hnc)
write.csv(mr_heterogeneity(dat_hnc),"./larynx_heterogeneity_results_ukbb.csv")

#Pleiotropy test
mr_pleiotropy_test(dat_hnc)
write.csv(mr_pleiotropy_test(dat_hnc),"./larynx_pleiotropy_results_ukbb.csv")

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

write.csv(mr_presso_global, "./larynx_ukbb_mr_presso_global.csv")

###############################################################################
#7. META-ANALYSIS FOR COMBINED ORAL AND OROPHARYNGEAL CANCER 
###############################################################################
GAMEON <- read.csv("gameon_hnc_results.csv", header=TRUE)
UKBiobank <- read.csv("hnc_ukbb_results.csv", header=TRUE)

#Separate exposure-method forest plots
for (exposure in GAMEON$mr_results_hnc.exposure){
  for (method in GAMEON$mr_results_hnc.method){
    GAMEON_beta <- GAMEON$mr_results_hnc.b [GAMEON$mr_results_hnc.exposure==exposure & GAMEON$mr_results_hnc.method==method]
    UKBiobank_beta <- UKBiobank$mr_results_hnc.b [UKBiobank$mr_results_hnc.exposure==exposure & UKBiobank$mr_results_hnc.method==method]
    TE <- c(GAMEON_beta, UKBiobank_beta) 
    GAMEON_se <- GAMEON$mr_results_hnc.se [GAMEON$mr_results_hnc.exposure==exposure & GAMEON$mr_results_hnc.method==method]
    UKBiobank_se <- UKBiobank$mr_results_hnc.se [UKBiobank$mr_results_hnc.exposure==exposure & UKBiobank$mr_results_hnc.method==method]
    seTE= c(GAMEON_se, UKBiobank_se )
    res <- metagen(TE, seTE, studlab= c("GAMEON", "UKBiobank"), sm="OR")
    pdf(paste0(exposure,"_",method,".pdf"))
    forest(res, leftcols="studlab", rightcols=c("effect", "ci")) 
    dev.off()
  }
}

#Combined exposure-method forest plots
GAMEON$studlab <- "GAMEON"
UKBiobank$studlab <- "UKBiobank"
meta <- rbind(GAMEON, UKBiobank)
meta$method <- meta$mr_results_hnc.method

for (exposure in meta$mr_results_hnc.exposure){
  exp <- meta[meta$mr_results_hnc.exposure==exposure,]
  res <- metagen(TE=exp$mr_results_hnc.b, seTE=exp$mr_results_hnc.se, studlab = exp$studlab,  byvar = exp$method, sm="OR", overall=F)
  png(paste0(exposure,".png"), width=600, height=800)
  forest(res, leftcols="studlab", rightcols=c("effect", "ci")) 
  dev.off()
}

#IVW only forest plots
GAMEON$studlab <- "GAMEON"
UKBiobank$studlab <- "UKBiobank"
meta <- rbind(GAMEON, UKBiobank)
meta$method <- meta$mr_results_hnc.method=="Inverse variance weighted"

for (exposure in meta$mr_results_hnc.exposure){
  exp <- meta[meta$mr_results_hnc.exposure==exposure,]
  res <- metagen(TE=exp$mr_results_hnc.b, seTE=exp$mr_results_hnc.se, studlab = exp$studlab, byvar = exp$method, sm="OR", overall=F)
  png(paste0(exposure,"ivw.png"), width=600, height=800)
  forest(res, leftcols="studlab", rightcols=c("effect", "ci")) 
  dev.off()
}

#IVW only forest plots PCSK9
GAMEON <- read.csv("gameon_pcsk9_results.csv", header=TRUE)
UKBiobank <- read.csv("oral_pharyngeal_ukbb_pcsk9_results.csv", header=TRUE)
GAMEON$studlab <- "GAMEON"
UKBiobank$studlab <- "UKBiobank"
meta <- rbind(GAMEON, UKBiobank)
meta$method <- meta$mr_results_hnc.method=="Inverse variance weighted"

for (exposure in meta$mr_results_hnc.exposure){
  exp <- meta[meta$mr_results_hnc.exposure==exposure,]
  res <- metagen(TE=exp$mr_results_hnc.b, seTE=exp$mr_results_hnc.se, studlab = exp$studlab, byvar = exp$method, sm="OR", overall=F)
  png(paste0(exposure,"ivw.png"), width=600, height=800)
  forest(res, leftcols="studlab", rightcols=c("effect", "ci")) 
  dev.off()
}

#IVW only forest plots LDLR
GAMEON <- read.csv("gameon_ldlr_results.csv", header=TRUE)
UKBiobank <- read.csv("oral_pharyngeal_ukbb_ldlr_results.csv", header=TRUE)
GAMEON$studlab <- "GAMEON"
UKBiobank$studlab <- "UKBiobank"
meta <- rbind(GAMEON, UKBiobank)
meta$method <- meta$mr_results_hnc.method=="Inverse variance weighted"

for (exposure in meta$mr_results_hnc.exposure){
  exp <- meta[meta$mr_results_hnc.exposure==exposure,]
  res <- metagen(TE=exp$mr_results_hnc.b, seTE=exp$mr_results_hnc.se, studlab = exp$studlab, byvar = exp$method, sm="OR", overall=F)
  png(paste0(exposure,"ivw.png"), width=600, height=800)
  forest(res, leftcols="studlab", rightcols=c("effect", "ci")) 
  dev.off()
}

