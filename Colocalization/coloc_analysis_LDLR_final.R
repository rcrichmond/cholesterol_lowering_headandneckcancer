########################## Execute coloc on chosen snps ##################

setwd("")


#clear ws 
rm(list = ls())

#Install and load packages 
install.packages("devtools")
devtools::install_github("MRCIEU/TwoSampleMR") #to update the package
devtools::install_github("MRCIEU/MRInstruments")
install.packages("plyr")
install.packages("ggplot2")
install.packages("xlsx")
install.packages("png")
install.packages("ggrepel")
install.packages("ggthemes")
install.packages("calibrate")
install.packages("vctrs")
BiocManager::install("snpStats")
install.packages("usethis")
devtools::install_github("mrcieu/ieugwasr")
install.packages("coloc")
install.packages("gwasglue")
install.packages("tidyverse")

BiocManager::install("biomaRt")
BiocManager::install("gwas")
BiocManager::install("rlang")
BiocManager::install("gwasgluev@4.0.2")
BiocManager::install("gwasvcf")
BiocManager::install("GenomicRanges")
BiocManager::install("gwasglue")

install.packages("gwasvcf")

BiocManager::install("Rhtslib") # use this one - 1
BiocManager::install("Rsamtools") # use this one - 2
BiocManager::install("VariantAnnotation") # use this one - 3

devtools::install_github("mrcieu/gwasvcf") # this one after 1,2,3 above
devtools::install_github("mrcieu/gwasglue") # This one after gwasvcf
devtools::install_github("mrcieu/gwasglue")

library(gwasglue)
library(data.table)
library(gwasvcf)
library(biomaRt)
library(Rhtslib)
library(GenomicRanges)
library(Matrix)
library(survival)
library(calibrate)
library(ggrepel)
library(ggthemes)
library(devtools)
library(TwoSampleMR)
library(MRInstruments)
library(plyr) 
library(ggplot2)
library(png)
library(data.table)
library(dplyr)
library(stringr)
library(rlang)
library(gwasglue)
library(snpStats)
library(usethis)
library(remotes)
library(arrow)
library(coloc)
library(tidyverse)
library("data.table")

#authorise use of MRBase tools if variables taken from here
ieugwasr::get_access_token()
ao<-available_outcomes() 

#Install the renviron package from github, it is large so you need an access token
# #creat access token 
# usethis::browse_github_pat()
# #add it to the renviron using instructions at this webpage https://happygitwithr.com/github-pat.html
# usethis::edit_r_environ()
# #check the access token is there 
# github_token()


# ####### COLOCALISATION OF PCSK9/ LDLR WITH LDL-C AND ORAL/OROPHARYNGEAL CANCER SUBTYPES #######

setwd("")

####################################################################
## Read in regional snps data #
####################################################################

exposure <- as.data.frame(fread("regional_exposure_rs6511720.txt"))

## Outcome - read in for overall OC/OPC ##

outcome <- as.data.frame(fread("regional_outcome_rs6511720.txt"))

# Edited SNP col- take only the snps which are in exposure
    
commonsnp <- outcome[outcome$SNP %in% exposure$SNP,]
dim(commonsnp)
comm_snp <- as.vector(commonsnp$SNP)


exposure <- exposure[exposure$SNP %in% comm_snp,]
outcome <- outcome[outcome$SNP %in% comm_snp,]
######################################################################################################
## Prep datasets for input into coloc.abf (coloc_analysis.R)
######################################################################################################

###Create dataset 1 from exposure data

## Variance = SE^2, create a varbeta column for exposure
exposure$varbeta <- ((exposure$se)^2)
exposure$samplesize<-188578

## create MAF exposure column ## 
exposure$MAF <- ifelse (exposure$EAF<=0.5, as.numeric(exposure$EAF),
                       ifelse(exposure$EAF>=0.5, as.numeric(1-(exposure$EAF)), NA))
is.na(exposure$MAF)

## Create dataset1 in format that can be used in coloc.abf function ## SNP, P, EAF, BETA, varbeta, type(quant or case control) and Number of people
dataset1 <- list(c(as.character(exposure$SNP)), c(as.numeric(exposure$P)), c(as.numeric(exposure$EAF)), c(as.numeric(exposure$MAF)), c(as.numeric(exposure$beta)),
                 c(as.numeric(exposure$varbeta)), c(as.numeric(exposure$se)), c("quant"), c(as.numeric(exposure$samplesize[[1]])))

names(dataset1) <- c("snp", "pvalues", "EAF", "MAF", "beta", "varbeta", "se", "type", "N")

lengths(dataset1)

## Create datset2 from outcome data ##

## create MAF outcome column ## 
outcome$MAF <- ifelse (outcome$EAF<=0.5, as.numeric(outcome$EAF),
                       ifelse(outcome$EAF>=0.5, as.numeric(1-(outcome$EAF)), NA))
is.na(outcome$MAF)

# Create varbeta column for outcome #

outcome$varbeta <- ((outcome$se)^2)

## Need to know number of cases and number of controls and get proportion of cases ## need it to put into dataset 2.

controls<- 6436
cases<- 6034
# s is the proportion of samples that are cases
proportion_of_cases <- cases/ (cases+controls)

dataset2 <- list(c(as.character(outcome$SNP)), c(as.numeric(outcome$P)), c(as.numeric(outcome$EAF)), c(as.numeric(outcome$MAF)), c(as.numeric(outcome$beta)),
                 c(as.numeric(outcome$varbeta)), c(as.numeric(outcome$se)), c("cc"), c(as.numeric(0.48)), c(as.numeric(outcome$N[[1]])))

names(dataset2) <- c("snp", "pvalues", "EAF", "MAF", "beta", "varbeta", "se", "type", "s", "N")

lengths(dataset2)


######################################################################################################
## Run Coloc 
######################################################################################################

coloc_results <- coloc.abf(dataset1=list(pvalues=dataset1$pvalues,N=dataset1$N, MAF=dataset1$MAF,type="quant"), 
                           dataset2=list(pvalues=dataset2$pvalues,N=dataset2$N, MAF=dataset2$MAF,type="cc",s=0.48))


coloc_results$summary

head(coloc_results$results)

coloc_results_summary <- as.data.frame(coloc_results$summary)

coloc_results_summary

coloc_results_all <- as.data.frame(coloc_results$results)

head(coloc_results_all)
dim(coloc_results_all)

write.table(coloc_results_summary, "coloc_results_summary_rs6511720_final.txt", quote=F, row.names=T, sep="\t")

write.table(coloc_results_all, "coloc_results_all_rs6511720_final.txt", quote=F, row.names=F, sep="\t")

# needs to be <500 to create plot 
temp <- coloc_to_gassocplot(coloc_results_all)

plot <- gassocplot::stack_assoc_plot(temp$markers, temp$z, temp$corr, traits=temp$traits)
