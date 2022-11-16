# Some packages being used
library(data.table)
library(dplyr)
library(tidyr)
library(cluster)
library(tidyverse)
library(ggplot2)

tpm_hipp = fread("D:/GTEx_SS/Data/GCT/gene_tpm_brain_hippocampus.gct", header = T)
samp_IDs = colnames(tpm_hipp)
ID_pattern = "GTEX-[a-zA-Z0-9]+|Name|Description"
new_IDs = str_match(samp_IDs,ID_pattern)
colnames(tpm_hipp) = new_IDs

tpm_amy = fread("D:/GTEx_SS/Data/GCT/gene_tpm_brain_amygdala.gct", header = T)
samp_IDs = colnames(tpm_amy)
ID_pattern = "GTEX-[a-zA-Z0-9]+|Name|Description"
new_IDs = str_match(samp_IDs,ID_pattern)
colnames(tpm_amy) = new_IDs

tpm_cor = fread("D:/GTEx_SS/Data/GCT/gene_tpm_brain_cortex.gct", header = T)
samp_IDs = colnames(tpm_cor)
ID_pattern = "GTEX-[a-zA-Z0-9]+|Name|Description"
new_IDs = str_match(samp_IDs,ID_pattern)
colnames(tpm_cor) = new_IDs

tpm_pfc = fread("D:/GTEx_SS/Data/GCT/gene_tpm_brain_frontal_cortex.gct", header = T)
samp_IDs = colnames(tpm_pfc)
ID_pattern = "GTEX-[a-zA-Z0-9]+|Name|Description"
new_IDs = str_match(samp_IDs,ID_pattern)
colnames(tpm_pfc) = new_IDs

tpm_sn = fread("D:/GTEx_SS/Data/GCT/gene_tpm_brain_substantia_nigra.gct", header = T)
samp_IDs = colnames(tpm_sn)
ID_pattern = "GTEX-[a-zA-Z0-9]+|Name|Description"
new_IDs = str_match(samp_IDs,ID_pattern)
colnames(tpm_sn) = new_IDs

tpm_hyp = fread("D:/GTEx_SS/Data/GCT/gene_tpm_brain_hypothalamus.gct", header = T)
samp_IDs = colnames(tpm_hyp)
ID_pattern = "GTEX-[a-zA-Z0-9]+|Name|Description"
new_IDs = str_match(samp_IDs,ID_pattern)
colnames(tpm_hyp) = new_IDs

tpm_na = fread("D:/GTEx_SS/Data/GCT/gene_tpm_brain_nucleus_accumbens.gct", header = T)
samp_IDs = colnames(tpm_na)
ID_pattern = "GTEX-[a-zA-Z0-9]+|Name|Description"
new_IDs = str_match(samp_IDs,ID_pattern)
colnames(tpm_na) = new_IDs



HIPP_file = "D:/GTEx_SS/Data/Raw/Brain_Hippocampus.v8.covariates.txt"
HIPP = fread(file = HIPP_file,header = T)
HIPP = HIPP[c(1:(nrow(HIPP)-3)),]

SN_file = "D:/GTEx_SS/Data/Raw/Brain_Substantia_nigra.v8.covariates.txt"
SN = fread(file = SN_file,header = T)
SN = SN[c(1:(nrow(SN)-3)),]

#Loading Co-variates for Amygdala
AM_file = "D:/GTEx_SS/Data/Raw/Brain_Amygdala.v8.covariates.txt"
AM = fread(file = AM_file,header = T)
AMG = AM[c(1:(nrow(AM)-3)),]

#Loading Co-variates for Cortex
COR_file = "D:/GTEx_SS/Data/Raw/Brain_Cortex.v8.covariates.txt"
COR = fread(file = COR_file,header = T)
COR = COR[c(1:(nrow(COR)-3)),]

#Loading Co-variates for pre frontal Cortex
PFC_file = "D:/GTEx_SS/Data/Raw/Brain_Frontal_Cortex_BA9.v8.covariates.txt"
PFC = fread(file = PFC_file,header = T)
PFC = PFC[c(1:(nrow(PFC)-3)),]

#Loading Co-variates for hypothalamus
HYP_file = "D:/GTEx_SS/Data/Raw/Brain_Hypothalamus.v8.covariates.txt"
HYP = fread(file = HYP_file,header = T)
HYP = HYP[c(1:(nrow(HYP)-3)),]

#Loading Co-variates for Nucleus Accumbens
NACC_file = "D:/GTEx_SS/Data/Raw/Brain_Cortex.v8.covariates.txt"
NACC= fread(file = NACC_file,header = T)
NACC = NACC[c(1:(nrow(NACC)-3)),]


reg_file = "D:/GTEx_SS/Data/Raw/Regulation of Gene Expression.csv"
gene_regs = read.csv(reg_file, header = T)
HIPP_genes = gene_regs %>% filter(Region.Code == 'HIPP') %>% select(Genes)
HIPP_genes = intersect(toupper(HIPP_genes$Genes),tpm_hipp$Description)

SN_genes = gene_regs %>% filter(Region.Code == 'SN') %>% select(Genes)
SN_genes = intersect(toupper(SN_genes$Genes),tpm_sn$Description)

# Amygdala
AM_genes = gene_regs %>% filter(Region.Code == 'AMG') %>% select(Genes)
AM_genes = intersect(toupper(AM_genes$Genes),tpm_amy$Description)

# Cortex
COR_genes = gene_regs %>% filter(Region.Code == 'COR') %>% select(Genes)
COR_genes = intersect(toupper(COR_genes$Genes),tpm_cor$Description)

# Prefrontal Cortex
PFC_genes = gene_regs %>% filter(Region.Code == 'PFC') %>% select(Genes)
PFC_genes = intersect(toupper(PFC_genes$Genes),tpm_pfc$Description)

# Hypothalamus
HYP_genes = gene_regs %>% filter(Region.Code == 'HYP') %>% select(Genes)
HYP_genes = intersect(toupper(HYP_genes$Genes),tpm_hyp$Description)

# Nucleus Accumbens
NACC_genes = gene_regs %>% filter(Region.Code == 'NACC') %>% select(Genes)
NACC_genes = intersect(toupper(NACC_genes$Genes),tpm_na$Description)


find_correlation = function(gene,brain_part,tpm){
  
  cols = intersect(colnames(brain_part),colnames(tpm))
  new_tpm = tpm%>% select("Description",cols) 
  gene_vals = new_tpm %>% filter(Description == gene)
  tpm_vals = t(gene_vals)
  tpm_vals = tpm_vals[-1]
  tpm_vals = as.numeric(tpm_vals)
  
  brain_part = brain_part %>% select("ID",cols)
  cor_vec = c()
  for(i in seq(1:(nrow(brain_part)))){
    pc_cov_vals = t(brain_part[i,])
    cov_pc = pc_cov_vals[1]
    pc_cov_vals = pc_cov_vals[-1]
    pc_cov_vals = as.numeric(pc_cov_vals)
    cor_vec[i] = cor(pc_cov_vals,tpm_vals)
  }
  names(cor_vec) = brain_part$ID
  return(cor_vec)
}

write_correlation = function(genes,cov_mat,part,tpm){
  
  rows = nrow(cov_mat)
  cols = length(genes)
  cor_df = data.frame(matrix(nrow = rows)) 
  for (i in genes){
    cor_df = cbind(cor_df,find_correlation(i,cov_mat,tpm))
    
  }
  cor_df = cor_df[,-1]
  colnames(cor_df) = genes
  file_name = paste('D:/GTEx_SS/Results/',part,'.csv')
  fwrite(cor_df, file = file_name)
  
  cor_df %>% slice_max()
}

write_correlation(HIPP_genes,cov_mat = HIPP, part = 'Hippocampus', tpm = tpm_hipp)
write_correlation(AM_genes, AMG, part = 'Amygdala', tpm_amy)
write_correlation(COR_genes, COR, part = 'Cortex', tpm_cor)
write_correlation(PFC_genes, PFC, part = 'PFC', tpm_pfc)
write_correlation(NACC_genes, NACC, part = 'Nucleus Accumbens', tpm_na)
write_correlation(HIPP_genes,cov_mat = HIPP, part = 'Substantia nigra', tpm_sn)
