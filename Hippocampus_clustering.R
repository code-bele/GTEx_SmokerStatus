# Some packages being used
library(data.table)
library(dplyr)
library(tidyr)
library(cluster)
library(tidyverse)
library(readxl)
library(FactoMineR)
library(factoextra)
library(cluster)

HIPP_file = "D:/GTEx_SS/Data/Raw/Brain_Hippocampus.v8.covariates.txt"
hipp = fread(file = HIPP_file,header = T)
HIPP = hipp[c(1:(nrow(hipp)-3)),]
HIPP = t(HIPP)
colnames(HIPP) = HIPP[1,]
HIPP = HIPP [-1,]
rows = row.names(HIPP)
HIPP = apply(HIPP,2,as.numeric)
row.names(HIPP) = rows

HIPP2 = as.data.frame(HIPP)


cluster_plot = function(covs){
  hipp2 = HIPP2 %>% select(covs)
  model = kmeans(hipp2,2,20)
  l = length(covs)
  save = paste0("D:/GTEx_SS/Results/Clustering_plots",
                "/hippocampus",as.character(l),".pdf")
  
  pdf(save)
  clusplot(hipp2,model$cluster,
           main = paste("hippocampus", "-",
                        as.character(l),"covariates"))
  dev.off()
  
  return(model)
}


model1 = cluster_plot(c(6:ncol(HIPP2)))
model2 = cluster_plot(c(6,7,8))
model3 = cluster_plot(c(6,7,8,9))


tpm_hipp_file = "D:/GTEx_SS/Data/GCT/gene_tpm_brain_hippocampus.gct"
tpm_hipp = fread(tpm_hipp_file, header =T)
samp_IDs = colnames(tpm_hipp)
ID_pattern = "GTEX-[a-zA-Z0-9]+|Name|Description"
new_IDs = str_match(samp_IDs,ID_pattern)
colnames(tpm_hipp) = new_IDs

sex = hipp %>% filter(ID == "sex")
gene_hipp = c("FOS", "CRHR1","CRHR2","NRG3","CREB1")
gene_reg = c("Up","Up","Up","Up","Up")

cluster_comparison = function(tpm, mode, gene, reg){
  cluster_1 = names(mode$cluster[mode$cluster==1])
  cluster_2 = names(mode$cluster[mode$cluster==2])
  result = c("Regulation","T test", "Mean of cluster 1", "Mean of cluster 2",
             "Wilcox test", "F test")
  for(i in gene){
    exp_1 = tpm_hipp %>% filter(Description == i) %>% 
      select(intersect(cluster_1,colnames(tpm_hipp)))
    exp_2 = tpm_hipp %>% filter(Description == i) %>% 
      select(intersect(cluster_2,colnames(tpm_hipp)))
    t = t.test(exp_1,exp_2)
    v = var.test(as.numeric(t(exp_1)[,1]),as.numeric(t(exp_2)[,1]))
    w = wilcox.test(as.numeric(t(exp_1)[,1]), as.numeric(t(exp_2)[,1]))
    loc = which(gene == i)
    genes = c(reg[loc], t$p.value, t$estimate, w$p.value, v$p.value)
    result = cbind(result, genes)
  }
  
  
  about_cluster = c("M:F cluster 1", "M:F cluster 2", 
                    "size cluster 1", "size of cluster 2",
                    "Cluster size ratio","")
  sex_ratio_1 = table(t(sex%>% select(intersect(cluster_1, colnames(sex)))))
  sex_ratio_2 = table(t(sex%>% select(intersect(cluster_2, colnames(sex)))))
  cl = c(sex_ratio_1[1]/sex_ratio_1[2], sex_ratio_2[1]/sex_ratio_2[2],
         length(cluster_1), 
         length(cluster_2),length(cluster_1)/length(cluster_2),NA)
  
  colnames(result) = c("Result", gene)
  result = cbind( result, about_cluster)
  result = cbind(result,Values =  cl)
  
  return(as.data.frame(result))
}



m2 = cluster_comparison(tpm_hipp,model2,gene = gene_hipp, reg = gene_reg)
m3 = cluster_comparison(tpm_hipp,model3,gene = gene_hipp, reg = gene_reg)

fwrite(as.data.frame(m2), 
       file = "D:/GTEx_SS/Results/Cluster_comparison_metrics/Hippocampus_model2_covs3.csv")

fwrite(as.data.frame(m3),
       file = "D:/GTEx_SS/Results/Cluster_comparison_metrics/Hippocampus_model3_covs4.csv")





