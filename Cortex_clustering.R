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

COR_file = "D:/GTEx_SS/Data/Raw/Brain_Cortex.v8.covariates.txt"
cor = fread(file = COR_file,header = T)
COR = cor[c(1:(nrow(cor)-3)),]
COR = t(COR)
colnames(COR) = COR[1,]
COR = COR[-1,]
rows = row.names(COR)
COR = apply(COR,2,as.numeric)
row.names(COR) = rows

cluster_plot = function(covs){
  cor = COR2 %>% select(covs)
  model = kmeans(cor,2,20)
  l = length(covs)
  save_as = paste0("D:/GTEx_SS/Results/Clustering_plots",
                   "/cortex",as.character(l),".pdf")
  
  pdf(save_as)
  clusplot(cor,model$cluster,
           main = paste("Cortex", "-",
                        as.character(l),"covariates"))
  dev.off()
  
  return(model)
}

COR2 = as.data.frame(COR)
model = cluster_plot(covs = c(6:ncol(COR2)))
model2 = cluster_plot(covs = c("InferredCov1", "InferredCov2",
                               "InferredCov4", "InferredCov11", 
                               "InferredCov18", "InferredCov21"))
model3 = cluster_plot(covs = c(6,7,9,16))
model4 = cluster_plot(covs = c("InferredCov1", "InferredCov2",
                               "InferredCov4"))

gene_cor = c("FOS","CHRNA4")
gene_reg = c("Up", "Up")

tpm_cor_file = "D:/GTEx_SS/Data/GCT/gene_tpm_brain_cortex.gct"
tpm_cor = fread(tpm_cor_file, header =T)
samp_IDs = colnames(tpm_cor)
ID_pattern = "GTEX-[a-zA-Z0-9]+|Name|Description"
new_IDs = str_match(samp_IDs,ID_pattern)
colnames(tpm_cor) = new_IDs

sex = cor %>% filter(ID == "sex")

cluster_comparison = function(tpm, mode, gene,exp){
  cluster_1 = names(mode$cluster[mode$cluster==1])
  cluster_2 = names(mode$cluster[mode$cluster==2])
  result = c("Regulation","T test", "Mean of cluster 1", "Mean of cluster 2",
             "Wilcox test", "F test")
  for(i in gene){
    exp_1 = tpm_cor %>% filter(Description == i) %>% 
      select(intersect(cluster_1,colnames(tpm_cor)))
    exp_2 = tpm_cor %>% filter(Description == i) %>% 
      select(intersect(cluster_2,colnames(tpm_cor)))
    t = t.test(exp_1,exp_2)
    v = var.test(as.numeric(t(exp_1)[,1]),as.numeric(t(exp_2)[,1]))
    w = wilcox.test(as.numeric(t(exp_1)[,1]), as.numeric(t(exp_2)[,1]))
    loc = which(gene == i)
    genes = c(exp[loc],t$p.value, t$estimate, w$p.value, v$p.value)
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

r1 = cluster_comparison(tpm_cor,model3,gene_cor, gene_reg)
r2 = cluster_comparison(tpm_cor, model2, gene_cor, gene_reg)
r3 = cluster_comparison(tpm_cor, model4, gene_cor, gene_reg)

fwrite(as.data.frame(r2), 
       file = "D:/GTEx_SS/Results/Cluster_comparison_metrics/Cortex_model2_covs4.csv")

fwrite(as.data.frame(r1),
       file = "D:/GTEx_SS/Results/Cluster_comparison_metrics/Cortex_model3_covs6.csv")

fwrite(as.data.frame(r3),
       file = "D:/GTEx_SS/Results/Cluster_comparison_metrics/Cortex_model4covs6.csv")



