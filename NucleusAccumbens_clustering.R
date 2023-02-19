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

NACC_file = "D:/GTEx_SS/Data/Raw/Brain_Nucleus_accumbens_basal_ganglia.v8.covariates.txt"
nacc = fread(file = NACC_file,header = T)
NACC = nacc[c(1:(nrow(nacc)-3)),]
NACC = t(NACC)
colnames(NACC) = NACC[1,]
NACC = NACC [-1,]
rows = row.names(NACC)
NACC = apply(NACC,2,as.numeric)
row.names(NACC) = rows

NACC2 = as.data.frame(NACC)
gene_nacc = c("CREB1", "FOSB","TNF")
gene_reg = c("Up","Up","Up")

cluster_plot = function(covs){
  nacc = NACC2 %>% select(covs)
  model = kmeans(nacc,2,20)
  l = length(covs)
  save_as = paste0("D:/GTEx_SS/Results/Clustering_plots",
                   "/Nucleus accumbens",as.character(l),".pdf")
  
  pdf(save_as)
  clusplot(nacc,model$cluster,
           main = paste("Nucleus accumbens", "-",
                        as.character(l),"covariates"))
  dev.off()
  
  return(model)
}



model = cluster_plot(covs = c(6:ncol(NACC2)))
model2 = cluster_plot(covs = c(6,7,9,8,14,28))
model3 = cluster_plot(covs = c(6,7,8,9))


tpm_nacc_file = "D:/GTEx_SS/Data/GCT/gene_tpm_brain_nucleus_accumbens.gct"
tpm_nacc = fread(tpm_nacc_file, header =T)
samp_IDs = colnames(tpm_nacc)
ID_pattern = "GTEX-[a-zA-Z0-9]+|Name|Description"
new_IDs = str_match(samp_IDs,ID_pattern)
colnames(tpm_nacc) = new_IDs

sex = nacc %>% filter(ID == "sex")

cluster_comparison = function( mode, gene, exp){
  cluster_1 = names(mode$cluster[mode$cluster==1])
  cluster_2 = names(mode$cluster[mode$cluster==2])
  result = c("Regulation","T test", "Mean of cluster 1", "Mean of cluster 2",
             "Wilcox test", "F test")
  
  for(i in gene){
    exp_1 = tpm_nacc %>% filter(Description == i) %>% 
      select(intersect(cluster_1,colnames(tpm_nacc)))
    exp_2 = tpm_nacc %>% filter(Description == i) %>% 
      select(intersect(cluster_2,colnames(tpm_nacc)))
    t = t.test(exp_1,exp_2)
    v = var.test(as.numeric(t(exp_1)[,1]),as.numeric(t(exp_2)[,1]))
    w = wilcox.test(as.numeric(t(exp_1)[,1]), as.numeric(t(exp_2)[,1]))
    loc = which(gene == i)
    g = c(exp[loc],t$p.value, t$estimate, w$p.value, v$p.value)
    result = cbind(result, g)
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

p1 = cluster_comparison(mode = model2, gene = gene_nacc, exp = gene_reg)
p2 = cluster_comparison(mode = model3, gene = gene_nacc, exp = gene_reg)

fwrite(as.data.frame(p2), 
       file = "D:/GTEx_SS/Results/Cluster_comparison_metrics/nucleus_acc_model3_covs4.csv")

fwrite(as.data.frame(p1),
       file = "D:/GTEx_SS/Results/Cluster_comparison_metrics/nucleus_acc_model2_covs6.csv")




