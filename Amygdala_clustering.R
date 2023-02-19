library(data.table)
library(dplyr)
library(tidyr)
library(cluster)
library(tidyverse)
library(readxl)
library(FactoMineR)
library(factoextra)
library(cluster)

AMG_file = "D:/GTEx_SS/Data/Raw/Brain_Amygdala.v8.covariates.txt"
amg = fread(file = AMG_file,header = T)
AMG = amg[c(1:(nrow(amg)-3)),]
AMG = t(AMG)
colnames(AMG) = AMG[1,]
AMG = AMG[-1,]
rows = row.names(AMG)
AMG = apply(AMG,2,as.numeric)
row.names(AMG) = rows

AMG2 = as.data.frame(AMG)

cluster_plot = function(covs){
  am = AMG2 %>% select(covs)
  model = kmeans(am,2,20)
  l = length(covs)
  save_as = paste0("D:/GTEx_SS/Results/Clustering_plots",
                   "/amygdala",as.character(l),".pdf")
  
  pdf(save_as)
  clusplot(am,model$cluster,
           main = paste("Amygdala", "-",
                        as.character(l),"covariates"))
  dev.off()
  
  return(model)
}


model = cluster_plot(covs = c(6:ncol(AMG2)))
model2 = cluster_plot(covs = c("InferredCov4", "InferredCov5",
                               "InferredCov6","InferredCov12", "InferredCov9"))
model3 = cluster_plot(covs = c(9,10,11,14))
gene_amy = c("FOS")
gene_reg = c("Up")

tpm_amg_file = "D:/GTEx_SS/Data/GCT/gene_tpm_brain_amygdala.gct"
tpm_amg = fread(tpm_amg_file, header =T)
samp_IDs = colnames(tpm_amg)
ID_pattern = "GTEX-[a-zA-Z0-9]+|Name|Description"
new_IDs = str_match(samp_IDs,ID_pattern)
colnames(tpm_amg) = new_IDs

sex = amg %>% filter(ID == "sex")

cluster_comparison = function( mode, gene, exp){
  cluster_1 = names(mode$cluster[mode$cluster==1])
  cluster_2 = names(mode$cluster[mode$cluster==2])
  result = c("Regulation","T test", "Mean of cluster 1", "Mean of cluster 2",
             "Wilcox test", "F test")
  
  
  for(i in gene){
    exp_1 = tpm_amg %>% filter(Description == i) %>% 
      select(intersect(cluster_1,colnames(tpm_amg)))
    exp_2 = tpm_amg %>% filter(Description == i) %>% 
      select(intersect(cluster_2,colnames(tpm_amg)))
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
  about = cbind(about_cluster, cl)
  result = cbind(about, result)
  colnames(result) = c("Cluster information", "Result",
                       "Cluster comparison metrics","Result")
  
  
  return(as.data.frame(result))
}

c1 = cluster_comparison(mode = model2, gene = gene_amy, exp = gene_reg)
View(c1)

fwrite(as.data.frame(c1),
       file = "D:/GTEx_SS/Results/Cluster_comparison_metrics/amygdala_model2_covs5.csv")

