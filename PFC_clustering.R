
PFC_file = "D:/GTEx_SS/Data/Raw/Brain_Frontal_Cortex_BA9.v8.covariates.txt"
pfc = fread(file = PFC_file,header = T)
PFC = pfc[c(1:(nrow(pfc)-3)),]
PFC = t(PFC)
colnames(PFC) = PFC[1,]
PFC = PFC [-1,]
rows = row.names(PFC)
PFC = apply(PFC,2,as.numeric)
row.names(PFC) = rows

PFC2 = as.data.frame(PFC)

cluster_plot = function(covs){
  pfc2 = PFC2 %>% select(covs)
  model = kmeans(pfc2,2,20)
  l = length(covs)
  save = paste0("D:/GTEx_SS/Results/Clustering_plots",
                   "/pre_frontal_cortex",as.character(l),".pdf")
  
  pdf(save)
  clusplot(pfc2,model$cluster,
           main = paste("Pre Frontal Cortex", "-",
                        as.character(l),"covariates"))
  dev.off()
  
  return(model)
}

model = cluster_plot(c(6:ncol(PFC2)))
model1 = cluster_plot(c(6,7,8))
model2 = cluster_plot(c(6,7,8,9,10))

tpm_pfc_file = "D:/GTEx_SS/Data/GCT/gene_tpm_brain_frontal_cortex.gct"
tpm_pfc = fread(tpm_pfc_file, header =T)
samp_IDs = colnames(tpm_pfc)
ID_pattern = "GTEX-[a-zA-Z0-9]+|Name|Description"
new_IDs = str_match(samp_IDs,ID_pattern)
colnames(tpm_pfc) = new_IDs

sex = pfc %>% filter(ID == "sex")
gene = c("DDN", "SIRT1", "FOS", "BDNF", "CRH")
regs = c("UP","DOWN","UP","UP","UP")

cluster_comparison = function(tpm, mode, genes, exp){
  cluster_1 = names(mode$cluster[mode$cluster==1])
  cluster_2 = names(mode$cluster[mode$cluster==2])
  result = c("Regulation", "T test", "Mean of cluster 1", "Mean of cluster 2",
             "Wilcox test", "F test")
  for(i in genes){
    exp_1 = tpm_pfc %>% filter(Description == i) %>% 
      select(intersect(cluster_1,colnames(tpm_pfc)))
    exp_2 = tpm_pfc %>% filter(Description == i) %>% 
      select(intersect(cluster_2,colnames(tpm_pfc)))
    t = t.test(exp_1,exp_2)
    v = var.test(as.numeric(t(exp_1)[,1]),as.numeric(t(exp_2)[,1]))
    w = wilcox.test(as.numeric(t(exp_1)[,1]), as.numeric(t(exp_2)[,1]))
    loc = which(genes == i)
    gen = c(exp[loc],t$p.value, t$estimate, w$p.value, v$p.value)
    result = cbind(result, gen)
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


k1 = cluster_comparison(tpm_pfc,model1,gene, exp =regs)
k2 = cluster_comparison(tpm_pfc,model2,gene, exp = regs)

fwrite(as.data.frame(k2), 
       file = "D:/GTEx_SS/Results/Cluster_comparison_metrics/PFC_model2_covs5.csv")

fwrite(as.data.frame(k1),
       file = "D:/GTEx_SS/Results/Cluster_comparison_metrics/PFC_model1_covs3.csv")






