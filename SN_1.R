# Some packages being used
library(data.table)
library(dplyr)
library(tidyr)
library(cluster)
library(tidyverse)
library(readxl)
library(FactoMineR)
library(factoextra)


# Loading the gene expression matrix
TPM_brain_file = "D:/GTEx_SS/Data/GTEx_Brain2.txt"
tpm_brain = fread(file = TPM_brain_file, header = T)
# Renaming columns to shorter form 
samp_IDs = colnames(tpm_brain)
ID_pattern = "GTEX-[a-zA-Z0-9]+|Name|Description"
new_IDs = str_match(samp_IDs,ID_pattern)
colnames(tpm_brain) = new_IDs

# Loading the co-variates matrix for Substantia Nigra
SN_file = "D:/GTEx_SS/Data/Raw/Brain_Substantia_nigra.v8.covariates.txt"
Sub_nigra = fread(file = SN_file,header = T)
# Taking transpose to obtain variables as columns
transpose_Sub_nigra = t(as.data.frame(Sub_nigra))
transpose_Sub_nigra = as.data.frame(transpose_Sub_nigra)
colnames(transpose_Sub_nigra) = as.character(transpose_Sub_nigra[1,])
transpose_Sub_nigra = transpose_Sub_nigra[-1,]

#loading the gene regulation matrix
reg_file = "D:/GTEx_SS/Data/Raw/Regulation of Gene Expression.csv"
gene_regs = read.csv(reg_file, header = T)
# Selecting Substantia Nigra
SN_genes = gene_regs %>% filter(Region.Code == 'SN') %>% select(Genes)
genes = toupper(SN_genes$Genes)

# Seperating the required columns from PEER matrix, PC's, PEERs, sexes
PC_matrix = apply(transpose_Sub_nigra[,c(1:5)],2,as.numeric)
row.names(PC_matrix) = row.names(transpose_Sub_nigra)

cov_matrix = apply(transpose_Sub_nigra[,c(6:20)],2,as.numeric)
row.names(cov_matrix) = row.names(transpose_Sub_nigra)

PC_Cov_matrix = apply(transpose_Sub_nigra[,c(1:20)],2,as.numeric)
row.names(PC_Cov_matrix) = row.names(transpose_Sub_nigra)

sex_matrix = transpose_Sub_nigra %>% select(sex) %>%
  mutate(sex= as.numeric(sex))
# ID's of males and females
males = row.names(sex_matrix %>% filter(sex == 1))
females = row.names(sex_matrix %>% filter(sex == 2))

# K means on all matrices
km_PC = kmeans(PC_matrix, 2, nstart= 20)
km_Cov = kmeans(cov_matrix,2,20)
km_PC_Cov = kmeans(PC_Cov_matrix,2,20)

# Plotting clusters and within cluster sum of squares (elbow method) for PC's, Cov's and PC + Cov's
fviz_nbclust(PC_matrix, kmeans, nstart=100, method = "wss")
fviz_cluster(km_PC, PC_matrix,palette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
                          geom = "point",
                          ellipse.type = "convex", 
                          ggtheme = theme_bw()
             )

fviz_nbclust(cov_matrix, kmeans, nstart=20, method = "wss")
fviz_cluster(km_Cov, cov_matrix,palette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)


fviz_nbclust(PC_Cov_matrix, kmeans, nstart=20, method = "wss") 
fviz_cluster(km_PC_Cov, PC_Cov_matrix,palette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)

# Sample run with 5 clusters (based on elbow plot) and plotting the same
km_PC_2 = kmeans(PC_matrix, 5, nstart= 20)      
fviz_nbclust(PC_matrix, kmeans, nstart=20, method = "wss") 
fviz_cluster(km_PC_2, PC_matrix,palette = c("#2E9FDF", "#00AFBB", "#E7B800",'Blue','Green'), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)

#Finding gender distribution and naming clusters as smokers and nonsmokers
gender_dist = function(km, males, females){
  cluster_1 = names(km$cluster[km$cluster==1])
  cluster_2 = names(km$cluster[km$cluster==2])
  
  males_in_cluster_1 = intersect(cluster_1,males)
  males_in_cluster_2 = intersect(cluster_2,males)
  females_in_cluster_1 = intersect(cluster_1,females)
  females_in_cluster_2 = intersect(cluster_2,females)
  
  males_ratio_1 = length(males_in_cluster_1) / length(cluster_1)
  males_ratio_2 = length(males_in_cluster_2) / length(cluster_2)
  females_ratio_1 = length(males_in_cluster_1) / length(cluster_1)
  females_ratio_2 = length(males_in_cluster_2) / length(cluster_2)
  
  if(males_ratio_1 > males_ratio_2){
    smokers = cluster_1
    non_smokers = cluster_2
  } else if(males_ratio_2 > males_ratio_1){
    smokers = cluster_2
    non_smokers = cluster_1
  }
  return(list(smoke = smokers, nosmoke = non_smokers))
}

PC_cluster = gender_dist(km_PC,males,females)
Cov_cluster = gender_dist(km_Cov,males,females)
PC_Cov_cluster = gender_dist(km_PC_Cov,males,females)

direction_test =  function(cluster, gene){
  goi_expression_s = tpm_brain %>% filter(Description == gene) %>% 
    select(intersect(cluster$smoke,colnames(tpm_brain)))
  s_mean = mean(as.numeric(goi_expression_s))
  s_median = median(as.numeric(goi_expression_s))
  goi_expression_ns = tpm_brain %>% filter(Description == gene) %>% 
    select(intersect(cluster$nosmoke,colnames(tpm_brain)))
  ns_mean = mean(as.numeric(goi_expression_ns))
  ns_median = median(as.numeric(goi_expression_ns))
  return(c(s_mean,s_median,ns_mean,ns_median))
}

direction_df = data.frame(Measures = c("Smoker_Mean", "Smoker_Median",
                             "Non_smoker_Mean", 
                             "Non_smoker_Median"))
for(i in genes){
  PC = direction_test(PC_cluster,i)
  Cov = direction_test(Cov_cluster,i)
  PC_Cov = direction_test(PC_Cov_cluster,i)
  df = data.frame(PC,Cov,PC_Cov)
  colnames(df) = c(paste('PC',i),paste('Cov',i),paste('PC_Cov',i))
  direction_df = cbind(direction_df,df)
}

direction_df = t(direction_df)
colnames(direction_df) = as.character(direction_df[1,])
direction_df = direction_df[-1,]
mean_diff = as.numeric(direction_df[,1]) - as.numeric(direction_df[,3])
median_diff = as.numeric(direction_df[,2]) - as.numeric(direction_df[,4])
direction_df = cbind(direction_df,mean_diff,median_diff)

# try kmeans on the tpm matrices brain part wise
# in directionality test, try comparison of means and median groupwise for tpm matrices for each brain part

