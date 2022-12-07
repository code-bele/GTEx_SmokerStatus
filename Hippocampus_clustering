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
HIPP = fread(file = HIPP_file,header = T)
HIPP = HIPP[c(1:(nrow(HIPP)-3)),]
HIPP = t(HIPP)
colnames(HIPP) = HIPP[1,]
HIPP = HIPP [-1,]
rows = row.names(HIPP)
HIPP = apply(HIPP,2,as.numeric)
row.names(HIPP) = rows


model = kmeans(HIPP, 2, 20)
clusplot(HIPP,model$cluster)
heatmap(HIPP)

# Selective covariates being modelled
HIPP2 = as.data.frame(HIPP)
HIPP2 = HIPP2 %>% select(InferredCov1, InferredCov2,
                     InferredCov4, InferredCov5, InferredCov7)
model3 = kmeans(HIPP2,2,20)
clusplot(HIPP2,model3$cluster)
