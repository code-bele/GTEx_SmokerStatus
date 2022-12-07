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

SN_file = "D:/GTEx_SS/Data/Raw/Brain_Substantia_nigra.v8.covariates.txt"
Sub_nigra = fread(file = SN_file,header = T)
SN = Sub_nigra[c(1:(nrow(Sub_nigra)-3)),]
SN = t(SN)
colnames(SN) = SN[1,]
SN = SN[-1,]
rows = row.names(SN)
SN = apply(SN,2,as.numeric)
row.names(SN) = rows

model = kmeans(SN, 2, 20)
fviz_cluster(model, SN,geom = "point", ellipse.type = "convex")
heatmap(SN)
clusplot(SN, model$cluster)
model2 = pam(SN,2)
clusplot(SN,model2$clustering)

#Selective model with few covariates
SN2 = as.data.frame(SN)
SN2 = SN2 %>% select(InferredCov5, InferredCov13,
                    InferredCov14, InferredCov11, InferredCov4)
model3 = kmeans(SN2,2,20)
clusplot(SN2,model3$cluster)
