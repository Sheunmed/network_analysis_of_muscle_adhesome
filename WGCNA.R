if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("preprocessCore", "impute" ))

library(dplyr)
library(GEOquery)
library(tidyverse) 
library(DESeq2)
library(data.table)
library(WGCNA)
library(CorLevelPlot)
library(ggplot2)
library(gridExtra) 
library(dynamicTreeCut)
options(stringsAsFactors = FALSE) 
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggupset)
library(limma) 


## tranpose data - swap genes to columns and samples to rows-----
# confirm number of genes present across all datasets with reduce(intersect) 

trans_111006 <- t(vsd_111006_df)
trans_111010 <- t(vsd_111010_df)
trans_111016 <- t(vsd_111016_df)
trans_113165 <- t(vsd_113165_df)
trans_120642 <- t(vsd_120642_df)
trans_126865 <- t(vsd_126865_df)
trans_144304 <- t(vsd_144304_df)
trans_151066 <- t(vsd_151066_df)
trans_159217 <- t(vsd_159217_df)
trans_163434 <- t(vsd_163434_df)
trans_164471 <- t(vsd_164471_df)
trans_165630 <- t(vsd_165630_df)
trans_174106 <- t(vsd_174106_df)
trans_175495 <- t(vsd_175495_df)
trans_226151 <- t(vsd_226151_df)
trans_235781 <- t(vsd_235781_df)
trans_242202 <- t(vsd_242202_df)

trans_list <- list(trans_111006, trans_111010, trans_111016, trans_113165, trans_120642, trans_126865,
     trans_144304, trans_151066, trans_159217, trans_163434, trans_164471, trans_165630,
     trans_174106, trans_175495, trans_226151, trans_235781, trans_242202) 

trans_list <- lapply(trans_list, na.omit)

inter_trans <- Reduce(intersect, lapply(trans_list, colnames)) 
str(inter_trans)



 


## choose a set of threshold powers -----
power_v <- c(seq(1,10, by=1), seq(12,30, by=2)) 

## soft threshold 



##----call the network topology analysis function-111006 determine soft threshold-----
soft_111006 <- pickSoftThreshold(trans_111006, powerVector = power_v,
                                networkType = "signed", verbose = 2)
soft_111006_data <- soft_111006$fitIndices

## plot the results to visualize the data
a1<- ggplot(soft_111006_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.85, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a2<- ggplot(soft_111006_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity')

## combines both graphs in one box
grid.arrange(a1, a2, nrow = 2)
## noting my soft threshold 
softthresh_111006 = 8

 

##-----soft threshold 111010----
soft_111010 <- pickSoftThreshold(trans_111010, powerVector = power_v,
                                 networkType = "signed", verbose = 2)
soft_111010_data <- soft_111010$fitIndices

## plot the results to visualize the data
a3<- ggplot(soft_111010_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a4<- ggplot(soft_111010_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity')

## combines both graphs in one box
grid.arrange(a3, a4, nrow = 2)
softthresh_111010 = 9




##--------soft threshold 111016----

soft_111016 <- pickSoftThreshold(trans_111016, powerVector = power_v,
                                 networkType = "signed", verbose = 2)
soft_111016_data <- soft_111016$fitIndices

## plot the results to visualize the data
a5<- ggplot(soft_111016_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a6<- ggplot(soft_111016_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity')

## combines both graphs in one box
grid.arrange(a5, a6, nrow = 2)
softthresh_111016 = 12 
 


##---- soft threshold 113165------
soft_113165 <- pickSoftThreshold(trans_113165, powerVector = power_v,
                                 networkType = "signed", verbose = 2) 
soft_113165_data <- soft_113165$fitIndices

## plot the results to visualize the data
a7<- ggplot(soft_113165_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a8<- ggplot(soft_113165_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity') 

## combines both graphs in one box
grid.arrange(a7, a8, nrow = 2)
softthresh_113165 = 14


##--- soft threshold 120642------ 
soft_120642 <- pickSoftThreshold(trans_120642, powerVector = power_v,
                                 networkType = "signed", verbose = 2) 
soft_120642_data <- soft_120642$fitIndices

## plot the results to visualize the data
a9<- ggplot(soft_120642_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a10<- ggplot(soft_120642_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity') 

## combines both graphs in one box
grid.arrange(a9, a10, nrow = 2)
softthresh_120642 = 10


##----soft threshold 126865------
soft_126865 <- pickSoftThreshold(trans_126865, powerVector = power_v,
                                 networkType = "signed", verbose = 2) 
soft_126865_data <- soft_126865$fitIndices

## plot the results to visualize the data
a11<- ggplot(soft_126865_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a12<- ggplot(soft_126865_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity') 

## combines both graphs in one box
grid.arrange(a11, a12, nrow = 2)
softthresh_126865 = 16


##---- soft threshold 144304------
soft_144304 <- pickSoftThreshold(trans_144304, powerVector = power_v,
                                 networkType = "signed", verbose = 2) 
soft_144304_data <- soft_144304$fitIndices

## plot the results to visualize the data
a13<- ggplot(soft_144304_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a14<- ggplot(soft_144304_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity') 

## combines both graphs in one box
grid.arrange(a13, a14, nrow = 2) 
softthresh_144304 = 14



##---- soft threshold 151066----
soft_151066 <- pickSoftThreshold(trans_151066, powerVector = power_v,
                                 networkType = "signed", verbose = 2) 
soft_151066_data <- soft_151066$fitIndices

## plot the results to visualize the data
a15<- ggplot(soft_151066_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a16<- ggplot(soft_151066_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity') 

## combines both graphs in one box
grid.arrange(a15, a16, nrow = 2)
softthresh_151066 = 14



##--- soft threshold 159217----- 
soft_159217 <- pickSoftThreshold(trans_159217, powerVector = power_v,
                                 networkType = "signed", verbose = 2) 
soft_159217_data <- soft_159217$fitIndices

## plot the results to visualize the data
a17<- ggplot(soft_159217_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a18<- ggplot(soft_159217_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity') 

## combines both graphs in one box
grid.arrange(a17, a18, nrow = 2)
softthresh_159217 = 14



##--- soft threshold 163434 ???-----
soft_163434 <- pickSoftThreshold(trans_163434, powerVector = power_v,
                                 networkType = "signed", verbose = 2) 
soft_163434_data <- soft_163434$fitIndices

## plot the results to visualize the data
a19<- ggplot(soft_163434_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a20<- ggplot(soft_163434_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity') 

## combines both graphs in one box
grid.arrange(a19, a20, nrow = 2) 
softthresh_163434


##--- s0ft threshold 164471 ---- 
soft_164471 <- pickSoftThreshold(trans_164471, powerVector = power_v,
                                 networkType = "signed", verbose = 2) 
soft_164471_data <- soft_164471$fitIndices

## plot the results to visualize the data
a21<- ggplot(soft_164471_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a22<- ggplot(soft_164471_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity') 

## combines both graphs in one box
grid.arrange(a21, a22, nrow = 2) 
softthresh_164471 = 14



##--- soft threshold 165630-----
soft_165630 <- pickSoftThreshold(trans_165630, powerVector = power_v,
                                 networkType = "signed", verbose = 2) 
soft_165630_data <- soft_165630$fitIndices

## plot the results to visualize the data
a23<- ggplot(soft_165630_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a24<- ggplot(soft_165630_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity') 

## combines both graphs in one box
grid.arrange(a23, a24, nrow = 2) 
softthresh_165630 = 16



##--- soft threshold 174106-----
soft_174106 <- pickSoftThreshold(trans_174106, powerVector = power_v,
                                 networkType = "signed", verbose = 2) 
soft_174106_data <- soft_174106$fitIndices

## plot the results to visualize the data
a25<- ggplot(soft_174106_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a26<- ggplot(soft_174106_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity') 

## combines both graphs in one box
grid.arrange(a25, a26, nrow = 2) 
softthresh_174106 = 10




##---soft threshold 175495----
soft_175495 <- pickSoftThreshold(trans_175495, powerVector = power_v,
                                 networkType = "signed", verbose = 2) 
soft_175495_data <- soft_175495$fitIndices

## plot the results to visualize the data
a27<- ggplot(soft_175495_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a28<- ggplot(soft_175495_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity') 

## combines both graphs in one box
grid.arrange(a27, a28, nrow = 2) 
softthresh_175495 = 12



##--- soft threshold 226151------
soft_226151 <- pickSoftThreshold(trans_226151, powerVector = power_v,
                                 networkType = "signed", verbose = 2) 
soft_226151_data <- soft_226151$fitIndices

## plot the results to visualize the data
a29<- ggplot(soft_226151_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a30<- ggplot(soft_226151_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity') 

## combines both graphs in one box
grid.arrange(a29, a30, nrow = 2) 
softthresh_226151 = 14



##---- soft threshold 235781-----
soft_235781 <- pickSoftThreshold(trans_235781, powerVector = power_v,
                                 networkType = "signed", verbose = 2) 
soft_235781_data <- soft_235781$fitIndices

## plot the results to visualize the data
a31<- ggplot(soft_235781_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a32<- ggplot(soft_235781_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity') 

## combines both graphs in one box
grid.arrange(a31, a32, nrow = 2) 
softthresh_235781 = 14






##--- soft threshold 242202----
soft_242202 <- pickSoftThreshold(trans_242202, powerVector = power_v,
                                 networkType = "signed", verbose = 2) 
soft_242202_data <- soft_242202$fitIndices

## plot the results to visualize the data
a33<- ggplot(soft_242202_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a34<- ggplot(soft_242202_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity') 

## combines both graphs in one box
grid.arrange(a33, a34, nrow = 2) 
softthresh_242202 = 12









##-- co-expression similarities and adjacencies-----
adj_111006 <- adjacency(trans_111006, type = "signed", power = softthresh_111006)
adj_111010 <- adjacency(trans_111010, type = "signed", power = softthresh_111010)
adj_111016 <- adjacency(trans_111016, type = "signed", power = softthresh_111016)
adj_113165 <- adjacency(trans_113165, type = "signed", power = softthresh_113165) 
adj_120642 <- adjacency(trans_120642, type = "signed", power = softthresh_120642)
adj_126865 <- adjacency(trans_126865, type = "signed", power = softthresh_126865) 
adj_144304 <- adjacency(trans_144304, type = "signed", power = softthresh_144304)
adj_151066 <- adjacency(trans_151066, type = "signed", power = softthresh_151066)
adj_159217 <- adjacency(trans_159217, type = "signed", power = softthresh_159217) 
adj_164471 <- adjacency(trans_164471, type = "signed", power = softthresh_164471) 
adj_165630 <- adjacency(trans_165630, type = "signed", power = softthresh_165630) 
adj_174106 <- adjacency(trans_174106, type = "signed", power = softthresh_174106) 
adj_175495 <- adjacency(trans_175495, type = "signed", power = softthresh_175495) 
adj_226151 <- adjacency(trans_226151, type = "signed", power = softthresh_226151) 
adj_235781 <- adjacency(trans_235781, type = "signed", power = softthresh_235781)
adj_242202 <- adjacency(trans_242202, type = "signed", power = softthresh_242202)



##--- Topological Overlap Matrix (TOM)-----
## Turn adjacency into topological overlap 111006
TOM_111006 <- TOMsimilarity(adj_111006, TOMType = "signed") 
rownames(TOM_111006) <- tom_gene_list$Gene_Symbol
colnames(TOM_111006) <- tom_gene_list$Gene_Symbol
## TOM111010
TOM_111010 <- TOMsimilarity(adj_111010, TOMType = "signed")
rownames(TOM_111010) <- tom_gene_list$Gene_Symbol

## TOM111016
TOM_111016 <- TOMsimilarity(adj_111016, TOMType = "signed")
rownames(TOM_111016) <- tom_gene_list$Gene_Symbol

## TOM113165
TOM_113165 <- TOMsimilarity(adj_113165, TOMType = "signed")
rownames(TOM_113165) <- tom_gene_list$Gene_Symbol

## TOM120642
TOM_120642 <- TOMsimilarity(adj_120642, TOMType = "signed")
rownames(TOM_120642) <- tom_gene_list$Gene_Symbol

## TOM126865
TOM_126865 <- TOMsimilarity(adj_126865, TOMType = "signed")
rownames(TOM_126865) <- tom_gene_list$Gene_Symbol

## TOM144304
TOM_144304 <- TOMsimilarity(adj_144304, TOMType = "signed")
rownames(TOM_144304) <- tom_gene_list$Gene_Symbol

## TOM151066
TOM_151066 <- TOMsimilarity(adj_151066, TOMType = "signed")
rownames(TOM_151066) <- tom_gene_list$Gene_Symbol

## TOM159217
TOM_159217 <- TOMsimilarity(adj_159217, TOMType = "signed")
rownames(TOM_159217) <- tom_gene_list$Gene_Symbol

## TOM164471
TOM_164471 <- TOMsimilarity(adj_164471, TOMType = "signed")
rownames(TOM_164471) <- tom_gene_list$Gene_Symbol

## TOM165630
TOM_165630 <- TOMsimilarity(adj_165630, TOMType = "signed")
rownames(TOM_165630) <- tom_gene_list$Gene_Symbol

## TOM174106
TOM_174106 <- TOMsimilarity(adj_174106, TOMType = "signed")
rownames(TOM_174106) <- tom_gene_list$Gene_Symbol

## TOM175495
TOM_175495 <- TOMsimilarity(adj_175495, TOMType = "signed")
rownames(TOM_175495) <- tom_gene_list$Gene_Symbol

## TOM226151
TOM_226151 <- TOMsimilarity(adj_226151, TOMType = "signed")
rownames(TOM_226151) <- tom_gene_list$Gene_Symbol

## TOM235781
TOM_235781 <- TOMsimilarity(adj_235781, TOMType = "signed")
rownames(TOM_235781) <- tom_gene_list$Gene_Symbol

## TOM242202
TOM_242202 <- TOMsimilarity(adj_242202, TOMType = "signed")
rownames(TOM_242202) <- tom_gene_list$Gene_Symbol


##-- calculating the consensus TOM for all TOM----
# creating a tom list to be parsed into the consensusTOM function

multiExpr <- list(TOM_111006 = list(data = TOM_111006), TOM_111010 = list(data = TOM_111010),
                 TOM_111016 = list(data = TOM_111016), TOM_113165 = list(data = TOM_113165),
                 TOM_120642 = list(data = TOM_120642), TOM_126865 = list(data = TOM_126865), 
                 TOM_144304 = list(data = TOM_144304), TOM_151066 = list(data = TOM_151066), 
                 TOM_159217 = list(data = TOM_159217), TOM_164471 = list(data = TOM_164471),
                 TOM_165630 = list(data = TOM_165630), TOM_174106 = list(data = TOM_174106),
                 TOM_175495 = list(data = TOM_175495), TOM_226151 = list(data = TOM_226151),
                 TOM_235781 = list(data = TOM_235781), TOM_242202 = list(data = TOM_242202))

## consensumTOM object
cons_mat_TOM <- consensusTOM(multiExpr, checkMissingData = T, networkType = "signed",
                        TOMType = "signed", networkCalibration = "full quantile", returnTOMs = T)

# $consensusTOM - A list containing consensus TOM for each block, stored as a distance structure
# note - this function has an inbuilt distance measurement

cons_mat_TOM$ 

  
  
   
  
  
## extracted the consensusTOM list from cons_mat_TOM -----
cons_TOM <- (cons_mat_TOM$consensusTOM)

## converting cons_TOM to a matrix for downstream clustering
# isolating distance values in cons_TOM
dist_cons_TOM <- cons_TOM[[1]]

# Create dissimilarity measure and convert dissimilarity vector to a matrix
dissTOM <- 1 - dist_cons_TOM

dissTOM <- as.matrix(dissTOM)

## clustering of TOMs consensus-----

# use heirarchical cluster function 
cons_Tree <- hclust(as.dist(dissTOM), method = "average") 

# We like large modules, so we set the minimum module size relatively high:
#minModuleSize = 30;
# Module identification using dynamic tree cut:
dissTOM_unmergedlabels <- cutreeDynamic(dendro = cons_Tree, distM = dissTOM,
                               deepSplit = 2, cutHeight = 0.995,
                               minClusterSize = 30,
                               pamRespectsDendro = FALSE) 
table(dissTOM_unmergedlabels)
dissTOM_unmergedColors = labels2colors(dissTOM_unmergedlabels)
table(dissTOM_unmergedColors) 

# Plot the dendrogram and module colors
sizeGrWindow(8,6)
plotDendroAndColors(cons_Tree, dissTOM_unmergedColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, 
                    main = "Gene dendrogram and module colors") 

# merge modules whose expression profiles are similar
# - created a multiset list of expression profile dataframes for module eingengenes calculation
multiExpr_trans <- list(trans_111006 = list(data = trans_111006), trans_111010 = list(data = trans_111010),
                  trans_111016 = list(data = trans_111016), trans_113165 = list(data = trans_113165),
                  trans_120642 = list(data = trans_120642), trans_126865 = list(data = trans_126865), 
                  trans_144304 = list(data = trans_144304), trans_151066 = list(data = trans_151066), 
                  trans_159217 = list(data = trans_159217), trans_164471 = list(data = trans_164471),
                  trans_165630 = list(data = trans_165630), trans_174106 = list(data = trans_174106),
                  trans_175495 = list(data = trans_175495), trans_226151 = list(data = trans_226151),
                  trans_235781 = list(data = trans_235781), trans_242202 = list(data = trans_242202))

# Calculate module eigengenes
unmergedMEs <- multiSetMEs(multiExpr_trans, colors = NULL, universalColors = dissTOM_unmergedColors)

# Calculate consensus dissimilarity of consensus module eigengenes
consMEDiss <- consensusMEDissimilarity(unmergedMEs) 

# Cluster consensus modules
consME_Tree <- hclust(as.dist(consMEDiss), method = "average")

# Plot the result
sizeGrWindow(7,6)
par(mfrow = c(1,1))
plot(consME_Tree, main = "Consensus clustering of consensus module eigengenes",
     xlab = "", sub = "")
abline(h=0.8, col = "red") 

# merging close modules close modules
merge <-  mergeCloseModules(multiExpr_trans, dissTOM_unmergedlabels,
                            cutHeight = 0.8, verbose = 3) 

## extract module colors and eigengenes needed for downstream analysis
# Numeric module labels
module_labels <- merge$colors
table(module_labels)
# Convert labels to colors
module_colors <- labels2colors(module_labels) 
table(module_colors)
# Eigengenes of the new merged modules:
cons_MEs <- merge$newMEs 

## visualization of merge on module colors -
sizeGrWindow(9,6)
plotDendroAndColors(cons_Tree, cbind(dissTOM_unmergedColors, module_colors),
                    c("Unmerged", "Merged"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)



## enrichment ------

# extract the genes present in each module
module_genes <- split(names(vsd_111006), module_labels) 

#create a list of the modules so I can iterate through for enrichment analysis
module_col_rep <- unique(module_labels) 

# extract list of genes to serve as "universe"
adhesome_gene_list <- ahhesome_genes$entrezgene_id
adhesome_gene_list <- na.omit(adhesome_gene_list)
adhesome_gene_list <- list(adhesome_gene_list)
# an empty list to store enrichment result
GO_enrich_res <- list() 

for (modules in module_col_rep) { 
  
  # calls a list of genes in each module
  mod_gene_list <- module_genes[[modules]] 
  
  # enrichment analysis using GO terms for Molecular Function (MF) 
  res_enrich <- enrichGO(gene = mod_gene_list, universe = adhesome_gene_list,
                         OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "MF", 
                         pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
  
  # Enrichment analysis using GO terms for Biological Process (BP)
  res_enrich_BP <- enrichGO(gene = mod_gene_list, universe = adhesome_gene_list,
                            OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", 
                            pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
  
  # Enrichment analysis using GO terms for Biological Process (CC)
  res_enrich_CC <- enrichGO(gene = mod_gene_list, universe = adhesome_gene_list,
                            OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "CC", 
                            pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
 
   # Store the results
  GO_enrich_res[[paste("ME", modules, "MF", sep = "_")]] <- res_enrich
  GO_enrich_res[[paste("ME", modules, "BP", sep = "_")]] <- res_enrich_BP
  GO_enrich_res[[paste("ME", modules, "CC", sep = "_")]] <- res_enrich_CC 
} 

#res_enrich_all <- enrichGO(gene = mod_gene_list, universe = adhesome_gene_list,
      #                    OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL", 
       #                   pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
#dotplot(res_enrich_all)
## visualizing gene ontology results 
dotplot(GO_enrich_res$ME_1_MF)
dotplot(GO_enrich_res$ME_1_BP)
dotplot(GO_enrich_res$ME_1_CC)
dotplot(GO_enrich_res$ME_2_MF)
dotplot(GO_enrich_res$ME_2_BP)
dotplot(GO_enrich_res$ME_2_CC)
dotplot(GO_enrich_res$ME_4_MF)
dotplot(GO_enrich_res$ME_4_BP)
dotplot(GO_enrich_res$ME_4_CC)
dotplot(GO_enrich_res$ME_5_MF)
dotplot(GO_enrich_res$ME_5_BP)
dotplot(GO_enrich_res$ME_5_CC)
dotplot(GO_enrich_res$ME_6_MF)
dotplot(GO_enrich_res$ME_6_BP)
dotplot(GO_enrich_res$ME_6_CC)
dotplot(GO_enrich_res$ME_7_MF)
dotplot(GO_enrich_res$ME_7_BP)
dotplot(GO_enrich_res$ME_7_CC)
dotplot(GO_enrich_res$ME_8_MF)
dotplot(GO_enrich_res$ME_8_BP)
dotplot(GO_enrich_res$ME_8_CC)
dotplot(GO_enrich_res$ME_9_MF)
dotplot(GO_enrich_res$ME_9_BP)
dotplot(GO_enrich_res$ME_9_CC)
dotplot(GO_enrich_res$ME_10_MF)
dotplot(GO_enrich_res$ME_10_BP) 
dotplot(GO_enrich_res$ME_10_CC)
dotplot(GO_enrich_res$ME_11_MF)
dotplot(GO_enrich_res$ME_11_BP) 
dotplot(GO_enrich_res$ME_11_CC) 
dotplot(GO_enrich_res$ME_12_MF)
dotplot(GO_enrich_res$ME_12_BP)
dotplot(GO_enrich_res$ME_12_CC)

## summary of GO result 
# Initialize a data frame to store the summary
go_summary <- data.frame()

# Loop through each module and summarize the results
for (module in names(GO_enrich_res)) {
  if (!is.null(GO_enrich_res[[module]])) {
    go_result <- as.data.frame(GO_enrich_res[[module]]@result)
    if (nrow(go_result) > 0) {
      go_result$Module <- module
      go_summary <- rbind(go_summary, go_result)
    }
  }
} 

# Select and rename columns for better readability
go_summary <- go_summary[, c("Module", "Description", "ID", "pvalue", "p.adjust", "qvalue", "geneID", "Count")]
colnames(go_summary) <- c("Module", "GO_Term", "Ontology", "PValue", "AdjustedPValue", "QValue", "GeneID", "GeneCount")

# Print the summary
print(go_summary)




## correlate modules to age group and disease condition - DE with limma -----
## identify significantly expressed modules
# extract module eigengene for each dataset. 
              
MEs_111006 <- t(merge[["newMEs"]][["trans_111006"]][["data"]])
MEs_111010 <- t(merge[["newMEs"]][["trans_111010"]][["data"]])
MEs_111016 <- t(merge[["newMEs"]][["trans_111016"]][["data"]])
MEs_113165 <- t(merge[["newMEs"]][["trans_113165"]][["data"]])
MEs_120642 <- t(merge[["newMEs"]][["trans_120642"]][["data"]])
MEs_126865 <- t(merge[["newMEs"]][["trans_126865"]][["data"]])
MEs_144304 <- t(merge[["newMEs"]][["trans_144304"]][["data"]])
MEs_151066 <- t(merge[["newMEs"]][["trans_151066"]][["data"]])
MEs_159217 <- t(merge[["newMEs"]][["trans_159217"]][["data"]])
MEs_164471 <- t(merge[["newMEs"]][["trans_164471"]][["data"]])
MEs_165630 <- t(merge[["newMEs"]][["trans_165630"]][["data"]])
MEs_174106 <- t(merge[["newMEs"]][["trans_174106"]][["data"]])
MEs_175495 <- t(merge[["newMEs"]][["trans_175495"]][["data"]])
MEs_226151 <- t(merge[["newMEs"]][["trans_226151"]][["data"]])
MEs_235781 <- t(merge[["newMEs"]][["trans_235781"]][["data"]])
MEs_242202 <- t(merge[["newMEs"]][["trans_242202"]][["data"]])

# convert variables to numeric values for DE analysis
# these numeric values will serve as design matrix  
# 111006
trait_111006 <- factor(cond_111006$sarcopenia, levels = c("yes", "no"))  
design_111006 <- model.matrix(~ 0 + trait_111006) 
colnames(design_111006) <- c("sarcopenia.yes", "sarcopenia.no") 
rownames(design_111006) <- rownames(cond_111006)
# DE limma
fit_111006 <- lmFit(MEs_111006, design_111006) 
contrast_111006 <- makeContrasts(sarcopenia.yesvssarcopenia.no = sarcopenia.yes - sarcopenia.no,
                                 levels = design_111006) 
fitcon_111006 <- contrasts.fit(fit_111006, contrast_111006) 
fitcon_111006 <- eBayes(fitcon_111006) 

# 111010
trait_111010 <- factor(cond_111010$sarcopenia, levels = c("yes", "no"))
design_111010 <- model.matrix(~ 0 + trait_111010)  
colnames(design_111010) <- c("sarcopenia.yes", "sarcopenia.no")   
# DE limma
fit_111010 <- lmFit(MEs_111010, design_111010) 
contrast_111010 <- makeContrasts(sarcopenia.yesvssarcopenia.no = sarcopenia.yes - sarcopenia.no,
                                 levels = design_111010) 
fitcon_111010 <- contrasts.fit(fit_111010, contrast_111010) 
fitcon_111010 <- eBayes(fitcon_111010) 

# 111016
trait_111016 <- factor(cond_111016$sarcopenia, levels = c("yes", "no"))
design_111016 <- model.matrix(~ 0 + trait_111016)
colnames(design_111016) <- c("sarcopenia.yes", "sarcopenia.no")
# DE limma
fit_111016 <- lmFit(MEs_111016, design_111016)  
contrast_111016 <- makeContrasts(sarcopenia.yesvssarcopenia.no = sarcopenia.yes - sarcopenia.no,
                                 levels = design_111016)  
fitcon_111016 <- contrasts.fit(fit_111016, contrast_111016) 
fitcon_111016 <- eBayes(fitcon_111016)


# 113165
# instruction from limma vignette
trait_113165 <- factor(cond_113165$age, levels = c("young", "old"))   
design_113165 <- model.matrix(~ 0 + trait_113165) 
colnames(design_113165) <- levels(trait_113165) 
# DE limma
fit_113165 <- lmFit(MEs_113165, design_113165)  
contrast_113165 <- makeContrasts(oldvsyoung = old - young, levels = design_113165) 
fitcon_113165 <- contrasts.fit(fit_113165, contrast_113165) 
fitcon_113165 <- eBayes(fitcon_113165)

# 120642
trait_120642 <- factor(cond_120642$diagnosis, 
                       levels = c("Healthy_Adult",
                                  "Critical_limb_ischemia",
                                  "Intermittent_claudicant"))
design_120642 <- model.matrix(~ 0 + trait_120642) 
colnames(design_120642) <- levels(trait_120642) 
# DE limma 
fit_120642 <- lmFit(MEs_120642, design_120642) 
contrast_120642 <- makeContrasts(Healthy_Adult - Critical_limb_ischemia,
                                 Healthy_Adult - Intermittent_claudicant,
                                 Intermittent_claudicant - Critical_limb_ischemia,
                                 levels = design_120642) 
fitcon_120642 <- contrasts.fit(fit_120642, contrast_120642)
fitcon_120642 <- eBayes(fitcon_120642) 

# 126865
# set factor and level
trait_126865 <- factor(cond_126865$group, levels = c("after", "before")) 
design_126865 <- model.matrix(~ 0 + trait_126865)
colnames(design_126865) <- levels(trait_126865) 
# DE limma
fit_126865 <- lmFit(MEs_126865, design_126865)  
contrast_126865 <- makeContrasts(beforevsafter = before-after, 
                                 levels = design_126865) 
fitcon_126865 <- contrasts.fit(fit_126865, contrast_126865)
fitcon_126865 <- eBayes(fitcon_126865)

# 144304
trait_144304 <- factor(cond_144304$muscle_status, levels = c("Young", "Frail", "Fit"))
design_144304 <- model.matrix(~ 0 + trait_144304) 
colnames(design_144304) <- c("Young", "Frail", "Fit")
# DE limma
fit_144304 <- lmFit(MEs_144304, design_144304) 
contrast_144304 <- makeContrasts(Frail - Young, Fit - Frail, Fit - Young,
                                 levels = design_144304)  
fitcon_144304 <- contrasts.fit(fit_144304, contrast_144304) 
fitcon_144304 <- eBayes(fitcon_144304) 

# 151066
trait_151066 <- factor(cond_151066$cohort, levels = c("Active", "Sedentary")) 
design_151066 <- model.matrix(~ 0 + trait_151066)
colnames(design_151066) <- levels(trait_151066) 
#DE limma
fit_151066 <- lmFit(MEs_151066, design_151066) 
contrast_151066 <- makeContrasts(SedentaryvsActive = Sedentary - Active,
                                 levels = design_151066)  
fitcon_151066 <- contrasts.fit(fit_151066, contrast_151066) 
fitcon_151066 <- eBayes(fitcon_151066)

# 159217
trait_159217 <- factor(cond_159217$age, levels = c("Older", "Young"))
design_159217 <- model.matrix(~ 0 + trait_159217) 
colnames(design_159217) <- c("Older", "Young") 
# DE limma
fit_159217 <- lmFit(MEs_159217, design_159217) 
contrast_159217 <- makeContrasts(OldervsYoung = Older-Young,
                                 levels = design_159217)
fitcon_159217 <- contrasts.fit(fit_159217, contrast_159217) 
fitcon_159217 <- eBayes(fitcon_159217) 

# 164471
cond_164471 <- cond_164471 %>%
  mutate(age = as.numeric(as.character(age)))%>%  # change to numeric - was factor before
  mutate(age_group = ifelse(age >= 20 & age <= 40, "young", "old")) #set age range newcolumn
# set factor and levels
trait_164471 <- factor(cond_164471$age_group, levels = c("old", "young"))
design_164471 <- model.matrix(~ 0 + trait_164471) 
colnames(design_164471) <- levels(trait_164471)
# DE limma
fit_164471 <- lmFit(MEs_164471, design_164471) 
contrast_164471 <- makeContrasts(oldvsyoung = old-young, 
                                 levels = design_164471) 
fitcon_164471 <- contrasts.fit(fit_164471, contrast_164471) 
fitcon_164471 <- eBayes(fitcon_164471)

# 165630
trait_165630 <- factor(cond_165630$title,
                       levels = c("endurance_trained", "resistance_trained", "sedentary"))
design_165630 <- model.matrix(~ 0 + trait_165630) 
colnames(design_165630) <- levels(trait_165630) 
# DE limma
fit_165630 <- lmFit(MEs_165630, design_165630) 
contrast_165630 <- makeContrasts(resistance_trained-endurance_trained, 
                                 sedentary-resistance_trained,
                                 sedentary-endurance_trained,
                                 levels = design_165630)  
fitcon_165630 <- contrasts.fit(fit_165630, contrast_165630) 
fitcon_165630 <- eBayes(fitcon_165630) 

# 174106
trait_174106 <- factor(cond_174106$lean_mass, levels = c("stable_low",
                                                         "sarcopenic",              
                                                         "gain", "stable_high"))   
design_174106 <- model.matrix(~ 0 + trait_174106) 
colnames(design_174106) <- levels(trait_174106) 
# DE limma
fit_174106 <- lmFit(MEs_174106, design_174106) 
contrast_174106 <- makeContrasts(stable_high-stable_low, gain-stable_low,
                                 gain-sarcopenic,
                                 stable_high-sarcopenic, gain-stable_high,
                                 stable_low-sarcopenic, 
                                 levels = design_174106)  
fitcon_174106 <- contrasts.fit(fit_174106, contrast_174106) 
fitcon_174106 <- eBayes(fitcon_174106) 


# 175495 
trait_175495 <- factor(cond_175495$age, levels = c("Young", "Old")) 
design_175495 <- model.matrix(~ 0 + trait_175495) 
colnames(design_175495) <- levels(trait_175495) 
# DE limma
fit_175495 <- lmFit(MEs_175495, design_175495) 
contrast_175495 <- makeContrasts(OldvsYoung = Old-Young,
                                 levels = design_175495)
fitcon_175495 <- contrasts.fit(fit_175495, contrast_175495) 
fitcon_175495 <- eBayes(fitcon_175495)

# 226151
trait_226151 <- factor(cond_226151$status, 
                       levels = c("healthy_aged", "pre_sarcopenia", "sarcopenia"))
design_226151 <- model.matrix(~ 0 + trait_226151) 
colnames(design_226151) <- levels(trait_226151) 
# DE limma
fit_226151 <- lmFit(MEs_226151, design_226151)  
contrast_226151 <- makeContrasts(sarcopenia-healthy_aged,
                                 sarcopenia-pre_sarcopenia, 
                                 pre_sarcopenia-healthy_aged,
                                 levels = design_226151) 
fitcon_226151 <- contrasts.fit(fit_226151, contrast_226151) 
fitcon_226151 <- eBayes(fitcon_226151)

# 235781
trait_235781 <- factor(cond_235781$age, levels = c("young", "old"))   
design_235781 <- model.matrix(~ 0 + trait_235781)  
colnames(design_235781) <- levels(trait_235781) 
# DE limma
fit_235781 <- lmFit(MEs_235781, design_235781)  
contrast_235781 <- makeContrasts(oldvsyoung = old-young,
                                 levels = design_235781) 
fitcon_235781 <- contrasts.fit(fit_235781, contrast_235781) 
fitcon_235781 <- eBayes(fitcon_235781) 

# 242202
trait_242202 <- factor(cond_242202$disease, levels = c("patient", "healthy")) 
design_242202 <- model.matrix(~ 0 + trait_242202) 
colnames(design_242202) <- levels(trait_242202)  
# DE limma
fit_242202 <- lmFit(MEs_242202, design_242202)  
contrast_242202 <- makeContrasts(healthyvspatient = healthy-patient,
                                 levels = design_242202) 
fitcon_242202 <- contrasts.fit(fit_242202, contrast_242202) 
fitcon_242202 <- eBayes(fitcon_242202)



## differentially expressed modules visualize -----

sigex_111006 <- topTable(fitcon_111006, number = 11)
## Visualize result of differentially expressed modules- bar plot
barplot(sigex_111006$logFC, names.arg = rownames(sigex_111006), las = 2, 
        col = ifelse(sigex_111006$logFC > 0, "navy", "brown"),
        ylab = "log Fold Change", main = "Differential Expression in sarcopenia yes vs no")
# Add a horizontal line at y=0
abline(h = 0, col = "black")

sigex_111010 <- topTable(fitcon_111010, number = 11)
## Visualize result of differentially expressed modules- bar plot
barplot(sigex_111010$logFC, names.arg = rownames(sigex_111010), las = 2, 
        col = ifelse(sigex_111010$logFC > 0, "navy", "brown"),
        ylab = "log Fold Change", main = "Differential Expression in sarcopenia yes vs no")
# Add a horizontal line at y=0
abline(h = 0, col = "black")

sigex_111016 <- topTable(fitcon_111016, number = 11)
## Visualize result of differentially expressed modules- bar plot
barplot(sigex_111016$logFC, names.arg = rownames(sigex_111016), las = 2, 
        col = ifelse(sigex_111016$logFC > 0, "navy", "brown"),
        ylab = "log Fold Change", main = "Differential Expression in sarcopenia yes vs no")
# Add a horizontal line at y=0
abline(h = 0, col = "black") 


sigex_113165 <- topTable(fitcon_113165, number = 11)
## Visualize result of differentially expressed modules- bar plot
barplot(sigex_113165$logFC, names.arg = rownames(sigex_113165), las = 2, 
        col = ifelse(sigex_113165$logFC > 0, "navy", "brown"),
        ylab = "log Fold Change", main = "Differential Expression in old vs young")
# Add a horizontal line at y=0
abline(h = 0, col = "black") 


sigex_120642 <- topTable(fitcon_120642, number = 11)
## Visualize result of differentially expressed modules- bar plot
barplot(sigex_120642$Healthy_Adult...Critical_limb_ischemia, names.arg = rownames(sigex_120642), las = 2, 
        col = ifelse(sigex_120642$Healthy_Adult...Critical_limb_ischemia > 0, "navy", "brown"),
        ylab = "log Fold Change", main = "Differential Expression in healthy adult vs critical limb ischemia")
# Add a horizontal line at y=0
abline(h = 0, col = "black")


sigex_126865 <- topTable(fitcon_126865, number = 11)
## Visualize result of differentially expressed modules- bar plot
barplot(sigex_126865$logFC, names.arg = rownames(sigex_126865), las = 2, 
        col = ifelse(sigex_126865$logFC > 0, "navy", "brown"),
        ylab = "log Fold Change", main = "Differential Expression in adults before vs after bedrest")
# Add a horizontal line at y=0
abline(h = 0, col = "black")

sigex_144304 <- topTable(fitcon_144304, number = 11)
## Visualize result of differentially expressed modules- bar plot
barplot(sigex_144304$Frail...Young, names.arg = rownames(sigex_144304), las = 2, 
        col = ifelse(sigex_144304$Frail...Young > 0, "navy", "brown"),
        ylab = "log Fold Change", main = "Differential Expression in frail adults vs young")
# Add a horizontal line at y=0
abline(h = 0, col = "black") 
barplot(sigex_144304$Fit...Frail, names.arg = rownames(sigex_144304), las = 2, 
        col = ifelse(sigex_144304$Fit...Frail > 0, "navy", "brown"),
        ylab = "log Fold Change", main = "Differential Expression in fit adults vs frail adults")
# Add a horizontal line at y=0
abline(h = 0, col = "black") 
barplot(sigex_144304$Fit...Young, names.arg = rownames(sigex_144304), las = 2, 
        col = ifelse(sigex_144304$Fit...Young > 0, "navy", "brown"),
        ylab = "log Fold Change", main = "Differential Expression in fit adults vs young")
# Add a horizontal line at y=0
abline(h = 0, col = "black")

sigex_151066 <- topTable(fitcon_151066, number = 11) 
## Visualize result of differentially expressed modules- bar plot
barplot(sigex_151066$logFC, names.arg = rownames(sigex_151066), las = 2, 
        col = ifelse(sigex_151066$logFC > 0, "navy", "brown"),
        ylab = "log Fold Change", main = "Differential Expression in sedentry vs active adults")
# Add a horizontal line at y=0
abline(h = 0, col = "black") 

sigex_159217 <- topTable(fitcon_159217, number = 11) 
## Visualize result of differentially expressed modules- bar plot
barplot(sigex_159217$logFC, names.arg = rownames(sigex_159217), las = 2, 
        col = ifelse(sigex_159217$logFC > 0, "navy", "brown"),
        ylab = "log Fold Change", main = "Differential Expression in older vs young")
# Add a horizontal line at y=0
abline(h = 0, col = "black") 

sigex_164471 <- topTable(fitcon_164471, number = 11) 
## Visualize result of differentially expressed modules- bar plot
barplot(sigex_164471$logFC, names.arg = rownames(sigex_164471), las = 2, 
        col = ifelse(sigex_164471$logFC > 0, "navy", "brown"),
        ylab = "log Fold Change", main = "Differential Expression in old vs young")
# Add a horizontal line at y=0
abline(h = 0, col = "black") 


sigex_165630 <- topTable(fitcon_165630, number = 11) 
## Visualize result of differentially expressed modules- bar plot
barplot(sigex_165630$resistance_trained...endurance_trained, names.arg = rownames(sigex_165630), las = 2, 
        col = ifelse(sigex_165630$resistance_trained...endurance_trained > 0, "navy", "brown"),
        ylab = "log Fold Change", main = "Differential Expression in RET vs ET")
# Add a horizontal line at y=0
abline(h = 0, col = "black") 

barplot(sigex_165630$sedentary...resistance_trained, names.arg = rownames(sigex_165630), las = 2, 
        col = ifelse(sigex_165630$sedentary...resistance_trained > 0, "navy", "brown"),
        ylab = "log Fold Change", main = "Differential Expression in sedentry adult vs RET")
# Add a horizontal line at y=0
abline(h = 0, col = "black")

barplot(sigex_165630$sedentary...endurance_trained, names.arg = rownames(sigex_165630), las = 2, 
        col = ifelse(sigex_165630$sedentary...endurance_trained > 0, "navy", "brown"),
        ylab = "log Fold Change", main = "Differential Expression in sedentry adult vs ET")
# Add a horizontal line at y=0
abline(h = 0, col = "black")

sigex_174106 <- topTable(fitcon_174106, number = 11) 
## Visualize result of differentially expressed modules- bar plot
barplot(sigex_174106$gain...sarcopenic, names.arg = rownames(sigex_174106), las = 2, 
        col = ifelse(sigex_174106$gain...sarcopenic > 0, "navy", "brown"),
        ylab = "log Fold Change", main = "Differential Expression in gain vs sarcopenic")
# Add a horizontal line at y=0
abline(h = 0, col = "black") 

barplot(sigex_174106$stable_high...sarcopenic, names.arg = rownames(sigex_174106), las = 2, 
        col = ifelse(sigex_174106$stable_high...sarcopenic > 0, "navy", "brown"),
        ylab = "log Fold Change", main = "Differential Expression in stable high vs sarcopenic")
# Add a horizontal line at y=0
abline(h = 0, col = "black")  

barplot(sigex_174106$stable_low...sarcopenic, names.arg = rownames(sigex_174106), las = 2, 
        col = ifelse(sigex_174106$stable_low...sarcopenic > 0, "navy", "brown"),
        ylab = "log Fold Change", main = "Differential Expression in stable low vs sarcopenic")
# Add a horizontal line at y=0
abline(h = 0, col = "black")

sigex_175495 <- topTable(fitcon_175495, number = 11)
## Visualize result of differentially expressed modules- bar plot
barplot(sigex_175495$logFC, names.arg = rownames(sigex_175495), las = 2, 
        col = ifelse(sigex_175495$logFC > 0, "navy", "brown"),
        ylab = "log Fold Change", main = "Differential Expression in old vs young")
# Add a horizontal line at y=0
abline(h = 0, col = "black") 

sigex_226151 <- topTable(fitcon_226151, number = 11)
## Visualize result of differentially expressed modules- bar plot
barplot(sigex_226151$sarcopenia...healthy_aged, names.arg = rownames(sigex_226151), las = 2, 
        col = ifelse(sigex_226151$sarcopenia...healthy_aged > 0, "navy", "brown"),
        ylab = "log Fold Change", main = "Differential Expression in sarcopenia vs healthy aged")
# Add a horizontal line at y=0
abline(h = 0, col = "black") 

sigex_235781 <- topTable(fitcon_235781, number = 11) 
## Visualize result of differentially expressed modules- bar plot
barplot(sigex_235781$logFC, names.arg = rownames(sigex_235781), las = 2, 
        col = ifelse(sigex_235781$logFC > 0, "navy", "brown"),
        ylab = "log Fold Change", main = "Differential Expression in old vs young")
# Add a horizontal line at y=0
abline(h = 0, col = "black")

sigex_242202 <- topTable(fitcon_242202, number = 11)
## Visualize result of differentially expressed modules for each dataset
# Create a bar plot
barplot(sigex_242202$logFC, names.arg = rownames(sigex_242202), las = 2, col = ifelse(sigex_242202$logFC > 0, "navy", "brown"),
        ylab = "log Fold Change", main = "Differential Expression of Module Eigengenes in healthy young vs old")
# Add a horizontal line at y=0
abline(h = 0, col = "black") 

summary(decideTests(fitcon_111016)) 
summary(decideTests(fitcon_174106)) 
summary(decideTests(fitcon_242202)) 



# module hub genes -----
hub_cons <- consensusKME(multiExpr, module_labels) 
# extract modules meta Z equal weights 
hub_gene_kME<- hub_cons %>%
  dplyr::select(c("ID", "meta.Z.equalWeights.kME1", "meta.Z.equalWeights.kME2", "meta.Z.equalWeights.kME4",
           "meta.Z.equalWeights.kME5", "meta.Z.equalWeights.kME6", "meta.Z.equalWeights.kME7", 
           "meta.Z.equalWeights.kME8", "meta.Z.equalWeights.kME9", "meta.Z.equalWeights.kME10",
           "meta.Z.equalWeights.kME11", "meta.Z.equalWeights.kME12"))



#create a df with gene id and gene symbol------
tom_gene_list <- vsd_113165_df %>%
  dplyr::select(1,5)
tom_gene_list$entrez_gene_id <- as.integer(rownames(adj_113165))  
tom_gene_list <- tom_gene_list %>%
  left_join(., ahhesome_genes, by = c('entrez_gene_id' = 'entrezgene_id'))
tom_gene_list <- tom_gene_list[, -c(1,2,5,6,7,8)] 
# filter by rownames of adj_111006 to match working gene number of 2382
tom_gene_list <- tom_gene_list[tom_gene_list$entrez_gene_id %in% rownames(adj_111006),]
#created unique to match rownames
unique_entrez_gene_id <- make.unique(as.character(tom_gene_list$entrez_gene_id))
rownames(tom_gene_list) <- unique_entrez_gene_id
# used intersect to find common genes and filtered tom gene list by intersect
inertom_gene_list <- intersect(rownames(tom_gene_list), rownames(adj_111006))
tom_gene_list <- tom_gene_list[inertom_gene_list, ]
# confirm the equal rownames
all(rownames(tom_gene_list) %in% rownames(adj_111006))
# make gene symbols the


