# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% QUTA & TEXAS RED OAKS : RUNNING DAPC %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This is an exploratory script, for running a DAPC on the QUTA and Texas Red Oak samples

library(adegenet)

# Folder containing different outputs
QUTA_TRO_wd <- '/RAID1/QUTA_TX_RedOaks/Genotyping/SNP_Calling/Output/'
setwd(QUTA_TRO_wd)
# SE-READS ONLY, THIRD CLUSTERING RUN ----
seOnly_Clust3_folder <- paste0(QUTA_TRO_wd,'Clust3_seOnly/pop_R98/')
# Read in genind file
QUTA_TRO_gen <- 
  read.genepop(paste0(seOnly_Clust3_folder, "populations.snps.gen"))
# Colors for different clusters
QUTA_colors <- brewer.pal(12, "Paired")

# DAPC ----
# %%%% K=3 ----
# # K-means clustering step--K=3
# QUTA_TRO_gen.grp <- find.clusters(QUTA_TRO_gen)
# # Retained PCs (here, there is no cost to retaining a lot of PCs, even greater than xlim): 150
# # Clusters (seeking to minimize the Bayesian Information Criterion value): 3
# 
# # PCA step
# QUTA_TRO_gen.dapc <- dapc(x=QUTA_TRO_gen, pop=QUTA_TRO_gen.grp$grp)
# # Retained PCs (specifying a large value can lead to overfitting/unstable membership probability): 80
# # Discriminant functions to retain (for less than 10 clusters, all eigenvalues can be retained): 3
# # List grouping of each individual (helpful for determining correct coloration)
# cbind(QUTA_TRO_gen.dapc$grp, rownames(QUTA_TRO_gen@tab))
# # Save the DAPC object to disk, so as to recreate in future analyses
# saveRDS(QUTA_TRO_gen.dapc,
#         "/home/akoontz/Documents/QUTA_TxRedOaks/Code/Clust3_DAPC_K3.Rdata")

# Read in the previously saved DAPC object
QUTA_TRO_gen.dapc <- 
  readRDS("/home/akoontz/Documents/QUTA_TxRedOaks/Code/Clust3_DAPC_K3.Rdata")
# List grouping of each individual (helpful for determining correct coloration)
cbind(QUTA_TRO_gen.dapc$grp, rownames(QUTA_TRO_gen@tab))
# Show DAPC as a scatterplot
scatter(QUTA_TRO_gen.dapc, scree.da=F, bg="white", pch=20, cell=0, cstar=0, solid=0.6, clab=0, legend=T,
        posi.leg=locator(n=1), cleg=1.0, cex=2, inset.solid=1,
        col = c('#B2DF8A','#1F78B4','#A6CEE3'),
        txt.leg = c('"emoryi"','"hypoleucoides"','"gravesii"'))
mtext("Quercus tardifolia/Texas Red Oaks: K=3", adj=0.09, line=1.3)

# %%%% K=4 ----
# K-means clustering step--K=3
# QUTA_TRO_gen.grp <- find.clusters(QUTA_TRO_gen)
# # Retained PCs (here, there is no cost to retaining a lot of PCs, even greater than xlim): 120
# # Clusters (seeking to minimize the Bayesian Information Criterion value): 4
# 
# # PCA step
# QUTA_TRO_gen.dapc <- dapc(x=QUTA_TRO_gen, pop=QUTA_TRO_gen.grp$grp)
# # Retained PCs (specifying a large value can lead to overfitting/unstable membership probability): 80
# # Discriminant functions to retain (for less than 10 clusters, all eigenvalues can be retained): 4
# # List grouping of each individual (helpful for determining correct coloration)
# cbind(QUTA_TRO_gen.dapc$grp, rownames(QUTA_TRO_gen@tab))
# # Save the DAPC object to disk, so as to recreate in future analyses
# saveRDS(QUTA_TRO_gen.dapc,
#         "/home/akoontz/Documents/QUTA_TxRedOaks/Code/Clust3_DAPC_K4.Rdata")

# Read in the previously saved DAPC object
QUTA_TRO_gen.dapc <- 
  readRDS("/home/akoontz/Documents/QUTA_TxRedOaks/Code/Clust3_DAPC_K4.Rdata")
# List grouping of each individual (helpful for determining correct coloration)
cbind(QUTA_TRO_gen.dapc$grp, rownames(QUTA_TRO_gen@tab))
# Show DAPC as a scatterplot
scatter(QUTA_TRO_gen.dapc, scree.da=F, bg="white", pch=20, cell=0, cstar=0, solid=0.6, clab=0, legend=T,
        posi.leg=locator(n=1), cleg=1.0, cex=2, inset.solid=1,
        col = QUTA_colors[1:4],
        txt.leg = c('"gravesii"','"graciliformis/emoryi','"hypoleucoides"','"canbyi"'))
mtext("Quercus tardifolia/Texas Red Oaks: K=4", adj=0.09, line=1.3)
