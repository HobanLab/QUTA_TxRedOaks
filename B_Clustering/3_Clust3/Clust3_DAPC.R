# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% QUTA & TEXAS RED OAKS : RUNNING DAPC %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This is an exploratory script, for running a DAPC on the QUTA and Texas Red Oak samples

library(adegenet)

# Folder containing different outputs
QUTA_TRO_wd <- '/RAID1/QUTA_TX_RedOaks/Genotyping/SNP_Calling/Output/'
setwd(QUTA_TRO_wd)
# Colors for different clusters
QUTA_colors <- brewer.pal(12, "Paired")

# %%%% CLUST3 %%%% ----
# SE-READS ONLY, THIRD CLUSTERING RUN
seOnly_Clust3_folder <- paste0(QUTA_TRO_wd,'Clust3_seOnly/pop_R98/')
# Read in genind file
QUTA_Clust3_gen <- 
  read.genepop(paste0(seOnly_Clust3_folder, "populations.snps.gen"))
# DAPC ----
# %%%% K=3 ----
# # K-means clustering step--K=3
# QUTA_Clust3_gen.grp <- find.clusters(QUTA_Clust3_gen)
# # Retained PCs (here, there is no cost to retaining a lot of PCs, even greater than xlim): 150
# # Clusters (seeking to minimize the Bayesian Information Criterion value): 3
# 
# # PCA step
# QUTA_Clust3_gen.dapc <- dapc(x=QUTA_Clust3_gen, pop=QUTA_Clust3_gen.grp$grp)
# # Retained PCs (specifying a large value can lead to overfitting/unstable membership probability): 80
# # Discriminant functions to retain (for less than 10 clusters, all eigenvalues can be retained): 3
# # List grouping of each individual (helpful for determining correct coloration)
# cbind(QUTA_Clust3_gen.dapc$grp, rownames(QUTA_Clust3_gen@tab))
# # Save the DAPC object to disk, so as to recreate in future analyses
# saveRDS(QUTA_Clust3_gen.dapc,
#         "/home/akoontz/Documents/QUTA_TxRedOaks/Code/Clust3_DAPC_K3.Rdata")

# Read in the previously saved DAPC object
QUTA_Clust3_gen.dapc <- 
  readRDS("/home/akoontz/Documents/QUTA_TxRedOaks/Code/Clust3_DAPC_K3.Rdata")
# List grouping of each individual (helpful for determining correct coloration)
cbind(QUTA_Clust3_gen.dapc$grp, rownames(QUTA_Clust3_gen@tab))
# Show DAPC as a scatterplot
scatter(QUTA_Clust3_gen.dapc, scree.da=F, bg="white", pch=20, cell=0, cstar=0, solid=0.6, clab=0, legend=T,
        posi.leg=locator(n=1), cleg=1.0, cex=2, inset.solid=1,
        col = c('#B2DF8A','#1F78B4','#A6CEE3'),
        txt.leg = c('"emoryi"','"hypoleucoides"','"gravesii"'))
mtext("Quercus tardifolia/Texas Red Oaks: K=3", adj=0.09, line=1.3)

# %%%% K=4 ----
# K-means clustering step--K=3
# QUTA_Clust3_gen.grp <- find.clusters(QUTA_Clust3_gen)
# # Retained PCs (here, there is no cost to retaining a lot of PCs, even greater than xlim): 120
# # Clusters (seeking to minimize the Bayesian Information Criterion value): 4
# 
# # PCA step
# QUTA_Clust3_gen.dapc <- dapc(x=QUTA_Clust3_gen, pop=QUTA_Clust3_gen.grp$grp)
# # Retained PCs (specifying a large value can lead to overfitting/unstable membership probability): 80
# # Discriminant functions to retain (for less than 10 clusters, all eigenvalues can be retained): 4
# # List grouping of each individual (helpful for determining correct coloration)
# cbind(QUTA_Clust3_gen.dapc$grp, rownames(QUTA_Clust3_gen@tab))
# # Save the DAPC object to disk, so as to recreate in future analyses
# saveRDS(QUTA_Clust3_gen.dapc,
#         "/home/akoontz/Documents/QUTA_TxRedOaks/Code/Clust3_DAPC_K4.Rdata")

# Read in the previously saved DAPC object
QUTA_Clust3_gen.dapc <- 
  readRDS("/home/akoontz/Documents/QUTA_TxRedOaks/Code/Clust3_DAPC_K4.Rdata")
# List grouping of each individual (helpful for determining correct coloration)
cbind(QUTA_Clust3_gen.dapc$grp, rownames(QUTA_Clust3_gen@tab))
# Show DAPC as a scatterplot
scatter(QUTA_Clust3_gen.dapc, scree.da=F, bg="white", pch=20, cell=0, cstar=0, solid=0.6, clab=0, legend=T,
        posi.leg=locator(n=1), cleg=1.0, cex=2, inset.solid=1,
        col = QUTA_colors[1:4],
        txt.leg = c('"gravesii"','"graciliformis/emoryi','"hypoleucoides"','"canbyi"'))
mtext("Quercus tardifolia/Texas Red Oaks: K=4", adj=0.09, line=1.3)

# %%%% SCY-HYPO-CAN-TAR %%%% ----
# SE-READS ONLY, SCY-HYPO-CAN-TAR RUN 
seOnly_ScyHypoCanTar_folder <- paste0(QUTA_TRO_wd,'ScyHypoCanTar_seOnly/pop_R97_WL10K/')
# Read in genind file
QUTA_ScyHypoCanTar_gen <- 
  read.genepop(paste0(seOnly_ScyHypoCanTar_folder, "populations.snps.gen"))
# # DAPC ----
# # %%%% K=2 ----
# # K-means clustering step--K=2
# QUTA_ScyHypoCanTar_gen.grp <- find.clusters(QUTA_ScyHypoCanTar_gen)
# # Retained PCs (here, there is no cost to retaining a lot of PCs, even greater than xlim): 100
# # Clusters (seeking to minimize the Bayesian Information Criterion value): 2
# 
# # PCA step
# QUTA_ScyHypoCanTar_gen.dapc <- dapc(x=QUTA_ScyHypoCanTar_gen, pop=QUTA_ScyHypoCanTar_gen.grp$grp)
# # Retained PCs (specifying a large value can lead to overfitting/unstable membership probability): 40
# # Discriminant functions to retain (for less than 10 clusters, all eigenvalues can be retained): 2
# # List grouping of each individual (helpful for determining correct coloration)
# cbind(QUTA_ScyHypoCanTar_gen.dapc$grp, rownames(QUTA_ScyHypoCanTar_gen@tab))
# # Save the DAPC object to disk, so as to recreate in future analyses
# saveRDS(QUTA_ScyHypoCanTar_gen.dapc,
#         "/home/akoontz/Documents/QUTA_TxRedOaks/Code/B_Clustering/DAPC_outputs/ScyHypoCanTar_DAPC_K2.Rdata")
# Read in the previously saved DAPC object
QUTA_ScyHypoCanTar_gen.dapc <- 
  readRDS("/home/akoontz/Documents/QUTA_TxRedOaks/Code/B_Clustering/DAPC_outputs/ScyHypoCanTar_DAPC_K2.Rdata")
# List grouping of each individual (helpful for determining correct coloration)
cbind(QUTA_ScyHypoCanTar_gen.dapc$grp, rownames(QUTA_ScyHypoCanTar_gen@tab))
# Show DAPC as a scatterplot
scatter(QUTA_ScyHypoCanTar_gen.dapc, scree.da=F, bg="white", pch=20, cell=0, cstar=0, solid=0.6, clab=0, legend=T,
        posi.leg=locator(n=1), cleg=1.0, cex=2, inset.solid=1,
        col = c('#B2DF8A','#1F78B4'),
        txt.leg = c('canbyi, miq., grav x hypo., tard','hypo., scy.'))
mtext("Quercus tardifolia, ScyHypoCanTar Dataset: K=2", adj=0.09, line=1.3)
