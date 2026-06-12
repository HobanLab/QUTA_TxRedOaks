# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% QUTA & TEXAS RED OAKS : RUNNING PCA %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This is an exploratory script, for running a PCA on the QUTA and Texas Red Oak samples. It uses
# the Clust3 dataset, and removes a few individuals from this dataset beforehand--robusta, sctophylla,
# others. This is to clarify the relations between tardifolia, gravesii, and hypoleucoides.
pacman::p_load(adegenet, dartR, vcfR, poppr, ggplot2, ggrepel, cowplot, viridis)

# Folder containing different outputs
QUTA_TRO_wd <- '/RAID1/QUTA_TX_RedOaks/Genotyping/SNP_Calling/Output/'
setwd(QUTA_TRO_wd)
# Specify image output directory
imageOut <- 
  '/home/akoontz/Documents/QUTA_TxRedOaks/Documentation/Images/ManuscriptDraft_IJPS/'

# ---- %%% 2026-01-30 UPDATE %%% ----
# %%% READ IN AND FORMAT DATA ----
# Read in the VCF and convert it to a genind
QUTA_Clust3_vcf <- read.vcfR('Clust3_seOnly/pop_R98/populations.snps.vcf')
QUTA_Clust3_genlight <- vcfR2genlight(QUTA_Clust3_vcf)
QUTA_Clust3_genlight <- gl.compliance.check(QUTA_Clust3_genlight)
# Get names and (cleaned) dets
QUTA_Clust3_names <- QUTA_Clust3_genlight@ind.names
QUTA_Clust3_dets <- 
  as.factor(read.csv('Clust3_seOnly/pop_R98/Clust3_NamesAndDets.csv', header=TRUE)[,2])
QUTA_Clust3_genlight@pop <- QUTA_Clust3_dets 

# %%% SUBSET INDIVIDUALS %%% ----
# Removing hybrids (miquihuanensis- autopista ; hypoxantha x gravesii, hypo. aff sid.e), 
# miquihuaensis, scytophylla, one oddball canbyi
QUTA_Clust3_indsRemove <- 
  c('SYST-MOR-0006099', 'SYST-MOR-0006654', 'SYST-MOR-0007328', 'SYST-MOR-0007329', 'SYST-MOR-0007280',
    'OAK-MOR-000896', 'OAK-MOR-000907', 'OAK-MOR-001178','OAK-MOR-001234', 'OAK-UMN-000180',
    'OAK-UMN-000268','SYST-MOR-0007305')
# Subset genlight object
QUTA_Clust3_genlight_Subset <- 
  QUTA_Clust3_genlight[!indNames(QUTA_Clust3_genlight) %in% QUTA_Clust3_indsRemove, ]

# %%% RUN PCA ----
# Run PCA on genlight object, selecting 2 axes
QUTA_Clust3_PCA <- glPca(QUTA_Clust3_genlight_Subset, nf = 2)
# Extract PCA scores
pca_df <- as.data.frame(QUTA_Clust3_PCA$scores)
pca_df$Taxa <- QUTA_Clust3_genlight_Subset@pop
# Shape helper: tardifolia vs others
pca_df$ShapeGroup <- ifelse(
  pca_df$Taxa == "tardifolia",
  "tardifolia",
  "other"
)
# Extract variances, used to build labels for PC axes
var_exp <- round(
  100 * QUTA_Clust3_PCA$eig / sum(QUTA_Clust3_PCA$eig),
  1
)
pc_labels <- paste0("PC", 1:length(var_exp), " (", var_exp, "%)")

# %%% PLOTTING WITH LABELS ----
# Specify a colorblind friendly palatte, for 11 different colors
palette_11 <- c(
  "#000000", "#56B4E9", "#D55E00", "#CC79A7", "#F0E442",
  "#0072B2", "#7FBF7B", "#009E73",
  "#882255", "#999999", "#E69F00"
)

# Select one representative per taxon for labeling
pca_labels <- pca_df %>%
  group_by(Taxa) %>%
  slice(1)

# Main PCA plot
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(
    aes(color = Taxa, shape = ShapeGroup),
    size = 4.5,
    alpha = 0.9
  ) +
  scale_color_manual(
    name = "Taxa",
    values = palette_11,
    guide = guide_legend(
      override.aes = list(
        shape = ifelse(
          levels(factor(pca_df$Taxa)) == "tardifolia",
          17,  # triangle
          16   # circle
        ),
        size = 4.5
      )
    )
  ) +
  scale_shape_manual(
    values = c(
      "tardifolia" = 17,
      "other"      = 16
    ),
    guide = "none"
  ) +
  labs(
    title = "Q. tardifolia and Texas Red Oaks (Dataset #2)",
    subtitle = paste0(nInd(QUTA_Clust3_genlight_Subset)," individuals, 9,954 loci, R98"),
    x = pc_labels[1],
    y = pc_labels[2]
  ) +
  geom_text_repel(
    data = pca_labels,
    aes(label = Taxa),
    size = 6,
    family = "mono",
    max.overlaps = Inf,
    nudge_x = 0.07 * (max(pca_df$PC1) - min(pca_df$PC1)),
    nudge_y = 0.07 * (max(pca_df$PC2) - min(pca_df$PC2)),
    force = 4,
    box.padding = 1.6,
    point.padding = 1.3,
    segment.size = 0.8,
    show.legend = FALSE
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(size = 22),
    plot.subtitle = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = "right"
  )

# %%% SCREE PLOT ----
scree_df <- data.frame(
  PC = factor(
    paste0("PC", 1:length(var_exp)),
    levels = paste0("PC", 1:length(var_exp))
  ),
  Variance = var_exp
)

scree_plot <- ggplot(scree_df, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "gray60") +
  labs(
    title = "Scree Plot",
    x = NULL,
    y = "Variance (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(size = 12, hjust = 0.5),
    axis.title.y = element_text(size = 11)
  )

# %%% COMBINE PCA AND INSET SCREE PLOT ----
final_plot <- ggdraw() +
  draw_plot(pca_plot) +
  draw_plot(
    scree_plot,
    x = 0.03, y = 0.05,
    width = 0.23, height = 0.23
  )
final_plot

# Generate the final plot, as TIFF
tiff(file = paste0(imageOut, "Fig3_PCA.tiff"), width = 1200, height = 795)
final_plot
dev.off()

# %%% ALL INDIVIDUALS %%% ----
# %%% RUN PCA ----
# Run PCA on genlight object, selecting 2 axes
QUTA_Clust3_PCA <- glPca(QUTA_Clust3_genlight, nf = 2)
# Extract PCA scores
pca_df <- as.data.frame(QUTA_Clust3_PCA$scores)
pca_df$Taxa <- QUTA_Clust3_genlight@pop
# Shape helper: tardifolia vs others
pca_df$ShapeGroup <- ifelse(
  pca_df$Taxa == "tardifolia",
  "tardifolia",
  "other"
)
# Extract variances, used to build labels for PC axes
var_exp <- round(
  100 * QUTA_Clust3_PCA$eig / sum(QUTA_Clust3_PCA$eig),
  1
)
pc_labels <- paste0("PC", 1:length(var_exp), " (", var_exp, "%)")

# %%% PLOTTING WITH LABELS ----
# Specify a colorblind friendly palatte, for 17 different colors...
palette_17 <- c(
  "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7",
  "#999999", "#882255", "#44AA99", "#DDCC77", "#332288",
  "#AA4499", "#117733", "#661100", "#88CCEE"
)
# Select one representative per population for labeling
pca_labels <- pca_df %>%
  group_by(Taxa) %>%
  slice(1)  # Take the first occurrence per population

# Main PCA plot
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(
    aes(color = Taxa, shape = ShapeGroup),
    size = 4.5,
    alpha = 0.9
  ) +
  scale_color_manual(
    name = "Taxa",
    values = palette_17,
    guide = guide_legend(
      override.aes = list(
        shape = ifelse(
          levels(factor(pca_df$Taxa)) == "tardifolia",
          17,  # triangle
          16   # circle
        ),
        size = 4.5
      )
    )
  ) +
  scale_shape_manual(
    values = c(
      "tardifolia" = 17,
      "other"      = 16
    ),
    guide = "none"
  ) +
  labs(
    title = "Q. tardifolia and Texas Red Oaks (Dataset #2, Complete)",
    subtitle = paste0(nInd(QUTA_Clust3_genlight)," individuals, 9,954 loci, R98"),
    x = pc_labels[1],
    y = pc_labels[2]
  ) +
  geom_text_repel(
    data = pca_labels,
    aes(label = Taxa),
    size = 5,
    family = "mono",
    max.overlaps = Inf,
    nudge_x = 0.06 * (max(pca_df$PC1) - min(pca_df$PC1)),
    nudge_y = 0.06 * (max(pca_df$PC2) - min(pca_df$PC2)),
    force = 4,
    box.padding = 1.6,
    point.padding = 1.3,
    segment.size = 0.8,
    show.legend = FALSE
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(size = 22),
    plot.subtitle = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = "right"
  )

# Scree plot (minimal, no x-axis labels)
scree_df <- data.frame(
  PC = factor(paste0("PC", 1:length(var_exp)), levels = paste0("PC", 1:length(var_exp))),
  Variance = var_exp
)
scree_plot <- ggplot(scree_df, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "gray60") +
  labs(
    title = "Scree Plot",
    x = NULL,
    y = "Variance (%)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(size = 10, hjust = 0.5),
    axis.title.y = element_text(size = 9)
  )

# Combine PCA and inset scree plot
final_plot <- ggdraw() +
  draw_plot(pca_plot) +
  draw_plot(scree_plot, x = 0.03, y = 0.05, width = 0.23, height = 0.23)
final_plot

# Generate the final plot, as TIFF
tiff(file = paste0(imageOut, "FigS7_PCA-Complete.tiff"), width = 1200, height = 795)
final_plot
dev.off()

# ---- %%% 2026-06-11 UPDATE : NEW HYBRIDS DATASET (GRAV-HYPO-TAR 2) %%% ----
# %%% READ IN AND FORMAT DATA ----
# Read in the VCF and convert it to a genind
QUTA_Clust3_GravHypoTar2_vcf <- read.vcfR('Clust3_seOnly/GravHypoTar2_R98_WL10K/populations.snps.vcf')
QUTA_Clust3_GravHypoTar2_genlight <- vcfR2genlight(QUTA_Clust3_GravHypoTar2_vcf)
QUTA_Clust3_GravHypoTar2_genlight <- gl.compliance.check(QUTA_Clust3_GravHypoTar2_genlight)
# Get names and (cleaned) dets
QUTA_Clust3_GravHypoTar2_names <- QUTA_Clust3_GravHypoTar2_genlight@ind.names
QUTA_Clust3_GravHypoTar2_dets <- 
  as.factor(read.csv('Clust3_seOnly/GravHypoTar2_R98_WL10K/Clust3_GravHypoTar2_NamesAndDets.csv', header=TRUE)[,2])
QUTA_Clust3_GravHypoTar2_genlight@pop <- QUTA_Clust3_GravHypoTar2_dets 

# %%% RUN PCA ----
# Run PCA on genlight object, selecting 2 axes
QUTA_Clust3_GravHypoTar2_PCA <- glPca(QUTA_Clust3_GravHypoTar2_genlight, nf = 2)
# Extract PCA scores
pca_df <- as.data.frame(QUTA_Clust3_GravHypoTar2_PCA$scores)
pca_df$Taxa <- QUTA_Clust3_GravHypoTar2_genlight@pop
# Extract variances, used to build labels for PC axes
var_exp <- round(
  100 * QUTA_Clust3_GravHypoTar2_PCA$eig / sum(QUTA_Clust3_GravHypoTar2_PCA$eig),
  1
)
pc_labels <- paste0("PC", 1:length(var_exp), " (", var_exp, "%)")

# %%% PLOTTING WITH LABELS ----
# Specify a colorblind friendly palatte, for 11 different colors
palette_11 <- c(
  "#000000", "#56B4E9", "#D55E00", "#CC79A7", "#F0E442",
  "#0072B2", "#7FBF7B", "#009E73",
  "#882255", "#999999", "#E69F00"
)

# Main PCA plot
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(
    aes(color = Taxa),
    size = 4.5,
    alpha = 0.9
  ) +
  scale_color_manual(
    name = "Taxa",
    values = palette_11,
    guide = guide_legend(
      override.aes = list(
        shape = ifelse(
          levels(factor(pca_df$Taxa)) == "tardifolia",
          17,  # triangle
          16   # circle
        ),
        size = 4.5
      )
    )
  ) +
  scale_shape_manual(
    values = c(
      "tardifolia" = 17,
      "other"      = 16
    ),
    guide = "none"
  ) +
  labs(
    title = "Hybrid vs. Parental Comparisons Dataset",
    subtitle = paste0(nInd(QUTA_Clust3_GravHypoTar2_genlight)," individuals, R98, 10,000 loci"),
    x = pc_labels[1],
    y = pc_labels[2]
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(size = 22),
    plot.subtitle = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = "right"
  )

# %%% SCREE PLOT ----
scree_df <- data.frame(
  PC = factor(
    paste0("PC", 1:length(var_exp)),
    levels = paste0("PC", 1:length(var_exp))
  ),
  Variance = var_exp
)

scree_plot <- ggplot(scree_df, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "gray60") +
  labs(
    title = "Scree Plot",
    x = NULL,
    y = "Variance (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(size = 12, hjust = 0.5),
    axis.title.y = element_text(size = 11)
  )

# %%% COMBINE PCA AND INSET SCREE PLOT ----
final_plot <- ggdraw() +
  draw_plot(pca_plot) +
  draw_plot(
    scree_plot,
    x = 0.68, y = 0.05,
    width = 0.23, height = 0.23
  )
final_plot

# Generate the final plot, as TIFF
tiff(file = paste0(imageOut, "revisedSubmission/FigS9_PCA.tiff"), width = 1200, height = 795)
final_plot
dev.off()
