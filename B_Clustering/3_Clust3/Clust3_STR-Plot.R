# %%%%%%%%%%%%%%%%%%%%%%
# %%% PLOT STRUCTURE %%%
# %%%%%%%%%%%%%%%%%%%%%%

# This script generates a STRUCTURE-like stacked barplot, with labels along the right hand side.
# It uses the Q matrix from Clust3, K=4 (main output; NOT Summary clusters). The CSV it reads in 
# is an adapted version of that Q. matrix, which includes updated dets for each sample. Samples
# within that CSV are also rearranged from the original order in the STRUCTURE run, to make 
# the presentation of the dets for samples simpler.

library(tidyverse)

# ---------------- VARIABLES ----------------
# Variable file path: directory containing Q matrix file to be read in.
QUTA_Clust3_inputCSV <- 
  "/home/akoontz/Documents/QUTA_TxRedOaks/Code/B_Clustering/3_Clust3/Clust3_K4_Qmat_CleanDets.csv"
# Variable file path: directory to save image to
QUTA_Clust3_outputPlotDir <- 
  "/home/akoontz/Documents/QUTA_TxRedOaks/Documentation/Images/ManuscriptDraft3/"
plot_title <- "STRUCTURE Results (K4)"  
plot_subtitle <- "R98, 84 individuals, 9,954 loci"

# Specify plotting colors
cluster_colors <- c(
  "#08306B", # dark blue
  "#6BAED6", # light blue
  "#A1D99B", # light green
  "#006D2C" # dark green
)

# Width and position of the classification guide on the right
segment_x <- 1.005          # x position of the vertical line segments
text_x    <- 1.008         # x position of the labels (just to the right)
# Specify text size for labels
labelSize <- 5
# Specify right hand margin
right_margin_mm <- 120 # size of free text placed between segments

# ---- SPECIFY CLASSIFICATION SEGMENTS ----
# ymin / ymax are row numbers (after ordering) defining each block
# label is the text shown for that block
class_segments <- tibble(
  ymin  = c(1, 16, 18.2, 24.2, 27, 36.5, 45, 53.2, 55.1, 64, 70.1, 76.1),
  ymax  = c(14, 16.5, 22.5, 25, 34.8, 43, 50, 53.8, 62, 68.8, 74.8, 83),
  label = c("Q. gravesii", "Q. gravesii", "Q. (emoryi x gravesii)", "Q. (emoryi x gravesii) blobbier",
            "Q. emoryi", "Q. graciliformis", "Q. canbyi", "Q. miquihuaensis", 
            "Q. hypoleucoides", "Q. scytophylla", "Q. (gravesii x hypoleucoides)", "Q. tardifolia")
)

# ---- SPECIFY FREE TEXT SEGMENTS ----
# Each row = one independent text block (no line drawn)
# ymin / ymax define the vertical span; text is centered between them
# This is used for annotations BETWEEN line segments
free_text_blocks <- tibble(
  ymin = c(14.8, 26, 44, 50.7, 52.0, 63, 84),
  ymax = c(14.8, 26, 44, 50.7, 52.0, 63, 84),
  label = c("Q. sp?", "Q. emoryi x ??", "Q. robusta", "Q. cf. miquihuanensis- autopista",
            "Q. hypoxantha", "Q. aff hypoleucoides aff. sideoxyla", "Q. hypoxantha x Q. gravesii")
)

# ---------------- READ & PREP DATA ----------------
df <- read.csv(QUTA_Clust3_inputCSV, check.names = FALSE)
# Specify column names
sample_col <- "Sample name"
class_col  <- "Updated Det (2026-01-13)"
# Identify Q columns automatically
q_cols <- grep("^Q", colnames(df), value = TRUE)
# Keep order as in file (top = first row, bottom = last row)
df <- df %>% mutate(row_id = row_number())
# Long format for plotting
df_long <- df %>%
  pivot_longer(cols = all_of(q_cols), names_to = "Cluster", values_to = "Q")

# ---------------- PLOT ----------------
p <- ggplot(df_long, aes(x = Q, y = factor(row_id, levels = unique(row_id)), fill = Cluster)) +
  geom_col(width = 0.9) +
  scale_fill_manual(values = cluster_colors) +
  scale_x_continuous(
    limits = c(0, 1.10),
    expand = c(0, 0),
    oob = scales::oob_keep
  ) +
  labs(x = "Q value", y = NULL, title = plot_title, subtitle = plot_subtitle) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title    = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 16),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5.5, right_margin_mm, 5.5, 5.5)
  ) +
  coord_cartesian(clip = "off")

# ---------------- ADD USER-DEFINED CLASSIFICATION SEGMENTS ----------------
# Vertical line segments
p <- p +
  geom_segment(
    data = class_segments,
    inherit.aes = FALSE,
    aes(x = segment_x, xend = segment_x, y = ymin - 0.5, yend = ymax + 0.5),
    linewidth = 0.8
  )

# Labels to the right of each line segment (centered vertically)
p <- p +
  geom_text(
    data = class_segments,
    inherit.aes = FALSE,
    aes(x = text_x, y = (ymin + ymax) / 2, label = label),
    hjust = 0,
    size = labelSize,
    fontface = "italic"
  )

# Free floating labels
p <- p +
  geom_text(
    data = free_text_blocks,
    inherit.aes = FALSE,
    aes(x = text_x, y = (ymin + ymax) / 2, label = label),
    hjust = 0,
    size = labelSize,
    fontface = "italic"
  )
print(p)

# Save the image to a PNG
png(file = paste0(QUTA_Clust3_outputPlotDir, "Fig2_STR_K4.png"), width = 1200, height = 795)
print(p)
dev.off()
