############################################################
### GENERALIZED DISTANCE OF MAXIMUM CORRELATION ANALYSIS SCRIPT
### Purpose: Estimate Distance of Maximum Correlation (DMC)
### Input: genind object, cost distance matrix
### Output: DMC plot, correlation curve, scatter plot
############################################################

### 01. LOAD PACKAGES ----
library(tidyverse)
library(adegenet)
library(graph4lg)  # required for mat_gen_dist, dist_max_corr, scatter_dist
library(readr)

### 02. USER INPUTS ----
# Define input paths
genind_file <- "file/path/to/genind_object.rds"        # RDS format
cost_dist_file <- "file/path/to/cost_distance_matrix.txt"

# Optional parameters
species_name <- "Your species name here"
output_prefix <- "DMC_result"  # used for saved files
bin_interval <- 100            # interv argument

# Color settings
point_color <- "black"
vline_color <- "black"

### 03. LOAD DATA ----
# Load genind object
genind_obj <- readRDS(genind_file)

# Compute genetic distance matrix (DPS)
mat_dps <- mat_gen_dist(x = genind_obj, dist = "DPS")

# Load cost distance matrix from Graphab
cost_df <- read_table(cost_dist_file)
cost_df[1:2] <- lapply(cost_df[1:2], factor, levels = unique(cost_df$Id1))
m_cost <- xtabs(Distance ~ Id1 + Id2, cost_df)
mat_cd <- as.matrix(as.dist(m_cost + t(m_cost)))
rownames(mat_cd) <- colnames(mat_cd) <- unique(cost_df$Id1)

### 04. CLEAN & ALIGN MATRICES ----
# Match population order between matrices
common_pops <- intersect(rownames(mat_cd), rownames(mat_dps))
mat_cd <- mat_cd[common_pops, common_pops]
mat_dps <- mat_dps[common_pops, common_pops]

### 05. RUN DMC ----
dmc <- dist_max_corr(
  mat_gd = mat_dps,
  mat_ld = mat_cd,
  interv = bin_interval,
  pts_col = point_color
)

dmc_dist <- dmc$`distance at which correlation reaches a maximum`
cat("✅ DMC identified at distance:", dmc_dist, "\n")

### 06. PLOT RESULTS ----

# Correlation curve
dmc_data <- data.frame(
  Distance = dmc[[3]],
  Correlation = dmc[[2]]
)

p_dmc <- ggplot(dmc_data, aes(x = Distance, y = Correlation)) +
  geom_line(color = "orange", size = 0.5) +
  geom_vline(xintercept = dmc_dist, linetype = "dashed", color = vline_color) +
  ggtitle(paste("DMC for", species_name)) +
  theme_minimal() +
  labs(x = "Distance threshold (cost units)", y = "Correlation Coefficient")

# Scatter plot of raw genetic vs cost distances
p_scatter <- scatter_dist(mat_gd = mat_dps, mat_ld = mat_cd, pts_col = point_color)

# Save plots
ggsave(paste0(output_prefix, "_DMC_lineplot.pdf"), p_dmc, width = 7, height = 5)
ggsave(paste0(output_prefix, "_DMC_scatterplot.pdf"), plot = p_scatter, width = 6, height = 5)

# Save DMC values
write.csv(dmc_data, paste0(output_prefix, "_DMC_curve.csv"), row.names = FALSE)
writeLines(as.character(dmc_dist), paste0(output_prefix, "_DMC_distance.txt"))

### 07. DONE ----
cat("✅ DMC analysis complete. Plots and results saved.\n")
