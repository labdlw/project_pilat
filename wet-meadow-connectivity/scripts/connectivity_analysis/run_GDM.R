############################################################
### GENERALIZED DISSIMILARITY MODE ANALYSIS SCRIPT
### Purpose: Run GDM to explore genetic-environment relationships
### Input: FST matrix, coordinates, resistance matrix, environmental rasters
### Output: Observed vs predicted plot, variable importance plot, GDM summary
############################################################

### 01. LOAD PACKAGES ----
library(gdm)
library(raster)
library(cluster)
library(caret)
library(tidyverse)

### 02. LOAD INPUT DATA ----
# Paths (edit for each species/project)
genetic_fst_file <- "file/path/to/fst_matrix.csv"
coords_file <- "file/path/to/occurrence_points.csv"
resistance_matrix_file <- "file/path/to/resistance_matrix.txt"
env_raster_dir <- "file/path/to/environmental_rasters/"

# Load FST matrix and convert to DPS
fst_matrix <- read.csv(genetic_fst_file, row.names = 1)
fst_matrix <- as.matrix(fst_matrix)
diag(fst_matrix) <- 0
fst_matrix[upper.tri(fst_matrix)] <- t(fst_matrix)[upper.tri(fst_matrix)]
genetic_dist <- fst_matrix / (1 - fst_matrix)
genetic_dist <- as.data.frame(genetic_dist)
genetic_dist$site <- rownames(genetic_dist)

# Load coordinates
occ_op_pop <- read.csv(coords_file)
colnames(occ_op_pop) <- c("site", "Long", "Lat")

# Load resistance distance matrix
resist_table <- read.table(resistance_matrix_file, header = TRUE)
resist_table[1:2] <- lapply(resist_table[1:2], factor, levels = unique(resist_table$Id1))
res_matrix <- xtabs(Distance ~ Id1 + Id2, resist_table)
res_matrix <- res_matrix + t(res_matrix)
resistance_dist <- as.matrix(as.dist(res_matrix))
rownames(resistance_dist) <- colnames(resistance_dist) <- unique(resist_table$Id1)

### 03. PROCESS ENVIRONMENTAL VARIABLES ----
# Load rasters (assumes files follow naming conventions)
env_stack <- stack(list.files(env_raster_dir, pattern = ".asc$", full.names = TRUE))
occ_values <- as.data.frame(extract(env_stack, occ_op_pop[, c("Long", "Lat")]))
# Remove highly correlated predictors
cor_mat <- cor(na.omit(occ_values))
hc <- findCorrelation(cor_mat, cutoff = 0.7)
reduced_data <- occ_values[, -hc]
reduced_data$site <- occ_op_pop$site
reduced_data$Long <- occ_op_pop$Long
reduced_data$Lat <- occ_op_pop$Lat

### 04. FORMAT DATA FOR GDM ----
# Merge genetic and spatial data
genetic_dist <- merge(genetic_dist, occ_op_pop, by = "site")

# Create site-pair table (environmental distances)
gdmTab_env <- formatsitepair(
  bioData = genetic_dist,
  bioFormat = 3,
  XColumn = "Long",
  YColumn = "Lat",
  predData = reduced_data,
  distPreds = as.matrix(daisy(reduced_data[, -which(names(reduced_data) %in% c("site", "Long", "Lat"))], metric = "gower")),
  siteColumn = "site"
)
gdmTab_env <- na.omit(gdmTab_env)

# Create site-pair table (resistance distances)
gdmTab_resist <- formatsitepair(
  bioData = genetic_dist,
  bioFormat = 3,
  XColumn = "Long",
  YColumn = "Lat",
  predData = reduced_data,
  distPreds = list(as.matrix(resistance_dist)),
  siteColumn = "site"
)
gdmTab_resist <- na.omit(gdmTab_resist)

### 05. RUN GDM MODELS ----
# Model with environmental variables
gdm_env <- gdm(data = gdmTab_env, geo = TRUE)

# Model with resistance matrix
gdm_resist <- gdm(data = gdmTab_resist, geo = TRUE)

# Summary and explained deviance
cat("Explained Deviance (Environment-only model):", round(gdm_env$explained, 3), "\n")
summary(gdm_env)

### 06. VARIABLE IMPORTANCE ----
gdm_imp <- gdm.varImp(gdmTab_env, geo = TRUE, nPerm = 500, parallel = TRUE, cores = 4)
imp_scores <- gdm_imp[[2]]
imp_sorted <- imp_scores[order(imp_scores$`Importance`, decreasing = TRUE), ]

# Plot and save
pdf("GDM_VariableImportance.pdf", width = 7, height = 5)
barplot(imp_sorted$Importance, names.arg = rownames(imp_sorted), las = 2,
        main = "GDM Variable Importance", col = "steelblue")
dev.off()

### 07. OBSERVED VS PREDICTED PLOT ----
gdm_pred <- predict(gdm_env, gdmTab_env)

pdf("GDM_Observed_vs_Predicted.pdf")
plot(gdmTab_env$distance,
     gdm_pred,
     xlab = "Observed Genetic Dissimilarity",
     ylab = "Predicted Dissimilarity",
     xlim = c(0, max(gdmTab_env$distance)),
     ylim = c(0, max(gdm_pred)),
     pch = 20,
     col = rgb(0, 0, 1, 0.5),
     main = "GDM Observed vs Predicted")
abline(0, 1, col = "red", lty = 2)
dev.off()

### 08. SAVE RESULTS ----
saveRDS(gdm_env, file = "gdm_env_model.rds")
write.csv(imp_sorted, "GDM_VariableImportance.csv", row.names = TRUE)
