############################################################
### GENERALIZED MANTEL & PARTIAL MANTEL ANALYSIS SCRIPT
### Purpose: Run Mantel and Partial Mantel tests for IBD, IBR, and IBE
### Input: VCF file, coordinate file, resistance, environmental and geo distances
### Output: Mantel coefficients, p-values, and plots
############################################################

### 01. LOAD PACKAGES ----
library(adegenet)
library(StAMPP)
library(vcfR)
library(cluster)
library(caret)
library(graph4lg)
library(raster)
library(readr)
library(vegan)
library(ggplot2)

### 02. USER INPUT ----
vcf_path <- "file/path/to/vcf_file.vcf"
coords_path <- "file/path/to/occurrence_coordinates.csv"
geo_dist_path <- "file/path/to/geographic_distance_matrix.csv"
cost_dist_path <- "file/path/to/cost_distance_matrix.txt"
env_raster_dir <- "file/path/to/environmental_rasters/"
species_name <- "Oenanthe peucedanifolia"

### 03. READ AND CLEAN GENETIC DATA ----
vcf <- read.vcfR(vcf_path)
genlight_obj <- vcfR2genlight(vcf)

# Derive population labels
site <- sapply(strsplit(indNames(genlight_obj), '_'), `[`, 1)
pop <- sapply(strsplit(indNames(genlight_obj), '_'), `[`, 3)
genlight_obj$pop <- factor(paste(site, pop, sep = "_"))

# Remove individuals with missing coordinates (customize if needed)
remove_inds <- c("p_op_4_9", "p_op_4_13", "p_op_4_4", "p_op_4_11", "p_op_4_12", "p_op_4_14", "p_op_4_15", "p_op_4_8", "p_op_4_3", "p_op_4_5")
genlight_obj <- genlight_obj[!indNames(genlight_obj) %in% remove_inds]

# Convert to FST matrix and transform
stampp <- stamppConvert(genlight_obj, type = "genlight")
fst_matrix <- stamppFst(stampp, nboots = 100, percent = 95, nclusters = 8)$Fsts
diag(fst_matrix) <- 0
fst_matrix[upper.tri(fst_matrix)] <- t(fst_matrix)[upper.tri(fst_matrix)]
genetic_dist <- as.matrix(fst_matrix / (1 - fst_matrix))

### 04. LOAD DISTANCE MATRICES ----
# Geographic
geo_dist <- read.csv(geo_dist_path, row.names = 1)
geo_dist <- as.matrix(geo_dist)

# Cost distance (Graphab output)
cost_df <- read_table(cost_dist_path)
cost_df[1:2] <- lapply(cost_df[1:2], factor, levels = unique(cost_df$Id1))
cost_matrix <- xtabs(Distance ~ Id1 + Id2, cost_df)
cost_matrix <- cost_matrix + t(cost_matrix)
cost_matrix <- as.matrix(cost_matrix)

# Environmental dissimilarity

######## fetch env ditance
Alti_RLA <- raster("Alti_RLA_HL_50_A.asc")
BDForet <- raster("BDForet_RLA_HL_50_A.asc")
CHELSEA_1 <- raster("CHELSEA_01_RLA_HL_50_A.asc")
CHELSEA_2 <- raster("CHELSEA_02_RLA_HL_50_A.asc")
CHELSEA_3 <- raster("CHELSEA_03_RLA_HL_50_A.asc")
CHELSEA_4 <- raster("CHELSEA_04_RLA_HL_50_A.asc")
CHELSEA_5 <- raster("CHELSEA_05_RLA_HL_50_A.asc")
CHELSEA_6 <- raster("CHELSEA_06_RLA_HL_50_A.asc")
CHELSEA_7 <- raster("CHELSEA_07_RLA_HL_50_A.asc")
CHELSEA_8 <- raster("CHELSEA_08_RLA_HL_50_A.asc")
CHELSEA_9 <- raster("CHELSEA_09_RLA_HL_50_A.asc")
CHELSEA_10 <- raster("CHELSEA_10_RLA_HL_50_A.asc")
CHELSEA_11 <- raster("CHELSEA_11_RLA_HL_50_A.asc")
CHELSEA_12 <- raster("CHELSEA_12_RLA_HL_50_A.asc")
CHELSEA_13 <- raster("CHELSEA_13_RLA_HL_50_A.asc")
CHELSEA_14 <- raster("CHELSEA_14_RLA_HL_50_A.asc")
CHELSEA_15 <- raster("CHELSEA_15_RLA_HL_50_A.asc")
CHELSEA_16 <- raster("CHELSEA_16_RLA_HL_50_A.asc")
CHELSEA_17 <- raster("CHELSEA_17_RLA_HL_50_A.asc")
CHELSEA_18 <- raster("CHELSEA_18_RLA_HL_50_A.asc")
CHELSEA_19 <- raster("CHELSEA_19_RLA_HL_50_A.asc")
dist_eau <- raster("Dist_eau_RLA_HL_50_A.asc")
Dist_route <- raster("Dist_routes_RLA_HL_50_A.asc")
Expo <- raster("Expo_RLA_HL_50_A.asc")
Occsol <- raster("Occsol_CESBIO_RLA_HL_50_A.asc")
ombrage <- raster("Ombrage_RLA_HL_50_A.asc")
pentes <- raster("Pentes_RLA_HL_50_A.asc")
ph <- raster("pH_RLA_HL_50_A.asc")
RPG <- raster("RPG_RLA_HL_50_A.asc")

env_raster_dir  <- stack(
  Alti_RLA,
  CHELSEA_1,
  CHELSEA_2,
  CHELSEA_3,
  CHELSEA_4,
  CHELSEA_5,
  CHELSEA_6,
  CHELSEA_7,
  CHELSEA_8,
  CHELSEA_9,
  CHELSEA_10,
  CHELSEA_11,
  CHELSEA_12,
  CHELSEA_13,
  CHELSEA_14,
  CHELSEA_15,
  CHELSEA_16,
  CHELSEA_17,
  CHELSEA_18,
  CHELSEA_19,
  dist_eau,
  Dist_route,
  Expo,
  ombrage,
  pentes,
  ph
)


raster_stack <- stack(list.files(env_raster_dir, pattern = ".asc$", full.names = TRUE))
occ_coords <- read.csv(coords_path)
env_values <- as.data.frame(extract(raster_stack, occ_coords[, c("Long", "Lat")]))
hc <- findCorrelation(cor(na.omit(env_values)), cutoff = 0.7)
reduced_env <- env_values[, -hc]
env_dist <- daisy(reduced_env, metric = "gower")
env_dist <- as.matrix(env_dist)

# Align all matrices by row/col names
shared_pops <- Reduce(intersect, list(rownames(genetic_dist), rownames(geo_dist), rownames(cost_matrix)))
genetic_dist <- genetic_dist[shared_pops, shared_pops]
geo_dist <- geo_dist[shared_pops, shared_pops]
cost_matrix <- cost_matrix[shared_pops, shared_pops]
env_dist <- env_dist[shared_pops, shared_pops]

### 05. MANTEL TESTS (GLOBAL) ----
cat("Global Mantel Tests\n")
cat("Genetic ~ Geographic:\n")
print(mantel.rtest(as.dist(genetic_dist), as.dist(geo_dist), nrepet = 9999))
cat("Genetic ~ Cost:\n")
print(mantel.rtest(as.dist(genetic_dist), as.dist(cost_matrix), nrepet = 9999))
cat("Genetic ~ Environment:\n")
print(mantel.rtest(as.dist(genetic_dist), as.dist(env_dist), nrepet = 9999))

### 06. PLOT GLOBAL CORRELATIONS ----
scatter_dist(genetic_dist, geo_dist, method = "lm") +
  labs(title = paste(species_name, "- Genetic vs Geographic"), x = "Geographic Distance", y = "FST/(1-FST)") +
  theme(plot.title = element_text(face = "italic"))

scatter_dist(genetic_dist, cost_matrix, method = "lm") +
  labs(title = paste(species_name, "- Genetic vs Cost Distance"), x = "Cost Distance", y = "FST/(1-FST)") +
  theme(plot.title = element_text(face = "italic"))

scatter_dist(genetic_dist, env_dist, method = "lm") +
  labs(title = paste(species_name, "- Genetic vs Environmental Distance"), x = "Environmental Dissimilarity", y = "FST/(1-FST)") +
  theme(plot.title = element_text(face = "italic"))

### 07. PARTIAL MANTEL TESTS ----
cat("\nPartial Mantel Tests (controlling for geographic distance)\n")
cat("Genetic ~ Cost | Geo:\n")
print(mantel.partial(as.dist(genetic_dist), as.dist(cost_matrix), as.dist(geo_dist), permutations = 9999))

cat("Genetic ~ Env | Geo:\n")
print(mantel.partial(as.dist(genetic_dist), as.dist(env_dist), as.dist(geo_dist), permutations = 9999))

### 08. OPTIONAL: SUBSET PER PLATEAU ----
# If you want to subset matrices by plateau (A, M, P), define site codes per region:
plateaus <- list(
  A = c("a_1", "a_2", "a_3", "a_4", "a_5", "a_6", "a_7", "a_8", "a_9", "a_10"),
  M = c("m_1", "m_2", "m_3", "m_4", "m_5", "m_6", "m_7", "m_8", "m_9", "m_10"),
  P = c("p_1", "p_2", "p_3", "p_5")
)

for (region in names(plateaus)) {
  pops <- plateaus[[region]]
  if (all(pops %in% rownames(genetic_dist))) {
    cat(paste0("\n--- Mantel Tests for plateau ", region, " ---\n"))
    gd <- as.dist(genetic_dist[pops, pops])
    gd_geo <- as.dist(geo_dist[pops, pops])
    gd_cost <- as.dist(cost_matrix[pops, pops])
    gd_env <- as.dist(env_dist[pops, pops])
    print(mantel.rtest(gd, gd_geo, nrepet = 9999))
    print(mantel.rtest(gd, gd_cost, nrepet = 9999))
    print(mantel.rtest(gd, gd_env, nrepet = 9999))
  }
}

### 09. DONE ----
cat("\n✅ Mantel analysis complete.\n")
