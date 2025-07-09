##########################################################################################################
######## PLS analysis using residuals between euclidean distance and the distance of cumulative cost #####
##########################################################################################################
setwd("/home/alle/Documents/Pilat_project/")
library(graph4lg)
library(sf)


proj_MM5 <- "PHPi_ModMoy5sp"

metrics_mm5 = get_graphab_metric(proj_name = proj_MM5) #metrics = 16445

#load in population coordinates not indivdiual
occ_lf <- read.csv("/home/alle/Documents/Pilat_project/process_all/processing_lf/denovom3_M1_n1/lf_filter_vcfs/samples_coord_scheme2_lambert93_by_POP.csv")

colnames(occ_lf)<- c("ID", "x", "y")
pointset <- graphab_pointset(proj_name = proj_name,
                             linkset = "Jeulien_cout",
                             pointset = occ_lf) #using coordinates from LF or OP populations for the PLS later on. 

lm_500 = glm(metrics_mm5$F_d388_p0.5_beta1_Graph_C_388 ~ metrics_mm5$F_d500_p0.5_beta1_Graph_E_500)
lm_1000 = glm(metrics_mm5$F_d1265_p0.5_beta1_Graph_C_1265 ~ metrics_mm5$F_d1000_p0.5_beta1_Graph_E_1000)
lm_2000 = glm(metrics_mm5$F_d4118_p0.5_beta1_Graph_C_4118 ~ metrics_mm5$F_d2000_p0.5_beta1_Graph_E_2000)
lm_3000 = glm(metrics_mm5$F_d8215_p0.5_beta1_Graph_C_8215 ~ metrics_mm5$F_d3000_p0.5_beta1_Graph_E_3000)
lm_5000 = glm(metrics_mm5$F_d19604_p0.5_beta1_Graph_19604 ~ metrics_mm5$F_d5000_p0.5_beta1_Graph_E_5000)
lm_10000 = glm(metrics_mm5$F_d63816_p0.5_beta1_Graph_C_63816 ~ metrics_mm5$F_d10000_p0.5_beta1_Graph_E_10000)

# these residuals are from the general model so they are not specific to any species.
residuals_500 = residuals(lm_500)
residuals_1000 = residuals(lm_1000)
residuals_2000 = residuals(lm_2000)
residuals_3000 = residuals(lm_3000)
residuals_5000 = residuals(lm_5000)
residuals_10000 = residuals(lm_10000)

pointset <- as.data.frame(pointset)
colnames(pointset)[1] <- "patch_ID"
# adding the patch ID column to the residuals data frames so that they can be merged together with the pointset of each species. 
residuals_500 = as.data.frame(residuals_500)
residuals_500$patch_ID = row.names(residuals_500)

residuals_1000 = as.data.frame(residuals_1000)
residuals_1000$patch_ID = row.names(residuals_1000)

residuals_2000 = as.data.frame(residuals_2000)
residuals_2000$patch_ID = row.names(residuals_2000)

residuals_3000 = as.data.frame(residuals_3000)
residuals_3000$patch_ID = row.names(residuals_3000)

residuals_5000 = as.data.frame(residuals_5000)
residuals_5000$patch_ID = row.names(residuals_5000)

residuals_10000 = as.data.frame(residuals_10000)
residuals_10000$patch_ID = row.names(residuals_10000)

#here we then merge this common column called patch_ID to the pointset which is extracted seprately for each species (different coordinates for different species)
pointset_res = merge(pointset, residuals_500, by = "patch_ID")
pointset_res = merge(pointset_res, residuals_1000, by = "patch_ID")
pointset_res = merge(pointset_res, residuals_2000, by = "patch_ID")
pointset_res = merge(pointset_res, residuals_3000, by = "patch_ID")
pointset_res = merge(pointset_res, residuals_5000, by = "patch_ID")
pointset_res = merge(pointset_res, residuals_10000, by = "patch_ID")
#here we remove irrelevent metrics, keeping only residual columns 
pointset_res_only <- pointset_res[,-1:-41]
pointset_res_only <- cbind(pointset_res_only, pointset_res$Capacity)
pointset_res_only <- cbind(pointset_res_only, pointset_res$Point_ID)
colnames(pointset_res_only)[8] <- "ID"
colnames(pointset_res_only)[7] <- "Capacity"

# load in genetic diversity data frame 
div_sh <- read.table("/home/alle/Documents/Pilat_project/process_all/processing_sh/resequenced_corrected_sh/scorzonere_diversite_populations_stacks.tsv", header = T)

div_lf <- read.csv("/home/alle/Documents/Pilat_project/lychnis_diversite_populations_stacks.csv")
colnames(div_lf)[1] <- "ID"

div_op <- read.csv("/home/alle/Downloads/oenanthe_diversite_populations_stacks.csv")
colnames(div_op)[1] <- "ID"

# here we add in the diversity measures and merge them by the column ID
div_connect <- merge(pointset_res_only, div_op[ c(1,21, 24, 15)], by = "ID")
write.table(div_connect, "/home/alle/Documents/Pilat_project/process_all/processing_op/diversity_connectivity_residuals_OP_table.tsv")

################################################################
vary <- c("Pi", "Exp_Het", "Fis")
varx <- c("Capacity", "residuals_500", "residuals_1000", "residuals_2000" ,  "residuals_3000" ,  "residuals_5000", "residuals_10000")
tab_pls2 <- div_connect[, c(vary, varx)]

ny <- length(vary)
nc <- ncol(tab_pls2) - ny
library(plsdepot)

model.pls <- plsreg2(tab_pls2[, c((ny+1):ncol(tab_pls2))], tab_pls2[, 1:ny], comps=nc, crosval=TRUE)

# Q2 prediction
Q2 <- model.pls$Q2
# R2 fit
R2 <- model.pls$expvar
# Variable weights
val_wgs2 <- model.pls$raw.wgs^2
# Variable correlation with components 
val_cor_xt <- model.pls$cor.xt

tab_coord <- rbind(data.frame(model.pls$cor.xt), data.frame(model.pls$cor.yt))
tab_coord$var <- row.names(tab_coord)
tab_coord$col <- "1"
tab_coord[which(tab_coord$var %in% vary), 'col'] <- "2"

library(ggplot2)
library(ggrepel)
biplot <- ggplot() +
  geom_hline(yintercept = 0, color = "#999999") +
  geom_vline(xintercept = 0, color = "#999999") +
  geom_segment(data = tab_coord, 
               aes(color = col, x = 0, y = 0, xend = t1, yend = t2),
               arrow = arrow(length = unit(0.1, "inches"))) +
  geom_text_repel(data = tab_coord, 
                  aes(color = col, label = var, x = t1, y = t2)) +
  geom_path(data = cercle, aes(x = x, y = y)) +
  theme_bw() +
  scale_color_manual(values = c("black", "darkblue")) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, .2)) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, .2)) +
  theme(legend.position = "none") +
  labs(x = paste0("t1 - R2x = ", round(model.pls$expvar[1, 1], digits = 2), 
                  " - R2y = ", round(model.pls$expvar[1, 3], digits = 2), 
                  #" - Q2 = ", round(model.pls$Q2[1, 1], digits = 2), " - ", round(model.pls$Q2[1, 2], digits = 2), " - ", round(model.pls$Q2[1, 3], digits = 2)), 
                  " - Q2 = ", round(model.pls$Q2[1, 3], digits = 2)), 
       y = paste0("t2 - R2x = ", round(model.pls$expvar[2, 1], digits = 2), 
                  " - R2y = ", round(model.pls$expvar[2, 3], digits = 2), 
                  " - Q2 = ", round(model.pls$Q2[2, 3], digits = 2))) +
  ggtitle("PLS of Scorzonera humilis")

biplot
