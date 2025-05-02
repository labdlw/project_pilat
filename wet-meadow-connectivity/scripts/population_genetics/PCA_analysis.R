# This script containes the code to run the PCA plot for the oenanthe species. the same was done for the other two. 

library(tibble)
library(vcfR)
library(adegenet)
library(FactoMineR)
library(factoextra)
#load vcf file 
vcf_op_262 <-
  read.vcfR(
    "/home/alle/Documents/Pilat_project/process_all/processing_op/corrected_op/denovo_op_m3_M5_n5/op_GQ20_minDP5_minmeanDP10_miss080_mac3_ind262.recode.vcf"
  )

#convert vcf to genind object
genind_vcf_op_262 <- vcfR2genind(vcf_op_262)


#add pop tab to genind object 
popmap_op <- read_table2("/home/alle/Documents/Pilat_project/process_all/processing_op/corrected_op/denovo_op_m3_M5_n5/popmap_op.txt", col_names = FALSE)
colnames(LQ_indv_op50) = "X1"
LQ_indv_op50 = as.data.frame(LQ_indv_op50)
popmap_op_262 = setdiff(popmap_op[,1],LQ_indv_op50)
popmap_op_262 = as.data.frame(popmap_op_262)
pop <- sapply(strsplit(popmap_op_262$X1, '_'),
              function(popmap_op_262){popmap_op_262[1]})
pop = as.factor(pop)
popmap_op_262$X2 = pop 
genind_vcf_op_262$pop = pop

#write the new popmap file with just 331 individuals out of 388
write_tsv(popmap_op_262, "/home/alle/Documents/Pilat_project/process_all/processing_op/corrected_op/denovo_op_m3_M5_n5/popmap_op_262.txt")

#scale genind for pca
genind_vcf_op_262_scaled = scaleGen(genind_vcf_op_262, NA.method =
                                          "mean")
pca_scheme1 <-
  dudi.pca(
    genind_vcf_op_262_scaled,
    cent = T,
    scale = F,
    scannf = F,
    nf = 10
  )
axis_all = pca_scheme1$eig * 100 / sum(pca_scheme1$eig)
barplot(axis_all[1:10], main = "PCA eigenvalues")
fviz_pca_ind(pca_scheme1)
ind = rownames(pca_scheme1)
#add site names to pca table
site <- rep(NA, length(ind))
site[grep("m", ind)] <- "Mornantais"
site[grep("p", ind)] <- "Pélussinois"
site[grep("a", ind)] <- "Annonéen"
pca_scheme1 = as.data.frame(pca_scheme1$li)
pca_scheme1 <- as.tibble(data.frame(pca_scheme1, site))

#plot the pca
ggplot(pca_scheme1, aes(x = Axis1, y = Axis2, col = site)) +
  geom_point(alpha = 0.6, size = 6)  +
  scale_color_manual(values = c("tomato", "skyblue", "orange")) +
  scale_shape_manual(values = c(1:10)) +
  labs(x = "PC 1 (5.3%)", y =  "PC 2 (1.3%)", color = "Sites") +
  geom_hline(yintercept = 0, colour = "grey") +
  geom_vline(xintercept = 0, colour = "grey") +
  scale_x_reverse() 

