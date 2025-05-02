#This script contains scripts to compute the FST between populations in oenanthe species. the same was calculated for the other two species. It also contains code to cpmpute genetic distance between individuals in all populations, called Nei genetic distance. 

library(StAMPP)
library(adegenet)
library(RColorBrewer)
library(gplots)


#load the vcf file 
vcf_op_262 <-
  read.vcfR(
    "/home/alle/Documents/Pilat_project/process_all/processing_op/corrected_op/denovo_op_m3_M5_n5/op_GQ20_minDP5_minmeanDP10_miss080_mac3_ind262.recode.vcf"
  )
#convert vcf to genlight
genlight_vcf_op_262 <- vcfR2genlight(vcf_op_262)

#add the populations per site, so 30 in total 
site <- sapply(strsplit(genlight_vcf_op_262@ind.names, '_'),
               function(genlight_vcf_op_262){genlight_vcf_op_262[1]})
pop <- sapply(strsplit(genlight_vcf_op_262@ind.names, '_'),
              function(genlight_vcf_op_262){genlight_vcf_op_262[3]})
popsite <- paste(site,pop,sep="_")
popsite = as.factor(popsite)
genlight_vcf_op_262$pop = popsite

#convert genlight to stammp object
stamp_op_262 = stamppConvert(genlight_vcf_op_262, type = "genlight")

##population fst
stamp_op_262_FST = stamppFst(stamp_op_262, nboots = 100, percent = 95, nclusters = 8)
stamp_op_262_FST = as.matrix(stamp_op_262_FST$Fsts)
diag(stamp_op_262_FST) <- 0
#copy lower tri to upper tri
stamp_op_262_FST[upper.tri(stamp_op_262_FST)]  <- t(stamp_op_262_FST)[upper.tri(stamp_op_262_FST)]
stamp_op_262_FST = as.data.frame(stamp_op_262_FST)
#order population names
stamp_op_262_FST=stamp_op_262_FST[order(row.names(stamp_op_262_FST)), order(colnames(stamp_op_262_FST))]
stamp_op_262_FST = as.matrix(stamp_op_262_FST)
#color the popsites accprding to the 3 sites
my_group <- as.numeric(as.factor(substr(rownames(stamp_op_262_FST), 1 , 1)))
colSide <- brewer.pal(9, "Set1")[my_group]
colSide <- c("tomato", "skyblue", "orange")[my_group]
#do heatmap
mtscaled <- scale(stamp_op_262_FST)
my_palette <- colorRampPalette(c("tomato", "skyblue", "orange"))(n = 1000)
heatmap.2(stamp_op_262_FST, Rowv=F, Colv="Rowv", scale='none', symm =T, trace="none",  RowSideColors=colSide)
heatmap.2(stamp_op_262_FST,  RowSideColors=colSide, trace="none")

#individual genetic distance
stamp_op_262_Nei=stamppNeisD(stamp_op_262, pop = FALSE)
stamp_op_262_Nei=as.matrix(stamp_op_262_Nei)
stamp_op_262_Nei=as.data.frame(stamp_op_262_Nei)
colnames(stamp_op_262_Nei)=row.names(stamp_op_262_Nei)
stamp_op_262_Nei=stamp_op_262_Nei[order(row.names(stamp_op_262_Nei)), order(colnames(stamp_op_262_Nei))]
diag(stamp_op_262_Nei) <- 0
#copy lower tri to upper tri
stamp_op_262_Nei[upper.tri(stamp_op_262_Nei)]  <- t(stamp_op_262_Nei)[upper.tri(stamp_op_262_Nei)]
#color the individuals accroding to the 3 sites 
my_group <- as.numeric(as.factor(substr(rownames(stamp_op_262_Nei), 1 , 1)))
colSide <- c("tomato", "skyblue", "orange")[my_group]
heatmap(stamp_op_262_Nei,RowSideColors=colSide)
stamp_op_262_Nei = as.matrix(stamp_op_262_Nei)
heatmap.2(stamp_op_262_Nei, Rowv=F, Colv="Rowv", scale='none', symm =T, trace="none",  RowSideColors=colSide)
heatmap.2(stamp_op_262_Nei,  RowSideColors=colSide, trace="none")
