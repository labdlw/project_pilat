
#To run the admixture analysis using snmf function, the vcf file was converted to structure file using the populations program from stacks.
#the structure file was then modified to replace all 0 (missing data) with -9 values. so that the function can read the file. 

library(LEA)
setwd("/home/alle/Documents/Pilat_project/process_all/processing_op/corrected_op/denovo_op_m3_M5_n5")
getwd()

struct2geno(input.file = "/home/alle/Documents/Pilat_project/process_all/processing_op/corrected_op/denovo_op_m3_M5_n5/final.txt", ploidy = 2, FORMAT = 2 ,extra.col = 2, extra.row = 1)



#to find the best K,  1 to 20 with 100 repetitions
#alpha = 1
obj.snmf_a1 = snmf("/media/alle/Ubuntu-Data1/snmf/op/final.txt.geno", K = 1:20, ploidy = 2, entropy = T,alpha = 1, project = "new", repetitions = 200, CPU = 4)

#alpha = 10
obj.snmf_a10 = snmf("/media/alle/Ubuntu-Data1/snmf/op/final.txt_a10.geno", K = 1:20, ploidy = 2, entropy = T,alpha = 10, project = "new", repetitions = 200, CPU = 4)

#alpha = 100
obj.snmf_a100 = snmf("/media/alle/Ubuntu-Data1/snmf/op/final.txt_a100.geno", K = 1:20, ploidy = 2, entropy = T,alpha = 100, project = "new", repetitions = 200, CPU = 4)

#alpha = 1000
obj.snmf_a1000 = snmf("/media/alle/Ubuntu-Data1/snmf/op/final.txt_a1000.geno", K = 1:20, ploidy = 2, entropy = T, alpha = 1000, project = "new", repetitions = 200, CPU = 4)



plot(obj.snmf_a1, col = "darkgreen", pch = 20, type = "o", main="α = 1",ylim=c(0.415, 0.482))
abline(h = seq(0.482,0.41, length.out = 11), lty = "dashed", col = "lightgray")

plot(obj.snmf_a10, col = "darkgreen", pch = 20,type = "o", main="α = 10",ylim=c(0.415, 0.482))
abline(h = seq(0.482,0.41, length.out = 11), lty = "dashed", col = "lightgray")
ce <- cross.entropy(obj.snmf_a10, K = 4)

plot(obj.snmf_a100, col = "darkgreen", pch = 20, type = "o", main="α = 100",ylim=c(0.415, 0.482))
abline(h = seq(0.482,0.41, length.out = 11), lty = "dashed", col = "lightgray")
ce <- cross.entropy(obj.snmf_a100, K = 3)

plot(obj.snmf_a1000, col = "darkgreen", pch = 20, type = "o", main="α = 1000",ylim=c(0.415, 0.482))
abline(h = seq(0.482,0.41, length.out = 11), lty = "dashed", col = "lightgray")



plot(obj.snmf_a10, col = "blue4")
ce <- cross.entropy(obj.snmf_a10, K = 4)

####after visualizing choose the K
best <- which.min(ce)
qmatrix = Q(obj.snmf_a10, K = 4, run = best)
popmap_op_262 = read.table("/home/alle/Documents/Pilat_project/process_all/processing_op/corrected_op/denovo_op_m3_M5_n5/popmap_op_262.txt", header = T)
barplot(t(qmatrix), col=RColorBrewer::brewer.pal(4,"Paired"), 
        border=NA, space=0, xlab="Individuals", 
        ylab="Admixture coefficients")

#Add population labels to the axis:
pops=popmap_op_262$X2
site <- sapply(strsplit(popmap_op_262$X1, '_'),
                            function(popmap_op_262){popmap_op_262[1]})
pop <- sapply(strsplit(popmap_op_262$X1, '_'),
              function(popmap_op_262){popmap_op_262[3]})
popsite <- paste(site,pop,sep="_")
popsite = as.factor(popsite)

my.colors = c("tomato", "skyblue", "darkblue", "orange")

q_df <- qmatrix %>% 
  as_tibble() %>% 
  # add the pops data for plotting
  mutate(individual = popmap_op_262$X1,
         site = site,
         population = popsite)

write.csv(q_df,  "/home/alle/Documents/Pilat_project/process_all/processing_op/snmf_plots/admix_k4_run86.csv")

library("readxl")
samples_coord = read_excel("/home/alle/Documents/Pilat_project/samples_all_coord.xlsx")
colnames(samples_coord)[1] <- "individual"
samples_coord_filtered <- subset(samples_coord, individual %in% q_df$individual)
samples_coord_filtered = samples_coord_filtered[,c(1,7,8)]


q_df_coord=merge(q_df, samples_coord_filtered, by="individual")
write.csv(q_df_coord,  "/home/alle/Documents/Pilat_project/process_all/processing_op/snmf_plots/admix_k4_run86_latitude.csv")

averages_df <- q_df_coord %>%
  group_by(population) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

write.csv(averages_df,  "/home/alle/Documents/Pilat_project/process_all/processing_op/snmf_plots/admix_k4_op_coord_AVERAGE.csv")




LEA::barchart(obj.snmf_a10, K = 4, run = best,
              border = NA, space = 0,
              col = my.colors,
              xlab = "Individuals",
              ylab = "Ancestry proportions",
              main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),
     labels = bp$order, las=1,
     cex.axis = .4)



admix_table=read.table("/home/alle/Documents/Pilat_project/process_all/processing_op/corrected_op/denovo_op_m3_M5_n5/final.txt.snmf/K4/run4/final.txt_r4.4.Q") 
admix_table$pop = popmap_op_262$X2
admix_table$ind = popmap_op_262$X1
write.csv(admix_table, "/home/alle/Documents/Pilat_project/process_all/processing_op/corrected_op/denovo_op_m3_M5_n5/admix_table262_op.csv")

#### K3 ####
admix_table3=read.table("/home/alle/Documents/Pilat_project/process_all/processing_op/corrected_op/denovo_op_m3_M5_n5/final.txt.snmf/K3/run4/final.txt_r4.3.Q") 
admix_table3$pop = popmap_op_262$X2
admix_table3$ind = popmap_op_262$X1
write.csv(admix_table3, "/home/alle/Documents/Pilat_project/process_all/processing_op/corrected_op/denovo_op_m3_M5_n5/admix_table262_opk3.csv")

table(popmap_op_262$X2)

qmatrix = as.data.frame(qmatrix)
pop =  popmap_op_262$X2
qmatrix_pop = as.tibble(data.frame(qmatrix, pop))
#ordered the Q values according to population name
qmatrix_pop_ordered <- q_df[order(q_df$population),]
barplot(t(qmatrix_pop_ordered[,c(1:4)]), col = c("orange","violet","lightgreen","tomato"), border = NA, space = 0,
        xlab = "Individuals", ylab = "Admixture coefficients")
admix_table4=read.table("/home/alle/Documents/Pilat_project/process_all/processing_op/corrected_op/denovo_op_m3_M5_n5/final.txt.snmf/K4/run4/final.txt_r4.4.Q") 
admix_table4$pop = popmap_op_262$X2
admix_table4$ind = popmap_op_262$X1
write.csv(admix_table4, "/home/alle/Documents/Pilat_project/process_all/processing_op/corrected_op/denovo_op_m3_M5_n5/admix_table262_op.csv")

####### add coordinates to the admix_table ###########################
library("readxl")
samples_coord = read_excel("/home/alle/Documents/Pilat_project/samples_all_coord.xlsx")
colnames(samples_coord)[1] <- "ind"
samples_coord_filtered <- subset(samples_coord, ind %in% admix_table$ind)
samples_coord_filtered = samples_coord_filtered[,c(1,7,8)]
### merge again coordinates iwth the admixture values
admix_table4_coord = merge(admix_table4, samples_coord_filtered, by = "ind")
admix_table3_coord = merge(admix_table3, samples_coord_filtered, by = "ind")

write.csv(admix_table3_coord, "/home/alle/Documents/Pilat_project/process_all/processing_op/corrected_op/admix_table3_coord.csv")
write.csv(admix_table4_coord, "/home/alle/Documents/Pilat_project/process_all/processing_op/corrected_op/admix_table4_coord.csv")
