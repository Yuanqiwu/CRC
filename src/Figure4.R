# ##############################################################################
#
##  function analysis
#
# ##############################################################################

library(RColorBrewer)
library(pheatmap)

function_combion <- read.csv(file = "Figure4.csv", header = TRUE, 
                             row.names = 2,check.names = FALSE,stringsAsFactors = FALSE)
# color line seting
row <- c(1,1,1,1,-0.250074656,-0.250074656)
#
function_combion <- rbind(function_combion,row)
function_combion$surperclass <- factor(function_combion$surperclass)

annotation_r <- data.frame(surperclass = function_combion[,2])

rownames(annotation_r) <- rownames(function_combion)

ann_colors <- list('surperclass' = c( 'Amine and Polyamine Degradation' = "#550A46",'Amino Acid Biosynthesis/Degradation' = "#902044",
                                      'Aromatic Compound Degradation'= "#CB414B", 'Carbohydrate Biosynthesis/Degradation' = "#D16F7C",
                                      'Carboxylate Degradation' =  "#E56B46", 'Cell Structure Biosynthesis' ="#e59346",
                                      'Cofactor, Prosthetic Group, Electron Carrier, and Vitamin Biosynthesis'= "#F4A862", 'Electron Transfer' = "#F6DB86",
                                      'Fatty Acid and Lipid Biosynthesis/Degradation' = "#e9f686", 'Fermentation' = "#DFE899", 'Inorganic Nutrient Metabolism' = "#A4D5A0",
                                      'Nucleoside and Nucleotide Biosynthesis/Degradation' = "#62BA9F", 'Secondary Metabolite Biosynthesis/Degradation' = "#3681AD",
                                      'TCA cycle' = "#5A4C97",'Others' = "#9D9B9C", "1" = "#FFFFFF"))


g <- pheatmap(as.data.frame(function_combion[,5:6]),color = colorRampPalette(c("navy", "white", "firebrick3"))(50), 
              cluster_cols = FALSE,fontsize = 6,  cellwidth=50, cellheight = 7, border_color = NA,
              cluster_rows = TRUE, clustering_method = "complete", scale = 'none',
              annotation_row  = annotation_r,   treeheight_row = 120,annotation_colors = ann_colors,
              display_numbers= matrix(ifelse(function_combion[3:4] < 0.01, "**", ifelse(function_combion[3:4] < 0.05, "*","")), 
                                      ncol = 2),number_color="black",fontsize_number=6)
function_heatmap <- as.data.frame(function_combion[,c(2,3,5,6)])
write.csv(function_heatmap, 'function_heamap_plot.csv')
