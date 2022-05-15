library(ComplexHeatmap)
taxanomyTop<-readxl::read_xlsx("Results/Top_10.xlsx", sheet = "Sheet1")
timeTop<-readxl::read_xlsx("Results/Top_10.xlsx", sheet = "Time")
groupTop<-readxl::read_xlsx("Results/Top_10.xlsx", sheet = "Group")
interactionTop<-readxl::read_xlsx("Results/Top_10.xlsx", sheet = "Time Group")

library(circlize)
col_fun = colorRamp2(c( 0,0.05,1), c( "red", "purple","blue"))
col_fun(seq(-3, 3))

ht_opt(heatmap_column_names_gp = gpar(fontface = "italic"), 
       heatmap_column_title_gp = gpar(fontsize = 10),
       legend_border = "black",
       heatmap_border = TRUE,
       annotation_border = TRUE
)
#ht_opt(RESET = TRUE)

mat1<-as.matrix(timeTop[,-c(1:2)])
rownames(mat1) = timeTop$Family
colnames(mat1) = c("NBMM","NBMM*","FZINBMM","FZINBMM*",     
                  "ZIGMM","ZIGMM*","ZIBR","ZIGMMA","ZIGMMA*","SplinectomeR")

mat_with_na1 = mat1
p1<-Heatmap(mat_with_na1, name = " p value", na_col = "white", col=col_fun,
        row_order = rownames(mat1),
        column_order = colnames(mat1),
        column_title = "Time Effect", 
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
        row_title="Taxa Family",
        row_title_side = "left",
        row_title_rot =90,
        row_title_gp = gpar(fontsize = 12, fontface = "bold"),
        rect_gp = gpar(col = "white", lwd = 0.5),
        column_split = c("Count", "Count", "Count","Count", "Count", "Count","RA", "RA", "RA","RA"),
        bottom_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:3),
                                                            labels = c("Count Methods", "RA Methods"), 
                                                            labels_gp = gpar(col = "black", fontsize = 12))))

mat2<-as.matrix(groupTop[,-c(1,2)])
rownames(mat2) = groupTop$Family
colnames(mat2) = c("NBMM","NBMM*","FZINBMM","FZINBMM*",     
                   "ZIGMM","ZIGMM*","ZIBR","ZIGMMA","ZIGMMA*","SplinectomeR")

mat_with_na2 = mat2
p2<-Heatmap(mat_with_na2, name = " p value", na_col = "white", col=col_fun,
            row_order = rownames(mat2),
            column_order = colnames(mat2),
            column_title = "Group Effect", 
            column_title_side = "top",
            column_title_gp = gpar(fontsize = 12, fontface = "bold"),
            row_title="Taxa Family",
            row_title_side = "left",
            row_title_rot =90,
            row_title_gp = gpar(fontsize = 12, fontface = "bold"),
            rect_gp = gpar(col = "white", lwd = 0.5),
            column_split = c("Count", "Count", "Count","Count", "Count", "Count","RA", "RA", "RA","RA"),
            bottom_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:3),
                                                                   labels = c("Count Methods", "RA Methods"), 
                                                                   labels_gp = gpar(col = "black", fontsize = 12))))

mat3<-as.matrix(interactionTop[,-c(1,2)])
rownames(mat3) = interactionTop$Family
colnames(mat3) = c("NBMM","NBMM*","FZINBMM","FZINBMM*",     
                   "ZIGMM","ZIGMM*","ZIBR","ZIGMMA","ZIGMMA*")

mat_with_na3 = mat3
p3<-Heatmap(mat_with_na3, name = " p value", na_col = "white", col=col_fun,
            row_order = rownames(mat3),
            column_order = colnames(mat3),
            column_title = "Time*Group Effect", 
            column_title_side = "top",
            column_title_gp = gpar(fontsize = 12, fontface = "bold"),
            row_title="Taxa Family",
            row_title_side = "left",
            row_title_rot =90,
            row_title_gp = gpar(fontsize = 12, fontface = "bold"),
            rect_gp = gpar(col = "white", lwd = 0.5),
            column_split = c("Count", "Count", "Count","Count", "Count", "Count","RA", "RA", "RA"),
            bottom_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:3),
                                                                   labels = c("Count Methods", "RA Methods"), 
                                                                   labels_gp = gpar(col = "black", fontsize = 12))))

ht_list =p1+p2+p3
png("Results/VRE_DA_top.png", width =14, height = 8, units = 'in', res = 300)
draw(ht_list, ht_gap = unit(1, "cm"))
dev.off()
