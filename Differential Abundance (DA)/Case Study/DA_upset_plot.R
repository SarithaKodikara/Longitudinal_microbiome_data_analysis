library(ComplexHeatmap)
library(UpSetR)
timeRes<-readxl::read_xlsx("Results/Time_All.xlsx")

lt = list(NBMM = timeRes[timeRes$NBMM=="s",]$OTU,
          NBMM_AR = timeRes[timeRes$NBMM_AR=="s",]$OTU,
          ZINBMM = timeRes[timeRes$ZINBMM=="s",]$OTU,
          ZINBMM_AR = timeRes[timeRes$ZINBMM_AR=="s",]$OTU,
          ZIGMM_count = timeRes[timeRes$ZIGMM_count=="s",]$OTU,
          ZIGMM_count_AR = timeRes[timeRes$ZIGMM_count_AR=="s",]$OTU,
          ZIBR = timeRes[timeRes$ZIBR=="s",]$OTU,
          ZIGMMA_RA = timeRes[timeRes$ZIGMMA_RA=="s",]$OTU,
          ZIGMMA_RA_AR = timeRes[timeRes$ZIGMMA_RA_AR=="s",]$OTU,
          SplinectomeR = timeRes[timeRes$SplinectomeR=="s",]$OTU)
names(lt) = c("NBMM", "NBMM*", "ZINBMM",
                     "ZINBMM*", "ZIGMM", "ZIGMM*",
                     "ZIBR", "ZIGMM", "ZIGMM*", "SplinectomeR")

m = make_comb_mat(lt)

png("Results/Plot1.png", width = 8, height = 8, units = 'in', res = 300)
set.seed(9810)
ht = draw(UpSet(m,  column_title="Time Effect", 
                bg_col = c("#F0F0FF", "#FFF0F0"), bg_pt_col = "#CCCCFF",
                top_annotation = HeatmapAnnotation(
                  "Number of\nCommon Taxa" = anno_barplot(comb_size(m), 
                                                          border = FALSE, 
                                                          gp = gpar(col = comb_degree(m), fill=comb_degree(m)), 
                                                          height = unit(5, "cm")
                  ), 
                  annotation_name_side = "left", 
                  annotation_name_rot = 0),
                right_annotation = rowAnnotation(
                  "Significant\nTaxa set size" = anno_barplot(set_size(m), 
                                                         border = FALSE, 
                                                         gp = gpar(fill = "grey", col="grey"), 
                                                         width = unit(3, "cm")
                  ),
                  `Input Type` = c("Count", "Count", "Count","Count", "Count", "Count","RA", "RA", "RA","RA"))))


od = column_order(ht)
cs = comb_size(m)
decorate_annotation("Number of\nCommon Taxa", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = "bottom", gp = gpar(fontsize = 8))
})
dev.off()



groupRes<-readxl::read_xlsx("Results/Group_All.xlsx")
lt = list(NBMM = groupRes[groupRes$NBMM=="s",]$OTU,
          NBMM_AR = groupRes[groupRes$NBMM_AR=="s",]$OTU,
          ZINBMM = groupRes[groupRes$ZINBMM=="s",]$OTU,
          ZINBMM_AR = groupRes[groupRes$ZINBMM_AR=="s",]$OTU,
          ZIGMM_count = groupRes[groupRes$ZIGMM_count=="s",]$OTU,
          ZIGMM_count_AR = groupRes[groupRes$ZIGMM_count_AR=="s",]$OTU,
          ZIBR = groupRes[groupRes$ZIBR=="s",]$OTU,
          ZIGMMA_RA = groupRes[groupRes$ZIGMMA_RA=="s",]$OTU,
          ZIGMMA_RA_AR = groupRes[groupRes$ZIGMMA_RA_AR=="s",]$OTU,
          SplinectomeR = groupRes[groupRes$SplinectomeR=="s",]$OTU)
names(lt) = c("NBMM", "NBMM*", "ZINBMM",
              "ZINBMM*", "ZIGMM", "ZIGMM*",
              "ZIBR", "ZIGMM", "ZIGMM*", "SplinectomeR")

m = make_comb_mat(lt)


png("Results/Plot2.png", width = 8, height = 8, units = 'in', res = 300)
set.seed(9810)
ht = draw(UpSet(m,  column_title="Group Effect", 
                bg_col = c("#F0F0FF", "#FFF0F0"), bg_pt_col = "#CCCCFF",
                top_annotation = HeatmapAnnotation(
                  "Number of\nCommon Taxa" = anno_barplot(comb_size(m), 
                                                          border = FALSE, 
                                                          gp = gpar(col = comb_degree(m), fill=comb_degree(m)), 
                                                          height = unit(5, "cm")
                  ), 
                  annotation_name_side = "left", 
                  annotation_name_rot = 0),
                right_annotation = rowAnnotation(
                  "Significant\nTaxa set size" = anno_barplot(set_size(m), 
                                                         border = FALSE, 
                                                         gp = gpar(fill = "grey", col="grey"), 
                                                         width = unit(3, "cm")
                  ),
                  `Input Type` = c("Count", "Count", "Count","Count", "Count", "Count","RA", "RA", "RA","RA"))))


od = column_order(ht)
cs = comb_size(m)
decorate_annotation("Number of\nCommon Taxa", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = "bottom", gp = gpar(fontsize = 8))
})
dev.off()


gTimeRes<-readxl::read_xlsx("Results/Group_Time_All.xlsx")
lt = list(NBMM = gTimeRes[gTimeRes$NBMM=="s",]$OTU,
          NBMM_AR = gTimeRes[gTimeRes$NBMM_AR=="s",]$OTU,
          ZINBMM = gTimeRes[gTimeRes$ZINBMM=="s",]$OTU,
          ZINBMM_AR = gTimeRes[gTimeRes$ZINBMM_AR=="s",]$OTU,
          ZIGMM_count = gTimeRes[gTimeRes$ZIGMM_count=="s",]$OTU,
          ZIGMM_count_AR = gTimeRes[gTimeRes$ZIGMM_count_AR=="s",]$OTU,
          ZIBR = gTimeRes[gTimeRes$ZIBR=="s",]$OTU,
          ZIGMMA_RA = gTimeRes[gTimeRes$ZIGMMA_RA=="s",]$OTU,
          ZIGMMA_RA_AR = gTimeRes[gTimeRes$ZIGMMA_RA_AR=="s",]$OTU)
names(lt) = c("NBMM", "NBMM*", "ZINBMM",
              "ZINBMM*", "ZIGMM", "ZIGMM*",
              "ZIBR", "ZIGMM", "ZIGMM*")

m = make_comb_mat(lt)


png("Results/Plot3.png", width = 8, height = 8, units = 'in', res = 300)
set.seed(9810)
ht = draw(UpSet(m,  column_title="Time and Group Interaction Effect", 
                bg_col = c("#F0F0FF", "#FFF0F0"), bg_pt_col = "#CCCCFF",
                top_annotation = HeatmapAnnotation(
                  "Number of\nCommon Taxa" = anno_barplot(comb_size(m), 
                                                          border = FALSE, 
                                                          gp = gpar(col = comb_degree(m), fill=comb_degree(m)), 
                                                          height = unit(5, "cm")
                  ), 
                  annotation_name_side = "left", 
                  annotation_name_rot = 0),
                right_annotation = rowAnnotation(
                  "Significant\nTaxa set size" = anno_barplot(set_size(m), 
                                                         border = FALSE, 
                                                         gp = gpar(fill = "grey", col="grey"), 
                                                         width = unit(3, "cm")
                  ),
                  `Input Type` = c("Count", "Count", "Count","Count", "Count", "Count","RA", "RA", "RA"))))


od = column_order(ht)
cs = comb_size(m)
decorate_annotation("Number of\nCommon Taxa", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = "bottom", gp = gpar(fontsize = 8))
})
dev.off()



