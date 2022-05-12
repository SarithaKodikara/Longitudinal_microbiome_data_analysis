library(ggplot2)
library(timeOmics)
library(tidyverse)
library(lmms)
library(magrittr)
library(reshape2)
library(cluster)
library(dtwclust)
library(tibble)
library(RColorBrewer)
library(patchwork)
set.seed(1234)
metafull_df <- read.csv("Data/Mu_etal_mouseVRE_MAPPING_FILE.txt",sep = "\t")

OTU_tsv <- read.csv("Data/Mu_etal_mouseVRE_FEATURE_TABLE_case_corrected.tsv",sep = "\t",skip = 1)
OTU_tsv <-cbind(OTU.ID=paste0("OTU_", seq(1,3574,1)),OTU_tsv)
colnames(OTU_tsv)[2]<-"Feature.ID"

taxa_clasification<-read.csv( "Data/taxonomy.tsv",sep = "\t")

taxa_details<-left_join(OTU_tsv[,1:2],taxa_clasification,
                        by='Feature.ID')

OTU_DF<-OTU_tsv[,-2]
meta_DF<-data.frame(Sample.ID=paste0("X",metafull_df$X.SampleID),
                    Time=metafull_df$day_of_experiment,
                    Subject.ID=metafull_df$host_subject_id)


metaData_0_5<-meta_DF %>%filter(.,Time<6)%>%
  filter(., Subject.ID!="not applicable")%>%
  mutate(.,New_Time=rep(1:4, each= 9))%>%
  mutate(.,Group=rep(0,36))%>%
  mutate(.,Subject_ID_new=paste0(Subject.ID,"_0"))
OTU_0_5<-OTU_DF[,c("OTU.ID",metaData_0_5$Sample.ID)]
librarySize_0_5<-data.frame(Sample.ID=colnames(OTU_0_5)[-1],
                            library.Size= colSums(OTU_0_5[,-1]))
meta_0_5<-left_join(metaData_0_5,librarySize_0_5,by="Sample.ID")


metaData_9_14<-meta_DF %>%filter(.,Time>8)%>%
  filter(., Subject.ID!="not applicable")%>%
  mutate(.,New_Time=rep(1:4, each= 9))%>%
  mutate(.,Group=rep(1,36))%>%
  mutate(.,Subject_ID_new=paste0(Subject.ID,"_1"))
OTU_9_14<-OTU_DF[,c("OTU.ID",metaData_9_14$Sample.ID)]
librarySize_9_14<-data.frame(Sample.ID=colnames(OTU_9_14)[-1],
                             library.Size= colSums(OTU_9_14[,-1]))
meta_9_14<-left_join(metaData_9_14,librarySize_9_14,by="Sample.ID")


metaData_new<-rbind(meta_0_5,meta_9_14)
rownames(metaData_new)<-metaData_new$Sample.ID
metaData_new<-metaData_new[,-1]
  
OTU_new<-left_join(OTU_0_5,OTU_9_14, by="OTU.ID")
#removinf low abundant OTUs (RA<0.01%)
keep.otu = OTU_new[which(rowSums(OTU_new[,-1])*100/(sum(rowSums(OTU_new[,-1]))) > 0.01),"OTU.ID"]
OTU_filt_count<-data.frame(t(OTU_new[OTU_new$OTU.ID %in%keep.otu,-1]))
colnames(OTU_filt_count)<-OTU_new[OTU_new$OTU.ID %in%keep.otu,1]

OTU_filt_RA<-OTU_filt_count/metaData_new$library.Size

taxanomy_new<-taxa_details[taxa_details$OTU.ID%in%colnames(OTU_filt_count),]

#lmms data
data_g1_g2<-left_join(rownames_to_column(metaData_new),
                      rownames_to_column(OTU_filt_RA),
                      by="rowname")
#Group0
data_g0<-data_g1_g2%>%filter(Group==0)

time<-data_g0$New_Time
con<-data_g0[,8:length(data_g0)]
row.names(con)<-data_g0$sample_name


lmms.output <- lmms::lmmSpline(data = con, time = time,
                               sampleID = rownames(con), deri = FALSE,
                               basis = "p-spline", numCores = 4, timePredict = 1:4,
                               keepModels = TRUE)
modelled.data <- t(slot(lmms.output, 'predSpline'))
modelled.data_f<- modelled.data[,colSums(modelled.data)!=0]

# centering and scaling
#scaleCentered.lmmData<- scale(modelled.data_f)

longlmms<-melt(modelled.data_f)
colnames(longlmms)<-c("Time","molecule", "LMMS value")

#PCA
pca.res <- mixOmics::pca(modelled.data_f, center = FALSE, scale =FALSE)
pca.res.cluster <- getCluster(pca.res) %>%
  dplyr::select(molecule, cluster) 
plotLong(pca.res, scale = FALSE, center = FALSE, 
         title = "PCA longitudinal clustering")
pca.res_g0 <- merge(longlmms, pca.res.cluster, by = 'molecule')
p <- ggplot(pca.res_g0, aes(y=`LMMS value`, x=Time, 
                            group=molecule,color=cluster)) + 
  scale_color_brewer(palette = "Set2")+ geom_line(size=0.6)
p0_pca<-p + facet_wrap(vars(cluster),drop = FALSE)+theme(legend.position = 'none')+ 
  ggtitle('PCA (Naive Phase)')
write.csv(pca.res.cluster , file="Results/pca_naive.csv")

#DTW
clust.pam <- tsclust(t(modelled.data_f),  k=4L, distance="dtw", centroid = "pam")
plot(clust.pam)
dtw.cluster<-data.frame(t(rbind(modelled.data_f[0,], cluster = clust.pam@cluster)))%>%
  rownames_to_column("molecule")
dtw.res_g0 <- merge(longlmms, dtw.cluster, by = 'molecule')
dtw.res_g0$cluster<-factor(dtw.res_g0$cluster)
p <- ggplot(dtw.res_g0, aes(y=`LMMS value`, x=Time, 
                            group=molecule,color=cluster)) + 
  scale_color_brewer(palette = "Set2")+ geom_line(size=0.6)
p0_dtw<-p + facet_wrap(vars(cluster))+theme(legend.position = 'none') +
ggtitle('DTW (Naive Phase)')
write.csv(dtw.cluster , file="Results/dtw_naive.csv")

##Hierarchical agglomerative clustering
temporal_dmat<-dist(t(modelled.data_f))
hc_full_cluster<-hclust(temporal_dmat)
clust_res<-cutree(hc_full_cluster,4)
agg.cluster<-data.frame(t(rbind(modelled.data_f[0,], cluster = clust_res)))%>%
  rownames_to_column("molecule")
agg.res_g0 <- merge(longlmms, agg.cluster, by = 'molecule')
agg.res_g0$cluster<-factor(agg.res_g0$cluster)
p <- ggplot(agg.res_g0, aes(y=`LMMS value`, x=Time, 
                            group=molecule,color=cluster)) + 
  scale_color_brewer(palette = "Set2")+ geom_line(size=0.6)
p0_agg<-p + facet_wrap(vars(cluster))+theme(legend.position = 'none') +
ggtitle('Agglomerative (Naive Phase)')
write.csv(agg.cluster , file="Results/agg_naive.csv")

##k-medoids clustering
clust_Kmedoids<-cluster::pam(temporal_dmat,k=4)
kmed.cluster<-data.frame(t(rbind(modelled.data_f[0,], cluster = clust_Kmedoids$clustering)))%>%
  rownames_to_column("molecule")
kmed.res_g0 <- merge(longlmms, kmed.cluster, by = 'molecule')
kmed.res_g0$cluster<-factor(kmed.res_g0$cluster)
p <- ggplot(kmed.res_g0, aes(y=`LMMS value`, x=Time, 
                            group=molecule,color=cluster)) + 
  scale_color_brewer(palette = "Set2")+ geom_line(size=0.6)
p0_kmed<-p + facet_wrap(vars(cluster))+theme(legend.position = 'none') +
ggtitle('Kmedoids (Naive Phase)')
write.csv(kmed.cluster , file="Results/kmed_naive.csv")


#Group1
data_g1<-data_g1_g2%>%filter(Group==1)

time<-data_g1$New_Time
con<-data_g1[,8:length(data_g1)]
row.names(con)<-data_g1$sample_name


lmms.output <- lmms::lmmSpline(data = con, time = time,
                               sampleID = rownames(con), deri = FALSE,
                               basis = "p-spline", numCores = 4, timePredict = 1:4,
                               keepModels = TRUE)
modelled.data <- t(slot(lmms.output, 'predSpline'))
modelled.data_f<- modelled.data[,colSums(modelled.data)!=0]
  
# Sentering and scaling
#scaleCentered.lmmData<- scale(modelled.data_f)

longlmms<-melt(modelled.data_f)
colnames(longlmms)<-c("Time","molecule", "LMMS value")

#PCA
pca.res <- mixOmics::pca(modelled.data_f, center = FALSE, scale =FALSE)
pca.res.cluster1 <- getCluster(pca.res) %>%
  dplyr::select(molecule, cluster) 
plotLong(pca.res, scale = FALSE, center = FALSE, 
         title = "PCA longitudinal clustering")
pca.res_g1 <- merge(longlmms, pca.res.cluster1, by = 'molecule')
p <- ggplot(pca.res_g1, aes(y=`LMMS value`, x=Time, 
                            group=molecule,color=cluster)) + 
  scale_color_brewer(palette = "Set2")+ geom_line(size=0.6)
p1_pca<-p + facet_wrap(vars(cluster),drop = FALSE)+theme(legend.position = 'none') +
ggtitle('PCA (VRE Phase)')
write.csv(pca.res.cluster1 , file="Results/pce_vre.csv")

#DTW
clust.pam <- tsclust(t(modelled.data_f),  k=4L, distance="dtw", centroid = "pam")
plot(clust.pam)
dtw.cluster1<-data.frame(t(rbind(modelled.data_f[0,], cluster = clust.pam@cluster)))%>%
  rownames_to_column("molecule")
dtw.res_g1 <- merge(longlmms, dtw.cluster1, by = 'molecule')
dtw.res_g1$cluster<-factor(dtw.res_g1$cluster)
p <- ggplot(dtw.res_g1, aes(y=`LMMS value`, x=Time, 
                            group=molecule,color=cluster)) + 
  scale_color_brewer(palette = "Set2")+ geom_line(size=0.6)
p1_dtw<-p + facet_wrap(vars(cluster))+theme(legend.position = 'none')+
  ggtitle('DTW (VRE Phase)')
write.csv(dtw.cluster1 , file="Results/dtw_vre.csv")

##Hierarchical agglomerative clustering
temporal_dmat<-dist(t(modelled.data_f))
hc_full_cluster<-hclust(temporal_dmat)
clust_res<-cutree(hc_full_cluster,4)
agg.cluster1<-data.frame(t(rbind(modelled.data_f[0,], cluster = clust_res)))%>%
  rownames_to_column("molecule")
agg.res_g1 <- merge(longlmms, agg.cluster1, by = 'molecule')
agg.res_g1$cluster<-factor(agg.res_g1$cluster)
p <- ggplot(agg.res_g1, aes(y=`LMMS value`, x=Time, 
                            group=molecule,color=cluster)) + 
  scale_color_brewer(palette = "Set2")+ geom_line(size=0.6)
p1_agg<-p + facet_wrap(vars(cluster))+theme(legend.position = 'none')+
  ggtitle('Agglomerative (VRE Phase)')
write.csv(agg.cluster1 , file="Results/agg_vre.csv")

##k-medoids clustering
clust_Kmedoids<-cluster::pam(temporal_dmat,k=4)
kmed.cluster1<-data.frame(t(rbind(modelled.data_f[0,], cluster = clust_Kmedoids$clustering)))%>%
  rownames_to_column("molecule")
kmed.res_g1 <- merge(longlmms, kmed.cluster1, by = 'molecule')
kmed.res_g1$cluster<-factor(kmed.res_g1$cluster)
p <- ggplot(kmed.res_g1, aes(y=`LMMS value`, x=Time, 
                            group=molecule,color=cluster)) + 
  scale_color_brewer(palette = "Set2")+ geom_line(size=0.6)
p1_kmed<-p + facet_wrap(vars(cluster))+theme(legend.position = 'none')+
  ggtitle('Kmedoids (VRE Phase)')
write.csv(kmed.cluster1 , file="Results/kmed_vre.csv")


png("Results/Plot4.png", width = 10, height = 10, units = 'in', res = 300)
(p0_pca+p1_pca)/(p0_dtw+p1_dtw)/(p0_agg+p1_agg)/(p0_kmed+p1_kmed)+ plot_annotation(tag_levels = 'A')

dev.off()



