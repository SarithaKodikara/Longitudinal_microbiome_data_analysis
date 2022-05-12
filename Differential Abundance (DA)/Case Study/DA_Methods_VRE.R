library(ggplot2)
library(magrittr)
library(dplyr)
library(tscount)
library(ZIBR)
library(NBZIMM)
library(splinectomeR)
library(epiR)
source("Utility_functions_DA.R")

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

OTU_new<-left_join(OTU_0_5,OTU_9_14, by="OTU.ID")
#removinf low abundant OTUs (RA<0.01%)
keep.otu = OTU_new[which(rowSums(OTU_new[,-1])*100/(sum(rowSums(OTU_new[,-1]))) > 0.01),"OTU.ID"]
OTU_filt_count<-data.frame(t(OTU_new[OTU_new$OTU.ID %in%keep.otu,-1]))
colnames(OTU_filt_count)<-OTU_new[OTU_new$OTU.ID %in%keep.otu,1]

OTU_filt_RA<-OTU_filt_count/metaData_new$library.Size

dominantTaxa<-names(sort(colSums(OTU_filt_RA), decreasing = TRUE)[1:10])

taxanomy_new<-taxa_details[taxa_details$OTU.ID%in%colnames(OTU_filt_count),]
write.csv(taxanomy_new, file="Results/taxanomy_filtered.csv")


####ZIGMM----

###RA- no AR ####
zigmm<-Gaussian_NB_methods(OTU_filt_RA,metaData_new, 
                           correlation_AR=FALSE, countData=FALSE, method="zigmm", nOTU=193)
write.csv(zigmm, file="Results/ZIGMM.csv")

###RA- with AR ####
zigmm_AR<-Gaussian_NB_methods(OTU_filt_RA,metaData_new, 
                              correlation_AR=TRUE, countData=FALSE, method="zigmm", nOTU=193)
write.csv(zigmm_AR, file="Results/ZIGMM_AR.csv")

###Count- no AR ####
zigmmCo<-Gaussian_NB_methods(OTU_filt_count,metaData_new, 
                             correlation_AR=FALSE, countData=TRUE, method="zigmm", nOTU=193)
write.csv(zigmmCo, file="Results/zigmmCo.csv")

###Count- with AR ####
zigmmCo_AR<-Gaussian_NB_methods(OTU_filt_count,metaData_new,
                                correlation_AR=TRUE, countData=TRUE, method="zigmm", nOTU=193)
write.csv(zigmmCo_AR, file="Results/zigmmCo_AR.csv")

####ZIBR----

zibr<-zibrRes(OTU_filt_RA,metaData_new,nOTU=193)
write.csv(zibr, file="Results/zibr.csv")

####NBMM----  

###Count- no AR ####
nbmm<-Gaussian_NB_methods(OTU_filt_count,metaData_new, 
                          correlation_AR=FALSE, countData=TRUE, method="nbmm", nOTU=193)
write.csv(nbmm, file="Results/nbmm.csv")

###Count- with AR ####
nbmm_AR<-Gaussian_NB_methods(OTU_filt_count,metaData_new, 
                             correlation_AR=TRUE, countData=TRUE, method="nbmm", nOTU=193)
write.csv(nbmm_AR, file="Results/nbmm_AR.csv")


####ZINBMM----

###Count- no AR ####
zinbmm<-Gaussian_NB_methods(OTU_filt_count,metaData_new, 
                            correlation_AR=FALSE, countData=TRUE, method="zinbmm", nOTU=193)
write.csv(zinbmm, file="Results/zinbmm.csv")

###Count- with AR ####
zinbmm_AR<-Gaussian_NB_methods(OTU_filt_count,metaData_new, 
                               correlation_AR=TRUE, countData=TRUE, method="zinbmm", nOTU=193)
write.csv(zinbmm_AR, file="Results/zinbmm_AR.csv")

####SplinectomeR----

SplinectomeR<-splineR(OTU_filt_RA,metaData_new, nOTU=193)
write.csv(SplinectomeR, file="Results/SplinectomeR.csv")