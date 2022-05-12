


#confusion table from timeOmics package
#https://github.com/abodein/timeOmics_frontiers/blob/e5a371bf84b0ec718d87fbb6888581e4a2a71c5c/Examples/Simulation.Rmd#L237

get_conf_table_bruteforce <- function(res.cluster, all.permutation.table, clusterLevels){
  CORRESPONDANCE <- NULL
  CONF.TABLE <- NULL
  ACC <- NULL
  res.cluster.tmp <- res.cluster %>% dplyr::select(cluster, first_cluster) %>%
    mutate(cluster = as.character(cluster) %>% as.numeric() %>% factor(levels = clusterLevels)) 
  for(i in 1:ncol(all.permutation.table)){
    correspondance <- list(new = c(all.permutation.table[,i]), 
                           old = c("c1", "c2", "c3", "c4")) %>% as.data.frame()
    conf <- suppressWarnings(left_join(res.cluster.tmp, correspondance, by = c("first_cluster"="old")))%>%
      dplyr::select(cluster, new) %>%
      mutate(new = factor(new, levels =  clusterLevels))
    acc <- accuracy(table(conf))
    #print(acc)
    if(is.null(ACC)){
      # first time
      ACC <- acc
      CONF.TABLE <- conf
      CORRESPONDANCE <- correspondance
    } 
    if (acc > ACC) {
      ACC <- acc
      CONF.TABLE <- conf
      CORRESPONDANCE <- correspondance
    } # else : pass
  }
  return(list(acc = ACC, conf.tab = CONF.TABLE))
}

#Clustering accuracy
accuracy <- function(conf.table){
  return(sum(diag(conf.table) / sum(conf.table)))
}

cluster_Methods_acc<-function(lmm.data){
  
  accuracy<-c()
  
  #PCA
  pca.res <- mixOmics::pca(lmm.data, center = FALSE, scale =FALSE)
  pca.res.cluster <- getCluster(pca.res) %>%
    dplyr::select(molecule, cluster) %>%
    mutate(first_cluster = molecule %>% str_split("\\.") %>% map_chr(~.x[1]))
  clusterL<-c(-2,-1,1,2)
  all.permutation.pca <- Deducer::perm(clusterL, duplicates = F) %>% t %>% 
    as.data.frame()
  conf.tab.pca <- get_conf_table_bruteforce(pca.res.cluster, all.permutation.pca, clusterL)
  ACC.pca <- conf.tab.pca$acc
  
  #DTW
  clust.pam <- tsclust(t(lmm.data),  k=4L, distance="dtw", centroid = "pam")
  #plot(clust.pam)
  dtw.cluster<-data.frame(t(rbind(lmm.data[0,], cluster = clust.pam@cluster)))%>%
    rownames_to_column("molecule") %>%
    mutate(first_cluster = molecule %>% str_split("\\.") %>% map_chr(~.x[1]))
  clusterL<-c(1,2,3,4)
  all.permutation <- Deducer::perm(clusterL, duplicates = F) %>% t %>% 
    as.data.frame()
  conf.tab.dtw <- get_conf_table_bruteforce(dtw.cluster, all.permutation, clusterL)
  ACC.dtw <- conf.tab.dtw$acc
  
  
  ##Hierarchical agglomerative clustering
  temporal_dmat<-dist(t(lmm.data))
  hc_full_cluster<-hclust(temporal_dmat)
  clust_res<-cutree(hc_full_cluster,4)
  agg.cluster<-data.frame(t(rbind(lmm.data[0,], cluster = clust_res)))%>%
    rownames_to_column("molecule") %>%
    mutate(first_cluster = molecule %>% str_split("\\.") %>% map_chr(~.x[1]))
  conf.tab.agg <- get_conf_table_bruteforce(agg.cluster, all.permutation, clusterL)
  ACC.agg <- conf.tab.agg$acc
  
  ##k-medoids clustering
  clust_Kmedoids<-cluster::pam(temporal_dmat,k=4)
  kmed.cluster<-data.frame(t(rbind(lmm.data[0,], cluster = clust_Kmedoids$clustering)))%>%
    rownames_to_column("molecule") %>%
    mutate(first_cluster = molecule %>% str_split("\\.") %>% map_chr(~.x[1]))
  conf.tab.kmed <- get_conf_table_bruteforce(kmed.cluster, all.permutation, clusterL)
  ACC.kmed <- conf.tab.kmed$acc
  
  accuracy<-c(PCA=ACC.pca,
              DTW= ACC.dtw, Agglomerative=ACC.agg, Kmedoids=ACC.kmed)
  return(accuracy)
}

plotFunction<-function(clusData, plotTitle, scale=FALSE, center=FALSE){
  # gather data
  if(scale & center){
    data.gathered <- clusData%>% scale(center = TRUE, scale = TRUE) %>% as.data.frame() %>% 
      rownames_to_column("time") %>%
      mutate(time = as.numeric(time)) %>%
      pivot_longer(names_to="feature", values_to = 'value', -time)%>%
      mutate(Cluster= substr(feature, start = 1, stop = 2)) 
  }
  else if (scale){
    data.gathered <- clusData%>% scale(center = FALSE, scale = TRUE) %>% as.data.frame() %>% 
      rownames_to_column("time") %>%
      mutate(time = as.numeric(time)) %>%
      pivot_longer(names_to="feature", values_to = 'value', -time)%>%
      mutate(Cluster= substr(feature, start = 1, stop = 2)) 
  }
  else if (center){
    data.gathered <- clusData%>% scale(center = TRUE, scale = FALSE) %>% as.data.frame() %>% 
      rownames_to_column("time") %>%
      mutate(time = as.numeric(time)) %>%
      pivot_longer(names_to="feature", values_to = 'value', -time)%>%
      mutate(Cluster= substr(feature, start = 1, stop = 2)) 
  }
  else{
    data.gathered <- clusData%>% as.data.frame() %>% 
      rownames_to_column("time") %>%
      mutate(time = as.numeric(time)) %>%
      pivot_longer(names_to="feature", values_to = 'value', -time)%>%
      mutate(Cluster= substr(feature, start = 1, stop = 2)) 
  }
  
  # plot profiles
  ggplot(data.gathered, 
         aes(x = time, y = value, group = feature,colour=Cluster)) + 
    geom_line() +
    theme_bw() + theme(legend.position="none")+ 
    scale_color_brewer(palette="Set2")+
    ggtitle(plotTitle) + ylab("Feature expression") +
    xlab("Time")+ 
    facet_wrap(~ Cluster)
} 

clusBoxplot<-function(clusRes, plotTitle){
  cl_long <- gather(clusRes, Method, Accuracy)
  
  p<-ggplot(cl_long, aes(x=Method, y=Accuracy, fill=Method)) + 
    geom_boxplot(na.rm = TRUE)+coord_flip()+
    ggtitle(plotTitle)+theme(axis.text=element_text(size=12),
                             axis.title=element_text(size=14,face="bold"),
                             legend.position = "none")+scale_y_continuous(limits = c(0.4, 1))
}  

