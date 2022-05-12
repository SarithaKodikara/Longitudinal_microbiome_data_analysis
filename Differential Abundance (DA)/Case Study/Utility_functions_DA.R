significance <- function(x) if_else(x<0.05, "s", "ns")

Gaussian_NB_methods<-function(df,metadf,correlation_AR=FALSE, countData=FALSE, method, nOTU=300){
  p_time_m1<-c()
  p_group_m1<-c()
  p_group_time_m1<-c()
  df_new<-cbind(metadf,df)
  
  options(warn=-1)#Suppress warnings 
  
  for(i in 1:nOTU){
    
    f<-if(countData) formula(paste(colnames(df[i]), "New_Time+Group+New_Time*Group+ offset(log(library.Size))", sep = "~")) else
      formula(paste(colnames(df[i]), "New_Time+Group+New_Time*Group", sep = "~"))
    args <- list(fixed = f, random = ~ 1 | Subject_ID_new, data = df_new)
    
    if (correlation_AR)
      args$correlation <- corAR1()
    
    if(method=="nbmm"){
      tryCatch( { fit <-  do.call('glmm.nb', args = args)
      fit.summary <- summary(fit)
      p_time_m1[i]=fit.summary$tTable[2,5]
      p_group_m1[i]=fit.summary$tTable[3,5]
      p_group_time_m1[i]=fit.summary$tTable[4,5]}, error = function(e) {
       p_time_m1[i]=NA
       p_group_m1[i]=NA
       p_group_time_m1[i]=NA})
    }
    else{
      if(sum(df[i]==0)>0){
        if(method=="zigmm"){
          fit <- do.call('lme.zig', args = args)
          fit.summary <- summary(fit)
          p_time_m1[i]=fit.summary$tTable[2,5]
          p_group_m1[i]=fit.summary$tTable[3,5]
          p_group_time_m1[i]=fit.summary$tTable[4,5]
        }
        else if(method=="zinbmm"){
          tryCatch( { fit <-  do.call('glmm.zinb', args = args)
          fit.summary <- summary(fit)
          p_time_m1[i]=fit.summary$tTable[2,5]
          p_group_m1[i]=fit.summary$tTable[3,5]
          p_group_time_m1[i]=fit.summary$tTable[4,5]}, error = function(e) {
            p_time_m1[i]=NA
            p_group_m1[i]=NA
            p_group_time_m1[i]=NA})
        }
      }
      else{
        p_time_m1[i]=NA
        p_group_m1[i]=NA
        p_group_time_m1[i]=NA
      }
    }
  } 
  p_value_DF<-data.frame(p_time=p_time_m1,
                p_group=p_group_m1,
                p_group_time=p_group_time_m1)
  p_value_DF<-cbind(OTU= colnames(df),p_value_DF,
                    mutate_each(p_value_DF,funs(significance)))
  return(p_value_DF)
}

zibrRes<-function(df_ra,metadf, nOTU=300){
  
  p_time_m1<-c()
  p_group_m1<-c()
  p_group_time_m1<-c()
  
  
  df_cov = data.frame(Time = metadf$New_Time, Group=metadf$Group, 
                      Time_Group=metadf$New_Time*metadf$Group)
  options(warn=-1)#Suppress warnings 

  zibr.res <-lapply(1:nOTU, function(i){
    zibr.fit <- zibr(logistic.cov = df_cov, 
                     beta.cov = df_cov, 
                     Y =   df_ra[,i] , 
                     subject.ind = metadf$Subject_ID_new,
                     time.ind =  metadf$New_Time)
    zibr.fit
  })
  
  for(i in 1:nOTU){
    #ZIBR----
    zibr.fit <- zibr.res[[i]]
    p_time_m1[i]=zibr.fit$beta.est.table[2,2]
    p_group_m1[i]=zibr.fit$beta.est.table[3,2]
    p_group_time_m1[i]=zibr.fit$beta.est.table[4,2]
  }
  
  p_value_DF<-data.frame(p_time=p_time_m1,
                     p_group=p_group_m1,
                     p_group_time=p_group_time_m1)
  
  p_value_DF<-cbind(OTU= colnames(df_ra),p_value_DF,
                  mutate_each(p_value_DF,funs(significance)))
  
  return(p_value_DF)
}

splineR<-function(df_ra,metadf,nOTU=300){
  
  p_time<-c()
  p_group<-c()
  df_new<-cbind(metadf,df_ra)
  
  options(warn=-1)#Suppress warnings 
  
  splinectomeR.Time <-lapply(1:nOTU, function(i){
    
    args <- list(data = df_new, x = 'New_Time', y = colnames(df_ra[i]),
                 cases = 'Subject_ID_new', perms = 999)
    
    fit <- do.call('trendyspliner', args = args)
    fit
  })
  splinectomeR.Group <-lapply(1:nOTU, function(i){
    
    args <- list(data = df_new, x = 'New_Time', y = colnames(df_ra[i]),
                 cases = 'Subject_ID_new', perms = 999, category = 'Group', groups = c('0','1'))
    
    fit <- do.call('permuspliner', args = args)
    fit
  })
  for(i in 1:nOTU){
    #ZIGMM----
    fit.Time <- splinectomeR.Time[[i]]
    p_time[i]=fit.Time$pval
    fit.Group <- splinectomeR.Group[[i]]
    p_group[i]=fit.Group$pval
  }
  
  p_value_DF<-data.frame(p_time=p_time,
                p_group=p_group)
 
  p_value_DF<-cbind(OTU= colnames(df_ra),p_value_DF,
                    mutate_each(p_value_DF,funs(significance)))
  
  return(p_value_DF)
}



