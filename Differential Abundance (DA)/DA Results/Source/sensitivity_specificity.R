

#Significant OTUs AT significance level 0.05
significance <- function(x) if_else(x<sig, 1, 0)

#### Function to calculate sensitivity and specificity results for ZIGMM, NBMM, FZINBMM####
se.sp<-function(df,metadf,t_effect,g_effect,g_t_effect,correlation_AR=FALSE, countData=FALSE, method){
  p_time_m1<-c()
  p_group_m1<-c()
  p_group_time_m1<-c()
  df_new<-cbind(metadf,df)
  
  options(warn=-1)#Suppress warnings 
  
  # Start the clock!
  ptm <- proc.time() 
  
  res.sum <-lapply(1:300, function(i){
    
    f<-if(countData) formula(paste(colnames(df[i]), "Time+Group+Time*Group+ offset(log(Library_size))", sep = "~")) else
      formula(paste(colnames(df[i]), "Time+Group+Time*Group", sep = "~"))
    args <- list(fixed = f, random = ~ 1 | Indiv, data = df_new)
    
    if (correlation_AR)
      args$correlation <- corAR1()
    
    if(method=="zigmm")
        fit <- do.call('lme.zig', args = args)
    if(method=="nbmm")
      fit <- do.call('glmm.nb', args = args)
    if(method=="zinbmm")
      fit <- do.call('glmm.zinb', args = args)
    summary(fit)
  })
  for(i in 1:300){
    #ZIGMM----
    fit.summary <- res.sum[[i]]
    p_time_m1[i]=fit.summary$tTable[2,5]
    p_group_m1[i]=fit.summary$tTable[3,5]
    p_group_time_m1[i]=fit.summary$tTable[4,5]
  }
  # Stop the clock
  print(proc.time() - ptm )
  runTime<-proc.time() - ptm 
  p<-data.frame(p_time=p_time_m1,
                      p_group=p_group_m1,
                      p_group_time=p_group_time_m1)
  p_adj<- p%>% apply(.,2,p.adjust,method = 'fdr') %>% data.frame()
  sig<-p_adj %>% mutate_each(.,funs(significance))
  #Confusion matrices
  conf_matrix_t<-table(factor(sig$p_time,levels = c(1,0)),
                       factor(t_effect,levels = c(1,0)))
  spe.t <- epi.tests(conf_matrix_t)$detail$sp$est
  sen.t<-epi.tests(conf_matrix_t)$detail$se$est
  
  conf_matrix_g<-table(factor(sig$p_group,levels = c(1,0)),
                       factor(g_effect,levels = c(1,0)))
  spe.g <- epi.tests(conf_matrix_g)$detail$sp$est
  sen.g<-epi.tests(conf_matrix_g)$detail$se$est
  
  conf_matrix_g_t<-table(factor(sig$p_group_time,levels = c(1,0)),
                         factor(g_t_effect,levels = c(1,0)))
  spe.g_t <- epi.tests(conf_matrix_g_t)$detail$sp$est
  sen.g_t<-epi.tests(conf_matrix_g_t)$detail$se$est
  
  se.sp_df<-data.frame(Time_sen=sen.t, Group_sen=sen.g, 
                       Group_Time_sen=sen.g_t,Time_spe=spe.t, Group_spe=spe.g, 
                       Group_Time_spe=spe.g_t)
  
  return(list(se.sp=se.sp_df, runT=runTime))
}

#### Function to calculate sensitivity and specificity results for splinectomeR####

se.sp_splineR<-function(df_ra,metadf,t_effect,g_effect){
  
  p_time<-c()
  p_group<-c()
  df_new<-cbind(metadf,df_ra)
  
  options(warn=-1)#Suppress warnings 
  
  # Start the clock!
  ptm <- proc.time()  
  
  splinectomeR.Time <-lapply(1:300, function(i){
    
    args <- list(data = df_new, x = 'Time', y = colnames(df_ra[i]),
                 cases = 'Indiv', perms = 999)
    
    fit <- do.call('trendyspliner', args = args)
    fit
  })
  splinectomeR.Group <-lapply(1:300, function(i){
    
    args <- list(data = df_new, x = 'Time', y = colnames(df_ra[i]),
                 cases = 'Indiv', perms = 999, category = 'Group', groups = c('0','1'))
    
    fit <- do.call('permuspliner', args = args)
    fit
  })
  for(i in 1:300){
    #ZIGMM----
    fit.Time <- splinectomeR.Time[[i]]
    p_time[i]=fit.Time$pval
    fit.Group <- splinectomeR.Group[[i]]
    p_group[i]=fit.Group$pval
  }
  # Stop the clock
  print(proc.time() - ptm )
  runTime<-proc.time() - ptm 
  p<-data.frame(p_time=p_time,
                p_group=p_group)
  p_adj<- p%>% apply(.,2,p.adjust,method = 'fdr') %>% data.frame()
  sig<-p_adj %>% mutate_each(.,funs(significance))
  #Confusion matrices
  conf_matrix_t<-table(factor(sig$p_time,levels = c(1,0)),
                       factor(t_effect,levels = c(1,0)))
  spe.t <- epi.tests(conf_matrix_t)$detail$sp$est
  sen.t<-epi.tests(conf_matrix_t)$detail$se$est
  
  conf_matrix_g<-table(factor(sig$p_group,levels = c(1,0)),
                       factor(g_effect,levels = c(1,0)))
  spe.g <- epi.tests(conf_matrix_g)$detail$sp$est
  sen.g<-epi.tests(conf_matrix_g)$detail$se$est
  
  
  se.sp_df<-data.frame(Time_sen=sen.t, Group_sen=sen.g,Time_spe=spe.t, Group_spe=spe.g)
  
  return(list(se.sp=se.sp_df, runT=runTime))
}

#### Function to calculate sensitivity and specificity results for ZIBR####
se.sp_zibr<-function(df_ra,metadf,t_effect,g_effect,g_t_effect){
  
  p_time_m1<-c()
  p_group_m1<-c()
  p_group_time_m1<-c()
  
  
  df_cov = data.frame(Time = metadf$Time, Group=metadf$Group, 
                      Time_Group=metadf$Time* metadf$Group)
  options(warn=-1)#Suppress warnings 
  
  # Start the clock!
  ptm <- proc.time()  
  
  zibr.res <-lapply(1:300, function(i){
    zibr.fit <- zibr(logistic.cov = df_cov, 
                     beta.cov = df_cov, 
                     Y =   df_ra[,i] , subject.ind = metadf$Indiv,
                     time.ind =  metadf$Time)
    zibr.fit
  })
  for(i in 1:300){
    #ZIBR----
    zibr.fit <- zibr.res[[i]]
    p_time_m1[i]=zibr.fit$beta.est.table[2,2]
    p_group_m1[i]=zibr.fit$beta.est.table[3,2]
    p_group_time_m1[i]=zibr.fit$beta.est.table[4,2]
  }
  # Stop the clock
  print(proc.time() - ptm )
  runTime<-proc.time() - ptm 
  zibr_p<-data.frame(p_time=p_time_m1,
                     p_group=p_group_m1,
                     p_group_time=p_group_time_m1)
  zibr_p_adj<- zibr_p%>% apply(.,2,p.adjust,method = 'fdr') %>% data.frame()
  zibr_sig<-zibr_p_adj %>% mutate_each(.,funs(significance))
  #Confusion matrices
  conf_matrix_t<-table(factor(zibr_sig$p_time,levels = c(1,0)),
                       factor(t_effect,levels = c(1,0)))
  spe.t <- epi.tests(conf_matrix_t)$detail$sp$est
  sen.t<-epi.tests(conf_matrix_t)$detail$se$est
  
  conf_matrix_g<-table(factor(zibr_sig$p_group,levels = c(1,0)),
                       factor(g_effect,levels = c(1,0)))
  spe.g <- epi.tests(conf_matrix_g)$detail$sp$est
  sen.g<-epi.tests(conf_matrix_g)$detail$se$est
  
  conf_matrix_g_t<-table(factor(zibr_sig$p_group_time,levels = c(1,0)),
                         factor(g_t_effect,levels = c(1,0)))
  spe.g_t <- epi.tests(conf_matrix_g_t)$detail$sp$est
  sen.g_t<-epi.tests(conf_matrix_g_t)$detail$se$est
  
  se.sp_df<-data.frame(Time_sen=sen.t, Group_sen=sen.g, 
                       Group_Time_sen=sen.g_t,Time_spe=spe.t, Group_spe=spe.g, 
                       Group_Time_spe=spe.g_t)
  
  return(list(se.sp=se.sp_df, runT=runTime))
}
