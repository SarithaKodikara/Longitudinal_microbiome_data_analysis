# Plotting the results from Scenario 1

boxplotFunction<-function(temp, l_lim=0, u_lim=1, plotTitle=""){
  sc_result<-c()
  for (i in 1:length(temp)) {
    tem<-unlist(str_split(temp[i], pattern=".csv"))[1]
    scRes<-read.csv(temp[i])
    scRes<-mutate(scRes,Method=rep(tem,50))
    sc_result<-dplyr::bind_rows(sc_result,scRes)
  }
  
  sc_long<-data.frame(Time=c(sc_result$Time_sen, sc_result$Time_spe),
                      Group=c(sc_result$Group_sen, sc_result$Group_spe),
                      Time_Group=c(sc_result$Group_Time_sen, sc_result$Group_Time_spe),
                      ErrorType=c(rep("Sensitivity",500),rep("Specificity",500)),
                      Method=rep(sc_result$Method,2))
  sc_long$Type<-ifelse((sc_long$Method=="ZIBR"|
                          sc_long$Method=="SplinectomeR"|sc_long$Method=="ZIGMM"|
                          sc_long$Method=="ZIGMMAr"),
                       "R.A","Count")
  sc_long$Method<-factor(sc_long$Method, levels=c("NBMM","NBMMAR","SplinectomeR","ZIBR","ZIGMM",
                                                  "ZIGMMAr","ZIGMMCo","ZIGMMCoAR","ZINBMM","ZINBMMAR"),
                         labels = c("NBMM","NBMM*","SplinectomeR`","ZIBR","ZIGMM",
                                    "ZIGMM*","ZIGMM","ZIGMM*","FZINBMM","FZINBMM*"))
  p1<-ggplot(sc_long, aes(x=Method, y=Time, fill=ErrorType)) + 
    geom_boxplot() +coord_flip()+ 
    scale_y_continuous(limits = c(l_lim, u_lim))+facet_grid(Type~., scales = "free_y")+ 
    ggtitle(plotTitle)+theme(axis.text=element_text(size=12),
                             axis.title=element_text(size=14,face="bold"))
  p2<-ggplot(sc_long, aes(x=Method, y=Group, fill=ErrorType)) + 
    geom_boxplot() +coord_flip()+ 
    scale_y_continuous(limits = c(l_lim, u_lim))+facet_grid(Type~., scales = "free_y")+ 
    ggtitle(plotTitle)+theme(axis.text=element_text(size=12),
                             axis.title=element_text(size=14,face="bold"))
  p3<-ggplot(sc_long, aes(x=Method, y=Time_Group, fill=ErrorType)) + 
    geom_boxplot(na.rm = TRUE)+coord_flip()+ 
    scale_y_continuous(limits = c(l_lim, u_lim))+facet_grid(Type~., scales = "free_y")+ 
    ggtitle(plotTitle)+theme(axis.text=element_text(size=12),
                             axis.title=element_text(size=14,face="bold"))
  
  return(list(p1,p2,p3))
}
