#!/bin/Rscript
###############################################################################
###################### RadarChart plots generation  ##########################
###############################################################################

#### Plot hexagono ####
# install.packages("fmsb")
library("fmsb"); library(dplyr)
test <- read.csv2("../../data/ML_out/two_class/optm/tables/results_test.csv",stringsAsFactors = F,sep=",")
colnames(test)[3] = "mean"
train <-  read.csv2("../../data/ML_out/two_class/optm/tables/results_train.csv",stringsAsFactors = F,sep=",")
# 
#### RadarChart una metrica todos los clasificadores y feature selections ####
plotradarchart<-function(tabla,metrica){
  tablita<-tabla[which(tabla$metric==metrica),]
  classifiers<-unique(tablita$classifier)
  mlinput<-unique(tablita$ml_input)
  # mlinput<-c("disgenet_expansion", "data_driven_expansion","disgenet_curated_guildify", "data_driven_guildify")
  tvalues<-c() ; bestselection<-c() ; bestvalue<-c()
  for(a in 1:length(classifiers)){
    # a<-1
    cuala<-which(tablita$classifier==classifiers[a])
    values<-c()
    for(b in 1:length(mlinput)){
      # b<-1
      cualb<-which(tablita$ml_input==mlinput[b])
      values<-c(values,as.numeric(tablita$mean[intersect(cuala,cualb)]))
    }
    mejorvalor<-which(values==max(values))
    ## An input works better than the rest ##
    if(length(mejorvalor)==1){
      bestselection<-c(bestselection,which(values==max(values)))
      bestvalue<-c(bestvalue,values[which(values==max(values))])
    }
    ## Several inputs are equally the best
    if(length(mejorvalor)>1){
      bestselection<-c(bestselection,7)
      bestvalue<-c(bestvalue,values[which(values==max(values))[1]])
    }
    tvalues<-cbind(tvalues,values)
  }
  ## Anadimos dos lineas con el valor minimo y el maximo con cada seleccion de genes ##
  tvalues <- as.data.frame(rbind(rep(0.9,6) , rep(0.5,6) , tvalues))
  colnames(tvalues)<-c(classifiers)
  rownames(tvalues)<-c("max","min",mlinput)
  ## Sacamos los datos
  ## Le damos los colores de Iria ##
  
  # colores<-c("#9AAAC3","#8EBCBF","#cfe2f3","#d2e7d6","black")
  colores<-c("#325486","#9AAAC3","#bea9de","#19787F","#8EBCBF","#FF9999","#ffbaba","#d2e7d6","#cfe2f3", "black")
  nombres<-c("COPD-related curated","COPD-related curated Omnipath","COPD-related entire","data-driven","data-driven OmniPath",
             "OmniPath intersection","OmniPath union", "data-driven GUILDify", "COPD-related curated GUILDify")
  # nombres<-c("COPD-related omnipath","data-driven omnipath", "disgenet guildify", "data-driven guildify")
  ## Hacemos el plot ##
  ## Labels indicando el valor mas alto de accuracy para el clasificador, coloreado por la seleccion de genes que lo consigue ##
  pdf(file=paste("../../data/ML_out/two_class/optm/tables/",name,"_",
                 metrica,".pdf", sep = ""),width = 12,height = 10)
    radarchart(tvalues,axistype = 2,pcol=colores,plwd=3,plty=1,
               cglcol="grey", axislabcol=colores[bestselection],
               paxislabels=round(bestvalue,3),
               cglwd=2,vlcex=1.2)
    legend(x=0.8, y=1.35, legend = nombres, bty = "n", pch=20 , col=colores , text.col = "black", cex=1.2, pt.cex=1.5)
  dev.off()
  ## Labels en las lineas indicando los valores de cada una de ellas ##
  pdf(file=paste("Test",metrica,"PlotValorLineas.pdf",sep="_"),width = 12,height = 10)
    radarchart(tvalues,axistype = 1,pcol=colores,plwd=3,plty=1,
               cglcol="grey", axislabcol="grey", caxislabels=round(seq(min(tvalues[2,]),max(tvalues[1,]),(max(tvalues[1,])-min(tvalues[2,]))/4),3),
               cglwd=2,vlcex=1.2)
    legend(x=0.7, y=1.3, legend = rownames(tvalues[-c(1,2),]), bty = "n", pch=20 , col=colores , text.col = "black", cex=1.2, pt.cex=1.5)
  # dev.off()
  ## Labels en ambas ##
  # pdf(file=paste("Test",metrica,"Plot.pdf",sep="_"),width = 12,height = 10)
  #   radarchart(tvalues,axistype = 3,pcol=colores,plwd=3,plty=1,
  #              cglcol="grey", axislabcol="grey", caxislabels=round(seq(min(tvalues[2,]),max(tvalues[1,]),(max(tvalues[1,])-min(tvalues[2,]))/4),3),
  #              paxislabels=nombres[bestselection],axislabels=c(1:5),
  #              cglwd=2,vlcex=1.2)
  #   legend(x=0.7, y=1.3, legend = rownames(tvalues[-c(1,2),]), bty = "n", pch=20 , col=colores , text.col = "black", cex=1.2, pt.cex=1.5)
  # dev.off()
}

# tabla<-test
# metrica<-"bal_accuracy"
name = "test";
plotradarchart(test, "mcc_normalized");
plotradarchart(test, "accuracy"); plotradarchart(test, "bal_accuracy")
name = "train";
plotradarchart(train, "mcc_normalized");
plotradarchart(train, "accuracy"); plotradarchart(train, "bal_accuracy")

#### RadarCharts, uno por clasificador, con las feature selections en las esquinas y las métricas representadas dentro ####
## Ponemos la tabla como toca ##
tabla<-test
classifiers2<-unique(tabla$classifier)
ltvalues<-list() ; allvals<-c()
metrics_1 <- c("bal_accuracy","accuracy","sens","spec")
metrics_2 <- c("macro_F1", "precision", "recall") # copd positive
metrics_3 <- c("mcc", "bal_accuracy", "accuracy") # ctrl positive
for(z in 1:length(classifiers2)){
  classifiers<-classifiers2[z]
  tablita<-tabla[which(tabla$classifier==classifiers),]
  mlinput<-unique(tablita$ml_input)
  metricas<-c("mcc", "bal_accuracy", "accuracy")

  tvalues<-c() ; bestselection<-c() ; bestvalue<-c()
  for(a in 1:length(mlinput)){
    # a<-1
    cuala<-which(tablita$ml_input==mlinput[a])
    values<-c()
    for(b in 1:length(metricas)){
      # b<-1
      cualb<-which(tablita$metric==metricas[b])
      values<-c(values,as.numeric(tablita$mean[intersect(cuala,cualb)]))
    }
    mejorvalor<-which(values==max(values))
    ## An input works better than the rest ##
    if(length(mejorvalor)==1){
      bestselection<-c(bestselection,which(values==max(values)))
      bestvalue<-c(bestvalue,values[which(values==max(values))])
    }
    ## Several inputs are equally the best
    if(length(mejorvalor)>1){
      bestselection<-c(bestselection,7)
      bestvalue<-c(bestvalue,values[which(values==max(values))[1]])
    }
    tvalues<-cbind(tvalues,values)
  }
  colnames(tvalues)<-c(mlinput)
  rownames(tvalues)<-metricas
  ltvalues[[classifiers]]$tvalues<-tvalues
  ltvalues[[classifiers]]$bestselection<-bestselection
  ltvalues[[classifiers]]$bestvalue<-bestvalue
  allvals<-c(allvals,as.numeric(tvalues))
}

## Le damos los colores de Iria ##
colores<-c("#aa6f73","#83adb5","#c7bbc9","#5e3c58")
## Add min and max values ##
for(a in 1:length(names(ltvalues))){
  # a<-1
  # ltvalues[[a]]$tvalues<-as.data.frame(rbind(rep(max(allvals),6) , rep(min(allvals),6) , ltvalues[[a]]$tvalues))
  ltvalues[[a]]$tvalues<-as.data.frame(rbind(rep(1,6) , rep(0.5,6) , ltvalues[[a]]$tvalues))
  rownames(ltvalues[[a]]$tvalues)[1:2]<-c("max","min")
}
## Plot ##
pdf(file="../../data/ML_out/two_class/optm/data_driven/metrics/Train_mcc_balacc_acc.pdf",width = 22,height = 14)
par(mfrow=c(2,3))
for(a in 1:length(names(ltvalues))){
  radarchart(ltvalues[[a]]$tvalues,axistype = 2,pcol=colores,plwd=3,plty=1,
             cglcol="grey", axislabcol=colores[ltvalues[[a]]$bestselection],
             paxislabels=round(ltvalues[[a]]$bestvalue,3),
             cglwd=2,vlcex=1.2)
  legend(x=0.7, y=1.3, legend = rownames(ltvalues[[a]]$tvalues[-c(1,2),]), bty = "n", pch=20 , col=colores , text.col = "black", cex=1.2, pt.cex=1.5)
  text(-1,1.2,names(ltvalues)[a],cex=4)
}
dev.off()
