#### Plot hexagono ####
# install.packages("fmsb")
library("fmsb")
test<-read.csv2("../../reports/tables/results_test_dd_dgn.csv",stringsAsFactors = F,sep=",")

## Ponemos la tabla como toca ##
testacc<-test[which(test$metric=="sens"),]
classifiers<-unique(testacc$classifier)
mlinput<-unique(testacc$ml_input)
tvalues<-c() ; bestselection<-c() ; bestvalue<-c()
for(a in 1:length(classifiers)){
  # a<-1
  cuala<-which(testacc$classifier==classifiers[a])
  values<-c()
  for(b in 1:length(mlinput)){
    # b<-1
    cualb<-which(testacc$ml_input==mlinput[b])
    values<-c(values,as.numeric(testacc$estimate[intersect(cuala,cualb)]))
  }
  bestselection<-c(bestselection,which(values==max(values)))
  bestvalue<-c(bestvalue,values[which(values==max(values))])
  tvalues<-cbind(tvalues,values)
}
## Anadimos dos lineas con el valor minimo y el maximo con cada seleccion de genes ##
tvalues <- as.data.frame(rbind(rep(max(tvalues),6) , rep(min(tvalues),6) , tvalues))
colnames(tvalues)<-c("RF","SVM-rad","SVM-poly","GLM","kNN","XGB")
rownames(tvalues)<-c("max","min",mlinput)

## Sacamos los datos 
## Le damos los colores de Iria ##
colores<-c("#005b96","#9AAAC3","#006666","#8EBCBF","#ed5555","#ffbaba")
nombres<-c("COPD-related","COPD-related expansion","data-driven","data-driven expansion","expansion intersection","expansion union")

## Hacemos el plot ##s
## Labels indicando el valor mas alto de accuracy para el clasificador, coloreado por la seleccion de genes que lo consigue ##
pdf(file="TestAccPlotValorMejorValorClasificador.pdf",width = 12,height = 10)
radarchart(tvalues,axistype = 2,pcol=colores,plwd=3,plty=1,
           cglcol="grey", axislabcol=colores[bestselection],
           paxislabels=round(bestvalue,3),
           cglwd=2,vlcex=1.2)
legend(x=0.7, y=1.3, legend = rownames(tvalues[-c(1,2),]), bty = "n", pch=20 , col=colores , text.col = "black", cex=1.2, pt.cex=1.5)
dev.off()

## Labels en las lineas indicando los valores de cada una de ellas ##
pdf(file="TestAccPlotValorLineas.pdf",width = 12,height = 10)
radarchart(tvalues,axistype = 1,pcol=colores,plwd=3,plty=1,
           cglcol="grey", axislabcol="grey", caxislabels=round(seq(min(tvalues[2,]),max(tvalues[1,]),(max(tvalues[1,])-min(tvalues[2,]))/4),3),
           cglwd=2,vlcex=1.2)
legend(x=0.7, y=1.3, legend = rownames(tvalues[-c(1,2),]), bty = "n", pch=20 , col=colores , text.col = "black", cex=1.2, pt.cex=1.5)
dev.off()

## Labels en ambas ##
pdf(file="TestAccPlot.pdf",width = 12,height = 10)
radarchart(tvalues,axistype = 3,pcol=colores,plwd=3,plty=1,
           cglcol="grey", axislabcol="grey", caxislabels=round(seq(min(tvalues[2,]),max(tvalues[1,]),(max(tvalues[1,])-min(tvalues[2,]))/4),3),
           paxislabels=nombres[bestselection],axislabels=c(1:5),
           cglwd=2,vlcex=1.2)
legend(x=0.7, y=1.3, legend = rownames(tvalues[-c(1,2),]), bty = "n", pch=20 , col=colores , text.col = "black", cex=1.2, pt.cex=1.5)
dev.off()
