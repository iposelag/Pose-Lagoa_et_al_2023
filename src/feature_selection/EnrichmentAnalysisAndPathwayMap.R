#!/bin/Rscript
###############################################################################
#################### Analisis de enriquecimiento ##############################
############################ con Reactome #####################################
###############################################################################

#### Analisis de enriquecimiento, con Reactome ####
if("Reactome"%in%list.files("Files/")==FALSE){dir.create("Files/Reactome")}
genesets<-list.files("Files/GeneSets/")
background<-read.csv2("Files/GeneSets/background.txt",stringsAsFactors = F,sep="\t",header=F)$V1
lreac<-list()
for(a in c(3,9,2,4)){
  tt<-read.csv2(paste("Files/GeneSets/",genesets[a],sep=""),stringsAsFactors = F,sep="\t",header=F)$V1
  ## enrichr ##
  enrich<-enrichr(tt,databases = "Reactome_2022")
  enrichresult<-enrich$Reactome_2022
  enrichresult<-enrichresult[,c(1:4,7:9)]
  enrichresult<-enrichresult[which(enrichresult$Adjusted.P.value<=0.05),]
  write.table(enrichresult,paste("Files/Reactome/",genesets[a],sep=""),quote=F,sep="\t",row.names = F)
  lreac[[gsub(".txt","",genesets[a])]]<-enrichresult
}

## Read Reactome pathways' information ##
reactomepathways<-read.csv2("Files/ReactomeFiles/ReactomePathways.txt",stringsAsFactors = F,sep="\t",header = F)
reactomepathways<-reactomepathways[grep("HSA",reactomepathways$V1),]
reacpathfactor<-reactomepathways$V2 ; names(reacpathfactor)<-reactomepathways$V1
## Read parent-child relationships and identify the pathways categories
padreshijos<-read.csv2("Files/ReactomeFiles/Red_reactome_padres.txt",stringsAsFactors = F,sep="\t")
padres<-setdiff(padreshijos$Padres,padreshijos$Hijos)
## For each pathway we are going to obtain the entire list of sons iteractively ##
parentschilds<-c()
for(a in 1:length(padres)){
  allchildren<-c()
  childs<-padreshijos$Hijos[which(padreshijos$Padres==padres[a])]
  while(length(childs)>0){
    allchildren<-c(allchildren,childs)
    childs2<-c(); for(b in 1:length(childs)){childs2<-c(childs2,padreshijos$Hijos[which(padreshijos$Padres==childs[b])])}
    childs<-unique(childs2)
  }
  parentschilds<-rbind(parentschilds,cbind(allchildren,padres[a]))
}
colnames(parentschilds)<-c("Hijos","Padres")

lreac2<-list()
for(a in 1:length(names(lreac))){
  thecategory<-c()
  terms<-gsub(".+ R-HSA","R-HSA",lreac[[a]]$Term)
  for(b in 1:length(terms)){
    cual<-which(parentschilds[,1]==terms[b])
    theparent<-as.character(reacpathfactor[parentschilds[cual,2]])
    if(length(theparent)>1){print(paste(b,"ojo! Pertenece a dos categorias!"))}
    if(length(cual)>0){thecategory<-c(thecategory,paste(theparent,collapse = ";"))}
    if(length(cual)==0){
      cual2<-which(parentschilds[,2]==terms[b])
      if(length(cual2)>0){thecategory<-c(thecategory,as.character(reacpathfactor[terms[b]]))}
      if(length(cual2)==0){thecategory<-c(thecategory,"")}
    }
  }
  lreac2[[names(lreac)[a]]]<-cbind(lreac[[a]],thecategory)
}

#### Create pathway-pathway connections ####
if("Networks"%in%list.files("")==FALSE){dir.create("Networks")}
net<-c()

## DEA ##
## @ @ ##
for(a in 1:(length(lreac2$dea$Term)-1)){
  one<-lreac2$dea[a,]
  onegenes<-strsplit(one$Genes,";")[[1]]
  ## Intra ##
  for(b in (a+1):length(lreac2$dea$Term)){
    two<-lreac2$dea[b,]
    twogenes<-strsplit(two$Genes,";")[[1]]
    net<-rbind(net,c(gsub(" R-HSA-.+","",lreac2$dea$Term[a]),gsub(" R-HSA-.+","",lreac2$dea$Term[b]),
                     length(intersect(onegenes,twogenes))/length(unique(c(onegenes,twogenes))),"DEA","Intra"))
  }
  ## Inter MinMax
  for(b in 1:length(lreac2$minmax$Term)){
    two<-lreac2$minmax[b,]
    twogenes<-strsplit(two$Genes,";")[[1]]
    net<-rbind(net,c(gsub(" R-HSA-.+","",lreac2$dea$Term[a]),gsub(" R-HSA-.+","",lreac2$minmax$Term[b]),
                     length(intersect(onegenes,twogenes))/length(unique(c(onegenes,twogenes))),"DEA-minmax","Inter"))
  }
  ## Inter DisGeNET
  for(b in 1:length(lreac2$disgenet$Term)){
    two<-lreac2$disgenet[b,]
    twogenes<-strsplit(two$Genes,";")[[1]]
    net<-rbind(net,c(gsub(" R-HSA-.+","",lreac2$dea$Term[a]),gsub(" R-HSA-.+","",lreac2$disgenet$Term[b]),
                     length(intersect(onegenes,twogenes))/length(unique(c(onegenes,twogenes))),"DEA-disgenet","Inter"))
  }
}

## minmax ##
## @@  @@ ##
for(a in 1:(length(lreac2$minmax$Term)-1)){
  one<-lreac2$minmax[a,]
  onegenes<-strsplit(one$Genes,";")[[1]]
  ## Intra
  for(b in (a+1):length(lreac2$minmax$Term)){
    two<-lreac2$minmax[b,]
    twogenes<-strsplit(two$Genes,";")[[1]]
    net<-rbind(net,c(gsub(" R-HSA-.+","",lreac2$minmax$Term[a]),gsub(" R-HSA-.+","",lreac2$minmax$Term[b]),
                     length(intersect(onegenes,twogenes))/length(unique(c(onegenes,twogenes))),"minmax","Intra"))
  }
  ## Inter DisGeNET
  for(b in 1:length(lreac2$disgenet$Term)){
    two<-lreac2$disgenet[b,]
    twogenes<-strsplit(two$Genes,";")[[1]]
    net<-rbind(net,c(gsub(" R-HSA-.+","",lreac2$minmax$Term[a]),gsub(" R-HSA-.+","",lreac2$disgenet$Term[b]),
                     length(intersect(onegenes,twogenes))/length(unique(c(onegenes,twogenes))),"minmax-disgenet","Inter"))
  }
}

## disgenet ##
## @@ @@ @@ ##
for(a in 1:(length(lreac2$disgenet$Term)-1)){
  # a<-1
  one<-lreac2$disgenet[a,]
  onegenes<-strsplit(one$Genes,";")[[1]]
  ## Intra
  for(b in (a+1):length(lreac2$disgenet$Term)){
    two<-lreac2$disgenet[b,]
    twogenes<-strsplit(two$Genes,";")[[1]]
    net<-rbind(net,c(gsub(" R-HSA-.+","",lreac2$disgenet$Term[a]),gsub(" R-HSA-.+","",lreac2$disgenet$Term[b]),
                     length(intersect(onegenes,twogenes))/length(unique(c(onegenes,twogenes))),"disgenet","Intra"))
  }
}

## Write the table with the network information, needed for the network representation ##
colnames(net)<-c("Pathway1","Pathway2","JaccardIndex","Interaction","Type")
write.table(net,"Networks/Network.txt",quote=F,sep="\t",row.names = F)

## Create a table with node information ##
tabla<-rbind(cbind(gsub(" R-HSA-.+","",lreac2$dea[,1]),lreac2$dea$thecategory,"DEA"),
             cbind(gsub(" R-HSA-.+","",lreac2$minmax[,1]),lreac2$minmax$thecategory,"minmax"),
             cbind(gsub(" R-HSA-.+","",lreac2$disgenet[,1]),lreac2$disgenet$thecategory,"disgenet"))
sharedpaths<-table(tabla[,1])[which(table(tabla[,1])>1)]
quitar<-c() ; newtab<-c()
for(a in 1:length(sharedpaths)){
  q<-which(tabla[,1]==names(sharedpaths[a]))
  quitar<-c(quitar,q)
  newtab<-rbind(newtab,c(tabla[q[1],1:2],paste(tabla[q,3],collapse = ";")))
}
tabla2<-rbind(tabla[-quitar,],newtab)
colnames(tabla2)<-c("Pathway","Category","Source")
coltab<-read.table("Files/ReactomeFiles/ReactomeCategoryColors.txt",stringsAsFactors = F,sep="\t",header=T)[,1]
colors<-c("#006347","#425D10","#89B524","#019477","#97D0A7",
          "#BCDA78","#AF9F01","#0091A0","#54B9C1","#BBDAF6",
          "#386373","#5173BD","#082F6D","#56187D","#9C58A1",
          "#EEA187","#EDA7C1","#D38895","#AC0123","#DD3F4E",
          "#DB2879","#DF4018","#EF9F26","#FFE001","#411C01",
          "#806726","#DEC990","#E2B53E","#000000")
names(colors)<-coltab
tabla2<-cbind(tabla2,as.character(colors[tabla2[,2]]))
tabla2[which(is.na(tabla2[,4])),4]<-"#FFFFFF"
colnames(tabla2)[4]<-"Colors"
## Write the table with the node characteristics (if it comes from DEA, minmax or DisGeNET) and the category and color of the reactome pathways ##
write.table(tabla2,"Networks/Tabla.txt",quote=F,sep="\t",row.names = F)
