
setwd("/Users/joannahench/Documents/毕设/mmp_test/PRISM_out/interface/cdh5-cdh13")
uni<-read.csv("cdh13.csv",header = T)

#沉默突变
sil<-subset(uni,select = c(AA.Mutation,Type),subset = (uni$Type=="Substitution - coding silent"))

#错义突变
mis<-subset(uni,select = c(AA.Mutation,Type),subset = (uni$Type=="Substitution - Missense"))


AA<-sil[,1]
#如果存在“p.”
del1<-vector()
for(i in 1:length(AA)){
  #用gsub函数把"p."去掉
  dela<-gsub("p.","",AA[i])
  del1<-c(del1,dela)
}
aat<-data.frame()
k<-1
for(i in 1:length(del1)){
  #用gsub函数把非数字的去掉
  t1<-gsub("[^0-9]","",del1[i])
  #用gsub函数把数字的去掉
  t2<-gsub("[0-9]"," ",del1[i])
  len<-length(as.character(unlist(strsplit(t2, split = " "))))
  t3<-as.character(unlist(strsplit(t2, split = " ")))[1]
  t4<-as.character(unlist(strsplit(t2, split = " ")))[len]
  aat[k,1]<-t3
  aat[k,2]<-t1
  aat[k,3]<-t4
  k<-k+1
}
colnames(aat)<-c("wild","pos","muta")
### 如果不转id
write.csv(aat,"res_cdh13_sil.csv",row.names = F)
### 如果要转id
setwd("/Users/joannahench/Documents/毕设/mmp_test/pdbprofiling")
pp<-read.csv("5ue3A.csv",header = T)
prof<-pp[,c(2,6)]
colnames(prof)<-c("pos","pos_pdb")
out<-merge(aat,prof,by="pos")
out1<-cbind(out[,2],out[,4],out[,3])
colnames(out1)<-c("wild","pos","muta")
setwd("/Users/joannahench/Documents/毕设/mmp_test/PRISM_out/interface/mmp1-mmp9")
write.csv(out1,"res_mmp9_sil.csv",row.names = F)




