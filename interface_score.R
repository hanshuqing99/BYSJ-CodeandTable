#处理PRISM——interface residue contacts
#计算概率
#Obs(M)/Obs(S)
#--------------
#Px(M)/Px(S)

#工作目录
setwd("/Users/joannahench/Documents/毕设/mmp_test/PRISM_out/interface/cdh5-cdh13")
#读入相互作用文件
d1<-read.table("cdh5-cdh13.txt")
##Target A
inter_A<-d1[-c(1,2),1]
#记录位点信息
itA<-data.frame()
k<-1
for(i in 1:length(inter_A)){
  pos<-as.character(unlist(strsplit(inter_A[i], split = "_")))[4]
  aa<-as.character(unlist(strsplit(inter_A[i], split = "_")))[3]
  chain<-as.character(unlist(strsplit(inter_A[i], split = "_")))[2]
  itA[k,1]<-chain
  itA[k,2]<-pos
  itA[k,3]<-aa
  k<-k+1
}
colnames(itA)<-c("chain","pos","aa")
itA<-unique(itA)


##Target B
inter_B<-d1[-c(1,2),3]
#记录位点信息
itB<-data.frame()
k<-1
for(i in 1:length(inter_B)){
  pos<-as.character(unlist(strsplit(inter_B[i], split = "_")))[4]
  aa<-as.character(unlist(strsplit(inter_B[i], split = "_")))[3]
  chain<-as.character(unlist(strsplit(inter_B[i], split = "_")))[2]
  itB[k,1]<-chain
  itB[k,2]<-pos
  itB[k,3]<-aa
  k<-k+1
}
colnames(itB)<-c("chain","pos","aa")
itB<-unique(itB)


#读入TA错义突变信息
m1<-read.csv("res_cdh13_mis.csv",header = T)
i1<-merge(m1,itA,by="pos")

#读入TB错义突变信息
m2<-read.csv("res_cdh5_mis.csv",header = T)
i2<-merge(m2,itB,by="pos")


resa<-data.frame()
q<-1
for(j in 1:nrow(i1)){
      resa[q,1]<-paste(i1[j,2],i1[j,4],i1[j,1],i1[j,3],";",sep="")
      q<-q+1
    }

#去重
resa<-unique(resa)
#输出
write.table(resa,"ta_mutalist.txt",sep='\t',col.names = F,row.names = F,quote = F)


resb<-data.frame()
q<-1
if(nrow(i2)>0){
  for(j in 1:nrow(i2)){
    resb[q,1]<-paste(i2[j,2],i2[j,4],i2[j,1],i2[j,3],";",sep="")
    q<-q+1
  }
}

#去重
resb<-unique(resb)
#输出
write.table(resb,"tb_mutalist.txt",sep='\t',col.names = F,row.names = F,quote = F)




#读入TA沉默突变信息
m1<-read.csv("res_cdh13_sil.csv",header = T)
i1<-merge(m1,itA,by="pos")
#读入TB沉默突变信息
m2<-read.csv("res_cdh5_sil.csv",header = T)
i2<-merge(m2,itB,by="pos")

resa_s<-data.frame()
q<-1
if(nrow(i1)>0){
  for(j in 1:nrow(i1)){
    resa_s[q,1]<-paste(i1[j,2],i1[j,4],i1[j,1],i1[j,3],";",sep="")
    q<-q+1
  }
}
#去重
resa_s<-unique(resa_s)

resb_s<-data.frame()
q<-1
if(nrow(i2)>0){
  for(j in 1:nrow(i2)){
    resb_s[q,1]<-paste(i2[j,2],i2[j,4],i2[j,1],i2[j,3],";",sep="")
    q<-q+1
  }
}
#去重
resb_s<-unique(resb_s)

## 记录突变
#观察——错义突变/沉默突变
obsM<-nrow(resa)+nrow(resb)
obsS<-nrow(resa_s)+nrow(resb_s)

#读入氨基酸突变概率表
setwd("/Users/joannahench/Documents/毕设/mmp_test")
#读入codon mis/sil table
table<-read.csv("sum.csv",header=T)
#计算TA的M和S
ma<-0
sa<-0
for(i in 1:nrow(itA)){
  ma<-ma+table[which(table$aa==itA[i,3]),2]
  sa<-sa+table[which(table$aa==itA[i,3]),3]
}

#计算TB的M和S
mb<-0
sb<-0
for(i in 1:nrow(itB)){
  mb<-mb+table[which(table$aa==itB[i,3]),2]
  sb<-sb+table[which(table$aa==itB[i,3]),3]
}



#估计——错义突变/沉默突变
pxM<-ma+mb
pxS<-sa+sb


if(obsS==0){
  obsS<-1
}
inter_socre<-(obsM/obsS)/(pxM/pxS)




