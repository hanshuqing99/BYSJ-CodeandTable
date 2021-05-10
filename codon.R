setwd("/Users/joannahench/Documents/research/毕设/mmp_test")
library(xlsx)
#读入codon table
codon<-read.xlsx("codon_table.xlsx",sheetIndex=1)
#存储20种氨基酸
t1<-levels(as.factor(codon[,1]))
#新表
sum<-data.frame()
m<-1
#对于每一个密码子进行循环
for(i in 1:nrow(codon)){

  mis<-0
  sil<-0
  for(j in 1:3){
    b1<-change(codon[i,2],j)
    for(k in 1:3){
      if(codon[which(codon$codon==b1[k]),1]==codon[i,1])
      {
        sil<-sil+1
        next
      }
      if(codon[which(codon$codon==b1[k]),1]=="stop"){
        next
      }
      else{
        mis<-mis+1
      }
    }
  }
  sum[m,1]<-codon[i,1]
  sum[m,2]<-codon[i,2]
  sum[m,3]<-mis
  sum[m,4]<-sil
  m<-m+1

}

write.csv(sum,"sum.csv")



change<-function(x,y){
  f0<-substring(x,1:3,1:3)
  strAA<-"ACGT"
  f1<-gsub(f0[y],"",strAA)
  f2<-substring(f1,1:3,1:3)
  out<-vector()
  if(y==1){
    out<-c(out,paste(f2[1],f0[2],f0[3],sep = ""),paste(f2[2],f0[2],f0[3],sep = ""),paste(f2[3],f0[2],f0[3],sep = ""))
  }
  if(y==2){
    out<-c(out,paste(f0[1],f2[1],f0[3],sep = ""),paste(f0[1],f2[2],f0[3],sep = ""),paste(f0[1],f2[3],f0[3],sep = ""))
  }
  if(y==3){
    out<-c(out,paste(f0[1],f0[2],f2[1],sep = ""),paste(f0[1],f0[2],f2[2],sep = ""),paste(f0[1],f0[2],f2[3],sep = ""))
  }
  return(out)
}



#读入codon mis/sil table
table<-read.csv("sum.csv",header=T)
#done！

