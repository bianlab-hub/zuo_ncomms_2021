library(ggplot2)
setwd("~/Documents/project/meiosis/RNA-seq/20210412/")
stages<-c('PreL1','PreL2','L1','L2','Z1','Z2')
for (stage in stages){
  a<-read.csv(paste(stage,'_100kb_fpkm.bedGraph',sep = ""),header = F,sep = '\t')
  b<-a[which(a$V1 == "1"),]
  d=loess.smooth(b$V3, b$V4, span = 1/100, degree = 1,
                 family ="gaussian", evaluation = length(b$V1))
  b$V4<-d$y 
  b$V5<-log10(b$V4/median(b$V4))
  b<-b[,-4]
  #b[which(b$V5 == 'NaN'),]$V5<-0
  b$V1<-paste('chr',b$V1,sep = "")
  #b<-na.omit(b$V5)
  write.table(b, paste(stage,'_100kb_chr1.bed',sep=""), sep='\t', col.names = F, row.names = FALSE, quote = FALSE)
  write.table(b, paste(stage,'_100kb_chr1.bedGraph',sep=""), sep='\t', col.names = F, row.names = FALSE, quote = FALSE)
}

for (stage in stages){
  a<-read.csv(paste(stage,'_100kb.eigs.cis.vecs.txt',sep=""),header = T,sep = '\t')
  b<-a[,c(1,2,3,5)]
  b<-na.omit(b)
  write.table(b, paste(stage,'_100kb_chr1_comp.bed',sep=""), sep='\t', col.names = F, row.names = FALSE, quote = FALSE)
  write.table(b, paste(stage,'_100kb_chr1_comp.bedGraph',sep=""), sep='\t', col.names = F, row.names = FALSE, quote = FALSE)
}

a<-read.csv('PreL1_chr1_100kb.bed',header = F,sep = '\t')
ggplot(a,aes(x=V8,y=V4))+geom_point(shape=19)+scale_y_continuous(limits = c(-1.2, 1))+scale_x_continuous(limits = c(-2, 1))
p=ggplot(a,aes(x=V8,y=V4))+geom_point(shape=19)+scale_y_continuous(limits = c(-1.2, 1))+scale_x_continuous(limits = c(-2, 1))
ggsave("PreL1.pdf", p)

a<-read.csv('PreL2_chr1_100kb.bed',header = F,sep = '\t')
ggplot(a,aes(x=V8,y=V4))+geom_point(shape=19)+scale_y_continuous(limits = c(-1.8, 1))+scale_x_continuous(limits = c(-2, 1))
p=ggplot(a,aes(x=V8,y=V4))+geom_point(shape=19)+scale_y_continuous(limits = c(-1.8, 1))+scale_x_continuous(limits = c(-2, 1))
ggsave("PreL2.pdf", p)

a<-read.csv('L1_chr1_100kb.bed',header = F,sep = '\t')
ggplot(a,aes(x=V8,y=V4))+geom_point(shape=19)+scale_y_continuous(limits = c(-1.2, 1.2))+scale_x_continuous(limits = c(-1.2, 1))
p=ggplot(a,aes(x=V8,y=V4))+geom_point(shape=19)+scale_y_continuous(limits = c(-1.2, 1.2))+scale_x_continuous(limits = c(-1.2, 1))
ggsave("L1.pdf", p)

a<-read.csv('L2_chr1_100kb.bed',header = F,sep = '\t')
ggplot(a,aes(x=V8,y=V4))+geom_point(shape=19)+scale_y_continuous(limits = c(-1, 1.1))+scale_x_continuous(limits = c(-1, 1))
p=ggplot(a,aes(x=V8,y=V4))+geom_point(shape=19)+scale_y_continuous(limits = c(-1, 1.1))+scale_x_continuous(limits = c(-1, 1))
ggsave("L2.pdf", p)

a<-read.csv('Z1_chr1_100kb.bed',header = F,sep = '\t')
ggplot(a,aes(x=V8,y=V4))+geom_point(shape=19)+scale_y_continuous(limits = c(-1.2, 1.2))+scale_x_continuous(limits = c(-1.2, 1))
p=ggplot(a,aes(x=V8,y=V4))+geom_point(shape=19)+scale_y_continuous(limits = c(-1.2, 1.2))+scale_x_continuous(limits = c(-1.2, 1))
ggsave("Z1.pdf", p)

a<-read.csv('Z2_chr1_100kb.bed',header = F,sep = '\t')
ggplot(a,aes(x=V8,y=V4))+geom_point(shape=19)+scale_y_continuous(limits = c(-1.5, 1))+scale_x_continuous(limits = c(-1.1, 1))
p=ggplot(a,aes(x=V8,y=V4))+geom_point(shape=19)+scale_y_continuous(limits = c(-1.5, 1))+scale_x_continuous(limits = c(-1.1, 1))
ggsave("Z2.pdf", p)
