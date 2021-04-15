# hicrep
library(strawr)
library(hicrep)
library(rhdf5)
library(Rcpp)
setwd("/Users/wzuo/Documents/project/Coding/Hi-C/hicrep")

mat1 = cool2matrix('spg1_500kb.cool', chr = 'chr1' )
mat2 = cool2matrix('spg2_500kb.cool', chr = 'chr1' )
h_value <- htrain(mat1, mat2, resol = 500000, lbr = 0, ubr = 5000000, range = 0:10)

a<-c('Sertoli1','Sertoli2','spg1','spg2','PreL1','PreL2','L1','L2','Z1','Z2','P1','P2','D1','D2','MII1','MII2','Sun1_1','Sun1_2')
for (i in c(1:17)){
  for (j in c((i+1):18)){
    all.scc <- list()
    for (k in paste("chr", as.character(1:19))){
      mat1 = cool2matrix(paste(a[i],'_500kb.cool',sep = ""), chr = k )
      mat2 = cool2matrix(paste(a[j],'_500kb.cool',sep = ""), chr = k )
      all.scc[[k]] = get.scc(mat1, mat2, resol = 500000, h = 1, lbr = 0, ubr = 5000000)
    }
    write.table(all.scc,paste(a[i],'_',a[j],'_500kb.csv',sep = ""),col.names = T,append = F,sep = '\t')
  }
}