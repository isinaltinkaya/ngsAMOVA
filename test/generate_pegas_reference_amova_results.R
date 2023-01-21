#!/usr/bin/env Rscript

### Replicate ngsAMOVA results with pegas

library(pegas)

d<-read.csv("test_s9_d1_1K.sfs.csv")

d=d[d$Method=="gle",]
m<-matrix(NA,9,9)
diag(m)<-0
dcol<-"Dij2"
m[lower.tri(m,diag=FALSE)]<- d[dcol][[1]]
dd.d<-as.dist(m)
dd.pops<-factor(c(rep("pop1",3),rep("pop2",3),rep("pop3",3)))
amv<-pegas::amova(dd.d ~ dd.pops,is.squared=TRUE)
sig2<-setNames(amv$varcomp$sigma2,rownames(amv$varcomp))
phi<-pegas::getPhi(sig2)


out <- paste(c(sprintf("%.d",amv$tab$df[3]),
  sprintf("%.6f",amv$tab$SSD[3]),
  sprintf("%.6f",amv$tab$MSD[3]),
  sprintf("%.d",amv$tab$df[1]),
  sprintf("%.6f",amv$tab$SSD[1]),
  sprintf("%.6f",amv$tab$MSD[1]),
  sprintf("%.d",amv$tab$df[2]),
  sprintf("%.6f",amv$tab$SSD[2]),
  sprintf("%.6f",amv$tab$MSD[2]),
  sprintf("%.6f",phi[1]),
  sprintf("%.6f",amv$varcoef[[1]]),
  sprintf("%.6f",amv$varcomp$sigma2[1]),
  sprintf("%.6f",amv$varcomp$sigma2[2])))


write.table(out,file="test_s9_d1_1K_mtd1_pegas.amova.csv", row.names=FALSE,  col.names=FALSE,  quote=FALSE, sep=',')



# 2 level AMOVA test


d<-read.csv("test_s9_d1_1K.sfs.csv")

d=d[d$Method=="gle",]
m<-matrix(NA,9,9)
diag(m)<-0
dcol<-"Dij2"
m[lower.tri(m,diag=FALSE)]<- d[dcol][[1]]
dd.d<-as.dist(m)
dd.pops<-factor(c(rep("pop1",3),rep("pop2",3),rep("pop3",3)))
dd.regs<-factor(c(rep("reg1",6),rep("reg2",3)))
amv2<-pegas::amova(dd.d ~ dd.regs/dd.pops,is.squared=TRUE)
sig2_2<-setNames(amv2$varcomp$sigma2,rownames(amv2$varcomp))
phi_2<-pegas::getPhi(sig2_2)



out <- paste(c(
  sprintf("%.d",amv2$tab$df[4]),
  sprintf("%.6f",amv2$tab$SSD[4]),
  sprintf("%.6f",amv2$tab$MSD[4]),
  sprintf("%.d",amv2$tab$df[1]),
  sprintf("%.6f",amv2$tab$SSD[1]),
  sprintf("%.6f",amv2$tab$MSD[1]),
  sprintf("%.d",amv2$tab$df[2]),
  sprintf("%.6f",amv2$tab$SSD[2]),
  sprintf("%.6f",amv2$tab$MSD[2]),
  sprintf("%.d",amv2$tab$df[3]),
  sprintf("%.6f",amv2$tab$SSD[3]),
  sprintf("%.6f",amv2$tab$MSD[3]),
  sprintf("%.6f",phi_2[1,1]),
  sprintf("%.6f",phi_2[2,2]),
  sprintf("%.6f",phi_2[1,2]),
  sprintf("%.6f",amv2$varcoef[[1]]),
  sprintf("%.6f",amv2$varcoef[[2]]),
  sprintf("%.6f",amv2$varcoef[[3]]),
  sprintf("%.6f",amv2$varcomp$sigma2[1]),
  sprintf("%.6f",amv2$varcomp$sigma2[2]),
  sprintf("%.6f",amv2$varcomp$sigma2[3])))

write.table(out,file="test_s9_d1_1K_mtd2_pegas.amova.csv", row.names=FALSE,  col.names=FALSE,  quote=FALSE, sep=',')




# 3 level AMOVA test


d<-read.csv("test_s9_d1_1K.sfs.csv")

d=d[d$Method=="gle",]
m<-matrix(NA,9,9)
diag(m)<-0
dcol<-"Dij2"
m[lower.tri(m,diag=FALSE)]<- d[dcol][[1]]
dd.d<-as.dist(m)
dd.subpops<-factor(c("subpop1","subpop2","subpop2","subpop3","subpop3","subpop4","subpop5","subpop6","subpop6"))
dd.pops<-factor(c(rep("pop1",3),rep("pop2",3),rep("pop3",3)))
dd.regs<-factor(c(rep("reg1",6),rep("reg2",3)))
amv3<-pegas::amova(dd.d ~ dd.regs/dd.pops/dd.subpops,is.squared=TRUE)
sig2_3<-setNames(amv3$varcomp$sigma2,rownames(amv3$varcomp))
phi_3<-pegas::getPhi(sig2_3)

amv3

out <- paste(c("pegas",sprintf("%.6f",c(amv2$tab$df[1],amv2$tab$SSD[1],amv2$tab$MSD[1],amv2$tab$df[2],amv2$tab$SSD[2],amv2$tab$MSD[2],amv2$tab$df[3],amv2$tab$SSD[3],amv2$tab$MSD[3],as.double(amv2$varcoef[[1]]),amv2$varcomp$sigma2[1],amv2$varcomp$sigma2[2],phi_2[1]))),collapse=",")
write.table(out,file="test_s9_d1_1K_mtd3_pegas.amova.csv", row.names=FALSE,  col.names=FALSE,  quote=FALSE, sep=',')





