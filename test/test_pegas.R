#!/usr/bin/env Rscript
### Comparing ngsAMOVA results with pegas
library(pegas)
args <- commandArgs(trailingOnly=TRUE)


d<-read.csv(args[1])

d$Pop1<-unlist(lapply(strsplit(as.character(d$Ind1),split="_"),"[[",1))
d$Pop2<-unlist(lapply(strsplit(as.character(d$Ind2),split="_"),"[[",1))
d=d[d$Method=="gle",]
m<-matrix(NA,9,9)
diag(m)<-0

#
# dcol<-"Dij"
# amv<-pegas::amova(dd.d ~ dd.pops,is.squared=FALSE)
#

dcol<-"Dij2"

m[lower.tri(m,diag=FALSE)]<- d[dcol][[1]]
dd.d<-as.dist(m)
dd.pops<-factor(c(rep("pop1",3),rep("pop2",3),rep("pop3",3)))

amv<-pegas::amova(dd.d ~ dd.pops,is.squared=TRUE)


sig2<-setNames(amv$varcomp$sigma2,rownames(amv$varcomp))
phi<-getPhi(sig2)


amvfile<-paste0(strsplit(args[1],".sfs.csv")[[1]],".amova.csv")

outfile<-paste0(args[1],"_testPegas.csv")

cat(paste0("\nPrinting to file: ",outfile))
cat("\n\n")
pegas_res<- data.frame("pegas",amv$tab$df[1],amv$tab$SSD[1],amv$tab$MSD[1],amv$tab$df[2],amv$tab$SSD[2],amv$tab$MSD[2],amv$tab$df[3],amv$tab$SSD[3],amv$tab$MSD[3],as.double(amv$varcoef[[1]]),amv$varcomp$sigma2[1],amv$varcomp$sigma2[2],phi[1])

colnames(pegas_res)<- c("Type","df1","SSD1","MSD1","df2","SSD2","MSD2","df3","SSD3","MSD3","varcoef","sigma2_1","sigma2_2","phi")

amv_res<-read.csv(amvfile,header=FALSE)
colnames(amv_res)<- c("Type","df1","SSD1","MSD1","df2","SSD2","MSD2","df3","SSD3","MSD3","varcoef","sigma2_1","sigma2_2","phi")


outtable=rbind(pegas_res,amv_res)

write.table(outtable,file=paste0(args[1],"_testPegas.csv"), row.names=FALSE,  col.names=FALSE,  quote=FALSE, sep=',')


