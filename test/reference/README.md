### test_s9_d1_1K_pegas.amova.csv

```
 ./ngsAMOVA -doEM 1 -doAMOVA 3 -doTest 0 -in test/test_s9_d1_1K.vcf -isSim 1 -minInd 2 -printMatrix 0  -doDist 1 -maxIter 100 -nThreads 0 -tole 1e-10  --hascolnames 1 -m test/metadata_with_header_1lvl.tsv 
```

```
#!/usr/bin/env Rscript
### Comparing ngsAMOVA results with pegas
library(pegas)
ff <- "test_s9_d1_1K.sfs.csv"
d<-read.csv(ff)
d$Pop1<-unlist(lapply(strsplit(as.character(d$Ind1),split="_"),"[[",1))
d$Pop2<-unlist(lapply(strsplit(as.character(d$Ind2),split="_"),"[[",1))
d=d[d$Method=="gle",]
m<-matrix(NA,9,9)
diag(m)<-0
dcol<-"Dij2"
m[lower.tri(m,diag=FALSE)]<- d[dcol][[1]]
dd.d<-as.dist(m)
dd.pops<-factor(c(rep("pop1",3),rep("pop2",3),rep("pop3",3)))
amv<-pegas::amova(dd.d ~ dd.pops,is.squared=TRUE)
sig2<-setNames(amv$varcomp$sigma2,rownames(amv$varcomp))
phi<-getPhi(sig2)
outfile<-"test_s9_d1_1K_pegas.amova.csv"
pegas_res<- data.frame("pegas",amv$tab$df[1],amv$tab$SSD[1],amv$tab$MSD[1],amv$tab$df[2],amv$tab$SSD[2],amv$tab$MSD[2],amv$tab$df[3],amv$tab$SSD[3],amv$tab$MSD[3],as.double(amv$varcoef[[1]]),amv$varcomp$sigma2[1],amv$varcomp$sigma2[2],phi[1])
colnames(pegas_res)<- c("Type","df1","SSD1","MSD1","df2","SSD2","MSD2","df3","SSD3","MSD3","varcoef","sigma2_1","sigma2_2","phi")
write.table(pegas_res,file=outfile, row.names=FALSE,  col.names=FALSE,  quote=FALSE, sep=',')
```