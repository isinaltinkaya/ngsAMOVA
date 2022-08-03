### Test


- test_mis, test_mis2

test_mis contains 14 sites: 9 missing for at least one individual and 5 that has data for both individuals.
test_mis2 contains only the 5 sites that has data for both individuals.

- test3.vcf

10 sites; (i1 missing sites: 2,6,8,10), (i2 missing sites: 8,10), (i3 missing sites: 4,8)
pairs; (i1,i2: 1,3,4,5,7,9) (i1,i3: 1,3,5,7,9) (i2,i3: 1,2,3,5,6,7,9)
onlyShared 1; (i1,i2: 1,3,5,7,9) (i1,i3: 1,3,5,7,9) (i2,i3: 1,3,5,7,9) 
onlyShared 0; (i1,i2: 1,3,4,5,7,9) (i1,i3: 1,3,5,7,9) (i2,i3: 1,2,3,5,6,7,9)
minInd 2 eq onlyShared 0; (i1,i2: 1,3,4,5,7,9) (i1,i3: 1,3,5,7,9) (i2,i3: 1,2,3,5,6,7,9)
minInd 3 eq onlyShared 1; (i1,i2: 1,3,5,7,9) (i1,i3: 1,3,5,7,9) (i2,i3: 1,3,5,7,9) 



- test_s9_d5_1K_pegas_amova_result.txt
```R
d<-read.csv("test_s9_d5_1K_test_minInd2_doAMOVA3.sfs.csv")
d$Pop1<-unlist(lapply(strsplit(d$Ind1,split="_"),"[[",1))
d$Pop2<-unlist(lapply(strsplit(d$Ind2,split="_"),"[[",1))
d=d[d$Method=="gt",]
m<-matrix(NA,9,9)
diag(m)<-0
dcol<-"Fij"
m[lower.tri(m,diag=FALSE)]<- d[dcol][[1]]
dd.d<-as.dist(m)
dd.pops<-factor(c(rep("pop1",3),rep("pop2",3),rep("pop3",3)))
(amv<-pegas::amova(dd.d ~ dd.pops,is.squared=FALSE))
pegas::write.pegas.amova(amv,file="test_s9_d5_1K_pegas_amova_result.txt")
```


