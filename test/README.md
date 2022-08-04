## Tests



### Test missing sites

- test_mis, test_mis2

test_mis contains 14 sites: 9 missing for at least one individual and 5 that has data for both individuals.
test_mis2 contains only the 5 sites that has data for both individuals.


### Comparing AMOVA results with pegas

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



### Testing minInd filter

- nSites where at least 4 individuals have data for site

nMissingInds should be < (nInd-minInd+1)
minInd=4, nInd=9
threshold=9-4+1=6

```sh
$ for l in `bcftools query -f '[%DP]\n'  -s pop1_ind1,pop1_ind2 test/test_s9_d1_1K.vcf`;do echo $l|grep -o 0|wc -l;done|grep 0|wc -l 
387                                          
```

387 sites are nonmissing for both inds in this pair


```sh
$ for l in `bcftools query -f '[%DP]\n' test/test_s9_d1_1K.vcf`;do if [ `echo ${l}|cut -c-2|grep -o 0|wc -l` -eq 0 ];then echo $l|cut -c3-|grep -o 0|wc -l;fi;done|awk '$0>5'|wc -l
4
```

Skip 4 sites for this individual pair because even though it is nonmissing for both inds in pair these sites cannot pass the minInd threshold. Therefore we should end up with 387-4=383 sites.

```sh
$ for l in `bcftools query -f '[%DP]\n' test/test_s9_d1_1K.vcf`;do if [ `echo ${l}|cut -c-2|grep -o 0|wc -l` -eq 0 ];then echo $l|cut -c3-|grep -o 0|wc -l;fi;done|awk '$0<6'|wc -l
383

$ cat amovaput.sfs.csv |grep pop1_ind1,pop1_ind2|grep "^gt"|cut -d, -f4-12|tr ',' '\n'|datamash sum 1
383

$ cat amovaput.sfs.csv |grep pop1_ind1,pop1_ind2|grep "^gle"|cut -d, -f4-12|tr ',' '\n'|datamash sum 1
383.000001
$ cat amovaput.sfs.csv |grep "pop1_ind1,pop1_ind2"|cut -f14 -d,|uniq
383
```


- nSites where all individuals have data for site

minInd=0, nInd=9


```sh
$ for l in `bcftools query -f '[%DP]\n'  test/test_s9_d1_1K.vcf`;do echo $l|grep -o 0|wc -l;done | awk '$0==0'|wc -l
16
```


- Use the site if it has data for both individuals for that specific pair at site

minInd 2  (or param not used)


```sh
$ ./ngsAMOVA -in test/test_s9_d1_1K.vcf -m test/metadata.tsv -doDist 2 -doAMOVA 3 -isSim 1  -minInd  2
$ for l in `bcftools query -f '[%DP]\n'  -s pop1_ind1,pop1_ind2 test/test_s9_d1_1K.vcf`;do echo $l|grep -o 0|wc -l;done | awk '$0==0'|wc -l
387
$ cat amovaput.sfs.csv |grep "pop1_ind1,pop1_ind2"|cut -f14 -d,|uniq
387
$ cat amovaput.sfs.csv |grep pop1_ind1,pop1_ind2|grep "^gt"|cut -d, -f4-12|tr ',' '\n'|datamash sum 1
387
$ cat amovaput.sfs.csv |grep pop1_ind1,pop1_ind2|grep "^gt"|cut -d, -f14
387
$ cat amovaput.sfs.csv |grep pop1_ind1,pop1_ind2|grep "^gle"|cut -d, -f14
387
$ for l in `bcftools query -f '[%DP]\n'  -s pop2_ind1,pop3_ind2 test/test_s9_d1_1K.vcf`;do echo $l|grep -o 0|wc -l;done | awk '$0==0'|wc -l
359
$ cat amovaput.sfs.csv |grep "pop2_ind1,pop3_ind2"|cut -f14 -d,|uniq
359
```


