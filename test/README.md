## Tests


### Prepare sites.txt for angsd comparison
```sh
cat test_s9_d1_1K.vcf |grep -v "#"|cut -f1,2,4,5|cut -d',' -f1 > test_s9_d1_1K_sites.txt
```


```sh
ANGSD=/maps/projects/lundbeck/scratch/pfs488/AMOVA/method/dev_2210/angsd/angsd
OUT1=testout1
OUT2=testout2
INVCF=/maps/projects/lundbeck/scratch/pfs488/AMOVA/method/dev_2210/ngsAMOVA/test/test_s9_d1_1K.vcf
INSITES=/maps/projects/lundbeck/scratch/pfs488/AMOVA/method/dev_2210/ngsAMOVA/test/test_s9_d1_1K_sites.txt

${ANGSD} sites index ${INSITES}

P1I1=${INVCF%.vcf}.pop1_ind1.vcf
P1I2=${INVCF%.vcf}.pop1_ind2.vcf
bcftools view -s pop1_ind1 ${INVCF} > ${P1I1}
bcftools view -s pop1_ind2 ${INVCF} > ${P1I2}

${ANGSD} -vcf-gl ${P1I1} -doSaf 5 -doMajorMinor 3 -sites ${INSITES} -out ${P1I1}
${ANGSD} -vcf-gl ${P1I2} -doSaf 5 -doMajorMinor 3 -sites ${INSITES} -out ${P1I2}

${ANGSD%angsd}/misc/realSFS -m 0 -cores 1 ${P1I1}.saf.idx ${P1I2}.saf.idx > pop1_ind1-pop2_ind2.sfs.txt
```

```sh
./ngsAMOVA -in test/test_s9_d1_1K.vcf -m test/metadata.tsv -doDist 1 -doAMOVA 3 -isSim 1  -out testput -tole 1e-10 -maxIter 100
cat testput.sfs.csv |grep "^gle,pop1_ind1,pop1_ind2" | cut -d, -f4-12 > testput
```


Difference between angsd-realSFS result and ngsAMOVA result

```R
> range(scan("testput",sep=",")-scan("test/pop1_ind1-pop1_ind2.realSFS.txt"))                                                                 
Read 9 items                                                                                                                                  
Read 9 items                                                                                                                                  
[1] -1e-06  1e-06
```


Create "good file" for test comparison

```sh
cut -d, -f4-12 testput.sfs.csv > test_s9_d1_1K_god.sfs
```

### Test missing sites

- test_mis, test_mis2

test_mis contains 14 sites: 9 missing for at least one individual and 5 that has data for both individuals.
test_mis2 contains only the 5 sites that has data for both individuals.


### Comparing AMOVA results with pegas


```
cd test;
../ngsAMOVA -in test_s9_d5_1K.vcf -m metadata.tsv -doDist 1 -minInd 2 -doAMOVA 1 -isSim 1 -out test_s9_d5_1K_Dij_minInd2_doAmova1
Rscript test_pegas.R test_s9_d5_1K_Dij_minInd2_doAmova1.sfs.csv
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




