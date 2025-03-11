# ngsAMOVA/test/validation

- Validate nj newick output

```r
Rscript validate_nj_newick_output.R ../../amovaput.newick ../../amovaput.distance_matrix.csv ../metadata_with_header_2lvl.tsv
```


- Test block bootstrapping implementation

```sh
make dev && ./ngsAMOVA --in-vcf test2_d0.5_snp_1Ksites.vcf --bcf-src 1 -dodist 3  -doem 1  -domajorminor 1 -dojgtm 1 --print-dm 1  --block-size 50000 --print-blocks 1 --nbootstraps 1 --seed 2 -o out2  -v 10 --print-bootstrap 1

bcftools view -H test2_d0.5_snp_1Ksites.vcf  -t c1:1-50000 > test2_d0.5_snp_1Ksites_block0.vcf
bcftools view -H test2_d0.5_snp_1Ksites.vcf  -t c1:50001-100000 > test2_d0.5_snp_1Ksites_block1.vcf
bcftools view -H test2_d0.5_snp_1Ksites.vcf  -t c1:100001-150000 > test2_d0.5_snp_1Ksites_block2.vcf
bcftools view -h test2_d0.5_snp_1Ksites.vcf  > test2_d0.5_snp_1Ksites_rep0.vcf

bcftools view -h test2_d0.5_snp_1Ksites_rep0.vcf  > delme1rep0hdr
paste delme1 delme1rep0nohdr |awk '{ $3=$4=""; print $0 }' OFS='\t'| sed 's/\t\t/\t/g' > delme1rep0sortednohdr
cat delme1rep0hdr delme1rep0sortednohdr > test2_d0.5_snp_1Ksites_rep0posfixed.vcf

make dev && ./ngsAMOVA --in-vcf test2_d0.5_snp_1Ksites_rep0posfixed.vcf --bcf-src 1 -dodist 1  -doem 1  -domajorminor 1 -dojgtm 1 --print-dm 1  --seed 2 -o out2_rep0 -v 10
diff -s <(cat out2.distance_matrix.txt|sed  1,826d ) <(cat out2_rep0.distance_matrix.txt | sed 1,46d )
```
