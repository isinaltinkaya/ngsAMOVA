# ngsAMOVA/test/reference


Generate the tab separated major minor allele file from the simulated bcf file from vcfgl, where the ancestral allelic state for all sites is A and derived state for all sites is C.

```sh
bcftools view -H sim_demes_v2-model1-1-rep13-d1.bcf |cut -f1,2,4,5|cut -d, -f1 > sim_demes_v2-model1-1-rep13-d1_ancder.tab
```

