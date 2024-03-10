# doc/internal.MD

This document contains internal documentation for the ngsAMOVA package.

## Individual names

If distance matrix input, reorder names in dmat_t based on metad_t names.

<!-- If VCF input, set pars->names using hdr->samples, reorder metadata names based on pars->names (??) -->

```mermaid
graph TD;
    MTD --- MTDNAMES
    MTD[Metadata file]
    DATAIN[Data input file]
    MTDNAMES[metad_t \n strArray names];
    INDM[Distance matrix file]
    DISTNAMES[dmat_t \n strArray names];
    NAMES_TO_USE[strArray names_to_use];
    VCFNAMES[vcfd_t \n hdr->samples];

    JOIN1[join]
    JOIN2[join]
    MTDNAMES -.-> JOIN1
    MTDNAMES -.-> JOIN2

    DATAIN --- INDM --- DISTNAMES
    DATAIN --- INVCF --- VCFNAMES
    DISTNAMES --- JOIN1 ---> NAMES_TO_USE
    VCFNAMES --- JOIN2 --->  NAMES_TO_USE


```
