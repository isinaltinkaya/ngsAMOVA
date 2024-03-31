# doc/internal.MD

This document contains internal documentation for the ngsAMOVA package.

## Individual names

If distance matrix input, reorder names in dmat_t based on metadata_t names.

<!-- If VCF input, set pars->names using hdr->samples, reorder metadata names based on pars->names (??) -->

```mermaid
graph TD;

    DATAIN[Main data input]

    subgraph INPUT
    INDM[Distance matrix file]
    INVCF[VCF file]
    DATAIN --- INDM --- DISTNAMES
    DATAIN --- INVCF --- VCFNAMES
    end

    MTD --- MTDNAMES
    MTD[Metadata file]
    MTDNAMES[metadata_t \n strArray names];
    DISTNAMES[dmat_t \n strArray names];
    NAMES_TO_USE[strArray names_to_use];
    VCFNAMES[vcfdata_t \n hdr->samples];


    MTDNAMES -.- DISTNAMES
    MTDNAMES -.- VCFNAMES

    DISTNAMES ---> NAMES_TO_USE
    VCFNAMES --->  NAMES_TO_USE

```

```
if input is distance matrix
    if metadata input is provided
        filter and reorder dmat_t based on metadata_t names
        dmat->names = ptrto metadata->names
        pars->names = ptrto metadata->names
    else
        dmat->names = allocated and filled from dmat input
        pars->names = ptrto dmat->names

else if input is VCF

    if metadata input is provided
        exclude samples from VCF that are not in metadata (set bool vcfd->includeSamples[bcf_hdr_nsamples(vcfd->hdr)]), no need to reorder
        pars->names = ptrto metadata->names
    else
        pars->names = allocated and filled from vcfd->hdr->samples

    if doDist
        dmat->names = allocated and filled from pars->names

```




```mermaid
%%{
  init: {
    'theme': 'default',
    'themeVariables': {
      'primaryColor': '#BB2528',
      'primaryTextColor': '#fff',
      'primaryBorderColor': '#7C0000',
      'lineColor': '#F8B229',
      'secondaryColor': '#006100',
      'tertiaryColor': '#fff'
    }
  }
}%%
graph TD

    subgraph "Input Files"
        subgraph " "
            IN_VCF("VCF File\n[--in-vcf]")
            IN_VCF --- GLDATA_T("Genotype likelihoods")
            IN_VCF --- gtdata("Genotypes")
        end
        IN_DM("Distance Matrix File\n[--in-dm]")
        IN_MTD("Metadata File\n[--in-mtd]")
        IN_JGTM("Joint Genotype Matrix File\n[--in-jgtm]")
    end

    JGTMAT_T("Joint Genotype Matrix")
    DMAT_T("Pairwise Distance Matrix")
    doEM("[-doEM]\nPerform EM Optimization")
    doJGTM_GL("[-doJGTM]\nGet Joint Genotype Matrix")

    gl_jgtmat_em --> doJGTM_GL

    doJGTM_GT("[-doJGTM]\nGet Joint Genotype Matrix")


    subgraph "Processing Commands"
        doDist("[-doDist]\nCalculate pairwise distances")
    end


    subgraph "Analysis Commands"
        doAMOVA("[-doAMOVA]\nPerform AMOVA Analysis")
        doPhylo("[-doPhylo INT+]\nPerform Neighbor Joining")
        doDxy("[-doDxy]\nPerform Dxy Estimation")
    end

  
    gtdata --> doJGTM
    JGTMAT_T --> doDist
    doDist --> DMAT_T 
    DMAT_T --> doAMOVA
    DMAT_T --> doPhylo
    DMAT_T --> doDxy
    IN_MTD --> doDxy
    IN_MTD --> doAMOVA

    style vcf fill:#f0e68c,stroke:#333,stroke-width:4px
    style GLDATA_T fill:#f0e68c,stroke:#333,stroke-width:4px
    style gtdata fill:#f0e68c,stroke:#333,stroke-width:4px
    style dm fill:#f0e68c,stroke:#333,stroke-width:4px
    style EM fill:#add8e6,stroke:#333,stroke-width:4px
    style doJGTM fill:#add8e6,stroke:#333,stroke-width:4px
    style doDist fill:#add8e6,stroke:#333,stroke-width:4px
    style doAMOVA fill:#98fb98,stroke:#333,stroke-width:4px
    style doPhylo fill:#98fb98,stroke:#333,stroke-width:4px
    style doDxy fill:#98fb98,stroke:#333,stroke-width:4px
```




Diagram of available analysis commands and their dependencies:



```mermaid
%%{
  init: {
    'theme': 'default',
    'themeVariables': {
      'primaryColor': '#BB2528',
      'primaryTextColor': '#fff',
      'primaryBorderColor': '#7C0000',
      'lineColor': '#F8B229',
      'secondaryColor': '#006100',
      'tertiaryColor': '#fff'
    }
  }
}%%
graph TD
    subgraph "Input Files"
        subgraph " "
            vcf("VCF File\n[--in-vcf]")
            vcf --- gldata("Genotype likelihoods")
            vcf --- gtdata("Genotypes")
        end
        vcf("VCF File\n[--in-vcf]")
        dm("Distance Matrix File\n[--in-dm]")
        metadata("Metadata File")
    end

    subgraph "Processing Commands"
        EM("[-doEM]\nEM Optimization")
        doJGTM("[-doJGTM]\nJoint Genotype Matrix")
    end

    subgraph "Analysis Commands"
        doDist("[-doDist]\nPairwise distance matrix")
        doAMOVA("[-doAMOVA]\nAMOVA Analysis")
        doPhylo("[-doPhylo INT+]\nNeighbor Joining")
        doDxy("[-doDxy]\nDxy Estimation")
    end

    gldata --> EM
    gtdata --> doJGTM
    EM --> doJGTM
    doJGTM --> doDist
    doDist --> doAMOVA
    doDist --> doPhylo
    doDist --> doDxy
    dm --> doAMOVA
    dm --> doPhylo
    dm --> doDxy
    metadata --> doDxy
    metadata --> doAMOVA

    style vcf fill:#f0e68c,stroke:#333,stroke-width:4px
    style gldata fill:#f0e68c,stroke:#333,stroke-width:4px
    style gtdata fill:#f0e68c,stroke:#333,stroke-width:4px
    style dm fill:#f0e68c,stroke:#333,stroke-width:4px
    style EM fill:#add8e6,stroke:#333,stroke-width:4px
    style doJGTM fill:#add8e6,stroke:#333,stroke-width:4px
    style doDist fill:#98fb98,stroke:#333,stroke-width:4px
    style doAMOVA fill:#98fb98,stroke:#333,stroke-width:4px
    style doPhylo fill:#98fb98,stroke:#333,stroke-width:4px
    style doDxy fill:#98fb98,stroke:#333,stroke-width:4px
```



## VCF input

```mermaid
%%{
  init: {
    'theme': 'default',
    'themeVariables': {
      'primaryColor': '#BB2528',
      'primaryTextColor': '#fff',
      'primaryBorderColor': '#7C0000',
      'lineColor': '#F8B229',
      'secondaryColor': '#006100',
      'subgraph_title' : 'font-size:14px;',
      'tertiaryColor': '#fff'
    }
  }
}%%

graph TD
    subgraph vcfdata["Data Input"]
        vcf("VCF File\n[--in-vcf]")
        vcf --- gldata("Genotype likelihoods")
        vcf --- gtdata("Genotypes")
    end
    metadata("Metadata File")

    subgraph "Processing Commands"
        EM("[-doEM]\nEM Optimization")
        doJGTM("[-doJGTM]\nJoint Genotype Matrix")
    end

    subgraph "Analysis Commands"
        doPhylo("[-doPhylo]\nNeighbor Joining")
        doDist("[-doDist]\nPairwise distance matrix")
        doDxy("[-doDxy]\nDxy Estimation")
        doAMOVA("[-doAMOVA]\nAMOVA Analysis")
    end

    gldata --> EM
    gtdata --> doJGTM
    EM --> doJGTM
    doJGTM --> doDist
    doDist --> doAMOVA
    doDist --> doPhylo
    doDist --> doDxy
    metadata -.-> doDxy
    metadata -.-> doAMOVA

    style vcf fill:#f0e68c,stroke:#333,stroke-width:4px
    style gldata fill:#f0e68c,stroke:#333,stroke-width:4px
    style gtdata fill:#f0e68c,stroke:#333,stroke-width:4px
    style metadata fill:#f0e68c,stroke:#333,stroke-width:4px
    style EM fill:#add8e6,stroke:#333,stroke-width:4px
    style doJGTM fill:#add8e6,stroke:#333,stroke-width:4px
    style doDist fill:#98fb98,stroke:#333,stroke-width:4px
    style doAMOVA fill:#98fb98,stroke:#333,stroke-width:4px
    style doPhylo fill:#98fb98,stroke:#333,stroke-width:4px
    style doDxy fill:#98fb98,stroke:#333,stroke-width:4px
```

## Distance matrix input

```mermaid
%%{
  init: {
    'theme': 'default',
    'themeVariables': {
      'primaryColor': '#BB2528',
      'primaryTextColor': '#fff',
      'primaryBorderColor': '#7C0000',
      'lineColor': '#F8B229',
      'secondaryColor': '#006100',
      'subgraph_title' : 'font-size:14px;',
      'tertiaryColor': '#fff'
    }
  }
}%%
graph TD
    subgraph dmat["Data Input"]
        dm("Distance Matrix File\n[--in-dm]")
    end
    metadata("Metadata File")

    subgraph "Analysis Commands"
        doPhylo("[-doPhylo]\nNeighbor Joining")
        doDxy("[-doDxy]\nDxy Estimation")
        doAMOVA("[-doAMOVA]\nAMOVA Analysis")
    end

    dm --> doAMOVA
    dm --> doPhylo
    dm --> doDxy
    metadata -.-> doDxy
    metadata -.-> doAMOVA

    style dm fill:#f0e68c,stroke:#333,stroke-width:4px
    style metadata fill:#f0e68c,stroke:#333,stroke-width:4px
    style doAMOVA fill:#98fb98,stroke:#333,stroke-width:4px
    style doPhylo fill:#98fb98,stroke:#333,stroke-width:4px
    style doDxy fill:#98fb98,stroke:#333,stroke-width:4px
```


## lifetime

vcfd: until we have jgtmat
jgtmat: until we have dmat
dmat: until the end of all analyses







## Analysis commands

<!-- 

--in-vcf <STR> : input VCF file
                returns: vcfdata, gldata

--in-dm <STR>  : input distance matrix file
                returns: dmat

-doJGTM <INT>  : get pairwise joint genotype matrix
                returns: jgtmat

-doDist <INT>  : estimate pairwise distance matrix
                requires: jgtmat
                returns: dmat

-doAMOVA <INT> : perform AMOVA analysis
                requires: dmat, metadata, formula

-doEM <INT>    : perform EM optimization
                requires: gldata
                returns: jgtmat

-doDxy <INT>   : estimate Dxy
                requires: dmat

-doPhylo <INT> : do neighbor-joining tree
                requires: dmat

-doMajorMinor <INT> : get major and minor alleles for each site 
                requires: vcfdata, alleles input file (optional)
                returns: alleles
-->


